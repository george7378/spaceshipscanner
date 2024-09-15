#define _USE_MATH_DEFINES
#include <windows.h>
#include <gdiplus.h>
#include <commctrl.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <urlmon.h>
#include <wininet.h>
#include <math.h>
#include <thread>

using namespace std;
using namespace Gdiplus;

#include "definitions.h"
#include "resource.h"

#pragma comment(lib, "Comctl32.lib")	//Common controls
#pragma comment(lib, "Gdiplus.lib")		//GDI+
#pragma comment(lib, "urlmon.lib")		//URLDownloadToFile
#pragma comment(lib, "wininet.lib")		//DeleteUrlCacheEntity

//Misc. functions
template <class T> wstring PrepareOutput(wstring in_prefix, T in_value, wstring in_unit, int in_precision)
{
	wstringstream stream;
	stream.setf(ios::fixed);
    stream << in_prefix << setprecision(in_precision) << in_value << in_unit;

    return stream.str();
}
string calendarDateToString(CalendarDate in_caldate)
{
	stringstream stream;
	stream << (unsigned(floor(in_caldate.day)) < 10 ? "0" : "") << unsigned(floor(in_caldate.day)) << "/" << (in_caldate.month < 10 ? "0" : "") << in_caldate.month << "/" << in_caldate.year << " ";
	
	double dayRemainder = 24*(in_caldate.day - floor(in_caldate.day));
	unsigned hour = unsigned(floor(dayRemainder));
	dayRemainder -= hour;
	dayRemainder *= 60;
	unsigned minute = unsigned(floor(dayRemainder));
	dayRemainder -= minute;
	dayRemainder *= 60;
	unsigned second = unsigned(floor(dayRemainder));
	stream << (hour < 10 ? "0" : "") << hour << ":" << (minute < 10 ? "0" : "") << minute << ":" << (second < 10 ? "0" : "") << second;

    return stream.str();
}
double wrapValue(double in_val, double in_max)
{
	double multiple = in_val/in_max;

	return (multiple - floor(multiple))*in_max;
}
double modifiedLongitude(double in_longitude)
{
	return in_longitude > 180 ? in_longitude - 360 : in_longitude;
}
double JulianDate(int in_day, int in_month, int in_year, double in_UT)
{
	//From 'Astronomical Algorithms'

    if (in_month <= 2){in_month += 12; in_year -= 1;}
	double A = floor(in_year/100.0);
	double B = 2 - A + floor(A/4);
    
    return floor(365.25*(in_year + 4716)) + floor(30.6001*(in_month + 1)) + in_day + B - 1524.5 + in_UT;
}
double LMSiderealTime(double in_JD, double in_longitude)
{
	//From 'Astronomical Algorithms'

	double T = (in_JD - 2451545)/36525;
	double GMST = 280.46061837 + 360.98564736629*(in_JD - 2451545) + 0.000387933*T*T - T*T*T/38710000;

	return wrapValue(GMST + in_longitude, 360);
}
CalendarDate CalendarFromJD(double in_JD)
{
	//From 'Astronomical Algorithms'

	CalendarDate out;
	ZeroMemory(&out, sizeof(CalendarDate));

	double Z = floor(in_JD + 0.5);
	double F = in_JD + 0.5 - Z;
	double A = Z;
	if (Z >= 2299161)
	{
		double Alpha = floor((Z - 1867216.25)/36524.25);
		A = Z + 1 + Alpha - floor(Alpha/4);
	}
	double B = A + 1524;
	double C = floor((B - 122.1)/365.25);
	double D = floor(365.25*C);
	double E = floor((B - D)/30.6001);

	out.day = B - D - floor(30.6001*E) + F;
	out.month = E < 14 ? unsigned(E - 1) : unsigned(E - 13);
	out.year = out.month > 2 ? int(C - 4716) : int(C - 4715);

	return out;
}
bool ActivateTLE(SingleTLE in_TLE)
{
	if (in_TLE.line2.length() < 69 || in_TLE.line3.length() < 69){return false;}

	string specialepoch = "20" + in_TLE.line2.substr(18, 2);																	//String containing epoch year
	string specialbstar = in_TLE.line2.substr(53, 1) + "0." + in_TLE.line2.substr(54, 5) + "e" + in_TLE.line2.substr(59, 2);	//String containing B-star drag term
	string specialmmd2 = in_TLE.line2.substr(33, 1) + "0" + in_TLE.line2.substr(34, 9);											//String containing first derivative of mean motion / 2
	string specialmmdd6 = in_TLE.line2.substr(44, 1) + "0." + in_TLE.line2.substr(45, 5) + "e" + in_TLE.line2.substr(50, 2);	//String containing second derivative of mean motion / 6

	double ejd = JulianDate(1, 1, stoi(specialepoch), 0) + stod(in_TLE.line2.substr(20, 12)) - 1;	//epochjd
	double bs = stod(specialbstar);																	//bstar
	double in = stod(in_TLE.line3.substr(8, 8));													//inclination
	double ra = stod(in_TLE.line3.substr(17, 8));													//raascendingnode
	double ec = stod("0." + in_TLE.line3.substr(26, 7));											//eccentricity
	double mm = stod(in_TLE.line3.substr(52, 11));													//meanmotion
	double mmd2 = stod(specialmmd2);																//mmDo2
	double mmdd6 = stod(specialmmdd6);																//mmDDo6
	double ap = stod(in_TLE.line3.substr(34, 8));													//argumentperigee
	double ma = stod(in_TLE.line3.substr(43, 8));													//meananomaly

	if (mm == 0 || ejd < 1){return false;}

	satName = in_TLE.line1;		//Satellite name
	epochjd = ejd;				//Julian date of epoch
	bstar = bs;					//B-star drag term
	inclination = in;			//Inclination
	raascendingnode = ra;		//Right ascension of ascending node
	eccentricity = ec;			//Eccentricity
	meanmotion = mm;			//Mean moton
	mmDo2 = mmd2;				//Mean motion dot over 2
	mmDDo6 = mmdd6;				//Mean motion dot dot over 6
	argumentperigee = ap;		//Argument of perigee
	meananomaly = ma;			//Mean anomaly
		
	inclinationRad = DEGRAD*inclination;
	raascendingnodeRad = DEGRAD*raascendingnode;
	argumentperigeeRad = DEGRAD*argumentperigee;
	meananomalyRad = DEGRAD*meananomaly;
	meanmotionRadPerMin = 2*M_PI*meanmotion/1440;
	mmDo2RadPerMinSq = 2*M_PI*mmDo2/1440/1440;
	mmDDo6RadPerMinCu = 2*M_PI*mmDDo6/1440/1440/1440;

	return true;
}
bool SaveSettings()
{
	ofstream prefs("prefs.txt", ios::out);

	if (prefs.is_open())
	{
		prefs << latitude << endl;
		prefs << longitude << endl;
		prefs << altitude;

		prefs.close();

		return true;
	}

	return false;
}
bool ReadSettings()
{
	ifstream prefs("prefs.txt", ios::in);

	if (prefs.is_open())
	{
		string lat, lon, alt;
		getline(prefs, lat);
		getline(prefs, lon);
		getline(prefs, alt);

		latitude = stod(lat);
		longitude = stod(lon);
		altitude = stod(alt);

		prefs.close();

		return true;
	}

	return false;
}

//Propagation functions
SunResults SolarPosition(double in_JD, double in_sidtimeobs)
{
	//From 'Astronomical Algorithms'

	SunResults out;
	ZeroMemory(&out, sizeof(SunResults));

	double T = (in_JD - 2451545)/36525;
	double L0 = 280.46645 + 36000.76983*T + 0.0003032*T*T;
	double M = 357.5291 + 35999.0503*T - 0.0001559*T*T - 0.00000048*T*T*T;
	double C = (1.9146 - 0.004817*T - 0.000014*T*T)*sin(DEGRAD*M) + (0.019993 - 0.000101*T)*sin(DEGRAD*2*M) + 0.00029*sin(DEGRAD*3*M);
	double L = L0 + C;
	double Eps = 23 + 26/60.0 + 21.448/3600 - (46.815*T + 0.00059*T*T - 0.001813*T*T*T)/3600;

	out.RA = RADDEG*atan2(cos(DEGRAD*Eps)*sin(DEGRAD*L), cos(DEGRAD*L));
	out.Dec = RADDEG*asin(sin(DEGRAD*Eps)*sin(DEGRAD*L));
	out.subSolarLat = out.Dec;
	out.subSolarLong = modifiedLongitude(wrapValue(out.RA - in_sidtimeobs + longitude, 360));

	return out;
}
PropagationResults SGP(double in_deltatMinutes)
{
	//From 'Space-Track report #3'

	PropagationResults out;
	ZeroMemory(&out, sizeof(PropagationResults));

	//Calculate constants based on TLE input
	double a1 = pow(ke/meanmotionRadPerMin, 2/3.0);
	double delta1 = 0.75*J2*ae*ae*(3*cos(inclinationRad)*cos(inclinationRad) - 1)/(a1*a1*pow(1 - eccentricity*eccentricity, 3/2.0));
	double a0 = a1*(1 - delta1/3 - delta1*delta1 - 134*delta1*delta1*delta1/81);
	double p0 = a0*(1 - eccentricity*eccentricity);
	double q0 = a0*(1 - eccentricity);
	double L0 = meananomalyRad + argumentperigeeRad + raascendingnodeRad;
	double minusDeriv = -1.5*J2*ae*ae*meanmotionRadPerMin*cos(inclinationRad)/(p0*p0);
	double omegaDeriv = 0.75*J2*ae*ae*meanmotionRadPerMin*(5*cos(inclinationRad)*cos(inclinationRad) - 1)/(p0*p0);

	double perigeeKM = RearthEq*(a0*(1 - eccentricity) - ae);
	double apogeeKM = RearthEq*(a0*(1 + eccentricity) - ae);

	//Secular effects of atmospheric drag
	double a = a0*pow(meanmotionRadPerMin/(meanmotionRadPerMin + 2*mmDo2RadPerMinSq*in_deltatMinutes
			   + 3*mmDDo6RadPerMinCu*in_deltatMinutes*in_deltatMinutes), 2/3.0);
	double e = a > q0 ? 1 - q0/a : 1e-6;
	double p = a*(1 - e*e);
	double minuss0 = raascendingnodeRad + minusDeriv*in_deltatMinutes;
	double omegas0 = argumentperigeeRad + omegaDeriv*in_deltatMinutes;
	double Ls = L0 + (meanmotionRadPerMin + omegaDeriv + minusDeriv)*in_deltatMinutes + mmDo2RadPerMinSq*in_deltatMinutes*in_deltatMinutes
				+ mmDDo6RadPerMinCu*in_deltatMinutes*in_deltatMinutes*in_deltatMinutes;

	//Long-period periodics
	double ayNSL = e*sin(omegas0) - 0.5*J3*ae*sin(inclinationRad)/(J2*p);
	double axNSL = e*cos(omegas0);
	double L = Ls - 0.25*J3*ae*axNSL*sin(inclinationRad)*((3 + 5*cos(inclinationRad))/(1 + cos(inclinationRad)))/(J2*p);
	double U = L - minuss0;

	double Eplusomegai = U;
	double Eplusomega;
	while (true)
	{
		Eplusomega = Eplusomegai + (U - ayNSL*cos(Eplusomegai) + axNSL*sin(Eplusomegai) - Eplusomegai)/(-ayNSL*sin(Eplusomegai)
					 - axNSL*cos(Eplusomegai) + 1);

		if (abs(Eplusomega - Eplusomegai) < 1e-10){break;}

		Eplusomegai = Eplusomega;
	}

	//Short-period preliminaries
	double ecosE = axNSL*cos(Eplusomega) + ayNSL*sin(Eplusomega);
	double esinE = axNSL*sin(Eplusomega) - ayNSL*cos(Eplusomega);
	double eLSq = axNSL*axNSL + ayNSL*ayNSL;
	double pL = a*(1 - eLSq);
	double r = a*(1 - ecosE);
	double rDot = ke*sqrt(a)*esinE/r;
	double rvDot = ke*sqrt(pL)/r;
	double sinu = a/r*(sin(Eplusomega) - ayNSL - axNSL*esinE/(1 + sqrt(1 - eLSq)));
	double cosu = a/r*(cos(Eplusomega) - axNSL + ayNSL*esinE/(1 + sqrt(1 - eLSq)));
	double u = atan2(sinu, cosu);

	//Short-period peturbations
	double rk = r + 0.25*J2*ae*ae*sin(inclinationRad)*sin(inclinationRad)*cos(2*u)/pL;
	double uk = u - 0.125*J2*ae*ae*(7*cos(inclinationRad)*cos(inclinationRad) - 1)*sin(2*u)/(pL*pL);
	double minusk = minuss0 + 0.75*J2*ae*ae*cos(inclinationRad)*sin(2*u)/(pL*pL);
	double ik = inclinationRad + 0.75*J2*ae*ae*sin(inclinationRad)*cos(inclinationRad)*cos(2*u)/(pL*pL);

	//Unit orientation vectors
	double Ux = -sin(minusk)*cos(ik)*sin(uk) + cos(minusk)*cos(uk);
	double Uy = cos(minusk)*cos(ik)*sin(uk) + sin(minusk)*cos(uk);
	double Uz = sin(ik)*sin(uk);
	
	double Vx = -sin(minusk)*cos(ik)*cos(uk) - cos(minusk)*sin(uk);
	double Vy = cos(minusk)*cos(ik)*cos(uk) - sin(minusk)*sin(uk);
	double Vz = sin(ik)*cos(uk);
	
	//Position and velocity
	double posX = rk*Ux;
	double posY = rk*Uy;
	double posZ = rk*Uz;

	double velX = rDot*Ux + rvDot*Vx;
	double velY = rDot*Uy + rvDot*Vy;
	double velZ = rDot*Uz + rvDot*Vz;
	
	//Store in the output structure
	out.posX = posX; out.posY = posY; out.posZ = posZ;
	out.velX = velX; out.velY = velY; out.velZ = velZ;
	out.perigeeKM = perigeeKM; out.apogeeKM = apogeeKM;

	return out;
}
PropagationResults SGP4(double in_deltatMinutes)
{
	//From 'Space-Track report #3'

	PropagationResults out;
	ZeroMemory(&out, sizeof(PropagationResults));

	//Calculate constants based on TLE input
	double a1 = pow(ke/meanmotionRadPerMin, 2/3.0);
	double delta1 = 1.5*k2*(3*cos(inclinationRad)*cos(inclinationRad) - 1)/(a1*a1*pow(1 - eccentricity*eccentricity, 3/2.0));
	double a0 = a1*(1 - delta1/3 - delta1*delta1 - 134*delta1*delta1*delta1/81);
	double delta0 = 1.5*k2*(3*cos(inclinationRad)*cos(inclinationRad) - 1)/(a0*a0*pow(1 - eccentricity*eccentricity, 3/2.0));
	double n0PrimePrime = meanmotionRadPerMin/(1 + delta0);
	double a0PrimePrime = a0/(1 - delta0);

	double perigeeKM = RearthEq*(a0PrimePrime*(1 - eccentricity)/ae - ae);
	double apogeeKM = RearthEq*(a0PrimePrime*(1 + eccentricity)/ae - ae);
	double s = 1.01222928;
	double q0ms24 = 1.88027916e-9;
	int TM = perigeeKM < 220 ? 0 : 1;

	if (perigeeKM >= 98 && perigeeKM <= 156)
	{
		double sStar = a0PrimePrime*(1 - eccentricity) - s + ae;
		q0ms24 = pow(pow(q0ms24, 1/4.0) + s - sStar, 4.0);
		s = sStar;
	}
	else if (perigeeKM < 98)
	{
		double sStar = 20/RearthEq + ae;
		q0ms24 = pow(pow(q0ms24, 1/4.0) + s - sStar, 4.0);
		s = sStar;
	}

	double theta = cos(inclinationRad);
	double xi = 1/(a0PrimePrime - s);
	double beta0 = sqrt(1 - eccentricity*eccentricity);
	double eta = a0PrimePrime*eccentricity*xi;
	
	double C2 = q0ms24*pow(xi, 4.0)*n0PrimePrime*pow(1 - eta*eta, -7/2.0)*(a0PrimePrime*(1 + 3*eta*eta/2 + 4*eccentricity*eta
				+ eccentricity*eta*eta*eta) + 1.5*k2*xi/(1 - eta*eta)*(-0.5 + 3*theta*theta/2)*(8 + 24*eta*eta + 3*pow(eta, 4.0)));
	double C1 = bstar*C2;
	double C3 = q0ms24*pow(xi, 5.0)*A30*n0PrimePrime*ae*sin(inclinationRad)/(k2*eccentricity);
	double C4 = 2*n0PrimePrime*q0ms24*pow(xi, 4.0)*a0PrimePrime*beta0*beta0*pow(1 - eta*eta, -7/2.0)*((2*eta*(1 + eccentricity*eta)
				+ eccentricity/2 + eta*eta*eta/2) - 2*k2*xi/(a0PrimePrime*(1 - eta*eta))*(3*(1 - 3*theta*theta)*(1 + 3*eta*eta/2
				- 2*eccentricity*eta - eccentricity*eta*eta*eta/2) + 0.75*(1 - theta*theta)*(2*eta*eta - eccentricity*eta
				- eccentricity*eta*eta*eta)*cos(2*argumentperigeeRad)));
	double C5 = 2*q0ms24*pow(xi, 4.0)*a0PrimePrime*beta0*beta0*pow(1 - eta*eta, -7/2.0)*(1 + 11*eta*(eta + eccentricity)/4
				+ eccentricity*eta*eta*eta);
	double D2 = 4*a0PrimePrime*xi*C1*C1;
	double D3 = 4*a0PrimePrime*xi*xi*(17*a0PrimePrime + s)*C1*C1*C1/3;
	double D4 = 2*a0PrimePrime*xi*xi*xi*(221*a0PrimePrime + 31*s)*pow(C1, 4.0)/3;

	//Secular effects of atmospheric drag
	double MDF = meananomalyRad + (1 + 3*k2*(-1 + 3*theta*theta)/(2*a0PrimePrime*a0PrimePrime*beta0*beta0*beta0)
				+ 3*k2*k2*(13 - 78*theta*theta + 137*pow(theta, 4.0))/(16*pow(a0PrimePrime, 4.0)*pow(beta0, 7.0)))*n0PrimePrime*in_deltatMinutes;
	double omegaDF = argumentperigeeRad + (-3*k2*(1 - 5*theta*theta)/(2*a0PrimePrime*a0PrimePrime*pow(beta0, 4.0)) + 3*k2*k2*(7
					 - 114*theta*theta + 395*pow(theta, 4.0))/(16*pow(a0PrimePrime, 4.0)*pow(beta0, 8.0)) + 5*k4*(3 - 36*theta*theta
					 + 49*pow(theta, 4.0))/(4*pow(a0PrimePrime, 4.0)*pow(beta0, 8.0)))*n0PrimePrime*in_deltatMinutes;
	double minusDF = raascendingnodeRad + (-3*k2*theta/(a0PrimePrime*a0PrimePrime*pow(beta0, 4.0)) + 3*k2*k2*(4*theta
					 - 19*theta*theta*theta)/(2*pow(a0PrimePrime, 4.0)*pow(beta0, 8.0)) + 5*k4*theta*(3
					 - 7*theta*theta)/(2*pow(a0PrimePrime, 4.0)*pow(beta0, 8.0)))*n0PrimePrime*in_deltatMinutes;
	double deltaomega = bstar*C3*cos(argumentperigeeRad)*in_deltatMinutes;
	double deltaM = -2/3.0*q0ms24*bstar*pow(xi, 4.0)*ae/(eccentricity*eta)*(pow(1 + eta*cos(MDF), 3.0) - pow(1 + eta*cos(meananomalyRad), 3.0));
	double Mp = MDF + TM*deltaomega + TM*deltaM;
	
	double omega = omegaDF - TM*deltaomega - TM*deltaM;
	double minus = minusDF - 21*n0PrimePrime*k2*theta*C1*in_deltatMinutes*in_deltatMinutes/(2*a0PrimePrime*a0PrimePrime*beta0*beta0);
	double e = eccentricity - bstar*C4*in_deltatMinutes - TM*bstar*C5*(sin(Mp) - sin(meananomalyRad));
	double a = a0PrimePrime*pow(1 - C1*in_deltatMinutes - TM*D2*in_deltatMinutes*in_deltatMinutes - TM*D3*in_deltatMinutes*in_deltatMinutes*in_deltatMinutes
			   - TM*D4*pow(in_deltatMinutes, 4.0), 2.0);
	double IL = Mp + omega + minus + n0PrimePrime*(3*C1*in_deltatMinutes*in_deltatMinutes/2
				+ TM*(D2 + 2*C1*C1)*in_deltatMinutes*in_deltatMinutes*in_deltatMinutes + TM*(3*D3 + 12*C1*D2 + 10*C1*C1*C1)*pow(in_deltatMinutes, 4.0)/4
				+ TM*(3*D4 + 12*C1*D3 + 6*D2*D2 + 30*C1*C1*D2 + 15*pow(C1, 4.0))*pow(in_deltatMinutes, 5.0)/5);
	double beta = sqrt(1 - e*e);
	double n = ke/pow(a, 3/2.0);

	//Long-period periodics
	double axN = e*cos(omega);
	double ILL = A30*sin(inclinationRad)*e*cos(omega)/(8*k2*a*beta*beta)*(3 + 5*theta)/(1 + theta);
	double ayNL = A30*sin(inclinationRad)/(4*k2*a*beta*beta);
	double ILT = IL + ILL;
	double ayN = e*sin(omega) + ayNL;
	double U = ILT - minus;
	
	double Eplusomegai = U;
	double Eplusomega;
	while (true)
	{
		Eplusomega = Eplusomegai + (U - ayN*cos(Eplusomegai) + axN*sin(Eplusomegai) - Eplusomegai)/(-ayN*sin(Eplusomegai)
					 - axN*cos(Eplusomegai) + 1);

		if (abs(Eplusomega - Eplusomegai) < 1e-10){break;}

		Eplusomegai = Eplusomega;
	}

	//Short-period preliminaries
	double ecosE = axN*cos(Eplusomega) + ayN*sin(Eplusomega);
	double esinE = axN*sin(Eplusomega) - ayN*cos(Eplusomega);
	double eL = sqrt(axN*axN + ayN*ayN);
	double pL = a*(1 - eL*eL);
	double r = a*(1 - ecosE);
	double rdot = ke*sqrt(a)*esinE/r;
	double rfdot = ke*sqrt(pL)/r;
	double cosu = a/r*(cos(Eplusomega) - axN + ayN*esinE/(1 + sqrt(1 - eL*eL)));
	double sinu = a/r*(sin(Eplusomega) - ayN - axN*esinE/(1 + sqrt(1 - eL*eL)));
	double u = atan2(sinu, cosu);
	
	double Changer = k2*(1 - theta*theta)*cos(2*u)/(2*pL);
	double Changeu = -k2*(7*theta*theta - 1)*sin(2*u)/(4*pL*pL);
	double Changeminus = 3*k2*theta*sin(2*u)/(2*pL*pL);
	double Changei = 3*k2*theta*sin(inclinationRad)*cos(2*u)/(2*pL*pL);
	double Changerdot = -k2*n*(1 - theta*theta)*sin(2*u)/pL;
	double Changerfdot = k2*n/pL*((1 - theta*theta)*cos(2*u) - 3*(1 - 3*theta*theta)/2);

	//Short-period periodics
	double rk = r*(1 - 3*k2*sqrt(1 - eL*eL)*(3*theta*theta - 1)/(2*pL*pL)) + Changer;
	double uk = u + Changeu;
	double minusk = minus + Changeminus;
	double ik = inclinationRad + Changei;
	double rdotk = rdot + Changerdot;
	double rfdotk = rfdot + Changerfdot;

	//Unit orientation vectors
	double Ux = -sin(minusk)*cos(ik)*sin(uk) + cos(minusk)*cos(uk);
	double Uy = cos(minusk)*cos(ik)*sin(uk) + sin(minusk)*cos(uk);
	double Uz = sin(ik)*sin(uk);
	
	double Vx = -sin(minusk)*cos(ik)*cos(uk) - cos(minusk)*sin(uk);
	double Vy = cos(minusk)*cos(ik)*cos(uk) - sin(minusk)*sin(uk);
	double Vz = sin(ik)*cos(uk);
	
	//Position and velocity
	double posX = rk*Ux;
	double posY = rk*Uy;
	double posZ = rk*Uz;

	double velX = rdotk*Ux + rfdotk*Vx;
	double velY = rdotk*Uy + rfdotk*Vy;
	double velZ = rdotk*Uz + rfdotk*Vz;

	//Store in the output structure
	out.posX = posX; out.posY = posY; out.posZ = posZ;
	out.velX = velX; out.velY = velY; out.velZ = velZ;
	out.perigeeKM = perigeeKM; out.apogeeKM = apogeeKM;

	return out;
}
PostpropResults PostPropFull(PropagationResults in_satresults, double in_UT, double in_sidtimeobs, double in_deltat)
{
	PostpropResults out;
	ZeroMemory(&out, sizeof(PostpropResults));
		
	//Sub-satellite latitude/longitude
	double R = sqrt(in_satresults.posX*in_satresults.posX + in_satresults.posY*in_satresults.posY);
	double satLongitude = modifiedLongitude(wrapValue(RADDEG*atan2(in_satresults.posY, in_satresults.posX) - in_sidtimeobs + longitude, 360));
	double satLatitudei = atan(in_satresults.posZ/R);
	double satLatitude, C;
	while (true)
	{
		C = 1/sqrt(1 - eSqEarth*sin(satLatitudei)*sin(satLatitudei));
		satLatitude = atan((in_satresults.posZ + C*eSqEarth*sin(satLatitudei))/R);

		if (abs(satLatitude - satLatitudei) < 1e-10){break;}

		satLatitudei = satLatitude;
	}
	satLatitude *= RADDEG;

	//Satellite altitude/velocity
	double sataltitude = RearthEq*(R/cos(DEGRAD*satLatitude) - C);
	double satvelocity = RearthEq*1440/86400.0*sqrt(in_satresults.velX*in_satresults.velX + in_satresults.velY*in_satresults.velY + in_satresults.velZ*in_satresults.velZ);

	//Satellite RA/Dec/range
	double Rg = pow(cos(DEGRAD*latitude)*cos(DEGRAD*latitude)/(RearthEq*RearthEq) + sin(DEGRAD*latitude)*sin(DEGRAD*latitude)/(RearthPo*RearthPo), -1/2.0) + altitude/1000;
	double ag = Rg*cos(DEGRAD*in_sidtimeobs)*cos(DEGRAD*latitude);
	double bg = Rg*sin(DEGRAD*in_sidtimeobs)*cos(DEGRAD*latitude);
	double cg = Rg*sin(DEGRAD*latitude);

	double xs = in_satresults.posX*RearthEq - ag;
	double ys = in_satresults.posY*RearthEq - bg;
	double zs = in_satresults.posZ*RearthEq - cg;

	double finalra = RADDEG*atan2(ys, xs);
	if (finalra < 0){finalra += 360;}
	double finalr = sqrt(xs*xs + ys*ys + zs*zs);
	double finaldec = RADDEG*asin(zs/finalr);

	//Satellite alt-az
	double hourangle = wrapValue(in_sidtimeobs - finalra, 360);
	double altitudeAngle = RADDEG*asin(sin(DEGRAD*latitude)*sin(DEGRAD*finaldec) + cos(DEGRAD*latitude)*cos(DEGRAD*finaldec)*cos(DEGRAD*hourangle));
	double azimuthAngle = 180 + RADDEG*atan2(sin(DEGRAD*hourangle), cos(DEGRAD*hourangle)*sin(DEGRAD*latitude) - tan(DEGRAD*finaldec)*cos(DEGRAD*latitude));
	if (azimuthAngle < 0){azimuthAngle += 360;}

	out.UT = in_UT; out.LST = in_sidtimeobs; out.DeltaT = in_deltat;
	out.RA = finalra; out.Dec = finaldec; out.Range = finalr; out.Alt = altitudeAngle; out.Az = azimuthAngle;
	out.Lat = satLatitude; out.Long = satLongitude; out.Altitude = sataltitude; out.Velocity = satvelocity;

	return out;
}
PostpropResults PostPropPassPredict(PropagationResults in_satresults, double in_sidtimeobs)
{
	PostpropResults out;
	ZeroMemory(&out, sizeof(PostpropResults));

	//Satellite RA/Dec/range
	double Rg = pow(cos(DEGRAD*latitude)*cos(DEGRAD*latitude)/(RearthEq*RearthEq) + sin(DEGRAD*latitude)*sin(DEGRAD*latitude)/(RearthPo*RearthPo), -1/2.0) + altitude/1000;
	double ag = Rg*cos(DEGRAD*in_sidtimeobs)*cos(DEGRAD*latitude);
	double bg = Rg*sin(DEGRAD*in_sidtimeobs)*cos(DEGRAD*latitude);
	double cg = Rg*sin(DEGRAD*latitude);

	double xs = in_satresults.posX*RearthEq - ag;
	double ys = in_satresults.posY*RearthEq - bg;
	double zs = in_satresults.posZ*RearthEq - cg;

	double finalra = RADDEG*atan2(ys, xs);
	if (finalra < 0){finalra += 360;}
	double finalr = sqrt(xs*xs + ys*ys + zs*zs);
	double finaldec = RADDEG*asin(zs/finalr);

	//Satellite alt-az
	double hourangle = wrapValue(in_sidtimeobs - finalra, 360);
	double altitudeAngle = RADDEG*asin(sin(DEGRAD*latitude)*sin(DEGRAD*finaldec) + cos(DEGRAD*latitude)*cos(DEGRAD*finaldec)*cos(DEGRAD*hourangle));
	double azimuthAngle = 180 + RADDEG*atan2(sin(DEGRAD*hourangle), cos(DEGRAD*hourangle)*sin(DEGRAD*latitude) - tan(DEGRAD*finaldec)*cos(DEGRAD*latitude));
	if (azimuthAngle < 0){azimuthAngle += 360;}

	out.Range = finalr; out.Alt = altitudeAngle; out.Az = azimuthAngle;

	return out;
}
PostpropResults PostPropGroundTrack(PropagationResults in_satresults, double in_sidtimeobs)
{
	PostpropResults out;
	ZeroMemory(&out, sizeof(PostpropResults));

	//Sub-satellite latitude/longitude
	double R = sqrt(in_satresults.posX*in_satresults.posX + in_satresults.posY*in_satresults.posY);
	double satLongitude = modifiedLongitude(wrapValue(RADDEG*atan2(in_satresults.posY, in_satresults.posX) - in_sidtimeobs + longitude, 360));
	double satLatitudei = atan(in_satresults.posZ/R);
	double satLatitude, C;
	while (true)
	{
		C = 1/sqrt(1 - eSqEarth*sin(satLatitudei)*sin(satLatitudei));
		satLatitude = atan((in_satresults.posZ + C*eSqEarth*sin(satLatitudei))/R);

		if (abs(satLatitude - satLatitudei) < 1e-10){break;}

		satLatitudei = satLatitude;
	}
	satLatitude *= RADDEG;

	out.Lat = satLatitude; out.Long = satLongitude;

	return out;
}

//Graphics functions
void update(HWND hwnd, HDC hdc)
{
	RECT rcClient;
	GetClientRect(hwnd, &rcClient);

	HDC hdcMem = CreateCompatibleDC(hdc);
	const int nMemDC = SaveDC(hdcMem);

	HBITMAP hBitmap = CreateCompatibleBitmap(hdc, rcClient.right - rcClient.left, rcClient.bottom - rcClient.top);
	SelectObject(hdcMem, hBitmap);

	Graphics graphics(hdcMem);
	
	Pen greenPen(Color(255, 0, 255, 0), 2);
	Pen grayPenDashed(Color(255, 150, 150, 150), 1);
	SolidBrush blackBrush(Color(255, 0, 0, 0));
	SolidBrush grayBrush(Color(255, 100, 100, 100));
	SolidBrush redBrush(Color(255, 255, 0, 0));
	SolidBrush yellowBrush(Color(255, 255, 255, 0));
	Font textfont(L"Ms Shell Dlg", 10);
	Font textfontLined(L"Ms Shell Dlg", 10, FontStyleUnderline);
	Image maptex(mapToUse.c_str());
	REAL equatorDashSpacing[2] = {1, 1};
	grayPenDashed.SetDashPattern(equatorDashSpacing, 2);

	//Draw the interface
	Rect maprect(int(rcClient.left + 10), int(rcClient.top + 10), int(rcClient.right - rcClient.left - 20), int(0.6*(rcClient.bottom - rcClient.top)));
	graphics.FillRectangle(&blackBrush, rcClient.left, rcClient.top, rcClient.right - rcClient.left, rcClient.bottom - rcClient.top);
	graphics.SetInterpolationMode(InterpolationModeBilinear);
	graphics.DrawImage(&maptex, maprect);
	graphics.DrawLine(&grayPenDashed, int(maprect.GetLeft()), int(maprect.GetTop() + 0.5*maprect.Height), int(maprect.GetRight()), int(maprect.GetTop() + 0.5*maprect.Height));
	graphics.DrawRectangle(&greenPen, int(maprect.GetLeft() - 1), int(maprect.GetTop() - 1), int(maprect.Width + 2), int(maprect.Height + 2));
	graphics.FillRectangle(&grayBrush, int(rcClient.left), int(rcClient.bottom - 20), int(rcClient.right - rcClient.left), 20);

	//Draw the observer location
	float observerMapPos[2] = {float(maprect.GetLeft() + maprect.Width/2.0 + maprect.Width/2.0*longitude/180 - 3), float(maprect.GetTop() + maprect.Height/2.0 - maprect.Height/2.0*latitude/90 - 3)};
	graphics.FillRectangle(&redBrush, observerMapPos[0], observerMapPos[1], 6.0f, 6.0f);

	//Draw the headings
	graphics.DrawString(L"Topocentric equatorial:", -1, &textfontLined, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + 10)), &yellowBrush);
	graphics.DrawString(L"Topocentric horizontal:", -1, &textfontLined, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + 3*textfont.GetHeight(&graphics) + 50)), &yellowBrush);
	graphics.DrawString(L"Orbital data:", -1, &textfontLined, PointF(float(rcClient.left + 410), float(maprect.GetBottom() + 10)), &yellowBrush);

	//Draw warning message
	if (epochjd < 1)
	{
		SolidBrush grayBrushAlpha(Color(180, 100, 100, 100));
		Font textfontLargeItalic(L"Ms Shell Dlg", 16, FontStyleItalic);

		graphics.FillRectangle(&grayBrushAlpha, int(rcClient.left + 80), int(rcClient.top + 70), 440, 60);
		graphics.DrawString(L"Select a satellite using the button below", -1, &textfontLargeItalic, PointF(float(rcClient.left + 100), float(rcClient.top + 86)), &yellowBrush);
	}

	//Draw the results
	if (displaySatResults)
	{
		displaySatResults = false;

		Pen whitePen(Color(255, 255, 255, 255), 1);
		Pen yellowPen(Color(255, 255, 255, 0), 1);
		Pen redPen(Color(255, 255, 0, 0), 1);
		SolidBrush whiteBrush(Color(255, 255, 255, 255));
		SolidBrush greenBrush(Color(255, 0, 255, 0));
		Image sunicon(L"sun.png");

		//Draw the solar terminator/Sun
		graphics.SetSmoothingMode(SmoothingModeAntiAlias);
		for (unsigned i = 1; i < terminatorPoints.size(); i++)
		{	
			PointF termPt1MapPos(float(maprect.GetLeft() + maprect.Width/2.0 + maprect.Width/2.0*terminatorPoints[i].Long/180), float(maprect.GetTop() + maprect.Height/2.0 - maprect.Height/2.0*terminatorPoints[i].Lat/90));
			PointF termPt2MapPos(float(maprect.GetLeft() + maprect.Width/2.0 + maprect.Width/2.0*terminatorPoints[i - 1].Long/180), float(maprect.GetTop() + maprect.Height/2.0 - maprect.Height/2.0*terminatorPoints[i - 1].Lat/90));
			int stepModulus = i % 4;
			graphics.DrawLine((stepModulus == 0 || stepModulus == 1) ? &yellowPen : &redPen, termPt1MapPos, termPt2MapPos);
		}
		graphics.SetSmoothingMode(SmoothingModeDefault);
		RectF sunIconRect(float(maprect.GetLeft() + maprect.Width/2.0 + maprect.Width/2.0*sunresult.subSolarLong/180 - 16), float(maprect.GetTop() + maprect.Height/2.0 - maprect.Height/2.0*sunresult.subSolarLat/90 - 16), 32, 32);
		graphics.DrawImage(&sunicon, sunIconRect);

		//Draw the ground track
		graphics.SetSmoothingMode(SmoothingModeAntiAlias);
		for (unsigned i = 1; i < groundTrackPoints.size(); i++)
		{
			PointF groundtrackPt1MapPos(float(maprect.GetLeft() + maprect.Width/2.0 + maprect.Width/2.0*groundTrackPoints[i].Long/180), float(maprect.GetTop() + maprect.Height/2.0 - maprect.Height/2.0*groundTrackPoints[i].Lat/90));
			PointF groundtrackPt2MapPos(float(maprect.GetLeft() + maprect.Width/2.0 + maprect.Width/2.0*groundTrackPoints[i - 1].Long/180), float(maprect.GetTop() + maprect.Height/2.0 - maprect.Height/2.0*groundTrackPoints[i - 1].Lat/90));
			if (fabsf(groundtrackPt1MapPos.X - groundtrackPt2MapPos.X) > maprect.Width/2.0){continue;}
			int stepModulus = i % 4;
			graphics.DrawLine((stepModulus == 0 || stepModulus == 1) ? &whitePen : &redPen, groundtrackPt1MapPos, groundtrackPt2MapPos);
		}
		graphics.SetSmoothingMode(SmoothingModeDefault);

		//Topocentric data
		graphics.DrawString(PrepareOutput(L"RA: ", postresult.RA, L" °", 2).c_str(), -1, &textfont, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + 20)), &yellowBrush);
		graphics.DrawString(PrepareOutput(L"Dec: ", postresult.Dec, L" °", 2).c_str(), -1, &textfont, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + textfont.GetHeight(&graphics) + 25)), &yellowBrush);
		graphics.DrawString(PrepareOutput(L"Range: ", postresult.Range, L" km", 1).c_str(), -1, &textfont, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + 2*textfont.GetHeight(&graphics) + 30)), &yellowBrush);

		graphics.DrawString(PrepareOutput(L"Altitude: ", postresult.Alt, L" °", 1).c_str(), -1, &textfont, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + 2*textfontLined.GetHeight(&graphics) + 3*textfont.GetHeight(&graphics) + 60)), &yellowBrush);
		graphics.DrawString(PrepareOutput(L"Azimuth: ", postresult.Az, L" °", 1).c_str(), -1, &textfont, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + 2*textfontLined.GetHeight(&graphics) + 4*textfont.GetHeight(&graphics) + 65)), &yellowBrush);
		if (postresult.Alt <= 0)
		{
			graphics.DrawString(L"Below horizon", -1, &textfont, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + 2*textfontLined.GetHeight(&graphics) + 5*textfont.GetHeight(&graphics) + 70)), &redBrush);
		}
		else
		{
			graphics.DrawString(L"Above horizon", -1, &textfont, PointF(float(rcClient.left + 220), float(maprect.GetBottom() + 2*textfontLined.GetHeight(&graphics) + 5*textfont.GetHeight(&graphics) + 70)), &greenBrush);
		}

		//Remaining data
		graphics.DrawString(PrepareOutput(L"Perigee: ", result.perigeeKM, L" km", 3).c_str(), -1, &textfont, PointF(float(rcClient.left + 410), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + 20)), &yellowBrush);
		graphics.DrawString(PrepareOutput(L"Apogee: ", result.apogeeKM, L" km", 3).c_str(), -1, &textfont, PointF(float(rcClient.left + 410), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + textfont.GetHeight(&graphics) + 25)), &yellowBrush);
		graphics.DrawString(PrepareOutput(L"Altitude: ", postresult.Altitude, L" km", 3).c_str(), -1, &textfont, PointF(float(rcClient.left + 410), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + 2*textfont.GetHeight(&graphics) + 30)), &yellowBrush);
		graphics.DrawString(PrepareOutput(L"Velocity: ", postresult.Velocity, L" km/s", 3).c_str(), -1, &textfont, PointF(float(rcClient.left + 410), float(maprect.GetBottom() + textfontLined.GetHeight(&graphics) + 3*textfont.GetHeight(&graphics) + 35)), &yellowBrush);
		
		graphics.DrawString(PrepareOutput(L"Sub-sat Latitude: ", postresult.Lat, L" °", 2).c_str(), -1, &textfont, PointF(float(rcClient.left + 410), float(maprect.GetBottom() + 2*textfontLined.GetHeight(&graphics) + 3*textfont.GetHeight(&graphics) + 60)), &whiteBrush);
		graphics.DrawString(PrepareOutput(L"Sub-sat Longitude: ", postresult.Long, L" °", 2).c_str(), -1, &textfont, PointF(float(rcClient.left + 410), float(maprect.GetBottom() + 2*textfontLined.GetHeight(&graphics) + 4*textfont.GetHeight(&graphics) + 65)), &whiteBrush);

		//Timing data
		graphics.DrawString(PrepareOutput(L"UTC: ", 24*postresult.UT, L" hr", 4).c_str(), -1, &textfont, PointF(float(rcClient.left + 10), float(rcClient.bottom - textfont.GetHeight(&graphics)/2.0 - 10)), &whiteBrush);
		graphics.DrawString(PrepareOutput(L"LST: ", 24*postresult.LST/360, L" hr", 4).c_str(), -1, &textfont, PointF(float(rcClient.left + 220), float(rcClient.bottom - textfont.GetHeight(&graphics)/2.0 - 10)), &whiteBrush);
		graphics.DrawString(PrepareOutput(L"TLE age: ", postresult.DeltaT, L" days", 2).c_str(), -1, &textfont, PointF(float(rcClient.left + 410), float(rcClient.bottom - textfont.GetHeight(&graphics)/2.0 - 10)), &whiteBrush);

		//Satellite blip on map
		float satBlipMapPos[2] = {float(maprect.GetLeft() + maprect.Width/2.0 + maprect.Width/2.0*postresult.Long/180 - 3), float(maprect.GetTop() + maprect.Height/2.0 - maprect.Height/2.0*postresult.Lat/90 - 3)};
		graphics.FillRectangle(&greenBrush, satBlipMapPos[0], satBlipMapPos[1], 6.0f, 6.0f);
		wstring satNameW(satName.begin(), satName.end());
		RectF satNameBoundRect;
		graphics.MeasureString(satNameW.c_str(), -1, &textfont, PointF(0, 0), StringFormat::GenericTypographic(), &satNameBoundRect);
		PointF satNameOrigin;
		if (postresult.Long < 0)
		{
			if (postresult.Lat > 0)
			{
				satNameOrigin = PointF(float(satBlipMapPos[0] + 8), float(satBlipMapPos[1])); //Top left
			}
			else
			{
				satNameOrigin = PointF(float(satBlipMapPos[0] + 8), float(satBlipMapPos[1] - satNameBoundRect.Height + 3)); //Bottom left
			}
		}
		else
		{
			if (postresult.Lat > 0)
			{
				satNameOrigin = PointF(float(satBlipMapPos[0] - satNameBoundRect.Width - 2), float(satBlipMapPos[1])); //Top right
			}
			else
			{
				satNameOrigin = PointF(float(satBlipMapPos[0] - satNameBoundRect.Width - 2), float(satBlipMapPos[1] - satNameBoundRect.Height + 3)); //Bottom right
			}
		}
		graphics.DrawString(satNameW.c_str(), -1, &textfont, satNameOrigin, StringFormat::GenericTypographic(), &greenBrush);
	}

	RECT rcClip;
	GetClipBox(hdc, &rcClip);
	BitBlt(hdc, rcClip.left, rcClip.top, rcClip.right - rcClip.left, rcClip.bottom - rcClip.top, hdcMem, rcClip.left, rcClip.top, SRCCOPY);

	RestoreDC(hdcMem, nMemDC);
	DeleteObject(hBitmap);
	DeleteDC(hdcMem);
}

//Windows functions
void on_satman_update(HWND hwnd)
{
	EnableWindow(hwnd, FALSE);

	bool failures = false;
	HWND edtStatusBox = GetDlgItem(hwnd, IDC_SATMAN_UPDATESTATUS);
	string curStatus = "Starting TLE update process:\r\n---";
	SetDlgItemText(hwnd, IDC_SATMAN_UPDATESTATUS, curStatus.c_str());

	for (unsigned i = 0; i < 8; i++)
	{
		curStatus += "\r\nAttempting download of " + tleDataFiles[i];
		SetDlgItemText(hwnd, IDC_SATMAN_UPDATESTATUS, curStatus.c_str());
		SendMessage(edtStatusBox, EM_LINESCROLL, 0, SendMessage(edtStatusBox, EM_GETLINECOUNT, 0, 0));

		string curURL = "http://celestrak.com/NORAD/elements/" + tleDataFiles[i];
		DeleteUrlCacheEntry(curURL.c_str());
		if (SUCCEEDED(URLDownloadToFile(NULL, curURL.c_str(), string("satdata/" + tleDataFiles[i]).c_str(), 0, NULL)))
		{
			curStatus += "\r\n- Started download of " + tleDataFiles[i];
		}
		else
		{
			failures = true;
			curStatus += "\r\n- FAILED to download " + tleDataFiles[i];
		}

		SetDlgItemText(hwnd, IDC_SATMAN_UPDATESTATUS, curStatus.c_str());
		SendMessage(edtStatusBox, EM_LINESCROLL, 0, SendMessage(edtStatusBox, EM_GETLINECOUNT, 0, 0));
	}

	curStatus += "\r\n---\r\nFinished TLE update process:\r\n";
	curStatus += failures ? "SOME UPDATES FAILED" : "ALL UPDATES SUCCESSFUL";
	SetDlgItemText(hwnd, IDC_SATMAN_UPDATESTATUS, curStatus.c_str());
	SendMessage(edtStatusBox, EM_LINESCROLL, 0, SendMessage(edtStatusBox, EM_GETLINECOUNT, 0, 0));

	EnableWindow(hwnd, TRUE);
}
void on_passes_calculate(HWND hwnd, double jd_begin,  double jd_end)
{
	EnableWindow(hwnd, FALSE);

	HWND pbrPassCalcProgress = GetDlgItem(hwnd, IDC_PASSES_PROGRESS);
	vector <PassPrediction> predictionResults;
	predictionResults.clear();
	SendMessage(lvwPassResults, LVM_DELETEALLITEMS, 0, 0);

	for (double coarseStep = jd_begin; coarseStep <= jd_end; coarseStep += 1/1440.0)
	{
		SendMessage(pbrPassCalcProgress, PBM_SETPOS, int(ceil(100*(coarseStep - jd_begin)/(jd_end - jd_begin))), 0);

		bool insideExistingPass = false;
		for (unsigned i = 0; i < predictionResults.size(); i++)
		{
			if (coarseStep >= predictionResults[i].riseTimeJD && coarseStep <= predictionResults[i].setTimeJD)
			{
				insideExistingPass = true;
				break;
			}
		}
		if (insideExistingPass){continue;}

		//Coarse propagation
		double coarseSiderealTimeObserver = LMSiderealTime(coarseStep, longitude);
		double coarseDeltaT = coarseStep - epochjd;
		PropagationResults coarsePropResult = SGP4(1440*coarseDeltaT);
		PostpropResults coarsePostPropResult = PostPropPassPredict(coarsePropResult, coarseSiderealTimeObserver);

		if (coarsePostPropResult.Alt > 0)
		{
			PassPrediction currentPass;
			ZeroMemory(&currentPass, sizeof(PassPrediction));

			double fineStep = coarseStep;
			PropagationResults finePropResult = coarsePropResult;
			PostpropResults finePostPropResult = coarsePostPropResult;
			currentPass.peakAltitude = finePostPropResult.Alt;

			//Fine backwards propagation
			while (finePostPropResult.Alt > 0 && fineStep >= jd_begin)
			{
				fineStep -= 1/86400.0;
				double fineSiderealTimeObserver = LMSiderealTime(fineStep, longitude);
				double fineDeltaT = fineStep - epochjd;
				finePropResult = SGP4(1440*fineDeltaT);
				finePostPropResult = PostPropPassPredict(finePropResult, fineSiderealTimeObserver);
				
				if (finePostPropResult.Alt > currentPass.peakAltitude)
				{
					currentPass.peakAltitude = finePostPropResult.Alt;
					currentPass.peakTimeJD = fineStep;
				}
			}
			currentPass.riseTimeJD = fineStep;
			currentPass.riseAzimuth = finePostPropResult.Az;

			//Fine forwards propagation
			fineStep = coarseStep;
			finePropResult = coarsePropResult;
			finePostPropResult = coarsePostPropResult;
			while (finePostPropResult.Alt > 0 && fineStep <= jd_end)
			{
				fineStep += 1/86400.0;
				double fineSiderealTimeObserver = LMSiderealTime(fineStep, longitude);
				double fineDeltaT = fineStep - epochjd;
				finePropResult = SGP4(1440*fineDeltaT);
				finePostPropResult = PostPropPassPredict(finePropResult, fineSiderealTimeObserver);
				
				if (finePostPropResult.Alt > currentPass.peakAltitude)
				{
					currentPass.peakAltitude = finePostPropResult.Alt;
					currentPass.peakTimeJD = fineStep;
				}
			}
			currentPass.setTimeJD = fineStep;
			currentPass.setAzimuth = finePostPropResult.Az;

			if (currentPass.peakAltitude >= 10)
			{
				predictionResults.push_back(currentPass);
			}
		}
	}

	if (predictionResults.size() == 0)
	{
		MessageBox(hwnd, "No satellite passes detected for chosen period.", "Information", MB_OK | MB_ICONINFORMATION);
	}
	else
	{
		char listViewTemp[256];
		for (unsigned j = 0; j < predictionResults.size(); j++)
		{
			LVITEM lvEntry;
			ZeroMemory(&lvEntry, sizeof(LVITEM));
			lvEntry.mask = LVIF_TEXT;
			lvEntry.cchTextMax = 256;
			lvEntry.iItem = j;
			strcpy_s(listViewTemp, calendarDateToString(CalendarFromJD(predictionResults[j].riseTimeJD)).c_str());
			lvEntry.pszText = listViewTemp;
			SendMessage(lvwPassResults, LVM_INSERTITEM, 0, (LPARAM)&lvEntry);
			sprintf_s(listViewTemp, "%d", int(floor(predictionResults[j].riseAzimuth + 0.5)));
			lvEntry.pszText = listViewTemp;
			lvEntry.iSubItem = 1;
			SendMessage(lvwPassResults, LVM_SETITEM, 0, (LPARAM)&lvEntry);
			strcpy_s(listViewTemp, calendarDateToString(CalendarFromJD(predictionResults[j].peakTimeJD)).c_str());
			lvEntry.pszText = listViewTemp;
			lvEntry.iSubItem = 2;
			SendMessage(lvwPassResults, LVM_SETITEM, 0, (LPARAM)&lvEntry);
			sprintf_s(listViewTemp, "%d", int(floor(predictionResults[j].peakAltitude + 0.5)));
			lvEntry.pszText = listViewTemp;
			lvEntry.iSubItem = 3;
			SendMessage(lvwPassResults, LVM_SETITEM, 0, (LPARAM)&lvEntry);
			strcpy_s(listViewTemp, calendarDateToString(CalendarFromJD(predictionResults[j].setTimeJD)).c_str());
			lvEntry.pszText = listViewTemp;
			lvEntry.iSubItem = 4;
			SendMessage(lvwPassResults, LVM_SETITEM, 0, (LPARAM)&lvEntry);
			sprintf_s(listViewTemp, "%d", int(floor(predictionResults[j].setAzimuth + 0.5)));
			lvEntry.pszText = listViewTemp;
			lvEntry.iSubItem = 5;
			SendMessage(lvwPassResults, LVM_SETITEM, 0, (LPARAM)&lvEntry);
		}
	}

	EnableWindow(hwnd, TRUE);
}
BOOL CALLBACK PrefDlgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch(msg)
	{
		case WM_INITDIALOG:
			{
				char strOutput[100];
				sprintf_s(strOutput, "%0.4f", latitude);
				SetDlgItemText(hwnd, IDC_PREF_LAT, strOutput);
				sprintf_s(strOutput, "%0.4f", longitude);
				SetDlgItemText(hwnd, IDC_PREF_LON, strOutput);
				sprintf_s(strOutput, "%0.0f", altitude);
				SetDlgItemText(hwnd, IDC_PREF_ALTASL, strOutput);
			}
			break;

		case WM_COMMAND:
			switch(LOWORD(wParam))
			{
				case IDC_PREF_UPDATE:
					{
						int lenLat = GetWindowTextLength(GetDlgItem(hwnd, IDC_PREF_LAT));
						int lenLon = GetWindowTextLength(GetDlgItem(hwnd, IDC_PREF_LON));
						int lenAlt = GetWindowTextLength(GetDlgItem(hwnd, IDC_PREF_ALTASL));
						char* bufLat = (char*)GlobalAlloc(GPTR, lenLat + 1);
						char* bufLon = (char*)GlobalAlloc(GPTR, lenLon + 1);
						char* bufAlt = (char*)GlobalAlloc(GPTR, lenAlt + 1);
						GetDlgItemText(hwnd, IDC_PREF_LAT, bufLat, lenLat + 1);
						GetDlgItemText(hwnd, IDC_PREF_LON, bufLon, lenLon + 1);
						GetDlgItemText(hwnd, IDC_PREF_ALTASL, bufAlt, lenAlt + 1);

						double latitudeTemp = stod(bufLat);
						double longitudeTemp = stod(bufLon);
						if (latitudeTemp > 90 || latitudeTemp < -90)
						{
							MessageBox(hwnd, "Latitude must be between +/- 90°", "Error", MB_OK | MB_ICONEXCLAMATION);
							break;
						}
						else if (longitudeTemp > 180 || longitudeTemp < -180)
						{
							MessageBox(hwnd, "Longitude must be between +/- 180°", "Error", MB_OK | MB_ICONEXCLAMATION);
							break;
						}

						latitude = latitudeTemp;
						longitude = longitudeTemp;
						altitude = stod(bufAlt);
					
						GlobalFree(bufLat); GlobalFree(bufLon); GlobalFree(bufAlt);

						if (SaveSettings())
						{
							MessageBox(hwnd, "Preferences successfully updated.", "Information", MB_OK | MB_ICONINFORMATION);
							EndDialog(hwnd, 0);
						}
						else
						{
							MessageBox(hwnd, "Unable to store preferences.", "Error", MB_OK | MB_ICONEXCLAMATION);
						}
					}
					break;
			}
			break;

		//case WM_CTLCOLORSTATIC:
		//	SetTextColor((HDC)wParam, RGB(220, 235, 120));
		//	SetBkMode((HDC)wParam, TRANSPARENT);
		//	return (INT_PTR)CreateSolidBrush(RGB(35, 125, 120));

		//case WM_CTLCOLORDLG:
		//	return (INT_PTR)CreateSolidBrush(RGB(35, 125, 120));

		case WM_CLOSE:
			EndDialog(hwnd, 0);
		    break;

		default:
			return FALSE;
	}

	return TRUE;
}
BOOL CALLBACK SatManDlgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch(msg)
	{
		case WM_INITDIALOG:
			{
				lboActiveSat = GetDlgItem(hwnd, IDC_SATMAN_ACTIVESAT);

				cboTLEFile = GetDlgItem(hwnd, IDC_SATMAN_TLEFILE);
				for (unsigned i = 0; i < 8; i++){SendMessage(cboTLEFile, CB_ADDSTRING, 0, (LPARAM)tleDataFiles[i].c_str());}
				SendMessage(cboTLEFile, CB_SETCURSEL, -1, 0);
				SetDlgItemText(hwnd, IDC_SATMAN_TLEFILE, string("Select a TLE source...").c_str());
			}
			break;

		case WM_COMMAND:
			switch(LOWORD(wParam))
			{
				case IDC_SATMAN_TLEFILE:
					switch(HIWORD(wParam))
					{
						case CBN_SELCHANGE:
							unsigned chosenTLESource = SendMessage(cboTLEFile, CB_GETCURSEL, 0, 0);
							ifstream tlesource("satdata/" + tleDataFiles[chosenTLESource], ios::in);
							SendMessage(lboActiveSat, LB_RESETCONTENT, 0, 0);
							tlesFromDataFile.clear();

							if (!tlesource.is_open())
							{
								MessageBox(hwnd, "Unable to load selected TLE source. Please update TLEs.", "Error", MB_OK | MB_ICONEXCLAMATION);
							}
							else
							{
								while (!tlesource.eof())
								{
									SingleTLE thisTLE;
									ZeroMemory(&thisTLE, sizeof(SingleTLE));

									getline(tlesource, thisTLE.line1);
									getline(tlesource, thisTLE.line2);
									getline(tlesource, thisTLE.line3);

									if (thisTLE.line1.size() > 0 && thisTLE.line2.size() > 0 && thisTLE.line3.size() > 0){tlesFromDataFile.push_back(thisTLE);}
								}
							}

							tlesource.close();
							
							for (unsigned i = 0; i < tlesFromDataFile.size(); i++) 
							{
								SendMessage(lboActiveSat, LB_ADDSTRING, 0, (LPARAM)tlesFromDataFile[i].line1.c_str());
							}
							SendMessage(lboActiveSat, LB_SETCURSEL, -1, 0);
							break;
					}
					break;

				case IDC_SATMAN_UPDATE:
					{
						thread t_satman_update(on_satman_update, hwnd);
						t_satman_update.detach();
					}
					break;

				case IDC_SATMAN_OK:
					int chosenTLE = SendMessage(lboActiveSat, LB_GETCURSEL, 0, 0);
					if (chosenTLE >= 0)
					{
						if (ActivateTLE(tlesFromDataFile[chosenTLE]))
						{
							if (1440/meanmotion > 225)
							{
								MessageBox(hwnd, "This satellite's orbital period is greater than 225 minutes, hence the SGP4 algorithm may not be completely accurate.", "Warning", MB_OK | MB_ICONEXCLAMATION);
							}

							terminatorCounter = 99;
							groundTrackCounter = 59;
							SetTimer(hwndMain, IDT_SATTIMER, 1000, NULL);
						}
						else
						{
							MessageBox(hwnd, "Unable to activate selected TLE. Please update the TLEs and try again.", "Error", MB_OK | MB_ICONEXCLAMATION);
						}
					}
					EndDialog(hwnd, 0);
					break;
			}
			break;

		case WM_CLOSE:
			EndDialog(hwnd, 0);
		    break;

		default:
			return FALSE;
	}

	return TRUE;
}
BOOL CALLBACK PassesDlgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch(msg)
	{
		case WM_INITDIALOG:
			{
				lvwPassResults = GetDlgItem(hwnd, IDC_PASSES_RESULTS);
				SendMessage(lvwPassResults, LVM_SETEXTENDEDLISTVIEWSTYLE, 0, (LPARAM)LVS_EX_FULLROWSELECT);
				LVCOLUMN lvCol;
				ZeroMemory(&lvCol, sizeof(LVCOLUMN));
				lvCol.mask = LVCF_WIDTH | LVCF_TEXT | LVCF_SUBITEM;
				lvCol.cx = 130;
				lvCol.pszText = "Rise time / UTC";
				SendMessage(lvwPassResults, LVM_INSERTCOLUMN, 0, (LPARAM)&lvCol);
				lvCol.cx = 70;
				lvCol.pszText = "Rise az / °";
				lvCol.iSubItem = 1;
				SendMessage(lvwPassResults, LVM_INSERTCOLUMN, 1, (LPARAM)&lvCol);
				lvCol.cx = 130;
				lvCol.pszText = "Peak time / UTC";
				lvCol.iSubItem = 2;
				SendMessage(lvwPassResults, LVM_INSERTCOLUMN, 2, (LPARAM)&lvCol);
				lvCol.cx = 70;
				lvCol.pszText = "Peak alt / °";
				lvCol.iSubItem = 3;
				SendMessage(lvwPassResults, LVM_INSERTCOLUMN, 3, (LPARAM)&lvCol);
				lvCol.cx = 130;
				lvCol.pszText = "Set time / UTC";
				lvCol.iSubItem = 4;
				SendMessage(lvwPassResults, LVM_INSERTCOLUMN, 4, (LPARAM)&lvCol);
				lvCol.cx = 70;
				lvCol.pszText = "Set az / °";
				lvCol.iSubItem = 5;
				SendMessage(lvwPassResults, LVM_INSERTCOLUMN, 5, (LPARAM)&lvCol);

				SYSTEMTIME st_passes_init;
				GetSystemTime(&st_passes_init);
				SendMessage(GetDlgItem(hwnd, IDC_PASSES_STARTDATE), DTM_SETSYSTEMTIME, GDT_VALID, (LPARAM)&st_passes_init);
				SendMessage(GetDlgItem(hwnd, IDC_PASSES_STARTTIME), DTM_SETSYSTEMTIME, GDT_VALID, (LPARAM)&st_passes_init);
				SendMessage(GetDlgItem(hwnd, IDC_PASSES_ENDDATE), DTM_SETSYSTEMTIME, GDT_VALID, (LPARAM)&st_passes_init);
				SendMessage(GetDlgItem(hwnd, IDC_PASSES_ENDTIME), DTM_SETSYSTEMTIME, GDT_VALID, (LPARAM)&st_passes_init);
			}
			break;

		case WM_COMMAND:
			switch(LOWORD(wParam))
			{
				case IDC_PASSES_CALCULATE:
					if (epochjd < 1)
					{
						MessageBox(hwnd, "Select a satellite before making pass predictions.", "Error", MB_OK | MB_ICONEXCLAMATION);
						break;
					}

					SYSTEMTIME st_predict_startdate, st_predict_starttime, st_predict_enddate, st_predict_endtime;
					SendMessage(GetDlgItem(hwnd, IDC_PASSES_STARTDATE), DTM_GETSYSTEMTIME, 0, (LPARAM)&st_predict_startdate);
					SendMessage(GetDlgItem(hwnd, IDC_PASSES_STARTTIME), DTM_GETSYSTEMTIME, 0, (LPARAM)&st_predict_starttime);
					SendMessage(GetDlgItem(hwnd, IDC_PASSES_ENDDATE), DTM_GETSYSTEMTIME, 0, (LPARAM)&st_predict_enddate);
					SendMessage(GetDlgItem(hwnd, IDC_PASSES_ENDTIME), DTM_GETSYSTEMTIME, 0, (LPARAM)&st_predict_endtime);

					double UT_predict_start = st_predict_starttime.wHour/24.0 + st_predict_starttime.wMinute/(60*24.0) + st_predict_starttime.wSecond/(60*60*24.0) + st_predict_starttime.wMilliseconds/(1000*60*60*24.0);
					double julianDt_predict_start = JulianDate(st_predict_startdate.wDay, st_predict_startdate.wMonth, st_predict_startdate.wYear, UT_predict_start);
					double UT_predict_end = st_predict_endtime.wHour/24.0 + st_predict_endtime.wMinute/(60*24.0) + st_predict_endtime.wSecond/(60*60*24.0) + st_predict_endtime.wMilliseconds/(1000*60*60*24.0);
					double julianDt_predict_end = JulianDate(st_predict_enddate.wDay, st_predict_enddate.wMonth, st_predict_enddate.wYear, UT_predict_end);
					if (julianDt_predict_end <= julianDt_predict_start)
					{
						MessageBox(hwnd, "Prediction range end must be after start.", "Error", MB_OK | MB_ICONEXCLAMATION);
						break;
					}

					thread t_passes_calculate(on_passes_calculate, hwnd, julianDt_predict_start, julianDt_predict_end);
					t_passes_calculate.detach();
					break;
			}
			break;

		case WM_CLOSE:
			EndDialog(hwnd, 0);
		    break;

		default:
			return FALSE;
	}

	return TRUE;
}
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	HDC hdc;
    PAINTSTRUCT ps;

    switch(msg)
    {
		case WM_CREATE:
			{
				RECT Client;
				GetClientRect(hwnd, &Client);
				CreateWindow("BUTTON", "Preferences", WS_VISIBLE | WS_CHILD | BS_FLAT, Client.left + 50, Client.bottom - 180, 120, 30, hwnd, (HMENU)IDC_PREF, GetModuleHandle(NULL), NULL);
				CreateWindow("BUTTON", "Satellites", WS_VISIBLE | WS_CHILD | BS_FLAT, Client.left + 50, Client.bottom - 140, 120, 30, hwnd, (HMENU)IDC_SATELLITES, GetModuleHandle(NULL), NULL);
				CreateWindow("BUTTON", "Passes", WS_VISIBLE | WS_CHILD | BS_FLAT, Client.left + 50, Client.bottom - 100, 120, 30, hwnd, (HMENU)IDC_PASSES, GetModuleHandle(NULL), NULL);

				if (!ReadSettings())
				{
					MessageBox(hwnd, "No preferences are stored. Please refer to the 'Preferences' menu.", "Warning", MB_OK | MB_ICONEXCLAMATION);
				}

				SYSTEMTIME st_init;
				GetSystemTime(&st_init);
				wstringstream mapStream;
				mapStream << "maps/map_" << unsigned(st_init.wMonth) << ".png";
				mapToUse = mapStream.str();
			}
			break;

		case WM_COMMAND:
			switch(LOWORD(wParam))
			{
				case IDC_PREF:
					DialogBox(GetModuleHandle(NULL), MAKEINTRESOURCE(IDD_PREFERENCES), hwnd, PrefDlgProc);
					break;

				case IDC_SATELLITES:
					DialogBox(GetModuleHandle(NULL), MAKEINTRESOURCE(IDD_SATMANAGER), hwnd, SatManDlgProc);
					break;

				case IDC_PASSES:
					DialogBox(GetModuleHandle(NULL), MAKEINTRESOURCE(IDD_PASSES), hwnd, PassesDlgProc);
					break;
			}
			break;

		case WM_TIMER:
			switch(LOWORD(wParam))
			{
				case IDT_SATTIMER:
					//Time calculations
					SYSTEMTIME st_update;
					GetSystemTime(&st_update);
					int year = st_update.wYear;
					int month = st_update.wMonth;
					int day = st_update.wDay;
					int hour = st_update.wHour;
					int minute = st_update.wMinute;
					int second = st_update.wSecond;
					int millisecond = st_update.wMilliseconds;

					double UT = hour/24.0 + minute/(60*24.0) + second/(60*60*24.0) + millisecond/(1000*60*60*24.0);
					double julianDt = JulianDate(day, month, year, UT);
					double siderealTimeObserver = LMSiderealTime(julianDt, longitude);

					//Solar terminator
					terminatorCounter += 1;
					if (terminatorCounter == 100)
					{
						terminatorCounter = 0;
						terminatorPoints.clear();
						sunresult = SolarPosition(julianDt, siderealTimeObserver);

						for (int i = -180; i <= 180; i++)
						{
							double terminatorSiderealTime = wrapValue(siderealTimeObserver + (i - longitude), 360);
							double terminatorHourAngle = wrapValue(terminatorSiderealTime - sunresult.RA, 360);
							double terminatorLatitude = RADDEG*atan(-cos(DEGRAD*terminatorHourAngle)/tan(DEGRAD*sunresult.Dec));

							terminatorPoints.push_back(MapCoordinates(terminatorLatitude, i));
						}
					}

					//Ground track
					groundTrackCounter += 1;
					if (groundTrackCounter == 60)
					{
						groundTrackCounter = 0;
						groundTrackPoints.clear();
						
						for (double gt_jd = julianDt - 1.1/(2*meanmotion); gt_jd <= julianDt + 1.1/(2*meanmotion); gt_jd += 1/2880.0)
						{
							double groundTrackSiderealTime = LMSiderealTime(gt_jd, longitude);
							PostpropResults groundTrackPoint = PostPropGroundTrack(SGP4(1440*(gt_jd - epochjd)), groundTrackSiderealTime);
							
							groundTrackPoints.push_back(MapCoordinates(groundTrackPoint.Lat, groundTrackPoint.Long));
						}
					}

					//Satellite propagation
					double deltat = julianDt - epochjd;
					result = SGP4(1440*deltat);
					postresult = PostPropFull(result, UT, siderealTimeObserver, deltat);
					
					displaySatResults = true;
					InvalidateRect(hwnd, NULL, TRUE);
					UpdateWindow(hwnd);
					break;
			}
			break;

		case WM_PAINT:
            hdc = BeginPaint(hwnd, &ps);
            update(hwnd, hdc);
            EndPaint(hwnd, &ps);
			break;

		case WM_ERASEBKGND:
			return FALSE;
			break;

        case WM_DESTROY:
			KillTimer(hwnd, IDT_SATTIMER);
            PostQuitMessage(0);
			break;

		case WM_CLOSE:
			KillTimer(hwnd, IDT_SATTIMER);
			PostQuitMessage(0);
			break;
    }

    return DefWindowProc(hwnd, msg, wParam, lParam);
}
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	INITCOMMONCONTROLSEX initex;
	initex.dwSize = sizeof(INITCOMMONCONTROLSEX);
	initex.dwICC  = ICC_PROGRESS_CLASS | ICC_DATE_CLASSES;
	InitCommonControlsEx(&initex);

	GdiplusStartupInput gdiplusStartupInput;
    ULONG_PTR gdiplusToken;
    GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

    //Step 1: Registering the Window Class
	WNDCLASSEX wc;
    ZeroMemory(&wc, sizeof(WNDCLASSEX));
    wc.cbSize = sizeof(WNDCLASSEX);
	wc.style = 0;
    wc.lpfnWndProc = WndProc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
	wc.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_MAINICON));
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)COLOR_WINDOW;
	wc.lpszMenuName = NULL;
    wc.lpszClassName = "WindowClass";
	wc.hIconSm = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_MAINICON));
    
	if(!RegisterClassEx(&wc))
	{
		MessageBox(NULL, "The window could not be registered!", "Error", MB_ICONEXCLAMATION | MB_OK);
		return 0;
	}

    //Step 2: Creating the Window
    hwndMain = CreateWindow("WindowClass", "SpaceShipScanner", WS_SPACESHIPSCANNER, GetSystemMetrics(SM_CXSCREEN)/2 - 300, GetSystemMetrics(SM_CYSCREEN)/2 - 315, 600, 630, NULL, NULL, hInstance, NULL);
	if(hwndMain == NULL)
	{
		MessageBox(NULL, "Window creation failed!", "Error", MB_ICONEXCLAMATION | MB_OK);
		return 0;
	}
	
	ShowWindow(hwndMain, nCmdShow);
	UpdateWindow(hwndMain);

    //Step 3: The Message Loop
	MSG Msg;
	while(GetMessage(&Msg, NULL, 0, 0))
    {
		TranslateMessage(&Msg);
		DispatchMessage(&Msg);
    }
	
	GdiplusShutdown(gdiplusToken);
	return Msg.wParam;
}
