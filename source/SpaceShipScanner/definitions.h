#ifndef DEFINITIONS_H
#define DEFINITIONS_H

//Constants
#define G			6.6726e-11
#define Mearth		5.9737e24
#define RearthEq	6378.1370
#define RearthPo	6356.7523
#define ke			0.743669161e-1
#define J2			1.082616e-3
#define J3			-0.253881e-5
#define J4			-1.65597e-6
#define DEGRAD		M_PI/180
#define RADDEG		180/M_PI
#define ae			1

//Constants relying on other constants
const double k2 = 0.5*J2*ae*ae;
const double A30 = -J3*ae*ae*ae;
const double k4 = -3*J4*ae*ae*ae*ae/8;
const double fEarth = 1 - RearthPo/RearthEq;
const double eSqEarth = 2*fEarth - fEarth*fEarth;

//Structures
struct PropagationResults
{
	double posX, posY, posZ;
	double velX, velY, velZ;
	double perigeeKM, apogeeKM;
};
struct PostpropResults
{
	double UT, LST, DeltaT;
	double RA, Dec, Range, Alt, Az;
	double Lat, Long, Altitude, Velocity;
};
struct SingleTLE
{
	string line1, line2, line3;
};
struct SunResults
{
	double RA, Dec;
	double subSolarLat, subSolarLong;
};
struct MapCoordinates
{
	double Lat, Long;

	MapCoordinates(double in_lat, double in_long)
	{Lat = in_lat; Long = in_long;};
};
struct PassPrediction
{
	double riseTimeJD, riseAzimuth;
	double peakTimeJD, peakAltitude;
	double setTimeJD, setAzimuth;
};
struct CalendarDate
{
	double day;
	unsigned month;
	int year;
};

//Active TLE parameters
string satName;
double epochjd, bstar, inclination, raascendingnode, eccentricity, meanmotion, mmDo2, mmDDo6, argumentperigee, meananomaly;
double inclinationRad, raascendingnodeRad, argumentperigeeRad, meananomalyRad, meanmotionRadPerMin, mmDo2RadPerMinSq, mmDDo6RadPerMinCu;

//Misc. parameters
const string tleDataFiles []  = {"stations.txt", "visual.txt", "weather.txt", "resource.txt", "iridium.txt", "science.txt", "amateur.txt", "cubesat.txt"};
double latitude = 0, longitude = 0, altitude = 0;
unsigned terminatorCounter = 99, groundTrackCounter = 59;
bool displaySatResults = false;
wstring mapToUse = L"maps/map_1.png";
HWND hwndMain = NULL, cboTLEFile = NULL, lboActiveSat = NULL, lvwPassResults = NULL;
vector <SingleTLE> tlesFromDataFile;
vector <MapCoordinates> terminatorPoints, groundTrackPoints;
PropagationResults result;
PostpropResults postresult;
SunResults sunresult;

#endif
