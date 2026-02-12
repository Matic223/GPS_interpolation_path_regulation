#ifndef GPS_INTERP_PATH_REG
#define GPS_INTERP_PATH_REG
//////////////////////////////////spremenjivke za preracune/////////////////////
#define pi 3.141592653589793
#define pi2 2*3.141592653589793
#define DEG2RAD (pi / 180.0)

#define g 8.91
#define Max_roll 30
extern double phi_cmd;
//extern int current_wp;  // start at first leg
extern double d1;
extern double d2;
extern double L1;
extern double lambda;
extern double SOG;  // ground speed
extern double COG; // course over ground [deg]
//extern int total_wp;

struct GPSdata{
    double lon[100];
    double lat[100];
    double alt[100];  
   };
struct xyz{
    double X[100];
    double Y[100]; 
    double Z[100]; 
   };
struct ENU{
    double E[100];
    double n[100]; 
    double U[100]; 
   };

   struct ZoneVals {
  double s;       // along-track projection
  double e_skalar;       // cross-track distance (>=0)
  double d1_zone;      // distance to A
  double d2_zone;      // distance to B
  double L1_zone; // your zone boundary distance
  double seg;     // length of segment AB (useful!)
  double Q1;
  double Q2;
  double s_raw;
  double L0;
  double d;
};

extern GPSdata gpsIZxyz,gpsIZxyzMatlab,gpsdata1;
extern GPSdata gpsIZxyzCIRC;
extern GPSdata gpstest;
extern xyz xyzInterp,xyzCIRC,xyzIzMatlab,xyzOut1;
extern ENU enuOut,enuInterp;
 // reference is always index 0-- za ecef-enu
extern double refx;
extern double refy;
extern double refz;

extern double lat0;
extern double lon0;
/*
// interpolacija pchip, makima
extern const int N;
extern const int NSAMP;
extern double t[];
// derivatives on ENU
extern double dE[], dN_[];
// build cubic coefficients on ENU
extern double aE[], bE[], cE[], dE0[];
extern double aN[], bN[], cN[], dN0[];
// sample curve in ENU
extern double tmin, tmax;
*/
void XYZtoWGS84(double *lon,double *lat,double *alt,double x,double y,double z);
void WGS84toXYZ(double *x,double *y,double *z,double lon, double lat, double alt);
void XYZ_interpolate(double *ptx, double *pty,double *ptz,double tocka1x,double tocka1y,double tocka1z,double tocka2x,double tocka2y,double tocka2z,int N);
void ECEFtoENU(double *e, double *n, double *u,
               double x, double y, double z,
               double x0, double y0, double z0,
               double lat0, double lon0);
void ENUtoECEF(double *x, double *y, double *z,
               double e, double n, double u,
               double x0, double y0, double z0,
               double lat0, double lon0);
void chord_length(const double *x, const double *y, int n, double *t);
void makima_derivatives(const double *t, const double *v, int n, double *d);
void pchip_slopes(const double *t, const double *v, int n, double *d);
double cubic_eval(double a, double b, double c, double d,double u);
void cubic_coeff(double v0, double v1, double d0, double d1,
                 double t0, double t1, double *a, double *b,
                 double *c, double *d);
void L0_guidence(double SOG,double COG,double plane_E,double plane_N,double A_E,double A_N,double B_E,double B_N);
void initReference();
ZoneVals zoneCalc(double plane_E,double plane_N,double A_E,double A_N,double B_E,double B_N,double SOG1);
static inline double clampd(double x, double lo, double hi);
double guidanceToPoint(double SOG, double COG_deg,
                       double plane_E, double plane_N,
                       double target_E, double target_N,
                       double L, double phiMaxRad);
int closestSegmentIndex(double plane_E, double plane_N, int total_wp);                      
#endif