//#include <Arduino.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include "Gps_interp_path_reg.h"

//////////////////////////////////spremenjivke za preracune/////////////////////
#define pi 3.141592653589793
#define pi2 2*3.141592653589793
#define DEG2RAD (pi / 180.0)
//////////////////////////////define za L0/////////////////////////////////////
#define g 9.81
#define Max_roll 30
double phi_cmd;
//int current_wp = 0;  // start at first leg
double d1;
double d2;
double L1;
double lambda;
double SOG;  // ground speed
double COG; // course over ground [deg]
//int total_wp;
////////////////////spremenljivke za gps ter interpolacijo///////////

GPSdata gpsIZxyz,gpsIZxyzMatlab,gpsdata1;
GPSdata gpsIZxyzCIRC;
GPSdata gpstest;
xyz xyzInterp,xyzCIRC,xyzIzMatlab,xyzOut1;
ENU enuOut,enuInterp;

double refx, refy, refz;
double lat0, lon0;
// reference is always index 0-- za ecef-enu
void initReference() {
    refx = xyzOut1.X[0];
    refy = xyzOut1.Y[0];
    refz = xyzOut1.Z[0];

    lat0 = gpsdata1.lat[0];
    lon0 = gpsdata1.lon[0];
}
//---------------------------------------------------------------------------
// http://www.navipedia.net/index.php/Ellipsoidal_and_Cartesian_Coordinates_Conversion
//---------------------------------------------------------------------------
// WGS84(a,b,h) = (long,lat,alt) [rad,rad,m]
// XYZ(x,y,z) [m]
//---------------------------------------------------------------------------
const double  _earth_a=6378137.00000;   // [m] WGS84 equator radius
const double  _earth_b=6356752.31414;   // [m] WGS84 epolar radius
const double  _earth_e=8.1819190842622e-2; //  WGS84 eccentricity
//const double  _earth_e=sqrt(1.0-((_earth_b/_earth_a)*(_earth_b/_earth_a)));
const double  _earth_ee=_earth_e*_earth_e;
//---------------------------------------------------------------------------
const double kmh=1.0/3.6;               // [km/h] -> [m/s]
//---------------------------------------------------------------------------

void XYZtoWGS84(double *lon,double *lat,double *alt,double x,double y,double z)
    {
    int i;
    double  a,b,h,l,n,db,s;
    //a=atan(y/x);
    a = atan2(y, x);  // Use atan2 for correct quadrant handling
    l=sqrt((x*x)+(y*y));
    
   // Handle near-polar cases (avoid division by zero)
    if (l < 1e-9) {
        b = (z > 0) ? (pi / 2) : (-pi / 2);  // Pole
        h = fabs(z) - _earth_b;
    } else {
        // Correct initial latitude estimate
        b = atan(z / (l * (1.0 - _earth_ee)));  // FIXED: Proper parentheses
        
        // Iterate to improve accuracy
        for (i = 0; i < 100; i++) {
            s = sin(b);
            db = b;
            n = _earth_a / sqrt(1.0 - (_earth_ee * s * s));  // Prime vertical radius
            h = (l / cos(b)) - n;
            b = atan(z / (l * (1.0 - (_earth_ee * n / (n + h)))));
            db = fabs(db - b);
            if (db < 1e-12) break;  // Convergence check
        }
    }
    
    if (b > 0.5 * pi) b -= pi2;  // Normalize latitude
    *lon = a;
    *lat = b;
    *alt = h;
}
  
//---------------------------------------------------------------------------
void WGS84toXYZ(double *x,double *y,double *z,double lon, double lat, double alt)  //double *abh
    {
    double  a,b,h,l,c,s;
    a = lon;  
    b = lat;
    h = alt;
    c=cos(b);
    s=sin(b);
    // WGS84 from eccentricity
    l=_earth_a/sqrt(1.0-(_earth_ee*s*s));
    *x=(l+h)*c*cos(a);
    *y=(l+h)*c*sin(a);
    *z=(((1.0-_earth_ee)*l)+h)*s;
    }

//---------------------------------------------------------------------------
void XYZ_interpolate(double *ptx, double *pty,double *ptz,double tocka1x,double tocka1y,double tocka1z,double tocka2x,double tocka2y,double tocka2z,int N)
    {
    const double  mz=_earth_a/_earth_b;
    const double _mz=_earth_b/_earth_a;
    double px[N],py[N],pz[N],r,r0,r1,t;
    
    // compute spherical radiuses of input points
    r0=sqrt((tocka1x*tocka1x)+(tocka1y*tocka1y)+(tocka1z*tocka1z*mz*mz)); //XYZ
    r1=sqrt((tocka2x*tocka2x)+(tocka2y*tocka2y)+(tocka2z*tocka2z*mz*mz)); //XYZ
    // linear interpolation
   for(int i = 0; i < N; i++) {
        double t = (double)i / (double)(N - 1);
        ptx[i] = tocka1x + (tocka2x - tocka1x) * t; // X
        pty[i] = tocka1y + (tocka2y - tocka1y) * t; // Y
        ptz[i] = tocka1z + (tocka2z - tocka1z) * t; // Z (NO SCALING)
    
   }
   
  }

  // Convert ECEF -> ENU
void ECEFtoENU(double *e, double *n, double *u,
               double x, double y, double z,
               double x0, double y0, double z0,
               double lat0, double lon0) 
{
    double dx = x - x0;
    double dy = y - y0;
    double dz = z - z0;

    double sinLat = sin(lat0);
    double cosLat = cos(lat0);
    double sinLon = sin(lon0);
    double cosLon = cos(lon0);

    *e = -sinLon * dx + cosLon * dy;
    *n = -sinLat * cosLon * dx - sinLat * sinLon * dy + cosLat * dz;
    *u =  cosLat * cosLon * dx + cosLat * sinLon * dy + sinLat * dz;
}


// Convert ENU -> ECEF
void ENUtoECEF(double *x, double *y, double *z,
               double e, double n, double u,
               double x0, double y0, double z0,
               double lat0, double lon0) 
{
    double sinLat = sin(lat0);
    double cosLat = cos(lat0);
    double sinLon = sin(lon0);
    double cosLon = cos(lon0);

    *x = x0 - sinLon * e - sinLat * cosLon * n + cosLat * cosLon * u;
    *y = y0 + cosLon * e - sinLat * sinLon * n + cosLat * sinLon * u;
    *z = z0 + cosLat * n + sinLat * u;
}

// chord length parameterization
void chord_length(const double *x, const double *y, int n, double *t) {
    t[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double dx = x[i] - x[i-1];
        double dy = y[i] - y[i-1];
        t[i] = t[i-1] + sqrt(dx*dx + dy*dy);
    }
}

// novi(2) MAKIMA derivative computation with straight boundary slopes
void makima_derivatives(const double *t, const double *v, int n, double *d) {
    double m[n-1];
    for (int j = 0; j < n-1; j++) {
        double h = t[j+1] - t[j];
        m[j] = (v[j+1] - v[j]) / h;
    }

    // extended slope array with padding
    double M[n+4];
    for (int j = 0; j < n-1; j++) M[j+2] = m[j];
    M[0] = M[1] = M[2];
    M[n+2] = M[n+1] = M[n] = M[n-2+2]; 
    M[n+3] = M[n+2];

    // compute derivatives
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            // enforce forward slope at start
            d[i] = (v[1] - v[0]) / (t[1] - t[0]);
        } else if (i == n-1) {
            // enforce backward slope at end
            d[i] = (v[n-1] - v[n-2]) / (t[n-1] - t[n-2]);
        } else {
            double m_im1 = M[i+1];
            double m_i   = M[i+2];
            double m_ip1 = M[i+3];
            double m_ip2 = M[i+4];

            double w1 = fabs(m_ip2 - m_ip1) + fabs(m_i - m_im1);
            double w2 = fabs(m_ip2 - m_i)   + fabs(m_i - m_im1);
            double w3 = fabs(m_ip1 - m_i)   + fabs(m_ip1 - m_im1);
            double w4 = fabs(m_ip1 - m_i)   + fabs(m_ip2 - m_i);

            double W = w1 + w2 + w3 + w4;

            if (W > 0.0)
                d[i] = (w1*m_im1 + w2*m_i + w3*m_ip1 + w4*m_ip2) / W;
            else
                d[i] = 0.25*(m_im1 + m_i + m_ip1 + m_ip2);
        }
    }
}

// PCHIP slope calculation
void pchip_slopes(const double *t, const double *v, int n, double *d) {
    double m[n-1];
    for (int j=0; j<n-1; j++) {
        double h = t[j+1] - t[j];
        m[j] = (v[j+1] - v[j]) / h;
    }

    // endpoint slopes
    d[0] = m[0];
    d[n-1] = m[n-2];

    for (int k=1; k<n-1; k++) {
        if (m[k-1]*m[k] > 0) {
            double w1 = 2*(t[k+1]-t[k]) + (t[k]-t[k-1]);
            double w2 = (t[k+1]-t[k]) + 2*(t[k]-t[k-1]);
            d[k] = (w1 + w2) / ( (w1/m[k-1]) + (w2/m[k]) );
        } else {
            d[k] = 0.0;
        }
    }
}


// evaluate cubic Hermite on one interval [t0,t1]
double cubic_eval(double a, double b, double c, double d,
                                double u) {
    return ((a*u + b)*u + c)*u + d;
}

// build cubic coefficients for one interval
void cubic_coeff(double v0, double v1, double d0, double d1,
                 double t0, double t1, double *a, double *b,
                 double *c, double *d) {
    double h = t1 - t0;
    double m = (v1 - v0)/h;
    *a = (d0 + d1 - 2*m)/(h*h);
    *b = (3*m - 2*d0 - d1)/h;
    *c = d0;
    *d = v0;
}
///////////////////////////////// END MAKIMA INTERPOLATION/////////////////////////////////////////////


void L0_guidence(double SOG,double COG,double plane_E,double plane_N,double A_E,double A_N,double B_E,double B_N){
double R;
double Va=SOG;
double Vg_vec[2];
double Vg_norm;
double L0;
double t_hat[2];
double r[2];
double s;
double Q[2];
double d_vect[2];
double d;
double gamma;
double alpha;


double T[2];
double dT[2];
double d_hat[2];
double v_hat[2];
double cross2D;
double dot2D;
double eta;
double alat;
  //MIN turn radius glede na omejotve ter hitrost
R=(Va*Va)/(g*tan(Max_roll*DEG2RAD));
L0=R;
double c0=B_E-A_E;
double c1=B_N-A_N;

// before building t_hat--clampamo vrednost protection bred NaN
//double dx = B_E - A_E, dy = B_N - A_N;
double seg = hypot(c0, c1);
if (seg < 1e-6) { phi_cmd = 0; return; }          // early safe return

t_hat[0] = (c0)/ seg ; //tangenta na pot (sqrt((c0*c0)+(c1*c1)))
t_hat[1] = (c1)/ seg ; //tangenta na pot (sqrt((c0*c0)+(c1*c1)))

 // % --- Ground velocity---
    Vg_vec[0] = Va*sin(COG*DEG2RAD);  //zamenjaqmo sin in cos zaradi cog // East component
    Vg_vec[1] = Va*cos(COG*DEG2RAD);   // North component
    
    Vg_norm = sqrt((Vg_vec[0]*Vg_vec[0])+(Vg_vec[1]*Vg_vec[1]));
// protection pred NaN
    if (Vg_norm < 0.5) { phi_cmd = 0; return; }       // low speed → no bank

    //--- Along-track projection ---
    r[0] = plane_E - A_E;   //vektor r od tocke A do letala
    r[1] = plane_N - A_N;

//C = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
    s = r[0]*t_hat[0]+r[1]*t_hat[1];   //skalarni produkt

    Q[0] = A_E + s * t_hat[0];   //closest point on segment
    Q[1] = A_N + s * t_hat[1];   //closest point on segment

    d_vect[0] = plane_E - Q[0];
    d_vect[1] = plane_N - Q[1];
    d = sqrt((d_vect[0]*d_vect[0])+(d_vect[1]*d_vect[1]));

   // % fernandez logivc
    // % --- Geometry ---
        gamma = atan2(B_N-A_N, B_E-A_E);
        alpha = atan2(plane_N-A_N, plane_E-A_E);
        lambda = gamma - alpha;
        d1 = hypot(plane_E-A_E, plane_N-A_N); //% A->aircraft
        d2 = hypot(plane_E-B_E, plane_N-B_N); //% B->aircraft

        L1 = sqrt((d*d) + (L0*L0));

         if (L1 < 1e-6) { phi_cmd = 0; return; }         // protect alat/L1

        T[0] = Q[0] + (L1 * t_hat[0]);  // % reference point offset
        T[1] = Q[1] + (L1 * t_hat[1]);  // % reference point offset


        // % --- Compute guidance commands ---
    dT[0] = T[0] - plane_E;
    dT[1] = T[1] - plane_N;

    double dTn = hypot(dT[0], dT[1]);
    if (dTn < 1e-6) { phi_cmd = 0; return; }          // already at T

    d_hat[0] = dT[0]/sqrt((dT[0]*dT[0]) + (dT[1]*dT[1]));
    d_hat[1] = dT[1]/sqrt((dT[0]*dT[0]) + (dT[1]*dT[1]));
    
    v_hat[0] = Vg_vec[0]/Vg_norm;
    v_hat[1] = Vg_vec[1]/Vg_norm;

    cross2D = v_hat[0]*d_hat[1] - v_hat[1]*d_hat[0];
    dot2D   = (v_hat[0]*d_hat[0]) + (v_hat[1]*d_hat[1]);       //dot(v_hat,d_hat);

    eta = atan2(cross2D,dot2D);

    if (eta > M_PI/2) eta = M_PI/2;
    if (eta < -M_PI/2) eta = -M_PI/2;
    
    alat = 2*(Vg_norm*Vg_norm) / L1 * sin(eta);    //  % classic + L0
    
    phi_cmd = atan(alat/g);

    if (!isfinite(phi_cmd)) { phi_cmd = 0; return; }

    // limit roll
    if (phi_cmd > Max_roll * DEG2RAD) phi_cmd = Max_roll * DEG2RAD;
    if (phi_cmd < -Max_roll * DEG2RAD) phi_cmd = -Max_roll * DEG2RAD;
    /*
    Serial.print("DBG: Vg="); Serial.print(Vg_norm,3);
    Serial.print(" d="); Serial.print(d,3);
    Serial.print(" L1="); Serial.print(L1,3);
    Serial.print(" eta(deg)="); Serial.print(eta * 180.0 / M_PI,3);
    Serial.print(" phi_cmd(deg)="); Serial.print(phi_cmd * 180.0 / M_PI,3);
    Serial.println();
    */
  }

////////////////////////end L0 guidence/////////////////////
double zoneCalc(double plane_E,double plane_N,double A_E,double A_N,double B_E,double B_N,double SOG) {
double values[5];
double s;
double e_skalar;
double e[2];
double c0=B_E-A_E;
double c1=B_N-A_N;
double t_hat[2];
double r[2];
double gamma;
double alpha;
double lambda;
double R;
double l0;
double d_vect[2];
double d;
double Q[2];
double L1_zone;
//MIN turn radius glede na omejotve ter hitrost
R=(SOG*SOG)/(g*tan(Max_roll*DEG2RAD));
l0=R;

// before building t_hat--clampamo vrednost protection bred NaN
//double dx = B_E - A_E, dy = B_N - A_N;
double seg = hypot(c0, c1);
if (seg < 1e-6) { phi_cmd = 0; return; }          // early safe return

t_hat[0] = (c0)/ seg ; //tangenta na pot (sqrt((c0*c0)+(c1*c1)))
t_hat[1] = (c1)/ seg ; //tangenta na pot (sqrt((c0*c0)+(c1*c1)))

 //--- Along-track projection ---
    r[0] = plane_E - A_E;   //vektor r od tocke A do letala
    r[1] = plane_N - A_N;

//C = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
    s = r[0]*t_hat[0]+r[1]*t_hat[1];   //skalarni produkt
    e[0]= r[0]+t_hat[0]*s;
    e[1]= r[1]+t_hat[1]*s;
    e_skalar=sqrt((e[0]*e[0])+(e[1]*e[1]));

    /////////ven iz L0 //////////////
    Q[0] = A_E + s * t_hat[0];   //closest point on segment
    Q[1] = A_N + s * t_hat[1];   //closest point on segment

    d_vect[0] = plane_E - Q[0];
    d_vect[1] = plane_N - Q[1];
    d = sqrt((d_vect[0]*d_vect[0])+(d_vect[1]*d_vect[1]));

   // % fernandez logivc
    // % --- Geometry ---
        gamma = atan2(B_N-A_N, B_E-A_E);
        alpha = atan2(plane_N-A_N, plane_E-A_E);
        lambda = gamma - alpha;
        d1 = hypot(plane_E-A_E, plane_N-A_N); //% A->aircraft
        d2 = hypot(plane_E-B_E, plane_N-B_N); //% B->aircraft

        L1_zone = sqrt((d*d) + (l0*l0));

      if (L1_zone < 1e-6) { 
        L1_zone = 0;
    }

values[0]=s;
values[1]=e_skalar;
values[2]=d1;
values[3]=d2;
values[4]=L1_zone;
return values[5];
}
