/* PROJECT 2: NUMERICAL INTEGRATION
 * This is a simple program which calculates a 6-dimensional integral
 * using Gaussian quadrature with Legendre and Laguerre polynomials,
 * and Monte Carlo metod with and without importance sampling.
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <random>
//#include <libiomp/omp.h>

#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10

using namespace std;


// FUNCTIONS
double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
double int_function_spherical(double r1,double r2,double theta1,double theta2,double phi1,double phi2);
double int_function_spherical_MC(double r1,double r2,double theta1,double theta2,double phi1,double phi2);
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double gammln(double xx);


// MAIN PROGRAM:
int main()
{
     #pragma omp parallel

     // GAUSS-LEGENDRE
     clock_t start_gauleg, finish_gauleg;               // declare start and final time for Gauss-Legendre
     start_gauleg = clock();                            // start the clock

     int N = 25;                                        // number of integration points
     double a = -3.0;                                   // lower limit
     double b = 3.0;                                    // upper limit

     double *x = new double [N];                        // array for mesh points
     double *w = new double [N];                        // array for weights
     gauleg(a,b,x,w,N);                                 // calculate mesh points and weights

     // Evaluate the integral with the Gauss-Legendre method
     double int_gauss = 0.;                             // initialize the sum

     //#pragma omp parallel for reduction(+:int_gauss) private(j,k,l,m,n)

     for (int i=0;i<N;i++){
         for (int j = 0;j<N;j++){
             for (int k = 0;k<N;k++){
                 for (int l = 0;l<N;l++){
                     for (int m = 0;m<N;m++){
                         for (int n = 0;n<N;n++){
                             int_gauss += w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
                  }}}}}
             }

     finish_gauleg = clock();                           // stop the clock


     // GAUSS-LAGUERRE-LEGENDRE COMBO
     clock_t start_gaulag, finish_gaulag;               // declare start and final time for Gauss-Laguerre-Legendre combo
     start_gaulag = clock();                            // start clock

     int N_GL = 25;                                     // number of integration points

     // r
     double alf = 2.0;
     double *x_r = new double [N_GL];                   // array for mesh points
     double *w_r = new double [N_GL];                   // array for weights
     gauss_laguerre(x_r,w_r,N_GL,alf);                  // calculate mesh points and weights


     // PHI
     double a_phi1 = 0;                                 // lower limit
     double b_phi1 = 2*M_PI;                            // upper limit
     double *x_phi = new double [N_GL];                 // array for mesh points
     double *w_phi = new double [N_GL];                 // array for weights
     gauleg(a_phi1,b_phi1,x_phi,w_phi,N_GL);            // calculate mesh points and weights


     // THETA
     double a_theta1 = 0;                               // lower limit
     double b_theta1 = M_PI;                            // upper limit
     double *x_theta = new double [N_GL];               // array for mesh points
     double *w_theta = new double [N_GL];               // array for weights
     gauleg(a_theta1,b_theta1,x_theta,w_theta,N_GL);    // calculate mesh points and weights


     // Solving the integral
     double int_spherical = 0.0;                        // initialize the sum

     //#pragma omp parallel for reduction(+:int_spherical) private(j,k,l,m,n)

     for (int i=1;i<=N_GL;i++){
         for (int j = 1;j<=N_GL;j++){
             for (int k = 0;k<N_GL;k++){
                 for (int l = 0;l<N_GL;l++){
                     for (int m = 0;m<N_GL;m++){
                         for (int n = 0;n<N_GL;n++){
                             int_spherical += w_r[i]*w_r[j]*w_theta[k]*w_theta[l]*w_phi[m]*w_phi[n]*int_function_spherical(x_r[i],x_r[j],x_theta[k],x_theta[l],x_phi[m],x_phi[n]);
                         }}}}}
             }

     int_spherical = int_spherical/1024.0;              // final calculation of integral
     finish_gaulag = clock();                           // stop the clock


     // MONTE-CARLO
     // Evaluate the integral with the a crude Monte-Carlo method
     clock_t start_MC, finish_MC;                       // declare start and final time for Monte Carlo
     start_MC = clock();                                // start the clock

     default_random_engine generator;                   // start random number generator
     uniform_real_distribution<double> distribution(0.0,1.0);   // pick type of random distribution

     int N_MC = 50000000;                               // number of calculations
     double *y = new double [N_MC];
     double fx;
     double MCint = 0;                                  // initialize the sum
     double MCintsqr2 = 0;                              // initialize standard deviation
     double length = 3;
     double jacobidet = pow((2*length),6);              // Jacobi determinant

     //#pragma omp parallel for reduction(+:MCint,MCintsqr2) private(j)

     for(int i=1;i<=N_MC;i++){
         for(int j=0;j<6;j++){
             y[j] = -length + 2*length*distribution(generator);
         }

         fx = int_function(y[0],y[1],y[2],y[3],y[4],y[5]);
         MCint += fx;
         MCintsqr2 += fx*fx;
     }

     MCint = jacobidet*MCint/((double) N_MC);           // finalize integral
     MCintsqr2 = jacobidet*MCintsqr2/((double) N_MC);   // finalize standard deviation
     double variance = MCintsqr2 - MCint*MCint;         // calculate variance

     finish_MC = clock();                               // stop the clock


     // Importance sampling
     clock_t start_MCi, finish_MCi;                     // declare start and final time for Monte Carlo with importance sampling
     start_MCi = clock();                               // start the clock

     int N_MCi = 1000000;                               // number of calculations

     double *z = new double [N_MCi];
     double fx_exp;
     double MCint_exp = 0;                              // initialize the sum
     double MCintsqr2_exp = 0;                          // initialize standard deviation
     double jacobidet_exp = 4*pow(acos(-1.),4.0)*1./16; // Jacobi determinant

     //#pragma omp parallel for reduction(+:MCint,MCintsqr2) private(j)

     for(int i=1;i<=N_MCi;i++){
         for(int j=0;j<6;j++){
             z[j] = distribution(generator);
         }

         fx_exp = int_function_spherical_MC(-0.25*log(1-z[0]),-0.25*log(1-z[1]),M_PI*z[2],M_PI*z[3],2*M_PI*z[4],2*M_PI*z[5]);
         MCint_exp += fx_exp;
         MCintsqr2_exp += fx_exp*fx_exp;
     }

     MCint_exp = jacobidet_exp*MCint_exp/((double) N_MCi); // finalize integral
     MCintsqr2_exp = jacobidet_exp*MCintsqr2_exp/((double) N_MCi); // finalize standard deviation
     double variance_exp = MCintsqr2_exp - MCint_exp*MCint_exp; // calculate variance

     finish_MCi = clock();                              // stop the clock


     // FINAL OUTPUT
     cout << "INPUT:" << endl;
     cout << "N = " << N << endl;
     cout << "N_GL = " << N_GL << endl;
     cout << "N_MC = " << N_MC << endl;
     cout << "N_MCi = " << N_MCi << endl;
     cout << "a = " << a << " (lower limit)" << endl;
     cout << "b = " << b << " (upper limit)" << endl;

     cout << endl << "RESULTS:" << endl;
     cout << "Gauss-Legendre "<< "\t" << setprecision(15)  << int_gauss << endl;
     cout << "Gauss-Laguerre " << "\t" << setprecision(15) << int_spherical << endl;
     cout << "Monte Carlo " << "\t" << setprecision(15) << MCint << " (variance = " << variance << ")"<< endl;
     cout << "Monte Carlo (imp.) " << "\t" << setprecision(15) << MCint_exp << " (variance = " << variance_exp << ")"<< endl;

     cout << endl << "Exact answer " << "\t" << 5*M_PI*M_PI/(16*16) << endl << endl;

     cout << "TIME USAGE:" << endl;
     cout << "Gauss-Legendre " << "\t" << ((float)(finish_gauleg - start_gauleg)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time
     cout << "Gauss-Laguerre " << "\t" << ((float)(finish_gaulag - start_gaulag)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time
     cout << "Monte Carlo " << "\t" << ((float)(finish_MC - start_MC)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time
     cout << "Monte Carlo (imp.) " << "\t" << ((float)(finish_MCi - start_MCi)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time


     // Clear memory
     delete [] x;
     delete [] w;
     delete [] x_r;
     delete [] w_r;
     delete [] x_theta;
     delete [] w_theta;
     delete [] x_phi;
     delete [] w_phi;
     delete [] y;
     delete [] z;

     return 0;
}
// End of main program


// DEFINITION OF FUNCTIONS:

// Integral (cartesian):
double int_function(double x1, double y1, double z1, double x2, double y2, double z2){
    double alpha = 2.0;
    double exp1 = -2*alpha*sqrt(x1*x1 + y1*y1 + z1*z1);
    double exp2 = -2*alpha*sqrt(x2*x2 + y2*y2 + z2*z2);
    double deno = sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));

    if(deno < pow(10.,-6.)){
        return 0;}
    else
        return exp(exp1+exp2)/deno;
    }

// Integral (spherical without r's):
double int_function_spherical(double r1,double r2,double theta1,double theta2,double phi1,double phi2){
    double cos_beta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double deno = r1*r1 + r2*r2 - 2*r1*r2*cos_beta;

    if(deno < pow(10.,-10.) || deno<0){
        return 0;}
    else
        return sin(theta1)*sin(theta2)/sqrt(deno);
}

// Integral (spherical with r's):
double int_function_spherical_MC(double r1,double r2,double theta1,double theta2,double phi1,double phi2){
    double cos_beta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double deno = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos_beta);

    if(deno < pow(10.,-6.)  || isnan(deno)){
        return 0;}
    else
        return r1*r1*r2*r2*sin(theta1)*sin(theta2)/deno;
}


/* GAUSS-LEGENDRE
 * The function gauleg() takes the lower and upper limits of integration x1, x2,
 * calculates and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
 * of length n of the Gauss--Legendre n--point quadrature formulae.
*/

void gauleg(double x1, double x2, double x[], double w[], int n){
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));
      do {
          p1 =1.0;
          p2 =0.0;

     for(j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /* p1 is now the desired Legrendre polynomial. Next compute ppp its derivative by standard relation
        * involving also p2, polynomial of one lower order. */

     pp = n * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

      /* Scale the root to the desired interval and put in its symmetric counterpart.
       * Compute the weight and its symmetric counterpart */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
}



// GAUSS-LAGUERRE
void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}



double gammln( double xx){
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
