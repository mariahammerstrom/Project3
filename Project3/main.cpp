//   This is a simple program which tests Gaussian quadrature using Legendre and Laguerre polynomials

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <random>

#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10

using namespace std;


// FUNCTIONS

double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
double int_function_spherical(double r1,double r2,double theta1,double theta2,double phi1,double phi2);
//double int_mc(double,double);
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double gammln(double);


int main()
{
     clock_t start, finish; // declare start and final time for Armadillo solver
     start = clock();

     int N = 10;
     double a = -1.0;
     double b = 1.0;

     // GAUSS-LEGENDRE
     // Mesh points weights and function values
     double *x = new double [N];
     double *w = new double [N];

     gauleg(a,b,x,w,N); // mesh points and weights

     // Evaluate the integral with the Gauss-Legendre method
     double int_gauss = 0.; // initialize the sum
     #pragma omp for reduction(+:int_gauss) private(i,j,k,l,m,n)
     for (int i=0;i<N;i++){
         for (int j = 0;j<N;j++){
             for (int k = 0;k<N;k++){
                 for (int l = 0;l<N;l++){
                     for (int m = 0;m<N;m++){
                         for (int n = 0;n<N;n++){
                             int_gauss += w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
                  }}}}}
             }



     // GAUSS-LAGUERRE
     // Mesh points weights and function values
     double *xgl = new double [N+1];
     double *wgl = new double [N+1];

     double alf = 1.0;
     gauss_laguerre(xgl,wgl,N,alf);

     // Evaluate the integral with the Gauss-Laguerre method
     double int_gausslag = 0.; // initialize the sum
     #pragma omp for reduction(+:int_gausslag) private(i,j,k,l,m,n)
     for (int i=0;i<N;i++){
         for (int j = 0;j<N;j++){
             for (int k = 0;k<N;k++){
                 for (int l = 0;l<N;l++){
                     for (int m = 0;m<N;m++){
                         for (int n = 0;n<N;n++){
                             int_gausslag += wgl[i]*wgl[j]*wgl[k]*wgl[l]*wgl[m]*wgl[n]*int_function_spherical(xgl[i],xgl[j],xgl[k],xgl[l],xgl[m],xgl[n]);
                  }}}}}
             }


     // MONTE-CARLO
     //double MCint = 0; // initialize the sum
     //double MCintsqr2 = 0.;
     //N = 100000;
     //double invers_period = 1./RAND_MAX; // initialise the random number generator
     //srand(time(NULL)); // This produces the so-called seed in MC jargon

     // Evaluate the integral with the a crude Monte-Carlo method

     //#pragma omp for reduction(+:MCint,MCintsqr2) private(i)
     default_random_engine generator;
     uniform_real_distribution<double> distribution(0.0,1.0);
     exponential_distribution<double> expdistribution(0.2928552*0.19276571); // Correct for length 3

     //int under = 0;
     /*for ( int i = 0; i <= N; i++){
         double x1 = distribution(generator); // old way: double(rand())*invers_period;
         double x2 = distribution(generator);
         double y1 = distribution(generator);
         double y2 = distribution(generator);
         double z1 = distribution(generator);
         double z2 = distribution(generator);
         double fx = int_function(x1,x2,y1,y2,z1,z2);*/

     // New Monte Carlo
     N = 100000;
     double *y = new double [N];
     double *z = new double [N];
     double fx;
     double fx_exp;
     double MCint = 0;
     double MCintsqr2 = 0;
     double length = 3;
     double jacobidet = pow((2*length),6);

     // importance sampling
     for(int i=1;i<=N;i++){
         for(int j=0;j<6;j++){
             //y[j] = -length + 2*length*rand()/RAND_MAX;
             y[j] = -length + 2*length*distribution(generator);
             z[j] = expdistribution(generator);
         }
         fx = int_function(y[0],y[1],y[2],y[3],y[4],y[5]);
         fx_exp = int_function_spherical(z[0],z[1],M_PI*z[2],M_PI*z[3],2*M_PI*z[4],2*M_PI*z[5]);;
         MCint += fx_exp;
         MCintsqr2 += fx_exp*fx_exp;
     }
     MCint = jacobidet*MCint/((double) N);
     MCintsqr2 = MCintsqr2/((double) N);
     double variance = MCintsqr2 - MCint*MCint;


     // FINAL OUTPUT
     cout << "INPUT:" << endl;
     cout << "N = " << N << endl;
     cout << "a = " << a << " (lower limit)" << endl;
     cout << "b = " << b << " (upper limit)" << endl;

     cout << endl << "RESULTS:" << endl;
     //cout << "Gauss-Legendre "<< "\t" << setprecision(15)  << int_gauss << endl;
     //cout << "Gauss-Laguerre " << "\t" << setprecision(15) << int_gausslag << endl;
     cout << "Monte Carlo " << "\t" << setprecision(15) << MCint << " (variance = " << variance << ")"<< endl;

     cout << endl << "Exact answer " << "\t" << 5*M_PI*M_PI/(16*16) << endl;

     finish = clock(); // final time
     cout << endl << "Total time: " << "\t" << ((finish - start)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time

     // Clear memory
     //delete [] x;
     //delete [] w;
     //delete [] xgl;
     //delete [] wgl;
     delete [] y;
     delete [] z;
     delete [] x;
     delete [] w;

     return 0;
}
// End of main program


// DEFINITION OF FUNCTIONS:

// INTEGRAL TO BE SOLVED
double int_function(double x1, double y1, double z1, double x2, double y2, double z2){
    double alpha = 2.0;
    double exp1 = -2*alpha*sqrt(x1*x1 + y1*y1 + z1*z1);
    double exp2 = -2*alpha*sqrt(x2*x2 + y2*y2 + z2*z2);
    double deno = sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));

    if(deno < pow(10.,-6.)){return 0;}
    else {return exp(exp1+exp2)/deno;}
    }


double int_function_spherical(double r1,double r2,double theta1,double theta2,double phi1,double phi2){
    double alpha = 2.0;
    double numerator = r1*r1*r2*r2;
    double cos_beta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double expr = -2*alpha*(r1+r2);
    double deno = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos_beta);

    if(deno <pow(10.,-6.)){return 0;}
    else {return numerator*exp(expr)/deno;}
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



double gammln( double xx)
{
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
