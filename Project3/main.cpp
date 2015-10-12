//   This is a simple program which tests Gaussian quadrature using Legendre and Laguerre polynomials

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
//#include <armadillo>
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10

using namespace std;


// FUNCTIONS

double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double gammln(double);


int main()
{
     int N = 20;
     double a = -5.0;
     double b = 5.0;


     // GAUSS-LEGENDRE
     // Mesh points weights and function values
     double *x = new double [N];
     double *w = new double [N];

     gauleg(a,b,x,w,N); // mesh points and weights

     // Evaluate the integral with the Gauss-Legendre method
     double int_gauss = 0.; // initialize the sum

     for (int i=0;i<N;i++){
                  for (int j = 0;j<N;j++){
                  for (int k = 0;k<N;k++){
                  for (int l = 0;l<N;l++){
                  for (int m = 0;m<N;m++){
                  for (int n = 0;n<N;n++){
                      int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
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

     for (int i=0;i<N;i++){
                  for (int j = 0;j<N;j++){
                  for (int k = 0;k<N;k++){
                  for (int l = 0;l<N;l++){
                  for (int m = 0;m<N;m++){
                  for (int n = 0;n<N;n++){
                      int_gausslag += wgl[i]*wgl[j]*wgl[k]*wgl[l]*wgl[m]*wgl[n]*int_function(xgl[i],xgl[j],xgl[k],xgl[l],xgl[m],xgl[n]);
                  }}}}}
             }


    /*
     // MONTE-CARLO
     double MCint, MCintsqr2, fx;
     MCint = MCintsqr2 = 0.;
     double invers_period = 1./RAND_MAX; // initialise the random number generator
     srand(time(NULL)); // This produces the so-called seed in MC jargon

     // Evaluate the integral with the a crude Monte-Carlo method
     for ( int i = 1; i <= n; i++){
         double x = double(rand())*invers_period;
         fx = int_function(x,alpha);
         MCint += fx;
         MCintsqr2 += fx*fx;
     }

     MCint = MCint/((double) n );
     MCintsqr2 = MCintsqr2/((double) n );
     double variance = MCintsqr2 - MCint*MCint;
     */


     // FINAL OUTPUT
     cout << "INPUT:" << endl;
     cout << "N = " << N << endl;
     cout << "a = " << a << " (lower limit)" << endl;
     cout << "b = " << b << " (upper limit)" << endl;


     cout << endl << "RESULTS:" << endl;
     cout << "Gauss-Legendre "<< "\t" << setprecision(15)  << int_gauss << endl;
     //cout << "Gauss-Legendre impr. " << "\t" << setprecision(15) << int_gausslegimproved << endl;
     cout << "Gauss-Laguerre " << "\t" << setprecision(15) << int_gausslag << endl;
     //cout << "Monte Carlo " << "\t" << MCint << " (variance = " << variance << ")"<< endl;

     cout << endl << "Exact answer " << "\t" << 5*M_PI*M_PI/(16*16) << endl;

     // Clear memory
     delete [] x;
     delete [] w;
     delete [] xgl;
     delete [] wgl;
     return 0;
}
// End of main program


// DEFINITION OF FUNCTIONS:

// INTEGRAL TO BE SOLVED
double int_function(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double alpha = 2.;
    double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));

    if(deno <pow(10.,-6.)){return 0;}
    else {return exp(exp1+exp2)/deno;}
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
