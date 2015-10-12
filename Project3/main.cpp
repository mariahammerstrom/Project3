//   This is a simple program which tests Gaussian quadrature using Legendre and Laguerre polynomials

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <armadillo>
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10

using namespace std;


// FUNCTIONS

double int_function(double x,double alpha);
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double gammln(double);


int main()
{
     int n = 100;
     double a = -50;
     double b = 50;
     double alpha = 2.0;


     // GAUSS-LEGENDRE
     // Mesh points weights and function values
     double *x = new double [n];
     double *w = new double [n];

     gauleg(a,b,x,w,n); // mesh points and weights

     // Evaluate the integral with the Gauss-Legendre method
     double int_gauss = 0.; // initialize the sum

     for ( int i = 0;  i < n; i++){
        int_gauss+=w[i]*int_function(x[i],alpha);
     }

     // Improved Gauss-Legendre, mapping of x \in [-1,1] to x \in [0, infinity)
     double *r = new double [n];
     double *s = new double [n];

     // Evaluate the integral with the improved Gauss-Legendre method. Mesh points with a tangent mapping.
     gauleg(-1.0, 1.0,x,w, n);
     double pi_4 = acos(-1.0)*0.25;

     for ( int i = 0;  i < n; i++){
       double xx=pi_4*(x[i]+1.0);
       r[i]= tan(xx);
       s[i]=pi_4/(cos(xx)*cos(xx))*w[i];
     }

     double int_gausslegimproved = 0.;
     for ( int i = 0;  i < n; i++){
       int_gausslegimproved += s[i]*int_function(r[i],alpha);
     }


     // GAUSS-LAGUERRE
     // Mesh points weights and function values
     double *xgl = new double [n+1];
     double *wgl = new double [n+1];

     double alf = 1.0;
     gauss_laguerre(xgl,wgl,n,alf);

     // Evaluate the integral with the Gauss-Laguerre method
     double int_gausslag = 0.; // initialize the sum
     for ( int i = 1;  i <= n; i++){
       int_gausslag += wgl[i];//*sin(xgl[i]);
     }



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



     // FINAL OUTPUT
     cout  << setiosflags(ios::showpoint | ios::uppercase);
     cout << "Gaussian-Legendre quad = "<< setw(20) << setprecision(15)  << int_gauss << endl;
     cout << "Gaussian-Legendre improved quad = " << setw(20) << setprecision(15) << int_gausslegimproved << endl;
     cout << "Gaussian-Laguerre quad = " << setw(20) << setprecision(15) << int_gausslag << endl;
     cout << "Correct answer = " << 5*3.14*3.14/(16*16) << endl;

     cout << "MC Integral = " << MCint << " Exact = " << M_PI << " Variance = " << variance << endl;

     // Clear memory
     delete [] x;
     delete [] w;
     delete [] xgl;
     delete [] wgl;
     delete [] s;
     delete [] r;
     return 0;
}
// End of main program


//  DEFINITION OF FUNCTIONS


// INTEGRAL TO BE SOLVED
double int_function(double x,double alpha)
{
  double value = exp(-2*alpha*x);
  return value;
}



/* GAUSS-LEGENDRE
 * The function gauleg() takes the lower and upper limits of integration x1, x2,
 * calculates and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
 * of length n of the Gauss--Legendre n--point quadrature formulae.
*/

void gauleg(double x1, double x2, double x[], double w[], int n)
{
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
