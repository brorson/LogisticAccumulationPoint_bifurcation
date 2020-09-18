#include <algorithm> 
#include "globals.h"
#include "f.h"
#include "J.h"


#define PREC 400

mpreal tol = 1e-100;

void bifurcation_calc(int, mpreal&, mpreal&);



//-------------------------------------------------
int main(int argc, char* argv[])
{

  mpreal::set_default_prec(PREC);
  mpreal x0;
  mpreal lam;
  int R;
  int ord;

  // Initialize calc
  ord = 4;
  x0 = 8.0e-01;
  lam = 3.4;


  // Open output file.  Give it name related to requested PREC.
  string filenamebase = "results_bifurcation_prec";
  int iprec = PREC;
  string sprec= to_string(iprec);
  string filename = filenamebase + sprec + ".txt";
  cout << "Opening filename = " << filename << endl;
  std::ofstream outfile;
  outfile.open (filename.c_str(), std::ofstream::out);

  // Set precision for mpfr output.  Divide by 3 since PREC is in bits
  // but I want decimal digits.
  outfile.precision(PREC/3);
  
  outfile << "PREC = " << PREC << ", tol = " << tol << "\n";

  
  // This just loops through each bifurcation point.
  for (int i=0; i<21; i++) {
    printf("=============================================\n");
    mpfr_printf("Calling with order = %d, lam = %10Rf, x0 = %10Rf \n", ord, lam, x0);
    R = ord+1;
    bifurcation_calc(R, lam, x0);
    mpfr_printf ("Computed order = %d, lam = %10Re, x0 = %10Rf \n", ord, lam, x0);

    outfile << "order = " << ord << ", lambda = " << lam << endl;
    outfile.flush();
    
    // Go to next bifurcation order
    ord = 2*ord;
  }

  // All done.
  outfile.close();
  
}


//-------------------------------------------------
void bifurcation_calc(int R, mpreal& lam, mpreal& x0) {
  mpreal::set_default_prec(PREC);  

  int i, j;
  Matrix<mpreal, Dynamic, 1> fn(R,1);
  DynamicSparseMatrix<mpreal> Jn(R, R);
  Matrix<mpreal, Dynamic, 1> delta(R,1);
  Matrix<mpreal, Dynamic, 1> xn(R, 1);
  SparseLU<SparseMatrix<mpreal>, COLAMDOrdering<int> >  solver;

  // Initialize J
  initJ(R, Jn);
  
  // Iterate create starting vector for Newton's method.
  mpreal ex = 3.3*log10(mpreal(R-1));   // log base 2 of R-1
  mpreal base = 0.21;
  mpreal dlam = pow(base, ex);
  lam = lam + dlam;
  printf("dlam = \n");
  mpfr_printf("%10Re \n", dlam);
  printf("starting lam = \n");
  mpfr_printf("%10Re \n", lam);
  for (i=0; i<10*R; i++) {
    x0 = lam*x0*(1-x0);
  }

  // Now iterate R-1 times more and collect x values
  for (i=0; i<R-1; i++) {
    x0 = lam*x0*(1-x0);
    xn(i) = x0;
  }
  xn(R-1) = lam;

  
// Turn this on for debugging only.  Currently is broken. 
#if 0
  // Now compute the difference between elements in this vector
  // and print out the smallest difference
  // Must do this operation on a copy
  Matrix<mpreal, Dynamic, 1> x1(R,1);
  std::copy(xn.begin(), xn.end()-1, std::begin(x1));
  std::sort(x1.begin(), x1.end());
  // Just do diff in place to save memory
  for (i=0; i<R-1; i++) {
    x1[i] = abs(x1[i]-x1[i+1]);
  }
  mpreal diff = 1.0;
  for (i=0; i<R-1; i++) {
    if (x1[i]<diff) {
      diff = x1[i];
    }
  }
  mpfr_printf("Minimum difference between zeros = %10Re \n", diff);  
#endif
  
  
  //--------------------------------------
  // Run a Newton's method loop
  for (int cnt=0; cnt<75; cnt++) {

    //printf("------------------------------------\n");
    
    // Compute f
    f(R, xn, fn);

    // Compute J
    J(R, xn, Jn);
    
    // Now compute delta = Newton step
    if (cnt == 0) {
      // configure solver on first iteration.
      solver.analyzePattern(Jn); 
    }
    solver.factorize(Jn);
    delta = solver.solve(fn); 

    // Take step
    for (i=0; i<R; i++) {    
      xn(i) = xn(i) - delta(i);
    }

    // Check if we're close enough to quit
    //mpfr_printf("norm(delta) = %10Re\n", delta.norm());
    if (delta.norm() < tol) {
      printf("  ---  Terminating after %d iterations because norm(delta) < tol.\n", cnt);
      lam = xn(R-1);
      x0 = xn(0);
      return;
    } else {
      //printf("Must take another loop iteration.\n");
    }
    
  } // for (int cnt=0; cnt<25; cnt++)

  printf("Convergence failed!  Throwing exception.\n");
  throw;

}
