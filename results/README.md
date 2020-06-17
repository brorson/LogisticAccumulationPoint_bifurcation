A calculation of the logistic map's accumulation point lambda_inf
=================================================================

The code in the parent directory implements a high-precision calculation of
the accumulation point of the logistic map using C++.  The C++ code
computes the accumulation point as follows:

1.  Using a starting point in the period 4 window, the right side of
the window (i.e. the 4 -> 8 bifurcation point) is computed using
Newton's method.  The system to be solved in covered in 
["Computing the bifurcation points and superstable orbits of the logistic map", 
I. S. Kotsireas1 and K. Karamanos2](https://carma.newcastle.edu.au/resources/jon/Preprints/Books/EMA/Exercises/For%20others/dhb-logistics%20B5.pdf)
specifically eqs (10) and (11).  Solving this system gives
the value of lambda_4, which is dumped into a file.

Note that eq (11) implies two possible values for lambda, one on the
right and one on the left of the window depending whether one takes +1
or -1 as the RHS.  This work attempts to solve for the right hand
side of the window, so this equation is modified appropriately.

2.  Then the program steps a little bit to the right from lambda4,
into the period 8 window.  It then runs the logistic iteration a
number of times to generate a good starting point for Newton's
method. 

3.  Newton's method is again used to compute the right side of the
window, giving the value of lambda_8, which is stored into the same
file. 

4.  The program runs steps 2 and 3 repeatedly, walking from one period
2^N window to the next, storing the values of lambda calculated at
each step.

5.  When the computation finishes, the file holds high-precision
values of lambda_4, lambda_8, lambda_16, ....  These values converge on
lambda_inf, the value of lambda at the accumulation point.  However,
convergence is slow.  Therefore, I use a Python program to read in the
lambda sequence, and apply Aitken's accleration method to squeeze as
many digits as I can out of the converging sequence.

---

This directory contains the results obtained so far.  The results are
stored in files named results_<precision>_<tolerance exp>.txt.  This
directory also contains the Python program compute_laminf.py.  The
Python program serves two purposes:

1.  It reads the lambda_n sequences from two files and enables me to
compare the two results.  The goal is to check that the lambda
sequences are the same to some number of digits when calculated using
different precision and tolerance settings.  This serves as a check
that the computation returns sane and hopefully correct results.

2.  The Python program also implements Aitken's method to accelerate
convergence of a sequence.

I generally hand-edit the Python program to choose which files to
read, which analysis to run, etc.

---

As of now, I match the value of lambda_inf from the OEIS up to
O(1e-27).  My digits which match the OEIS are:

3.56994567187094490184200515

Stuart Brorson 6.17.2020.



