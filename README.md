# Amplitudemaker



This program calculates the amplitude spectrum of your fitness landscape
and is able to fit a generalized LK-landscape to the experimental data

Please prepare the data in the form 1   0   0   ... 1   2.400049
with the binary genome  1   0   0   ... 1 and fitness 2.400049.


It is possible to call the program with parameters landscape-filename and
genome-length to simplify the work with shell scripts. Example: to calculate
the amplitude spectrum of landscape "foo.dat" with genome length "bar"
you can call ./amplitudemaker foo.dat bar

For details on amplitude spectra, see
http://arxiv.org/pdf/1301.1923


For the calculations and fitting, the GSL libraries are used: http://www.gnu.org/software/gsl/
For plotting, gnuplot is used via system call: http://www.gnuplot.info/

I/O is not implemented very elegantly, but it works.
