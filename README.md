# Amplitudemaker



This program calculates the amplitude spectrum of your fitness landscape
and is able to fit a generalized LK-landscape to the experimental data

Please prepare the data in the form 1   0   0   ... 1   2.400049
with the binary genome  1   0   0   ... 1 and fitness 2.400049.


It is possible to call the programm with parameters landscapefilename and
genomelength to simplify the work with shell scripts. Example: to calculate
the amplitude spectrum of landscape "foo.dat" with genome length "bar"
you can call ./amplitudemaker foo.dat bar

For details on aplitude spectra, see
http://arxiv.org/pdf/1301.1923
