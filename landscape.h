/*This class introduces a method to work with experimental fitness landscapes
given in a data file in the format
sequence          fitness  
1   0   1   ...   1.987
and brings the adjacency matrix as well as other helpful properties

(cc) Johannes Neidhart, 2012*/

#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

using namespace std;

class Clandscape {
    private: 
        vector < vector < int> > linekey;
        int n, L;
        string name;
        int delta(int, int);
        int hamming( vector < int >, vector < int > );
    public: 
        /*Clandscape gets a vector containing the fitness values, in the order of the
        linenumbers of the data file, a vector of vectors, containing the correspondance
        of the line number and the sequence. Sequence length and number of sequences are
        also given.
        */

        double binomial (int, int);
        Clandscape (string, int , int);
        int distance (int, int);
        vector < double > fitness;
        void fitit( vector < double >, vector < double >, vector < double > &);
        void fit(vector <double>, vector < double > &, vector < double > &, double);
        vector < vector < int > > adjacency;
        vector < double > spectrum();
        vector < vector < int > > hypercube;

};

#endif
