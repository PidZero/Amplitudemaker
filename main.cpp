/*
This program can read in data of fitness landscapes and
gives the normalised amplitude spectrum.

compile with g++ main.cpp landscape.cpp -lgsl -lblas -lgslcblas

Needs GSL for calculating the spectrum and fitting and Gnuplot for plotting

(cc) Johannes Neidhart, 2015
*/

#include <iostream>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "landscape.h"
#include <string.h>

using namespace std;
void save_vector(vector < double > vec, string s){
    stringstream strstr;
    string str;
    strstr.str("");
    strstr<<s;
    str=strstr.str();
    ofstream myfile;
    myfile.open(str.c_str());

    for(size_t i = 0; i < vec.size(); i++){
        myfile << i+1 << "\t\t" << vec.at(i) << endl;
    }
    myfile.close();
}
void save_vector(vector < double > vec, vector < double > vec2, vector < double > vec3, string s){
    stringstream strstr;
    string str;
    strstr.str("");
    strstr<<s;
    str=strstr.str();
    ofstream myfile;
    myfile.open(str.c_str());

    for(size_t i = 0; i < vec.size(); i++){
        myfile << i+1 << "\t\t" << vec.at(i) << "\t\t" << vec2.at(i)<<  "\t\t" << vec3.at(i)<< endl;
    }
    myfile.close();
}
void interactive(){
    system("clear");
    cout<<endl;
    cout<<endl;
    cout<<"\t\t\t  oooooooooooooooooooooooo"<<endl;
    cout<<"\t\t\t |   Amplitude Maker 0.6  |"<<endl;
    cout<<"\t\t\t | (cc) Johannes Neidhart |"<<endl;
    cout<<"\t\t\t |           2015         |"<<endl;
    cout<<"\t\t\t  ooooooooooooooooooooooo"<<endl;
    cout<<endl;
    cout<<endl;
    cout<< "This program calculates the amplitude spectrum of your fitness landscape"<<endl;
    cout<< "and is able to fit a generalized LK-landscape to the experimental data"<<endl;
    cout<<endl;
    cout<<"Please prepare the data in the form 1   0   0   ... 1   2.400049"<<endl;
    cout<<"with the binary genome  1   0   0   ... 1 and fitness 2.400049."<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"It is possible to call the programm with parameters landscapefilename and"<<endl;
    cout<<"genomelength to simplify the work with shell scripts. Example: to calculate"<<endl;
    cout<<"the amplitude spectrum of landscape \"foo.dat\" with genome length \"bar\""<<endl;
    cout<<"you can call ./amplitudemaker foo.dat bar"<<endl;
    cout<<endl;

    //name and location of the landscape data file
    cout<<"Please enter the name of the data file: ";
    
    string name;
    cin>>name;
    cout<<endl;

    int L;
    cout<<"Please enter the genome length: L = ";
    cin>>L;
    cout<<endl;
    if(cin.fail()){
        cout<<"L has to be an integer."<<endl;
        exit(0);
    }
   
    //define the number of sequences
    int n=pow(2,L);
    

    Clandscape ls(name, n, L);

    vector < double > spec (L);
    spec = ls.spectrum();
    cout<<"The normalized amplitude Spectrum is:"<<endl<<endl;
    vector < double > plotspec;
    for(int i=1;i<L+1;i++){
        cout<< "p = " << i << "\t B_p = " << spec.at(i)<<endl;
        plotspec.push_back(spec.at(i));
    }
    save_vector(plotspec, "tmp.dat");

    cout<<endl;
    cout<<"What would you like to do?"<<endl;
    cout<<"p: plot only;\ts: save as..;\tf: fit;\tc: change fitting tolerance\tq: Quit"<<endl;

    char hin;
    bool repeat = true;
    bool fit_yet = false;
    stringstream strstr;
    double tol = 0;
    vector < double > fit;
    vector < double > A;
    string h_str;
   
    while(repeat){  
        cin>>h_str;
        if(h_str.size()==1){
            hin = h_str.at(0);
            if(hin == 'q' or hin == 's' or hin =='c' or hin == 'f' or hin == 'p' or hin == '\n'){
                switch(hin){
                    case 'q': 
                    {
                        repeat = false;
                        break;
                    }
                    case 'p':
                    {
                        if(fit_yet){
                            system("gnuplot -e \"set terminal wxt persist; \
                            set logscale y; set xrange [1:*]; \
                            set xlabel 'p'; set ylabel 'B_p'; \
                            pl './tmp.dat' u 1:2 title 'Spectrum', 'tmp_fit.dat' u 1:2 w lines title 'Fit'\" & > /dev/null"); 
                        }else{
                            system("gnuplot -e \"set terminal wxt persist; \
                                    set logscale y; \
                                    set xrange [1:*]; \
                                    set xlabel 'p'; \
                                    set ylabel 'B_p'; \
                                    pl './tmp.dat' u 1:2\" & > /dev/null"); 
                        }
                        break;
                    }
                    case 's':
                    {
                        cout<<"Where to save? ";
                        cin>>name;
                        if(fit_yet){
                            save_vector(plotspec, fit, A, name);
                        }else{
                            save_vector(plotspec, name);
                        }
                        break;
                    }
                    case 'f':
                    {
                        fit.resize(spec.size()-1, 0);
                        A.resize(spec.size());
                        ls.fit(spec, A, fit, tol );

                        save_vector(A, "tmp_A.dat");
                        save_vector(fit, "tmp_fit.dat");

                        cout<<endl;
                        cout<<endl;
                        cout<<endl;
                        cout<<"The normalized amplitude Spectrum, fit and A_i are:"<<endl;
                        for(int i=0;i<L;i++){
                            cout<<"p = "<<i+1<<"\tB_p: "<< plotspec.at(i) << "\tfit: "<< fit.at(i) <<"\tA_i " << A.at(i)<<endl;
                        }
                        fit_yet = true;
                        break;
                    }
                    case 'c':
                    {
                        cout<<endl;
                        cout<<"Standard tolerance is A_i >= 0."<<endl;
                        cout<<"New tolerance: ";
                        cin>>tol;
                        cout<<tol;
                        cout<<endl;
                    }
                }
                }else{
                    cout<<endl;
                    cout<<hin<<" is not a valid option."<<endl;
                    cout<<"What would you like to do?"<<endl;
                    cout<<"p: plot only;\ts: save as..;\tf: fit;\tc: change fitting tolerance\tq: Quit"<<endl;
                }
            }else{
                cout<<endl;
                cout<<"Please enter one character only!"<<endl;
                cout<<endl;
            }
                
    } 
    system("clear");
    system("rm tmp.dat");
    if(fit_yet){
        system("rm tmp_A.dat");
        system("rm tmp_fit.dat");
    }

}


int main(int argc, char* argv[]){
    if(argc == 1){
        interactive();
    }else{
        if(argc == 3){
            string name;
            name = argv[1];
            string L_str;
            L_str = argv[2];
            bool isint=true;
            for(size_t i = 0; i < L_str.size(); i++){
                if(not isdigit(L_str.at(i))){
                    isint = false;
                    break;
                }
            }
            int L;
            if(isint){
                L = atoi(argv[2]);
                int n=pow(2,L);

                Clandscape ls(name, n, L);

                vector < double > spec (L);
                spec = ls.spectrum();
                cout<<"The normalized amplitude Spectrum is:"<<endl<<endl;
                for(int i=1;i<L+1;i++){
                    cout<< "p = " << i << "\t B_p = " << spec.at(i)<<endl;
                }
            }else{
                cout<<"The second parameter has to be integer."<<endl;
                exit(0);
            }
        }else{
            cout<<"Wrong number of parameters."<<endl;
        }
    }
    return(0);
}
