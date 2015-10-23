#include "landscape.h"
Clandscape::Clandscape (string a, int c, int d){

    /*initialise the private variables*/
    name = a;
    n=c;
    L=d;


    /*initialize the landscape and
    read out the data file*/
    fitness.clear();
    linekey.clear();
    fitness.resize(n);
    linekey.resize(n);
    string line;
    ifstream myfile; 
    myfile.exceptions( ifstream::failbit | ifstream::badbit );
    try{
        myfile.open(name.c_str());
    }
    catch(ifstream::failure e){
        cerr<<endl<<"ERROR: Could not open "<<name.c_str()<<"; No such file."<<endl<<endl;
        myfile.close();
    }
    double data;
    for(int i=0;i<n;i++){
        getline(myfile,line);
        istringstream is(line);
        for(int j=0; j<L+1; j++){
            is>>data;
            if(j<L){
                linekey.at(i).push_back((int)data);
            }else{
                fitness.at(i)=data;
            }
        }
    }
    myfile.close();


    /*calculate adjacency matrix and store it in adjacency*/

    int h;
    adjacency.resize(n);

    for(int i = 0; i<n; i++){
        for(int j=0; j<n; j++){
            h=0;
            for(int k=0; k<L; k++){
                h+=1-delta(linekey.at(i).at(k),linekey.at(j).at(k));
            }
            if(h==1){
                adjacency.at(i).push_back(1);
            }else{
                adjacency.at(i).push_back(0);
            }
        }
    }
}
int Clandscape::delta(int a, int b) {
    if(a==b){return(1);}
    else{return(0);}
}

int Clandscape::distance(int a, int b) {
    int h=0;
    for(int i=0;i<L;i++){
        h+=1-delta(linekey.at(a).at(i),linekey.at(b).at(i));
    }
    return(h);
}

double Clandscape::binomial(int n, int k){
    if(k>n){return(0);}else{
    if(k==0){return(1);}else{
    if(k==n){return(1);}else{
    double h=1;
    for(double i=1;i<=(double)k;i++){
        h*=(double(n-(k-i))/i);
    }
    return(h);
    }}}
}

void Clandscape::fit(vector <double> spec, vector < double > &A_o , vector < double > &fit, double tol){
    vector < double > sp;
    vector < double > picks;
    A_o.clear();
    A_o.resize(spec.size()-1, 0);
    bool repeat;
    repeat = true;
    for(size_t i = 1; i < spec.size(); i++){
        sp.push_back(spec.at(i));
        picks.push_back(i);
    }
    size_t p_max = sp.size();
    double chisq;
    double tol_q = tol / 20.;
    while(repeat){

        fit.clear();
        fit.resize(spec.size()-1, 0);
        gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(p_max, picks.size());
        gsl_matrix *X = gsl_matrix_alloc(p_max, picks.size());
        gsl_vector *A = gsl_vector_calloc(picks.size());
        gsl_vector *s = gsl_vector_alloc(p_max);
        gsl_matrix *cov = gsl_matrix_alloc(picks.size(), picks.size());

        repeat = false;
        for(size_t i = 0; i < p_max; i++){
            int p = i+1;
            for(size_t j = 0; j < picks.size(); j++){
                gsl_matrix_set(X, i, j,
                    pow(2., -picks.at(j))*binomial(picks.at(j), p)
                    );
            }
            gsl_vector_set(s, i, sp.at(i));
        }
        gsl_multifit_linear(X, s, A, cov, &chisq, work);
        for(size_t i = 0; i<picks.size(); i++){
            if(gsl_vector_get(A, i) < tol_q){
                repeat = true;
                picks.erase(picks.begin()+i);
            }
        }
        if(tol_q < tol){
            tol_q += tol / 20.;
            repeat = true;
        }
        if(repeat){
            gsl_matrix_free(X);
            gsl_matrix_free(cov);
            gsl_vector_free(A);
            gsl_vector_free(s);
            gsl_multifit_linear_free(work);
        }
        if(not repeat){
            for(size_t i = 0; i < p_max; i++){
                int p = i+1;
                for(size_t j = 0; j < picks.size(); j++){
                    fit.at(i) +=  pow(2., -(double)picks.at(j))*binomial(picks.at(j), p)*gsl_vector_get(A, j);
                }
            }
            for(size_t i = 0; i<picks.size(); i++){
                A_o.at(picks.at(i)-1) = gsl_vector_get(A, i);
            }
            gsl_matrix_free(X);
            gsl_matrix_free(cov);
            gsl_vector_free(A);
            gsl_vector_free(s);
            gsl_multifit_linear_free(work);
        }

    }
}


//Calculate the spectrum
vector < double > Clandscape::spectrum(){
    gsl_vector *fitvec=gsl_vector_alloc(n);
    gsl_vector *specgsl = gsl_vector_alloc(n);
    vector < double > spec;
    vector < double > specout;
    vector < vector < double > > evalkey; //stores at 0 the 
    evalkey.resize(L); //eigenvalue and then the positions where the 
    //same are 
    double epsilon = 0.00001;  //precision of eigenvalue comparison
    double sum;
    gsl_matrix *m=gsl_matrix_alloc(n,n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            gsl_matrix_set(m,i,j,+(double)adjacency[i][j]);
            //cout<<i<<"\t"<<j<<"\t"<<adjacency.at(i).at(j)<<endl;
        }
        gsl_matrix_set(m,i,i,(double)adjacency[i][i]-(double)L);
        gsl_vector_set(fitvec,i,fitness.at(i));
    }
/*      
    Calculate the eigenvectors with the GSL
    */
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
                
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
    //Calculate the eigenvectors and eigenvalues: 
    gsl_eigen_symmv (m, eval, evec, w);
    //Sort eigenvectors and eigenvalues simultaneously:
    gsl_eigen_symmv_sort (eval, evec, 
                                    GSL_EIGEN_SORT_VAL_DESC);
    /*
    evec contains the eigenvectors as columns, so evec is transposed, multyplied 
    by 1.0 and fitvec, while 0.0 is added. The result is written in specgsl.
    */
    gsl_eigen_symmv_free (w);
    gsl_blas_dgemv( CblasTrans, 1.0, evec, fitvec, 0.0, specgsl );
    spec.clear();
    for(int i=0; i<n; i++){
        spec.push_back(gsl_vector_get(specgsl,i));
    }
    sum=0;
    for(int i=1; i<n;i++){
        sum+=spec.at(i)*spec.at(i);
    }
    int h=0;
    int i=0;
    int j;
    specout.resize(L+1,0);
    while(h<L+1){
        j=i;
        while(i<n 
            and gsl_vector_get(eval,i)>(gsl_vector_get(eval,j)-epsilon)
            and gsl_vector_get(eval,i)<(gsl_vector_get(eval,j)+epsilon)){
            specout.at(h)+=spec.at(i)*spec.at(i)/sum;

            i++;
        }
        h++;
    }
    gsl_matrix_free(m);
    gsl_matrix_free(evec);
    gsl_vector_free(fitvec);
    gsl_vector_free(eval);
    gsl_vector_free(specgsl);
    
    return(specout);
}

