#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;
#include "Matrix.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "Hamiltonian.h"


#include "random"


int main(int argc, char *argv[]) {

    string ex_string_original =argv[0];

    string ex_string;
    //ex_string.substr(ex_string_original.length()-5);
    ex_string=ex_string_original.substr (2);
    cout<<"'"<<ex_string<<"'"<<endl;



    if(ex_string=="MoireBands"){
        string model_inputfile = argv[1];

        if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }

        Parameters Parameters_;
        Parameters_.Initialize(model_inputfile);

        Coordinates Coordinates_(Parameters_.Grid_L1, Parameters_.Grid_L2, 2);
        Hamiltonian Hamiltonian_(Parameters_, Coordinates_);


        string file_bands_out="Bands_energy.txt";
        ofstream FileBandsOut(file_bands_out.c_str());
        FileBandsOut<<"#index kx_value ky_value E0(k)  E1(k)   E2(k) ....."<<endl;


        int n1, n2;
        Mat_1_intpair k_path;
        k_path.clear();
        pair_int temp_pair;
        int L1_,L2_;
        L1_=Parameters_.BZ_L1;
        L2_=Parameters_.BZ_L2;

        //K+' to K-
        n1=int(2*L1_/3);
        n2=int(L2_/3);
        while(n2>=int(-L2_/3)){
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2--;
            n1=int(2*L1_/L2_)*n2;
        }

        //K- to K+
        n1=int(-2*L1_/3);
        n2=int(-L2_/3);
        n2--;n1++;
        while(n1<=int(-L1_/3)){
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2--;
            n1++;
        }

        //K+ to K+'
        n1=int(-L1_/3);
        n2=int(-2*L2_/3);
        n2=n2+2; //in principle it should be n2=n2+1, n1=n1+1
        n1=n1+2;
        while(n1<=int(2*L1_/3)){
            temp_pair.first = n1;
            temp_pair.second = n2;
            k_path.push_back(temp_pair);
            n2=n2+2;  //in principle it should be n2=n2+1, n1=n1+1
            n1=n1+2;
        }



        for(int index=0;index<k_path.size();index++){
        n1=k_path[index].first;
        n2=k_path[index].second;


        Hamiltonian_.kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3)*L1_))  +  n2*(1.0/(sqrt(3)*L2_)));
        Hamiltonian_.ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
        Hamiltonian_.HTBCreate();
        char Dflag='V';
        Hamiltonian_.Diagonalize(Dflag);

        cout << k_path.size()<<"  "<<index <<"  "<<Hamiltonian_.Ham_.n_col()<<endl;
        FileBandsOut<<index<<"  "<<Hamiltonian_.kx_<<"  "<<Hamiltonian_.ky_<<"   ";
        for(int band=0;band<Hamiltonian_.Ham_.n_col();band++){
            FileBandsOut<<Hamiltonian_.eigs_[band]<<"  ";
        }
        FileBandsOut<<endl;
        }



        Hamiltonian_.Get_Wannier_function(0);

    }






    cout << "--------THE END--------" << endl;
} // main
