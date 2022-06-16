#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);
//zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);

class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__, Coordinates&  Coordinates__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__)

    {
        Initialize();
    }


    void Initialize();    //::DONE
    void Hoppings();        //::DONE
    double GetCLEnergy();    //::DONE
    void InteractionsCreate();   //::DONE
    void Check_Hermiticity();  //::DONE
    void HTBCreate();   //::DONE
    double chemicalpotential(double muin,double Particles);    //::DONE
    void Get_Wannier_function(int band);

    double TotalDensity();   //::DONE
    double E_QM();   //::DONE
    double NIEnergy(double kx_val, double ky_val);

    void Diagonalize(char option);   //::DONE
    void copy_eigs(int i);  //::DONE

    int convert_jm_to_int(string jm_val);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    int ns_, l1_, l2_;
    double kx_, ky_;
    double k_plusx, k_minusx, k_plusy, k_minusy;

    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_,sx_,sy_,sz_;

    //real space  effective H params
    int L1_eff, L2_eff;
    Mat_2_Complex_doub Tij;
    Mat_2_Complex_doub Uij;

};


void Hamiltonian::Initialize(){

    ns_=Parameters_.ns;
    l1_=Parameters_.Grid_L1;
    l2_=Parameters_.Grid_L2;

    int space=ns_*2;

    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    eigs_saved_.resize(space);

    k_plusx = (-1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    k_plusy = (-1.0/3.0)*(2.0*PI/Parameters_.a_moire);
    k_minusx = (-1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    k_minusy = (1.0/3.0)*(2.0*PI/Parameters_.a_moire);


    //real space  effective H params
    L1_eff=12;L2_eff=12;
    Tij.resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Tij[i].resize(L1_eff*L2_eff);
    }
    Uij.resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Uij[i].resize(L1_eff*L2_eff);
    }

} // ----------

double Hamiltonian::TotalDensity(){

    double n1=0.0;
    /*
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    */
    return n1;

} // ----------



double Hamiltonian::E_QM(){

    return 0.0;

} // ----------

double Hamiltonian::NIEnergy(double kx_val, double ky_val){

    double energy_;
    //energy_ = -1.0*(Parameters_.RedPlanckConst*Parameters_.RedPlanckConst*(kx_val*kx_val  + ky_val*ky_val))*(0.5/Parameters_.MStar);
    energy_ = -1.0*(((3.809842*1000)/Parameters_.MStar)*(kx_val*kx_val  + ky_val*ky_val));

    return energy_;
}

double Hamiltonian::GetCLEnergy(){

    return 0.0;

} // ----------


int Hamiltonian::convert_jm_to_int(string jm_val){

    int val;
    if(jm_val=="3by2_m3by2"){val=0;}
    if(jm_val=="3by2_3by2"){val=1;}
    if(jm_val=="3by2_m1by2"){val=2;}
    if(jm_val=="3by2_1by2"){val=3;}
    if(jm_val=="1by2_m1by2"){val=4;}
    if(jm_val=="1by2_1by2"){val=5;}
    return val;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<Ham_.n_row();i++) {
        for(int j=0;j<Ham_.n_row();j++) {
            if(
                    abs(Ham_(i,j) - conj(Ham_(j,i)))>0.00001
                    ) {
                cout<<Ham_(i,j)<<endl;
                cout<<conj(Ham_(j,i))<<endl;

            }
            assert(
                        abs(Ham_(i,j) - conj(Ham_(j,i)))<0.00001
                        ); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}





void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}


void Hamiltonian::HTBCreate(){


    Ham_.resize(ns_*2,ns_*2);
    double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);

    int Bottom_, Top_;
    Bottom_=0;Top_=1;

    //l1_/2,l2_/2 is the k-point

    double kx_local, ky_local;

    int row, col;
    int i1_neigh, i2_neigh;
    for(int i1=0;i1<l1_;i1++){
        for(int i2=0;i2<l2_;i2++){
            kx_local = kx_ + (-(l1_/2)+i1)*(b1x_) + (-(l2_/2)+i2)*(b2x_);
            ky_local = ky_ + (-(l1_/2)+i1)*(b1y_) + (-(l2_/2)+i2)*(b2y_);
            for(int orb=0;orb<2;orb++){
                row=Coordinates_.Nbasis(i1, i2, orb);
                if(orb==Bottom_){

                    //1
                    col = row;
                    Ham_(row,col) += NIEnergy(kx_local - k_plusx, ky_local - k_plusy) + (0.5*Parameters_.Vz_);

                    //2 i.e +/- b1
                    i1_neigh = i1 + 1;
                    i2_neigh = i2;
                    if(i1_neigh<l1_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(-iota_complex*Parameters_.Psi_param);
                    }
                    i1_neigh = i1 - 1;
                    i2_neigh = i2;
                    if(i1_neigh>=0){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(iota_complex*Parameters_.Psi_param);
                    }


                    //3 i.e +/- b3
                    i1_neigh = i1 - 1;
                    i2_neigh = i2 + 1;
                    if(i1_neigh>=0 && i2_neigh<l2_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(-iota_complex*Parameters_.Psi_param);
                    }
                    i1_neigh = i1 + 1;
                    i2_neigh = i2 - 1;
                    if(i1_neigh<l1_ && i2_neigh>=0){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(iota_complex*Parameters_.Psi_param);
                    }


                    //4 i.e +/- b5
                    i1_neigh = i1;
                    i2_neigh = i2 - 1;
                    if(i2_neigh>=0){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(-iota_complex*Parameters_.Psi_param);
                    }
                    i1_neigh = i1;
                    i2_neigh = i2 + 1;
                    if(i2_neigh<l2_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(iota_complex*Parameters_.Psi_param);
                    }


                    //5
                    col = Coordinates_.Nbasis(i1, i2, Top_);
                    Ham_(row,col) += Parameters_.omega_param;

                    //6
                    i1_neigh = i1;
                    i2_neigh = i2+1;
                    if(i2_neigh<l2_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, Top_);
                        Ham_(row,col) += Parameters_.omega_param;
                    }

                    //7
                    i1_neigh = i1-1;
                    i2_neigh = i2+1;
                    if(i1_neigh>=0  && i2_neigh<l2_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, Top_);
                        Ham_(row,col) += Parameters_.omega_param;
                    }


                }
                else{//i.e. orb=Top_


                    //1
                    col = row;
                    Ham_(row,col) += NIEnergy(kx_local - k_minusx, ky_local - k_minusy) - (0.5*Parameters_.Vz_);

                    //2 i.e +/- b1
                    i1_neigh = i1 + 1;
                    i2_neigh = i2;
                    if(i1_neigh<l1_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(iota_complex*Parameters_.Psi_param);
                    }
                    i1_neigh = i1 - 1;
                    i2_neigh = i2;
                    if(i1_neigh>=0){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(-iota_complex*Parameters_.Psi_param);
                    }


                    //3 i.e +/- b3
                    i1_neigh = i1 - 1;
                    i2_neigh = i2 + 1;
                    if(i1_neigh>=0 && i2_neigh<l2_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(iota_complex*Parameters_.Psi_param);
                    }
                    i1_neigh = i1 + 1;
                    i2_neigh = i2 - 1;
                    if(i1_neigh<l1_ && i2_neigh>=0){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(-iota_complex*Parameters_.Psi_param);
                    }


                    //4 i.e +/- b5
                    i1_neigh = i1;
                    i2_neigh = i2 - 1;
                    if(i2_neigh>=0){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(iota_complex*Parameters_.Psi_param);
                    }
                    i1_neigh = i1;
                    i2_neigh = i2 + 1;
                    if(i2_neigh<l2_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, orb);
                        Ham_(row,col) += Parameters_.V_param*exp(-iota_complex*Parameters_.Psi_param);
                    }


                    //5
                    col = Coordinates_.Nbasis(i1, i2, Bottom_);
                    Ham_(row,col) += Parameters_.omega_param;

                    //6
                    i1_neigh = i1;
                    i2_neigh = i2-1;
                    if(i2_neigh>=0){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, Bottom_);
                        Ham_(row,col) += Parameters_.omega_param;
                    }

                    //7
                    i1_neigh = i1+1;
                    i2_neigh = i2-1;
                    if(i2_neigh>=0  && i1_neigh<l1_){//OBC
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, Bottom_);
                        Ham_(row,col) += Parameters_.omega_param;
                    }

                }

            }
        }
    }



} // ----------



void Hamiltonian::Get_Wannier_function(int band){


    double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);

    double kx_local, ky_local;
    complex<double> temp_factor;
    complex<double> checknorm;
    int eigen_no, comp;
    double dis_x,dis_y;


    double rx_min, ry_min, rx_max, ry_max, d_rx, d_ry;
    rx_min=-5.0*Parameters_.a_moire;
    ry_min=-5.0*Parameters_.a_moire;
    rx_max=5.0*Parameters_.a_moire;
    ry_max=5.0*Parameters_.a_moire;
    int r_ind;
    int space_slices=200;
    d_rx=(rx_max-rx_min)/(space_slices);
    d_ry=(ry_max-ry_min)/(space_slices);


    double qx_min, qy_min, qx_max, qy_max, d_qx, d_qy;
    qx_min=-0.5*PI;//Parameters_.a_moire;
    qy_min=-0.5*PI;//Parameters_.a_moire;
    qx_max=0.5*PI;//Parameters_.a_moire;
    qy_max=0.5*PI;//Parameters_.a_moire;
    int q_slices=200;
    d_qx=(qx_max-qx_min)/(q_slices);
    d_qy=(qy_max-qy_min)/(q_slices);
    double eta_q=0.001;
    double d_sqr=360000;
    double screening=0.0;


    double q_max, d_q,  d_theta;
    q_max=0.25*PI;
    d_q=(q_max)/(q_slices);
    int theta_slices=200;
    d_theta = (2.0*PI)/(theta_slices);
   
    int L1_,L2_;
    L1_=Parameters_.BZ_L1;
    L2_=Parameters_.BZ_L2;

    int Bottom_, Top_;
    Bottom_=0;Top_=1;
    Mat_2_Complex_doub Wnr_state_, Psi_state_;
    Wnr_state_.resize(2);Psi_state_.resize(2);
    for(int type=Bottom_;type<=Top_;type++){
        Wnr_state_[type].resize(space_slices*space_slices);
        Psi_state_[type].resize(space_slices*space_slices);
    }


    Mat_2_Complex_doub Mq_;
    Mq_.resize(q_slices);
    for(int q1=0;q1<q_slices;q1++){
        Mq_[q1].resize(q_slices);
    }

    Mat_2_Complex_doub Mq_SphC;
    Mq_SphC.resize(q_slices);
    for(int q_i=0;q_i<q_slices;q_i++){
        Mq_SphC[q_i].resize(theta_slices);
    }



    Mat_2_Complex_doub Vq_;
    Vq_.resize(q_slices);
    for(int q1=0;q1<q_slices;q1++){
        Vq_[q1].resize(q_slices);
    }


    int center_=(L1_eff/2) + L1_eff*(L2_eff/2);
    int center_neigh;
    //Tij[center_][center_p1]=0.0;
    eigen_no=(2*l1_*l2_)-1-band;

     temp_factor=0.0;
    for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){

            cout<<"doing "<<n1<<"  "<<n2<<endl;

            kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3)*L1_))  +  n2*(1.0/(sqrt(3)*L2_)));
            ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
            HTBCreate();
            char Dflag='V';
            Diagonalize(Dflag);

            //------------------



	    //------------------
/*
	    int k_n=n1+L1_*n2;
	   //Print C_{t,b}(k+n1G1 + n2G2) for fixed k, varying n1,n2
	   string file_Ck_out="States_C_k" + to_string(k_n) +".txt";
           ofstream FileCkOut(file_Ck_out.c_str());
           FileCkOut<<"# n1 ,n2(brillioun zone) = "<<n1<<" "<<n2<<endl;
           FileCkOut<<"#eigen_no(Reciprocal space vectors)  Cvalue"<<endl;

	   for(int i1=0;i1<l1_;i1++){ //on Reciprocal lattice
               for(int i2=0;i2<l2_;i2++){
                int comp_n = Coordinates_.Nbasis(i1, i2, 0);	
		FileCkOut<<i1<<"  "<<i2<<"   "<<Ham_(comp_n,eigen_no).real()<<"  "<<Ham_(comp_n,eigen_no).imag()<<endl;
	   }
	FileCkOut<<endl;
	}
*/			
	   //------------------


            //Hopping-real space-------
            for(int r1_=0;r1_<L1_eff;r1_++){
                for(int r2_=0;r2_<L2_eff;r2_++){
                    center_neigh = (r1_) + L1_eff*(r2_);
                    dis_x = ((sqrt(3.0)/2.0)*(r1_-(L1_eff/2)) +  (sqrt(3.0)/2.0)*(r2_-(L2_eff/2)))*Parameters_.a_moire;
                    dis_y = (-0.5*(r1_-(L1_eff/2)) + 0.5*(r2_-(L2_eff/2)))*Parameters_.a_moire;
                    Tij[center_][center_neigh]+=(1.0/(L1_*L2_))*exp(iota_complex*( kx_*(dis_x) +  ky_*(dis_y) ))*eigs_[eigen_no];
                }
            }
            //--------------------

            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    r_ind=rx_ind + (space_slices)*ry_ind;
                    for(int orb=0;orb<2;orb++){
                        Psi_state_[orb][r_ind]=0.0;
                        for(int i1=0;i1<l1_;i1++){ //on Reciprocal lattice
                            for(int i2=0;i2<l2_;i2++){
                                kx_local = kx_ + (-(l1_/2)+i1)*(b1x_) + (-(l2_/2)+i2)*(b2x_);
                                ky_local = ky_ + (-(l1_/2)+i1)*(b1y_) + (-(l2_/2)+i2)*(b2y_);
                                comp = Coordinates_.Nbasis(i1, i2, orb);

                                Psi_state_[orb][r_ind] += Ham_(comp,eigen_no)*exp(iota_complex*( kx_local*(rx_min + rx_ind*d_rx) +  ky_local*(ry_min + ry_ind*d_ry) ));
                            }
                        }

                        if( orb==0 && (abs(rx_min + rx_ind*d_rx)<=0.00000001 &&
                                       abs(ry_min + ry_ind*d_ry)<=0.00000001) ){
                            temp_factor = Psi_state_[orb][r_ind];
                        }
                    }

                }
            }

            //choosing phase
            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    r_ind=rx_ind + (space_slices)*ry_ind;
                    for(int orb=0;orb<2;orb++){
                        Psi_state_[orb][r_ind] *= conj(temp_factor)/abs(temp_factor);
                    }
                }
            }

            //            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
            //                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
            //                    r_ind=rx_ind + (space_slices)*ry_ind;
            //                    for(int orb=0;orb<2;orb++){
            //                     checknorm +=   Psi_state_[orb][r_ind]*conj(Psi_state_[orb][r_ind]);
            //                    }
            //                }
            //            }

            //-------------------


            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    r_ind=rx_ind + (space_slices)*ry_ind;
                    for(int type=Bottom_;type<=Top_;type++){
                        Wnr_state_[type][r_ind] += (1.0/sqrt(L1_*L2_))*(Psi_state_[type][r_ind]);

                    }
                }
            }



        }
    }



    checknorm=0.0;
    for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
        for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
            r_ind=rx_ind + (space_slices)*ry_ind;
            checknorm += d_rx*d_ry*(Wnr_state_[0][r_ind]*conj(Wnr_state_[0][r_ind]) + Wnr_state_[1][r_ind]*conj(Wnr_state_[1][r_ind]));
        }
    }

    for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
        for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
            r_ind=rx_ind + (space_slices)*ry_ind;
            Wnr_state_[0][r_ind] = Wnr_state_[0][r_ind]*(1.0/sqrt(abs(checknorm)));
            Wnr_state_[1][r_ind] = Wnr_state_[1][r_ind]*(1.0/sqrt(abs(checknorm)));
        }
    }


    checknorm=0.0;
    string file_Wnr_out="Wannier_functions_band" + to_string(band) +".txt";
    ofstream FileWNROut(file_Wnr_out.c_str());
    FileWNROut<<"#index rx ry W_bottom W_top  ....."<<endl;
    for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
        for(int ry_ind=0;ry_ind<space_slices;ry_ind++){

            complex<double> exp_kpr, exp_kmr;
            
	    kx_=(2.0*PI/Parameters_.a_moire)*(((-1.0*L1_)/3.0)*(1.0/(sqrt(3)*L1_))  +  ((-2.0*L2_)/3.0)*(1.0/(sqrt(3)*L2_)));
            ky_=(2.0*PI/Parameters_.a_moire)*(((-1.0*L1_)/3.0) *(-1.0/(L1_))  +  ((-2.0*L2_)/3.0)*(1.0/(L2_)));
	    exp_kpr=exp (-1.0*iota_complex*( (kx_*(rx_min + rx_ind*d_rx) ) + (ky_*(ry_min + ry_ind*d_ry) ) ) );


	    kx_=(2.0*PI/Parameters_.a_moire)*(((-2.0*L1_)/3.0)*(1.0/(sqrt(3)*L1_))  +  ((-1.0*L2_)/3.0)*(1.0/(sqrt(3)*L2_)));
            ky_=(2.0*PI/Parameters_.a_moire)*(((-2.0*L1_)/3.0) *(-1.0/(L1_))  +  ((-1.0*L2_)/3.0)*(1.0/(L2_)));
            exp_kmr=exp (-1.0*iota_complex*( (kx_*(rx_min + rx_ind*d_rx) ) + (ky_*(ry_min + ry_ind*d_ry) )  ));



            r_ind=rx_ind + (space_slices)*ry_ind;
            FileWNROut<<r_ind<<"  "<<(rx_min + rx_ind*d_rx)/(Parameters_.a_moire)<<"   "<<(ry_min + ry_ind*d_ry)/(Parameters_.a_moire)<<"   "<<abs(Wnr_state_[0][r_ind])<<"   "<<abs(Wnr_state_[1][r_ind])<< "   "<<arg((exp_kpr*Wnr_state_[0][r_ind])/abs(Wnr_state_[0][r_ind]))<<"   "<< arg((exp_kmr*Wnr_state_[1][r_ind])/abs(Wnr_state_[1][r_ind]))<<endl;
            checknorm += d_rx*d_ry*(Wnr_state_[0][r_ind]*conj(Wnr_state_[0][r_ind]) + Wnr_state_[1][r_ind]*conj(Wnr_state_[1][r_ind]));
        }
        FileWNROut<<endl;

    }

    FileWNROut<<"#a_moire = "<<Parameters_.a_moire<<endl;
    FileWNROut<<"#norm = "<<checknorm.real()<<"  "<<checknorm.imag()<<endl;


    //cout<<"Tij[0][0+a1] = "<<Tij[center_][center_p1]<<endl;

    cout<<"--------------------------------------Tij[center][neigh]--------------------------------------"<<endl;
    for(int r2_=0;r2_<L2_eff;r2_++){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Tij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;





    //For Uij--------------------------------------------------------------------//

/*
  double q_max, d_q,  d_theta;
    q_max=2.0*PI;
    d_q=(q_max)/(q_slices);
    int theta_slices=200;
    d_theta = (2.0*PI)/(theta_slices);
*/
 
double q_val, theta_val;    
double qx_,qy_;
for(int q_ind=0;q_ind<q_slices;q_ind++){
 q_val = q_ind*d_q;
for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
 theta_val = theta_ind*d_theta;

qx_=q_val*cos(theta_val);
qy_=q_val*sin(theta_val);


      Mq_SphC[q_ind][theta_ind]=0.0;
      for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
      for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
      r_ind=rx_ind + (space_slices)*ry_ind;
      Mq_SphC[q_ind][theta_ind] += d_rx*d_ry*(
      Wnr_state_[0][r_ind]*conj(Wnr_state_[0][r_ind]) +
      Wnr_state_[1][r_ind]*conj(Wnr_state_[1][r_ind]) )
                            *exp(iota_complex*(qx_*((rx_min + rx_ind*d_rx)) + qy_*((ry_min + ry_ind*d_ry))));
                }
            }


}
}



    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;

            Mq_[qx_ind][qy_ind]=0.0;
            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    r_ind=rx_ind + (space_slices)*ry_ind;

                    Mq_[qx_ind][qy_ind] += d_rx*d_ry*(
                                Wnr_state_[0][r_ind]*conj(Wnr_state_[0][r_ind]) +
                            Wnr_state_[1][r_ind]*conj(Wnr_state_[1][r_ind]) )
                            *exp(iota_complex*(qx_*((rx_min + rx_ind*d_rx)) + qy_*((ry_min + ry_ind*d_ry))));
                }
            }
            //Mq_[qx_ind][qy_ind]=1.0;
        }
    }


/*
    double rx_, ry_;
    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;

            Vq_[qx_ind][qy_ind]=0.0;
            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    rx_ = rx_min + rx_ind*d_rx;
                    ry_ = ry_min + ry_ind*d_ry;

                    Vq_[qx_ind][qy_ind] += d_rx*d_ry*((14.3952*1000)/(Parameters_.eps_DE))*
                            ((1.0/(sqrt(rx_*rx_ + ry_*ry_ )+eta_q)) - screening*(1.0/ (sqrt(rx_*rx_ + ry_*ry_ + d_sqr)) ))
                            *exp(iota_complex*(qx_*((rx_)) + qy_*((ry_))));
                }
            }
        }
    }

*/

    double V_;
    for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            dis_x = ((sqrt(3.0)/2.0)*(r1_-(L1_eff/2)) +  (sqrt(3.0)/2.0)*(r2_-(L2_eff/2)))*Parameters_.a_moire;
            dis_y = (-0.5*(r1_-(L1_eff/2)) + 0.5*(r2_-(L2_eff/2)))*Parameters_.a_moire;

            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;
            Uij[center_][center_neigh]=0.0;
		for(int q_ind=0;q_ind<q_slices;q_ind++){
		 q_val = q_ind*d_q;
		for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
		 theta_val = theta_ind*d_theta;
		qx_=q_val*cos(theta_val);
		qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.3952*1000)/(Parameters_.eps_DE);
                    Uij[center_][center_neigh]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*abs(Mq_SphC[q_ind][theta_ind])*abs(Mq_SphC[q_ind][theta_ind]))
                            *exp(iota_complex*( qx_*(dis_x) +  qy_*(dis_y) ));
                }
            }
        }
    }





    cout<<endl;
    cout<<"--------------------------------------Uij[center][neigh]--------------------------------------"<<endl;
    for(int r2_=0;r2_<L2_eff;r2_++){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Uij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;

    string Mq_file="Mq.txt";
    ofstream MqFILE(Mq_file.c_str());

    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;
            MqFILE<<qx_<<"   "<<qy_<<"   "<<Mq_[qx_ind][qy_ind].real()<<"   "<<Mq_[qx_ind][qy_ind].imag()<<endl;
        }
        MqFILE<<endl;
    }


    string MqSphC_file="MqSphC.txt";
    ofstream MqSphCFILE(MqSphC_file.c_str());

    for(int q_ind=0;q_ind<q_slices;q_ind++){
	q_val = q_ind*d_q;
	for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
	theta_val = theta_ind*d_theta;         
            MqSphCFILE<<q_val<<"   "<<theta_val<<"   "<<Mq_SphC[q_ind][theta_ind].real()<<"   "<<Mq_SphC[q_ind][theta_ind].imag()<<endl;
        }
        MqSphCFILE<<endl;
    }




/*
    string Vq_file="Vq.txt";
    ofstream VqFILE(Vq_file.c_str());

    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;
            VqFILE<<qx_<<"   "<<qy_<<"   "<<Vq_[qx_ind][qy_ind].real()<<"   "<<Vq_[qx_ind][qy_ind].imag()<<endl;
        }
        VqFILE<<endl;
    }

    string Vr_file="Vr.txt";
    ofstream VrFILE(Vr_file.c_str());

    for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
        for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
            rx_ = rx_min + rx_ind*d_rx;
            ry_ = ry_min + ry_ind*d_ry;

            VrFILE<<rx_<<"   "<<ry_<<"   "<< ((14.3952*1000)/(Parameters_.eps_DE))*((1.0/(sqrt(rx_*rx_ + ry_*ry_ )+eta_q)) - screening*(1.0/ (sqrt(rx_*rx_ + ry_*ry_ + d_sqr)) )) <<endl;
        }
        VrFILE<<endl;
    }

*/
    
    //---------------------------------------------------------------------------


}

void Hamiltonian::Hoppings(){

} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


#endif
