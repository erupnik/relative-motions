#include "bundle.h"


/* Cost functions */
class cResidualError;//pose and 3d points as parameters
class cResidualOnPose;//only pose as parameter
class cPoseConstraint;//constrant on pose
class cResidualOn3ViewPoseDecomp;
class cResidualOn2ViewPoseDecomp;
class cResidualOn3ViewsPoseDecompLAB;

/*    1, -Wz, Wy,
     Wz, 1, -Wx,
    -Wy, Wx, 1; */

template <typename T>
bool cResidualOn3ViewsPoseDecomp::operator() (const T* const C0, const T* const W0,
                                              const T* const C1, const T* const W1,
                                              const T* const C2, const T* const W2, T* Residual) const
{
    Eigen::Matrix<T,18,1> c_wr;
  

    //calculate c 
    // C = 1.0/lambda * alpha.inverse() * (c - beta);
    // c  = lambda * alpha * C + beta 

    c_wr[0] = T(m_beta[0]) + T(m_L) * ( T(m_alpha(0,0))*C0[0] + T(m_alpha(0,1))*C0[1] + T(m_alpha(0,2))*C0[2] );
    c_wr[1] = T(m_beta[1]) + T(m_L) * ( T(m_alpha(1,0))*C0[0] + T(m_alpha(1,1))*C0[1] + T(m_alpha(1,2))*C0[2] );
    c_wr[2] = T(m_beta[2]) + T(m_L) * ( T(m_alpha(2,0))*C0[0] + T(m_alpha(2,1))*C0[1] + T(m_alpha(2,2))*C0[2] );



    c_wr[6] = T(m_beta[0]) + m_L * ( T(m_alpha(0,0))*C1[0] + T(m_alpha(0,1))*C1[1] + T(m_alpha(0,2))*C1[2] );
    c_wr[7] = T(m_beta[1]) + m_L * ( T(m_alpha(1,0))*C1[0] + T(m_alpha(1,1))*C1[1] + T(m_alpha(1,2))*C1[2] );
    c_wr[8] = T(m_beta[2]) + m_L * ( T(m_alpha(2,0))*C1[0] + T(m_alpha(2,1))*C1[1] + T(m_alpha(2,2))*C1[2] );



    c_wr[12] = T(m_beta[0]) + m_L * ( T(m_alpha(0,0))*C2[0] + T(m_alpha(0,1))*C2[1] + T(m_alpha(0,2))*C2[2] );
    c_wr[13] = T(m_beta[1]) + m_L * ( T(m_alpha(1,0))*C2[0] + T(m_alpha(1,1))*C2[1] + T(m_alpha(1,2))*C2[2] );
    c_wr[14] = T(m_beta[2]) + m_L * ( T(m_alpha(2,0))*C2[0] + T(m_alpha(2,1))*C2[1] + T(m_alpha(2,2))*C2[2] );
   

    Eigen::Matrix<T, 3, 3> Id = Eigen::Matrix<T, 3, 3>::Identity();
   
    //calculate wr 
    //r = alpha * R  
    // Wr = r0.inv * alpha * R0 * (WR + I) - I //


    //wr0 
    Eigen::Matrix<T, 3, 3> WR0plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR0plusI(0,1) = -W0[2];
    WR0plusI(1,0) = W0[2];
    WR0plusI(0,2) = W0[1];
    WR0plusI(2,0) = -W0[1];
    WR0plusI(2,1) = W0[0];
    WR0plusI(1,2) = -W0[0];

    Eigen::Matrix<T,3,3> r0_inv = m_rV[0].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R00    = m_R0V[0].cast<T>(); 
    Eigen::Matrix<T,3,3> wr0 = r0_inv * m_alpha * R00 * WR0plusI - Id;


    //wr1 
    Eigen::Matrix<T, 3, 3> WR1plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR1plusI(0,1) = -W1[2];
    WR1plusI(1,0) = W1[2];
    WR1plusI(0,2) = W1[1];
    WR1plusI(2,0) = -W1[1];
    WR1plusI(2,1) = W1[0];
    WR1plusI(1,2) = -W1[0];

    Eigen::Matrix<T,3,3> r1_inv = m_rV[1].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R10    = m_R0V[1].cast<T>(); 
    Eigen::Matrix<T,3,3> wr1 = r1_inv * m_alpha * R10 * WR1plusI - Id;

    //wr2 
    Eigen::Matrix<T, 3, 3> WR2plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR2plusI(0,1) = -W2[2];
    WR2plusI(1,0) = W2[2];
    WR2plusI(0,2) = W2[1];
    WR2plusI(2,0) = -W2[1];
    WR2plusI(2,1) = W2[0];
    WR2plusI(1,2) = -W2[0];

    Eigen::Matrix<T,3,3> r2_inv = m_rV[2].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R20    = m_R0V[2].cast<T>(); 
    Eigen::Matrix<T,3,3> wr2 = r2_inv * m_alpha * R20 * WR2plusI - Id;

    //update c_wr
    c_wr[3] = wr0(2,1);
    c_wr[4] = wr0(0,2);
    c_wr[5] = wr0(1,0);

    c_wr[9]  = wr1(2,1);
    c_wr[10] = wr1(0,2);
    c_wr[11] = wr1(1,0);


    c_wr[15] = wr2(2,1);
    c_wr[16] = wr2(0,2);
    c_wr[17] = wr2(1,0);

    /*for (int i=0; i<18; i++)
        std::cout << c_wr[i] << " \n" ;
    std::cout << "\n";
*/ 

    //res = w (l * C - cste)
    Eigen::Matrix<T,18,1> res_c_wr = 1.0/m_tot_res * (m_Li * c_wr - m_Cstei);
    
   /* for (int i=0; i<18; i++)
        std::cout << res_c_wr[i] << " \n" ;
    std::cout << "\n";*/ 
 
    Residual[0] = T(std::sqrt(m_Wi(0))) * res_c_wr(0);
    Residual[1] = T(std::sqrt(m_Wi(1))) * res_c_wr(1);
    Residual[2] = T(std::sqrt(m_Wi(2))) * res_c_wr(2);

    Residual[3] = T(std::sqrt(m_Wi(3))) * res_c_wr(3);
    Residual[4] = T(std::sqrt(m_Wi(4))) * res_c_wr(4);
    Residual[5] = T(std::sqrt(m_Wi(5))) * res_c_wr(5);

    Residual[6] = T(std::sqrt(m_Wi(6))) * res_c_wr(6);
    Residual[7] = T(std::sqrt(m_Wi(7))) * res_c_wr(7);
    Residual[8] = T(std::sqrt(m_Wi(8))) * res_c_wr(8);

    Residual[9] = T(std::sqrt(m_Wi(9))) * res_c_wr(9);
    Residual[10] = T(std::sqrt(m_Wi(10))) * res_c_wr(10);
    Residual[11] = T(std::sqrt(m_Wi(11))) * res_c_wr(11);

    Residual[12] = T(std::sqrt(m_Wi(12))) * res_c_wr(12);
    Residual[13] = T(std::sqrt(m_Wi(13))) * res_c_wr(13);
    Residual[14] = T(std::sqrt(m_Wi(14))) * res_c_wr(14);
    Residual[15] = T(std::sqrt(m_Wi(15))) * res_c_wr(15);
    Residual[16] = T(std::sqrt(m_Wi(16))) * res_c_wr(16);
    Residual[17] = T(std::sqrt(m_Wi(17))) * res_c_wr(17);

    /*for (int i=0; i<18; i++)
        std::cout << Residual[i] << " \n" ;
    std::cout << "\n";
    getchar();*/


    return true;
}

CostFunction* cResidualOn3ViewsPoseDecomp::Create(const Mat3d alpha, const Vec3d beta, const double L, 
                                                 const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                                 const VecXd Wi, const MatXd Li, const VecXd Cstei,const double total_res)
{
    return  (new AutoDiffCostFunction<cResidualOn3ViewsPoseDecomp,18,3,3,3,3,3,3> (new cResidualOn3ViewsPoseDecomp(alpha,beta,L,
                                                                                                                rV,R0V,
                                                                                                                Wi,Li,Cstei,
                                                                                                                total_res)));
}


template <typename T>
bool cResidualOn2ViewsPoseDecomp::operator() (const T* const C0, const T* const W0,
                                              const T* const C1, const T* const W1, T* Residual) const
{
    Eigen::Matrix<T,12,1> c_wr;
   
    //calculate c 
    // C = 1.0/lambda * alpha.inverse() * (c - beta);
    // c  = lambda * alpha * C + beta 

    c_wr[0] = m_beta[0] + m_lambda[0] * ( T(m_alpha(0,0))*C0[0] + T(m_alpha(0,1))*C0[1] + T(m_alpha(0,2))*C0[2] );
    c_wr[1] = m_beta[1] + m_lambda[0] * ( T(m_alpha(1,0))*C0[0] + T(m_alpha(1,1))*C0[1] + T(m_alpha(1,2))*C0[2] );
    c_wr[2] = m_beta[2] + m_lambda[0] * ( T(m_alpha(2,0))*C0[0] + T(m_alpha(2,1))*C0[1] + T(m_alpha(2,2))*C0[2] );



    c_wr[6] = m_beta[0] + m_lambda[0] * ( T(m_alpha(0,0))*C1[0] + T(m_alpha(0,1))*C1[1] + T(m_alpha(0,2))*C1[2] );
    c_wr[7] = m_beta[1] + m_lambda[0] * ( T(m_alpha(1,0))*C1[0] + T(m_alpha(1,1))*C1[1] + T(m_alpha(1,2))*C1[2] );
    c_wr[8] = m_beta[2] + m_lambda[0] * ( T(m_alpha(2,0))*C1[0] + T(m_alpha(2,1))*C1[1] + T(m_alpha(2,2))*C1[2] );
   

    Eigen::Matrix<T, 3, 3> Id = Eigen::Matrix<T, 3, 3>::Identity();
   
    //calculate wr 
    //r = alpha * R  
    // Wr = r0.inv * alpha * R0 * (WR + I) - I //


    //wr0 
    Eigen::Matrix<T, 3, 3> WR0plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR0plusI(0,1) = -W0[2];
    WR0plusI(1,0) = W0[2];
    WR0plusI(0,2) = W0[1];
    WR0plusI(2,0) = -W0[1];
    WR0plusI(2,1) = W0[0];
    WR0plusI(1,2) = -W0[0];

    Eigen::Matrix<T,3,3> r0_inv = m_rV[0].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R00    = m_R0V[0].cast<T>(); 
    Eigen::Matrix<T,3,3> wr0 = r0_inv * m_alpha * R00 * WR0plusI - Id;


    //wr1 
    Eigen::Matrix<T, 3, 3> WR1plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR1plusI(0,1) = -W1[2];
    WR1plusI(1,0) = W1[2];
    WR1plusI(0,2) = W1[1];
    WR1plusI(2,0) = -W1[1];
    WR1plusI(2,1) = W1[0];
    WR1plusI(1,2) = -W1[0];

    Eigen::Matrix<T,3,3> r1_inv = m_rV[1].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R10    = m_R0V[1].cast<T>(); 
    Eigen::Matrix<T,3,3> wr1 = r1_inv * m_alpha * R10 * WR1plusI - Id;

    //update c_wr
    c_wr[3] = wr0(2,1);
    c_wr[4] = wr0(0,2);
    c_wr[5] = wr0(1,0);

    c_wr[9]  = wr1(2,1);
    c_wr[10] = wr1(0,2);
    c_wr[11] = wr1(1,0);

    //res = w (l * C - cste)
    Eigen::Matrix<T,12,1> res_c_wr =  1.0/m_tot_res *  (m_Li * c_wr - m_Cstei);
    
    Residual[0] = T(std::sqrt(m_Wi(0))) * res_c_wr(0);
    Residual[1] = T(std::sqrt(m_Wi(1))) * res_c_wr(1);
    Residual[2] = T(std::sqrt(m_Wi(2))) * res_c_wr(2);

    Residual[3] = T(std::sqrt(m_Wi(3))) * res_c_wr(3);
    Residual[4] = T(std::sqrt(m_Wi(4))) * res_c_wr(4);
    Residual[5] = T(std::sqrt(m_Wi(5))) * res_c_wr(5);

    Residual[6] = T(std::sqrt(m_Wi(6))) * res_c_wr(6);
    Residual[7] = T(std::sqrt(m_Wi(7))) * res_c_wr(7);
    Residual[8] = T(std::sqrt(m_Wi(8))) * res_c_wr(8);

    Residual[9] = T(std::sqrt(m_Wi(9))) * res_c_wr(9);
    Residual[10] = T(std::sqrt(m_Wi(10))) * res_c_wr(10);
    Residual[11] = T(std::sqrt(m_Wi(11))) * res_c_wr(11);



    return true;
}

CostFunction* cResidualOn2ViewsPoseDecomp::Create(const Mat3d alpha, const double* beta, const double* lambda, 
                                                 const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                                 const VecXd Wi, const MatXd Li, const VecXd Cstei,const double total_res)
{
    return  (new AutoDiffCostFunction<cResidualOn2ViewsPoseDecomp,12,3,3,3,3> (new cResidualOn2ViewsPoseDecomp(alpha,beta,lambda,
                                                                                                                rV,R0V,
                                                                                                                Wi,Li,Cstei,
                                                                                                                total_res)));
}

/* r           = alpha           * R
 * r0 (I + Wr) = alpha0 (I + Wa) * R0 (I + Wr) 
 * so,
 * I + Wr = r0inv * alpha0 (I + Wa) * R0 (I + Wr) 
 * Wr = r0inv * alpha0 (I + Wa) * R0 (I + Wr) - I 
 *
 * */ 
template <typename T>
bool cResidualOn3ViewsPoseDecompLAB::operator() (const T* const C0, const T* const W0,
                                              const T* const C1, const T* const W1,
                                              const T* const C2, const T* const W2, 
                                              const T* const alpha_beta_L, T* Residual) const
{
    Eigen::Matrix<T,18,1> c_wr;

    const T* const  Wa   = alpha_beta_L;
    const T* const beta = alpha_beta_L+3;
    const T* const L    = alpha_beta_L+6;

    // c= Lambda * alpha0 (I+Wa) * C + beta 
    Eigen::Matrix<T, 3, 3> WaPlusI = Eigen::Matrix<T, 3, 3>::Identity();
    WaPlusI(0,1) = -Wa[2];
    WaPlusI(1,0) = Wa[2];
    WaPlusI(0,2) = Wa[1];
    WaPlusI(2,0) = -Wa[1];
    WaPlusI(2,1) = Wa[0];
    WaPlusI(1,2) = -Wa[0];

    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C0m(C0); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C1m(C1); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C2m(C2); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > betam(beta); 

    Eigen::Matrix<T,3,1> c0 = L[0] * m_alpha0 * WaPlusI * C0m + betam;
    Eigen::Matrix<T,3,1> c1 = L[0] * m_alpha0 * WaPlusI * C1m + betam;
    Eigen::Matrix<T,3,1> c2 = L[0] * m_alpha0 * WaPlusI * C2m + betam;

    // Wr = r0_inv * alpha0 (I+Wa)  * R0    (I+WR) -I  
    //Wr0 
    Eigen::Matrix<T, 3, 3> WR0plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR0plusI(0,1) = -W0[2];
    WR0plusI(1,0) = W0[2];
    WR0plusI(0,2) = W0[1];
    WR0plusI(2,0) = -W0[1];
    WR0plusI(2,1) = W0[0];
    WR0plusI(1,2) = -W0[0];

    Eigen::Matrix<T,3,3> r0_inv = m_rV[0].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R00    = m_R0V[0].cast<T>(); 
    Eigen::Matrix<T,3,3> wr0 = r0_inv * m_alpha0 * WaPlusI * R00 * WR0plusI - Eigen::Matrix<T, 3, 3>::Identity();

    //Wr1
    Eigen::Matrix<T, 3, 3> WR1plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR1plusI(0,1) = -W1[2];
    WR1plusI(1,0) = W1[2];
    WR1plusI(0,2) = W1[1];
    WR1plusI(2,0) = -W1[1];
    WR1plusI(2,1) = W1[0];
    WR1plusI(1,2) = -W1[0];

    Eigen::Matrix<T,3,3> r1_inv = m_rV[1].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R01    = m_R0V[1].cast<T>(); 
    Eigen::Matrix<T,3,3> wr1 = r1_inv * m_alpha0 * WaPlusI * R01 * WR1plusI - Eigen::Matrix<T, 3, 3>::Identity();

    //Wr2 
    Eigen::Matrix<T, 3, 3> WR2plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR2plusI(0,1) = -W2[2];
    WR2plusI(1,0) = W2[2];
    WR2plusI(0,2) = W2[1];
    WR2plusI(2,0) = -W2[1];
    WR2plusI(2,1) = W2[0];
    WR2plusI(1,2) = -W2[0];

    Eigen::Matrix<T,3,3> r2_inv = m_rV[2].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R02    = m_R0V[2].cast<T>(); 
    Eigen::Matrix<T,3,3> wr2 = r2_inv * m_alpha0 * WaPlusI * R02 * WR2plusI - Eigen::Matrix<T, 3, 3>::Identity();

   
    //update the c_wr vector 
    c_wr[0] = c0[0];
    c_wr[1] = c0[1];
    c_wr[2] = c0[2];
    c_wr[3] = wr0(2,1);
    c_wr[4] = wr0(0,2);
    c_wr[5] = wr0(1,0);

    c_wr[6] = c1[0];
    c_wr[7] = c1[1];
    c_wr[8] = c1[2];
    c_wr[9] = wr1(2,1);
    c_wr[10] = wr1(0,2);
    c_wr[11] = wr1(1,0);

    c_wr[12] = c2[0];
    c_wr[13] = c2[1];
    c_wr[14] = c2[2];
    c_wr[15] = wr2(2,1);
    c_wr[16] = wr2(0,2);
    c_wr[17] = wr2(1,0);

/*    std::cout << "Xpred\n";
    for (int i=0; i<18; i++)
        std::cout << c_wr[i] << " ";
    std::cout << "\n";
getchar();*/ 

    //res = w (l * C - cste)
    Eigen::Matrix<T,18,1> res_c_wr =  1.0/m_tot_res * (m_Li * c_wr - m_Cstei);


    Residual[0] = T(std::sqrt(m_Wi(0))) * res_c_wr(0);
    Residual[1] = T(std::sqrt(m_Wi(1))) * res_c_wr(1);
    Residual[2] = T(std::sqrt(m_Wi(2))) * res_c_wr(2);

    Residual[3] = T(std::sqrt(m_Wi(3))) * res_c_wr(3);
    Residual[4] = T(std::sqrt(m_Wi(4))) * res_c_wr(4);
    Residual[5] = T(std::sqrt(m_Wi(5))) * res_c_wr(5);

    Residual[6] = T(std::sqrt(m_Wi(6))) * res_c_wr(6);
    Residual[7] = T(std::sqrt(m_Wi(7))) * res_c_wr(7);
    Residual[8] = T(std::sqrt(m_Wi(8))) * res_c_wr(8);

    Residual[9] = T(std::sqrt(m_Wi(9))) * res_c_wr(9);
    Residual[10] = T(std::sqrt(m_Wi(10))) * res_c_wr(10);
    Residual[11] = T(std::sqrt(m_Wi(11))) * res_c_wr(11);

    Residual[12] = T(std::sqrt(m_Wi(12))) * res_c_wr(12);
    Residual[13] = T(std::sqrt(m_Wi(13))) * res_c_wr(13);
    Residual[14] = T(std::sqrt(m_Wi(14))) * res_c_wr(14);
    Residual[15] = T(std::sqrt(m_Wi(15))) * res_c_wr(15);
    Residual[16] = T(std::sqrt(m_Wi(16))) * res_c_wr(16);
    Residual[17] = T(std::sqrt(m_Wi(17))) * res_c_wr(17);

    return true;
}

CostFunction* cResidualOn3ViewsPoseDecompLAB::Create(const Mat3d alpha0, 
                                                     const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                                     const VecXd Wi, const MatXd Li, const VecXd Cstei,const double total_res)
{

    return  (new AutoDiffCostFunction<cResidualOn3ViewsPoseDecompLAB,18,3,3,3,3,3,3,7> (new cResidualOn3ViewsPoseDecompLAB(alpha0,                                                                                                                rV,R0V,
                                                                                                                Wi,Li,Cstei,
                                                                                                                total_res)));
}

template <typename T>
bool cResidualOn2ViewsPoseDecompLAB::operator() (const T* const C0, const T* const W0,
                                              const T* const C1, const T* const W1,
                                              const T* const alpha_beta_L, T* Residual) const
{
    Eigen::Matrix<T,12,1> c_wr;

    const T* const Wa   = alpha_beta_L;
    const T* const beta = alpha_beta_L+3;
    const T* const L    = alpha_beta_L+6;

    // c= Lambda * alpha0 (I+Wa) * C + beta 
    Eigen::Matrix<T, 3, 3> WaPlusI = Eigen::Matrix<T, 3, 3>::Identity();
    WaPlusI(0,1) = -Wa[2];
    WaPlusI(1,0) = Wa[2];
    WaPlusI(0,2) = Wa[1];
    WaPlusI(2,0) = -Wa[1];
    WaPlusI(2,1) = Wa[0];
    WaPlusI(1,2) = -Wa[0];

    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C0m(C0); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C1m(C1); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > betam(beta); 

    Eigen::Matrix<T,3,1> c0 = L[0] * m_alpha0 * WaPlusI * C0m + betam;
    Eigen::Matrix<T,3,1> c1 = L[0] * m_alpha0 * WaPlusI * C1m + betam;

    // Wr = r0_inv * alpha0 (I+Wa) (I+WR) -I 
   
    //Wr0 
    Eigen::Matrix<T, 3, 3> WR0plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR0plusI(0,1) = -W0[2];
    WR0plusI(1,0) = W0[2];
    WR0plusI(0,2) = W0[1];
    WR0plusI(2,0) = -W0[1];
    WR0plusI(2,1) = W0[0];
    WR0plusI(1,2) = -W0[0];

    Eigen::Matrix<T,3,3> r0_inv = m_rV[0].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R00    = m_R0V[0].cast<T>(); 
    Eigen::Matrix<T,3,3> wr0 = r0_inv * m_alpha0 * WaPlusI * R00 * WR0plusI - Eigen::Matrix<T, 3, 3>::Identity();

    //Wr1
    Eigen::Matrix<T, 3, 3> WR1plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR1plusI(0,1) = -W1[2];
    WR1plusI(1,0) = W1[2];
    WR1plusI(0,2) = W1[1];
    WR1plusI(2,0) = -W1[1];
    WR1plusI(2,1) = W1[0];
    WR1plusI(1,2) = -W1[0];

    Eigen::Matrix<T,3,3> r1_inv = m_rV[1].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R01    = m_R0V[1].cast<T>(); 
    Eigen::Matrix<T,3,3> wr1 = r1_inv * m_alpha0 * WaPlusI * R01 * WR1plusI - Eigen::Matrix<T, 3, 3>::Identity();

    //update the c_wr vector 
    c_wr[0] = c0[0];
    c_wr[1] = c0[1];
    c_wr[2] = c0[2];
    c_wr[3] = wr0(2,1);
    c_wr[4] = wr0(0,2);
    c_wr[5] = wr0(1,0);

    c_wr[6] = c1[0];
    c_wr[7] = c1[1];
    c_wr[8] = c1[2];
    c_wr[9] = wr1(2,1);
    c_wr[10] = wr1(0,2);
    c_wr[11] = wr1(1,0);

    //res = w (l * C - cste)
    Eigen::Matrix<T,12,1> res_c_wr = 1.0/m_tot_res * (m_Li * c_wr - m_Cstei);


    Residual[0] = T(std::sqrt(m_Wi(0))) * res_c_wr(0);
    Residual[1] = T(std::sqrt(m_Wi(1))) * res_c_wr(1);
    Residual[2] = T(std::sqrt(m_Wi(2))) * res_c_wr(2);

    Residual[3] = T(std::sqrt(m_Wi(3))) * res_c_wr(3);
    Residual[4] = T(std::sqrt(m_Wi(4))) * res_c_wr(4);
    Residual[5] = T(std::sqrt(m_Wi(5))) * res_c_wr(5);

    Residual[6] = T(std::sqrt(m_Wi(6))) * res_c_wr(6);
    Residual[7] = T(std::sqrt(m_Wi(7))) * res_c_wr(7);
    Residual[8] = T(std::sqrt(m_Wi(8))) * res_c_wr(8);

    Residual[9] = T(std::sqrt(m_Wi(9))) * res_c_wr(9);
    Residual[10] = T(std::sqrt(m_Wi(10))) * res_c_wr(10);
    Residual[11] = T(std::sqrt(m_Wi(11))) * res_c_wr(11);



    return true;
}

CostFunction* cResidualOn2ViewsPoseDecompLAB::Create(const Mat3d alpha0, 
                                                     const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                                     const VecXd Wi, const MatXd Li, const VecXd Cstei,const double total_res)
{

    return  (new AutoDiffCostFunction<cResidualOn2ViewsPoseDecompLAB,12,3,3,3,3,7> (new cResidualOn2ViewsPoseDecompLAB(alpha0,                                                                                                                rV,R0V,
                                                                                                                Wi,Li,Cstei,
                                                                                                                total_res)));
}

template <typename T>
bool cResidualOn3ViewsPoseBasicLAB::operator()(const T* const C0, const T* const W0,
                                               const T* const C1, const T* const W1,
                                               const T* const C2, const T* const W2, 
                                               const T* const alpha_beta_L, T* Residual) const
{
    Eigen::Matrix<T,18,1> c_wr;

    const T* const  Wa   = alpha_beta_L;
    const T* const beta = alpha_beta_L+3;
    const T* const L    = alpha_beta_L+6;

    // c= Lambda * alpha0 (I+Wa) * C + beta 
    Eigen::Matrix<T, 3, 3> WaPlusI = Eigen::Matrix<T, 3, 3>::Identity();
    WaPlusI(0,1) = -Wa[2];
    WaPlusI(1,0) = Wa[2];
    WaPlusI(0,2) = Wa[1];
    WaPlusI(2,0) = -Wa[1];
    WaPlusI(2,1) = Wa[0];
    WaPlusI(1,2) = -Wa[0];

    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C0m(C0); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C1m(C1); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C2m(C2); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > betam(beta); 

    Eigen::Matrix<T,3,1> c0 = L[0] * m_alpha0 * WaPlusI * C0m + betam;
    Eigen::Matrix<T,3,1> c1 = L[0] * m_alpha0 * WaPlusI * C1m + betam;
    Eigen::Matrix<T,3,1> c2 = L[0] * m_alpha0 * WaPlusI * C2m + betam;

    // Wr = r0_inv * alpha0 (I+Wa)  * R0    (I+WR) -I  
    //Wr0 
    Eigen::Matrix<T, 3, 3> WR0plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR0plusI(0,1) = -W0[2];
    WR0plusI(1,0) = W0[2];
    WR0plusI(0,2) = W0[1];
    WR0plusI(2,0) = -W0[1];
    WR0plusI(2,1) = W0[0];
    WR0plusI(1,2) = -W0[0];

    Eigen::Matrix<T,3,3> r0_inv = m_rV[0].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R00    = m_R0V[0].cast<T>(); 
    Eigen::Matrix<T,3,3> wr0 = r0_inv * m_alpha0 * WaPlusI * R00 * WR0plusI - Eigen::Matrix<T, 3, 3>::Identity();

    //Wr1
    Eigen::Matrix<T, 3, 3> WR1plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR1plusI(0,1) = -W1[2];
    WR1plusI(1,0) = W1[2];
    WR1plusI(0,2) = W1[1];
    WR1plusI(2,0) = -W1[1];
    WR1plusI(2,1) = W1[0];
    WR1plusI(1,2) = -W1[0];

    Eigen::Matrix<T,3,3> r1_inv = m_rV[1].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R01    = m_R0V[1].cast<T>(); 
    Eigen::Matrix<T,3,3> wr1 = r1_inv * m_alpha0 * WaPlusI * R01 * WR1plusI - Eigen::Matrix<T, 3, 3>::Identity();

    //Wr2 
    Eigen::Matrix<T, 3, 3> WR2plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR2plusI(0,1) = -W2[2];
    WR2plusI(1,0) = W2[2];
    WR2plusI(0,2) = W2[1];
    WR2plusI(2,0) = -W2[1];
    WR2plusI(2,1) = W2[0];
    WR2plusI(1,2) = -W2[0];

    Eigen::Matrix<T,3,3> r2_inv = m_rV[2].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R02    = m_R0V[2].cast<T>(); 
    Eigen::Matrix<T,3,3> wr2 = r2_inv * m_alpha0 * WaPlusI * R02 * WR2plusI - Eigen::Matrix<T, 3, 3>::Identity();

   
    //update the c_wr vector 
    c_wr[0] = c0[0];
    c_wr[1] = c0[1];
    c_wr[2] = c0[2];
    c_wr[3] = wr0(2,1);
    c_wr[4] = wr0(0,2);
    c_wr[5] = wr0(1,0);

    c_wr[6] = c1[0];
    c_wr[7] = c1[1];
    c_wr[8] = c1[2];
    c_wr[9] = wr1(2,1);
    c_wr[10] = wr1(0,2);
    c_wr[11] = wr1(1,0);

    c_wr[12] = c2[0];
    c_wr[13] = c2[1];
    c_wr[14] = c2[2];
    c_wr[15] = wr2(2,1);
    c_wr[16] = wr2(0,2);
    c_wr[17] = wr2(1,0);

    //- check whether the input Wr is the original Wr or the one corrected by the local bundle
    //if it is the original then Wr = 0 
    //- print residuals to make sure equations are OK


    Residual[0]  = (m_cV[0][0] - c_wr[0])/mPdsSqrt_c; 
    Residual[1]  = (m_cV[0][1] - c_wr[1])/mPdsSqrt_c; 
    Residual[2]  = (m_cV[0][2] - c_wr[2])/mPdsSqrt_c; 

    Residual[3]  = (0.0 - c_wr[3])/mPdsSqrt_w;
    Residual[4]  = (0.0 - c_wr[4])/mPdsSqrt_w;
    Residual[5]  = (0.0 - c_wr[5])/mPdsSqrt_w;

    Residual[6]  = (m_cV[1][0] - c_wr[6])/mPdsSqrt_c; 
    Residual[7]  = (m_cV[1][1] - c_wr[7])/mPdsSqrt_c;
    Residual[8]  = (m_cV[1][2] - c_wr[8])/mPdsSqrt_c;

    Residual[9]  = (0.0 - c_wr[9])/mPdsSqrt_w; 
    Residual[10] = (0.0 - c_wr[10])/mPdsSqrt_w;
    Residual[11] = (0.0 - c_wr[11])/mPdsSqrt_w;

    Residual[12] = (m_cV[2][0] - c_wr[12])/mPdsSqrt_c;
    Residual[13] = (m_cV[2][1] - c_wr[13])/mPdsSqrt_c;
    Residual[14] = (m_cV[2][2] - c_wr[14])/mPdsSqrt_c;

    Residual[15] = (0.0 - c_wr[15])/mPdsSqrt_w;
    Residual[16] = (0.0 - c_wr[16])/mPdsSqrt_w;
    Residual[17] = (0.0 - c_wr[17])/mPdsSqrt_w;

    /*for (int i=0;i<18; i++)
        std::cout << Residual[i] << " " ;
    std::cout << "\n";

    getchar();*/


    return true;
}

CostFunction* cResidualOn3ViewsPoseBasicLAB::Create( const Mat3d alpha0, 
                                                     const std::vector<double*> cV,
                                                     const std::vector<Mat3d> rV, 
                                                     const std::vector<Mat3d> R0V,
                                                     const double             Pds_c,
                                                     const double             Pds_w)
{

    return  (new AutoDiffCostFunction<cResidualOn3ViewsPoseBasicLAB,18,3,3,3,3,3,3,7> 
                  (new cResidualOn3ViewsPoseBasicLAB(alpha0,cV,rV,R0V,Pds_c,Pds_w)));
}


template <typename T>
bool cResidualOn2ViewsPoseBasicLAB::operator()(const T* const C0, const T* const W0,
                                               const T* const C1, const T* const W1,
                                               const T* const alpha_beta_L, T* Residual) const
{
    Eigen::Matrix<T,12,1> c_wr;

    const T* const  Wa   = alpha_beta_L;
    const T* const beta = alpha_beta_L+3;
    const T* const L    = alpha_beta_L+6;

    // c= Lambda * alpha0 (I+Wa) * C + beta 
    Eigen::Matrix<T, 3, 3> WaPlusI = Eigen::Matrix<T, 3, 3>::Identity();
    WaPlusI(0,1) = -Wa[2];
    WaPlusI(1,0) = Wa[2];
    WaPlusI(0,2) = Wa[1];
    WaPlusI(2,0) = -Wa[1];
    WaPlusI(2,1) = Wa[0];
    WaPlusI(1,2) = -Wa[0];

    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C0m(C0); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > C1m(C1); 
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > betam(beta); 

    Eigen::Matrix<T,3,1> c0 = L[0] * m_alpha0 * WaPlusI * C0m + betam;
    Eigen::Matrix<T,3,1> c1 = L[0] * m_alpha0 * WaPlusI * C1m + betam;

    // Wr = r0_inv * alpha0 (I+Wa)  * R0    (I+WR) -I  
    //Wr0 
    Eigen::Matrix<T, 3, 3> WR0plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR0plusI(0,1) = -W0[2];
    WR0plusI(1,0) = W0[2];
    WR0plusI(0,2) = W0[1];
    WR0plusI(2,0) = -W0[1];
    WR0plusI(2,1) = W0[0];
    WR0plusI(1,2) = -W0[0];

    Eigen::Matrix<T,3,3> r0_inv = m_rV[0].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R00    = m_R0V[0].cast<T>(); 
    Eigen::Matrix<T,3,3> wr0 = r0_inv * m_alpha0 * WaPlusI * R00 * WR0plusI - Eigen::Matrix<T, 3, 3>::Identity();

    //Wr1
    Eigen::Matrix<T, 3, 3> WR1plusI =  Eigen::Matrix<T, 3, 3>::Identity();
    WR1plusI(0,1) = -W1[2];
    WR1plusI(1,0) = W1[2];
    WR1plusI(0,2) = W1[1];
    WR1plusI(2,0) = -W1[1];
    WR1plusI(2,1) = W1[0];
    WR1plusI(1,2) = -W1[0];

    Eigen::Matrix<T,3,3> r1_inv = m_rV[1].inverse().cast<T>(); 
    Eigen::Matrix<T,3,3> R01    = m_R0V[1].cast<T>(); 
    Eigen::Matrix<T,3,3> wr1 = r1_inv * m_alpha0 * WaPlusI * R01 * WR1plusI - Eigen::Matrix<T, 3, 3>::Identity();

    
    //update the c_wr vector 
    c_wr[0] = c0[0];
    c_wr[1] = c0[1];
    c_wr[2] = c0[2];
    c_wr[3] = wr0(2,1);
    c_wr[4] = wr0(0,2);
    c_wr[5] = wr0(1,0);

    c_wr[6] = c1[0];
    c_wr[7] = c1[1];
    c_wr[8] = c1[2];
    c_wr[9] = wr1(2,1);
    c_wr[10] = wr1(0,2);
    c_wr[11] = wr1(1,0);


    //- check whether the input Wr is the original Wr or the one corrected by the local bundle
    //if it is the original then Wr = 0 
    //- print residuals to make sure equations are OK


    Residual[0]  = (m_cV[0][0] - c_wr[0])/mPdsSqrt_c; 
    Residual[1]  = (m_cV[0][1] - c_wr[1])/mPdsSqrt_c; 
    Residual[2]  = (m_cV[0][2] - c_wr[2])/mPdsSqrt_c; 

    Residual[3]  = (0.0 - c_wr[3])/mPdsSqrt_w;
    Residual[4]  = (0.0 - c_wr[4])/mPdsSqrt_w;
    Residual[5]  = (0.0 - c_wr[5])/mPdsSqrt_w;

    Residual[6]  = (m_cV[1][0] - c_wr[6])/mPdsSqrt_c; 
    Residual[7]  = (m_cV[1][1] - c_wr[7])/mPdsSqrt_c;
    Residual[8]  = (m_cV[1][2] - c_wr[8])/mPdsSqrt_c;

    Residual[9]  = (0.0 - c_wr[9])/mPdsSqrt_w; 
    Residual[10] = (0.0 - c_wr[10])/mPdsSqrt_w;
    Residual[11] = (0.0 - c_wr[11])/mPdsSqrt_w;


    /*for (int i=0;i<12; i++)
        std::cout << Residual[i] << " " ;
    std::cout << "\n";

    getchar();*/


    return true;

}
CostFunction* cResidualOn2ViewsPoseBasicLAB::Create( const Mat3d alpha0, 
                                                     const std::vector<double*> cV,
                                                     const std::vector<Mat3d> rV, 
                                                     const std::vector<Mat3d> R0V,
                                                     const double             Pds_c,
                                                     const double             Pds_w)
{

    return  (new AutoDiffCostFunction<cResidualOn2ViewsPoseBasicLAB,12,3,3,3,3,7> 
            (new cResidualOn2ViewsPoseBasicLAB(alpha0,cV,rV,R0V,Pds_c,Pds_w)));
}




template <typename T>
bool cResidualError::operator()(const T* const aW,
                                const T* const aC,
                                const T* const aPt3d,
                                T* Residual) const
{
    // - move 3d point
    // - rotate by the current R
    // - rotate by the small rotation (parameter in estimation)
    // - project to bundle

    const T & Wx = aW[0];
    const T & Wy = aW[1];
    const T & Wz = aW[2];

    // Vector P->Cam
    T XPC = aPt3d[0]-aC[0];
    T YPC = aPt3d[1]-aC[1];
    T ZPC = aPt3d[2]-aC[2];
    //std::cout << "XYZ PC=" << XPC << " " << YPC << " " << ZPC << "<=>\n";
    // Coordinate of points in  camera coordinate system, do not integrate "tiny" rotation
    T  XCam0 = mRot0(0,0)*XPC +  mRot0(1,0)*YPC +  mRot0(2,0)*ZPC;
    T  YCam0 = mRot0(0,1)*XPC +  mRot0(1,1)*YPC +  mRot0(2,1)*ZPC;
    T  ZCam0 = mRot0(0,2)*XPC +  mRot0(1,2)*YPC +  mRot0(2,2)*ZPC;

    // Now "tiny" rotation
    //  Wx      X      Wy * Z - Wz * Y
    //  Wy  ^   Y  =   Wz * X - Wx * Z
    //  Wz      Z      Wx * Y - Wy * X

    /*//  apply small rotation P =  P0 + W ^ P0
    T  XCam = XCam0 + Wy * ZCam0 - Wz * YCam0;
    T  YCam = YCam0 + Wz * XCam0 - Wx * ZCam0;
    T  ZCam = ZCam0 + Wx * YCam0 - Wy * XCam0;
    //std::cout << "XYZCam=" << XCam << " " << YCam << " " << ZCam << "<=>\n";
*/
    //  apply small rotation P =  P0 + W^T ^ P0
    T  XCam = XCam0 +  Wz * YCam0 - Wy * ZCam0;
    T  YCam = YCam0 +  Wx * ZCam0 - Wz * XCam0;
    T  ZCam = ZCam0 +  Wy * XCam0 - Wx * YCam0;

    // project to bundle
    T xPi =  XCam/ZCam;
    T yPi =  YCam/ZCam;
    //std::cout << "xPi yPi=" << xPi << " " << yPi << "<=>\n";

    //double focal = 30975.0;
    //double px = 26460.0;
    //double py = 17004.0;
    // prediction - observation
    //Residual[0] = ((xPi - mPtBundle(0))*focal + px )/mPdsSqrt;
    //Residual[1] = ((yPi - mPtBundle(1))*focal + py )/mPdsSqrt;
    Residual[0] = (xPi - mPtBundle(0))*mFoc /mPdsSqrt;
    Residual[1] = (yPi - mPtBundle(1))*mFoc /mPdsSqrt;
   
   
    /*std::cout << "************ " << Residual[0] << " " << Residual[1] << "\n";
    std::cout << mPtBundle[0]*726.287+354.64968 << " " << mPtBundle[1]*726.287+186.4656 << "\n";
    std::cout << xPi*726.287+354.64968 << " " << yPi*726.287+186.4656 
              << ", 3D=" << aPt3d[0] << " " << aPt3d[1] << " " << aPt3d[2] << "\n";
    getchar();*/
/*
    T xpred = (xPi * 5652.66392) + 2614.60004014115202;
    T ypred = (yPi * 5652.66392) + 1886.05247082382 ;

    double xobs = (mPtBundle(0) * 5652.66392) + 2614.60004014115202;
    double yobs = (mPtBundle(1) * 5652.66392) + 1886.05247082382 ;

    std::cout << "Res in px \nx=" << xpred - xobs << "\ny=" << ypred - yobs << "\n";
    std::cout << "Res in bundle\n" << Residual[0] << "\n" << Residual[1] << "\n";
    
    std::cout << "coords= " << xobs << " " << yobs << "\n";
    std::cout << "C= " << aC[0] << " " << aC[1] << " " << aC[2] << "\n";
    getchar();*/
/*
          <PP>2614.60004014115202 1886.0524708238238</PP>
          <F>5652.66392026270023</F>
 * 
 * */

    /*std::cout << "eee\n W " << Wx << " " << Wy << " " << Wz << "\n"
              << "C "   << aC[0] << " " << aC[1] << " " << aC[2] << "\n"
              << "Pt3 " << aPt3d[0] << " " << aPt3d[1] << " " << aPt3d[2] << "\n"
              << "Rot\n" << mRot0 << "\n";

    std::cout << Residual[0] << " " << Residual[1] << ", bunP=" << xPi << "," << yPi
              << ", bunO=" << mPtBundle.transpose() <<  "-" << mPtBundle(0,0) << " " << mPtBundle(1,0) <<  "\n";
    */
    return true;

}

CostFunction * cResidualError::Create(const Mat3d mRotCur,const Vec2d PtBundle,const double PdsSqrt,const double Foc)
{
    return  (new AutoDiffCostFunction<cResidualError,2,3,3,3> (new cResidualError(mRotCur,PtBundle,PdsSqrt,Foc)));
}
