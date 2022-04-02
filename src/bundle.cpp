#include "bundle.h"


/* Cost functions */
class cResidualError;//pose and 3d points as parameters
class cResidualOnPose;//only pose as parameter
class cResidualOnViewPose; //local to global view transformation with covariance propagation
class cPoseConstraint;//constrant on pose

template <typename T>
bool cResidualOnViewPose::operator()(const T* const C, const T* const W, T* Residual) const
{
    
     /*    1, -Wz, Wy,
           Wz, 1, -Wx,
          -Wy, Wx, 1; */
    //Iw =  inv(r0) @ M @ R0 @ (I + WMat) => output is also skew

    //lambda * alpha * C  + beta - c = e
    T res_cx = m_lambda * (m_alpha(0,0)*C[0] +  m_alpha(0,1)*C[1] +  m_alpha(0,2)*C[2]) + m_beta[0] - m_c0[0];
    T res_cy = m_lambda * (m_alpha(1,0)*C[0] +  m_alpha(1,1)*C[1] +  m_alpha(1,2)*C[2]) + m_beta[1] - m_c0[1];
    T res_cz = m_lambda * (m_alpha(2,0)*C[0] +  m_alpha(2,1)*C[1] +  m_alpha(2,2)*C[2]) + m_beta[2] - m_c0[2];
    //std::cout << res_cx << " " << res_cy << " " << res_cz << "\n";

    // r = alpha * R , global to local
    // r0 (w+I) = alpha R0 (W+I) 
    //     w+I = r0**-1 alpha * R0 * (I+W)
    //
    //w   = r0**-1 alpha * R0 * (I+W) -I => r0**-1 alpha * R0 * (I+W) -I -w = e  (here w=0) 
    T res_wx = m_r0_inv_alpha_R0(2,0)*(-W[2]) + m_r0_inv_alpha_R0(2,1)*1 + m_r0_inv_alpha_R0(2,2)*W[0]; //wx (2,1)
    T res_wy = m_r0_inv_alpha_R0(0,0)*W[1] + m_r0_inv_alpha_R0(0,1)*(-W[0]) + m_r0_inv_alpha_R0(0,2)*1; //wy (0,2)
    T res_wz = m_r0_inv_alpha_R0(1,0)*1 + m_r0_inv_alpha_R0(1,1)*(W[2]) + m_r0_inv_alpha_R0(1,2)*(-W[1]); //wz (1,0)

    //weigh by the hessian = 1/sigma 
    Residual[0] = m_cov(0,0)*res_cx + m_cov(0,1)*res_cy + m_cov(0,2)*res_cz + m_cov(0,3)*res_wx + m_cov(0,4)*res_wy + m_cov(0,5)*res_wz; 
    Residual[1] = m_cov(1,0)*res_cx + m_cov(1,1)*res_cy + m_cov(1,2)*res_cz + m_cov(1,3)*res_wx + m_cov(1,4)*res_wy + m_cov(1,5)*res_wz;
    Residual[2] = m_cov(2,0)*res_cx + m_cov(2,1)*res_cy + m_cov(2,2)*res_cz + m_cov(2,3)*res_wx + m_cov(2,4)*res_wy + m_cov(2,5)*res_wz;
    Residual[3] = m_cov(3,0)*res_cx + m_cov(3,1)*res_cy + m_cov(3,2)*res_cz + m_cov(3,3)*res_wx + m_cov(3,4)*res_wy + m_cov(3,5)*res_wz;
    Residual[4] = m_cov(4,0)*res_cx + m_cov(4,1)*res_cy + m_cov(4,2)*res_cz + m_cov(4,3)*res_wx + m_cov(4,4)*res_wy + m_cov(4,5)*res_wz;
    Residual[5] = m_cov(5,0)*res_cx + m_cov(5,1)*res_cy + m_cov(5,2)*res_cz + m_cov(5,3)*res_wx + m_cov(5,4)*res_wy + m_cov(5,5)*res_wz;

    /*std::cout << "res pure: " << res_cx << " " << res_cy << " " << res_cz << " " 
                              << res_wx << " " << res_wy << " " << res_wz << "\n";
    std::cout << "Residuals : " << Residual[0] << " " << Residual[1] << " " << Residual[2] << " " 
                                << Residual[3] << " " << Residual[4] << " " << Residual[5] << "\n";*/

    return true;
}

CostFunction * cResidualOnViewPose::Create(const Mat3d alpha, const Vec3d beta, const double lambda,
                                           const Mat3d r0, const Vec3d c0,const Mat3d R0, const Mat6d covariance)
{
    return  (new AutoDiffCostFunction<cResidualOnViewPose,6,3,3> (new cResidualOnViewPose(alpha,beta,lambda,r0,c0,R0,covariance)));
}

template <typename T>
bool cResidualOnPose::operator()(const T* const aW,
                                const T* const aC,
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
    T XPC = mPt3d[0]-aC[0];
    T YPC = mPt3d[1]-aC[1];
    T ZPC = mPt3d[2]-aC[2];

    // Coordinate of points in  camera coordinate system, do not integrate "tiny" rotation
    T  XCam0 = mRot0(0,0)*XPC +  mRot0(1,0)*YPC +  mRot0(2,0)*ZPC;
    T  YCam0 = mRot0(0,1)*XPC +  mRot0(1,1)*YPC +  mRot0(2,1)*ZPC;
    T  ZCam0 = mRot0(0,2)*XPC +  mRot0(1,2)*YPC +  mRot0(2,2)*ZPC;

    // Now "tiny" rotation
    //  Wx      X      Wy * Z - Wz * Y
    //  Wy  ^   Y  =   Wz * X - Wx * Z
    //  Wz      Z      Wx * Y - Wy * X

    //  apply small rotation P =  P0 + W ^ P0
    T  XCam = XCam0 + Wy * ZCam0 - Wz * YCam0;
    T  YCam = YCam0 + Wz * XCam0 - Wx * ZCam0;
    T  ZCam = ZCam0 + Wx * YCam0 - Wy * XCam0;


    // project to bundle
    T xPi =  XCam/ZCam;
    T yPi =  YCam/ZCam;

    // prediction - observation
    Residual[0] = (xPi - mPtBundle(0,0))/mPdsSqrt;
    Residual[1] = (yPi - mPtBundle(1,0))/mPdsSqrt;



    return true;
}

CostFunction * cResidualOnPose::Create(const Mat3d mRotCur,const Vec2d PtBundle,const double* Pt3d,const double PdsSqrt)
{
    return  (new AutoDiffCostFunction<cResidualOnPose,2,3,3> (new cResidualOnPose(mRotCur,PtBundle,Pt3d,PdsSqrt)));
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
    //std::cout << "XYZCam0=" << XCam0 << " " << YCam0 << " " << ZCam0 << "<=>\n";

    // Now "tiny" rotation
    //  Wx      X      Wy * Z - Wz * Y
    //  Wy  ^   Y  =   Wz * X - Wx * Z
    //  Wz      Z      Wx * Y - Wy * X

    //  apply small rotation P =  P0 + W ^ P0
    T  XCam = XCam0 + Wy * ZCam0 - Wz * YCam0;
    T  YCam = YCam0 + Wz * XCam0 - Wx * ZCam0;
    T  ZCam = ZCam0 + Wx * YCam0 - Wy * XCam0;
    //std::cout << "XYZCam=" << XCam << " " << YCam << " " << ZCam << "<=>\n";


    // project to bundle
    T xPi =  XCam/ZCam;
    T yPi =  YCam/ZCam;
    //std::cout << "xPi yPi=" << xPi << " " << yPi << "<=>\n";

    // prediction - observation
    Residual[0] = (xPi - mPtBundle(0,0))/mPdsSqrt;
    Residual[1] = (yPi - mPtBundle(1,0))/mPdsSqrt;


    /*std::cout << "eee\n W " << Wx << " " << Wy << " " << Wz << "\n"
              << "C "   << aC[0] << " " << aC[1] << " " << aC[2] << "\n"
              << "Pt3 " << aPt3d[0] << " " << aPt3d[1] << " " << aPt3d[2] << "\n"
              << "Rot\n" << mRot0 << "\n";

    std::cout << Residual[0] << " " << Residual[1] << ", bunP=" << xPi << "," << yPi
              << ", bunO=" << mPtBundle.transpose() <<  "-" << mPtBundle(0,0) << " " << mPtBundle(1,0) <<  "\n";
    */
    return true;

}

CostFunction * cResidualError::Create(const Mat3d mRotCur,const Vec2d PtBundle,const double PdsSqrt)
{
    return  (new AutoDiffCostFunction<cResidualError,2,3,3,3> (new cResidualError(mRotCur,PtBundle,PdsSqrt)));
}
