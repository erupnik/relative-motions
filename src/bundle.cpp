#include "bundle.h"


/* Cost functions */
class cResidualError;//pose and 3d points as parameters
class cResidualOnPose;//only pose as parameter
class cPoseConstraint;//constrant on pose

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
    std::cout << "XYZ PC=" << XPC << " " << YPC << " " << ZPC << "<=>\n";
    // Coordinate of points in  camera coordinate system, do not integrate "tiny" rotation
    T  XCam0 = mRot0(0,0)*XPC +  mRot0(1,0)*YPC +  mRot0(2,0)*ZPC;
    T  YCam0 = mRot0(0,1)*XPC +  mRot0(1,1)*YPC +  mRot0(2,1)*ZPC;
    T  ZCam0 = mRot0(0,2)*XPC +  mRot0(1,2)*YPC +  mRot0(2,2)*ZPC;
    std::cout << "XYZCam0=" << XCam0 << " " << YCam0 << " " << ZCam0 << "<=>\n";

    // Now "tiny" rotation
    //  Wx      X      Wy * Z - Wz * Y
    //  Wy  ^   Y  =   Wz * X - Wx * Z
    //  Wz      Z      Wx * Y - Wy * X

    //  apply small rotation P =  P0 + W ^ P0
    T  XCam = XCam0 + Wy * ZCam0 - Wz * YCam0;
    T  YCam = YCam0 + Wz * XCam0 - Wx * ZCam0;
    T  ZCam = ZCam0 + Wx * YCam0 - Wy * XCam0;
    std::cout << "XYZCam=" << XCam << " " << YCam << " " << ZCam << "<=>\n";


    // project to bundle
    T xPi =  XCam/ZCam;
    T yPi =  YCam/ZCam;
    std::cout << "xPi yPi=" << xPi << " " << yPi << "<=>\n";

    // prediction - observation
    Residual[0] = (xPi - mPtBundle(0,0))/mPdsSqrt;
    Residual[1] = (yPi - mPtBundle(1,0))/mPdsSqrt;


    std::cout << "eee\n W " << Wx << " " << Wy << " " << Wz << "\n"
              << "C "   << aC[0] << " " << aC[1] << " " << aC[2] << "\n"
              << "Pt3 " << aPt3d[0] << " " << aPt3d[1] << " " << aPt3d[2] << "\n"
              << "Rot\n" << mRot0 << "\n";

    std::cout << Residual[0] << " " << Residual[1] << ", bunP=" << xPi << "," << yPi
              << ", bunO=" << mPtBundle.transpose() <<  "-" << mPtBundle(0,0) << " " << mPtBundle(1,0) <<  "\n";
    return true;

}

CostFunction * cResidualError::Create(const Mat3d mRotCur,const Vec2d PtBundle,const double PdsSqrt)
{
    return  (new AutoDiffCostFunction<cResidualError,2,3,3,3> (new cResidualError(mRotCur,PtBundle,PdsSqrt)));
}
