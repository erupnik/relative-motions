#ifndef _BUNDLE_H
#define _BUNDLE_H

#include "all.h"

/* Ceres */
#include "ceres/ceres.h"
#include "ceres/types.h"
#include "ceres/rotation.h"

/* Eigen */
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

/* Logging */
#include "gflags/gflags.h"
#include "glog/logging.h"

using namespace ceres;


/*  pose_pds is the weight matrix (square root) */
class cPoseConstraint
{
    public:
      cPoseConstraint(const double* init_pose, double pose_pds) :
        init_pose(init_pose), pose_pds(pose_pds) {}

      ~cPoseConstraint(){}

      template <typename T>
      bool operator()(const T* const obs_pose,
                      T* Residual) const {
            Residual[0] = (obs_pose[0] - init_pose[0])/pose_pds;
            Residual[1] = (obs_pose[1] - init_pose[1])/pose_pds;
            Residual[2] = (obs_pose[2] - init_pose[2])/pose_pds;
            return true;
      }


      static CostFunction * Create(const double* pose_value, const double pose_pds){
          return  (new AutoDiffCostFunction<cPoseConstraint,3,3> (new cPoseConstraint(
            pose_value,pose_pds)));
      }

    private:

      const double*  init_pose;
      const double   pose_pds;
};

/* under dev */
class cResidualOnViewPoseAffFree
{
    public: 
        cResidualOnViewPoseAffFree(const Mat3d alpha0, const Mat3d r0, const Vec3d c0, const Mat3d R0, const Mat6d covariance) :
                                               m_alpha0(alpha0),
                                               m_r0(r0),
                                               m_c0(c0),
                                               m_R0(R0),
                                               m_cov(covariance) {}
        ~cResidualOnViewPoseAffFree(){}


        template <typename T>
            bool operator()(const T* const C, const T* const W, 
                            const T* const Walpha, const T* const beta, const T* const lambda, 
                            T* Residual) const;
                                
        static CostFunction* Create(const Mat3d alpha0,
                                    const Mat3d r0, const Vec3d c0, const Mat3d R0, const Mat6d covariance);

    private:
        const Mat3d m_alpha0;

        const Mat3d m_r0;
        const Vec3d m_c0;

        const Mat3d m_R0;

        const Mat6d m_cov;
};

class cResidualOnViewPose
{
    public: 
        cResidualOnViewPose(const Mat3d alpha, const double* beta, const double* lambda,
                            const Mat3d r0, const Vec3d c0, const Mat3d R0, const Mat6d covariance) :
                                               m_alpha(alpha),
                                               m_beta(beta),
                                               m_lambda(lambda),
                                               m_r0(r0),
                                               //m_r0_inv_alpha_R0(r0.inverse() * alpha * R0),
                                               m_alpha_R0(alpha * R0),
                                               m_c0(c0),
                                               m_R0(R0),
                                               m_cov(covariance) {}
        ~cResidualOnViewPose(){}

        template <typename T>
            bool operator()(const T* const C, const T* const W, T* Residual) const;

        static CostFunction* Create(const Mat3d alpha, const double* beta, const double* lambda,
                                    const Mat3d r0, const Vec3d c0, const Mat3d R0, const Mat6d covariance);
        //static CostFunction* Create(const Mat3d alpha, const Vec3d beta, const double lambda,
        //                           const Mat3d r0, const Vec3d c0, const Mat3d R0, const Mat6d covariance);

    private:
        const Mat3d m_alpha;
        const double* m_beta;
        const double* m_lambda;

        const Mat3d m_r0;
        const Mat3d m_alpha_R0;
        const Vec3d m_c0;

        const Mat3d m_R0;

        const Mat6d m_cov;
};

class cResidualOnPose
{
    public :
        cResidualOnPose(const Mat3d aRot0,const Vec2d aPtBundle,const double* aPt3d,const double aPdsSqrt) :
            mRot0(aRot0), mPtBundle(aPtBundle), mPt3d(aPt3d),mPdsSqrt(aPdsSqrt) {}
        ~cResidualOnPose(){}

        // Parameters: pose
        template <typename T>
        bool operator()(const T* const aW,
                        const T* const aC,
                        T* Residual) const;
        static CostFunction * Create(const Mat3d mRotCur,
                                     const Vec2d PtBundle,
                                     const double* Pt3d,
                                     const double PdsSqrt);

    private:

        const Mat3d   mRot0;
        const Vec2d   mPtBundle;
        const double *mPt3d;
        const double  mPdsSqrt;

};

class cResidualError
{
    public :
        cResidualError(const Mat3d aRot0,const Vec2d aPtBundle,const double aPdsSqrt) :
            mRot0(aRot0), mPtBundle(aPtBundle), mPdsSqrt(aPdsSqrt) {}
        ~cResidualError(){}

        // Parameters: pose and 3d points
        template <typename T>
        bool operator()(const T* const aW,
                        const T* const aC,
                        const T* const aPt3,
                        T* Residual) const;
        static CostFunction * Create(const Mat3d mRotCur,
                                     const Vec2d PtBundle,
                                     const double PdsSqrt);

    private:

        const Mat3d   mRot0;
        const Vec2d   mPtBundle;
        const double  mPdsSqrt;

};



#endif //_BUNDLE_H
