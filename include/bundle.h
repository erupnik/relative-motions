#ifndef _BUNDLE_H
#define _BUNDLE_H

#include "all.h"

/* Ceres */
#include "ceres/ceres.h"
#include "ceres/types.h"
#include "ceres/rotation.h"

/* Eigen */
#include <Eigen/Sparse>

/* Logging */
#include "gflags/gflags.h"
#include "glog/logging.h"

using namespace ceres;


/* Mat3d aPds is the weight matrix (square root) */
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

class cResidualError
{
    public :
        cResidualError(Mat3d aRot0,Vec2d aPtBundle,double aPdsSqrt) :
            mRot0(aRot0), mPtBundle(aPtBundle), mPdsSqrt(aPdsSqrt) {}
        ~cResidualError(){}

        template <typename T>
        bool operator()(const T* const aW,
                        const T* const aC,
                        const T* const aPt3,
                        T* Residual) const;

        static CostFunction * Create(const Mat3d mRotCur,const Vec2d PtBundle,const double PdsSqrt);
        /*{
          return  (new AutoDiffCostFunction<cResidualError,2,3,3,3> (new cResidualError(mRotCur,PtBundle,PdsSqrt)));
        }*/

    private:

        Mat3d   mRot0;
        Vec2d   mPtBundle;
        double  mPdsSqrt;

};

#endif //_BUNDLE_H
