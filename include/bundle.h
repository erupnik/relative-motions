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

class cCoordConstraint
{
    public:
        cCoordConstraint(const double* init_val, double val_pds):
            init_val(init_val), val_pds(val_pds) {}
        ~cCoordConstraint() {}

        template <typename T>
        bool operator()(const T* const obs_val, T* Residual) const {
              Residual[0] = (obs_val[0] - init_val[0])/val_pds;
              return true;
        }
        static CostFunction * Create(const double* coord_value, const double coord_pds){
          return  (new AutoDiffCostFunction<cCoordConstraint,1,3> (new cCoordConstraint(
            coord_value,coord_pds)));
      }


    private:
        const double* init_val;
        const double  val_pds;
};

/* Residuals weighted by the covariance
 * (with decomposition of the quadratic cost into a sum of linear terms)
 * - unknowns C, WR
 * - constants alpha, beta, lambda (similarity trafo) and eveyrthing related to covariances */
class cResidualOn3ViewsPoseDecomp
{
    public:
        cResidualOn3ViewsPoseDecomp(const Mat3d alpha, const Vec3d beta, const double L, 
                                   const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                   const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res) :
                                  m_beta(beta),
                                  m_L(L),
                                  m_rV(rV),
                                  m_R0V(R0V),
                                  m_Wi(Wi),
                                  m_Li(Li),
                                  m_Cstei(Cstei),
                                  m_tot_res(total_res) {}

        ~cResidualOn3ViewsPoseDecomp(){}

        template <typename T>
            bool operator()(const T* const C0, const T* const W0,
                            const T* const C1, const T* const W1,
                            const T* const C2, const T* const W2, T* Residual) const;
                                
        static CostFunction* Create(const Mat3d alpha, const Vec3d beta, const double L, 
                                    const std::vector<Mat3d> r0V, const std::vector<Mat3d> R0V,
                                    const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res);


            private:
        /* Constants */
        const Mat3d m_alpha;
        const Vec3d m_beta;
        const double m_L;
        const std::vector<Mat3d> m_rV;
        const std::vector<Mat3d> m_R0V;

        /* Covariance-related */
        const VecXd m_Wi;
        const MatXd m_Li;
        const VecXd m_Cstei;

        const double m_tot_res;

};
class cResidualOn2ViewsPoseDecomp
{
    public:
        cResidualOn2ViewsPoseDecomp(const Mat3d alpha, const double* beta, const double* lambda, 
                                   const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                   const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res) :
                                  m_alpha(alpha),
                                  m_beta(beta),
                                  m_lambda(lambda),
                                  m_rV(rV),
                                  m_R0V(R0V),
                                  m_Wi(Wi),
                                  m_Li(Li),
                                  m_Cstei(Cstei),
                                  m_tot_res(total_res) {}

        ~cResidualOn2ViewsPoseDecomp(){}

        template <typename T>
            bool operator()(const T* const C0, const T* const W0,
                            const T* const C1, const T* const W1, T* Residual) const;
                                
        static CostFunction* Create(const Mat3d alpha, const double* beta, const double* lambda, 
                                    const std::vector<Mat3d> r0V, const std::vector<Mat3d> R0V,
                                    const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res);


            private:
        /* Constants */
        const Mat3d m_alpha;
        const double* m_beta;
        const double* m_lambda;
        const std::vector<Mat3d> m_rV;
        const std::vector<Mat3d> m_R0V;

        /* Covariance-related */
        const VecXd m_Wi;
        const MatXd m_Li;
        const VecXd m_Cstei;
        
        const double m_tot_res;

};

class cResidualOn3ViewsPoseDecompLAB
{
    public:
        cResidualOn3ViewsPoseDecompLAB(const Mat3d alpha0, const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                   const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res) :
                                  m_alpha0(alpha0),
                                  m_rV(rV),
                                  m_R0V(R0V),
                                  m_Wi(Wi),
                                  m_Li(Li),
                                  m_Cstei(Cstei),
                                  m_tot_res(total_res) {}

        ~cResidualOn3ViewsPoseDecompLAB(){}

        template <typename T>
            bool operator()(const T* const C0, const T* const W0,
                            const T* const C1, const T* const W1,
                            const T* const C2, const T* const W2, 
                            const T* const alpha_beta_L, T* Residual) const;
                                
        static CostFunction* Create(const Mat3d alpha0, 
                                    const std::vector<Mat3d> r0V, const std::vector<Mat3d> R0V,
                                    const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res);


    private:
        /* Constants */
        const Mat3d m_alpha0;
        const std::vector<Mat3d> m_rV;
        const std::vector<Mat3d> m_R0V;

        /* Covariance-related */
        const VecXd m_Wi;
        const MatXd m_Li;
        const VecXd m_Cstei;
        
        const double m_tot_res;

};

class cResidualOn2ViewsPoseDecompLAB
{
    public:
        cResidualOn2ViewsPoseDecompLAB(const Mat3d alpha0, 
                                   const std::vector<Mat3d> rV, const std::vector<Mat3d> R0V,
                                   const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res) :
                                  m_alpha0(alpha0),
                                  m_rV(rV),
                                  m_R0V(R0V),
                                  m_Wi(Wi),
                                  m_Li(Li),
                                  m_Cstei(Cstei),
                                  m_tot_res(total_res) {}

        ~cResidualOn2ViewsPoseDecompLAB(){}

        template <typename T>
            bool operator()(const T* const C0, const T* const W0,
                            const T* const C1, const T* const W1, 
                            const T* const alpha_beta_L, T* Residual) const;
                                
        static CostFunction* Create(const Mat3d alpha0, 
                                    const std::vector<Mat3d> r0V, const std::vector<Mat3d> R0V,
                                    const VecXd Wi, const MatXd Li, const VecXd Cstei, const double total_res);


            private:
        /* Constants */
        const Mat3d m_alpha0;
        const std::vector<Mat3d> m_rV;
        const std::vector<Mat3d> m_R0V;

        /* Covariance-related */
        const VecXd m_Wi;
        const MatXd m_Li;
        const VecXd m_Cstei;

        const double m_tot_res;
};

/* Triplets: Basic adjustment without covariance propagation 
 * * Wr=0 */
class cResidualOn3ViewsPoseBasicLAB
{
    public:
        cResidualOn3ViewsPoseBasicLAB(const Mat3d alpha0, 
                                    const std::vector<double*> cV, 
                                    const std::vector<Mat3d> rV, 
                                    const std::vector<Mat3d> R0V,
                                    const double Pds_c,const double Pds_w) :
                                  m_alpha0(alpha0),
                                  m_cV(cV),
                                  m_rV(rV),
                                  m_R0V(R0V),
                                  mPdsSqrt_c(Pds_c),
                                  mPdsSqrt_w(Pds_w){}

        ~cResidualOn3ViewsPoseBasicLAB(){}

        template <typename T>
            bool operator()(const T* const C0, const T* const W0,
                            const T* const C1, const T* const W1,
                            const T* const C2, const T* const W2, 
                            const T* const alpha_beta_L, T* Residual) const;
                                
        static CostFunction* Create(const Mat3d alpha0, 
                                    const std::vector<double*> c0V,
                                    const std::vector<Mat3d> r0V, 
                                    const std::vector<Mat3d> R0V,
                                    const double             Pds_c,
                                    const double             Pds_w);


    private:
        /* Constants */
        const Mat3d m_alpha0;
        const std::vector<double*> m_cV;
        const std::vector<Mat3d> m_rV;
        const std::vector<Mat3d> m_R0V;
        
        const double mPdsSqrt_c;
        const double mPdsSqrt_w;

   
};

/* Pairs: Basic adjustment without covariance propagation 
 * * Wr=0 */
class cResidualOn2ViewsPoseBasicLAB
{
    public:
        cResidualOn2ViewsPoseBasicLAB(const Mat3d alpha0, 
                                    const std::vector<double*> cV, 
                                    const std::vector<Mat3d> rV, 
                                    const std::vector<Mat3d> R0V,
                                    const double Pds_c,const double Pds_w) :
                                  m_alpha0(alpha0),
                                  m_cV(cV),
                                  m_rV(rV),
                                  m_R0V(R0V),
                                  mPdsSqrt_c(Pds_c),
                                  mPdsSqrt_w(Pds_w){}

        ~cResidualOn2ViewsPoseBasicLAB(){}

        template <typename T>
            bool operator()(const T* const C0, const T* const W0,
                            const T* const C1, const T* const W1,
                            const T* const alpha_beta_L, T* Residual) const;
                                
        static CostFunction* Create(const Mat3d alpha0, 
                                    const std::vector<double*> c0V,
                                    const std::vector<Mat3d> r0V, 
                                    const std::vector<Mat3d> R0V,
                                    const double             Pds_c,
                                    const double             Pds_w);


    private:
        /* Constants */
        const Mat3d m_alpha0;
        const std::vector<double*> m_cV;
        const std::vector<Mat3d> m_rV;
        const std::vector<Mat3d> m_R0V;
        
        const double mPdsSqrt_c;
        const double mPdsSqrt_w;

   
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
        cResidualError(const Mat3d aRot0,const Vec2d aPtBundle,const double aPdsSqrt,const double Foc) :
            mRot0(aRot0), mPtBundle(aPtBundle), mPdsSqrt(aPdsSqrt), mFoc(Foc) {}
        ~cResidualError(){}

        // Parameters: pose and 3d points
        template <typename T>
        bool operator()(const T* const aW,
                        const T* const aC,
                        const T* const aPt3,
                        T* Residual) const;
        static CostFunction * Create(const Mat3d mRotCur,
                                     const Vec2d PtBundle,
                                     const double PdsSqrt,
                                     const double Foc);

    private:

        const Mat3d   mRot0;
        const Vec2d   mPtBundle;
        const double  mPdsSqrt;
        
        const double mFoc;
};



#endif //_BUNDLE_H
