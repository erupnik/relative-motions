#include "all.h"
#include "cov_in_motions.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>

/*
  TODO
  - DONE + BuildProblem added - add optimization of triplets in cAppCovInMotion::OptimizeRelMotions()
  - DONE and verif ceres cov with your Hessian
  - DONE verify that your global similitude corresponds to the good trafo direction
  - DONE make sure the order in Hi and gi correponds to the good poses
  - DONE finish propagation test
  - DONE fix ReadSimGlob to update new variables
  - check inputs:
      + edges:
      + similarities: computed between edges and Ori
      + global poses: from input Ori
  - parallellize

 NOTES:
  - we are doing the final refinement so the pairs/triplets have to be without gross errors
    otherwise the global H will be contaminated with errors
  -

*/

/*   LOCAL BUNDLE ADJUSTEMENT per relative motion

      Fn(x) = I(Proj(pose(x)))  => collinearity 3D-2D
      I: pp, f (C)         3D-2D
      Proj: Z              3D-3D
                 I+proj  => intrinsic parameters!!!!
      pose: {r,t}        => extrinsic parameters!!!!!


      res = (Fn(P) - p)
      LOSS = res**2

      delta = dRt     pose update 6x1

      J_ = d_res/d_delta =
           d_Fn/d_delta = d_I /d_(Xq,Yq) @  d_Proj/d_(Xc,Yc,Zc) @ d_pose/d_(r,t)
                        = J_p2d_p3q      @  J_p3q_p3d           @ J_p3d_rt

      sizes:
        res:            N x 1
        J_:             N x 6    -> 1 pose, N observations
        J_:             N x 6M   -> M poses, N observations
        (J_t W J_)^-1:  6M x 6M  -> M poses, N observations
        J_t W res:      6M x 1   -> M poses

      Solution:
      delta_ = = H_^-1 @ g_ = (J_t W J_)^-1 J_t W res

      Update:
      [r+,t+] = exp(delta_)^T * [r | t]

*/

/*   GLOBAL BUNDLE ADJUSTEMENT all images with covariance propagation

     global to local transformation:
      pose = s * alpha * Pose + beta

    // reminder:
    // - pose and Pose are extrinsic parameters
    // - {s, alpha, beta} describe transformation from absolute 3D to relative 3D
    // - P and P_ are 3D pts in relative and absolute frames respectively
     res = ( Fn(P) - p )                                 //local
         =   I(Proj(          pose(P)           )) -p    //local
         =   I(Proj( s * alpha * Pose(P_) + beta   )) -p    //

         = I(Proj( simil( Pose(P_) ) )) => global   ->   local     ->   projection  -> calib
                                          (p3d^g)    s   (p3d^l)  proj    (p3q)     I  (p2d)

    derivatives (2months later),
         (alpha,beta) are constants:
           res' =  J_p2d_p3q @ J_p3q_p3d^l @ J_p3d^l_rt @ J_rt_RT
         (alpha,beta) are NOT constants:
           res' =  J_p2d_p3q @ J_p3q_p3d^l @ J_p3d^l_rt @ J_rt_simil @ J_simil_RT

    J = d_res/d_delta =
             d_Fn / d_delta = d_I/d_(Xq,Yq) @               //
                              d_Proj/d_(Xc,Yc,Zc) @         //
                              d_pose/d_(r,t) @              //============
                              d_Pose/d_(R,T)                //
                            = J_p2d_p3q @ J_p3q_p3d @ J_p3d_rt @ J_rt_RT
                            = J_ @ J_rt_RT


    sizes :
       res:          N x 1
       J:            N x 6    -> 1 pose, N observations
       J:            N x 6M   -> M poses, N observations
       (Jt W J)^-1:  6M x 6M  -> M poses, N observations
       Jt W res:     6M x 1   -> M poses
       J_rt_RT:      6M x 6M

    reformulate :
       hesian          H = (Jt W J) = (J_ @ J_rt_RT)^t @ W @ (J_ @ J_rt_RT)
                         = J_rt_RT^t @ (J_^t @ W @ J_) @ J_rt_RT)
                         = J_rt_RT^t @  H_   @ J_rt_RT
       gradient vector g = Jt W res = J_rt_RT^t @ (J_^t @ W @ res)
                         = J_rt_RT^t @ g_

       J_rt_R =  alpha * skew_symmetric(1,1,1)
              =  |  a11 a12  a13 |     |  0 -1   1 |
                 |  a21 a22  a23 |  @  |  1  0  -1 |
                 |  a31 a32  a33 |     | -1  1   0 |
       J_rt_T = s * alpha
              = s *  |  a11 a12  a13 |
                     |  a21 a22  a23 |
                     |  a31 a32  a33 |

      J_rt_RT = [ J_rt_R | J_rt_T ]
              = [ alpha @ [1,1,1]x   | s * alpha ]
              = [ |  a11 a12  a13 |                     |  a11 a12  a13 |
                  |  a21 a22  a23 | @ [1,1,1]x  |   s * |  a21 a22  a23 |
                  |  a31 a32  a33 |                     |  a31 a32  a33 | ]


    Solution:
       delta =  H^-1 @ g = (Jt W J)^-1 Jt W res
       size:
       delta: 6M x 1 =  (6M x 6M) @ (6M x 1)

    Update:
       [R+,T+] = exp(delta)^T * [R | T]

         ==============
    H - the big H will be a sum of the small Hi
                 image1                  image 2                        image 3 ....
          dTx, dTy, dTz, wx, wy wz    dTx, dTy, dTz, wx, wy wz  ....
      dTx
      dTy
      dTz
      wx
      wy
      wz
      .
      .
      .



         ==============

     LOSS = (x a x + b x + d) = x (Jt W J) x + p^t J x  + p^t * p
     minimize LOSS :
         LOSS'= ax + b = 0  => x = a^-1 @ -b






*/

double RandUnif_0_1()
{
   return ((double) rand() / RAND_MAX );
}

double RandUnif_C()
{
   return (RandUnif_0_1()-0.5) * 2.0;
}

/*********** POSES *******************/
Vec3d& cPose::C_()
{
    throw "C_() not implemented in the derived cPose class.";
    return mC_;
}

void cPose::Show() const
{
    std::cout << "============ " << mName << " ============" << "\n";
    std::cout << "R0\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        for (int aK2=0; aK2<3; aK2++)
            std::cout << "\t" << mR(aK1,aK2) << " ";

        std::cout << "\n" ;
    }
    std::cout << "Omega\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        std::cout << "\t" << mOmega[aK1] << " ";
    }
   /* std::cout << "\nCovariance of Omega\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        for (int aK2=0; aK2<3; aK2++)
            std::cout << "\t" << mCov_Omega[3*aK1+aK2] << " ";

        std::cout << "\n" ;
    }*/
    std::cout << "C \t" << mC[0] << " " << mC[1] << " " << mC[2] << "\n";
    /*std::cout << "Covariance of C\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        for (int aK2=0; aK2<3; aK2++)
            std::cout << "\t" << mCov_C[3*aK1+aK2] << " ";

        std::cout << "\n" ;
    }*/
}

void cPoseBasic::Show() const
{
    std::cout << "============ " << mName << " ============" << "\n";
    std::cout << "R\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        for (int aK2=0; aK2<3; aK2++)
            std::cout << "\t" << mR(aK1,aK2) << " ";

        std::cout << "\n" ;
    }

    std::cout << "\nC \t" << mC[0] << " " << mC[1] << " " << mC[2] << "\n";


}

/*********** VIEWS  *******************/
template <typename T, typename U>
bool cNviewPoseT<T,U>::propagate_cov()
{
    if (_COV_PROP)
      return true;

    Mat3d J_rt_T = affine_trafo.lambda * affine_trafo.alpha;

    Mat3d skew_sym;
    skew_sym << 0, -1, 1,
                1, 0, -1,
               -1, 1, 0;
    Mat3d J_rt_R = affine_trafo.alpha * skew_sym;

    T J_rt_RT = T::Zero();
    J_rt_RT.block(0,0,3,3) = J_rt_T;
    J_rt_RT.block(3,3,3,3) = J_rt_R;
    if (J_rt_RT.cols()==12)
    {
      J_rt_RT.block(6,6,3,3) = J_rt_T;
      J_rt_RT.block(9,9,3,3) = J_rt_R;
    }


    T JtHJ = J_rt_RT.transpose() * mHg->H_() * J_rt_RT;

    mHg->H_() = JtHJ;
    mHg->g_() = J_rt_RT.transpose() * mHg->g_();

    _COV_PROP = true;
    return _COV_PROP;
}

bool cNviewPoseX::propagate_cov()
{
    if (_COV_PROP)
      return true;

    Mat3d J_rt_T = affine_trafo.lambda * affine_trafo.alpha;

    //Mat3d skew_sym;
    /*skew_sym << 0, -1, 1,
                1, 0, -1,
               -1, 1, 0;*/
    /*skew_sym << 1, 0, 0,
                0, 1, 0,
                0, 0, 1;
    Mat3d J_rt_R = affine_trafo.alpha * skew_sym;*/
    Mat3d J_rt_R = J_rt_T;

    int rc = _GAUGE_FIRST_CAM_FIX ? 6*(NbView()-1) : 6*NbView();

    MatXd J_rt_RT = MatXd::Zero(rc,rc);
    J_rt_RT.block(0,0,3,3) = J_rt_T;
    J_rt_RT.block(3,3,3,3) = J_rt_R;
    if (J_rt_RT.cols()==12)
    {
      J_rt_RT.block(6,6,3,3) = J_rt_T;
      J_rt_RT.block(9,9,3,3) = J_rt_R;
    }
    else if (J_rt_RT.cols()==18)
    {

      J_rt_RT.block(6,6,3,3) = J_rt_T;
      J_rt_RT.block(9,9,3,3) = J_rt_R;

      J_rt_RT.block(12,12,3,3) = J_rt_T;
      J_rt_RT.block(15,15,3,3) = J_rt_R;

    }

    //std::cout << "J_rt_RT " << J_rt_RT.size() << ", mHg->H_() " << mHg->H_().size()
    //          << ", mHg->g_() " << mHg->g_().size() << "\n";
    MatXd JtHJ = J_rt_RT.transpose() * mHg->H_() * J_rt_RT;

    mHg->H_() = JtHJ;
    mHg->g_() = J_rt_RT.transpose() * mHg->g_();

    _COV_PROP = true;
    return _COV_PROP;
}

bool cAppCovInMotion::ReadFeatures()
{
    FILE* fptr = fopen(mfeats_file.c_str(), "r");
    if (fptr == NULL) {
      return false;
    };

    int OK=1;

    int aNbAll=0;
    while (!std::feof(fptr)  && !ferror(fptr))
    {

        int aNbPt;
        fscanf(fptr, "%i ", &aNbPt);
        //std::cout << aNbPt << "\n";

        std::string         aViewName="";
        std::vector<Vec2d>  aFeat2d;
        double            * aFeat3d = new double [3];

        //3d image coordinates
        for (int aK=0; aK<aNbPt; aK++)
        {
            //image name
            char aName[50];
            fscanf(fptr, "%s ", aName);
            //std::cout << aName << "\n";

            if (aK!=0)
              aViewName+="-";

            aViewName += aName;

            //pts coordinate
            Vec2d aPt;
            OK = fscanf(fptr, "%lf %lf", &aPt[0], &aPt[1]);
            aFeat2d.push_back(aPt);
            //std::cout <<  aPt << "\n";

        }
        //std::cout << "  " << OK << " " << aViewName << "\n";

        //3d coordinates
        OK = fscanf(fptr, "%lf %lf %lf", &(aFeat3d[0]),&(aFeat3d[1]),&(aFeat3d[2]));
        // std::cout << OK << " " << aFeat3d[0] << " " << aFeat3d[1] << " " << aFeat3d[2] << "\n";

        if (!DicBoolFind(mFeat3dMap,aViewName) && (OK>0))
            mFeat3dMap[aViewName] = new std::vector<double*>();

        if ((OK>0))
        {
            mFeat3dMap[aViewName]->push_back(aFeat3d);
            /*std::cout << aViewName << " " << aFeat3d[0] << " " << aFeat3d[1] << " " << aFeat3d[2] << "\n";
            std::vector<double*>* aTest = mFeat3dMap[aViewName];
            std::cout << " " << aTest->at(aTest->size()-1)[0] << " " << aTest->at(aTest->size()-1)[1] << " " << aTest->at(aTest->size()-1)[2] << "\n";
            getchar();*/
        }


        if (!DicBoolFind(mFeatViewMap,aViewName) && (OK>0))
        {
            mFeatViewMap[aViewName] = new std::vector<std::vector<Vec2d>>();
        }

        if ((OK>0))
          mFeatViewMap[aViewName]->push_back(aFeat2d);

            //ELISE_ASSERT((aNb==3) || (aNb==1),"Could not read 3 or 1 values");
    }
    fclose(fptr);


    //print out
    if (0)
    {
      for (auto aV : mFeatViewMap)
      {
          std::cout << aV.first << " ";

          for (auto aPt : *(aV.second))
          {
              for (auto aPtV : aPt)
              {
                  std::cout << aPtV[0] << " " << aPtV[1] << " ";
              }
              std::cout << "\n";
          }
      }
    }
    if (0)
    {
      for (auto aVPt : mFeat3dMap)
      {
          std::cout << aVPt.first << "\n";

          for (auto aPt : *(aVPt.second))
          {
              std::cout << "\t" << *aPt << " " << aPt[0] << " " << aPt[1] << " " << aPt[2] << "\n";
              getchar();
          }
      }

    }

    return true;

}

/*********** MANAGER CLASS  *******************/
bool cAppCovInMotion::ReadViews()
{
  //std::srand(std::time(nullptr));

  FILE* fptr = fopen(mviews_file.c_str(), "r");
  if (fptr == NULL) {
    return false;
  };

  int OK=1;

  while (!std::feof(fptr)  && !ferror(fptr))
  {

      int aNbV;
      fscanf(fptr, "%i", &aNbV);
      //std::cout << aNbV << "\n";
      //getchar();

      std::string              aViewName ="";
      std::vector<std::string> aPoseNameV;
      for (int aK=0; aK<aNbV; aK++)
      {
          char aName[50];
          fscanf(fptr, "%s", aName);
          aPoseNameV.push_back(aName);

          if (aK!=0)
            aViewName+="-";

          aViewName += aName;
          //std::cout <<  " " << aName << " aK=" << aK;
      }
      //std::cout << " aViewName " << aViewName << "\n";

      //first view
      double * aCId = new double[3]{0,0,0};
      //double * aCId = new double[3]{double(std::rand()) / RAND_MAX *0.1,
      //                              double(std::rand()) / RAND_MAX *0.1,
      //                              double(std::rand()) / RAND_MAX *0.1};
      //cPose aPose1(Mat3d::Identity(),aCId,aPoseNameV[0]);
      cPose* aPose1_ = new cPose(Mat3d::Identity(),aCId,aPoseNameV[0]);

      // second view
      Mat3d aRot21;
      double * aC21 = new double [3];
      for (int aK1=0; aK1<3; aK1++)
      {
          for (int aK2=0; aK2<3; aK2++)
          {
              OK = fscanf(fptr, "%lf", &aRot21(aK1,aK2));
          }
      }
      for (int aK1=0; aK1<3; aK1++)
      {
          OK = fscanf(fptr, "%lf", &aC21[aK1]);
      }


      if (OK>0)
      {
        if (aNbV==2) // two-view
        {
            //cPose aPose21(aRot21,aC21,aPoseNameV[1]);
            cPose* aPose21_ = new cPose(aRot21,aC21,aPoseNameV[1]);

            mAllViewMap_[aViewName] = new cNviewPoseX
                     (aPose1_,aPose21_,NULL,new cHessianGradientX(Mat6d::Zero(),Vec6d::Zero()));

        }
        else if (aNbV==3) //three-view?
        {
          // third view
          Mat3d aRot31;
          double * aC31 = new double [3];

          for (int aK1=0; aK1<3; aK1++)
          {
              for (int aK2=0; aK2<3; aK2++)
              {
                  OK = fscanf(fptr, "%lf", &aRot31(aK1,aK2));
                  //std::cout << aRot21(aK1,aK2) << " + ";
              }
          }
          for (int aK1=0; aK1<3; aK1++)
          {
              OK = fscanf(fptr, "%lf", &aC31[aK1]);
              //std::cout << aC31[aK1] << " TTT  ";
          }

          //cPose aPose21(aRot21,aC21,aPoseNameV[1]);
          //cPose aPose31(aRot31,aC31,aPoseNameV[2]);
          cPose* aPose21_ = new cPose(aRot21,aC21,aPoseNameV[1]);
          cPose* aPose31_ = new cPose(aRot31,aC31,aPoseNameV[2]);

          mAllViewMap_[aViewName] = new cNviewPoseX
                   (aPose1_,aPose21_,aPose31_,new cHessianGradientX(Mat12d::Zero(),Vec12d::Zero()));

        }
        else
        {
            std::cout << "Too many views. We currently support up to 3 views." << "\n";
            return 0;
        }
      }
    }

    return true;
}

bool cAppCovInMotion::ReadRotTrS(FILE* fptr,Mat3d& alpha,Vec3d& beta,double& s)
{
    bool OK;
    //rotation
    for (int aK1=0; aK1<3; aK1++)
    {
        for (int aK2=0; aK2<3; aK2++)
        {
            OK = fscanf(fptr, "%lf", &alpha(aK1,aK2));
            //std::cout << alpha(aK1,aK2) << " " << aK1 << " ";
        }
    }

    //translation
    for (int aK=0; aK<3; aK++)
    {
        OK = fscanf(fptr, "%lf", &beta(aK));
        //std::cout << beta(aK) << " ";
    }

    //scale
    OK = fscanf(fptr, "%lf", &s);
    //std::cout << s << " " << OK << " s";

    return OK;
}

bool cAppCovInMotion::ReadSimGlobal()
{
    std::cout << "ReadSimGlobal\n";
    //getchar();

    FILE* fptr = fopen(msimil_file.c_str(), "r");
    if (fptr == NULL) {
      return false;
    };

    int OK=1;

    while (!std::feof(fptr)  && !ferror(fptr))
    {

        int aNbV;
        fscanf(fptr, "%i", &aNbV);
        //std::cout << aNbV << "\n";
        //getchar();

        std::string              aViewName ="";
        std::vector<std::string> aPoseNameV;
        for (int aK=0; aK<aNbV; aK++)
        {
            char aName[50];
            fscanf(fptr, "%s", aName);
            aPoseNameV.push_back(aName);

            if (aK!=0)
              aViewName+="-";

            aViewName += aName;
            //std::cout <<  " " << aName << " aK=" << aK;
        }
        //std::cout <<  aViewName << "\n";

        if (DicBoolFind(mAllViewMap_,aViewName))
        {
            double& L = mAllViewMap_[aViewName]->lambda();
            Mat3d&  alpha = mAllViewMap_[aViewName]->alpha();
            Vec3d&  beta  = mAllViewMap_[aViewName]->beta();

            if (! ReadRotTrS(fptr,alpha,beta,L))
            {
                std::cout << "ERROR reading global similitude in cAppCovInMotion::ReadSimGlobal for "
                          << aViewName << "\n" ;
                return false;
            }

            //m2ViewMap_[aViewName]->PrintAlpha();
            //m2ViewMap_[aViewName]->PrintBeta();
            //getchar();

        }
        else
        {
            double i;
            for (int aK=0; aK<13; aK++)
            {
                fscanf(fptr, "%lf", &i);
            }
        }

      }
      return true;
}

bool cAppCovInMotion::ReadGlobalPoses()
{
    std::cout << "ReadGlobalPoses\n";
    //getchar();

    FILE* fptr = fopen(mglob_p_file.c_str(), "r");
    if (fptr == NULL) {
      return false;
    };

    int OK=1;
    while (!std::feof(fptr)  && !ferror(fptr))
    {
        Mat3d aR;
        double aT[3];

        char aName[50];
        OK = fscanf(fptr, "%s", aName);
        std::string PoseName = std::string(aName);

        for (int aK1=0; aK1<3; aK1++)
        {
            for (int aK2=0; aK2<3; aK2++)
            {
                OK = fscanf(fptr, "%lf", &aR(aK1,aK2));
                //std::cout << aR(aK1,aK2) << " " << aK1 << " ";
            }
        }
        for (int aK1=0; aK1<3; aK1++)
        {
            OK = fscanf(fptr, "%lf", &aT[aK1]);
            //std::cout << aT[aK1] << " " << aK1 << " ";
        }

        mGlobalPoses[PoseName] = new cPose(aR,aT,PoseName);
        //mGlobalPoses[PoseName]->Show();
        //getchar();

    }

    return true;
}

void cAppCovInMotion::PrintAllViews()
{
      std::cout << "Pairs and triplets:\n";
      for (auto Pose : mAllViewMap_)
      {
          Pose.second->View(0).Show();
          Pose.second->View(1).Show();
          if (Pose.second->NbView()==3)
            Pose.second->View(2).Show();

      }
}

void cAppCovInMotion::PrintAllPts3d()
{
    for (auto aVPt : mFeat3dMap)
    {
        std::cout << aVPt.first << "\n";

        for (auto aPt : *(aVPt.second))
        {
            std::cout << "\t" << *aPt << " " << aPt[0] << " " << aPt[1] << " " << aPt[2] << "\n";
            //getchar();
        }
    }
}


bool cAppCovInMotion::BuildProblem_(cNviewPoseX*& views,std::string views_name)
{
  /*bool gauge_first_cam_fix = false;
  bool gauge_base_fix = false;
  bool gauge_rappel_poses = true;*/

  // save initial structure
  //std::string aIniStruc = views_name + "_init.ply";
  //WriteToPLYFile(aIniStruc, views_name);

  //create ceres problem
  ceres::Problem * aProblem = new ceres::Problem;;
  ceres::Solver::Summary        aSummary;

  //covariance
  Covariance::Options aCovOpts;
  //aCovOpts.null_space_rank = 10;
  //aCovOpts.algorithm_type = DENSE_SVD;
  //aCovOpts.min_reciprocal_condition_number = 1.95146e-300;
  Covariance covariance(aCovOpts);
  std::vector<std::pair<const double*, const double*> > aCovBlocks;

  // get features
  std::vector<std::vector<Vec2d>>* aFeat2d = mFeatViewMap[views_name];

  int aNbPts = 0;
  if (aFeat2d)
    aNbPts = aFeat2d->size();
  else
  {
    VLOG(1) << "====>" << aNbPts << " features for " << views_name << "\n";
    return false;
  }

  // get points 3D
  std::vector<double*>* aFeat3d = mFeat3dMap[views_name];
  if (int(aFeat3d->size()) != aNbPts )
  {
      std::cout << "Inconsistency in image observations and 3D points. int(aPtVec->size()) != aNbPts";
      return false;
  }

  //create cost function per observation and AddResidualBlock
  std::vector<double*> param_bloc_CR;

  int NbCam = (&(views->View(2))==NULL) ? 2 : 3;
  for (int aCam=0; aCam<NbCam; aCam++)
  {
      //get EO
      double*  aC = views->View(aCam).C();
      double   aCPds = _T_PDS;

      Mat3d   aR0 = views->View(aCam).R();
      double*  aW = views->View(aCam).Omega();
      double   aWPds = _ROT_PDS;

      //residuals on features
      //std::cout << "features nb=" << aNbPts << "\n";
      for (int aK=0; aK<aNbPts; aK++)
      {
          //residuals on features
          if (1)
          {
              CostFunction * aCostF = cResidualError::Create( aR0,
                                                              aFeat2d->at(aK).at(aCam),
                                                              _FEAT_PDS);
              LossFunction * aLossF = new HuberLoss(1.00);
              aProblem->AddResidualBlock(aCostF,aLossF,
                                         aW, aC, aFeat3d->at(aK));
          }
          else
          {
              CostFunction * aCostF = cResidualOnPose::Create( aR0,
                                                               aFeat2d->at(aK).at(aCam),
                                                               aFeat3d->at(aK),
                                                               _FEAT_PDS);
              LossFunction * aLossF = new HuberLoss(1.0);
              aProblem->AddResidualBlock(aCostF,aLossF,
                                         aW, aC);
          }
      }

      // ===== BEGIN FIX GAUGE TO AVOID RANK DEFICIENCY
      //lock the first camera pose
      if (_GAUGE_FIRST_CAM_FIX)
      {
        if (aCam==0)
        {
          //perspective center
          aProblem->SetParameterBlockConstant(aC);
          //rotation
          aProblem->SetParameterBlockConstant(aW);

          //std::cout << "_GAUGE_FIRST_CAM_FIX" << "\n";
        }
      }
      //lock the base_x of the second camera
      if (_GAUGE_BASE_FIX)
      {
        if (aCam==1)
        {
          ceres::SubsetParameterization *constant_base_parameterization = NULL;
          std::vector<int> constant_translation;
          constant_translation.push_back(0);//x component

          constant_base_parameterization =
              new ceres::SubsetParameterization(3, constant_translation);

          aProblem->SetParameterization(aC, constant_base_parameterization);

          //std::cout << "_GAUGE_BASE_FIX" << "\n";

        }
      }
      // ===== END FIX GAUGE


      // "rappel" on the poses
      if (_GAUGE_RAPPEL_POSES)
      {
        if (_GAUGE_FIRST_CAM_FIX && aCam==0)
        {
            std::cout << "rappel pose cam0?\n" ;
            getchar();

        } // do nothing
        else
        {
            //perspective center
            CostFunction * aCostC = cPoseConstraint::Create(aC,aCPds);
            aProblem->AddResidualBlock(aCostC,NULL,(aC));
            param_bloc_CR.push_back(aC);

            //small rotation
            CostFunction * aCostR = cPoseConstraint::Create(aW,aWPds);
            aProblem->AddResidualBlock(aCostR,NULL,(aW));
            param_bloc_CR.push_back(aW);

            if (GET_COVARIANCES)
            {
                aCovBlocks.push_back(std::make_pair(aC,aC));
                aCovBlocks.push_back(std::make_pair(aW,aW));

            }
            //std::cout << "_GAUGE_RAPPEL_POSES" << "\n";
        }
      }

  }

    //solve the least squares problem
    ceres::Solver::Options aOpts;
    SetCeresOptions(aOpts);

    //ceres::Solve(aOpts,aProblem,&aSummary);
    //std::cout << aSummary.FullReport() << "\n";

    // show first camera
    //views->View(0).Show();

    //get covariances
    if (GET_COVARIANCES)
    {
      CHECK(covariance.Compute(aCovBlocks, aProblem));
      for (int aCam=1; aCam<NbCam; aCam++)
      {
          covariance.GetCovarianceBlock(views->View(aCam).C(), views->View(aCam).C(), views->View(aCam).Cov_C());
          covariance.GetCovarianceBlock(views->View(aCam).Omega(),views->View(aCam).Omega(), views->View(aCam).Cov_Omega());

          // show nth camera
          //views->View(aCam).Show();
      }
    }

    //get the jacobians
    double cost;// = aSummary.final_cost;//
    Problem::EvaluateOptions aEvalOpts;

    aEvalOpts.parameter_blocks = param_bloc_CR;
    std::vector<double> res;
    std::vector<double> JtWres;
    CRSMatrix *aJ = new CRSMatrix;

    bool eval_success = aProblem->Evaluate(aEvalOpts, &cost, &res, &JtWres, aJ);


    //hessian and gradient vector

    Eigen::SparseMatrix<double> A = MapJacToSparseM(aJ);
    Eigen::SparseMatrix<double> JtWJ = (A.transpose()*A);
    //std::cout << "Jacobian:\n" << A << "\n";

    //map to MatXd
    int JtWJ__size = _GAUGE_FIRST_CAM_FIX ? 6*(views->NbView()-1) : 6*views->NbView();
    //std::cout << "JtWJ__size " << JtWJ__size << " JtWJ " << JtWJ.rows() << " " <<  JtWJ.cols() << "\n";
    MatXd JtWJ_ = MatXd::Zero(JtWJ__size,JtWJ__size);

    MapJacToEigMat(JtWJ,JtWJ_,((_GAUGE_BASE_FIX) ? 1 : 0));//offset 1 bc tx is the constant base

    if (_GAUGE_BASE_FIX)
    {
      double grad_x_const [] = {0};
      JtWres.insert(JtWres.begin(),grad_x_const,grad_x_const+1);
    }

    Eigen::VectorXd my_vect = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(JtWres.data(), JtWres.size());//fixed-size to dynamic
    views->Hg_() = *(new cHessianGradientX(JtWJ_,my_vect));

    bool checkRank = false;
    if (checkRank)
    {

        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(views->Hg_().H_());
        std::cout << "\nThe rank of A is " << lu_decomp.rank()
                  << ", det=" << views->Hg_().H_().determinant() << "\n";

    }

    if (VLOG_IS_ON(2))
    {
        VLOG(2) << "Bundle adjustment SUCCESS: " << eval_success;
        VLOG(2) << "\n" << views_name
                << " _GAUGE_BASE_FIX=" << _GAUGE_BASE_FIX << ", "
                << " _GAUGE_FIRST_CAM_FIX=" <<  _GAUGE_FIRST_CAM_FIX << " "
                << "\nHHHHHHHHHHHHHHHHHHHHHHHh\n" << JtWJ_;
        //views->Hg_().printH();
        VLOG(2) << "\nggggggggggggggggggggggggggggg\n" << my_vect;
        //views->Hg_().printG();

        if (!eval_success)
            getchar();
    }
    //print residuals vector
    if (0)
    {
        std::cout << "residuals vector \n";
        for (int i=0; i<res.size(); i++)
        {
            std::cout << res[i] << ", ";
        }
        std::cout << "\n";
    }
    //covariance check
    if (GET_COVARIANCES)
    {
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
      solver.compute(JtWJ);
      Eigen::SparseMatrix<double> I(aJ->num_cols,aJ->num_cols);
      I.setIdentity();
      Eigen::SparseMatrix<double> JtWJ_inv = solver.solve(I);

      //print
      if (1)
      {
        std::cout << "hessian\n";
        for (int aK1=0; aK1<aJ->num_cols; aK1++)
        {
          for (int aK2=0; aK2<aJ->num_cols; aK2++)
          {
              std::cout << JtWJ.coeffRef(aK1,aK2) << " ";
          }
          std::cout << "\n";
        }
        std::cout << "covariance\n";
        for (int aK1=0; aK1<aJ->num_cols; aK1++)
        {
          for (int aK2=0; aK2<aJ->num_cols; aK2++)
          {
              std::cout << JtWJ_inv.coeffRef(aK1,aK2) << " ";
          }
          std::cout << "\n";
        }
        getchar();
      }
    }

    //save final structure
    //std::string aFinStruc = views_name + "_final.ply";
    //WriteToPLYFile(aFinStruc, views_name);

  return true;
}


void cAppCovInMotion::MapJacToEigMat(Eigen::SparseMatrix<double>& J, MatXd& JMat, int offset)
{
    for (int i = 0; i < J.rows(); i++)
    {
      for (int k = 0; k < J.cols(); k++)
      {
        JMat(i+offset, k+offset) = J.coeffRef(i,k);
      }
    }
}

Eigen::SparseMatrix<double> cAppCovInMotion::MapJacToSparseM(CRSMatrix* J)
{
    Eigen::SparseMatrix<double> A(J->num_rows, J->num_cols);
    //std::cout << J->num_rows << " " << J->num_cols << "\n";

    int row = 0;
    int el = 0;
    for (int i = 1; i < J->rows.size(); i++)
    {
      int NumInRow = J->rows[i] - J->rows[i - 1];
      for (int k = 0; k < NumInRow; k++)
      {
        int col = J->cols[el];

        A.coeffRef(row, col) = J->values[el++];
        //std::cout << A.coeffRef(row, col) << " ";
      }
      //std::cout << "\n" << i << "," << NumInRow << "\n";
      row++;
    }

    return A;
}

void cAppCovInMotion::SetMinimizer(ceres::Solver::Options& aSolOpt)
{
  //S - the reduced camera matrix / the Shur complement;
  //uses the SHURR trick; solves S as a dense matrix with Cholesky factorization; i
  //for problems up to several hundreds of cameras
  aSolOpt.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;
  //uses the SHURR trick; solves S as a sparse matrix with Cholesky factorization;
  //aSolOpt.linear_solver_type = ceres::SPARSE_SCHUR;
  //applies Preconditioned Conjugate Gradients to S; implements inexact step algorithm;
  //choose a precondition, e.g. CLUSTER_JACOBI,CLUSTER_TRIDIAGONAL that exploits camera-point visibility structure
  //mSolopt.linear_solver_type = ceres::ITERATIVE_SHUR;

  //relaxes the requirement to decrease the obj function at each iter step;
  //may turn very efficient in the long term;
  aSolOpt.use_nonmonotonic_steps = false;

  aSolOpt.max_num_iterations = 1;
  aSolOpt.minimizer_progress_to_stdout = true;
  /*aSolOpt.num_threads = 20;*/

  aSolOpt.use_inner_iterations = false;

}

void cAppCovInMotion::SetCeresOptions(ceres::Solver::Options& aSolOpt)
{
  SetMinimizer(aSolOpt);
  //SetOrdering();


}

bool cAppCovInMotion::OptimizeRelMotions()
{
    int Cnt = 0;

    for (auto a3v : mAllViewMap_)
    {
        std::cout << "=============== " << Cnt << " " << a3v.first << "===========================\n";
        if (!BuildProblem_(a3v.second,a3v.first))
		    {
                std::cout << " not added in the block because no observations" << "\n";
		        mAllViewMap_.erase(a3v.first);
		    }
        else
            Cnt++;
    }

    return true;

}

bool cAppCovInMotion::OptimizeGlobally()
{
      std::cout << "\nOptimize globally " <<  "\n";
    /*
      one iteration
        1- Propagate covariances
        2- Accumulate H, g
        3- Solve for delta
        4- Update global poses
		5- Save new poses
    */

    //stores poses' names and their index
    std::map<std::string,int> poses_in_motions;

    //1- Propagate covariances

    int pose_idx=0;
    //N-view
    for (auto a3v : mAllViewMap_)
    {
        std::cout << "=============== "  << a3v.first << "===========================\n";

        if (a3v.second->propagate_cov())
        {

          if (!DicBoolFind(poses_in_motions,a3v.second->View(0).Name()))
            poses_in_motions[a3v.second->View(0).Name()] = pose_idx++;

          if (!DicBoolFind(poses_in_motions,a3v.second->View(1).Name()))
            poses_in_motions[a3v.second->View(1).Name()] = pose_idx++;

          if (a3v.second->NbView()==3)
          {
              if (!DicBoolFind(poses_in_motions,a3v.second->View(2).Name()))
                  poses_in_motions[a3v.second->View(2).Name()] = pose_idx++;

          }

        }
    }

    if (VLOG_IS_ON(1))
    {
        for (auto pose : poses_in_motions)
            std::cout << "\n" << pose.first << " " << pose.second;
    }


    //2- Accumulate H, g
    if (poses_in_motions.size()>0)
    {
        int Hg_size = (poses_in_motions.size())*6;
        std::cout << "\nGlobal hessian size: " << Hg_size << "\n";

        Eigen::MatrixXd H_global = Eigen::MatrixXd::Zero(Hg_size,Hg_size);
        Eigen::VectorXd g_global = Eigen::VectorXd::Zero(Hg_size);

        for (auto a3v : mAllViewMap_)
        {
            int NumView = a3v.second->NbView();
            if (a3v.second->IS_COV_PROP())
            {
 	   		int start_cam_id = (_GAUGE_FIRST_CAM_FIX) ? 1 : 0;
                for (int vi=start_cam_id; vi<NumView; vi++)
                {
                    int glob_id = 6* poses_in_motions[a3v.second->View(vi).Name()];
                    int loc_id = (vi-start_cam_id)*6;
 	   			//std::cout << "VIEW3 " << glob_id << " " << loc_id << "\n";

                    //std::cout << a3v.second->View(vi).Name() << " "  << glob_id <<  "\n";
                    //std::cout << "loc_id=" << a3v.second->Hg_().g_().rows() << "/" << loc_id << "\n";

                    H_global.block(glob_id,glob_id,6,6) += a3v.second->Hg_().H_().block(loc_id,loc_id,6,6);
                    g_global.block(glob_id,0,6,1) += a3v.second->Hg_().g_().block(loc_id,0,6,1);

                    //std::cout << "H_global current\n" << H_global << "\n";

                  //a3v.second->Hg_().printH();

                }
            }
        }

        std::cout << "\nH_global\n" << H_global << "\n";
        VLOG(1) << "\nH_global\n" << H_global << "\n";
        VLOG(1) << "\ng_global\n" << g_global << "\n";

        Eigen::MatrixXd x = H_global.fullPivLu().solve(g_global);

        double relative_error = (H_global*x - g_global).norm() / g_global.norm();
        VLOG(1) << "\nThe realtive error is:\n" << relative_error << "\n";
        VLOG(1) << "\nThe determinant is:\n" << H_global.determinant() << "\n";
        VLOG(1) << "\nX:\n" << x << "\n";


        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(H_global);
        VLOG(1) << "\nThe rank of A is " << lu_decomp.rank() << "\n";

        /* Update poses and save */
        //mGlobalPoses[PoseName] = new cPose(aR,aT,PoseName)
        //std::string save_eo_name = "global_poses_new.txt";
        //std::ofstream eo_file(save_eo_name);

        for (auto pose : poses_in_motions)
        {
            int id_start = 6*pose.second;
            Vec3d dT(x(id_start,0),
                     x(id_start+1,0),
                     x(id_start+2,0));

            Mat3d dW;
            dW << 1,              -x(id_start+5,0), x(id_start+4,0),
                  x(id_start+5,0), 1,              -x(id_start+3,0),
                 -x(id_start+4,0), x(id_start+3,0), 1;


            cPose * pose0 = mGlobalPoses[pose.first];
            pose0->Show();

            std::cout << "dT=" << dT << ",\ndW" << dW << "\n";


            Mat3d newR = pose0->R()*dW;
            double * newT = pose0->C();
            newT[0] += dT[0];
            newT[1] += dT[1];
            newT[2] += dT[2];
            std::cout << "pose0_C=\n" << pose0->C()[0] << " "
                                   << pose0->C()[1] << " "
                                   << pose0->C()[2] << "\n";
            pose0->Show();
         //             << ", \nnewT=" << newT[0] << " " << newT[1] << " " << newT[2] << "\n";

            VLOG(1) << pose.first << " " << newR(0,0) << " " << newR(0,1) << " " << newR(0,2)
                                  << " " << newR(1,0) << " " << newR(1,1) << " " << newR(1,2)
                                  << " " << newR(2,0) << " " << newR(2,1) << " " << newR(2,2)
                                  << " " << pose0->C()[0]+dT[0]
                                  << " " << pose0->C()[1]+dT[1]
                                  << " " << pose0->C()[2]+dT[2] << "\n";
        }
        //eo_file.close();

    }
    else
        LOG(INFO) << "0 poses to propagate.";

    return true;
}

void cAppCovInMotion::WriteToPLYFile(const std::string& filename, const std::string& viewName)
{

  std::vector<double*>* aFeat3d = mFeat3dMap[viewName];

  if (aFeat3d)
  {
    int NbPts = int((aFeat3d)->size());
    std::ofstream of(filename.c_str());

    of << "ply"
     << '\n' << "format ascii 1.0"
     << '\n' << "element vertex " << 2 + NbPts
     << '\n' << "property float x"
     << '\n' << "property float y"
     << '\n' << "property float z"
     << '\n' << "property uchar red"
     << '\n' << "property uchar green"
     << '\n' << "property uchar blue"
     << '\n' << "end_header" << std::endl;

     // camera perspective centers
     if (DicBoolFind(mAllViewMap_,viewName))
     {

         double*  aC1  = mAllViewMap_[viewName]->View(0).C();
         double*  aC2  = mAllViewMap_[viewName]->View(1).C();

         of << aC1[0] << ' ' << aC1[1] << ' ' << aC1[2]
               << " 0 255 0" << '\n';

         of << aC2[0] << ' ' << aC2[1] << ' ' << aC2[2]
               << " 0 255 0" << '\n';

         if (mAllViewMap_[viewName]->NbView()==3)
         {
            double*  aC3  = mAllViewMap_[viewName]->View(2).C();
            of << aC3[0] << ' ' << aC3[1] << ' ' << aC3[2]
               << " 0 255 0" << '\n';
         }


     }


      // 3D structure
      for (auto aPt : (*aFeat3d))
      {
         of << aPt[0] << ' ' << aPt[1] << ' ' << aPt[2] << ' '
         << " 255 255 255" << '\n';
      }
      of.close();

  }
}

cAppCovInMotion::cAppCovInMotion(const std::string& avfile,
                                 const std::string& affile,
                                 const std::string& asimfile,
                                 const std::string& aglobpfile,
                                 const bool get_covariances) :
  mviews_file(avfile),
  mfeats_file(affile),
  msimil_file(asimfile),
  mglob_p_file(aglobpfile),
  GET_COVARIANCES(get_covariances)
{
    ReadFeatures();
    ReadViews();

    if (0)
      PrintAllViews();
    if (0)
      PrintAllPts3d();

    /* Get covariances per motion */
    OptimizeRelMotions();

    getchar();

    if (msimil_file!="")
      ReadSimGlobal();

    if (mglob_p_file!="")
      ReadGlobalPoses();

    /* Run global bundle adjustment using "relative" covariances */
    OptimizeGlobally();

}

int cov_in_motions_main(std::string views,
                        std::string tracks,
                        std::string similarities,
                        std::string global_poses,
                        bool get_covariances)
{



  cAppCovInMotion anAp(views,tracks,similarities,global_poses,get_covariances);
  //cAppCovInMotion anAp(aArg1,aArg2,aArg3,aArg4);


  return 1;
}
