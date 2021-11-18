#include "all.h"
#include "cov_in_motions.h"
#include <Eigen/IterativeLinearSolvers>

/*
  TODO
  - DONE + BuildProblem added - add optimization of triplets in cAppCovInMotion::OptimizeRelMotions()
  - DONE and verif ceres cov with your Hessian
  - DONE verify that your global similitude corresponds to the good trafo direction
  - make sure the order in Hi and gi correponds to the good poses
  - finish propagation test
  -
*/

/*   LOCAL BUNDLE ADJUSTEMENT per relative motion

      Fn(x) = I(Proj(pose(x)))  => collinearity 3D-2D
      I: pp, f (C)         3D-2D
      Proj: Z              3D-3D
      pose: R,t (rx+t,C)


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


     res = ( Fn(P) - p )                                 //local
         =   I(Proj(          pose(P)           )) -p    //local
         =   I(Proj(  alpha * Pose(P) + beta    )) -p    //global

     J = d_res/d_delta =
         d_Fn / d_delta = d_I/d_(Xq,Yq) @
                          d_Proj/d_(Xc,Yc,Zc) @
                          d_pose/d_(r,t) @
                          d_Pose/d_(R,T)
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
    H
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
    std::cout << "\nCovariance of Omega\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        for (int aK2=0; aK2<3; aK2++)
            std::cout << "\t" << mCov_Omega[3*aK1+aK2] << " ";

        std::cout << "\n" ;
    }
    std::cout << "C \t" << mC[0] << " " << mC[1] << " " << mC[2] << "\n";
    std::cout << "Covariance of C\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        for (int aK2=0; aK2<3; aK2++)
            std::cout << "\t" << mCov_C[3*aK1+aK2] << " ";

        std::cout << "\n" ;
    }
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

void cHessianGradient2::printH() const
{
  std::cout << "========================" << "\n";

  for (int aK1=0; aK1<6; aK1++)
  {
      for (int aK2=0; aK2<6; aK2++)
          std::cout << "\t" << mH(aK1,aK2) << " ";

      std::cout << "\n" ;
  }

}

void cHessianGradient3::printH() const
{
  std::cout << "========================" << "\n";

  for (int aK1=0; aK1<12; aK1++)
  {
      for (int aK2=0; aK2<12; aK2++)
          std::cout << "\t" << mH(aK1,aK2) << " ";

      std::cout << "\n" ;
  }

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

bool cAppCovInMotion::ReadViews()
{
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
      cPose aPose1(Mat3d::Identity(),aCId,aPoseNameV[0]);

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
            cPose aPose21(aRot21,aC21,aPoseNameV[1]);

            m2ViewMap[aViewName] = new c2viewPose(aPose1,aPose21);
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

          cPose aPose21(aRot21,aC21,aPoseNameV[1]);
          cPose aPose31(aRot31,aC31,aPoseNameV[2]);
          m3ViewMap[aViewName] = new c3viewPose(aPose1,aPose21,aPose31);

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
            std::cout << alpha(aK1,aK2) << " " << aK1 << " ";
        }
    }
    std::cout << OK << " a \n";
    //translation
    for (int aK=0; aK<3; aK++)
    {
        OK = fscanf(fptr, "%lf", &beta(aK));
        std::cout << beta(aK) << " ";
    }
    std::cout << OK << " b \n";
    //scale
    OK = fscanf(fptr, "%lf", &s);
    std::cout << s << " " << OK << " s";

    return OK;
}

bool cAppCovInMotion::ReadSimGlobal()
{
    std::cout << "ReadSimGlobal\n";
    getchar();

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
            std::cout <<  " " << aName << " aK=" << aK;
        }
        std::cout <<  aViewName << "\n";

        if (DicBoolFind(m2ViewMap,aViewName))
        {
            double& L = m2ViewMap[aViewName]->lambda();
            Mat3d&  alpha = m2ViewMap[aViewName]->alpha();
            Vec3d&  beta  = m2ViewMap[aViewName]->beta();

            if (! ReadRotTrS(fptr,alpha,beta,L))
            {
                std::cout << "ERROR reading global similitude in cAppCovInMotion::ReadSimGlobal for "
                          << aViewName << "\n" ;
                return false;
            }

            m2ViewMap[aViewName]->PrintAlpha();
            m2ViewMap[aViewName]->PrintBeta();
            getchar();

        }
        else if (DicBoolFind(m3ViewMap,aViewName))
        {
          double& L = m3ViewMap[aViewName]->lambda();
          Mat3d&  alpha = m3ViewMap[aViewName]->alpha();
          Vec3d&  beta  = m3ViewMap[aViewName]->beta();

          if (! ReadRotTrS(fptr,alpha,beta,L))
          {
              std::cout << "ERROR reading global similitude in cAppCovInMotion::ReadSimGlobal for "
                        << aViewName << "\n" ;
              return false;
          }

          m3ViewMap[aViewName]->PrintAlpha();
          m3ViewMap[aViewName]->PrintBeta();
          getchar();


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
    getchar();

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
                std::cout << aR(aK1,aK2) << " " << aK1 << " ";
            }
        }
        for (int aK1=0; aK1<3; aK1++)
        {
            OK = fscanf(fptr, "%lf", &aT[aK1]);
            std::cout << aT[aK1] << " " << aK1 << " ";
        }

        mGlobalPoses[PoseName] = new cPose(aR,aT,PoseName);
        //mGlobalPoses[PoseName]->Show();
        //getchar();

    }

    return true;
}

void cAppCovInMotion::PrintAllViews()
{
      std::cout << "Pairs:\n";
      for (auto Pose : m2ViewMap)
      {
          Pose.second->View(0).Show();
          Pose.second->View(1).Show();

      }
      std::cout << "Triplets:\n";
      for (auto Pose : m3ViewMap)
      {
          Pose.second->View(0).Show();
          Pose.second->View(1).Show();
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

bool cAppCovInMotion::BuildProblem(const std::string& aNameCam,int aNbCam)
{
  //retrieve the motion (pair or triplet)
  cNviewPose* view_i = 0;
  if (aNbCam==2)
    view_i = m2ViewMap[aNameCam];
  if (aNbCam==3)
    view_i = m3ViewMap[aNameCam];

  // save initial structure
  std::string aIniStruc = aNameCam + "_init.ply";
  WriteToPLYFile(aIniStruc, aNameCam);

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
  std::vector<std::vector<Vec2d>>* aFeat2d = mFeatViewMap[aNameCam];
  int aNbPts = aFeat2d->size();
  std::cout << aNameCam << " " << aNbPts << "\n";


  // get points 3D
  std::vector<double*>* aFeat3d = mFeat3dMap[aNameCam];
  if (int(aFeat3d->size()) != aNbPts )
  {
      std::cout << "Inconsistency in image observations and 3D points. int(aPtVec->size()) != aNbPts";
      return false;
  }

  //create cost function per observation and AddResidualBlock
  //for each camera
  ///std::vector<ResidualBlockId> aPoseResBlockId;
  std::vector<double*> param_bloc_CR;
  //Eigen::SparseMatrix<double> S_inv(aNbPts*2*3+(aNbCam-1)*6, aNbPts*2*3+(aNbCam-1)*6);

  for (int aCam=0; aCam<aNbCam; aCam++)
  {
      //get EO
      double*  aC = view_i->View(aCam).C();
      double   aCPds = _T_PDS;

      Mat3d   aR0 = view_i->View(aCam).R();
      double*  aW = view_i->View(aCam).Omega();
      double   aWPds = _ROT_PDS;

      //residuals on features
      std::cout << "features nb=" << aNbPts << "\n";
      for (int aK=0; aK<aNbPts; aK++)
      {
          //residuals on features
          if (1)
          {
              CostFunction * aCostF = cResidualError::Create( aR0,
                                                              aFeat2d->at(aK).at(aCam),
                                                              _FEAT_PDS);
              LossFunction * aLossF = new HuberLoss(1.0);
              aProblem->AddResidualBlock(aCostF,aLossF,
                                         aW, aC, aFeat3d->at(aK));


              //weight matrix
              //int pos_in_mat = aCam*aNbPts*2+2*aK + ((aCam>1) ? ((aCam-1)*6) : 0) ;
              //S_inv.coeffRef(pos_in_mat, pos_in_mat) = 1/(_FEAT_PDS);
              //S_inv.coeffRef(pos_in_mat+1, pos_in_mat+1) = 1/(_FEAT_PDS);
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
      if (aCam==0)
      {
        //perspective center
        aProblem->SetParameterBlockConstant(aC);
        //rotation
        aProblem->SetParameterBlockConstant(aW);
      }
      //lock the base_x of the second camera
      ceres::SubsetParameterization *constant_base_parameterization = NULL;
      if (aCam==1)
      {
          std::vector<int> constant_translation;
          constant_translation.push_back(0);//x component

          constant_base_parameterization =
              new ceres::SubsetParameterization(3, constant_translation);

          aProblem->SetParameterization(aC, constant_base_parameterization);
        }
      // ===== END FIX GAUGE


      // "rappel" on the poses (except for the 1st one which is constant)
      if (aCam!=0)
      {
          //perspective center
          CostFunction * aCostC = cPoseConstraint::Create(aC,aCPds);
          //ResidualBlockId aIdC =
          aProblem->AddResidualBlock(aCostC,NULL,(aC));
          //aPoseResBlockId.push_back(aIdC);
          param_bloc_CR.push_back(aC);

          //small rotation
          CostFunction * aCostR = cPoseConstraint::Create(aW,aWPds);
          //ResidualBlockId aIdR =
          aProblem->AddResidualBlock(aCostR,NULL,(aW));
          //aPoseResBlockId.push_back(aIdR);
          param_bloc_CR.push_back(aW);


          //weight matrix
          /*int pos_in_mat = (aCam+1)*aNbPts*2 + ((aCam>1) ? ((aCam-1)*6) : 0);

          S_inv.coeffRef(pos_in_mat, pos_in_mat) = 1/(_T_PDS);
          S_inv.coeffRef(pos_in_mat+1, pos_in_mat+1) = 1/(_T_PDS);
          S_inv.coeffRef(pos_in_mat+2, pos_in_mat+2) = 1/(_T_PDS);
          S_inv.coeffRef(pos_in_mat+3, pos_in_mat+3) = 1/(_ROT_PDS);
          S_inv.coeffRef(pos_in_mat+4, pos_in_mat+4) = 1/(_ROT_PDS);
          S_inv.coeffRef(pos_in_mat+5, pos_in_mat+5) = 1/(_ROT_PDS);*/

      }

      //request covariance computation
      if (aCam!=0)
      {
        aCovBlocks.push_back(std::make_pair(aC,aC));
        aCovBlocks.push_back(std::make_pair(aW,aW));
      }
    }

    //solve the least squares problem
    ceres::Solver::Options aOpts;
    SetCeresOptions(aOpts);

    ceres::Solve(aOpts,aProblem,&aSummary);
    std::cout << aSummary.FullReport() << "\n";

    // show first camera
    view_i->View(0).Show();

    //get covariances
    CHECK(covariance.Compute(aCovBlocks, aProblem));
    for (int aCam=1; aCam<aNbCam; aCam++)
    {

        covariance.GetCovarianceBlock(view_i->View(aCam).C(), view_i->View(aCam).C(), view_i->View(aCam).Cov_C());
        covariance.GetCovarianceBlock(view_i->View(aCam).Omega(),view_i->View(aCam).Omega(), view_i->View(aCam).Cov_Omega());

        // show nth camera
        view_i->View(aCam).Show();

    }

    //get the jacobians - testing of the cov extraced above
    double cost = aSummary.final_cost;// 0.0;
    Problem::EvaluateOptions aEvalOpts;
    //aEvalOpts.residual_blocks = aPoseResBlockId;
    aEvalOpts.parameter_blocks = param_bloc_CR;
    std::vector<double> evaluated_residuals;
    std::vector<double> gradient_vector;
    CRSMatrix *aJ = new CRSMatrix;

    //aProblem->Evaluate(Problem::EvaluateOptions(), &cost, &evaluated_residuals, &gradient_vector, aJ);
    aProblem->Evaluate(aEvalOpts, &cost, &evaluated_residuals, &gradient_vector, aJ);

    std::cout << "Rows: " << aJ->num_rows << ", cols: " << aJ->num_cols << "\n";

    //save final structure
    std::string aFinStruc = aNameCam + "_final.ply";
    WriteToPLYFile(aFinStruc, aNameCam);


    //hessian on all residual blocks
    //Eigen::MappedSparseMatrix which will let you wrap a CRS matrix directly
    Eigen::SparseMatrix<double> A = MapJacToSparseM(aJ);
    Eigen::SparseMatrix<double> JtJ = (A.transpose()*A); //

    //print the hessian
    if (1)
    {
      std::cout << "hessian\n";
      for (int aK1=0; aK1<aJ->num_cols; aK1++)
      {
        for (int aK2=0; aK2<aJ->num_cols; aK2++)
        {
            std::cout << JtJ.coeffRef(aK1,aK2) << ", ";
        }
        std::cout << "\n";
      }
    }

    //print gradient vector
    if (1)
    {
        std::cout << "gradient vector \n";
        for (int i=0; i<gradient_vector.size(); i++)
        {
            std::cout << gradient_vector[i] << ", ";
        }
        std::cout << "\n";
    }
    //print residuals vector
    if (0)
    {
        std::cout << "residuals vector \n";
        for (int i=0; i<evaluated_residuals.size(); i++)
        {
            std::cout << evaluated_residuals[i] << ", ";
        }
        std::cout << "\n";
    }
    //covariance check
    if (0)
    {
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
      solver.compute(JtJ);
      Eigen::SparseMatrix<double> I(aJ->num_cols,aJ->num_cols);
      I.setIdentity();
      Eigen::SparseMatrix<double> JtJ_inv = solver.solve(I);

      //print
      if (1)
      {
        std::cout << "covariance\n";
        for (int aK1=0; aK1<aJ->num_cols; aK1++)
        {
          for (int aK2=0; aK2<aJ->num_cols; aK2++)
          {
              std::cout << JtJ_inv.coeffRef(aK1,aK2) << " ";
          }
          std::cout << "\n";
        }
      }
    }

    return true;
}

Eigen::SparseMatrix<double> cAppCovInMotion::MapJacToSparseM(CRSMatrix* J)
{
    Eigen::SparseMatrix<double> A(J->num_rows, J->num_cols);
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
  aSolOpt.use_nonmonotonic_steps = true;

  aSolOpt.max_num_iterations = 1;
  aSolOpt.minimizer_progress_to_stdout = true;
  /*aSolOpt.num_threads = 20;*/

  aSolOpt.use_inner_iterations = true;

}

void cAppCovInMotion::SetCeresOptions(ceres::Solver::Options& aSolOpt)
{
  SetMinimizer(aSolOpt);
  //SetOrdering();


}

bool cAppCovInMotion::OptimizeRelMotions()
{

    //2-view movements
    for (auto a2v : m2ViewMap)
    {
        if (!BuildProblem(a2v.first,2))
          return false;
    }

    for (auto a3v : m3ViewMap)
    {
        //std::cout << a3v.first << "\n";
        if (!BuildProblem(a3v.first,3))
          return false;
    }

    return true;

}

bool cAppCovInMotion::OptimizeGlobally()
{
    /*
        r=r0(I+w)=r0+r0w, r0 is scalar so we can forget it

        given that, let's assume
        r=r0*w

        x = aX + b

        x_ = x + dx
        residual (x_ - x) = dx

        dx = d(aX + b) = a* dX



    */

    return true;
}

void cAppCovInMotion::WriteToPLYFile(const std::string& filename, const std::string& viewName)
{
    int NbPts = int((mFeat3dMap[viewName])->size());
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
     if (DicBoolFind(m2ViewMap,viewName))
     {
        double*  aC1  = m2ViewMap[viewName]->View(0).C();
        double*  aC2  = m2ViewMap[viewName]->View(1).C();

        of << aC1[0] << ' ' << aC1[1] << ' ' << aC1[2]
           << " 0 255 0" << '\n';

        of << aC2[0] << ' ' << aC2[1] << ' ' << aC2[2]
              << " 0 255 0" << '\n';

     }
     else if (DicBoolFind(m3ViewMap,viewName))
     {
         double*  aC1  = m3ViewMap[viewName]->View(0).C();
         double*  aC2  = m3ViewMap[viewName]->View(1).C();
         double*  aC3  = m3ViewMap[viewName]->View(2).C();

         of << aC1[0] << ' ' << aC1[1] << ' ' << aC1[2]
               << " 0 255 0" << '\n';

         of << aC2[0] << ' ' << aC2[1] << ' ' << aC2[2]
               << " 0 255 0" << '\n';

         of << aC3[0] << ' ' << aC3[1] << ' ' << aC3[2]
               << " 0 255 0" << '\n';
     }


      // 3D structure
      std::vector<double*>* aFeat3d = mFeat3dMap[viewName];

      for (auto aPt : (*aFeat3d))
      {
          of << aPt[0] << ' ' << aPt[1] << ' ' << aPt[2] << ' '
          << " 255 255 255" << '\n';

      }

      of.close();

}

cAppCovInMotion::cAppCovInMotion(const std::string& avfile,
                                 const std::string& affile,
                                 const std::string& asimfile,
                                 const std::string& aglobpfile) :
  mviews_file(avfile),
  mfeats_file(affile),
  msimil_file(asimfile),
  mglob_p_file(aglobpfile)
{
    ReadFeatures();
    ReadViews();

    if (1)
      PrintAllViews();
    if (0)
      PrintAllPts3d();

    /* Get covariances per motion */
    OptimizeRelMotions();

    if (msimil_file!="")
      ReadSimGlobal();

    if (mglob_p_file!="")
      ReadGlobalPoses();

    /* Run global bundle adjustment using "relative" covariances */
    OptimizeGlobally();

}

int cov_in_motions_main(int argc, char** argv)
{

  std::string aArg1 = std::string(argv[1]);
  std::string aArg2 = std::string(argv[2]);

  std::cout << "argc=" << argc << "\n";
  std::string aArg3 = "";
  if (argc>3)
    aArg3 = std::string(argv[3]);
  std::string aArg4 = "";
  if (argc>4)
    aArg4 = std::string(argv[4]);

  cAppCovInMotion anAp(aArg1,aArg2,aArg3,aArg4);


  return 1;
}
