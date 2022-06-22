#include "cov_in_motions.h"

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

    old der
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

void cPose::ShowImmutable() const
{
    std::cout << "============ " << mName << " immutable ============" << "\n";
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
        std::cout << "\t" << mOmega_immutable[aK1] << " ";
    }
  
    std::cout << "C \t" << mC_immutable[0] << " " << mC_immutable[1] << " " << mC_immutable[2] << "\n";

    std::cout << "Delta Omega / delta C\n";
    for (int aK1=0; aK1<3; aK1++)
    {
        std::cout << "\t" << mOmega[aK1] - mOmega_immutable[aK1] << " \t";
    }
  
    std::cout << mC[0] - mC_immutable[0] << " " << mC[1] - mC_immutable[1] << " " << mC[2] - mC_immutable[2] << "\n";

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
/*template <typename T, typename U>
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

bool cNviewPoseX::propagate_cov(Mat3d& Rg0,Mat3d& Rg1,Mat3d& Rg2)
{
    if (_COV_PROP)
      return true;

    Mat3d id_skew_sym;
    id_skew_sym << 1, -1, 1,
                   1, 1, -1,
                   -1, 1, 1;

    // R = alpha.inverse() * r;
    // alpha * R = r 
    //C = 1.0/lambda * alpha.inverse() * (c - beta);
    //lambda * alpha * C  + beta = c   
    Mat3d J_rt_T = affine_trafo.lambda * affine_trafo.alpha;
    // r = alpha * R , global to local
    // r0 (w+I) = alpha R0 (W+I) 
    //     w+I = r0**-1 alpha * R0 * (I+W)
    //     w   = r0**-1 alpha * R0 * (I+W) -I
    //dw/dW = r0**-1 * alpha * R0* (Id+skew) 
    Mat3d J_rt_R0 = (mView1->R()).inverse() * affine_trafo.alpha * Rg0 * id_skew_sym ;
    Mat3d J_rt_R1 = (mView21->R()).inverse() * affine_trafo.alpha * Rg1 *id_skew_sym ;
    Mat3d J_rt_R2;
    if (NbView()==3) J_rt_R2 = (mView31->R()).inverse() * affine_trafo.alpha * Rg2 * id_skew_sym;


    int rc = _GAUGE_FIRST_CAM_FIX ? 6*(NbView()-1) : 6*NbView();
    MatXd J_rt_RT = MatXd::Zero(rc,rc);
    J_rt_RT.block(0,0,3,3) = J_rt_T;
    J_rt_RT.block(3,3,3,3) = _GAUGE_FIRST_CAM_FIX ? J_rt_R1 : J_rt_R0;
    if (J_rt_RT.cols()==12)
    {
      J_rt_RT.block(6,6,3,3) = J_rt_T;
      J_rt_RT.block(9,9,3,3) = _GAUGE_FIRST_CAM_FIX ? J_rt_R2 : J_rt_R1;
    }
    else if (J_rt_RT.cols()==18)
    {

      J_rt_RT.block(6,6,3,3) = J_rt_T;
      J_rt_RT.block(9,9,3,3) = J_rt_R1;

      J_rt_RT.block(12,12,3,3) = J_rt_T;
      J_rt_RT.block(15,15,3,3) = J_rt_R2;

    }

    //std::cout << "J_rt_RT " << J_rt_RT.size() << ", mHg->H_() " << mHg->H_().size()
    //          << ", mHg->g_() " << mHg->g_().size() << "\n";
    MatXd JtHJ = J_rt_RT.transpose() * mHg->H_() * J_rt_RT;

    mHg->H_() = JtHJ;
    mHg->g_() = J_rt_RT.transpose() * mHg->g_();

    _COV_PROP = true;
    return _COV_PROP;
}*/

bool cAppCovInMotion::ReadFeatures(std::string tracks_file)
{
    FILE* fptr = fopen(tracks_file.c_str(), "r");
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
              aViewName+=_DELIMITER_;

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

    std::cout << " #" << mFeatViewMap.size() << " feature tracks\n";

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

    return EXIT_SUCCESS;

}

/*********** MANAGER CLASS  *******************/

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

void cAppCovInMotion::WriteLocalCovs()
{
    for(auto vi : mTriSet->mAllViewMap) 
    {
        std::string view_name = vi.first;

        //save covariance diagonal
        std::string covR_file_name = "covR" + m_lba_opts._KEY + ".txt";
        std::fstream covR_file;
        covR_file.open(covR_file_name.c_str(), std::ios_base::app);
        std::string covC_file_name = "covC" + m_lba_opts._KEY + ".txt";
        std::fstream covC_file;
        covC_file.open(covC_file_name.c_str(), std::ios_base::app);
       
       
        covR_file << vi.second->NbView() << " "; 
        covC_file << vi.second->NbView() << " "; 
        std::vector<std::string> vname_decom = mTriSet->DecompViewNames(view_name);
        for (auto vn : vname_decom)
        {
            covR_file << vn << " " ; 
            covC_file << vn << " " ; 
        }

        MatXd H = vi.second->Hg_().H_();
        for (int aK1=0; aK1<H.rows(); aK1++)
        { 
            for (int aK2=0; aK2<H.cols(); aK2++)
            {
                if (aK1==aK2)
                {
                    int id_block = std::floor(double(aK1)/3);
                 
                    if (id_block % 2 == 0) // if even -> C
                        covC_file << H(aK1,aK2) << " ";
                    else // if odd -> R
                        covR_file << H(aK1,aK2) << " ";
                }
            }
        }
        covR_file << "\n";
        covR_file.close();
        covC_file << "\n";
        covC_file.close();

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

      Mat3d   aR0 = views->View(aCam).R();
      double*  aW = views->View(aCam).Omega();

      //residuals on features
      //std::cout << "features nb=" << aNbPts << "\n";
      for (int aK=0; aK<aNbPts; aK++)
      {
          //residuals on features
          if (1)
          {
              CostFunction * aCostF = cResidualError::Create( aR0,
                                                              aFeat2d->at(aK).at(aCam),
                                                              m_lba_opts._FEAT_PDS);
              LossFunction * aLossF = new HuberLoss(m_lba_opts._HUBER);
              aProblem->AddResidualBlock(aCostF,NULL,
                                         aW, aC, aFeat3d->at(aK));
          }
          else
          {
              CostFunction * aCostF = cResidualOnPose::Create( aR0,
                                                               aFeat2d->at(aK).at(aCam),
                                                               aFeat3d->at(aK),
                                                               m_lba_opts._FEAT_PDS);
              LossFunction * aLossF = new HuberLoss(m_lba_opts._HUBER);
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
	    ceres::SubsetManifold* constant_base_manifold = nullptr;
            std::vector<int> constant_translation;
            constant_translation.push_back(0);//x component
	    
	    constant_base_manifold = new
		    ceres::SubsetManifold(3, constant_translation);
          
	    aProblem->SetManifold(aC, constant_base_manifold);
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
            double*  aC_immu = views->View(aCam).C_immutable();
            CostFunction * aCostC = cPoseConstraint::Create(aC_immu,m_lba_opts._C_PDS);
            aProblem->AddResidualBlock(aCostC,NULL,(aC));
            param_bloc_CR.push_back(aC);

            //small rotation
            double*  aW_immu = views->View(aCam).Omega_immutable();
            CostFunction * aCostR = cPoseConstraint::Create(aW_immu,m_lba_opts._ROT_PDS);
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
    //ceres::Solver::Options aOpts;
    //SetCeresOptions(aOpts);

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
    double cost=0;// = aSummary.final_cost;//
    Problem::EvaluateOptions aEvalOpts;
    aEvalOpts.apply_loss_function = true;

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

    // return if trace is over the threshold
    if (JtWJ_.trace() > m_lba_opts._TRACE_H)
    {
        std::vector<std::string> vname_decom = mTriSet->DecompViewNames(views_name);

        std::cout << " -Too large value of Hessian trace for ";
        for (auto i_n : vname_decom)
            std::cout << i_n << " ";  

        std::cout << "(" << JtWJ_.trace() <<  ")\n";
        
        return false;

    }

    if (_GAUGE_BASE_FIX)
    {
      double grad_x_const [] = {0};
      JtWres.insert(JtWres.begin(),grad_x_const,grad_x_const+1);
    }

    Eigen::VectorXd my_vect = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(JtWres.data(), JtWres.size());//fixed-size to dynamic
   

    //take square root of the hessian 
    views->Hg_() = *(new cHessianGradientX(JtWJ_.sqrt(),my_vect));
    //views->Hg_() = *(new cHessianGradientX(JtWJ_,my_vect));

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
                << "\nHHHHHHHHHHHHHHHHHHHHHHHh\n" << JtWJ_.sqrt();
        //views->Hg_().printH();
        //VLOG(2) << "\nggggggggggggggggggggggggggggg\n" << my_vect;
        //views->Hg_().printG();

        if (!eval_success)
            getchar();
    }
    //print residuals vector
    if (0)
    {
        //double* cmut = views->View(1).C();
        //double* cimmut = views->View(1).C_immutable();
        //VLOG(2) << "C =" << cmut[0] << " " << cmut[1] << " " << cmut[2] << 
        //           ", imm=" << cimmut[0] << " " << cimmut[1] << " " << cimmut[2] << "\n";
        
        VLOG(2) << "residuals vector \n";
        for (int i=0; i<res.size(); i++)
        {
            VLOG(2) << res[i] << ", ";
        }
        VLOG(2) << "\n";
    }
    //covariance check
    if (GET_COVARIANCES)
    {
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
      solver.compute(JtWJ);
      Eigen::SparseMatrix<double> I(aJ->num_cols,aJ->num_cols);
      I.setIdentity();
      Eigen::SparseMatrix<double> JtWJ_inv = solver.solve(I);

      //save covariance diagonal
      std::string covR_file_name = "covR" + m_lba_opts._KEY + ".txt";
      std::fstream covR_file;
      covR_file.open(covR_file_name.c_str(), std::ios_base::app);
      std::string covC_file_name = "covC" + m_lba_opts._KEY + ".txt";
      std::fstream covC_file;
      covC_file.open(covC_file_name.c_str(), std::ios_base::app);


      covR_file << views->NbView() << " "; 
      covC_file << views->NbView() << " "; 
      std::vector<std::string> vname_decom = mTriSet->DecompViewNames(views_name);
      for (auto vn : vname_decom)
      {
          covR_file << vn << " " ; 
          covC_file << vn << " " ; 
      }

      VLOG(2) << views_name << " ";
      MatXd JtWJ_sqrt = JtWJ_.sqrt();
      for (int aK1=0; aK1<JtWJ_sqrt.rows(); aK1++)
      {
          for (int aK2=0; aK2<JtWJ_sqrt.cols(); aK2++)
          {
              if (aK1==aK2)
              {
                  int id_block = std::floor(double(aK1)/3);
               
                  if (id_block % 2 == 0) // if even -> C
                      covC_file << JtWJ_sqrt(aK1,aK2) << " ";
                  else // if odd -> R
                      covR_file << JtWJ_sqrt(aK1,aK2) << " ";

                 VLOG(2) << JtWJ_sqrt(aK1,aK2) << " ";

              }
          }
      }
      VLOG(2) << "\n";
      covR_file << "\n";
      covR_file.close();
      covC_file << "\n";
      covC_file.close();

      //save the gradient vector  my_vect
      std::string gradR_file_name = "gradR" + m_lba_opts._KEY + ".txt";
      std::fstream gradR_file;
      gradR_file.open(gradR_file_name.c_str(), std::ios_base::app);
      std::string gradC_file_name = "gradC" + m_lba_opts._KEY + ".txt";
      std::fstream gradC_file;
      gradC_file.open(gradC_file_name.c_str(), std::ios_base::app);


      gradR_file << views->NbView() << " "; 
      gradC_file << views->NbView() << " "; 
      for (auto vn : vname_decom)
      {
          gradR_file << vn << " " ; 
          gradC_file << vn << " " ; 
      }

      VLOG(2) << "gradients \n";
      for (int aK=0; aK<my_vect.size(); aK++)
      {
          int id_block = std::floor(double(aK)/3);

          if (id_block % 2 == 0) // if even -> C
              gradC_file << my_vect[aK] << " ";
          else // if odd -> R
              gradR_file << my_vect[aK] << " ";
                 
          VLOG(2) << my_vect[aK] << " ";

      }
      VLOG(2) << "\n";

      gradR_file << "\n";
      gradR_file.close();
      gradC_file << "\n";
      gradC_file.close();


      //print hessian 
      if (0)
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
Eigen::SparseMatrix<double> cAppCovInMotion::MatToSparseM(Eigen::MatrixXd& M)
{
    int r = M.rows();
    int c =  M.cols();
    Eigen::SparseMatrix<double> A(r, c);

    for (int i = 0; i < r; i++)
    {
      for (int k = 0; k < c; k++)
      {
        A.coeffRef(i, k) = M(i,k);
        //std::cout << A.coeffRef(i, k) << " " ;
      }
    }

    return A;

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

void cAppCovInMotion::SetMinimizerLocal(ceres::Solver::Options& aSolOpt)
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

  aSolOpt.max_num_iterations = 100;
  aSolOpt.minimizer_progress_to_stdout = true;
  aSolOpt.num_threads = _PROC_COUNT;

  aSolOpt.use_inner_iterations = m_lba_opts._INNER_ITER;

}

//final, global least-squares simultaneously on all poses  
void cAppCovInMotion::SetMinimizerGlobal(ceres::Solver::Options& aSolOpt)
{
  //S - the reduced camera matrix / the Shur complement;
  //uses the SHURR trick; solves S as a dense matrix with Cholesky factorization; i
  //for problems up to several hundreds of cameras
  aSolOpt.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;
  //uses the SHURR trick; solves S as a sparse matrix with Cholesky factorization;
  //aSolOpt.linear_solver_type = ceres::SPARSE_SCHUR;
  //applies Preconditioned Conjugate Gradients to S; implements inexact step algorithm;
  //choose a precondition, e.g. CLUSTER_JACOBI,CLUSTER_TRIDIAGONAL that exploits camera-point visibility structure
  //aSolOpt.linear_solver_type = ceres::ITERATIVE_SCHUR;

  //relaxes the requirement to decrease the obj function at each iter step;
  //may turn very efficient in the long term;
  aSolOpt.use_nonmonotonic_steps = false;

  aSolOpt.max_num_iterations = 100;
  aSolOpt.minimizer_progress_to_stdout = true;
  aSolOpt.num_threads = _PROC_COUNT;

  aSolOpt.use_inner_iterations = m_gba_opts._INNER_ITER;

}


void cAppCovInMotion::InitCovariances(cNviewPoseX*& views,std::string views_name)
{
    int NbCam = (&(views->View(2))==NULL) ? 2 : 3;
    
    // Hessian is an Id matrix, no covariances 
    MatXd H = MatXd::Identity(NbCam*6,NbCam*6);
    VecXd g = VecXd::Zero(NbCam*6); 
    
    views->Hg_() = *(new cHessianGradientX(H,g));
    std::cout << "init ";
}

bool cAppCovInMotion::OptimizeRelMotions()
{
    std::cout << "Local view optimization\n";

    
    std::vector<std::string> views_to_del_by_key;

    #pragma omp parallel num_threads(_PROC_COUNT) 
    {
        #pragma omp for
        for(int vi=0; vi<mTriSet->mAllViewMap.size(); vi++) 
        {
            auto vi_it = mTriSet->mAllViewMap.begin();
            std::advance(vi_it, vi);

            if (m_lba_opts._RUN_PROP)
            {
                if (! BuildProblem_((*vi_it).second,(*vi_it).first))
                    views_to_del_by_key.push_back((*vi_it).first);
            }
            else 
                InitCovariances((*vi_it).second,(*vi_it).first);
        }

    }

    //delete views with no features 
    for (auto vi : views_to_del_by_key)
        mTriSet->mAllViewMap.erase(vi);

    std::cout << " # " << views_to_del_by_key.size() << " view removed\n";

    // save local covariances to file if requested 
    if (m_lba_opts._WRITE_COV)
        WriteLocalCovs();

    return EXIT_SUCCESS;

}

bool cAppCovInMotion::OptimizeRelMotionsGlobally()
{

    std::cout << "Global view optimization\n";
    
    ceres::Problem * aProblem = new ceres::Problem;;
    ceres::Solver::Summary        aSummary;


    //add equation for each view and camera 
    for(auto view_i : mTriSet->mAllViewMap)
    {

        //similitude
        Mat3d alpha0 = view_i.second->alpha0();
        double* Walpha = view_i.second->Walpha();
        double* beta  = view_i.second->beta();
        double* lambda = view_i.second->lambda();


        int NbCam = (&(view_i.second->View(2))==NULL) ? 2 : 3;
        for (int aCam=0; aCam<NbCam; aCam++)
        {
            //local pose 
            Vec3d c (view_i.second->View(aCam).C()[0], 
                     view_i.second->View(aCam).C()[1],
                     view_i.second->View(aCam).C()[2]);
            Mat3d r =  view_i.second->View(aCam).R();
            //double*  aW = view_i->View(aCam).Omega(); to consider including it

            //initial global pose
            cPose*& pose = mTriSet->Pose(view_i.second->View(aCam).Name());
            //constant 
            Mat3d R0 = pose->R();
            //parameters 
            double* C = pose->C();
            double* W = pose->Omega();
           
            //covariance 
            int loc = aCam*6;
            Mat6d cov = view_i.second->Hg_().H_().block(loc,loc,6,6);

    
            /*CostFunction * aCost = 
                cResidualOnViewPose::Create(alpha0,beta,lambda,r,c,R0,cov);

            LossFunction * aLoss = new HuberLoss(m_gba_opts._HUBER_S);
            
            aProblem->AddResidualBlock(aCost,aLoss,C,W);*/

            CostFunction * aCost = 
                cResidualOnViewPoseAffFree::Create(alpha0,r,c,R0,cov);

            LossFunction * aLoss = new HuberLoss(m_gba_opts._HUBER_S);
            aProblem->AddResidualBlock(aCost,aLoss,C,W,Walpha,beta,lambda);
            
            //flag as refined 
            pose->SetRefined();
        }
    }

    //add constraint on initial global poses 
    std::vector<ceres::ResidualBlockId> res_bloc_CR;
    for (auto pose_i : mTriSet->mGlobalPoses)
    {
        if (pose_i.second->IsRefined())
        {
            //parameters 
            double* C = pose_i.second->C();
            double* C_immutable = pose_i.second->C_immutable();
            double* W = pose_i.second->Omega();
            double* W_immutable = pose_i.second->Omega_immutable();
        
            LossFunction * aLossW = new HuberLoss(m_gba_opts._HUBER_P);
            LossFunction * aLossC = new HuberLoss(m_gba_opts._HUBER_P);
            
            //perspective center
            CostFunction * CostC = cPoseConstraint::Create(C_immutable,m_gba_opts._C_PDS);
            ceres::ResidualBlockId res_C_id = aProblem->AddResidualBlock(CostC,aLossC,(C));
            res_bloc_CR.push_back(res_C_id);

            //small rotation
            CostFunction * CostR = cPoseConstraint::Create(W_immutable,m_gba_opts._ROT_PDS);
            ceres::ResidualBlockId res_W_id = aProblem->AddResidualBlock(CostR,aLossW,(W));
            res_bloc_CR.push_back(res_W_id);
        }
    }
    //TODO add constraint on initial similarity 
    //add immutable parameters


    //solve the least squares problem
    ceres::Solver::Options Opts;
    SetMinimizerGlobal(Opts);

    ceres::Solve(Opts,aProblem,&aSummary);
    std::cout << aSummary.FullReport() << "\n";

    if (VLOG_IS_ON(1))
    {
        //Print residuals 
        //FILE* file_residual = fopen("residuals.csv", "w+");

        double cost=aSummary.final_cost;
        Problem::EvaluateOptions eval_opts;
        eval_opts.apply_loss_function = false;
 
        eval_opts.residual_blocks = res_bloc_CR;
        std::vector<double> residuals;
 
        bool eval_success = aProblem->Evaluate(eval_opts, &cost, &residuals, nullptr, nullptr);

        double res_total=0;
        double res_max=0;
        int cnt=0;
        for (auto res : residuals)
        {
            res_total+=std::abs(res);
            if (std::abs(res) > res_max)
                res_max = std::abs(res);


            //fprintf(file_residual, "%7.9f \n", res);
            //fflush(file_residual); 

            VLOG(1) << " " << res << ", scaled=" << ((cnt % 2) ? (res * m_gba_opts._C_PDS) : (res * m_gba_opts._ROT_PDS))  << " " 
                           << ((cnt % 2) ? (m_gba_opts._C_PDS) : (m_gba_opts._ROT_PDS)) << "\n";
            
            cnt++;
        }
        //fclose(file_residual);

        VLOG(1) << "\nAverage residual: " << res_total/residuals.size() 
                <<  ", Max residual: " << res_max << "\n";
    }

    /* Update the similarity transformation */
    mTriSet->UpdateAllAffine();

    return EXIT_SUCCESS;
}

//obsolete
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
/*
    //stores poses' names and their index
    std::map<std::string,int> poses_in_motions;

    //1- Propagate covariances

    int pose_idx=0;
    Mat3d Rg3Dummy = Mat3d::Identity();

    //N-view
    for (auto a3v : mTriSet->mAllViewMap)
    {
        //std::cout << "=============== "  << a3v.first << "===========================\n";

        if (a3v.second->propagate_cov(mTriSet->mGlobalPoses[poses_in_motions,a3v.second->View(0).Name()]->R(),
                                      mTriSet->mGlobalPoses[poses_in_motions,a3v.second->View(1).Name()]->R(),
                                      (a3v.second->NbView()==3) ? 
                                      mTriSet->mGlobalPoses[poses_in_motions,a3v.second->View(2).Name()]->R() : Rg3Dummy)); 
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

        for (auto a3v : mTriSet->mAllViewMap)
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

        VLOG(1) << "\nH_global\n" << H_global << "\n";
        VLOG(1) << "\ng_global\n" << g_global << "\n";

        Eigen::MatrixXd x = H_global.fullPivLu().solve(-g_global);
    
        double relative_error = (H_global*x - g_global).norm() / g_global.norm();
        VLOG(1) << "\nThe realtive error is:\n" << relative_error << "\n";
        VLOG(1) << "\nThe determinant is:\n" << H_global.determinant() << "\n";
        VLOG(1) << "\nX:\n" << x << "\n";


        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(H_global);
        VLOG(1) << "\nThe rank of A is " << lu_decomp.rank() << "\n";

        // Update poses and save //
        std::string out_p_file = "output_poses.txt";
        std::fstream eo_file;
        eo_file.open(out_p_file.c_str(), std::istream::out);

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


            cPose * pose0 = mTriSet->mGlobalPoses[pose.first];

            Mat3d newR = pose0->R()*dW;
            pose0->Show();
            std::cout << "dT=" << dT << "\n";
            eo_file << pose.first << " " << newR(0,0) << " " << newR(0,1) << " " << newR(0,2)
                                  << " " << newR(1,0) << " " << newR(1,1) << " " << newR(1,2)
                                  << " " << newR(2,0) << " " << newR(2,1) << " " << newR(2,2)
                                  << " " << pose0->C()[0]+dT[0]
                                  << " " << pose0->C()[1]+dT[1]
                                  << " " << pose0->C()[2]+dT[2] << "\n";
        }
        eo_file.close();

    }
    else
        LOG(INFO) << "0 poses to propagate.";
*/
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
     if (DicBoolFind(mTriSet->mAllViewMap,viewName))
     {

         double*  aC1  = mTriSet->mAllViewMap[viewName]->View(0).C();
         double*  aC2  = mTriSet->mAllViewMap[viewName]->View(1).C();

         of << aC1[0] << ' ' << aC1[1] << ' ' << aC1[2]
               << " 0 255 0" << '\n';

         of << aC2[0] << ' ' << aC2[1] << ' ' << aC2[2]
               << " 0 255 0" << '\n';

         if (mTriSet->mAllViewMap[viewName]->NbView()==3)
         {
            double*  aC3  = mTriSet->mAllViewMap[viewName]->View(2).C();
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

cAppCovInMotion::~cAppCovInMotion()
{
    delete mTriSet;
}

cAppCovInMotion::cAppCovInMotion(const InputFiles& inputs,
                    const LocalBundleOptions& lba_opts,
                    const GlobalBundleOptions& gba_opts,
                    const bool do_extra_cov_check) :
    GET_COVARIANCES(do_extra_cov_check)
{
    std::cout << _PROC_COUNT << "\n";

    try 
    {
        // initialize bundle adjustment params 
        m_lba_opts = lba_opts;
        m_gba_opts = gba_opts;


        // initialize the structure containing local motions 
        mTriSet = new cTripletSet(inputs.views_file,
                              inputs.similarities_file,
                              inputs.global_poses_file);

        // read input data 
        if ( mTriSet->ReadViews())
            throw "Couldn't read the views";
 
        if ( mTriSet->ReadSimGlobal())
            throw "Couldn't read the similarities";
 
        if ( mTriSet->ReadGlobalPoses())
            throw "Couldn't read the initial global poses";
 
        if ( ReadFeatures(inputs.tracks_file))
            throw "Couldn't read the features";
 
        if (0)
          mTriSet->PrintAllViews();
        if (0)
          PrintAllPts3d();
        if (0)
          mTriSet->PrintAllPoses();

 
        /* Get covariances per motion */
        OptimizeRelMotions();
 
        /* Run global bundle adjustment using "relative" covariances */
        OptimizeRelMotionsGlobally();
        //OptimizeRelMotionsGlobally();
 
        mTriSet->SaveGlobalPoses(inputs.output_poses_file);
        
        if (0)
          mTriSet->PrintAllPosesDelta();

    }
    catch (std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
    }
    catch(...) 
    {
        std::cerr << "Exception of unknown type!\n";
    }
}

int cov_in_motions_main(InputFiles& inputs,
                        LocalBundleOptions& lba_opts,
                        GlobalBundleOptions& gba_opts,
                        bool do_extra_cov_check)
{

    cAppCovInMotion refinement(inputs,lba_opts,gba_opts,
                               do_extra_cov_check);

    return EXIT_SUCCESS;
}
