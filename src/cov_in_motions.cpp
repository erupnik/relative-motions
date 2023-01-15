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
  
    std::cout << "\nC \t" << mC_immutable[0] << " " << mC_immutable[1] << " " << mC_immutable[2] << "\n";

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
void cNviewPoseX::decompose_H()
{
    /* Get covariance matrix A and vector B */
    MatXd H = mHg->H_();
    VecXd g = mHg->g_();

    /* Get X0 = {c1x,c1y,c1z,w1x,w1y,w1z,
     *           c2x,c2y,c2z,w2x,w2y,w2z,
     *(optional) c3x,c3y,c3z,w3x,w3y,w3z} */ 
    int sz = g.size();

    VecXd X0 = VecXd::Zero(sz); 

    for (int v=0; v<mNbV; v++)
    {
        X0[v*6] = this->View(v).C()[0]; 
        X0[v*6+1] = this->View(v).C()[1];
        X0[v*6+2] = this->View(v).C()[2];
        X0[v*6+3] = this->View(v).Omega()[0]; 
        X0[v*6+4] = this->View(v).Omega()[1];
        X0[v*6+5] = this->View(v).Omega()[2];
    }
    //std::cout << "X0\n" << X0 << "\n";
    /* Retrieve residual vector */
    VecXd res_vec = g + H * X0;

    /* Solve again for X0 */
    VecXd Sol = H.ldlt().solve(res_vec); 
    //std::cout << "Sol\n" << Sol << "\n";

    /* Compute eigenvectors and eigenvalues */
    int NbEl = 6*mNbV;
    mWi = VecXd::Zero(NbEl);
    mLi = MatXd::Zero(NbEl,NbEl);
    mCstei = VecXd::Zero(NbEl);

    //for (int v=0; v<mNbV; v++)
    //{

        /*Eigen::EigenSolver<MatXd> es(H.block(6*v,6*v,6,6), true);

        MatXd EV = es.eigenvectors().real().transpose();
 
        mLi.block(6*v,6*v,6,6)  = EV;
        mCstei.block(6*v,0,6,1)       = EV * Sol.block(6*v,0,6,1);
        mWi.block(6*v,0,6,1)          = es.eigenvalues().real();

*/
        /*for (int unk=0 ; unk<2; unk++)
        {
            Eigen::EigenSolver<MatXd> es(H.block(6*v+3*unk,6*v+3*unk,3,3), true);

            MatXd EV = es.eigenvectors().real().transpose();
 
            mLi.block(6*v+3*unk,6*v+3*unk,3,3)  = EV;
            mCstei.block(6*v+3*unk,0,3,1)       = EV * Sol.block(6*v+3*unk,0,3,1);
            mWi.block(6*v+3*unk,0,3,1)          = es.eigenvalues().real();

        }*/
        

    //}

/*    Eigen::EigenSolver<MatXd> es(H, true);

    //MatXd EV = es.eigenvectors().transpose();
    MatXd EV = es.eigenvectors().real().transpose();
 
    mLi     = EV;
    mCstei  = EV * Sol;
    mWi     = es.eigenvalues().real();
   */ 


   Eigen::SelfAdjointEigenSolver<MatXd> es;
   es.compute(H);

   mLi = es.eigenvectors().transpose() ;
   mWi = es.eigenvalues();
   mCstei = mLi * Sol;

   //zero-out negative eigenvalues 
   //(sometimes due to numerical instability values can be neg)
   if (mWi[0] < 0)
   {
       mWi[0]=0.0;
       VLOG(1) << "=== Zero-out negative eigenvalue " << mWi[0] << "\n" ;
   }
/*
SelfAdjointEigenSolver<Matrix4f> es;
Matrix4f X = Matrix4f::Random(4,4);
Matrix4f A = X + X.transpose();
es.compute(A);
cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << endl;
es.compute(A + Matrix4f::Identity(4,4)); // re-use es to compute eigenvalues of A+I
cout << "The eigenvalues of A+I are: " << es.eigenvalues().transpose() << endl;
*/

/*
   std::cout << "mH\n" <<  H << "\n";
    std::cout << "mLi\n" <<  mLi << "\n";
    std::cout << "mW\n" << mWi  << "\n";
    std::cout << "Cstei\n" <<  mCstei << "\n";
getchar(); */
   
}

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
   


    // save initial structure
    //std::string aIniStruc = views_name + "_init.ply";
    //WriteToPLYFile(aIniStruc, views_name);
 
    //create ceres problem
    ceres::Problem * aProblem = new ceres::Problem;;
    ceres::Solver::Summary        aSummary;
 
    //covariance
    Covariance::Options aCovOpts;
    //aCovOpts.null_space_rank = -1;
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
 
    //std::cout << views_name << " " << aNbPts << "\n";
 
    // get points 3D
    std::vector<double*>* aFeat3d = mFeat3dMap[views_name];
    if (int(aFeat3d->size()) != aNbPts )
    {
        std::cout << "Inconsistency in image observations and 3D points. int(aPtVec->size()) != aNbPts";
        return false;
    }
 
    //create cost function per observation and AddResidualBlock
    std::vector<double*> param_bloc_CRP;
    //std::vector<double*> param_bloc_PT3D;
    //std::vector<ceres::ResidualBlockId> res_bloc_id;
 
    int NbCam = (&(views->View(2))==NULL) ? 2 : 3;
    for (int aCam=0; aCam<NbCam; aCam++)
    {
        double*  aC = views->View(aCam).C();
        double*  aW = views->View(aCam).Omega();
        
        //keep reference to parameters to later retrieve their H,g
        param_bloc_CRP.push_back(aC);
        param_bloc_CRP.push_back(aW);
    }

    int num_im_obs=0;
    for (int aCam=0; aCam<NbCam; aCam++)
    {
        //get EO
        double*  aC = views->View(aCam).C();
 
        Mat3d   aR0 = views->View(aCam).R();
        double*  aW = views->View(aCam).Omega();
 
        double PdsAtten = CostAttenuate(aNbPts,m_lba_opts._NB_LIAIS);

        //residuals on features
        //std::cout << "****CAM" << aCam << " features nb=" << aNbPts << " " << PdsAtten << "\n";
        for (int aK=0; aK<aNbPts; aK++)
        {
            //residuals on features
            if (1)
            {
                //compute the residual first 
                cResidualError eval    (aR0,
                                        aFeat2d->at(aK).at(aCam),
                                        1.0,
                                        m_lba_opts._FOCAL);
                double res_non_pond =0;
                eval(aW, aC, aFeat3d->at(aK),&res_non_pond);
                //if (res_non_pond*m_lba_opts._FOCAL <3)

                if (res_non_pond<m_lba_opts._MAX_ERR)//viabon and PVA with "3"
                {
                    //std::cout << "==>" << res_non_pond*5560 << " " ;
                 
                    CostFunction * aCostF = cResidualError::Create( aR0,
                                                                    aFeat2d->at(aK).at(aCam),
                                                                    PdsAtten,
                                                                    //m_lba_opts._FEAT_PDS,
                                                                    //m_lba_opts._FOCAL); rounding errors lead to
                                                                    //negative eigenvalues 
                                                                    1.0);
                    LossFunction * aLossF = new HuberLoss(m_lba_opts._HUBER);
                    aProblem->AddResidualBlock(aCostF,aLossF,
                                               aW, aC, aFeat3d->at(aK));

                    num_im_obs++;
                    
                }
               // else
                //    std::cout << "PONDDD ERR=" << res_non_pond << "\n";
                //res_bloc_id.push_back(res_PC_id);
                //if (aCam==0) //add the 3D points only once
                //  param_bloc_CRP.push_back(aFeat3d->at(aK));
            }
            
            
        }  


        // "rappel" on the poses
        if (_GAUGE_RAPPEL_POSES)
        {
          
            //perspective center
            double*  aC_immu = views->View(aCam).C_immutable();
            CostFunction * aCostC = cPoseConstraint::Create(aC_immu,m_lba_opts._C_PDS);
            aProblem->AddResidualBlock(aCostC,NULL,(aC));
            //param_bloc_CR.push_back(aC);
            
            //small rotation
            double*  aW_immu = views->View(aCam).Omega_immutable();
            CostFunction * aCostR = cPoseConstraint::Create(aW_immu,m_lba_opts._ROT_PDS);
            aProblem->AddResidualBlock(aCostR,NULL,(aW));
            //param_bloc_CR.push_back(aW);
           
        }
 
    }
    //std::cout << "PONDERATION inlier/all=" << double(num_im_obs)/(aNbPts*3)*100.0 << "\n";
    
    // check that there is a min number of image observations 
    if (num_im_obs < m_lba_opts._MIN_NUM_OBS)
    {
        VLOG(1) << "Not enough image observations in " << views_name << "\n" ;
        return false;
    }

    // define between which parameters the covs should be computed 
    for (int aCam0=0; aCam0<NbCam; aCam0++)
    {
        //get EO
        double*  aC0 = views->View(aCam0).C();
        double*  aW0 = views->View(aCam0).Omega();
        
        for (int aCam1=0; aCam1<NbCam; aCam1++)
        {
            //std::cout << aCam0 << " " <<  aCam1 << "\n";
            //get EO
            double*  aC1 = views->View(aCam1).C();
            double*  aW1 = views->View(aCam1).Omega();

            aCovBlocks.push_back(std::make_pair(aC0,aC1));
            aCovBlocks.push_back(std::make_pair(aC0,aW1));
            aCovBlocks.push_back(std::make_pair(aW0,aW1));
            aCovBlocks.push_back(std::make_pair(aW0,aC1));
        }
    }


    //solve the least squares problem
    ceres::Solver::Options aOpts;
    SetMinimizerLocal(aOpts);
    
    // Set ordering (necessary for Schur)
    ceres::ParameterBlockOrdering* ordering =
    new ceres::ParameterBlockOrdering;
    // The points come before the cameras.
    for (int aK=0; aK<aNbPts; aK++)
    {
        ordering->AddElementToGroup(aFeat3d->at(aK), 0);
    }
    for (int aCam=0; aCam<NbCam; aCam++)
    {
        ordering->AddElementToGroup(views->View(aCam).C(), 1);
        ordering->AddElementToGroup(views->View(aCam).Omega(), 1);
        
    }
    aOpts.linear_solver_ordering.reset(ordering);


    ceres::Solve(aOpts,aProblem,&aSummary);
    //std::cout << aSummary.FullReport() << "\n";
    
    // show first camera
    //views->View(0).Show();

    //get covariances
    MatXd CovMat = MatXd::Zero(NbCam*6,NbCam*6);
    MatXd HMat = MatXd::Zero(NbCam*6,NbCam*6);


    int NbVar = std::pow(NbCam,2); // =4 for NbCam=2 and 9 for NbCam=3
    double cov_c[NbVar][18]; //18: 3x3 for Cov(Ci,Cj) and 3x3 for Cov(Ci,Wj)
    double cov_w[NbVar][18];
    

    //if rank deficiency in jacobian, ignore this view 
    bool RANK_OK = covariance.Compute(aCovBlocks, aProblem);


    if (!RANK_OK)
    {
        VLOG(1) << "====> RANK DEFICIT IN JACOBIAN IN " << views_name << "\n";
        return false;
    }
    else //- retrieve the covariances &  fill the CovMat matrix
    {
        for (int aCam0=0; aCam0<NbCam; aCam0++)
        {
            for (int aCam1=0; aCam1<NbCam; aCam1++)
            {
             
                  covariance.GetCovarianceBlock(views->View(aCam0).C(), views->View(aCam1).C(),     &cov_c[aCam0*NbCam+aCam1][0]);
                  covariance.GetCovarianceBlock(views->View(aCam0).C(), views->View(aCam1).Omega(), &cov_c[aCam0*NbCam+aCam1][9]);
                  
                  covariance.GetCovarianceBlock(views->View(aCam0).Omega(),views->View(aCam1).Omega(), &cov_w[aCam0*NbCam+aCam1][0]);
                  covariance.GetCovarianceBlock(views->View(aCam0).Omega(),views->View(aCam1).C(),     &cov_w[aCam0*NbCam+aCam1][9]);
            }
        }
        
        for (int aCam0=0; aCam0<NbCam; aCam0++)
        {
            for (int aCam1=0; aCam1<NbCam; aCam1++)
            {
     
                //c-c 
                int cmpt=0;
                for (int i=0; i<3; i++)
                   for (int j=0; j<3; j++)
                   { 
                       CovMat(6*aCam0+i,6*aCam1+j) = cov_c[aCam0*NbCam+aCam1][cmpt];
                       cmpt++;
                    }
                //cw 
                for (int i=0; i<3; i++)
                   for (int j=0; j<3; j++)
                   {
                       CovMat(6*aCam0+i,6*aCam1+3+j) = cov_c[aCam0*NbCam+aCam1][cmpt];
                       cmpt++;
                    }
                //ww
                cmpt=0;
                for (int i=0; i<3; i++)
                   for (int j=0; j<3; j++)
                   {
                       CovMat(3+6*aCam0+i,3+6*aCam1+j) = cov_w[aCam0*NbCam+aCam1][cmpt];
                       cmpt++;
                    }
                //wc 
                for (int i=0; i<3; i++)
                   for (int j=0; j<3; j++)
                   {
                       CovMat(3+6*aCam0+i,3+6*aCam1-3+j) = cov_w[aCam0*NbCam+aCam1][cmpt];
                       cmpt++;
                    }
            }
        }
    }

    // - compute the hessian HMat = CovMat^-1
    Eigen::ConjugateGradient<MatXd > solver;
    solver.compute(CovMat);
    Eigen::SparseMatrix<double> I(NbCam*6,NbCam*6);
    I.setIdentity();
    Eigen::SparseMatrix<double> CovMat_inv = solver.solve(I);
    MapJacToEigMat(CovMat_inv,HMat,0);

    //make symmetric again, otherwise negative eigenvalues 
    //MatXd HMatSym = 0.5* (HMat + HMat.transpose());


    //get the residuals (vector B) 
    double cost=aSummary.final_cost;//
    Problem::EvaluateOptions aEvalOpts;
    aEvalOpts.apply_loss_function = true;

    aEvalOpts.parameter_blocks = param_bloc_CRP;
    std::vector<double> res;
    //std::vector<double> JtWres;
    //CRSMatrix *aJ = new CRSMatrix;

    bool eval_success = aProblem->Evaluate(aEvalOpts, &cost, &res, NULL, NULL);
    double res_total=0.0;
    for (auto r : res)
    {
        res_total+=std::abs(r);
        //std::cout << "res=" << r << "\n";
    }


    //add 1 to avois div by 0
    res_total += 1.0;
    //std::cout << "\ntotal res=" << res_total << "\n";
    //res_total *= m_lba_opts._FOCAL; 

    views->TotalRes() =  1.0;//aNbPts/res_total; //res_total/aNbPts res 0.69


    //Eigen::VectorXd my_vect = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(JtWres.data(), JtWres.size());//fixed-size to dynamic
    //Eigen::VectorXd my_vect_poses = my_vect.block(0,0,6*views->NbView(),1);
    VecXd my_vect_poses = VecXd::Zero(6*views->NbView());

    //store  hessian and vector B 
    views->Hg_() = *(new cHessianGradientX(HMat,my_vect_poses));

    //std::cout << "HMat\n" << HMat << "\n";
    //std::cout << "HMatSym\n" << HMatSym << "\n";  
    //getchar();

    //decompose the quadratic form into a sum of linear terms 
    views->decompose_H();

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
                << " _GAUGE_FIRST_CAM_FIX=" <<  _GAUGE_FIRST_CAM_FIX << " ";
             //   << "\nHHHHHHHHHHHHHHHHHHHHHHHh\n" << JtWJ_.sqrt();
       
    }
    //print residuals vector
    if (0)
    {
        VLOG(2) << "residuals vector \n";
        for (int i=0; i<res.size(); i++)
        {
            VLOG(2) << res[i] << ", ";
        }
        VLOG(2) << "\n";
    }
    

    if (0)
    {
      //save covariance diagonal
      std::string covR_file_name = "covRsqrt" + m_lba_opts._KEY + ".txt";
      std::fstream covR_file;
      covR_file.open(covR_file_name.c_str(), std::ios_base::app);
      std::string covC_file_name = "covCsqrt" + m_lba_opts._KEY + ".txt";
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

   
      VLOG(2) << "\n";

      gradR_file << "\n";
      gradR_file.close();
      gradC_file << "\n";
      gradC_file.close();

    }
  
    //delete ordering;
    delete aProblem;
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
  //aSolOpt.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;
  //uses the SHURR trick; solves S as a sparse matrix with Cholesky factorization;
  aSolOpt.linear_solver_type = ceres::DENSE_SCHUR;// ceres::SPARSE_SCHUR;
  //applies Preconditioned Conjugate Gradients to S; implements inexact step algorithm;
  //choose a precondition, e.g. CLUSTER_JACOBI,CLUSTER_TRIDIAGONAL that exploits camera-point visibility structure
  //aSolOpt.linear_solver_type = ceres::ITERATIVE_SCHUR;
  //aSolOpt.use_explicit_schur_complement = true;

  //relaxes the requirement to decrease the obj function at each iter step;
  //may turn very efficient in the long term;
  aSolOpt.use_nonmonotonic_steps = false;

  aSolOpt.max_num_iterations = 10;
  aSolOpt.minimizer_progress_to_stdout =  false;//true;
  aSolOpt.num_threads = m_lba_opts._PROC_COUNT;

  aSolOpt.use_inner_iterations = m_lba_opts._INNER_ITER;

}

//final, global least-squares simultaneously on all poses  
void cAppCovInMotion::SetMinimizerGlobal(ceres::Solver::Options& aSolOpt)
{
  //S - the reduced camera matrix / the Shur complement;
  //uses the SHURR trick; solves S as a dense matrix with Cholesky factorization; i
  //for problems up to several hundreds of cameras
  //aSolOpt.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;//ceres::DENSE_SCHUR;
  //aSolOpt.linear_solver_type = ceres::CGNR;
  //uses the SHURR trick; solves S as a sparse matrix with Cholesky factorization;
  aSolOpt.linear_solver_type = ceres::SPARSE_SCHUR;
  //aSolOpt.linear_solver_type = ceres::DENSE_SCHUR;
  //applies Preconditioned Conjugate Gradients to S; implements inexact step algorithm;
  //choose a precondition, e.g. CLUSTER_JACOBI,CLUSTER_TRIDIAGONAL that exploits camera-point visibility structure
  //aSolOpt.linear_solver_type = ceres::ITERATIVE_SCHUR;

  //relaxes the requirement to decrease the obj function at each iter step;
  //may turn very efficient in the long term;
  aSolOpt.use_nonmonotonic_steps = false;

  /*aSolOpt.function_tolerance = 0.0;
  aSolOpt.gradient_tolerance = 0.0;
  aSolOpt.parameter_tolerance = 0.0;*/

  aSolOpt.max_num_iterations = 200;
  aSolOpt.minimizer_progress_to_stdout = true;
  aSolOpt.num_threads = m_gba_opts._PROC_COUNT;

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

    
    //std::vector<std::string> views_to_del_by_key;

    //#pragma omp parallel num_threads(1) 
    #pragma omp parallel num_threads(m_lba_opts._PROC_COUNT) 
    {
        #pragma omp for
        for(int vi=0; vi<mTriSet->mAllViewMap.size(); vi++) 
        {
            auto vi_it = mTriSet->mAllViewMap.begin();
            std::advance(vi_it, vi);

            if (m_lba_opts._RUN_PROP)
            {
                //std::cout << (*vi_it).first << "\n";
                if (! BuildProblem_((*vi_it).second,(*vi_it).first))
                {
                    (*vi_it).second->SetOutlier();
                    std::cout << "ELIMINATED!" << "\n";
                    //views_to_del_by_key.push_back((*vi_it).first);
                }

            }
            else 
                InitCovariances((*vi_it).second,(*vi_it).first);
        }

    }

    //count eliminated triplets 
    int NbTripElim=0;
    for(auto vi : mTriSet->mAllViewMap) 
    {
        if (! (vi.second)->IsInlier())
            NbTripElim++;
    }
    std::cout << "Nb eliminated triplets: " << NbTripElim << "\n"; 
    //delete views with no features or views with rank deficient jacobians  
    //std::cout << " # " << views_to_del_by_key.size() << " view removed\n";
    /*{

    //for (auto vi : views_to_del_by_key)
    for (int vi=0; vi<int(views_to_del_by_key.size()); vi++)
    {
        std::cout << "-vi " << (views_to_del_by_key[vi]) ;
        std::map<std::string, cNviewPoseX*>::iterator itr = mTriSet->mAllViewMap.find(views_to_del_by_key[vi]);
        std::cout << "- " << itr->first << "\n" ;
        if (itr != mTriSet->mAllViewMap.end())
        {
        //   delete itr->second;
           mTriSet->mAllViewMap.erase(itr);
           //mTriSet->mAllViewMap.erase(vi);
           std::cout << "+\n";
        }
        else
           std::cout << "not found \n";

    }
    

    }*/

/*
 * std::map<std::string, Texture*>::iterator itr = textureMap.find("some/path.png");
   if (itr != textureMap.end())
   {
    // found it - delete it
    delete itr->second;
    textureMap.erase(itr);
   }
 *
 *
 * */

    std::cout << "***************\n";

    // save local covariances to file if requested 
    if (m_lba_opts._WRITE_COV)
        WriteLocalCovs();

    //print out the total res
    //for (auto v : mTriSet->mAllViewMap)
    //    v.second->PrintTotalRes();
    //print out number of features per motion
    /*for (auto v :  mFeatViewMap)
        std::cout << v.second->size() << " ";
    std::cout << "\n";*/ 

   // getchar();
    return EXIT_SUCCESS;

}

bool cAppCovInMotion::OptimizeRelMotionsGloballyAllViewsWithDecomp()
{

    std::cout << "Global view optimization with quadratic form decomposition\n";
    
    ceres::Problem * aProblem = new ceres::Problem;;
    ceres::Solver::Summary        aSummary;


    //add equation for each view and camera 
    for(auto view_i : mTriSet->mAllViewMap)
    {
        if (view_i.second->IsInlier())
        {
            //std::cout << view_i.first << "\n";
         
            //similitude
            Mat3d alpha0 = view_i.second->alpha0();
            //double* Walpha = view_i.second->Walpha();
            //double* beta  = view_i.second->beta();
            //double* lambda = view_i.second->lambda();
            double* alpha_beta_l = view_i.second->alpha_beta_l();
         
            int NbCam = (&(view_i.second->View(2))==NULL) ? 2 : 3;
         
            //collect unknowns 
            std::vector<double*>  cVec; //relative center (only for basic adj) 
            std::vector<Mat3d>    rVec; //relative rotations 
            std::vector<Mat3d>    RVec; //global rotations 
            
         
            for (int aCam=0; aCam<NbCam; aCam++)
            {
                //local pose 
                rVec.push_back(view_i.second->View(aCam).R());
                cVec.push_back(view_i.second->View(aCam).C());
         
         
                //initial global pose
                cPose*& pose = mTriSet->Pose(view_i.second->View(aCam).Name());
         
                //parameters / constants  
                RVec.push_back(pose->R());
              
            }
         
            //covariance 
            VecXd Wi = view_i.second->Wi();
            MatXd Li = view_i.second->Li();
            VecXd Cstei = view_i.second->Cstei();
         
            double total_res = view_i.second->TotalRes();
         
         
            cPose*& pose0 = mTriSet->Pose(view_i.second->View(0).Name());
            cPose*& pose1 = mTriSet->Pose(view_i.second->View(1).Name());
            
            double * C0 = pose0->C();
            double * W0 = pose0->Omega();
            
            double * C1 = pose1->C();
            double * W1 = pose1->Omega();
         
            if (NbCam==2)
            {
                if (m_gba_opts._PROPAGATE)
                {
                    /* Unknown poses, constant similarity */ 
                    /*CostFunction * aCost = cResidualOn2ViewsPoseDecomp::Create(alpha0,beta,lambda,
                                                          rVec,RVec,
                                                          Wi,Li,Cstei) ;
                    aProblem->AddResidualBlock(aCost,NULL,
                                         C0,W0,C1,W1);*/
                 
                    CostFunction * aCost = cResidualOn2ViewsPoseDecompLAB::Create(alpha0,
                                                          rVec,RVec,
                                                          Wi,Li,Cstei,total_res) ;
              
                    LossFunction * aLoss = NULL; //new HuberLoss(m_gba_opts._HUBER_S);
         
                    aProblem->AddResidualBlock(aCost,aLoss,
                                         C0,W0,C1,W1,alpha_beta_l);
                }
                else
                {
                    CostFunction* Cost = cResidualOn3ViewsPoseBasicLAB::Create(alpha0,cVec,rVec,RVec,
                                                                                m_gba_opts._C_PDS,
                                                                                m_gba_opts._ROT_PDS);
                    LossFunction * Loss = NULL;
                    aProblem->AddResidualBlock(Cost,Loss,C0,W0,C1,W1,alpha_beta_l);
                                                                               
                }
            }
         
         
            if (NbCam==3)
            {
                
                cPose*& pose2 = mTriSet->Pose(view_i.second->View(2).Name());
             
                double * C2 = pose2->C();
                double * W2 = pose2->Omega();
            
                if (m_gba_opts._PROPAGATE)
                {
                    /* Unknown poses, constant similarity */ 
                    /*Vec3d beta {alpha_beta_l[3],alpha_beta_l[4],alpha_beta_l[5]};
                    double L = alpha_beta_l[6];
                    CostFunction * aCost = cResidualOn3ViewsPoseDecomp::Create(alpha0,beta,L,
                                                          rVec,RVec,
                                                          Wi,Li,Cstei) ;
                    aProblem->AddResidualBlock(aCost,NULL,
                                         C0,W0,C1,W1,C2,W2); */
              
                    CostFunction * aCost = cResidualOn3ViewsPoseDecompLAB::Create(alpha0,
                                                          rVec,RVec,
                                                          Wi,Li,Cstei,
                                                          total_res) ; 
              
                    LossFunction * aLoss = NULL; //new HuberLoss(m_gba_opts._HUBER_S);
                    aProblem->AddResidualBlock(aCost,aLoss,
                                         C0,W0,C1,W1,C2,W2,alpha_beta_l);
                }
                else
                {
                    CostFunction* Cost = cResidualOn3ViewsPoseBasicLAB::Create(alpha0,cVec,rVec,RVec,
                                                                              m_gba_opts._C_PDS,
                                                                              m_gba_opts._ROT_PDS);
                    LossFunction * Loss = NULL;
                    aProblem->AddResidualBlock(Cost,Loss,C0,W0,C1,W1,C2,W2,alpha_beta_l);
         
                }
             
               
            }
           
            for (int aCam=0; aCam<NbCam; aCam++)
                mTriSet->Pose(view_i.second->View(aCam).Name())->SetRefined();

        }
    }

    //add constraint on initial global poses 
    if (m_gba_opts._CONSTRAIN_GPOSE)
    {
        for (auto pose_i : mTriSet->mGlobalPoses)
        {
            if (pose_i.second->IsRefined())
            {
                //parameters 
                double* C = pose_i.second->C();
                double* C_immutable = pose_i.second->C_immutable();
                double* W = pose_i.second->Omega();
                double* W_immutable = pose_i.second->Omega_immutable();
            
                //LossFunction * aLossW = new HuberLoss(m_gba_opts._HUBER_P);
                //LossFunction * aLossC = new HuberLoss(m_gba_opts._HUBER_P);
                
                //perspective center
                CostFunction * CostC = cPoseConstraint::Create(C_immutable,m_gba_opts._C_PDS);
                ceres::ResidualBlockId res_C_id = aProblem->AddResidualBlock(CostC,NULL,(C));
 
                //small rotation
                CostFunction * CostR = cPoseConstraint::Create(W_immutable,m_gba_opts._ROT_PDS);
                ceres::ResidualBlockId res_W_id = aProblem->AddResidualBlock(CostR,NULL,(W));
            }
        }
    }
    
    //solve the least squares problem
    ceres::Solver::Options Opts;
    SetMinimizerGlobal(Opts);

    // Set ordering (necessary for Schur)
    ceres::ParameterBlockOrdering* ordering =
    new ceres::ParameterBlockOrdering;
    // The 3d similarity come before the cameras.
    for(auto view_i : mTriSet->mAllViewMap)
    {
        if (view_i.second->IsInlier())
        {
            //similitude
            double* alpha_beta_l = view_i.second->alpha_beta_l();
         
            ordering->AddElementToGroup(alpha_beta_l, 0);
        }
    }
    for (auto pose_i : mTriSet->mGlobalPoses)
    {
        if (pose_i.second->IsRefined())
        {
            //parameters 
            double* C = pose_i.second->C();
            double* W = pose_i.second->Omega();

            ordering->AddElementToGroup(C, 1);
            ordering->AddElementToGroup(W, 1);
        }
    }
    Opts.linear_solver_ordering.reset(ordering);

/*
 *  // Set ordering (necessary for Schur)
    ceres::ParameterBlockOrdering* ordering =
    new ceres::ParameterBlockOrdering;
    // The points come before the cameras.
    for (int aK=0; aK<aNbPts; aK++)
    {
        ordering->AddElementToGroup(aFeat3d->at(aK), 0);
    }
    for (int aCam=0; aCam<NbCam; aCam++)
    {
        ordering->AddElementToGroup(views->View(aCam).C(), 1);
        ordering->AddElementToGroup(views->View(aCam).Omega(), 1);
        
    }
    aOpts.linear_solver_ordering.reset(ordering);
 * */


    ceres::Solve(Opts,aProblem,&aSummary);
    std::cout << aSummary.FullReport() << "\n";

    if (VLOG_IS_ON(1))
    {
        //Print residuals 
        //FILE* file_residual = fopen("residuals.csv", "w+");

        double cost=aSummary.final_cost;
        Problem::EvaluateOptions eval_opts;
        eval_opts.apply_loss_function = false;
 
        //eval_opts.residual_blocks = res_bloc_CR;
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
    //mTriSet->UpdateAllAffine();

    delete aProblem;

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
    std::cout << m_lba_opts._PROC_COUNT << " "<< m_gba_opts._PROC_COUNT << "\n";

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

        if (0)
            mTriSet->WriteGlobalPFromRelPAndSim("pose_predictions.txt");
 
        /* Get covariances per motion */
        if (m_lba_opts._RUN_PROP)
            OptimizeRelMotions();
 
        /* Run global bundle adjustment using "relative" covariances */
        OptimizeRelMotionsGloballyAllViewsWithDecomp();
        //OptimizeRelMotionsGloballyWithDecomp();
        //OptimizeRelMotionsGlobally();
//getchar(); 
        mTriSet->SaveGlobalPoses(inputs.output_poses_file);
        
        if (0)
            mTriSet->PrintAllPosesDelta();
        if (1)
            mTriSet->PrintSumDelta();

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
