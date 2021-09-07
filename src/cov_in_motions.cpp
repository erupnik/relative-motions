#include "all.h"
#include "cov_in_motions.h"

/* TODO
-
  */
double _FEAT_PDS=1;

double RandUnif_0_1()
{
   return ((double) rand() / RAND_MAX );
}

double RandUnif_C()
{
   return (RandUnif_0_1()-0.5) * 2.0;
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
          double aC31[3];

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
              //std::cout << aC21(aK1,0) << " TTT  ";
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

void cAppCovInMotion::BuildProblem(int aNbCam)
{
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
bool cAppCovInMotion::Optimize()
{

    //2-view movements
    for (auto a2v : m2ViewMap)
    {
      // save initial structure
      std::string aIniStruc = a2v.first + "_init.ply";
      WriteToPLYFile(aIniStruc, a2v.first);

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
      std::vector<std::vector<Vec2d>>* aFeat2d = mFeatViewMap[a2v.first];
      int aNbPts = aFeat2d->size();
      std::cout << a2v.first << " " << aNbPts << "\n";

      // get points 3D
      std::vector<double*>* aFeat3d = mFeat3dMap[a2v.first];
      if (int(aFeat3d->size()) != aNbPts )
      {
          std::cout << "Inconsistency in image observations and 3D points. int(aPtVec->size()) != aNbPts";
          return false;
      }

      //create cost function per observation and AddResidualBlock
      //for each camera
      int aNbCam=2;
      std::vector<ResidualBlockId> aPoseResBlockId;
      for (int aCam=0; aCam<aNbCam; aCam++)
      {
          //get EO
          double*  aC = a2v.second->View(aCam).C();
          double   aCPds = 0.1;

          Mat3d   aR0 = a2v.second->View(aCam).R();
          double*  aW = a2v.second->View(aCam).Omega();
          double   aWPds = 0.0017*2;

          //residuals on features
          for (int aK=0; aK<aNbPts; aK++)
          {
              //residuals on features
              CostFunction * aCostF = cResidualError::Create( aR0, aFeat2d->at(aK).at(aCam),_FEAT_PDS);
              LossFunction * aLossF = new HuberLoss(1.0);
              aProblem->AddResidualBlock(aCostF,aLossF,
                                         aW, aC, aFeat3d->at(aK));
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
          /*ceres::SubsetParameterization *constant_base_parameterization = NULL;
          if (aCam==1)
          {
              std::vector<int> constant_translation;
              constant_translation.push_back(0);//x component

              constant_base_parameterization =
                  new ceres::SubsetParameterization(3, constant_translation);

              aProblem->SetParameterization(aC, constant_base_parameterization);
            }*/
          // ===== END FIX GAUGE


          // "rappel" on the poses (except for the 1st one which is constant)
          if (aCam!=0)
          {
              //perspective center
              ResidualBlockId aIdC = aProblem->AddResidualBlock(cPoseConstraint::Create(aC,aCPds),NULL,(aC));
              aPoseResBlockId.push_back(aIdC);
              //small rotation
              ResidualBlockId aIdR = aProblem->AddResidualBlock(cPoseConstraint::Create(aW,aWPds),NULL,(aW));
              aPoseResBlockId.push_back(aIdR);
          }

          //request covariance computation
          if (aCam==1)
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

        //get the jacobians
        double cost = aSummary.final_cost;// 0.0;
        Problem::EvaluateOptions aEvalOpts;
        aEvalOpts.residual_blocks = aPoseResBlockId;
        std::vector<double> evaluated_residuals;
        CRSMatrix *aJ = new CRSMatrix;

        //aProblem->Evaluate(aEvalOpts, &cost, &evaluated_residuals, nullptr, aJ);
        aProblem->EvaluateResidualBlockAssumingParametersUnchanged(aPoseResBlockId,true,&cost,&evaluated_residuals,aJ);
        //void Problem::GetParameterBlocksForResidualBlock(const ResidualBlockId residual_block, vector<double*> *parameter_blocks)const

        std::cout << aJ->num_rows << " " << aJ->num_cols << "\n";
        getchar();

        //Eigen::MappedSparseMatrix which will let you wrap a CRS matrix directly
        Eigen::SparseMatrix<double> A(aJ->num_rows, aJ->num_cols);
      	int row = 0;
      	int el = 0;
      	for (int i = 1; i < aJ->rows.size(); i++)
      	{
      		int NumInRow = aJ->rows[i] - aJ->rows[i - 1];
      		for (int k = 0; k < NumInRow; k++)
      		{
      			int col = aJ->cols[el];
      			A.coeffRef(row, col) = aJ->values[el++];
            std::cout << A.coeffRef(row, col) << " ";
      		}
          std::cout << "\n";
      		row++;
      	}
        for( auto i=0; i<evaluated_residuals.size() ; i++ )
            std::cout << i << ": " << evaluated_residuals[i] << "\n";

        //get covariances
        CHECK(covariance.Compute(aCovBlocks, aProblem));
        for (int aCam=1; aCam<aNbCam; aCam++)
        {

            covariance.GetCovarianceBlock(a2v.second->View(aCam).C(), a2v.second->View(aCam).C(), a2v.second->View(aCam).Cov_C());
            covariance.GetCovarianceBlock(a2v.second->View(aCam).Omega(), a2v.second->View(aCam).Omega(), a2v.second->View(aCam).Cov_Omega());

            a2v.second->View(aCam).Show();

        }



        //save final structure
        std::string aFinStruc = a2v.first + "_final.ply";
        WriteToPLYFile(aFinStruc, a2v.first);

    }

    for (auto a3v : m3ViewMap)
    {
        //std::cout << a3v.first << "\n";
    }

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

cAppCovInMotion::cAppCovInMotion(const std::string& avfile,const std::string& affile) :
  mviews_file(avfile),
  mfeats_file(affile)
{
    ReadFeatures();
    ReadViews();

    if (0)
      PrintAllViews();
    if (0)
      PrintAllPts3d();

    Optimize();
}

int cov_in_motions_main(int argc, char** argv)
{

  std::string aArg1 = std::string(argv[1]);
  std::string aArg2 = std::string(argv[2]);

  cAppCovInMotion anAp(aArg1,aArg2);


  return 1;
}
