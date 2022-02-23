#include "triplet_utils.h"

DEFINE_string(views_file,
              "",
              "The file containing relative motions");
DEFINE_string(similarities_file,
              "",
              "The file containing the 3D similarities between the relative and global poses");
DEFINE_string(global_poses_file,
              "",
              "The file containing the inital global poses");
DEFINE_string(output_poses_file,
              "output_poses.txt",
              "Output file with inlier (inlier_out_poses.txt) or outlier poses (outlier_out_poses.txt) or all of it output_poses.txt");
DEFINE_bool(filter,
              false,
              "Filter relative motions");
DEFINE_string(filtered_view_file,
              "",
              "The file with inlier and outlier relative motions");


bool cTripletSet::ReadViews()
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
                  
    return EXIT_SUCCESS;
}

bool cTripletSet::ReadSimGlobal()
{
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
    
    return EXIT_SUCCESS;
}

bool cTripletSet::ReadRotTrS(FILE* fptr,Mat3d& alpha,Vec3d& beta,double& s)
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

bool cTripletSet::ReadGlobalPoses()
{
    FILE* fptr = fopen(mglob_p_file.c_str(), "r");
    if (fptr == NULL) {
      return false;
    };

    int OK=1;
    while (!std::feof(fptr)  && !ferror(fptr))
    {
        Mat3d aR;
        double * aT = new double [3];

        char aName[20];
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
        if (OK>0)
        {
            for (int aK1=0; aK1<3; aK1++)
            {
                OK = fscanf(fptr, "%lf", &aT[aK1]);
                std::cout << aT[aK1] << " " << aK1 << " ";
            }

            mGlobalPoses[PoseName] = new cPose(aR,aT,PoseName);
            //mGlobalPoses[PoseName]->Show();
            //getchar();
        }
        
    }
    
    return EXIT_SUCCESS;
}



int main(int argc,char** argv)
{

   FLAGS_log_dir = "./";

   google::InitGoogleLogging(argv[0]);
   GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

   LOG(INFO) << "triplet_utils --views_file=" << FLAGS_views_file
            << " --similarities_file=" << FLAGS_similarities_file
            << " --global_poses_file=" << FLAGS_global_poses_file
            << " --output_poses_file=" << FLAGS_output_poses_file
            << " --filter=" << FLAGS_filter
            << " --filtered_file=" << FLAGS_filtered_view_file;
     
    cTripletSet TriSet(FLAGS_views_file,FLAGS_similarities_file,FLAGS_global_poses_file);

    std::cout << "reussite!" << std::endl;
    /*cAppCovInMotion anAp(FLAGS_views_file,"",
                         FLAGS_similarities_file,
                         FLAGS_global_poses_file,
                         FLAGS_output_poses_file,
                         false);
  */ 
    return EXIT_SUCCESS;
}
