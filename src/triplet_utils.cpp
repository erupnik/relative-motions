#include "triplet_utils.h"

extern int main(int argc,char** argv);
extern double euclid(Vec3d );
static Mat3d NearestRotation(Mat3d&);

//std::string _DELIMITER_ = "-zyx-";

std::vector<std::string> cTripletSet::DecompViewNames(std::string& name)
{
    std::vector<std::string> decomposed_names;

    //first name 
    std::size_t found_delim = name.find(_DELIMITER_);

    if (found_delim!=std::string::npos)
        decomposed_names.push_back(name.substr(0,found_delim));

    //second name
    found_delim = name.find(_DELIMITER_);
    std::size_t found_delim_last = name.rfind(_DELIMITER_);

    //image pair
    //
    if (found_delim == found_delim_last) 
    {
        decomposed_names.push_back(name.substr(found_delim+_DELIMITER_.size()));
    }
    else 
    {
        int length = found_delim_last - found_delim - _DELIMITER_.size();
        decomposed_names.push_back(name.substr(found_delim+_DELIMITER_.size(), length));
        decomposed_names.push_back(name.substr(found_delim_last+_DELIMITER_.size()));

    }

    return decomposed_names;
}

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
            char aName[60];
            fscanf(fptr, "%s", aName);
            aPoseNameV.push_back(aName);
 
            if (aK!=0)
              aViewName+=_DELIMITER_;
 
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
 
                mAllViewMap[aViewName] = new cNviewPoseX
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
 
              mAllViewMap[aViewName] = new cNviewPoseX
                       (aPose1_,aPose21_,aPose31_,new cHessianGradientX(Mat12d::Zero(),Vec12d::Zero()));
 
            }
            else
            {
                std::cout << "Too many views. We currently support up to 3 views." << "\n";
                return 0;
            }
          }
        }
    
    fclose(fptr);

    std::cout << " #" << mAllViewMap.size() << " views\n";

    return EXIT_SUCCESS;
}

bool cTripletSet::ReadSimGlobal()
{

    FILE* fptr = fopen(msimil_file.c_str(), "r");
    if (fptr == NULL) {
      return false;
    };

    int Num=0;
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
            char aName[60];
            fscanf(fptr, "%s", aName);
            aPoseNameV.push_back(aName);

            if (aK!=0)
              aViewName+=_DELIMITER_;

            aViewName += aName;
            //std::cout <<  " " << aName << " aK=" << aK;
        }
        //std::cout <<  aViewName << "\n";

        if (DicBoolFind(mAllViewMap,aViewName))
        {
            //double& L = mAllViewMap[aViewName]->lambda();
            double* L = mAllViewMap[aViewName]->lambda();
            Mat3d&  alpha = mAllViewMap[aViewName]->alpha0();
            double*  beta  = mAllViewMap[aViewName]->beta();

            if (! ReadRotTrS(fptr,alpha,beta,L))
            {
                std::cout << "ERROR reading global similitude in cAppCovInMotion::ReadSimGlobal for "
                          << aViewName << "\n" ;
                return false;
            }

            /*mAllViewMap[aViewName]->PrintAlpha();
            mAllViewMap[aViewName]->PrintBeta();
            mAllViewMap[aViewName]->PrintLambda();
            getchar();*/

            mAllViewMap[aViewName]->SetInit();
            Num++;
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
  
    //Remove rel motions without similarities
    auto v = mAllViewMap.begin();
    while (v != mAllViewMap.end()) 
    {
        if(! (*v).second->Init())
        {
            VLOG(2) << "No similarity for this motion: " << (*v).first << "\n";    
            v = mAllViewMap.erase(v);

        }
        else
            ++v;
    }

    fclose(fptr);

    std::cout << " #" << Num << " similarity transformations\n";

    return EXIT_SUCCESS;
}
bool cTripletSet::ReadRotTrS(FILE* fptr,Mat3d& alpha,double*& beta,double*& s)
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
        OK = fscanf(fptr, "%lf", &beta[aK]);
        //std::cout << beta(aK) << " ";
    }

    //scale
    OK = fscanf(fptr, "%lf", &s[0]);
    //std::cout << s << " " << OK << " s";

    return OK;
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

        char aName[60];
        OK = fscanf(fptr, "%s", aName);
        std::string PoseName = std::string(aName);
        //std::cout << PoseName << "\n";


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
                //std::cout << aT[aK1] << " " << aK1 << " ";
            }

            mGlobalPoses[PoseName] = new cPose(aR,aT,PoseName);
            //mGlobalPoses[PoseName]->Show();
            //getchar();
        }

        if (!OK)
            return EXIT_SUCCESS;
    }
   
    fclose(fptr);

    std::cout << " #" << mGlobalPoses.size() << " inital global poses\n";
    
    return EXIT_SUCCESS;
}

void cTripletSet::PrintAllViews()
{
      std::cout << "Pairs and triplets:\n";
      for (auto Pose : mAllViewMap)
      {
          Pose.second->View(0).Show();
          Pose.second->View(1).Show();
          if (Pose.second->NbView()==3)
            Pose.second->View(2).Show();

      }
}
void cTripletSet::PrintAllPosesDelta()
{
    std::cout << "Global poses:\n";
    for (auto pose : mGlobalPoses)
    {
        std::cout << pose.first << "\n";
        pose.second->ShowImmutable();
    }
}

void cTripletSet::PrintAllPoses()
{
    std::cout << "Global poses:\n";
    for (auto pose : mGlobalPoses)
    {
        std::cout << pose.first << "\n";
        pose.second->Show();
    }
}

void cTripletSet::UpdateAllAffine()
{
    for (auto view_i : mAllViewMap)
    {
        AffineFromLocGlob(view_i.second);
    }
}

static double SetToMatDirecte(Eigen::MatrixXd& Mat)
{
   if (Mat.determinant() >0)
     return 1;

   for (int y=0; y<Mat.cols() ; y++)
       Mat(0,y) *= -1;
   return -1;
}


static Mat3d NearestRotation(Mat3d& rot)
{
    Eigen::MatrixXd U, V;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(rot,Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd SVec = svd.singularValues();
    U = svd.matrixU();
    V = svd.matrixV();

    Mat3d S = Mat3d::Zero(3,3);
    S(0,0) = SVec(0);
    S(1,1) = SVec(1);
    S(2,2) = SVec(2);

    double sign = SetToMatDirecte(U) * SetToMatDirecte(V);
    if (sign < 0)
        S(0,0) *= -1;
    V.transposeInPlace();


    Mat3d Mul;
    for (int i=0; i<3; i++)
        Mul(i,i) = (S(i,i) > 0) ? 1 : -1;
    

    return U * Mul * V;


}

void cTripletSet::AffineFromLocGlob(cNviewPoseX* view_i)
{
    bool TVIEW = (view_i->NbView()==3) ? 1 : 0;
    double NbDiv = 1/double(view_i->NbView());

    Mat3d alpha_ = Mat3d::Zero(3,3);
    Vec3d beta_  = Vec3d::Zero(3);


    //current global poses
    cPose*& Pose1 = this->Pose(view_i->View(0).Name());
    cPose*& Pose2 = this->Pose(view_i->View(1).Name());

    
    Vec3d C1 (Pose1->C()[0], Pose1->C()[1], Pose1->C()[2]);
    Vec3d C2 (Pose2->C()[0], Pose2->C()[1], Pose2->C()[2]);
    Vec3d C3 (TVIEW ? Pose(view_i->View(2).Name())->C()[0] : 0, 
              TVIEW ? Pose(view_i->View(2).Name())->C()[1] : 0, 
              TVIEW ? Pose(view_i->View(2).Name())->C()[2] : 0);


    //current local poses 
    Vec3d c1 (view_i->View(0).C()[0],
              view_i->View(0).C()[1],
              view_i->View(0).C()[2]);
    Vec3d c2 (view_i->View(1).C()[0], 
              view_i->View(1).C()[1],
              view_i->View(1).C()[2]);
    Vec3d c3 ( TVIEW ? view_i->View(2).C()[0] : 0, 
               TVIEW ? view_i->View(2).C()[1] : 0,
               TVIEW ? view_i->View(2).C()[2] : 0);

    

    //lambda 
    double L1 = euclid(c1 - c2) / euclid(C1 - C2);
    double L2 = euclid(c1 - c3) / euclid(C1 - C3);
    view_i->lambda()[0] = (L1 + L2) / ((TVIEW) ? 2 : 1);

    //rotation  Rk = rk_i * Ri^-1
    Mat3d Rk_i = view_i->View(0).R() * Pose1->R().transpose();
    Mat3d Rk_j = view_i->View(1).R() * Pose2->R().transpose();
    Mat3d MZ = Mat3d::Zero(3,3);
    Mat3d Rk_m = (TVIEW) ? (view_i->View(2).R() * (Pose(view_i->View(2).Name()))->R().transpose()) : MZ;
   

    Mat3d Rk_accum = (Rk_i+Rk_j+Rk_m)*NbDiv;
    view_i->alpha0() = NearestRotation(Rk_accum);

    //translation Ck = ck_i - lambda * Rk * Ci
    Vec3d RkC_1 = view_i->alpha0() * C1;
    Vec3d RkC_2 = view_i->alpha0() * C2;
    Vec3d RkC_3 = view_i->alpha0() * C3;
    Vec3d Ck_1 = c1 - view_i->lambda()[0] * RkC_1;
    Vec3d Ck_2 = c2 - view_i->lambda()[0] * RkC_2; 
    Vec3d Ck_3 = c3 - view_i->lambda()[0] * RkC_3; 

    //view_i->beta() = (Ck_1+Ck_2+Ck_3)*NbDiv;
    Vec3d beta = (Ck_1+Ck_2+Ck_3)*NbDiv;
    view_i->beta()[0] = beta(0);
    view_i->beta()[1] = beta(1);
    view_i->beta()[2] = beta(2);

}

void cTripletSet::LocalToGlobal(cNviewPoseX* pose,const int& view, Mat3d& R,Vec3d& C)
{
    //pose->Show();

    Mat3d r = pose->View(view).R();
    const double * cptr = pose->View(view).C();
    Vec3d c;
    c << cptr[0], cptr[1], cptr[2];

    Mat3d alpha = pose->alpha0();
    Vec3d beta;
    beta << pose->beta()[0], pose->beta()[1], pose->beta()[2];
    double lambda = pose->lambda()[0];

    R = alpha.inverse() * r;
    C = 1.0/lambda * alpha.inverse() * (c - beta);

    //lambda * alpha * C + beta = c 
}

void cTripletSet::SaveGlobalPoses(const std::string& filename)
{
    std::fstream eo_file;
    eo_file.open(filename.c_str(), std::istream::out);
    eo_file << std::fixed << std::setprecision(8) ;

    for (auto pose : mGlobalPoses)
    {
        if (pose.second->IsRefined())
        {
            cPose * pose0 = mGlobalPoses[pose.first];
            
            Mat3d dW;
            dW << 1,                 -pose0->Omega()[2], pose0->Omega()[1],
                  pose0->Omega()[2],    1,              -pose0->Omega()[0],
                 -pose0->Omega()[1], pose0->Omega()[0],             1;
         
         
            Mat3d newR = pose0->R()*dW;
            //pose0->Show();
         
            eo_file << pose.first << " " << newR(0,0) << " " << newR(0,1) << " " << newR(0,2)
                                  << " " << newR(1,0) << " " << newR(1,1) << " " << newR(1,2)
                                  << " " << newR(2,0) << " " << newR(2,1) << " " << newR(2,2)
                                  << " " << pose0->C()[0]
                                  << " " << pose0->C()[1]
                                  << " " << pose0->C()[2] << "\n";
        }
    }
    eo_file.close();

}

void cTripletSet::WriteGlobalPFromRelPAndSim(const std::string& filename)
{
    std::fstream eo_file;
    eo_file.open(filename.c_str(), std::istream::out);

    int cmpt=0;
    for (auto a3v : mAllViewMap)
    {
        //std::cout << a3v.first << "\n";
        int NumView = a3v.second->NbView();
        for (int view=0; view<NumView; view++)
        {
            std::string im_name =  std::to_string(cmpt) + "_" + a3v.second->View(view).Name();

            Mat3d R;
            Vec3d T;
            LocalToGlobal(a3v.second,view,R,T);

            eo_file << im_name << " " << R(0,0) << " " << R(0,1) << " " << R(0,2)
                               << " " << R(1,0) << " " << R(1,1) << " " << R(1,2)
                               << " " << R(2,0) << " " << R(2,1) << " " << R(2,2)
                               << " " << T[0]
                               << " " << T[1]
                               << " " << T[2] << "\n";

            //std::cout << im_name << "\n";
            //std::cout << "alpha=" << alpha << "\n beta=" << beta << "\n s=" << lambda << "\n";
            //std::cout << "r=" << r << "\nt" << c << "\n";
            //std::cout << "R=" << R << "\nT=" << T << "\n";
            std::cout << "cp DummyIm.tif " << im_name << "\n";

        }
        cmpt++;
        //getchar();

    }
    eo_file.close();
}
void cFilterTrip::CalcConservStats(std::vector<Vec3d>& data, 
                                   std::vector<int>& inlier_id,
                                   std::vector<int>& outliers_id)
{
    //compute geometric mean
    //compute residuals 
    //get 0.25 quantile as the new mean 
    ////compute residuals with the new mean 
    //compute std_dev and remove everything above 1sigma 

    double geom_cx=0,geom_cy=0,geom_cz=0;
    std::vector<double> residuals;

    //geometric mean 
    int num=0;
    //compute geometric center
    for (auto d : data)
    {
        geom_cx += d[0];
        geom_cy += d[1];
        geom_cz += d[2];
        num++;

    }

    geom_cx /= num;
    geom_cy /= num;
    geom_cz /= num;

    //residuals
    Vec3d geom_xyz;
    geom_xyz << geom_cx, geom_cy, geom_cz;
    residuals = CalcVecRes(data,geom_xyz);

    //compute quantile 
    double qperc = 0.25;
    double vmed = FindQuantile(residuals,qperc);

    std::vector<double>::iterator it_vmed = std::find (residuals.begin(), residuals.end(), vmed);
    double dist_vmed = std::distance(residuals.begin(),it_vmed);

    Vec3d geom_med_xyz;
    geom_med_xyz << data[dist_vmed][0], data[dist_vmed][1], data[dist_vmed][2];

    //compute residuals wrt to the new "mean" 
    residuals.clear();
    residuals = CalcVecRes(data,geom_med_xyz);
    
    //compute quantile on new residuals 
    double qperc2 = 2*qperc;
    double vmed2 = FindQuantile(residuals,qperc2);


    //standard deviation
    double std_dev=0;
    for (auto r : residuals)
        std_dev += r;
    
    std_dev /= num;
    std_dev = std::sqrt(std_dev);


    //keep motions with residuals <1sigma
    int cnt=0;
    for (auto r : residuals)
    {
        double res = std::sqrt(r);
        if (res< (std::sqrt(vmed2)+std_dev))
        {
            //std::cout << "inlier=" << res << "\n";
            inlier_id.push_back(cnt);
        }
        else
        {
            outliers_id.push_back(cnt);
            //std::cout << "   (!outlier=" << res << "!) \n";
        }
        cnt++;
    }

    std::cout << "new mean/std_dev= " << std::sqrt(vmed2) << "/" << std_dev << "\n";

}
std::vector<double> cFilterTrip::CalcVecRes(std::vector<Vec3d>& data,Vec3d& point)
{
    std::vector<double> residuals;

    for (auto d : data)
    {
        double res = std::pow(d[0] - point[0],2) + 
                     std::pow(d[1] - point[1],2) +
                     std::pow(d[2] - point[2],2);
        residuals.push_back(res);
        //std::cout << res << "\n";
    }
    return  residuals;
}

/* Calculate stats using rotation and translation */
bool cFilterTrip::CalcStats(std::vector<Vec3d>& dataT,std::vector<Mat3d>& dataR,
                            std::vector<int>& inlier_id,std::vector<int>& outlier_id)
{
    if (dataT.size() != dataR.size())
        return EXIT_FAILURE;

    int num =  dataT.size();

    //pick a random pose 
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_int_distribution<> distrib(0, num);
    int random = distrib(gen);

    //compute residuals wrt the random pose 
    std::vector<double> residuals;
    std::vector<double> residuals_cpy;
    for (int d=0; d<num; d++)
    {
        if (d!=random)
            residuals.push_back(DistanceRot(dataR[random],dataT[random],dataR[d],dataT[d]));
    }
    //get residual corresponding to 25% quantile 
    residuals_cpy = residuals;
    double qperc = 0.75;
    double res_quant = FindQuantile(residuals_cpy,qperc);
    
    //find the pose correponding to the qunatile
    std::vector<double>::iterator it_quant = std::find (residuals.begin(), residuals.end(), res_quant);
    double pos_quant = std::distance(residuals.begin(),it_quant);

    Vec3d quantT = dataT[pos_quant] ;
    Mat3d quantR = dataR[pos_quant];

    //compute residuals wrt to the new pose (recovered as 25%quantile)
    residuals.clear();
    residuals_cpy.clear();
    double vmin=1e10, vmax=-1e10;
    for (int d=0; d<num; d++)
    {
        double res = DistanceRot(quantR,quantT,dataR[d],dataT[d]);
        if (res>vmax) vmax = res;
        if (res<vmin) vmin = res;

        residuals.push_back(res);
    }

    //compute quantile on new residuals 
    //residuals_cpy = residuals;
    double qperc2 = 1*qperc;
    double res_quant2 = FindQuantile(residuals_cpy,qperc2);

    //standard deviation
    double std_dev=0;
    for (auto r : residuals)
        std_dev += (r-res_quant2)*(r-res_quant2);
    
    std_dev /= num;
    std_dev = std::sqrt(std_dev);

    std::cout << "quant/std_dev=" << res_quant2 << " " << std_dev 
              << ", min/max=" << vmin << " " << vmax << "\n"; 

    //keep motions with residuals <1sigma
    int cnt=0;
    for (auto r : residuals)
    {
        double res = r;
        if (res< (res_quant2+std_dev))
        {
            //std::cout << "inlier=" << res << "\n";
            inlier_id.push_back(cnt);
        }
        else
        {
            outlier_id.push_back(cnt);
            //std::cout << "   (!outlier=" << res << "!) \n";
        }
        cnt++;
    }

    //getchar();

    return EXIT_SUCCESS;
}

double cFilterTrip::FindQuantile(std::vector<double>& v, double p)
{
    std::sort (v.begin(), v.end());

    std::vector<double>::iterator start = v.begin();
    std::vector<double>::iterator end = v.end();
    
    const std::size_t pos = p * std::distance(start, end);

    std::vector<double>::iterator quantile = start;
    std::advance(quantile,pos);

    std::nth_element(start,quantile,end);

    //std::cout << "median=" << *quantile << " " <<  "\n";

    /*for (auto it : v)
        std::cout << it << " \n"; */
    return (*quantile);
}

double euclid(Vec3d pt)
{
    return std::sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]);
}

double cFilterTrip::scal (const Vec3d& p1,const Vec3d& p2)
{
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

double cFilterTrip::DistBase(Vec3d B1,Vec3d B2)
{
      if (scal(B1,B2) < 0) B2 = - B2;
      double D1 = euclid(B1); 
      double D2 = euclid(B2); 

      if (D1 > D2)
         B1 = B1 * (D2/D1);
      else
         B2 = B2 * (D1/D2);

      return euclid(B1-B2);
}

double cFilterTrip::DistanceRot(const Mat3d& R1,const Vec3d& C1, const Mat3d& R2,const Vec3d& C2)
{
      Mat3d Dif = R1 - R2;
      double DistRot = Dif.norm();
      double DistTr =  DistBase(C1,C2) * 0.3 ;

      //std::cout << "distTr=" << DistTr << ", DistR=" << DistRot << "\n";
      return DistTr + DistRot;
}
void cFilterTrip::SaveViews(cTripletSet& Tri,
                            const std::list<std::string>& views,
                            std::string filename)
{
    std::fstream in_file;
    in_file.open(filename.c_str(), std::istream::out);
   
    for(auto motion : views)
    {
        std::vector<std::string> names = Tri.DecompViewNames(motion);
        
        cNviewPoseX* NViewP = Tri.mAllViewMap[motion];
        
        in_file << NViewP->NbView() << " ";
        for (int vi=0; vi<NViewP->NbView(); vi++)
            in_file << names[vi] << " "; 


        for (int vi=1; vi<NViewP->NbView(); vi++)
        {
            
            Mat3d R =  NViewP->View(vi).R();
            double  * C = NViewP->View(vi).C();

            in_file << R(0,0) << " " << R(0,1) << " " << R(0,2) << " "
                    << R(1,0) << " " << R(1,1) << " " << R(1,2) << " " 
                    << R(2,0) << " " << R(2,1) << " " << R(2,2) << " "
                    << C[0] << " " 
                    << C[1] << " "
                    << C[2] << " ";

        
        }
        in_file << "\n";

    }
    in_file.close();
}

cFilterTrip::cFilterTrip(std::string views,
                         std::string simil,
                         bool do_only_pred,
                         std::string outgpose,
                         std::string gpose,
                         std::string filtered_views)
{
    std::cout << "Relative motion filter" << "\n";

    ONLY_PREDICTION = do_only_pred;

    cTripletSet aTriSet(views,simil,gpose);
    aTriSet.ReadViews();
    aTriSet.ReadSimGlobal();
    //aTriSet.PrintAllViews();

    if (ONLY_PREDICTION)
    {
        aTriSet.WriteGlobalPFromRelPAndSim(outgpose);
        getchar();
    }
    else //do filter
    {

        //iterate over all motions
        std::map<std::string,std::vector<Vec3d>>          map_pose_predC;
        std::map<std::string,std::vector<Mat3d>>          map_pose_predR;
        //std::vector<Mat3d>        pose_predR; for now I'll identifty outliers by looking at C only
        std::map<std::string,std::vector<std::string>>  map_pose_relmotion;
        std::list<std::string> list_rmotion_inliers;
        std::list<std::string> list_rmotion_outliers;
       
        //iterate over relative motions and compute global poses using relative poses
        for (auto a3 : aTriSet.mAllViewMap)
        {
            int NumView = a3.second->NbView();
            for (int view=0; view<NumView; view++)
            {
 
                Mat3d R;
                Vec3d T;
                aTriSet.LocalToGlobal(a3.second,view,R,T);
 
                //std::cout << a3.second->View(view).Name() << " " << T << "\n";
                if (!DicBoolFind(map_pose_predC,a3.second->View(view).Name()))
                {
                    map_pose_predC[a3.second->View(view).Name()] = std::vector<Vec3d>();
                    map_pose_predC[a3.second->View(view).Name()].push_back(T);
 
                    map_pose_predR[a3.second->View(view).Name()] = std::vector<Mat3d>();
                    map_pose_predR[a3.second->View(view).Name()].push_back(R);
 
                    map_pose_relmotion[a3.second->View(view).Name()] = std::vector<std::string>();
                    map_pose_relmotion[a3.second->View(view).Name()].push_back(a3.first);
                }
                else
                {
                    map_pose_predC[a3.second->View(view).Name()].push_back(T);
                    map_pose_predR[a3.second->View(view).Name()].push_back(R);
                    map_pose_relmotion[a3.second->View(view).Name()].push_back(a3.first);
                }
            }
        }
 
        //collect all outliers first 
        for (auto p : map_pose_predC)
        {
            std::vector<std::string> rel_poses_i = map_pose_relmotion[p.first];
            
            std::vector<int> inliers;
            std::vector<int> outliers;
 
            CalcStats(p.second,map_pose_predR[p.first],inliers,outliers);
            if (outliers.size())
            {
 
                for(auto out : outliers)
                {
                    std::list<std::string>::iterator f_it = std::find(list_rmotion_outliers.begin(),
                              list_rmotion_outliers.end(), rel_poses_i[out]);
 
                    std::cout << "rel_poses_i[in]=" << rel_poses_i[out] << ", out=" << out << "\n";
                    if ( list_rmotion_outliers.size()==0)
                    {
                        list_rmotion_outliers.push_back(rel_poses_i[out]); 
                    }
                    else if(f_it == list_rmotion_outliers.end())
                    {
                        list_rmotion_outliers.push_back(rel_poses_i[out]); 
                    }
 
                }
            }
        }
 
 
        //identify inliers and save
        int cnt_inliers=0, cnt_all=0;
        for (auto p : map_pose_predC)
        {
            std::cout << p.first << " ";
            std::vector<std::string> rel_poses_i = map_pose_relmotion[p.first];
            
            std::vector<int> inliers;
            std::vector<int> outliers;
 
            //CalcConservStats(p.second,inliers,outliers);
            CalcStats(p.second,map_pose_predR[p.first],inliers,outliers);
 
 
            if (inliers.size())
            {
                cnt_all += p.second.size();
 
                std::cout << " inliers=" << inliers.size() << "/" << p.second.size() << "\n";
 
                for (auto in : inliers)
                {
                    std::list<std::string>::iterator f_it = std::find(list_rmotion_inliers.begin(),
                              list_rmotion_inliers.end(), rel_poses_i[in]);
                    std::list<std::string>::iterator fo_it = std::find(list_rmotion_outliers.begin(),
                              list_rmotion_outliers.end(), rel_poses_i[in]);
 
 
                    if (fo_it==list_rmotion_outliers.end())
                    {
                        cnt_inliers++;
                        
                        std::cout << "rel_poses_i[in]=" << rel_poses_i[in] << ", in=" << in << "\n";
                        if (list_rmotion_inliers.size()==0) 
                        {
                            list_rmotion_inliers.push_back(rel_poses_i[in]); 
                        }
                        else if(f_it == list_rmotion_inliers.end())
                        {
                            list_rmotion_inliers.push_back(rel_poses_i[in]); 
                        }
                    }
                }
            }
        }
        std::cout << "Filter acceptance rate: " << double(cnt_inliers)/(cnt_all)*100 << "\n";
 
        SaveViews(aTriSet,list_rmotion_inliers,"inlier-"+filtered_views);
        SaveViews(aTriSet,list_rmotion_outliers,"outlier-"+filtered_views);


    }

}
