#ifndef _COV_IN_MOT_H
#define _COV_IN_MOT_H

#include "bundle.h"

class cPose;
class c2viewPose;
class c3viewPose;
class cAppCovInMotion;

template <class TVal,class TCont>
bool DicBoolFind(const TCont & aCont,const TVal & aVal)
{
   return aCont.find(aVal) != aCont.end();
}

class cPose //any pose, relative or absolute
{
    public:
        cPose(Mat3d R,double* C,std::string& Name) :
             mR(R), mC(C), mName(Name), mOmega(new double[3]),
             mCov_C(new double [9]), mCov_Omega(new double [9]) {
             memset(mOmega, 0, sizeof(int)*3);
             memset(mCov_C, 0, sizeof(int)*9);
             memset(mCov_Omega, 0, sizeof(int)*9);
           }
        ~cPose(){}

        double*   C() {return mC;}
        Mat3d&    R() {return mR;}
        double*   Omega() {return mOmega;}

        double*   Cov_C() {return mCov_C;}
        double*   Cov_Omega() {return mCov_Omega;}

        const std::string& Name() {return mName;}

        void Show() const;

    private:

        double * mC;
        Mat3d mR;
        double * mOmega; //tiny rotation

        double * mCov_C;
        double * mCov_Omega;

        //int         mId; //Id that will link with mKP and bundler EO init ?
        std::string mName;
};

 
class c2viewPose
{
  public:
    c2viewPose(cPose v1,cPose v2) : mView1(v1), mView21(v2) {}
    cPose& View1() {return mView1;}
    cPose& View21() {return mView21;}
    cPose& View(int NbV) {
         if (NbV==0) return mView1; else return mView21; }

  private:
    cPose                            mView1; // Identity
    cPose                            mView21;
};

class c3viewPose
{
  public:
    c3viewPose(cPose v1,cPose v2,cPose v3) : mView1(v1), mView21(v2), mView31(v3) {}
    cPose& View1() {return  mView1;}
    cPose& View21() {return  mView21;}
    cPose& View31() {return  mView31;}
    cPose& View(int NbV) {
         if (NbV==0) return mView1; else if (NbV==1) return mView21; else return mView31; }

  private:
    cPose mView1; // Identity
    cPose mView21;
    cPose mView31;

};

class cAppCovInMotion
{
  public:
    cAppCovInMotion(const std::string& ,const std::string& );

  private:
    bool ReadFeatures();
    bool ReadViews();
    void PrintAllViews();
    void PrintAllPts3d();

    void BuildProblem(int aNbCam);
    void SetCeresOptions(ceres::Solver::Options&);
    void SetMinimizer(ceres::Solver::Options&);
    bool Optimize();

    void WriteToPLYFile(const std::string& filename,
                        const std::string& viewName);

    std::string mviews_file;
    std::string mfeats_file;

    std::map<std::string,c2viewPose*> m2ViewMap;
    std::map<std::string,c3viewPose*> m3ViewMap;

    /* view id is the string; the vector stores features */
    std::map<std::string,std::vector<std::vector<Vec2d>>* > mFeatViewMap;
    std::map<std::string,std::vector<double*>*> mFeat3dMap;
};


#endif //_COV_IN_MOT_H
