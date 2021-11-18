#ifndef _COV_IN_MOT_H
#define _COV_IN_MOT_H

#include "bundle.h"


class cPose;
class cNviewPose;
class c2viewPose;
class c3viewPose;
class cAppCovInMotion;

template <class TVal,class TCont>
bool DicBoolFind(const TCont & aCont,const TVal & aVal)
{
   return aCont.find(aVal) != aCont.end();
}

struct sAffine
{
  double lambda;
  Mat3d alpha;
  Vec3d beta;
};

class cPoseGen
{
  public:
    virtual Mat3d&    R()  =0;
    virtual double*   C()  =0;
    virtual Vec3d&    C_()  =0;
    virtual void      Show() const =0;

    virtual double*   Omega() =0;
    virtual double*   Cov_C() =0;
    virtual double*   Cov_Omega() =0;


};
class cPoseBasic : public cPoseGen //basic pose on eigen matrices
{
    public:
        cPoseBasic(Mat3d R,Vec3d C,std::string Name) :
                    mR(R), mC(C), mName(Name) {}
        ~cPoseBasic(){}

        Vec3d&    C_() {return mC;}
        double*   C() {throw "Not implemented!"; return 0;}
        Mat3d&    R() {return mR;}

        double*   Omega() {throw "Not implemented!"; return 0;}
        double*   Cov_C() {throw "Not implemented!"; return 0;}
        double*   Cov_Omega() {throw "Not implemented!"; return 0;}

        void Show() const;

    private:
        Vec3d mC;
        Mat3d mR;

        std::string mName;
};

class cPose : public cPoseGen //any pose, relative or absolute
{
    public:
        cPose(Mat3d R,double* C,std::string& Name) :
             mR(R), mC(C), mName(Name), mOmega(new double[3]),
             mCov_C(new double [9]), mCov_Omega(new double [9]) {
             memset(mOmega, 0, sizeof(int)*3);
             memset(mCov_C, 0, sizeof(int)*9);
             memset(mCov_Omega, 0, sizeof(int)*9);
             mC_ << mC[0], mC[1], mC[2];
           }
        ~cPose(){}

        Vec3d&    C_();//not implemented
        double*   C() {return mC;}
        Mat3d&    R() {return mR;}
        double*   Omega() {return mOmega;}

        double*   Cov_C() {return mCov_C;}
        double*   Cov_Omega() {return mCov_Omega;}

        const std::string& Name() {return mName;}

        void Show() const;

    private:
        double*  mC;
        Vec3d    mC_;//not used
        Mat3d    mR;
        double * mOmega; //tiny rotation

        double * mCov_C;
        double * mCov_Omega;

        std::string mName;
};

class cHessianGradient
{
  public:

      virtual void printH() const =0;
      virtual void printG() const =0;
};

//two-view hessiand and gradient
class cHessianGradient2 : public cHessianGradient
{
  public:
    cHessianGradient2(Mat6d hessian,Vec6d gradient) :
                mH(hessian), mg(gradient) {}
    void printH() const;
    void printG() const {std::cout << mg;}

  private:
    Mat6d mH;
    Vec6d mg;
};

//three-view hessiand and gradient
class cHessianGradient3 : public cHessianGradient
{

  public:
    cHessianGradient3(Mat12d hessian,Vec12d gradient) :
                mH(hessian), mg(gradient) {}
    void printH() const;
    void printG() const {std::cout << mg;}

  private:
    Mat12d mH;
    Vec12d mg;

};

class cNviewPose
{
  public:
    virtual cPoseGen& View(int)  =0;
    //virtual cHessianGradient& Hg_() {return *Hg;}
    virtual cHessianGradient& Hg_() =0;

    double& lambda() {return affine_trafo.lambda;}
    Mat3d&  alpha()  {return affine_trafo.alpha;}
    Vec3d&  beta()   {return affine_trafo.beta;}

    void PrintAlpha() {std::cout << affine_trafo.alpha << "\n";}
    void PrintBeta() {std::cout << affine_trafo.beta << "\n";}
    void PrintLambda() {std::cout << affine_trafo.lambda << "\n";}

  private:
    sAffine affine_trafo;

};

class c2viewPose_ : public cNviewPose
{
  public:
    /*c2viewPose_(cPoseGen* v1,cPoseGen* v2) :
                mView1(v1), mView21(v2)
                { mHg = new  cHessianGradient2(Mat6d::Zero(),Vec6d::Zero());}*/
    c2viewPose_(cPoseGen* v1,cPoseGen* v2,cHessianGradient* aHg) :
                            mView1(v1), mView21(v2), mHg(aHg) {}

    cPoseGen& View(int NbV)  {
         if (NbV==0) return *mView1; else return *mView21; }


    cHessianGradient& Hg_() {return *mHg;}

  private:
    cPoseGen*                           mView1;
    cPoseGen*                           mView21;

    cHessianGradient* mHg;

};
class c2viewPose : public cNviewPose
{
  public:
    c2viewPose(cPose v1,cPose v2) : mView1(v1), mView21(v2) {}
    cPose& View1() {return mView1;}
    cPose& View21() {return mView21;}
    cPose& View(int NbV)  {
         if (NbV==0) return mView1; else return mView21; }

    cHessianGradient& Hg_() {return *Hg;}

  private:
    cPose                            mView1; // Identity
    cPose                            mView21;

    cHessianGradient* Hg;
    /*Mat12d a_i;
    Vec12d b_i;*/

};

class c3viewPose_ : public cNviewPose
{
  public:
    c3viewPose_(cPoseGen* v1,cPoseGen* v2, cPoseGen* v3,cHessianGradient* aHg) :
                    mView1(v1), mView21(v2), mView31(v3), mHg(aHg) {}

    cPoseGen& View(int NbV)  {
         if (NbV==0) return *mView1; else if (NbV==1) return *mView21; else return *mView31; }

    cHessianGradient& Hg_() {return *mHg;}

  private:
    cPoseGen*                           mView1;
    cPoseGen*                           mView21;
    cPoseGen*                           mView31;

    cHessianGradient* mHg;
};

class c3viewPose : public cNviewPose
{
  public:
    c3viewPose(cPose v1,cPose v2,cPose v3) : mView1(v1), mView21(v2), mView31(v3) {}
    cPose& View1() {return  mView1;}
    cPose& View21() {return  mView21;}
    cPose& View31() {return  mView31;}
    cPose& View(int NbV) {
         if (NbV==0) return mView1; else if (NbV==1) return mView21; else return mView31; }

    cHessianGradient& Hg_() {return *Hg;}

  private:
    cPose mView1; // Identity
    cPose mView21;
    cPose mView31;

    cHessianGradient* Hg;
    /*Mat18d a_i;
    Vec18d b_i;*/

    //Mat18d A_i;
    //Vec18d B_i;

};

class cAppCovInMotion
{
  public:
    cAppCovInMotion(const std::string& ,
                    const std::string&,
                    const std::string&,
                    const std::string& );

  private:
    bool ReadFeatures();
    bool ReadViews();
    bool ReadSimGlobal();
    bool ReadGlobalPoses();
    bool ReadRotTrS(FILE* fptr,Mat3d& alpha,Vec3d& beta,double& s);

    void PrintAllViews();
    void PrintAllPts3d();

    bool BuildProblem(const std::string&,int);
    void SetCeresOptions(ceres::Solver::Options&);
    void SetMinimizer(ceres::Solver::Options&);
    bool OptimizeRelMotions();
    bool OptimizeGlobally();

    void WriteToPLYFile(const std::string& filename,
                        const std::string& viewName);

    Eigen::SparseMatrix<double> MapJacToSparseM(CRSMatrix*);

    std::string mviews_file;
    std::string mfeats_file;
    std::string msimil_file;
    std::string mglob_p_file;

    /* relative motions and global poses */
    std::map<std::string,c2viewPose*> m2ViewMap;//ok
    std::map<std::string,c3viewPose*> m3ViewMap;//ok
    std::map<std::string,cPose*>      mGlobalPoses;//ok

    /* view id is the string; the vector stores features */
    std::map<std::string,std::vector<std::vector<Vec2d>>* > mFeatViewMap;//ok
    std::map<std::string,std::vector<double*>*> mFeat3dMap;//ok
};


#endif //_COV_IN_MOT_H
