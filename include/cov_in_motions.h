#ifndef _COV_IN_MOT_H
#define _COV_IN_MOT_H

#include "all.h"
#include "options.h"
#include "triplet_utils.h"
#include "bundle.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>

class cPose;
class cNviewPose;
class cAppCovInMotion;
class cTripletSet;

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
    virtual double*   C_immutable()  =0;
    virtual Vec3d&    C_()  =0;
    virtual void      Show() const =0;
    virtual const std::string& Name() =0;

    virtual double*   Omega() =0;
    virtual double*   Omega_immutable() =0;
    virtual double*   Cov_C() =0;
    virtual double*   Cov_Omega() =0;

    virtual ~cPoseGen(){}


};
class cPoseBasic : public cPoseGen //basic pose on eigen matrices
{
    public:
        cPoseBasic(Mat3d R,Vec3d C,std::string Name) :
                    mR(R), mC(C), mName(Name) {}
        ~cPoseBasic(){}

        Vec3d&    C_() {return mC;}
        double*   C()  {throw "Not implemented!"; return 0;}
        double*   C_immutable()  {throw "Not implemented!"; return 0;}
        Mat3d&    R()  {return mR;}

        double*   Omega() {throw "Not implemented!"; return 0;}
        double*   Omega_immutable() {throw "Not implemented!"; return 0;}
        double*   Cov_C() {throw "Not implemented!"; return 0;}
        double*   Cov_Omega() {throw "Not implemented!"; return 0;}

        const std::string& Name() {return mName;}
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
             mR(R), mC(C), mC_immutable(new double[3]), mName(Name), 
             mOmega(new double[3]), mOmega_immutable(new double[3]),
             mCov_C(new double [9]), mCov_Omega(new double [9]),
             FLAG_Refined(false) {
             mC_immutable[0] = mC[0];//copy
             mC_immutable[1] = mC[1];//copy 
             mC_immutable[2] = mC[2];//copy 
             memset(mOmega, 0, sizeof(int)*3);
             memset(mOmega_immutable, 0, sizeof(int)*3);
             memset(mCov_C, 0, sizeof(int)*9);
             memset(mCov_Omega, 0, sizeof(int)*9);
             mC_ << mC[0], mC[1], mC[2];
             mOmega[0] = 0;//_ETA* double(std::rand()) / RAND_MAX;
             mOmega[1] = 0;//_ETA*double(std::rand()) / RAND_MAX;
             mOmega[2] = 0;//_ETA*double(std::rand()) / RAND_MAX;
             mOmega_immutable[0] = 0;//_ETA* double(std::rand()) / RAND_MAX;
             mOmega_immutable[1] = 0;//_ETA*double(std::rand()) / RAND_MAX;
             mOmega_immutable[2] = 0;//_ETA*double(std::rand()) / RAND_MAX;
             
           }
        ~cPose(){ 
            delete mC;
            delete mC_immutable;
            delete mOmega; 
            delete mOmega_immutable; 
            delete mCov_C; 
            delete mCov_Omega; }

        Vec3d&    C_();//not implemented
        double*   C() {return mC;}
        double*   C_immutable() {return mC_immutable;}
        Mat3d&    R() {return mR;}
        double*   Omega() {return mOmega;}
        double*   Omega_immutable() {return mOmega_immutable;}

        double*   Cov_C() {return mCov_C;}
        double*   Cov_Omega() {return mCov_Omega;}

        const std::string& Name() {return mName;}

        bool IsRefined() { return FLAG_Refined; }
        void SetRefined() { FLAG_Refined=true; }

        void Show() const;

    private:
        double*  mC;
        double*  mC_immutable;
        Vec3d    mC_;//not used
        Mat3d    mR;
        double * mOmega; //tiny rotation
        double * mOmega_immutable; //tiny rotation

        double * mCov_C;
        double * mCov_Omega;

        bool FLAG_Refined;
        std::string mName;
};

class cHessianGradient
{
  public:

      virtual void printH() const =0;
      virtual void printG() const =0;
};


template <typename T,typename U>
class cHessianGradientN : public cHessianGradient
{

  public:
    cHessianGradientN(T hessian,U gradient) :
                mH(hessian), mg(gradient) {}
    void printH() const  {std::cout << mH << "\n";}
    void printG() const {std::cout << mg << "\n";}

    T& H_() {return mH;}
    U& g_() {return mg;}

  private:
    T mH;
    U mg;

};

class cHessianGradientX : public cHessianGradient
{

  public:
    cHessianGradientX(MatXd hessian,VecXd gradient) :
                mH(hessian), mg(gradient) {}
    void printH() const  {std::cout << mH << "\n";}
    void printG() const {std::cout << mg << "\n";}

    MatXd& H_() {return mH;}
    VecXd& g_() {return mg;}

  private:
    MatXd mH;
    VecXd mg;

};

class cNviewPose
{
  public:
    virtual cPoseGen& View(int)  =0;
    //virtual cHessianGradient& Hg_() {return *Hg;}
    virtual cHessianGradient& Hg_() =0;
    virtual int NbView() const =0;

    double& lambda() {return affine_trafo.lambda;}
    Mat3d&  alpha()  {return affine_trafo.alpha;}
    Vec3d&  beta()   {return affine_trafo.beta;}

    void PrintAlpha() {std::cout << affine_trafo.alpha << "\n";}
    void PrintBeta() {std::cout << affine_trafo.beta << "\n";}
    void PrintLambda() {std::cout << affine_trafo.lambda << "\n";}

  private:
    sAffine affine_trafo;

};

template <typename T, typename U>
class cNviewPoseT
{
  public:
    cNviewPoseT<T,U>(cPoseGen* v1,
                     cPoseGen* v2,
                     cPoseGen* v3,
                     cHessianGradientN<T,U>* aHg) :
                            mView1(v1),
                            mView21(v2),
                            mView31(v3),
                            mHg(aHg),
                            mNbV( (&(*mView31)==NULL) ? 2 : 3 ),
                            _COV_PROP(false) {}

    cPoseGen& View(int NbV)  {
         if (NbV==0) return *mView1; else if (NbV==1) return *mView21; else return *mView31; }

    bool propagate_cov();
    bool IS_COV_PROP() {return _COV_PROP;}

    cHessianGradientN<T,U>& Hg_() {return *mHg;};
    int NbView() const {return mNbV;};

    double& lambda() {return affine_trafo.lambda;}
    Mat3d&  alpha()  {return affine_trafo.alpha;}
    Vec3d&  beta()   {return affine_trafo.beta;}

    void PrintAlpha() {std::cout << affine_trafo.alpha << "\n";}
    void PrintBeta() {std::cout << affine_trafo.beta << "\n";}
    void PrintLambda() {std::cout << affine_trafo.lambda << "\n";}

  private:
    sAffine affine_trafo;

    cPoseGen*                            mView1; // Identity
    cPoseGen*                            mView21;
    cPoseGen*                            mView31;
    cHessianGradientN<T,U>*              mHg;

    int                                  mNbV;
    bool                                 _COV_PROP;

};

class cNviewPoseX
{
  public:
    cNviewPoseX     (cPoseGen* v1,
                     cPoseGen* v2,
                     cPoseGen* v3,
                     cHessianGradientX* aHg) :
                            mView1(v1),
                            mView21(v2),
                            mView31(v3),
                            mHg(aHg),
                            mNbV( (&(*mView31)==NULL) ? 2 : 3 ),
                            _COV_PROP(false),
                            _FLAG_OUTLIER(false),
                            _INIT(false) {}

    cPoseGen& View(int NbV)  {
         if (NbV==0) return *mView1; else if (NbV==1) return *mView21; else return *mView31; }

    bool propagate_cov(Mat3d&,Mat3d&,Mat3d& );
    bool IS_COV_PROP() {return _COV_PROP;}

    cHessianGradientX& Hg_() {return *mHg;};
    int NbView() const {return mNbV;};

    double& lambda() {return affine_trafo.lambda;}
    Mat3d&  alpha()  {return affine_trafo.alpha;}
    Vec3d&  beta()   {return affine_trafo.beta;}

    void PrintAlpha() {std::cout << affine_trafo.alpha << "\n";}
    void PrintBeta() {std::cout << affine_trafo.beta << "\n";}
    void PrintLambda() {std::cout << affine_trafo.lambda << "\n";}

    bool Init() {return _INIT;}
    void SetInit() {_INIT=true;}
    void SetOutlier() {_FLAG_OUTLIER=true;}
    void Show() {mView1->Show(); mView21->Show(); if (mNbV==3) mView31->Show();}

  private:
    sAffine affine_trafo;

    cPoseGen*                            mView1; // Identity
    cPoseGen*                            mView21;
    cPoseGen*                            mView31;
    cHessianGradientX*                   mHg;

    int                                  mNbV;
    bool                                 _COV_PROP;

    bool                                 _FLAG_OUTLIER;
    bool                                 _INIT;

};

class cAppCovInMotion
{
  public:
    /*cAppCovInMotion(const std::string& ,
                    const std::string&,
                    const std::string&,
                    const std::string&,
                    const std::string&,
                    const bool,
                    const std::string& );*/
    cAppCovInMotion(const InputFiles&,
                    const LocalBundleOptions& ,
                    const GlobalBundleOptions& ,
                    const bool );
    ~cAppCovInMotion();

    bool ReadFeatures(std::string);
    void PrintAllPts3d();


    bool BuildProblem_(cNviewPoseX*&,std::string);
    void InitCovariances(cNviewPoseX*&,std::string);

    void SetCeresOptions(ceres::Solver::Options&);
    void SetMinimizer(ceres::Solver::Options&);
    void SetMinimizerGlobal(ceres::Solver::Options&);
    bool OptimizeRelMotions();
    bool OptimizeRelMotionsGlobally();
    bool OptimizeGlobally();

    void WriteLocalCovs();
    void WriteToPLYFile(const std::string& filename,
                        const std::string& viewName);

  private:
    void MapJacToEigMat(Eigen::SparseMatrix<double>&, MatXd&, int offset);
    Eigen::SparseMatrix<double> MapJacToSparseM(CRSMatrix*);
    Eigen::SparseMatrix<double> MatToSparseM(Eigen::MatrixXd&);

    LocalBundleOptions  m_lba_opts;
    GlobalBundleOptions m_gba_opts;

    /*std::string mviews_file;
    std::string mfeats_file;
    std::string msimil_file;
    std::string mglob_p_file;
    std::string mout_p_file;*/

    bool        GET_COVARIANCES;

    /* relative motions and global poses */
    cTripletSet                        * mTriSet;

    /* view id is the string; the vector stores features */
    std::map<std::string,std::vector<std::vector<Vec2d>>* > mFeatViewMap;//ok
    std::map<std::string,std::vector<double*>*> mFeat3dMap;//ok
};


#endif //_COV_IN_MOT_H
