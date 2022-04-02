#ifndef _TRI_UTILS
#define _TRI_UTILS

#include "all.h"
#include "cov_in_motions.h"
#include <random>

class cPose;
class cNviewPoseX;

class cTripletSet
{
    public: 
        cTripletSet(std::string views,
                    std::string simil,
                    std::string gpose="") : 
            mviews_file(views),  
            msimil_file(simil),
            mglob_p_file(gpose) { std::cout << "Initialization\n"; }
        ~cTripletSet(){};

        bool ReadViews();
        bool ReadSimGlobal();
        bool ReadGlobalPoses();
        void PrintAllViews();
        void PrintAllPoses();

        void LocalToGlobal(cNviewPoseX*,const int& view,Mat3d&,Vec3d&);
        void AffineFromLocGlob(cNviewPoseX*); 
        void UpdateAllAffine();

        void WriteGlobalPFromRelPAndSim(const std::string& );
        void SaveGlobalPoses(const std::string&);
        
        cPose*& Pose(std::string name) {return mGlobalPoses[name];}        

    private:
        friend class cAppCovInMotion;
        friend class cFilterTrip;

        bool ReadRotTrS(FILE* fptr,Mat3d& alpha,Vec3d& beta,double& s);
        std::vector<std::string> DecompViewNames(std::string&);

        std::string mviews_file;
        std::string msimil_file;
        std::string mglob_p_file;
 
        /* relative motions and global poses */
        std::map<std::string,cNviewPoseX*>   mAllViewMap;//ok
 
        std::map<std::string,cPose*>      mGlobalPoses;//ok
};

class cFilterTrip
{
    public:
        cFilterTrip(std::string views,
                    std::string simil,
                    bool do_only_pred,
                    std::string outgpose="",
                    std::string gpose="",
                    std::string filterd_views="");

    private:
        bool CalcStats(std::vector<Vec3d>&, std::vector<Mat3d>&,
                       std::vector<int>&,std::vector<int>&);
        void CalcConservStats(std::vector<Vec3d>&, std::vector<int>&,std::vector<int>&);
        double DistanceRot(const Mat3d& R1,const Vec3d& C1, const Mat3d& R2,const Vec3d& C2);
        double DistBase(Vec3d B1,Vec3d aB2);
        double scal (const Vec3d& p1,const Vec3d& p2);

        void SaveViews(cTripletSet&,const std::list<std::string>&,std::string);

        std::vector<double> CalcVecRes(std::vector<Vec3d>&,Vec3d&);
        double FindQuantile(std::vector<double>& v,double p);

        bool ONLY_PREDICTION;
};

#endif //_TRI_UTILS
