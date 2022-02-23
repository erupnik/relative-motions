#ifndef _TRI_UTILS
#define _TRI_UTILS

#include "all.h"
#include "cov_in_motions.h"

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
            mglob_p_file(gpose) {}
        ~cTripletSet(){};

        bool ReadViews();
        bool ReadSimGlobal();
        bool ReadGlobalPoses();

 
        void PrintAllViews();

    private:
        friend class cAppCovInMotion;

        bool ReadRotTrS(FILE* fptr,Mat3d& alpha,Vec3d& beta,double& s);
        
        std::string mviews_file;
        std::string msimil_file;
        std::string mglob_p_file;
 
        /* relative motions and global poses */
        std::map<std::string,cNviewPoseX*>   mAllViewMap_;//ok
 
        std::map<std::string,cPose*>      mGlobalPoses;//ok
};
/*
class cFilterTrip
{
    public:
        cFilterTrip();
}*/

#endif //_TRI_UTILS
