#ifndef _OPTIONS
#define _OPTIONS

#include "all.h"
#include <random>

struct GlobalBundleOptions
{
    public:
        // square root PDS for rotations
        double _ROT_PDS=0.01 ;

        // square root PDS for perspective centers
        double _C_PDS=0.5 ;

        // huber loss threshold for similarities 
        // (in the unit of the local motion frame)
        double _HUBER_S=0.001;

        // huber loss threshold for poses (soft constraint)
        // (in the units of the global poses)
        double _HUBER_P=0.001;

        //use inner iterations
        bool _INNER_ITER=false;

    private:
};

// Local bundle adjustment
// stereopolis small:
//    _FEAT_PDS=100 
//    _ROT_PDS=1 
//    _T_PDS=1 
//    _HUBER_S=0.1
//    _ROT_PDS_g=0.01 ;
//    _T_PDS_g=0.5 ;
//     _HUBER_g=0.001;
//
struct LocalBundleOptions 
{
    public:
        // turn on/off the covariance propagation
        bool _RUN_PROP = true;

        // squared root PDS for normalized features (pix/F)
        double _FEAT_PDS=1000;
        
        // square root PDS for initial rotations 
        double _ROT_PDS=1 ;
        
        // square root PDS for initial perspective centers
        double _C_PDS=1 ;
        
        // huber loss threshold
        // (in the unit of the normalized features)
        double _HUBER=1; 
      
        // maximum Hessian trace value, 
        // motions above this value are not considered
        double _TRACE_H = 10e10;

        // print covariances to a file 
        bool _WRITE_COV = false;

        // covariance file key
        std::string _KEY = "_log";

        //use inner iterations
        bool _INNER_ITER=false;

    private:

};

struct InputFiles
{
    public:
        // file with feature tracks 
        std::string tracks_file;

        // file with views (relative motions)
        std::string views_file;

        // file with local to global frame similarities 
        std::string similarities_file;

        // file with initial global poses 
        std::string global_poses_file;

        // output, refined poses file
        std::string output_poses_file = "output_poses.txt";

    private:
};

#endif //_OPTIONS
