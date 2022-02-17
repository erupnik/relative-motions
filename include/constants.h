#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

/* Project-wide contants */

// Local bundle adjustment
const double _FEAT_PDS=1.00;
const double _ROT_PDS=1.0;//0.17*2;
const double _T_PDS=2.0;//0.1;

const double _ETA_up=0.1;
const double _ETA=0.001;

const bool _GAUGE_FIRST_CAM_FIX=false;
const bool _GAUGE_BASE_FIX=false;
const bool _GAUGE_RAPPEL_POSES=true;


/*
 *  15 - no fix                     => where 3DOF go?
 *  10 - fix first cam (6 params)   => where 2DOF go?
 *  9 - fix base x,  fix 1 cam (7)  => where 2DOF go?
 *  14 - fix base x (1)             => where 3DOF go?
 * */

#endif // _CONSTANTS_H_
