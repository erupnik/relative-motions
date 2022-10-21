#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <all.h>
#include <thread>

/* Project-wide contants */

const double _ETA_up=0.1;
const double _ETA=0.001;

const bool _GAUGE_FIRST_CAM_FIX=false;
const bool _GAUGE_BASE_FIX=false;
const bool _GAUGE_RAPPEL_POSES=true;//only first pose and x of the second pose 


const std::string _DELIMITER_="-zyx-";

#endif // _CONSTANTS_H_
