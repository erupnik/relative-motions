#ifndef _ALL_H
#define _ALL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>


#include "Eigen/Core"
#include "Eigen/Dense"

typedef Eigen::Matrix<double,18,18> Mat18d;
typedef Eigen::Matrix<double,12,12> Mat12d;
typedef Eigen::Matrix<double,6,6> Mat6d;
typedef Eigen::Matrix<double,3,3> Mat3d;
typedef Eigen::Matrix<double,18,1> Vec18d;
typedef Eigen::Matrix<double,12,1> Vec12d;
typedef Eigen::Matrix<double,6,1> Vec6d;
typedef Eigen::Matrix<double,3,1> Vec3d;
typedef Eigen::Matrix<double,2,1> Vec2d;

#include <constants.h>

#endif //_ALL_H
