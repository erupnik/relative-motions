#ifndef _ALL_H
#define _ALL_H

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <ctime>

#include <stdio.h>
#include <omp.h>

#include "Eigen/Core"
#include "Eigen/Dense"

#include <boost/program_options.hpp>

typedef Eigen::MatrixXd MatXd;
typedef Eigen::Matrix<double,18,18> Mat18d;
typedef Eigen::Matrix<double,12,12> Mat12d;
typedef Eigen::Matrix<double,11,11> Mat11d;
typedef Eigen::Matrix<double,6,6> Mat6d;
typedef Eigen::Matrix<double,5,5> Mat5d;
typedef Eigen::Matrix<double,3,3> Mat3d;
typedef Eigen::VectorXd VecXd;
typedef Eigen::Matrix<double,18,1> Vec18d;
typedef Eigen::Matrix<double,12,1> Vec12d;
typedef Eigen::Matrix<double,11,1> Vec11d;
typedef Eigen::Matrix<double,6,1> Vec6d;
typedef Eigen::Matrix<double,6,1> Vec5d;
typedef Eigen::Matrix<double,3,1> Vec3d;
typedef Eigen::Matrix<double,2,1> Vec2d;

#include <constants.h>

#include "glog/logging.h"

#endif //_ALL_H
