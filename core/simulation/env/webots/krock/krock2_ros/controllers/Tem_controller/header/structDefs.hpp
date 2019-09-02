#ifndef STRUCTDEFS_HPP
#define STRUCTDEFS_HPP

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include <stdlib.h>


struct gaitParams {
	Eigen::Matrix<double, 2, 1> spineCPGscaling;
  	Eigen::Matrix<double, 2, 1> Duty;
    Eigen::Matrix<double, 4, 1> phShifts;
	Eigen::Matrix<double, 3, 4> midStance;
	Eigen::Matrix<double, 1, 4> ellipse_a, ellipse_b;
	Eigen::Matrix<double, 1, 4> swing_height, swing_width;
	
    
    Eigen::Matrix<double, 3, 4> nSurf;
    Eigen::Matrix<double, 3, 4> nLO, nTD;

	Eigen::Matrix<double, 3, 2>  bezierScaling;
	double tSclSwing;
 	
 	Eigen::Matrix<double, 4, 4> qNULL;


};






#endif
