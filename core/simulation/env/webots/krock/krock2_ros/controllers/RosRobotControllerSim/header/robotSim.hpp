#ifndef ROBOTSIM_HPP
#define ROBOTSIM_HPP


#define OPTIMIZATION
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include <sys/time.h>
#include "joystick.h"
#include "Utils/utils.hpp"
#include "Utils/misc_math.hpp"
#include <fstream>
#include <iostream>

//#include <webots/Servo.hpp>

#include "optimization_tools.hpp" 


#include <webots/Supervisor.hpp>
#include <webots/Robot.hpp>
#include <webots/Gyro.hpp>
#include <webots/Accelerometer.hpp>
#include <webots/Compass.hpp>
#include <webots/GPS.hpp>
#include <webots/Motor.hpp>
#include <webots/TouchSensor.hpp>
#include <webots/Node.hpp>
#include <webots/Field.hpp>
#include <webots/PositionSensor.hpp>



#define N_TOUCH_SENSORS    4
const int maxSpeed = 45;    // mm/s


extern int IS_SIMULATION;
extern int IS_OPTIMIZATION;
extern int USE_JOYSTICK;
extern int SWIM;
extern int SPINE_COMPENSATION_WITH_FEEDBACK;
extern int USE_REFLEXES;
extern int USE_TAIL;
extern int USE_IMU;
extern int USED_MOTORS[30];
extern int AUTO_RUN, AUTO_RUN_STATE;
extern int LOG_DATA;
extern int JOYSTICK_RECORDING_REPLAYING;
extern int Ntrunk, Nlegs, Ntail, Nneck, NUM_MOTORS;


using namespace std;
//using namespace Eigen;
//using namespace Robot;
//using Eigen::Matrix Matrix;



class RobotSim : public webots::Supervisor{


  public:
    //================== public variables ===============================================
    //webots::Servo *servo[NUM_MOTORS], *linservo[3];
    vector<webots::Motor*> rm_motor;
    vector<webots::PositionSensor*> ps_motor;
    webots::Gyro *gyro;
    webots::Accelerometer *acc;
    webots::Compass *compass;
    webots::GPS *gps_fgird, *gps_hgird;
    webots::Camera *camera;
    webots::TouchSensor *touch_sensor[N_TOUCH_SENSORS];
    webots::TouchSensor *touch_sensor_spine[12];
    webots::Node *fgirdle, *FL_marker, *FR_marker, *HL_marker, *HR_marker, *CoM_marker, *CoM_marker_proj, *roboDef, *tsdefFL, *tsdefFR, *tsdefHL, *tsdefHR;
    const double *compassData, *gpsDataFgird, *gpsDataHgird, *gyroData, *accData, *ts_fl, *ts_fr, *ts_hl, *ts_hr, *rotMat, *posFL, *posFR, *posHL, *posHR;
    webots::Field *roboRot, *roboPos, *CoM_marker_pos, *CoM_marker_pos_proj;
    double gamma, tcyc, tcyc_an;
    unsigned char camImg[230400];
    double t_total;

    webots::Node *supportPolyDEF;
    
    //================== public functions ===============================================
    RobotSim(int TIME_STEP); // constructor
    void setAngles(double*,int*);
    void torques(double*, int*);
    void getPositionTorques(double *d_posture, double *d_torques);
    void GetPosition(double *gpsData1, double *gpsData2);
    void InitIMU();
    void ReadIMUAttitude(double *rotmat);
    void ReadTouchSensors(double *ts_data);
    void GetCamera();
    void InitCamera();
    void killSimulation();
    void ColorMarkers(double *logic, double trans, double *col1);
    void setPositionRotation(double *p, double *r);
    void setPositionOfRobot(double *p);
    void setCoMmarker(double *p);
    void GetCompass(double *data_i);
    void GetFeetGPS(double *FL_feet_gpos, double *FR_feet_gpos, double *HL_feet_gpos, double *HR_feet_gpos);
    void setServoMaxForce(double *force);
    void GetTailForce(double *tailForce);
    void setSupportPoly(MatrixXd globalPoly, MatrixXd stance, double height);
  private:
    //================== private variables ===============================================




};


#endif
