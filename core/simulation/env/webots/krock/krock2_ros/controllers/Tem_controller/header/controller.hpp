#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP


#define OPTIMIZATION2
//#define EIGEN_DONT_ALIGN_STATICALLY
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "eigen3/Eigen/SVD"
#include "joystick.h"
#include "Utils/utils.hpp"
#include "Utils/misc_math.hpp"
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string.h>
#include <thread>  
#include <mutex> 
#include "header/structDefs.hpp"


/////////////////////////////////
#include "header/optimization_tools.hpp"


#define my_pi          3.141592653589793

#define N_anTr      1000
#define JOYSTICK_DEVNAME "/dev/input/js0"
/*
#define USE_JOYSTICK   1
#define JOYSTICK_TYPE 1     //PS1 - 1, PS3 - 2, Logitech - 3
#define SWIM   0
#define IS_SIMULATION 1
#define SPINE_COMPENSATION_WITH_FEEDBACK 0
*/

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



extern std::mutex mpc_in_mtx, mpc_out_mtx, mpc_sync_mtx;
extern std::mutex fdo_in_mtx, fdo_out_mtx, fdo_sync_mtx;

using namespace Eigen;

enum{WALKING, CRAWLING, STANDING, POSING, SWIMMING, ANIMAL_WALKING, ANIMAL_SWIMMING, ANIMAL_AQUASTEP, INITIAL, NONPERIODIC};


class Controller{


  public:
    //================== public variables ===============================================
    Joystick js;
    int state;

    double turning_curvature, crawling_turning_curvature, freq_swim, freq_walk, freq_aquastep;
    Matrix<double, Dynamic, 1> fbck_torques, fbck_angles, fbck_current, joint_angles;
    double q_neck;
    Vector2d q_tail;
    Matrix<double,12, 1> force, force_filt;
    Matrix<double,3, 1> posture;
    Matrix<double, 4, 1> legs_stance, legs_contact, legs_contact_control, legs_out_of_reach;
    double compassData[3], gpsData[3], gyroData[3], accData[3];
    Matrix<double, 3, Dynamic> global_joint_pos;
    Matrix<double, 3, 4> global_feet_pos;
    Matrix<double, 3, 4> feet_force, feet_force_reference, feet_force_filt, fc_feet_force, feet_force_from_current, nnEstForce;
    Matrix<double, 4, 4> rawTCompOptoforce;
    Matrix<double, 12, 1> rawForce;
    vector<Matrix4d> HJs, HJs_g, HJt, HJt_g;
    vector<Matrix4d> HJfl, HJfr, HJhl, HJhr, HJfl_g, HJfr_g, HJhl_g, HJhr_g;
    vector<MatrixXd> legJacob;
    Matrix<double, 4, 4> Fgird, Hgird, Fgird0, Hgird0, Fgird_mean, Hgird_mean, Hhead;
    double hgirdlePhi;
    Vector4d phaseTransitionStanceSwing, phaseTransitionSwingStance, legs_stance_old;
    double IG;


    // needed for RUN
    Matrix<double, 11, 1> amps0, amps, phases0, phases, offsets0, offsets, phases_swim, offsets_swim, amps_swim, qs_posing, qsr;
    Matrix<double, Dynamic, 1> trunk_kin, neck_kin, tail_kin;
    Matrix<double, 5, 3> FL_kin, FR_kin, HL_kin, HR_kin;
    Matrix<double, 4, 1> stancePhase, swingPhase, legPhase;

    // gait parameters (currently in use), walking parameters, crawling parameters
    gaitParams GP, WP, CP, WPlow;

    //joystick manipulation
    double stick_r_lim, joy_walk_max_freq, joy_walk_min_freq, joy_walk_speed_change_filt, disable_crab_walk_lim;
    double posing_xFH, posing_yF, posing_zF, posing_head;
    double crab_rotAng, ellipse_small_axis;
    double posing_joy_x1_rate, posing_joy_y1_rate, posing_joy_y2_rate, posing_head_rate, posing_head_limit;
    double joy_swim_max_offset, joy_swim_max_freq, joy_swim_speed_change_filt, joy_swim_turning_dead_zone, aquastep_pitch_offset;

    //cpg
    double wbck, wfwd, wlat, ph_lagbck, ph_lagfwd0, ph_lagfwd, ph_laglat, nu_cpg, a_cpg, cpg_offset, phase_lag_change_rate, phase_lag_min, phase_lag_max, R_cpg_gain;

    

    //joystick
    double joy_x1, joy_x2, joy_y1, joy_y2, joy_x3, joy_y3;
    int joy_l1, joy_l2, joy_l3, joy_r1, joy_r2, joy_r3, joy_sel, joy_start, joy_bD, joy_bL, joy_bR, joy_bU, joy_aD, joy_aL, joy_aR, joy_aU;
    double joy_lsr, joy_rsr;
    double joy_lsphi, joy_rsphi;



    //reflexes 
    Vector2d reflON;
    double stumbleForceLimit, stumbleTimeout, stumbleFilterConstant;
    Vector3d stumbleImpulse;
    Vector2d stumblePhaseLimits;

    double extendForceLimit, extendTimeout, extendFilterConstant;
    Vector3d extendImpulse;
    Vector2d extendPhaseLimits;

    // old
    int stu_reflex_active[4], ext_reflex_active[4];
    double extRefForceLim;
    double extRefOnFilter;
    double extRefOffFilter;
    double extRefSpike;
    double extRefTimeout, stuRefTimeout;
    double stuRefForceLimX, stuRefForceLimZ;
    double stuRefOnFilter;
    double stuRefOffFilter;
    double stuRefDx;
    double stuRefDz;




    Matrix<double, 3, 3> fbck_fgirdRotMat, fgird2ForRotMat, forRotMat, compassRotMat;
    //Matrix<double, 4, 4> Fgird2RobotHmat, Hgird2RobotHmat;

    double attData[3];
    Vector3d globalPosFL, globalPosFR, globalPosHL, globalPosHR;
    Vector3d gpsPos;
    Matrix3d rotMatFgird;

    Vector2d FFParam;
    Vector3d forRPY;
    // girdle trajectories
    Matrix<double, 6, 2> girdleTraj;
    Matrix<double, 3, 2> girdleVelocity, girdleVelocity_filtered;
    Matrix<double, 3, 2> forVelocity, forVelocity_filtered;
    Matrix<double, 1, 2> girdleAngularVelocity, girdleAngularVelocity_filtered;
    Matrix<double, 1, 2> forAngularVelocity, forAngularVelocity_filtered;
    Matrix<double, 6, 2> forTraj, forTraj_posing0;
    Matrix<double, 3, 300> forTrajHist;
    std::vector<Matrix<double, 3, 4>> feetTrajHist;
    Matrix<double, Dynamic, 1> q0_trunk_from_spline;
    Matrix<double, 2, 2> girdleRealOrientation;
    Matrix<double, 2, 4> legStanceTraj;
    double walkingDirection, bodyAngularVelocity;
    double fgird_heading_angle;
    double t, dt;
    double walking_forward_velocity, walking_angular_velocity, walking_rotational_velocity;
    Vector3d realCoM;

    // walking done properly
    Matrix<double, 3, 4> ikinRef, feetReference, feetReference_reflex, feetReference_forceCorrections, feetFeedbackPosition, feetMod;
    Matrix<double, 3, 4> feetLocations, feetLocationsOdo, feetVelocities, feetLocationsLocal;
    double posingFeetHeight;
    
    Vector2d girdleCpgOutput;
    double *joysticRecordings;
    int sizeOfJoystickRecording;

    // pipe crawling

    enum CRSTATE{initial, selectLegs, moveBody, swingLeg, loadingLeg, unloadingLeg};
    CRSTATE crstate;



    //================== public functions ===============================================
    Controller(double time_step); // constructor
    void setTimeStep(double time_step);
    bool runStep();
    void readJoystick();
    bool updateState();
    void getAngles(double *table);
    void getTorques(double *table);
    void forwardKinematics();
    std::vector<Matrix<double, 4, 4>> legKinematics(Vector4d q, int leg);
    Matrix<double, 3, 4> Jacob(Vector4d q, int leg);
    MatrixXd forceEstimation();
    Matrix<double, 3, 4> contactForceController(Matrix<double, 3, 4> feetForceRef);
    Matrix<double, 3, 4> contactNormalForceController(Matrix<double, 3, 4> feetForceRef, Matrix<double, 3, 4> feetForceEst, Matrix<double, 3, 4> nVec);
    void updateRobotState(double *d_posture, double *d_torque);
    void getForce(double *force);
    void getSensors(int acc, double *accData_i, int gyro, double *gyroData_i, int compass, double *compassData_i, int gps, double *gpsData_i);
    void getAttitude(double *rotmat);
    void getRPY(double *rpy);
    void globalKinematics(double gpsPosition[3], double rotMat[9]);
    void getGlobalCoM(double pos[3]);
    Matrix<double, 3, 4>  getGlobalSupportPoly();
    void getGPS(double *gps, double *gps_feet1, double *gps_feet2, double *gps_feet3, double *gps_feet4);

    #ifdef OPTIMIZATION
        void optimizationInit();
        void optimizationEnd();
        int optimizationShouldWeStopIt(double timestep);
        std::vector<double> params;
        std::vector<double> settings;
        std::vector<std::string> params_names;
        std::vector<std::string> settings_names;
        double rolling, pitching, yawing, applied_torques, applied_torques_l, applied_torques_s;
    #endif
    
    // playing with torque
    bool DirectTorqueSetup();

    // estimation
    void getAcceleration(double acc[3]);
    void getFeetPosition(double *fl, double *fr, double *hl, double *hr);
    MatrixXd getReferenceFootsteps();
    void GetCompass(double data[3]);
    Vector3d getCoM();
    void torquePID();


    // static walking
    Vector3d getSwingTrajectory(Vector3d initPoint, Vector3d middlePoint, Vector3d finalPoint, double phase, int leg);
    bool moveSwingLeg(int leg);
    bool moveSwingLegManualPhase(int leg, double swing_phase);
    void moveAlongVector(int leg, Vector3d nVec, double dist);
    //void moveBody(Vector3d bodyVelocity, double headingAngularVelocity);
    void moveGirdle(Vector3d girdleVelocity, double girdleAngularVelocity, int girdleNum);
    void walkingStateMachine();
    void crawlingStateMachine();
    Vector3d trunkForwardKinematics(Vector3d fgird_pos, MatrixXd q_trunk);
    void trunkInverseKinematics();
    void getWalkingSpeedFrequency();

    // MPC 
    //std::vector<Matrix<double,3,4>> predictedFootsteps;
    std::vector<Matrix<double,3,4>> supportPolys, supportPolysIMU;
    std::thread mpc_thread;
    Matrix<double,3,4> followCoMReference(MatrixXd comref);

    //FDO
    std::thread fdo_thread;
    FDO fdo;

    MatrixXd calcForceDistribution(MatrixXd NS, MatrixXd maxForce, MatrixXd ddx, MatrixXd dw, MatrixXd fg);

    void sendStatusToMatlab(int downsampling, const char *IPadr);



  private:
    //================== private variables ===============================================
    //vector<qpOASES::QProblemB> ikinQpSolver; 
    vector<qpOASES::QProblem> ikinQpSolver; 
    int useAnDF, useAnDH, useAnSP, ikin_maxIter;
    double ikin_tol;
    double T_trans, T_trans0, T_stand;
    Vector3d ikin_constr_penalty;
    // joystick objects
    js_event event;
    bool speed_or_frequency;
    double maxSpeed, maxFrequency;




    Matrix<double, 3, 4> stanceStart, stanceEstEnd;

    // constraints on joint angles
    Matrix<double, 2, 4> constrFL, constrHL, constrFR, constrHR; 
    double constrS;

    Matrix<double, 3, 4> workspaceEllipsoid, workspaceEllipsoidCenter;


    Matrix<double, 3, 1> rFL, rFR, rHL, rHR, pFL, pFR, pHL, pHR, tmp3, rFL_posing, rFR_posing, rHL_posing, rHR_posing;
    Matrix<double, 4, 1> reCon, rsCon, vmCon;
    double forceFilterExt, forceFilterStu;
    double z_offset;
    vector<double> angSignsCorrIkin2Webots, angShiftsCorrIkin2Webots;
    vector<double> angSignsCorrIkin2Robot, angShiftsCorrIkin2Robot;


    Matrix<double, Dynamic, 1> angles, torques;


    // initial conditions
    Vector4d q0FL, q0FR, q0HL, q0HR, qFL, qFR, qHL, qHR;
    Matrix<double, Dynamic, 1> init_angles, q_trunk;


    // auxiliary variables
    Matrix<double, 11, 1> qs;
    Matrix<double, 4,4> CinvF, CinvH, MF, MH;
    double max_dist;
    Vector2d lamF, lamH;





    // dynamics
    double Tf1; // trajectory filtering
    double Tfilt_angles;
    bool is_flipped;




    // trajectories
    //double nfl, nfr, nhl, nhr, rFLy, rFRy, rHLy, rHRy, fxfl, fxfr, fxhl, fxhr;
    

    double xGain, yGain, Troll_posture, Tpitch_posture;
    

    Matrix<double, 11, 1> spine_gains0;
    Matrix<double, 4, 1> legs_height0;
    double erRefOnFilt, erRefOffFilt, er_timeout, er_duration, er_flow_limit;
    double er_slipping_Tfilt, er_slipping_threshold;


    //swimming
    Matrix<double, 16, 1> swim_legs_pos;

    Matrix<double, 3, 4> feetGPS;

    Matrix<double, 4, 23> masses;
    Matrix<double, 3, 1> CoM;
    


    // force sensors and friction cones
    Matrix<double, 4, 1> feet_friction_angles, feet_is_slipping;
    
    // MPC
    MPC mpc_com;
    Matrix<double, 2, 1> mpcComRef;
    //Matrix<double, Dynamic, Dynamic> predLegPhase;

    // FDO and crawling
    double mu_wall, mu_ground, mu_conservative_scaling, forceTreshold, maxContactForce, fdo_dt, rawForcesFilter, leg_loading_filter, leg_approach_speed;
    double contactF_pid_Kp, contactF_pid_Ki, contactF_pid_Kd, contactF_pid_Td, contactF_pid_limit;
    Matrix<double,Dynamic,4> crawlingSequence;
    bool crawling_continue;
    //================== private functions ===============================================

    Vector4d iKinNullIDLS(int leg, Vector3d pref, Vector4d qref, Vector4d q0, Vector2d lam, Matrix4d M, 
                            double max_dist, int maxIter, double tol, MatrixXd constr, Vector3d constr_penalty);
    Vector4d iKinQpOases(int leg, Vector3d pref, Vector4d qref, Vector4d q0, Vector2d lam, Matrix4d M, 
                            double max_dist, int maxIter, double tol, MatrixXd constr, Vector3d constr_penalty);
    void joystickRecord();
    void joystickReplay();
    void joystickAuto();
    bool getParameters();
    void walkingTrajectories();
    void swimFun();
    void contactDetection();
    void postureControl();
    void legExtension();
    void stumbleReflex();
    MatrixXd transformation_Ikin_Webots(MatrixXd joint_angles_in, int direction, int shifts);
    MatrixXd transformation_Ikin_Robot(MatrixXd joint_angles_in, int direction, int shifts);
    MatrixXd flipKinematics(MatrixXd joint_angles_in);
    void swimFunCPG();
    MatrixXd inverseKinematicsController(MatrixXd pRef);
    void slippingFromFrictionCone();
    void legPhaseDynamics(double dt);
    void girdleTrajectories(double v, double w, Vector3d V, double W, double fgirdArc, Vector2d spineCPGscaling);
    void girdleVelocities();
    void legSwingTrajectory();
    void legStanceTrajectory();
    void getLegJacobians();
    std::vector<Matrix<double,3,4>> predictTrajectories(int N, double time_step, MatrixXd *predLegPhase);
    void runComMPC();
    void girdleOscillations(Vector2d spineCPGscaling);
    void getSwingStanceEndPoints();
    MatrixXd standingTransition(MatrixXd pRef);
    void posingManipulation();
    void pipePathPlanner();
    void getForceFromCurrent();
    void mlpContactForceEstimationDummy();

    


    
};




#endif
