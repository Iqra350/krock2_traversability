#include "controller.hpp"
#include <iostream>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>



using namespace std;
using namespace Eigen;



// ---------------------------------------------- FUNCTIONS ----------------------------------------------------------
/* Controller constructor */
Controller :: Controller(double time_step)
{

    //############################# RESIZE DYNAMIC MATRICES ##########################################
    fbck_torques.resize(NUM_MOTORS,1);
    fbck_angles.resize(NUM_MOTORS,1);
    joint_angles.resize(NUM_MOTORS,1);
    fbck_current.resize(NUM_MOTORS,1);

    global_joint_pos.resize(3, NUM_MOTORS);

    angles.resize(NUM_MOTORS,1);
    torques.resize(NUM_MOTORS,1);


    init_angles.resize(NUM_MOTORS,1);


    


    q0_trunk_from_spline.resize(Ntrunk, 1);
    q_trunk.resize(Ntrunk, 1);

    trunk_kin.resize(Ntrunk+1, 1);

    if(Nneck>0){
        neck_kin.resize(Nneck, 1);
    }
    
    if(Ntail>0){
        tail_kin.resize(Ntail, 1);
        HJt.resize(Ntail+1);
        HJt_g.resize(Ntail+1);
    }   
    

    HJs.resize(Ntrunk+1);
    HJfl.resize(5);
    HJfr.resize(5);
    HJhl.resize(5);
    HJhr.resize(5);
    HJs_g.resize(Ntrunk+1);
    HJfl_g.resize(5);
    HJfr_g.resize(5);
    HJhl_g.resize(5);
    HJhr_g.resize(5);
    

    ikinQpSolver.resize(4);
    legJacob.resize(4);


    
    //############################ initialize MPC controller ##########################################
    setTimeStep(time_step);
    mpc_com.readParameters("config/MPC/mpc.config");
    if(mpc_com.enabled){
        cout << "setConstraintMatrices"<<endl;
        mpc_com.setConstraintMatrices();
        cout << "constructQP"<<endl;
        mpc_com.constructQP();
        cout << "printMatrices"<<endl;
        mpc_com.printMatrices("MPC_TEST.txt");
        cout << "initQPSolver"<<endl;
        mpc_com.initQPSolver();
    }
    angSignsCorrIkin2Webots.reserve(NUM_MOTORS);
    angShiftsCorrIkin2Webots.reserve(NUM_MOTORS);
    angSignsCorrIkin2Robot.reserve(NUM_MOTORS);
    angShiftsCorrIkin2Robot.reserve(NUM_MOTORS);

    


    // other variables
    supportPolys.resize(mpc_com.N);
    supportPolysIMU.resize(mpc_com.N);


    //############################ PARAMETERS ########################################################
    getParameters();
    GP=WP;
    cout<<"I successfully read the parameters"<<endl;
    state=INITIAL;
    T_trans=T_trans0;
    T_stand=T_trans0/2;


    t=0;
    if(USE_JOYSTICK){
        js.load(JOYSTICK_DEVNAME);

        if(js.is_ready()) {
            printf("Joystick Ready\n");
            updateState();
        }
    }

    forRPY.setZero();
    is_flipped=false;
    turning_curvature=0;
    //########## STANDING POSITION #################

    force_filt << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    qs << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    // initial conditions
    qFL.setZero();
    qFR.setZero();
    qHL.setZero();
    qHR.setZero();
    init_angles.setZero();
    init_angles*=0;

    legs_out_of_reach.setZero();




    feetReference = GP.midStance;



    legs_stance << 1,1,1,1;
    swingPhase << 0, 0, 0, 0;





    joint_angles.setZero();
    feet_force<<0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    feet_force_filt=feet_force;
    fbck_angles.setZero();
    freq_walk=0;

    rotMatFgird=Matrix3d::Identity();

    girdleVelocity_filtered=MatrixXd::Zero(3,2);
    girdleAngularVelocity_filtered=MatrixXd::Zero(1,2);
    forVelocity_filtered=MatrixXd::Zero(3,2);
    forAngularVelocity_filtered=MatrixXd::Zero(1,2);

    feetReference_forceCorrections.setZero();

    // initialize trajectory container
    for(int i=0; i<300; i++){
        forTrajHist(0,i)=-i/250.*IG;
        forTrajHist(1,i)=0;
        forTrajHist(2,i)=0;
    }
    crawling_turning_curvature=0;
    crawling_continue = false;
    

    //############################ initialize FDO controller ##########################################
    fdo.initSolver(masses.sum(), mu_wall, fdo_dt);
    legs_contact_control << 0,0,0,0;

    crstate = initial;
    //#############################Joystick dynamics########################################


/*
    nfl=0; nfr=0; nhl=0; nhr=0;
    fxfl=0; fxfr=0; fxhl=0; fxhr=0;*/

    posing_xFH=0; posing_yF=0; posing_zF=0; posing_head=0;
    crab_rotAng=0;

    //############################# REFLEXES ########################################
    for(int i=0; i<4; i++){
        stu_reflex_active[i]=0;
        ext_reflex_active[i]=0;
    }



    fgird2ForRotMat=MatrixXd::Identity(3,3);


    //############################# LEG PHASES ########################################
    for(int i=0; i<4; i++){
        legPhase(i)=GP.phShifts(i)<0?GP.phShifts(i)+1:GP.phShifts(i);
    }
    legs_contact.setZero();
    
    girdleTraj << 0,0,0,0,0,0, 0,0,0,0,0,0;
    girdleTraj(3,0)=-IG; girdleTraj(3,1)=-IG;
    //girdleRealOrientation << 0,0, 0,0;
    forTraj << 0,0,0, 0,0,0, 0,0,0, 0,0,0;
    forTraj(3,0)=-IG; forTraj(3,1)=-IG;
    fgird_heading_angle=0;
}

/* Sets time step */
void
Controller :: setTimeStep(double time_step)
{
    dt=time_step;
    return;
}

/* Updates angles for all joints - MAIN FUNCTION */
bool
Controller :: runStep()
{

    if(!updateState()){
        return false;
    }
    
    if(state==INITIAL){
        angles=init_angles;
        return true;
    }


    joint_angles=angles;
    //joint_angles=fbck_angles;
    forwardKinematics();
    realCoM=getCoM();

    mlpContactForceEstimationDummy();
    //############################ TIME #########################################################
    if(state==WALKING || state==CRAWLING  || state==SWIMMING){
        t+=dt;
    }

    //state=WALKING;
    //GP=CP;


    //--------------------------- joystick input-----------------------------
    if(state==WALKING || state==CRAWLING){
        walkingDirection=joy_lsphi;
        getWalkingSpeedFrequency();
        double r_elly=GP.ellipse_a(0)*GP.ellipse_b(0)/sqrt(GP.ellipse_a(0)*GP.ellipse_a(0)*sin(walkingDirection)*sin(walkingDirection)+
                                                  GP.ellipse_b(0)*GP.ellipse_b(0)*cos(walkingDirection)*cos(walkingDirection));
        walking_forward_velocity*=r_elly/GP.ellipse_a(0);


        if(joy_r1==0){
            turning_curvature=pt1((joy_x2), turning_curvature, 0.2, dt);
            walking_angular_velocity = (-turning_curvature)*walking_forward_velocity*2;
            walking_rotational_velocity=pt1(0, walking_rotational_velocity, 0.2, dt);
        }
        else{
            walking_rotational_velocity=pt1(joy_x2*0.3, walking_rotational_velocity, 0.2, dt);
        }
        
    }
    //-----------------------------------------------------------------------
    getSwingStanceEndPoints();


    /*if(state==WALKING && t>0){
        state=CRAWLING;
        crawling_continue=true;
        GP=CP;
    }*/
    if(state==WALKING){
        walkingStateMachine();
    }
    if(state==CRAWLING){
        crawlingStateMachine();
    }
   
    if(state==POSING){
        posingManipulation();
    }

    // reflexes
    feetReference_reflex = MatrixXd::Zero(3,4);
    stumbleReflex();
    legExtension();

    // combine all the factors for feet reference
    for(int i=0; i<4; i++){
        ikinRef.block<3,1>(0,i)=feetReference.block<3,1>(0,i) + 
                                1*feetReference_forceCorrections.block<3,1>(0,i)+
                                1*feetReference_reflex.block<3,1>(0,i);
    }
    
    MatrixXd tempLegAngles(4,4);
    ikinRef = standingTransition(ikinRef);

    tempLegAngles = inverseKinematicsController(ikinRef); 

    angles.block<4,1>(0,0) =tempLegAngles.block<4,1>(0,0);
    angles.block<4,1>(4,0) =tempLegAngles.block<4,1>(0,1);
    angles.block<4,1>(8,0) =tempLegAngles.block<4,1>(0,2);
    angles.block<4,1>(12,0)=tempLegAngles.block<4,1>(0,3);
    angles.block<2,1>(16,0)=q_trunk; 

    if(Nneck>0){
       angles(18)=q_neck; 
    }
    if(Ntail>0){
        angles.block<2,1>(19,0)=q_tail;
    }
    
    joint_angles=angles;


 




    if(LOG_DATA){
        static ofstream timeLOG("/data/thorvat/webotsLogs/timeLOG.txt");
        timeLOG << t << endl;
        
        static ofstream crstateLOG("/data/thorvat/webotsLogs/crstateLOG.txt");
        crstateLOG << crstate << endl;

        static ofstream fkinLog("/data/thorvat/webotsLogs/fkinLog.txt");
        for(int i=0;i<5;i++){
            for(int j=0;j<3;j++)
                fkinLog<<HJfl_g[i](j, 3)<<"\t";
        }
        for(int i=0;i<5;i++){
            for(int j=0;j<3;j++)
                fkinLog<<HJfr_g[i](j, 3)<<"\t";
        }
        for(int i=0;i<5;i++){
            for(int j=0;j<3;j++)
                fkinLog<<HJhl_g[i](j, 3)<<"\t";
        }
        for(int i=0;i<5;i++){
            for(int j=0;j<3;j++)
                fkinLog<<HJhr_g[i](j, 3)<<"\t";
        }
        for(int j=0;j<3;j++)
                fkinLog<<Fgird(j, 3)<<"\t";
        for(int i=0;i<Ntrunk;i++)
            for(int j=0;j<3;j++)
                fkinLog<<HJs_g[i](j, 3)<<"\t";
        for(int j=0;j<3;j++)
            fkinLog<<Hgird(j, 3)<<"\t";
        fkinLog<<endl;

        static ofstream forceLOG("/data/thorvat/webotsLogs/forceLOG.txt");
        forceLOG << feet_force.block<3,1>(0,0).transpose() << "\t";
        forceLOG << feet_force.block<3,1>(0,1).transpose() << "\t";
        forceLOG << feet_force.block<3,1>(0,2).transpose() << "\t";
        forceLOG << feet_force.block<3,1>(0,3).transpose() << "\t";
        forceLOG << feet_force_reference.block<3,1>(0,0).transpose() << "\t";
        forceLOG << feet_force_reference.block<3,1>(0,1).transpose() << "\t";
        forceLOG << feet_force_reference.block<3,1>(0,2).transpose() << "\t";
        forceLOG << feet_force_reference.block<3,1>(0,3).transpose() << endl;

        static ofstream imuLOG("/data/thorvat/webotsLogs/imuLOG.txt");
        imuLOG << forRPY.transpose() << endl;

        static ofstream positionLog("/data/thorvat/webotsLogs/positionLog.txt");
        static ofstream torqueLog("/data/thorvat/webotsLogs/torqueLog.txt");

        positionLog << fbck_angles.transpose() << endl;
        torqueLog << fbck_torques.transpose() << endl;

        static ofstream feetReferenceLog("/data/thorvat/webotsLogs/feetReferenceLog.txt");
        for(int i=0; i<4; i++){
            feetReferenceLog << ikinRef.block<3,1>(0,i).transpose() << "\t";
        }
        for(int i=0; i<4; i++){
            feetReferenceLog << feetLocationsLocal.block<3,1>(0,i).transpose() << "\t";
        }
        for(int i=0; i<4; i++){
            feetReferenceLog << feetReference.block<3,1>(0,i).transpose() << "\t";
        }
        for(int i=0; i<4; i++){
            feetReferenceLog << feetReference_forceCorrections.block<3,1>(0,i).transpose() << "\t";
        }
        feetReferenceLog << endl;

    }   



    return true;
}



/* Writes angles calculated by runStep function into a table - interface with Pleurobot class */
void
Controller :: getAngles(double *table)
{

        static MatrixXd angles_old=init_angles;
        static MatrixXd angles2=transformation_Ikin_Robot(init_angles,1,1);
        static MatrixXd angles2_old=transformation_Ikin_Robot(init_angles,1,1);


        is_flipped = abs(forRPY(0))>my_pi/2;
        static bool is_flipped_old = false;
        if(is_flipped!=is_flipped_old){
            T_trans=T_trans0;
        }
        is_flipped_old=is_flipped;

        if(is_flipped){
            angles=flipKinematics(angles);
        }

        angles=pt1_vec(angles, angles_old, T_trans, dt);
        angles_old=angles;


        if(IS_SIMULATION){
            angles2=transformation_Ikin_Webots(angles,1, 1);
        }

        if(!IS_SIMULATION){
            angles2=transformation_Ikin_Robot(angles,1,1);
        }


        angles2=pt1_vec(angles2, angles2_old, Tfilt_angles, dt);
        angles2_old=angles2;

        for(int i=0; i<NUM_MOTORS; i++){
            table[i]=angles2(i);
        }
        static ofstream anglesLOG("/data/thorvat/webotsLogs/anglesLOG.txt");
        anglesLOG << angles2.transpose() << endl;

    return;
}

/* Writes torques calculated by runStep function into a table - interface with Pleurobot class */

void
Controller :: getTorques(double *table)
{
        static MatrixXd torques2(NUM_MOTORS,1);
        //static MatrixXd torques(27,1);

        static MatrixXd torques2_old=MatrixXd::Zero(NUM_MOTORS, 1);
        static MatrixXd torques_old=MatrixXd::Zero(NUM_MOTORS, 1);




        if(IS_SIMULATION){
            torques2=transformation_Ikin_Webots(torques, 1, 0);
        }

        if(!IS_SIMULATION){
            torques2=transformation_Ikin_Robot(torques, 1, 0);
        }


        for(int i=0; i<NUM_MOTORS; i++){
            table[i]=torques2(i);
        }


        // CONSTRAINTS
    if(!IS_SIMULATION){
        for(int i=0; i<NUM_MOTORS; i++){
            table[i]=table[i]>7?7:table[i];
            table[i]=table[i]<-7?-7:table[i];
        }
    }

    return;
}


/* Estimates foot contact for stance and swing (obstacle) phases and for different reflexes*/
void
Controller :: contactDetection()
{
    double limF, limH, forceFilter=0.1;
    Vector3d tmpForce;
    reCon << 0, 0, 0, 0;
    rsCon << 0, 0, 0, 0;
    vmCon << 0, 0, 0, 0;

    // leg extension reflex
    static MatrixXd force_filt1=MatrixXd::Zero(12, 1);
    limF=10;
    limH=10;
    force_filt1=pt1_vec(force, force_filt1, forceFilter, dt);

    if(force_filt1.block<3,1>(0,0).norm()>limF)
        reCon(0)=1;
    if(force_filt1.block<3,1>(3,0).norm()>limF)
        reCon(1)=1;
    if(force_filt1.block<3,1>(6,0).norm()>limH)
        reCon(2)=1;
    if(force_filt1.block<3,1>(9,0).norm()>limH)
        reCon(3)=1;


    // stumble reflex
    static MatrixXd force_filt2=MatrixXd::Zero(12, 1);
    limF=30;
    limH=50;
    force_filt2=pt1_vec(force, force_filt2, forceFilter, dt);

    if(force_filt2.block<3,1>(0,0).norm()>limF)
        rsCon(0)=1;
    if(force_filt2.block<3,1>(3,0).norm()>limF)
        rsCon(1)=1;
    if(force_filt2.block<3,1>(6,0).norm()>limH)
        rsCon(2)=1;
    if(force_filt2.block<3,1>(9,0).norm()>limH)
        rsCon(3)=1;





    // posture control
    //forceFilter=0.15;
    static double ffl3=0, ffr3=0, fhl3=0, fhr3=0;

    ffl3=pt1(force.block<3,1>(0,0).norm(), ffl3, forceFilter, dt);
    ffr3=pt1(force.block<3,1>(3,0).norm(), ffr3, forceFilter, dt);
    fhl3=pt1(force.block<3,1>(6,0).norm(), fhl3, forceFilter, dt);
    fhr3=pt1(force.block<3,1>(9,0).norm(), fhr3, forceFilter, dt);
}


/* Update robot positions and fbck_torques */
void
Controller :: updateRobotState(double *d_posture, double *d_current)
{
    for(int i=0; i<NUM_MOTORS; i++){
        fbck_angles(i)=d_posture[i];
        fbck_current(i)=d_current[i];
    }

    if(IS_SIMULATION){
        fbck_angles=transformation_Ikin_Webots(fbck_angles, -1, 1);
        fbck_torques=transformation_Ikin_Webots(fbck_current, -1, 0);
    }
    else{
        fbck_angles=transformation_Ikin_Robot(fbck_angles, -1, 1);
        fbck_torques=transformation_Ikin_Robot(fbck_current, -1, 0);
    }

}

/* Get sensor data */
void
Controller :: getSensors(int acc, double *accData_i, int gyro, double *gyroData_i, int compass, double *compassData_i, int gps, double *gpsData_i)
{

    if(acc){
        accData[0]=accData_i[0];
        accData[1]=accData_i[1];
        accData[2]=accData_i[2];
    }
    if(gyro){
        gyroData[0]=gyroData_i[0];
        gyroData[1]=gyroData_i[1];
        gyroData[2]=gyroData_i[2];
    }
    if(compass){
        compassData[0]=compassData_i[0];
        compassData[1]=compassData_i[1];
        compassData[2]=compassData_i[2];
    }
    if(gps){
        gpsData[0]=gpsData_i[0];
        gpsData[1]=gpsData_i[1];
        gpsData[2]=gpsData_i[2];
    }

    #ifdef OPTIMIZATION
    rolling+=gyroData[0]*gyroData[0]*dt;
    pitching+=gyroData[2]*gyroData[2]*dt;
    yawing+=gyroData[1]*gyroData[1]*dt;
    if(t>0.7){
    for(int i=0; i<27; i++){
        applied_torques+=fbck_torques(i)*fbck_torques(i)*dt;
        if(i<11){
            applied_torques_s+=fbck_torques(i)*fbck_torques(i)*dt;
        }
        else{
            applied_torques_l+=fbck_torques(i)*fbck_torques(i)*dt;
        }
    }
    }
    #endif
}


/* Get acc data */
void
Controller :: getAcceleration(double acc[3])
{
    accData[0]=acc[0];
    accData[1]=acc[1];
    accData[2]=acc[2];
}


/* Get gps */
void
Controller :: getGPS(double *gps, double *gps_feet1, double *gps_feet2, double *gps_feet3, double *gps_feet4){
    for(int i=0; i<3; i++){
        feetGPS(i,0)=gps_feet1[i];
        feetGPS(i,1)=gps_feet2[i];
        feetGPS(i,2)=gps_feet3[i];
        feetGPS(i,3)=gps_feet4[i];
        gpsPos(i)=gps[i];
    }

}


void
Controller :: getGlobalCoM(double pos[3]){
    Vector3d tmp;

    tmp=forRotMat*realCoM;

    pos[0]=tmp(0)+gpsPos(0);
    pos[1]=tmp(2)+gpsPos(1);
    pos[2]=-tmp(1)+gpsPos(2);


}


Matrix<double, 3, 4>
Controller :: getGlobalSupportPoly(){
    Vector3d tmp;
    MatrixXd globalPoly(3,4);
    Vector4d tmp_vector;

    for(int i=0; i<4; i++){
        tmp_vector << feetMod.block<3,1>(0,i), 1;
        if(i<2){
            //tmp_vector.block<3,1>(0,0)=AngleAxisd(forTraj(2,1)-girdleTraj(2,1), Vector3d::UnitZ())*tmp_vector.block<3,1>(0,0); // girdle frame to girdle
            tmp_vector.block<3,1>(0,0)=Fgird.block<3,3>(0,0)*tmp_vector.block<3,1>(0,0);
        }
        else{
            //tmp_vector.block<3,1>(0,0)=AngleAxisd(forTraj(5,1)-girdleTraj(5,1), Vector3d::UnitZ())*tmp_vector.block<3,1>(0,0); // girdle frame to girdle
            tmp_vector.block<3,1>(0,0)=Hgird.block<3,3>(0,0)*tmp_vector.block<3,1>(0,0);
        }

        tmp=forRotMat*(supportPolys[0].block<3,1>(0,i)+tmp_vector.block<3,1>(0,0));
        //tmp=(supportPolys[0].block<3,1>(0,i)+tmp_vector.block<3,1>(0,0));
        //tmp=forRotMat*(supportPolys[0].block<3,1>(0,i)+feetMod.block<3,1>(0,i));

        globalPoly(0,i)=tmp(0)+gpsPos(0);
        globalPoly(1,i)=tmp(2)+gpsPos(1);
        globalPoly(2,i)=-tmp(1)+gpsPos(2);
    }


    return globalPoly;
}

/* Get IMU data */
void
Controller :: getAttitude(double *rotmat)
{
    for(int i=0; i<3; i++){
        fbck_fgirdRotMat(0,i)=rotmat[i];
        fbck_fgirdRotMat(1,i)=rotmat[i+3];
        fbck_fgirdRotMat(2,i)=rotmat[i+6];
    }


    forRotMat = fgird2ForRotMat.inverse()*fbck_fgirdRotMat;

    forRPY = forRotMat.eulerAngles(2,1,0);



    forRPY(0)=    atan2( forRotMat(2,1),forRotMat(2,2));
    forRPY(1)=    atan2( -forRotMat(2,0),sqrt(forRotMat(2,1)*forRotMat(2,1)+forRotMat(2,2)*forRotMat(2,2)));
    forRPY(2)=    atan2( forRotMat(1,0),forRotMat(1,1));

}

void
Controller :: getRPY(double *rpy)
{

    forRPY(0)=  rpy[0];  
    forRPY(1)=  rpy[1];  
    forRPY(2)=  rpy[2];  

}


void
Controller :: sendStatusToMatlab(int downsampling, const char *IPadr)
{
    static unsigned long int cnt=0;
    cnt++;
    if(cnt%downsampling){
        return;
    }

    std::vector<float>data;

    // angles
    for(int i=0; i<NUM_MOTORS; i++){
        data.push_back((float)angles(i));
    }

    //feedback angles
    for(int i=0; i<NUM_MOTORS; i++){
        data.push_back((float)fbck_angles(i));
    }

    // current
    for(int i=0; i<NUM_MOTORS; i++){
        data.push_back((float)fbck_current(i));
    }

    // raw optoforce
    for(int i=0; i<12; i++){
        data.push_back((float)rawForce(i));
    }

    // feet force - optoforce
    for(int i=0; i<12; i++){
        data.push_back((float)feet_force(i));
    }

    // feet force
    for(int i=0; i<12; i++){
        data.push_back((float)nnEstForce(i));
    }
    
    // IMU
    for(int i=0; i<3; i++){
        data.push_back((float)forRPY(i));
    }
    // kinematics
    for(unsigned int i=0; i<HJfl_g.size(); i++){
        for(unsigned int j=0; j<3; j++)
            data.push_back((float)HJfl_g[i](j, 3));
    }
    for(unsigned int i=0; i<HJfr_g.size(); i++){
        for(unsigned int j=0; j<3; j++)
            data.push_back((float)HJfr_g[i](j, 3));
    }
    for(unsigned int i=0; i<HJhl_g.size(); i++){
        for(unsigned int j=0; j<3; j++)
            data.push_back((float)HJhl_g[i](j, 3));
    }
    for(unsigned int i=0; i<HJhr_g.size(); i++){
        for(unsigned int j=0; j<3; j++)
            data.push_back((float)HJhr_g[i](j, 3));
    }
    for(unsigned int i=0; i<HJs_g.size(); i++){
        for(unsigned int j=0; j<3; j++)
            data.push_back((float)HJs_g[i](j, 3));
    }
    for(unsigned int i=0; i<HJt_g.size(); i++){
        for(unsigned int j=0; j<3; j++)
            data.push_back((float)HJt_g[i](j, 3));
    }




    //cout << data.size() << "\t" << sizeof(float) << endl;
    sendUDP(data.data(), data.size()*sizeof(float), IPadr, 8472);



}



