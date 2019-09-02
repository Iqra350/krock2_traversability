
#include "robotSim.hpp"






#define USE_CAMERA 0
    



extern double get_timestamp();



/* RobotSim constructor */
RobotSim :: RobotSim(int TIME_STEP)
{
    cout << "RobotSim constructor"<<endl;

    /*
    gyro=new webots::Gyro("head_gyro");                    
    gyro->enable(TIME_STEP);                                

    acc=new webots::Accelerometer("head_acc");
    acc->enable(TIME_STEP);
    compass=new webots::Compass("compass_fgirdle");
    compass->enable(TIME_STEP);
    
    gyroData=gyro->getValues();
    accData=acc->getValues();
    compassData=compass->getValues();
    */

    gps_fgird=new webots::GPS("gps_fgirdle");
    gps_fgird->enable(TIME_STEP);
    gps_hgird=new webots::GPS("gps_hgirdle");
    gps_hgird->enable(TIME_STEP);

    gpsDataFgird=gps_fgird->getValues();
    gpsDataHgird=gps_hgird->getValues();




    roboDef=getFromDef("ROBOT"); 
    
    roboPos=roboDef->getField("translation");
    roboRot=roboDef->getField("rotation");


    //webots::Field fldTranslation = roboDef.getField("translation");
    // INITIALIZE MARKERS
/*
    FL_marker=getFromDef("FL_MARKER"); 
    FR_marker=getFromDef("FR_MARKER"); 
    HL_marker=getFromDef("HL_MARKER"); 
    HR_marker=getFromDef("HR_MARKER"); 
*/
    CoM_marker=getFromDef("COM_MARKER"); 
    CoM_marker_pos=CoM_marker->getField("translation");

    CoM_marker_proj=getFromDef("COM_MARKER_PROJ"); 
    CoM_marker_pos_proj=CoM_marker_proj->getField("translation");

    supportPolyDEF=getFromDef("SUPPORT_POLY"); 
    
    vector<string> RM_NAMES = 
    {
        "FLpitch_motor", "FLyaw_motor", "FLroll_motor", "FLknee_motor",
        "FRpitch_motor", "FRyaw_motor", "FRroll_motor", "FRknee_motor",
        "HLpitch_motor", "HLyaw_motor", "HLroll_motor", "HLknee_motor",
        "HRpitch_motor", "HRyaw_motor", "HRroll_motor", "HRknee_motor",
        "spine1_motor", "spine2_motor",
        "neck1_motor", "tail1_motor", "tail2_motor"
    };

    vector<string> PS_NAMES = 
    {
        "FLpitch_sensor", "FLyaw_sensor", "FLroll_sensor", "FLknee_sensor",
        "FRpitch_sensor", "FRyaw_sensor", "FRroll_sensor", "FRknee_sensor",
        "HLpitch_sensor", "HLyaw_sensor", "HLroll_sensor", "HLknee_sensor",
        "HRpitch_sensor", "HRyaw_sensor", "HRroll_sensor", "HRknee_sensor",
        "spine1_sensor", "spine2_sensor",
        "neck1_sensor", "tail1_sensor", "tail2_sensor"
    };



    const char *TOUCH_SENSOR_NAMES[N_TOUCH_SENSORS] =
    {
    "fl_touch", "fr_touch", "hl_touch", "hr_touch",
    };


    
    rm_motor.resize(NUM_MOTORS);
    ps_motor.resize(NUM_MOTORS);
    // get the servos
    cout << "connecting to motors" << endl;
    for(int i=0;i<NUM_MOTORS;i++)
    {   
        //rm_motor[i].getMotor(RM_NAMES[i]);
        rm_motor[i] = webots::Robot::getMotor(RM_NAMES[i]);
        ps_motor[i] = webots::Robot::getPositionSensor(PS_NAMES[i]);
        ps_motor[i]->enable(TIME_STEP);
        rm_motor[i]->enableTorqueFeedback(TIME_STEP);
    }
    cout << "motors collected" << endl;

  



    for(int i=0;i<N_TOUCH_SENSORS;i++){
        touch_sensor[i] = new webots::TouchSensor(TOUCH_SENSOR_NAMES[i]);
        touch_sensor[i]->enable(TIME_STEP);
    }

    tsdefFL=getFromDef("TS_FL"); 
    tsdefFR=getFromDef("TS_FR"); 
    tsdefHL=getFromDef("TS_HL"); 
    tsdefHR=getFromDef("TS_HR"); 
    

}

/* Sets the position of all servos to a table of angles in radians */
void
RobotSim :: setAngles(double *table, int *ids)
{
        for(int i=0; i<NUM_MOTORS;  i++){
            if(ids[i]){
                rm_motor[i]->setPosition(table[i]);
            }
        }
}

/* Sets the torques of all servos to a table of torques */
void
RobotSim :: torques(double *table, int *ids)
{
        for(int i=0; i<NUM_MOTORS;  i++){
            if(ids[i]){
                //servo[i]->setForce(table[i]);
            }
        }
}






/* Reads positions, torques   */
void
RobotSim::getPositionTorques(double *d_posture, double *d_torques)
{
    for(int i=0;i<NUM_MOTORS;i++)
    {        
        d_posture[i]=ps_motor[i]->getValue(); 
        d_torques[i]=rm_motor[i]->getForceFeedback(); 
    }
}

/* Reads positions, torques and IMU data */
void
RobotSim::GetPosition(double *gpsDataFgird_i, double *gpsDataHgird_i)
{
    gpsDataFgird_i[0]=gpsDataFgird[0];
    gpsDataFgird_i[1]=gpsDataFgird[1];
    gpsDataFgird_i[2]=gpsDataFgird[2];

    gpsDataHgird_i[0]=gpsDataHgird[0];
    gpsDataHgird_i[1]=gpsDataHgird[1];
    gpsDataHgird_i[2]=gpsDataHgird[2];
}
/* Initializes IMU */
void
RobotSim::InitIMU()
{
    
    fgirdle=getFromDef("GPS_FGIRDLE"); 
    rotMat = fgirdle->getOrientation();
}



/* Reads positions, torques and IMU data */
void
RobotSim::ReadIMUAttitude(double *rotmat)
{
    

    double rotmat_tmp[9];
    rotMat = fgirdle->getOrientation();

    for(int i=0; i<9; i++){
        rotmat[i]=rotMat[i];
        rotmat_tmp[i]=rotMat[i];
    }
    rotmat[3]=rotmat_tmp[6];
    rotmat[4]=rotmat_tmp[7];
    rotmat[5]=rotmat_tmp[8];
    rotmat[6]=rotmat_tmp[3];
    rotmat[7]=rotmat_tmp[4];
    rotmat[8]=rotmat_tmp[5];

    for(int i=0; i<9; i++){
        rotmat_tmp[i]=rotmat[i];
    }

    rotmat[1]=rotmat_tmp[2];
    rotmat[4]=rotmat_tmp[5];
    rotmat[7]=rotmat_tmp[8];
    rotmat[2]=rotmat_tmp[1];
    rotmat[5]=rotmat_tmp[4];
    rotmat[8]=rotmat_tmp[7];

    rotmat[1]*=-1;
    rotmat[3]*=-1;
    rotmat[5]*=-1;
    rotmat[7]*=-1;

}
/* Reads touch sensor data */
void
RobotSim::ReadTouchSensors(double *ts_data)
{

    ts_fl=touch_sensor[0]->getValues();
    ts_fr=touch_sensor[1]->getValues();
    ts_hl=touch_sensor[2]->getValues();
    ts_hr=touch_sensor[3]->getValues();
    /*for(int i=0; i<3; i++){
        ts_data[i]=ts_fl[i];
        ts_data[i+3]=ts_fr[i];
        ts_data[i+6]=ts_hl[i];
        ts_data[i+9]=ts_hr[i];
    }*/

    ts_data[0]=ts_fl[0];
    ts_data[1]=-ts_fl[2];
    ts_data[2]=ts_fl[1];

    ts_data[0+3]=ts_fr[0];
    ts_data[1+3]=-ts_fr[2];
    ts_data[2+3]=ts_fr[1];

    ts_data[0+6]=ts_hl[0];
    ts_data[1+6]=-ts_hl[2];
    ts_data[2+6]=ts_hl[1];

    ts_data[0+9]=ts_hr[0];
    ts_data[1+9]=-ts_hr[2];
    ts_data[2+9]=ts_hr[1];
}


/* Quits simulation */
void
RobotSim::killSimulation()
{
    simulationQuit(0);
}





/* Changes color of the feet to indiccate stuff */
void
RobotSim::ColorMarkers(double *logic, double trans, double *col1)
{
    webots::Field *transparency;
    webots::Field *color;
    double col2[3]={0.8, 0.8, 0.8}, col[12];
    double p[3];

    for(int i=0; i<4; i++){
        if(logic[i]){
            col[0+i*3]=col1[0+i*3];
            col[1+i*3]=col1[1+i*3];
            col[2+i*3]=col1[2+i*3];
        }
        else{
            col[0+i*3]=col2[0];
            col[1+i*3]=col2[1];
            col[2+i*3]=col2[2];
        }
    }

    //FL
    transparency = FL_marker->getField("transparency");
    transparency->setSFFloat(trans);

    p[0]=col[0];
    p[1]=col[1];
    p[2]=col[2];
    color = FL_marker->getField("diffuseColor");
    color->setSFColor((const double*) p);

    //FR
    transparency = FR_marker->getField("transparency");
    transparency->setSFFloat(trans);

    p[0]=col[0+3];
    p[1]=col[1+3];
    p[2]=col[2+3];
    color = FR_marker->getField("diffuseColor");
    color->setSFColor((const double*) p);

    //HL
    transparency = HL_marker->getField("transparency");
    transparency->setSFFloat(trans);

    p[0]=col[0+6];
    p[1]=col[1+6];
    p[2]=col[2+6];
    color = HL_marker->getField("diffuseColor");
    color->setSFColor((const double*) p);

    //HR
    transparency = HR_marker->getField("transparency");
    transparency->setSFFloat(trans);

    p[0]=col[0+9];
    p[1]=col[1+9];
    p[2]=col[2+9];
    color = HR_marker->getField("diffuseColor");
    color->setSFColor((const double*) p);
}

/* Sets position ond orientation of the Robot */
void 
RobotSim::setPositionRotation(double *p, double *r)
{
    roboPos->setSFVec3f((const double*)p);
    roboRot->setSFRotation(r);
}

/* Sets position ond orientation of the Robot */
void 
RobotSim::setPositionOfRobot(double *p)
{
    roboPos->setSFVec3f((const double*)p);
}


/* Sets position of CoM marker */
void 
RobotSim::setCoMmarker(double *p)
{
    CoM_marker_pos->setSFVec3f((const double*)p);
    //p[1]=0.0011;
    //CoM_marker_pos_proj->setSFVec3f((const double*)p);
}


/* Reads compass data */
void
RobotSim::GetCompass(double *data)
{
    data[0]=compassData[0];
    data[1]=compassData[1];
    data[2]=compassData[2];
}



/* Reads compass data */
void
RobotSim::GetFeetGPS(double *FL_feet_gpos, double *FR_feet_gpos, double *HL_feet_gpos, double *HR_feet_gpos)
{
    static const double *FL_feet_gposC=tsdefFL->getPosition();
    static const double *FR_feet_gposC=tsdefFR->getPosition();
    static const double *HL_feet_gposC=tsdefHL->getPosition();
    static const double *HR_feet_gposC=tsdefFR->getPosition();
   // cout<<FL_feet_gpos[0]<<FL_feet_gpos[1]<<FL_feet_gpos[2]<<endl;
    
    for(int i=0; i<3; i++){
        FL_feet_gpos[i]=FL_feet_gposC[i];
        FR_feet_gpos[i]=FR_feet_gposC[i];
        HL_feet_gpos[i]=HL_feet_gposC[i];
        HR_feet_gpos[i]=HR_feet_gposC[i];
    }
}

/* Draws support polygon */
void
RobotSim::setSupportPoly(MatrixXd globalPoly, MatrixXd stance, double height){

    webots::Field *supportPolyCoord;
    webots::Field *supportPolycoordIndex;
    webots::Node *Coordinate;

    Coordinate=getFromDef("poly_coordinate"); 

    supportPolyCoord=Coordinate->getField("point");
    supportPolycoordIndex=supportPolyDEF->getField("coordIndex");


    //cout << globalPoly << endl << endl;


    const double p1[3]={globalPoly(0,0), height, globalPoly(2,0)};
    const double p2[3]={globalPoly(0,1), height, globalPoly(2,1)};
    const double p3[3]={globalPoly(0,3), height, globalPoly(2,3)};
    const double p4[3]={globalPoly(0,2), height, globalPoly(2,2)};

    
    if(globalPoly(0,0)>-100 || globalPoly(0,0)<100){
        

        supportPolyCoord->setMFVec3f(0, p1);
        supportPolyCoord->setMFVec3f(1, p2);
        supportPolyCoord->setMFVec3f(2, p3);
        supportPolyCoord->setMFVec3f(3, p4);


        int j=0;
        if(stance(2)){
            supportPolycoordIndex->setMFInt32(j, 3);
            j++;
        }
        if(stance(3)){
            supportPolycoordIndex->setMFInt32(j, 2);
            j++;
        }
        if(stance(1)){
            supportPolycoordIndex->setMFInt32(j, 1);
            j++;
        }
        if(stance(0)){
            supportPolycoordIndex->setMFInt32(j, 0);
            j++;
        }
        supportPolycoordIndex->setMFInt32(j, -1);

    }
    
}