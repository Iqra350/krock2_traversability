#include "controller.hpp"

using namespace std;
using namespace Eigen;


MatrixXd
Controller :: transformation_Ikin_Webots(MatrixXd joint_angles_in, int direction, int shifts)
{

    shifts=shifts==0?0:1;

    if(direction<0){
        //WEBOTS 2 ROBOT
        for(int i=0;i<NUM_MOTORS;i++){
            joint_angles_in(i)=(joint_angles_in(i)-angShiftsCorrIkin2Webots[i]*shifts*my_pi/180.)*angSignsCorrIkin2Webots[i];
        }
    }
    else if(direction>0){
        //ROBOT 2 WEBOTS
        for(int i=0;i<NUM_MOTORS;i++){
            joint_angles_in(i)=joint_angles_in(i)*angSignsCorrIkin2Webots[i] + angShiftsCorrIkin2Webots[i]*shifts*my_pi/180.;
        }
    }

    return joint_angles_in;
}

MatrixXd
Controller :: transformation_Ikin_Robot(MatrixXd joint_angles_in, int direction, int shifts)
{

    shifts=shifts==0?0:1;

    if(direction<0){
        //ROBOT 2 IKIN
        for(int i=0;i<NUM_MOTORS;i++){
            joint_angles_in(i)=(joint_angles_in(i)-angShiftsCorrIkin2Robot[i]*shifts*my_pi/180.)*angSignsCorrIkin2Robot[i];
        }
    }
    else if(direction>0){
        //IKIN 2 ROBOT
        for(int i=0;i<NUM_MOTORS;i++){
            joint_angles_in(i)=joint_angles_in(i)*angSignsCorrIkin2Robot[i] + angShiftsCorrIkin2Robot[i]*shifts*my_pi/180.;
        }
    }

    return joint_angles_in;
}


MatrixXd
Controller :: flipKinematics(MatrixXd joint_angles_in)
{
    MatrixXd flipped_angles=joint_angles_in;
    // flip left-right
    flipped_angles.block<4,1>(0,0)=joint_angles_in.block<4,1>(4,0);
    flipped_angles.block<4,1>(4,0)=joint_angles_in.block<4,1>(0,0);
    flipped_angles.block<4,1>(8,0)=joint_angles_in.block<4,1>(12,0);
    flipped_angles.block<4,1>(12,0)=joint_angles_in.block<4,1>(8,0);

    // correct axes
    for(int i=0; i<4; i++){
        // pitch
        //flipped_angles(i*4,0)*=-1;
        // yaw
        flipped_angles(i*4+1,0)*=-1;
        // roll
        flipped_angles(i*4+2,0)*=-1;
        // knee
        flipped_angles(i*4+3,0)*=-1;
    }

    // flip spine
    flipped_angles.block(16,0,NUM_MOTORS-16,1)*=-1;

    // flip roll
/*    flipped_angles(2,0)*=-1;
    flipped_angles(2,1)*=-1;
    flipped_angles(2,2)*=-1;
    flipped_angles(2,3)*=-1;*/

    return flipped_angles;


}
