#include "controller.hpp"

using namespace std;
using namespace Eigen;



/* Modifies trajectories to correct posture (roll and pitch angles) */
void
Controller :: postureControl()
{

}

/* Leg extension reflex */
void
Controller :: legExtension()
{
    static MatrixXd reflexCorrection = MatrixXd::Zero(3,4);
    static double reflexTimeout[4]={0};

    // check if the reflex is disabled
    if(reflON(1)==0){
        return;
    }
    for(int leg=0; leg<4; leg++){

        // timeout
        
        
        // get rid of swing legs
        if(legs_stance(leg)==0){
            reflexCorrection.block<3,1>(0,leg) << 0,0,0;
            reflexTimeout[leg]=0;
            continue;
        }

        // no stance contact
        if(nnEstForce.block<3,1>(0,leg).norm()<extendForceLimit && state==WALKING){

            reflexTimeout[leg] +=dt;


            if(reflexTimeout[leg]>extendTimeout){
                reflexCorrection.block<3,1>(0,leg)+=extendImpulse*dt;

                cout << "EX: " << leg << " corr: " << reflexCorrection.block<3,1>(0,leg).transpose() << endl;
            }


        }

        

        // apply reflex correction
        feetReference_reflex.block<3,1>(0,leg) += AngleAxisd(forTraj(2+3*(leg/2),1)-girdleTraj(2+3*(leg/2),1), Vector3d::UnitZ())*reflexCorrection.block<3,1>(0,leg);

        // decay the correction if the contact is detected
        if(legPhase(leg)<GP.Duty(leg/2)){
            //reflexCorrection.block<3,1>(0,leg) -=reflexCorrection.block<3,1>(0,leg)/(1-legPhase(leg))*dt;
            reflexCorrection.block<3,1>(0,leg)=pt1_vec(MatrixXd::Zero(3,1), reflexCorrection.block<3,1>(0,leg), extendFilterConstant, dt);
        }
        else{
            reflexCorrection.block<3,1>(0,leg) << 0,0,0;
        }
        

    }
}

/* Stumble reflex */
void
Controller :: stumbleReflex()
{
    static MatrixXd reflexCorrection = MatrixXd::Zero(3,4);
    static double reflexTimeout[4]={0};

    // check if the reflex is disabled
    if(reflON(0)==0){
        return;
    }

    for(int leg=0; leg<4; leg++){

        // timeout
        reflexTimeout[leg]-=dt;
        if(reflexTimeout[leg]<0){
            reflexTimeout[leg]=0;
        }

        if(leg==1){
            //cout << legs_stance(leg) << " first lim " << (legPhase(leg)-GP.Duty(leg/2))/(1-GP.Duty(leg/2)) << " sec lim " << (legPhase(leg)-GP.Duty(leg/2))/(1-GP.Duty(leg/2)) << " force " << nnEstForce.block<3,1>(0,leg).norm() << endl;
        }


        // get rid of stance legs
        if(legs_stance(leg)==1 || (legPhase(leg)-GP.Duty(leg/2))/(1-GP.Duty(leg/2))<stumblePhaseLimits(0) || (legPhase(leg)-GP.Duty(leg/2))/(1-GP.Duty(leg/2))>stumblePhaseLimits(1)){
            reflexCorrection.block<3,1>(0,leg) << 0,0,0;
            continue;
        }

        // swing contact detected

        if(nnEstForce.block<3,1>(0,leg).norm()>stumbleForceLimit && state==WALKING && reflexTimeout[leg]==0){

            if(abs(nnEstForce(1,leg))>0.5*stumbleForceLimit){
                reflexTimeout[leg] = stumbleTimeout;
                reflexCorrection.block<3,1>(0,leg)+=stumbleImpulse;
            }
            cout << "ST: " << leg << endl;
        }

        

        // apply reflex correction
        feetReference_reflex.block<3,1>(0,leg) += AngleAxisd(forTraj(2+3*(leg/2),1)-girdleTraj(2+3*(leg/2),1), Vector3d::UnitZ())*reflexCorrection.block<3,1>(0,leg);

        // decay the correction
        if(legPhase(leg)>GP.Duty(leg/2)){
            //reflexCorrection.block<3,1>(0,leg) -=reflexCorrection.block<3,1>(0,leg)/(1-legPhase(leg))*dt;
            reflexCorrection.block<3,1>(0,leg)=pt1_vec(MatrixXd::Zero(3,1), reflexCorrection.block<3,1>(0,leg), stumbleFilterConstant, dt);
        }
        else{
            reflexCorrection.block<3,1>(0,leg) << 0,0,0;
        }
        

    }

    //cout << feetReference_reflex << endl << endl;

    //feetReference.block<3,1>(0,leg) = AngleAxisd(forTraj(2+3*(leg/2),1)-girdleTraj(2+3*(leg/2),1), Vector3d::UnitZ())*trajPoint;
}

