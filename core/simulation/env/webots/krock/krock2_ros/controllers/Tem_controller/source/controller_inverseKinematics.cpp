#include "controller.hpp"

using namespace std;
using namespace Eigen;

extern int IS_SIMULATION, USE_JOYSTICK, SWIM, SPINE_COMPENSATION_WITH_FEEDBACK, USE_REFLEXES, LOG_DATA;

// ---------------------------------------------- SETUP SPINE OPTIMIZATION ----------------------------------------------------------
typedef dlib::matrix<double,0,1> column_vector;


class optimizer_spineOscillations
{
    public:
        optimizer_spineOscillations (
            Controller *ctrl_input, Vector2d old_solution_input
        )
        {
            ctrl=ctrl_input;
            old_solution=old_solution_input;
        }

        double operator() ( const column_vector& x) const
        //double operator() ( const double& x) const
        {
            
            MatrixXd q_trunk(Ntrunk,1);
            q_trunk=ctrl->q0_trunk_from_spline*x(0) + MatrixXd::Constant(Ntrunk,1,x(1)); 
            // rotate front girdle fore x(1)
            Vector3d fgird_pos; 
            fgird_pos(0)=ctrl->girdleTraj(0,1);
            fgird_pos(1)=ctrl->girdleTraj(1,1);
            fgird_pos(2)=ctrl->girdleTraj(2,1);

            Vector3d hgird_pos=ctrl->trunkForwardKinematics(fgird_pos, q_trunk);
            
            // find minimal distance of the hind girdle form the stored trajectory
            double mindist=9999, dist; 
            int indx=0;
            for(int i=50; i<300; i++){
                dist=sqrt((ctrl->forTrajHist(0,i)-hgird_pos(0))*(ctrl->forTrajHist(0,i)-hgird_pos(0)) + 
                          (ctrl->forTrajHist(1,i)-hgird_pos(1))*(ctrl->forTrajHist(1,i)-hgird_pos(1)));
                if(dist<mindist){
                    mindist=dist;
                    indx=i;
                }
            }

            // cost function
            double f1, f2, f3, f4;
            
            // distance from the curve
            f1 = 300*mindist;

            // intrusivness of the solution -> scaling close to 1, offset close to 0
            f2 = (x(0)-1)*(x(0)-1) + x(1)*x(1); 

            // hind girdle angle
            //f3 = 0*(hgird_pos(2) - ctrl->girdleCpgOutput(1))*(hgird_pos(2) - ctrl->girdleCpgOutput(1));
            double hgirdAngRef;
            hgirdAngRef=ctrl->forTraj(Ntrunk,1) + ctrl->girdleCpgOutput(1);
            f3 = 20*(hgird_pos(2) - hgirdAngRef)*(hgird_pos(2) - hgirdAngRef);
            
            // distance from the old solution
            f4= 30*((x(0)-old_solution(0))*(x(0)-old_solution(0)) + (x(1)-old_solution(1))*(x(1)-old_solution(1)));
                

            return f1 + f2 + f3 + f4;
        }

    private:
        Controller *ctrl;
        Vector2d old_solution;
};


/* Calculates forward kinematics for a single leg */
vector<Matrix4d>
Controller :: legKinematics(Vector4d q, int leg)
{

    static vector<Matrix4d> H(5);

    if(leg==0){     //FL
        H[0] <<  1, 0, 0, FL_kin(0,0),
                0, cos(q(0)), -sin(q(0)), FL_kin(0,1),
                0, sin(q(0)), cos(q(0)), FL_kin(0,2),
                0, 0, 0, 1;
        H[1] <<  cos(q(1)), -sin(q(1)), 0, FL_kin(1,0),
                sin(q(1)), cos(q(1)), 0, FL_kin(1,1),
                0, 0, 1, FL_kin(1,2),
                0, 0, 0, 1;
        H[2] <<  cos(q(2)), 0, sin(q(2)), FL_kin(2,0),
                0, 1, 0, FL_kin(2,1),
                -sin(q(2)), 0, cos(q(2)), FL_kin(2,2),
                0, 0, 0, 1;
        H[3] <<  cos(q(3)), -sin(q(3)), 0, FL_kin(3,0),
                sin(q(3)), cos(q(3)), 0, FL_kin(3,1),
                0, 0, 1, FL_kin(3,2),
                0, 0, 0, 1;
        H[4] <<  1, 0, 0, FL_kin(4,0),
                0, 1, 0, FL_kin(4,1),
                0, 0, 1, FL_kin(4,2),
                0, 0, 0, 1;
    }
    else if (leg==1){       //FR
        H[0] <<  1, 0, 0, FR_kin(0,0),
                0, cos(q(0)), -sin(q(0)), FR_kin(0,1),
                0, sin(q(0)), cos(q(0)), FR_kin(0,2),
                0, 0, 0, 1;
        H[1] <<  cos(q(1)), -sin(q(1)), 0, FR_kin(1,0),
                sin(q(1)), cos(q(1)), 0, FR_kin(1,1),
                0, 0, 1, FR_kin(1,2),
                0, 0, 0, 1;
        H[2] <<  cos(q(2)), 0, sin(q(2)), FR_kin(2,0),
                0, 1, 0, FR_kin(2,1),
                -sin(q(2)), 0, cos(q(2)), FR_kin(2,2),
                0, 0, 0, 1;
        H[3] <<  cos(q(3)), -sin(q(3)), 0, FR_kin(3,0),
                sin(q(3)), cos(q(3)), 0, FR_kin(3,1),
                0, 0, 1, FR_kin(3,2),
                0, 0, 0, 1;
        H[4] <<  1, 0, 0, FR_kin(4,0),
                0, 1, 0, FR_kin(4,1),
                0, 0, 1, FR_kin(4,2),
                0, 0, 0, 1;
    }
    else if (leg==2){       //HL
        H[0] <<  1, 0, 0, HL_kin(0,0),
                0, cos(q(0)), -sin(q(0)), HL_kin(0,1),
                0, sin(q(0)), cos(q(0)), HL_kin(0,2),
                0, 0, 0, 1;
        H[1] <<  cos(q(1)), -sin(q(1)), 0, HL_kin(1,0),
                sin(q(1)), cos(q(1)), 0, HL_kin(1,1),
                0, 0, 1, HL_kin(1,2),
                0, 0, 0, 1;
        H[2] <<  cos(q(2)), 0, sin(q(2)), HL_kin(2,0),
                0, 1, 0, HL_kin(2,1),
                -sin(q(2)), 0, cos(q(2)), HL_kin(2,2),
                0, 0, 0, 1;
        H[3] <<  cos(q(3)), -sin(q(3)), 0, HL_kin(3,0),
                sin(q(3)), cos(q(3)), 0, HL_kin(3,1),
                0, 0, 1, HL_kin(3,2),
                0, 0, 0, 1;
        H[4] <<  1, 0, 0, HL_kin(4,0),
                0, 1, 0, HL_kin(4,1),
                0, 0, 1, HL_kin(4,2),
                0, 0, 0, 1;
    }
    else if (leg==3){       //HR
        H[0] <<  1, 0, 0, HR_kin(0,0),
                0, cos(q(0)), -sin(q(0)), HR_kin(0,1),
                0, sin(q(0)), cos(q(0)), HR_kin(0,2),
                0, 0, 0, 1;
        H[1] <<  cos(q(1)), -sin(q(1)), 0, HR_kin(1,0),
                sin(q(1)), cos(q(1)), 0, HR_kin(1,1),
                0, 0, 1, HR_kin(1,2),
                0, 0, 0, 1;
        H[2] <<  cos(q(2)), 0, sin(q(2)), HR_kin(2,0),
                0, 1, 0, HR_kin(2,1),
                -sin(q(2)), 0, cos(q(2)), HR_kin(2,2),
                0, 0, 0, 1;
        H[3] <<  cos(q(3)), -sin(q(3)), 0, HR_kin(3,0),
                sin(q(3)), cos(q(3)), 0, HR_kin(3,1),
                0, 0, 1, HR_kin(3,2),
                0, 0, 0, 1;
        H[4] <<  1, 0, 0, HR_kin(4,0),
                0, 1, 0, HR_kin(4,1),
                0, 0, 1, HR_kin(4,2),
                0, 0, 0, 1;
    }


    //return H[0]*H[1]*H[2]*H[3]*H[4];
    return H;
}

/* Returns jacobian matrix of one leg */
Matrix<double, 3, 4>
Controller :: Jacob(Vector4d q, int leg)
{
    static MatrixXd J(3, 4);
    static double S1, S2, S3, S4, C1, C2, C3, C4;
    S1 = sin(q(0));
    S2 = sin(q(1));
    S3 = sin(q(2));
    S4 = sin(q(3));
    C1 = cos(q(0));
    C2 = cos(q(1));
    C3 = cos(q(2));
    C4 = cos(q(3));


    if(leg==0){     //FL
        J(0) = 0;
        J(1) = 0.132*C4*(C1*S3 - 1.0*C3*S1*S2) - 0.1474*C2*S1 - 0.132*C2*S1*S4;
        J(2) = 0.1474*C1*C2 + 0.132*C4*(S1*S3 + C1*C3*S2) + 0.132*C1*C2*S4;
        J(3) = - 0.1474*C2 - 0.132*C2*S4 - 0.132*C3*C4*S2;
        J(4) = -0.0022*C1*(67.0*S2 + 60.0*S2*S4 - 60.0*C2*C3*C4);
        J(5) = -0.0022*S1*(67.0*S2 + 60.0*S2*S4 - 60.0*C2*C3*C4);
        J(6) = -0.132*C2*C4*S3;
        J(7) = 0.132*C4*(C3*S1 - 1.0*C1*S2*S3);
        J(8) = -0.132*C4*(C1*C3 + S1*S2*S3);
        J(9) = - 0.132*C4*S2 - 0.132*C2*C3*S4;
        J(10) = 0.132*C1*C2*C4 - 0.132*S4*(S1*S3 + C1*C3*S2);
        J(11) = 0.132*S4*(C1*S3 - 1.0*C3*S1*S2) + 0.132*C2*C4*S1;
    }
    else if (leg==1){       //FR
        J(0) = 0;
        J(1) = 0.1474*C2*S1 + 0.132*C4*(C1*S3 - 1.0*C3*S1*S2) - 0.132*C2*S1*S4;
        J(2) = 0.132*C4*(S1*S3 + C1*C3*S2) - 0.1474*C1*C2 + 0.132*C1*C2*S4;
        J(3) = 0.1474*C2 - 0.132*C2*S4 - 0.132*C3*C4*S2;
        J(4) = 0.0022*C1*(67.0*S2 - 60.0*S2*S4 + 60.0*C2*C3*C4);
        J(5) = 0.0022*S1*(67.0*S2 - 60.0*S2*S4 + 60.0*C2*C3*C4);
        J(6) = -0.132*C2*C4*S3;
        J(7) = 0.132*C4*(C3*S1 - 1.0*C1*S2*S3);
        J(8) = -0.132*C4*(C1*C3 + S1*S2*S3);
        J(9) = - 0.132*C4*S2 - 0.132*C2*C3*S4;
        J(10) = 0.132*C1*C2*C4 - 0.132*S4*(S1*S3 + C1*C3*S2);
        J(11) = 0.132*S4*(C1*S3 - 1.0*C3*S1*S2) + 0.132*C2*C4*S1;
    }
    else if (leg==2){       //HL
        J(0) = 0;
        J(1) = 0.132*C2*S1*S4 - 0.132*C4*(C1*S3 - 1.0*C3*S1*S2) - 0.1474*C2*S1;
        J(2) = 0.1474*C1*C2 - 0.132*C4*(S1*S3 + C1*C3*S2) - 0.132*C1*C2*S4;
        J(3) = 0.132*C2*S4 - 0.1474*C2 + 0.132*C3*C4*S2;
        J(4) = -0.0022*C1*(67.0*S2 - 60.0*S2*S4 + 60.0*C2*C3*C4);
        J(5) = -0.0022*S1*(67.0*S2 - 60.0*S2*S4 + 60.0*C2*C3*C4);
        J(6) = 0.132*C2*C4*S3;
        J(7) = -0.132*C4*(C3*S1 - 1.0*C1*S2*S3);
        J(8) = 0.132*C4*(C1*C3 + S1*S2*S3);
        J(9) = 0.132*C4*S2 + 0.132*C2*C3*S4;
        J(10) = 0.132*S4*(S1*S3 + C1*C3*S2) - 0.132*C1*C2*C4;
        J(11) = - 0.132*S4*(C1*S3 - 1.0*C3*S1*S2) - 0.132*C2*C4*S1;
    }
    else if (leg==3){       //HR
        J(0) = 0;
        J(1) = 0.1474*C2*S1 - 0.132*C4*(C1*S3 - 1.0*C3*S1*S2) + 0.132*C2*S1*S4;
        J(2) = - 0.1474*C1*C2 - 0.132*C4*(S1*S3 + C1*C3*S2) - 0.132*C1*C2*S4;
        J(3) = 0.1474*C2 + 0.132*C2*S4 + 0.132*C3*C4*S2;
        J(4) = 0.0022*C1*(67.0*S2 + 60.0*S2*S4 - 60.0*C2*C3*C4);
        J(5) = 0.0022*S1*(67.0*S2 + 60.0*S2*S4 - 60.0*C2*C3*C4);
        J(6) = 0.132*C2*C4*S3;
        J(7) = -0.132*C4*(C3*S1 - 1.0*C1*S2*S3);
        J(8) = 0.132*C4*(C1*C3 + S1*S2*S3);
        J(9) = 0.132*C4*S2 + 0.132*C2*C3*S4;
        J(10) = 0.132*S4*(S1*S3 + C1*C3*S2) - 0.132*C1*C2*C4;
        J(11) = - 0.132*S4*(C1*S3 - 1.0*C3*S1*S2) - 0.132*C2*C4*S1;
    }



    return J;
}


/* Gets entire forward kinematics */
void 
Controller :: forwardKinematics()
{

    static MatrixXd feetLocations_old(3,4);
    


    // Front girdle
    Fgird.setZero();
    Fgird.block<3,3>(0, 0)=MatrixXd::Identity(3,3);
    Fgird(3,3)=1;

    //Trunk
    HJs[0]<<cos(joint_angles(16)), -sin(joint_angles(16)), 0, -trunk_kin(0),
                sin(joint_angles(16)), cos(joint_angles(16)), 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
    HJs_g[0]=HJs[0];  
    for(int i=1; i<Ntrunk; i++){
        HJs[i]<<cos(joint_angles(16+i)), -sin(joint_angles(16+i)), 0, -trunk_kin(i),
                sin(joint_angles(16+i)), cos(joint_angles(16+i)), 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
        HJs_g[i]=HJs_g[i-1]*HJs[i];
    }
    HJs[Ntrunk]<<1, 0, 0, -trunk_kin(Ntrunk),
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
    HJs_g[Ntrunk]=HJs_g[Ntrunk-1]*HJs[Ntrunk];
    Hgird=HJs_g[Ntrunk];

    // LEGS
    HJfl=legKinematics(joint_angles.block<4,1>(0,0), 0);
    HJfr=legKinematics(joint_angles.block<4,1>(4,0), 1);
    HJhl=legKinematics(joint_angles.block<4,1>(8,0), 2);
    HJhr=legKinematics(joint_angles.block<4,1>(12,0), 3);

    // FRONT LEFT LEG
    HJfl_g[0]=Fgird*HJfl[0];
    for(int i=1;i<5;i++){
        HJfl_g[i]=HJfl_g[i-1]*HJfl[i];
    }
    // FRONT RIGHT LEG
    HJfr_g[0]=Fgird*HJfr[0];
    for(int i=1;i<5;i++){
        HJfr_g[i]=HJfr_g[i-1]*HJfr[i];
    }

    // HIND LEFT LEG
    HJhl_g[0]=Hgird*HJhl[0];
    for(int i=1;i<5;i++){
        HJhl_g[i]=HJhl_g[i-1]*HJhl[i];
    }

    // HIND RIGHT LEG
    HJhr_g[0]=Hgird*HJhr[0];
    for(int i=1;i<5;i++){
        HJhr_g[i]=HJhr_g[i-1]*HJhr[i];
    }

    

    // TRANSFORM EVERYTHING INTO ROBOT's COO FRAME (passing through both girdles)
    static Matrix4d HTransform, HTranslate, HRotate;
   
    hgirdlePhi = atan2(Hgird(1, 3)-Fgird(1, 3), -Hgird(0, 3)+Fgird(0, 3));

    HTransform.block(0,0,3,3) << cos(hgirdlePhi), -sin(hgirdlePhi), 0,    sin(hgirdlePhi), cos(hgirdlePhi), 0,    0, 0, 1;
    HTransform.block(0,3,3,1) = -Fgird.block<3,1>(0,3);
    HTransform.block(3,0,1,4) << 0, 0, 0, 1;

    HTranslate=Matrix4d::Identity();
    HTranslate.block<3,1>(0,3) = -Fgird.block<3,1>(0,3);
    HRotate=Matrix4d::Identity();
    HRotate.block<3,3>(0,0) << cos(hgirdlePhi), -sin(hgirdlePhi), 0,    sin(hgirdlePhi), cos(hgirdlePhi), 0,    0, 0, 1;

    HTransform=HTransform.inverse().eval();
    HTransform=HRotate*HTranslate;
    for(int i=0;i<5;i++){
        HJfl_g[i]=HTransform*HJfl_g[i];
        HJfr_g[i]=HTransform*HJfr_g[i];
        HJhl_g[i]=HTransform*HJhl_g[i];
        HJhr_g[i]=HTransform*HJhr_g[i];
    }

    for(int i=0;i<Ntrunk;i++){
        HJs_g[i]=HTransform*HJs_g[i];
    }

    Hgird=HTransform*Hgird;
    Fgird=HTransform*Fgird;


    // HEAD
    if(Nneck>0){
        Hhead<<cos(-joint_angles(18)), -sin(-joint_angles(18)), 0, neck_kin(0),
                sin(-joint_angles(18)), cos(-joint_angles(18)), 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
        Hhead=Fgird*Hhead;
    }
    

    // TAIL
    if(Ntail>0){
        HJt[0]<<cos(joint_angles(19)), -sin(joint_angles(19)), 0, -tail_kin(0),
                    sin(joint_angles(19)), cos(joint_angles(19)), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
        HJt_g[0]=Hgird*HJt[0]; 

        HJt[1]<<cos(joint_angles(20)), -sin(joint_angles(20)), 0, -tail_kin(1),
                    sin(joint_angles(20)), cos(joint_angles(20)), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1;
        HJt_g[1]=HJt_g[0]*HJt[1];         
    }          

    // get fgird rotation in respect to the 1-frame of reference
    fgird2ForRotMat=AngleAxisd(hgirdlePhi, Vector3d::UnitZ());


    feetLocations.block<3,1>(0,0)=HJfl_g[4].block<3,1>(0,3);
    feetLocations.block<3,1>(0,1)=HJfr_g[4].block<3,1>(0,3);
    feetLocations.block<3,1>(0,2)=HJhl_g[4].block<3,1>(0,3);
    feetLocations.block<3,1>(0,3)=HJhr_g[4].block<3,1>(0,3);


    feetVelocities=feetLocations - feetLocations_old;

    feetLocations_old=feetLocations;


    // local feet locations
    
    feetLocationsLocal.block<3,1>(0,0)=(HJfl[0]*HJfl[1]*HJfl[2]*HJfl[3]*HJfl[4]).block<3, 1>(0, 3);
    feetLocationsLocal.block<3,1>(0,1)=(HJfr[0]*HJfr[1]*HJfr[2]*HJfr[3]*HJfr[4]).block<3, 1>(0, 3);
    feetLocationsLocal.block<3,1>(0,2)=(HJhl[0]*HJhl[1]*HJhl[2]*HJhl[3]*HJhl[4]).block<3, 1>(0, 3);
    feetLocationsLocal.block<3,1>(0,3)=(HJhr[0]*HJhr[1]*HJhr[2]*HJhr[3]*HJhr[4]).block<3, 1>(0, 3);
    
    

}

/* Gets trunk forward kinematics - returns position of Hind girdle for the given position of Front girdle and spine angles*/
Vector3d
Controller :: trunkForwardKinematics(Vector3d fgird_pos, MatrixXd q_trunk_in)
{
    MatrixXd spine_joints(3,Ntrunk);
    //first joint
    spine_joints(0,0)=fgird_pos(0) + cos(fgird_pos(2))*(-trunk_kin(0));
    spine_joints(1,0)=fgird_pos(1) + sin(fgird_pos(2))*(-trunk_kin(0));
    spine_joints(2,0)=fgird_pos(2);

    //spine
    for(int i=1; i<Ntrunk; i++){
        spine_joints(2,i)=q_trunk_in(i-1)+spine_joints(2,i-1);
        spine_joints(0,i)=spine_joints(0,i-1) + cos( spine_joints(2,i) ) * (-trunk_kin(i));
        spine_joints(1,i)=spine_joints(1,i-1) + sin( spine_joints(2,i) ) * (-trunk_kin(i));
    }

    //hind girdle
    Vector3d hgird_pos;
    hgird_pos(2)=spine_joints(2,Ntrunk-1) + q_trunk_in(Ntrunk-1);
    hgird_pos(0)=spine_joints(0,Ntrunk-1) + cos( hgird_pos(2) ) * (-trunk_kin(Ntrunk));
    hgird_pos(1)=spine_joints(1,Ntrunk-1) + sin( hgird_pos(2) ) * (-trunk_kin(Ntrunk));

    return hgird_pos;

}

/* Solves spine inverse kinematics */
void
Controller :: trunkInverseKinematics()
{
    // ------------------------------------------ RUNNING GIRDLE CPG AND SPLINING STUFF -------------------------------
    // spline parameter
    static MatrixXd t_spline(1, Ntrunk+2);
    static MatrixXd p_spline(2, t_spline.size());
    static bool is_init;
    static Vector2d old_solution;
    if(!is_init)
    {
        old_solution << 1,0;
        // spline parameter
        double spinecumsum=trunk_kin(0);
        t_spline(0)=0;
        for(int i=0; i<Ntrunk; i++)
        {
            t_spline(i+1)=spinecumsum/IG;
            spinecumsum+=trunk_kin(i+1);
        }
        t_spline(Ntrunk+1)=spinecumsum/IG;
    }
    // spline points and vectors
    static Vector2d p0_spline, p1_spline, m0_spline, m1_spline;
    p0_spline=forTraj.block(0,0,2,1);
    p1_spline=forTraj.block(3,0,2,1); 

    //girdleCpgOutput*=(1-abs(gamma));
    m0_spline(0)=-cos(forTraj(2,0) + girdleCpgOutput(0));
    m0_spline(1)=-sin(forTraj(2,0) + girdleCpgOutput(0));
    m1_spline(0)=-cos(forTraj(5,0) + girdleCpgOutput(1));
    m1_spline(1)=-sin(forTraj(5,0) + girdleCpgOutput(1));

    // get basic angles from spline points
    static double const scl=0.4;
    p_spline=hermiteSpline(p0_spline, p1_spline, scl*m0_spline, scl*m1_spline, t_spline);
    for(int i=0; i<Ntrunk; i++){
        p0_spline=p_spline.block(0, i+1, 2, 1) - p_spline.block(0, i  , 2, 1);
        p1_spline=p_spline.block(0, i+2, 2, 1) - p_spline.block(0, i+1, 2, 1);
        q0_trunk_from_spline(i)=copysign(           
                    SafeAcos(p0_spline.dot(p1_spline)/p0_spline.norm()/p1_spline.norm()),
                    p0_spline(0)*p1_spline(1) - p1_spline(0)*p0_spline(1)
                    );
    }


    girdleTraj(0,1)=forTraj(0,1);
    girdleTraj(1,1)=forTraj(1,1);
    girdleTraj(2,1)=forTraj(2,1) + girdleCpgOutput(0);

    // ----------------------------------------------------- OPTIMIZATION ---------------------------------------------
    // --------------- find_optimal_parameters --------------------     
    static column_vector starting_point(2), x_lower(2), x_upper(2);
    static bool is_init_opt=false;

    if(!is_init_opt){
        is_init_opt=true;
        starting_point=1,0;
        x_lower=0.7, -1;
        x_upper=1.3, 1;
    }
    dlib::find_optimal_parameters(
        0.3,                                                        //double initial_search_radius,
        0.00001,                                                     //double eps,
        200,                                                        //const unsigned int max_f_evals,
        starting_point,                                             //matrix<double,0,1>& x,
        x_lower,                                                //const matrix<double,0,1>& x_lower,
        x_upper,                                                //const matrix<double,0,1>& x_upper,
        optimizer_spineOscillations(this, old_solution)                     //const funct& f
    );
    old_solution(0)=starting_point(0);
    old_solution(1)=starting_point(1);
    q_trunk=q0_trunk_from_spline*starting_point(0) + MatrixXd::Constant(Ntrunk,1,starting_point(1));

    for(int i=0; i<q_trunk.size(); i++){
        q_trunk(i) = q_trunk(i) > constrS ? constrS : q_trunk(i);
        q_trunk(i) = q_trunk(i) <-constrS ?-constrS : q_trunk(i);
    }

}
/* Calculates global position of all the joints for a given GPS position and front girdle orientation */
/*void 
Controller :: globalKinematics(double gpsPosition[3], double rotMat[9])
{
    gpsPos(0)=gpsPosition[0];
    gpsPos(1)=-gpsPosition[2];
    gpsPos(2)=gpsPosition[1];

   rotMatFgird = fgird2ForRotMat.inverse()*compassRotMat;    

    for(int i=0; i<11; i++){
        global_joint_pos.block<3,1>(0,i) =rotMatFgird*HJs_g[i].block<3,1>(0, 3) + gpsPos;
    }
    for(int i=0; i<4; i++){
        global_joint_pos.block<3,1>(0,11+i) = rotMatFgird*HJfl_g[i].block<3,1>(0, 3) + gpsPos;
        global_joint_pos.block<3,1>(0,15+i) = rotMatFgird*HJfr_g[i].block<3,1>(0, 3) + gpsPos;
        global_joint_pos.block<3,1>(0,19+i) = rotMatFgird*HJhl_g[i].block<3,1>(0, 3) + gpsPos;
        global_joint_pos.block<3,1>(0,23+i) = rotMatFgird*HJhr_g[i].block<3,1>(0, 3) + gpsPos;
    }
    global_feet_pos.block<3,1>(0,0) = rotMatFgird*HJfl_g[4].block<3,1>(0, 3) + gpsPos;
    global_feet_pos.block<3,1>(0,1) = rotMatFgird*HJfr_g[4].block<3,1>(0, 3) + gpsPos;
    global_feet_pos.block<3,1>(0,2) = rotMatFgird*HJhl_g[4].block<3,1>(0, 3) + gpsPos;
    global_feet_pos.block<3,1>(0,3) = rotMatFgird*HJhr_g[4].block<3,1>(0, 3) + gpsPos;

    //cout<<Fgird.block<3,3>(0,0)*rotMatFgird<<endl<<endl;
    static Matrix3d tmp33;
    tmp33=rotMatFgird;
    //cout<<Fgird.block<3,3>(0,0)<<endl<<endl;
    //FoROrientation=atan2(tmp33(1,0), tmp33(0,0));

} */

/* get global feet position from webots */
void
Controller :: getFeetPosition(double *fl, double *fr, double *hl, double *hr)
{
    globalPosFL(0)=fl[0];
    globalPosFL(1)=-fl[2];
    globalPosFL(2)=fl[1];

    globalPosFR(0)=fr[0];
    globalPosFR(1)=-fr[2];
    globalPosFR(2)=fr[1];

    globalPosHL(0)=hl[0];
    globalPosHL(1)=-hl[2];
    globalPosHL(2)=hl[1];

    globalPosHR(0)=hr[0];
    globalPosHR(1)=-hr[2];
    globalPosHR(2)=hr[1];

}

/* Calculates angles from IKIN and compensates for spine_kin */
MatrixXd
Controller :: inverseKinematicsController(MatrixXd pRef)
{   
    static MatrixXd qSol = MatrixXd::Zero(4,4);

    // check if reference is outside the constraint ellipsoid
    legs_out_of_reach.setZero();
    for(int i=0; i<4; i++){
        double ellyTmp =    pow((pRef(0,i)-workspaceEllipsoidCenter(0,i)) / workspaceEllipsoid(0,i),  2) + 
                            pow((pRef(1,i)-workspaceEllipsoidCenter(1,i)) / workspaceEllipsoid(1,i),  2) + 
                            pow((pRef(2,i)-workspaceEllipsoidCenter(2,i)) / workspaceEllipsoid(2,i),  2);

        if(ellyTmp>1){
            pRef.block<3,1>(0,i) = projectOntoEllipsoid(pRef.block<3,1>(0,i), workspaceEllipsoid.block<3,1>(0,i), workspaceEllipsoidCenter.block<3,1>(0,i));
            legs_out_of_reach(i)=1;
        }
    }
    
    // ----------------------------------- SOLVE INVERSE KINEMATICS -----------------------------------------
    //double t1=get_timestamp();
    qSol.block<4,1>(0,0)=iKinQpOases(0, pRef.block<3,1>(0,0), GP.qNULL.block<4,1>(0,0), qSol.block<4,1>(0,0), lamF, MF, max_dist, ikin_maxIter, ikin_tol, constrFL, ikin_constr_penalty);
    qSol.block<4,1>(0,1)=iKinQpOases(1, pRef.block<3,1>(0,1), GP.qNULL.block<4,1>(0,1), qSol.block<4,1>(0,1), lamF, MF, max_dist, ikin_maxIter, ikin_tol, constrFR, ikin_constr_penalty);
    qSol.block<4,1>(0,2)=iKinQpOases(2, pRef.block<3,1>(0,2), GP.qNULL.block<4,1>(0,2), qSol.block<4,1>(0,2), lamH, MH, max_dist, ikin_maxIter, ikin_tol, constrHL, ikin_constr_penalty);
    qSol.block<4,1>(0,3)=iKinQpOases(3, pRef.block<3,1>(0,3), GP.qNULL.block<4,1>(0,3), qSol.block<4,1>(0,3), lamH, MH, max_dist, ikin_maxIter, ikin_tol, constrHR, ikin_constr_penalty);
    //double t2=get_timestamp();

    //cout << (t2-t1)*1000 << endl;
    return qSol;
}


/* Solves inverse kinematics for one leg by using Damped Inverse Jacobian method */
Vector4d
Controller :: iKinNullIDLS(int leg, Vector3d pref, Vector4d qref, Vector4d q0, Vector2d lam, Matrix4d M, 
                            double max_dist, int maxIter, double tol, MatrixXd constr, Vector3d constr_penalty)
{
    static MatrixXd J(3,4), tmp33(3,3), tmp44(4,4), H(4,4), D(4,4), U(3,3), V(4,4), S(3,1), Jnull(4,4), C(4,4);
    static Vector4d dqref, q, dq;
    static Vector3d dpref, p0;
    static double norm_dp;
    static bool constr_violated=false;
    static Matrix4d P=MatrixXd::Constant(4,4,constr_penalty(0)/constr_penalty(1));
    vector<Matrix4d> H_vec(5);

    constr_penalty(0)/=constr_penalty(1);
    //P=MatrixXd::Zero(1,4);
    C=MatrixXd::Zero(4,4);
  /*  for(int i=0; i<4;i++){
        P(leg,i)=constr_penalty(0);
    }*/
    // BASE CONTROLLER
    q=q0;
    for(int k=0; k<maxIter; k++){
        
        // get current position
        H_vec=legKinematics(q, leg);
        H=H_vec[0]*H_vec[1]*H_vec[2]*H_vec[3]*H_vec[4];

        p0=H.block<3, 1>(0, 3);


        J=Jacob(q, leg);

        dqref=qref-q;

        dpref=pref-p0;

        norm_dp=sqrt(dpref(0)*dpref(0)+dpref(1)*dpref(1)+dpref(2)*dpref(2));
        if(norm_dp>max_dist){
            dpref=max_dist/norm_dp*dpref;
        }
        
        if(norm_dp<tol && !constr_violated ){
            break;
        }
        
        
        tmp44 = lam(0)*Matrix4d::Identity() + J.transpose()*J;
        dq = tmp44.inverse() * (J.transpose()*dpref  - C*P.block<1,4>(leg, 0).transpose());
        q+=dq;

        // check and penalize limits
        constr_violated=false;
        C=MatrixXd::Zero(4,4);
        for(int i=0; i<4; i++){
            if(q(i)<constr(0, i)){
                C(i,i)=-1;
                constr_violated=true;

            }
            else if(q(i)>constr(1, i)){
                C(i,i)=1;
                constr_violated=true;

            }
            if(C(i,i)){
                P(leg, i)*=constr_penalty(1);
            }
            else if (P(leg, i)>constr_penalty(0)){
                P(leg, i)/=constr_penalty(1);
            }
        }
        if(!C(0,0) && !C(1,1) && !C(2,2) && !C(3,3)){
            constr_violated==false;
            tol*=constr_penalty(2);
            

        }


    }

    // NULL SPACE POSTURE CONTROLLER
    // compute SVD decomposition of Jacobian
    J=Jacob(q, leg);
    JacobiSVD<MatrixXd> svd(J, ComputeFullU | ComputeFullV);
    S=svd.singularValues();
    U=svd.matrixU();
    V=svd.matrixV();

    D=Matrix4d::Zero();
    D(3,3)=1;
    for(int i=2; i<0;i++){
        if(abs(S(i))<0.001)
            D(i, i)=1;
    }
    Jnull=V*D*V.transpose();


    tmp44=Jnull.transpose()*M.transpose()*Jnull+lam(1)*Matrix4d::Identity();
    dq=tmp44.inverse()*Jnull.transpose()*M.transpose()*dqref;

    q+=dq;
    return q;
}

/* Solves inverse kinematics for one leg by using Damped Inverse Jacobian method */
Vector4d
Controller :: iKinQpOases(int leg, Vector3d pref, Vector4d qref, Vector4d q0, Vector2d lam, Matrix4d M, 
                            double max_dist, int maxIter, double tol, MatrixXd constr, Vector3d constr_penalty)
{
    static MatrixXd J(3,4), H(4,4), D(4,4), U(3,3), V(4,4), S(3,1), Jnull(4,4);
    static MatrixXd H_qp(4,4), h_qp(4,1), lb(4,1), ub(4,1), A_qp(2,4), lbA(2,1), ubA(2,1);
    static qpOASES::real_t H_qpoases[4*4], h_qpoases[4], lb_qpoases[4], ub_qpoases[4], A_qpoases[2*4], lbA_qpoases[2], ubA_qpoases[2];
    static Vector4d dqref, dq, q0_first;
    static Vector3d dpref, p0;
    static double norm_dp;
    static vector<bool> is_init={false, false, false, false};
    vector<Matrix4d> H_vec(5);  
    q0_first = q0;
    for(int k=0; k<maxIter; k++){
        // get current position
        H_vec=legKinematics(q0, leg);
        H=H_vec[0]*H_vec[1]*H_vec[2]*H_vec[3]*H_vec[4];
        p0=H.block<3, 1>(0, 3);

        // get Jacobians
        J=Jacob(q0, leg);
        JacobiSVD<MatrixXd> svd(J, ComputeFullU | ComputeFullV);
        S=svd.singularValues();
        U=svd.matrixU();
        V=svd.matrixV();
        D=Matrix4d::Zero();
        D(3,3)=1;
        /*for(int i=0; i<3;i++){
            if(abs(S(i))<0.01)
                D(i, i)=1;
        }*/

        Jnull=V*D*V.transpose();
        // get reference velocities
        dqref=qref-q0;
        dpref=pref-p0;

        // limit maximum linear velocity
        norm_dp=sqrt(dpref(0)*dpref(0)+dpref(1)*dpref(1)+dpref(2)*dpref(2));
        if(norm_dp>max_dist){
            dpref=max_dist/norm_dp*dpref;
        }
        if(dpref.norm()<tol){
            break;
        }
        
        //=========================================== QP ================================================
        // construct qp matrices

        H_qp=J.transpose()*J + M + lam(0)*MatrixXd::Identity(4,4);
        h_qp=-J.transpose()*dpref - M*Jnull*dqref;
        lb=constr.block<1,4>(0,0).transpose() - q0;
        ub=constr.block<1,4>(1,0).transpose() - q0;


        Map<MatrixXd>( &H_qpoases[0], H_qp.cols(), H_qp.rows() ) = H_qp.transpose();
        Map<MatrixXd>( &h_qpoases[0], h_qp.rows(), h_qp.cols() ) = h_qp;
        Map<MatrixXd>( &lb_qpoases[0], lb.rows(), lb.cols() ) = lb;
        Map<MatrixXd>( &ub_qpoases[0], ub.rows(), ub.cols() ) = ub;

        A_qp << 1,1,0,0,
                1,-1,0,0;
        lbA <<  constr(0,0)*0 -13500*my_pi/180. -(q0(0)+q0(1)), 
                constr(0,1)*0 -13500*my_pi/180. -(q0(0)-q0(1));                  
        ubA <<  constr(1,0)*0 +13500*my_pi/180. -(q0(0)+q0(1)), 
                constr(1,1)*0 +13500*my_pi/180. -(q0(0)-q0(1));   
        Map<MatrixXd>( &A_qpoases[0], A_qp.cols(), A_qp.rows() ) = A_qp.transpose();
        Map<MatrixXd>( &lbA_qpoases[0], lbA.rows(), lbA.cols() ) = lbA;
        Map<MatrixXd>( &ubA_qpoases[0], ubA.rows(), ubA.cols() ) = ubA;
        // init qp solver
        qpOASES::int_t nWSR = 300;
        if(!is_init[leg]){
            //qpOASES::QProblemB tmpSolver(4);
            qpOASES::QProblem tmpSolver(4, 2);
            ikinQpSolver[leg]=tmpSolver;
            qpOASES::Options myOptions;
            //myOptions.setToMPC();
            myOptions.printLevel = qpOASES::PL_NONE;
            ikinQpSolver[leg].setOptions(myOptions);
            
        }
        
        //ikinQpSolver[leg].init( H_qpoases,h_qpoases,lb_qpoases,ub_qpoases, nWSR, 0 );
        ikinQpSolver[leg].init( H_qpoases, h_qpoases, A_qpoases, lb_qpoases, ub_qpoases, lbA_qpoases, ubA_qpoases, nWSR, 0 );
        //ikinQpSolver[leg].hotstart( H_qpoases, h_qpoases, A_qpoases, lb_qpoases, ub_qpoases, lbA_qpoases, ubA_qpoases, nWSR, 0 );

        // get solution
        qpOASES::real_t sol_qpoases[4];
        ikinQpSolver[leg].getPrimalSolution(sol_qpoases);
        dq(0)=sol_qpoases[0];
        dq(1)=sol_qpoases[1];
        dq(2)=sol_qpoases[2];
        dq(3)=sol_qpoases[3];

        q0+=dq;
        //M/=1.1;
        if(k==(maxIter-1)){
            //q0=q0_first;
        }
    }

    return q0;
}


void
Controller :: getLegJacobians(){
    legJacob[0]=Jacob(joint_angles.block<4,1>(0,0), 0);
    legJacob[1]=Jacob(joint_angles.block<4,1>(4,0), 1);
    legJacob[2]=Jacob(joint_angles.block<4,1>(8,0), 2);
    legJacob[3]=Jacob(joint_angles.block<4,1>(12,0), 3);

}

Vector3d
Controller :: getCoM()
{

    static Matrix4d mvec=MatrixXd::Identity(4,4);
    static double total_mass;

    CoM<<0,0,0;
    total_mass=0;


    //head
    mvec.block<3,1>(0,3)=masses.block<3,1>(1,0);
    CoM+=masses(0,0)*(Hhead*mvec).block<3,1>(0,3);
    total_mass+=masses(0,0);
    


    //FGIRD
    mvec.block<3,1>(0,3)=masses.block<3,1>(1,1);
    CoM+=masses(0,1)*(Fgird*mvec).block<3,1>(0,3);
    total_mass+=masses(0,1);

    //TRUNK
    mvec.block<3,1>(0,3)=masses.block<3,1>(1,2);
    CoM+=masses(0,2)*(HJs_g[0]*mvec).block<3,1>(0,3);
    total_mass+=masses(0,2);

    //HGIRD
    mvec.block<3,1>(0,3)=masses.block<3,1>(1,3);
    CoM+=masses(0,3)*(Hgird*mvec).block<3,1>(0,3);
    total_mass+=masses(0,3);

    //FL
    for(int i=0; i<Nlegs; i++){
        mvec.block<3,1>(0,3)=masses.block<3,1>(1,Ntrunk+2+i);
        CoM+=masses(0,Ntrunk+2+i)*(HJfl_g[i]*mvec).block<3,1>(0,3);
        total_mass+=masses(0,Ntrunk+2+i);
    }

    //FR
    for(int i=0; i<Nlegs; i++){
        mvec.block<3,1>(0,3)=masses.block<3,1>(1,Ntrunk+2+i+Nlegs);
        CoM+=masses(0,Ntrunk+2+i+Nlegs)*(HJfr_g[i]*mvec).block<3,1>(0,3);
        total_mass+=masses(0,Ntrunk+2+i+Nlegs);
    }

    //HL
    for(int i=0; i<Nlegs; i++){
        mvec.block<3,1>(0,3)=masses.block<3,1>(1,Ntrunk+2+i+2*Nlegs);
        CoM+=masses(0,Ntrunk+2+i+2*Nlegs)*(HJhl_g[i]*mvec).block<3,1>(0,3);
        total_mass+=masses(0,Ntrunk+2+i+2*Nlegs);
    }

    //HR
    for(int i=0; i<Nlegs; i++){
        mvec.block<3,1>(0,3)=masses.block<3,1>(1,Ntrunk+2+i+3*Nlegs);
        CoM+=masses(0,Ntrunk+2+i+3*Nlegs)*(HJhr_g[i]*mvec).block<3,1>(0,3);
        total_mass+=masses(0,Ntrunk+2+i+3*Nlegs);
    }

    //TAIL
    for(int i=0; i<Ntail; i++){
        mvec.block<3,1>(0,3)=masses.block<3,1>(1,Ntrunk+2+i+4*Nlegs);
        CoM+=masses(0,Ntrunk+2+i+4*Nlegs)*(HJt_g[i]*mvec).block<3,1>(0,3);
        total_mass+=masses(0,Ntrunk+2+i+4*Nlegs);
    }

    CoM/=total_mass;

    return CoM;
}


/* Follows the reference for COM */
Matrix<double, 3, 4>
Controller :: followCoMReference(MatrixXd comRef){
    static Vector2d comReal, comEst, comU;
    static Vector3d comTmp;
    static bool is_init=false;
    static PID PIDx(1,0,0,dt), PIDy;
    MatrixXd feetMod=MatrixXd::Zero(3,4);

    if(is_init==false){
        comU << 0,0;
        PIDy=PIDx;
        is_init=true;
    }

    // get com estimation
    comTmp=getCoM();

    /*
    if(t>0.5)
    {
        comTmp= AngleAxisd(forRPY(1), Vector3d::UnitY())*
                AngleAxisd(forRPY(0), Vector3d::UnitX())*
                comTmp;
    };
    */


    comEst=comTmp.block(0,0,2,1);
    comReal = comEst + comU*0;


    // PID
    comU(0)=PIDx.calc(comRef(0) - comReal(0));
    comU(1)=PIDy.calc(comRef(1) - comReal(1));

    //comU(0)+=(comRef(0) - comReal(0));
    //comU(1)+=(comRef(1) - comReal(1));



    feetMod << -comU(0), -comU(0), -comU(0), -comU(0),
                -comU(1), -comU(1), -comU(1), -comU(1),
                0, 0, 0, 0;
    //cout << feetMod << endl;

    feetMod.block<3,1>(0,0)=Fgird.block<3,3>(0,0).inverse()*feetMod.block<3,1>(0,0);
    feetMod.block<3,1>(0,1)=Fgird.block<3,3>(0,0).inverse()*feetMod.block<3,1>(0,1);
    feetMod.block<3,1>(0,2)=Hgird.block<3,3>(0,0).inverse()*feetMod.block<3,1>(0,2);
    feetMod.block<3,1>(0,3)=Hgird.block<3,3>(0,0).inverse()*feetMod.block<3,1>(0,3);

    /*
    for(int k=0; k<4; k++){
        if(!legs_stance(k)){
            feetMod.block<3,1>(0,k)*=0;
        }
    }*/



    /*if(LOG_DATA){
        static ofstream shiftLog("shiftLog.txt");
        shiftLog << comU(0) << "\t" << comU(1) <<endl;

    }*/

    return feetMod;


}