 #include "controller.hpp"

using namespace std;
using namespace Eigen;

extern int IS_SIMULATION, IS_PLEUROBOT, USE_JOYSTICK, JOYSTICK_TYPE, SWIM, SPINE_COMPENSATION_WITH_FEEDBACK, USE_REFLEXES;

typedef dlib::matrix<double,0,1> column_vector;

double getCircleFittingCost(MatrixXd points_input, column_vector x, Vector2d fgirdLoc, Vector2d hgirdLoc){
    Matrix<double, 2, Dynamic> p1, p2;
    p1 = points_input.block(0,0,   2,points_input.cols()/2);
    p2 = points_input.block(0,points_input.cols()/2,    2, points_input.cols()/2);

    int N=p1.cols();

    // get distances to the center
    MatrixXd R1(1,N), R2(1,N);
    for(int i=0; i<N; i++){
        R1(i) = sqrt((x(0) - p1(0,i))*(x(0) - p1(0,i)) + (x(1) - p1(1,i))*(x(1) - p1(1,i)));
        R2(i) = sqrt((x(0) - p2(0,i))*(x(0) - p2(0,i)) + (x(1) - p2(1,i))*(x(1) - p2(1,i)));
    }
    //double R = 0.5*(R1.mean() + R2.mean());
    double R = sqrt ( (fgirdLoc(0)-x(0))*(fgirdLoc(0)-x(0)) + (fgirdLoc(1)-x(1))*(fgirdLoc(1)-x(1)) );
    double Rh = sqrt ( (hgirdLoc(0)-x(0))*(hgirdLoc(0)-x(0)) + (hgirdLoc(1)-x(1))*(hgirdLoc(1)-x(1)) );

    double f=0;
    double weights[4] = {1, 0.5, 1, 0.5};
    for(int i=0; i<N; i++){
        f += weights[i]*pow(  abs(R1(i) - R) - abs(R2(i) - R) , 2 );
    }
    f += 0*pow(  Rh - R, 2 );
    //f = (    (R*MatrixXd::Ones(1,N)-R1) - (R2-R*MatrixXd::Ones(1,N))    ).norm();  

    return f;
}


class optimizer_circleFitting
{
    public:
        optimizer_circleFitting (Controller *ctrl_input, MatrixXd points_input, Vector2d fgirdLocIn, Vector2d hgirdLocIn)
        {
            ctrl=ctrl_input;
            p1 = points_input.block(0,0,   2,points_input.cols()/2);
            p2 = points_input.block(0,points_input.cols()/2,    2, points_input.cols()/2);
            fgirdLoc=fgirdLocIn;
            hgirdLoc=hgirdLocIn;
        }

        double operator() ( const column_vector& x) const
        {
            int N=p1.cols();

            // get distances to the center
            MatrixXd R1(1,N), R2(1,N);
            for(int i=0; i<N; i++){
                R1(i) = sqrt((x(0) - p1(0,i))*(x(0) - p1(0,i)) + (x(1) - p1(1,i))*(x(1) - p1(1,i)));
                R2(i) = sqrt((x(0) - p2(0,i))*(x(0) - p2(0,i)) + (x(1) - p2(1,i))*(x(1) - p2(1,i)));
            }
            //double R = 0.5*(R1.mean() + R2.mean());
            double R = sqrt ( (fgirdLoc(0)-x(0))*(fgirdLoc(0)-x(0)) + (fgirdLoc(1)-x(1))*(fgirdLoc(1)-x(1)) );
            double Rh = sqrt ( (hgirdLoc(0)-x(0))*(hgirdLoc(0)-x(0)) + (hgirdLoc(1)-x(1))*(hgirdLoc(1)-x(1)) );

            double f=0;
            double weights[4] = {1, 0.5, 1, 0.5};
            for(int i=0; i<N; i++){
                f += weights[i]*pow(  abs(R1(i) - R) - abs(R2(i) - R) , 2 );
            }
            f += 0*pow(  Rh - R, 2 );

            //double f = ((R1-R*MatrixXd::Ones(1,N)).cwiseAbs() - (R2-R*MatrixXd::Ones(1,N)).cwiseAbs()).norm();

            double e=0;
            double d;
            if(R1.mean() > R2.mean()){
                for(int i=0; i<N; i++){
                    d = R1(i) - R;
                    if(d<0)
                        e=e-d;
                    d = R-R2(i);
                    if(d<0)
                        e=e-d;
                }
            }
            else{
                for(int i=0; i<N; i++){
                    d = R2(i) - R;
                    if(d<0)
                        e=e-d;
                    d = R-R1(i);
                    if(d<0)
                        e=e-d;
                }
            }

           
            f=f+e*100; 


            return f;
        }

    private:
        Controller *ctrl;
        Matrix<double, 2, Dynamic> p1, p2;
        Vector2d fgirdLoc;
        Vector2d hgirdLoc;
};





void
Controller :: GetCompass(double data[3])
{	
	compassData[0]=data[0];  
	compassData[1]=data[1];
	compassData[2]=data[2];

}

/* Choose appropriate frequency based on the desired velocity and range of the legs */
void
Controller :: getWalkingSpeedFrequency(){
    double front_range=2*GP.ellipse_a(0);
    double hind_range=2*GP.ellipse_a(2);
    front_range*=1-abs(turning_curvature)/3.;
    hind_range*=1-abs(turning_curvature)/3.;

    if(state == CRAWLING){
        joy_lsr=1;
    }

    if(!speed_or_frequency){    // SPEED PRIORITY
        walking_forward_velocity=maxSpeed*min(joy_lsr + joy_rsr*joy_r1,1.0);
        freq_walk=walking_forward_velocity/min(front_range, hind_range)*GP.Duty(0);
    }
    else{                   // FREQUENCY PRIORITY
        freq_walk=maxFrequency*min(joy_lsr + joy_rsr*joy_r1,1.0);
        walking_forward_velocity=maxFrequency*joy_lsr*min(front_range, hind_range)/GP.Duty(0);
    }
}

/* Choose appropriate stance starting point to get a maximum possible stance length on the ellipse */
void
Controller :: getSwingStanceEndPoints()
{

    // radius and angle of side stepping 
    double front_angle=atan2(forVelocity_filtered(1,0), forVelocity_filtered(0,0));
    double hind_angle=atan2(forVelocity_filtered(1,1), forVelocity_filtered(0,1));

    if(abs(forVelocity_filtered(1,1)) + abs(forVelocity_filtered(0,1)) < 1e-3){
        hind_angle=0;
    }


    //front_angle=0;
    //hind_angle=0;
    double r_front=GP.ellipse_a(0)*GP.ellipse_b(0)/sqrt(GP.ellipse_a(0)*GP.ellipse_a(0)*sin(front_angle)*sin(front_angle)+
                                                  GP.ellipse_b(0)*GP.ellipse_b(0)*cos(front_angle)*cos(front_angle));
    double r_hind=GP.ellipse_a(2)*GP.ellipse_b(2)/sqrt(GP.ellipse_a(2)*GP.ellipse_a(2)*sin(hind_angle)*sin(hind_angle)+
                                                  GP.ellipse_b(2)*GP.ellipse_b(2)*cos(hind_angle)*cos(hind_angle));
    
    r_front*=1-abs(turning_curvature)/3.;
    r_hind*=1-abs(turning_curvature)/3.;
   

    for(int i=0; i<2; i++){
        stanceStart(0,i)=r_front*cos(front_angle)+GP.midStance(0,i);
        stanceStart(0,i+2)=r_hind*cos(hind_angle)+GP.midStance(0,i+2);
        stanceStart(1,i)=r_front*sin(front_angle)+GP.midStance(1,i);
        stanceStart(1,i+2)=r_hind*sin(hind_angle)+GP.midStance(1,i+2);
        stanceStart(2,i)=GP.midStance(2,i);
        stanceStart(2,i+2)=GP.midStance(2,i+2);
    }
    for(int i=0; i<2; i++){
        stanceEstEnd(0,i)=r_front*cos(front_angle+my_pi)+GP.midStance(0,i);
        stanceEstEnd(0,i+2)=r_hind*cos(hind_angle+my_pi)+GP.midStance(0,i+2);
        stanceEstEnd(1,i)=r_front*sin(front_angle+my_pi)+GP.midStance(1,i);
        stanceEstEnd(1,i+2)=r_hind*sin(hind_angle+my_pi)+GP.midStance(1,i+2);
        stanceEstEnd(2,i)=GP.midStance(2,i);
        stanceEstEnd(2,i+2)=GP.midStance(2,i+2);
    }


    // distort in case of turning
    if(state == WALKING){
        if(freq_walk>0.1){
            for(int i=0; i<2; i++){
                stanceStart(0,i)  +=0.5*GP.Duty(0)/freq_walk * (-forAngularVelocity_filtered(0)) * stanceStart(1,i);
                stanceStart(0,i+2)+=0.5*GP.Duty(1)/freq_walk * (-forAngularVelocity_filtered(1)) * stanceStart(1,i+2);
            } 
        }
    }
    
    /*if(state == CRAWLING){
        if(freq_walk>0.1){
            for(int i=0; i<2; i++){
                stanceStart(0,i)  +=0.05*GP.Duty(0) * (-crawling_turning_curvature) * stanceStart(1,i);
                stanceStart(0,i+2)+=0.05*GP.Duty(1) * (-crawling_turning_curvature) * stanceStart(1,i+2);
            } 
        }
    }*/


    

}

/* Updates phase of the legs */
void
Controller :: legPhaseDynamics(double dt)
{   
    
    double sumterm;
    double leg_cpg_weight=0.3;
    
    static MatrixXd legCPGtheta=MatrixXd::Ones(4,2)*0.5*2*my_pi; 

    
    double freq_walk_adaptive;
    //===================================== LEG CPG ===================================================
	for(int i=0; i<4; i++){
        if(state==STANDING && legs_stance(i)==0){
            freq_walk_adaptive = 1;
        }
        else{
            freq_walk_adaptive = freq_walk;
        }
		sumterm=0;
        for(int j=0; j<4; j++){
            sumterm+=leg_cpg_weight*(  (legCPGtheta(j,0) - GP.phShifts(j)*2*my_pi) - (legCPGtheta(i,0) - GP.phShifts(i)*2*my_pi)  )*dt;
        }
        legCPGtheta(i,1)=legCPGtheta(i,0) + 2*my_pi*freq_walk*dt + sumterm;

        legPhase(i)=fmod(legCPGtheta(i,1)/2./my_pi, 1);
        

    }
    legCPGtheta.block<4,1>(0,0)=legCPGtheta.block<4,1>(0,1);    //remember old value
    //============================== DETECT STANCE AND SWING TRANSITIONS =============================
    for(int i=0; i<4; i++){
        if(legPhase(i)<GP.Duty(i/2))
        {
            legs_stance(i)=1;
        }
        else{
            legs_stance(i)=0;
        }

        if(legs_stance_old(i)==1 && legs_stance(i)==0)
            phaseTransitionStanceSwing(i)=1;
        else
            phaseTransitionStanceSwing(i)=0;
        if(legs_stance_old(i)==0 && legs_stance(i)==1)
            phaseTransitionSwingStance(i)=1;
        else
            phaseTransitionSwingStance(i)=0;
    }
    legs_stance_old=legs_stance;
    //================================================================================================

}


/* Returns a point for a corresponding phase on the swing trajectory defined by a besiere spline passing through
initial, middle and final point */
Vector3d 
Controller :: getSwingTrajectory(Vector3d initPoint, Vector3d middlePoint, Vector3d finalPoint, double phase, int leg)
{   


    // arange time
    double tB0=0;
    double tB1=GP.tSclSwing;
    double tB2=1;
    double tB3=1;

    if(GP.tSclSwing>0){
        phase = pow(1-phase,3)*tB0 + 3*phase*pow(1-phase,2)*tB1 + 3*pow(phase,2)*(1-phase)*tB2 + pow(phase,3)*tB3;
    }



    // initialize beziere variables                                
    MatrixXd B1(3, 4), B2(3, 4);
    Vector3d Bcurve(3, 1), tmp1, tmp2;
    double tB;

    tmp1=finalPoint-initPoint;
    tmp1/=tmp1.norm();

    tmp2=-finalPoint+initPoint;
    tmp2/=tmp2.norm();

    B1.block<3, 1>(0,0)=initPoint;
    B1.block<3, 1>(0,1)=initPoint + GP.nLO.block<3,1>(0,leg)*GP.bezierScaling(0,leg/2);
    B1.block<3, 1>(0,2)=middlePoint + tmp2*GP.bezierScaling(1,leg/2);
    B1.block<3, 1>(0,3)=middlePoint;

    B2.block<3, 1>(0,0)=middlePoint;
    B2.block<3, 1>(0,1)=middlePoint + tmp1*GP.bezierScaling(1,leg/2);
    B2.block<3, 1>(0,2)=finalPoint + GP.nTD.block<3,1>(0,leg)*GP.bezierScaling(2,leg/2);
    B2.block<3, 1>(0,3)=finalPoint;






    // get swing points
    if(phase<0.5){
        tB=phase*2;

        Bcurve= pow(1-tB, 3)*B1.block<3, 1>(0,0) +
                3*pow(1-tB, 2)*tB* B1.block<3, 1>(0,1)  +
                3*pow(tB, 2)*(1-tB)* B1.block<3, 1>(0,2)  +
                pow(tB, 3)*B1.block<3, 1>(0,3); 
    }
    else{
        tB=phase*2-1;

        Bcurve= pow(1-tB, 3)*B2.block<3, 1>(0,0) +
                3*pow(1-tB, 2)*tB* B2.block<3, 1>(0,1)  +
                3*pow(tB, 2)*(1-tB)* B2.block<3, 1>(0,2)  +
                pow(tB, 3)*B2.block<3, 1>(0,3);
    }
    //velFR=(Bcurve - rFR)/dt;
    //rFR=rFR + velFR*dt;

    //legs_stance(1)=0;



    return Bcurve;

} 

/* Takes care of a swing phase of the swinging leg */
bool
Controller :: moveSwingLeg(int leg)
{
    static Vector3d trajPoint;
    static MatrixXd initPoint(3,4), middlePoint(3,4), finalPoint(3,4);
    static bool is_init=false;
    if(!is_init){
        initPoint.block<3,1>(0,0)=(Fgird.inverse()*HJfl_g[4]).block<3,1>(0,3);
        initPoint.block<3,1>(0,1)=(Fgird.inverse()*HJfr_g[4]).block<3,1>(0,3);
        initPoint.block<3,1>(0,2)=(Hgird.inverse()*HJhl_g[4]).block<3,1>(0,3);
        initPoint.block<3,1>(0,3)=(Hgird.inverse()*HJhr_g[4]).block<3,1>(0,3);
        is_init=true;
    }

    if(phaseTransitionStanceSwing(leg)){
            initPoint.block<3,1>(0,leg)=AngleAxisd(-forTraj(2+3*(leg/2),1)+girdleTraj(2+3*(leg/2),1), Vector3d::UnitZ())*feetReference.block<3,1>(0,leg);
    }

    middlePoint.block<3,1>(0,leg)=GP.midStance.block<3,1>(0,leg);
    if(leg%2){
        middlePoint(1,leg)-=GP.swing_width(leg);
    }
    else{
        middlePoint(1,leg)+=GP.swing_width(leg);
    }
    
    middlePoint(2,leg)+=GP.swing_height(leg);

    finalPoint=stanceStart;

    legs_stance(leg)=0;

    middlePoint(0,leg)=(initPoint(0,leg)+finalPoint(0,leg))/2.;
    
    swingPhase(leg)=(legPhase(leg)-GP.Duty(leg/2))/(1.-GP.Duty(leg/2));
    trajPoint=getSwingTrajectory(initPoint.block<3,1>(0,leg), middlePoint.block<3,1>(0,leg), 
                                        finalPoint.block<3,1>(0,leg), swingPhase(leg), leg);
    
    // translate to cancel the shift of hind girdle due to the spine bending 
    if(leg>1){
        trajPoint(0) += -IG-Hgird(0, 3);
    }


    feetReference.block<3,1>(0,leg) = AngleAxisd(forTraj(2+3*(leg/2),1)-girdleTraj(2+3*(leg/2),1), Vector3d::UnitZ())*trajPoint;
    
    if(swingPhase(leg)>=1){
        return true;
    }
    return false;
    
}



/* Takes care of a swing phase of the swinging leg */
bool
Controller :: moveSwingLegManualPhase(int leg, double swing_phase)
{
    static Vector3d trajPoint;
    static MatrixXd initPoint(3,4), middlePoint(3,4), finalPoint(3,4);
    static bool is_init=false;
    if(!is_init){
        initPoint.block<3,1>(0,0)=(Fgird.inverse()*HJfl_g[4]).block<3,1>(0,3);
        initPoint.block<3,1>(0,1)=(Fgird.inverse()*HJfr_g[4]).block<3,1>(0,3);
        initPoint.block<3,1>(0,2)=(Hgird.inverse()*HJhl_g[4]).block<3,1>(0,3);
        initPoint.block<3,1>(0,3)=(Hgird.inverse()*HJhr_g[4]).block<3,1>(0,3);
        is_init=true;
    }

    if(swing_phase==0){
            initPoint.block<3,1>(0,leg)=AngleAxisd(-forTraj(2+3*(leg/2),1)+girdleTraj(2+3*(leg/2),1), Vector3d::UnitZ())*feetReference.block<3,1>(0,leg);
            //initPoint.block<3,1>(0,leg) << 0, GP.midStance(1,leg), 0;
    }

    middlePoint.block<3,1>(0,leg)=GP.midStance.block<3,1>(0,leg);
    if(leg%2){
        middlePoint(1,leg)-=GP.swing_width(leg);
    }
    else{
        middlePoint(1,leg)+=GP.swing_width(leg);
    }
    
    middlePoint(2,leg)+=GP.swing_height(leg);

    finalPoint=stanceStart;

    legs_stance(leg)=0;

    middlePoint(0,leg)=(initPoint(0,leg)+finalPoint(0,leg))/2.;
    
    trajPoint=getSwingTrajectory(initPoint.block<3,1>(0,leg), middlePoint.block<3,1>(0,leg), 
                                        finalPoint.block<3,1>(0,leg), swing_phase, leg);
    
    // translate to cancel the shift of hind girdle due to the spine bending 
    if(leg>1){
        trajPoint(0) += -IG-Hgird(0, 3);
    }

    feetReference.block<3,1>(0,leg) = AngleAxisd(forTraj(2+3*(leg/2),1)-girdleTraj(2+3*(leg/2),1), Vector3d::UnitZ())*trajPoint;

    if(swing_phase>=1){
        legs_stance(leg)=1;
        return true;
    }
    return false;
    
}



/* Moves leg along arbitrary vector */
void 
Controller :: moveAlongVector(int leg, Vector3d nVec, double vel)
{
    feetReference.block<3,1>(0,leg) += vel*nVec*dt;
}


/* Path planning for crawling through pipes */
void 
Controller :: pipePathPlanner()
{
    static PIDvec pidPath;
    bool isUpdated=false;
    static column_vector starting_point(2), x_lower(2), x_upper(2);
    static bool is_init=false;
    if(!is_init){
        is_init=true;
        //pidPath.init(0.4, 0.0, 0.05, 0.05, 1, dt);
        pidPath.init(1.2, 0.0, 0.05, 0.05, 1, dt);
        x_lower=-300, -300;
        x_upper=300, 300;

        
    }

    // --------------- FIT THE CIRCLE --------------------
    

    static MatrixXd points(2,8), points_old(2,8);
    static Vector2d xc;
    static double R;
    Vector2d fg = forTraj.block<2,1>(0,1);
    static Vector2d fg_old; fg_old << 0,0;
    Vector2d hg = forTraj.block<2,1>(3,1);
    // left
    points.block<2,1>(0,0) = feetTrajHist[feetTrajHist.size()-1].block<2,1>(0,0);
    points.block<2,1>(0,1) = feetTrajHist[feetTrajHist.size()-2].block<2,1>(0,0);
    points.block<2,1>(0,2) = feetTrajHist[feetTrajHist.size()-1].block<2,1>(0,2);
    points.block<2,1>(0,3) = feetTrajHist[feetTrajHist.size()-2].block<2,1>(0,2);
    // right
    points.block<2,1>(0,4) = feetTrajHist[feetTrajHist.size()-1].block<2,1>(0,1);
    points.block<2,1>(0,5) = feetTrajHist[feetTrajHist.size()-2].block<2,1>(0,1);
    points.block<2,1>(0,6) = feetTrajHist[feetTrajHist.size()-1].block<2,1>(0,3);
    points.block<2,1>(0,7) = feetTrajHist[feetTrajHist.size()-2].block<2,1>(0,3);

    if((points-points_old).norm()>0.00001){
        isUpdated = true;
        MatrixXd startingPoints(4, 2);
        startingPoints <<    50.0 + fg(0),  50.0 + fg(1), 
                            -50.0 + fg(0),  50.0 + fg(1), 
                            -50.0 + fg(0), -50.0 + fg(1), 
                             50.0 + fg(0), -50.0 + fg(1);
        startingPoints*=0.5;
        double cost = 0, cost_old = 99999999999999999;

        for(int i=0; i<4; i++){
            starting_point=startingPoints(i,0), startingPoints(i,1);

            dlib::find_min_using_approximate_derivatives(dlib::lbfgs_search_strategy(10),
                                           dlib::objective_delta_stop_strategy(1e-7),
                                           optimizer_circleFitting(this, points, fg, hg), starting_point, -1);

            cost = getCircleFittingCost(points, starting_point, fg, hg);
            cout << cost << "\t";
            if(cost < cost_old){
                cost_old = cost;
                xc << starting_point(0), starting_point(1);
            }
        }
        cout << endl;
        

        MatrixXd R1(1,4), R2(1,4);
        for(int i=0; i<4; i++){
            R1(i) = sqrt((xc(0) - points(0,i))*(xc(0) - points(0,i)) + (xc(1) - points(1,i))*(xc(1) - points(1,i)));
            R2(i) = sqrt((xc(0) - points(0,4+i))*(xc(0) - points(0,4+i)) + (xc(1) - points(1,4+i))*(xc(1) - points(1,4+i)));
        }
        //R = 0.5*(R1.mean() + R2.mean());
        R = sqrt ( (fg(0)-xc(0))*(fg(0)-xc(0)) + (fg(1)-xc(1))*(fg(1)-xc(1)) );
        static ofstream circfitlog("/data/thorvat/webotsLogs/circfitlog.txt");
        circfitlog << xc.transpose() << "\t" << R << "\t"<< points.block<1,8>(0,0) << "\t" << points.block<1,8>(1,0) <<endl;
        cout << "xc\t" << xc.transpose() << "\t\t R\t" << R << endl;
        cout << endl <<points.block<1,8>(0,0) << "\t" << points.block<1,8>(1,0) << endl << endl;
    }
    points_old = points;


    // DEFINE THE PATH 
    double lhdist = 0.2;
    

    // find closest point to the front girdle on the center line
    double phifg = atan2(fg(1) - xc(1), fg(0) - xc(0));
    //double phihg = atan2(hg(1) - xc(1), hg(0) - xc(0));

    Vector2d plh1, plh2, v1, v2;
    plh1 << R*cos(phifg+0.01) + xc(0) , R*sin(phifg+0.01) + xc(1);
    plh2 << R*cos(phifg-0.01) + xc(0) , R*sin(phifg-0.01) + xc(1);

    v1 = fg-hg;
    v2 = plh1 - plh2;
    // look ahead

    double philh;
    if( v1.dot(v2) > 0 ){ // POSSIBLE PROBLEM WITH WRAPPING (plh1-fg_old).norm() > (plh2-fg_old).norm()
        philh = phifg + lhdist / R;
    }
    else{
        philh = phifg - lhdist / R;
    }


    Vector2d plh;
    plh << R*cos(philh) + xc(0) , R*sin(philh) + xc(1);

    static double angref_old=0;
    static int Npi=0;
    if(isUpdated){
        Npi = 0;
    }
    double angref = atan2(plh(1)-fg(1), plh(0)-fg(0));
    

    

    /*if(angref-angref_old> 6){
        Npi-=2;
        angref+=-2*my_pi;
    }
    else if(angref - angref_old< -6){
        Npi+=2;
        angref+= 2*my_pi;
    }

    if(angref-angref_old> 3){
        Npi--;
        angref+=-my_pi;
    }
    else if(angref - angref_old< -3){
        Npi++;
        angref+= my_pi;
    }*/

    

    angref_old = angref;

    
    double e = angref - forTraj(2,1);
    Npi = round(e / my_pi); 
    e = e - Npi*my_pi;
    static double e_old = 0;
    e = pt1(e, e_old, 0.3, dt);
    e_old=e;

    walking_angular_velocity = pidPath.calc(e);
    

    /*double wlks = walking_angular_velocity;
    v1 << cos(forTraj(2,1)), sin(forTraj(2,1));
    v2 = xc-fg;

    if( v1(0)*v2(1)-v1(1)*v2(0) >= 0){
        walking_angular_velocity = walking_forward_velocity / R ;  // sign info is missing
    }
    else{
        walking_angular_velocity = -walking_forward_velocity / R ;  // sign info is missing
    }
    cout << wlks << "\t\t||\t\t" << walking_angular_velocity << endl;*/

    
    crawling_turning_curvature = walking_angular_velocity / walking_forward_velocity;
    static ofstream trajpidlog("/data/thorvat/webotsLogs/trajpidlog.txt");
    trajpidlog << angref<< "\t" << e << "\t" << walking_angular_velocity << "\t" << forTraj.block<6,1>(0,1).transpose() << "\t" <<plh.transpose() <<"\t" << fg.transpose()<< endl;
     

     fg_old=fg;
    

}

MatrixXd 
Controller :: standingTransition(MatrixXd pRef)
{
    
    static MatrixXd pRef_old(3,4);
    static bool is_init=false;

    if(!is_init){
        is_init=true;
        pRef_old=pRef;
        return pRef;
    }


    // average height of stance legs
    double stanceHeight=0;
    int cnt=0;

    if(state == STANDING || state == POSING){
        for(int i=0; i<4; i++){
            if(legs_stance(i)){
                stanceHeight+=pRef(2,i);
                cnt++;
            }
        }
        if(cnt>0){
            stanceHeight/=cnt;
        }
        stanceHeight = GP.midStance(2,0);

        // get new reference
        for(int i=0; i<4; i++){
            if(legs_stance(i)==0){
                pRef(2,i)=stanceHeight;
            }
        }
    }
    if(state==POSING){
        pRef(2,0)=posingFeetHeight;
        pRef(2,1)=posingFeetHeight;
    }
    // filter
    pRef=pt1_vec(pRef, pRef_old, T_stand, dt);
    pRef_old=pRef;

    return pRef;


}


