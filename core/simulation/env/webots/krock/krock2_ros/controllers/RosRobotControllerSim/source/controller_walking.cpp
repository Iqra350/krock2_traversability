#include "controller.hpp"

using namespace std;
using namespace Eigen;

extern int IS_SIMULATION, IS_PLEUROBOT, USE_JOYSTICK, JOYSTICK_TYPE, SWIM, SPINE_COMPENSATION_WITH_FEEDBACK, USE_REFLEXES, LOG_DATA, USE_TAIL;



// ---------------------------------------------- SOLVE MPC THREAD ---------------------------------------------------------

void 
solveMPC(MPC *mpc_com, MatrixXd *sol)
{
 	double t1, t2, t3, dt;
 	//static ofstream mpcTimeLog("./data/mpcTimeLog.txt");
	while(true){
		t1=get_timestamp();


		// update QP arguments with mutex protection
		mpc_in_mtx.lock();
		mpc_com->updateQpArg();
		mpc_in_mtx.unlock();

		// solve QP
		//mpc_com->printMatrices("MPC_TEST.txt");

		// sync
		//mpc_sync_mtx.lock();
		mpc_com->solveQP();

		
		//mpc_sync_mtx.unlock();

		// update results with mutex protection
		mpc_out_mtx.lock();
		*sol = mpc_com->sol;
		mpc_out_mtx.unlock();






		t2=get_timestamp();

		
		if(mpc_com->dt > (t2-t1)){
			usleep((mpc_com->dt - (t2-t1))*1000000);
		}


		//cout << "MPC time: " << (t2-t1)*1000 <<" ms"<<endl;
		//mpcTimeLog << (t2-t1)*1000 << endl;
		

	}	
}




// ---------------------------------------------- FUNCTIONS ----------------------------------------------------------

/* handles walking and all subcontrollers */
void 
Controller :: walkingStateMachine()
{
	// get leg phases (accounts for GP.Duty cycles and leg phase offsets)
	legPhaseDynamics(dt);

	// move swing leg
	for(int i=0; i<4; i++){
		if(legs_stance(i)==0){
			moveSwingLeg(i);
		}
	}

	// get girdle trajectories  
	// double v, double w, Vector3d V, double W, double fgirdArc, Vector2d spineCPGscaling
	girdleTrajectories(walking_forward_velocity, walking_angular_velocity, MatrixXd::Zero(3,1), walking_rotational_velocity, 0, GP.spineCPGscaling);
	
	// udpate girdle velocities
	girdleVelocities();



	// RUN MPC
	//runComMPC();
}

/* Translates and rotates a girdle by using a correspoinging leg pair while in stance */
void 
Controller :: moveGirdle(Vector3d girdleVelocity, double girdleAngularVelocity, int girdleNum)
{
	static Transform<double,3,Affine> Girdle_transformation;
	//Girdle_transformation = AngleAxisd(-girdleAngularVelocity*dt, Vector3d::UnitZ()) * Translation3d(-girdleVelocity(0)*dt, -girdleVelocity(1)*dt, -girdleVelocity(2)*dt);
	Girdle_transformation = AngleAxisd(-girdleAngularVelocity*dt, Vector3d::UnitZ()) * Translation3d(-girdleVelocity*dt);

	for(int i=girdleNum*2; i<2*(1+girdleNum); i++){
		//if(legs_stance(i)){
			feetReference.block<3,1>(0,i) = Girdle_transformation*feetReference.block<3,1>(0,i);
		//}
	}

}

/* calculates girdle trajctories / velocities for given input commands - linear and angular velocity of the front girdle */
void
Controller :: girdleTrajectories(double v, double w, Vector3d V, double W, double fgirdArc, Vector2d spineCPGscaling)
{	

	// INITIALIZE
	MatrixXd tmp_traj(3,300);

	static double mindisp=0.003333333333;

	// ------------------------------------------ REMEMBERING OLD STATES -------------------------------
	forTraj.block<6,1>(0,0)=forTraj.block<6,1>(0,1);
    girdleTraj.block<6,1>(0,0)=girdleTraj.block<6,1>(0,1);

    // ------------------------------------------ FORWARD LOCOMOTION -------------------------------
    // move frame of reference
    forTraj(0,1)=forTraj(0,0)+v*cos(forTraj(2,0))*dt;
    forTraj(1,1)=forTraj(1,0)+v*sin(forTraj(2,0))*dt;
    forTraj(2,1)=forTraj(2,0)+w*dt;

    // remember FOR trajectory - update if displacement is greater than mindelta
    static double disp, hist_disp;
    if(state!=POSING){
	    disp =sqrt((forTraj(0,1)-forTrajHist(0,0))*(forTraj(0,1)-forTrajHist(0,0))+(forTraj(1,1)-forTrajHist(1,0))*(forTraj(1,1)-forTrajHist(1,0)));
	    hist_disp =sqrt((forTraj(0,1)-forTrajHist(0,299))*(forTraj(0,1)-forTrajHist(0,299))+(forTraj(1,1)-forTrajHist(1,299))*(forTraj(1,1)-forTrajHist(1,299)));
	    if(disp>mindisp && hist_disp>1.2*IG)
	    {
	    	tmp_traj=forTrajHist.block<3,299>(0,0);
	        forTrajHist.block<3,299>(0,1)=tmp_traj;
	    	forTrajHist.block<3,1>(0,0)=forTraj.block<3,1>(0,0);
	    }
	}
 	


    

 	// -------------------------------------- FRONT GIRDLE ARC -----------------------------------------
	// get the centerpoint
	Vector2d centerPoint;

	forTraj(2,1) += 2*fgirdArc*dt;
	girdleTraj(2,1) += 2*fgirdArc*dt;
	
	

	centerPoint << forTraj(3,1), forTraj(4,1);
	Affine2d A = Translation<double,2>(centerPoint)* Rotation2D<double>(fgirdArc*dt)  * Translation<double,2>(-centerPoint);
	forTraj.block<2,1>(0,1) = A * forTraj.block<2,1>(0,1);
	girdleTraj.block<2,1>(0,1) = A * girdleTraj.block<2,1>(0,1);
	

    


	// -------------------------------------- GIRDLE OSCILLATIONS -----------------------------------------
    girdleOscillations(spineCPGscaling);

    // -------------------------------------- SPINE INVERSE KINEMATICS -----------------------------------------
    trunkInverseKinematics();

    // --------------- update kinematics -------------------- 
    girdleTraj.block<3,1>(3,1)=trunkForwardKinematics(girdleTraj.block<3,1>(0,1), q_trunk);
    // find minimal distance of the hind girdle form the stored trajectory
	double mindist=9999, dist; 
	static int indx;
	indx=0;
	for(int i=50; i<300; i++){
		dist=sqrt((forTrajHist(0,i)-girdleTraj(3))*(forTrajHist(0,i)-girdleTraj(3)) + 
				  (forTrajHist(1,i)-girdleTraj(4))*(forTrajHist(1,i)-girdleTraj(4)));
		if(dist<mindist){
			mindist=dist;
			indx=i;
		}
	}

    // get a new hgird position and orientation
    forTraj(3,1)=girdleTraj(3,1);
    forTraj(4,1)=girdleTraj(4,1);
    forTraj(5,1)=forTrajHist(2,indx);


    // -------------------------------------- WHOLE BODY TRANSLATION --------------------------------------
    // get centerline from forTraj
    double centerLine = 0.5*(forTraj(2,1) + forTraj(5,1));

    // rotate whole body translation vector to align with the centerline
    V = AngleAxisd(centerLine, Vector3d::UnitZ()) * V;

    // move girdles and frames of reference
    forTraj(0,1) +=V(0)*dt;
    forTraj(1,1) +=V(1)*dt;
    forTraj(3,1) +=V(0)*dt;
    forTraj(4,1) +=V(1)*dt;

    girdleTraj(0,1) +=V(0)*dt;
    girdleTraj(1,1) +=V(1)*dt;
    girdleTraj(3,1) +=V(0)*dt;
    girdleTraj(4,1) +=V(1)*dt;

    // update history
    if(state!=POSING){
	    for(int i=0; i<300; i++){
	    	forTrajHist(0,i) +=V(0)*dt;
		    forTrajHist(1,i) +=V(1)*dt;
		}
	}

    // -------------------------------------- WHOLE BODY ROTATION -----------------------------------------
	// get the centerpoint
	centerPoint << 0.5*(forTraj(0,1) + forTraj(3,1)), 0.5*(forTraj(1,1) + forTraj(4,1));
	A = Translation<double,2>(centerPoint)* Rotation2D<double>(W*dt)  * Translation<double,2>(-centerPoint);

	forTraj.block<2,1>(0,1) = A * forTraj.block<2,1>(0,1);
	forTraj.block<2,1>(3,1) = A * forTraj.block<2,1>(3,1);
	forTraj(2,1) += W*dt;
	forTraj(5,1) += W*dt;
	girdleTraj.block<2,1>(0,1) = A * girdleTraj.block<2,1>(0,1);
	girdleTraj.block<2,1>(3,1) = A * girdleTraj.block<2,1>(3,1);
	girdleTraj(2,1) += W*dt;
	girdleTraj(5,1) += W*dt;

	// update history
    for(int i=0; i<300; i++){
    	forTrajHist.block<2,1>(0,i) = A*forTrajHist.block<2,1>(0,i);
    	forTrajHist(2,i) += W*dt;
	}



    // neck and tail
    if(Nneck>0){
    	q_neck=-(forTraj(2,1)-girdleTraj(2,1));
    }
    if(Ntail>0){
    	for(int i=0; i<2; i++){
			q_tail(i)=1.5*(forTraj(5,1)-girdleTraj(5,1))/2;
		}
    }
    

	// -------------------------------------- TRACK FOOTSTEPS ODOMETRY-----------------------------------------
    // link 1-frame with 2-frame representation
    double hgirdlePhiOdo;
    hgirdlePhiOdo = atan2(forTraj(1,1)-forTraj(4,1), forTraj(0,1)-forTraj(3,1));

    for(int i=0; i<4; i++){
        feetLocationsOdo.block<3,1>(0,i) = AngleAxisd(hgirdlePhiOdo, Vector3d::UnitZ())*feetLocations.block<3,1>(0,i);
        feetLocationsOdo(0,i)+=forTraj(0,1);
        feetLocationsOdo(1,i)+=forTraj(1,1);
    }


	if(LOG_DATA){
		static ofstream girdleTrajLog("/data/thorvat/webotsLogs/girdleTrajLog.txt");
		//girdleTrajLog << forTraj.block<6,1>(0,0).transpose() << "\t" <<forTrajHist(0,0) <<"\t"<<forTrajHist(0,299) <<"\t" << mindist<<endl;
		girdleTrajLog << girdleTraj.block<6,1>(0,0).transpose()<<"\t"<<girdleTraj.block<6,1>(0,1).transpose()<<"\t";
		girdleTrajLog << forTraj.block<6,1>(0,0).transpose()<<"\t"<<forTraj.block<6,1>(0,1).transpose()<<"\t";

		for(int i=0; i<4; i++){
			girdleTrajLog << feetLocationsOdo.block<3,1>(0,i).transpose()<<"\t";
		}
		girdleTrajLog << endl;
	}
	
}


/* predicts future trajectories over prediction horizon for MPC or anything else*/
std::vector<Matrix<double,3,4>> 
Controller :: predictTrajectories(int N, double time_step, MatrixXd *predLegPhase_in)
{

	std::vector<Matrix<double,3,4>> predictedFootsteps;
	predictedFootsteps.resize(N);

	double omega=2*my_pi*freq_walk;
	MatrixXd predLegPhase(N,4);
	predLegPhase.block(0,0,1,4)=legPhase.transpose();

	// initial (current) points (all in girdle frame reference)
	predictedFootsteps[0].block<3,1>(0,0)=AngleAxisd(-forTraj(2,1)+girdleTraj(2,1), Vector3d::UnitZ())*feetReference.block<3,1>(0,0);
	predictedFootsteps[0].block<3,1>(0,1)=AngleAxisd(-forTraj(2,1)+girdleTraj(2,1), Vector3d::UnitZ())*feetReference.block<3,1>(0,1);
	predictedFootsteps[0].block<3,1>(0,2)=AngleAxisd(-forTraj(5,1)+girdleTraj(5,1), Vector3d::UnitZ())*feetReference.block<3,1>(0,2);
	predictedFootsteps[0].block<3,1>(0,3)=AngleAxisd(-forTraj(5,1)+girdleTraj(5,1), Vector3d::UnitZ())*feetReference.block<3,1>(0,3);

	// starting point for stance phase
	MatrixXd stanceStartingPoint(3,4);

	stanceStartingPoint=stanceStart;
	Transform<double,3,Affine> Girdle_transformation;


	for(int j=0; j<4; j++){
		if(predLegPhase(0,j)>GP.Duty(j/2)){
			predictedFootsteps[0].block<3,1>(0,j) << -9999,
							                   -9999,
			   				                   -9999;
		}
	}


	
	// LOOK INTO THE FUTURE
	for(int i=1; i<N; i++){

		// run legPhaseDynamics and determine if swing or stance
		predLegPhase.block(i,0,1,4)=predLegPhase.block(i-1,0,1,4)+freq_walk*time_step*MatrixXd::Ones(1,4);
		for(int j=0; j<4; j++){
			predLegPhase(i,j)=predLegPhase(i,j)>1 ? (0+(predLegPhase(i,j)-1)) : predLegPhase(i,j);
		}

		// when in stance move the leg with current speed profiles
		for(int j=0; j<4; j++){
			
			// ignore if in swing phase										
			if(predLegPhase(i,j)>GP.Duty(j/2)){
				predictedFootsteps[i].block<3,1>(0,j) << -9999,
								                   -9999,
				   				                   -9999;
               	continue;
			}
			// if comming back from the swing phase (landing)
			if(predictedFootsteps[i-1](2,j)==-9999){
				predictedFootsteps[i].block<3,1>(0,j)=stanceStartingPoint.block<3,1>(0,j);
				continue;
			}

			// otherwise, move leg with a current speed
			Girdle_transformation = AngleAxisd(-forAngularVelocity_filtered(j/2)*time_step, Vector3d::UnitZ()) * 
												Translation3d((-forVelocity_filtered.block<3,1>(0,j/2))*time_step);

			predictedFootsteps[i].block<3,1>(0,j)=Girdle_transformation*predictedFootsteps[i-1].block<3,1>(0,j);
			
		}
		
	}




	*predLegPhase_in=predLegPhase;
	return predictedFootsteps;
}

/* MPC calls */
void 
Controller :: runComMPC()
{

	// initialize system
	static bool is_init=false;
	static VectorXd x(mpc_com.n), u(mpc_com.m);
	static Vector3d com;
	static MatrixXd A=mpc_com.A, B=mpc_com.B;

	



	// make prediction for future leg trajectories (N timesteps)
	std::vector<Matrix<double,3,4>> predictedFootsteps;
	MatrixXd predLegPhase(mpc_com.N, 4);
	predictedFootsteps=predictTrajectories(mpc_com.N, mpc_com.dt, &predLegPhase);

	static ofstream predictedFootstepsLog("supportPolysLog.txt");
	// put everyting into robot's coordinate frame
	Vector4d tmp_vector;
	for(int i=0; i<mpc_com.N; i++){
		for(int j=0; j<4; j++){
			tmp_vector << predictedFootsteps[i].block<3,1>(0,j), 1;
			if(j<2){
				tmp_vector.block<3,1>(0,0)=AngleAxisd(forTraj(2,1)-girdleTraj(2,1), Vector3d::UnitZ())*tmp_vector.block<3,1>(0,0); // girdle frame to girdle
				tmp_vector=Fgird*tmp_vector;
				//tmp_vector=tmp_vector+Fgird.block<4,1>(0,3);
			}
			else{
				tmp_vector.block<3,1>(0,0)=AngleAxisd(forTraj(5,1)-girdleTraj(5,1), Vector3d::UnitZ())*tmp_vector.block<3,1>(0,0); // girdle frame to girdle
				tmp_vector=Hgird*tmp_vector;
				//tmp_vector=tmp_vector+Hgird.block<4,1>(0,3);
			}
			supportPolys[i].block<3,1>(0,j)=tmp_vector.block<3,1>(0,0);
			supportPolysIMU[i].block<3,1>(0,j)=	AngleAxisd(forRPY(1), Vector3d::UnitY())*
												AngleAxisd(forRPY(0), Vector3d::UnitX())*
												supportPolys[i].block<3,1>(0,j);
		}
		//if(LOG_DATA){
		//	predictedFootstepsLog << supportPolys[i] << endl;
		//}
	}

	// get future reference from predicted footsteps
	MatrixXd xref(mpc_com.N*6,1);
	Vector2d polyCenter;
	int cnt;
	for(int i=0; i<mpc_com.N; i++){
		polyCenter << 0,0;
		cnt=0;
		for(int j=0; j<4; j++){
			if(supportPolys[i](0,j)>-5){
				polyCenter = polyCenter + supportPolys[i].block(0,j,2,1);
				cnt++;
			}
		}
		polyCenter=polyCenter/cnt;
		xref.block(i*6, 0, 6, 1) << polyCenter(0)*0-0.212,
									0,
									0,
									polyCenter(1)*0,
									0,
									0;

	}


	static MatrixXd sol_container=MatrixXd::Zero(mpc_com.N*(mpc_com.n+mpc_com.m+mpc_com.n_soft),1);
	if(!is_init){
		is_init=true;
		com=getCoM();
		x << com(0), 0, 0, com(1), 0, 0;
		sol_container.block(0,0,6,1)=x;
		mpc_com.updateConstraints_SupportPolygons(supportPolys);
		mpc_com.updateReference(xref);
		mpc_com.updateState(x);
	}
	
	


	// ----------------------------------- UPDATE constraints, reference and state with MUTEX protection -------------------------------
	// rotate weights
	MatrixXd Qx_new(mpc_com.n, mpc_com.n), Qf_new(mpc_com.n, mpc_com.n), Ru_new(mpc_com.m, mpc_com.m);
	Qx_new=mpc_com.Qx0;
	Qf_new=mpc_com.Qf0;
	Ru_new=mpc_com.Ru0;

	// rotate velocities
	double mphi=SafeAcos(abs(cos(walkingDirection)))*2/my_pi;
	Qx_new(1,1)=(1-mphi)*mpc_com.Qx0(1,1) + (  mphi)*mpc_com.Qx0(3,3);
	Qx_new(3,3)=(  mphi)*mpc_com.Qx0(1,1) + (1-mphi)*mpc_com.Qx0(3,3);
	Qf_new(1,1)=(1-mphi)*mpc_com.Qf0(1,1) + (  mphi)*mpc_com.Qf0(3,3);
	Qf_new(3,3)=(  mphi)*mpc_com.Qf0(1,1) + (1-mphi)*mpc_com.Qf0(3,3);
	// rotate accelerations
	Qx_new(2,2)=(1-mphi)*mpc_com.Qx0(2,2) + (  mphi)*mpc_com.Qx0(4,4);
	Qx_new(4,4)=(  mphi)*mpc_com.Qx0(2,2) + (1-mphi)*mpc_com.Qx0(4,4);
	Qf_new(2,2)=(1-mphi)*mpc_com.Qf0(2,2) + (  mphi)*mpc_com.Qf0(4,4);
	Qf_new(4,4)=(  mphi)*mpc_com.Qf0(2,2) + (1-mphi)*mpc_com.Qf0(4,4);


	mpc_in_mtx.lock();
		//updateConstraints_Box(x_min_max_in, u_min_max_in);
		mpc_com.updateWeights(Qx_new, Qf_new, Ru_new);
		mpc_com.updateConstraints_SupportPolygons(supportPolys);
		mpc_com.updateReference(xref);
		mpc_com.updateState(x);
	mpc_in_mtx.unlock();

	
	// ----------------------------------- start thread -------------------------------
	static bool mpc_thread_started=false;
	
	static MatrixXd xsol(mpc_com.N*(mpc_com.n),1);
	static MatrixXd usol(mpc_com.N*(mpc_com.m),1);
	static MatrixXd esol(mpc_com.N*(mpc_com.n_soft),1);
	if(!mpc_thread_started && t>3){
		
		if(mpc_com.enabled){
			cout << "STARTING THREAD" << endl;
			mpc_thread=std::thread(solveMPC, &mpc_com, &sol_container);
			mpc_thread.detach();
		}
		
		mpc_thread_started=true;

	}


	// ----------------------------------- get solutions with mutex protection -------------------------------
	mpc_out_mtx.lock();
		xsol=sol_container.block(0,0,mpc_com.N*mpc_com.n,1);
		usol=sol_container.block(mpc_com.N*mpc_com.n,0,mpc_com.N*mpc_com.m,1);
		esol=sol_container.block(mpc_com.N*(mpc_com.n+mpc_com.m),0,mpc_com.N*mpc_com.n_soft,1);
	mpc_out_mtx.unlock();



	// ----------------------------------- RUN SYSTEM ---------------------------------------
	//sync
	static int syncCounter=0;
	syncCounter++;
	syncCounter%=5;
	//if(!syncCounter){
	//	mpc_sync_mtx.unlock();
	//}
	u=usol.block(0,0,mpc_com.m, 1)*dt/mpc_com.dt;
	x=xsol.block(0,0,mpc_com.n, 1);
	//cout << u << endl << endl;
	//x=A*x+B*u;
	//cout << "u: "<<usol.block(0,0,mpc_com.m, 1).transpose() << endl << endl;

	//	cout << xsol.transpose() <<endl; //<<"\t"<<"\t"<<"\t"<<u.transpose()<< endl;


	// ----------------------------------- LOGGING  ---------------------------------------
	//if(LOG_DATA){
	//	static ofstream mpcReferenceLog("mpcReferenceLog.txt");
	//	mpcReferenceLog << xref.transpose() << endl;
	//	static ofstream mpc_solution("mpc_solution.txt");
	//	mpc_solution << xsol.transpose() << "\t" << usol.transpose() << "\t" << esol.transpose() << endl;
	//}
	mpcComRef << x(0), x(3);
	//mpcComRef << -0.22, 0;


}


/* Calculate girdle oscillations to follow leg movements */
void
Controller :: girdleOscillations(Vector2d spineCPGscaling)
{


	int N=60;
	double pred_dt=0.05;
	// predict trajectories
	std::vector<Matrix<double,3,4>> predictedFootsteps;
	predictedFootsteps.resize(N);

	static MatrixXd predLegPhase(N,4);
	predictedFootsteps=predictTrajectories(N, pred_dt, &predLegPhase);
	static MatrixXd scaled_phases(4,N);
	scaled_phases=predLegPhase.transpose();


	//------------------------------- MESSY FOOTSTEP ANGLE PREDICTION -------------------------------
	double footX, footY, tmp;
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<4; j++)
		{	
			if(predLegPhase(i,j)<0)
					predLegPhase(i,j)+=1;

			if(predLegPhase(i, j)<GP.Duty(j/2))
			{
				tmp=(predLegPhase(i,j))/(GP.Duty(j/2));	

				footX=(1-tmp)*stanceStart(0,0) + tmp*stanceEstEnd(0,0);
				footY=(1-tmp)*stanceStart(1,0) + tmp*stanceEstEnd(1,0);
			}
			else
			{
				
				tmp=(predLegPhase(i,j)-GP.Duty(j/2))/(1-GP.Duty(j/2));	  

				footX=(1-tmp)*stanceEstEnd(0,0) + tmp*stanceStart(0,0);
				footY=(1-tmp)*stanceEstEnd(1,0) + tmp*stanceStart(1,0);
			}

			scaled_phases(j,i)=atan2(footX, footY);
		}
	}

	scaled_phases.block(0,0,1,N)*=-1;
	scaled_phases.block(2,0,1,N)*=-1;



	// get reference for both girdles
	static MatrixXd girdleRefFromLegs(2,N);
	girdleRefFromLegs.block(0,0,1,N) = spineCPGscaling(0)*(scaled_phases.block(0,0,1,N) + scaled_phases.block(1,0,1,N));
	girdleRefFromLegs.block(1,0,1,N) = spineCPGscaling(1)*(scaled_phases.block(2,0,1,N) + scaled_phases.block(3,0,1,N));



	// init cpg
	double a=1;
	Vector2d R=girdleRefFromLegs.rowwise().maxCoeff();
	static MatrixXd theta=MatrixXd::Zero(2,2), r=MatrixXd::Zero(2,2);
	Vector2d dr, dtheta;



	//----------------------------- DFT PHASE ESTIMATION ------------------------------
	double re, im;
	static MatrixXd phase_est=MatrixXd::Zero(2,2);
	double phase_diff;
	static Vector2d phase_correction;
	for(int i=0; i<2; i++){
		re=0; im=0;
		for(int k=0; k<N; k++){
			re = re + girdleRefFromLegs(i,k)*cos(2*my_pi*freq_walk*k*pred_dt);
			im = im - girdleRefFromLegs(i,k)*sin(2*my_pi*freq_walk*k*pred_dt);
		}
		phase_est(i, 0)=phase_est(i, 1);
		phase_est(i, 1)=atan2(im, re);


		phase_correction(i)=phase_est(i,0) + round((theta(i,0) - phase_est(i,0))/(2*my_pi))*2*my_pi   - theta(i,0);
	}
	
	

	//------------------------------------- RUN CPG -------------------------------------

	dtheta = 2*my_pi*freq_walk*MatrixXd::Ones(2,1) + 2*phase_correction;

	// update radius
	R=R*(1-abs(turning_curvature))*abs(cos(walkingDirection));
	dr = a*(R-r.block(0,0,2,1));

	// euler integration
	theta.block(0,1,2,1)=theta.block(0,0,2,1) + dt*dtheta;
	r.block(0,1,2,1)=r.block(0,0,2,1) + dt*dr;

	// angle output
	girdleCpgOutput(0)=r(0,0)*cos(theta(0,0));
	girdleCpgOutput(1)=r(1,0)*cos(theta(1,0));

	// udpate old values
	theta.block(0,0,2,1)=theta.block(0,1,2,1);
	r.block(0,0,2,1)=r.block(0,1,2,1);


	return;

}

/* calculates velocities of both girdles needed in inverse kinematics */
void 
Controller :: girdleVelocities()
{

	// frame of reference velocities
	forVelocity(0,0)=(forTraj(0,1)-forTraj(0,0))/dt;
    forVelocity(1,0)=(forTraj(1,1)-forTraj(1,0))/dt;
    forVelocity(2,0)=0;
    forVelocity(0,1)=(forTraj(3,1)-forTraj(3,0))/dt;
    forVelocity(1,1)=(forTraj(4,1)-forTraj(4,0))/dt;
    forVelocity(2,1)=0;
    forAngularVelocity(0)=(forTraj(2,1)-forTraj(2,0))/dt;
    forAngularVelocity(1)=(forTraj(5,1)-forTraj(5,0))/dt;

    forVelocity.block<3,1>(0,0)=AngleAxisd(walkingDirection, Vector3d::UnitZ())*
    							AngleAxisd(-forTraj(2), Vector3d::UnitZ())*forVelocity.block<3,1>(0,0);
	forVelocity.block<3,1>(0,1)=AngleAxisd(walkingDirection, Vector3d::UnitZ())*
    							AngleAxisd(-forTraj(5), Vector3d::UnitZ())*forVelocity.block<3,1>(0,1);

	// get front girdle velocities in the body (girdle) frame
    girdleVelocity(0,0)=(girdleTraj(0,1)-girdleTraj(0,0))/dt;
    girdleVelocity(1,0)=(girdleTraj(1,1)-girdleTraj(1,0))/dt;
    girdleVelocity(2,0)=0;
    girdleVelocity.block<3,1>(0,0)=AngleAxisd(walkingDirection, Vector3d::UnitZ())*AngleAxisd(-girdleTraj(2,1), Vector3d::UnitZ())*girdleVelocity.block<3,1>(0,0);
    girdleAngularVelocity(0)=(girdleTraj(2,1)-girdleTraj(2,0))/dt;

    // move front girdle
	moveGirdle(girdleVelocity.block<3,1>(0,0), girdleAngularVelocity(0), 0);

	// get hind girdle velocities in the body (girdle) frame
	girdleVelocity(0,1)=(girdleTraj(3,1)-girdleTraj(3,0))/dt;
    girdleVelocity(1,1)=(girdleTraj(4,1)-girdleTraj(4,0))/dt;
    girdleVelocity(2,1)=0;
    girdleVelocity.block<3,1>(0,1)=AngleAxisd(walkingDirection, Vector3d::UnitZ())*AngleAxisd(-girdleTraj(5,1), Vector3d::UnitZ())*girdleVelocity.block<3,1>(0,1);
    girdleAngularVelocity(1)=(girdleTraj(5,1)-girdleTraj(5,0))/dt;

    // move hind girdle
	moveGirdle(girdleVelocity.block<3,1>(0,1),girdleAngularVelocity(1), 1);

	girdleVelocity_filtered=pt1_vec(girdleVelocity, girdleVelocity_filtered, Tf1, dt);
	girdleAngularVelocity_filtered=pt1_vec(girdleAngularVelocity, girdleAngularVelocity_filtered, Tf1, dt);
	forVelocity_filtered=pt1_vec(forVelocity, forVelocity_filtered, Tf1, dt);
	forAngularVelocity_filtered=pt1_vec(forAngularVelocity, forAngularVelocity_filtered, Tf1, dt);


}


/* high level walking controller for locomotion through pipes */
void 
Controller :: crawlingStateMachine()
{
	

	static bool runLegDynamics=true;
	static int leg=0;
	static vector<int> legs;
	static int cnt=100;
	static Vector4d loadingLegShouldStop;
	static int crawlingSequenceSelector=-1;


	static MatrixXd maxF=MatrixXd::Zero(4,1);
	static Vector4d t_loading = Vector4d::Zero(), swing_phase = Vector4d::Zero();
	static Vector4d legsDone = Vector4d::Zero();
	Vector3d fg, dw, ddx, fg_effective;
	fg << 0,0,-9.81;
	dw << 0,0,0; 
	ddx << 0,0,0;

	fg_effective(0)=fg(2)*cos(forRPY(1))*mu_ground + fg(2)*sin(-forRPY(1));
	fg_effective(1)=0;
	fg_effective(2)=0;
	if(!crawling_continue){
		swing_phase << 0.5, 0.5, 0.5, 0.5;
	}
	//===================== initial state =============================
	if(crstate==initial){
		loadingLegShouldStop << 0,0,0,0;
		swing_phase(0)+=freq_walk*dt;
		for(int i=0; i<4; i++){
			if(swing_phase(0)<=1){
				moveSwingLegManualPhase(i, swing_phase(0));
			}
			if(!crawling_continue){
				continue;
			}
			else if(!legs_contact(i)){
				if(abs(feetReference(1,i))>abs(WP.midStance(1,i))){
					GP.nSurf.block<3,1>(0,i)=WP.nSurf.block<3,1>(0,i);
				}
				moveAlongVector(i, -GP.nSurf.block<3,1>(0,i), leg_approach_speed/4);
				
			}
			else{
				legs_contact_control(i)=1;
				maxF(i) = maxF(i) + dt*leg_loading_filter;
				//maxF(i)=pt1(maxContactForce, maxF(i), leg_loading_filter, dt);
			}	
		}

		if(swing_phase(0)>1 && legs_contact_control.sum()==4){
			t_loading(0) +=dt;
		}

		 if(maxF(0) >= maxContactForce){
		 	maxF(0) = maxContactForce;
			t_loading(0)=0;
			crstate=moveBody;
			swing_phase(0)=0;
			legs_stance << 1,1,1,1;
		}
	}	

	

	//===================== moving sequence selector =============================
	if(crstate==selectLegs){
		crawlingSequenceSelector++;

		if(crawlingSequenceSelector>=crawlingSequence.rows()){
			crstate=moveBody;
			crawlingSequenceSelector = -1;

			// remember feet positions
			feetTrajHist.push_back(feetLocationsOdo);
			if(feetTrajHist.size()>5){
				feetTrajHist.erase(feetTrajHist.begin(),feetTrajHist.begin()+1);
			}
			for(unsigned int i=0; i<feetTrajHist.size(); i++){
				cout << feetTrajHist[i] << endl << endl;
			}
			cout << endl << endl;
		}
		else{
			crstate=unloadingLeg;
		}
	}

	// get currently active legs
	legs.clear();
	for(int i=0; i<4; i++){
		if(crawlingSequenceSelector>=0 && crawlingSequence(crawlingSequenceSelector, i)){
			legs.push_back(i);
		}
	}

	if(crstate==unloadingLeg){
		for(int i=0; i<legs.size(); i++){
			maxF(legs[i]) = feet_force_reference.block<3,1>(0,legs[i]).cwiseAbs().maxCoeff();
		}
	}
	
	
	

	// ============================== COMMON STATES FOR ALL THE LEGS =======================================
	for(unsigned int k=0; k<legs.size(); k++){
		leg = legs[k];

		//===================== unloading leg =============================
		if(crstate==unloadingLeg){
			//cout << "unloadingLeg " << leg<< "\t\t tloading " << t_loading.transpose() << "\t\t legs done " << legsDone.transpose() << endl;
			legs_contact_control(leg)=1;


			/*maxF(leg)=pt1(0, maxF(leg), leg_loading_filter, dt);
			t_loading(leg) +=dt;
			// check if the leg is done with the current state
			if(t_loading(leg) > 3*leg_loading_filter){
				legsDone(leg)=1;
				GP.nSurf.block<3,1>(0,leg)=CP.nSurf.block<3,1>(0,leg);
			}*/

			// decrease max force 
			maxF(leg) = maxF(leg) - dt * leg_loading_filter;
			//check if the leg is done with the current state
			if(maxF(leg)<=0){
				maxF(leg)=0;
				legsDone(leg)=1;
				GP.nSurf.block<3,1>(0,leg)=CP.nSurf.block<3,1>(0,leg);
			}

			// switch the state if all the legs finished with the current state
			if(legsDone.sum() == legs.size()){
				t_loading.setZero();
				swing_phase.setZero();
				legsDone.setZero();
				GP.nSurf.block<3,1>(0,leg)=CP.nSurf.block<3,1>(0,leg);
				crstate=swingLeg;
			}
		}

		//===================== run phase =============================
		if(crstate==swingLeg){
			//cout << "swingLeg "<< leg<< endl;
			legs_contact_control(leg)=0;
			maxF(leg)=0;
			//===================== move leg =============================
			if(!legsDone(leg)){
				legsDone(leg)=moveSwingLegManualPhase(leg, swing_phase(leg));
				swing_phase(leg)+=freq_walk*dt;
			}
			//===================== no contact =============================
			if(!legs_contact(leg)){
				// if a leg cannot reach a wall, change objective to search for a ground plane
				if(abs(feetReference(1,leg))>abs(WP.midStance(1,leg))+0.05){
					GP.nSurf.block<3,1>(0,leg)=WP.nSurf.block<3,1>(0,leg);
				}
				// reach out
				moveAlongVector(leg, -GP.nSurf.block<3,1>(0,leg), leg_approach_speed);
			}

			// if all the finished legs are in contact, switch the state
			if(legs_contact.cwiseProduct(legsDone).sum() == legs.size()){
				crstate=loadingLeg;
				legsDone.setZero();
			}

		}	

		//===================== loading leg =============================
		if(crstate==loadingLeg){
			//cout << "loadingLeg " << leg<< endl;
			legs_contact_control(leg)=1;

			/*maxF(leg)=pt1(maxContactForce, maxF(leg), leg_loading_filter, dt);

			t_loading(leg) +=dt;

			// check if the leg is done with the current state
			if(t_loading(leg) > 3*leg_loading_filter){
				legsDone(leg)=1;
			}*/


			maxF(leg) = maxF(leg) + dt * leg_loading_filter;
			if(maxF(leg) >= maxContactForce || loadingLegShouldStop(leg)){
				maxF(leg) = maxContactForce;
				legsDone(leg)=1;
			}

			// switch the state if all the legs finished with the current state
			if(legsDone.sum() == legs.size())
			{
				t_loading.setZero();
				legsDone.setZero();
				crstate=selectLegs;
			}

		}
	}

	//===================== move body =============================
	if(crstate==moveBody){
		//cout << "moveBody"<< endl;
		static double t_moveBody=0;
		legs_stance << 1,1,1,1;
		// set forces to max
		for(int i=0; i<4; i++){
			maxF(i)=maxContactForce;
		}
		
		
		walking_angular_velocity=0;
		/*if(feetTrajHist.size()<=1){
			feetTrajHist.push_back(feetLocationsOdo);
		}*/
		if(feetTrajHist.size()>1){
			pipePathPlanner();
		}
		//pipePathPlanner();
		
		

		// get girdle trajectories  
		// double v, double w, Vector3d V, double W, double fgirdArc, Vector2d spineCPGscaling
		girdleTrajectories(walking_forward_velocity, walking_angular_velocity, MatrixXd::Zero(3,1), 0, 0, MatrixXd::Zero(2,1));
		//moveBothGirdles(MatrixXd::Zero(3,1), -walking_angular_velocity);

		// udpate girdle velocities
		girdleVelocities();

		// check if need to swing a leg
		t_moveBody+=dt;
		if(t_moveBody > (1/freq_walk*GP.Duty(0))){
			t_moveBody=0;
			crstate=selectLegs;
		}


		// control ROLL
		for(int i=0; i<4; i++){
			feetReference.block<3,1>(0,i) = AngleAxisd(0.1 * forRPY(0), Vector3d::UnitX()) * feetReference.block<3,1>(0,i);
		}



	}
	else{
		// double v, double w, Vector3d V, double W, double fgirdArc, Vector2d spineCPGscaling
		girdleTrajectories(0, 0, MatrixXd::Zero(3,1), 0, 0, MatrixXd::Zero(2,1));
		//moveBothGirdles(MatrixXd::Zero(3,1), 0);

		// udpate girdle velocities
		girdleVelocities();
	}

	// ============================== GET FORCE DISTRIBUTION ===============================
	

	dw << 0,0,0;

	feet_force_reference=calcForceDistribution(GP.nSurf, maxF, ddx, dw, fg_effective);
	for(int i=0; i<4; i++){
		if(abs(feet_force_reference(1,i))>(maxF(i) + 5)){
			loadingLegShouldStop(i)=1;
		}
		else{
			loadingLegShouldStop(i)=0;
		}
	}
	

	// ============================== CONTACT FORCE CONTROL ===============================

	feetReference_forceCorrections=contactNormalForceController(-feet_force_reference, -feet_force, GP.nSurf);


	//cout << "state: " << crstate << "\t" << "\t contact: " << legs_contact.transpose()<<"\t legsDone: " << legsDone.transpose() << endl;
	//cout << feetReference_forceCorrections << endl;


	



}