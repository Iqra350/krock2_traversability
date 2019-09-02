#include "controller.hpp"

using namespace std;
using namespace Eigen;



extern int LOG_DATA;


template <class MatT>
Eigen::Matrix<typename MatT::Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime>
pseudoInverse(const MatT &mat, typename MatT::Scalar tolerance = typename MatT::Scalar{1e-4}) // choose appropriately
{
    typedef typename MatT::Scalar Scalar;
    auto svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    const auto &singularValues = svd.singularValues();
    Eigen::Matrix<Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime> singularValuesInv(mat.cols(), mat.rows());
    singularValuesInv.setZero();
    for (unsigned int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > tolerance)
        {
            singularValuesInv(i, i) = Scalar{1} / singularValues(i);
        }
        else
        {
            singularValuesInv(i, i) = Scalar{0};
        }
    }
    return svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
}





// ---------------------------------------------- SOLVE FDO THREAD ---------------------------------------------------------
void
FDOThread(FDO *fdo, MatrixXd  *sol)
{
    double t1, t2, t3, dt;
    while(true){
        t1=get_timestamp();

        // update arguments with mutex protection
        fdo_in_mtx.lock();
            fdo->updateArgs();
        fdo_in_mtx.unlock();

        fdo->solveNLP();

        
        // update results with mutex protection
        fdo_out_mtx.lock();
            *sol = fdo->getSolution();
            //cout << *sol << endl;
        fdo_out_mtx.unlock();
        
        t2=get_timestamp();

        if(fdo->dt > (t2-t1)){
            usleep((fdo->dt - (t2-t1))*1000000);
        }
        //cout << "FDO time: " << (t2-t1)*1000 <<" ms"<<endl;
    }   
}

/* Get force sensors */
void 
Controller :: getForce(double *force)
{
	static MatrixXd feet_force_old=MatrixXd::Zero(3, 4);
    MatrixXd ftmp(12,1);
    for(int i=0; i<12; i++){
        ftmp(i)=force[i];
    }

    feet_force.block<3,1>(0,0)=HJfl_g[4].block<3,3>(0,0)*ftmp.block<3,1>(0,0);
    feet_force.block<3,1>(0,1)=HJfr_g[4].block<3,3>(0,0)*ftmp.block<3,1>(3,0);
    feet_force.block<3,1>(0,2)=HJhl_g[4].block<3,3>(0,0)*ftmp.block<3,1>(6,0);
    feet_force.block<3,1>(0,3)=HJhr_g[4].block<3,3>(0,0)*ftmp.block<3,1>(9,0);

    // filter
    feet_force_filt=pt1_vec(feet_force, feet_force_old, rawForcesFilter, dt);
	feet_force_old=feet_force_filt;


	// detect contact
	for(int i=0; i<4; i++){
		if(feet_force.block<3,1>(0,i).norm()>forceTreshold){
			legs_contact(i)=1;
		}
		else{
			legs_contact(i)=0;
		}
	}

}


/* Control the contact forces by adjusting the feet reference position */
Matrix<double, 3, 4> 
Controller :: contactNormalForceController(Matrix<double, 3, 4> feetForceRef, Matrix<double, 3, 4> feetForceEst, Matrix<double, 3, 4> nVec)
{
	MatrixXd feetRefCorrections(3,4);
	static MatrixXd feetRefCorrections_old=MatrixXd::Zero(3,4);
	feetRefCorrections.setZero();

	static MatrixXd e(4,1), tmp(1,1); 
	static MatrixXd eint = MatrixXd::Zero(4,1), u_old = MatrixXd::Zero(4,1), u = MatrixXd::Zero(4,1);
	

	for(int i=0; i<4; i++){

		tmp = nVec.block<3,1>(0,i).transpose() * (feetForceRef.block<3,1>(0,i) - feetForceEst.block<3,1>(0,i));
		e(i) = tmp(0);

		eint(i) = eint(i) + e(i)*dt;

		u(i) = 1/(1+contactF_pid_Kd/dt) * ( contactF_pid_Kp * e(i) + contactF_pid_Ki * eint(i) + contactF_pid_Kd/dt * u_old(i)  );

		u_old(i) = u(i);

	}
	//=================================== run controllers ==================================================
	

	for(int i=0; i<4; i++){
		if(legs_contact_control(i)){
			feetRefCorrections.block<3,1>(0,i)=nVec.block<3,1>(0,i)*u(i);
			feetRefCorrections_old.block<3,1>(0,i)=feetRefCorrections.block<3,1>(0,i);
		}
		else{
			// prevent jumps in feetRefCorrections once force control switches off
			feetRefCorrections.block<3,1>(0,i) = pt1_vec(feetRefCorrections.block<3,1>(0,i), feetRefCorrections_old.block<3,1>(0,i), 0.5, dt);
			feetRefCorrections_old.block<3,1>(0,i)=feetRefCorrections.block<3,1>(0,i);
			eint(i)=0;
			u_old(i)=0;
		}
	}

	return feetRefCorrections;
	
}


MatrixXd
Controller::calcForceDistribution(MatrixXd NS, MatrixXd maxForce, MatrixXd ddx, MatrixXd dw, MatrixXd fg)
{
	MatrixXd XC(3,4), forceDistribReference(3,4);
    

    //================= FEET LOCATIONS IN RESPECT TO THE COM =====================

    // com
    Vector3d com=getCoM();
    for(int i=0; i<4; i++){
    	XC.block<3,1>(0,i)=feetLocations.block<3,1>(0,i) - com;
    }


    static MatrixXd fdo_sol(3,4);
    // update parameters
/*    cout << "XC\t" << XC << endl << endl;
    cout << "NS\t" << NS << endl << endl;
    cout << "maxForce\t" << maxForce << endl << endl;
    cout << "fg\t" << fg << endl << endl;
    cout << "ddx\t" << ddx << endl << endl;
    cout << "dw\t" << dw << endl << endl;*/
    fdo_in_mtx.lock();
        fdo.updateParameters(XC, NS, maxForce, fg, ddx, dw);
    fdo_in_mtx.unlock();
    
    // start thread
    static bool fdo_thread_started=false;
    if(!fdo_thread_started){
        cout << "STARTING FDO THREAD" << endl;
        fdo_thread=std::thread(FDOThread, &fdo, &fdo_sol);
        fdo_thread.detach();

       fdo_thread_started=true;
    }

    // get solution
    fdo_out_mtx.lock();
        forceDistribReference=fdo_sol;
    fdo_out_mtx.unlock();

    return forceDistribReference;

}

/* Get force from servo current */
void
Controller :: getForceFromCurrent()
{	
	for(int leg=0; leg<4; leg++){
		MatrixXd J=Jacob(fbck_angles.block<4,1>(0+4*leg,0), leg);
		feet_force_from_current.block<3,1>(0,leg)= pseudoInverse(J.transpose())*fbck_current.block<4,1>(0+4*leg,0);
	}
}



/* Use neural network to estimate contact force */
void
Controller :: mlpContactForceEstimationDummy()
{
	nnEstForce=feet_force;
}


