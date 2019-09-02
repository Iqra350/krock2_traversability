#ifndef OPTIMIZATION_TOOLS_HPP
#define OPTIMIZATION_TOOLS_HPP

#include <casadi/casadi.hpp>
#include <fstream>
#include <pthread.h>
#include <vector>
#include <string.h>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <dlib-19.0/dlib/optimization.h>
#include <dlib-19.0/dlib/optimization/find_optimal_parameters.h>
#include "qpOASES-3.2.0/include/qpOASES.hpp"
#include <mutex> 


using namespace Eigen;



using namespace std;


class MPC{

	public:
    //================== public variables ===============================================
	int n, m, ne;
	bool enabled;
	int N;
	double dt;
	int n_state_constraints, n_input_constraints, n_poly, n_soft;
	int Esoft_pattern[50];
	// system matrices
	Matrix<double, Dynamic, Dynamic> A, B;	
	//cost matrices
	Matrix<double, Dynamic, Dynamic> Qx, Ru, Qf, Psoft, Qx0, Ru0, Qf0;
	//constraint matrices
	Matrix<double, Dynamic, Dynamic> Fx, fx, Mu, mu, Esoft;
	Matrix<double, Dynamic, 2> u_min_max, x_min_max, u_min_max0, x_min_max0;

	// variables
	Matrix<double, Dynamic, 1> xref;

	// solutions
	Matrix<double, Dynamic,1> sol, xsol, usol, esol;
	qpOASES::SQProblem oasesSolver;

    //================== public functions ===============================================
	MPC();
	void constructQP();
	void setConstraintMatrices();
	void readParameters(const char *filename);
	void printMatrices(const char *filename);
	void initQPSolver();
	void updateConstraints_SupportPolygons(std::vector<Matrix<double,3,4>> supportPolys);
	void updateConstraints_Box(MatrixXd x_min_max_new, MatrixXd u_min_max_new);
	void updateReference(MatrixXd xref_in);
	void updateWeights(MatrixXd Qx_new, MatrixXd Qf_new, MatrixXd Ru_new);
	void updateState(MatrixXd x0);
	void solveQP();
	void updateQpArg();

	


    private:
    //================== private variables ===============================================
    bool USE_CASADI, use_sparse;
	Matrix<double, Dynamic, Dynamic> H_qp, h_qp, G_qp, g_qp, T_qp, t_qp;	   
	Matrix<double, Dynamic, Dynamic> A_qp, lba, uba;	 
	casadi::Function qp_solver;  	
	casadi::DM cas_H_qp, cas_h_qp, cas_G_qp, cas_g_qp, cas_T_qp, cas_t_qp;	
	casadi::DM cas_A_qp, cas_lba, cas_uba;	
	vector<double> vec_H_qp, vec_h_qp, vec_G_qp, vec_g_qp, vec_T_qp, vec_t_qp, vec_A_qp, vec_lba, vec_uba;
	casadi::DMDict qp_arg;
	double Pmargin;

	qpOASES::real_t *H_qpoases, *h_qpoases, *A_qpoases, *lba_qpoases, *uba_qpoases;	
	qpOASES::SymSparseMat *Hs_qpoases;	
	qpOASES::SparseMatrix *As_qpoases;	


	//================== private functions ===============================================
	



};




void solveSomeQP();


void forceDistributionIPOPT(MatrixXd XC, MatrixXd NS, MatrixXd contact, Vector3d fg, Vector3d ddx, Vector3d dw, double m, double mu);
//void setupMPC(MPC *mpc);


class FDO{

	public:
		FDO(){};
		void initSolver(double m, double mu, double dt_in);
		void updateParameters(MatrixXd XC, MatrixXd NS, MatrixXd maxForce, Vector3d fg, Vector3d ddx, Vector3d dw);
		void updateArgs();
		void solveNLP();
		Matrix<double, 3, 4>  getSolution();

		double dt;

	private:
		casadi::Function solver;
		casadi::DMDict arg, arg_tmp;
		casadi::DMDict res;


};



#endif
