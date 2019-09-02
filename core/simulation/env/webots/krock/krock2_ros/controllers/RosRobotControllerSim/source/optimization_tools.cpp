#include "optimization_tools.hpp" 
#include "utils.hpp"

std::mutex mpc_in_mtx, mpc_out_mtx, mpc_sync_mtx;
std::mutex fdo_in_mtx, fdo_out_mtx, fdo_sync_mtx;

/*  Constructor */
MPC :: MPC()
{
	n=0;
	m=0;
	n_soft=0;
	N=1;
}

/* Read config file with MPC parameters (weights, constraints) */
void 
MPC :: readParameters(const char *filename)
{
	// prep file reading
	stringstream stringstream_file;
	ifstream file_mpc;
	file_mpc.open(filename);
	readFileWithLineSkipping(file_mpc, stringstream_file);

	//---------------------------------------- SYSTEM INFO ---------------------------------------------------
	stringstream_file >> enabled;
	stringstream_file >> USE_CASADI;
	stringstream_file >> use_sparse;
	stringstream_file >> n;
	stringstream_file >> m;
	stringstream_file >> n_soft;

	// A & B matrices
	A.resize(n,n);
	B.resize(n,m);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			stringstream_file >> A(i,j);
		}
	}
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			stringstream_file >> B(i,j);
		}
	}

	//---------------------------------------- HORIZON --------------------------------------------------- 
	stringstream_file >> dt;

	stringstream_file >> N;

	//---------------------------------------- COST MATRICES --------------------------------------------------- 
	Qx=MatrixXd::Zero(n,n);
	Qf=MatrixXd::Zero(n,n);
	Ru=MatrixXd::Zero(m,m);
	Qx0=MatrixXd::Zero(n,n);
	Qf0=MatrixXd::Zero(n,n);
	Ru0=MatrixXd::Zero(m,m);
	for(int i=0; i<n; i++){
		stringstream_file >> Qx(i,i);
		Qx0(i,i)=Qx(i,i);
	}
	for(int i=0; i<n; i++){
		stringstream_file >> Qf(i,i);
		Qf0(i,i)=Qf(i,i);
	}
	for(int i=0; i<m; i++){
		stringstream_file >> Ru(i,i);
		Ru0(i,i)=Ru(i,i);
	}


	//---------------------------------------- BOX CONSTRAINTS ---------------------------------------------------
	x_min_max.resize(n,2); x_min_max0.resize(n,2);
	for(int i=0; i<n; i++){
		stringstream_file >> x_min_max(i,0);
		stringstream_file >> x_min_max(i,1);
	}
	x_min_max0=x_min_max;

	u_min_max.resize(n,2); u_min_max0.resize(n,2);
	for(int i=0; i<m; i++){
		stringstream_file >> u_min_max(i,0);
		stringstream_file >> u_min_max(i,1);
	}
	u_min_max0=u_min_max;
	//------------------------------------- ADITIONAL and SOFT CONSTRAINTS -------------------------------------------
	stringstream_file >> n_poly;
	stringstream_file >> n_soft;

	for(int i=0; i<(2*n+n_poly); i++){
		stringstream_file >> Esoft_pattern[i];
	}

	Psoft.resize(n_soft, n_soft);
	Psoft=MatrixXd::Zero(n_soft, n_soft);
	for(int i=0; i<n_soft; i++){
		stringstream_file >> Psoft(i,i);
	}
	stringstream_file >> Pmargin;
}

/*  Set constraint matrices */
void
MPC :: setConstraintMatrices()
{
	//Fx
	Fx.resize(2*n+n_poly, n);
	Fx=MatrixXd::Zero(2*n+n_poly, n);
	for(int i=0; i<n; i++){
		Fx(2*i,i)=-1;
		Fx(2*i+1,i)=1;
	}

	//fx
	fx.resize(2*n+n_poly, 1);
	fx=MatrixXd::Zero(2*n+n_poly, 1);
	for(int i=0; i<n; i++){
		fx(2*i,0)=-x_min_max(i,0);
		fx(2*i+1,0)=x_min_max(i,1);
	}

	//Mu
	Mu.resize(2*m,m);
	Mu=MatrixXd::Zero(2*m,m);
	for(int i=0; i<m; i++){
		Mu(2*i,i)=-1;
		Mu(2*i+1,i)=1;
	}

	//mu
	mu.resize(2*m,1);
	mu=MatrixXd::Zero(2*m,1);
	for(int i=0; i<m; i++){
		mu(2*i,0)=-u_min_max(i,0);
		mu(2*i+1,0)=u_min_max(i,1);
	}

	// embed Esoft
	Esoft.resize(2*n+n_poly, n_soft);
	Esoft=MatrixXd::Zero(2*n+n_poly, n_soft);
	for(int j=0, i=0; i<2*n+n_poly; i++){
		if(Esoft_pattern[i]==1){
			Esoft(i, j)=-1;
			j++;
		}
	}
}

/*  Construct QP matrices from the system and constraints */
void
MPC :: constructQP()
{


	//------------------------ INEQUALITY MATRIX-------------------------------
	const int Frows=Fx.rows();
	const int Fcols=Fx.cols();
	const int Mrows=Mu.rows();
	const int Mcols=Mu.cols();
	const int Erows=Esoft.rows();
	const int Ecols=Esoft.cols();
	MatrixXd F_poly=MatrixXd::Zero(n_poly, n);
	F_poly.block(0,0,n_poly,1)=MatrixXd::Ones(n_poly,1);
	F_poly.block(0,3,n_poly,1)=MatrixXd::Ones(n_poly,1);
	// put hard constraints into the G_qp matrix (diagonal elements)
	// G_qp
	G_qp.resize((Frows+Mrows+n_soft)*N,(Fcols+Mcols+n_soft)*N);
	G_qp=MatrixXd::Zero((Frows+Mrows+n_soft)*N,(Fcols+Mcols+n_soft)*N);
	for(int i=0; i<N; i++){
		G_qp.block(i*Frows, i*Fcols, Frows, Fcols)=Fx;
		G_qp.block(N*Frows+i*Mrows, N*Fcols+i*Mcols, Mrows, Mcols)=Mu;

		if(n_soft>0){
			G_qp.block(N*(Frows+Mrows)+i*n_soft, N*(Fcols+Mcols)+i*n_soft, n_soft, n_soft)=-MatrixXd::Identity(n_soft,n_soft);
			G_qp.block(i*Frows, N*(Fcols+Mcols)+i*n_soft, Erows, Ecols)=Esoft;
			G_qp.block(i*Frows+2*n, i*Fcols, n_poly, n)=F_poly;
		}
	}

	// g_qp
	g_qp.resize(G_qp.rows(),1);
	g_qp=MatrixXd::Zero(G_qp.rows(),1);
	for(int i=0; i<N; i++){
		g_qp.block(i*Frows, 0, Frows, 1)=fx;
		g_qp.block(N*Frows+i*Mrows, 0, Mrows, 1)=mu;
	}

	




	//------------------------ EQUALITY MATRIX-------------------------------
	// T_qp
	T_qp.resize(n*N,(n+m)*N);
	T_qp=MatrixXd::Zero(n*N,(n+m+n_soft)*N);
	for(int i=0; i<N; i++){
		T_qp.block(i*n, i*n, n, n)=MatrixXd::Identity(n,n);
		T_qp.block(i*n, N*n+i*m, n, m)=-B;
	}
	for(int i=0; i<(N-1); i++){
		T_qp.block((i+1)*n, i*n, n, n)=-A;
	}

	//t_qp
	t_qp.resize(n*N,n);
	t_qp=MatrixXd::Zero(n*N,n);
	t_qp.block(0,0, n, n)=A;



	//------------------------ COST MATRIX-------------------------------
	//H_qp
	H_qp.resize((n+m+n_soft)*N,(n+m+n_soft)*N);
	H_qp=MatrixXd::Zero((n+m+n_soft)*N,(n+m+n_soft)*N);
	for(int i=0; i<N; i++){
		H_qp.block(i*n, i*n, n, n)=Qx;
		H_qp.block(N*n+i*m, N*n+i*m, m, m)=Ru;
		if(n_soft>0){
			H_qp.block(N*(n+m)+i*n_soft,N*(n+m)+i*n_soft, n_soft, n_soft)=Psoft;//*pow(10,i);
		}
	}
	H_qp.block((N-1)*n, (N-1)*n, n, n)=Qf;

	//h_qp & xref
	h_qp.resize((n+m+n_soft)*N,1);
	h_qp=MatrixXd::Zero((n+m+n_soft)*N,1);
	xref.resize(n*N,1);
	//xref=MatrixXd::Zero(n*N,1);
	xref=MatrixXd::Ones(n*N,1);
	for(int i=0; i<N; i++){
		h_qp.block(i*n, 0, n, 1)=-H_qp.block(i*n, i*n, n, n)*xref.block(i*n,0,n,1);
	}


	//------------------------ VERTCAT MATRICES-------------------------------
	MatrixXd CoM(6,1);
	CoM << 1, 2, 3, 4, 5, 6;
	A_qp.resize(G_qp.rows()+T_qp.rows(), G_qp.cols());
	A_qp << G_qp,
			T_qp;

	lba.resize(A_qp.rows(), 1);
	lba << -1.0e20*MatrixXd::Ones(G_qp.rows(),1),
			t_qp*CoM;

	uba.resize(A_qp.rows(), 1);
	uba << g_qp,
			t_qp*CoM;


	//----------------------- INIT QP SOLVER ------------------------------
	//initQPSolver();
}

/*  Initializes QP solver */
void
MPC :: initQPSolver()
{

	// EIGEN MATRICES
	MatrixXd H2(2,2);
	H2=2*MatrixXd::Identity(2,2);            
	
	//============================================== CASADI INITIALIZATION ===============================================
	if(USE_CASADI)
	{
		vec_H_qp.resize(H_qp.size());
		vec_h_qp.resize(h_qp.size());
		vec_G_qp.resize(G_qp.size());
		vec_g_qp.resize(g_qp.size());
		vec_T_qp.resize(T_qp.size());
		vec_t_qp.resize(t_qp.size());

		vec_A_qp.resize(A_qp.size());
		vec_lba.resize(lba.size());
		vec_uba.resize(uba.size());



		// map to std::vector
		Map<MatrixXd>( &vec_H_qp[0], H_qp.rows(), H_qp.cols() ) = H_qp;
		Map<MatrixXd>( &vec_h_qp[0], h_qp.rows(), h_qp.cols() ) = h_qp;
		Map<MatrixXd>( &vec_G_qp[0], G_qp.rows(), G_qp.cols() ) = G_qp;
		Map<MatrixXd>( &vec_g_qp[0], g_qp.rows(), g_qp.cols() ) = g_qp;
		Map<MatrixXd>( &vec_T_qp[0], T_qp.rows(), T_qp.cols() ) = T_qp;
		Map<MatrixXd>( &vec_t_qp[0], t_qp.rows(), t_qp.cols() ) = t_qp;

		Map<MatrixXd>( &vec_A_qp[0], A_qp.rows(), A_qp.cols() ) = A_qp;
		Map<MatrixXd>( &vec_lba[0], lba.rows(), lba.cols() ) = lba;
		Map<MatrixXd>( &vec_uba[0], uba.rows(), uba.cols() ) = uba;

		// construct casadi matrices
		cas_H_qp=casadi::Matrix<double>(vec_H_qp);
		cas_h_qp=casadi::Matrix<double>(vec_h_qp);
		cas_G_qp=casadi::Matrix<double>(vec_G_qp);
		cas_g_qp=casadi::Matrix<double>(vec_g_qp);
		cas_T_qp=casadi::Matrix<double>(vec_T_qp);
		cas_t_qp=casadi::Matrix<double>(vec_t_qp);

		cas_A_qp=casadi::Matrix<double>(vec_A_qp);
		cas_lba=casadi::Matrix<double>(vec_lba);
		cas_uba=casadi::Matrix<double>(vec_uba);

		// reshape casadi matrices to the correct size
		cas_H_qp=cas_H_qp.reshape(cas_H_qp, H_qp.rows(), H_qp.cols());
		cas_h_qp=cas_h_qp.reshape(cas_h_qp, h_qp.rows(), h_qp.cols());
		cas_G_qp=cas_G_qp.reshape(cas_G_qp, G_qp.rows(), G_qp.cols());
		cas_g_qp=cas_g_qp.reshape(cas_g_qp, g_qp.rows(), g_qp.cols());
		cas_T_qp=cas_T_qp.reshape(cas_T_qp, T_qp.rows(), T_qp.cols());
		cas_t_qp=cas_t_qp.reshape(cas_t_qp, t_qp.rows(), t_qp.cols());

		cas_A_qp=cas_A_qp.reshape(cas_A_qp, A_qp.rows(), A_qp.cols());
		cas_lba=cas_lba.reshape(cas_lba, lba.rows(), lba.cols());
		cas_uba=cas_uba.reshape(cas_uba, uba.rows(), uba.cols());


		// create QP solver
		casadi::SpDict qp={};
		qp["h"]=cas_H_qp.sparsity();
		qp["a"]=cas_A_qp.sparsity();
		casadi::Dict opts;
		opts["printLevel"] = "none";

		// my MPC config
		opts["terminationTolerance"] = 1.0e-6;
		opts["numRegularisationSteps"] = 1;
		opts["numRefinementSteps"] = 0;
		opts["enableNZCTests"] = false;
		opts["enableFlippingBounds"] = false;
		opts["enableRegularisation"] = true;
		opts["enableRamping"] = false;
		opts["enableEqualities"] = true;
		//opts["sparse"] = true;

		//qp_solver = casadi::qpsol("solver", "qpoases", qp, opts);
	}
	else
	{
	//============================================== QPOASES INITIALIZATION ===============================================
		H_qpoases=new qpOASES::real_t[H_qp.rows()*H_qp.cols()];
		h_qpoases=new qpOASES::real_t[h_qp.rows()*h_qp.cols()];
		A_qpoases=new qpOASES::real_t[A_qp.rows()*A_qp.cols()];
		lba_qpoases=new qpOASES::real_t[lba.rows()*lba.cols()];
		uba_qpoases=new qpOASES::real_t[uba.rows()*uba.cols()];

		Map<MatrixXd>( &H_qpoases[0], H_qp.cols(), H_qp.rows() ) = H_qp.transpose();
		Map<MatrixXd>( &h_qpoases[0], h_qp.rows(), h_qp.cols() ) = h_qp;
		Map<MatrixXd>( &A_qpoases[0], A_qp.cols(), A_qp.rows() ) = A_qp.transpose();		
		Map<MatrixXd>( &lba_qpoases[0], lba.rows(), lba.cols() ) = lba;
		Map<MatrixXd>( &uba_qpoases[0], uba.rows(), uba.cols() ) = uba;

		cout << "CCCCCCCCCCCCCC: 1" << endl;
		cout << A_qp << endl;
		Hs_qpoases = new qpOASES::SymSparseMat(H_qp.rows(), H_qp.cols(), H_qp.rows(), &H_qpoases[0]);
		Hs_qpoases->createDiagInfo();
		cout << "CCCCCCCCCCCCCC: 2" << endl;
		As_qpoases = new qpOASES::SparseMatrix(A_qp.rows(), A_qp.cols(), A_qp.cols(), &A_qpoases[0]);
		cout << "CCCCCCCCCCCCCC: 3" << endl;

		qpOASES::SQProblem tmpSolver(H_qp.rows(), A_qp.rows());
		oasesSolver=tmpSolver;


		qpOASES::Options myOptions;


		//myOptions.setToDefault();
		//myOptions.enableRamping                 =  qpOASES::BT_FALSE;
		//myOptions.enableFarBounds               =  qpOASES::BT_TRUE;
		//myOptions.enableFlippingBounds          =  qpOASES::BT_FALSE;
		//myOptions.enableRegularisation          =  qpOASES::BT_TRUE;
		//myOptions.enableNZCTests                =  qpOASES::BT_FALSE;
		//myOptions.enableDriftCorrection         =  0;
		//myOptions.enableEqualities              =  qpOASES::BT_TRUE;


		//myOptions.terminationTolerance=1.0e9 * 2.221e-16;





		//myOptions.setToReliable();
		myOptions.setToMPC();
		//myOptions.enableEqualities              =  qpOASES::BT_FALSE;
		//myOptions.terminationTolerance = 1.0e7 * 2.221e-16;
		//myOptions.boundTolerance       =  1.0e5 * 2.221e-16;
		//myOptions.initialStatusBounds       =  qpOASES::ST_LOWER;


		//myOptions.enableRamping         =  qpOASES::BT_TRUE;






		myOptions.printLevel = qpOASES::PL_NONE;
		oasesSolver.setOptions(myOptions);



		qpOASES::int_t nWSR = 300;
		if(!use_sparse)
			oasesSolver.init( H_qpoases,h_qpoases,A_qpoases,NULL,NULL,lba_qpoases,uba_qpoases, nWSR, 0 );
		else
			oasesSolver.init( Hs_qpoases,h_qpoases,As_qpoases,NULL,NULL,lba_qpoases,uba_qpoases, nWSR, 0 );



	}
  	//SET TO MPC option modifications
	/*
	setToDefault( );

	enableRamping                 =  BT_FALSE;
	enableFarBounds               =  BT_TRUE;
	enableFlippingBounds          =  BT_FALSE;
	enableRegularisation          =  BT_TRUE;
	enableNZCTests                =  BT_FALSE;
	enableDriftCorrection         =  0;
	enableEqualities              =  BT_TRUE;

	terminationTolerance          =  1.0e9 * EPS;
	
	initialStatusBounds           =  ST_INACTIVE;
	numRegularisationSteps        =  2;
	numRefinementSteps            =  0;


	#ifdef __USE_SINGLE_PRECISION__
	const real_t EPS = 1.193e-07;
	#else
	const real_t EPS = 2.221e-16;
	#endif 

	*/

}


/*  Prints QP matrices into a file */
void
MPC :: printMatrices(const char *filename)
{

	ofstream file(filename);
	file << "#-------------------------------SYSTEM-----------------------------------" << endl;
	file << "# A" << endl;
	file << A << endl << endl;
	file << "# B" << endl;
	file << B << endl << endl;

	file << "#-------------------------------COST-----------------------------------" << endl;
	file << "# Qx" << endl;
	file << Qx << endl << endl;
	file << "# H_qp" << endl;
	file << H_qp << endl << endl;
	file << "# h_qp" << endl;
	file << h_qp << endl << endl;

	file << "#-------------------------------CONSTRAINTS-----------------------------------" << endl;
	file << "# G_qp" << endl;
	file << G_qp << endl << endl;

	file << "# g_qp" << endl;
	file << g_qp << endl << endl;

	file << "# T_qp" << endl;
	file << T_qp << endl << endl;

	file << "# t_qp" << endl;
	file << t_qp << endl << endl;

	file << "#-------------------------------QPOASES-----------------------------------" << endl;
	file << "# A_qp" << endl;
	file << A_qp << endl << endl;

	file << "# lba" << endl;
	file << lba << endl << endl;

	file << "# uba" << endl;
	file << uba << endl << endl;
}

/* Inserts new values of support polygon (soft) constraints into inequality matrices */
void 
MPC :: updateConstraints_SupportPolygons(std::vector<Matrix<double,3,4>> supportPolys)
{
	//------------------------ UPDATE INEQUALITY MATRIX-------------------------------
	const int Frows=Fx.rows();
	const int Fcols=Fx.cols();

	MatrixXd F_poly=MatrixXd::Zero(n_poly, n), supportPolysContainer=MatrixXd::Zero(2, n_poly+1), f_poly=MatrixXd::Zero(n_poly,1);
	// put hard constraints into the G_qp matrix (diagonal elements)
	// G_qp
	int active_indx;
	Vector2d w_poly;
	for(int i=0; i<N; i++){

		// reorder support polygon points and remove swing points
		active_indx=0;
		if(supportPolys[i](2,0)>-5){
			supportPolysContainer.block(0,active_indx,2,1)=supportPolys[i].block(0,0,2,1);  	//fl
			active_indx++;
		}
		if(supportPolys[i](2,1)>-5){
			supportPolysContainer.block(0,active_indx,2,1)=supportPolys[i].block(0,1,2,1);  	//fr
			active_indx++;
		}
		if(supportPolys[i](2,3)>-5){
			supportPolysContainer.block(0,active_indx,2,1)=supportPolys[i].block(0,3,2,1);  	//hr
			active_indx++;
		}
		if(supportPolys[i](2,2)>-5){
			supportPolysContainer.block(0,active_indx,2,1)=supportPolys[i].block(0,2,2,1);  	//hl
			active_indx++;
		}
		// add first as last point - close support polygon
		supportPolysContainer.block(0,active_indx,2,1)=supportPolysContainer.block(0,0,2,1);


		// create polytopic constraints from the support polygon

		for(int k=0; k<active_indx; k++){
			// calc direction vector
			w_poly(0)=supportPolysContainer(1,k)-supportPolysContainer(1,k+1);
			w_poly(1)=supportPolysContainer(0,k+1)-supportPolysContainer(0,k);
			//F_poly(k,0)=supportPolysContainer(1,k)-supportPolysContainer(1,k+1);
			//F_poly(k,3)=supportPolysContainer(0,k+1)-supportPolysContainer(0,k);
			// normalize dir vector
			w_poly=w_poly/w_poly.norm();
			//F_poly(k,0)=F_poly(k,0)/sqrt(F_poly(k,0)*F_poly(k,0)+F_poly(k,3)*F_poly(k,3));
			//F_poly(k,3)=F_poly(k,3)/sqrt(F_poly(k,0)*F_poly(k,0)+F_poly(k,3)*F_poly(k,3));
			F_poly(k,0)=w_poly(0);
			F_poly(k,3)=w_poly(1);
			// calc offset
			f_poly(k,0)=F_poly(k,0)*supportPolysContainer(0,k) + F_poly(k,3)*supportPolysContainer(1,k) - Pmargin*i/(N-1);

		}


		// insert F_poly into G_qp matrix
		G_qp.block(i*Frows+2*n, i*Fcols, n_poly, n)=F_poly;

		// insert f_poly into g_qp matrix
		g_qp.block(i*Frows+2*n, 0, n_poly, 1)=f_poly;
	}
}


/* Inserts new values of support polygon (soft) constraints into inequality matrices */
void 
MPC :: updateConstraints_Box(MatrixXd x_min_max_new, MatrixXd u_min_max_new)
{
	x_min_max=x_min_max_new;
	u_min_max=u_min_max_new;

	//------------------------ UPDATE INEQUALITY MATRIX-------------------------------
	const int Frows=Fx.rows();
	const int Fcols=Fx.cols();
	const int Mrows=Mu.rows();
	const int Mcols=Mu.cols();

	//fx
	for(int i=0; i<n; i++){
		fx(2*i,0)=-x_min_max(i,0);
		fx(2*i+1,0)=x_min_max(i,1);
	}

	//mu
	for(int i=0; i<m; i++){
		mu(2*i,0)=-u_min_max(i,0);
		mu(2*i+1,0)=u_min_max(i,1);
	}

	// g_qp
	for(int i=0; i<N; i++){
		g_qp.block(i*Frows, 0, Frows, 1)=fx;
		g_qp.block(N*Frows+i*Mrows, 0, Mrows, 1)=mu;
	}

}

/* Inserts new values of reference over prediction horizon */
void 
MPC :: updateReference(MatrixXd xref_in)
{
	xref=xref_in;
	for(int i=0; i<N; i++){
		h_qp.block(i*n, 0, n, 1)=-H_qp.block(i*n, i*n, n, n)*xref.block(i*n,0,n,1);
	}



	
}

/* Updates current state */
void 
MPC :: updateState(MatrixXd x0)
{	

	lba << -casadi::inf*MatrixXd::Ones(G_qp.rows(),1),
			t_qp*x0;

	uba << g_qp,
			t_qp*x0;

		
}

/* Updates cost weight matrices */
void
MPC :: updateWeights(MatrixXd Qx_new, MatrixXd Qf_new, MatrixXd Ru_new)
{
	Qx=Qx_new;
	Qf=Qf_new;
	Ru=Ru_new;
	//------------------------ COST MATRIX-------------------------------
	//H_qp
	for(int i=0; i<N; i++){
		H_qp.block(i*n, i*n, n, n)=Qx;
		H_qp.block(N*n+i*m, N*n+i*m, m, m)=Ru;
	}
	H_qp.block((N-1)*n, (N-1)*n, n, n)=Qf;

	



}
/* Updates the arguments of the qp - not needed but allowes better mutex control */
void 
MPC :: updateQpArg()
{	

	// combine constraints into A_qp matrix
	A_qp << G_qp,
			T_qp;


	//------------------------ CONVERSION-------------------------------
	if(USE_CASADI)
	{		
		// map to std::vector
		Map<MatrixXd>( &vec_A_qp[0], A_qp.rows(), A_qp.cols() ) = A_qp;
		Map<MatrixXd>( &vec_h_qp[0], h_qp.rows(), h_qp.cols() ) = h_qp;
		Map<MatrixXd>( &vec_H_qp[0], H_qp.rows(), H_qp.cols() ) = H_qp;
		Map<MatrixXd>( &vec_lba[0], lba.rows(), lba.cols() ) = lba;
		Map<MatrixXd>( &vec_uba[0], uba.rows(), uba.cols() ) = uba;


		// construct casadi matrices
		cas_A_qp=casadi::Matrix<double>(vec_A_qp);
		cas_h_qp=casadi::Matrix<double>(vec_h_qp);
		cas_H_qp=casadi::Matrix<double>(vec_H_qp);
		cas_lba=casadi::Matrix<double>(vec_lba);
		cas_uba=casadi::Matrix<double>(vec_uba);

		// reshape casadi matrices to the correct size
		cas_A_qp=cas_A_qp.reshape(cas_A_qp, A_qp.rows(), A_qp.cols());
		cas_lba=cas_lba.reshape(cas_lba, lba.rows(), lba.cols());
		cas_uba=cas_uba.reshape(cas_uba, uba.rows(), uba.cols());	
		cas_h_qp=cas_h_qp.reshape(cas_h_qp, h_qp.rows(), h_qp.cols());
		cas_H_qp=cas_H_qp.reshape(cas_H_qp, H_qp.rows(), H_qp.cols());


		qp_arg={{"h",cas_H_qp}, {"g",cas_h_qp}, {"a", cas_A_qp}, {"lba", cas_lba}, {"uba", cas_uba}};
	}	
	else
	{
		Map<MatrixXd>( &H_qpoases[0], H_qp.cols(), H_qp.rows() ) = H_qp.transpose();
		Map<MatrixXd>( &h_qpoases[0], h_qp.rows(), h_qp.cols() ) = h_qp;
		Map<MatrixXd>( &A_qpoases[0], A_qp.cols(), A_qp.rows() ) = A_qp.transpose();		
		Map<MatrixXd>( &lba_qpoases[0], lba.rows(), lba.cols() ) = lba;
		Map<MatrixXd>( &uba_qpoases[0], uba.rows(), uba.cols() ) = uba;

		Hs_qpoases = new qpOASES::SymSparseMat(H_qp.rows(), H_qp.cols(), H_qp.rows(), &H_qpoases[0]);
		Hs_qpoases->createDiagInfo();
		As_qpoases = new qpOASES::SparseMatrix(A_qp.rows(), A_qp.cols(), A_qp.cols(), &A_qpoases[0]);
	}

	//oasesSolver.init(vec_H_qp.data(),&vec_h_qp[0],&vec_A_qp[0],0,0,&vec_lba[0],&vec_uba[0], 10);
	

}

/* Solves formulated QP problem */
void
MPC :: solveQP()
{
	
	if(USE_CASADI)
	{	
		casadi::DMDict res = qp_solver(qp_arg);

		vector<double> sol_std(res.at("x"));

		sol = MatrixXd::Map(sol_std.data(), N*(n+m+n_soft),1); 
	}
	else
	{
		qpOASES::int_t nWSR = 1000;
		qpOASES::real_t sol_qpoases[H_qp.rows()];
		qpOASES::real_t cputime = 0;
		qpOASES::real_t cputime_init = 0;
		qpOASES::returnValue retval;
		
		qpOASES::real_t tic, toc;
		tic = qpOASES::getCPUtime();
		if(!use_sparse)
			retval=oasesSolver.hotstart( H_qpoases,h_qpoases,A_qpoases,NULL,NULL,lba_qpoases,uba_qpoases, nWSR, 0 );
		else
			retval=oasesSolver.hotstart( Hs_qpoases,h_qpoases,As_qpoases,NULL,NULL,lba_qpoases,uba_qpoases, nWSR, 0 );
		
		if(qpOASES::getSimpleStatus(retval)!=0){
			oasesSolver.reset();
			if(!use_sparse)
				oasesSolver.init( H_qpoases,h_qpoases,A_qpoases,NULL,NULL,lba_qpoases,uba_qpoases, nWSR, 0 );
			else
				oasesSolver.init( Hs_qpoases,h_qpoases,As_qpoases,NULL,NULL,lba_qpoases,uba_qpoases, nWSR, 0 );
		}
		toc = qpOASES::getCPUtime();
		cout << (toc-tic)*1000 <<" ms" << endl;

		static ofstream mpcTimeLog("./data/mpcTimeLog.txt");
		mpcTimeLog << (toc-tic)*1000 << endl;


		oasesSolver.getPrimalSolution(sol_qpoases);
		sol = MatrixXd::Map(&sol_qpoases[0], N*(n+m+n_soft),1);

		
		

		


	}

	xsol=sol.block(0,0,N*n,1);
	usol=sol.block(N*n,0,N*m,1);
	esol=sol.block(N*(n+m),0,N*n_soft,1);
}




void
FDO :: initSolver(double m, double mu_wall, double dt_in)
{
	dt=dt_in;
	// ============================= FORMULATE THE PROBLEM =============================================
	// Variables (forces)
	vector<casadi::SX> f;
	for(int i=0; i<4; i++){
		f.push_back(casadi::SX::sym("f" + std::to_string(i), 3, 1));
	}

	// Parameters
	vector<casadi::SX> param, xc, ns, ts1, ts2;
	for(int i=0; i<4; i++){
		xc.push_back(casadi::SX::sym("xc" + std::to_string(i), 3, 1));	// positions of legs in core coo frame
		ns.push_back(casadi::SX::sym("ns" + std::to_string(i), 3, 1));	// surface normals
		ts1.push_back(casadi::SX::sym("ts1" + std::to_string(i), 3, 1));	// surface tangent1
		ts2.push_back(casadi::SX::sym("ts2" + std::to_string(i), 3, 1));	// surface tangent2
		
		param.push_back(xc[i]);	// positions of legs in core coo frame
		param.push_back(ns[i]);	// surface normals
		param.push_back(ts1[i]);	// surface tangent1
		param.push_back(ts2[i]);	// surface tangent2
	}
	casadi::SX fg_sym=casadi::SX::sym("fg", 3, 1);
	param.push_back(fg_sym);
	casadi::SX ddx_sym=casadi::SX::sym("ddx", 3, 1);
	param.push_back(ddx_sym);
	casadi::SX dw_sym=casadi::SX::sym("dw", 3, 1);
	param.push_back(dw_sym);

	// Objective function
	//casadi::SX costFun = dot(f[0],f[0]) + dot(f[1],f[1]) + dot(f[2],f[2]) + dot(f[3],f[3]);
	casadi::SX forces_core = f[0]+f[1]+f[2]+f[3]+m*(ddx_sym + fg_sym);
						
	casadi::SX torques_core = cross(xc[0], f[0])+cross(xc[1], f[1])+cross(xc[2], f[2])+cross(xc[3], f[3]) + dw_sym;
	
	casadi::SX costFun = 	dot(forces_core,forces_core) + 
							dot(torques_core, torques_core) +
							(dot(f[0],f[0]) + dot(f[1],f[1]) + dot(f[2],f[2]) + dot(f[3],f[3]))*0.001;

	// Constraints
	vector<casadi::SX> constr;

	// force on the core
	//constr.push_back(f[0]+f[1]+f[2]+f[3]);

	// friction cones
	mu_wall*=sqrt(0.5);
	for(int i=0; i<4; i++){
		constr.push_back(mu_wall*dot(ns[i], f[i]) + dot(ts1[i], f[i]));
		constr.push_back(mu_wall*dot(ns[i], f[i]) + dot(ts2[i], f[i]));

		constr.push_back(-mu_wall*dot(ns[i], f[i]) + dot(ts1[i], f[i]));
		constr.push_back(-mu_wall*dot(ns[i], f[i]) + dot(ts2[i], f[i]));
	}

	// NLP solver
	casadi::SXDict nlp = {{"x", vertcat(f)}, {"p", vertcat(param)}, {"f", costFun}, {"g", vertcat(constr)}};	

	casadi::Dict opts;
	opts["print_time"] = false;
	opts["ipopt.print_level"] = 0;
	//opts["verbose_init"] = false;
	//opts["verbose"] = false;
	//cout << casadi::nlpsol_options("ipopt") << endl;
	
	solver = casadi::nlpsol("solver", "ipopt", nlp, opts);


}

void
FDO :: updateParameters(MatrixXd XC, MatrixXd NS, MatrixXd maxForce, Vector3d fg, Vector3d ddx, Vector3d dw)
{
	// =============================  GET SURFACE TANGENTS =============================================
	MatrixXd TS1(3,4), TS2(3,4);
	TS1.setZero();
	TS2.setZero();
	for(int i=0; i<4; i++){
		if(abs(NS(0,i))>0.001){ // x not zero
			TS1(1,i)=1; 
			TS1(2,i)=1;
			TS1(0,i)=( -NS(1,i)*TS1(1,i) - NS(2,i)*TS1(2,i) )/NS(0,i);
		}
		else if(abs(NS(1,i))>0.001){
			TS1(0,i)=1; 
			TS1(2,i)=1;
			TS1(1,i)=( -NS(0,i)*TS1(0,i) - NS(2,i)*TS1(2,i) )/NS(1,i);
		}
		else{
			TS1(0,i)=1; 
			TS1(1,i)=1;
			TS1(2,i)=( -NS(0,i)*TS1(0,i) - NS(1,i)*TS1(1,i) )/NS(2,i);
		}
		TS1.block<3,1>(0,i)=TS1.block<3,1>(0,i)/TS1.block<3,1>(0,i).norm();

		TS2.block<3,1>(0,i)=TS1.block<3,1>(0,i).cross(NS.block<3,1>(0,i));
		TS2.block<3,1>(0,i)=TS2.block<3,1>(0,i)/TS2.block<3,1>(0,i).norm();
	}


	// ============================= UPDATE VALUES (PARAMETERS AND BOUNDS) =============================================
	// initial guess
	vector<double> x0;
	for(int i=0; i<12; i++){
		x0.push_back(0);
	}

	// Parameter values
	vector<double> p0;
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++)
			p0.push_back(XC(j,i));
		for(int j=0; j<3; j++)
			p0.push_back(NS(j,i));
		for(int j=0; j<3; j++)
			p0.push_back(TS1(j,i));
		for(int j=0; j<3; j++)
			p0.push_back(TS2(j,i));		
	}
	for(int j=0; j<3; j++)
		p0.push_back(fg(j));
	for(int j=0; j<3; j++)
		p0.push_back(ddx(j));
	for(int j=0; j<3; j++)
		p0.push_back(dw(j));


	// Variable bounds 
  	vector<double> lbx, ubx;
  	for(int i=0; i<4; i++){
  		for(int j=0; j<3; j++){
			lbx.push_back(-maxForce(i));
			ubx.push_back(maxForce(i));
  		}
  	}


  	// Constraint bounds
  	vector<double> lbg, ubg;

	// force on the core
	//for(int i=0; i<3; i++){
	//	lbg.push_back(-m*fg(i));
	//	ubg.push_back(-m*fg(i));
	//}

	// friction cones
	for(int i=0; i<4; i++){
		lbg.push_back(0);
		lbg.push_back(0);
		ubg.push_back(casadi::inf);
		ubg.push_back(casadi::inf);

		lbg.push_back(-casadi::inf);
		lbg.push_back(-casadi::inf);
		ubg.push_back(0);
		ubg.push_back(0);
	}

	arg_tmp = {
		{"x0", x0},
		{"p", p0},
		{"lbx", lbx},
        {"ubx", ubx},
        {"lbg", lbg},
        {"ubg", ubg}
    };  
}

void
FDO :: updateArgs()
{
	// Pack all arguments
	arg = arg_tmp;
}

void
FDO :: solveNLP()
{
	 res = solver(arg);
}

Matrix<double, 3, 4> 
FDO :: getSolution()
{
	casadi::DM xopt=res["x"];
	MatrixXd xoptmat(3,4);
	for(int i=0; i<4; i++){
		xoptmat.block<3,1>(0,i) << (double)xopt(0+i*3), (double)xopt(1+i*3), (double)xopt(2+i*3);
	}

	return xoptmat;
}

















void
forceDistributionIPOPT(MatrixXd XC, MatrixXd NS, MatrixXd contact, Vector3d fg, Vector3d ddx, Vector3d dw, double m, double mu)
{
	static bool is_init;
	static casadi::Function solver;
	// ============================= FORMULATE THE PROBLEM =============================================
	if(!is_init){
		is_init=true;
		// Variables (forces)
		vector<casadi::SX> f;
		for(int i=0; i<4; i++){
			f.push_back(casadi::SX::sym("f" + std::to_string(i), 3, 1));
		}

		// Parameters
		vector<casadi::SX> param, xc, ns, ts1, ts2;
		for(int i=0; i<4; i++){
			xc.push_back(casadi::SX::sym("xc" + std::to_string(i), 3, 1));	// positions of legs in core coo frame
			ns.push_back(casadi::SX::sym("ns" + std::to_string(i), 3, 1));	// surface normals
			ts1.push_back(casadi::SX::sym("ts1" + std::to_string(i), 3, 1));	// surface tangent1
			ts2.push_back(casadi::SX::sym("ts2" + std::to_string(i), 3, 1));	// surface tangent2
			
			param.push_back(xc[i]);	// positions of legs in core coo frame
			param.push_back(ns[i]);	// surface normals
			param.push_back(ts1[i]);	// surface tangent1
			param.push_back(ts2[i]);	// surface tangent2
		}
		casadi::SX fg_sym=casadi::SX::sym("fg", 3, 1);
		param.push_back(fg_sym);
		casadi::SX ddx_sym=casadi::SX::sym("ddx", 3, 1);
		param.push_back(ddx_sym);
		casadi::SX dw_sym=casadi::SX::sym("dw", 3, 1);
		param.push_back(dw_sym);

		// Objective function
		//casadi::SX costFun = dot(f[0],f[0]) + dot(f[1],f[1]) + dot(f[2],f[2]) + dot(f[3],f[3]);
		casadi::SX forces_core = f[0]+f[1]+f[2]+f[3]+m*(ddx_sym + fg_sym);
							
		casadi::SX torques_core = cross(xc[0], f[0])+cross(xc[1], f[1])+cross(xc[2], f[2])+cross(xc[3], f[3]) + dw_sym;
		
		casadi::SX costFun = 	dot(forces_core,forces_core) + 
								dot(torques_core, torques_core) +
								(dot(f[0],f[0]) + dot(f[1],f[1]) + dot(f[2],f[2]) + dot(f[3],f[3]))*0.001;

		// Constraints
		vector<casadi::SX> constr;

		// force on the core
		//constr.push_back(f[0]+f[1]+f[2]+f[3]);

		// friction cones
		mu*=sqrt(0.5);
		for(int i=0; i<4; i++){
			constr.push_back(mu*dot(ns[i], f[i]) + dot(ts1[i], f[i]));
			constr.push_back(mu*dot(ns[i], f[i]) + dot(ts2[i], f[i]));

			constr.push_back(-mu*dot(ns[i], f[i]) + dot(ts1[i], f[i]));
			constr.push_back(-mu*dot(ns[i], f[i]) + dot(ts2[i], f[i]));
		}

		// NLP solver
		casadi::SXDict nlp = {{"x", vertcat(f)}, {"p", vertcat(param)}, {"f", costFun}, {"g", vertcat(constr)}};	

		casadi::Dict opts;
		opts["print_time"] = false;
		opts["ipopt.print_level"] = 0;
		//opts["verbose_init"] = false;
		//opts["verbose"] = false;
		//cout << casadi::nlpsol_options("ipopt") << endl;
		
		solver = casadi::nlpsol("solver", "ipopt", nlp, opts);
	}

	// =============================  GET SURFACE TANGENTS =============================================
	MatrixXd TS1(3,4), TS2(3,4);
	TS1.setZero();
	TS2.setZero();
	for(int i=0; i<4; i++){
		if(abs(NS(0,i))>0.001){ // x not zero
			TS1(1,i)=1; 
			TS1(2,i)=1;
			TS1(0,i)=( -NS(1,i)*TS1(1,i) - NS(2,i)*TS1(2,i) )/NS(0,i);
		}
		else if(abs(NS(1,i))>0.001){
			TS1(0,i)=1; 
			TS1(2,i)=1;
			TS1(1,i)=( -NS(0,i)*TS1(0,i) - NS(2,i)*TS1(2,i) )/NS(1,i);
		}
		else{
			TS1(0,i)=1; 
			TS1(1,i)=1;
			TS1(2,i)=( -NS(0,i)*TS1(0,i) - NS(1,i)*TS1(1,i) )/NS(2,i);
		}
		TS1.block<3,1>(0,i)=TS1.block<3,1>(0,i)/TS1.block<3,1>(0,i).norm();

		TS2.block<3,1>(0,i)=TS1.block<3,1>(0,i).cross(NS.block<3,1>(0,i));
		TS2.block<3,1>(0,i)=TS2.block<3,1>(0,i)/TS2.block<3,1>(0,i).norm();
	}


	// ============================= UPDATE VALUES (PARAMETERS AND BOUNDS) =============================================
	// initial guess
	vector<double> x0;
	for(int i=0; i<12; i++){
		x0.push_back(0);
	}

	// Parameter values
	vector<double> p0;
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++)
			p0.push_back(XC(j,i));
		for(int j=0; j<3; j++)
			p0.push_back(NS(j,i));
		for(int j=0; j<3; j++)
			p0.push_back(TS1(j,i));
		for(int j=0; j<3; j++)
			p0.push_back(TS2(j,i));		
	}
	for(int j=0; j<3; j++)
		p0.push_back(fg(j));
	for(int j=0; j<3; j++)
		p0.push_back(ddx(j));
	for(int j=0; j<3; j++)
		p0.push_back(dw(j));


	// Variable bounds 
  	vector<double> lbx, ubx;
  	for(int i=0; i<4; i++){
  		for(int i=0; i<3; i++){
  			if(contact(i)){
  				lbx.push_back(-casadi::inf);
  				ubx.push_back(casadi::inf);
  			}
  			else{
  				lbx.push_back(0);
  				ubx.push_back(0);
  			}
  		}
  	}


  	// Constraint bounds
  	vector<double> lbg, ubg;

	// force on the core
	//for(int i=0; i<3; i++){
	//	lbg.push_back(-m*fg(i));
	//	ubg.push_back(-m*fg(i));
	//}

	// friction cones
	for(int i=0; i<4; i++){
		lbg.push_back(0);
		lbg.push_back(0);
		ubg.push_back(casadi::inf);
		ubg.push_back(casadi::inf);

		lbg.push_back(-casadi::inf);
		lbg.push_back(-casadi::inf);
		ubg.push_back(0);
		ubg.push_back(0);
	}



	// Get the optimal solution
	casadi::DMDict arg = {
				{"x0", x0},
				{"p", p0},
				{"lbx", lbx},
	            {"ubx", ubx},
	            {"lbg", lbg},
	            {"ubg", ubg}
	        };  


	casadi::DMDict res = solver(arg);
	cout << res["x"] << endl;
}
