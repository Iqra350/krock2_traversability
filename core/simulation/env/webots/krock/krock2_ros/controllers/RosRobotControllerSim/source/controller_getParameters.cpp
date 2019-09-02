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

   
gaitParams 
readGaitParameters(stringstream& stringstream_file){
    int i, j;
    gaitParams GPtmp;
    //spineCPGscaling
    stringstream_file>>GPtmp.spineCPGscaling(0);
    stringstream_file>>GPtmp.spineCPGscaling(1);

     // duty cycles for legs (front and hind)
    for(i=0; i<2; stringstream_file>>GPtmp.Duty(i), i++);      
    cout << "Duty\t" << GPtmp.Duty << endl << endl;
    // phase shift between legs
    for(i=0; i<4; stringstream_file>>GPtmp.phShifts(i), i++);  
    cout << "phShifts\t" << GPtmp.phShifts << endl << endl;
    // trajectories MID STANCE
    for(i=0; i<3; i++){
        stringstream_file>>GPtmp.midStance(i, 0);
    }
    GPtmp.midStance.block<3,1>(0,1)=GPtmp.midStance.block<3,1>(0,0); GPtmp.midStance(1,1)*=-1;
    for(i=0; i<3; i++){
        stringstream_file>>GPtmp.midStance(i, 2);
    }  
    GPtmp.midStance.block<3,1>(0,3)=GPtmp.midStance.block<3,1>(0,2); GPtmp.midStance(1,3)*=-1;
    cout << "midStance\t" << GPtmp.midStance << endl << endl;

    // ellipse props
    stringstream_file >> GPtmp.ellipse_a(0); 
    GPtmp.ellipse_a(1)=GPtmp.ellipse_a(0);
    GPtmp.ellipse_a(2)=GPtmp.ellipse_a(0);
    GPtmp.ellipse_a(3)=GPtmp.ellipse_a(0);


    stringstream_file >> GPtmp.ellipse_b(0); 
    GPtmp.ellipse_b(1)=GPtmp.ellipse_b(0);
    GPtmp.ellipse_b(2)=GPtmp.ellipse_b(0);
    GPtmp.ellipse_b(3)=GPtmp.ellipse_b(0);

    cout << "ellipse_a\t" << GPtmp.ellipse_a << endl << endl;
    cout << "ellipse_b\t" << GPtmp.ellipse_b << endl << endl;

    //cout << "============================MIDSTANCE ====================================" << endl; 
    //cout << midStance << endl;
    //cout << "============================ELLIPSE ====================================" << endl; 
    //cout << ellipse_a << " b: "<<ellipse_b << endl;
    // swing height FRONT
    stringstream_file >> GPtmp.swing_height(0);
    GPtmp.swing_height(1)=GPtmp.swing_height(0); 
    // swing height HIND
    stringstream_file >> GPtmp.swing_height(2); 
    GPtmp.swing_height(3)=GPtmp.swing_height(2);

    // swing width FRONT
    stringstream_file >> GPtmp.swing_width(0);
    GPtmp.swing_width(1)=GPtmp.swing_width(0); 
    // swing width HIND
    stringstream_file >> GPtmp.swing_width(2); 
    GPtmp.swing_width(3)=GPtmp.swing_width(2);

    cout << "swing_height\t" << GPtmp.swing_height << endl << endl;
    cout << "swing_width\t" << GPtmp.swing_width << endl << endl;
    // nLO FRONT
    for(i=0; i<3; i++){
        stringstream_file >> GPtmp.nLO(i,0);
    }
    GPtmp.nLO.block<3,1>(0,0)=GPtmp.nLO.block<3,1>(0,0)/GPtmp.nLO.block<3,1>(0,0).norm();
    GPtmp.nLO.block<3,1>(0,1)=GPtmp.nLO.block<3,1>(0,0);
    GPtmp.nLO(1,1)*=-1;


    // nLO HIND
    for(i=0; i<3; i++){
        stringstream_file >> GPtmp.nLO(i,2);
    }
    GPtmp.nLO.block<3,1>(0,2)=GPtmp.nLO.block<3,1>(0,2)/GPtmp.nLO.block<3,1>(0,2).norm();
    GPtmp.nLO.block<3,1>(0,3)=GPtmp.nLO.block<3,1>(0,2);
    GPtmp.nLO(1,3)*=-1;

    cout << "nLO\t" << GPtmp.nLO << endl << endl;
     // nTD FRONT
    for(i=0; i<3; i++){
        stringstream_file >> GPtmp.nTD(i,0);
    }
    GPtmp.nTD.block<3,1>(0,0)=GPtmp.nTD.block<3,1>(0,0)/GPtmp.nTD.block<3,1>(0,0).norm();
    GPtmp.nTD.block<3,1>(0,1)=GPtmp.nTD.block<3,1>(0,0);
    GPtmp.nTD(1,1)*=-1;


    // nTD HIND
    for(i=0; i<3; i++){
        stringstream_file >> GPtmp.nTD(i,2);
    }
    GPtmp.nTD.block<3,1>(0,2)=GPtmp.nTD.block<3,1>(0,2)/GPtmp.nTD.block<3,1>(0,2).norm();
    GPtmp.nTD.block<3,1>(0,3)=GPtmp.nTD.block<3,1>(0,2);
    GPtmp.nTD(1,3)*=-1;

    cout << "nTD\t" << GPtmp.nTD << endl << endl;

    // some stupid angles for splines
    for(i=0; i<3; i++){
        stringstream_file>>GPtmp.bezierScaling(i,0);
        stringstream_file>>GPtmp.bezierScaling(i,1);
    }                                                    

    cout << "bezierScaling\t" << GPtmp.bezierScaling << endl << endl;

    //tSclSwing
    stringstream_file >> GPtmp.tSclSwing;

    cout << "tSclSwing\t" << GPtmp.tSclSwing << endl << endl;
    //nSurf
    for(i=0; i<3; i++){
        for(j=0; j<4; j++){
            stringstream_file >> GPtmp.nSurf(i,j);
        }
    }
    for(i=0; i<4; i++){
        GPtmp.nSurf.block<3,1>(0,i)/=GPtmp.nSurf.block<3,1>(0,i).norm();
    }

    cout << "nSurf\t" << GPtmp.nSurf << endl << endl;
    //qNULL
    for(i=0; i<4; i++){
        for(j=0; j<4; j++){
            stringstream_file >> GPtmp.qNULL(i,j);
        }
    }
    GPtmp.qNULL*=my_pi/180.;
    cout << "qNULL\t" << GPtmp.qNULL << endl << endl;
    return GPtmp;

}




/* Reads parameters from parameters_ikin.txt file */
bool
Controller :: getParameters()
{
    int i, j;
    stringstream stringstream_file;
    ifstream file_kinematics;
    ifstream file_ikinController;
    ifstream file_reflex_posture;
    ifstream file_cpg;
    ifstream file_joystick;
    ifstream file_contactForceControllerConfig;
    ifstream file_masses;
    ifstream file_torquesCtrl;
    ifstream file_joysticRec;
    
    int numoflines;
    



    file_kinematics.open("config/kinematics.config");
    file_ikinController.open("config/ikinController.config");
    file_torquesCtrl.open("config/torquesCtrl.config");
    file_joysticRec.open("data/joysticRec.txt");
    file_reflex_posture.open("config/reflex_posture.config");
    file_cpg.open("config/cpg.config");
    file_joystick.open("config/joystick.config");
    file_contactForceControllerConfig.open("config/contactForceControllerConfig.config");
        
    file_masses.open("config/masses.config");
    if(file_masses.is_open()){
        cout << "reading file_masses" << endl;
        stringstream_file.str(string());
        readFileWithLineSkipping(file_masses, stringstream_file);
        for(i=0; i<22; i++){
            stringstream_file >> masses(0, i);  // mass
            stringstream_file >> masses(1, i);  // coordinates
            stringstream_file >> masses(2, i);
            stringstream_file >> masses(3, i);
            //cout<<i<<endl;
        }
        cout << "MASSES: " << endl;
        cout<<masses.transpose()<<endl;
    }

    if(file_reflex_posture.is_open()){
        cout << "reading file_reflex_posture" << endl;
        stringstream_file.str(string());
        readFileWithLineSkipping(file_reflex_posture, stringstream_file);
        
        // stumble
        stringstream_file >> reflON(0);
        stringstream_file >> reflON(1);
        stringstream_file >> stumbleForceLimit;
        stringstream_file >> stumbleTimeout;
        stringstream_file >> stumblePhaseLimits(0);
        stringstream_file >> stumblePhaseLimits(1);
        stringstream_file >> stumbleImpulse(0);
        stringstream_file >> stumbleImpulse(1);
        stringstream_file >> stumbleImpulse(2);
        stringstream_file >> stumbleFilterConstant;

        // extend
        stringstream_file >> extendForceLimit;
        stringstream_file >> extendTimeout;
        stringstream_file >> extendPhaseLimits(0);
        stringstream_file >> extendPhaseLimits(1);
        stringstream_file >> extendImpulse(0);
        stringstream_file >> extendImpulse(1);
        stringstream_file >> extendImpulse(2);
        stringstream_file >> extendFilterConstant;

        
    }


    if(file_kinematics.is_open()) {
        cout << "reading file_kinematics" << endl;
        stringstream_file.str(std::string());
        readFileWithLineSkipping(file_kinematics, stringstream_file);

        //spine
        for(i=0; i<Ntrunk+1; i++){    //cout<<"spine"<<spine_kin.transpose()<<endl;
            stringstream_file>>trunk_kin(i);
        }

        for(i=0; i<Nneck; i++){    //cout<<"spine"<<spine_kin.transpose()<<endl;
            stringstream_file>>neck_kin(i);
        }

        for(i=0; i<Ntail; i++){    //cout<<"spine"<<spine_kin.transpose()<<endl;
            stringstream_file>>tail_kin(i);
        }

        // IG
        IG=0;
        for(int i=0; i<Ntrunk+1; i++){
            IG+=trunk_kin(i);
        }
        cout<<"trunk_kin\n"<<trunk_kin<<endl;
        // legs
        for(i=0;i<Nlegs+1;i++){
            for(j=0;j<3;j++){
                stringstream_file>>FL_kin(i, j); 
            }
        }    
        cout<<"FL_kin\n"<<FL_kin<<endl;
        for(i=0;i<Nlegs+1;i++){
            for(j=0;j<3;j++){
                stringstream_file>>FR_kin(i, j);
            }
        }   
        cout<<"FR_kin\n"<<FR_kin<<endl;
        for(i=0;i<Nlegs+1;i++){
            for(j=0;j<3;j++){
                stringstream_file>>HL_kin(i, j);
            }
        }   
        cout<<"HL_kin\n"<<HL_kin<<endl;
        for(i=0;i<Nlegs+1;i++){
            for(j=0;j<3;j++){
                stringstream_file>>HR_kin(i, j);
            }
        }       
        cout<<"HR_kin\n"<<HR_kin<<endl;
        for(i=0; i<NUM_MOTORS; stringstream_file>>angSignsCorrIkin2Webots[i], i++);
        for(i=0; i<NUM_MOTORS; stringstream_file>>angShiftsCorrIkin2Webots[i], i++);

        for(i=0; i<NUM_MOTORS; stringstream_file>>angSignsCorrIkin2Robot[i], i++);
        for(i=0; i<NUM_MOTORS; stringstream_file>>angShiftsCorrIkin2Robot[i], i++);

    }
    
    

    if(file_ikinController.is_open()) {
        cout << "reading file_ikinController" << endl;
        stringstream_file.str(string());
        readFileWithLineSkipping(file_ikinController, stringstream_file);
        //speed_or_frequency
        stringstream_file>>speed_or_frequency;

        //maxSpeed
        stringstream_file>>maxSpeed;

        //maxFrequency
        stringstream_file>>maxFrequency;

        stringstream_file>>T_trans0;
        T_trans=T_trans0;

        stringstream_file>>Tfilt_angles;

        // trajectory filtering
        stringstream_file>>Tf1;


        // READ GAIT PARAMETERS
        // WALKING
        WP = readGaitParameters(stringstream_file);

        WPlow = readGaitParameters(stringstream_file);

        // CRAWLING
        int crawlingSequence_N;
        stringstream_file>>crawlingSequence_N;

        crawlingSequence.resize(crawlingSequence_N,4);
        for(int i=0; i<crawlingSequence_N; i++){
            stringstream_file>>crawlingSequence(i,0);
            stringstream_file>>crawlingSequence(i,1);
            stringstream_file>>crawlingSequence(i,2);
            stringstream_file>>crawlingSequence(i,3);
        }
        CP = readGaitParameters(stringstream_file);



        //================================ ikin related ================================================
        // constraints
        for(i=0; i<4; stringstream_file>>constrFL(0, i), constrFL(0, i)*=my_pi/180., i++);
        for(i=0; i<4; stringstream_file>>constrFL(1, i), constrFL(1, i)*=my_pi/180., i++);   //cout<<constrFL.transpose()<<endl;
        for(i=0; i<4; stringstream_file>>constrFR(0, i), constrFR(0, i)*=my_pi/180., i++);
        for(i=0; i<4; stringstream_file>>constrFR(1, i), constrFR(1, i)*=my_pi/180., i++);   //cout<<constrFR.transpose()<<endl;
        for(i=0; i<4; stringstream_file>>constrHL(0, i), constrHL(0, i)*=my_pi/180., i++);
        for(i=0; i<4; stringstream_file>>constrHL(1, i), constrHL(1, i)*=my_pi/180., i++);   //cout<<constrHL.transpose()<<endl;
        for(i=0; i<4; stringstream_file>>constrHR(0, i), constrHR(0, i)*=my_pi/180., i++);
        for(i=0; i<4; stringstream_file>>constrHR(1, i), constrHR(1, i)*=my_pi/180., i++);   //cout<<constrHR.transpose()<<endl;
        stringstream_file >> constrS;
        constrS*=my_pi/180.;
        
        for(i=0; i<3; i++){
            stringstream_file >> workspaceEllipsoid(i,0);
        }
        for(i=0; i<3; i++){
            stringstream_file >> workspaceEllipsoid(i,2);
        }
        workspaceEllipsoid.block<3,1>(0,1)=workspaceEllipsoid.block<3,1>(0,0);
        workspaceEllipsoid.block<3,1>(0,3)=workspaceEllipsoid.block<3,1>(0,2);

        for(i=0; i<3; i++){
            stringstream_file >> workspaceEllipsoidCenter(i,0);
        }
        for(i=0; i<3; i++){
            stringstream_file >> workspaceEllipsoidCenter(i,2);
        }
        workspaceEllipsoidCenter.block<3,1>(0,1)=workspaceEllipsoidCenter.block<3,1>(0,0);
        workspaceEllipsoidCenter.block<3,1>(0,3)=workspaceEllipsoidCenter.block<3,1>(0,2);
        workspaceEllipsoidCenter(1,1)*=-1;
        workspaceEllipsoidCenter(1,3)*=-1;

        // lam, M, Cinv and max_dist for iKinDLS

        stringstream_file >> lamF(0);
        stringstream_file >> lamF(1);
        stringstream_file >> lamH(0);
        stringstream_file >> lamH(1);

        for(i=0; i<4; stringstream_file>>MF(0, i), i++);
        for(i=0; i<4; stringstream_file>>MF(1, i), i++);
        for(i=0; i<4; stringstream_file>>MF(2, i), i++);
        for(i=0; i<4; stringstream_file>>MF(3, i), i++);   //cout<<MF<<endl;
        for(i=0; i<4; stringstream_file>>MH(0, i), i++);
        for(i=0; i<4; stringstream_file>>MH(1, i), i++);
        for(i=0; i<4; stringstream_file>>MH(2, i), i++);
        for(i=0; i<4; stringstream_file>>MH(3, i), i++);   //cout<<MH<<endl;


        CinvF=lamF(0)/2.*Matrix4d::Identity()+MF.transpose()*MF;
        CinvF=CinvF.inverse().eval();

        CinvH=lamH(0)/2.*Matrix4d::Identity()+MH.transpose()*MH;
        CinvH=CinvH.inverse().eval();


        stringstream_file >> max_dist;

        stringstream_file >> ikin_tol; //cout<<ikin_tol<<endl;

        stringstream_file >> ikin_maxIter; //cout<<ikin_maxIter<<endl;

        stringstream_file >> ikin_constr_penalty(0);
        stringstream_file >> ikin_constr_penalty(1);
        stringstream_file >> ikin_constr_penalty(2);

        for(i=0; i<Nneck+Ntrunk+Ntail; stringstream_file>>amps_swim(i), i++);    //cout<<amps_swim.transpose()<<endl;

        for(i=0; i<Nneck+Ntrunk+Ntail; stringstream_file>>phases_swim(i), i++);    //cout<<phases_swim.transpose()<<endl;


        



        for(i=0; i<4*Nlegs; stringstream_file>>swim_legs_pos(i), i++);
        cout<<swim_legs_pos<<endl;


    }
    
    
    if(file_contactForceControllerConfig.is_open()){ //cout<<"FILE_FORCES_OPENED"<<endl;
        cout << "reading file_forces" << endl;
        stringstream_file.str(string());
        readFileWithLineSkipping(file_contactForceControllerConfig, stringstream_file);
        stringstream_file >> mu_wall;
        stringstream_file >> mu_ground;
        stringstream_file >> mu_conservative_scaling;
        mu_wall = mu_wall * (1-mu_conservative_scaling/100.);
        mu_ground = mu_ground * (1+mu_conservative_scaling/100.);


        stringstream_file >> forceTreshold;
        stringstream_file >> maxContactForce;
        stringstream_file >> fdo_dt;
        //======================= FILTERS ======================================
        stringstream_file >> rawForcesFilter;
        stringstream_file >> leg_loading_filter;
        stringstream_file >> leg_approach_speed;

        //==================== CONTACT FORCES PID =================================
        stringstream_file >> contactF_pid_Kp;
        stringstream_file >> contactF_pid_Ki;
        stringstream_file >> contactF_pid_Kd;
        stringstream_file >> contactF_pid_Td;

        stringstream_file >> contactF_pid_limit;
    }



    if(file_joystick.is_open()){
        cout << "reading file_joystick" << endl;
        stringstream_file.str(string());
        readFileWithLineSkipping(file_joystick, stringstream_file);
        stringstream_file >> stick_r_lim;
        stringstream_file >> joy_walk_max_freq;
        stringstream_file >> joy_walk_min_freq;
        stringstream_file >> joy_walk_speed_change_filt;


        stringstream_file >> disable_crab_walk_lim;
        stringstream_file >> ellipse_small_axis;
        stringstream_file >> posing_joy_y1_rate;
        stringstream_file >> posing_joy_x1_rate;
        stringstream_file >> posing_joy_y2_rate;
        stringstream_file >> posing_head_rate;
        stringstream_file >> posing_head_limit;
        stringstream_file >> joy_swim_max_freq;
        stringstream_file >> joy_swim_speed_change_filt;
        stringstream_file >> joy_swim_max_offset;
        stringstream_file >> phase_lag_change_rate;
        stringstream_file >> phase_lag_min;
        stringstream_file >> phase_lag_max;
        stringstream_file >> joy_swim_turning_dead_zone; joy_swim_turning_dead_zone*=my_pi/180;

    }



    
    







    if(file_joysticRec.is_open()) {
        cout << "reading file_joysticRec" << endl;
        numoflines=readFileWithLineSkipping(file_joysticRec, stringstream_file);
        //cout << "NUMOFLINES: " << numoflines << endl;
        sizeOfJoystickRecording=numoflines;
        if(numoflines>0){
            joysticRecordings=new double[sizeOfJoystickRecording*4];
            for(i=0; i<sizeOfJoystickRecording; i++){
                for(j=0; j<4; j++){
                    stringstream_file >> joysticRecordings[i*4+j];
                }
            }
        }
        
    }

    return true;
    
}
