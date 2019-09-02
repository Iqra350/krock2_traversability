#include "OptiCtrl.hpp"


#include <sstream>



#include <webots/Robot.hpp>
#include <webots/Supervisor.hpp>
#include <webots/Servo.hpp>
#include <webots/Emitter.hpp>
#include <webots/Receiver.hpp>

#define nParameters 4
#define nSettings 1
#define nFitness 3



#ifdef OPTIMIZATION




OptiCtrl :: OptiCtrl(){
	t_total=0;


	// init fitness
	fitfun[0]=99999;
	fitfun[1]=99999;
	fitfun[2]=99999;
};

/* Initializes optimization */
void
OptiCtrl :: optimizationInit(Controller *controller, PleurobotSim *pleuro){
    // get pointer to the webots opti stuff
    
    if(opti){
        // get names of parameters and settings
        for(int i=0; i<nParameters; i++)
        {
        string str="P";
            stringstream s;
                    s << "P";
                    s << i; cout<<s.str()<<endl;
            params_names.push_back(string(s.str().c_str()));
        }
        for(int i=0; i<nSettings; i++)
        {
            stringstream s("");
                    s << "S";
                    s << i;
            settings_names.push_back(s.str());
        }

        // get parameters
        Task::Parameter p;
        for(int i=0; i<params_names.size(); i++)
        {
            if (opti.Parameter(params_names[i],p))
                {
                params.push_back(p.value());
            }
            else
            {
                cerr << "optimization parameter ' " << params_names[i] << " ' not found! " << endl;
            }
        }

        // get settings
        string param;
        double val;

        for(int i=0; i<settings_names.size(); i++)
        {
            if (opti.Setting(settings_names[i], param))
                {
                    stringstream s(param);
                    s >> val;
                settings.push_back(val); cout<<settings[i]<<endl;
            }
            else
            {
                cerr << "optimization setting ' " << settings_names[i] << " ' not found! " << endl;
            }
        }

        /////////////////////////////////// modify stuff ////////////////////////////////////////////

        // GAIT PARAMETERS
        double xoff, yoff, duty, phase_off0, phase_offDiff, stride_len;
        xoff=params[0];
        yoff=params[1];
        duty=params[2];
        phase_off0=params[3];
        stride_len=controller->trPointsFL0(0,0) - controller->trPointsFL0(0,2);
        stride_len=0.345;

        // gait kinematics
        controller->trPointsFL0(0,0)=stride_len/2+xoff;
        controller->trPointsFL0(0,1)=xoff;
        controller->trPointsFL0(0,2)=-stride_len/2+xoff;

        controller->trPointsFR0(0,0)=stride_len/2+xoff;
        controller->trPointsFR0(0,1)=xoff;
        controller->trPointsFR0(0,2)=-stride_len/2+xoff;

        controller->trPointsHL0(0,0)=stride_len/2+yoff;
        controller->trPointsHL0(0,1)=yoff;
        controller->trPointsHL0(0,2)=-stride_len/2+yoff;

        controller->trPointsHR0(0,0)=stride_len/2+yoff;
        controller->trPointsHR0(0,1)=yoff;
        controller->trPointsHR0(0,2)=-stride_len/2+yoff;

        controller->trPointsFL=controller->trPointsFL0;
        controller->trPointsFR=controller->trPointsFR0;
        controller->trPointsHL=controller->trPointsHL0;
        controller->trPointsHR=controller->trPointsHR0;

        // duty
        controller->Duty(0)=duty;
        controller->Duty(1)=duty;

        // phases
        
        double front_hind_footplacement_offset=0.16;
        double IG=-0.4368;



        phase_offDiff=duty/stride_len*( front_hind_footplacement_offset + IG - (xoff-yoff) );
        

        if(phase_offDiff<-1)
            phase_offDiff+=1;
        else if (phase_offDiff>1)
            phase_offDiff-=1;
        





        double phShifts[4];
        phShifts[0]=phase_off0;
        phShifts[1]=phShifts[0]+0.5;
        phShifts[2]=phase_off0+phase_offDiff;
        phShifts[3]=phShifts[2]+0.5;



        for(int i=0; i<4; i++){
            //controller->legPhase(i)=phShifts[i]<0?phShifts[i]+1:phShifts[i];
            controller->legPhase(i)=fmod(phShifts[i],1);
            controller->phShifts(i)=phShifts[i];
        }
        

        cout <<  "stride_len= "<< stride_len << endl;
        cout <<  "xoff= "<< xoff << endl;
        cout <<  "yoff= "<< yoff << endl;
        cout <<  "duty= "<< duty << endl;
        cout <<  "phd0= "<< phase_off0 << endl;
        cout <<  "phdiff= "<< phase_offDiff << endl;
        //cout <<  "delta= "<< -IG-stride_len+xoff-yoff+duty*stride_len*phase_offDiff << endl;
        cout <<  "delta= "<< phase_offDiff/duty*stride_len + (xoff-yoff) - IG << endl;
        cout << "leg=" << controller->trPointsFL0<<endl;
/*
        double front_hind_footplacement_offset=0.16;
        double IG=-0.4368;
        phase_offDiff=duty/stride_len*( (xoff-yoff) - IG - front_hind_footplacement_offset);
*/

        ///////////////////////////////////////////////////////////////////////////////////////////
    }
}


/* Ends optimization */
void
OptiCtrl :: optimizationEnd(PleurobotSim *pleuro){
    vector<double> fit;
    vector<string> fitness_names;
    map<string, double> fitness;

    // fitness variables
    for(int i=0; i<nFitness; i++){
    	fit.push_back(fitfun[i]);
    }
    for(int i=0; i<nFitness; i++)
    {
        stringstream s("");
                s << "F";
                s << i;
        fitness_names.push_back(s.str());
    }


    // pack them
    for(int i=0; i<fit.size(); i++)
	{
		fitness[fitness_names[i]] = fit[i];
	}

    // Send the actual response
	if (opti.Setting("optiextractor"))
	{
		cout << "Not sending response, because in optiextractor mode! Enjoy the show!" << endl;
	}
	else
	{
		opti.Respond(fitness);
	}




    // cout if not in optimization
	if(opti.Setting("optiextractor") || !opti)
	{
	//	cout << "Fitness: " << fitfun[0] <<"\t"<< fitfun[1] <<"\t"<< fitfun[2] <<"\t"<< fitfun[0]*fitfun[1]*fitfun[2]<< endl;
	}



	// quit simulation
	if (opti && !opti.Setting("optiextractor"))
		pleuro->killSimulation();

}

/* Checks stopping criterion */
int
OptiCtrl :: optimizationShouldWeStopIt(double timestep){
    static double t_stop=0;
    t_stop+=timestep;
    //optimization::Webots &opti = optimization::Webots::Instance();
    int crit=0;
    if(opti){
        if(t_stop>settings[0]){
            crit=1;
        }

    }

    return crit;
}

/* Updates internal time and runs state machine */
int
OptiCtrl :: OptiStep(double dt, Controller *controller, PleurobotSim *pleuro){


	static int terrains=0;
	static int reruns=0;
	static double t_run=0;
	static double maxReruns=3, init_pos[3], final_pos[3];
	//timings
	static double T1max=34, T2max=34, T3max=34, t_timeout=4;
	bool flipped, init_pos_set=false;
	t_total+=dt;
	t_run+=dt;
	


	// get initial position
	if(t_run>t_timeout && !init_pos_set){
		init_pos_set=true;
		pleuro->GetPosition(init_pos, init_pos);
	}


	static double roll, pitch, yaw, Iroll=0, Ipitch=0;

	roll=atan2(controller->rotMatFgird(2,1), controller->rotMatFgird(1,1));
	pitch=atan2(controller->rotMatFgird(0,2), controller->rotMatFgird(2,2));
	yaw=atan2(controller->rotMatFgird(1,0), controller->rotMatFgird(0,0));


	if(abs(roll)>1.5 || abs(pitch)>1.5){
        
		//flipped=1;
	}

	Iroll+=roll*roll*dt;
	Ipitch+=pitch*pitch*dt;


	// normalize to flat terrain (vel of 0.1750m/s)
	if(t_run>settings[0] || flipped){
		if(!flipped){
			pleuro->GetPosition(final_pos, final_pos);
			fitfun[0]=final_pos[0]-init_pos[0];	
			fitfun[1]=-Iroll;
			fitfun[2]=-Ipitch;	
		}
		else{
			flipped=0;
			fitfun[0]=-999999;
			fitfun[1]=-999999;
			fitfun[2]=-999999;
		}
		optimizationEnd(pleuro);
	}

	return 1;
}













































#endif
