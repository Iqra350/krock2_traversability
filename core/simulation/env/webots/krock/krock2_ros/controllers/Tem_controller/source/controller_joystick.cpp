 #include "controller.hpp" 


/* Operators on gaitParameters */
gaitParams
combineGaitParameters(gaitParams *gp1, gaitParams *gp2, double x)
{
    gaitParams gp=*gp1;

    gp.spineCPGscaling =    (1-x)*gp1->spineCPGscaling +      x*gp2->spineCPGscaling;
    gp.Duty =               (1-x)*gp1->Duty +                 x*gp2->Duty;
    gp.phShifts =           (1-x)*gp1->phShifts +             x*gp2->phShifts;
    gp.midStance =          (1-x)*gp1->midStance +            x*gp2->midStance;
    gp.ellipse_a =          (1-x)*gp1->ellipse_a +            x*gp2->ellipse_a;
    gp.ellipse_b =          (1-x)*gp1->ellipse_b +            x*gp2->ellipse_b;
    gp.swing_height =       (1-x)*gp1->swing_height +         x*gp2->swing_height;
    gp.swing_width =        (1-x)*gp1->swing_width +          x*gp2->swing_width;
    gp.nSurf =              (1-x)*gp1->nSurf +                x*gp2->nSurf;
    gp.nLO =                (1-x)*gp1->nLO +                  x*gp2->nLO;
    gp.nTD =                (1-x)*gp1->nTD +                  x*gp2->nTD;
    gp.bezierScaling =      (1-x)*gp1->bezierScaling +        x*gp2->bezierScaling;
    gp.tSclSwing =          (1-x)*gp1->tSclSwing +            x*gp2->tSclSwing;
    gp.qNULL =              (1-x)*gp1->qNULL +                x*gp2->qNULL;




    return gp;



}

/* Reads joystick */
void
Controller :: readJoystick()
{
    js.update();

    // PS3 CONTROLLER
    joy_sel=js.buttons[0];
    joy_l3=js.buttons[1];
    joy_r3=js.buttons[2];
    joy_start=js.buttons[3];
    joy_l2=js.buttons[8];
    joy_r2=js.buttons[9];
    joy_l1=js.buttons[10];
    joy_r1=js.buttons[11];

    joy_aU=js.buttons[4];
    joy_aR=js.buttons[5];
    joy_aD=js.buttons[6];
    joy_aL=js.buttons[7];

    joy_bU=js.buttons[12];
    joy_bR=js.buttons[13];
    joy_bD=js.buttons[14];
    joy_bL=js.buttons[15];

    joy_x1=js.axes[0];
    joy_y1=-js.axes[1];
    joy_x2=js.axes[2];
    joy_y2=-js.axes[3];
    joy_x3=js.axes[4];
    joy_y3=-js.axes[5];


    if(AUTO_RUN){
        joystickAuto();
    }
    if(JOYSTICK_RECORDING_REPLAYING==1){
        joystickRecord();
    }
    else if(JOYSTICK_RECORDING_REPLAYING==2){
        joystickReplay();
    }


    joy_lsr=sqrt(joy_x1*joy_x1+joy_y1*joy_y1);
    joy_lsr=joy_lsr>1?1:joy_lsr;
    if(joy_lsr<0.1){
        joy_lsr=0;
    }
    if(joy_lsr>0.2){
        joy_lsphi=atan2(-joy_x1, joy_y1);
    }
    else{
        joy_lsphi=0;
    }


    joy_rsr=sqrt(joy_x2*joy_x2+joy_y2*joy_y2);
    joy_rsr=joy_rsr>1?1:joy_rsr;
    if(joy_rsr<0.1){
        joy_rsr=0;
    }
    if(joy_rsr>0.2){
        joy_rsphi=atan2(joy_y2, joy_x2)-my_pi/2;
    }
    else{
        joy_rsphi=0;
    }
    /*for(int i=0; i<16; i++){
        cout << "B" << i << ":" << js.buttons[i] << "\t\t";
    }
    cout << endl;*/

}

/* Reads joystick inputs and modifies trajectories */
bool
Controller :: updateState()
{

    if(USE_JOYSTICK){
        readJoystick();
    }
    else{
        joy_bU=0;
        joy_bR=0;
        joy_bD=0;
        joy_bL=0;
        joy_l2=0;
        joy_r2=0;
        joy_l1=0;
        joy_r1=0;
        joy_sel=0;
        joy_start=0;
        joy_l3=0;
        joy_r3=0;

        joy_bU=0;
        joy_bR=0;
        joy_bD=0;
        joy_bL=0;

        joy_x1=0;
        joy_y1=0;
        joy_x2=0;
        joy_y2=0;
        joy_x3=0;
        joy_y3=0;

        joy_lsr=1;
        joy_rsr=1;
        joy_lsphi=0;
        joy_rsphi=0;
    }

    phases=phases0;
    //============================== EMERGENCY STOP ==============================
    if(joy_start == 1){
        return false;
    }

    static double gpx = 0, gpx_old=0;
    if(state == WALKING || state == STANDING || state == POSING){
        gpx = gpx + (joy_l2 - joy_r2)*dt/2.;
        gpx = gpx<0 ? 0 : gpx;
        gpx = gpx>1 ? 1 : gpx;

        if(gpx_old!=gpx){
            GP=combineGaitParameters(&WP, &WPlow, gpx);
            gpx_old=gpx;
        }
    }
    
    
    //=============================== POSING ===================================
    if(joy_sel==0 && joy_l1==1 && state != POSING){
        state=POSING;
        forTraj_posing0 = forTraj;
        posing_head=0;
        //T_trans=T_trans0;
        T_stand = T_trans0/2;
    }
    if(joy_l1==0 && state==POSING){
        state=STANDING;
        T_stand = T_trans0/2;
    }

    //=============================== CRAWLING ===================================
    static int joy_bL_old=0;

    if(joy_sel==0 && joy_bL==1 && joy_bL_old == 0 && state == CRAWLING){
        crawling_continue = true;
    }
    joy_bL_old = joy_bL;
    if(joy_sel==0 && joy_bL==1 && state != CRAWLING){
        state=CRAWLING;
        GP=CP;
        T_trans=T_trans0;
        crawling_continue = false;
        crstate = initial;
    }
    
    
    //=============================== STANDING ===================================
    else if(joy_sel==0 && joy_bR==1 && !(state==WALKING || state == STANDING)){
        state=STANDING;
        //GP=WP;
        GP=combineGaitParameters(&WP, &WPlow, gpx);
        T_trans=T_trans0;
        T_stand = T_trans0/2;
    }

    //=============================== WALKING ===================================
    if(((state == STANDING && (   joy_lsr>0.2  || (joy_rsr*joy_r1)>0.2 )    ) || !USE_JOYSTICK)){
        state=WALKING;
        //GP=WP;
        GP=combineGaitParameters(&WP, &WPlow, gpx);
        T_stand = T_trans0/2;
    }

    //=============================== STANDING ===================================
    else if(state == WALKING && !(joy_lsr>0.2 || (joy_rsr*joy_r1)>0.2 )){
        state=STANDING;
        //GP=WP;
        GP=combineGaitParameters(&WP, &WPlow, gpx);
        T_stand = T_trans0/2;
    }

    //=============================== INITIAL ===================================
    if(joy_sel==0 && joy_bD==1 && !(state==INITIAL)){
        state=INITIAL;
        T_trans=T_trans0;
    }

    //=============================== SWIMMING ===================================
    if(joy_sel==1 && joy_bR==1 && state!=SWIMMING){
        state=SWIMMING;
        T_trans=T_trans0;
        freq_swim=0;
        cpg_offset=0;
    }





    //=============================== DECAY TRANSITION FILTERING CONSTANT ===================================
    T_trans=pt1(0, T_trans, 1, dt);
    T_stand=pt1(0, T_stand, 1, dt);

    //=============================== JOYSTICK INTERACTION ===================================
    if(AUTO_RUN){
        joy_lsr=1;
        if(t<0.1){
            state=AUTO_RUN_STATE;   
        }
        
    }

    
    return true;
}

/* Joystick manipulation of different parameters/variables */
void
Controller :: posingManipulation()
{   

    double centerLine = 0.5*(forTraj(2,1) + forTraj(5,1));
    double centerLine2 = 0.5*(forTraj_posing0(2,1) + forTraj_posing0(5,1));

    // rotating body
    double posing_rotational_velocity=0;
    posing_rotational_velocity = (-joy_x2*0.2 + forTraj_posing0(5,1) - forTraj(5,1))*3;

    // front girdle position
    double fgirdArc = (-joy_x1*0.4 + forTraj(5,1) - forTraj(2,1))*2;

    // forward - backward
    Vector3d transVel = MatrixXd::Zero(3,1);
    Vector3d proj0, proj;
    proj0 = AngleAxisd(-forTraj(5,1), Vector3d::UnitZ()) * forTraj_posing0.block<3,1>(3,1);
    proj =  AngleAxisd(-forTraj(5,1), Vector3d::UnitZ()) * forTraj.block<3,1>(3,1);


    transVel(0) = (joy_y1*0.1 + proj0(0) - proj(0))*2;
    transVel = AngleAxisd(forTraj(5,1)-centerLine, Vector3d::UnitZ()) * transVel;

    // front girdle up - down
    static double posingFgirdHeightMod = 0;
    static double normalHeight = GP.midStance(2,0);
    //if(posingFgirdHeightMod == 0){
    //    normalHeight = 0.5*(feetReference(2,0) + feetReference(2,1));
    //}
    posingFgirdHeightMod = -joy_y2*0.1;


    posingFeetHeight = normalHeight + posingFgirdHeightMod;
    


    // double v, double w, Vector3d V, double W, double fgirdArc, Vector2d spineCPGscaling
    girdleTrajectories(0, 0, transVel, posing_rotational_velocity, fgirdArc, MatrixXd::Zero(2,1)); 
    girdleVelocities();


}

void
Controller :: joystickRecord()
{
	static ofstream joysticRec("./data/joysticRec.txt");
	joysticRec << joy_x1 <<"\t" << joy_y1 << "\t";
	joysticRec << joy_x2 <<"\t" << joy_y2 << endl;

}

void
Controller :: joystickReplay()
{	
	static int linecount=0;
	if(linecount>=sizeOfJoystickRecording)
		return;
	else{
		joy_x1=joysticRecordings[4*linecount+0];
		joy_y1=joysticRecordings[4*linecount+1];
		joy_x2=joysticRecordings[4*linecount+2];
		joy_y2=joysticRecordings[4*linecount+3];
		linecount++;
	}
}

void 
Controller :: joystickAuto(){
    
    static bool is_init=false, is_enabled=false, is_oscillating=false;
    static int num_commands;
    static double in_array[10][5];
    static double joy_freq;

    // ------------------- INITIALIZE --------------------------
    if(is_init==false){
        is_init=true;
        stringstream stringstream_file;
        ifstream file_joystickAuto;

        // init reading
        file_joystickAuto.open("config/joystickAuto.config");
        stringstream_file.str(std::string());
        readFileWithLineSkipping(file_joystickAuto, stringstream_file);

        // read stuff
        stringstream_file >> is_enabled;
        stringstream_file >> num_commands;


        for(int i=0; i<num_commands; i++){
            for(int j=0; j<5; j++){
                stringstream_file >> in_array[i][j];
            }
        }

        stringstream_file >> is_oscillating;
        stringstream_file >> joy_freq;

    }

    // ------------------- RUN --------------------------

    if(is_enabled){
        if(t<0.1){
            state=AUTO_RUN_STATE;
        }
        int seg=0;
        for(int i=0; i<num_commands; i++){
            if(in_array[i][0] <= t){
                seg=i;
            }
        }




        joy_x1=in_array[seg][1];
        joy_y1=in_array[seg][2];
        joy_x2=in_array[seg][3];
        joy_y2=in_array[seg][4];



        if(is_oscillating){
            if( t<= in_array[0][0]){
                joy_x1=0;
                joy_y1=1;
                joy_x2=cos(2*my_pi*joy_freq*t);;
                joy_y2=0;
            }
            else{
                joy_x1=cos(2*my_pi*joy_freq*t);
                joy_y1=sin(2*my_pi*joy_freq*t);
                joy_x2=0;
                joy_y2=0;
            }
            
        }
    }


}