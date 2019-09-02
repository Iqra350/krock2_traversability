#include "controller.hpp"
#include <iostream>
#include <deque>

#define CAM_WIDTH 320
#define CAM_HEIGHT 240


#define IMG_W (160*2)
#define IMG_H (120*2)
#define CAM_PYRAMID 1


#define INTEGRATION_TIME 2

using namespace cv;
using namespace Eigen;
using namespace std;

pthread_t thread_imgProc;
pthread_mutex_t imgIsReadyMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  imgIsReadyCond   = PTHREAD_COND_INITIALIZER;

double movavg(double u, double *vec, int N){

    double sum=0;
    for(int i=0;i<N-1;i++){
        vec[i]=vec[i+1];
        sum+=vec[i];
    }
    vec[N-1]=u;
    sum+=u;
    sum=sum/N;
    return sum;



}

void * ControllerImageProcThread(void *ctrl){
    cout<<"it's in"<<endl;
    Controller *ctrl2=(Controller*)(ctrl);
    while(1){
        pthread_cond_wait( &imgIsReadyCond, &imgIsReadyMutex );
        ctrl2->optFlowLK();
        ctrl2->getCamera();
     //   ctrl2->visualFeedback();
        pthread_mutex_unlock( &imgIsReadyMutex );
    }
}

void
Controller :: getCamera(){
    static int cnt=0;
  /*  static VideoWriter writer("testVideo.avi", 
               CV_FOURCC('D','I','V','X'),
               25,
               Size(CAM_WIDTH, CAM_HEIGHT));*/
    
    
    if(!cnt++){
     //   namedWindow("MyWindow", CV_WINDOW_AUTOSIZE);
        cout<<"window created"<<endl;
    }
    //cout<<"I'm here"<<endl;
    //cout<<*imgArray<<endl;
    camImg=Mat(CAM_HEIGHT, CAM_WIDTH, CV_8UC3, imgArray);
   // imshow("MyWindow", camImg);
    //writer.write(camImg);
    //cout<<"WRITING"<<endl<<endl;
    //cout<<*imgPtr<<endl;
    waitKey(1);
    //cout<<"I'm here22"<<endl;
}

void
Controller :: initImageProcessing(unsigned char *srcPtr){

    imgArray=srcPtr;
    pthread_create( &thread_imgProc, NULL, ControllerImageProcThread, this);

}

void
Controller :: optFlowFarneback(){
    static Mat src, src_old, flow, xy[2], channel[3], magnitude, angle, frame2;
    static int cnt=0;
    cout<<"FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW START"<<endl<<endl;

    src=Mat(CAM_HEIGHT*CAM_PYRAMID, CAM_WIDTH*CAM_PYRAMID, CV_8UC3, imgArray);
    cvtColor( src, src, CV_RGB2GRAY );
    if(CAM_PYRAMID>1){
        pyrDown(src, src, Size(CAM_WIDTH, CAM_HEIGHT));
    }
    if(!cnt++){
        src_old=src;
      //  namedWindow("flow", CV_WINDOW_AUTOSIZE);
        cout<<"window created"<<endl;
    }


    Sobel(src, src, 3, 1, 1, 1, 1, 0, BORDER_DEFAULT );
    
    calcOpticalFlowFarneback(src_old, src, flow, 0.5, 1, 1, 1, 1, 1., 0);
    
    //calcOpticalFlowFarneback(src_old, src, flow, 0.6, 4, 5, 5, 3, 5., 0);
    split(flow, xy);
    cartToPolar(xy[0], xy[1], magnitude, angle, true);
    angle=angle/2;
    magnitude=magnitude*255;
    magnitude.convertTo(magnitude, CV_8UC1);
    angle.convertTo(angle, CV_8UC1);

    channel[0]=angle;
    channel[2]=magnitude;
    channel[1]=Mat::ones(CAM_HEIGHT, CAM_WIDTH, CV_8UC1)*255;

    angle.convertTo(angle, CV_32F);

    cv2eigen(xy[0],imX);
    cv2eigen(xy[1],imY);








    //cout<<b<<endl;


    merge(channel, 3, frame2);
    cvtColor(frame2, frame2, CV_HSV2RGB, 0);
   // imshow("flow", frame2);
    src_old=src;
    waitKey(1);

    cout<<"FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW END"<<endl<<endl;
}

void
Controller :: visualFeedback(){
    static ofstream oflow("oflow.txt");
    static int cnt=0;
    double leftFlowX, leftFlowY, rightFlowX, rightFlowY, Q1X, Q2X, Q3X, Q4X, Q1Y, Q2Y, Q3Y, Q4Y;
    static double leftFlowX_old=0, leftFlowY_old=0, rightFlowX_old=0, rightFlowY_old=0, Q1X_old=0, Q2X_old=0, Q3X_old=0, Q4X_old=0, Q1Y_old=0, Q2Y_old=0, Q3Y_old=0, Q4Y_old=0, gamma_old=0;
    const int width_off = 100, height_off=0;
    static double w_yaw=0;
    const int N=15;
    static double Q1Xvec[N]={0}, Q2Xvec[N]={0}, Q3Xvec[N]={0}, Q4Xvec[N]={0}, Q1Yvec[N]={0}, Q2Yvec[N]={0}, Q3Yvec[N]={0}, Q4Yvec[N]={0};

    leftFlowX=imX.block<CAM_HEIGHT, CAM_WIDTH/2>(0,0).sum();
    leftFlowY=imY.block<CAM_HEIGHT, CAM_WIDTH/2>(0,0).sum();
    rightFlowX=imX.block<CAM_HEIGHT, CAM_WIDTH/2>(0,CAM_WIDTH/2).sum();
    rightFlowY=imY.block<CAM_HEIGHT, CAM_WIDTH/2>(0,CAM_WIDTH/2).sum();

    Q1X=imX.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(0,CAM_WIDTH/2+width_off).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;
    Q2X=imX.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(0,0).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;
    Q3X=imX.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(CAM_HEIGHT/2+height_off,0).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;
    Q4X=imX.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(CAM_HEIGHT/2+height_off,CAM_WIDTH/2+width_off).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;
    Q1Y=imY.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(0,CAM_WIDTH/2+width_off).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;
    Q2Y=imY.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(0,0).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;
    Q3Y=imY.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(CAM_HEIGHT/2+height_off,0).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;
    Q4Y=imY.block<CAM_HEIGHT/2-height_off, CAM_WIDTH/2-width_off>(CAM_HEIGHT/2+height_off,CAM_WIDTH/2+width_off).sum()/(double)((CAM_HEIGHT/2.-height_off)*(CAM_WIDTH/2.-width_off));;

    // filtering
    double img_filt=3;
    double gamma_filt=0.1;
    double gyro_filt=0.01;

    leftFlowX = pt1(leftFlowX, leftFlowX_old, img_filt, dt);
    leftFlowY = pt1(leftFlowY, leftFlowY_old, img_filt, dt);
    rightFlowX = pt1(rightFlowX, rightFlowX_old, img_filt, dt);
    rightFlowY = pt1(rightFlowY, rightFlowY_old, img_filt, dt);
/*
    Q1X = pt1(Q1X, Q1X_old, img_filt, dt);
    Q2X = pt1(Q2X, Q2X_old, img_filt, dt);
    Q3X = pt1(Q3X, Q3X_old, img_filt, dt);
    Q4X = pt1(Q4X, Q4X_old, img_filt, dt);
    Q1Y = pt1(Q1Y, Q1Y_old, img_filt, dt);
    Q2Y = pt1(Q2Y, Q2Y_old, img_filt, dt);
    Q3Y = pt1(Q3Y, Q3Y_old, img_filt, dt);
    Q4Y = pt1(Q4Y, Q4Y_old, img_filt, dt);
*/


  Q1X = movavg(Q1X, Q1Xvec, N);
    Q2X = movavg(Q2X, Q2Xvec, N);
    Q3X = movavg(Q3X, Q3Xvec, N);
    Q4X = movavg(Q4X, Q4Xvec, N);

    Q1Y = movavg(Q1Y, Q1Yvec, N);
    Q2Y = movavg(Q2Y, Q2Yvec, N);
    Q3Y = movavg(Q3Y, Q3Yvec, N);
    Q4Y = movavg(Q4Y, Q4Yvec, N);



    w_yaw=pt1(gyroData[1], w_yaw, gyro_filt, dt);

    //cout<<gamma<<"   "<<gamma_old<<endl;

    //gamma = -0.03 * (Q1X + Q2X + Q3X + Q4X) ;
    gamma = -3 * (Q1X + Q2X)  -  2*(Q1Y  -  Q2Y) + w_yaw*0.0;
    gamma=gamma*1.0;


    //cout<<gamma<<endl;

    gamma=pt1(gamma, gamma_old, gamma_filt, dt);

    //gamma=my_pi/2*sin(0.5*t);
    //cout<<leftFlowX<<"\t"<<leftFlowY<<"\t"<<rightFlowX<<"\t"<<rightFlowY<<endl;

    oflow<<Q1X<<"\t"<<Q2X<<"\t"<<Q3X<<"\t"<<Q4X<<"\t"<<Q1Y<<"\t"<<Q2Y<<"\t"<<Q3Y<<"\t"<<Q4Y<<"\t"<<gamma<<endl;
    // update old values
    Q1X_old=Q1X;
    Q2X_old=Q2X;
    Q3X_old=Q3X;
    Q4X_old=Q4X;
    Q1Y_old=Q1Y;
    Q2Y_old=Q2Y;
    Q3Y_old=Q3Y;
    Q4Y_old=Q4Y;
    gamma_old=gamma;

    //cout<<tcyc<<endl;
    tcyc=1.0;




}

void
Controller :: optFlowLK(){
    
    static int MA_WINDOW_SIZE = (int)round(MA_WINDOW_SIZE_TIME/(frame_time_step/1000.));
    static deque<double> hist_contX(MA_WINDOW_SIZE, 0);
    static deque<double> hist_contY(MA_WINDOW_SIZE, 0);
    
    static char rawWindow[] = "Raw Video";
    static char opticalFlowWindow[] = "Optical Flow Window";

    static int cnt=0;
    
    static Mat frame, frame2, grayFrames, grayFramesOrig, rgbFrames, prevGrayFrame;
    static Mat opticalFlow = Mat(IMG_H,  IMG_W, CV_32FC3);
    
    static vector<Point2f> points1;
    static vector<Point2f> points2;

    static Point2f diff;

    static vector<uchar> status;
    static vector<float> err;

    static RNG rng(12345);
    static Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255),rng.uniform(0, 255));
				    
    static int i, k;
    static TermCriteria termcrit(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, LK_flow_term_crit_max_count, LK_flow_term_crit_eps);
    static Size winSize(LK_flow_win_size[0], LK_flow_win_size[1]);			//subPixWinSize(10, 10), 
    static bool needToInit=true;
    
    static double angle;
    static double t=0, dt=0, t_reset=0;
    static double t0, t1, t2;
    static int wait_t;
    
    // Lambda Matrix
    static Mat lambda( 2, 4, CV_32FC1 );
    // Input Quadilateral or Image plane coordinates
    static Point2f inputQuad[4]; 
    // Output Quadilateral or World plane coordinates
    static Point2f outputQuad[4];
    static ofstream flowVectors("flowVectors.txt");
    static ofstream filtflow("filtflow.txt");
    static double meanFlowX, meanFlowY;
    
    
    // ESTIMATION

    
    
  //  cout<<"FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW START"<<endl<<endl;


    if(!cnt++){
	namedWindow(rawWindow, CV_WINDOW_AUTOSIZE);
	//namedWindow(opticalFlowWindow, CV_WINDOW_AUTOSIZE);
	/////////////////////////////////// PERSPECTIVE TRANSFORMATION ////////////////////////////////////////////////////////
	inputQuad[0] = Point2f( 0-warp_in_width,warp_in_height/CAM_PYRAMID );
	inputQuad[1] = Point2f( IMG_W-1+warp_in_width,warp_in_height/CAM_PYRAMID);
	inputQuad[2] = Point2f( IMG_W-1,IMG_H-1);
	inputQuad[3] = Point2f( 0,IMG_H-1  );  
	
	outputQuad[0] = Point2f( 0,0 );
	outputQuad[1] = Point2f( IMG_W-1,0);
	outputQuad[2] = Point2f( IMG_W-(IMG_W/2-crop_W/2),IMG_H-1);
	outputQuad[3] = Point2f( IMG_W/2-crop_W/2,IMG_H-1  );

	// Get the Perspective Transform Matrix i.e. lambda 
	lambda = getPerspectiveTransform( inputQuad, outputQuad );
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
        
        cout<<"window created"<<endl;
    }

    t+=dt;
    t_reset+=dt;
    t0=get_timestamp();
    
    ////////////////////////////   PREPROCESSING    //////////////////////////////////////////////////////////////
    frame=Mat(CAM_HEIGHT, CAM_WIDTH, CV_8UC3, imgArray);
    
    if(CAM_PYRAMID>1){
        pyrDown(frame, frame, Size(IMG_W, IMG_H));
    }
    
    warpPerspective(frame,frame2,lambda,frame.size() );
    frame=frame2;
    frame.copyTo(rgbFrames);
    cvtColor(rgbFrames, grayFramesOrig, CV_BGR2GRAY);
    //grayFrames=grayFramesOrig(Rect(0, 150/CAM_PYRAMID, IMG_W, IMG_H-150/CAM_PYRAMID-1)).clone();
    grayFrames=grayFramesOrig(Rect(IMG_W/2-crop_W/2+1, crop_H0, crop_W-2, IMG_H-crop_H0-2)).clone();
    //grayFrames=grayFramesOrig;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	meanFlowX=0;
	meanFlowY=0;
    if (needToInit || (cnt<3)) {
	goodFeaturesToTrack(grayFrames, points1, track_feat_max_count, track_feat_quality, track_feat_min_dist, Mat(), 
			    track_feat_block_size, track_feat_use_Harris, track_feat_Harris_k);
	needToInit = false;

    } 
    else if (!points2.empty()) {
	
	calcOpticalFlowPyrLK(prevGrayFrame, grayFrames, points2, points1,
	    status, err, winSize, LK_flow_max_level, termcrit, LK_flow_flags, LK_flow_min_eig_treshhold);
	
	for (i = k = 0; i < points2.size(); i++) {
	  //  cout<<"doin something"<<endl;
	 //   if(points1[i].y<130)
	//	continue;
	/*    cout << "Optical Flow Difference... X is "
		<< int(points1[i].x - points2[i].x) << "\t Y is "
		<< int(points1[i].y - points2[i].y) << "\t\t" << i
		<< "\n";*/
	    flowVectors << points1[i].x - points2[i].x << "\t" << points1[i].y - points2[i].y<< "\t" ;
	    meanFlowX+=points1[i].x - points2[i].x;
	    meanFlowY+=points1[i].y - points2[i].y;
	    //cout<<points1[i].y<<endl;

	    if ((points1[i].x - points2[i].x) > 0) {
		line(grayFrames, points1[i], points2[i], Scalar(0, 0, 255),
		1, 1, 0);

		circle(grayFrames, points1[i], 2, Scalar(255, 0, 0), 1, 1,
		0);

		line(opticalFlow, points1[i], points2[i], Scalar(0, 0, 255),
		1, 1, 0);
		circle(opticalFlow, points1[i], 1, Scalar(255, 0, 0), 1, 1,
		0);
	    } 
	    else {
		line(grayFrames, points1[i], points2[i], Scalar(0, 255, 0),
		1, 1, 0);

		circle(grayFrames, points1[i], 2, Scalar(255, 0, 0), 1, 1,
		0);

		line(opticalFlow, points1[i], points2[i], Scalar(0, 255, 0),
		1, 1, 0);
		circle(opticalFlow, points1[i], 1, Scalar(255, 0, 0), 1, 1,
		0);
	    }
	    points1[k++] = points1[i];

	}
	flowVectors << endl;
	meanFlowX/=(double)points2.size();
	meanFlowY/=(double)points2.size();
	
	goodFeaturesToTrack(grayFrames, points1, track_feat_max_count, track_feat_quality, track_feat_min_dist, Mat(), 
			    track_feat_block_size, track_feat_use_Harris, track_feat_Harris_k);

	
	
	
	filtFlowX+=meanFlowX/(double)MA_WINDOW_SIZE;
	filtFlowY+=meanFlowY/(double)MA_WINDOW_SIZE;
	
	filtFlowX-=hist_contX.back();
	filtFlowY-=hist_contY.back();

	hist_contX.pop_back();
	hist_contY.pop_back();
	hist_contX.push_front(meanFlowX/(double)MA_WINDOW_SIZE);
	hist_contY.push_front(meanFlowY/(double)MA_WINDOW_SIZE);
	filtflow<<filtFlowX<<"\t"<<filtFlowY<<"\t"<<meanFlowX<<"\t"<<meanFlowY<<endl;
	
	if(meanFlowY==0){
	    cout<<meanFlowY<<"\t"<<points2.size()<<endl;
	}
	
	
    }
    swap(points2, points1);
    points1.clear();
    grayFrames.copyTo(prevGrayFrame);

    
    imshow(rawWindow, grayFrames);
    
  //  imshow(opticalFlowWindow, opticalFlow);
    
    /////////////////////////////////////////// ESTIMATE VELOCITY /////////////////////////////////////////////////////

    //static double filtFlowX=0, filtFlowY=0;

    
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    t1=get_timestamp();     /// END TIME
    
    
    if(frame_time_step/1000.>(t1-t0)){
            usleep((frame_time_step/1000.- (t1-t0))*1000000 );
    }
        
        
    wait_t=frame_time_step-(int)((t1-t0)*1000);
    wait_t=wait_t>1?wait_t:1;
    //if(waitKey(wait_t) >= 0) break;


    t2=get_timestamp();     /// MEASURE TIME
  //  cout<<(t1-t0)*1000<<"\t"<<(t2-t0)*1000<<endl;
    dt=(t2-t0);

  //  cout<<"FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW END"<<endl<<endl;
}

