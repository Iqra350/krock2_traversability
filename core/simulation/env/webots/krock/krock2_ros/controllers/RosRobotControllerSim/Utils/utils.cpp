#include "utils.hpp"

/* Read current time */
double
get_real_time()
{
    struct timespec real_time;

    if (clock_gettime(CLOCK_REALTIME, &real_time) == -1 )
    {
        perror("clock gettime");
        exit( EXIT_FAILURE );
    }

    return((double)real_time.tv_sec + (double)(real_time.tv_nsec/BILLION));
}

/* Read current time - this one is being used*/
double
get_timestamp()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  (now.tv_usec + (unsigned long long)now.tv_sec * 1000000)/1000000.;
}



/* Sends an UDP package "128.178.148.59" */
int
sendUDP(void *data, int len, const char *IP, int UDP_PORT)
{

    static int sock = socket(AF_INET, SOCK_DGRAM, 0);

    struct sockaddr_in saddr = {0,};

    saddr.sin_family = AF_INET;
    saddr.sin_addr.s_addr = inet_addr(IP);
    saddr.sin_port = htons(UDP_PORT);


    //double t=0;




    //  double data[27 * 4];
    //  for(int i=0;i<27*4;i++){
    //      data[i]=t;

    //  }
    //  t+=0.001;
       // sendto(sock, (char *)data, len*sizeof(double), 0, (struct sockaddr *)&saddr, sizeof(saddr));
        sendto(sock, (char *)data, len, 0, (struct sockaddr *)&saddr, sizeof(saddr));
     //       double data2[3]={rand(),2,3};
     //      printf("%d\n", sendto(sock, (char *)data2, sizeof(data2), 0, (struct sockaddr *)&saddr, sizeof(saddr)));
       //     usleep(8e3);


 // create socket
/*
    const int port = 8472;
    Socket socket;
    if ( !socket.Open( port ) )
    {
        printf( "failed to create socket!\n" );
        return false;
    }efined reference to `JOYSTICK_TYPE'


    // send a packet

    const char data[] = "hello world!";
    socket.Send( Address(128,178,148,59,port), data, sizeof( data ) );
*/

    return 0;
}

/* Reads file, skipps commented lines starting with "//" and stores the data into an array */
int
readFileWithLineSkipping(ifstream& inputfile, stringstream& file){
    string line;
    file.str(std::string());
    int linenum=0;
    while (!inputfile.eof()){
        getline(inputfile,line);

        //line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));
        if (line.length() == 0 || (line[0] == '#') || (line[0] == ';')){
            //continue;
        }
        else
            file << line << "\n";
            linenum++;
        }
    return linenum-1;
}



/* sleep */
//void 
//sleepBoost(float sTime){
//    boost::this_thread::sleep(boost::posix_time::microseconds(sTime*1000000));
//}



