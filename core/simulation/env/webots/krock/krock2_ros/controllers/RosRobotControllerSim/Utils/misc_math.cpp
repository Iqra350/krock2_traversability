#include "misc_math.hpp"





/* PT1 filter for scalars */
double
pt1(double u, double y, double T, double Ts)
{
    y=exp(-Ts/T)*y+(1-exp(-Ts/T))*u;
    return y;
}

/* PT1 filter for vectors */
MatrixXd
pt1_vec(MatrixXd u, MatrixXd y, double T, double Ts)
{
    y=exp(-Ts/T)*y+(1-exp(-Ts/T))*u;
    return y;
}

/* rotation matrix around x axis */
Matrix3d
rotx(double ang)
{
    Matrix3d mat;
    mat <<  1, 0, 0, 
            0, cos(ang), -sin(ang), 
            0, sin(ang), cos(ang);
    return mat;
}

/* rotation matrix around y axis */
Matrix3d
roty(double ang)
{
    Matrix3d mat;
    mat <<  cos(ang), 0, sin(ang),
            0, 1, 0,
            -sin(ang), 0, cos(ang);
    return mat;
}

/* rotation matrix around z axis */
Matrix3d
rotz(double ang)
{
    Matrix3d mat;
    mat <<  cos(ang), -sin(ang), 0,
            sin(ang), cos(ang), 0,
            0,0,1;

    return mat;
}



/* Hermite spline */
MatrixXd
hermiteSpline(Vector2d p0, Vector2d p1, Vector2d m0, Vector2d m1, MatrixXd t)
{
    MatrixXd p(2,t.size());
    for(int i=0; i<2; i++){
        p.block(i,0,1,t.size())=( 2*t.array().pow(3) - 3*t.array().pow(2) + 1)*p0(i) + 
                                (   t.array().pow(3) - 2*t.array().pow(2) + t.array())*m0(i) + 
                                (-2*t.array().pow(3) + 3*t.array().pow(2))*p1(i) + 
                                (   t.array().pow(3) -   t.array().pow(2))*m1(i);
    }
    return p;
}



/* Save acos */
double SafeAcos (double x)
{
	if (x < -1.0) x = -1.0 ;
	else if (x > 1.0) x = 1.0 ;
	return acos (x) ;
}




/* --------------------------- PID class --------------------------------*/
PID :: PID(double Kpin, double Kiin, double Kdin, double dt){
  Kp=Kpin;
  Ki=Kiin;
  Kd=Kdin;
  Ts=dt;

  ek1=0;
  ek2=0;
  uk1=0;

  a=Kp+Ki*Ts/2.+Kd/Ts;
  b=-Kp+Ki*Ts/2.-2*Kd/Ts;
  c=Kd/Ts;
}

double
PID :: calc(double e)
{
  double u;
  u=uk1 + a*e + b*ek1 + c*ek2;

  uk1=u;
  ek2=ek1;
  ek1=e;

  return u;
}

PID& PID::operator=( const PID& other ) 
{
      Kp = other.Kp;
      Ki = other.Ki;
      Kd = other.Kd;
      Ts = other.Ts;
      a = other.a;
      b = other.b;
      c = other.c;
      ek1 = other.ek1;
      ek2 = other.ek2;
      uk1 = other.uk1;
      return *this;
}





/* --------------------------- PIDvec class --------------------------------*/


void
PIDvec :: init(MatrixXd Kpin, MatrixXd Kiin, MatrixXd Kdin, MatrixXd Tdin, double lim, double dt){
  Kp=Kpin;
  Ki=Kiin;
  Kd=Kdin;
  Td=Tdin;
  Ts=dt;

  ek1.resize(Kp.rows(), Kp.cols());
  ek2.resize(Kp.rows(), Kp.cols());
  uk1.resize(Kp.rows(), Kp.cols());
  uk2.resize(Kp.rows(), Kp.cols());

  ek1.setZero();
  ek2.setZero();
  uk1.setZero();
  uk2.setZero();

  //a=Kp+Ki*Ts/2.+Kd/Ts;
  //b=-Kp+Ki*Ts/2.-2*Kd/Ts;
  //c=Kd/Ts;
  a=Kd+Kp.cwiseProduct(Td);
  b=-2*Kd - 2*Kp.cwiseProduct(Td) + Kp*Ts + Ki.cwiseProduct(Td)*Ts;
  c=Kd + Kp.cwiseProduct(Td) - Kp*Ts - Ki.cwiseProduct(Td)*Ts + Ki*Ts*Ts;


  d=Td; 
  f=2*Td - Ts*MatrixXd::Ones(Kp.rows(), Kp.cols());
  g=Ts*MatrixXd::Ones(Kp.rows(), Kp.cols()) - Td;

  a=a.cwiseQuotient(d);
  b=b.cwiseQuotient(d);
  c=c.cwiseQuotient(d);
  f=f.cwiseQuotient(d);
  g=g.cwiseQuotient(d);

  limit=lim;
}

void
PIDvec :: init(double Kpin, double Kiin, double Kdin, double Tdin, double lim, double dt){
  Kp.resize(1,1);
  Ki.resize(1,1);
  Kd.resize(1,1);
  Td.resize(1,1);

  Kp(0)=Kpin;
  Ki(0)=Kiin;
  Kd(0)=Kdin;
  Td(0)=Tdin;
  Ts=dt;
  

  ek1.resize(1,1);
  ek2.resize(1,1);
  uk1.resize(1,1);
  uk2.resize(1,1);

  ek1.setZero();
  ek2.setZero();
  uk1.setZero();
  uk2.setZero();

  //a=Kp+Ki*Ts/2.+Kd/Ts;
  //b=-Kp+Ki*Ts/2.-2*Kd/Ts;
  //c=Kd/Ts;
  a=Kd+Kp.cwiseProduct(Td);
  b=-2*Kd - 2*Kp.cwiseProduct(Td) + Kp*Ts + Ki.cwiseProduct(Td)*Ts;
  c=Kd + Kp.cwiseProduct(Td) - Kp*Ts - Ki.cwiseProduct(Td)*Ts + Ki*Ts*Ts;


  d=Td; 
  f=2*Td - Ts*MatrixXd::Ones(Kp.rows(), Kp.cols());
  g=Ts*MatrixXd::Ones(Kp.rows(), Kp.cols()) - Td;

  a=a.cwiseQuotient(d);
  b=b.cwiseQuotient(d);
  c=c.cwiseQuotient(d);
  f=f.cwiseQuotient(d);
  g=g.cwiseQuotient(d);

  limit=lim;
}

MatrixXd
PIDvec :: calc(MatrixXd e)
{
  MatrixXd u = uk1.cwiseProduct(f) + uk2.cwiseProduct(g) + e.cwiseProduct(a) + ek1.cwiseProduct(b) + ek2.cwiseProduct(c);

  for(int i=0; i<u.size(); i++){
    u(i) = u(i) > limit ? limit : u(i);
    u(i) = u(i) <-limit ?-limit : u(i);
  }

  uk2=uk1;
  uk1=u;
  ek2=ek1;
  ek1=e;

  return u;
}

double
PIDvec :: calc(double e)
{
  MatrixXd u = uk1.cwiseProduct(f) + uk2.cwiseProduct(g) + e*a + ek1.cwiseProduct(b) + ek2.cwiseProduct(c);

  for(int i=0; i<u.size(); i++){
    u(i) = u(i) > limit ? limit : u(i);
    u(i) = u(i) <-limit ?-limit : u(i);
  }




  uk2=uk1;
  uk1=u;
  ek2=ek1;
  ek1(0)=e;
  
  return u(0);
}

PIDvec& PIDvec::operator=( const PIDvec& other ) 
{
      Kp = other.Kp;
      Ki = other.Ki;
      Kd = other.Kd;
      Ts = other.Ts;
      a = other.a;
      b = other.b;
      c = other.c;
      ek1 = other.ek1;
      ek2 = other.ek2;
      uk1 = other.uk1;
      return *this;
}

void 
PIDvec::flush()
{
  uk1.setZero();
  uk2.setZero();
  ek2.setZero();
  ek1.setZero();
}



Vector3d 
projectOntoEllipsoid(Vector3d p, Vector3d e, Vector3d c)
{
  Vector3d proj;
  // center the point
  p=p-c;

  // get angles of the point
  double phi = atan2(p(1)*e(0),p(0)*e(1));
  double theta = atan2(p(2)/e(2), sqrt(pow(p(0)/e(0), 2) + pow(p(1)/e(1), 2)));

  proj(0) = e(0)*cos(theta)*cos(phi) + c(0);
  proj(1) = e(1)*cos(theta)*sin(phi) + c(1);
  proj(2) = e(2)*sin(theta) + c(2);

  return proj;
}


/* --------------------------- movAVG class --------------------------------*/
movAVG :: movAVG(int dim, int span)
{
  N=span;
  n=dim;
  old.resize(dim, N);
  old.setZero();
}

MatrixXd
movAVG :: calc(MatrixXd in)
{
  MatrixXd tmp(n, N-1);

  tmp=old.block(0,0,n,N-1);
  old.block(0,1,n,N-1)=tmp;

  old.block(0,0,n,1)=in;

/*
  MatrixXd out(n,1);
  for(int k=0; k<n; k++){
      out(k,0)=old.block(k,0,1,N).sum()/N;
  }*/

  return old.rowwise().sum()/N;
  //return out;
}