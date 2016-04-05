#include <Fepic/CustomEigen>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

// space
typedef Matrix<double, Dynamic,1,0,3,1>              Vector;
typedef Matrix<double, Dynamic,Dynamic,RowMajor,3,3> Tensor;
const double pi  = 3.141592653589793;
const double pi2 = pi*pi;

double pho(Vector const& X, int tag);
double gama(Vector const& X, double t, int tag);
double muu(int tag);
Vector force(Vector const& X, double t, int tag);   //density*gravity (force/vol)
Vector u_exact(Vector const& X, double t, int tag);
double p_exact(Vector const& X, double t, int tag);
Vector z_exact(Vector const& X, double t, int tag);
Tensor grad_u_exact(Vector const& X, double t, int tag);
Vector grad_p_exact(Vector const& X, double t, int tag);
Vector traction(Vector const& X, Vector const& normal, double t, int tag);
Vector u_initial(Vector const& X, int tag);
double p_initial(Vector const& X, int tag);
Vector z_initial(Vector const& X, int tag);
Vector solid_normal(Vector const& X, double t, int tag);
Tensor feature_proj(Vector const& X, double t, int tag);
Vector gravity(Vector const& X);


// gota estática 2d/////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
  if (tag == 16)
  {
    return 1.0e3;
  }
  else
  {
    return 8.0e2;
  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.02;
}

double muu(int tag)
{
  if (tag == 16)
  {
    return 1.0e-3;
  }
  else
  {
    return 3.0e-1;
  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
  f(1) = 10;
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double un = 998*9.8*2.85e-4*2.85e-4/(3*1e-3);
  Vector v(Vector::Zero(X.size()));

  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  return dxU;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  if (tag == 16)
  {
    return 400.0;
  }
  else
  {
    return 0.0;
  }
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const tet = 10. * (pi/180.);
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

#endif


// solid asc 2d/////////////////////////////////////////////////////////////
#if (true)

double pho(Vector const& X, int tag)
{
//  if (tag == 15)
//  {
    return 1.0;///1e4;
//  }
//  else
//  {
//    return 0.0;
//  }
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0.0;
}

double beta_diss()
{
  return 0.0;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.5;
}

double muu(int tag)
{
//  if (tag == 15)
//  {
    return 0.1;
//  }
//  else
//  {
//    return 0.0;
//  }
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
//  if (tag == 15)
//  {
    f(1) = 10.0*1.0;//*1e4;//*1e3;
//
//  else
//  {
//    f(1) = 0.0;  //-8e-4*1e4;
//  }
  return f;
}

Vector gravity(Vector const& X){
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
  //if (tag == 15)
  //{
    f(1) = 10.0;//-8e-4;  //*1e3;
  //}
  //else
  //{
  //  f(1) = 0.0;  //-8e-4*1e4;
  //}
  return f;
}


Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double un = 998*9.8*2.85e-4*2.85e-4/(3*1e-3);
  //Vector v(Vector::Ones(X.size()));
  Vector v(Vector::Zero(X.size()));

  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  return dxU;
}

Vector z_exact(Vector const& X, double t, int tag)
{
  Vector v(Vector::Zero(3)); //v << 0, 0, -1000;
  //Vector v(Vector::Ones(3));
  return v;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  
  return 0.0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

Vector z_initial(Vector const& X, int tag)
{
  return z_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (tag == 2) {N(1) = 1;};
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

#endif


// canal /////////////////////////////////////////////////////////////
#if (false)

double pho(Vector const& X, int tag)
{
  return 1.0;
}

double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 0;
}

double muu(int tag)
{
  return  1.0;
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
  return f;
}

Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double umax = 1;
  Vector v(Vector::Zero(X.size()));

  if (tag == 3 || tag == 2)
  //if (tag == 2)
  {
    v(0) = umax*y*(2-y);
    v(1) = 0;
  }

  return v;
}

Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  return dxU;
}

double p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0;
}

Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(X.size());

  return dxP;
}

Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));
  //T(0) = -p_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));
  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));
  T = (- p_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;
  return T;
}

Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}

double p_initial(Vector const& X, int tag)
{
  return p_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const tet = 10. * (pi/180.);
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //N(1) = 1;
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  if (tag == 3)
  {
    f(0,0) = 1;
  }
//  else if (tag == 1)
//  {
//    f(0,0) = 1;
//  }
  return f;
}

#endif
