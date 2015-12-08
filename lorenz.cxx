#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

void func (const double& x, const double& y, const double& z, double* f) {
  double a = 10, b = 28, c = 8/3;
  f[0] = a*(y-x);
  f[1] = x*(b-z)-y;
  f[2] = x*y-c*z;
}

int main() {
  double r[3];
  r[0] = 1, r[1] = 1, r[2] = 1;		// starting values
  
  double k1[3], k2[3], k3[3], k4[3];	// coefficients
  
  double const t_start = 0; 		// starting time
  double const t_end = 100; 		// end time
  double const dt = 0.001; 		// step size
  int const N = (t_end - t_start)/dt;
  
  ofstream out("d0001");
  out << t_start << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << endl;
  for (int i = 0; i<N; i++) {
    func (r[0], r[1], r[2], k1);
    func (r[0] + 0.5*k1[0]*dt, r[1] + 0.5*k1[1]*dt, r[2] + 0.5*k1[2]*dt, k2);
    func (r[0] + 0.5*k2[0]*dt, r[1] + 0.5*k2[1]*dt, r[2] + 0.5*k2[2]*dt, k3);
    func (r[0] + k3[0]*dt, r[1] + k3[1]*dt, r[2] + k3[2]*dt, k4);
    
    r[0] = r[0] + (1.0/6*k1[0] + 1.0/3*k2[0] + 1.0/3*k3[0] + 1.0/6*k4[0])*dt;
    r[1] = r[1] + (1.0/6*k1[1] + 1.0/3*k2[1] + 1.0/3*k3[1] + 1.0/6*k4[1])*dt;
    r[2] = r[2] + (1.0/6*k1[2] + 1.0/3*k2[2] + 1.0/3*k3[2] + 1.0/6*k4[2])*dt;
    out << i*dt << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << endl;
  }
  out.close();
  
  
  return 0;
}