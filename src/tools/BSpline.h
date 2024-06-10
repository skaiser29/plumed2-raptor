#ifndef BSPLINE_H
#define BSPLINE_H

namespace PLMD {

class BSpline
{
 public:

  BSpline(int, double, double, int);
  ~BSpline();

  double Compute_Basis0(double, double, double, double);
  double Compute_Basis(double);

  bool do_write_debug; // Output debugging info
  int order;  // Polynomial order
  int num_control_points; // Number of (potentially) non-zero basis functions
  int num_bins;
  double x_min, x_max, dx;

  int num_knots;
  double * knots;

  double * basis;
  double * basis_old;
  double * basis_deriv;
  double bsplinenorm; 
};
}
#endif
