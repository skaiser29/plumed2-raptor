/* -------------------------------------------------------------------------- */

#include "BSpline.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define TOL_ZERO 1e-8
#define BSPLINE_SCALE 0.99
#define BSPLINE_FDX 1e-4

/* -------------------------------------------------------------------------- */
namespace PLMD {

BSpline::BSpline(int arg_order, double arg_min, double arg_max, int arg_num_bins)
{
  do_write_debug = false;

  // Settings from fit.inp
  order     = arg_order;
  x_min     = arg_min;
  x_max     = arg_max;
  num_bins  = arg_num_bins;

  basis = basis_old = basis_deriv = NULL;
  knots = NULL;

  // Initialize
  num_knots = num_bins + 1 + order * 2;
  //num_knots = num_bins + 1 + order ;
  num_control_points = num_bins + order;

  dx = (x_max - x_min) / static_cast<double>(num_bins);

  knots = (double *) malloc(sizeof(double)*(num_knots+2));
  int pos = 0;
  for(int i=0; i<order; i++)    knots[pos++] = x_min;
  for(int i=0; i<num_bins; i++) knots[pos++] = x_min + static_cast<double>(i) * dx;
  //for(int i=0; i<=order; i++)   knots[pos++] = x_max;
  for(int i=0; i<=order+2; i++)   knots[pos++] = x_max;
  
  basis = (double *) malloc(sizeof(double)*num_knots);
  if(!basis) printf("BSpline: cannot malloc for basis");
  basis_old = (double *) malloc(sizeof(double)*num_knots);
  if(!basis_old) printf("BSpline: cannot malloc for basis_old");
  basis_deriv = (double *) malloc(sizeof(double)*num_knots);
  if(!basis_deriv) printf("BSpline: cannot malloc for basis_deriv");
}

BSpline::~BSpline()
{
  free(knots);
  free(basis);
  free(basis_old);
  free(basis_deriv);
}

double BSpline::Compute_Basis(double x)
{
  // fprintf(stdout,"Inside Compute_Basis with x= %f\n",x);

  int istart, istop;
  istart = 0;
  istop = num_knots;

  memset(basis, 0, sizeof(double)*num_knots);

  double num1, denom1, ratio1, term1, b01;
  double num2, denom2, ratio2, term2, b02;

  // Compute linear basis functions
  for(int i=istart; i<istop; i++) {
  //for(int i=istart; i<istop-2; i++) {

    // Compute first term
    num1 = x - knots[i];
    denom1 = knots[i+1] - knots[i];
    if(fabs(denom1) > TOL_ZERO) ratio1 = num1 / denom1;
    else ratio1 = 1.0;
    b01 = Compute_Basis0(x, knots[i], knots[i+1], x_max);
    term1 = ratio1 * b01;

    // Compute second term
    num2 = knots[i+2] - x;
    denom2 = knots[i+2] - knots[i+1];
    if(fabs(denom2) > TOL_ZERO) ratio2 = num2 / denom2;
    else ratio2 = 1.0;
    b02 = Compute_Basis0(x, knots[i+1], knots[i+2], x_max);
    term2 = ratio2 * b02;
    
    basis[i] = term1 + term2;
    
    // fprintf(stdout,"Computing basis function: i= %i  basis= %f\n",i,basis[i]);
  }


  double bsplinenorm=0.0;
  if(order == 1) return 0;

  // Evaluate higher order basis functions
  for(int i=2; i<=order; i++) {
    memcpy(basis_old, basis, sizeof(double)*num_knots);
    memset(basis, 0, sizeof(double)*num_knots);

    istop--;

    bsplinenorm = 0.0;
    for(int j=istart; j<istop; j++) {
      num1 = x - knots[j];
      denom1 = knots[j+i] - knots[j];
      if(fabs(denom1) > TOL_ZERO) ratio1 = num1 / denom1;
      else ratio1 = 0.0;
      term1 = ratio1 * basis_old[j];

      num2 = knots[j+i+1] - x;
      denom2 = knots[j+i+1] - knots[j+1];
      if(fabs(denom2) > TOL_ZERO) ratio2 = num2 / denom2;
      else ratio2 = 0.0;
      term2 = ratio2 * basis_old[j+1];

      basis[j] = term1 + term2;

      bsplinenorm+= basis[j];
//if(i==order and j==istop-1) printf(" 0:  x=%lf, bsplinenorm=%lf term1=%lf term2=%lf\n",x,bsplinenorm,term1,term2);


      // fprintf(stdout,"num1= %f  denom1= %f  ratio1= %f  term1= %f\n",num1,denom1,ratio1,term1);
      // fprintf(stdout,"num2= %f  denom2= %f  ratio2= %f  term2= %f\n",num2,denom2,ratio2,term2);
      // fprintf(stdout,"Computing basis functions of order= %i : i= %i  basis= %f  norm= %f\n",
      // 	      i,j,basis[j],norm);
    }
  }

//if(i==order) printf("basis[%d]=%lf bsplinenorm=%lf\n",j,basis[j],bsplinenorm);
//printf(" 1:  x=%lf, bsplinenorm=%lf\n",x,bsplinenorm);

  // Evaluate derivatives of basis functions
  memset(basis_deriv, 0, sizeof(double)*num_knots);
  for(int i=istart; i<istop; i++) {
    num1 = order;
    denom1 = knots[i+order] - knots[i];
    if(fabs(denom1) > TOL_ZERO) ratio1 = num1 / denom1;
    else ratio1 = 0.0;
    term1 = ratio1 * basis_old[i];
//printf(" 2:  x=%lf, bsplinenorm=%lf\n",x,bsplinenorm);

    num2 = order;
    denom2 = knots[i+order+1] - knots[i+1];
    if(fabs(denom2) > TOL_ZERO) ratio2 = num2 / denom2;
    else ratio2 = 0.0;
    term2 = ratio2 * basis_old[i+1];

    basis_deriv[i] = term1 - term2;

    // fprintf(stdout,"i= %i  ratio1= %f  term1= %f  num2= %f  denom2= %f  ratio2= %f  term2= %f  basis_deriv= %f\n",
    // 	    i,ratio1,term1,num2,denom2,ratio2,term2,basis_deriv[i]);
  }

//printf(" 3:  x=%lf, bsplinenorm=%lf\n\n",x,bsplinenorm);
  //if(fabs(1.0 - bsplinenorm) > TOL_ZERO) return 1;
  return bsplinenorm;
}

double BSpline::Compute_Basis0(double xx, double lo, double hi, double max)
{
  double x = xx;

  // If point is close to largest knot, then shift to slightly to right of largest knot
  if(fabs(x - max) < TOL_ZERO) x = max + TOL_ZERO * TOL_ZERO;

  if(do_write_debug) {
    if(fabs(hi-lo) < TOL_ZERO)  fprintf(stdout,"Test 1 (hi < lo ) is true.\n");
    if(fabs(x-hi) < TOL_ZERO)   fprintf(stdout,"Test 2 ( x - hi ) is true.\n");
    if(fabs(hi-max) < TOL_ZERO) fprintf(stdout,"Test 3 (hi - max) is true.\n");
    if((lo <= x) && (x < hi)) fprintf(stdout,"Test 4 (lo <= x < hi) is true.\n");
    fprintf(stdout,"xx= %f  x= %f  lo= %f  hi= %f  max= %f\n",xx,x,lo,hi,max);
  }

  // If hi and lo are the same knot
  if(fabs(hi-lo) < TOL_ZERO) 
    if((fabs(x-hi) < TOL_ZERO) && (hi-max < TOL_ZERO)) return 1.0;

  // if hi and lo are different knots
  if((lo <= x) && (x < hi)) return 1.0;

  return 0.0;
}
}
