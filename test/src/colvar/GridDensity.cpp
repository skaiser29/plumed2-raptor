/*
 * Author: Chenghan Li
 *
 * As I allow use switching function whose dmax is larger than `binsize` here,
 * number of grid points associated with one atom is no longer fixed 8 as that
 * in GridNN. Thus, `evaluateDensity` used in GridNN needs generalization to 
 * handle dynamical number of mapping points.
 */
#include "Colvar.h"
#include "tools/Vector.h"
#include "tools/SwitchingFunction.h"
#include "ActionRegister.h"


namespace PLMD {
namespace colvar {

class GridDensity: public Colvar {
   bool do_pbc, do_mpi;
   SwitchingFunction sf;
   unsigned natoms;
   Vector rx,ry,rz;
   double xl,xh,yl,yh,zl,zh, dmax;
   double binsize;
   int nx, ny, nz;
   unsigned ngrids;

public:
   explicit GridDensity(const ActionOptions&);
   void calculate();
   static void registerKeywords(Keywords& keys);

};

PLUMED_REGISTER_ACTION(GridDensity, "GRID_DENSITY")

void GridDensity::registerKeywords(Keywords& keys) {
   Colvar::registerKeywords(keys);
   keys.add("atoms", "CENTER", "Center of the grid");
   keys.add("atoms", "ATOMS", "Atoms to be mapped onto grids");
   keys.add("optional", "ROTX", "Rotation vector in x axis");
   keys.add("optional", "ROTY", "Rotation vector in y axis");
   keys.add("optional", "ROTZ", "Rotation vector in z axis");
   keys.add("compulsory", "RANGES", 
            "Ranges of the grid in the order of lx hx ly hy lz hz");
   keys.add("compulsory", "BINSIZE", "Binsize of the grid");
   keys.add("compulsory", "SWITCH", 
            "Switching function to map atoms onto grids");
   keys.addFlag("SERIAL", false, "Do not use mpi");

   keys.addOutputComponent("gd", "default", 
         "result densities flattened in row-major order");
}

GridDensity::GridDensity(const ActionOptions&ao) :
PLUMED_COLVAR_INIT(ao),
do_pbc(true), do_mpi(false)
{
   bool nopbc = !do_pbc, nompi = !do_mpi;
   parseFlag("NOPBC", nopbc); do_pbc = !nopbc;
   parseFlag("SERIAL", nompi); do_mpi = !nompi;
   if(do_pbc) log.printf("  PBC will be considered\n");
   else       log.printf("  PBC will not be considered\n");
   if(do_mpi) log.printf("  MPI will be used\n");
   else       log.printf("  MPI will not be used\n");

   /*--- parse ATOMS ---*/
   std::vector<AtomNumber> center, atoms;
   parseAtomList("CENTER", center);
   if(center.size() != 1) error("one and only one atom expected in CENTER");
   parseAtomList("ATOMS", atoms);
   natoms = atoms.size();
   log.printf("  number of atoms used: %d\n", natoms);
   // put center at the beginning of atoms to be requested
   atoms.insert(atoms.begin(), center[0]);

   /*--- parse RANGES ---*/
   std::vector<double> ranges;
   parseVector("RANGES", ranges);
   if(ranges.size() != 6) error("expected 6 numbers in RANGES");
   xl = ranges[0]; xh = ranges[1];
   yl = ranges[2]; yh = ranges[3];
   zl = ranges[4]; zh = ranges[5];
   if(xl >= xh) error("x lower boundary is higher than the upper one");
   if(yl >= yh) error("y lower boundary is higher than the upper one");
   if(zl >= zh) error("z lower boundary is higher than the upper one");
   log.printf("  ranges of the grid:\n");
   log.printf("  from %f to %f in x\n", xl, xh);
   log.printf("  from %f to %f in y\n", yl, yh);
   log.printf("  from %f to %f in z\n", zl, zh);

   /*--- parse BINSIZE ---*/
   parse("BINSIZE", binsize);
   if(binsize < 0) error("binsize should be positive");
   log.printf("  binsize of grid: %f\n", binsize);
   nx = ny = nz = 0;
   for(double t = xl; t < xh + binsize; t += binsize) nx++;
   for(double t = yl; t < yh + binsize; t += binsize) ny++;
   for(double t = zl; t < zh + binsize; t += binsize) nz++;
   ngrids = nx * ny * nz;
   log.printf("  number of grid poinst: %u x %u x %u\n", nx, ny, nz, ngrids);

   /*--- finishes grid definition ---*/
   std::vector<double> rx_, ry_, rz_;
   parseVector("ROTX", rx_);
   parseVector("ROTY", ry_);
   parseVector("ROTZ", rz_);
   if(rx_.empty() and ry_.empty() and rz_.empty()) {
      double r[3] = {1.0, 0.0, 0.0};
      rx_ = std::vector<double>(r, r + 3);
      r[0] = 0.0; r[1] = 1.0; r[2] = 0.0;
      ry_=std::vector<double>(r, r + 3);
      r[0] = 0.0; r[1] = 0.0; r[2] = 1.0;
      rz_=std::vector<double>(r, r + 3);
   }
   if(rx_.size() != 3 || ry_.size() != 3 || rz_.size() != 3)
      error("rotation vectors should be 3-dimensional");
   rx = Vector(rx_[0], rx_[1], rx_[2]);
   ry = Vector(ry_[0], ry_[1], ry_[2]);
   rz = Vector(rz_[0], rz_[1], rz_[2]);
   double norm = rx.modulo2();
   if(norm < 1e-30) error("ROTX should be non-zero");
   rx /= sqrt(norm);
   norm = ry.modulo2();
   if(norm < 1e-30) error("ROTY should be non-zero");
   ry /= sqrt(norm);
   norm = rz.modulo2();
   if(norm < 1e-30) error("ROTZ should be non-zero");
   rz /= sqrt(norm);
   if(fabs(dotProduct(rx,ry)) > 1e-5) error("ROTX and ROTY are not orthogonal");
   if(fabs(dotProduct(ry,rz)) > 1e-5) error("ROTY and ROTZ are not orthogonal");
   if(fabs(dotProduct(rx,rz)) > 1e-5) error("ROTZ and ROTX are not orthogonal");
   log.printf("  the rotation matrix is :\n");
   log<<"  "<<rx[0]<<" "<<rx[1]<<" "<<rx[2]<<"\n";
   log<<"  "<<ry[0]<<" "<<ry[1]<<" "<<ry[2]<<"\n";
   log<<"  "<<rz[0]<<" "<<rz[1]<<" "<<rz[2]<<"\n";
   log.printf("\n");


   /*--- parse SWITCH ---*/
   std::string sf_dscpt, errors;
   parse("SWITCH", sf_dscpt);
   if(sf_dscpt.length() > 0) {
      sf.set(sf_dscpt, errors);
      if(errors.length() != 0) 
         error("problem reading SWITCH Keyword : " + errors);
   } else 
      error("SWITCH should be specified");
   dmax = sf.get_dmax();

   /*--- add components ---*/
   for(unsigned i = 1; i <= ngrids; ++i) {
      std::string num; Tools::convert(i, num);
      addComponentWithDerivatives("gd-" + num);
      componentIsNotPeriodic("gd-" + num);
   }

   checkRead();
   requestAtoms(atoms);
}

void GridDensity::calculate() {

   // vector of density on each grid point on this rank
   std::vector<double> densities;
   densities.resize(ngrids);
   // vector of tuple (atom_index, grid_index, [derivx, derivy, derivz])
   // where [derivx, derivy, derivz] = d grid_(grid_index) / d r_(atom_index)
   std::vector<std::tuple<unsigned, unsigned, Vector> > atom_grid_derivs;
   // d grid / d r_atom for each grid and each atom
   Matrix<Vector> dgrid_dr; 
   if(do_mpi) {
      dgrid_dr.resize(ngrids, getNumberOfAtoms());
   }


   unsigned stride = comm.Get_size();
   unsigned rank = comm.Get_rank();
   if(!do_mpi) { rank = 0; stride = 1; }
   for(unsigned iatom = rank; iatom < natoms; iatom += stride) {

      /*--- calculte atom position in grid coordinate system ---*/
      Vector distance = do_pbc?
         pbcDistance(getPosition(0), getPosition(iatom+1)):
         delta(getPosition(0), getPosition(iatom+1));
      double tx = dotProduct(distance, rx);
      double ty = dotProduct(distance, ry);
      double tz = dotProduct(distance, rz);

      /*--- determine which grids this atom should be mapped onto ---*/
      // grids[i][j][k] = (xl + i * binsize, yl + j * binsize, zl + k * binsize)
      // grids[i][j][k][0] \in [tx - dmax, tx + dmax] <=>
      // i >= (tx - dmax - xl) / binsize && 
      // i <= (tx + dmax - xl) / binsize <=>
      // i >= std::ceil( (tx - dmax - xl) / binsize) &&
      // i <= std::floor((tx + dmax - xl) / binsize)
      int il = std::ceil( (tx - dmax - xl) / binsize);
      int ih = std::floor((tx + dmax - xl) / binsize);
      int jl = std::ceil( (ty - dmax - yl) / binsize);
      int jh = std::floor((ty + dmax - yl) / binsize);
      int kl = std::ceil( (tz - dmax - zl) / binsize);
      int kh = std::floor((tz + dmax - zl) / binsize);
      il = std::max(0, il); ih = std::min(ih, nx - 1);
      jl = std::max(0, jl); jh = std::min(jh, ny - 1);
      kl = std::max(0, kl); kh = std::min(kh, nz - 1);

      /*--- loop on those grids to compute density and derivatives ---*/
      Value* v;
      for(int i = il; i <= ih; ++i) {
         double d = tx - xl - i * binsize;
         double dfx, fx = sf.calculate(d, dfx); dfx *= d;
         for(int j = jl; j <= jh; ++j) {
            d = ty - yl - j * binsize;
            double dfy, fy = sf.calculate(d, dfy); dfy *= d;
            for(int k = kl; k <= kh; ++k) {
               d = tz - zl - k * binsize;
               double dfz, fz = sf.calculate(d, dfz); dfz *= d;

               unsigned igrid = nz * ny * i + nz * j + k;
               densities[igrid] += fx * fy * fz;
               v = getPntrToComponent(igrid);
               Vector deriv, td(dfx * fy * fz, fx * dfy * fz, fx * fy * dfz);
               // rotate td into original xyz coordinate system
               for(int b = 0; b < 3; b++)
                  deriv[b] = td[0]*rx[b] + td[1]*ry[b] + td[2]*rz[b];
               if(do_mpi) {
                  dgrid_dr(igrid, 0) -= deriv;
                  atom_grid_derivs.emplace_back(iatom+1, igrid, deriv);
               } else {


                  setAtomsDerivatives(v, iatom + 1, deriv);


                  setAtomsDerivatives(v, 0, -deriv);
               }
            }
         }
      }

   } // end of loop over atoms

   if(do_mpi) {

      /*--- collect info from each rank ---*/
      for(const auto& ele : atom_grid_derivs) 
         dgrid_dr(std::get<1>(ele), std::get<0>(ele)) += std::get<2>(ele);
         comm.Sum(densities);
         comm.Sum(dgrid_dr);

      /*--- set the output and derivatives ---*/
      for(unsigned igrid = 0; igrid < ngrids; ++igrid) {
         Value* v = getPntrToComponent(igrid);
         v->set(densities[igrid]);
         for(unsigned jatom = 0; jatom < getNumberOfAtoms(); ++jatom)
            setAtomsDerivatives(v, jatom, dgrid_dr(igrid, jatom));
      }

   } else {
      for(unsigned igrid = 0; igrid < ngrids; ++igrid) 
         getPntrToComponent(igrid)->set(densities[igrid]);
   }

}

}
}
