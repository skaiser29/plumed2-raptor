/*
 * Author: Chenghan Li
 * Created: Feb 4
 * Description: Implement a NN CV based on a set of grids
 */

// ActionWithGrid is a good starting point but it seems no rotation is supported
// Modify that behaviour and use that to replace my ad hoc implementation if I 
// have time...

// read a grid input file or define grid by giving ranges and bins?

// how to parallel?
// parallel by atom in density evaluation O(Natom/Nrank)
// parallel by grid for forward propogation O(Ngrid/Nrank)

// how to compute forces on atoms?
// build a vector of tuple of 
// (atom_index, grid_index, Vector (derivx, derivy, derivz) )
// deriv[atom_index] += dcv_dgrid[grid_index] * (derivx, derivy, derivz)
// as derivx, derivy, derivz are pre-computed, this should be fast
// ----------
// deriv[0] += dcv_dgrid[grid_index] * center_derivs[grid_index]

#ifndef NDEBUG
#define __DEBUG__
#endif

#include "Colvar.h"
#include "tools/Vector.h"
#include "tools/OpenMP.h"
#include "ActionRegister.h"
#include "tools/Matrix.h"

#include <string>
#include <fstream>
#include <sstream>

namespace PLMD {

namespace colvar {

class GridNN : public Colvar {
    bool do_pbc, is_sparse;
    unsigned natoms;
    Vector rx,ry,rz;
    double xl,xh,yl,yh,zl,zh;
    double binsize;
    unsigned nx, ny, nz;
    typedef double (GridNN::*d2d_mf)(double,double&) const;
    d2d_mf fsw;
    std::vector<double (*)(double,double&)> actfs;
    std::vector<Matrix<double>> weights;
    std::vector<std::vector<double>> biases;
    std::vector<std::vector<double>> dactfs;
    double scale;

/// updated by evaluateDensity:
/// contains d grid_j / d r_atom_i
    std::vector<
       std::tuple<unsigned, unsigned, Vector>
       > atom_grid_derivs;
/// d grid_i / d r_atom_center
    std::vector<Vector> center_derivs;
/// input layer or refered as density
    std::vector<double> input_layer;

#ifdef __DEBUG__
    std::ofstream fs_debug;
#endif

public:
    explicit GridNN(const ActionOptions&);
    void calculate();

/// collection of activation function
    static double sigmoid(double x, double& df);

/// collection of switching function acting as kernels
    double rational3(double x, double& df) const;
    //double drational3(double x) const;
    double rational5(double x, double& df) const;
    //double drational5(double x) const;

    void evaluateDensity();
/// next_layer = f_act( multi(w[:-2],this_layer)+w[-1] )
    void forwardOneLayer(
          unsigned this_layer_index, const std::vector<double>& this_layer,
          std::vector<double>& next_layer, bool usempi);
    void conv3d();
    static void registerKeywords(Keywords& keys);

#ifdef __DEBUG__
    ~GridNN();
#endif
};

PLUMED_REGISTER_ACTION(GridNN, "GRIDNN")

void GridNN::registerKeywords( Keywords& keys) {
    Colvar::registerKeywords(keys);
    keys.add("atoms","ATOM", "Center of the grid");
    keys.add("atoms","SPECIES", "Atoms to be predicted interested property of");
    keys.add("optional","ROTX","Rotation vector in x axis");
    keys.add("optional","ROTY","Rotation vector in y axis");
    keys.add("optional","ROTZ","Rotation vector in z axis");
    keys.add("compulsory","RANGES","Ranges of the grid in the order of lx hx ly hy lz hz");
    keys.add("compulsory","BINSIZE","Binsize of the grid");
    keys.add("compulsory", "WEIGHTS", "File names for trained weights where the last line should be the bias");
    keys.add("compulsory", "ACTIVATIONS", "Activation functions for each layer");
    keys.add("compulsory", "SWITCH", "Switching function to map atoms onto grid");
    keys.add("optional","SCALE","Scaling factor of the final cv value");
    keys.addFlag("NOTSPARSE",false,"Atoms in the channel are not a minority of all atoms");
}

GridNN::GridNN(const ActionOptions&ao) :
PLUMED_COLVAR_INIT(ao),
do_pbc(true), is_sparse(true), scale(1.0)
#ifdef __DEBUG__
    ,fs_debug("debug.out")
#endif
{
    bool nopbc=!do_pbc, notsparse=!is_sparse;
    parseFlag("NOPBC",nopbc); parseFlag("NOTSPARSE",notsparse);
    do_pbc=!nopbc; is_sparse = !notsparse;
    parse("SCALE", scale);
    log.printf("  the final cv value will be multiplied by %f\n", scale);

    std::vector<AtomNumber> atom,species;
    parseAtomList("ATOM",atom);
    parseAtomList("SPECIES",species);
    natoms = species.size();
    log.printf("  number of atoms used: %d\n", natoms);
    if(atom.size()!=1) error("one and only one atom expected for ATOM");
    species.insert(species.begin(), atom[0]);

    //range def
    std::vector<double> ranges;
    parseVector("RANGES",ranges);
    if(ranges.size()!=6) error("expected 6 numbers in RANGES");
    xl=ranges[0]; xh=ranges[1]; 
    yl=ranges[2]; yh=ranges[3]; 
    zl=ranges[4]; zh=ranges[5];
    if(xl>=xh) error("lower boundary should be lower than the upper one in x");
    if(yl>=yh) error("lower boundary should be lower than the upper one in y");
    if(zl>=zh) error("lower boundary should be lower than the upper one in z");
    log.printf("  using grid ranging as follows (after rotation if any):\n");
    log.printf("  from %f to %f in x\n",xl,xh);
    log.printf("  from %f to %f in y\n",yl,yh);
    log.printf("  from %f to %f in z\n",zl,zh);
    std::string sw;
    parse("SWITCH", sw);
    if(sw == "rational3") {
       fsw = &GridNN::rational3;
       //dfsw = &GridNN::drational3;
    } else if(sw == "rational5") {
       fsw = &GridNN::rational5;
       //dfsw = &GridNN::drational5;
    } else error("unknown switching function " + sw);

    parse("BINSIZE", binsize);
    if(binsize < 0) error("binsize should be positive");
    log.printf("  binsize of grid: %f\n", binsize);
    nx = ny = nz = 0;
    for(double t = xl; t < xh + binsize; t += binsize) nx++;
    for(double t = yl; t < yh + binsize; t += binsize) ny++;
    for(double t = zl; t < zh + binsize; t += binsize) nz++;
    unsigned dim = nx * ny * nz;
    input_layer.resize(dim);
    center_derivs.resize(dim);
    log.printf("  number of grid points: %ux%ux%u = %u\n", nx, ny, nz, dim);
    std::vector<std::string> fnames, funcs;
    parseVector("WEIGHTS", fnames);
    parseVector("ACTIVATIONS", funcs);
    if(fnames.size() != funcs.size()) 
       error("number of weight files does not match number of activation functions");
    for(std::string fname : fnames) {
       std::ifstream fs(fname);
       if(!fs.is_open()) error("cannot open file " + fname);
       std::string line; unsigned line_count = 0;
       std::vector<double> content, bias;
       unsigned ncol;
       while(getline(fs, line)) {
          std::stringstream ss(line);
          double t;
          if(line_count == dim)
             while(ss >> t) bias.push_back(t);
          else
             while(ss >> t) content.push_back(t);
          if(!line_count) ncol = content.size();
          line_count++;
       }
       if(line_count != dim + 1) error("wrong number of lines in "+fname);
       if(content.size() != (line_count-1) * ncol || bias.size() != ncol) 
          error("incomplete weight file "+fname);
       biases.push_back(bias);
       Matrix<double> weight; weight.resize(line_count-1, ncol);
       weight = content;
       weights.push_back(weight);
       fs.close();

       // number of columns + 1 is the number of line of next weight
       dim = ncol;

       log.printf("  input weight file %s read, dims of weights = %ldx%ld"
             " dim of bias = %ld\n", 
             fname.c_str(), 
             weights.back().nrows(), weights.back().ncols(), 
             biases.back().size());
    }
    if(dim != 1) error("wrong number of columns in "+fnames.back());
    for(std::string func: funcs) {
       if(func == "sigmoid") {
          actfs.push_back(&GridNN::sigmoid);
          //dactfs.push_back(&GridNN::dsigmoid);
       } else {
          error("unkown activation function "+func);
       }
    }
    dactfs.resize(weights.size());

    addValueWithDerivatives(); setNotPeriodic();
    if(do_pbc) log.printf("  using periodic boundary conditions\n");

    //rot
    std::vector<double> rx_, ry_, rz_;
    parseVector("ROTX", rx_);
    parseVector("ROTY", ry_);
    parseVector("ROTZ", rz_);
    if(rx_.empty() and ry_.empty() and rz_.empty()) {
        double r[3] = {1.0,0.0,0.0};
        rx_=std::vector<double>(r, r + 3);
        r[0] = 0.0; r[1] = 1.0; r[2] = 0.0;
        ry_=std::vector<double>(r, r + 3);
        r[0] = 0.0; r[1] = 0.0; r[2] = 1.0;
        rz_=std::vector<double>(r, r + 3);
    }
    if(rx_.size()!=3||ry_.size()!=3||rz_.size()!=3) 
       error("rotation vectors should be 3-dimensional");
    rx = Vector(rx_[0], rx_[1], rx_[2]);
    ry = Vector(ry_[0], ry_[1], ry_[2]);
    rz = Vector(rz_[0], rz_[1], rz_[2]);
    double norm = rx.modulo2();
    if(norm<1e-30) error("ROTX should be non-zero");
    rx /= sqrt(norm);
    norm = ry.modulo2();
    if(norm<1e-30) error("ROTY should be non-zero");
    ry /= sqrt(norm);
    norm = rz.modulo2();
    if(norm<1e-30) error("ROTZ should be non-zero");
    rz /= sqrt(norm);

    if(fabs(dotProduct(rx,ry))>1e-5) error("ROTX and ROTY are not orthogonal");
    if(fabs(dotProduct(ry,rz))>1e-5) error("ROTY and ROTZ are not orthogonal");
    if(fabs(dotProduct(rz,rx))>1e-5) error("ROTZ and ROTX are not orthogonal");
    log.printf("  the rotation matrix is :\n");
    log<<"  "<<rx[0]<<" "<<rx[1]<<" "<<rx[2]<<"\n";
    log<<"  "<<ry[0]<<" "<<ry[1]<<" "<<ry[2]<<"\n";
    log<<"  "<<rz[0]<<" "<<rz[1]<<" "<<rz[2]<<"\n";
    log.printf("\n");

    checkRead();
    requestAtoms(species);
}

#ifdef __DEBUG__
GridNN::~GridNN() {
   fs_debug.close();
}
#endif

double GridNN::sigmoid(double x, double& df) {
   double t = 1.0/(1.0 + exp(-x));
   df = t * (1.0 -t);
   return t;
}

//double GridNN::dsigmoid(double x) const{
//   double t = sigmoid(x);
//   return t * (1.0 - t);
//}

double GridNN::rational3(double x, double& df) const{
   double t = fabs(x) / binsize - 0.5;
   double t2 = t * t;
   df = static_cast<int>((x > 0.0) - (x < 0.0)) * 
      (6 * t2 - 1.5) / binsize;
   return (t*(2 * t2 - 1.5) + 0.5);
}

//double GridNN::drational3(double x) const{
//   double t = fabs(x) / binsize - 0.5;
//   return (6*t*t -1.5) / binsize;
//}

double GridNN::rational5(double x, double& df) const{
   (void) x;
   df = 0.0;
   return 0.0;
}

//double GridNN::drational5(double x) const{
//   return 0.0;
//}

void GridNN::evaluateDensity() {
   input_layer = std::vector<double>(nx*ny*nz, 0.0);
   atom_grid_derivs.clear();
   std::fill(center_derivs.begin(), center_derivs.end(), Vector());
   unsigned stride = comm.Get_size();
   unsigned rank = comm.Get_rank();
   //std::vector<double> local_density(nx*ny*nz, 0.0);
   //std::vector<std::tuple<unsigned, unsigned, Vector>> local_atom_grid_derivs;
   for(unsigned i = rank; i < natoms; i += stride) {
      Vector distance = do_pbc?
         pbcDistance(getPosition(0), getPosition(i+1)):
         delta(getPosition(0), getPosition(i+1));
      // distance in grid coordinates:
      double tx = dotProduct(distance, rx);
      double ty = dotProduct(distance, ry);
      double tz = dotProduct(distance, rz);
      //determine which grid an atom should be in
      if(tx >= xh + binsize || ty >= yh + binsize || 
            tz >= zh + binsize || tx <= xl - binsize ||
            ty <= yl - binsize || tz <= zl - binsize) continue;
      unsigned ix = (unsigned) ((tx-xl)/binsize);
      unsigned iy = (unsigned) ((ty-yl)/binsize);
      unsigned iz = (unsigned) ((tz-zl)/binsize);

      for(int kx = 0; kx < 2; ++kx) {
         if(ix + kx > nx - 1) continue;
         double diff = tx - xl - (ix + kx) * binsize;
         if(diff < -binsize) continue;
         //double fx = (this->*fsw)(diff), dfx = (this->*dfsw)(diff);
         double dfx, fx = (this->*fsw)(diff, dfx);
         for(int ky = 0; ky < 2; ++ky) {
            if(iy + ky > ny - 1) continue;
            double diff = ty - yl - (iy + ky) * binsize;
            if(diff < -binsize) continue;
            //double fy = (this->*fsw)(diff), dfy = (this->*dfsw)(diff);
            double dfy, fy = (this->*fsw)(diff, dfy);
            for(int kz = 0; kz < 2; ++kz) {
               if(iz + kz > nz - 1) continue;
               double diff = tz - zl - (iz + kz) * binsize;
               if(diff < -binsize) continue;
               //double fz = (this->*fsw)(diff), dfz = (this->*dfsw)(diff);
               double dfz, fz = (this->*fsw)(diff, dfz);

               unsigned igrid = nz*ny*(ix+kx)+nz*(iy+ky)+(iz+kz);
               input_layer[igrid] += fx*fy*fz;
               Vector deriv, td(dfx*fy*fz,fx*dfy*fz,fx*fy*dfz);
               for(int b = 0; b < 3; b++)
                  deriv[b] = td[0]*rx[b] + td[1]*ry[b] + td[2]*rz[b];
               //local_atom_grid_derivs.push_back(make_tuple(i+1, igrid, deriv));
               atom_grid_derivs.emplace_back(i+1, igrid, deriv);
               center_derivs[igrid] -= deriv;

#ifdef __DEBUG__
               // check atom_grid_derivs
               if(i == 99928) {
                  printf("for atom 99929:\n");
                  printf("distance = %f %f %f\n", distance[0], distance[1], distance[2]);
                  printf("tx = %f ty = %f tz = %f\n",tx, ty, tz);
                  printf("ix  = %u iy  = %u iz  = %u\n",ix, iy, iz);
                  printf("ix' = %u iy' = %u iz' = %u\n",ix+kx, iy+ky, iz+kz);
                  printf("diffx = %f diffy = %f diffz = %f\n",tx-xl-(ix+kx)*binsize,ty-yl-(iy+ky)*binsize,tz-zl-(iz+kz)*binsize);
                  printf(" fx = %f  fy = %f  fz = %f\ndfx = %f dfy = %f dfz = %f\n",
                        fx, fy, fz, dfx, dfy, dfz);
                  printf("atom_grid_derivs = %u %u (%f %f %f)\n",
                        i+1, igrid, deriv[0], deriv[1], deriv[2]);
               }
#endif

            }
         }
      }
   }
   //mpi sum the local density
   comm.Sum(input_layer);
   //mpi sum the center_derivs
#ifdef __DEBUG__
   printf("rank %u: center_derivs[5355] = %f %f %f\n", comm.Get_rank(),
         center_derivs[5355][0], center_derivs[5355][1], center_derivs[5355][2]);
#endif
   comm.Sum(center_derivs);

   //do not mpi concatenate local atom_grid_derivs
   //instead, compute deriv[iatom] on each rank and 
   //mpi sum the final resulting deriv in compute
}

void GridNN::forwardOneLayer(
      unsigned this_layer_index, const std::vector<double>& this_layer,
      std::vector<double>& next_layer, bool usempi) 
{
   //next_layer = f_act( multi(w[:-2],this_layer)+w[-1] )
   //where w = weights[indx]; b = biased[indx]
   unsigned indx = this_layer_index;
   const Matrix<double>& w = weights[indx];
   const std::vector<double>& b = biases[indx];
   next_layer.resize(b.size(), 0.0);
   if(!usempi) 
      mult(this_layer, w, next_layer);
   else {
      //parallel by this_layer
      unsigned stride = comm.Get_size();
      unsigned rank = comm.Get_rank();
      unsigned ntasks = this_layer.size();
      //std::fill(next_layer.begin(), next_layer.end(), 0.0);
      for(unsigned i = rank; i < ntasks; i += stride) 
         for(unsigned j = 0; j < next_layer.size(); ++j) 
            next_layer[j] += w(i,j) * this_layer[i];
      comm.Sum(next_layer);
   }

   dactfs[indx].resize(next_layer.size());
   for(unsigned i = 0; i < next_layer.size(); i++) {
      next_layer[i]
         = actfs[indx](next_layer[i] + b[i], dactfs[indx][i]);
   }
}

void GridNN::calculate() {

#ifdef  __DEBUG__
   fs_debug << "grid:\n";
   for(unsigned i = 0; i < nx; i++)
      for(unsigned j = 0; j < ny; j++)
         for(unsigned k = 0; k < nz; k++)
            fs_debug << xl+i*binsize <<' '<< 
               yl+j*binsize <<' '<< zl+k*binsize <<'\n';

   fs_debug << "\n";
   for(unsigned i = 0; i < weights.size(); ++i) {
      fs_debug << "weight" << i << ":\n";
      //matrixOut(fs_debug, weights[i]);
      for(unsigned j = 0; j < weights[i].nrows(); ++j) {
         for(unsigned k = 0; k < weights[i].ncols(); ++k)
            fs_debug << weights[i](j,k) <<" ";
         fs_debug << '\n';
      }
      fs_debug << '\n';
   }
#endif

   //prepare the input layer
   dactfs.clear();
   evaluateDensity();

#ifdef __DEBUG__
   fs_debug << "density:\n";
   for(unsigned i = 0; i < nx; i++)
      for(unsigned j = 0; j < ny; j++)
         for(unsigned k = 0; k < nz; k++) {
            unsigned igrid = nz*ny*i+nz*j+k;
            fs_debug << input_layer[igrid] <<'\n';
         }
#endif

   //foward propagation
   std::vector<double> resulting_layer;
   for(unsigned i = 0; i < weights.size(); i++) {
      if(!i)
         forwardOneLayer(i, input_layer, resulting_layer, true);
      else
         forwardOneLayer(i, input_layer, resulting_layer, false);
      if(i != weights.size() - 1)
         input_layer = resulting_layer;
   }
   plumed_dbg_assert(resulting_layer.size() == 1);

   setValue(scale * resulting_layer[0]);

   //decompose forces onto atoms
   std::vector<Vector> deriv(getNumberOfAtoms());
   std::vector<double> dcv_dgrid(1, scale);
   //backward propagation
   for(unsigned ilayerc = 0; ilayerc < weights.size(); ilayerc++) {
      unsigned ilayer = weights.size() - 1 - ilayerc;
      unsigned nnode = weights[ilayer].nrows();
      unsigned next_nnode = (ilayer == weights.size() - 1) ? 1 :
         weights[ilayer+1].nrows();
      std::vector<double> prev_layer(nnode, 0.0);
      for(unsigned i = 0; i < nnode; ++i) 
         for(unsigned j = 0; j < next_nnode; ++j)
            prev_layer[i] += (dcv_dgrid[j] * dactfs[ilayer][j] * weights[ilayer](i,j));
      dcv_dgrid = prev_layer;
   }
   plumed_dbg_assert(dcv_dgrid.size() == nx*ny*nz);
#ifdef __DEBUG__
   //check the dcv_dgrid
   fs_debug << "\ndcv_dgrid:\n";
   for(double d : dcv_dgrid)
      fs_debug << d <<'\n';
#endif
   std::vector<unsigned> atom_indexes, union_atom_indexes;
   for(auto ele : atom_grid_derivs) {
      deriv[std::get<0>(ele)] += dcv_dgrid[std::get<1>(ele)]*std::get<2>(ele);
      atom_indexes.push_back(std::get<0>(ele));
   }
   unsigned stride = comm.Get_size();
   unsigned rank = comm.Get_rank(), ngrid = nx * ny * nz;

   if(is_sparse) {
      //if atoms in the channel is a minority, directly mpi sum deriv wastes
      //a lot of band width, it should be better to mk a reduced_deriv and mpi sum it 
      
      //allgatherv union_atom_index
      std::vector<int> counts(stride), displ(stride);
      int count = static_cast<int>(atom_indexes.size());
      comm.Allgather(count, counts);
      displ[0] = 0;
      for(size_t i = 1; i < stride; ++i) displ[i] = displ[i-1] + counts[i-1];
      union_atom_indexes.resize(displ[stride-1]+counts[stride-1]);
      if(union_atom_indexes.empty()) return;
      comm.Allgatherv(
            atom_indexes, union_atom_indexes, counts.data(), displ.data());
      //comm.Allgatherv(&atom_indexes[0], 
      //      &union_atom_indexes[0], &counts[0], &displ[0]);
   #ifdef __DEBUG__
      unsigned natoms_ = union_atom_indexes.size();
      std::unique(union_atom_indexes.begin(), union_atom_indexes.end());
   #endif
      plumed_dbg_assert(natoms_ == union_atom_indexes.size());
      //mk a reduced deriv, mpi sum it
      std::vector<Vector> reduced_deriv;
      for(unsigned iatom : union_atom_indexes) 
         reduced_deriv.push_back(deriv[iatom]);
      comm.Sum(reduced_deriv);
      //cp back the reduced deriv to deriv
      for(unsigned i = 0; i < union_atom_indexes.size(); ++i)
         deriv[union_atom_indexes[i]] = reduced_deriv[i];
   }

   //do the deriv of center atom
   for(unsigned igrid = rank; igrid < ngrid; igrid += stride)
      deriv[0] += dcv_dgrid[igrid] * center_derivs[igrid];

   if(is_sparse) {
      comm.Sum(deriv[0]);
   } else {
      // if atoms in the channel is not a minority of all atoms
      // i.e. deriv is not sparse then just mpi sum the whole deriv
      comm.Sum(deriv);
   }

   for(unsigned i = 0; i < getNumberOfAtoms(); ++i) 
#ifdef __DEBUG__
   {
#endif
      setAtomsDerivatives(i, deriv[i]);
#ifdef __DEBUG__
      if(deriv[i].modulo2()>1e-50)
         printf("deriv[%u] = %f %f %f\n", i,
               deriv[i][0], deriv[i][1], deriv[i][2]);
   }
#endif

}


}
}

#ifdef __DEBUG__
#undef __DEBUG__
#endif
