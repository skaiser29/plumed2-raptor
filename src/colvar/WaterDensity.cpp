/*
 * Author: Chenghan Li
 * Created: Sept 10
 * Last Modified: Sept 13
 * Description: Implement Water Density (unsatisfied with the poor performance in \
 * previous implementation)
 */

#include "Colvar.h"
#include "tools/Vector.h"
#include "tools/NeighborListMPI.h"
#include "tools/OpenMP.h"
#include "ActionRegister.h"

#include <string>

namespace PLMD {

class NeighborListMPI;

namespace colvar {

class WaterDensity : public Colvar {
    bool do_pbc;
    NeighborListMPI* nl;
    bool invalidateList;
    bool firsttime;
    std::vector<double> rx,ry,rz, sigmas;
    double nl_cut,xl,xh,yl,yh,zl,zh;
    unsigned nl_st;

public:
    explicit WaterDensity(const ActionOptions&);
    ~WaterDensity();
    void calculate();
    void prepare();
// return fsw(distance)
    double pairing(double distance, double &dfunc, double sigma, double lowb, double highb) const;
// return \prod_i fsw(distance[i]) i.e. the density
    double pairing(const Vector& distance, Vector&dfunc, const std::vector<double>&sigmas) const;
//double pairing(const Vector& distance, Vector&dfunc, const std::vector<double>& sigmas,bool p) const; //@@@
    static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(WaterDensity, "WATERDENSITY")

void WaterDensity::registerKeywords( Keywords& keys) {
    Colvar::registerKeywords(keys);
    //keys.addFlag("NOPBC",false,"Do not consider PBC");
    keys.add("atoms","ATOM", "The reference atom of the box");
    keys.add("atoms","SPECIES", "Atoms to be calculated density of");
    keys.add("compulsory","NL_CUTOFF","The cutoff for the neighbour list");
    keys.add("compulsory","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
    keys.add("optional","ROTX","Rotation vector in x axis");
    keys.add("optional","ROTY","Rotation vector in y axis");
    keys.add("optional","ROTZ","Rotation vector in z axis");
    keys.add("compulsory","BOX","the lower boundary relative to the coordinate of the atom");
    keys.add("compulsory","SIGMAS","the width of the function to be used for density estimation");
}

WaterDensity::WaterDensity(const ActionOptions&ao) :
PLUMED_COLVAR_INIT(ao),
do_pbc(true),
invalidateList(true),
firsttime(true)
{
    bool nopbc=!do_pbc;
    parseFlag("NOPBC",nopbc);
    do_pbc=!nopbc;

    std::vector<AtomNumber> atom,species;
    parseAtomList("ATOM",atom);
    parseAtomList("SPECIES",species);
    log.printf("  number of atoms used: %d\n",species.size());
    if(atom.size()!=1) error("only one atom expected for ATOM");

    //sigmas
    parseVector("SIGMAS",sigmas);
    if(sigmas.size()!=3) error("expected 3 numbers for SIGMAS");
    if(sigmas[0]<=0.0 or sigmas[1]<=0.0 or sigmas[2]<=0.0)
        error("SIGMAS should be explicitly specified and non-negative");
    log.printf("  sigmas in three dimensions: %f %f %f\n",
            sigmas[0],sigmas[1],sigmas[2]);

    //box def
    std::vector<double> solbox;
    parseVector("BOX",solbox);
    if(solbox.size()!=6) error("expected 6 numbers in BOX");
    xl=solbox[0]; xh=solbox[1]; yl=solbox[2]; yh=solbox[3]; zl=solbox[4]; zh=solbox[5];
    if(xl>=xh) error("lower boundary should be lower than the upper bound in x");
    if(yl>=yh) error("lower boundary should be lower than the upper bound in y");
    if(zl>=zh) error("lower boundary should be lower than the upper bound in z");
    log.printf("  using a box ranging:\n");
    log.printf("  from %f to %f in x\n",xl,xh);
    log.printf("  from %f to %f in y\n",yl,yh);
    log.printf("  from %f to %f in z\n",zl,zh);

    //neighbor list stuff
    parse("NL_CUTOFF",nl_cut);
    if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
    if(nl_cut<=fabs(xl)+2*sigmas[0] or nl_cut<=fabs(xh)+2*sigmas[0] or 
            nl_cut<=fabs(yl)+2*sigmas[1] or nl_cut<=fabs(yh)+2*sigmas[1] or
            nl_cut<=fabs(zl)+2*sigmas[2] or nl_cut<=fabs(zh)+2*sigmas[2])
        error("NL_CUTOFF should be larger than box size plus 2*sigma");
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");

    addValueWithDerivatives(); setNotPeriodic();
    nl= new NeighborListMPI(atom,species,false,do_pbc,getPbc(),nl_cut,nl_st,&comm);

    requestAtoms(nl->getFullAtomList());

    if(do_pbc) log.printf("  using periodic boundary conditions\n");
    log.printf("  using neighbor lists with\n");
    log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);

    //rot
    parseVector("ROTX", rx);
    parseVector("ROTY", ry);
    parseVector("ROTZ", rz);
    if(rx.size()==0 and ry.size()==0 and rz.size()==0) {
        double r[3] = {1.0,0.0,0.0};
        rx=std::vector<double>(r, r + 3);
        r[0] = 0.0; r[1] = 1.0; r[2] = 0.0;
        ry=std::vector<double>(r, r + 3);
        r[0] = 0.0; r[1] = 0.0; r[2] = 1.0;
        rz=std::vector<double>(r, r + 3);
    }
    if(rx.size()!=3||ry.size()!=3||rz.size()!=3) error("rotation vectors should be 3-dimensional");
    double norm=0;
    for(unsigned i=0;i<3;++i) norm+=rx[i]*rx[i];
    if(norm<1e-30) error("ROTX should be non-zero");
    for(unsigned i=0;i<3;++i) rx[i]/=sqrt(norm);
    norm=0;
    for(unsigned i=0;i<3;++i) norm+=ry[i]*ry[i];
    if(norm<1e-30) error("ROTY should be non-zero");
    for(unsigned i=0;i<3;++i) ry[i]/=sqrt(norm);
    norm=0;
    for(unsigned i=0;i<3;++i) norm+=rz[i]*rz[i];
    if(norm<1e-30) error("ROTZ should be non-zero");
    for(unsigned i=0;i<3;++i) rz[i]/=sqrt(norm);

    double inner=0;
    for(unsigned i=0;i<3;++i) inner+=rx[i]*ry[i];
    if(fabs(inner)>1e-5) error("ROTX and ROTY are not orthogonal");
    inner=0;
    for(unsigned i=0;i<3;++i) inner+=ry[i]*rz[i];
    if(fabs(inner)>1e-5) error("ROTY and ROTZ are not orthogonal");
    inner=0;
    for(unsigned i=0;i<3;++i) inner+=rz[i]*rx[i];
    if(fabs(inner)>1e-5) error("ROTZ and ROTX are not orthogonal");
    log.printf("  the rotation matrix is :\n");
    log<<"  "<<rx[0]<<" "<<rx[1]<<" "<<rx[2]<<"\n";
    log<<"  "<<ry[0]<<" "<<ry[1]<<" "<<ry[2]<<"\n";
    log<<"  "<<rz[0]<<" "<<rz[1]<<" "<<rz[2]<<"\n";
    log.printf("\n");

    checkRead();
}

WaterDensity::~WaterDensity() {
    delete nl;
}

void WaterDensity::prepare() {
  if(nl->getStride()>0){
    if(firsttime || (getStep()%nl->getStride()==0)){
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    }else{
      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

//calculator
void WaterDensity::calculate() {
    double density=0.0;
    std::vector<Vector> deriv(getNumberOfAtoms());
    if(nl->getStride()>0 and invalidateList) {
        nl->update(getPositions());
//if(comm.Get_rank()==1) printf("# of all pairs = %d\n",nl->size());
//if(comm.Get_rank()==1)  {
//for(int i=0; i<nl->size(); ++i) 
//printf("%d ", getAbsoluteIndex(nl->getClosePair(i).second).serial()); //@@@
//printf("\n");
//}
    }

    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    unsigned nt=OpenMP::getNumThreads();
    const unsigned nn=nl->size();

    if(nt*stride*10>nn) nt=nn/stride/10;
    if(nt==0) nt=1;

#pragma omp parallel num_threads(nt)
    {
        std::vector<Vector> omp_deriv(getPositions().size());
#pragma omp for reduction(+:density) nowait
        for(unsigned int i=rank;i<nn;i+=stride) {
            Vector distance, tmp(0.0,0.0,0.0), dfunc;
            unsigned i0=nl->getClosePair(i).first;
            unsigned i1=nl->getClosePair(i).second;

            if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;

            if(do_pbc) {
                distance=pbcDistance(getPosition(i0),getPosition(i1));
            } else {
                distance=delta(getPosition(i0),getPosition(i1));
            }
            //rotate distance before use!!
            for(unsigned i=0;i<3;++i) {
                tmp[0]+=distance[i]*rx[i];
                tmp[1]+=distance[i]*ry[i];
                tmp[2]+=distance[i]*rz[i];
            }
            distance = tmp;
            //dfunc will be \partial density/\partial distance
//int index=getAbsoluteIndex(i1).serial(); //@@@
//if(index==70)
            density += pairing(distance,dfunc, sigmas);
//else
//            density += pairing(distance,dfunc, sigmas, false);

//if(index==70 or index==145 or index==469 or index==574 or index==1849) {
//printf("%d%18.12f from rank %d\n",index,pairing(distance,dfunc, sigmas),comm.Get_rank()); //@@@
//}
            Vector dd; 
            //rotate dfunc back to get \partial density/\partial r_{i1}
            for(unsigned int i=0; i<3; ++i) {
                dd[i]=dfunc[0]*rx[i]+dfunc[1]*ry[i]+dfunc[2]*rz[i];
            }

            if(nt>1) {
                omp_deriv[i0]-=dd;
                omp_deriv[i1]+=dd;
            } else {
                deriv[i0]-=dd;
                deriv[i1]+=dd;
            }
        }
#pragma omp critical
        //reduce to get deriv from omp_deriv
        if(nt>1) {
            for(int i=0;i<getPositions().size();++i) deriv[i]+=omp_deriv[i];
        }

    } // end of pragma omp parallel
    
    comm.Sum(density); //get density from all ranks
    if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());

    for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]);
    setValue(density);

//if(comm.Get_rank()==1) //@@@
//for(int i=0; i<deriv.size();++i) {
//int index=getAbsoluteIndex(i).serial();
//if(index==70 or index==145 or index==469 or index==574 or index==1849) {
//printf("deriv %d%18.12f%18.12f%18.12f\n",index,deriv[i][0],deriv[i][1],deriv[i][2]); //@@@
//printf("pos %d%18.12f%18.12f%18.12f\n",index,getPosition(i)[0],getPosition(i)[1],getPosition(i)[2]);
//}
//}
}

double WaterDensity::pairing(double distance, double &dfunc, double sigma, double lowb, double highb) const {
    double cv,ref,middle,xs,xs3,xs5,xs6,xs12,f;
    ref=(highb-lowb)/2.0;
    middle=ref+lowb;
    cv=distance-middle; //no need for pbc considerations here as distance is pbc'ed
    xs=(fabs(cv)-ref)/sigma+1.0; //this is the x in Yuxing's code

    if(xs<0) { dfunc=0.0; return 1.0;}
    else {
        if(fabs(xs-1.0)<1e-5) f=0.5;
        else {
            xs3=xs*xs*xs; xs6=xs3*xs3; xs12=xs6*xs6;
            f=(1-xs6)/(1-xs12);
        }
        xs5=xs*xs*xs3;
        dfunc=-6/sigma*xs5/(1+xs6)/(1+xs6);
    }
    if(cv<0) dfunc=-dfunc;
    return f;
}

double WaterDensity::pairing(const Vector& distance, Vector&dfunc, const std::vector<double>& sigmas) const {
    double fx,fy,fz, dfx, dfy, dfz;
    fx=pairing(distance[0], dfx, sigmas[0], xl, xh);
    fy=pairing(distance[1], dfy, sigmas[1], yl, yh);
    fz=pairing(distance[2], dfz, sigmas[2], zl, zh);
    dfunc[0]=dfx*fy*fz;
    dfunc[1]=fx*dfy*fz;
    dfunc[2]=fx*fy*dfz;
    return fx*fy*fz;
}

//double WaterDensity::pairing(const Vector& distance, Vector&dfunc, const std::vector<double>& sigmas,bool p) const {
//    double fx,fy,fz, dfx, dfy, dfz;
//    fx=pairing(distance[0], dfx, sigmas[0], xl, xh);
//    fy=pairing(distance[1], dfy, sigmas[1], yl, yh);
//    fz=pairing(distance[2], dfz, sigmas[2], zl, zh);
//    dfunc[0]=dfx*fy*fz;
//    dfunc[1]=fx*dfy*fz;
//    dfunc[2]=fx*fy*dfz;
//    if(p) printf("x=%f,y=%f,z=%f,fx=%f,fy=%f,fz=%f,dfx=%f,dfy=%f,dfz=%f\n",distance[0],distance[1],distance[2],fx,fy,fz,dfx,dfy,dfz); //@@@
//    return fx*fy*fz;
//}

}
}

