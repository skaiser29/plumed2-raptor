/*
    created: Sept 13 2017
    last modified: Sept 13 2017
    author: Chenghan Li

    Implementattion of the Proton Indicator in Lin's paper. 
    No box derivatives included, not suitable for NPT ensemble
*/
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "../tools/OpenMP.h"
#include "mpi.h"

#include<vector>
#include<cmath>
#include<algorithm>
#include<numeric>

using namespace std;

namespace PLMD {
namespace vatom{

class ProtonIndicator : public ActionWithVirtualAtom {
//private:
    bool do_pbc,firsttime; //firsttime=if the first time of this class is used
    int do_serial;
    unsigned natomh,natomx,donorIndex;
    Vector donorPos;
    vector<int> wx; //number of hydrogens bonded when protonated , i.e. M in Lin's paper
    //the idea of using reducedAtomList is abandoned because building it in MPI
    // will shuffle the list and thus requires a reorder
    //vector<AtomNumber> fullAtomList, reducedAtomList;
    //Vector lastpos; //store the last pos of PI helping determine donor atom
    Vector lastpos; //store the last pos of DONOR atom NOT PI!
    double r0,rlist,rho0, thres;
    //unsigned int getDonorIndex();
    double updateDonorIndex();
    //store the indexes of possible acceptors in fullatomlist
    //store the indexes of candidate hydrogens for bonding H of donor in fullatomlist
    //i.e. X and H within rlist of donor
    vector<unsigned> acceptorList, bondedHList, candidateHList;
    void setDonorPos();
    double gSW(double rho, double &dgsw) const;  //switching function g(x(\rho))
    void plainsort(vector<double>& metric, vector<Vector>& list);
    vector<unsigned> firstShell; //store the indexes in fullAtomList of the first sol
                                //shell oxygen
    vector<double> rH2_;
public:
    explicit ProtonIndicator(const ActionOptions&ao);
    //void prepare();
    void calculate();
    static void registerKeywords( Keywords& keys );
    //void updateAtomList();
    void updateNList();

//@@@
vector<AtomNumber> fullAtomList;
};

PLUMED_REGISTER_ACTION(ProtonIndicator,"PROTONINDICATOR")

void ProtonIndicator::registerKeywords(Keywords& keys) {
    ActionWithVirtualAtom::registerKeywords(keys);
    keys.add("atoms","ATOMH","Hydrogen atoms");
    keys.add("atoms","ATOMX","Heavy atoms");
    keys.add("atoms","INIT","Initial serial number of the donor");
    keys.add("compulsory","WX","The weights for heavy atoms, typically should be negative");
    keys.add("optional","R0","r_\\text{DH}^0");
    keys.add("optional","RLIST","r_\\text{LIST}");
    //keys.add("optional","NSTRIDE","The actual atom list is updated every this steps");
    keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances ");
    keys.add("optional","PARALLEL","level of parallelization");
    keys.add("optional","FIRSTSHELL",
             "oxygens whose rho is over this value will be regarded as first shell atoms");
}

ProtonIndicator::ProtonIndicator(const ActionOptions&ao):
    Action(ao),
    ActionWithVirtualAtom(ao),
    do_pbc(true),r0(1.0),rlist(14.0), thres(0.5),
    firsttime(true),do_serial(1)//, nstride(2)
{
    /*--- parse ---*/
    vector<AtomNumber> atomH,atomX,initDonorSerial;
    parseAtomList("ATOMH",atomH);
    parseAtomList("ATOMX",atomX);
    parseAtomList("INIT",initDonorSerial);
    parseVector("WX",wx);
    bool nopbc;
    parseFlag("NOPBC",nopbc);
    int do_para=1;
    parse("PARALLEL",do_para);
    do_pbc=!nopbc;
    parse("R0",r0);parse("RLIST",rlist);parse("FIRSTSHELL",thres);
    //parse("NSTRIDE",nstride);
    checkRead();

    /*--- error ---*/
    natomh= static_cast<unsigned int>(atomH.size());
    natomx= static_cast<unsigned int>(atomX.size());
    if(natomh==0) error("Hydrogen atoms should be specified");
    if(atomX.empty()) error("Heavy atoms should be specified");
    if(wx.size()!=atomX.size()) error("Number of WX should be the same as number of ATOMX");
    for(unsigned int i=0; i<natomx; i++) {
        if(wx[i]<=0)
            error("WX should be positive numbers");
    }
    if(initDonorSerial.size()!=1) error("Only one atom expected in INIT");
    if(find(atomX.begin(),atomX.end(),initDonorSerial[0]) == atomX.end()) {
        error("INIT is not a member of ATOMX");
    }
    if(r0<=0.0) error("R0 should be positive");
    if(rlist<=0.0) error("RLIST should be positive");
    if(r0/rlist>0.4) error("Too small RLIST");
    do_serial=2-do_para;
    if(do_serial!=0 and do_serial!=1 and do_serial!=2)
        error("Available PARALLEL options: 0, 1 and 2");
    //if(nstride<=0) error("NSTRIDE should be positive");

    /*--- past processing ---*/
    atomH.insert(atomH.end(),atomX.begin(),atomX.end()); //concatenate
    donorIndex = static_cast<unsigned int>(
              find(atomH.begin(), atomH.end(), initDonorSerial[0])
                     - atomH.begin());
    rho0=r0/rlist;

    /*--- write log ---*/
    log.printf("  %d hydrogen atoms and %d heavy atoms involved\n",
    natomh, static_cast<int>(atomX.size()));
    log<<"  serial number of hydrogen atoms:\n  ";
    for(unsigned int i=0; i<natomh; ++i) {
        log.printf("%d ", atomH[i].serial());
    }
    log<<"\n  serial number of heavy atoms:\n  ";
    for(unsigned int i=0; i<natomx; ++i) {
        log.printf("%d ", atomH[i+natomh].serial());
    }
    log<<"\n  serial number of initial donor: ";
    log.printf("%d ", initDonorSerial[0].serial());
    log <<"\n  weights of heavy atoms:\n  ";
    for(unsigned int i=0; i<natomx; ++i) {
        log<<wx[i]<<" ";
    }
    log<<"\n";
    log<< "  r0: " << r0 << " rlist: " << rlist <<"\n";
    log<<"  atoms whose rho is over "<< thres <<" are the first shell atoms\n";
    if(!do_pbc){
      log<<"  WARNING: PBC will be ignored; could cause severe problem\n";
    } else {
      log<<"  PBC will be respected\n";
    }
    if(do_serial==2) {
        log<<"  do the calculation in serial\n";
    }
    else if(do_serial==1) {
        log<<"  do the neighbor list in parallel\n";
    }
    else {
        log<<"  do the neighbor list and caculation in parallel\n";
    }
    log<<"  WARNING: charge and mass of this virtual atom are set to be 0 and 1 respectively\n";
    //log<<"  WARNING: PBC was not tested extensively, be careful...\n";
    //initialize indexes, allocate positions, masses, ...
    requestAtoms(atomH); //atomH now is the fullAtomList

//@@@
fullAtomList = atomH;
}

void ProtonIndicator::setDonorPos() {
    donorPos=getPosition(donorIndex);
}

//void ProtonIndicator::prepare() {
//    if(firsttime|| (getStep()%nstride==0)) {
//        requestAtoms(fullAtomList);
//        invalidateList=true; //in calculate(), true makes updateAtomList involved
//    } else {
//        requestAtoms(reducedAtomList);
//        invalidateList=false;
//        if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps");
//    }
//    if(getExchangeStep()) firsttime=true;
//}

//update the reducedAtomList to make it validate
//void ProtonIndicator::updateAtomList() {
//    //this function always sees fullAtomList
//}

//update the acceptorList and candidateHList
void ProtonIndicator::updateNList() {
    acceptorList.clear();
    candidateHList.clear();
    unsigned nt,stride,rank, natom=getNumberOfAtoms();
    double rlist2 = rlist * rlist;
    vector<unsigned> local_h, local_a; //lists on each rank
    if(do_serial==2) {
        nt = 1;
        stride = 1;
        rank = 0;
    } else {
        nt = OpenMP::getNumThreads();
        stride = comm.Get_size();
        rank = comm.Get_rank();
        if(nt*stride*10>natom) nt=natom/stride/10;
        if(nt==0) nt=1;
    }
#pragma omp parallel num_threads(nt)
    {
        vector<unsigned> omp_h, omp_a;
#pragma omp for nowait
        for (unsigned i = rank; i < getNumberOfAtoms(); i+=stride) {
            Vector distance;
            if (do_pbc)
                distance = pbcDistance(donorPos, getPosition(i));
            else
                distance = delta(donorPos, getPosition(i));
            double value=modulo2(distance);
            if(value<=rlist2) {
                if(nt>1)
                    if(i<natomh) omp_h.push_back(i);
                    else omp_a.push_back(i);
                else
                    if(do_serial==2) {
                        if(i<natomh)
                            candidateHList.push_back(i);
                        else
                            acceptorList.push_back(i);
                    } else {
                        if(i<natomh) local_h.push_back(i);
                        else local_a.push_back(i);
                    }

            }
        }
#pragma omp critical
        if(nt>1) {
            local_a.insert(local_a.end(),omp_a.begin(),omp_a.end());
            local_h.insert(local_h.end(),omp_h.begin(),omp_h.end());
        }
    } //end of omp parallel
    if(do_serial!=2) {
        int* counts_h=(int*)malloc(stride* sizeof(int));
        int* displ_h=(int*)malloc(stride*sizeof(int));
        int* counts_a=(int*)malloc(stride* sizeof(int));
        int* displ_a=(int*)malloc(stride*sizeof(int));
        for(int i=0;i<stride;++i) {
            if(rank==i) {
                counts_a[i]=local_a.size();
                counts_h[i]=local_h.size();
            }
            //MPI_Bcast(&counts_a[i],1,MPI_INT,i,comm->Get_world().Get_comm());
            //MPI_Bcast(&counts_h[i],1,MPI_INT,i,comm->Get_world().Get_comm());
            comm.Bcast(counts_a[i],i);
            comm.Bcast(counts_h[i],i);
        }
        displ_a[0]=0; displ_h[0]=0;
        for(int i=1;i<stride;++i) {
            displ_a[i]=displ_a[i-1]+counts_a[i-1];
            displ_h[i]=displ_h[i-1]+counts_h[i-1];
        }
        int total_a=displ_a[stride-1]+counts_a[stride-1]; //total # of a among all ranks
        int total_h=displ_h[stride-1]+counts_h[stride-1]; //total # of h
        unsigned* list_a=(unsigned*)malloc(sizeof(unsigned)*total_a);
        unsigned* list_h=(unsigned*)malloc(sizeof(unsigned)*total_h);
        MPI_Allgatherv(&local_a[0],local_a.size(),MPI_UNSIGNED,&list_a[0],
                counts_a, displ_a, MPI_UNSIGNED, comm.Get_world().Get_comm()
                );
        MPI_Allgatherv(&local_h[0],local_h.size(),MPI_UNSIGNED,&list_h[0],
                counts_h, displ_h, MPI_UNSIGNED, comm.Get_world().Get_comm()
                );
        acceptorList=vector<unsigned>(list_a,list_a+total_a);
        candidateHList=vector<unsigned>(list_h,list_h+total_h);
        free(counts_h); free(displ_h);
        free(counts_a); free(displ_a);
        free(list_a); free(list_h);
    }
//    int maxh=wx[donorIndex-natomh];
//    for(unsigned i=0;i<candidateHList.size();++i) {
//        if(i<maxh)
//            bondedHList.push_back(candidateHList[i]);
//        else
//
//    }
}

double ProtonIndicator::updateDonorIndex() {
    //return the distance between this donor and the excess proton
    //firstShell implementation & threshold
    if(firsttime) {
        firsttime = false;
        return 0.0;
    }
    else {
        if(firstShell.empty()) return 0.0;
        double bondLength=0.0;
        Vector protonPos;

//@@@
//printf("#0 step %ld from rank %d \n",getStep(),comm.Get_rank());

        for(int m=0; m<wx[donorIndex-natomh]; ++m) {
            if(rH2_[m]>bondLength) {
                bondLength=rH2_[m];
                protonPos=getPosition(bondedHList[m]);
            }
        }

//@@@
//printf("#1 step %ld from rank %d donorIndex=%d donorM=%d bondedHList=",getStep(),comm.Get_rank(),fullAtomList[donorIndex].serial(),wx[donorIndex]);
//for(int m=0;m<wx[donorIndex-natom];++m) {
//   printf("%d ",fullAtomList[bondedHList[m]].serial());
//}
//printf("\n");

        bool hop=false;
        double d2;
        for(unsigned j=0; j<firstShell.size(); ++j) {
            Vector distance=do_pbc?pbcDistance(protonPos,getPosition(firstShell[j])):
                                    delta(protonPos,getPosition(firstShell[j]));
            d2=distance.modulo2();
            if(d2<bondLength) {
                bondLength=d2;
                hop=true; donorIndex=firstShell[j];
            }
        }

//@@@
//printf("#1 step %ld from rank %d donorIndex=%d\n",getStep(),comm.Get_rank(),fullAtomList[donorIndex].serial());

        if(!hop) d2=bondLength;
        return d2;
    }
}



//double ProtonIndicator::updateDonorIndex() {
//    //return the distance between lastpos and this donor
//    if(firsttime)  {
//        //if this is the firsttime, do nothing
//        firsttime = false;
//        //return donorIndex;
//        return 0.0;
//    } else {
//        unsigned nacceptor=acceptorList.size();
//        unsigned rank,stride;
////        //I assume nacceptor is not very large so no need for multithreading
//        //not parallel for O(N_neighbor) operation
//        //if(true) {
//        if(do_serial) {
//            rank=0; stride=1;
//        } 
//        else {
//            rank=comm.Get_rank();
//            stride=comm.Get_size();
//            if(stride*10>nacceptor) {
//                rank=0; stride=1;
//            }
//        }
//        double mind2_local;
//        Vector distance;
//        if(do_pbc) 
//            distance=pbcDistance(lastpos,getPosition(donorIndex));
//        else
//            distance=delta(lastpos,getPosition(donorIndex));
//        mind2_local=modulo2(distance);
//        unsigned donorIndex_local=donorIndex;
//        for(unsigned i=rank; i<nacceptor; i+=stride) {
//            if(do_pbc)
//                distance=pbcDistance(lastpos,getPosition(acceptorList[i]));
//            else
//                distance=delta(lastpos,getPosition(acceptorList[i]));
//            double value=modulo2(distance);
//            if(value<mind2_local) {
//                mind2_local=value; donorIndex_local=acceptorList[i];
//            }
//        }
//        if(stride==1) {
//            donorIndex=donorIndex_local;
//            return mind2_local;
//        }
//        //else stride!=1
//        double* mind2_list=(double*) malloc(sizeof(double)*stride);
//        for(unsigned i=0;i<stride;i++) {
//            if(rank==i) 
//                mind2_list[i]=mind2_local;
//            comm.Bcast(mind2_list[i],i);
//        }
//        double mind2=mind2_list[0];
//        unsigned rankMinMin=0;
//        for(unsigned i=1;i<stride;i++) {
//            if(mind2_list[i]<mind2) {
//                mind2=mind2_list[i];
//                rankMinMin=i;
//            }
//        }
//        if(rank==rankMinMin) {
//            donorIndex=donorIndex_local;
//        }
//        comm.Bcast(donorIndex,rankMinMin);
//        free(mind2_list);
//        return mind2;
//    }
//}

double ProtonIndicator::gSW(double rho,double &dgsw) const {
    double x=1.0-(rho-rho0)/(0.5-rho0);
    if(x>=1) {
        dgsw=0.0; return 0;
    }
    else if(x<=0.0) {
        dgsw=0.0; return 1.0;
    }
    dgsw=-30*x+60;
    dgsw=(dgsw*x-30)*x*x/(rho0-0.5);
    double g=-6*x+15;
    return((g*x-10)*x*x*x+1.0);
}

//a template maybe better?
//any sort algorihtm in src/tools?
void ProtonIndicator::plainsort(vector<double>& metric, vector<Vector>& list) {
    unsigned size=metric.size();
    plumed_assert(size==list.size());
    plumed_assert(size=bondedHList.size());
    for(unsigned i=1;i<size;++i) 
        for(unsigned j=0;j<size-i;++j) {
            if(metric[j]>metric[j+1]) {
                double tmp=metric[j];
                metric[j]=metric[j+1]; metric[j+1]=tmp;
                Vector ttmp=list[j];
                list[j]=list[j+1]; list[j+1]=ttmp;
                unsigned tttmp=bondedHList[j];
                bondedHList[j]=bondedHList[j+1]; bondedHList[j+1]=tttmp;
            }
        }
}

void ProtonIndicator::calculate() {
    Vector pos;
    //rho=\rho_{mj}, x=x(\rho_{mj}), gi=g_\text{I}
    double rho, gI=0.0;
    double tmp;
    unsigned natom=getNumberOfAtoms();
    vector<Tensor> deriv(natom);

    updateDonorIndex(); //donorIndex is the index in fullatomlist

//@@@
#ifndef NDEBUG
//if(comm.Get_rank()==getStep()%4) {
log.printf("step = %ld donorIndex = %d Pos = %f %f %f\n",getStep(), fullAtomList[donorIndex].serial(), getPosition(donorIndex)[0], getPosition(donorIndex)[1],getPosition(donorIndex)[2]);
//printf("step = %ld donorIndex = %d Pos = %f %f %f\n",getStep(), fullAtomList[donorIndex].serial(), getPosition(donorIndex)[0], getPosition(donorIndex)[1],getPosition(donorIndex)[2]);
//}
#endif

    setDonorPos();
    updateNList(); //candidateH and acceptors obtained

//@@@
//if(comm.Get_rank()==(getStep()+1)%4) {
//printf("acceptor: (%d)\n",acceptorList.size());
//for(int i=0;i<acceptorList.size();++i) {
//printf("%d ",fullAtomList[acceptorList[i]].serial());
//}
//printf("\n");
//}

//@@@
//if(comm.Get_rank()==(getStep()+1)%4) {
//printf("candidateHList: (%d)\n",candidateHList.size());
//for(int i=0;i<candidateHList.size();++i) {
//printf("%d ",fullAtomList[candidateHList[i]].serial());
//}
//printf("\n");
//}

    //substract every position by r_D
    unsigned nh=candidateHList.size(), nx=acceptorList.size();
    vector<Vector> rX(nx); //X_{X_i} - X_\text{D}
    vector<double> rX2(nx); // (X_{X_i} - X_\text{D})^2
    int *nxcounts,*nxdispl;
    int nxbin,nxmax,nxlast_bin;
    unsigned stride,rank;
    if(do_serial) {
        stride=1;rank=0;
        nxcounts=NULL;
        nxdispl=NULL;
    } 
    else {
        stride=comm.Get_size();
        rank=comm.Get_rank();
        nxbin=nx/stride;
        nxmax=(rank==stride-1)?nx:(rank+1)*nxbin;
        nxlast_bin=nx-(stride-1)*nxbin;
        nxcounts=(int*) malloc(stride* sizeof(int));
        nxdispl=(int*) malloc(stride* sizeof(int));
        for(unsigned i=0;i<stride-1;++i) nxcounts[i]= static_cast<int>(nxbin);
        nxcounts[stride-1]=nxlast_bin;
        nxdispl[0]=0;
        for(int i=1;i<stride;++i) {
            nxdispl[i]=nxdispl[i-1]+nxcounts[i-1];
        }
    }
    //if(do_serial or nx/stride<50) {
        for(unsigned i=0; i<nx; ++i) {
            if(do_pbc) {
                rX[i] = pbcDistance(donorPos, getPosition(acceptorList[i]));
            }
            else {
                rX[i] = delta(donorPos, getPosition(acceptorList[i]));
            }
            rX2[i]=(modulo2(rX[i]));
        }
    //} else { //do parallel
    //    vector<Vector> rX_local(nxbin); 
    //    vector<double> rX2_local(nxbin);
    //    if(nxbin!=nxlast_bin and rank==stride-1) {
    //        rX_local.resize(nxlast_bin);
    //        rX2_local.resize(nxlast_bin);
    //    }
    //    for(unsigned j=rank*nxbin;j<nxmax;j+=1) {
    //        unsigned jj=j-rank*nxbin;
    //        if(do_pbc) {
    //            rX_local[jj] = pbcDistance(donorPos,
    //                    getPosition(acceptorList[j]));
    //        }
    //        else {
    //            rX_local[jj] = delta(donorPos, getPosition(acceptorList[j]));
    //        }
    //        rX2_local[jj]=(modulo2(rX_local[jj]));
    //    }
    //    comm.Allgatherv(rX_local,rX,nxcounts,nxdispl);
    //    comm.Allgatherv(rX2_local,rX2,nxcounts,nxdispl);
    //    rX_local.clear();
    //    rX2_local.clear();
    //}

//@@@
//if(comm.Get_rank()==(getStep()+1)%4) {
//for(unsigned i=0;i<nx;++i)
//    printf("rX[%d],rX2= %f %f %f, %f\n",fullAtomList[acceptorList[i]].serial(),rX[i][0],rX[i][1],rX[i][2],rX2[i]);
//}

    //find the bonded hydrogens to the donor
    int donorM=wx[donorIndex-natomh]; // M value of the donor

//@@@
//printf("donorM=%d\n",donorM);

    vector<Vector> rH(donorM);
    bondedHList.clear();
    //vector<unsigned> bondedHList(donorM);
    vector<double> rH2(donorM);
    for(unsigned i=0; i<nh; ++i) {
        if(i<donorM) {
            rH[i]=do_pbc?
                pbcDistance(donorPos,getPosition(candidateHList[i])):
                delta(donorPos,getPosition(candidateHList[i]));
            bondedHList.push_back(candidateHList[i]);
            rH2[i]=modulo2(rH[i]);
            if(i==donorM-1) plainsort(rH2,rH);
        } else {
            Vector distance=do_pbc?
                pbcDistance(donorPos,getPosition(candidateHList[i])):
                delta(donorPos,getPosition(candidateHList[i]));
            double value=modulo2(distance);
            if(value<rH2.back()) {
                rH.back()=distance; rH2.back()=value; 
                bondedHList.back()=candidateHList[i];
                plainsort(rH2,rH);
            }
        }
    }

//@@@
//if(comm.Get_rank()==(getStep()+1)%4) {
//printf("bondedHList: ");
//for(int i=0;i<donorM;i++) {
//printf("%d ",fullAtomList[bondedHList[i]].serial());
//}
//printf("\n");
//}


    /*--- compute pos and derivatives ---*/
    //vector<double> maxrho(nx,0.0); //max rho for each acceptor
    vector<int> isfirstshell(nx,0);
    //pos
    //g_\text{I} requires O(donorM*nx), parallelize!
    unsigned nloop=donorM*nx; //nx=J
    vector<double> g(nloop),dg(nloop),rhos(nloop);
    int *counts,*displ;
    int bin,max,last_bin;
    if(do_serial) {
        //stride=1;rank=0;
        counts=NULL;
        displ=NULL;
    } else {
        //stride=comm.Get_size();
        //rank=comm.Get_rank();
        bin=nloop/stride;
        max=(rank==stride-1)?nloop:(rank+1)*bin;
        last_bin=nloop-(stride-1)*bin;
        counts=(int*) malloc(stride* sizeof(int));
        displ=(int*) malloc(stride* sizeof(int));
        for(unsigned i=0;i<stride-1;++i) counts[i]= static_cast<int>(bin);
        counts[stride-1]=last_bin;
        displ[0]=0;
        for(int i=1;i<stride;++i) {
            displ[i]=displ[i-1]+counts[i-1];
        }
    }
    if(do_serial or nloop/stride<10) {
        for(unsigned i=0;i<nloop;++i) {
            unsigned m=i%donorM, j=i/donorM;
            if(acceptorList[j]!=donorIndex) {
                double rho=dotProduct(rH[m],rX[j])/rX2[j];
                rhos[i]=rho;
                //if(rho>maxrho[j]) maxrho[j]=rho;
                if(rho>thres) isfirstshell[j]=1;
                pos+=((g[i]=gSW(rho,dg[i]))*rX[j]);
                gI+=g[i];
            }
        }
    }
    else { //do parallel


//@@@
//printf("bin=%d on rank %d\n",rank==stride-1?last_bin:bin,rank);

        vector<double> g_local(bin),dg_local(bin), rho_local(bin);
        //vector<double> maxrho_local(nx,0.0);
        if(bin!=last_bin and rank==stride-1) {
            g_local.resize(last_bin);
            dg_local.resize(last_bin);
            rho_local.resize(last_bin);
        }

//@@@
//printf("rank%d reached #1\n",rank);
        
        for(unsigned i=rank*bin;i<max;i+=1) {
            unsigned m=i%donorM, j=i/donorM;
            if(acceptorList[j]!=donorIndex) {
                double rho=dotProduct(rH[m],rX[j])/rX2[j];
                unsigned ii=i-rank*bin;
                rho_local[ii]=rho;
                //if(rho>maxrho_local[j]) maxrho_local[j]=rho;
                if(rho>thres) isfirstshell[j]=1;
                pos+=((g_local[ii]=gSW(rho,dg_local[ii]))*rX[j]);

//@@@
//if(i-rank*bin==459 or ii==460 or ii==461 )
//unsigned serial=fullAtomList[acceptorList[j]].serial();
//if(serial==460)
//printf("rank=%d rho=%f g_local[%d]=%f dg_local[%d]=%f\n",rank,rho,ii,g_local[ii],ii,dg_local[ii]);
//if(fabs(g_local[ii]-0.000035)<1e-7) {
//printf("ii=%d i=%d m=%d j=%d serial=%d\n",ii,i,m,j,serial);
//}
    
                gI+=g_local[ii];
            }
        }

//@@@
//printf("rank%d reached #2\n",rank);

        comm.Allgatherv(g_local,g,counts,displ);
//double* g_list=(double*) malloc(nloop*sizeof(double));

//@@@
//if(rank==0) printf("!!!!!! %f %f %f\n",g_local[6],g_local[7],g_local[8]);
//if(rank==1) {
//for(int i=0;i<stride;++i) {
//printf("displ[%d]=%d counts[%d]=%d\n",i,displ[i],i,counts[i]);
//}
//}


//MPI_Allgatherv(&g_local[0],g_local.size(),MPI_DOUBLE,&g_list[0],counts,displ,MPI_DOUBLE,comm.Get_world().Get_comm());
//g=vector<double>(g_list,g_list+nloop);

//@@@
//printf("rank%d reached #3\n",rank);

        comm.Allgatherv(dg_local,dg,counts,displ);
        comm.Allgatherv(rho_local,rhos,counts,displ);

        //MPI_Allreduce(&maxrho_local[0],&maxrho[0],nx,MPI_DOUBLE,
        //              MPI_MAX,comm.Get_world().Get_comm());
        comm.Sum(isfirstshell);

//@@@
////if(getStep()==2){ error("");}
//if(comm.Get_rank()==(getStep()+1)%4) 
//if(getStep()==0){
////for(unsigned i=0;i<g_local.size();++i)
//for(unsigned i=0;i<nx;++i)
//    //printf("g_local[%d],dg_local= %f, %f\n",i,g_local[i],dg_local[i]);
////    printf("g[%d],dg= %f %f %f, %f %f %f\n",fullAtomList[acceptorList[i]].serial(),g[i*3],g[i*3+1],g[i*3+2],dg[i*3],dg[i*3+1],dg[i*3+2]);
//    printf("g[%d],dg= %f %f %f, %f %f %f\n",acceptorList[i]-natomh,g_local[i*3],g_local[i*3+1],g_local[i*3+2],dg_local[i*3],dg_local[i*3+1],dg_local[i*3+2]);
//}


//@@@
//printf("rank%d reached #4\n",rank);

        comm.Sum(gI);

//@@@
//printf("rank%d reached #5 gI=%f\n",rank,gI);

        comm.Sum(pos);
        g_local.clear();
        dg_local.clear();
        rho_local.clear();
        //maxrho_local.clear();

//@@@
//printf("rank%d reached #6 pos=%f %f %f\n",rank,pos[0],pos[1],pos[2]);
//comm.Barrier();


    }//end of computing g and dg

//@@@
//if(comm.Get_rank()==getStep())
//for(int i=0;i<donorM*nx;++i) {
//printf("%f,%f ",g[i],dg[i]);
//if(i%donorM==2)
//printf("\n");
//}

//@@@
//if(getStep()==0 and do_serial){
////for(unsigned i=0;i<g_local.size();++i)
//for(unsigned i=0;i<nx;++i)
//    //printf("g_local[%d],dg_local= %f, %f\n",i,g_local[i],dg_local[i]);
//    printf("g[%d],dg= %f %f %f, %f %f %f\n",fullAtomList[acceptorList[i]].serial(),g[i*3],g[i*3+1],g[i*3+2],dg[i*3],dg[i*3+1],dg[i*3+2]);
//    //printf("g[%d],dg= %f %f %f, %f %f %f\n",acceptorList[i]-natomh,g_local[i*3],g_local[i*3+1],g_local[i*3+2],dg_local[i*3],dg_local[i*3+1],dg_local[i*3+2]);
//}

//@@@
//printf("rank%d reached #7 donorPos=%f %f %f\n",rank,donorPos[0],donorPos[1],donorPos[2]);

//@@@
//printf("pre-pos=%f %f %f\n",pos[0],pos[1],pos[2]);

    //construct first shell
    for(unsigned j=0;j<nx;++j) {
        //if(maxrho[j]>thres) {
        if(isfirstshell[j]) {
            firstShell.push_back(acceptorList[j]);
        }
    }

//@@@
//log.printf("firstShell: ");
//for(unsigned j=0;j<firstShell.size();++j) {
//log<< fullAtomList[firstShell[j]].serial() << " ";
//}

    gI++;
    Vector rI=pos/gI;
    pos=rI+donorPos;

//@@@
//printf("rank%d reached #8 pos=%f %f %f\n",rank,pos[0],pos[1],pos[2]);

//@@@
//printf("rank %d reaches #1 in step %ld\n",comm.Get_rank(),getStep());

    //end of pos
    //deriv
    //deriv(X_{A_j})
    if(do_serial or 10*stride>nx) {
        for(unsigned j=0;j<nx;++j) {
            unsigned global_j=acceptorList[j];
            if(global_j==donorIndex) continue;
            for (int m = 0; m < donorM; ++m) {
                unsigned jm=donorM*j+m;
                deriv[global_j]+=
                        extProduct(
                        dg[jm]*(rH[m]-2*rhos[jm]*rX[j]),
                        (rX[j]-rI)/rX2[j])
                        +
                        g[jm]*Tensor::identity();
            }
            deriv[global_j]/=gI;
        }
    }
    else {  //do in parallel
        for(unsigned j=rank*nxbin;j<nxmax;j+=1) {
            unsigned global_j=acceptorList[j];
            if(global_j==donorIndex) continue;
            for (int m = 0; m < donorM; ++m) {
                unsigned jm=donorM*j+m;
                deriv[global_j]+=
                        extProduct(
                        dg[jm]*(rH[m]-2*rhos[jm]*rX[j]),
                        (rX[j]-rI)/rX2[j])
                        +
                        g[jm]*Tensor::identity();
            }
            deriv[global_j]/=gI;
        }
        comm.Sum(deriv);
    }
    //deriv(X_\text{D})
    if(do_serial or 10*stride>nloop) {
        deriv[donorIndex]=Tensor::identity();
        for(unsigned j=0;j<nx;++j) {
//
//@@@
//if(comm.Get_rank()==1) {
//printf("step%ld j=%d deriv[donorIndex]={{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",
//        getStep(),j,deriv[donorIndex][0][0],deriv[donorIndex][0][1],
//        deriv[donorIndex][0][2],deriv[donorIndex][1][0],
//        deriv[donorIndex][1][1],deriv[donorIndex][1][2],
//        deriv[donorIndex][2][0],deriv[donorIndex][2][1],
//        deriv[donorIndex][2][2]);
//}
//
//@@@
//if(j<2) 
//printf("rank %d reaches #2 in step %ld\n",comm.Get_rank(),getStep());

            if(acceptorList[j]==donorIndex) continue;
            for(unsigned m=0;m<donorM;++m) {
                unsigned jm=donorM*j+m;
//
//@@@
//if(j==8)
//printf("step%ld m=%d rank=%d nh=%d rhos[jm]=%f rX[j]={%f,%f,%f} bondedHList[m]=%d rH[bondedHList[m]={%f,%f,%f} rH2[bondedHList[m]]=%f rX2[j]=%f rI={%f,%f,%f} dg[jm]=%f\n",
//        getStep(),m,comm.Get_rank(),nh,
//        rhos[jm],rX[j][0],rX[j][1],rX[j][2],bondedHList[m],
//        rH[bondedHList[m]][0],rH[bondedHList[m]][1],rH[bondedHList[m]][2],
//        rH2[bondedHList[m]],
//        rX2[j],rI[0],rI[1],rI[2],dg[jm]);
//
                deriv[donorIndex]+=extProduct(
                        ((2*rhos[jm]-1)*rX[j]-rH[m])/rX2[j],
                        (rX[j]-rI)*dg[jm]);
            }
        }
        deriv[donorIndex]/=gI;
    } else {
        for(unsigned jm=rank*bin;jm<max;++jm) {
            unsigned m=jm%donorM, j=jm/donorM;
            if(acceptorList[j]==donorIndex) continue;
            deriv[donorIndex]+=extProduct(
                    ((2*rhos[jm]-1)*rX[j]-rH[m])/rX2[j],
                    (rX[j]-rI)*dg[jm]);
        }
        comm.Sum(deriv[donorIndex]);
        deriv[donorIndex]=(deriv[donorIndex]+Tensor::identity())/gI;
    }
    //deriv(X_{H_m})
    if(do_serial or 10*stride>nx) {
        for (unsigned m = 0; m < donorM; ++m) {
            unsigned global_m=bondedHList[m];
            for (unsigned j = 0; j < nx; ++j) {
                if(acceptorList[j]==donorIndex) continue;
                deriv[global_m]+=
                        extProduct(
                                rX[j]/rX2[j],
                                (rX[j]-rI)*dg[donorM*j+m]
                        );
            }
            deriv[global_m]/=gI;
        }
    } else {
        for (unsigned m = 0; m < donorM; ++m) {
            unsigned global_m=bondedHList[m];
            for (unsigned j = rank; j < nx; j+=stride) {
                if(acceptorList[j]==donorIndex) continue;
                deriv[global_m]+=
                        extProduct(
                                rX[j]/rX2[j],
                                (rX[j]-rI)*dg[donorM*j+m]
                        );
            }
            comm.Sum(deriv[global_m]);
            deriv[global_m]/=gI;
        }
    }
    
//@@@
//if(comm.Get_rank()==1) {
//for(unsigned j=0;j<nx;++j) {
//printf("rX2[%d]=%f\n",fullAtomList[acceptorList[j]].serial(),rX2[j]);
//}
//}

//@@@
//if(comm.Get_rank()==1) {
//for(unsigned j=0;j<nx;j++) {
//printf("rhos[%d]=%f,%f,%f\n",fullAtomList[acceptorList[j]].serial(),
//        rhos[3*j],rhos[3*j+1],rhos[3*j+2]);
//}
//}

//@@@
//printf("step%d rank%d: nx=%d donorIndex=%d gI=%f pos={%f,%f,%f} deriv[donorIndex]={{%f %f %f},{%f %f %f},{%f %f %f}}\n",
//        getStep(),comm.Get_rank(),nx,donorIndex,gI,
//        pos[0],pos[1],pos[2],
//        deriv[donorIndex][0][0],deriv[donorIndex][0][1],
//        deriv[donorIndex][0][2],deriv[donorIndex][1][0],
//        deriv[donorIndex][1][1],deriv[donorIndex][1][2],
//        deriv[donorIndex][2][0],deriv[donorIndex][2][1],
//        deriv[donorIndex][2][2]);
//
//@@@
//for(unsigned j=0;j<nx;++j) {
//unsigned jj=acceptorList[j];
//unsigned aj=fullAtomList[jj].serial();
//if(aj==718 or aj==2548 or aj==2557)
//printf("step%d rank%d: nx=%d j=%d gI=%f pos={%f,%f,%f} deriv[m]={{%f %f %f},{%f %f %f},{%f %f %f}}\n",
//        getStep(),comm.Get_rank(),nx,acceptorList[j]-natomh+1,gI,
//        pos[0],pos[1],pos[2],
//        deriv[jj][0][0],deriv[jj][0][1],
//        deriv[jj][0][2],deriv[jj][1][0],
//        deriv[jj][1][1],deriv[jj][1][2],
//        deriv[jj][2][0],deriv[jj][2][1],
//        deriv[jj][2][2]);
//}

//@@@
//for(unsigned m=0;m<donorM;++m) {
//unsigned mm=bondedHList[m];
//printf("step%d rank%d: nx=%d m=%d gI=%f pos={%f,%f,%f} deriv[j]={{%f %f %f},{%f %f %f},{%f %f %f}}\n",
//        getStep(),comm.Get_rank(),nx,mm,gI,
//        pos[0],pos[1],pos[2],
//        deriv[mm][0][0],deriv[mm][0][1],
//        deriv[mm][0][2],deriv[mm][1][0],
//        deriv[mm][1][1],deriv[mm][1][2],
//        deriv[mm][2][0],deriv[mm][2][1],
//        deriv[mm][2][2]);
//}

    setPosition(pos);
    setMass(1.0);
    setCharge(0.0);
    setAtomsDerivatives(deriv);

    //lastpos = pos;
    lastpos = donorPos;
    rH2_=rH2;

    if(counts) free(counts);
    if(displ) free(displ);
    if(nxcounts) free(nxcounts);
    if(nxdispl) free(nxdispl);
}

}
}
