/*
    created: May 22 2017
    last modified: Oct 5 2017
    author: Chenghan Li

    Implementattion of the CEC in QC's paper. 
    Only one group of pair is supported. Atoms in PAIRS are considered to be in the same molecule(eg. suitable for the tetra histidines in M2 but not suitable for proton transfer from one histidine to another.
    PBC is used by default.
    Number of hydrogen atoms and heavy atoms should not exceed the range of unsigned.
    No box derivatives included, not suitable for NPT ensemble
*/
// TODO use std::swap to replace tedious swap in this code
// TODO multiple pairs support
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include<vector>
#include<cmath>
#include<algorithm>
#include<numeric>

using namespace std;

namespace PLMD {
namespace vatom{

class CEC : public ActionWithVirtualAtom {
//private:
    bool do_pbc;
    bool do_list;
    unsigned natomh,natomnp,natomp,natom; // I hope unsigned can hold the number of hydrogens
    vector<double> wx; // weights for AtomX; not compact
    vector<double> r0,d0; // store r0,d0 for each AtomX; not compact
    int nn; //nn=15 by default
    Vector lastpos;
    double rlist,rbond; //rlist is for finding heavy atoms
                        //rbond is for defining bonded H
                        //I assume rbond is the same for diff. ATOMX
    double rwarn, rerror; //the same assumption as above
    Vector initPos();
    unsigned getTopo(const vector<unsigned>& hindexes,
          const vector<unsigned>& oindexes,vector<long int>& topo);
    void updateList();
    vector<unsigned> heavyList, bondedHList;
    int stride;
    bool firsttime, do_updateTopo;

//@@@
//vector<AtomNumber> fullAtomList;

public:
    explicit CEC(const ActionOptions&ao);
    void calculate();
    static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(CEC,"CEC")

void CEC::registerKeywords(Keywords& keys) {
    ActionWithVirtualAtom::registerKeywords(keys);
    keys.add("atoms","ATOMH","Hydrogen atoms");
    keys.add("atoms","ATOMX","Heavy atoms");
    keys.add("atoms","D","The d parameter of the switching function for each ATOMX");
    keys.add("atoms","R0","The r_0 parameter of the switching function for each ATOMX");
    keys.add("atoms","WX","The weights for heavy atoms, typically should be negative");
    keys.add("optional","NN","The exponent used in the continuous maximum function (in the denominator)");
    keys.add("optional","RLIST","The cutoff for finding heavy atoms");
    keys.add("optional","RBOND","The distance for defining bonded hydrogen atoms");
    keys.add("hidden","RWARN","The distance for generating warning");
    keys.add("hidden","RERROR","The distance for generating error");
    keys.add("optional","STRIDE","Update lists every this step(s)");
    keys.add("atoms","PAIRS","Heavy Atoms invloved in the correction");
    keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances ");
    keys.addFlag("UPDATETOPO",false,"update topology in calculation, only useful when water cluster is larger than half box");
    keys.addFlag("NOLIST",false,"do not update lists");
    //keys.addFlag("NOUNWRAP",false,
    //             "do not unwrap coordinates when doing calculations");
}

CEC::CEC(const ActionOptions&ao):
    Action(ao),
    ActionWithVirtualAtom(ao),
    do_pbc(true), //unwrap(true),//firsttime(true),
    //r0(0.0),d0(0.0),
    nn(15),rlist(8.0),rbond(1.65),stride(2),rerror(1.2),
    firsttime(true),rwarn(1.45),do_updateTopo(false),do_list(true)
{
    //parse
    vector<AtomNumber> atomH,atomX,atomP;
    //vector<vector<AtomNumber>> atomXs; // compact version of atomX
    //vector<vector<double>> wxs, r0s, d0s; // compact version of wx, r0, d0
    //unsigned nsets = 0;
    parseAtomList("ATOMH",atomH);
    parseAtomList("ATOMX",atomX);
    //for(unsigned i = 0; ; ++i) {
    //   vector<AtomNumber> tempAtomX;
    //   if(!parseNumberedVector("ATOMX", i, tempAtomX)) {
    //      printf("i = %ld\n", i);
    //      nsets = i;
    //      break;
    //   }
    //   if(tempAtomX.empty())  error("empty ATOMX");
    //   atomXs.push_back(tempAtomX);
    //}
    //readAtomsLikeKeyword("ATOMX", nsets, atomX);
    //if(!nsets) error("Heavy atoms should be specified");
    //printf("nsets = %ld\n", nsets);
    parseAtomList("PAIRS",atomP);
    parseVector("WX",wx);
    //for(unsigned i = 0; ; ++i) {
    //   double tempwx;
    //   if(!parseNumbered("WX", i, tempwx)) {
    //      if(i != nsets) 
    //         error("Number of WX does not match number of ATOMX");
    //      break;
    //   }
    //   wxs.push_back(tempwx);
    //}
    parseVector("D",d0);
    //for(unsigned i = 0; ; ++i) {
    //   double tempd;
    //   if(!parseNumbered("D", i, tempd)) {
    //      if(i != nsets) 
    //         error("Number of D does not match number of ATOMX");
    //      break;
    //   }
    //   d0s.push_back(tempd);
    //}
    parseVector("R0",r0);
    //for(unsigned i = 0; ; ++i) {
    //   double tempr;
    //   if(!parseNumbered("D", i, tempr)) {
    //      if(i != nsets) 
    //         error("Number of R does not match number of ATOMX");
    //      break;
    //   }
    //   r0s.push_back(tempr);
    //}
    parse("NN",nn);
    parse("RLIST",rlist); parse("RBOND",rbond);
    parse("STRIDE",stride);
    parse("RWARN",rwarn); parse("RERROR",rerror);
    parseFlag("UPDATETOPO",do_updateTopo);
    bool nopbc;//, nounwrap;
    parseFlag("NOPBC",nopbc);
    bool no_list;
    parseFlag("NOLIST",no_list);
    do_list = !no_list;
    if(!do_list) do_updateTopo = false;
    //parseFlag("NOUNWRAP",nounwrap);
    do_pbc=!nopbc; //unwrap=!nounwrap;
    //lastpos=Vector(0.0,0.0,0.0);
    checkRead();
    //error
    natomh= static_cast<unsigned int>(atomH.size());
    if(natomh==0) error("Hydrogen atoms should be specified");
    //unpress wxs, r0s, d0s;
    //for(unsigned i = 0; i < nsets; ++i) {
    //   for(unsigned ii = 0; ii < wxs[i].size(); ++ii) {
    //      wx.push_back(wxs[i][ii]); r0.push_back(r0s[i][ii]);
    //      d0.push_back(d0s[i][ii]);
    //   }
    //}
    if(atomX.empty()) error("Heavy atoms should be specified");
    //if(r0<=0.0||d0<=0.0) error("R0 and D should be explicitly specified and positive");
    if(nn<2) error("Please try some larger NN value for good approximation of the maximum");
    if(!atomP.empty()) { //check whether pairs is included in atomX
        for(vector<AtomNumber>::iterator itp=atomP.begin();\
        itp!=atomP.end();++itp) {
            if(find(atomX.begin(),atomX.end(),(*itp))==atomX.end()) {
                error("Every PAIRS atom should be included in ATOMX");
            }
        }
    }
    if(wx.size()!=atomX.size()) error("Number of WX should be the same as number of ATOMX");
    if(r0.size()!=atomX.size()) error("Number of R0 should be the same as number of ATOMX");
    if(d0.size()!=atomX.size()) error("Number of D should be the same as number of ATOMX");
    double total=natomh-1;
    for(unsigned i=0;i<wx.size();++i) {
        total+=wx[i];
    }
    if(fabs(total)>1e-6) error("Number of excess proton should be exactly one");
    if(do_list && rlist<=0) error("RLIST should be positive");
    if(do_list && rbond<=0) error("RBOND should be positive");
    if(do_list && stride<=0) error("STRIDE should be positive");
    //past processing
    natomnp= static_cast<unsigned int>(atomX.size() - atomP.size());
    for(unsigned i=0;i<atomP.size();i++) { //iterate atomX, trying to move all the pair atoms to the bottom of atomX
        vector<double>::iterator itw=wx.begin(), itr=r0.begin(), itd=d0.begin();
        for(vector<AtomNumber>::iterator itx=atomX.begin();itx!=atomX.end()-i;++itx) {
            if(find(atomP.begin(),atomP.end(),(*itx))!=atomP.end()) { //if atomj is a pair atom, then we consider it to be ""heavier" than atomj+1
                AtomNumber an_tmp=*itx; double wx_tmp=*itw, d_tmp=*itd, r_tmp=*itr;
                (*itx)=*(atomX.end()-i-1); (*itw)=*(wx.end()-i-1);
                (*itr)=*(r0.end()-i-1); (*itd)=*(d0.end()-i-1);
                *(atomX.end()-i-1)=an_tmp; *(wx.end()-i-1)=wx_tmp;
                *(r0.end()-i-1)=r_tmp; *(d0.end()-i-1)=d_tmp;
                break;
            }
            ++itw; ++itr; ++itd;
        }
    }
    atomH.insert(atomH.end(),atomX.begin(),atomX.end()); //concatenate 
    vector<AtomNumber>::iterator atomH_end=atomH.begin()+natomh;
    vector<AtomNumber>::iterator atomNP_end=atomH_end+natomnp;
    natomnp+=natomh; //now natomnp is the number of non pair atoms including H
    natomp=atomP.size();
    natom=natomp+natomnp;
    for(unsigned i = 0; i < natomh; ++i) bondedHList.push_back(i);
    for(unsigned i = natomh; i < natom; ++i) heavyList.push_back(i);
    //write log
    log.printf("  %d hydrogen atoms and %d heavy atoms involved\n",\
    natomh,atomX.size());
    log<<"  serial number of hydrogen atoms:\n  ";
    for(vector<AtomNumber>::iterator it=atomH.begin();it!=atomH_end;++it) {
        log.printf("%d ",(*it).serial());
    }
    log<<"\n  serial number of heavy atoms:\n  ";
    for(vector<AtomNumber>::iterator it=atomH_end;it!=atomH.end();++it) {
        log.printf("%d ",(*it).serial());
    }
    if(!atomP.empty()) {
        log<<"\n  serial number of pair atoms:\n  ";
        for(vector<AtomNumber>::iterator it=atomNP_end;it!=atomH.end();\
        ++it) {
            log.printf("%d ",(*it).serial());
        }
    }
    else
        log<<"\n  no pair atoms specified, work in single proton mode";
    log <<"\n  weights of heavy atoms:\n  ";
    for(vector<double>::iterator it=wx.begin();it!=wx.end();++it)    {
        log<<*it<<" ";
    }
    log <<"\n  r parameters of heavy atoms:\n  ";
    for(vector<double>::iterator it=r0.begin();it!=r0.end();++it)    {
        log<<*it<<" ";
    }
    log <<"\n  d parameters of heavy atoms:\n  ";
    for(vector<double>::iterator it=d0.begin();it!=d0.end();++it)    {
        log<<*it<<" ";
    }
    if(do_list) {
      log<<"\n  using "<< rlist << " to find heavy atoms";
      log<<"\n  hydrogen within "<< rbond << " of heavy atoms will be regarded as bonded";
      log<<"\n  lists will be updated every " << stride << " step(s)";
    }
    
    log<<"\n";
    if(do_pbc){
        log<<"  PBC will be considered\n";
    } else {
        log<<"  PBC will be ignored\n";
    }
    if(do_updateTopo){
        log<<"  Topology will be updated every step\n";
    } else {
        log<<"  Do not consider topology in calculation\n";
    }
    log<<"  charge and mass of the CEC are set to be 0 and 1 in this implementation\n";
    if(do_list) {
      log<<"  dist between nonbonded heavy atom and hydrogen lower than "
         << rwarn << " will give warnings\n";
      log<<"  dist between nonbonded heavy atom and hydrogen lower than "
         << rerror << " will stop the program\n";
    }
    //log<<"  WARNING: PBC was not tested extensively, be careful...\n";
    //initialize indexes, allocate positions, masses, ...
    requestAtoms(atomH);

//printf("atomH.size() = %d\n",atomH.size()); // @@@

//@@@
//fullAtomList = atomH;

}

//getTopo() for given H's and heavy's
// return the index of the excess hydrogen atom in hindexes
// the *indexes contain the indexes in fullAtomList
unsigned CEC::getTopo(const vector<unsigned>& hindexes,
        const vector<unsigned>& oindexes, vector<long int>& topo) {
    topo=vector<long int>(hindexes.size(),-1);
    vector<int> counts;
    vector<unsigned> index_pair; // index_pair[i] = index in oindexes
    for(unsigned i=0;i<oindexes.size();++i) {
        unsigned io = oindexes[i];
        if(io<natomnp) {
            counts.push_back(int(0.5-wx[io-natomh]));
        } else {
            index_pair.push_back(i);
            counts.push_back(int(0.5-wx[io-natomh]*natomp));
        }
    }
    if(!index_pair.empty() && index_pair.size()!=natomp) 
        plumed_merror("ERROR: pair atom missing");

    //loop for O:
    //  find the closest unassigned H until its count is zero !!!! PAIR
    for(unsigned j=0;j<oindexes.size();++j) { 
        while(counts[j]) {
            unsigned closest_h=0;
            double min_d2=rlist*rlist;
            for(unsigned i=0;i<hindexes.size();++i) {
                if(topo[i]!=-1) continue;
                double d2;
                //distance between two atoms
                if(oindexes[j]<natomnp) {
                   Vector distance = do_pbc?
                       pbcDistance(getPosition(hindexes[i]),getPosition(oindexes[j])):
                       delta(getPosition(hindexes[i]),getPosition(oindexes[j]));
                   d2=distance.modulo2();
                   if(d2<min_d2) {
                       min_d2=d2; closest_h=i;
                   }
                }
                //distance between one atom and one molecule
                else {
                   for(unsigned k=0;k<index_pair.size();++k) {
                      Vector distance=do_pbc?
                       pbcDistance(getPosition(hindexes[i]),
                             getPosition(oindexes[index_pair[k]])):
                       delta(getPosition(hindexes[i]),
                             getPosition(oindexes[index_pair[k]]));
                      d2=distance.modulo2();
                      if(d2<min_d2) {
                         min_d2=d2; closest_h=i;
                      }
                   }
                }
            }
            topo[closest_h]=j;
            //modify counts
            unsigned jo=oindexes[j];
            if(jo<natomnp) {
                counts[j]--;
            } else {
                for(unsigned k=0;k<natomp;++k) {
                    counts[index_pair[k]]--;
                }
            }
        }
    }

    // find the -1 in topo and assign the closest O to that
    unsigned rc=0;
    unsigned count_rc=0;
    plumed_assert(topo.size()==hindexes.size());
    for(unsigned i=0;i<topo.size();++i) {
        if(topo[i]==-1) {
            rc=i; count_rc++;
        }
    }
    plumed_assert(count_rc==1);
    double min_d2 = rlist * rlist;
    unsigned closet_o=0;
    for(unsigned j=0;j<oindexes.size();++j) {
        Vector distance=do_pbc?
                 pbcDistance(getPosition(hindexes[rc]),getPosition(oindexes[j])):
                 delta(getPosition(hindexes[rc]),getPosition(oindexes[j]));
        double d2 = distance.modulo2();
        if(d2<min_d2) {
            min_d2 = d2; closet_o = j;
        }
    }
    topo[rc] = closet_o;
    return rc;
}

Vector CEC::initPos() {
    if(!do_list)
      return Vector(0.0,0.0,0.0);
    vector<unsigned> hindexes(natomh),oindexes(natom-natomh);
    vector<long int> topo;
    for(unsigned i=0;i<natomh;++i) {
        hindexes[i]=i;
    }
    for(unsigned i=natomh;i<natom;++i) {
        oindexes[i-natomh]=i;
    }
    unsigned rc = getTopo(hindexes,oindexes,topo);
    return getPosition(rc);
}

//void CEC::updateList(vector<unsigned>& heavyList, vector<unsigned>& bondedHList) {
void CEC::updateList() {
    //unsigned natom = getNumberOfAtoms();
    double rlist2 =  rlist*rlist;
    double rbond2 = rbond*rbond;
    double total=0; // total weights of heavy atoms

    heavyList.clear();
    for(unsigned j=natomh; j<natom; ++j) {
        Vector distance=do_pbc?pbcDistance(lastpos,getPosition(j)):
            delta(lastpos,getPosition(j));
        double d2=distance.modulo2();
        if(d2<=rlist2) {
            heavyList.push_back(j);
            total += wx[j-natomh];
        }
    }
    unsigned np_heavyList=0; // # of pair in heavyList
    for(unsigned i=0; i<heavyList.size(); ++i) {
        if(heavyList[i]>=natomnp) np_heavyList++;
    }
    if(np_heavyList!=natom-natomnp) {
        log<<"WARNING: Pair atoms are separated by cutoff at step "<<getStep()
            <<"\n";
        // if pair molecules are separated by cutoff
        // then erase all the pairs from heavyList
        for(vector<unsigned>::iterator it=heavyList.begin();
                it!=heavyList.end();) {
            if(*it>=natomnp) {
                total -= wx[(*it)-natomh];
                it=heavyList.erase(it);
            }
            else it++;
        }
    }
    unsigned nh = int(-total+0.5); //expected number of H
    if(fabs(nh+total)>1e-6) plumed_merror("ERROR: total is not an integer");
    nh++;

    vector<int> isbonded(natomh,0);
    vector<double> mind2(natomh,rbond2); //min_{all heavy} d2(heavy,H) ???
    for(unsigned jj=0; jj<heavyList.size(); ++jj) {
        for(unsigned i=0;i<natomh;++i) {
            Vector distance=do_pbc?
                pbcDistance(getPosition(heavyList[jj]),getPosition(i)):
                delta(getPosition(heavyList[jj]),getPosition(i));
            double d2=distance.modulo2();
            if(d2<=rbond2) {
                isbonded[i]++;
                if(d2<mind2[i]) mind2[i]=d2;
            }
        }
    }

    //find the closest nh H's in natomh H's
    unsigned count=0, index=0;
    vector<double> bondedHd2;
    bondedHList.clear();
    for(unsigned i=0; i<natomh; ++i) {
        if(isbonded[i]==0) continue;
        bondedHList.push_back(i);
        bondedHd2.push_back(mind2[i]);
        count++;
        if(count==nh) {
            index = i;
            //sort the bondedHList according to bondedHd2 
            for(unsigned j=0; j<nh; j++) {
                for(unsigned k=0; k<nh-j-1; k++) {
                    if(bondedHd2[k]>bondedHd2[k+1]) {
                        unsigned ind = bondedHList[k];
                        double s = bondedHd2[k];
                        bondedHList[k] = bondedHList[k+1];
                        bondedHd2[k] = bondedHd2[k+1];
                        bondedHList[k+1] = ind;
                        bondedHd2[k+1] = s;
                    }
                }
            }
            break;
        }
    }
    if(count!=nh) {
        std::string counts;
        Tools::convert(count-nh+1, counts);
        std::string msg="ERROR: the number of excess charge is not 1 but "+
           counts;
        //plumed_merror("ERROR: the number of excess charge is not 1");
        plumed_merror(msg);
    }
    double minbondedHd2_back=rbond2+1;
    unsigned argminbondedHd2_back=0;
    for(unsigned i=index+1; i<natomh; ++i) {
        if(mind2[i]<bondedHd2.back()) {
            double &back = bondedHd2.back();
            if(back<minbondedHd2_back) { 
                minbondedHd2_back=back;
                argminbondedHd2_back = bondedHList.back();
            }
            bondedHList.back() = i;
            back = mind2[i];

            // O(N) reorder in bondedHList and bondedHd2 
            for(unsigned j=nh-1; j>0; --j) {
                if(bondedHd2[j]<bondedHd2[j-1]) {
                    unsigned ind = bondedHList[j];
                    double s = bondedHd2[j];
                    bondedHList[j] = bondedHList[j-1];
                    bondedHd2[j] = bondedHd2[j-1];
                    bondedHList[j-1] = ind;
                    bondedHd2[j-1] = s;
                }
                else break;
            }
        }
        else { //the first nh H's are the best at this moment
           if(mind2[i]<minbondedHd2_back) {
              minbondedHd2_back=mind2[i];
              argminbondedHd2_back=i;
           }
        }
    }
    if(minbondedHd2_back <= rerror*rerror) {
        std::string errstr = "ERROR: partial decomposition happens at step ";
        std::string tmpstr;
        Tools::convert(getStep(), tmpstr);
        errstr = errstr + tmpstr + " for atom ";
        Tools::convert(argminbondedHd2_back, tmpstr);
        errstr = errstr + tmpstr + "\n";
        plumed_merror(errstr);
    }
    if(minbondedHd2_back <= rwarn*rwarn) {
        log<<"WARNING: partial decomposition may happen at step "<< getStep()
            << " for atom "<<argminbondedHd2_back<< "\n";
        log<<"         " << sqrt(minbondedHd2_back) <<" <= "<< rwarn<<"\n";
    }
}

void CEC::calculate() {

//printf("0: %ld\n",getStep()); //@@@

    Vector pos;
    Vector vdij; //\vec{r}^{H_i}-\vec{r}^{X_j}
    double dij,fij,dfij; //dij=norm of vdij, fij=f_{sw}(dij),dfij=f'_{sw}(dij)
    double tmp;
    vector<Tensor> deriv(natom);
    if(firsttime) {
        lastpos=initPos();
        log<<"Info: Initial guess for CEC is "
           <<lastpos[0] <<' '<< lastpos[1] <<' '<< lastpos[2] <<'\n';
        if(do_list) updateList();
        firsttime = false;
    }
    else if(getStep()%stride==0) {
        if(do_list) updateList();
    }
    vector<long int> topo;

//printf("1: %ld\n",getStep()); //@@@


    unsigned nh = bondedHList.size(), nx = heavyList.size();
    vector<Vector> relpos(nh+nx);
    unsigned np=0;
    for(unsigned j=0;j<nx;j++) {
        if(heavyList[j]>=natomnp)  np++;
    }
    if(np!=0&&np!=natomp) plumed_merror("ERROR: np!=natomp");
    unsigned nnp = nh+nx-np, n = nh+nx;

//printf("natomh = %d natomnp = %d natom = %d nh = %d nx = %d n = %d natomp = %d\n",natomh, natomnp, natom, nh, nx, n, natomp); //@@@

    //pbc all the heavy by lastpos
    for(unsigned j=0;j<nx;++j) {
        if(do_list) {
          Vector dist=do_pbc?
              pbcDistance(lastpos,getPosition(heavyList[j])):
              delta(lastpos,getPosition(heavyList[j]));
          relpos[j+nh] = lastpos + dist;
        } else {
          // if no nb list is used, just use the original positions
          relpos[j+nh] = getPosition(heavyList[j]);
        }
    }
    //pbc all the H by lastpos or by heavy
    if(do_updateTopo) { 
        getTopo(bondedHList,heavyList,topo);
        for(unsigned i=0;i<nh;++i) {
            unsigned iref = nh + topo[i];
            Vector dist=do_pbc?
                pbcDistance(relpos[iref],getPosition(bondedHList[i])):
                delta(relpos[iref],getPosition(bondedHList[i]));
            relpos[i]=relpos[iref]+dist;
        }
    } else {
        for(unsigned i=0;i<nh;++i) {
            if(do_list) {
               Vector dist=do_pbc?
                   pbcDistance(lastpos,getPosition(bondedHList[i])):
                   delta(lastpos,getPosition(bondedHList[i]));
               relpos[i] = lastpos + dist;
            } else {
               // if no nb list is used, just use the original positions
               relpos[i] = getPosition(bondedHList[i]);
            }
        }
    }

    vector<double> wx_reduce, r0_reduce, d0_reduce;
    for(unsigned j=0;j<nx;++j) {
        wx_reduce.push_back(wx[heavyList[j]-natomh]);
        r0_reduce.push_back(r0[heavyList[j]-natomh]);
        d0_reduce.push_back(d0[heavyList[j]-natomh]);
    }

    vector<double> fnnm1(np*nh); //fnnm1(i*np+j) stores f_sw(d_HiXj)^(nn-1) fnnm1 means f^(nn minus 1)
    vector<double> fnn(np*nh); //fnn(i*np+j) stores f_sw(d_HiXj)^nn fnn means f^nn
    vector<double> dfsw(np*nh); //dfsw(i*np+j) stores df_sw(d_HiXj)/d_HiXj
    vector<double> tnn(np,0.0); //tnn[j] stores \sum_k f_sw(d_HkXj)^nn
    vector<double> tnnp1(np,0.0); //tnnp1[j] stores \sum_k f_sw(d_HkXj)^nn+1
    vector<double> mi(np); //mi[j] stores m_j
    //compute pos and derivatives
    for(unsigned j=nh;j<n;++j) { //loop for heavy atoms
        pos+=(wx_reduce[j-nh])*relpos[j];

//if(fullAtomList[heavyList[j-nh]].serial()==868) { // @@@
//   if(fabs(wx_reduce[j-nh])>1e-6)
//   log.printf("deriv += wx[%d] (%f) * identity\n",j-nh,wx_reduce[j-nh]);
//}

        deriv[heavyList[j-nh]]=wx_reduce[j-nh]*Tensor::identity();
    }
    for(unsigned i=0;i<nh;++i) {  //loop for hydrogens
        pos+=relpos[i];
        deriv[bondedHList[i]]=Tensor::identity();
        for(unsigned j=nh;j<n;++j) { //loop for heavy atoms
            //vdij=do_pbc?pbcDistance(getPosition(j),getPosition(i)):
            //     delta(getPosition(j),getPosition(i));
            vdij=do_pbc?pbcDistance(relpos[j],relpos[i]):
               delta(relpos[j],relpos[i]);
            dij=vdij.modulo();
            //tmp=exp((dij-r0)/d0);fij=1.0/(1.0+tmp);dfij=-pow(fij,2)*tmp/d0;
            if((tmp = (dij-r0_reduce[j-nh])/d0_reduce[j-nh]) <= 600) {
               tmp=exp(tmp);
               fij=1.0/(1.0+tmp);dfij=-pow(fij,2)*tmp/d0_reduce[j-nh];
            } else {
               fij = 0.0; dfij = 0.0;
            }
            pos-=fij*vdij;
            deriv[bondedHList[i]]-=fij*Tensor::identity();
            deriv[bondedHList[i]]-=dfij*extProduct(vdij,vdij)/dij;

//if(fullAtomList[heavyList[j-nh]].serial()==868) { // @@@
//   if(fabs(fij)>1e-6)
//   log.printf("deriv += f_%u_%u (%f) * identity\n",fullAtomList[bondedHList[i]].serial(), fullAtomList[heavyList[j-nh]].serial(), fij);
//   if(fabs(dfij/dij)>1e-6)
//   log.printf("deriv += df_%u_%u (%f) / d_%u_%u (%f) * vdij x vdij\n",fullAtomList[bondedHList[i]].serial(), fullAtomList[heavyList[j-nh]].serial(), dfij,fullAtomList[bondedHList[i]].serial(), fullAtomList[heavyList[j-nh]].serial(), dij);
//}

            deriv[heavyList[j-nh]]+=fij*Tensor::identity();
            deriv[heavyList[j-nh]]+=dfij*extProduct(vdij,vdij)/dij;
            //for pair correction
            if(j>=nnp) { //if atom_j is a pair atom
                fnnm1[i*np+j-nnp]=pow(fij,nn-1);
                fnn[i*np+j-nnp]=(fnnm1[i*np+j-nnp])*fij;
                tnn[j-nnp]+=fnn[i*np+j-nnp];
                tnnp1[j-nnp]+=fnn[i*np+j-nnp]*fij;
                dfsw[i*np+j-nnp]=dfij/dij;
            }
        }
    }
    for(unsigned j=0;j<np;++j) { //loop in pairs
        if(tnn[j]<1e-300) 
            mi[j]=0.0;
        else
            mi[j]=tnnp1[j]/tnn[j];
    }
    for(unsigned i=nnp;i<n;++i) { //correction term for pos
        for(unsigned j=nnp;j<n;++j) {
            vdij=relpos[j]-relpos[i];
            pos+=1.0/np*mi[i-nnp]*vdij;
        }
    }

//printf("-2: %ld\n",getStep()); //@@@

    double mi_sum=accumulate(mi.begin(),mi.end(),0.0);
    for(unsigned i=0;i<np;++i) { //loop for pair atom
        deriv[heavyList[i+nnp-nh]]+=(mi_sum/np-mi[i])*Tensor::identity();
        Vector vtmp(0.0,0.0,0.0);
        for(unsigned q=0;q<nh;++q) { //loop for hydrogen
            unsigned ind=q*np+i;
            vtmp+=dfsw[ind]*(fnn[ind]*(nn+1)-fnnm1[ind]*nn*mi[i])*
                (relpos[i+nnp]-relpos[q]);
        }
        for(unsigned l=0;l<np;++l) { //loop for pair
            if(tnn[i]>1e-300) {
                Vector dist=relpos[l+nnp]-relpos[i+nnp];
                deriv[heavyList[i+nnp-nh]]+=
                        extProduct(vtmp,dist)/tnn[i]/np;
            }
        }
    }
    for(unsigned i=0;i<nh;i++) { //loop for hydrogens
        for(unsigned l=0;l<np;++l) { //loop for pair atoms
            for(unsigned q=0;q<np;++q) { //loop for pair atoms
                unsigned ind=i*np+q;
                if(tnn[q]>1e-300) {
                    Vector dist1 = relpos[i] -relpos[q+nnp];
                    Vector dist2 = relpos[l+nnp] -relpos[q+nnp];
                    deriv[bondedHList[i]]+=1.0/np/tnn[q]*dfsw[ind]*
                              ((nn+1)*fnn[ind]-nn*mi[q]*fnnm1[ind])*extProduct(dist1,dist2);
                }
            }
        }
    }
    lastpos=pos;

    setPosition(pos);
    setMass(1.0);
    setCharge(0.0);
    setAtomsDerivatives(deriv);

//printf("-1: %ld\n",getStep()); //@@@

//printf("pos = %f %f %f\n", pos[0], pos[1], pos[2]); // @@@

}

}
}
