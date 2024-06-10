/*
    created: Dec 11 2017
    author: Chenghan Li

    description: approximation to MS-RMD CEC w/o running MS-RMD
*/
#include "RMDCEC.h"

namespace PLMD {
namespace vatom {

PLUMED_REGISTER_ACTION(RMDCEC,"RMDCEC")

void RMDCEC::registerKeywords(Keywords& keys) {
    ActionWithVirtualAtom::registerKeywords(keys);
    keys.add("atoms","ATOMH","hydrogen atoms");
    keys.add("atoms","ATOMX","heavy atoms");
    keys.add("compulsory","WX","max number of protons each heavy atom can accept");
    keys.add("compulsory","MOLID","molecular id for heavy atoms");
    //keys.add("optional","KAPPA","The exponential parameter for the switching function");
    //keys.add("optional","RLOWER","The inner cutoff for the switching function");
    //keys.add("optional","RUPPER","The outter cutoff for the switching function");
    keys.add("compulsory","SWITCH","switching function to be used");
    keys.add("optional","SWITCHROO","switching function to be used");
    //keys.add("atoms","PAIRS","Heavy Atoms invloved in the correction");
    keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances ");
    keys.add("optional","ATOMH_CHARGES","proton's charge when bonding to heavy atoms");
    keys.add("optional","ATOMX_CHARGES","heavy atoms' charges when bonded");
    keys.add("optional","MAXNH","max number of protons one molecule has when no excess charge present");
    //keys.add("optional","A2","coefficient for gmat^2");
    //keys.add("optional","A3","coefficient for gmat^3");
    keys.addFlag("USEMASK",false,"use the mask to accelerate calculation");
    keys.addFlag("FOURTHSHELL",false,"include the fourth solvation shell");
    keys.add("optional","MOLINFO","a file containing WX MOLID ATOMH_CHARGES ATOMX_CHARGES MAXNH");
    for(unsigned i = 0; i < MAXNUMSPECIAL; ++i) {
      std::string str; Tools::convert(i, str);
      keys.add("optional","SWITCH"+str,"switching function to be used");
      keys.add("optional","SWITCHROO"+str,"switching function to be used");
    }
    keys.add("optional", "RATIO", "ratio of max water O-H bond length over dmax");
    keys.add("optional", "INIT_TOPO", "initial topology of the system");
    keys.addFlag("PRINT_REACTION", false, "if print reaction in plumed out file");
    keys.addFlag("PRINT_DELTA", false, "if print three delta values");
    keys.addFlag("PRINT_SHELLS", false, "if print out CEC contributions from each shell");
    keys.addFlag("PRINT_O_H", false, "if print out CEC contributions from oxygens and hydrogens");
    keys.addFlag("PRINT_ACCEPTORS", false, "if print out three most probable acceptors");
    keys.addFlag("PRINT_NSTATES", false, "print out N number of states with leading c_i^2");


keys.add("optional","H","for debug");//@@@
keys.add("optional","O","for debug");//@@@

#ifdef __FAKE_VATOM__
    keys.add("compulsory", "MOLIDS2WATCH", 
          "Molecule indexes in this action to be watched");
#endif
}

unsigned RMDCEC::getSerial(unsigned i) const {
   return fullAtomList[oindexes[i]].serial();
}

#ifdef __DOFITTING__
RMDCEC::~RMDCEC() {
   fs_out.close();  // for fitting @@@
   //fs_delta.close();  // for fitting @@@
}
#endif

RMDCEC::RMDCEC(const ActionOptions&ao):
    Action(ao),
    ActionWithVirtualAtom(ao),
    do_pbc(true), 
    //kappa(4.39419), r1(1.5), r2(1.8),
    use_mask(false), do_fourth(false), firsttime(true),
    ratio(0.75),
    print_reaction(false),
    print_nstates(0)
#ifdef __DOFITTING__
    //, fs_out("/project/gavoth/chhli/fitrmd/0_glu_in_water/rmdcec/fit_paras/ci_plumed.dat") // for fitting @@@
    , fs_out("ci_plumed.dat") // for fitting @@@
    //, fs_delta("delta.csv") // for fitting @@@
#endif
{
    //parse
h_serial=-1; o_serial=-1;
parse("H",h_serial);//@@@
parse("O",o_serial);//@@@
    
    vector<AtomNumber> atomH,atomX,atomP;
    parseAtomList("ATOMH",atomH);
    parseAtomList("ATOMX",atomX);
    //parseAtomList("PAIRS",atomP);
    string fn_molinfo;
    parse("MOLINFO", fn_molinfo);
    if(fn_molinfo.length() != 0) {
       //TODO use IFile class to subrogate this ad hoc solution
       log.printf("  input molinfo file:\n  %s\n", fn_molinfo.c_str());
       ifstream fs_molinfo(fn_molinfo);
       if(!fs_molinfo.is_open())
          error("cannot open file "+fn_molinfo);
       string line; unsigned line_count = 0;
       while(getline(fs_molinfo, line)) {
          line_count++;
          int wx_, max_proton_;
          unsigned  molid_;
          double proton_charge, heavy_charge;
          stringstream ss(line);
          ss >> wx_ >> molid_ >> proton_charge >> heavy_charge >> max_proton_;
          wx.push_back(wx_); molid.push_back(molid_); 
          proton_charges.push_back(proton_charge); 
          heavy_charges.push_back(heavy_charge);
          max_proton.push_back(max_proton_);
       }
       fs_molinfo.close();
    } else {
      parseVector("WX",wx);
      parseVector("MOLID",molid);
      parseVector("ATOMH_CHARGES",proton_charges);
      parseVector("ATOMX_CHARGES",heavy_charges);
      parseVector("MAXNH",max_proton);
    }
    //double A2 = 1.0, A3 = 1.0;
    //parse("A2", A2); parse("A3", A3);
    //A = array<double, 3> ({1.0, A2, A3});
    //parse("KAPPA", kappa);
    //parse("RLOWER", r1); parse("RUPPER", r2);
    //sf_indexes.resize(atomX.size(), atomX.size());
    /*--- parse switch for delta ---*/
    std::string sw, errors; parse("SWITCH", sw);
    if(sw.length()) {
       SwitchingFunction sf_;
       //sf.set(sw, errors);
       sf_.set(sw, errors); sf_delta_library.push_back(sf_);
       if(errors.length()) error("problem reading SWITCH keyword : " + errors);
       //for(unsigned i = 0; i < atomX.size(); ++i) 
       //  for(unsigned j = 0; j < atomX.size(); ++j)
       //     sf_indexes[i][j] = 0;
    }
    else error("SWITCH is required");
    for(unsigned i = 0; i < atomX.size() && i < MAXNUMSPECIAL; ++i) {
       std::string sw, errors, num;
       Tools::convert(i, num);
       parse("SWITCH"+num, sw);
       if(sw.length()) {
          SwitchingFunction sf_, isf_;
          sf_.set(sw, errors); sf_delta_library.push_back(sf_);
          if(errors.length()) 
             error("problem reading SWITCH"+num+" keyword : "+errors);
          // inverse of sf_:
          isf_ = sf_; isf_.set_d0(- sf_.get_d0());
          sf_delta_library.push_back(isf_);
          //for(unsigned j = 0; j < atomX.size(); ++j) {
          //   // In this way, the sf between special atoms are overwritten 
          //   // again and again. It is fine as long as those sf are not useful
          //   // sf[i][j] = sf_; sf[j][i] = isf_;
          //   sf_indexes[i][j] = sf_delta_library.size() - 2; 
          //   sf_indexes[j][i] = sf_delta_library.size() - 1;
          //}
          map2sf_delta_index[i + 1] = sf_delta_library.size() - 2;
          map2sf_delta_index[-static_cast<long>(i + 1)] = sf_delta_library.size() - 1;
       }
    }
    /*--- parse switch for roo ---*/
    sw.clear();
    parse("SWITCHROO", sw);
    if(sw.length()) {
       SwitchingFunction sf_;
       //sf.set(sw, errors);
       sf_.set(sw, errors); sf_roo_library.push_back(sf_);
       if(errors.length()) error("problem reading SWITCHROO keyword : " + errors);
#ifdef __DEBUG__
       log.printf("just to make sure I did not get here for r_OO\n");
#endif
       //for(unsigned i = 0; i < atomX.size(); ++i) 
       //  for(unsigned j = 0; j < atomX.size(); ++j)
       //     sf_indexes[i][j] = 0;
    }
    else {
       SwitchingFunction sf_;
#ifdef __DEBUG__
       log.printf("just to make sure I got here for r_OO\n");
#endif
       sf_.set("CONSTANT", errors); sf_roo_library.push_back(sf_);
    }
    for(unsigned i = 0; i < atomX.size() && i < MAXNUMSPECIAL; ++i) {
       std::string sw, errors, num;
       Tools::convert(i, num);
       parse("SWITCHROO"+num, sw);
       if(sw.length()) {
          SwitchingFunction sf_, isf_;
          sf_.set(sw, errors); sf_roo_library.push_back(sf_);
          if(errors.length()) 
             error("problem reading SWITCHROO"+num+" keyword : "+errors);
          //// inverse of sf_:
          //isf_ = sf_; isf_.set_d0(- sf_.get_d0());
          //sf_roo_library.push_back(isf_);
          map2sf_roo_index[i + 1] = sf_roo_library.size() - 1;
          //map2sf_roo_index[-static_cast<long>(i + 1)] = sf_roo_library.size() - 1;
       }
    }
    parseFlag("PRINT_REACTION", print_reaction);
    parseFlag("PRINT_DELTA", print_delta);
    parseFlag("PRINT_SHELLS", print_shells);
    parseFlag("PRINT_O_H", print_o_h);
    parseFlag("PRINT_ACCEPTORS", print_acceptors);
    parse("PRINT_NSTATES", print_nstates);
    bool nopbc;//, nounwrap;
    parseFlag("NOPBC",nopbc);
    do_pbc=!nopbc; //unwrap=!nounwrap;
    parseFlag("USEMASK",use_mask);
    parseFlag("FOURTHSHELL",do_fourth);
    //error
    natomh= static_cast<unsigned int>(atomH.size());
    if(natomh==0) error("Hydrogen atoms should be specified");
    if(atomX.empty()) error("Heavy atoms should be specified");
    //if(!atomP.empty()) { //check whether pairs is included in atomX
    //    for(vector<AtomNumber>::iterator itp=atomP.begin();\
    //    itp!=atomP.end();++itp) {
    //        if(find(atomX.begin(),atomX.end(),(*itp))==atomX.end()) {
    //            error("Every PAIRS atom should be included in ATOMX");
    //        }
    //    }
    //}
    if(wx.size()!=atomX.size()) {
       log.printf("WX : %d, ATOMX : %d\n", wx.size(), atomX.size());
       error("Number of WX should be the same as number of ATOMX");
    }
    if(molid.size()!=atomX.size()) {
       log.printf("MOLID : %d, ATOMX : %d\n", molid.size(), atomX.size());
       error("Number of MOLID should be the same as number of ATOMX");
    }
    if(std::find(molid.begin(), molid.end(), 0) == molid.end())
       error("MOLID should starts from 0");
    int molid_max = *std::max_element(molid.begin(), molid.end());
    for(int i = 1; i < molid_max ; ++i) {
      if(std::find(molid.begin(), molid.end(), i) == molid.end())
         error("MOLID should be consecutive numbers");
    }
    if(proton_charges.empty()) 
       proton_charges = vector<double>(atomX.size(), 0.44);
    if(heavy_charges.empty())
       heavy_charges = vector<double>(atomX.size(), -0.32);
    if(proton_charges.size() != atomX.size())  {
       log.printf("ATOMH_CHARGES : %d, ATOMX : %d\n", 
             proton_charges.size(), atomX.size());
       error("Number of ATOMH_CHARGES should be the same as number of ATOMX");
    }
    if(heavy_charges.size() != atomX.size()) {
       log.printf("ATOMX_CHARGES : %d, ATOMX : %d\n", 
             heavy_charges.size(), atomX.size());
       error("Number of ATOMX_CHARGES should be the same as number of ATOMH");
    }
    //past processing
    bondedHList = vector<vector<unsigned>>(atomX.size());
    heavyList = vector<vector<unsigned>>(molid_max + 1);
    for(auto it = molid.begin(); it != molid.end(); ++it) 
       heavyList[*it].push_back(it - molid.begin());
    if(max_proton.empty()) {
       // if user does not refer max_proton, 
       // use \sum_{j, j\in molid[i]} wx[j] 
       for(auto it = heavyList.begin(); it != heavyList.end(); ++it) {
          int s = 0;
          for(unsigned j : *it) s += wx[j];
          max_proton.push_back(s-1);
       }
       log.printf("WARNING: MAXNH not found, using default\n");
    }
    else {
       if(max_proton.size() != atomX.size() 
             && max_proton.size() != heavyList.size()) {
         error("MAXNH should match the number of molecules or heavy atoms");
       }
       if(max_proton.size() == atomX.size()) {
         // check whether max_proton are the same for atoms in the same mol
         auto tmp = max_proton; max_proton.clear();
         for(auto it = heavyList.begin(); it != heavyList.end(); ++it) {
            int s = tmp[(*it)[0]];
            max_proton.push_back(s);
            for(unsigned j : *it) 
               if(tmp[j] != s) {
                  string str;
                  Tools::convert(static_cast<unsigned>(it-heavyList.begin()),str);
                  error("Inconsistent MAXNH for molecule " + str);
               }
         }
       }
    }
    //int total = 0;
    //for(auto l : heavyList) total += max_proton[l[0]];
    int total = std::accumulate(max_proton.begin(), max_proton.end(), 0);
    if(total + 1!= int(natomh)) { 
       log.printf("total max_proton = %d, natomh = %d\n", total, natomh);
       error("Number of excess proton should be exactly one");
    }
    // assume dmax and RATIO * dmax are the max possible O-H bond lengths
    // for hydronium and water respectively
    vector<double> dmaxs;
    //for(unsigned i = 0; i < atomX.size(); ++i) 
    //   for(unsigned j = 0; j < atomX.size(); ++j)
    //      dmaxs.push_back(sf_delta_library[sf_indexes[i][j]].get_dmax());
    for(const auto& sf : sf_delta_library) dmaxs.push_back(sf.get_dmax());
    for(const auto& sf : sf_roo_library) dmaxs.push_back(sf.get_dmax());
    //rlisto = sf.get_dmax() * (2 + 0.75 * 4 + 3);
    //if(do_fourth) rlisto += sf.get_dmax() * (0.75 *2 + 1);
    parse("RATIO", ratio);
    max_dmax2 = *max_element(dmaxs.begin(), dmaxs.end());
    rlisto = max_dmax2 * (2 + ratio * 4 + 3);
    if(do_fourth) rlisto += max_dmax2 * (ratio * 2 + 1);
    rlisto = rlisto * rlisto;
    max_dmax2 = max_dmax2 * max_dmax2;
    atomH.insert(atomH.end(),atomX.begin(),atomX.end()); //concatenate 
    //vector<AtomNumber>::iterator atomH_end=atomH.begin()+natomh;
    auto atomH_end=atomH.begin()+natomh;
    //vector<AtomNumber>::iterator atomNP_end=atomH_end+natomnp;
    //natomnp+=natomh; //now natomnp is the number of non pair atoms including H
    //natomp=atomP.size();
    natom = atomH.size(); //natom=natomp+natomnp;
    for(unsigned i = 0; i < natomh; ++i) hindexes.push_back(i);
    for(unsigned i = natomh; i < natom; ++i) oindexes.push_back(i);
    //write log
    log.printf("  %d hydrogen atoms and %d heavy atoms involved\n",\
    natomh,atomX.size());
    log<<"  serial number of hydrogen atoms:\n  ";
    for(auto it=atomH.begin();it!=atomH_end;++it) {
        log.printf("%d ",it->serial()); 
        if((it-atomH.begin() + 1) % 20 == 0) log.printf("\n  ");
    }
    log<<"\n  serial number of heavy atoms:\n  ";
    for(auto it=atomH_end;it!=atomH.end();++it) {
        log.printf("%d ",it->serial());
        if((it-atomH_end + 1) % 20 == 0) log.printf("\n  ");
    }
    //if(!atomP.empty()) {
    //    log<<"\n  serial number of pair atoms:\n  ";
    //    for(vector<AtomNumber>::iterator it=atomNP_end;it!=atomH.end();\
    //    ++it) {
    //        log.printf("%d ",(*it).serial());
    //    }
    //}
    //else
    //    log<<"\n  no pair atoms specified, work in single proton mode";
    log <<"\n  max num of protons for heavy atoms:\n  ";
    for(auto it=wx.begin();it!=wx.end();++it)    {
        log<<*it<<" ";
        if((it-wx.begin() + 1) % 20 == 0) log.printf("\n  ");
    }
    log <<"\n  charges of heavy atoms:\n  ";
    for(auto it=heavy_charges.begin();it!=heavy_charges.end();++it) {
        log<<*it<<" ";
        if((it-heavy_charges.begin() + 1) % 20 == 0) log.printf("\n  ");
    }
    log <<"\n  charges of protons:\n  ";
    for(auto it=proton_charges.begin();it!=proton_charges.end();++it) {
        log<<*it<<" ";
        if((it-proton_charges.begin() + 1) % 20 == 0) log.printf("\n  ");
    }
    log <<"\n  molecular id for heavy atoms:\n  ";
    for(auto it=molid.begin();it!=molid.end();++it)    {
        log<<*it<<" ";
        if((it-molid.begin() + 1) % 20 == 0) log.printf("\n  ");
    }
    log <<"\n  max num of protons for molecules:\n  ";
    for(auto it=heavyList.begin();it!=heavyList.end();++it)    {
        log<<max_proton[(*it)[0]]<<" ";
        if((it-heavyList.begin() + 1) % 20 == 0) log.printf("\n  ");
    }
    log<<"\n";
//#ifdef __PRINTF_DEBUG__
//    for(unsigned i = 0; i < atomX.size(); ++i)
//       for(unsigned j = 0; j < atomX.size(); ++j)
//         log.printf("  sf[%d][%d] info: r0 = %s\n", i, j,
//                                (sf_delta_library[sf_[i][j]].description()).c_str() );
//
//#else
    log.printf("  switching function of delta between water: r0 = %s\n",
                                        (sf_delta_library[0].description()).c_str() );
    for(const auto& p : map2sf_delta_index) {
       if(p.first > 0)
          log.printf("  switching function of delta "
                "for atom serial %ld to water: r0 = %s\n",
                p.first, sf_delta_library[p.second].description().c_str());
       else
          log.printf("  switching function of delta "
                "for water to atom serial %ld: r0 = %s\n",
                -p.first, sf_delta_library[p.second].description().c_str());
    }
    log.printf("  switching function of r_OO between water: r0 = %s\n",
                                        (sf_roo_library[0].description()).c_str() );
    for(const auto& p : map2sf_roo_index) {
       log.printf("  switching function of r_OO "
             "for atom serial %ld to water: r0 = %s\n",
             p.first, sf_roo_library[p.second].description().c_str());
    }
//#endif
    string fn_topo;
    parse("INIT_TOPO", fn_topo);
    guess_initial_topo = fn_topo.empty();
    if(guess_initial_topo) 
       log.printf("  initial topology not provided,"
             " will construct at the first step");
    else {
       fs_topo.open(fn_topo);
       if(!fs_topo.is_open()) error("Cannot open " + fn_topo);
       log.printf("  will use initial topology given in %s\n",
             fn_topo.c_str());
    }
    if(do_pbc){
        log<<"  PBC will be considered\n";
    } else {
        log<<"  PBC will be ignored\n";
    }
    if(use_mask){
        log<<"  masks will be used to accelerate calculations\n";
        log<<"  only heavy atoms within " << std::sqrt(rlisto)
           << " of refpos will contribute to CEC\n";
    }

    if(print_delta) {
       log<<"WARNING: PRINT_DELTA only makes sense for proton in water\n";
    }
    if(print_shells) {
       log<<"WARNING: PRINT_SHELLS does not treat reaction rings properly currently\n";
       log.printf("shells: #step refpos.x refpos.y refpos.z");
       log.printf(" shell0.x shell0.y shell0.z");
       log.printf(" shell1.x shell1.y shell1.z");
       log.printf(" shell2.x shell2.y shell2.z");
       log.printf(" shell3.x shell3.y shell3.z\n");
    }
    if(print_o_h) {
       log.printf("O_H: #step refpos.x refpos.y refpos.z");
       log.printf(" O.x O.y O.z");
       log.printf(" H.x H.y H.z\n");
    }
    if(print_acceptors) {
       log<<"WARNING: PRINT_ACCEPTORS only makes sense for proton in water\n";
       log.printf("acceptors: #step Ox Oy Oz deltax deltay deltaz\n");
    }
    if(print_nstates) {
       log<<"WARNING: PRINT_NSTATES does not treat reaction rings properly currently\n";
       log.printf("nstates: # will print serial of leading %u states\n", 
             print_nstates);
    }
    //log<<"  second shell will be weighted by "<< A2 << '\n';
    //log<<"  third shell will be weighted by "<< A3 << '\n';
    //if(do_updateTopo){
    //    log<<"  Topology will be updated every step\n";
    //} else {
    //    log<<"  Do not consider topology in calculation\n";
    //}
    log<<"  charge and mass of RMDCEC are set to be 0 and 1 in this implementation\n";
    //log<<"  dist between nonbonded heavy atom and hydrogen lower than "
    //   << rwarn << " will give warnings\n";
    //log<<"  dist between nonbonded heavy atom and hydrogen lower than "
    //   << rerror << " will stop the program\n";
    //initialize indexes, allocate positions, masses, ...
    requestAtoms(atomH);

//@@@
fullAtomList = atomH;

#ifdef __FAKE_VATOM__
    parseVector("MOLIDS2WATCH", molids_2b_watched);
    if(molids_2b_watched.empty()) error("MOLIDS2WATCH not specified");
    if(molids_2b_watched.size() > 3) error("cannot specify more than 3 MOLIDS2WATCH");
    log.printf("  molids to be watched:");
    for(unsigned molid : molids_2b_watched) log.printf(" %u", molid);
    log.printf("\n");
    ci2_2b_watched = Vector(0.0, 0.0, 0.0);
#endif

    checkRead();

}

inline
SwitchingFunction RMDCEC::sfij_delta(unsigned i, unsigned j) const {
   auto p = map2sf_delta_index.find(i + 1);
   if(p != map2sf_delta_index.end()) {
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
#endif
#ifdef __DEBUG_07312019__
if((i == 2 && j == 4) or (i == 4 && j == 2)) {
#endif

#ifdef __DEBUG__

      printf("atom index %u and %d on rank %d, sf_delta[%u]: ",
            i, j, comm.Get_rank(), p->second);
      std::fflush(stdout);
      printf("%s\n", 
            sf_delta_library[p->second].description().c_str());
#endif

#ifdef __DEBUG_07312019__
}
#endif
#ifdef __DEBUG_LOW__
}                              
#endif

      return sf_delta_library[p->second];
   } 
   p = map2sf_delta_index.find(-static_cast<long>(j + 1));
   if(p != map2sf_delta_index.end()) {
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
#endif
#ifdef __DEBUG_07312019__
if((i == 2 && j == 4) or (i == 4 && j == 2)) {
#endif

#ifdef __DEBUG__
      printf("atom index %u and %d on rank %d, sf_delta[%u]: ",
            i, j, comm.Get_rank(), p->second);
      std::fflush(stdout);
      printf("%s\n", 
            sf_delta_library[p->second].description().c_str());
#endif
#ifdef __DEBUG_07312019__
}
#endif
#ifdef __DEBUG_LOW__
}                               
#endif
      return sf_delta_library[p->second];
   }
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
#endif
#ifdef __DEBUG_07312019__
if((i == 2 && j == 4) or (i == 4 && j == 2)) {
#endif

#ifdef __DEBUG__
   printf("atom index %u and %d on rank %d, sf_delta[%u]: ",
         i, j, comm.Get_rank(), 0);
   std::fflush(stdout);
   printf("%s\n", 
         sf_delta_library[0].description().c_str());
#endif

#ifdef __DEBUG_07312019__
}
#endif
#ifdef __DEBUG_LOW__
}                                
#endif
   return sf_delta_library[0];
}

inline
SwitchingFunction RMDCEC::sfij_roo(unsigned i, unsigned j) const {
   auto p = map2sf_roo_index.find(i + 1);
   if(p != map2sf_roo_index.end()) {
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
#endif
#ifdef __DEBUG__
      printf("atom index %u and %d on rank %d, sf_roo[%u]: ",
            i, j, comm.Get_rank(), p->second);
      std::fflush(stdout);
      printf("%s\n", 
            sf_roo_library[p->second].description().c_str());
#endif
#ifdef __DEBUG_LOW__
}                             
#endif
      return sf_roo_library[p->second];
   } 
   p = map2sf_roo_index.find(j + 1);
   if(p != map2sf_roo_index.end()) {
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
#endif
#ifdef __DEBUG__
      printf("atom index %u and %d on rank %d, sf_roo[%u]: ",
            i, j, comm.Get_rank(), p->second);
      std::fflush(stdout);
      printf("%s\n", 
            sf_roo_library[p->second].description().c_str());
#endif
#ifdef __DEBUG_LOW__
}                            
#endif
      return sf_roo_library[p->second];
   }
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
#endif
#ifdef __DEBUG__
   printf("atom index %u and %d on rank %d, sf_roo[%u]: ",
         i, j, comm.Get_rank(), 0);
   std::fflush(stdout);
   printf("%s\n", 
         sf_roo_library[0].description().c_str());
#endif
#ifdef __DEBUG_LOW__
}                              
#endif
   return sf_roo_library[0];
}

// return the index of the excess hydrogen atom in hindexes
// build bondedOList and bondedHList
// the *indexes contain the indexes in fullAtomList
// this should only used in the first time
unsigned RMDCEC::readTopo()
{
   bondedOList = vector<long int>(hindexes.size(),-1);
   string line;
   while(getline(fs_topo, line)) {
      unsigned i, j;
      stringstream ss(line);
      ss >> i >> j;
      bondedOList.at(i) = j;
   }
   fs_topo.close();

   for(unsigned i = 0; i < bondedOList.size(); ++i) {
      if(bondedOList[i] == -1) {
         string index; Tools::convert(i, index);
         plumed_merror("atom H of index " + index + 
               " in this action is not bonded to any heavy atom");
      }
   }

   // construct bondedHList using bondedOList
   for(auto it = bondedOList.begin(); it != bondedOList.end(); ++it) 
      bondedHList[*it].push_back(it - bondedOList.begin());

   // find the pivot heavy atom, i.e. the atom that has WX protons
   for(unsigned i = 0; i < bondedHList.size(); ++i) {

      if(static_cast<int>(bondedHList[i].size()) != max_proton[i]) {

         if(static_cast<int>(bondedHList[i].size()) == wx[i]) {
            // Note that it is possible that there are more than one pivot atom
            // for instance, HSP
            // we just use the first one we meet in this loop
            // Note that the pivot heavy atom may have multiple H atoms
            // we just return the last one
            return bondedHList[i].back();
         }
         else { 
            string index; Tools::convert(i, index);
            plumed_merror("ill bonding topology for atom X of index " 
                  + index + "in this action");
         }

      }

   }

   plumed_merror("should not reach here");
   return 0;
}

// return the index of the excess hydrogen atom in hindexes
// build bondedOList and bondedHList
// the *indexes contain the indexes in fullAtomList
// this should only used in the first time
unsigned RMDCEC::getTopo()
{

    if(!guess_initial_topo) return readTopo();

    // I will use bondedOList to construct bondedHList
    // bondedOList[i] stores the index of bonded heavy atom of hydrogen i
    bondedOList=vector<long int>(hindexes.size(),-1);
    vector<int> counts;
#ifdef __WOPAIR__
    vector<unsigned> index_pair; // index_pair[i] stores index in oindexes
    for(unsigned i=0;i<oindexes.size();++i) {
        unsigned io = oindexes[i];
        if(true) { //if(io<natomnp) {
            counts.push_back(int(0.5+wx[io-natomh]));
        } else {
            index_pair.push_back(i);
            counts.push_back(int(0.5+wx[io-natomh]*natomp));
        }
    }
#else
    //counts = max_proton;
    //num_proton = vector<int>(molid.size());
#endif
    //if(!index_pair.empty() && index_pair.size()!=natomp) 
    //    plumed_merror("ERROR: pair atom missing");

#ifdef __WOPAIR__
    //loop for O:
    //  find the closest unassigned H until its count is zero 
    for(unsigned j=0;j<oindexes.size();++j) { 
        while(counts[j]) {
            unsigned closest_h=0;
            double min_d2=std::numeric_limits<double>::max();
            for(unsigned i=0;i<hindexes.size();++i) {
                if(bondedOList[i]!=-1) continue;
                double d2;
                //distance between two atoms
                if(true) { //if(oindexes[j]<natomnp) {
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
            bondedOList[closest_h]=j;
            //modify counts
            unsigned jo=oindexes[j];
            if(true) { //if(jo<natomnp) {
                counts[j]--;
            } else {
                for(unsigned k=0;k<natomp;++k) {
                    counts[index_pair[k]]--;
                }
            }
        }
    }
#else
    // for every molecule, find the closest max_proton protons
    for(unsigned mol_id = 0; mol_id < heavyList.size(); ++mol_id) {
       vector<tuple<double,unsigned,unsigned>> closest;
       for(unsigned h : hindexes) {
          if(bondedOList[h] != -1) continue;
          closest.push_back( dist2MolH(mol_id,h) );
       }
       std::partial_sort(closest.begin(), 
             closest.begin() + max_proton[mol_id], closest.end());
       for(int i = 0; i < max_proton[mol_id]; ++i)
          bondedOList[std::get<2>(closest[i])] = std::get<1>(closest[i]);
    }
#endif

#ifdef __DEBUG_07312019__
if(comm.Get_rank() == 0)
for(unsigned i = 0; i < bondedOList.size(); ++i) {
   printf("bondedOList[%u] = %ld\n", i, bondedOList[i]);
   std::fflush(stdout);
}
#endif

    // find the -1 in bondedOList and assign the closest O to that
    unsigned rc=0;
    unsigned count_rc=0;
    //plumed_assert(bondedOList.size()==hindexes.size());
    for(unsigned i=0;i<bondedOList.size();++i) {
        if(bondedOList[i]==-1) {
            rc=i; count_rc++;
        }
    }
    plumed_assert(count_rc==1);
    double min_d2 = std::numeric_limits<double>::max();
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
    bondedOList[rc] = closet_o;

    //construct bondedHList using bondedOList
    for(auto it = bondedOList.begin(); it != bondedOList.end(); ++it) 
       bondedHList[*it].push_back(it - bondedOList.begin());

#ifndef __WOPAIR__
    for(unsigned o : bondedOList) {
       //if( (o == closet_o && bondedHList[o].size() > wx[o]) ||
       //    (o != closet_o && bondedHList[o].size() >= wx[o]) ) {
       if(static_cast<int>(bondedHList[o].size()) > wx[o]) {
          string str1, str2;
          Tools::convert(fullAtomList[oindexes[o]].serial(),str1);
          Tools::convert(getStep(), str2);
          plumed_merror("ERROR: atom "+str1+" is overbonded at step "+str2);
       }
    }
#endif
    
    return rc;
}

Vector RMDCEC::initPos() {
    //set bondedOList as class member for use in next step //vector<long int> bondedOList;
    unsigned rc = getTopo();//,bondedOList);
    pivot = bondedOList[rc];
    return getPosition(rc);
}

void RMDCEC::normDistOH(unsigned i, unsigned j, Vector& result) const
{
   result = do_pbc ? pbcDistance(
            getPosition(hindexes[j]),
            getPosition(oindexes[i])) :
      delta(
            getPosition(hindexes[j]),
            getPosition(oindexes[i]));
   result /= std::sqrt(result.modulo2());
}

Vector RMDCEC::normDistOH(unsigned i, unsigned j) const
{
   Vector result = do_pbc ? pbcDistance(
              getPosition(hindexes[j]),
              getPosition(oindexes[i])) :
        delta(
              getPosition(hindexes[j]),
              getPosition(oindexes[i]));
   return (result / std::sqrt(result.modulo2()));
}

Vector RMDCEC::normDistOO(unsigned i, unsigned j) const
{
   Vector result = do_pbc ? pbcDistance(
              getPosition(oindexes[j]),
              getPosition(oindexes[i])) :
        delta(
              getPosition(oindexes[j]),
              getPosition(oindexes[i]));
   return (result / std::sqrt(result.modulo2()));
}

double RMDCEC::dist2OH(unsigned i, unsigned j) const
{
   Vector distance = do_pbc ?
          pbcDistance(
                getPosition(oindexes[i]),
                getPosition(hindexes[j])) :
          delta(
                getPosition(oindexes[i]),
                getPosition(hindexes[j]));
   return distance.modulo2();
}

double RMDCEC::dist2OO(unsigned i, unsigned j) const
{
   Vector distance = do_pbc ?
          pbcDistance(
                getPosition(oindexes[i]),
                getPosition(oindexes[j])) :
          delta(
                getPosition(oindexes[i]),
                getPosition(oindexes[j]));
   return distance.modulo2();
}

tuple<double,unsigned,unsigned> 
RMDCEC::dist2MolH(unsigned mol_id, unsigned h) const
{
   vector<pair<double,unsigned>> t;
   for(unsigned j : heavyList[mol_id]) 
      t.push_back(make_pair( dist2OH(j,h), j ));
   auto tt = min_element(t.begin(), t.end());
   return make_tuple(tt->first, tt->second, h);
}

double RMDCEC::calculateDelta(unsigned i, unsigned j, unsigned k) const 
{
   //donor i acceptor j hydrogen k
   double d1 = dist2OH(i, k), d2 = dist2OH(j, k);
   return std::sqrt(d2) - std::sqrt(d1);
}

bool RMDCEC::shouldSkip(double d1, double d2, double d0) const
{
   if(d2 < d1) return false;
   double t = d1 + d2 - d0;
   return t * t > 4 * d1 * d2;
}

void RMDCEC::buildReactionTree() {

   // use masks to accelerate calculations
   
   //vector<bool> omasks(oindexes.size(), true);
   //if(use_mask) {
   //   for(auto io = oindexes.begin(); io != oindexes.end(); ++io) {
   //      if(relpos[*io].modulo2() > rlisto)
   //         omasks[io - oindexes.begin()] = false;
   //   }
   //}
   vector<unsigned> masked_oindexes;
   if(use_mask) {
      for(auto io = oindexes.begin(); io != oindexes.end(); ++io) {
         if(relpos[*io].modulo2() <= rlisto)
            masked_oindexes.push_back(io - oindexes.begin());
      }
   } else {
      // this should be faster than std::itoa
      for(unsigned i = 0; i < oindexes.size(); ++i)
         masked_oindexes.push_back(i);
   }

//printf("0: \n");//@@@

   reactionTree.clear(); 
   gmat_delta.clear(); dgmat_delta.clear();
   gmat_roo.clear(); dgmat_roo.clear();
   unsigned pivot_molid = molid[pivot];
   vector<unsigned> ilist = heavyList[pivot_molid];
   for(unsigned i : ilist) {
      //maybe better to use iterator here but I want to make minor changes
      for(unsigned j_=0; j_ < masked_oindexes.size(); ++j_) {
         unsigned j = masked_oindexes[j_];
      //for(unsigned j=0; j < oindexes.size(); ++j) {
      //   if(!omasks[j]) continue;
         //if(bondedHList[j].size() >= int(wx[j]+1.0) ) continue;
         if(static_cast<int>(bondedHList[j].size()) >= wx[j] ) continue;
         if(molid[j] == pivot_molid) continue;

//printf("1: \n");//@@@

         for(unsigned hi : bondedHList[i]) {

//printf("2: \n");//@@@

            //FIXME this can be optimized by maintaining a bondedHLength2[i]
            // bondedHLength2 can be constructed in the loop of first j
            // then it is not changing for other j
            double d1 = dist2OH(i, hi), d2 = dist2OH(j, hi);

//printf("i = %d hi = %d j = %d d1 = %f d2 = %f\n",
//      fullAtomList[i].serial(),hi,j,d1,d2);//@@@

            //if(shouldSkip(d1, d2, sf.get_dmax2())) continue;
            auto sf_delta = sfij_delta(i, j);
            if(
                  shouldSkip(d1, d2, sf_delta.get_dmax2())
                  ) continue;
            //if(shouldSkip(d1, d2, max_dmax2)) continue;

//printf("3: \n");//@@@

            auto chain = vector<array<unsigned,3>>(
//                  1, array<unsigned,3>({pivot, j, hi}) );
                  1, array<unsigned,3>({i, j, hi}) );
            double d = std::sqrt(d2) - std::sqrt(d1);
            //double df, f = sf.calculate(d, df); df *= d;

//if(getSerial(i)==19 and getSerial(j) == 1688) { // @@@
//   printf("i = %d j = %d r0 = %f d0 = %f\n",i,j,sf[i][j].get_r0(),sf[i][j].get_d0());
//}

            double df, f = sf_delta.calculate(d, df); 
#ifdef __DEBUG__
#ifdef __DEBUG_07312019__
if((i == 2 && j == 4) or (i == 4 && j == 2)) {
#endif
            printf("d = %f f = %f df = %f \n", d, f, df);
#ifdef __DEBUG_07312019__
}
#endif
#endif
            if(sf_delta.divided()) df *= d;
            auto g_delta = vector<double>(1, f);
            auto dg_delta = vector<double>(1, df);
            gmat_delta.push_back(g_delta); dgmat_delta.push_back(dg_delta);
            reactionTree.push_back(chain);
            double roo = std::sqrt(dist2OO(i, j));
            auto sf_roo  = sfij_roo(i, j);
            f = sf_roo.calculate(roo, df);
            if(sf_roo.divided()) df *= roo;
            vector<double> g_roo, dg_roo;
            if(d >= abs(sf_delta.get_d0())) {
               g_roo = vector<double>(1, f);
               dg_roo = vector<double>(1, df);
            } else {
               g_roo = vector<double>(1, 1.0 / f);
               dg_roo = vector<double>(1, - df / f / f);
            }
            gmat_roo.push_back(g_roo); dgmat_roo.push_back(dg_roo);

            // start 2nd shell
            //bondedHList[j].push_back(hi); FIXME
            unsigned j_molid = molid[j];
            vector<unsigned> jlist = heavyList[j_molid];
            for(unsigned jb : jlist) {
               for(unsigned k_ = 0; k_ < masked_oindexes.size(); ++k_) {
                  unsigned k = masked_oindexes[k_];
               //for(unsigned k = 0; k < oindexes.size(); ++k) {
               //   if(!omasks[k]) continue;
                  //if(bondedHList[k].size() >= int(wx[k]+1.0) ) continue;
                  if(static_cast<int>(
                           (bondedHList[k].size())) >= wx[k] ) continue;
                  if(molid[k] == j_molid || molid[k] == pivot_molid) continue;
                  for(unsigned hj : bondedHList[jb]) {
                     double d1 = dist2OH(jb, hj), d2 = dist2OH(k, hj);
                     //if(shouldSkip(d1, d2, sf.get_dmax2())) continue;
                     auto sf = sfij_delta(jb, k);
                     if(shouldSkip(d1, d2, sf.get_dmax2())) continue;
                     //if(shouldSkip(d1, d2, max_dmax2)) continue;
                     chain.push_back(array<unsigned,3>({jb,k,hj}));
                     double d = std::sqrt(d2) - std::sqrt(d1);
                     //double df, f = sf.calculate(d, df); df *= d;
                     double df, f = sf.calculate(d, df); 
#ifdef __DEBUG__
            printf("   d = %f f = %f df = %f \n", d, f, df);
#endif
                     if(sf.divided()) df *= d;
                     g_delta.push_back(f); dg_delta.push_back(df);
                     reactionTree.push_back(chain);
                     gmat_delta.push_back(g_delta); dgmat_delta.push_back(dg_delta);
                     double roo = std::sqrt(dist2OO(k, jb));
                     auto sf_roo  = sfij_roo(k, jb);
                     f = sf_roo.calculate(roo, df);
#ifdef __DEBUG__
            printf("   roo = %f f = %f df = %f \n", roo, f, df);
#endif
                     if(sf_roo.divided()) df *= roo;
                     if(d >= abs(sf_delta.get_d0())) {
                        g_roo.push_back(f);
                        dg_roo.push_back(df);
                     } else {
                        g_roo.push_back(1.0 / f);
                        dg_roo.push_back(- df / f / f);
                     }
                     gmat_roo.push_back(g_roo); dgmat_roo.push_back(dg_roo);
                     
                     // start 3rd shell
                     //bondedHList[k].push_back(hj); FIXME
                     unsigned k_molid = molid[k];
                     vector<unsigned> klist = heavyList[k_molid];
                     for(unsigned kb : klist) {
                        for(unsigned l_=0; l_ < masked_oindexes.size(); ++l_) {
                           unsigned l = masked_oindexes[l_];
                        //   if(!omasks[l]) continue;
                           //if(bondedHList[l].size() >= int(wx[l]+1.0) ) continue;
                           if(static_cast<int>(
                                    bondedHList[l].size()) >= wx[l]) continue;
                           if(molid[l] == k_molid || 
                                 molid[l] == j_molid || molid[l] == pivot_molid)
                              continue;
                           for(unsigned hl : bondedHList[kb]) {
                              double d1 = dist2OH(kb, hl), d2 = dist2OH(l, hl);
                              //if(shouldSkip(d1, d2, sf.get_dmax2())) continue;
                              auto sf = sfij_delta(kb, l);
                              if(shouldSkip(d1,d2,sf.get_dmax2())) 
                              //if(shouldSkip(d1,d2,max_dmax2)) 
                                                                       continue;
                              chain.push_back(array<unsigned,3>({kb,l,hl}));
                              double d = std::sqrt(d2) - std::sqrt(d1);
                              //double df, f = sf.calculate(d,df); df *= d;
                              double df, f = sf.calculate(d,df); 
#ifdef __DEBUG__
            printf("      d = %f f = %f df = %f \n", d, f, df);
#endif
                              if(sf.divided()) df *= d;
                              g_delta.push_back(f); dg_delta.push_back(df);
                              reactionTree.push_back(chain);
                              gmat_delta.push_back(g_delta); dgmat_delta.push_back(dg_delta);
                              double roo = std::sqrt(dist2OO(l, kb));
                              auto sf_roo  = sfij_roo(l, kb);
                              f = sf_roo.calculate(roo, df);
#ifdef __DEBUG__
            printf("      roo = %f f = %f df = %f \n", roo, f, df);
#endif
                              if(sf_roo.divided()) df *= roo;
                              if(d >= abs(sf_delta.get_d0())) {
                                 g_roo.push_back(f);
                                 dg_roo.push_back(df);
                              } else {
                                 g_roo.push_back(1.0 / f);
                                 dg_roo.push_back(- df / f / f);
                              }
                              gmat_roo.push_back(g_roo); dgmat_roo.push_back(dg_roo);


                              //if there is fourth shell, then it will go here
                              if(do_fourth) {
                                 unsigned l_molid = molid[l];
                                 vector<unsigned> llist = heavyList[l_molid];
                                 for(unsigned lb : llist) {
                                    for(unsigned m_=0; m_ < masked_oindexes.size(); ++m_) {
                                       unsigned m = masked_oindexes[m_];
                                       if(static_cast<int>(
                                                bondedHList[m].size()) >= wx[m]) continue;
                                       if(molid[m] == l_molid || molid[m] == k_molid ||
                                             molid[l] == j_molid || molid[l] == pivot_molid)
                                          continue;
                                       for(unsigned hm : bondedHList[lb]) {
                                          double d1 = dist2OH(lb, hm), d2 = dist2OH(m, hm);
                                          //if(shouldSkip(d1, d2, sf.get_dmax2())) continue;
                                          auto sf = sfij_delta(lb, m);
                                          if(shouldSkip(d1, d2, sf.get_dmax2())) continue;
                                          //if(shouldSkip(d1, d2, max_dmax2)) continue;
                                          chain.push_back(array<unsigned,3>({lb,m,hm}));
                                          double d = std::sqrt(d2) - std::sqrt(d1);
                                          //double df, f = sf.calculate(d, df); df *= d;
                                          double df, f = sf.calculate(d, df); 
                                          if(sf.divided()) df *= d;
                                          g_delta.push_back(f); dg_delta.push_back(df);
                                          reactionTree.push_back(chain);
                                          gmat_delta.push_back(g_delta); dgmat_delta.push_back(dg_delta);
                                          double roo = std::sqrt(dist2OO(m, lb));
                                          auto sf_roo  = sfij_roo(m, lb);
                                          f = sf_roo.calculate(roo, df);
                                          if(sf_roo.divided()) df *= roo;
                                          if(d >= abs(sf_delta.get_d0())) {
                                             g_roo.push_back(f);
                                             dg_roo.push_back(df);
                                          } else {
                                             g_roo.push_back(1.0 / f);
                                             dg_roo.push_back(- df / f / f);
                                          }
                                          gmat_roo.push_back(g_roo); dgmat_roo.push_back(dg_roo);
                                          //if there is fifth shell, then it will go here
                                          {;}
                                          chain.pop_back(); 
                                          g_delta.pop_back(); dg_delta.pop_back();
                                          g_roo.pop_back(); dg_roo.pop_back();
                                          
                                       } //hm
                                    } //m
                                 } //lb
                              }


                              chain.pop_back(); 
                              g_delta.pop_back(); dg_delta.pop_back();
                              g_roo.pop_back(); dg_roo.pop_back();
                              
                           } //hl
                        } //l
                     } //kb
                     chain.pop_back();
                     g_delta.pop_back(); dg_delta.pop_back();
                     g_roo.pop_back(); dg_roo.pop_back();
                  }
               }
            }
            chain.pop_back();
            g_delta.pop_back(); dg_delta.pop_back();
            g_roo.pop_back(); dg_roo.pop_back();
         }
      }
   }
}

void RMDCEC::printNstates(unsigned nstates) {

   if(nstates == 0) {
      return;
   } else if(nstates == 1) {
   }

   // get the indexes of states with n leading ci2
   // this should be in <algorithm> ??!!
   priority_queue<pair<double,unsigned>> q;
   for(unsigned i = 0; i < ci_2.size(); ++i) {
      q.emplace(ci_2[i], i);
   }
   vector<unsigned> states;
   for(unsigned i = 0; i < nstates; ++i) {
      if(q.empty()) break;
      states.push_back(q.top().second);
      q.pop();
   }

   vector<unsigned> serials2print;

   for(auto it = states.rbegin(); it != states.rend(); ++it) {
      unsigned istate = *it;

      // if this state is pivot
      if(istate + 1 == ci_2.size()) {
         serials2print.push_back(fullAtomList[oindexes[pivot]].serial());
         // if this is the largest ci state
         if(it + 1 == states.rend()){
            for(auto& h : bondedHList[pivot]) {
               serials2print.push_back(fullAtomList[h].serial());
            }
         }
         // if this is not the largest ci state
         else {
            // get the reaction that transfers the pivot proton to largest ci_2 state
            auto& r = reactionTree[states[0]][0];
            for(auto& h : bondedHList[pivot]) {
               if(h == r[2]) continue;
               serials2print.push_back(fullAtomList[h].serial());
            }
         }
      }
      // if this state is not pivot
      else {
         // get the last step of PT reaction chain
         auto& r = reactionTree[istate].back();
         unsigned A = r[1];
         unsigned H = r[2];
         serials2print.push_back(fullAtomList[oindexes[A]].serial());
         for(auto& h : bondedHList[A]) {
            serials2print.push_back(fullAtomList[h].serial());
         }

         // if this is the largest ci state
         if(it + 1 == states.rend()){
            serials2print.push_back(fullAtomList[hindexes[H]].serial());
         }
         // if this is not the largest ci state
         else {
            ;
         }
      }
   }

   log.printf("nstates:");
   for(auto& s : serials2print) {
      log.printf(" %u", s);
   }
   log.printf("\n");
}

void RMDCEC::calculate() {

#ifdef __FAKE_VATOM__
    ci2_2b_watched = Vector(0.0, 0.0, 0.0);
#endif

    Vector pos;
    //Vector vdij; //\vec{r}^{H_i}-\vec{r}^{X_j}
    //double dij,fij,dfij; //dij=norm of vdij, fij=f_{sw}(dij),dfij=f'_{sw}(dij)
    vector<Tensor> deriv(natom);
    if(firsttime) {
        initPos(); //getTopo was done in initPos()
        log<<"Info: refpos for RMDCEC is "
           <<getPosition(oindexes[pivot])[0] 
           <<' '<< getPosition(oindexes[pivot])[1] 
           <<' '<< getPosition(oindexes[pivot])[2] <<'\n';
        //updateList(); //update the oindexes and fullBondedHList
        firsttime = false;
    }
    // BY NOW, I have molid bondedOList bondedHList heavyList 
    // I also have pivot
    
    //else if(getStep()%stride==0) {
    //    updateList();
    //}

#ifdef __DEBUG_07312019__
if(comm.Get_rank() == 0)
for(unsigned i = 0; i < 8; ++i) {
   printf("bondedHList[%u] = ", i);
   std::fflush(stdout);
   for(const auto & h : bondedHList[i]) {
      printf("%u ", h);
      std::fflush(stdout);
   }
   printf("\n");
   std::fflush(stdout);
}
#endif

//modify pivot to see whether an alternative pivot will cause a large error
//@@@ do pivot bondeHList
unsigned proton = 0;
if(h_serial!=static_cast<unsigned>(-1) && o_serial!=static_cast<unsigned>(-1)) {
for(auto hi=hindexes.begin(); hi!=hindexes.end(); ++hi) {
   //find the changing proton
   if(fullAtomList[*hi].serial() == h_serial) { proton = hi - hindexes.begin(); break;}
}
// erase the proton in bondedHList[pivot]
bondedHList[pivot].erase(std::find(bondedHList[pivot].begin(), bondedHList[pivot].end(), proton));
for(auto oi=oindexes.begin(); oi!=oindexes.end(); ++oi) {
   // change pivot to this serial
   if(fullAtomList[*oi].serial() == o_serial) {
      pivot = oi - oindexes.begin();
      bondedHList[pivot].push_back(proton);
   }
}
}


    unsigned nh = hindexes.size(), nx = oindexes.size();
    Vector refpos = getPosition(pivot + nh);
    relpos = vector<Vector>(nh+nx);
    //unsigned np=0;
    for(unsigned j=0;j<nx;j++) {
        //if(oindexes[j]>=natomnp)  np++;
    }
    //if(np!=0&&np!=natomp) plumed_merror("ERROR: np!=natomp");
    //unsigned nnp = nh+nx-np, n = nh+nx;
    //pbc all the heavy by refpos
    for(unsigned j=0;j<nx;++j) {
        Vector dist=do_pbc?
            pbcDistance(refpos,getPosition(oindexes[j])):
            delta(refpos,getPosition(oindexes[j]));
        //relpos[j+nh] = refpos + dist;
        relpos[j+nh] = dist;
    }
    //pbc all the H by refpos
    for(unsigned i=0;i<nh;++i) {
        Vector dist=do_pbc?
            pbcDistance(refpos,getPosition(hindexes[i])):
            delta(refpos,getPosition(hindexes[i]));
        //relpos[i] = refpos + dist;
        relpos[i] = dist;
    }

    buildReactionTree();


#ifdef __NOISY__
if(comm.Get_rank()==0) {
   printf("%ld states found in step %ld\n",
         gmat_delta.size(),getStep()); //@@@
}
#endif
#ifdef __DOFITTING__
vector<double> ci_2_on_atom(nx); // for fitting @@@
fs_out << getStep() << ' '; // for fitting @@@
#endif

    pos = refpos;

    /*--- calculate the ci_2 and its derivs ---*/
//    vector<double> ci_2;
    ci_2.clear();
    //\[pd] ci_2 / \[pd] delta before normed
    vector<vector<double>> pc_pdelta; 
    //\[pd] ci_2 / \[pd] r_oo before normed
    vector<vector<double>> pc_proo;
    for(unsigned i = 0; i < gmat_delta.size(); ++i) {
       //double ci = A[gmat[i].size() - 1];
       //for(double g : gmat[i]) {
       //   ci *=  g;
       //}
       vector<double> pd_i;
       vector<double> proo_i;
       auto & gmat_delta_i = gmat_delta[i];
       auto & gmat_roo_i = gmat_roo[i];
       auto & dgmat_delta_i = dgmat_delta[i];
       auto & dgmat_roo_i = dgmat_roo[i];
       double ci_delta = 
          std::accumulate(
             gmat_delta_i.begin(), gmat_delta_i.end(), 
             1.0, std::multiplies<double>());
       double ci_roo = 
          std::accumulate(
             gmat_roo_i.begin(), gmat_roo_i.end(), 
             1.0, std::multiplies<double>());
       // loop for nodes in this reaction chain:
       for(unsigned jnode = 0; jnode < gmat_delta_i.size(); ++jnode) {
          // and so on ...
          //ci *= gmat_delta_i[jnode] * gmat_roo_i[jnode];
          double pd = dgmat_delta_i[jnode] * ci_roo;
          pd = std::accumulate(
                gmat_delta_i.begin(), gmat_delta_i.begin() + jnode,
                pd, std::multiplies<double>());
          pd = std::accumulate(
                gmat_delta_i.begin() + +jnode + 1, gmat_delta_i.end(), 
                pd, std::multiplies<double>());
          pd_i.push_back(pd);
          double proo = dgmat_roo_i[jnode] * ci_delta;
          proo = std::accumulate(
                gmat_roo_i.begin(), gmat_roo_i.begin() + jnode,
                proo, std::multiplies<double>());
          proo = std::accumulate(
                gmat_roo_i.begin() + +jnode + 1, gmat_roo_i.end(), 
                proo, std::multiplies<double>());
          proo_i.push_back(proo);
       }
       ci_2.push_back(ci_delta * ci_roo);
       pc_pdelta.push_back(pd_i);
       pc_proo.push_back(proo_i);
    }

    ci_2.push_back(1.0);
    double ctotal = std::accumulate(ci_2.begin(), ci_2.end(), 0.0);
    for(double & c : ci_2) c /= ctotal;

    // BY NOW pc_pdelta does not have pivot

    /*--- compute the positions and derivatives ---*/
    // NOTICE: DO NOT FORGET PIVOT WHEN CALCULATING DERIVATIVES
    
    /*** compute d ctotal / d r ***/

    // use mapping to save memory, i-th atom is the mapping[i]-th 
    // atom in reactionTree also the mapping[i]-th atom in 
    // pctotal_pr
    // mapped[i] indicates whether i-th atom is in the reactionTree
    // I do not use a vector of pair because 1) I do not want to call
    // pair constructor natom times 2) vector of bool saves memory a little
    vector<unsigned> mapping(natom);
    vector<bool> mapped(natom);
    // stores the index of atoms in reactionTree without duplication
    // i == mapping[reactionAtoms[i]]
    vector<unsigned> reactionAtoms;
    // count the atom in reactionTree
    // final count should equal to the number of true in mapped
    unsigned count = 0;
    // construct the reactionAtoms mapping and mapped:
    // essentially record the first appearance of all atoms in reactionTree
    for(auto ri =reactionTree.begin(); ri !=reactionTree.end(); ++ri) {
       for(auto r = ri->begin(); r != ri->end(); ++r) {
          unsigned i = (*r)[0], j = (*r)[1], h = (*r)[2];
          if(!mapped[oindexes[i]]) { // if oindexes[i] first appears
             mapped[oindexes[i]] = true; 
             mapping[oindexes[i]] = count;
             count ++; reactionAtoms.push_back(oindexes[i]);
          }
          if(!mapped[oindexes[j]]) { // if oindexes[j] first appears
             mapped[oindexes[j]] = true; 
             mapping[oindexes[j]] = count;
             count ++; reactionAtoms.push_back(oindexes[j]);
          }
          if(!mapped[hindexes[h]]) { // if hindexes[h] first appears
             mapped[hindexes[h]] = true; 
             mapping[hindexes[h]] = count;
             count ++; reactionAtoms.push_back(hindexes[h]);
          }
       }

#ifdef __NOISY__
if(comm.Get_rank()==0)
for(auto r : *ri) {
   printf("%d -> %d through %d | ",fullAtomList[oindexes[r[0]]].serial(),
         fullAtomList[oindexes[r[1]]].serial(),
         fullAtomList[hindexes[r[2]]].serial());
       }
if(comm.Get_rank()==0) printf("ci = %e (%e)\n",ci_2[ri-reactionTree.begin()],
      ci_2[ri-reactionTree.begin()] * ctotal);//@@@
#endif 

    }

#ifdef __NOISY__
if(comm.Get_rank()==0) printf("\n"); //@@@
#endif

    //\[pd] c / \[pd] \vec{r} 
    vector<vector<Vector>> pc_pr(
          gmat_delta.size(), vector<Vector>(reactionAtoms.size()));
    for(auto ri =reactionTree.begin(); ri !=reactionTree.end(); ++ri) {
       unsigned istate = ri - reactionTree.begin();
       auto & pci_pdelta = pc_pdelta[istate];
       auto & pci_proo = pc_proo[istate];
       for(auto r = ri->begin(); r != ri->end(); ++r) {
          unsigned jnode = r - ri->begin();
          double pci_pdeltaijh = pci_pdelta[jnode];
          double pci_prooij = pci_proo[jnode];
          unsigned i = (*r)[0], j = (*r)[1], h = (*r)[2];
          //d \vec{c} / d vec{r} = d \vec{c} / d delta * d delta / d \vec{r}
          //                     + d \vec{c} / d roo * d roo / d \vec{r}
          //
          //delta_ijh = rjh - rih
          //d delta_ijh / d vec{r}_i = - nih
          //d delta_ijh / d vec{r}_j =   njh
          //d delta_ijh / d vec{r}_h =   nih - njh
          //
          //d roo_ij / d vec{r}_i = nij
          //d roo_ij / d vec{r}_j = -nij

          // d \vec{c} / d delta * d delta / d \vec{r}:
          Vector nih = normDistOH(i, h);
          Vector njh = normDistOH(j, h);
          pc_pr[istate][mapping[oindexes[i]]] += (-nih) * pci_pdeltaijh;
          pc_pr[istate][mapping[oindexes[j]]] += njh * pci_pdeltaijh;
          pc_pr[istate][mapping[hindexes[h]]] += (nih - njh) * pci_pdeltaijh;
          // d \vec{c} / d roo * d roo / d \vec{r}:
          Vector nij = normDistOO(i, j);
          pc_pr[istate][mapping[oindexes[i]]] += nij * pci_prooij;
          pc_pr[istate][mapping[oindexes[j]]] += nij * (-pci_prooij);
          
          // for clarity, I moved this part outside this block
          //pctotal_pr[mapping[oindexes[i]]] += (-nih) * pci_pdeltaijh;
          //pctotal_pr[mapping[oindexes[j]]] += njh * pci_pdeltaijh;
          //pctotal_pr[mapping[hindexes[h]]] += (nih - njh) * pci_pdeltaijh;
       }
    }
    // the resulting pc_pr has dim of nstates x reactionAtoms.size() x 3
    //\[pd] ctotal / \[pd] \vec{r}_j before normed
    // mapping required to access this
    vector<Vector> pctotal_pr(reactionAtoms.size());
    for(const auto & pci_pr : pc_pr) {
       for(unsigned j = 0; j < reactionAtoms.size(); ++j) 
          pctotal_pr[j] += pci_pr[j];
    }
    // the resulting pctotal_pr has dim of reactionAtoms.size() x 3

    /*** compute COC and zero-th order of derivs on the fly ***/

    vector<Vector> r_coc(gmat_delta.size() + 1);
    for(auto r=reactionTree.begin(); r!=reactionTree.end(); ++r) {
       unsigned oi = r->back()[1];
       unsigned istate = r - reactionTree.begin();
       auto this_heavyList = heavyList[molid[oi]];

#ifdef __DOFITTING__
ci_2_on_atom[oi] += ci_2[istate]; // for fitting @@@
#endif

#ifdef __FAKE_VATOM__
for(auto molid_iter = molids_2b_watched.begin(); 
      molid_iter != molids_2b_watched.end(); molid_iter++) {
   if(molid[oi] == *molid_iter) {
      ci2_2b_watched[ static_cast<unsigned>(
            molid_iter - molids_2b_watched.begin()
            ) ] += ci_2[istate];
   }
}
#endif

       // temporarily push the proton to bondedHList[oi]
       bondedHList[oi].push_back(r->back()[2]);
       // temporarily erase the proton from the very parent
       unsigned op = (*r)[0][0];
       bondedHList[op].erase(
         std::find(bondedHList[op].begin(), bondedHList[op].end(), (*r)[0][2]));
       auto & t = r_coc[r - reactionTree.begin()];
       double qtot4check = 0.0;
#ifdef __PRINTF_DEBUG__
if(comm.Get_rank() == 0) {
   printf("COC atoms of mol %u:", molid[oi]);
   std::fflush(stdout);
}
#endif
       for(unsigned i : this_heavyList) {
         if(!bondedHList[i].empty()) {
           t += heavy_charges[i] * relpos[oindexes[i]];
           qtot4check += heavy_charges[i];
#ifdef __PRINTF_DEBUG__
if(comm.Get_rank() == 0) {
   printf("heavy atom %u ", getAbsoluteIndex(i).serial());
   std::fflush(stdout);
}
#endif
#ifndef __FAKE_VATOM__
           deriv[oindexes[i]]+=ci_2[istate]*heavy_charges[i]*Tensor::identity();
#endif
         }
         for(unsigned h : bondedHList[i]) {
            t += proton_charges[i] * relpos[hindexes[h]];
            qtot4check += proton_charges[i];
#ifdef __PRINTF_DEBUG__
if(comm.Get_rank() == 0) {
   printf("proton %u ", getAbsoluteIndex(h).serial());
   std::fflush(stdout);
}
#endif
#ifndef __FAKE_VATOM__
            deriv[hindexes[h]]
               += ci_2[istate]*proton_charges[i] * Tensor::identity();
#endif

//if(fullAtomList[hindexes[h]].serial()==980 && comm.Get_rank()==0) { // @@@
//   log.printf("coeff for 980 = %e x %e = %e\n", ci_2[istate], proton_charges[i], ci_2[istate]*proton_charges[i]);
//}

         }
       }
       if(fabs(qtot4check - 1.0) > 1e-4) {
          log.printf("WARNING: Total excess charge of mol %u = %lf at step %u\n",
                molid[oi], qtot4check, getStep());
       }
       //r_coc[r - reactionTree.begin()] = t;
       // erase the proton 
       bondedHList[oi].pop_back();
       // push the proton back
       bondedHList[op].push_back( (*r)[0][2] );
#ifdef __PRINTF_DEBUG__
if(comm.Get_rank() == 0) {
   printf("\n");
   std::fflush(stdout);
}
#endif
    }
    auto pivot_heavyList = heavyList[molid[pivot]];
    for(unsigned i : pivot_heavyList) {
       if(!bondedHList[i].empty()) {
         r_coc[gmat_delta.size()] += heavy_charges[i] * relpos[oindexes[i]];
#ifndef __FAKE_VATOM__
         deriv[oindexes[pivot]] += heavy_charges[i]/ctotal * Tensor::identity();
#endif
       }

#ifdef __PRINTF_DEBUG__
#ifndef __DEBUG_LOW__
printf("Now pivot_coc = %f %f %f\n", 
    r_coc[gmat_delta.size()][0], r_coc[gmat_delta.size()][1], r_coc[gmat_delta.size()][2]); // @@@
#endif
#endif

       for(unsigned h : bondedHList[i]) {
          r_coc[gmat_delta.size()] += proton_charges[i] * relpos[hindexes[h]];

#ifdef __PRINTF_DEBUG__
#ifndef __DEBUG_LOW__
printf("Now pivot_coc = %f %f %f\n", 
    r_coc[gmat_delta.size()][0], r_coc[gmat_delta.size()][1], r_coc[gmat_delta.size()][2]); // @@@
#endif
#endif

#ifndef __FAKE_VATOM__
          deriv[hindexes[h]] += proton_charges[i]/ctotal * Tensor::identity();
#endif
       }
    }
    for(auto coc_iter = r_coc.begin(); coc_iter != r_coc.end(); ++coc_iter) {

#ifdef __PRINTF_DEBUG__
if(comm.Get_rank()==0) {
printf("COC = %f %f %f ci_2 = %f\n",(*coc_iter)[0],(*coc_iter)[1],(*coc_iter)[2],
      ci_2[coc_iter - r_coc.begin()]); // @@@
}
#endif

       pos += (*coc_iter) * ci_2[coc_iter - r_coc.begin()];
    }

#ifdef __DOFITTING__
ci_2_on_atom[pivot] = 1 / ctotal; //@@@ for fitting
for(unsigned i = 0; i < nx; ++i) { // for fitting @@@
   if(ci_2_on_atom[i] > 0.0) 
      fs_out << fullAtomList[oindexes[i]].serial() << ' '
         << std::sqrt(ci_2_on_atom[i]) << ' ';
}
fs_out << '\n';
firsttime = true;
#endif

#ifdef __FAKE_VATOM__
for(auto molid_iter = molids_2b_watched.begin(); 
      molid_iter != molids_2b_watched.end(); molid_iter++) {
   if(molid[pivot] == *molid_iter) {
      ci2_2b_watched[ static_cast<unsigned>(
            molid_iter - molids_2b_watched.begin()
            ) ] += 1 / ctotal;
   }
}
#endif

    /*** compute d c/ctotal / d r on the fly and finialize the calculation ***/
    for(unsigned istate = 0; istate < gmat_delta.size(); istate++) {
       // d c/ctotal / d r =
       // ( d c / d r - (d ctotal / d r) * (c / ctotal) ) / ctotal
       for(unsigned ira = 0; ira < reactionAtoms.size(); ira++) {
#ifdef __FAKE_VATOM__

for(unsigned i_molid = 0; i_molid < molids_2b_watched.size(); i_molid++) {
   if(molid[reactionTree[istate].back()[1]] == molids_2b_watched[i_molid]) {
      for(unsigned alpha = 0; alpha < 3; ++alpha) {
         // d c^2_{i_molid} / d r_{ira, alpha}:
         deriv[reactionAtoms[ira]](alpha, i_molid) += 
             ( pc_pr[istate][ira][alpha] - pctotal_pr[ira][alpha] * ci_2[istate] ) / ctotal;
      }
   }
}

#else

          deriv[reactionAtoms[ira]] +=
#ifndef __TRANSPOSE__
             Tensor(
             ( pc_pr[istate][ira] - pctotal_pr[ira] * ci_2[istate] ) / ctotal,
             r_coc[istate] );
#else
             Tensor(
             r_coc[istate],
             ( pc_pr[istate][ira] - pctotal_pr[ira] * ci_2[istate] ) / ctotal);
#endif

#endif
       }
    }
    // ...  and the pivot state
    // pc_pr[gamt.size()] is not defined, otherwise this part can be
    // integrated into the above
    for(unsigned ira = 0; ira < reactionAtoms.size(); ira++) {
#ifdef __FAKE_VATOM__

for(unsigned i_molid = 0; i_molid < molids_2b_watched.size(); i_molid++) {
   if(molid[pivot] == molids_2b_watched[i_molid]) {
      for(unsigned alpha = 0; alpha < 3; ++alpha) {
         deriv[reactionAtoms[ira]](alpha, i_molid) += 
             - pctotal_pr[ira][alpha] / ctotal / ctotal;
      }
   }
}

#else

       deriv[reactionAtoms[ira]] += 
#ifndef __TRANSPOSE__
          Tensor(
                - pctotal_pr[ira] / ctotal / ctotal,
                r_coc[gmat_delta.size()] );
#else
          Tensor(
                r_coc[gmat.size()],
                - pctotal_pr[ira] / ctotal / ctotal);
#endif

#endif
    }

    // print out some fancy or not too fancy stuff
    if(print_delta) {
       // WARNING: this only makes sense for proton in water
       vector<double> deltas;
       for(const auto& r : reactionTree) {
          if(r.size() > 1) continue;
          if(r[0][0] == pivot) {
             // collect delta if the donor atom of 
             // the first reaction in chain r is the pivot 
             const auto& h = r[0][2];
             const auto& A = r[0][1];
             deltas.push_back(calculateDelta(pivot, A, h));
          } 
       }
       std::sort(deltas.begin(), deltas.end());
       log.printf("deltas: %ld", getStep());
       for(const auto& d : deltas) {
          log.printf(" %lf", d);
       }
       log.printf("\n");
    }
    if(print_shells) {
       vector<double> cec_1st(3, 0.0);
       vector<double> cec_2nd(3, 0.0);
       vector<double> cec_3rd(3, 0.0);
       unsigned nreactions = reactionTree.size();
       for(unsigned istate = 0; istate < nreactions; ++istate) {
          const auto& r = reactionTree[istate];
          double ci2 = ci_2[istate];
          const auto& coc = r_coc[istate];
          if(r.size() == 1) {
             for(int k = 0; k < 3; ++k) cec_1st[k] += ci2 * coc[k];
          } else if(r.size() == 2) {
             for(int k = 0; k < 3; ++k) cec_2nd[k] += ci2 * coc[k];
          } else if(r.size() == 3) {
             for(int k = 0; k < 3; ++k) cec_3rd[k] += ci2 * coc[k];
          } else {
             plumed_merror("PRINT_SHELLS not supporting fourth shell yet\n");
          }
       }
       log.printf("shells: %ld", getStep());
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", refpos[k]);
       }
       for(int k = 0; k < 3; ++k) {
          // pivot appears as the last one in ci_2 and r_coc
          // this is the H3O core, i.e. 0th shell
          log.printf(" %lf", ci_2[nreactions] * r_coc[nreactions][k]);
       }
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", cec_1st[k]);
       }
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", cec_2nd[k]);
       }
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", cec_3rd[k]);
       }
       log.printf("\n");
    }
    if(print_o_h) {
       // WARNING: by doing the following, I assume no reaction cycles
       // see computing r_coc above for more rigorous treatment by updating
       // bondedHList along the chain
       Vector cec_o(0.0, 0.0, 0.0);
       Vector cec_h(0.0, 0.0, 0.0);
       unsigned nreactions = reactionTree.size();
       for(unsigned istate = 0; istate < nreactions; ++istate) {
          const auto& r = reactionTree[istate];
          double ci2 = ci_2[istate];
          const auto& A = r.back()[1];
          const auto& H = r.back()[2];
          const auto& this_heavyList = heavyList[molid[A]];
          for(const auto& x : this_heavyList) {
             if(!bondedHList[x].empty()) {
                cec_o += ci2 * heavy_charges[x] * relpos[oindexes[x]];
                cec_h += ci2 * proton_charges[x] * relpos[hindexes[H]];
             }
             for(const auto& h : bondedHList[x]) {
                cec_h += ci2 * proton_charges[x] * relpos[hindexes[h]];
             }
          }
       } // end of loop for all acceptors
       // and also do for the H3O indeed
       const auto& pivot_heavyList = heavyList[molid[pivot]];
       for(unsigned x : pivot_heavyList) {
          if(!bondedHList[x].empty()) {
             cec_o += ci_2[nreactions] * heavy_charges[x] * relpos[oindexes[x]];
          }
          for(const auto& h : bondedHList[x]) {
             cec_h += ci_2[nreactions] * proton_charges[x] * relpos[hindexes[h]];
          }
       }
       log.printf("O_H: %ld", getStep());
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", refpos[k]);
       }
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", cec_o[k]);
       }
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", cec_h[k]);
       }
       log.printf("\n");
    }
    if(print_acceptors) {
       // WARNING: this only makes sense for proton in water
       plumed_assert(bondedHList[pivot].size() == 3);
       vector<double> rohs;
       for(int k = 0; k < 3; ++k) {
           rohs.push_back(std::sqrt(dist2OH(pivot, bondedHList[pivot][k])));
       }
       vector<pair<double,unsigned>> deltas_serials;
       for(unsigned i = 0; i < oindexes.size(); ++i) {
          if(i == pivot) continue;
          vector<double> deltas;
          for(int k = 0; k < 3; ++k) {
              deltas.push_back(
                 std::sqrt(dist2OH(i, bondedHList[pivot][k])) - rohs[k]);
          }
          deltas_serials.emplace_back(
                *std::min_element(deltas.begin(), deltas.end()), getSerial(i));
       } // end of loop of all oxygens
       // find the smallest three
       std::partial_sort( deltas_serials.begin(), 
                          deltas_serials.begin()+3, 
                          deltas_serials.end() );
       log.printf("acceptors: %ld", getStep());
       for(int k = 0; k < 3; ++k) {
          log.printf(" %ld", deltas_serials[k].second);
       }
       for(int k = 0; k < 3; ++k) {
          log.printf(" %lf", deltas_serials[k].first);
       }
       log.printf("\n");
    }



// check derivs @@@
#ifdef __NOISY__
log.printf("atoms in reactionTree: ");
for(unsigned i = 0; i < reactionAtoms.size(); ++i) {
   log.printf("%d ",fullAtomList[reactionAtoms[i]].serial());
}
log.printf("\n");
#endif
// end of checking derivs @@@
  
    /*--- print n leading states atom serials ---*/
    printNstates(print_nstates);

    /*--- update bondedHList and pivot ---*/
    // save bondedHList before updating
    bondedHList_ = bondedHList;

    //pop_back to get rid of pivot to adapt to the following code
    ci_2.pop_back();
    // get the largest ci_2 and reaction
    auto max_ci_iter = std::max_element(ci_2.begin(), ci_2.end());
#ifdef __SECONDLARGESTCI__
    if(*max_ci_iter < 1.0 / ctotal)  
    // if max smaller than currect pivot then max is the second largest
#else
    if(*max_ci_iter > 1.0 / ctotal) 
#endif
    {
      auto max_ci_reaction = reactionTree[max_ci_iter - ci_2.begin()];
      //plumed_assert(max_ci_reaction.size() == 1);
      if(max_ci_reaction.size() != 1) {
         log.printf("WARNING: max_ci_reaction.size() = %d at step %ld\n",
               max_ci_reaction.size(), getStep());
      }
      // push the H to the new o, erase the H in the old o
      for(auto one_reaction : max_ci_reaction) {
         pivot = one_reaction[0];
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
   printf("trying to delete hydrogen %u from heavy atom %u\n",
         one_reaction[2], pivot);
   if( std::find( bondedHList[pivot].begin(), 
   bondedHList[pivot].end(), one_reaction[2]) == bondedHList[pivot].end()) {
      printf("but hydrogen not found !!!!\n");
   }
   std::fflush(stdout);
}
#endif
         bondedHList[pivot].erase(
               std::find(
                  bondedHList[pivot].begin(), 
                  bondedHList[pivot].end(), one_reaction[2]));
                  //max_ci_reaction.back()[2]));
         //pivot = max_ci_reaction.back()[1];
         pivot = one_reaction[1];
         //bondedHList[pivot].push_back(max_ci_reaction.back()[2]);
         bondedHList[pivot].push_back(one_reaction[2]);
#ifdef __DEBUG_LOW__
if(comm.Get_rank() == 0) {
   printf("added hydrogen %u to heavy atom %u\n",
         one_reaction[2], pivot);
   std::fflush(stdout);
}
#endif
      }

//log.printf("step = %ld pivot = %d next pivot = %d through %d\n",
//      getStep(),
//      fullAtomList[oindexes[max_ci_reaction[0][0]]].serial(),
//      fullAtomList[oindexes[max_ci_reaction.back()[1]]].serial(),
//      fullAtomList[hindexes[max_ci_reaction.back()[2]]].serial());//@@@

      if(print_reaction) {
         log.printf("step = %ld ",getStep());
         for(auto r : max_ci_reaction) {
            log.printf("| %d -> %d through %d ",
                  fullAtomList[oindexes[r[0]]].serial(),
                  fullAtomList[oindexes[r[1]]].serial(),
                  fullAtomList[hindexes[r[2]]].serial());
                }
         log.printf("\n"); 
      }

    }
else {  
#ifdef __SECONDLARGESTCI__
   // max is larger than the currect pivot then
   // find the largest in {pivot}U{ci_2}/{max}
   double max_ci2 = 1.0 / ctotal;
   bool pt_happen = false;
   auto prev_max_ci_iter = max_ci_iter;
   auto max_ci_iter = ci_2.begin();
   for(auto it = ci_2.begin(); it != ci_2.end(); ++it) {
      if(it == prev_max_ci_iter) continue;
      if(*it > max_ci2) {
         pt_happen = true; 
         max_ci2 = *it; max_ci_iter = it;
      }
   }
   if(pt_happen) {
      auto max_ci_reaction = reactionTree[max_ci_iter - ci_2.begin()];
      if(max_ci_reaction.size() != 1) {
         log.printf("WARNING: max_ci_reaction.size() = %d at step %ld\n",
               max_ci_reaction.size(), getStep());
      }
      // push the H to the new o, erase the H in the old o
      for(auto one_reaction : max_ci_reaction) {
         pivot = oone_reaction[0];
         bondedHList[pivot].erase(
               std::find(
                  bondedHList[pivot].begin(), 
                  bondedHList[pivot].end(), one_reaction[2]));
                  //max_ci_reaction.back()[2]));
         //pivot = max_ci_reaction.back()[1];
         pivot = one_reaction[1];
         //bondedHList[pivot].push_back(max_ci_reaction.back()[2]);
         bondedHList[pivot].push_back(one_reaction[2]);
      }

//log.printf("step = %ld pivot = %d next pivot = %d through %d\n",
//      getStep(),
//      fullAtomList[oindexes[max_ci_reaction[0][0]]].serial(),
//      fullAtomList[oindexes[max_ci_reaction.back()[1]]].serial(),
//      fullAtomList[hindexes[max_ci_reaction.back()[2]]].serial());//@@@
log.printf("step = %ld ",getStep());
for(auto r : max_ci_reaction) {
   log.printf("| %d -> %d through %d ", 
         fullAtomList[oindexes[r[0]]].serial(),
         fullAtomList[oindexes[r[1]]].serial(),
         fullAtomList[hindexes[r[2]]].serial());
       }
log.printf("\n"); //@@@

   }

#ifdef __NOISY__
   else { //@@@
//log.printf("step = %ld pivot = %d next pivot = %d\n",
//      getStep(),
//      fullAtomList[oindexes[pivot]].serial(),
//      fullAtomList[oindexes[pivot]].serial());  //@@@
log.printf("step = %ld no reaction\n", getStep()); //@@@
   }
#endif

#else

#ifdef __NOISY__
//log.printf("step = %ld pivot = %d next pivot = %d\n",
//      getStep(),
//      fullAtomList[oindexes[pivot]].serial(),
//      fullAtomList[oindexes[pivot]].serial());  //@@@
log.printf("step = %ld no reaction\n", getStep()); //@@@
#endif

#endif
}


#ifdef __FAKE_VATOM__
    setPosition(ci2_2b_watched);
#else
    setPosition(pos);
#endif
    // I need this because I pop_back() before updating bondedHList
    ci_2.push_back(1.0 / ctotal);
    setMass(1.0);
    setCharge(0.0);
    setAtomsDerivatives(deriv);

#ifdef __DOFITTING__
bondedHList = vector<vector<unsigned>>(nx); // for fitting @@@
#endif

}

}
}
#ifdef __DOFITTING__
#undef __DOFITTING__
#endif
#ifdef __SECONDLARGESTCI__
#undef __SECONDLARGESTCI__
#endif
#ifdef __NOISY__
#undef __NOISY__
#endif
#ifdef __WOPAIR__
#undef __WOPAIR__
#endif
#ifdef __TRANSPOSE__
#undef __TRANSPOSE__
#endif
#ifdef __FAKE_VATOM__
#undef __FAKE_VATOM__
#endif
#ifdef __DEBUG__
#undef __DEBUG__
#endif
#ifdef __DEBUG_LOW__
#undef __DEBUG_LOW__
#endif
