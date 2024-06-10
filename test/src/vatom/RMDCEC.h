#ifndef __PLUMED_vatom_RMDCEC_h
#define __PLUMED_vatom_RMDCEC_h
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/SwitchingFunction.h"

#include <numeric>
#include <tuple>
#include <limits>
#include <unordered_map>
#include <queue>

#include "tools/File.h" //@@@ for fitting
#include <fstream>
#include <sstream>

// set vatom coordinates to atom's ci^2 by referring atom index
//#define __FAKE_VATOM__
//#define __DEBUG__
//#define __DEBUG_LOW__
//#define __DEBUG_07312019__
//
//#define __DOFITTING__
//#define __SECONDLARGESTCI__
//#define __NOISY__
//#define __WOPAIR__
//#define __PRINTF_DEBUG__
//#define __TRANSPOSE__

//TODO USE TRUE DFS !!!
//TODO use EXTRA_HEAVY0 to indicate the atoms in COC but not reactive in mol0
//use EXTRA_CHARGES0 to indicate the charges of EXTRA_HEAVY0
//TODO use makewhole to treat pbc
//TODO use REFMOL to use bonded atom in that mol as refpos
//TODO use ABSCHARGE to use abs value of charges

using namespace std;

namespace PLMD {
namespace vatom{

class RMDCEC : public ActionWithVirtualAtom {
//private:
/// Max number of special mols
    static const unsigned MAXNUMSPECIAL = 100;
/// Max value of dmax2 over all sf
    double max_dmax2;
/// Whether to use pbc condition
    bool do_pbc;
/// Whether to use masks to accelerate computation
    bool use_mask;
/// May not be useful 
    unsigned natomh,natomnp,natomp,natom; 
/// Cutoff distance square for oxygens
    double rlisto;
/// Construct these two vectors for historical reasons
/// List of indexes of hydrogen atoms
    vector<unsigned> hindexes;
/// List of indexes of heavy atoms
    vector<unsigned> oindexes;
/// Index of donor in oindexes
    unsigned pivot; 
/// Number of protons acceptable for AtomX; should be 3 for water
    vector<int> wx; 
/// BondedHList for each heavy atom; changing
    vector<vector<unsigned>> bondedHList; 
/// BondedHList for each heavy atom before updating
    vector<vector<unsigned>> bondedHList_; 
/// Bonded heavy atom for each hydrogen atom; changing
    vector<long int> bondedOList; 
/// Molecule id for each heavy atom; not changing
    vector<unsigned> molid; 
/// Heavy atom list for each molecule; not changing
    vector<vector<unsigned>> heavyList; 
/// Max number of protons each molecule has when no excess charge present
    vector<int> max_proton;
///// Number of protons each molecule has currently; changing
//    vector<int> num_proton;
/// Reaction tree
    vector<vector<array<unsigned,3>>> reactionTree; 
/// The weight matrix of delta part
    vector<vector<double>> gmat_delta;
/// The derivs of weight matrix of delta part
    vector<vector<double>> dgmat_delta;
/// The weight matrix of rOO part
    vector<vector<double>> gmat_roo;
/// The derivs of weight matrix of roo part
    vector<vector<double>> dgmat_roo;
///// The weight for g, g^2 and g^3
//    array<double,3> A;
/// Stores the relative positions compared to refpos
    vector<Vector> relpos;
/// Gives the i-th heavy atom's charge when bonded
    vector<double> heavy_charges;
/// Gives the charge of proton when bonding to i-th heavy atom
    vector<double> proton_charges;
/// Whether to do the fourth shell
    bool do_fourth;
/// Stores REFMOL
    unsigned refmol;
    //SwitchingFunction sf;
///// Switching function index to be used
//    Matrix<unsigned> sf_indexes;
/// Map from special atom index to sf_delta index in sf_delta_library
    unordered_map<long, unsigned> map2sf_delta_index;
/// Map from special atom index to sf_roo index in sf_roo_library
    unordered_map<long, unsigned> map2sf_roo_index;
/// Return the switching function for delta used by atom (i, j)
    SwitchingFunction sfij_delta(unsigned i, unsigned j) const;
/// Return the switching function for r_OiOj
    SwitchingFunction sfij_roo(unsigned i, unsigned j) const;
/// All sf used in this action
    vector<SwitchingFunction> sf_delta_library, sf_roo_library;
    bool firsttime;
    Vector initPos();
    unsigned getTopo();
    void updateList();
/// DFS state search
    void buildReactionTree();
/// Calculate normed vector (r_oindexes[i] - r_hindexes[j]) / norm
    void normDistOH(unsigned i, unsigned j, Vector& result) const;
    Vector normDistOH(unsigned i, unsigned j) const;
/// Calculate normed vector (r_oindexes[i] - r_oindexes[j]) / norm
    Vector normDistOO(unsigned i, unsigned j) const;
/// Calculate distance^2 between oindexes[i] and hindexes[j]
    double dist2OH(unsigned i, unsigned j) const;
/// Calculate distance^2 between oindexes[i] and oindexes[j]
    double dist2OO(unsigned i, unsigned j) const;
    double calculateDelta(unsigned i, unsigned j, unsigned k) const;
/// Check whether √d2 - √d1 > √d0 
    bool shouldSkip(double d1, double d2, double d0) const;
/// Calculate the distance between hydrogen and mol
    tuple<double,unsigned,unsigned> dist2MolH(unsigned molid, unsigned h) const;
/// Ratio of max water O-H bond length over dmax of sf
    double ratio;
/// File stream for initial topo if provided
    ifstream fs_topo;
/// Whether to guess or read initial topo
    bool guess_initial_topo;
/// Read initial topology
    unsigned readTopo();
/// Print reaction
    bool print_reaction; 
/// Print out three delta values
    bool print_delta;
/// Print out decomposition of CEC into solvation shells
    bool print_shells;
/// Print out decomposition of CEC into H part and O part
    bool print_o_h;
/// Print out three most probable acceptors
    bool print_acceptors;
/// Number of states to print out with leading c_i^2
    unsigned print_nstates;
/// Function to do the above
    void printNstates(unsigned);
/// ci^2 to be passed to other actions
    vector<double> ci_2;

#ifdef __FAKE_VATOM__
/// MOLIDS whose ci^2 are to be set in resulting pos
    vector<unsigned> molids_2b_watched;
/// ci^2 of molecules in molids_2b_watched
    Vector ci2_2b_watched;
#endif

#ifdef __DOFITTING__
~RMDCEC(); // for fitting @@@
ofstream fs_out;// output ci's for fitting @@@
//ofstream fs_delta;// output delta for fitting @@@
#endif

//@@@
unsigned getSerial(unsigned i) const;
vector<AtomNumber> fullAtomList;
unsigned h_serial, o_serial;

public:
    explicit RMDCEC(const ActionOptions&ao);
    void calculate();
    static void registerKeywords( Keywords& keys );
    // API needed for WeightedEDS
    const std::vector<AtomNumber> & getAbsoluteIndexes() const {
       return ActionAtomistic::getAbsoluteIndexes();
    }
    AtomNumber getAbsoluteIndex(int i) const {
       return ActionAtomistic::getAbsoluteIndex(i);
    }
    const vector<vector<array<unsigned,3>>> & getReactionTree() const {
       return reactionTree;
    }
    const vector<unsigned>& getHIndexes() const {
       return hindexes;
    }
    const vector<unsigned>& getOIndexes() const {
       return oindexes;
    }
    vector<vector<unsigned>>& getBondedHList() {
       return bondedHList_;
    }
    const vector<double>& getCi2() const {
       return ci_2;
    }
    unsigned getNAtomH() const {
       return natomh;
    }
};

}
}

#endif
