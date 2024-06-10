/*
 *    author: Chenghan Li
 */
#include "BiasAtomistic.h"
#include "core/ActionSet.h"
#include "tools/SwitchingFunction.h"
#include "tools/File.h"

#include "vatom/RMDCEC.h"

// I already used std in RMDCEC.h
//using namespace std;

namespace PLMD {
namespace biasatomistic {

//#define __TEST_EDS__
//#define __DEBUG__

class WeightedEDS : public BiasAtomistic
{
   vatom::RMDCEC* cec;
   SwitchingFunction sf;

   vector<int> rpowers;
   vector<double> centers, couplings;
   vector<Value*> moment_pts;
   IFile eds_in;
   
   // TODO cell list or at least neighor list

   unsigned natomh;
/// j-th atom in this action is the indexes[j]-th atom in cec
/// size == number of atoms here
   vector<unsigned> indexes_;

/// update bondedHList given a reaction_chain
   template<class T>
   void forwardReaction(
     vector<vector<unsigned>>& bondedHList, T reaction_chain);
/// recover bondedHList to the state before appling reaction_chain
   template<class T>
   void backwardReaction(
     vector<vector<unsigned>>& bondedHList,T reaction_chain);
public:
   explicit WeightedEDS(const ActionOptions&);
   void calculate();
   static void registerKeywords(Keywords& keys);

};

PLUMED_REGISTER_ACTION(WeightedEDS,"WEDS")

void WeightedEDS::registerKeywords(Keywords& keys) {
   Bias::registerKeywords(keys);
   // do not use ARG otherwise ActionWithArguments::parseArgument will
   // be called in its constructor, then we'll lose the opportunity to 
   // get a pointer to RMDCEC
   keys.add("compulsory", "RMDCEC", "the input RMDCEC action");
   keys.add("optional", "SWITCH", "the switching function used in EDS");
   keys.add("optional", "NN", "nn of the switching function used in EDS");
   keys.add("optional", "MM", "mm of the switching function used in EDS");
   keys.add("optional", "R_0", "r_0 of the switching function used in EDS");
   keys.add("optional", "D_0", "d_0 of the switching function used in EDS");
   keys.add("compulsory", "RPOWERS", "moment numbers in EDS");
   keys.add("compulsory", "CENTERS", "reference centers in EDS");
   keys.add("optional", "COUPLINGS", "couplings in EDS");
   keys.add("optional", "IN_RESTART", 
         "EDS restart file; ignored if COUPLINGS is given");
   keys.add("optional", "ARGS", 
         "cv names of coordination number in EDS input file");
   keys.add("atoms", "ATOMH", "input hydrogens");
   keys.add("atoms", "ATOMX", "input heavy atoms");
}


WeightedEDS::WeightedEDS (const ActionOptions& ao) : 
   Action(ao),
   BiasAtomistic(ao)
{
   vector<AtomNumber> atomh, atomx;
   parseAtomList("ATOMH", atomh);
   if(atomh.empty()) plumed_merror("ATOMH not given");
   parseAtomList("ATOMX", atomx);
   if(atomx.empty()) plumed_merror("ATOMX not given");

   parseVector("RPOWERS", rpowers);
   if(rpowers.empty()) plumed_merror("RPOWERS not given");
   parseVector("CENTERS", centers);
   if(centers.empty()) plumed_merror("CENTERS not given");
   if(rpowers.size() != centers.size())
      plumed_merror("RPOWERS and CENTERS do not have same size");
   parseVector("COUPLINGS", couplings);
   if(couplings.empty()) {
      string eds_fn;
      parse("IN_RESTART", eds_fn);
      if(eds_fn.empty()) plumed_merror("COUPLINGS nor IN_RESTART not given");
      vector<string> cvnames;
      parseVector("ARGS", cvnames);
      if(cvnames.empty()) plumed_merror("ARGS not given");
      if(rpowers.size() != cvnames.size())
         plumed_merror("RPOWERS and ARGS do not have same size");

      // things I do not use but have to read
      eds_in.open(eds_fn);
      double junk;
      if(eds_in.FieldExist("kbt")) {
        eds_in.scanField("kbt", junk);
      } else { error("No field 'kbt' in restart file"); }
      if(eds_in.FieldExist("update_period")) {
        eds_in.scanField("update_period", junk);
      } else { error("No field 'update_period' in restart file"); }
      if(eds_in.FieldExist("adaptive")) {
        //note, no version of scanField for boolean
        eds_in.scanField("adaptive", junk);
      } else { error("No field 'adaptive' in restart file"); }
      if(eds_in.FieldExist("seed")) {
        eds_in.scanField("seed", junk);
      } else { error("No field 'seed' in restart file"); }

      // read couplings from eds restart file
      couplings.resize(centers.size());
      double time;
      while(eds_in.scanField("time", time)) {
         for(unsigned k = 0; k < rpowers.size(); ++k) {
            const auto& cv_name = cvnames[k];
            double junk;
            eds_in.scanField(cv_name + "_center", junk);
            eds_in.scanField(cv_name + "_set", junk);
            eds_in.scanField(cv_name + "_target", junk);
            eds_in.scanField(cv_name + "_coupling", couplings[k]);
            eds_in.scanField(cv_name + "_maxrange", junk);
            eds_in.scanField(cv_name + "_maxgrad", junk);
            eds_in.scanField(cv_name + "_accum", junk);
            eds_in.scanField(cv_name + "_mean", junk);
            eds_in.scanField(cv_name + "_std", junk);
         }
         eds_in.scanField();
      }
      eds_in.close();
   }
   if(rpowers.size() != couplings.size())
      plumed_merror("RPOWERS and COUPLINGS do not have same size");

   vector<string> rmdcec_lab;
   parseVector("RMDCEC", rmdcec_lab);
   if(rmdcec_lab.empty()) plumed_merror("RMDCEC not given");
   if(rmdcec_lab.size() != 1) plumed_merror("more than 1 RMDCEC is given");
   cec = plumed.getActionSet().selectWithLabel<vatom::RMDCEC*>(rmdcec_lab[0]);
   if(!cec) plumed_merror("input RMDCEC is not a RMDCEC action");


   // check if atomh and atomx belong to RMDCEC atoms
   // and construct indexes_;
   const auto & cec_atoms = cec->getAbsoluteIndexes();
   natomh = atomh.size();
   // now atomh has all atoms
   atomh.insert(atomh.end(), atomx.begin(), atomx.end());
#ifndef __TEST_EDS__
   for(unsigned i =0 ; i < atomh.size(); ++i) {
      const auto& a = atomh[i];
      const auto& p = std::find(cec_atoms.begin(), cec_atoms.end(), a);
      if(p == cec_atoms.end()) {
         plumed_merror(
           "atom " + to_string(a.serial()) + " cannot be found in RMDCEC");
      }
      indexes_.push_back(p - cec_atoms.begin());
   }
#endif

   requestAtoms(atomh);
   // add dependency of RMDCEC
   // this cannot be done before requestAtoms because requestAtoms clears
   // all dependencies
   addDependency(cec);

   // read sf, exactly the same as COORDINATIONNUMBER
   string sw, errors; parse("SWITCH",sw);
   if(sw.length() > 0) {
      sf.set(sw, errors);
      if(errors.length() != 0) error("problem reading SWITCH keyword : " + errors);
   } else { 
      double r_0=-1.0, d_0; int nn, mm = 0;
      parse("NN", nn); parse("MM", mm);
      parse("R_0", r_0); parse("D_0", d_0);
      if(r_0 < 0.0) error("you must set a value for R_0");
      sf.set(nn, mm, r_0, d_0);
   }

   checkRead();

   // TODO logs
   log.printf(
     "  coordination of central atom and those within %s\n",
     (sf.description()).c_str() );
   log.printf("  moment orders:\n  ");
   for(unsigned k = 0; k < couplings.size(); ++k)
      log.printf("%u ", rpowers[k]);
   log.printf("\n");
   log.printf("  coupling constants:\n  ");
   for(unsigned k = 0; k < couplings.size(); ++k)
      log.printf("%lf ", couplings[k]);
   log.printf("\n");
   log.printf("  reference moment values:\n  ");
   for(unsigned k = 0; k < couplings.size(); ++k)
      log.printf("%lf ", centers[k]);
   log.printf("\n");

   for(unsigned i : rpowers) {
      addComponent("moment-" + to_string(i));
      componentIsNotPeriodic("moment-" + to_string(i));
      moment_pts.push_back(getPntrToComponent("moment-" + to_string(i)));
   }
}

template<class T>
void WeightedEDS::forwardReaction(
  vector<vector<unsigned>>& bondedHList, T reaction_chain) {

   unsigned pivot = reaction_chain[0][0];
   for(const auto & one_reaction : reaction_chain) {
      bondedHList[pivot].erase(
        std::find( bondedHList[pivot].begin(), 
                   bondedHList[pivot].end(), 
                   one_reaction[2] ) );
      pivot = one_reaction[1];
      bondedHList[pivot].push_back(one_reaction[2]);
   }

}

template<class T>
void WeightedEDS::backwardReaction(
  vector<vector<unsigned>>& bondedHList, T reaction_chain) {

   unsigned pivot = reaction_chain.back()[1];
   for(auto one_reaction = reaction_chain.rbegin();
            one_reaction != reaction_chain.rend();
            one_reaction++ ) {
      bondedHList[pivot].erase(
        std::find( bondedHList[pivot].begin(),
                   bondedHList[pivot].end(),
                   (*one_reaction)[2] ) );
      pivot = (*one_reaction)[0];
      bondedHList[pivot].push_back((*one_reaction)[2]);
   }
}

void WeightedEDS::calculate() {

   auto& forces = modifyForces();

   double ene = 0;

#ifdef __TEST_EDS__
   // let me first check if I can reproduce EDS O-H bias force
   vector<double> moments(centers.size(), 0.0);
   for(unsigned i = 0; i < natomh; ++i)
      for(unsigned j = natomh; j < getNumberOfAtoms(); ++j) {
         // rj - ri
         Vector vij = pbcDistance(getPosition(i), getPosition(j));
         double rij = vij.modulo();
         if(rij < 1.2) continue;
         if(rij > 6.8940444834) continue;
         vij /= rij;
         double dfsw, fsw = sf.calculate(rij, dfsw);
         dfsw *= rij;
         for(unsigned k = 0; k < centers.size(); k++) {
            double power = 1.0, dpower = 0.0;
            if(rpowers[k] > 0)  {
               power = pow(rij, rpowers[k]);
               dpower = rpowers[k] * pow(rij, rpowers[k] - 1);
            }
            ene += couplings[k] * power * fsw;
            moments[k] += power * fsw;
            double de_dr = couplings[k] * (power * dfsw + dpower * fsw);
            forces[i] += vij * de_dr;
            forces[j] -= vij * de_dr;
         }
      }

   ene /= natomh;
   for(auto& f : forces) f /= natomh;
   for(unsigned k = 0; k < centers.size(); ++k) {
      ene -= couplings[k] * centers[k];
      moments[k] /= natomh;
   }

log.printf("step=%u 0-th moment=%lf 1-st moment=%lf\n",
      getStep(), moments[0], moments[1]
      );

#else

   // retrieve the bondedHList in RMDCEC before updating
   auto& bondedHList = cec->getBondedHList();
   const auto& reactionTree = cec->getReactionTree();
//   const auto& hindexes = cec->getHIndexes();
//   const auto& oindexes = cec->getOIndexes();
   const auto& ci_2 = cec->getCi2();
   const unsigned natomh_cec = cec->getNAtomH();

   unsigned rank = comm.Get_rank();
   unsigned stride = comm.Get_size();

   vector<double> moments(centers.size(), 0.0);
//   for(unsigned ireact = 0; ireact < ci_2.size(); ++ireact) {
   for(unsigned ireact = rank; ireact < ci_2.size(); ireact += stride) {

      // I put the pivot ci_2 at the end of ci_2
      // so no updating bondedHList is needed for that
      if(ireact != ci_2.size() - 1)
         forwardReaction(bondedHList, reactionTree[ireact]);

      double ci2 = ci_2.at(ireact);
#ifdef __DEBUG__
      printf("state %u ci2 = %lf\n", ireact, ci2);
#endif

      for(unsigned j = natomh; j < getNumberOfAtoms(); ++j) {
         // atom j is the jj-th oxygen in RMDCEC
         unsigned jj = indexes_[j] - natomh_cec;
         for(unsigned i = 0; i < natomh; ++i) {
            // atom i is the ii-th hydrogen in RMDCEC
            unsigned ii = indexes_[i];
            // skip if i is bonded to j
#ifdef __DEBUG__
            if( std::find(
                  bondedHList[jj].begin(), bondedHList[jj].end(), ii) 
                != bondedHList[jj].end()) {
               printf("state %u skip %u %u\n", 
                 ireact,
                 getAbsoluteIndex(j).serial(),
                 getAbsoluteIndex(i).serial() );
               continue;
            }
#else
            if( std::find(
                  bondedHList[jj].begin(), bondedHList[jj].end(), ii) 
                != bondedHList[jj].end()) continue;
#endif
            // rj - ri
            Vector vij = pbcDistance(getPosition(i), getPosition(j));
            double rij = vij.modulo();
            vij /= rij;
            double dfsw, fsw = sf.calculate(rij, dfsw);
            dfsw *= rij;
            for(unsigned k = 0; k < centers.size(); k++) {
               double power = 1.0, dpower = 0.0;
               if(rpowers[k] > 0)  {
                  power = pow(rij, rpowers[k]);
                  dpower = rpowers[k] * pow(rij, rpowers[k] - 1);
               }
               ene += couplings[k] * power * fsw * ci2;
               moments[k] += power * fsw * ci2;
               double de_dr = couplings[k] * (power * dfsw + dpower * fsw);
               forces[i] += vij * de_dr * ci2;
               forces[j] -= vij * de_dr * ci2;
            } // end of lopp moments
         } // end of loop hydrogen
      } // end of loop oxygen

      // recover the bondedHList
      if(ireact != ci_2.size() - 1)
         backwardReaction(bondedHList, reactionTree[ireact]);
   } // end of loop states

   comm.Sum(moments);
   comm.Sum(forces);
   comm.Sum(ene);

   ene /= natomh;
   for(auto& f : forces) f /= natomh;
   for(unsigned k = 0; k < centers.size(); ++k) {
      ene -= couplings[k] * centers[k];
      moments[k] /= natomh;
   }

#endif

   setBias(ene);
   for(unsigned k = 0; k < moment_pts.size(); ++k) {
      moment_pts[k]->set(moments[k]);
   }
}

}   
}
