/*
 * author: Chenghan Li (lch004218@gmail.com)
 *
 */
#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "KSPBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

void KSPBase::registerKeywords(Keywords& keys) {
   ActionWithInputMatrix::registerKeywords(keys);
   keys.addFlag("NOINVERSE",false,"Do not take the inverse of the weights");
   keys.addFlag("NOSYMMETRIZE",false,"Do not symmetrize the input matrix");
   //keys.addFlag("NOCHECK",false,"Do not check the (k+1)-th path");
   //keys.addFlag("NOWARNINGS",false,"Do not give warnings");
   //keys.addFlag("NOHALTING",false,"Do not halt program")
   keys.add("compulsory","K","number of shortest path");
   keys.add("atoms","SOURCE","source of KSP");
   keys.add("atoms","TARGET","target of KSP");
}

KSPBase::KSPBase(const ActionOptions&ao):
   Action(ao),
   ActionWithInputMatrix(ao),
   inverse_weight(true),
   symmetrize(true),
   npaths(0)
{
   usespecies = false;
   parse("K",npaths);
   if(npaths == 0) error("should specify K");
   log.printf("  will compute the first %u paths\n",npaths);

   bool noinverse_weight, nosymmetrize;//, nowarnings;
   //bool nohalting, nocheck;
   parseFlag("NOINVERSE", noinverse_weight);
   parseFlag("NOSYMMETRIZE", nosymmetrize);
   //parseFlag("NOWARNINGS", nowarnings);
   //parseFlag("NOHALTING", nohalting);
   //parseFlag("NOCHECK", nocheck);
   inverse_weight = !noinverse_weight;
   symmetrize = !nosymmetrize;
   //warnings = !nowarnings;
   //halting = !nohalting; check = !nocheck;
   //if(compute_all) check = false;
   //if(!check) { warnings = false; halting = false; }
   if(inverse_weight) log.printf("  will use the inverse of weights\n");
   else  log.printf("  does not use the inverse\n");
   if(symmetrize) log.printf("  will symmetrize the adjacency matrix\n");
   else  log.printf("  does not symmetrize the adjacency matrix\n");
   //if(check) log.printf("  will check the (k+1)-th path\n");
   //else log.printf("  does not check the (k+1)-th path\n");
   //if(warnings) log.printf("  will give warnings if ");
   //else log.printf("");
   parseAtomList("SOURCE", start);
   parseAtomList("TARGET", end);
   if(start.empty()) error("no atoms provided for SOURCE");
   if(end.empty()) error("no atoms provided for TARGET");
   if( getAdjacencyVessel() ){
      //if( getNumberOfNodeTypes()!=1 ) 
      //   error("should only be running KSP with one base multicolvar in function");
      if( !getAdjacencyVessel()->isSymmetric() and !symmetrize) 
         error("asymmetric graph needs symmetrization");
      if( getAdjacencyVessel()->isSymmetric() and symmetrize) 
         log.printf("  WARNING: symmetric graph does not need symmetrization\n");
   }
   if( keywords.exists("MATRIX") ){
      std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
   }
}

void KSPBase::calculate() 
{
   //use boost to get KSP
   //ksp_nedges, ksp, ksp_edge_list, ksp_arch_lengths are built
   buildKSP();
}

}
}

#endif
