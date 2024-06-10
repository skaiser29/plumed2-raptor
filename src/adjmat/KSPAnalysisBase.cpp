#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "KSPAnalysisBase.h"

namespace PLMD {
namespace adjmat {

void KSPAnalysisBase::registerKeywords( Keywords &keys )
{
   MultiColvarBase::registerKeywords( keys );
   keys.add("compulsory", "KSP", "the action that computes k-shortest paths");
}

KSPAnalysisBase::KSPAnalysisBase(const ActionOptions& ao):
   Action(ao),
   MultiColvarBase(ao),
   kspbase(NULL)
{
   //TODO what is the point of this?
   matsums = true; usespecies = true;
   std::vector<AtomNumber> fatoms;
   if(!parseMultiColvarAtomList("KSP",-1,fatoms)) error("unable to read KSP");
   if(mybasemulticolvars.size()!=1) error("should be only one KSP inputted");
   kspbase = dynamic_cast<KSPBase*>(mybasemulticolvars[0]);
   if(!kspbase) error("input KSP is not a k-shortest path object");
}


}
}

#endif
