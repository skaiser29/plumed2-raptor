#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "MSTAnalysisBase.h"

namespace PLMD {
namespace adjmat {

void MSTAnalysisBase::registerKeywords( Keywords &keys )
{
   MultiColvarBase::registerKeywords( keys );
   keys.add("compulsory", "MST", "the action that builds minimum spanning tree");
}

MSTAnalysisBase::MSTAnalysisBase(const ActionOptions& ao):
   Action(ao),
   MultiColvarBase(ao),
   mstbase(NULL)
{
   //TODO what is the point of this?
   matsums = true; usespecies = true;
   std::vector<AtomNumber> fatoms;
   if(!parseMultiColvarAtomList("MST",-1,fatoms)) error("unable to read MST");
   if(mybasemulticolvars.size()!=1) error("should be only one MST inputted");
   mstbase = dynamic_cast<MSTBase*>(mybasemulticolvars[0]);
   if(!mstbase) error("input MST is not a minimum spanning tree");
}


}
}

#endif
