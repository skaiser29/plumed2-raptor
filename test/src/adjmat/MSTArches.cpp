#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "MSTAnalysisBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

class MSTArches : public MSTAnalysisBase {
   //max number of components
   unsigned max_num_comps;
public:
   static void registerKeywords(Keywords& keys);
   explicit MSTArches(const ActionOptions&);
   void calculate();
   //double compute(const unsigned&, multicolvar::AtomValuePack&) const;
};

PLUMED_REGISTER_ACTION(MSTArches,"MSTARCHES")

void MSTArches::registerKeywords(Keywords& keys)
{
   MSTAnalysisBase::registerKeywords(keys);
   componentsAreNotOptional(keys);
   keys.addOutputComponent("edge","default","edge length in the MST");
   keys.add("optional","MAXNUM","max number of components");
}

MSTArches::MSTArches(const ActionOptions&ao):
   Action(ao),
   MSTAnalysisBase(ao)
{
   max_num_comps = mstbase->getNumberOfNodes() - 1;
   unsigned tmp_max_num_comps = 0;
   parse("MAXNUM",tmp_max_num_comps);
   if(tmp_max_num_comps !=0 && tmp_max_num_comps < max_num_comps) 
      max_num_comps = tmp_max_num_comps;
   for(unsigned i = 0; i < max_num_comps; ++i) {
      std::string num; Tools::convert(i + 1, num);
      addComponentWithDerivatives("edge-"+num);
      componentIsNotPeriodic("edge-"+num);
      getPntrToComponent(i)
         ->resizeDerivatives(mstbase->getNumberOfDerivatives());
   }
}

//double MSTArches::compute(
//      const unsigned& tind, multicolvar::AtomValuePack& myatoms) const
void MSTArches::calculate()
{
   std::vector<double> arches = mstbase->retrieveMSTArchLengths();

printf("arches has a size of %ld\n", arches.size()); //@@@

   for(unsigned i = 0; i < max_num_comps; ++i) {
      getPntrToComponent(i)->set(arches[i]);
   }

   //TODO derivatives


}

}
}


#endif
