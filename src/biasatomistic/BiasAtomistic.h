#ifndef __PLUMED_biasatomistic_BiasAtomistic_h
#define __PLUMED_biasatomistic_BiasAtomistic_h

#include "bias/Bias.h"
#include "core/ActionAtomistic.h"

namespace PLMD {
namespace biasatomistic {

class BiasAtomistic :
   public ActionAtomistic,
   public bias::Bias
{
   // for fixing the confusion from multiple inheritance
   void lockRequests();
   void unlockRequests();
   void calculateNumericalDerivatives( ActionWithValue* a = NULL);

public:
   explicit BiasAtomistic(const ActionOptions&);
   static void registerKeywords(Keywords&) {}
};

inline void BiasAtomistic::lockRequests() {
   ActionWithArguments::lockRequests();
   ActionAtomistic::lockRequests();
}

inline void BiasAtomistic::unlockRequests() {
   ActionWithArguments::unlockRequests();
   ActionAtomistic::unlockRequests();
}

inline void BiasAtomistic::calculateNumericalDerivatives(ActionWithValue* a) {
   (void) a;
   plumed_merror(
         "BiasAtomistic::calculateNumericalDerivatives is not implemented"
         );
}

}
}

#endif
