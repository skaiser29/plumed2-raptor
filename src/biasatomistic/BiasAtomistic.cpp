#include "BiasAtomistic.h"

namespace PLMD {
namespace biasatomistic {

BiasAtomistic::BiasAtomistic(const ActionOptions & ao) :
   Action(ao),
   ActionAtomistic(ao),
   Bias(ao) 
{
}

//void BiasAtomistic::apply() {
//   // NOTE `forces' declared in ActionAtomicistic is the reduced vector
//   // of atomic forces of atoms that used in the action
//   // apply() and applyForce() of ActionAtomistic will both be called
//   // in PlumedMain
//   // ActionWithVirtualAtom::apply() propagates the force on vatom into atoms 
//   // used in the action, i.e. computes `forces'
//   // then ActionAtomistic::applyForce() add `forces' to full atomic forces
//   // vector `atoms.forces'
//   
//   // I will do forces = sum_i c_i^2 F_i in WEDS 
//   // so I do not need this
//}

}
}
