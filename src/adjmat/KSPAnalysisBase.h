#ifdef __PLUMED_HAS_BOOST_GRAPH

#ifndef __PLUMED_adjmat_KSPAnalysisBase_H
#define __PLUMED_adjmat_KSPAnalysisBase_H

#include "KSPBase.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class KSPAnalysisBase : public multicolvar::MultiColvarBase {
protected:
   KSPBase* kspbase;
public:
   static void registerKeywords( Keywords& keys );
   explicit KSPAnalysisBase(const ActionOptions&);
   bool isPeriodic();
   virtual double compute(
         const unsigned& tindex, multicolvar::AtomValuePack& myatoms) const
   { plumed_error(); }
   unsigned getNumberOfQuantities() const;
   unsigned getNumberOfDerivatives() const;
   const std::vector<std::pair<unsigned,unsigned>>&
      retrieveSPEdgeList(unsigned k) const;
   double getPathLength(unsigned k) const;
   double retrieveConnectionValue(
         const unsigned& i, const unsigned& j, std::vector<double>& vals) const;
   bool ifSymmetrize() const;
   bool ifInverseWeight() const;
   AdjacencyMatrixVessel* getAdjacencyVessel() const;
};

inline bool KSPAnalysisBase::isPeriodic()
{
   return mybasemulticolvars[0]->isPeriodic();
}

inline
unsigned KSPAnalysisBase::getNumberOfQuantities() const
{
   return kspbase->getNumberOfQuantities();
}

inline
unsigned KSPAnalysisBase::getNumberOfDerivatives() const
{
   return kspbase->getNumberOfDerivatives();
}

inline
double KSPAnalysisBase::retrieveConnectionValue(
      const unsigned& i, const unsigned& j, std::vector<double>& vals) const
{
   return kspbase->retrieveConnectionValue(i, j, vals);
}

inline
bool KSPAnalysisBase::ifInverseWeight() const
{
   return kspbase->ifInverseWeight();
}

inline
bool KSPAnalysisBase::ifSymmetrize() const
{
   return kspbase->ifSymmetrize();
}

inline
AdjacencyMatrixVessel* KSPAnalysisBase::getAdjacencyVessel() const {
   return kspbase->getAdjacencyVessel();
}

inline
double KSPAnalysisBase::getPathLength(unsigned k) const
{
   return kspbase->getPathLengths()[k];
}

inline
const std::vector<std::pair<unsigned,unsigned>>&
      KSPAnalysisBase::retrieveSPEdgeList(unsigned k) const
{
   return kspbase->retrieveSPEdgeList(k);
}


}
}

#endif

#endif
