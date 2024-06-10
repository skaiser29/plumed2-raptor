#ifdef __PLUMED_HAS_BOOST_GRAPH

#ifndef __PLUMED_adjmat_MSTAnalysisBase_H
#define __PLUMED_adjmat_MSTAnalysisBase_H

#include "MSTBase.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class MSTAnalysisBase : public multicolvar::MultiColvarBase {
protected:
   MSTBase* mstbase;
public:
   static void registerKeywords( Keywords& keys );
   explicit MSTAnalysisBase(const ActionOptions&);
   bool isPeriodic();
   virtual double compute(
         const unsigned& tindex, multicolvar::AtomValuePack& myatoms) const
   { plumed_error(); }
   unsigned getNumberOfQuantities() const;
   unsigned getNumberOfDerivatives() const;
   const std::vector<std::pair<unsigned,unsigned>>&
      retrieveMSTEdgeList() const;
   const std::vector<double>& retrieveMSTArchLengths() const;
   const std::pair<unsigned,unsigned> getMSTEdge(const int i) const;
   double getMSTArchLength(const int i) const;
   double retrieveConnectionValue(
         const unsigned& i, const unsigned& j, std::vector<double>& vals) const;
   bool ifSymmetrize() const;
   bool ifInverseWeight() const;
   AdjacencyMatrixVessel* getAdjacencyVessel() const;
};

inline bool MSTAnalysisBase::isPeriodic()
{
   return mybasemulticolvars[0]->isPeriodic();
}

inline
unsigned MSTAnalysisBase::getNumberOfQuantities() const
{
   return mstbase->getNumberOfQuantities();
}

inline
unsigned MSTAnalysisBase::getNumberOfDerivatives() const
{
   return mstbase->getNumberOfDerivatives();
}

inline
const std::vector<std::pair<unsigned,unsigned>>&
      MSTAnalysisBase::retrieveMSTEdgeList() const
{
   return mstbase->retrieveMSTEdgeList();
}

inline
const std::vector<double>& MSTAnalysisBase::retrieveMSTArchLengths() const
{
   return mstbase->retrieveMSTArchLengths();
}

inline
const std::pair<unsigned,unsigned> 
             MSTAnalysisBase::getMSTEdge(const int i) const
{
   return mstbase->retrieveMSTEdgeList()[i];
}

inline
double MSTAnalysisBase::getMSTArchLength(const int i) const
{
   return mstbase->retrieveMSTArchLengths()[i];
}

inline
double MSTAnalysisBase::retrieveConnectionValue(
      const unsigned& i, const unsigned& j, std::vector<double>& vals) const
{
   return mstbase->retrieveConnectionValue(i, j, vals);
}

inline
bool MSTAnalysisBase::ifInverseWeight() const
{
   return mstbase->ifInverseWeight();
}

inline
bool MSTAnalysisBase::ifSymmetrize() const
{
   return mstbase->ifSymmetrize();
}

inline
AdjacencyMatrixVessel* MSTAnalysisBase::getAdjacencyVessel() const {
   return mstbase->getAdjacencyVessel();
}

}
}

#endif

#endif
