/*
 * author: Chenghan Li
 */
#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "KSPAnalysisBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

class KSPDistribution : public KSPAnalysisBase {
   unsigned nderivs;
public:
   static void registerKeywords(Keywords& keys);
   explicit KSPDistribution(const ActionOptions&);
   void calculate();
   //double compute(const unsigned&, multicolvar::AtomValuePack&) const;
   void performTask(
     const unsigned& taskid, const unsigned& current, MultiValue& myvals) const;
};

PLUMED_REGISTER_ACTION(KSPDistribution,"KSPDISTRIBUTION")

void KSPDistribution::registerKeywords(Keywords& keys)
{
   KSPAnalysisBase::registerKeywords(keys);
   //TODO gives error when using switching function based cv to filter atoms and use MAX
   keys.use("SUM"); keys.use("MIN"); keys.use("MEAN");
   keys.use("LOG_MIN");
   keys.use("LOWEST");
   //componentsAreNotOptional(keys);
   //keys.addOutputComponent("edge","default","edge length in the KSP");
   //keys.add("optional","MAXNUM","max number of components");
}

KSPDistribution::KSPDistribution(const ActionOptions&ao):
   Action(ao),
   KSPAnalysisBase(ao),
   nderivs(0)
{
   for(unsigned i = 0; i < kspbase->getK(); ++i) {
      addTaskToList(i); taskFlags[i] = 1;
   }
   std::vector<AtomNumber> fake_atoms; setupMultiColvarBase(fake_atoms);
}

void KSPDistribution::calculate()
{
   nderivs = getNumberOfDerivatives();
   deactivateAllTasks();

   //not affected by filter or epsilon
   //if(comm.Get_rank()==0) printf("nderivs = %u\n",nderivs);//@@@

   for(unsigned i = 0; i < kspbase->getTrueNPaths(); ++i)
      taskFlags[i] = 1;
   lockContributors();
   runAllTasks();

}

void KSPDistribution::performTask(
     const unsigned& taskid, const unsigned& current, MultiValue& myvals) const
{
   myvals.setValue(0, 1.0);
   myvals.setValue(
         1, getPathLength(current)
         );

   //derivatives
   if(!doNotCalculateDerivatives()) {
      for(auto edge : retrieveSPEdgeList(current))
      {
         unsigned i1, i2;
         i1 = edge.first;
         i2 = edge.second;
         std::vector<double> vals(getNumberOfQuantities());
         std::vector<double> convals(getNumberOfQuantities());
         retrieveConnectionValue(i1, i2, vals);
         retrieveConnectionValue(i2, i1, convals);

         //wheter or not symmetrize, tvals are always needed
         MultiValue tvals(getNumberOfQuantities(), getNumberOfDerivatives());
         MultiValue contvals(getNumberOfQuantities(), getNumberOfDerivatives());
         unsigned ele = getAdjacencyVessel()
            ->getStoreIndexFromMatrixIndices(i1, i2);
         getAdjacencyVessel()->retrieveDerivatives(ele, false, tvals);

//if(current==0) printf("tvals.getNumberActive() = %d\n",tvals.getNumberActive());//@@@
//printf("tvals.getNumberActive() = %d\n",tvals.getNumberActive());//@@@


         if(!ifInverseWeight()) {
            /*--- no inverse ---*/
            // tvals is \partial a_ij / \partial atom_k
            if(ifSymmetrize()) {
               //ifsymmetrize derivs come from both aij and aji
               ele = getAdjacencyVessel()->getStoreIndexFromMatrixIndices(i2, i1);
               getAdjacencyVessel()->retrieveDerivatives(ele, false, contvals);
               for(unsigned k = 0; k < tvals.getNumberActive(); ++k) {
                  unsigned kat = tvals.getActiveIndex(k);
                  myvals.addDerivative(
                        1, kat, 0.5 * (vals[0] * tvals.getDerivative(1, kat) +
                        vals[1] * tvals.getDerivative(0, kat)));
               }
               for(unsigned k = 0; k < contvals.getNumberActive(); ++k) {
                  unsigned kat = contvals.getActiveIndex(k);
                  myvals.addDerivative(
                        1, kat, 0.5 * (convals[0] * contvals.getDerivative(1, kat) +
                        convals[1] * contvals.getDerivative(0, kat)));
               }
            } else {
               //else derivs only come from aij
               for(unsigned k = 0; k < tvals.getNumberActive(); ++k) {
                  unsigned kat = tvals.getActiveIndex(k);
                  myvals.addDerivative(
                        1, kat, vals[0] * tvals.getDerivative(1, kat) +
                        vals[1] * tvals.getDerivative(0, kat));
               }
            }
         } else {
            /*--- weights are inversed, derivs need a bit more care ---*/
            if(ifSymmetrize()) {
               //ifsymmetrize derivs come from both aij and aji
               ele = getAdjacencyVessel()->getStoreIndexFromMatrixIndices(i2, i1);
               getAdjacencyVessel()->retrieveDerivatives(ele, false, contvals);
               for(unsigned k = 0; k < tvals.getNumberActive(); ++k) {
                  unsigned kat = tvals.getActiveIndex(k);
                  myvals.addDerivative(
                        1, kat, -0.5 * (vals[0]/vals[1] * tvals.getDerivative(1, kat)
                           + tvals.getDerivative(0, kat))/vals[1]);
               }
               for(unsigned k = 0; k < contvals.getNumberActive(); ++k) {
                  unsigned kat = contvals.getActiveIndex(k);
                  myvals.addDerivative(
                        1, kat, -0.5 * (
                           convals[0]/convals[1] * contvals.getDerivative(1, kat)
                           + contvals.getDerivative(0, kat) ) / convals[1]);
               }
            } else {
               //else derivs only come from aij
               for(unsigned k = 0; k < tvals.getNumberActive(); ++k) {
                  unsigned kat = tvals.getActiveIndex(k);
                  myvals.addDerivative(
                        1, kat, -(vals[0]/vals[1] * tvals.getDerivative(1, kat)
                           + tvals.getDerivative(0, kat))/vals[1]);
               }
            }
         }
      }
   }

//printf("%f == %f?\n",vals[0]/vals[1], mstbase->retrieveMSTArchLengths()[current]);//@@@
}

}
}


#endif
