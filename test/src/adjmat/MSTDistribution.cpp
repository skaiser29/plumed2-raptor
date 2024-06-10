/*
 * author: Chenghan Li
 */
#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "MSTAnalysisBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

class MSTDistribution : public MSTAnalysisBase {
   unsigned nderivs;
public:
   static void registerKeywords(Keywords& keys);
   explicit MSTDistribution(const ActionOptions&);
   void calculate();
   //double compute(const unsigned&, multicolvar::AtomValuePack&) const;
   void performTask(
     const unsigned& taskid, const unsigned& current, MultiValue& myvals) const;
};

PLUMED_REGISTER_ACTION(MSTDistribution,"MSTDISTRIBUTION")

void MSTDistribution::registerKeywords(Keywords& keys)
{
   MSTAnalysisBase::registerKeywords(keys);
   //TODO gives error when using switching function based cv to filter atoms and use MAX
   keys.use("SUM"); keys.use("MAX"); keys.use("MEAN");
   //componentsAreNotOptional(keys);
   //keys.addOutputComponent("edge","default","edge length in the MST");
   //keys.add("optional","MAXNUM","max number of components");
}

MSTDistribution::MSTDistribution(const ActionOptions&ao):
   Action(ao),
   MSTAnalysisBase(ao),
   nderivs(0)
{
   //this part is now in MinimumSpanningTree
   //max_num_comps = mstbase->getNumberOfNodes() - 1;
   //unsigned tmp_max_num_comps = 0;
   //parse("MAXNUM",tmp_max_num_comps);
   //if(tmp_max_num_comps !=0 && tmp_max_num_comps < max_num_comps) 
   //   max_num_comps = tmp_max_num_comps;
   //for(unsigned i = 0; i < max_num_comps; ++i) {
   //   std::string num; Tools::convert(i + 1, num);
   //   addComponentWithDerivatives("edge-"+num);
   //   componentIsNotPeriodic("edge-"+num);
   //   getPntrToComponent(i)
   //      ->resizeDerivatives(mstbase->getNumberOfDerivatives());
   //}
   for(unsigned i = 0; i < mstbase->getMSTTrueSize(); ++i) {
      addTaskToList(i); taskFlags[i] = 1;
   }
   std::vector<AtomNumber> fake_atoms; setupMultiColvarBase(fake_atoms);
}

//double MSTDistribution::compute(
//      const unsigned& tind, multicolvar::AtomValuePack& myatoms) const
void MSTDistribution::calculate()
{
   nderivs = getNumberOfDerivatives();

//not affected by filter or epsilon
if(comm.Get_rank()==0) printf("nderivs = %u\n",nderivs);//@@@

//   for(unsigned i = 0; i < max_num_comps; ++i) {
//      getPntrToComponent(i)->set(arches[i]);
//   }

   //do we need this???
   lockContributors();
   runAllTasks();

}

void MSTDistribution::performTask(
     const unsigned& taskid, const unsigned& current, MultiValue& myvals) const
{
   myvals.setValue(0, 1.0);
   myvals.setValue(
         //1, mstbase->retrieveMSTArchLengths()[current]
         1, getMSTArchLength(current)
         );

   //derivatives
   if(!doNotCalculateDerivatives()) {
      unsigned i1, i2;
      //i1 = mstbase->retrieveMSTEdgeList()[current].first;
      //i2 = mstbase->retrieveMSTEdgeList()[current].second;
      i1 = getMSTEdge(current).first;
      i2 = getMSTEdge(current).second;
      //std::vector<double> vals(mstbase->getNumberOfQuantities());
      std::vector<double> vals(getNumberOfQuantities());
      std::vector<double> convals(getNumberOfQuantities());
      //mstbase->retrieveConnectionValue(i1, i2, vals);
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
//
//if(i1==3035 && i2==2544) printf("a_3035_2544 involves %u tvals.getDerivative(1,kat) = %f\n",kat,tvals.getDerivative(1,kat)); //@@@
//
               myvals.addDerivative(
                     1, kat, -(vals[0]/vals[1] * tvals.getDerivative(1, kat)
                        + tvals.getDerivative(0, kat))/vals[1]);
            }
         }
      }
   }

//printf("%f == %f?\n",vals[0]/vals[1], mstbase->retrieveMSTArchLengths()[current]);//@@@
}

}
}


#endif
