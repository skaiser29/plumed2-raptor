/*
 * author: Chenghan Li
 */

#ifndef __PLUMED_adjmat_KSPBase_H
#define __PLUMED_adjmat_KSPBase_H

#include "ActionWithInputMatrix.h"
#include "multicolvar/AtomValuePack.h"
#include "AdjacencyMatrixVessel.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/yen_ksp.hpp>
#include <boost/graph/custom_dijkstra_call.hpp>

namespace PLMD {
namespace adjmat {

class KSPBase : public ActionWithInputMatrix {
public:
   //whether other class can access these type defs ???
   const double iepsilon = 1 / epsilon; 
   typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, double>> Graph;
   typedef boost::graph_traits<Graph>::edge_descriptor Edge;
   typedef std::pair<unsigned,unsigned> unsigned_pair;
   static void registerKeywords( Keywords& keys );
   explicit KSPBase(const ActionOptions&);
   void calculate();
   // do not apply forces here, instead do the forces in KSPAnalysis
   void apply() {}
   bool ifInverseWeight() const;
   bool ifSymmetrize() const;
   const std::list< std::pair< double,std::list<Edge> > >& retrieveKSP() const;
   const std::vector<std::vector<unsigned_pair>>& retrieveKSPEdgeList() const;
   const std::vector<unsigned_pair>& retrieveSPEdgeList(unsigned k) const;
   //const std::vector<std::vector<double>>& retrieveKSPArchLengths() const;
   //const std::vector<double>& retrieveSPArchLengths(unsigned k) const;
   const std::vector<unsigned>& getKSPSize() const;
   unsigned getSPSize(unsigned k) const;
   // Get the number of paths found; could be less than npaths*Nsource*Ntarget
   unsigned getTrueNPaths() const;
   unsigned getK() const;
   const std::vector<double>& getPathLengths() const;

protected:
   //whether to take the inverse, symmetrize the matrix
   //whether to give warnings or halt the program
   bool inverse_weight, symmetrize;//, warnings, halting, check;
   std::vector<AtomNumber> start, end;
   virtual void buildKSP()=0;
   //number of edges in resulting ksp
   std::vector<unsigned> ksp_nedges;
   //stores the resulting ksp
   std::list< std::pair< double,std::list<Edge> > > ksp;
   //stores the edge list in ksp
   std::vector<std::vector<unsigned_pair>> ksp_edge_list;
   //stores the arch lengths in ksp
   //std::vector<std::vector<double>> ksp_arch_lengths;
   //stores the path lengths
   // the K in KSP
   unsigned npaths; 
   // path lengths
   std::vector<double> ksp_path_lengths;
};

inline
bool KSPBase::ifInverseWeight() const 
{
   return inverse_weight;
}

inline
bool KSPBase::ifSymmetrize() const
{
   return symmetrize;
}

inline
const std::list<std::pair<double,std::list<KSPBase::Edge>>>& 
KSPBase::retrieveKSP() const
{
   return ksp;
}

//inline
//const std::vector<std::vector<double>>& KSPBase::retrieveKSPArchLengths() const
//{
//   return ksp_arch_lengths;
//}

//inline
//const std::vector<std::vector<double>>& 
//KSPBase::retrieveSPArchLengths(unsigned k) const
//{
//   return ksp_arch_lengths[k];
//}

inline
const std::vector<std::vector<KSPBase::unsigned_pair>>& 
                             KSPBase::retrieveKSPEdgeList() const
{
   return ksp_edge_list;
}

inline
const std::vector<KSPBase::unsigned_pair>& 
KSPBase::retrieveSPEdgeList(unsigned k) const
{
   return ksp_edge_list[k];
}

inline
const std::vector<unsigned>& KSPBase::getKSPSize() const
{
   return ksp_nedges;
}

inline
unsigned KSPBase::getSPSize(unsigned k) const
{
   return ksp_nedges[k];
}

inline
unsigned KSPBase::getTrueNPaths() const 
{
   return ksp.size();
}

inline
const std::vector<double>& KSPBase::getPathLengths() const
{
   return ksp_path_lengths;
}

inline
unsigned KSPBase::getK() const
{
   return npaths;
}

}
}

#endif

