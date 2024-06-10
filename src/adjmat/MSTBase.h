/*
 * author: Chenghan Li
 */
#ifdef __PLUMED_HAS_BOOST_GRAPH

#ifndef __PLUMED_adjmat_MSTBase_H
#define __PLUMED_adjmat_MSTBase_H

#include "ActionWithInputMatrix.h"
#include "multicolvar/AtomValuePack.h"
#include "AdjacencyMatrixVessel.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

namespace PLMD {
namespace adjmat {

class MSTBase : public ActionWithInputMatrix {
public:
   //whether other class can access these type defs ???
   const double iepsilon = 1 / epsilon; 
   typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, double>> Graph;
   typedef boost::graph_traits<Graph>::edge_descriptor Edge;
   static void registerKeywords( Keywords& keys );
   explicit MSTBase(const ActionOptions&);
   void calculate();
   // do not apply forces here, instead do the forces in MSTAnalysis
   void apply() {}
   bool ifInverseWeight() const;
   bool ifSymmetrize() const;
   const std::vector<Edge>& retrieveMST() const;
   const std::vector<std::pair<unsigned,unsigned>>& retrieveMSTEdgeList() const;
   const std::vector<double>& retrieveMSTArchLengths() const;
   unsigned getMSTSize() const;
   unsigned getMSTTrueSize() const;

protected:
   //whether to take the inverse, symmetrize the matrix
   //whether to give warnings
   bool inverse_weight, symmetrize, warnings;
   virtual void buildMST()=0;
   //number of edges in resulting mst
   unsigned mst_nedges;
   //stores the resulting mst
   std::vector<Edge> mst;
   //stores the edge list in mst
   std::vector<std::pair<unsigned,unsigned>> mst_edge_list;
   //stores the arch lengths in mst
   std::vector<double> mst_arch_lengths;
   //max number of components
   unsigned max_num;
};

inline
bool MSTBase::ifInverseWeight() const 
{
   return inverse_weight;
}

inline
bool MSTBase::ifSymmetrize() const
{
   return symmetrize;
}

inline
const std::vector<MSTBase::Edge>& MSTBase::retrieveMST() const
{
   return mst;
}

inline
const std::vector<double>& MSTBase::retrieveMSTArchLengths() const
{
   return mst_arch_lengths;
}

inline
const std::vector<std::pair<unsigned,unsigned>>& 
                             MSTBase::retrieveMSTEdgeList() const
{
   return mst_edge_list;
}

inline
unsigned MSTBase::getMSTSize() const
{
   return mst_nedges;
}

inline
unsigned MSTBase::getMSTTrueSize() const
{
   return max_num;
}

}
}

#endif

#endif
