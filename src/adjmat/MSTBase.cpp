/*
 * author: Chenghan Li
 */
#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "MSTBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

void MSTBase::registerKeywords(Keywords& keys) {
   ActionWithInputMatrix::registerKeywords(keys);
   keys.addFlag("NOINVERSE",false,"Do not take the inverse of the weights");
   keys.addFlag("NOSYMMETRIZE",false,"Do not symmetrize the input matrix");
   keys.addFlag("NOWARNINGS",false,"Do not give warnings");
   keys.add("optional","MAXNUM","max number of components");
}

MSTBase::MSTBase(const ActionOptions&ao):
   Action(ao),
   ActionWithInputMatrix(ao),
   inverse_weight(true),
   symmetrize(true),
   warnings(true)
{
   usespecies = false;
   max_num = getNumberOfNodes() - 1;
   unsigned tmp_max_num = 0;
   parse("MAXNUM",tmp_max_num);
   if(tmp_max_num != 0 && tmp_max_num < max_num) max_num = tmp_max_num;

   bool noinverse_weight, nosymmetrize, nowarnings;
   parseFlag("NOINVERSE", noinverse_weight);
   parseFlag("NOSYMMETRIZE", nosymmetrize);
   parseFlag("NOWARNINGS", nowarnings);
   inverse_weight = !noinverse_weight;
   symmetrize = !nosymmetrize;
   warnings = !nowarnings;
   if( getAdjacencyVessel() ){
      //if( getNumberOfNodeTypes()!=1 ) 
      //   error("should only be running MST with one base multicolvar in function");
      if( !getAdjacencyVessel()->isSymmetric() and !symmetrize) 
         error("asymmetric graph needs symmetrization");
      if( getAdjacencyVessel()->isSymmetric() and symmetrize) 
         log.printf("  WARNING: symmetric graph does not need symmetrization\n");
   }
   if( keywords.exists("MATRIX") ){
      std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
   }
}

void MSTBase::calculate() 
{
   //use boost to get MST
   //mst_nedges, mst, mst_edge_list, mst_arch_lengths are built

//printf("MSTBase::calcualte called on rank %d\n",comm.Get_rank());//@@@

   buildMST();
}

}
}

#endif
