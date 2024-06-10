/*
 * author: Chenghan Li
 */
#ifdef __PLUMED_HAS_BOOST_GRAPH

#include "AdjacencyMatrixBase.h"
#include "MSTBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

class MinimumSpanningTree : public MSTBase {
   void preprocess();
   void buildMST();
   DynamicList<unsigned> active_elements;
   Matrix<double> thematrix, matrix_inversed;
   std::vector<std::pair<unsigned,unsigned>> edge_list;
   //if edges in this list are in the resulting MST, give warning
   std::vector<std::pair<unsigned,unsigned>> warning_edge_list;
   std::vector<double> processed_archlengths;
   unsigned nedges;
public:
   static void registerKeywords( Keywords& keys );
   explicit MinimumSpanningTree(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(MinimumSpanningTree,"MINIMUMSPANNINGTREE")

void MinimumSpanningTree::registerKeywords( Keywords& keys )
{
   MSTBase::registerKeywords(keys);
   componentsAreNotOptional(keys); //√√√
   keys.addOutputComponent("edge","default","edge length in the MST");
   //these things do not work; need another class MSTAnalysis
   //keys.use("MEAN"); keys.use("SUM");

//printf("epsilon = %e\n",epsilon); //@@@
}

MinimumSpanningTree::MinimumSpanningTree(const ActionOptions&ao):
   Action(ao),
   MSTBase(ao)//,
   //thematrix(getNumberOfNodes(), getNumberOfNodes()),
   //matrix_inversed(getNumberOfNodes(), getNumberOfNodes())
{
   //TODO handle with symmetry
   //think: how to reduce edge_list to n*(n-1)/2 after symmetrization??
   //resize the edge_list to ensure it can hold all the edges
   //some of the edges will be omitted because they are too small
   //so edge_list is typically not full
   //they are omitted not because their derivs are small but because they
   //are unlikely to be in the MST
   //in an asymmetric graph, A->B is kept in edge_list while B->A could be omitted
   //I will just keep A->B in edge_list but give warning when A->B is in resulting MST
   edge_list.resize(getNumberOfNodes() * (getNumberOfNodes() - 1));

   for(unsigned i=0; i<max_num; ++i){ //√√√
      std::string num; Tools::convert(i+1,num);
      addComponentWithDerivatives("edge-"+num);
      componentIsNotPeriodic("edge-"+num);
      getPntrToComponent(i)->resizeDerivatives( getNumberOfDerivatives() );
   }

   unsigned ntasks = getNumberOfNodes() * (getNumberOfNodes() - 1);
   if(getAdjacencyVessel()->isSymmetric()) ntasks /= 2;
   for(unsigned i = 0; i < ntasks; ++i) active_elements.addIndexToList(i);

//
//printf("getNumberOfNodes() = %d, getNumberOfDerivatives() = %d\n", getNumberOfNodes(), getNumberOfDerivatives()); //@@@
//
}

void MinimumSpanningTree::preprocess() {
   processed_archlengths.clear();
   warning_edge_list.clear();
   //loop in edge list to reduce computational cost
   if(inverse_weight) {
      for(unsigned i = 0; i < nedges; ++i) {
         unsigned i1 = edge_list[i].first;
         unsigned i2 = edge_list[i].second;
         std::vector<double> vals(getNumberOfQuantities());
         unsigned ele = getAdjacencyVessel()
                  ->getStoreIndexFromMatrixIndices(i1, i2);
         getAdjacencyVessel()->retrieveValueWithIndex(ele, false, vals);
         processed_archlengths.push_back(vals[0]/vals[1]);

//@@@
if(comm.Get_rank()==0) {
if(vals[0]<0.7) printf("i1 = %d i2 = %d vals[0] = %f vals[1] = %f\n",getAbsoluteIndexOfCentralAtom(i1).serial(),getAbsoluteIndexOfCentralAtom(i2).serial(),vals[0],vals[1]);
//printf("i1 = %d i2 = %d vals[0] = %f vals[1] = %f\n",getAbsoluteIndexOfCentralAtom(i1).serial(),getAbsoluteIndexOfCentralAtom(i2).serial(),vals[0],vals[1]);
if((i1 == 2544 and i2==809) or (i1 == 809 and i2==2544)) {
   printf("2544vs809: vals[0] = %e vals[1] = %e\n",vals[0],vals[1]);
}
if((i1 == 10 and i2==57) or (i1 == 57 and i2==10)) {
   printf("10vs57: vals[0] = %e vals[1] = %e\n",vals[0],vals[1]);
}
if((i1 == 10 and i2==41) or (i1 == 41 and i2==10)) {
   printf("10vs41: vals[0] = %e vals[1] = %e\n",vals[0],vals[1]);
}
if(i1 == 57 and i2==41) {
   printf("57vs41: vals[0] = %e vals[1] = %e\n",vals[0],vals[1]);
}
if(i1 == 2698 and i2==1868) {
   printf("2698vs1868: vals[0] = %e vals[1] = %e\n",vals[0],vals[1]);
}
if((i1 == 3035 and i2 == 2544) or (i2 == 3035 and i1 == 2544)) {
   printf("3035vs2544: vals[0] = %e vals[1] = %e\n",vals[0],vals[1]);
   printf("%u %u\n",getAbsoluteIndexOfCentralAtom(3035).serial(),getAbsoluteIndexOfCentralAtom(2544).serial());
}
}
//if(getAdjacencyVessel()->getMatrixAction()->mybasemulticolvars.size() && vals[0] > 0.1) {
//if(vals[0] > 0.1) {
//   printf("%dvs%d: vals[0] = %e vals[1] = %e\n",i1,i2,vals[0],vals[1]);
//}
//@@@

      }
      if(symmetrize) {
         //inverse_weight && symmetrize
         for(unsigned i = 0; i < nedges; ++i) {
            std::pair<unsigned,unsigned> edge = edge_list[i];
            std::pair<unsigned,unsigned> con_edge(edge.second, edge.first);
            std::vector<std::pair<unsigned,unsigned>>::iterator it = 
               std::find(edge_list.begin() + i, edge_list.end(), con_edge);
            if(it != edge_list.end()) {
               //if the conjugated part not found
               processed_archlengths[i] *= 0.5;
               processed_archlengths[i] += 0.5 * iepsilon;
               if(warnings) warning_edge_list.push_back(edge);
            } else {
               //else average the conjugated part
               unsigned j = it - edge_list.begin();
               processed_archlengths[i] = 
                  (processed_archlengths[i] + processed_archlengths[j]) / 2.0;
               processed_archlengths[j] = processed_archlengths[i];
            }
         }
      } 
   } else {
      for(unsigned i = 0; i < nedges; ++i) {
         unsigned i1 = edge_list[i].first;
         unsigned i2 = edge_list[i].second;
         std::vector<double> vals(getNumberOfQuantities());
         processed_archlengths.push_back(retrieveConnectionValue(i1, i2, vals));
      }
      if(symmetrize) {
         //!inverse_weight && symmetrize
         for(unsigned i = 0; i < nedges; ++i) {
            std::pair<unsigned,unsigned> edge = edge_list[i];
            std::pair<unsigned,unsigned> con_edge(edge.second, edge.first);
            std::vector<std::pair<unsigned,unsigned>>::iterator it = 
               std::find(edge_list.begin() + i, edge_list.end(), con_edge);
            if(it != edge_list.end()) {
               //if the conjugated part not found
               processed_archlengths[i] *= 0.5;
            } else {
               //else average the conjugated part
               unsigned j = it - edge_list.begin();
               processed_archlengths[i] = 
                  (processed_archlengths[i] + processed_archlengths[j]) / 2.0;
               processed_archlengths[j] = processed_archlengths[i];
            }
         }
      }
   }

//   if(inverse_weight) {
//      //loop in full matrix instead of active_elements to avoid a linear search
//      for(int i=0; i<getNumberOfNodes(); ++i) {
//         for(int j=i + 1; j<getNumberOfNodes(); ++j) {
//            if(!symmetrize) { //no symmetrize requested so input is symmetric
//               if(thematrix[i][j] < epsilon) 
//                  matrix_inversed[i][j] = matrix_inversed[j][i] = iepsilon;
//               else
//                  matrix_inversed[i][j] = 
//                     matrix_inversed[j][i] = 1.0 / thematrix[i][j];
//
//if(i==36 && j==43 && comm.Get_rank()==0) printf("thematrix[%d][%d]=%f\n",i,j,thematrix[i][j]);//@@@
//
//            } else {
//               //if both inverse and symmetrized are requested
//               //this kind of thing could happen:
//               //E(A->B) >= epsilon but E(B->A) < epsilon
//               //then the retrieved edge list will not contain E(B->A)
//               //I will push A->B to the warning_edge_list 
//               //if the resulting MST contains any of the warning_edge_list,
//               //I will give warnings 
//               warning_edge_list.clear();
//               bool f1 = (thematrix[i][j] < epsilon);
//               bool f2 = (thematrix[j][i] < epsilon);
//               if(f1 && f2) 
//                  matrix_inversed[i][j] = matrix_inversed[j][i] = iepsilon;
//               else if(!f1 && !f2)
//                  matrix_inversed[i][j] = matrix_inversed[j][i] 
//                     = 0.5 / thematrix[i][j] + 0.5 / thematrix[j][i];
//               else if(!f1 && f2) {
//                  matrix_inversed[i][j] = matrix_inversed[j][i] 
//                     = 0.5 / thematrix[i][j] + 0.5 * iepsilon; 
//                  warning_edge_list.push_back(
//                        std::pair<unsigned,unsigned>(i,j)
//                        );
//               }
//               else if(f1 && !f2) {
//                  matrix_inversed[i][j] = matrix_inversed[j][i] 
//                     = 0.5 / thematrix[j][i] + 0.5 * iepsilon; 
//                  warning_edge_list.push_back(
//                        std::pair<unsigned,unsigned>(j,i)
//                        );
//               }
//            }
//         }
//      }
//   }

//   if(symmetrize) {
//      if(inverse_weight && warnings) {
//         //if both inverse and symmetrized are requested
//         //this kind of thing could happen:
//         //E(A->B) >= epsilon but E(B->A) < epsilon
//         //then the retrieved edge list will not contain E(B->A)
//         //I will push A->B to the warning_edge_list 
//         //if the resulting MST contains any of the warning_edge_list,
//         //I will give warnings 
//      }
//   }
}

void MinimumSpanningTree::buildMST() {

//printf("MinimumSpanningTree::buildMST called on rank %d\n",comm.Get_rank());//@@@

   //getAdjacencyVessel()->retrieveMatrix(active_elements, thematrix);
   //get the edge list
   nedges = 0;
   // in AdjacencyMatrixVessel.cpp: retrieveEdgeList returns edges involving
   // stored values; nedges is affected by epsilon
   getAdjacencyVessel()->retrieveEdgeList(nedges, edge_list);

//if(comm.Get_rank()==0) { //@@@
//FILE* fp = fopen("mat.dat","w");
//for(int i=0; i<getNumberOfNodes(); ++i) {
//   for(int j=0; j<getNumberOfNodes()-1; ++j) {
//      fprintf(fp,"%f,",thematrix[i][j]);
//   }
//   fprintf(fp,"%f\n",thematrix[i][getNumberOfNodes()-1]);
//}
//fclose(fp);
//fp = fopen("invmat.dat","w");
//for(int i=0; i<getNumberOfNodes(); ++i) {
//   for(int j=0; j<getNumberOfNodes()-1; ++j) {
//      fprintf(fp,"%f,",matrix_inversed[i][j]);
//   }
//   fprintf(fp,"%f\n",matrix_inversed[i][getNumberOfNodes()-1]);
//}
//fclose(fp);
//}

   preprocess();

//if(comm.Get_rank()==0) printf("matrix_inversed[%d][%d] = %f\n",6,59,matrix_inversed[6][59]); //@@@
if(comm.Get_rank()==0) printf("nedges = %d\n", nedges); //@@@

   if(!nedges) plumed_merror("no edges found");
//   std::vector<double> weights;
//   for(unsigned i = 0; i < nedges; ++i) {
//      std::pair<unsigned,unsigned> edge = edge_list[i];
////
////if(comm.Get_rank()==0)
////printf("matrix_iversed[%d][%d] = %f\n",edge.first,edge.second,matrix_inversed[edge.first][edge.second]); //@@@
////
//      weights.push_back(matrix_inversed[edge.first][edge.second]);
//   }
   
//   //exploit active_elements instead of loop in edge_list; orders are the same??
//   for(unsigned i = 0; i < active_elements.getNumberActive(); ++i) {
//      unsigned j, k;
//      getAdjacencyVessel()->getMatrixIndices(active_elements[i], j, k);
//
//printf("matrix_iversed[%d][%d] = %f\n",j,k,matrix_inversed[j][k]); //@@@
//
//      weights.push_back(matrix_inversed[j][k]);
//   }
   //boost MST
//   Graph g(
//         edge_list.begin(), edge_list.begin() + nedges, weights.begin(), getNumberOfNodes()
//       );
   Graph g(
         edge_list.begin(), edge_list.begin() + nedges, 
         processed_archlengths.begin(), getNumberOfNodes()
       );
   mst.clear();
   boost::property_map<Graph, boost::edge_weight_t>::type 
                                            weight = get(boost::edge_weight, g);
   boost::kruskal_minimum_spanning_tree(g, std::back_inserter(mst));

   mst_nedges = mst.size();
   mst_edge_list.clear(); mst_arch_lengths.clear();

//if(comm.Get_rank()==0) {
//   printf("step %u\n", getStep());
//   for(std::vector<std::pair<unsigned,unsigned>>::iterator it = edge_list.begin();
//         it != edge_list.begin() + nedges; ++it) {
//      //printf("%d <--> %d: %f\n",it->first,it->second,matrix_inversed[it->first][it->second]);
//      printf("%d <--> %d\n",it->first,it->second);
//   }
//} //@@@


//printf("mst got on rank %d\n",comm.Get_rank()); //@@@

   //for(unsigned i=0;i<getNumberOfComponents();++i){
   //   getPntrToComponent(i)->set( 1.0 );
   //}
   
   //reverse iter for descending sort
   for(std::vector<Edge>::reverse_iterator it = mst.rbegin();
         it != mst.rend(); ++it)
   {
      int i1 = boost::source(*it, g), i2 = boost::target(*it, g);
//@@@
//printf("it goes to %d on rank %d\n",static_cast<unsigned>(it-mst.begin()),comm.Get_rank()); //@@@
if(comm.Get_rank()==0) {
int i1 = boost::source(*it, g), i2 = boost::target(*it, g);
//printf("%d-th edge: %d <--> %d: %f matrix_inversed[%d][%d] = %f,%f\n",static_cast<unsigned>(it-mst.rbegin()),i1,i2,weight[*it],i1,i2,matrix_inversed[i1][i2],matrix_inversed[i2][i1]);
//printf("%d-th edge: %d <--> %d : %f\n",static_cast<unsigned>(it-mst.rbegin()),i1,i2,weight[*it]);
printf("%d-th edge: %d <--> %d : %f\n",static_cast<unsigned>(it-mst.rbegin()),getAbsoluteIndexOfCentralAtom(i1).serial(),getAbsoluteIndexOfCentralAtom(i2).serial(),weight[*it]);
}
//@@@

      unsigned icomp = it - mst.rbegin();
      if(icomp >= max_num) break;
      getPntrToComponent( //√√√
            static_cast<unsigned>(icomp)
            //)->set(matrix_inversed[i1][i2]);
            )->set(weight[*it]);
      mst_edge_list.push_back(std::pair<unsigned,unsigned>(i1, i2));
      mst_arch_lengths.push_back(weight[*it]);
   }

   if(warnings && !warning_edge_list.empty()) {
      for(unsigned i = 0; i < mst_edge_list.size(); ++i) {
         std::vector<std::pair<unsigned,unsigned>>::iterator it = 
            std::find(
            warning_edge_list.begin(), warning_edge_list.end(), mst_edge_list[i]
                  );
         if(it != warning_edge_list.end()) {
            log.printf("  WARNING: Edge(%d,%d) = %f is in MST\n",
                  it->first, it->second, 
                  mst_arch_lengths[it - warning_edge_list.begin()]);
         }
      }
   }

//@@@
//if(comm.Get_rank()==0) {
//   for(unsigned i = 0; i < mst_arch_lengths.size(); ++i) {
//      printf("MST: arch-%d = %f\n", i + 1, mst_arch_lengths[i]);
//   }
//}
//@@@


////number of multicolvars used in MATRIX
//printf("getNumberOfNodeTypes() = %d on rank %d\n",getNumberOfNodeTypes(), comm.Get_rank());//@@@
////printf("getBaseMultiColvar(0) = %d on rank %d\n",getBaseMultiColvar(0), comm.Get_rank()); //@@@
////number of atoms used in the first multicolvar used in MATRIX
////printf("getNumberOfAtomsInGroup(0) = %d on rank %d\n",getNumberOfAtomsInGroup(0), comm.Get_rank());//@@@
////number of components in this multicolvar
//printf("getNumberOfComponents() = %d on rank %d\n",getNumberOfComponents(), comm.Get_rank());//@@@
////dynamic value: n * (n - 1) for asymmetric graph; n * (n - 1)/2 for symmetric 
////where n is the number of atoms after filtering
//printf("getNumberOfStoredValues() = %d on rank %d\n",getAdjacencyVessel()->getNumberOfStoredValues(), comm.Get_rank());//@@@

}

}
}

#endif
