/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "NeighborListMPI.h"
#include "Vector.h"
#include "Pbc.h"
#include "AtomNumber.h"
#include "Tools.h"
#include "Communicator.h"
#include "mpi.h"
#include "OpenMP.h"
#include <vector>
#include <algorithm>

namespace PLMD{
using namespace std;

NeighborListMPI::NeighborListMPI(const vector<AtomNumber>& list0, const vector<AtomNumber>& list1,
                           const bool& do_pair, const bool& do_pbc, const Pbc& pbc,
                           const double& distance, const unsigned& stride, PLMD::Communicator* comm):
                           reduced(false),
                           do_pair_(do_pair), do_pbc_(do_pbc), pbc_(&pbc),
                           distance_(distance), stride_(stride), comm_(comm)
{
// store full list of atoms needed
 //fullatomlist_=list1;
 //fullatomlist_.push_front(list0);
 //nlist0_=1;
 fullatomlist_=list0;
 fullatomlist_.insert(fullatomlist_.end(),list1.begin(),list1.end());
 nlist0_=list0.size();
 nlist1_=list1.size();
 twolists_=true;
 if(!do_pair){
  nallpairs_=nlist0_*nlist1_;
 }else{
  plumed_assert(nlist0_==nlist1_);
  nallpairs_=nlist0_;
 }
 initialize();
 lastupdate_=0;
 
 //the following define the mpi type for pair<unsigned,unsigned>
 const int nitems=2;
 const int blocklengths[2]={1,1};
 MPI_Datatype types[2]={MPI_UNSIGNED,MPI_UNSIGNED};
 MPI_Aint offsets[2];
 offsets[0]=offsetof(unsigned_pair,first);
 offsets[1]=offsetof(unsigned_pair,second);
 MPI_Type_create_struct(nitems, blocklengths, offsets, &types[0],&mpi_unsigned_pair);
 MPI_Type_commit(&mpi_unsigned_pair);
}

NeighborListMPI::NeighborListMPI(const vector<AtomNumber>& list0, const bool& do_pbc,
                           const Pbc& pbc, const double& distance,
                           const unsigned& stride, PLMD::Communicator* comm): 
                           reduced(false),
                           do_pbc_(do_pbc), pbc_(&pbc),
                           distance_(distance), stride_(stride), comm_(comm){
 fullatomlist_=list0;
 nlist0_=list0.size();
 twolists_=false;
 nallpairs_=nlist0_*(nlist0_-1)/2;
 initialize();
 lastupdate_=0;
 //the following define the mpi type for pair<unsigned,unsigned>
 const int nitems=2;
 const int blocklengths[2]={1,1};
 MPI_Datatype types[2]={MPI_UNSIGNED,MPI_UNSIGNED};
 MPI_Aint offsets[2];
 offsets[0]=offsetof(unsigned_pair,first);
 offsets[1]=offsetof(unsigned_pair,second);
 MPI_Type_create_struct(nitems, blocklengths, offsets, &types[0],&mpi_unsigned_pair);
 MPI_Type_commit(&mpi_unsigned_pair);
}

void NeighborListMPI::initialize() {
 neighbors_.clear();
 for(unsigned int i=0;i<nallpairs_;++i){
   neighbors_.push_back(getIndexPair(i));
 }
}

vector<AtomNumber>& NeighborListMPI::getFullAtomList() {
 return fullatomlist_;
}

pair<unsigned,unsigned> NeighborListMPI::getIndexPair(unsigned ipair) {
 pair<unsigned,unsigned> index;
 if(twolists_ && do_pair_){
  index=pair<unsigned,unsigned>(ipair,ipair+nlist0_);
 }else if (twolists_ && !do_pair_){
  index=pair<unsigned,unsigned>(ipair/nlist1_,ipair%nlist1_+nlist0_);
 }else if (!twolists_){
  unsigned ii = nallpairs_-1-ipair;
  unsigned  K = unsigned(floor((sqrt(double(8*ii+1))+1)/2));
  unsigned jj = ii-K*(K-1)/2;
  index=pair<unsigned,unsigned>(nlist0_-1-K,nlist0_-1-jj);
 }
 return index;
}

void NeighborListMPI::update(const vector<Vector>& positions) {
 neighbors_.clear();
 vector<unsigned_pair> neighbors_local;
 //unsigned nt=OpenMP::getNumThreads();
 unsigned stride=comm_->Get_size();
 unsigned rank=comm_->Get_rank();
 const double d2=distance_*distance_;
// check if positions array has the correct length 
 plumed_assert(positions.size()==fullatomlist_.size());
// delete the multithreading for thread safety considerations
//#pragma omp parallel num_threads(nt)
//{
// std::vector<unsigned_pair> omp_neighbors;
//#pragma omp for nowait
 for(unsigned int i=rank;i<nallpairs_;i+=stride){
   pair<unsigned,unsigned> index=getIndexPair(i);
   unsigned index0=index.first;
   unsigned index1=index.second;
   Vector distance;
   if(do_pbc_){
    distance=pbc_->distance(positions[index0],positions[index1]);
   } else {
    distance=delta(positions[index0],positions[index1]);
   }
   double value=modulo2(distance);
   //if(value<=d2) {omp_neighbors.push_back(index);} 
   if(value<=d2) {neighbors_local.push_back(index);} 
 }
//#pragma omp critical
// neighbors_local.insert(
//         neighbors_local.end(),omp_neighbors.begin(),omp_neighbors.end());
//}

//printf("found %d pairs on rank %d\n",neighbors_local.size(),comm_->Get_rank());
 //collect neighbors_local from all ranks
 //TODO use plumed wrapped gatherv instead
 int* counts=(int*)malloc(stride*sizeof(int));
 int* displ=(int*)malloc(stride*sizeof(int));
 for(int i=0;i<stride;++i) {
     if(rank==i)
         counts[i]=neighbors_local.size();
     MPI_Bcast(&counts[i],1,MPI_INT,i,comm_->Get_world().Get_comm());
 }
 displ[0]=0;
 for(int i=1;i<stride;++i) {
     displ[i]=displ[i-1]+counts[i-1];
 }
 int total=displ[stride-1]+counts[stride-1]; //total # of pairs
 unsigned_pair* neighbors_list=
    (unsigned_pair*) malloc(sizeof(unsigned_pair)*total);
    MPI_Allgatherv(&neighbors_local[0],neighbors_local.size(),mpi_unsigned_pair,
            &neighbors_list[0],counts,displ,mpi_unsigned_pair,comm_->Get_world().Get_comm());
 neighbors_=std::vector<unsigned_pair>(neighbors_list,neighbors_list+total);

//printf("neighbors_local has %d atoms, neighbors_list has %d atoms on rank %d\n",neighbors_local.size(),total, comm_->Get_rank());//@@@

 free(counts); free(displ); free(neighbors_list);
 setRequestList();
}

void NeighborListMPI::setRequestList() {
 requestlist_.clear();
 for(unsigned int i=0;i<size();++i){
  requestlist_.push_back(fullatomlist_[neighbors_[i].first]);
  requestlist_.push_back(fullatomlist_[neighbors_[i].second]);
 }
 Tools::removeDuplicates(requestlist_);
 reduced=false;
}

vector<AtomNumber>& NeighborListMPI::getReducedAtomList() {
 if(!reduced)for(unsigned int i=0;i<size();++i){
  unsigned newindex0=0,newindex1=0;
  // neighbors_ now stores the indexes in fullatomlist_
  // we need to modify it so that it stores the indexes in requestlist_
  AtomNumber index0=fullatomlist_[neighbors_[i].first];
  AtomNumber index1=fullatomlist_[neighbors_[i].second];
// I exploit the fact that requestlist_ is an ordered vector
  vector<AtomNumber>::iterator p;
  p = std::find(requestlist_.begin(), requestlist_.end(), index0); plumed_assert(p!=requestlist_.end()); newindex0=p-requestlist_.begin();
  p = std::find(requestlist_.begin(), requestlist_.end(), index1); plumed_assert(p!=requestlist_.end()); newindex1=p-requestlist_.begin();
  neighbors_[i]=pair<unsigned,unsigned>(newindex0,newindex1);
 }
 reduced=true;
 return requestlist_;
}

unsigned NeighborListMPI::getStride() const {
 return stride_;
}

unsigned NeighborListMPI::getLastUpdate() const {
 return lastupdate_;
}

void NeighborListMPI::setLastUpdate(unsigned step) {
 lastupdate_=step;
}

unsigned NeighborListMPI::size() const {
 return neighbors_.size();
}

 pair<unsigned,unsigned> NeighborListMPI::getClosePair(unsigned i) const {
 return neighbors_[i];
}

vector<unsigned> NeighborListMPI::getNeighbors(unsigned index) {
 vector<unsigned> neighbors;
 for(unsigned int i=0;i<size();++i){
  if(neighbors_[i].first==index)  neighbors.push_back(neighbors_[i].second);
  if(neighbors_[i].second==index) neighbors.push_back(neighbors_[i].first);
 }
 return neighbors;
}

}
