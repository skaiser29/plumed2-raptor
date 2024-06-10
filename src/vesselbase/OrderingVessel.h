/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_vesselbase_OrderingVessel_h
#define __PLUMED_vesselbase_OrderingVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "ValueVessel.h"
#include "StoreDataVessel.h"
#include "core/Value.h"

namespace PLMD {
namespace vesselbase {

class OrderingVessel : public ValueVessel {
private:
  StoreDataVessel* mydata;
public:
  static void registerKeywords( Keywords& keys );
  explicit OrderingVessel( const VesselOptions& da );
  void resize() override;
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const override {}
  void finish( const std::vector<double>& buffer ) override;
  virtual bool compare( const double&, const double& )=0;
};

}
}
#endif

