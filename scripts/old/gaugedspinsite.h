#ifndef __ITENSOR_SPINHI_H
#define __ITENSOR_SPINHI_H
#include "itensor/mps/siteset.h"

namespace itensor {

class GaugedSpinSite;

using GaugedSpin = BasicSiteSet<GaugedSpinSite>;

class GaugedSpinSite    {
    Index s;
    public:

    CustomSpinSite() { }

    CustomSpinSite(Index I) : s(I) { }
    
    CustomSpinSite(Args const& args = Args::global())
      }
