#include "itensor/all.h"
#include "ladder.h"

using namespace itensor;
using Schwinger = MixedSiteSet<LadderSite, SpinHalfSite>;
int schwinger(const char *inputfile, Args const& args = Args::global());
