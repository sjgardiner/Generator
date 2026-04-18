#ifdef __CINT__
#include "Framework/Conventions/GBuild.h"

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;

#ifdef __GENIE_MARLEY_ENABLED__
#pragma link C++ class genie::MarleyInterface;
#pragma link C++ class genie::MarleyGenerator;
#pragma link C++ class genie::MarleyPXSec;
#pragma link C++ class genie::MarleyDeExcitation;
#endif

#endif
