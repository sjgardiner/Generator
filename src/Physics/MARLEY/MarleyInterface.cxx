#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// GENIE includes
#include "Physics/MARLEY/MarleyInterface.h"

// MARLEY includes
#include "marley/RootJSONConfig.hh"

using namespace genie;

//____________________________________________________________________________
MarleyInterface::MarleyInterface() : Algorithm( "genie::MarleyInterface" )
{

}
//____________________________________________________________________________
MarleyInterface::MarleyInterface(string config) :
Algorithm( "genie::MarleyInterface", config )
{

}
//____________________________________________________________________________
MarleyInterface::~MarleyInterface()
{

}
//____________________________________________________________________________
void MarleyInterface::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void MarleyInterface::Configure( string config )
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void MarleyInterface::LoadConfig(void)
{
  // Get configuration file name for initializing MARLEY
  std::string config_file_name;
  GetParam( "ConfigFileName", config_file_name ) ;

  // Initialize a new marley::Generator object
  marley::RootJSONConfig jc( config_file_name );
  fMarleyGenerator = jc.create_generator();
}
//____________________________________________________________________________
marley::Generator* MarleyInterface::GetMarleyGenerator() const {
  return &fMarleyGenerator;
}
//____________________________________________________________________________

#endif // __GENIE_MARLEY_ENABLED__
