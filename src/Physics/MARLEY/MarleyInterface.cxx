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

// Standard library includes
#include <cstdlib>

// GENIE includes
#include "Physics/MARLEY/MarleyInterface.h"
#include "Framework/Numerical/RandomGen.h"

// MARLEY includes
#include "marley/JSONConfig.hh"

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
  std::string full_path = std::getenv( "GENIE" );
  full_path += "/data/evgen/marley/" + config_file_name;
  marley::JSONConfig jc( full_path );
  fMarleyGenerator = jc.create_generator();

  // Seed the MARLEY random number generator using the GENIE seed
  genie::RandomGen* rnd = RandomGen::Instance();
  long int genie_seed = rnd->GetSeed();
  uint_fast64_t marley_seed = static_cast< uint_fast64_t >( genie_seed );
  fMarleyGenerator.reseed( marley_seed );
}
//____________________________________________________________________________
marley::Generator* MarleyInterface::GetMarleyGenerator() const {
  return &fMarleyGenerator;
}
//____________________________________________________________________________

#endif // __GENIE_MARLEY_ENABLED__
