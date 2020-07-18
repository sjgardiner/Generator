#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*!

\class    genie::MarleyInterface

\brief    Algorithm that provides a lightweight wrapper for a marley::Generator
          object.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  July 18, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _MARLEY_INTERFACE_H_
#define _MARLEY_INTERFACE_H_

// Standard library includes
#include <memory>

// GENIE includes
#include "Framework/Algorithm/Algorithm.h"

// MARLEY includes
#include "marley/Generator.hh"

namespace genie {

class MarleyInterface : public Algorithm {

public:

  MarleyInterface();
  MarleyInterface(string config);
  virtual ~MarleyInterface();

  // Override the Algorithm::Configure methods to load configuration data to
  // private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  marley::Generator* GetMarleyGenerator();

private:

  void LoadConfig (void);

  marley::Generator fMarleyGenerator;

};

} // genie namespace

#endif // _MARLEY_INTERFACE_H

#endif // __GENIE_MARLEY_ENABLED__
