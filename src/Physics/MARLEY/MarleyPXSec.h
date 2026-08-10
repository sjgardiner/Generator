#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__
//____________________________________________________________________________
/*!

\class    genie::MarleyPXSec

\brief    Computes differential cross sections by interfacing with
          the external low-energy neutrino interaction generator MARLEY.\n
          This is a concrete implementation of the XSecAlgorithmI interface.

\ref      TODO: add refs.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  July 18, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _MARLEY_CROSS_SECTION_H_
#define _MARLEY_CROSS_SECTION_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/MARLEY/MarleyInterface.h"

namespace genie {

class MarleyPXSec : public XSecAlgorithmI {

public:
  MarleyPXSec();
  MarleyPXSec(string config);
  virtual ~MarleyPXSec();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction* i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction* i) const;
  bool   ValidProcess    (const Interaction* i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfig (void);

  const MarleyInterface* fMARLEY;

};

}       // genie namespace

#endif // _MARLEY_CROSS_SECTION_H_
#endif // __GENIE_MARLEY_ENABLED__
