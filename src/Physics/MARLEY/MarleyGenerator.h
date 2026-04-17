#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*!

\class    genie::MarleyGenerator

\brief    Simulate events using an interface to the external MARLEY
          generator for low-energy neutrino interactions

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  July 18, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MARLEY_GENERATOR_H
#define _MARLEY_GENERATOR_H

#include "TLorentzVector.h"

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Physics/MARLEY/MarleyInterface.h"

namespace marley {
  class Particle;
}

namespace genie {

class Interaction;

class MarleyGenerator : public EventRecordVisitorI {

public :

  MarleyGenerator();
  MarleyGenerator(string config);
  ~MarleyGenerator();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord* event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry& config);
  void Configure(string config);

  void AddMarleyParticle( GHepRecord* event,
    const HepMC3::GenParticle& part, int mom_index,
    GHepStatus_t status, const TLorentzVector& v4 ) const;

private:

  void LoadConfig(void);

  const genie::MarleyInterface* fMARLEY;

};

}      // genie namespace
#endif // _MARLEY_GENERATOR_H
#endif // __GENIE_MARLEY_ENABLED__
