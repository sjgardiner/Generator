#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*!

\class    genie::MarleyDeExcitation

\brief    Interface to the MARLEY nuclear de-excitation model

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  April 27, 2026

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MARLEY_DEEXCITATION_H
#define _MARLEY_DEEXCITATION_H

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/MARLEY/MarleyInterface.h"

namespace genie {

class MarleyDeExcitation : public EventRecordVisitorI {

public :

  MarleyDeExcitation();
  MarleyDeExcitation(string config);
  ~MarleyDeExcitation();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord* event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry& config);
  void Configure(string config);

private:

  void LoadConfig(void);

  const genie::MarleyInterface* fMARLEY;

};

}      // genie namespace
#endif // _MARLEY_DEEXCITATION_H
#endif // __GENIE_MARLEY_ENABLED__
