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

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Physics/MARLEY/MarleyInterface.h"
#include "Physics/MARLEY/MarleyPXSec.h"

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
MarleyPXSec::MarleyPXSec() : XSecAlgorithmI("genie::MarleyPXSec")
{

}
//____________________________________________________________________________
MarleyPXSec::MarleyPXSec(string config) :
  XSecAlgorithmI("genie::MarleyPXSec", config)
{

}
//____________________________________________________________________________
MarleyPXSec::~MarleyPXSec()
{

}
//____________________________________________________________________________
double MarleyPXSec::XSec(const Interaction* /*interaction*/,
  KinePhaseSpace_t /*kps*/) const
{
  // Dummy method to fulfill the requirements of the genie::XSecAlgorithmI
  // interface. All real differential cross sections happen "under the hood"
  // within MARLEY itself.
  return 0.;
}
//____________________________________________________________________________
double MarleyPXSec::Integral(const Interaction* in) const
{
  // PDG codes for the nuclear target and probe
  const InitialState& init_state = in->InitState();
  int tgt_pdg = init_state.TgtPdg();
  int probe_pdg = init_state.ProbePdg();

  TLorentzVector* temp_probe_p4 = init_state.GetProbeP4( kRfLab );
  TLorentzVector probe_p4 = *temp_probe_p4;
  delete temp_probe_p4;

  // Probe kinetic energy (MeV)
  double probe_KE = ( probe_p4.E() - probe_p4.M() ) / genie::units::MeV;

  marley::Generator* marley_gen = fMARLEY->GetMarleyGenerator();
  double tot_xsec = marley_gen->total_xs( probe_pdg, probe_KE, tgt_pdg );

  // Convert from MARLEY's units (MeV^(-2)) to GENIE's (GeV^(-2))
  tot_xsec /= genie::units::MeV * genie::units::MeV;

  return tot_xsec;
}
//____________________________________________________________________________
bool MarleyPXSec::ValidProcess(const Interaction * interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState& init_state = interaction->InitState();
  const ProcessInfo&  proc_info  = interaction->ProcInfo();

  if ( !proc_info.IsQuasiElastic() ) return false;

  int nuc = init_state.Tgt().HitNucPdg();
  int nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool prcok = proc_info.IsWeakCC() && ( (isP && isnub) || (isN && isnu) );
  if ( !prcok ) return false;

  return true;
}
//____________________________________________________________________________
void MarleyPXSec::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void MarleyPXSec::Configure(string config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void MarleyPXSec::LoadConfig(void)
{
  fMARLEY = dynamic_cast<const MarleyInterface*>( this->SubAlg("MarleyAlg") );
  assert( fMARLEY );
}
//____________________________________________________________________________

#endif // __GENIE_MARLEY_ENABLED__
