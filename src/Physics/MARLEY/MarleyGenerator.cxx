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
#include <array>

// GENIE includes
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/EventGen/HepMC3Converter.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/MARLEY/MarleyGenerator.h"
#include "Physics/MARLEY/MarleyInterface.h"

#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

namespace {
  constexpr int HEPMC3_FINAL_STATE_STATUS = 1;
}

//___________________________________________________________________________
MarleyGenerator::MarleyGenerator() :
  EventRecordVisitorI( "genie::MarleyGenerator" )
{

}
//___________________________________________________________________________
MarleyGenerator::MarleyGenerator(string config) :
  EventRecordVisitorI( "genie::MarleyGenerator", config )
{

}
//___________________________________________________________________________
MarleyGenerator::~MarleyGenerator()
{

}
//___________________________________________________________________________
void MarleyGenerator::ProcessEventRecord(GHepRecord* event) const
{
  genie::Interaction* inter = event->Summary();

  const InitialState& init_state = inter->InitState();
  int probe_pdg = init_state.ProbePdg();
  int tgt_pdg = init_state.TgtPdg();

  TLorentzVector* temp_probe_p4 = init_state.GetProbeP4( kRfLab );
  TLorentzVector probe_p4 = *temp_probe_p4;
  delete temp_probe_p4;

  // Probe kinetic energy (MeV)
  double probe_KE = ( probe_p4.E() - probe_p4.M() ) / genie::units::MeV;

  // Build a unit vector in the probe's direction of motion
  TVector3 dir = probe_p4.Vect().Unit();
  std::array< double, 3 > probe_dir = { dir.X(), dir.Y(), dir.Z() };

  // Get the 4-position of the interaction vertex
  TLorentzVector v4( *event->Probe()->X4() );

  // Create a new MARLEY event for the given initial state
  // TODO: seed the MARLEY generator object appropriately
  marley::Generator* marley_gen = fMARLEY->GetMarleyGenerator();
  std::shared_ptr< HepMC3::GenEvent > marley_evt
    = marley_gen->create_event( probe_pdg, probe_KE, tgt_pdg, probe_dir );

  const int probe_idx = event->ProbePosition();
  const int tgt_idx = event->TargetNucleusPosition();

  // Add the final-state particles to the event record
  for ( const auto& part : marley_evt->particles() ) {

    // Skip everything except the final-state particles
    // TODO: revisit this so we can store the full MARLEY event history
    if ( part->status() != HEPMC3_FINAL_STATE_STATUS ) continue;

    int pdg = part->pid();
    int mom_index = tgt_idx;
    if ( genie::pdg::IsLepton(pdg) ) mom_index = probe_idx;

    this->AddMarleyParticle( event, *part, mom_index, kIStStableFinalState, v4 );
  }

  // Set the final-state lepton's polarization
  genie::utils::SetPrimaryLeptonPolarization( event );
}
//___________________________________________________________________________
void MarleyGenerator::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyGenerator::Configure(string config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyGenerator::LoadConfig(void)
{
  fMARLEY = dynamic_cast<const MarleyInterface*>( this->SubAlg("MarleyAlg") );
}
//___________________________________________________________________________
void MarleyGenerator::AddMarleyParticle( GHepRecord* event,
  const HepMC3::GenParticle& part, int mom_index,
  GHepStatus_t status, const TLorentzVector& v4 ) const
{
  // Get the particle's 4-momentum, and convert to using GeV instead of MeV
  const HepMC3::FourVector& mom4 = part.momentum();
  TLorentzVector p4( mom4.px() * genie::units::MeV,
    mom4.py() * genie::units::MeV, mom4.pz() * genie::units::MeV,
    mom4.e() * genie::units::MeV );

  event->AddParticle( part.pid(), status, mom_index,
    -1, -1, -1, p4, v4 );
}

#endif // __GENIE_MARLEY_ENABLED__
