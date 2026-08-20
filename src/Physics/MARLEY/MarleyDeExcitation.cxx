#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_MARLEY_ENABLED__

//____________________________________________________________________________
/*
 Copyright (c) 2003-2026, The GENIE Collaboration
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
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/HepMC3Converter.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/MARLEY/MarleyDeExcitation.h"
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
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/NucleusDecayer.hh"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
MarleyDeExcitation::MarleyDeExcitation() :
  EventRecordVisitorI( "genie::MarleyDeExcitation" )
{

}
//___________________________________________________________________________
MarleyDeExcitation::MarleyDeExcitation(string config) :
  EventRecordVisitorI( "genie::MarleyDeExcitation", config )
{

}
//___________________________________________________________________________
MarleyDeExcitation::~MarleyDeExcitation()
{

}
//___________________________________________________________________________
void MarleyDeExcitation::ProcessEventRecord(GHepRecord* event) const
{
  static genie::HepMC3Converter hepmc3_conv;
  static marley::NucleusDecayer marley_deex;

  // Convert the input GHepRecord* to EventRecord* so that it plays well
  // with the HepMC3 converter interface
  auto* ev_rec = dynamic_cast< genie::EventRecord* >( event );

  // Convert the input event from GHEP to HepMC3 format (MARLEY uses the latter)
  std::shared_ptr< HepMC3::GenEvent > marley_event
    = hepmc3_conv.ConvertToHepMC3( *ev_rec );

  // Run the de-excitation model on the event
  try {
    auto* marley_gen = fMARLEY->GetMarleyGenerator();
    marley_deex.process_event( *marley_event, *marley_gen );
  }
  // If MARLEY runs into a problem, convert its exception into
  // a GENIE EVGThreadException and rethrow while backing up the
  // event generation thread
  catch ( const marley::Error& err ) {
    ev_rec->EventFlags()->SetBitNumber( genie::kDecayErr, true );
    genie::exceptions::EVGThreadException exception;
    exception.SetReason( err.what() );
    exception.SwitchOnStepBack();
    exception.SetReturnStep( 0 );
    throw exception;
  }

  // MARLEY will append new particles in the event, so start the loop
  // over the new additions by skipping all pre-existing particles
  // NOTE: This relies on the HepMC3 particles being listed in the same
  // order as in the original GHepRecord. If that changes, then
  // we will need to revisit this code.
  int num_original_genie_particles = event->GetEntries();
  const auto& marley_particles = marley_event->particles();
  int num_marley_particles = marley_particles.size();

  // NucleusDecayer::process_event() restores the event's original momentum
  // unit before returning (it saves old_p4_unit up front and calls
  // event.set_units(old_p4_unit, ...) at the end), and that call rescales
  // every particle currently in the event -- including the ones MARLEY just
  // added while working internally in MeV -- to match. Since GENIE created
  // this event in GeV, by this point it should already be back in GeV.
  // Check the event's actual declared unit rather than assuming MeV: a
  // blind MeV->GeV conversion here double-converts already-correct GeV
  // values (confirmed empirically -- see debug trace from event 20).
  double p4_scale = 1.0;
  if ( marley_event->momentum_unit() == HepMC3::Units::MEV ) {
    p4_scale = genie::units::MeV;
  }

  for ( int p = num_original_genie_particles; p < num_marley_particles; ++p ) {
    // Convert the HepMC3 status code to the native GENIE one
    auto& mar_part = marley_particles.at( p );
    genie::GHepStatus_t status
      = hepmc3_conv.GetGHepParticleStatus( mar_part->status() );

    // Get the particle 4-momentum, already in GENIE's native GeV units
    const HepMC3::FourVector& mom4 = mar_part->momentum();
    TLorentzVector p4( mom4.px() * p4_scale, mom4.py() * p4_scale,
      mom4.pz() * p4_scale, mom4.e() * p4_scale );

    // Get the parent particle via the production vertex. We will use
    // it to assign a 4-position.
    int mommy_id = mar_part->production_vertex()->particles_in().front()->id();
    int mommy_index = mommy_id - 1;
    GHepParticle* mommy = event->Particle( mommy_index );

    event->AddParticle( mar_part->pid(), status, mommy_index, -1, -1, -1,
      p4, *mommy->X4() );
  }

}
//___________________________________________________________________________
void MarleyDeExcitation::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyDeExcitation::Configure(string config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//___________________________________________________________________________
void MarleyDeExcitation::LoadConfig(void)
{
  fMARLEY = dynamic_cast<const MarleyInterface*>( this->SubAlg("MarleyAlg") );
}

#endif // __GENIE_MARLEY_ENABLED__
