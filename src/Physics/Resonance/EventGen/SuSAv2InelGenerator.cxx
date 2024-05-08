//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Resonance/EventGen/SuSAv2InelGenerator.h"

namespace {

  enum HadronicFSChoice { Unknown, DeltaRES, OtherRES, DIS };

}

//___________________________________________________________________________
genie::SuSAv2InelGenerator::SuSAv2InelGenerator()
  : genie::KineGeneratorWithCache( "genie::SuSAv2InelGenerator" )
{
}

//___________________________________________________________________________
genie::SuSAv2InelGenerator::SuSAv2InelGenerator( std::string config )
  : genie::KineGeneratorWithCache( "genie::SuSAv2InelGenerator", config )
{
}

//___________________________________________________________________________
genie::SuSAv2InelGenerator::~SuSAv2InelGenerator()
{
}

//___________________________________________________________________________
void genie::SuSAv2InelGenerator::ProcessEventRecord( GHepRecord* evrec ) const
{
  // Access the cross section algorithm for running thread (presumably the
  // SuSAv2 inelastic one)
  genie::RunningThreadInfo* rt_info = RunningThreadInfo::Instance();
  const genie::EventGeneratorI* evg = rt_info->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  // Make the remnant nucleus after removing the struck nucleon (already set)
  this->MakeNuclearRemnant( evrec );

  // Generate the outgoing lepton
  this->SelectLeptonKinematics( evrec );

  // Generate the outgoing hadrons
  this->MakeHadronicFinalState( evrec );
}

//___________________________________________________________________________
void genie::SuSAv2InelGenerator::SelectLeptonKinematics( GHepRecord* evrec )
  const
{
  // Get access to the interaction from within the event record
  genie::Interaction* interaction = evrec->Summary();
  interaction->SetBit( genie::kISkipProcessChk );

  // Get access to the random number generators
  genie::RandomGen* rnd = RandomGen::Instance();

  // Retrieve the probe 4-momentum in the lab frame. Delete the cloned version
  // to avoid a memory leak.
  TLorentzVector* temp_p4v = interaction->InitState().GetProbeP4();
  TLorentzVector p4v( *temp_p4v );
  delete temp_p4v;

  // Get the lab-frame probe energy from the input interaction
  double Ev = p4v.E();

  // Get the final-state lepton mass from the input interaction
  double ml = interaction->FSPrimLepton()->Mass();

  // Set the kinematic limits for sampling
  // TODO: revisit these and assign better values if needed
  double W_min = genie::constants::kNeutronMass + genie::constants::kPionMass;
   double W_max =genie::constants::kNeutronMass + Ev;// - Tl - ml;

  double cth_min = -1.;
  double cth_max =  1.;

  double Tl_min = 0.;
  double Tl_max =Ev - ml;

  // Calculate the maximum differential cross section or retrieve it
  double xsec_max = this->MaxXSec( evrec );

  double W, cth, Tl, Q2, xsec;
  TLorentzVector p4l;
  unsigned int iter = 0;
  bool accept = false;
  while ( true ) {
    ++iter;
    if( iter > genie::controls::kRjMaxIterations ) {
      LOG( "SuSAv2Inel", pWARN )
        << "*** Could not select a valid (W, cth, Tl) tuple after "
        << iter << " iterations";
      evrec->EventFlags()->SetBitNumber( genie::kKineGenErr, true );
      genie::exceptions::EVGThreadException exception;
      exception.SetReason( "Couldn't select kinematics" );
      exception.SwitchOnFastForward();
      throw exception;
    }

    // Currently sampled values of the kinematic variables of interest

    cth = cth_min + ( cth_max - cth_min ) * rnd->RndKine().Rndm();
    Tl = Tl_min + ( Tl_max - Tl_min ) * rnd->RndKine().Rndm();
    W = W_min + ( W_max - W_min ) * rnd->RndKine().Rndm();

    if ( W > genie::constants::kNeutronMass + Ev - Tl - ml ) continue;



    LOG( "SuSAv2Inel", pDEBUG ) << "Trying: W = " << W
      << ", cth = " << cth << ", Tl = " << Tl;

    // Update the interaction object with the sampled values
    interaction->KinePtr()->SetW( W );
    interaction->KinePtr()->SetKV( genie::kKVctl, cth );
    interaction->KinePtr()->SetKV( genie::kKVTl, Tl );

    // Compute the differential cross section for the current kinematics
   double xsec_ini = fXSecModel->XSec( interaction, genie::kPSTlctl );
     xsec=xsec_ini*W/pow(genie::constants::kNeutronMass,2);
    LOG( "SuSAv2Inel", pDEBUG ) << "xsec = " << xsec;

    // Do the usual check for rejection sampling of the kinematics
    double y = xsec_max * rnd->RndKine().Rndm();
    accept = (y < xsec);

    LOG( "SuSAv2Inel", pDEBUG ) << "xsec_max = " << xsec_max << ", y = " << y;

    // Double-check that the precomputed maximum differential cross section
    // wasn't exceeded by the current one (otherwise we have a rejection
    // sampling problem)
    this->AssertXSecLimits( interaction, xsec, xsec_max );

    // If the generated kinematics are accepted, then double-check the Q^2
    // limit before breaking out of the rejection sampling loop. Skip the
    // check if we're not going to accept the variables thrown to save time.
    if ( accept ) {

      // Select an azimuthal angle for the outgoing lepton (uniformly
      // distributed)
      double phi = 2. * genie::constants::kPi * rnd->RndKine().Rndm();

      // Build the outgoing lepton 4-momentum in the lab frame
      double El = ml + Tl;
      double pl = std::sqrt( std::max(0., El*El - ml*ml) );
      double sth = std::sqrt( std::max(0., 1. - cth*cth) );

      double pxl = pl * sth * std::cos( phi );
      double pyl = pl * sth * std::sin( phi );
      double pzl = pl * cth;

      // Rotate lepton 3-momentum vector so that the prior z direction is along
      // the probe direction
      TVector3 unit_probe_dir = p4v.Vect().Unit();
      TVector3 p3l( pxl, pyl, pzl );
      p3l.RotateUz( unit_probe_dir );

      p4l = TLorentzVector( p3l, El );

      // Compute the 4-momentum transfer and Q^2
      TLorentzVector q4 = p4v - p4l;
      Q2 = -1. * q4.Mag2();

      // Reset the 'trust' bits before exiting the loop
      interaction->ResetBit( kISkipProcessChk );
      interaction->ResetBit( kISkipKinematicChk );

      LOG( "SuSAv2Inel", pINFO ) << "Selected: W = " << W
        << ", cth = " << cth << ", Tl = " << Tl;

      break;
    }
  } // iterations of the rejection sampling loop

  // Store the differential cross section for the selected kinematics
  evrec->SetDiffXSec( xsec, kPSTlctl );

  // Compute x and y from W and Q^2
  // NOTE: The values are currently approximate because we assume an on-shell
  // initial nucleon at rest (the hadronic system has not been generated yet
  // since the SuSAv2 differential cross section is inclusive).
  // TODO: revisit and do a better job later
  double x = -1.;
  double y = -1.;
  double M = genie::constants::kNucleonMass;
  genie::utils::kinematics::WQ2toXY( Ev, M, W, Q2, x, y );

  // Store selected kinematics and clear running values
  interaction->KinePtr()->SetQ2( Q2, true );
  interaction->KinePtr()->SetW( W,  true );
  interaction->KinePtr()->Setx( x,  true );
  interaction->KinePtr()->Sety( y,  true );
  interaction->KinePtr()->SetFSLeptonP4( p4l );
  interaction->KinePtr()->ClearRunningValues();

  // Create a particle object for the outgoing lepton. Make it the daughter of
  // the probe and assign it the same 4-position in the event record.
  int pdgc = interaction->FSPrimLepton()->PdgCode();
  int mom_index = evrec->ProbePosition();
  TLorentzVector pos4v( *evrec->Probe()->X4() );
  evrec->AddParticle( pdgc, kIStStableFinalState, mom_index, -1, -1, -1,
    p4l, pos4v );
}

//___________________________________________________________________________
void genie::SuSAv2InelGenerator::MakeNuclearRemnant( genie::GHepRecord* evrec )
  const
{
  // Add the remnant nucleus (= initial nucleus - struck nucleon) to the event
  // record
  GHepParticle* target = evrec->TargetNucleus();
  GHepParticle* hit_nuc = evrec->HitNucleon();

  int Z = target->Z();
  int A = target->A();

  // Remove the struck nucleon
  A -= 1;
  if ( hit_nuc->Pdg() == kPdgProton ) Z -= 1;

  int remnant_pdg = genie::pdg::IonPdgCode( A, Z );

  const TLorentzVector& p4_hit_nuc = *( hit_nuc->P4() );
  const TLorentzVector& p4tgt = *( target->P4() );

  const TLorentzVector p4 = p4tgt - p4_hit_nuc;
  const TLorentzVector v4( 0., 0., 0., 0. );

  int mom_idx = evrec->TargetNucleusPosition();

  evrec->AddParticle( remnant_pdg, kIStStableFinalState, mom_idx, -1, -1,
    -1, p4, v4 );
}
//___________________________________________________________________________
void genie::SuSAv2InelGenerator::Configure( const genie::Registry& config )
{
  genie::Algorithm::Configure( config );
  this->LoadConfig();
}

//____________________________________________________________________________
void genie::SuSAv2InelGenerator::Configure( std::string config )
{
  genie::Algorithm::Configure( config );
  this->LoadConfig();
}

//____________________________________________________________________________
void genie::SuSAv2InelGenerator::LoadConfig()
{
  // Safety factor for the maximum differential cross section
  this->GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.25 );

  // Minimum energy for which max xsec would be cached, forcing explicit
  // calculation for lower eneries
  this->GetParamDef( "Cache-MinEnergy", fEMin, 0.5 );

  // Load Wcut used in DIS/RES join scheme
  this->GetParam( "Wcut", fWcut );

  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  this->GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. );
  assert( fMaxXSecDiffTolerance >= 0. );

  // Get sub-algorithm that handles preparation of RES final states
  fRESHadronGenerator = dynamic_cast< const EventRecordVisitorI* >(
    this->SubAlg( "RESHadronGenerator" )
  );
  assert( fRESHadronGenerator );

  // Get sub-algorithm that handles preparation of DIS final states
  fDISHadronGenerator = dynamic_cast< const EventRecordVisitorI* >(
    this->SubAlg( "DISHadronGenerator" )
  );
  assert( fDISHadronGenerator );
}

//____________________________________________________________________________
double genie::SuSAv2InelGenerator::ComputeMaxXSec(
  const genie::Interaction* interaction ) const
{
  // Computes the maximum differential cross section in the requested phase
  // space. This method overloads the KineGeneratorWithCache::ComputeMaxXSec()
  // method. The value is cached at a circular cache branch for retrieval
  // during subsequent event generation. The computed max differential cross
  // section does not need to be the exact maximum. The number used in the
  // rejection method will be scaled up by a safety factor. But this needs to
  // be fast -- do not use a very fine grid.

  double Ev = interaction->InitState().ProbeE( genie::kRfLab );
  double ml = interaction->FSPrimLepton()->Mass();

  // Start with a guess of where the maximum xsec lies in phase space
  genie::PDGLibrary* pdg_library = genie::PDGLibrary::Instance();
  double mDelta = pdg_library->Find( genie::kPdgP33m1232_DeltaP )->Mass();

  double W = mDelta;
  if (Ev<0.52){W=0.468*Ev+0.9595;}
  double Tl = 0.9516*Ev - 0.3228;
  double cth = std::min( 1., 0.5395*std::exp(-1.984*Ev) + 0.9956 );

  interaction->KinePtr()->SetW( W );
  interaction->KinePtr()->SetKV( genie::kKVTl, Tl );
  interaction->KinePtr()->SetKV( genie::kKVctl, cth );

  double max_xsec = fXSecModel->XSec( interaction, kPSTlctl );
  LOG( "SuSAv2Inel", pDEBUG ) << "xsec(W= " << W << ", Tl= " << Tl
    << ", cth = " << cth << ") = " << max_xsec;

  // Step around the guess in all directions with a decreasing step size,
  // move on each iteration towards the maximum value
  double step_W = 0.5;
  double step_Tl = 0.5;
  double step_cth = 0.5;
  const int MAX_STEPS = 1000;
  const int MIN_STEP_W = 0.01;
  int step_iter = 0;

  while ( step_iter < MAX_STEPS && step_W > MIN_STEP_W ) {
    std::vector< double > Ws = { W - step_W, W, W + step_W };
    std::vector< double > Tls = { Tl - step_Tl, Tl, Tl + step_Tl };
    std::vector< double > cths = { cth - step_cth, cth, cth + step_cth };

    std::vector< genie::Interaction > interaction_vec;
    for ( const double& cur_W : Ws ) {
      for ( const double& cur_Tl : Tls ) {
        for ( const double& cur_cth : cths ) {
          interaction_vec.emplace_back( *interaction );
          auto& last_interaction = interaction_vec.back();
          last_interaction.KinePtr()->SetW( cur_W );
          last_interaction.KinePtr()->SetKV( genie::kKVTl, cur_Tl );
          last_interaction.KinePtr()->SetKV( genie::kKVctl, cur_cth );
        }
      }
    }

    const genie::Interaction* max_inter = interaction;
    for ( const auto& test_inter : interaction_vec ) {
      double xsec = fXSecModel->XSec( &test_inter, kPSTlctl );
      if ( xsec > max_xsec ) {
        LOG( "SuSAv2Inel", pDEBUG ) << "xsec(W= " << W << ", Tl= " << Tl
          << ", cth = " << cth << ") = " << xsec;
        max_inter = &test_inter;
        max_xsec = xsec;
        W = test_inter.Kine().W();
        Tl = test_inter.Kine().GetKV( genie::kKVTl );
        cth = test_inter.Kine().GetKV( genie::kKVctl );
      }
    }

    // If we found a phase-space point with a higher differential cross
    // section, then just continue the loop. If the current choice is still
    // best, decrease the step sizes for the next attempt.
    if ( max_inter == interaction ) {
      step_W /= 5.;
      step_Tl /= 5.;
      step_cth /= 5.;
    }
  }

  LOG( "SuSAv2Inel", pDEBUG ) << "max_xsec(W= " << W << ", Tl= " << Tl
    << ", cth = " << cth << ") = " << max_xsec;

  // Apply safety factor
  max_xsec *= fSafetyFactor;

  return max_xsec;
}
//___________________________________________________________________________
void genie::SuSAv2InelGenerator::MakeHadronicFinalState( GHepRecord* evrec )
  const
{
  HadronicFSChoice my_choice = HadronicFSChoice::Unknown;

  // The SuSAv2 calculation is for the inclusive cross section, so we adopt
  // an ad hoc procedure for deciding whether to generate a final state
  // based on GENIE's native RES or DIS treatment.
  genie::Interaction* interaction = evrec->Summary();
  double Q2 = interaction->Kine().Q2( true );
  double W = interaction->Kine().W( true );
  double Ev = interaction->InitState().ProbeE( genie::kRfLab );

  // If Q^2 or W is above one of these thresholds, unconditionally treat the
  // event as DIS
  // TODO: remove hard-coded kinematic limits in this function (make them
  // configurable via XML)
  // TODO: also adjust to respect existing W cut in CommonParam.xml
  if ( Q2 > 3.0 || W > 2.1 ) {
    my_choice = HadronicFSChoice::DIS;
  }
  // If W is below this threshold, assume pure RES through the Delta resonance
  else if ( W < 1.4 ) {
    my_choice = DeltaRES;
    interaction->ExclTagPtr()->SetResonance( genie::EResonance::kP33_1232 );
    fRESHadronGenerator->ProcessEventRecord( evrec );
  }
  else {
    // We are in the intermediate region of W, so consider a mixture of RES and
    // DIS interactions
    double ratio_other_RES = 0.085 * Ev + 0.045;
    double ratio_DIS = std::max( 0., -0.053 * Ev + 0.5183 );
    double ratio_Delta_RES = std::max( 0., 1. - ratio_other_RES - ratio_DIS );
    double sum_ratios = ratio_other_RES + ratio_DIS + ratio_Delta_RES;

    // Get access to the random number generators
    genie::RandomGen* rnd = RandomGen::Instance();

    // Choose a channel based on the relative probabilities
    double rand = sum_ratios * rnd->RndKine().Rndm();
    if ( rand <= ratio_Delta_RES ) {
      my_choice = HadronicFSChoice::DeltaRES;
    }
    else if ( rand <= ratio_Delta_RES + ratio_DIS ) {
      my_choice = HadronicFSChoice::DIS;
    }
    else {
      my_choice = HadronicFSChoice::OtherRES;
    }
  }

  // Delegate handling of the hadronic final state to the appropriate
  // sub-algorithm. In the case of "other RES," first choose the resonance to
  // use
  switch ( my_choice ) {
    case HadronicFSChoice::DIS: {

      // Switch the process information, etc. so that this event is labeled as
      // DIS (default for this workflow is RES)
      this->SwitchToDISEvent( evrec );

      // Delegate further event handling to the chosen DIS generator
      fDISHadronGenerator->ProcessEventRecord( evrec );
      break;
    }
    case HadronicFSChoice::DeltaRES: {
      interaction->ExclTagPtr()->SetResonance( genie::EResonance::kP33_1232 );
      fRESHadronGenerator->ProcessEventRecord( evrec );
      break;
    }
    case HadronicFSChoice::OtherRES: {
      // TODO: implement this part!
      // FIXME: first attempt
      double ratio_other_RES = 0.085 * Ev + 0.045;
      // Get access to the random number generators
      genie::RandomGen* rnd = RandomGen::Instance();
      genie::Resonance_t other_res;
      int other_res_norm = int( ratio_other_RES / 0.1);
      int other_res_choice = int( (other_res_norm + 1) * rnd->RndKine().Rndm());
      switch (other_res_choice) {
        case 0: other_res = genie::EResonance::kP11_1440; break;
        case 1: other_res = genie::EResonance::kP33_1600; break;
        case 2: other_res = genie::EResonance::kS11_1650; break;
        case 3: other_res = genie::EResonance::kD15_1675; break;
        case 4: other_res = genie::EResonance::kD33_1700; break;
        case 5: other_res = genie::EResonance::kP11_1710; break;
        case 6: other_res = genie::EResonance::kP31_1910; break;
        case 7: other_res = genie::EResonance::kP33_1920; break;
        case 8: other_res = genie::EResonance::kF37_1950; break;
        default:{
          LOG( "SuSAv2Inel", pERROR )
            << "*** Failed to select valid *other* hadronic final state";
          evrec->EventFlags()->SetBitNumber( genie::kHadroSysGenErr, true );
          genie::exceptions::EVGThreadException exception;
          exception.SetReason( "Couldn't select *other* hadronic final state" );
          exception.SwitchOnFastForward();
          throw exception;
        }
      }
      interaction->ExclTagPtr()->SetResonance( other_res );
      fRESHadronGenerator->ProcessEventRecord( evrec );
      break;
    }
    default: {
      LOG( "SuSAv2Inel", pERROR )
        << "*** Failed to select valid hadronic final state";
      evrec->EventFlags()->SetBitNumber( genie::kHadroSysGenErr, true );
      genie::exceptions::EVGThreadException exception;
      exception.SetReason( "Couldn't select hadronic final state" );
      exception.SwitchOnFastForward();
      throw exception;
    }
  }

}
//___________________________________________________________________________
void genie::SuSAv2InelGenerator::SwitchToDISEvent( GHepRecord* evrec ) const
{
  // Relabel the event (default for this generator is RES) to DIS in the process
  // info. Keep the interaction type (CC, etc.) unchanged
  genie::Interaction* interaction = evrec->Summary();
  genie::ProcessInfo* proc = interaction->ProcInfoPtr();
  genie::InteractionType_t inter_type = proc->InteractionTypeId();
  proc->Set( genie::kScDeepInelastic, inter_type );

  // Set the hit quark. The recipe depends on the probe and interaction type
  int probe_pdg = interaction->InitState().ProbePdg();

  bool is_cc = proc->IsWeakCC();
  bool is_nc_or_em = proc->IsWeakNC();

  // For CC events, charge conservation ensures that only one valence quark
  // species can participate, so just pick the relevant one and call it good.
  int hit_q_pdg = kPdgDQuark;
  genie::Target* tgt = interaction->InitState().TgtPtr();

  if ( is_cc ) {
    // For antineutrinos, we want the other option
    if ( probe_pdg < 0 ) hit_q_pdg = kPdgUQuark;
  }
  else if ( !is_nc_or_em ) {
    LOG( "SuSAv2Inel", pWARN )
      << "*** Could not select a valid hit quark for the interaction";
    evrec->EventFlags()->SetBitNumber( genie::kKineGenErr, true );
    genie::exceptions::EVGThreadException exception;
    exception.SetReason( "Couldn't select hit quark" );
    exception.SwitchOnFastForward();
    throw exception;
  }
  else {
    // For NC and EM processes, set the hit quark considering only valence
    // quarks in the struck nucleon with equal probabilities
    // TODO: do something better!
    int hit_nuc_pdg = tgt->HitNucPdg();

    // Set the up quark probability for a neutron first
    double up_prob = 1. / 3.;
    if ( hit_nuc_pdg == genie::kPdgProton ) {
      up_prob = 2. / 3.;
    }

    // Throw a random number to pick an up or down quark
    genie::RandomGen* rnd = RandomGen::Instance();
    double rand = rnd->RndKine().Rndm();

    // Already set as the default value above
    //hit_q_pdg = kPdgDQuark;
    if ( rand <= up_prob ) hit_q_pdg = kPdgUQuark;
  }

  // Set the hit quark PDG code in the interaction
  tgt->SetHitQrkPdg( hit_q_pdg );

  // We only consider valence quarks for now
  // TODO: do something better
  tgt->SetHitSeaQrk( false );
}
