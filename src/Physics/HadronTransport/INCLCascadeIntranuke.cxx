#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

//---------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

// GENIE
#include "INCLCascadeIntranuke.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Algorithm/AlgConfigPool.h"

// INCL++
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"


#include "G4INCLPauliBlocking.hh"

#include "G4INCLCrossSections.hh"
#include "G4INCLDecayAvatar.hh"

#include "G4INCLPhaseSpaceGenerator.hh"

#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"

#include "G4INCLINuclearPotential.hh"

#include "G4INCLCoulombDistortion.hh"

#include "G4INCLClustering.hh"

#include "G4INCLIntersection.hh"

#include "G4INCLBinaryCollisionAvatar.hh"

#include "G4INCLCascadeAction.hh"
#include "G4INCLAvatarDumpAction.hh"

#include "G4INCLClusterDecay.hh"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

#include "Physics/NuclearState/NucleusGenI.h"



// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif

//#ifdef HAS_BOOST_DATE_TIME
//#include <boost/date_time/posix_time/posix_time.hpp>
//namespace bpt = boost::posix_time;
//#endif

#ifdef HAS_BOOST_TIMER
#include <boost/timer/timer.hpp>
namespace bt = boost::timer;
#endif


// --------------------------------------Include for GENIE---------------------
// GENIE

#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include "Physics/HadronTransport/G4INCLGENIEAvatar.h"
#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"

// ROOT
#include "TSystem.h"

using namespace genie;
using namespace genie::utils;
using namespace G4INCL;
using std::ostringstream;
using namespace std;

INCLCascadeIntranuke::INCLCascadeIntranuke() :
  EventRecordVisitorI("genie::INCLCascadeIntranuke"),
  theINCLConfig(0), theINCLModel(0), theDeExcitation(0)
{
  LOG("INCLCascadeIntranuke", pDEBUG)
    << "default ctor";
}

//______________________________________________________________________________
INCLCascadeIntranuke::INCLCascadeIntranuke(string config) :
  EventRecordVisitorI("genie::INCLCascadeIntranuke", config),
  theINCLConfig(0), theINCLModel(0), theDeExcitation(0)
{
  LOG("INCLCascadeIntranuke", pDEBUG)
    << "ctor from config " << config;
}

//______________________________________________________________________________
INCLCascadeIntranuke::~INCLCascadeIntranuke()
{

  // Config is owned by model once handed over
  if ( theINCLConfig   ) { theINCLConfig=0;   }
  if ( theINCLModel    ) { delete theINCLModel;    theINCLModel=0;    }
  if ( theDeExcitation ) { delete theDeExcitation; theDeExcitation=0; }

}

//______________________________________________________________________________
void INCLCascadeIntranuke::LoadConfig(void)
{
  fResonanceDecayer = nullptr;
  fResonanceDecayer = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("Decayer"));
  assert(fResonanceDecayer);

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  theConfig = incl_nucleus->getConfig();
  if(theConfig->getINCLXXDataFilePath().size() == 0){ // if the inclxx data file path is not setup, then initialize the genie::NucleusGenINCL
    AlgFactory * algf = AlgFactory::Instance();
    auto incl = dynamic_cast<const NucleusGenI*> (algf->GetAlgorithm("genie::NucleusGenINCL", "Default")); // initialize INCL configuration for hadron model, will never use it.
    assert(incl);
  }

}

//______________________________________________________________________________
bool INCLCascadeIntranuke::AddDataPathFlags(size_t& nflags, char** flags) {
  return false;

}

//______________________________________________________________________________

bool INCLCascadeIntranuke::LookForAndAddValidPath(std::vector<std::string>& datapaths,
    size_t defaultIndx,
    const char* optString,
    size_t& nflags, char** flags) {

  return false;
}

//______________________________________________________________________________
// a helper function to add final state from INCL hadron nucleus scattering to GENIE GHEP Event Record.
void INCLCascadeIntranuke::AddINCLParticle(int i, G4INCL::EventInfo &result, GHepRecord * evrec, int first_mother, int second_mother) const {
  double EKin  = result.EKin[i];
  double px  = result.px[i];
  double py  = result.py[i];
  double pz  = result.pz[i];
  short origin = result.origin[i];
  int pdg = result.PDGCode[i];
  double mass = ((px*px + py*py + pz*pz) - EKin * EKin) / (2.0 * EKin);
  double E = EKin + mass;

  TLorentzVector p4(px/1000., py/1000., pz/1000., E/1000.);
  TLorentzVector x4(0,0,0,0);
  if(pdg > 1000){
    if(pdg != 2212 && pdg != 2112 && pdg != 3122){
      int S = pdg / 1000000;
      int pdg_no_s = pdg % 1000000;
      int A = pdg_no_s % 1000;
      int Z = pdg_no_s / 1000;
      pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
      pdg = this->INCLPDG_to_GHEPPDG(pdg, A, Z, S);
    }
  }
  evrec->AddParticle(pdg, kIStStableFinalState, first_mother, second_mother, -1, -1, p4, x4);
}
//______________________________________________________________________________
int INCLCascadeIntranuke::doCascade(GHepRecord * evrec) const {

  // do hadron nucleus cascade
  // re-implement the interface from inclxx/main/src/INCLCascade.cc
  int tpos = evrec->TargetNucleusPosition();
  GHepParticle * target = evrec->Particle(tpos);
  GHepParticle * pprobe = evrec->Probe();
  G4INCL::ParticleType theType = this->PDG_to_INCLType(pprobe->Pdg());
  G4INCL::ParticleSpecies  theSpecies(theType);
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  // setup INCL config
  theConfig = incl_nucleus->getConfig();
  theConfig->setProjectileSpecies(theSpecies);
  theConfig->setProjectileKineticEnergy((pprobe->E() - pprobe->Mass())*1000.);
  theConfig->setTargetA(target->A());
  theConfig->setTargetZ(target->Z());
  theConfig->setTargetS(0);

  // initialize INCL model
  theINCLModel = new G4INCL::INCL(theConfig);
  G4INCL::EventInfo result;
  // process INCL event
  result = theINCLModel->processEvent();
  if( result.transparent ) {
    evrec->AddParticle(pprobe->Pdg(), kIStStableFinalState, 0, -1, -1, -1, *pprobe->P4(),*pprobe->X4());
    evrec->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1, -1, -1, -1, *target->P4(),*target->X4());
  }
  else{

    int n_outgoing = 0;
    // add INCL final states from cascade to GHEP Event Record
    for(int i=0; i<result.nParticles; ++i) {
      INCLCascadeIntranuke::AddINCLParticle(i, result, evrec, 0, 1);
      n_outgoing++;
    }

    // add nuclear remnant from cascade to GHEP Event Record
    // do we always have only one remnant?
    double Rem_px = result.pxRem[0]/1000.;
    double Rem_py = result.pyRem[0]/1000.;
    double Rem_pz = result.pzRem[0]/1000.;
    double Rem_Kin = result.EKinRem[0]/1000.;
    double Rem_EStar = result.EStarRem[0]/1000.;
    int A = result.ARem[0];
    int Z = result.ZRem[0];
    int S = result.SRem[0];
    int pdg = genie::pdg::IonPdgCode( A , Z, S, 0 );
    TParticlePDG * prem = PDGLibrary::Instance()->Find(pdg);
    int PreDeExPDG = pdg;
    double Rem_E = 0;
    if(!prem){
      PreDeExPDG = kPdgHadronicBlob;
      double Rem_p2 = std::sqrt(Rem_px*Rem_px + Rem_py*Rem_py + Rem_pz*Rem_pz);
      double Rem_mass = (Rem_p2*Rem_p2 - Rem_Kin*Rem_Kin) / ( 2.0 * Rem_Kin);
      Rem_E = Rem_Kin + Rem_mass;
    }
    else{
      double Rem_MStar = Rem_EStar + prem->Mass();
      Rem_E  = std::sqrt(Rem_px*Rem_px + Rem_py*Rem_py + Rem_pz*Rem_pz + Rem_MStar*Rem_MStar);
    }
    TLorentzVector p4mom(Rem_px, Rem_py, Rem_pz, Rem_E);
    TLorentzVector p4posi(0,0,0,0);
    evrec->AddParticle(PreDeExPDG, kIStPreDeExNuclearRemnant, 1, -1, -1, -1, p4mom, p4posi);

    // running de-excitation model

    switch(theConfig->getDeExcitationType()){
      case G4INCL::DeExcitationABLAXX:
        {
          std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLAblaInterface>(theConfig);
          theDeExcitation->deExcite(&result);
          break;
        }
      case G4INCL::DeExcitationABLA07:
        {
          std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<ABLA07CXX::Abla07Interface>(theConfig);
          theDeExcitation->deExcite(&result);
          break;
        }
      case G4INCL::DeExcitationGEMINIXX:
        {
          std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLGEMINIXXInterface>(theConfig);
          theDeExcitation->deExcite(&result);
          break;
        }
      case G4INCL::DeExcitationNone:
        {
          // Skip the de-excitation step entirely and end processing early
          return 0;
          break;
        }
      default:
        {
          exit(1);
          break;
        }
    }

    // add final states from de-excitation to GHEP Event Record
    for(int i=n_outgoing; i<result.nParticles; ++i) {
      INCLCascadeIntranuke::AddINCLParticle(i, result, evrec, n_outgoing + 2);
      if(i == (result.nParticles - 1) && ( pdg > 1000000000 && pdg%10000/10 > 4)){
        evrec->Particle(i)->SetStatus(kIStFinalStateNuclearRemnant);
      }
    }
  }
  return 0;
}

void INCLCascadeIntranuke::ProcessEventRecord(GHepRecord * evrec)  const {
  LOG("INCLCascadeIntranuke", pINFO) << "Start with this event";

  fGMode = evrec->EventGenerationMode();
  if(fGMode == kGMdHadronNucleus || fGMode == kGMdPhotonNucleus){
    this->doCascade(evrec);
    return;
  }

  this->PreparePrimaryVertex(evrec);
  //this->DecayResonance(evrec);

  // LOG("INCLCascadeIntranuke", pNOTICE) << "is resonacne : " << evrec->Summary()->ProcInfo().IsResonant();

  // FIXME: start with the NC process
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  incl_target = incl_nucleus->getNuclues();
  theConfig = incl_nucleus->getConfig();
  propagationModel =  incl_nucleus->getPropagationModel();
  incl_target->setParticleNucleusCollision();
  std::unique_ptr<G4INCL::FinalState> finalState(new FinalState);
  // primarylepton = evrec->FinalStatePrimaryLepton();
  prob = evrec->Probe();

  this->preCascade();

  double currentTime = 0.0;
  // double temfin;

  tempFinalState.clear();
  stepFinalState.clear();
  istep = 0;

  double maximumTime = 29.8 * std::pow(incl_target->getA(), 0.16);
  propagationModel->setStoppingTime(maximumTime);
  propagationModel->setNucleus(incl_target);
  propagationModel->generateAllAvatars();

  this->fillFinalState(evrec, finalState.get());

  /*
     int idx = 0;
     double TLab;
     temfin = 29.8 * std::pow(incl_target->getA(), 0.16);
     if(TLab > 2000.)
     temfin *= (5.8E4-TLab)/5.6E4;
     */

  //    double maximumTime = temfin;
  //    propagationModel->setStoppingTime(maximumTime);
  //    propagationModel->setNucleus(incl_target);
  //    propagationModel->generateAllAvatars();

  // stopping time:
  // INCL don't have the stopping time for neutrino.
  // we can calculate the longest stopping time for all the daughters from primary interaction.

  //LOG("INCLCascadeIntranuke", pWARN) << incl_target->print();
  incl_target->applyFinalState(finalState.get());
  //LOG("INCLCascadeIntranuke", pWARN) << incl_target->print();

  // LOG("INCLCascadeIntranuke", pWARN) << incl_target->print();
  //    incl_target->getStore()->getBook().incrementCascading();   // FIXME
  incl_target->getStore()->getBook().incrementAcceptedCollisions();
  int step = 0;


  while(step < 10000000 && continueCascade()){
    IAvatar *avatar = propagationModel->propagate(finalState.get());
    finalState->reset();
    if(avatar == 0) break; // No more avatars in the avatar list.
    G4INCL::ParticleList mother_list = avatar->getParticles();
    backup_mother.clear();
    for(G4INCL::ParticleIter imom =  mother_list.begin(); imom != mother_list.end(); imom++){
      this->fillStep(*imom, backup_mother, -1, -1);
    }
    avatar->fillFinalState(finalState.get());
    // Must fill event record before incl nucleus applyFinalState
    // applyFinalState will delete destroyed particles.
    LOG("INCLCascadeIntranuke", pNOTICE) << finalState->print();
    this->fillEventRecord(finalState.get(), mother_list, evrec, propagationModel->getCurrentTime(), avatar->getType());
    incl_target->applyFinalState(finalState.get());
    LOG("INCLCascadeIntranuke", pNOTICE) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
    delete avatar;
    step++;
  }


  primarylepton = evrec->FinalStatePrimaryLepton();
  // std::cout << incl_target->print() << std::endl;
  // put the nuclear remnant in the event record
  //
  LOG("INCLCascadeIntranuke", pWARN) << "cascade step: " << step;

  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  double debug_A = incl_target->getA();
  double debug_Z = incl_target->getZ();
  this->postCascade(evrec, finalState.get());
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  TObjArrayIter piter(evrec);
  GHepParticle * p = nullptr;
  int stable_finalstate = 0;
  while ( (p = (GHepParticle *) piter.Next() ) ) {
    if(p->Status() == kIStStableFinalState){
      // FIXME: for NC QE, it is only neutron and proton
      stable_finalstate++;
    }
  }

  const int n_outgoing = theEventInfo.nParticles;

  LOG("INCLCascadeIntranuke", pWARN) << "stable_finalstate and n_outgoing: " << stable_finalstate << "  " << n_outgoing;

  // LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Before De-Excitation : " << stable_finalstate - theEventInfo.nParticles;
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants   : " << theEventInfo.nRemnants;
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants A : " << theEventInfo.ARem[0];
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants Z : " << theEventInfo.ZRem[0];
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nRemnants S : " << theEventInfo.SRem[0];
  // LOG("INCLCascadeIntranuke", pNOTICE) << "nParticles   : " << theEventInfo.nParticles;
  // for(int i = 0; i < theEventInfo.nParticles; i++){
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles PDG: " << theEventInfo.PDGCode[i];
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles px : " << theEventInfo.px[i];
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles py : " << theEventInfo.py[i];
  //   LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles pz : " << theEventInfo.pz[i];
  // }


  double Rem_p2 = theEventInfo.pxRem[0]*theEventInfo.pxRem[0]
    + theEventInfo.pyRem[0]*theEventInfo.pyRem[0]
    + theEventInfo.pzRem[0]*theEventInfo.pzRem[0];
  double Rem_Kin = theEventInfo.EKinRem[0];
  double Rem_M = (Rem_p2 - Rem_Kin*Rem_Kin) / ( 2.0 * Rem_Kin);

  // GHepParticle *p4l = evrec->FinalStatePrimaryLepton();

  double Rem_px = theEventInfo.pxRem[0]/1000.;
  double Rem_py = theEventInfo.pyRem[0]/1000.;
  double Rem_pz = theEventInfo.pzRem[0]/1000.;
  double Rem_E =  std::sqrt(Rem_p2 + Rem_M*Rem_M)/1000.;

  TLorentzVector p4mom(Rem_px, Rem_py, Rem_pz, Rem_E);

  TLorentzVector p4posi(0,0,0,0);
  // LOG("INCLCascadeIntranuke", pWARN) << incl_target->print();
  int A = incl_target->getA();
  int Z = incl_target->getZ();
  int S = incl_target->getS(); // INCL and ABLA support hypernuclei
  int pdg = genie::pdg::IonPdgCode( A , Z, S, 0 );

  LOG("INCLCascadeIntranuke", pWARN) << "remnant A Z S: " <<  A << " " << Z << " " << S;
  LOG("INCLCascadeIntranuke", pWARN) << "remnant pdg: " <<  pdg;
  TParticlePDG * prem = PDGLibrary::Instance()->Find(pdg);
  int PreDeExPDG = pdg;
  if(!prem){
    PreDeExPDG = kPdgHadronicBlob;
  }
  int rem_mom_id = (evrec->Particle(3)->FirstDaughter() == -1) ? 3 : evrec->Particle(3)->FirstDaughter();
  evrec->AddParticle(PreDeExPDG, kIStPreDeExNuclearRemnant, rem_mom_id, -1, -1, -1, p4mom, p4posi);

  switch(theConfig->getDeExcitationType()){
    case G4INCL::DeExcitationABLAXX:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLAblaInterface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    case G4INCL::DeExcitationABLA07:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<ABLA07CXX::Abla07Interface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    case G4INCL::DeExcitationGEMINIXX:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLGEMINIXXInterface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    case G4INCL::DeExcitationNone:
      {
        // If de-excitations are disabled, then skip this step entirely
        // and return early
        return;
        break;
      }
    default:
      {
        exit(1);
        break;
      }
  }

  // particle from de-excitation
  const int remnant_id = evrec->GetEntries() - 1;
  if(theEventInfo.nParticles > n_outgoing){
    for(int i = n_outgoing; i < theEventInfo.nParticles; i++){
      LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles PDG: " << theEventInfo.PDGCode[i];
      LOG("INCLCascadeIntranuke", pNOTICE) << "Final State Particles PDG: " << this->INCLPDG_to_GHEPPDG(theEventInfo.PDGCode[i], theEventInfo.A[i], theEventInfo.Z[i], theEventInfo.S[i]);
      int depdg = this->INCLPDG_to_GHEPPDG(theEventInfo.PDGCode[i], theEventInfo.A[i], theEventInfo.Z[i], theEventInfo.S[i]);
      TParticlePDG * p = PDGLibrary::Instance()->Find(depdg);

      LOG("INCLCascadeIntranuke", pNOTICE) << depdg;

      double p2 = theEventInfo.px[i]*theEventInfo.px[i]
        + theEventInfo.py[i]*theEventInfo.py[i]
        + theEventInfo.pz[i]*theEventInfo.pz[i];
      double EKin = theEventInfo.EKin[i];
      double mass = (p2 - EKin*EKin) / (2.0 * EKin);
      double M = p->Mass();
      if(std::fabs(mass/1000. - M) > 0.05){
        LOG("INCLCascadeIntranuke", pERROR) << " particle from de-excitation is unphysical for (" << depdg << "), mass = " << mass/1000. << "(" << M <<")";
      }
      //double E = sqrt(p2/1000000. + M*M);
      double E = (EKin + mass) / 1000.;
      TLorentzVector p4mom(theEventInfo.px[i] / 1000.,
          theEventInfo.py[i] / 1000.,
          theEventInfo.pz[i] / 1000.,
          E);
      TLorentzVector p4posi(0,0,0,0);

      EGHepStatus ptype;
      if(theEventInfo.A[i] > 4){
        ptype = kIStFinalStateNuclearRemnant;
      } else {
        ptype = kIStStableFinalState;
      }
      evrec->AddParticle(depdg, ptype, remnant_id, -1, -1, -1, p4mom, p4posi);
    }
  }
  else{
    TLorentzVector p4posi(0,0,0,0);
    int A = incl_target->getA();
    int Z = incl_target->getZ();
    int S = incl_target->getS();
    int pdg = genie::pdg::IonPdgCode( A , Z, std::abs(S), 0 );
    TParticlePDG * p = PDGLibrary::Instance()->Find(pdg);
    double M = p->Mass();
    double Rem_E = sqrt(Rem_p2/1000000. + M*M);
    TLorentzVector p4mom(Rem_px, Rem_py, Rem_pz, Rem_E);
    LOG("INCLCascadeIntranuke", pNOTICE) << pdg;
    evrec->AddParticle(pdg, kIStFinalStateNuclearRemnant, remnant_id, -1, -1, -1, p4mom, p4posi);
  }


  // print interaction in each step
  for(std::map<int, std::vector<INCLRecord>>::iterator it = stepFinalState.begin(); it != stepFinalState.end(); it++){
    LOG("INCLCascadeIntranuke", pINFO) << "step : " << it->first;
    for(auto ip = it->second.begin(); ip != it->second.end(); ip++){
      std::string fs_particle;
      switch (ip->mother_index) {
        case -2: fs_particle = "Modified  particle "; break;
        case -3: fs_particle = "Outgoing  particle "; break;
        case -4: fs_particle = "Destoried particle "; break;
        case -5: fs_particle = "Created   particle "; break;
        case -6: fs_particle = "Mother    particle "; break;
        default:;
      }
      LOG("INCLCascadeIntranuke", pINFO)  << fs_particle  << "  global id : " << ip->global_index << "  "
        << "pdg id : " << ip->pdgid << "  "
        << "mother id : " << ip->mother_index << "  "
        << "local id : " << ip->local_index;
    }
  }
  bool is_conservation = this->BaryonNumberConservation(evrec);
  if(!is_conservation){
    exit(1);
  }
  //evrec->Print(std::cout);
  //if(debug_Z == 6 && debug_A == 6){
  //  exit(1);
  //}

  LOG("INCLCascadeIntranuke", pINFO) << "Done with this event";
}

bool INCLCascadeIntranuke::CanRescatter(const GHepParticle * p) const {

  return false;
}

void INCLCascadeIntranuke::TransportHadrons(GHepRecord * evrec) const {


}

//______________________________________________________________________________
int INCLCascadeIntranuke::pdgcpiontoA(int pdgc) const {

  if      ( pdgc == 2212 || pdgc == 2112 ) return 1;
  else if ( pdgc ==  211 || pdgc == -211 || pdgc == 111 ) return 0;
  return 0;  // return something

}

//______________________________________________________________________________
int INCLCascadeIntranuke::pdgcpiontoZ(int pdgc) const {

  if      ( pdgc == 2212 || pdgc == 211 ) return 1;
  else if ( pdgc == 2112 || pdgc == 111 ) return 0;
  else if ( pdgc == -211 ) return -1;
  return 0; // return something

}

//______________________________________________________________________________
bool INCLCascadeIntranuke::NeedsRescattering(const GHepParticle * p) const {

  // checks whether the particle should be rescattered
  assert(p);
  // attempt to rescatter anything marked as 'hadron in the nucleus'
  return ( p->Status() == kIStHadronInTheNucleus );

}

//______________________________________________________________________________
void INCLCascadeIntranuke::Configure(const Registry & config) {

  LOG("INCLCascadeIntranuke", pDEBUG)
    << "Configure from Registry: '" << config.Name() << "'\n"
    << config;

  Algorithm::Configure(config);
  this->LoadConfig();

}

//___________________________________________________________________________
void INCLCascadeIntranuke::Configure(string param_set) {

  LOG("INCLCascadeIntranuke", pDEBUG)
    << "Configure from param_set name: " << param_set;

  Algorithm::Configure(param_set);
  this->LoadConfig();

}

bool INCLCascadeIntranuke::continueCascade() const{

  bool continueCascade_ = true;


  if(propagationModel->getCurrentTime() > propagationModel->getStoppingTime()){
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop time : " << propagationModel->getCurrentTime() << " : " << propagationModel->getStoppingTime();

  }
  if(incl_target->getStore()->getBook().getCascading()==0 &&
      incl_target->getStore()->getIncomingParticles().empty()){
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop cascading ";
  }
  int minRemnantSize = 4;
  if(incl_target->getA() <= minRemnantSize) {
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop min size ";
  }

  if(incl_target->getTryCompoundNucleus()) {
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop Compound ";
  }


  return continueCascade_;

}

void INCLCascadeIntranuke::fillEventRecord(G4INCL::FinalState *fs, G4INCL::ParticleList mother_list, GHepRecord * evrec, double time, G4INCL::AvatarType avaType) const {
  std::vector<INCLRecord> stepParticleList;
  stepParticleList.clear();
  istep++;
  LOG("INCLCascadeIntranuke", pDEBUG) << "istep : " << istep;
  for(ParticleIter iter=mother_list.begin(); iter!=mother_list.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -6, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -6, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "mother list ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode();
  }
  ParticleList modified = fs->getModifiedParticles();
  for(ParticleIter iter=modified.begin(); iter!=modified.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -2, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -2, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Modified ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode();
  }
  ParticleList outgoing = fs->getOutgoingParticles();
  for(ParticleIter iter=outgoing.begin(); iter!=outgoing.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -3, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -3, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Outgoing ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode();
  }
  ParticleList destroyed = fs->getDestroyedParticles();
  for(ParticleIter iter=destroyed.begin();  iter!=destroyed.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -4, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -4, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Destroyed ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode();
  }
  ParticleList created = fs->getCreatedParticles();
  for(ParticleIter iter=created.begin(); iter!=created.end(); ++iter){
    this->fillStep(*iter, stepParticleList, -5, time);
    ParticleSpecies ptype = (*iter)->getSpecies();
    stepFinalState[istep].emplace_back((*iter)->getID(), ptype.getPDGCode(), -5, 0);
    LOG("INCLCascadeIntranuke", pDEBUG) << "Created ID : " << (*iter)->getID() << " pdg : " << ptype.getPDGCode();
  }

  int num_partiles = tempFinalState.size();
  int num_temp = stepParticleList.size();
  if(created.size() == 0 && outgoing.size() == 0 && modified.size() == 0) return;
  int index_ = num_partiles;
  if(mother_list.size() == 1){
    // FIXME: channel with mother size = 1 could be reflection, outgoing and decay
    // Try to match the mother in step to the event record history
    // if the mother of this step is not in the event record history and it is a reflection
    // avatar, we will not put this step into event record history
    int mother_position = -1;
    for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
      if(ip->mother_index != -6) continue;
      for(auto it = tempFinalState.begin(); it != tempFinalState.end(); ++it){
        if(it->global_index == ip->global_index){
          mother_position = it->local_index;
        }
      }
    }

    if(mother_position != -1){

      if(avaType == G4INCL::SurfaceAvatarType){
        //
        for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
          if(ip->mother_index == -6 || ip->mother_index == -4) continue;
          LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
          LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
          // put this particle in event record
          tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->p4mom, ip->p4posi);
          evrec->Particle(mother_position)->SetRescatterCode(int(avaType));

          // get the pdg in genie style
          int pdg = ip->pdgid;
          LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
          if(ip->theType == G4INCL::Composite){
            if(pdg != 2212 && pdg != 2112 && pdg != 3122){
              int S = pdg / 1000000;
              int pdg_no_s = pdg % 1000000;
              int A = pdg_no_s % 1000;
              int Z = pdg_no_s / 1000;
              pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
              pdg = this->INCLPDG_to_GHEPPDG(pdg, A, Z, S);
            }
          }
          LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
          // get the type of the particle
          EGHepStatus ptype;

          if(ip->mother_index == -3){
            //ptype = kIStPreDecayResonantState;
            ptype = kIStStableFinalState;
          }
          else{
            ptype = kIStHadronInTheNucleus;
          }

          GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
          evrec->AddParticle(p);
        }
      }
      else if(avaType == G4INCL::DecayAvatarType){

        /* If the decay was Pauli-blocked, make sure the propagation model
         * generates a new decay avatar on the next call to propagate().
         *
         * \bug{Note that we don't generate new decay avatars for deltas that
         * could not satisfy energy conservation. This is in keeping with
         * INCL4.6, but doesn't seem to make much sense to me (DM), as energy
         * conservation can be impossible to satisfy due to weird local-energy
         * conditions, for example, that evolve with time.}
         */

        /* decay channel without daughters, skip it
         * a decay must have more than 1 daughters
         */

        if((created.size() + outgoing.size() + modified.size()) != 1){
          evrec->Particle(mother_position)->SetRescatterCode(int(avaType));
          evrec->Particle(mother_position)->SetStatus(kIStDecayedState);
          for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
            if(ip->mother_index == -6 || ip->mother_index == -4) continue;
            LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
            LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
            // put this particle in event record
            tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->p4mom, ip->p4posi);

            // get the pdg in genie style
            int pdg = ip->pdgid;
            LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
            if(ip->theType == G4INCL::Composite){
              if(pdg != 2212 && pdg != 2112 && pdg != 3122){
                int S = pdg / 1000000;
                int pdg_no_s = pdg % 1000000;
                int A = pdg_no_s % 1000;
                int Z = pdg_no_s / 1000;
                pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
              }
            }
            LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
            // get the type of the particle
            EGHepStatus ptype;

            if(ip->mother_index == -3){
              //ptype = kIStPreDecayResonantState;
              ptype = kIStStableFinalState;
            }
            else{
              ptype = kIStHadronInTheNucleus;
            }

            GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
            evrec->AddParticle(p);
          }
        }
      }
      else if(avaType == G4INCL::UnknownParticle){

      }
    }
    else if(avaType == G4INCL::DecayAvatarType){
      LOG("INCLCascadeIntranuke", pWARN) << "From decayMe: ";
      for(auto imom = backup_mother.begin(); imom != backup_mother.end(); ++imom){
        LOG("INCLCascadeIntranuke", pWARN) << "Mother pdgid: " << imom->pdgid;
        if(!(imom->mother_index == -7 && imom->theType == G4INCL::Composite)){
          LOG("INCLCascadeIntranuke", pFATAL) << "Is not a decay of a composite: " << imom->pdgid;
          exit(1);
        }
        int A = imom->pdgid%1000;
        int Z = (imom->pdgid/1000)%1000;
        int S = (imom->pdgid/1000000)%1000;
        int pdg = genie::pdg::IonPdgCode(A,Z,S,0);
        pdg = this->INCLPDG_to_GHEPPDG(pdg, A, Z, S);
        GHepParticle p(pdg, kIStPreDeExNuclearRemnant, 3, -1, -1, -1, imom->p4mom, imom->p4posi);
        tempFinalState.emplace_back(imom->global_index, imom->pdgid, 3, index_++, imom->p4mom, imom->p4posi);
        evrec->AddParticle(p);
      }
      //evrec->Print(std::cout);
      mother_position = index_ - 1;
      LOG("INCLCascadeIntranuke", pWARN) << "the index of cluster : " << mother_position;
      for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
        if(ip->mother_index == -6 || ip->mother_index == -4) continue;
        LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
        LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
        LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
        LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
        // put this particle in event record
        tempFinalState.emplace_back(ip->global_index, ip->pdgid, mother_position, index_++, ip->p4mom, ip->p4posi);
        evrec->Particle(mother_position)->SetRescatterCode(int(avaType));

        // get the pdg in genie style
        int pdg = ip->pdgid;
        LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
        if(ip->theType == G4INCL::Composite){
          if(pdg != 2212 && pdg != 2112 && pdg != 3122){
            int S = pdg / 1000000;
            int pdg_no_s = pdg % 1000000;
            int A = pdg_no_s % 1000;
            int Z = pdg_no_s / 1000;
            pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
            pdg = this->INCLPDG_to_GHEPPDG(pdg, A, Z, S);
          }
        }
        LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
        // get the type of the particle
        EGHepStatus ptype;

        if(ip->mother_index == -3){
          //ptype = kIStPreDecayResonantState;
          ptype = kIStStableFinalState;
        }
        else{
          ptype = kIStHadronInTheNucleus;
        }

        GHepParticle p(pdg, ptype, mother_position, -1, -1, -1, ip->p4mom, ip->p4posi);
        evrec->AddParticle(p);
      }
      //evrec->Print(std::cout);
      //exit(1);
    }
  }
  else if(mother_list.size() == 2){
    int first_mother_position = -1;
    int second_mother_position = -1;
    // find parents for nn-collision
    for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
      if(ip->mother_index != -6) continue;
      int mother_position = -1;
      bool is_spectator = true;
      for(auto it = tempFinalState.begin(); it != tempFinalState.end(); ++it){
        if(it->global_index == ip->global_index){
          mother_position = it->local_index;
          is_spectator = false;
        }
      }
      if(is_spectator){
        for(auto imom = backup_mother.begin(); imom != backup_mother.end(); ++imom){
          if(imom->global_index == ip->global_index){
            tempFinalState.emplace_back(ip->global_index, imom->pdgid, -1, index_++, imom->p4mom, imom->p4posi);
            GHepParticle p(imom->pdgid, kIStSpectator, -1, -1, -1, -1, imom->p4mom, imom->p4posi);
            evrec->AddParticle(p);
          }
        }
        mother_position = index_ - 1;
      }
      // assign the mother index
      if(first_mother_position == -1)
        first_mother_position = mother_position;
      else
        second_mother_position = mother_position;
    }

    if(first_mother_position > second_mother_position){
      int temp_position = first_mother_position;
      first_mother_position = second_mother_position;
      second_mother_position = temp_position;
    }
    LOG("INCLCascadeIntranuke", pNOTICE) << first_mother_position;
    LOG("INCLCascadeIntranuke", pNOTICE) << second_mother_position;
    evrec->Particle(first_mother_position)->SetRescatterCode(int(avaType));
    evrec->Particle(second_mother_position)->SetRescatterCode(int(avaType));

    for(auto ip = stepParticleList.begin(); ip != stepParticleList.end(); ++ip){
      if(ip->mother_index == -6 || ip->mother_index == -4) continue;
      LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
      LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
      // put this particle in event record
      tempFinalState.emplace_back(ip->global_index, ip->pdgid, first_mother_position, index_++, ip->p4mom, ip->p4posi);

      // get the pdg in genie style
      int pdg = ip->pdgid;
      LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << ip->pdgid;
      if(ip->theType == G4INCL::Composite){
        if(pdg != 2212 && pdg != 2112 && pdg != 3122){
          int S = pdg / 1000000;
          int pdg_no_s = pdg % 1000000;
          int A = pdg_no_s % 1000;
          int Z = pdg_no_s / 1000;
          pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
        }
      }
      LOG("INCLCascadeIntranuke", pINFO) << "PDG : " << pdg;
      // get the type of the particle
      EGHepStatus ptype;
      ptype = kIStHadronInTheNucleus;
      GHepParticle p(pdg, ptype, first_mother_position, second_mother_position, -1, -1, ip->p4mom, ip->p4posi);
      evrec->AddParticle(p);
    }
    // update the daughter's indices for spectator
    if(evrec->Particle(first_mother_position)->FirstDaughter() == -1){
      evrec->Particle(first_mother_position)->SetFirstDaughter(evrec->Particle(second_mother_position)->FirstDaughter());
      evrec->Particle(first_mother_position)->SetLastDaughter(evrec->Particle(second_mother_position)->LastDaughter());
    }
    else if(evrec->Particle(second_mother_position)->FirstDaughter() == -1){
      evrec->Particle(second_mother_position)->SetFirstDaughter(evrec->Particle(first_mother_position)->FirstDaughter());
      evrec->Particle(second_mother_position)->SetLastDaughter(evrec->Particle( first_mother_position)->LastDaughter());
    }

  }
  else{
    LOG("INCLCascadeIntranuke", pNOTICE) << "wrong mother size : " << mother_list.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "created: " << created.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "outgoing: " << outgoing.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "modified: " << modified.size();
    LOG("INCLCascadeIntranuke", pNOTICE) << "destroyed: " << destroyed.size();
    exit(1);
  }


}

void INCLCascadeIntranuke::fillStep(G4INCL::Particle *par, std::vector<INCLRecord> &stepList, int type, double time) const {
  ParticleSpecies ptype = par->getSpecies();
  TLorentzVector p4mom;
  p4mom.SetPx(par->getMomentum().getX() / 1000.);
  p4mom.SetPy(par->getMomentum().getY() / 1000.);
  p4mom.SetPz(par->getMomentum().getZ() / 1000.);
  p4mom.SetE(par->getEnergy() / 1000.);
  TLorentzVector p4posi;
  p4posi.SetX(par->getPosition().getX());
  p4posi.SetY(par->getPosition().getY());
  p4posi.SetZ(par->getPosition().getZ());
  p4posi.SetT(time);
  stepList.emplace_back(par->getID(), ptype.getPDGCode(), type, 0, p4mom, p4posi, ptype.theType);
}

void INCLCascadeIntranuke::postCascade(GHepRecord * evrec, G4INCL::FinalState * finalState) const {
  // Fill in the event information
  theEventInfo.stoppingTime = propagationModel->getCurrentTime();

  // The event bias
  theEventInfo.eventBias = (Double_t) Particle::getTotalBias();

  // Check if the nucleus contains strange particles
  theEventInfo.sigmasInside = incl_target->containsSigma();
  theEventInfo.antikaonsInside = incl_target->containsAntiKaon();
  theEventInfo.lambdasInside = incl_target->containsLambda();
  theEventInfo.kaonsInside = incl_target->containsKaon();

  LOG("INCLCascadeIntranuke", pNOTICE) << "n Sigma and anti-Kaon " << incl_target->containsSigma() << "  " << incl_target->containsAntiKaon() << "  " << incl_target->containsLambda() << "  " << incl_target->containsKaon();

  // FIXME: It seems the INCL will transfer Sigma and anti-Kaon
  // into Lambda via interacting with the first nucleon in Store
  // list
  // - However, de-excitation
  // Capture antiKaons and Sigmas and produce Lambda instead
  theEventInfo.absorbedStrangeParticle = this->decayInsideStrangeParticles(evrec, finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  // Emit strange particles still inside the nucleus
  this->emitInsideStrangeParticles(evrec, finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  theEventInfo.emitKaon = this->emitInsideKaon(evrec, finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  theEventInfo.emitLambda = this->emitInsideLambda(evrec, finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  // Check if the nucleus contains deltas
  theEventInfo.deltasInside = incl_target->containsDeltas();

  // Take care of any remaining deltas

  if(theEventInfo.deltasInside){
    G4INCL::ParticleList const &inside = incl_target->getStore()->getOutgoingParticles();
  }
  //theEventInfo.forcedDeltasOutside = incl_target->decayOutgoingDeltas();   //FIXME: leave resonances to pythia
  theEventInfo.forcedDeltasInside = this->decayInsideDeltas(evrec, finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  // Take care of any remaining etas, omegas, neutral Sigmas and/or neutral kaons
  double timeThreshold=theConfig->getDecayTimeThreshold();
  theEventInfo.forcedPionResonancesOutside = this->decayOutgoingPionResonances(timeThreshold, evrec, finalState);  //FIXME: leave resonances to pythia
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  this->decayOutgoingSigmaZero(timeThreshold, evrec, finalState); //FIXME: leave resonances to pythia
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  this->decayOutgoingNeutralKaon(evrec, finalState); //FIXME: leave resonances to pythia
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  //this->emitInsidePions(evrec, finalState); // FIXME: is it valid?

  // Apply Coulomb distortion, if appropriate
  // Note that this will apply Coulomb distortion also on pions emitted by
  // unphysical remnants (see decayInsideDeltas). This is at variance with
  // what INCL4.6 does, but these events are (should be!) so rare that
  // whatever we do doesn't (shouldn't!) make any noticeable difference.
  // CoulombDistortion::distortOut(incl_target->getStore()->getOutgoingParticles(), incl_target);

  // If the normal cascade predicted complete fusion, use the tabulated
  // masses to compute the excitation energy, the recoil, etc.
  //if(incl_target->getStore()->getOutgoingParticles().size()==0
  //    && (!incl_target->getProjectileRemnant()
  //	|| incl_target->getProjectileRemnant()->getParticles().size()==0)) {
  if( false && (!incl_target->getProjectileRemnant()
        || incl_target->getProjectileRemnant()->getParticles().size()==0)) {

    LOG("INCLCascadeIntranuke", pINFO) << "Cascade resulted in complete fusion, using realistic fusion kinematics";

    incl_target->useFusionKinematics();

    if(incl_target->getExcitationEnergy()<0.) {
      // Complete fusion is energetically impossible, return a transparent
      LOG("INCLCascadeIntranuke", pINFO) << "Complete-fusion kinematics yields negative excitation energy, returning a transparent!";
      theEventInfo.transparent = true;
      return;
    }

  }
  else {

    // Set the excitation energy
    incl_target->setExcitationEnergy(incl_target->computeExcitationEnergy());

    // Make a projectile pre-fragment out of the geometrical and dynamical
    // spectators
    //theEventInfo.nUnmergedSpectators = makeProjectileRemnant();

    // Compute recoil momentum, energy and spin of the nucleus
    int minRemnantSize = 4;
    if(incl_target->getA()==1 && minRemnantSize>1) {
      LOG("INCLCascadeIntranuke", pINFO) << "Computing one-nucleon recoil kinematics. We should never be here nowadays, cascade should stop earlier than this.";
    }
    incl_target->computeRecoilKinematics();

    // Make room for the remnant recoil by rescaling the energies of the
    // outgoing particles. FIXME
    if(incl_target->hasRemnant()) this->rescaleOutgoingForRecoil();
  }

  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  theEventInfo.clusterDecay = this->decayOutgoingClusters(evrec, finalState) || this->decayMe(evrec, finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  incl_target->fillEventInfo(&theEventInfo);

}

bool INCLCascadeIntranuke::preCascade() const {
  // Reset theEventInfo
  theEventInfo.reset();
  EventInfo::eventNumber++;

  // Fill in the event information
  // neutrino projectile
  // theEventInfo.projectileType = projectileSpecies.theType;
  theEventInfo.Ap = 0;
  theEventInfo.Zp = 0;
  theEventInfo.Sp = 0;
  theEventInfo.Ep = prob->P4()->E() * 1000.;

  theEventInfo.At = incl_target->getA();
  theEventInfo.Zt = incl_target->getZ();
  theEventInfo.St = incl_target->getS();

  theEventInfo.transparent = false;
  theEventInfo.impactParameter = 0.;

  theEventInfo.effectiveImpactParameter = 0.;

  return true;

}

void INCLCascadeIntranuke::rescaleOutgoingForRecoil() const {
  TLorentzVector *p4prob = prob->P4();
  TLorentzVector *p4lep  = primarylepton->P4();
  G4INCL::ThreeVector transferQ((p4prob->Px() - p4lep->Px()) * 1000.,
      (p4prob->Py() - p4lep->Py()) * 1000.,
      (p4prob->Pz() - p4lep->Pz()) * 1000.);

  RecoilCMFunctor theRecoilFunctor(incl_target, theEventInfo, transferQ);

  // Apply the root-finding algorithm
  const G4INCL::RootFinder::Solution theSolution = G4INCL::RootFinder::solve(&theRecoilFunctor, 1.0);
  if(theSolution.success) {
    theRecoilFunctor(theSolution.x); // Apply the solution
  } else {
    LOG("INCLCascadeIntranuke", pINFO) << "Couldn't accommodate remnant recoil while satisfying energy conservation, root-finding algorithm failed.";
  }
}

int INCLCascadeIntranuke::INCLPDG_to_GHEPPDG(int pdg, int A, int Z, int S) const{
  //  TParticlePDG * p = PDGLibrary::Instance()->Find(pdg);
  TDatabasePDG * fDatabasePDG = TDatabasePDG::Instance();
  TParticlePDG * p = fDatabasePDG->GetParticle(pdg);
  if(!p){
    if(A != 0){
      int ion_pdg =  genie::pdg::IonPdgCode( A , Z, std::abs(S), 0 );
      TParticlePDG *ion = fDatabasePDG->GetParticle(ion_pdg);
      if(S != 0 && !ion){
        PDGLibrary *pdg_library = PDGLibrary::Instance();
        pdg_library->AddHypernucleus(ion_pdg);
      }
      return ion_pdg;
    }
    else{
      LOG("INCLCascadeIntranuke", pERROR) << "Particle is not identified: pdg = " << pdg << "; A, Z, S => "
        << A << " " << Z << " " << S;
      exit(1);
    }
  }
  else{
    return pdg;
  }
}

void INCLCascadeIntranuke::fillFinalState(GHepRecord * evrec, G4INCL::FinalState * finalState) const{

  // neutrino interaction info
  const ProcessInfo & proc_info = evrec->Summary()->ProcInfo();

  // FIXME: start with the NC process
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();

  double TLab;
  temfin = 29.8 * std::pow(incl_target->getA(), 0.16);

  TObjArrayIter piter(evrec);
  GHepParticle * p = nullptr;
  std::vector<G4INCL::GENIEParticleRecord> eventRecord;
  eventRecord.clear();
  while ( (p = (GHepParticle *) piter.Next() ) ) {
    // the code of the particles in primary neutrino interaction
    G4INCL::GENIERecordCode recordCode;
    if(eventRecord.size() == evrec->ProbePosition())
      recordCode = G4INCL::kProbe;
    else if(eventRecord.size() == evrec->HitNucleonPosition())
      recordCode = G4INCL::kHitNucleon;
    else if(eventRecord.size() == evrec->FinalStatePrimaryLeptonPosition())
      recordCode = G4INCL::kFinalStateLepton;
    else if(eventRecord.size() == evrec->TargetNucleusPosition())
      recordCode = G4INCL::kTarget;
    else if(eventRecord.size() == evrec->RemnantNucleusPosition())
      recordCode = G4INCL::kRemnant;
    else
      recordCode = G4INCL::kUnknown;

    eventRecord.emplace_back(p, int(proc_info.ScatteringTypeId()), recordCode);
  }

  NuclearModel_t nucl_model = incl_nucleus->getHybridModel();
  LOG("INCLCascadeIntranuke", pWARN) << "Nuclear model " << NuclearModel::AsString(nucl_model);

  bool isHybridModel;
  if(nucl_model == kNucmINCL){
    isHybridModel = false;
  }
  else if(nucl_model == kNucmUndefined){
    exit(1);
  }
  else{
    isHybridModel = true;
  }

  LOG("INCLCascadeIntranuke", pWARN) << "Nuclear model " << NuclearModel::AsString(nucl_model) << "  :  " << isHybridModel;

  std::shared_ptr<G4INCL::IAvatar> avatar;
  if(proc_info.IsMEC()){
    avatar = std::make_shared<G4INCL::GENIEAvatar>(0, (incl_nucleus->getHitNNCluster()).get(), incl_nucleus->getNuclues(), &eventRecord, isHybridModel);
    avatar->fillFinalState(finalState);
  }
  else{
    avatar = std::make_shared<G4INCL::GENIEAvatar>(0, incl_nucleus->getHitParticle(), incl_nucleus->getNuclues(), &eventRecord, isHybridModel);
    avatar->fillFinalState(finalState);
  }

  if(finalState->getValidity() != G4INCL::ValidFS){
    LOG("INCLCascadeIntranuke", pWARN)
      << "Enforcing energy conservation: failed! ";
    evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Couldn't select kinematics");
    exception.SwitchOnFastForward();
    throw exception;
  }

  // update the event record after INCL postInteraction
  //
  //

  TObjArrayIter piter1(evrec);
  GHepParticle * p1 = nullptr;
  auto er = eventRecord.begin();
  int idx =0;
  while ( (p1 = (GHepParticle *) piter1.Next() ) ) {
    TLorentzVector *p4 = p1->P4();
    p4->SetPx(er->P3().getX()/1000.);
    p4->SetPy(er->P3().getY()/1000.);
    p4->SetPz(er->P3().getZ()/1000.);
    p4->SetE(std::sqrt(er->P3().mag2() + er->Mass()*er->Mass())/1000.);
    tempFinalState.emplace_back(er->ID(), er->Pdg(), er->FirstMother(), idx++);
    LOG("INCLCascadeIntranuke", pWARN) << er->ID();
    er++;
  }
  evrec->Print(std::cout);

  //LOG("INCLCascadeIntranuke", pWARN) << finalState->print();

  // put the out-going particle into event record

  ParticleList outgoing = finalState->getOutgoingParticles();

  for(ParticleIter iter=outgoing.begin(); iter!=outgoing.end(); ++iter){
    int outp_mother_idx = -1;
    int tmp_idx_ = -1;
    int pdg = 0;
    GHepParticle * p1 = nullptr;
    for(auto er = eventRecord.begin(); er != eventRecord.end(); er++){
      tmp_idx_++;
      if((*iter)->getID() == er->ID()){
        outp_mother_idx = tmp_idx_;
        pdg = er->Pdg();
      }
    }

    GHepParticle p(pdg, kIStStableFinalState, outp_mother_idx, -1, -1, -1,
        TLorentzVector((*iter)->getMomentum().getX() / 1000,
          (*iter)->getMomentum().getY() / 1000,
          (*iter)->getMomentum().getZ() / 1000,
          (*iter)->getEnergy() / 1000),
        TLorentzVector((*iter)->getPosition().getX(),
          (*iter)->getPosition().getY(),
          (*iter)->getPosition().getZ(),
          0)
        );
    evrec->AddParticle(p);
    tempFinalState.emplace_back((*iter)->getID(), pdg, outp_mother_idx, idx++);
  }
  //evrec->Print(std::cout);
  //LOG("INCLCascadeIntranuke", pWARN) << finalState->print();

  if(!outgoing.empty()){

    std::unique_ptr<G4INCL::FinalState> tmpfs(new FinalState);
    G4INCL::ParticleList modified  = finalState->getModifiedParticles();
    G4INCL::ParticleList created   = finalState->getCreatedParticles();
    G4INCL::ParticleList Destroyed = finalState->getDestroyedParticles();

    for(ParticleIter i = modified.begin(), e = modified.end(); i!=e; ++i){
      tmpfs->addModifiedParticle((*i));
    }
    for(ParticleIter i = Destroyed.begin(), e = Destroyed.end(); i!=e; ++i){
      tmpfs->addDestroyedParticle((*i));
    }
    for(ParticleIter i = created.begin(), e = created.end(); i!=e; ++i){
      if(!((*i)->isOutOfWell())){
        tmpfs->addCreatedParticle((*i));
      }
    }
    for(ParticleIter i = outgoing.begin(), e = outgoing.end(); i!=e; ++i){
      tmpfs->addOutgoingParticle((*i));
    }

    finalState->reset();
    modified  = tmpfs->getModifiedParticles();
    created   = tmpfs->getCreatedParticles();
    Destroyed = tmpfs->getDestroyedParticles();
    outgoing =  tmpfs->getOutgoingParticles();

    for(ParticleIter i = modified.begin(), e = modified.end(); i!=e; ++i){
      finalState->addModifiedParticle((*i));
    }
    for(ParticleIter i = Destroyed.begin(), e = Destroyed.end(); i!=e; ++i){
      finalState->addDestroyedParticle((*i));
    }
    for(ParticleIter i = created.begin(), e = created.end(); i!=e; ++i){
      finalState->addCreatedParticle((*i));
    }
    for(ParticleIter i = outgoing.begin(), e = outgoing.end(); i!=e; ++i){
      finalState->addOutgoingParticle((*i));
    }
  }
  //LOG("INCLCascadeIntranuke", pWARN) << finalState->print();

  return;

}

G4INCL::ParticleType INCLCascadeIntranuke::PDG_to_INCLType(int pdg) const {
  switch(pdg){
    case 2212: return G4INCL::Proton;
    case 2112: return G4INCL::Neutron;
    case 211: return G4INCL::PiPlus;
    case -211: return G4INCL::PiMinus;
    case 111: return G4INCL::PiZero;
    case 2224: return G4INCL::DeltaPlusPlus;
    case 2214: return G4INCL::DeltaPlus;
    case 2114: return G4INCL::DeltaZero;
    case 1114: return G4INCL::DeltaMinus;
    case 221: return G4INCL::Eta;
    case 223: return G4INCL::Omega;
    case 331: return G4INCL::EtaPrime;
    case 22: return G4INCL::Photon;
    case 3122: return G4INCL::Lambda;
    case 3222: return G4INCL::SigmaPlus;
    case 3212: return G4INCL::SigmaZero;
    case 3112: return G4INCL::SigmaMinus;
    case -2212: return G4INCL::antiProton;
    case 3312: return G4INCL::XiMinus;
    case 3322: return G4INCL::XiZero;
    case -2112: return G4INCL::antiNeutron;
    case -3122: return G4INCL::antiLambda;
    case -3222: return G4INCL::antiSigmaPlus;
    case -3212: return G4INCL::antiSigmaZero;
    case -3112: return G4INCL::antiSigmaMinus;
    case -3312: return G4INCL::antiXiMinus;
    case -3322: return G4INCL::antiXiZero;
    case 321: return G4INCL::KPlus;
    case 311: return G4INCL::KZero;
    case -311: return G4INCL::KZeroBar;
    case -321: return G4INCL::KMinus;
    case 310: return G4INCL::KShort;
    case 130: return G4INCL::KLong;
    default:
              return G4INCL::UnknownParticle;
  }

}

void INCLCascadeIntranuke::PreparePrimaryVertex(GHepRecord * event_rec) const{
  // 1. for CCQE
  int inucl = -1;
  inucl = event_rec->RemnantNucleusPosition();
  GHepParticle * nucl = event_rec->Particle(inucl);
  nucl->SetStatus(kIStIntermediateState);

  //exit(1);
}

void INCLCascadeIntranuke::DecayResonance(GHepRecord *evrec) const{
  if(evrec->Summary()->ProcInfo().IsResonant()){
    bool decay_flag = true;
    while(decay_flag){
      TObjArrayIter piter(evrec);
      GHepParticle * p = nullptr;
      decay_flag = false;
      while ( (p = (GHepParticle *) piter.Next() ) ) {
        if(p->Status() == kIStPreDecayResonantState && p->FirstDaughter() == -1 && PDG_to_INCLType(p->Pdg()) == G4INCL::UnknownParticle){
          decay_flag = true;
        }
      }
      if(decay_flag){
        fResonanceDecayer->ProcessEventRecord(evrec);
      }
    }
  }
}

bool INCLCascadeIntranuke::BaryonNumberConservation(GHepRecord *evrec) const {
  int final_C = 0, final_A = 0;
  Target *tgt = evrec->Summary()->InitStatePtr()->TgtPtr();
  int inital_C = evrec->Probe()->Charge() + tgt->Charge();
  int target_A = tgt->A();
  LOG("INCLCascadeIntranuke", pINFO) << "Inital states charge and A: " << inital_C << " " << target_A;
  TObjArrayIter piter(evrec);
  GHepParticle * p = nullptr;
  while ( (p = (GHepParticle *) piter.Next() ) ) {
    // the code of the particles in primary neutrino interaction
    std::string particle_class = std::string(PDGLibrary::Instance()->Find(p->Pdg())->ParticleClass());
    if(p->Status() == kIStStableFinalState || p->Status() == kIStFinalStateNuclearRemnant){
      if(genie::pdg::IsNeutronOrProton(p->Pdg())){
        final_A++;
      }
      else if(genie::pdg::IsBaryonResonance(p->Pdg())){
        final_A++;
      }
      else if(particle_class.find("Baryon") != std::string::npos){
        if(p->Pdg() < 0){
          final_A--;
        }
        else{
          final_A++;
        }
      }
      else if(genie::pdg::IsIon(p->Pdg())){
        final_A+=p->A();
      }
      final_C += p->Charge()/3;
      LOG("INCLCascadeIntranuke", pINFO) << "PDG and Charge: " << p->Pdg() << "  " << p->Charge() << "  " << p->A() << "  " << final_A << "  " << final_C;
    }
  }
  LOG("INCLCascadeIntranuke", pINFO) << "Final states charge and A: " << final_C << " " << final_A;
  if(inital_C != final_C || target_A != final_A){
    evrec->Print(std::cout);
    LOG("INCLCascadeIntranuke", pFATAL) << "Violating charge number and baryon number!";
    LOG("INCLCascadeIntranuke", pFATAL) << "Inital states charge and A: " << inital_C << " " << target_A;
    LOG("INCLCascadeIntranuke", pFATAL) << "Final states charge and A: " << final_C << " " << final_A;
    exit(1);
    return false;
  }
  return true;

}

#include "INCLPostCascade.icc"

#endif // __GENIE_INCL_ENABLED__
