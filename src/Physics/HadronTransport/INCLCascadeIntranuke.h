#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _INCLCascadeIntranuke_H_
#define _INCLCascadeIntranuke_H_

#include <string>
using std::string;

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Conventions/GMode.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include <TLorentzVector.h>

namespace G4INCL {
  class Config;
  class INCL;
  class IDeExcitation;
}

namespace genie {

  class GHepParticle;

  class INCLCascadeIntranuke: public EventRecordVisitorI {

    public :
      INCLCascadeIntranuke();
      INCLCascadeIntranuke(std::string config);
      ~INCLCascadeIntranuke();

      int pdgcpiontoA(int pdgc) const;
      int pdgcpiontoZ(int pdgc) const;

      // implement the EventRecordVisitorI interface
      // also the LoadConfig interface

      void Configure (const Registry & config);
      void Configure (string param_set);

      virtual void ProcessEventRecord(GHepRecord * event_rec) const;

    protected:
      virtual void LoadConfig (void);

      bool CanRescatter(const GHepParticle * p) const;
      bool IsInNucleus(const GHepParticle * p) const;
      void TransportHadrons(GHepRecord * evrec) const;
      int  doCascade(GHepRecord * event_rec) const;
      void AddINCLParticle(int i, G4INCL::EventInfo &result, GHepRecord * event_rec, int first_mother = -1, int second_mother = -1) const;
      bool NeedsRescattering(const GHepParticle * p) const;

      bool AddDataPathFlags(size_t& nflags, char** flags);
      bool LookForAndAddValidPath(std::vector<std::string>& datapaths,
          size_t defaultIndx,
          const char* optString,
          size_t& nflags, char** flags);


      mutable int            fRemnA;         ///< remnant nucleus A
      mutable int            fRemnZ;         ///< remnant nucleus Z
      mutable TLorentzVector fRemnP4;        ///< P4 of remnant system
      mutable GEvGenMode_t   fGMode;
      mutable G4INCL::Config        *theINCLConfig;
      mutable G4INCL::INCL          *theINCLModel;
      mutable G4INCL::IDeExcitation *theDeExcitation;
    private:

      // GENIE method
      mutable GHepParticle *prob;
      mutable GHepParticle *primarylepton;
      G4INCL::ParticleType PDG_to_INCLType(int pdg) const;
      const EventRecordVisitorI * fResonanceDecayer;

      // INCL can handle delta resonances, but it can't
      // handle higher order resonances, 
      // Will decay the kIStPreDecayResonantState except delta
      // resonances
      void PreparePrimaryVertex(GHepRecord * event_rec) const;
      void DecayResonance(GHepRecord * event_rec) const;

      bool BaryonNumberConservation(GHepRecord * event_rec) const;

      // INCL method
      void postCascade(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;
      //void postCascadeEventRecord(GHepRecord * event_rec, G4INCL::FinalState * finalState, int pre, int post) const;
      bool preCascade() const;

      bool continueCascade() const;
      mutable G4INCL::Nucleus  *incl_target;
      mutable G4INCL::Particle *hit_nucleon;
      mutable G4INCL::Particle *hit_cluster;
      mutable G4INCL::Config * theConfig;
      mutable G4INCL::StandardPropagationModel * propagationModel;
      mutable G4INCL::EventInfo theEventInfo;

      mutable double temfin;
      void fillFinalState(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;
      //void fillFinalStateNCEL(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;

      // FIXME: put the G4INCL::Nucleus::<post cascade func> in here 
      // to get the event record
      bool decayInsideStrangeParticles(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;
      // Emit strange particles still inside the nucleus
      //  \brief Force emission of all strange particles inside the nucleus.
      void emitInsideStrangeParticles(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;
      /// \brief Force emission of all Lambda (desexitation code with strangeness not implanted yet)
      int  emitInsideLambda(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;
      /// \brief Force emission of all Kaon inside the nucleus
      bool emitInsideKaon(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;

      /** \brief Force the decay of deltas inside the nucleus.
       * 
       * \return true if any delta was forced to decay.
       */

      bool decayInsideDeltas(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;
      /// \brief Force emission of all pions inside the nucleus.
      void emitInsidePions(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;

      /** \brief Force the decay of unstable outgoing clusters.
       *
       * \return true if any cluster was forced to decay.
       */
      bool decayOutgoingClusters(GHepRecord * event_rec, G4INCL::FinalState * finalState) const;

      /** \brief Force the transformation of outgoing Neutral Kaon into propation eigenstate.
       * \return true if any kaon was forced to decay.
       */
      bool decayOutgoingNeutralKaon(GHepRecord * evrec, G4INCL::FinalState * finalState) const;

      /** \brief Force the decay of outgoing Neutral Sigma.
       * \return true if any Sigma was forced to decay.
       */
      bool decayOutgoingSigmaZero(double timeThreshold, GHepRecord * evrec, G4INCL::FinalState * finalState) const;

      /** \brief Force the decay of outgoing PionResonances (eta/omega).
       * \return true if any eta was forced to decay.
       */

      bool decayOutgoingPionResonances(double timeThreshold, GHepRecord * evrec, G4INCL::FinalState * finalState) const;

      /** \brief Force the phase-space decay of the Nucleus.
       *
       * Only applied if Z==0 or N==0.
       *
       * \return true if the nucleus was forced to decay.
       */
      bool decayMe(GHepRecord * event_rec, G4INCL::FinalState * finalState) const ;


      struct INCLRecord{
        int global_index;        // Each particles in INCLXX will have a unique ID, it is a global index for every simulation run.
        int pdgid;	       // PDG ID of particles in INCLXX
        int mother_index;        // mother index of particles in each event
        int local_index;         // local index of particles in each event
        TLorentzVector p4mom;
        TLorentzVector p4posi;
        G4INCL::ParticleType theType;   // INCL Particle type

        INCLRecord(int g_id, int p_id, int m_id, int l_id):
          global_index(g_id),
          pdgid(p_id),
          mother_index(m_id),
          local_index(l_id){}
        INCLRecord(int g_id, int p_id, int m_id, int l_id, TLorentzVector mom, TLorentzVector posi):
          global_index(g_id),
          pdgid(p_id),
          mother_index(m_id),
          local_index(l_id),
          p4mom(mom),
          p4posi(posi){}
        INCLRecord(int g_id, int p_id, int m_id, int l_id, TLorentzVector mom, TLorentzVector posi, G4INCL::ParticleType pType):
          global_index(g_id),
          pdgid(p_id),
          mother_index(m_id),
          local_index(l_id),
          p4mom(mom),
          p4posi(posi), 
          theType(pType){}
      };



      mutable std::vector<INCLRecord> tempFinalState;
      mutable std::map<int, std::vector<INCLRecord>> stepFinalState;
      mutable std::vector<INCLRecord> backup_mother; // this cantainer is used to store the initial momentum and position of mother particles in binary collision
      mutable int istep;

      void fillEventRecord(G4INCL::FinalState *fs, G4INCL::ParticleList mother_list, GHepRecord * evrec, double time, G4INCL::AvatarType avaType) const;
      void fillStep(G4INCL::Particle *par, std::vector<INCLRecord> &stepList, int type, double time) const;

      int INCLPDG_to_GHEPPDG(int pdg, int A, int Z, int S) const;

      /// \brief Class to adjust remnant recoil
      class RecoilFunctor : public G4INCL::RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleEnergies
           *
           * The constructor sets the private class members.
           */
          RecoilFunctor(G4INCL::Nucleus * const n, const G4INCL::EventInfo &ei) :
            RootFunctor(0., 1E6),
            nucleus(n),
            outgoingParticles(n->getStore()->getOutgoingParticles()),
            theEventInfo(ei) {
              for(G4INCL::ParticleIter p=outgoingParticles.begin(), e=outgoingParticles.end(); p!=e; ++p) {
                particleMomenta.push_back((*p)->getMomentum());
                particleKineticEnergies.push_back((*p)->getKineticEnergy());
              }
              G4INCL::ProjectileRemnant * const aPR = n->getProjectileRemnant();
              if(aPR && aPR->getA()>0) {
                particleMomenta.push_back(aPR->getMomentum());
                particleKineticEnergies.push_back(aPR->getKineticEnergy());
                outgoingParticles.push_back(aPR);
              }
            }
          virtual ~RecoilFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle energies
           * \return the energy-conservation violation
           */
          double operator()(const double x) const {
            scaleParticleEnergies(x);
            return nucleus->getConservationBalance(theEventInfo,true).energy;
          }

          /// \brief Clean up after root finding
          void cleanUp(const bool success) const {
            if(!success)
              scaleParticleEnergies(1.);
          }

        private:
          /// \brief Pointer to the nucleus
          G4INCL::Nucleus *nucleus;
          /// \brief List of final-state particles.
          G4INCL::ParticleList outgoingParticles;
          // \brief Reference to the EventInfo object
          G4INCL::EventInfo const &theEventInfo;
          /// \brief Initial momenta of the outgoing particles
          std::list<G4INCL::ThreeVector> particleMomenta;
          /// \brief Initial kinetic energies of the outgoing particles
          std::list<double> particleKineticEnergies;

          /** \brief Scale the kinetic energies of the outgoing particles.
           *
           * \param rescale scale factor
           */
          void scaleParticleEnergies(const double rescale) const {
            // Rescale the energies (and the momenta) of the outgoing particles.
            G4INCL::ThreeVector pBalance = nucleus->getIncomingMomentum();
            std::list<G4INCL::ThreeVector>::const_iterator iP = particleMomenta.begin();
            std::list<double>::const_iterator iE = particleKineticEnergies.begin();
            for( G4INCL::ParticleIter i = outgoingParticles.begin(), e = outgoingParticles.end(); i!=e; ++i, ++iP, ++iE)
            {
              const double mass = (*i)->getMass();
              const double newKineticEnergy = (*iE) * rescale;

              (*i)->setMomentum(*iP);
              (*i)->setEnergy(mass + newKineticEnergy);
              (*i)->adjustMomentumFromEnergy();

              pBalance -= (*i)->getMomentum();
            }

            nucleus->setMomentum(pBalance);
            const double remnantMass = G4INCL::ParticleTable::getTableMass(nucleus->getA(),nucleus->getZ(),nucleus->getS()) + nucleus->getExcitationEnergy();
            const double pRem2 = pBalance.mag2();
            const double recoilEnergy = pRem2/
              (std::sqrt(pRem2+remnantMass*remnantMass) + remnantMass);
            nucleus->setEnergy(remnantMass + recoilEnergy);
          }
      };

      /// \brief Class to adjust remnant recoil in the reaction CM system
      class RecoilCMFunctor : public G4INCL::RootFunctor {
        public:
          /** \brief Prepare for calling the () operator and scaleParticleEnergies
           *
           * The constructor sets the private class members.
           */
          RecoilCMFunctor(G4INCL::Nucleus * const n, const G4INCL::EventInfo &ei, const G4INCL::ThreeVector mom) :
            RootFunctor(0., 1E6),
            nucleus(n),
            transferQ(mom),
            theIncomingMomentum(nucleus->getIncomingMomentum()),
            outgoingParticles(n->getStore()->getOutgoingParticles()),
            theEventInfo(ei) {
              if(theIncomingMomentum.mag() == 0.){ //change the condition
                thePTBoostVector = {0., 0., 0.};
                //INCL_WARN("PTBoostVector at rest is zero" << '\n');      
              }
              else{
                thePTBoostVector = nucleus->getIncomingMomentum()/(nucleus->getInitialEnergy()); //D
                //INCL_WARN("PTBoostVector" << '\n');
              }
              for(G4INCL::ParticleIter p=outgoingParticles.begin(), e=outgoingParticles.end(); p!=e; ++p) {
                (*p)->boost(thePTBoostVector);
                particleCMMomenta.push_back((*p)->getMomentum());
              }
              G4INCL::ProjectileRemnant * const aPR = n->getProjectileRemnant();
              if(aPR && aPR->getA()>0) {
                aPR->boost(thePTBoostVector);
                particleCMMomenta.push_back(aPR->getMomentum());
                outgoingParticles.push_back(aPR);
              }
            }
          virtual ~RecoilCMFunctor() {}

          /** \brief Compute the energy-conservation violation.
           *
           * \param x scale factor for the particle energies
           * \return the energy-conservation violation
           */
          double operator()(const double x) const {
            scaleParticleCMMomenta(x);
            return nucleus->getConservationBalance(theEventInfo,true).energy;
          }

          /// \brief Clean up after root finding
          void cleanUp(const bool success) const {
            if(!success)
              scaleParticleCMMomenta(1.);
          }

        private:
          /// \brief Pointer to the nucleus
          G4INCL::Nucleus *nucleus;
          /// \brief Projectile-target CM boost vector
          G4INCL::ThreeVector thePTBoostVector;
          /// \brief Incoming momentum
          G4INCL::ThreeVector theIncomingMomentum;
          /// \brief List of final-state particles.
          G4INCL::ParticleList outgoingParticles;
          /// \brief Reference to the EventInfo object
          G4INCL::EventInfo const &theEventInfo;
          /// \brief Initial CM momenta of the outgoing particles
          std::list<G4INCL::ThreeVector> particleCMMomenta;
          /// \brief Transfer momentum
          G4INCL::ThreeVector transferQ;

          /** \brief Scale the kinetic energies of the outgoing particles.
           *
           * \param rescale scale factor
           */
          void scaleParticleCMMomenta(const double rescale) const {
            // Rescale the CM momenta of the outgoing particles.
            // theIncomingMomentum is zero; I put the transfer momentun 
            // as the start of cascade.
            G4INCL::ThreeVector remnantMomentum = theIncomingMomentum;
            remnantMomentum += transferQ;
            // LOG("INCLCascadeIntranuke", pINFO) << remnantMomentum.print();
            std::list<G4INCL::ThreeVector>::const_iterator iP = particleCMMomenta.begin();
            for( G4INCL::ParticleIter i = outgoingParticles.begin(), e = outgoingParticles.end(); i!=e; ++i, ++iP)
            {
              (*i)->setMomentum(*iP * rescale);
              (*i)->adjustEnergyFromMomentum();
              (*i)->boost(-thePTBoostVector);

              remnantMomentum -= (*i)->getMomentum();
            }

            nucleus->setMomentum(remnantMomentum);
            const double remnantMass = G4INCL::ParticleTable::getTableMass(nucleus->getA(),nucleus->getZ(),nucleus->getS()) + nucleus->getExcitationEnergy();
            const double pRem2 = remnantMomentum.mag2();
            const double recoilEnergy = pRem2/
              (std::sqrt(pRem2+remnantMass*remnantMass) + remnantMass);
            nucleus->setEnergy(remnantMass + recoilEnergy);
          }
      };

      /** \brief Rescale the energies of the outgoing particles.
       *
       * Allow for the remnant recoil energy by rescaling the energy (and
       * momenta) of the outgoing particles.
       */
      void rescaleOutgoingForRecoil() const;

  };

}

#endif
#endif // __GENIE_INCL_ENABLED__
