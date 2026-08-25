//____________________________________________________________________________
/*!

  \class    genie::INCLNucleus

  \brief    INCLXX nuclear model. Implements the NuclearModelI 
  interface.

  \ref      

  \author   Liang Liu, (liangliu@fnal.gov)

  \created  Oct. 2024

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________


#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _INCL_NUCLEUS_H_
#define _INCL_NUCLEUS_H_

// ROOT
#include "TVector3.h"


// INCL++
// For configuration
#include "G4INCLConfig.hh"
#include "G4INCLNucleus.hh"
// GENIE
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Interaction/Target.h"
#include "Physics/NuclearState/NuclearModel.h"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLCascadeAction.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"
#include "G4INCLRootFinder.hh"
#include "G4INCLParticleTable.hh"


namespace genie {

  class INCLNucleus {

    public: 
      static INCLNucleus * Instance (void);
      void initialize(const Target * tgt);
      void reset(const Target * tgt);
      void configure();

      TVector3 getHitNucleonPosition();
      TVector3 getHitNucleonMomentum();
      double   getHitNucleonEnergy();
      double   getHitNucleonMass();
      double   getMass();
      double   getRemovalEnergy();
      G4INCL::Nucleus * getNuclues();
      G4INCL::StandardPropagationModel * getPropagationModel();

      void setHitParticle(const int pdg, TVector3 &posi);
      void setHitNNCluster(const int pdg1, const int pdg2, TVector3 &posi);
      G4INCL::Particle *getHitParticle();
      std::shared_ptr<G4INCL::Cluster>  getHitNNCluster();

      G4INCL::Config *getConfig(){return theConfig_;}

//      bool     nextNucleonIndex(int pdg){
//	if(pdg::IsProton(pdg)){
//	  nucleon_index_++;
//	  if(nucleon_index_ >= nucleus_->getZ())
//	    return false;
//	}
//	else{
//	  nucleon_index_++;
//	  if(nucleon_index_ >= nucleus_->getA())
//	    return false;
//	}
//	hitNucleon_ = nucleus_->getStore()->getParticles().at(nucleon_index_);
//	return true;
//      }
//      void     initNucleonIndex(int pdg){
//	if(pdg::IsProton(pdg)){
//	  nucleon_index_ = 0;
//	}
//	else{
//	  nucleon_index_ = nucleus_->getZ();
//	}
//	hitNucleon_ = nucleus_->getStore()->getParticles().at(nucleon_index_);
//      }

      double getMaxUniverseRadius() {return maxUniverseRadius_;}

      void setINCLXXDataFilePath(std::string str){ INCLXXDataFilePath_ = str; }
      void setABLAXXDataFilePath(std::string str){ ablaxxDataFilePath_ = str; }
      void setABLA07DataFilePath(std::string str){ abla07DataFilePath_ = str; }
      void setGEMINIXXDataFilePath(std::string str){ geminixxDataFilePath_ = str; }
      void setDeExcitationType(G4INCL::DeExcitationType deExType){ deExcitationType_ = deExType; }
      void setPotentialType(G4INCL::PotentialType type) { potentialType_ = type; }
      void setPauliType(G4INCL::PauliType type) { pauliType_ = type; }
      void setPauliString(std::string str) {pauliString_ = str;}

      void setLocalEnergyBBType(G4INCL::LocalEnergyType type) {localEnergyTypeBB_ = type;}
      void setLocalEnergyPiType(G4INCL::LocalEnergyType type) {localEnergyTypePi_ = type;}
      void setHadronizationTime(const double t) { hadronizationTime_=t; }


      bool isRPValid(double r, double p);

      void ResamplingHitNucleon();
      TVector3 ResamplingVertex(const int pdg);

      void setHybridModel(NuclearModel_t model){
        model_type_ = model;
      }

      NuclearModel_t getHybridModel(){
        return model_type_;
      }

    private:
      INCLNucleus();
      ~INCLNucleus();
      void initUniverseRadius(const int A, const int Z);
      G4INCL::Particle* getNucleon(const int pdg);
      std::shared_ptr<G4INCL::Cluster> getNNCluster(const int pdg1, const int pdg2);
      static INCLNucleus *fInstance;

      //TVector3 v3_; // position of initial nucleon 
      //TVector3 p3_; // fermi momentum of initial nucleon
      double  energy_; // off-shell energy of initial nucleon
      G4INCL::Config *theConfig_;
      G4INCL::Nucleus *nucleus_;
      G4INCL::Particle *hitNucleon_;
      // NN cluster for MEC channel
      std::shared_ptr<G4INCL::Cluster>  clusterNN_;
//      const int maxClusterMass = 2;
//      G4INCL::Particle *selectedParticles[maxClusterMass];



      G4INCL::StandardPropagationModel *propagationModel_;
      G4INCL::CascadeAction *cascadeAction_;
      const G4INCL::NuclearDensity *theDensityForLepton;
      const G4INCL::NuclearDensity *theDensity;
      const G4INCL::NuclearPotential::INuclearPotential *thePotential;

      int nucleon_index_;
      int cluster_index1_;
      int cluster_index2_;

      double maxUniverseRadius_;
      // double maxInteractionDistance_;
      double minRemnantSize_;
      double hadronizationTime_;

      std::string INCLXXDataFilePath_;
      std::string abla07DataFilePath_;
      std::string ablaxxDataFilePath_;
      std::string geminixxDataFilePath_;

      G4INCL::DeExcitationType deExcitationType_;
      G4INCL::PotentialType potentialType_;
      G4INCL::PauliType pauliType_;
      std::string pauliString_;

      G4INCL::LocalEnergyType localEnergyTypeBB_;
      G4INCL::LocalEnergyType localEnergyTypePi_;

      NuclearModel_t model_type_;
  };

}         // genie namespace
#endif    // _INCL_NUCLEUS_H_
#endif // __GENIE_INCL_ENABLED__
