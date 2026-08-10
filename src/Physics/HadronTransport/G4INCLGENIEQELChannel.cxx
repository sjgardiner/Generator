#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__
#include "G4INCLGENIEQELChannel.h"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  GENIEQELChannel::GENIEQELChannel(Particle *p, Nucleus *n, std::vector<GENIEParticleRecord> *eventRecord)
    :hitParticle(p), theNucleus(n),
    genie_evtrec(eventRecord)
  {
  }

  GENIEQELChannel::~GENIEQELChannel()
  {
  }

  void GENIEQELChannel::fillFinalState(FinalState *fs)
  {
    std::vector<GENIEParticleRecord>::iterator ip;
    for(ip = genie_evtrec->begin(); ip != genie_evtrec->end(); ip++){
      if(ip->Status() == 14){
        ip->setID(int(hitParticle->getID()));
        hitParticle->setType(ip->Type());
        hitParticle->setMomentum(ip->P3());
        std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << "  " << hitParticle->getPosition().print() << std::endl;
        hitParticle->setPosition(ip->X3());
        std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << "  " << hitParticle->getPosition().print() << std::endl;
        std::cout << "DEBUG: " << __FILE__ << ":" << __LINE__ << " position radius difference: " << hitParticle->getPosition().mag() - ip->X3().mag() << std::endl;
        hitParticle->adjustEnergyFromMomentum();
      }
    }
    fs->addModifiedParticle(hitParticle);
  }
}

#endif // __GENIE_INCL_ENABLED__
