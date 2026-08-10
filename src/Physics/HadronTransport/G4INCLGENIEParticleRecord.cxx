#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"


namespace G4INCL{
  GENIEParticleRecord::GENIEParticleRecord(genie::GHepParticle *p, int scType, GENIERecordCode recordCode){
    fPdgCode        = p->Pdg();
    fCharge         = p->Charge()/3;
    // if(p->Status() == genie::kIStHadronInTheNucleus) fStatus = 1;
    // else fStatus = 0;
    fStatus         = p->Status();
    fFirstMother    = p->FirstMother();
    fLastMother     = p->LastMother();
    fFirstDaughter  = p->FirstDaughter();
    fLastDaughter   = p->LastDaughter();
    fScatteringType = scType;
    fPType          = this->PDG_to_INCLType(p->Pdg());
    fMass           = p->P4()->M()*1000.;
    fP3.set(p->P4()->Px()*1000., p->P4()->Py()*1000., p->P4()->Pz()*1000.);
    fX3.set(p->X4()->X(), p->X4()->Y(), p->X4()->Z());
    fID = -1;
    fGenieRecordCode = recordCode;


  }

  GENIEParticleRecord::~GENIEParticleRecord(){

  }
  void GENIEParticleRecord::setID(int id){
    fID = id;
  }
}
#endif // __GENIE_INCL_ENABLED__
