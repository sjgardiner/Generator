//____________________________________________________________________________
/*!

\class    genie::HepMC3Converter

\brief    Converts genie::EventRecord objects into HepMC3::GenEvent objects
          following the NuHepMC standard (https://github.com/NuHepMC/Spec)
	  Also allows the inverse direction, 
	  HepMC3::GenEvent to genie::EventRecord.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 3, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HEPMC3_CONVERTER_H_
#define _HEPMC3_CONVERTER_H_

#include "Framework/Conventions/GBuild.h"
#include "Framework/GHEP/GHepStatus.h"
#include "HepMC3/Attribute.h"
#include <cmath>
#ifdef __GENIE_HEPMC3_INTERFACE_ENABLED__
// Forward-declare needed HepMC3 classes here
namespace HepMC3 {
  class GenEvent;
  class GenRunInfo;
  class Attribute;

// Class that stores doubles with more precision.
class PreciseDoubleAttribute : public Attribute {

  public:

    
  PreciseDoubleAttribute(): Attribute(), m_val(0.0) {}

  /** @brief Constructor initializing attribute value */
  PreciseDoubleAttribute(double val): Attribute(), m_val(val) {}

  /** @brief Implementation of Attribute::from_string */
  bool from_string(const std::string &att)  override {
    m_val = atof( att.c_str() );
    set_is_parsed(true);
    return true;
  }

  /** @brief Implementation of Attribute::to_string */
  bool to_string(std::string &att) const  override {
    std::ostringstream oss;
    oss << std::setprecision(16)
  	  << m_val;
    att = oss.str();
    return true;
  }

  /** @brief get the value associated to this Attribute. */
  double value() const {
    return m_val;
  }

  /** @brief set the value associated to this Attribute. */
  void set_value(const double& d) {
    m_val = d;
    set_is_parsed(true);
  }

  private:

  double m_val; ///< Attribute value

};

// Vector variant of PreciseDoubleAttribute.
class PreciseVectorDoubleAttribute : public Attribute {
  public:

  /** @brief Default constructor */
  PreciseVectorDoubleAttribute():Attribute(),m_val() {}

  /** @brief Constructor initializing attribute value */
  PreciseVectorDoubleAttribute(std::vector<double> val):Attribute(),m_val(val) {}

  /** @brief Implementation of Attribute::from_string */
  bool from_string(const std::string &att) override {
    double  datafoo;
    m_val.clear();
    std::stringstream datastream(att);
    while (datastream >> datafoo) m_val.emplace_back(datafoo);
    set_is_parsed(true);
    return true;
  }

  /** @brief Implementation of Attribute::to_string */
  bool to_string(std::string &att) const  override {
    att.clear();
    for (const auto& a:  m_val) {
  	if (att.length()) att+=" ";  
  	//att+=std::to_string(a); // insufficient precision
  	std::ostringstream oss;
  	oss << std::setprecision(16) << a;
  	att += oss.str();
    }
    return true;
  }

  /** @brief get the value associated to this Attribute. */
  std::vector<double> value() const {
    return m_val;
  }

  /** @brief set the value associated to this Attribute. */
  void set_value(const std::vector<double>& i) {
    m_val = i;
    set_is_parsed(true);
  }

  private:
  std::vector<double> m_val; ///< Attribute value
};
}

namespace genie {

class EventRecord;
class GHepParticle;
class Interaction;

class HepMC3Converter {

public:

  HepMC3Converter(void);

  std::shared_ptr< HepMC3::GenEvent > ConvertToHepMC3(
    const genie::EventRecord& gevrec );

  std::shared_ptr< genie::EventRecord > RetrieveGHEP(
    const HepMC3::GenEvent& evt );

  int GetNuHepMCParticleStatus( const genie::GHepParticle* gpart,
    const genie::EventRecord& gevrec ) const;

  int GetNuHepMCProcessID( const genie::Interaction& inter ) const;

  genie::GHepStatus_t GetGHepParticleStatus( int nuhepmc_status ) const;

protected:

  void StoreInteraction( const genie::Interaction& inter,
    HepMC3::GenEvent& evt );

  genie::Interaction* RetrieveInteraction( const HepMC3::GenEvent& evt );

  void PrepareRunInfo( const genie::EventRecord* gevrec );

  std::shared_ptr< HepMC3::GenRunInfo > fRunInfo;

private:
  
  bool fTuneLoaded;

};

} // genie namespace

#endif // __GENIE_HEPMC3_INTERFACE_ENABLED__
#endif // _HEPMC3_CONVERTER_H_
