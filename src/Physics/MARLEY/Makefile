#
# Makefile for GENIE Neutrino Generator
#
# Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

PACKAGE           = Physics/MARLEY
PACKAGE_ABBREV    = PhMARLEY
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

# Only need these here ... not in every package

ifeq ($(strip $(GOPT_ENABLE_MARLEY)),YES)
  # extra flags, include paths, and libraries to link to
  CXXFLAGS               += $(MARLEY_FLAGS)
  CPP_INCLUDES           += $(MARLEY_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(MARLEY_INCLUDES)
  EXTRA_EXT_LIBS         += $(MARLEY_LIBRARIES)
else
  $(info $(PACKAGE) not build against MARLEY)
endif

# MARLEY's own headers (marley/NuclearReaction.hh) unconditionally include
# HepMC3/*.h, so this package needs the HepMC3 include path too, regardless
# of whether any GENIE code in this package uses HepMC3 directly.
ifeq ($(strip $(GOPT_ENABLE_HEPMC3_INTERFACE)),YES)
  CXXFLAGS               += $(HEPMC3_FLAGS)
  CPP_INCLUDES           += $(HEPMC3_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(HEPMC3_INCLUDES)
  EXTRA_EXT_LIBS         += $(HEPMC3_LIBRARIES)
else
  $(info $(PACKAGE) not built against HepMC3)
endif

FORCE:
