# Detect OS and Architecture
OS := $(shell uname -s)
ifeq ($(findstring CYGWIN, $(OS)),CYGWIN)
    OS := Windows
endif

# Check for external gsl repository needed on Windows
ifeq ($(OS), Windows)
  GSL_REPO = $(wildcard ../gsl)
  ifeq ($(GSL_REPO),)
    $(error GSL source code not found. Run 'git clone https://github.com/rtsoliday/gsl.git' next to the elegant repository)
  endif
endif

# Check for external SDDS repository
SDDS_REPO = $(wildcard ../SDDS)
ifeq ($(SDDS_REPO),)
  $(error SDDS source code not found. Run 'git clone https://github.com/rtsoliday/SDDS.git' next to the elegant repository)
endif

include Makefile.rules

DIRS = $(GSL_REPO)
DIRS += $(SDDS_REPO)/meschach
DIRS += $(SDDS_REPO)/zlib
DIRS += $(SDDS_REPO)/lzma
DIRS += $(SDDS_REPO)/mdblib
DIRS += $(SDDS_REPO)/mdbmth
DIRS += $(SDDS_REPO)/rpns/code
DIRS += $(SDDS_REPO)/namelist
DIRS += $(SDDS_REPO)/SDDSlib
DIRS += $(SDDS_REPO)/fftpack
DIRS += $(SDDS_REPO)/matlib
DIRS += $(SDDS_REPO)/mdbcommon
ifneq ($(MPI_CC),)
DIRS += $(SDDS_REPO)/pgapack
endif
DIRS += physics
DIRS += xraylib
DIRS += src
DIRS += src/elegantTools
DIRS += src/sddsbrightness

.PHONY: all $(DIRS) clean distclean

all: $(DIRS)

ifneq ($(GSL_REPO),)
  $(GSL_REPO):
	$(MAKE) -C $@ -f Makefile.MSVC all
endif
$(SDDS_REPO)/meschach:
	$(MAKE) -C $@
$(SDDS_REPO)/zlib: $(SDDS_REPO)/meschach
	$(MAKE) -C $@
$(SDDS_REPO)/lzma: $(SDDS_REPO)/zlib
	$(MAKE) -C $@
$(SDDS_REPO)/mdblib: $(SDDS_REPO)/lzma
	$(MAKE) -C $@
$(SDDS_REPO)/mdbmth: $(SDDS_REPO)/mdblib
	$(MAKE) -C $@
$(SDDS_REPO)/rpns/code: $(SDDS_REPO)/mdbmth $(GSL_REPO)
	$(MAKE) -C $@
$(SDDS_REPO)/namelist: $(SDDS_REPO)/rpns/code
	$(MAKE) -C $@
ifeq ($(MPI_CC),)
$(SDDS_REPO)/SDDSlib: $(SDDS_REPO)/namelist
	$(MAKE) -C $@
else
$(SDDS_REPO)/SDDSlib: $(SDDS_REPO)/namelist
	$(MAKE) -C $@
	$(MAKE) -C $@ -f Makefile.mpi
endif
$(SDDS_REPO)/fftpack: $(SDDS_REPO)/SDDSlib
	$(MAKE) -C $@
$(SDDS_REPO)/matlib: $(SDDS_REPO)/fftpack
	$(MAKE) -C $@
$(SDDS_REPO)/mdbcommon: $(SDDS_REPO)/matlib
	$(MAKE) -C $@
ifneq ($(MPI_CC),)
$(SDDS_REPO)/pgapack: $(SDDS_REPO)/mdbcommon
	$(MAKE) -C $@
endif
physics: $(SDDS_REPO)/mdbcommon
	$(MAKE) -C $@
xraylib: physics
	$(MAKE) -C $@
ifeq ($(MPI_CC),)
src: xraylib
	$(MAKE) -C $@
else
src: xraylib
	$(MAKE) -C $@
	$(MAKE) -C $@ -f Makefile.mpi
endif
src/elegantTools: src
	$(MAKE) -C $@
src/sddsbrightness: src/elegantTools
	$(MAKE) -C $@

clean:
	$(MAKE) -C physics clean
	$(MAKE) -C xraylib clean
	$(MAKE) -C src clean
	$(MAKE) -C src/elegantTools clean
	$(MAKE) -C src/sddsbrightness clean

distclean: clean
	rm -rf bin/$(OS)-$(ARCH)
	rm -rf lib/$(OS)-$(ARCH)
