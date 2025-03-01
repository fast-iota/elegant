SDDS_REPO = $(wildcard ../../../SDDS)
ifeq ($(SDDS_REPO),)
  $(error SDDS source code not found. Run 'git clone https://github.com/rtsoliday/SDDS.git' next to the SDDS-Python-module repository)
endif

DD=../
include ../../Makefile.rules

SCRIPTS = brightnessEnvelope

PROD = sddsbrightness sddsurgent sddsfluxcurve sddsws sddsbs sddsecon sddscompton sddsxra sddsxrf xrltool xrltest

sddsbrightness_SRC = sddsbrightness.c usbC.c brighteC.c huntC.c
sddsurgent_SRC = sddsurgent.c
sddsfluxcurve_SRC = sddsfluxcurve.c

CFLAGS +=  -I$(SDDS_REPO)/include -I../../physics -I../../xraylib -I$(OBJ_DIR)

ifeq ($(OS), Linux)
  PROD += sddssyncflux
  CFLAGS += -fopenmp -DUSE_GSL
  LDFLAGS := -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) -L../../lib/$(OS)-$(ARCH) -fopenmp $(LDFLAGS)
  PROD_SYS_LIBS := $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS = -lxrl -loagphy -lmdbcommon -lSDDS1 -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Darwin)
  PROD += sddssyncflux
  CFLAGS += -DUSE_GSL
  LDFLAGS := -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) -L../../lib/$(OS)-$(ARCH) $(LDFLAGS)
  PROD_SYS_LIBS := $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS = -lxrl -loagphy -lmdbcommon -lSDDS1 -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Windows)
  CFLAGS += -I$(SDDS_REPO)/lzma -wd4244
  PROD_LIBS = xrl.lib oagphy.lib mdbcommon.lib SDDS1.lib rpnlib.lib mdbmth.lib mdblib.lib lzma.lib z.lib
  LIB_LINK_DIRS = -LIBPATH:$(SDDS_REPO)/lib/$(OS)-$(ARCH) -LIBPATH:../../lib/$(OS)-$(ARCH)
endif

include ../../Makefile.build

$(OBJ_DIR)/sddsbrightness$(EXEEXT): $(sddsbrightness_OBJS) $(PROD_DEPS) $(OBJ_DIR)/usb.$(OBJEXT) $(OBJ_DIR)/hunt.$(OBJEXT) $(OBJ_DIR)/brighte.$(OBJEXT)
	$(F77) $(RPATH) $(OUTPUTEXE) $(sddsbrightness_OBJS) $(OBJ_DIR)/usb.$(OBJEXT) $(OBJ_DIR)/hunt.$(OBJEXT) $(OBJ_DIR)/brighte.$(OBJEXT) $(LDFLAGS) $(LIB_LINK_DIRS) $(PROD_LIBS) $(TIFF_LIBS) $(PROD_SYS_LIBS)
	cp -f $@ $(BIN_DIR)/

$(OBJ_DIR)/sddsurgent$(EXEEXT): $(sddsurgent_OBJS) $(PROD_DEPS) $(OBJ_DIR)/us.$(OBJEXT) $(OBJ_DIR)/hunt.$(OBJEXT) $(OBJ_DIR)/brighte.$(OBJEXT) $(OBJ_DIR)/urgent.$(OBJEXT)
	$(F77) $(RPATH) $(OUTPUTEXE) $(sddsurgent_OBJS) $(OBJ_DIR)/us.$(OBJEXT) $(OBJ_DIR)/hunt.$(OBJEXT) $(OBJ_DIR)/brighte.$(OBJEXT) $(OBJ_DIR)/urgent.$(OBJEXT) $(LDFLAGS) $(LIB_LINK_DIRS) $(PROD_LIBS) $(TIFF_LIBS) $(PROD_SYS_LIBS)
	cp -f $@ $(BIN_DIR)/

$(OBJ_DIR)/sddsfluxcurve$(EXEEXT): $(sddsfluxcurve_OBJS) $(PROD_DEPS) $(OBJ_DIR)/us.$(OBJEXT) $(OBJ_DIR)/hunt.$(OBJEXT) $(OBJ_DIR)/brighte.$(OBJEXT)
	$(F77) $(RPATH) $(OUTPUTEXE) $(sddsfluxcurve_OBJS) $(OBJ_DIR)/us.$(OBJEXT) $(OBJ_DIR)/hunt.$(OBJEXT) $(OBJ_DIR)/brighte.$(OBJEXT) $(LDFLAGS) $(LIB_LINK_DIRS) $(PROD_LIBS) $(TIFF_LIBS) $(PROD_SYS_LIBS)
	cp -f $@ $(BIN_DIR)/

$(OBJ_DIR)/usb.$(OBJEXT): usb.f
	$(F77) usb.f -c -o $(OBJ_DIR)/usb.$(OBJEXT)	

$(OBJ_DIR)/us.$(OBJEXT): us.f
	$(F77) us.f -c -o $(OBJ_DIR)/us.$(OBJEXT)	

$(OBJ_DIR)/hunt.$(OBJEXT): hunt.f
	$(F77) hunt.f -c -o $(OBJ_DIR)/hunt.$(OBJEXT)	

$(OBJ_DIR)/brighte.$(OBJEXT): brighte.f
	$(F77) brighte.f -c -o $(OBJ_DIR)/brighte.$(OBJEXT)	

$(OBJ_DIR)/urgent.$(OBJEXT): urgent.f
	$(F77) urgent.f -c -o $(OBJ_DIR)/urgent.$(OBJEXT)	

