include make.inc
#################################################
### Defining Compiling variables              ###
#################################################
#-Wno-parentheses
CXX           = c++ -Wno-deprecated -Wall -std=c++17
#-m32  -Wno-deprecated
LD            = g++ -Wno-deprecated -Wall -std=c++17
#-m32  -Wno-deprecated 
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

SHAREDLIBFLAGS =
SHAREDCXXFLAGS =-I./CommonUtils 
SUBDIRS = CommonUtils/

TAUSPINERFLAG = 

override CXXFLAGS += 
ALLCXXFLAGS = $(CXXFLAGS) -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lRooFitCore -lGenVector -lTMVA -lPostscript -lMatrix -lPhysics -lGui -LObj -lMinuit2 -lThread -lTreePlayer -lRIO -lMathCore -rdynamic -Wl,--rpath $(ROOTSYS)/lib

# -ldl -lpthread

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

TARGETS       =	Ntuple_Controller \
		NtupleReader \
		Parameters \
		Selection_Base \
		Selection_Factory \
		Selection \
		Example \
		Tables \
		Plots \
		Tools \
		SkimConfig \
		TauDecay \
		DataMCType \
		PDGInfo \
		HistoConfig \
                FakeFactors \
		DoubleEventRemoval \
		ImpactParameter \
		BDTClassification \
		Logger

CINTTARGETS = 



# now include conditional options set by flags in make.inc
ifdef USE_TauSpinner
	SUBDIRS += CommonUtils/TauPolSoftware/TauDecaysInterface/
	SHAREDLIBFLAGS +=	-L./TauSpiner/lhapdf/lib -lLHAPDF \
				-L./TauSpiner/TAUOLA/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaTauSpinner \
				-L./TauSpiner/hepmc/install/lib -lHepMC
        SHAREDCXXFLAGS +=	-I./TauSpiner/hepmc/install/include \
				-I./TauSpiner/TAUOLA/include/Tauola/ \
				-I./TauSpiner/TAUOLA/TauSpinner/include/TauSpinner \
				-I./TauSpiner/TAUOLA/include/TauSpinner/ \
				-I./TauSpiner/TAUOLA/TauSpinner/examples/ \
				-I./TauSpiner/TAUOLA/include/ \
				-I./TauSpiner/lhapdf/include/ \
				-I./TauSpiner/pythia8/176/include/
	DEFS += -DUSE_TauSpinner=1
endif

ifdef USE_SVfit
#	TARGETS += SVfitProvider DataStorage SVFitStorage
#	CINTTARGETS += SVFitObject
#	SHAREDLIBFLAGS += -L./CommonUtils/lib -lSVfit
#	DEFS += -DUSE_SVfit=1
endif

ifdef USE_cherepanov
	TARGETS += MyTest ZTauHTauH SingleMuSkim ZTauMuTauH TTBar ZTauTau ZMuTau HTauTau TemplateStudyMT
	DEFS += -DUSE_cherepanov=1
endif


ifdef USE_goe
	TARGETS += MyTest
	DEFS += -DUSE_goe=1
endif

ifdef USE_lebihan
	TARGETS += 
	DEFS += -DUSE_lebihan=1
endif

ifdef USE_gbourgat
	TARGETS += HCPTauTau HCPMuTau HCPPiTau
	DEFS += -DUSE_gbourgat=1
endif

ifdef USE_cgrimault
	TARGETS += 
	DEFS += -DUSE_cgrimault=1
endif

ifdef USE_msessini
	TARGETS += HCPMuTau
	DEFS += -DUSE_msessini=1
endif



HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cxx)  Analysis.cxx
OBJS          = $(SRCS:.cxx=.o)

PROGRAM       = Analysis.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"

VPATH = utilities:i386_linux

ifdef USE_cherepanov
	vpath %.cxx cherepanov
	vpath %.h cherepanov
endif
ifdef USE_goe
	vpath %.cxx goe
	vpath %.h goe
endif
ifdef USE_lebihan
	vpath %.cxx lebihan
	vpath %.h lebihan
endif


ifdef USE_gbourgat
	vpath %.cxx gbourgat 
	vpath %.h gbourgat
endif

ifdef USE_cgrimault
	vpath %.cxx cgrimault 
	vpath %.h cgrimault
endif

ifdef USE_msessini
	vpath %.cxx msessini 
	vpath %.h msessini
endif


vpath %.cxx CommonFiles
vpath %.h CommonFiles

#	@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o ../$(PROGRAM)

##########################################################
###                                                    ###
##########################################################

$(OBJS): %.o : %.cxx
	$(CXX) $(ALLCXXFLAGS) $(DEFS) $< -o i386_linux/$@ 

.PHONY: clean cleanall cleandf all dataformats install sharedlib 

install: Analysis.exe


clean:
	@rm i386_linux/*.o
	@rm Analysis.exe

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@rm i386_linux/*.o
	@rm Analysis.exe

all: install

sharedlib:
	@cd CommonUtils/TauPolSoftware/TauDecaysInterface/; make; cp bin/libTDCILib.so ../../lib; cd ../../../;
	$(install)


