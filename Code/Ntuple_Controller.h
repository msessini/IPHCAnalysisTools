//Ntuple_Controller.h HEADER FILE

#ifndef Ntuple_Controller_h
#define Ntuple_Controller_h


// Root include files
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TSystem.h"

// Include files (C & C++ libraries)
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <utility>      // std::pair
#include <tuple>
#include <functional>

#include "NtupleReader.h"

#include "HistoConfig.h"
#ifdef USE_TauSpinner
#include "TauSpinerInterface.h"
#endif
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay.h"

#ifdef USE_SVfit
#include "DataFormats/SVFitObject.h"
#include "SVFitStorage.h"
#include "SVfitProvider.h"
#endif


#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"

#include "PileUp.h"
//#include "TauTriggerSFs/interface/TauTriggerSFs2017.h"
//#include "TauTriggerSFs/interface/SFProvider.h"
#include "TauIDSFs/interface/TauIDSFTool.h"
//#include "TauIDSFTool.h"
#include "RecoilCorrector.h"
#include "MEtSys.h"

#include "RooWorkspace.h"
#include "RooFunctor.h"
#include <memory>
#include <unistd.h>

// Rochester muon momentum correction
//#include "CommonFiles/rochcor2012jan22.h"

// small struct needed to allow sorting indices by some value
struct sortIdxByValue {
  bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second > right.second;
  }
};

///////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//*
//*   Class: Ntuple_Controller
//*   
//*   Purpose: The purpose of this class is to provide a interface to the
//*            Ntuple
//*
//*   Designed by: Vladimir Cherepanov
//*
//*
//*****************************************************************************
///////////////////////////////////////////////////////////////////////////////


class Ntuple_Controller{
 private:
  NtupleReader *Ntp;
  TFile *newfile;
  TTree *SkimmedTree;
  int nbytes;
  int jentry;
  int nb;
  bool copyTree;

  int currentEvent;

  bool cannotObtainHiggsMass; // avoid repeated printing of warning when running locally

  int EmbedID;

  // Ntuple Access Functions
  virtual void Branch_Setup(TString B_Name, int type);
  virtual void Branch_Setup(){}

  // Functions to configure objects
  virtual void ConfigureObjects(); 
  void doElectrons();
  void doPhotons();
  void doJets();
  void doMuons();
  void doTaus();
  void doMET();
  unsigned int ObjEvent;

  // helper functions for internal calculations
  void printMCDecayChain(unsigned int par, unsigned int level = 0, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false);

  // Object Variables
  std::vector<TLorentzVector> electrons_default;
  std::vector<TLorentzVector> photons_default;
  std::vector<TLorentzVector> jets_default;
  std::vector<TLorentzVector> muons_default;
  std::vector<TLorentzVector> taus_default;
  TLorentzVector              met_default;
  std::vector<TLorentzVector> electrons;
  std::vector<TLorentzVector> photons;
  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> muons;
  std::vector<TLorentzVector> taus;
  TLorentzVector              met;

  // TString flags for object corrections
  TString tauCorrection;
  TString muonCorrection;
  TString elecCorrection;
  TString jetCorrection;

  // Systematic controls variables
  int theSys;
  HistoConfig HConfig;

  // Interfaces
#ifdef USE_TauSpinner  
  TauSpinerInterface TauSpinerInt;
#endif
  HistoConfig HistoC;

  // Fit Variables
  LorentzVectorParticle               theTau;
  std::vector<LorentzVectorParticle>  daughter;
  double                              LC_chi2;
  double                              ndof;
  bool                                fitStatus;
  bool                                isInit;

  // muon correction related objects
  //  rochcor2012*   rmcor;
  std::vector<TLorentzVector> Muon_corrected_p4;
  void           CorrectMuonP4();
  bool           Muon_isCorrected;

  // helpers for SVFit
#ifdef USE_SVfit
  // create SVFitObject from standard muon and standard tau_h
  void runAndSaveSVFit_MuTauh(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, unsigned tauIdx, double scaleMu, double scaleTau, bool save = true);
  // create SVFitObject from standard muon and fully reconstructed 3prong tau
  void runAndSaveSVFit_MuTau3p(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, double scaleMu, double scaleTau, bool save = true);
  // create SVFitObject from standard tau_h and standard tau_h
  void runAndSaveSVFit_TauhTauh(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned tauIdx1, unsigned tauIdx2, double scaleTau1, double scaleTau2, bool save = false);
#endif

 public:
  // Constructor
  Ntuple_Controller(std::vector<TString> RootFiles, TString SysTree);

  // Destructor
  virtual ~Ntuple_Controller() ;

  // Event initializer
  void InitEvent();

  //TauSpiner function
  double TauSpinerGet(int SpinType,char* CPstate);
  void TauSpinerSetSignal(int signalcharge){
#ifdef USE_TauSpinner
    TauSpinerInt.SetTauSignalCharge(signalcharge);
#endif
  }
  
  enum beamspot{BS_x0,BS_y0,BS_z0,BS_sigmaZ,BS_dxdz,BS_dydz,BS_BeamWidthX,NBS_par};
  enum TrackQuality {
    undefQuality = -1, loose = 0, tight = 1, highPurity = 2,
    confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6,
    qualitySize = 7
  };
  enum TrackPar{i_qoverp = 0, i_lambda, i_phi, i_dxy,i_dsz};

  enum GenParticleFlag{Bit_isPrompt=0,
		       Bit_isDecayedLeptonHadron,
		       Bit_isTauDecayProduct,
		       Bit_isPromptTauDecayProduct,
		       Bit_isDirectTauDecayProduct,
		       Bit_isDirectPromptTauDecayProduct,
		       Bit_isDirectHadronDecayProduct,
		       Bit_isHardProcess,
		       Bit_fromHardProcess,
		       Bit_isHardProcessTauDecayProduct,
		       Bit_isDirectHardProcessTauDecayProduct,
		       Bit_fromHardProcessBeforeFSR,
		       Bit_isFirstCopy,
		       Bit_isLastCopy,
		       Bit_isLastCopyBeforeFSR,
		       Bit_isVBFParton};




  enum TauQualityBitMask{Bit_byLooseCombinedIsolationDeltaBetaCorr3Hits=0,
			 Bit_byMediumCombinedIsolationDeltaBetaCorr3Hits,
			 Bit_byTightCombinedIsolationDeltaBetaCorr3Hits,
			 Bit_againstMuonLoose3,
			 Bit_againstMuonTight3,
			 Bit_againstElectronVLooseMVA6,
			 Bit_againstElectronLooseMVA6,
			 Bit_againstElectronMediumMVA6,
			 Bit_againstElectronTightMVA6,
			 Bit_againstElectronVTightMVA6,
			 Bit_byVLooseIsolationMVArun2v1DBoldDMwLT,
			 Bit_byLooseIsolationMVArun2v1DBoldDMwLT,
			 Bit_byMediumIsolationMVArun2v1DBoldDMwLT,
			 Bit_byTightIsolationMVArun2v1DBoldDMwLT,
			 Bit_byVTightIsolationMVArun2v1DBoldDMwLT,
			 Bit_byVLooseIsolationMVArun2v1DBnewDMwLT,
			 Bit_byLooseIsolationMVArun2v1DBnewDMwLT,
			 Bit_byMediumIsolationMVArun2v1DBnewDMwLT,
			 Bit_byTightIsolationMVArun2v1DBnewDMwLT,
			 Bit_byVTightIsolationMVArun2v1DBnewDMwLT,
			 Bit_byLooseIsolationMVArun2v1DBdR03oldDMwLT,
			 Bit_byMediumIsolationMVArun2v1DBdR03oldDMwLT,
			 Bit_byTightIsolationMVArun2v1DBdR03oldDMwLT,
			 Bit_byVTightIsolationMVArun2v1DBdR03oldDMwLT,
			 Bit_byVLooseIsolationMVArun2017v1DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byLooseIsolationMVArun2017v1DBoldDMwLT2017,  //FRA syncApr2018
			 Bit_byMediumIsolationMVArun2017v1DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byTightIsolationMVArun2017v1DBoldDMwLT2017,  //FRA syncApr2018
			 Bit_byVTightIsolationMVArun2017v1DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byVLooseIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byLooseIsolationMVArun2017v2DBoldDMwLT2017,  //FRA syncApr2018
			 Bit_byMediumIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byTightIsolationMVArun2017v2DBoldDMwLT2017,  //FRA syncApr2018
			 Bit_byVTightIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byVVTightIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			 Bit_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, //FRA syncApr2018
			 Bit_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017,  //FRA syncApr2018
			 Bit_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017, //FRA syncApr2018
			 Bit_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017,  //FRA syncApr2018
			 Bit_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, //FRA syncApr2018
			 Bit_byVVVLooseDeepTau2017v2p1VSjet,
			 Bit_byVVLooseDeepTau2017v2p1VSjet, 
			 Bit_byVLooseDeepTau2017v2p1VSjet,  
			 Bit_byLooseDeepTau2017v2p1VSjet,   
			 Bit_byMediumDeepTau2017v2p1VSjet,  
			 Bit_byTightDeepTau2017v2p1VSjet,   
			 Bit_byVTightDeepTau2017v2p1VSjet,  
			 Bit_byVVTightDeepTau2017v2p1VSjet, 
			 Bit_byVVVLooseDeepTau2017v2p1VSe,  
			 Bit_byVVLooseDeepTau2017v2p1VSe, 
			 Bit_byVLooseDeepTau2017v2p1VSe,   
			 Bit_byLooseDeepTau2017v2p1VSe,	
			 Bit_byMediumDeepTau2017v2p1VSe,   
			 Bit_byTightDeepTau2017v2p1VSe,	
			 Bit_byVTightDeepTau2017v2p1VSe,   
			 Bit_byVVTightDeepTau2017v2p1VSe,   
			 Bit_byVLooseDeepTau2017v2p1VSmu, 
			 Bit_byLooseDeepTau2017v2p1VSmu, 
			 Bit_byMediumDeepTau2017v2p1VSmu, 
			 Bit_byTightDeepTau2017v2p1VSmu};





  enum MuonQualityBitMask{Bit_MuonLoose=0,
			  Bit_MuonSoft,
			  Bit_MuonMedium,
			  Bit_MuonTight,
			  Bit_MuonHighPt,
			  Bit_MuonTight_noVtx};



  enum particleType {
    MUON = 0,
    ELECTRON = 1,
    TAU =2
  };
  float m_MVAEleIDCuts[2][2][3] ;

  enum pairType {
    MuHad  = 0,
    EHad   = 1,
    HadHad = 2,
    MuMu   = 3,
    EE     = 4,
    EMu    = 5,
    EEPrompt = 6, // prompt Z->ee/mumu decays
    MuMuPrompt = 7,
    Other  = 8 // for e.g. h->bb
  };

  enum eleMVAIDWP {
    EMVATight = 0, // 80% eff
    EMVALoose = 1  // 90% eff
  };

  enum muIDWP {
    MuLoose  = 0,
    MuSoft   = 1,
    MuMedium = 2,
    MuTight  = 3,
    MuHighPt = 4
  };

  enum aeleWP {
    aeleVLoose = 0,
    aeleLoose  = 1,
    aeleMedium = 2,
    aeleTight  = 3,
    aeleVTight = 4
  };

  enum amuWP {
    amuLoose = 0,
    amuTight = 1
  };

  typedef std::vector<float> tauPair_t; // pt1 - iso1 - idx1 - pt2 - iso2 - idx2 - idxoriginalPair
  typedef std::tuple <float, float, int, float, float, int, int> tauPair_tuple; // pt1 - iso1 - idx1 - pt2 - iso2 - idx2 - idxoriginalPair
  // access to SVFit
#ifdef USE_SVfit
  SVFitObject* getSVFitResult_MuTauh(SVFitStorage& svFitStor, TString metType, unsigned muIdx, unsigned tauIdx, unsigned rerunEvery = 5000, TString suffix = "", double scaleMu = 1 , double scaleTau = 1);
  SVFitObject* getSVFitResult_MuTau3p(SVFitStorage& svFitStor, TString metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, TString suffix = "", double scaleMu = 1, double scaleTau = 1);
  SVFitObject* getSVFitResult_TauhTauh(SVFitStorage& svFitStor, TString metType, unsigned tauIdx1, unsigned tauIdx2, unsigned rerunEvery  = 5000 , TString suffix  ="" , double scaleTau1  =1 , double scaleTau2  =1 );


#endif


  // Ntuple Access Functions
  virtual Int_t Get_Entries();
  virtual void Get_Event(int _jentry);
  virtual Int_t Get_EventIndex();
  virtual TString Get_File_Name();

  //Ntuple Cloning Functions
  virtual void CloneTree(TString n);
  virtual void SaveCloneTree();
  inline void AddEventToCloneTree(){if(copyTree)SkimmedTree->Fill();}

  // Systematic controls
  enum    Systematic {Default=0,NSystematics};

  int     SetupSystematics(TString sys_);
  void    SetSysID(int sysid){theSys=sysid;}


  // Data/MC switch and thin
  TH1F* hLLRCounters;
  void ThinTree();

  // Set object corrections to be applied
  void SetTauCorrections(TString tauCorr){tauCorrection = tauCorr;}
  void SetMuonCorrections(TString muonCorr){muonCorrection = muonCorr;}
  void SetElecCorrections(TString elecCorr){elecCorrection = elecCorr;}
  void SetJetCorrections(TString jetCorr){jetCorrection = jetCorr;}
  // corresponding getters
  const TString& GetTauCorrections() const {return tauCorrection;}
  const TString& GetMuonCorrections() const {return muonCorrection;}
  const TString& GetElecCorrections() const {return elecCorrection;}
  const TString& GetJetCorrections() const {return jetCorrection;}

  // Information from input Ntuple path name
  TString GetInputDatasetName();
  TString GetInputPublishDataName();
  int getSampleHiggsMass();
  int readHiggsMassFromString(TString input);

  // resonance mass
  int getHiggsSampleMassFromGenInfo();
  double getResonanceMassFromGenInfo(bool useZ0 = true, bool useHiggs0 = true, bool useW = true);

  // Physics Variable Get Functions
  Long64_t GetMCID();
  // Event Variables
  ULong64_t EventNumber()   {return Ntp->_evt;}
  ULong64_t RunNumber()     {return Ntp->_runNumber;}
  Int_t Lumi()		    {return Ntp->_lumi;}
  Int_t Year()		    {return Ntp->_theYear;}
  bool isData()		    {return Ntp->_isData;}
  int tauIndex()	    {return Ntp->_tauIndex;}
  int tauGenMatch()	    {return Ntp->_tauGenMatch;}	
  float tauDM()	            {return Ntp->_tauDM;}
  Double_t tauPt()	    {return Ntp->_tauPt;}
  Double_t tauEta()         {return Ntp->_tauEta;}
  Double_t tauPhi()	    {return Ntp->_tauPhi;}
  Double_t tauE()           {return Ntp->_tauE;}
  float tauIPx()	    {return Ntp->_tauIPx;}
  float tauIPy()            {return Ntp->_tauIPy;}
  float tauIPz()            {return Ntp->_tauIPz;}
  float tauIPsignificance() {return Ntp->_tauIPsignificance;}
  int muIndex()	            {return Ntp->_muIndex;}
  int muGenMatch()	    {return Ntp->_muGenMatch;}
  Double_t muPt()           {return Ntp->_muPt;}
  Double_t muEta()	    {return Ntp->_muEta;}
  Double_t muPhi() 	    {return Ntp->_muPhi;}
  Double_t muE()	    {return Ntp->_muE;}
  double muIso()	    {return Ntp->_muIso;}
  float muIPx()		    {return Ntp->_muIPx;}
  float muIPy()		    {return Ntp->_muIPy;}
  float muIPz()		    {return Ntp->_muIPz;}
  float muIPsignificance()  {return Ntp->_muIPsignificance;}
  bool isOSpair()	    {return Ntp->_isOSpair;}
  bool isIso()		    {return Ntp->_isIso;}
  double pairvisMass()      {return Ntp->_pairvisMass;}
  int Njets()		    {return Ntp->_Njets;}
  int Nbjets()		    {return Ntp->_Nbjets;}
  double leadingjetPt()	    {return Ntp->_leadingjetPt;}
  double trailingjetPt()    {return Ntp->_trailingjetPt;}
  double dijetMass()   	    {return Ntp->_dijetMass;}
  double dijetPt()	    {return Ntp->_dijetPt;}
  double dijetdeltaEta()    {return Ntp->_dijetdeltaEta;}
  double ditauPt()	    {return Ntp->_ditauPt;}
  double ditauMass()	    {return Ntp->_ditauMass;}
  float muMETmt()	    {return Ntp->_muMETmt;}
  float PUPPImet()	    {return Ntp->_PUPPImet;}
  float PUPPImetphi() 	    {return Ntp->_PUPPImetphi;}
  float PUPPIMETCovXX()	    {return Ntp->_PUPPIMETCov00;}
  float PUPPIMETCovXY()     {return Ntp->_PUPPIMETCov10;}
  float PUPPIMETCovYY()     {return Ntp->_PUPPIMETCov11;}
  double pvPhiCP()	    {return Ntp->_pvPhiCP;}
  double dpPhiCP()	    {return Ntp->_dpPhiCP;}
  double wEven()	    {return Ntp->_wEven;}
  double wOdd()		    {return Ntp->_wOdd;}
  double wMM()		    {return Ntp->_wMM;}
  double wPrefiring()	    {return Ntp->_wPrefiring;}
  double wPrefiringUp()	    {return Ntp->_wPrefiringUp;}
  double wPrefiringDown()   {return Ntp->_wPrefiringDown;}
  double wIDvsJet()	    {return Ntp->_wIDvsJet;}
  double wIDvsJetUp()	    {return Ntp->_wIDvsJetUp;}
  double wIDvsJetDown()	    {return Ntp->_wIDvsJetDown;}
  double wIDvsEle()	    {return Ntp->_wIDvsEle;}
  double wIDvsEleUp()	    {return Ntp->_wIDvsEleUp;}
  double wIDvsEleDown()	    {return Ntp->_wIDvsEleDown;}
  double wIDvsMu()	    {return Ntp->_wIDvsMu;}
  double wIDvsMuUp()	    {return Ntp->_wIDvsMuUp;}
  double wIDvsMuDown()	    {return Ntp->_wIDvsMuDown;}
  double wTrg()	 	    {return Ntp->_wTrg;}
  double wTrgUp()	    {return Ntp->_wTrgUp;}
  double wTrgDown()	    {return Ntp->_wTrgDown;}
  double wIDMu()	    {return Ntp->_wIDMu;}
  double wTrkMu()	    {return Ntp->_wTrkMu;}
  double wPU()	  	    {return Ntp->_wPU;}
  double wZpT()		    {return Ntp->_wZpT;}
  double wZpTUp()	    {return Ntp->_wZpTUp;}
  double wZpTDown()	    {return Ntp->_wZpTDown;}
  double wToppT()	    {return Ntp->_wToppT;}
  double wToppTUp()	    {return Ntp->_wToppTUp;}
  double wToppTDown()	    {return Ntp->_wToppTDown;}
  double wBtag()	    {return Ntp->_wBtag;}
  double wBtagUp()	    {return Ntp->_wBtagUp;}
  double wBtagDown()	    {return Ntp->_wBtagDown;}
  double wPSISRUp()	    {return Ntp->_wPSISRUp;}
  double wPSISRDown()	    {return Ntp->_wPSISRDown;}
  double wPSFSRUp()	    {return Ntp->_wPSFSRUp;}
  double wPSFSRDown()	    {return Ntp->_wPSFSRDown;}
  double wScaleUp()	    {return Ntp->_wScaleUp;}
  double wScaleDown()	    {return Ntp->_wScaleDown;}
  double wMC()	  	    {return Ntp->_wMC;}
  double wSignal()	    {return Ntp->_wSignal;}
  float pvx()	 	    {return Ntp->_pvx;}
  float pvy()		    {return Ntp->_pvy;}
  float pvz()		    {return Ntp->_pvz;}
  float pvCov00()	    {return Ntp->_pvCov00;}
  float pvCov11()	    {return Ntp->_pvCov11;}
  float pvCov22() 	    {return Ntp->_pvCov22;}
  float pvCov01()           {return Ntp->_pvCov01;}
  float pvCov02()           {return Ntp->_pvCov02;}
  float pvCov12()           {return Ntp->_pvCov12;}
  Int_t Id()		    {return Ntp->_Id;}
  Int_t trueId()            {return Ntp->_trueId;}
  int Npartons()	    {return Ntp->_Npartons;}
  bool isZ()		    {return Ntp->_isZ;}
  bool isW()		    {return Ntp->_isW;}
  bool isH()		    {return Ntp->_isH;}
  bool isSignal()	    {return Ntp->_isSignal;}
  bool isQCD()		    {return Ntp->_isQCD;}
  bool isVV()		    {return Ntp->_isVV;}
  bool isTTbar()	    {return Ntp->_isTTbar;}
  bool isSingleTop()	    {return Ntp->_isSingleTop;}

};

#endif

