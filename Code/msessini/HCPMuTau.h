#ifndef HCPMuTau_h
#define HCPMuTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "boost/functional/hash.hpp"
//#include "SVFitStorage.h"
#include "TVector3.h"
//#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Objects.h"
//#include "PUReweight.h"
#include "FakeFactors.h"
#include <TMVA/Reader.h>

#include "RooWorkspace.h"
#include "RooFunctor.h"
#include <memory>

class HCPMuTau : public Selection {

 public:
  HCPMuTau(TString Name_, TString id_, char* Channel_);
  virtual ~HCPMuTau();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {Trigger=0,
	     //Id_and_Kin, 
	     /* NPairsFound, */
             //ZTTMC=0,
             //METFilters,
	     Id_and_Kin,
	     //genmatch,
	     TausIsolation, 
	     //Tau2Isolation,
	     AgainstEleMu,
	     LeptonVeto,
	     PairCharge, 
	     PairMass,
	     //MTM,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  char* Channel;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  //PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  //cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  
  //DataMCCorrections DataMC_Corr;
  //tauTrigSFreader tauTrgSF;
  //Int_t year;
  
  //ClassicSVfit svfitAlgo1;
  //ClassicSVfit svfitAlgo2;
  //  SVFitStorage svfitstorage;
  TH2D *ff_fracs_qcd_;
  TH2D *ff_fracs_wjets_;
  TH2D *ff_fracs_qcd_ss_;
  TH2D *ff_fracs_wjets_ss_;
  TH2D *ff_fracs_qcd_aiso_;
  TH2D *ff_fracs_wjets_aiso_;
  TH2D *ff_fracs_qcd_highmt_;
  TH2D *ff_fracs_wjets_highmt_;
  std::shared_ptr<RooWorkspace> ff_ws_;
  FakeFactors* _FF;
  
 private:
  // Selection Variables and Histos
  
  std::vector<TH1D> polarimetricAcopAngleEven;
  std::vector<TH1D> polarimetricAcopAngleOdd;
  std::vector<TH1D> polarimetricAcopAngleMM;
  std::vector<TH1D> AcopAngleEven;
  std::vector<TH1D> AcopAngleOdd;
  std::vector<TH1D> AcopAngleMM;
  std::vector<TH1D> genpolarimetricAcopAngleEven;
  std::vector<TH1D> genpolarimetricAcopAngleOdd;
  std::vector<TH1D> genpolarimetricAcopAngleMM;
  std::vector<TH1D> genAcopAngleEven;
  std::vector<TH1D> genAcopAngleOdd;
  std::vector<TH1D> genAcopAngleMM;
  std::vector<TH1D> pullPVx;
  std::vector<TH1D> pullPVy;
  std::vector<TH1D> pullPVz;
  std::vector<TH1D> pullTauSVx;
  std::vector<TH1D> pullTauSVy;
  std::vector<TH1D> pullTauSVz;
  std::vector<TH1D> pullTauE;
  std::vector<TH1D> pullTauPt;
  std::vector<TH1D> pullTauPhi;
  std::vector<TH1D> pullTauEta;
  std::vector<TH1D> pullMuonE;
  std::vector<TH1D> pullMuonPt;
  std::vector<TH1D> pullMuonPhi;
  std::vector<TH1D> pullMuonEta;  
};

#endif
