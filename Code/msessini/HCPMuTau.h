#ifndef HCPMuTau_h
#define HCPMuTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "boost/functional/hash.hpp"
//#include "SVFitStorage.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
//#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "ReferenceScaleFactors.h"
#include "ScaleFactor.h"
#include "Objects.h"
//#include "PUReweight.h"
#include "tauTrigSFreader.h"
#include "DataMCCorrections.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "RecoilCorrector.h"
#include "MEtSys.h"
#include "BDTClassification.h"
#include "FakeFactors.h"

#include "RooWorkspace.h"
#include "RooFunctor.h"
#include <memory>

class HCPMuTau : public Selection {

 public:
  HCPMuTau(TString Name_, TString id_, char* Channel_, char* CPstate_);
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
  ReferenceScaleFactors *RSF;
  char* Channel;
  char* CPstate;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  //PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  //cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  
  //DataMCCorrections DataMC_Corr;
  //tauTrigSFreader tauTrgSF;
  //Int_t year;
  BDTClassification *BDT;
  
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
  
  std::vector<TH1D> polarimetricAcopAngle;
  
};

#endif
