#ifndef SingleMuSkim_h
#define SingleMuSkim_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "Objects.h"
class SingleMuSkim : public Selection {

 public:
  SingleMuSkim(TString Name_, TString id_);
  virtual ~SingleMuSkim();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,PrimeVtx,ntaus,nmuons,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables and Histos

  //------------- Truth Gen ---------
  std::vector<TH1D> TauDecayMode;

  std::vector<TH1D> OSPairMass;
  std::vector<TH1D> SSPairMass;
  std::vector<TH1D> MissingTEnergy;
  std::vector<TH1D> isPairCandOS;
  std::vector<TH1D> TauPT;
  std::vector<TH1D> JetPT;
  

};
#endif
