#ifndef FakeFactors_h
#define FakeFactors_h

#include <TH2D.h>
#include <string>
#include <TMVA/Reader.h>
#include <RooFunctor.h>
#include <RooWorkspace.h>
#include <ImpactParameter.h>
#include "TLorentzVector.h"


class FakeFactors {
 public:
  FakeFactors(Int_t theYear, TH2D* ff_fracs_qcd, TH2D* ff_fracs_wjets, TH2D* ff_fracs_qcd_ss, TH2D* ff_fracs_wjets_ss, TH2D* ff_fracs_qcd_aiso, TH2D* ff_fracs_wjets_aiso, TH2D* ff_fracs_qcd_highmt, TH2D* ff_fracs_wjets_highmt, std::shared_ptr<RooWorkspace> ff_ws_);
  ~FakeFactors(){};
  void Initialize(TLorentzVector taup4, TLorentzVector mup4, TLorentzVector metp4, int tauDM, int Njets, double dijetMass, double muMETmt, double muIso, float ipsig, bool isOS, bool isIso);
  std::map<std::string, double> GetFakeFactors(TString sysType);

 private:
  std::shared_ptr<RooWorkspace> w_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;
  std::map<std::string, double> fake_factors_;

  TH2D *ff_fracs_qcd_;
  TH2D *ff_fracs_wjets_;
  TH2D *ff_fracs_qcd_ss_;
  TH2D *ff_fracs_wjets_ss_;
  TH2D *ff_fracs_qcd_aiso_;
  TH2D *ff_fracs_wjets_aiso_;
  TH2D *ff_fracs_qcd_highmt_;
  TH2D *ff_fracs_wjets_highmt_;

  TString xml_file, year_;
  TMVA::Reader *reader_;

  float met_, pt_1_, pt_2_, mva_dm_2_, mt_1_, m_vis_, pt_tt_, mjj_, n_jets_;

  std::vector<std::string> systs_mvadm_;
  std::vector<double> args_, args_qcd_, args_w_, args_ttbar_;
  //
};

#endif
