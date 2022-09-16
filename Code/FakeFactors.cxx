#include <FakeFactors.h>

FakeFactors::FakeFactors(Int_t theYear, TH2D* ff_fracs_qcd, TH2D* ff_fracs_wjets, TH2D* ff_fracs_qcd_ss, TH2D* ff_fracs_wjets_ss, TH2D* ff_fracs_qcd_aiso, TH2D* ff_fracs_wjets_aiso, TH2D* ff_fracs_qcd_highmt, TH2D* ff_fracs_wjets_highmt, std::shared_ptr<RooWorkspace> ff_ws_) {

  auto year_ = std::to_string(theYear);
  ff_fracs_qcd_ = ff_fracs_qcd;
  ff_fracs_wjets_ = ff_fracs_wjets;
  ff_fracs_qcd_ss_ = ff_fracs_qcd_ss;
  ff_fracs_wjets_ss_ = ff_fracs_wjets_ss;
  ff_fracs_qcd_aiso_ = ff_fracs_qcd_aiso;
  ff_fracs_wjets_aiso_ = ff_fracs_wjets_aiso;
  ff_fracs_qcd_highmt_ = ff_fracs_qcd_highmt;
  ff_fracs_wjets_highmt_ = ff_fracs_wjets_highmt;

  systs_mvadm_ = {"","_wjets_syst_up","_wjets_syst_down","_wjets_met_up","_wjets_met_down","_wjets_l_pt_up","_wjets_l_pt_down","_wjets_stat_unc1_njet0_mvadm10_up","_wjets_stat_unc2_njet0_mvadm10_up","_wjets_stat_unc1_njet0_mvadm10_down","_wjets_stat_unc2_njet0_mvadm10_down","_wjets_stat_unc1_njet1_mvadm10_up","_wjets_stat_unc2_njet1_mvadm10_up","_wjets_stat_unc1_njet1_mvadm10_down","_wjets_stat_unc2_njet1_mvadm10_down","_wjets_stat_unc1_njet2_mvadm10_up","_wjets_stat_unc2_njet2_mvadm10_up","_wjets_stat_unc1_njet2_mvadm10_down","_wjets_stat_unc2_njet2_mvadm10_down","_qcd_syst_up","_qcd_syst_down","_qcd_met_up","_qcd_met_down","_qcd_l_pt_up","_qcd_l_pt_down","_qcd_stat_unc1_njet0_mvadm10_up","_qcd_stat_unc2_njet0_mvadm10_up","_qcd_stat_unc1_njet0_mvadm10_down","_qcd_stat_unc2_njet0_mvadm10_down","_qcd_stat_unc1_njet1_mvadm10_up","_qcd_stat_unc2_njet1_mvadm10_up","_qcd_stat_unc1_njet1_mvadm10_down","_qcd_stat_unc2_njet1_mvadm10_down","_qcd_stat_unc1_njet2_mvadm10_up","_qcd_stat_unc2_njet2_mvadm10_up","_qcd_stat_unc1_njet2_mvadm10_down","_qcd_stat_unc2_njet2_mvadm10_down","_ttbar_syst_up","_ttbar_syst_down","_ttbar_met_up","_ttbar_met_down"};

  for(auto s : systs_mvadm_) {
    fns_["ff_lt_medium_mvadmbins"+s] = std::shared_ptr<RooFunctor>(
								   ff_ws_->function(("ff_mt_medium_mvadmbins"+s).c_str())->functor(ff_ws_->argSet("pt,mvadm,ipsig,njets,m_pt,os,met_var_qcd,met_var_w,mt,m_iso,pass_single,mvis,WpT,wjets_frac,qcd_frac,ttbar_frac")));
  }

  fns_["ff_lt_medium_mvadmbins_qcd"] = std::shared_ptr<RooFunctor>(
            ff_ws_->function("ff_mt_medium_mvadmbins_qcd")->functor(ff_ws_->argSet("pt,mvadm,ipsig,njets,m_pt,os,met_var_qcd,m_iso,pass_single")));
  fns_["ff_lt_medium_mvadmbins_wjets"] = std::shared_ptr<RooFunctor>(
            ff_ws_->function("ff_mt_medium_mvadmbins_wjets")->functor(ff_ws_->argSet("pt,mvadm,ipsig,njets,m_pt,met_var_w,mt,pass_single,mvis,WpT")));
}

void FakeFactors::Initialize(TLorentzVector taup4, TLorentzVector mup4, TLorentzVector metp4, int tauDM, int Njets, double dijetMass, double muMETmt, double muIso, float ipsig, bool isOS, bool isIso) {

  TLorentzVector metp4w(metp4.Px() + mup4.Px(), metp4.Py() + mup4.Py(), 0., mup4.Pt());
  //
  pt_tt_ = sqrt(pow((taup4.Px() + mup4.Px() + metp4.Px()),2) + pow((taup4.Py() + mup4.Py() + metp4.Py()),2));
  pt_1_ = mup4.Pt();
  pt_2_ = taup4.Pt();
  met_ = metp4.Pt();
  m_vis_ = (taup4 + mup4).M();
  n_jets_ = Njets;
  mjj_ = dijetMass;
  mva_dm_2_ = tauDM;
  mt_1_ = muMETmt;
  //
  double iso_1_ = muIso;
  double met_var_qcd = (metp4.Pt()/taup4.Pt())*cos(metp4.DeltaPhi(taup4));
  double met_var_w = (metp4w.Pt()/taup4.Pt())*cos(metp4w.DeltaPhi(taup4));
  double WpT = metp4w.Pt();
  //
  double singlemupt = 25.;
  if(year_ == "2016") singlemupt = 23.;
  double pass_single = 1.;
  if(mup4.Pt() < singlemupt) pass_single = 0.;
  //
  // load MVA scroes reader for fractions
  reader_ = new TMVA::Reader();
  reader_->AddVariable("pt_tt", &pt_tt_);
  reader_->AddVariable("pt_1", &pt_1_);
  reader_->AddVariable("pt_2", &pt_2_);
  reader_->AddVariable("met", &met_);
  reader_->AddVariable("m_vis", &m_vis_);
  reader_->AddVariable("n_jets", &n_jets_);
  reader_->AddVariable("mjj", &mjj_);
  reader_->AddVariable("mva_dm_2", &mva_dm_2_);
  reader_->AddVariable("mt_1", &mt_1_);
  xml_file="/opt/sbg/cms/safe1/cms/msessini/IPHCAnalysisTools/Code/CommonFiles/FakeFactors/fractions_2018_mt.xml";
  reader_->BookMVA("BDT method", xml_file);
  //
  std::vector<float> scores = reader_->EvaluateMulticlass("BDT method");
  double qcd_score = scores[1];
  double w_score = scores[0];
  //
  double w_frac = ff_fracs_wjets_->GetBinContent(ff_fracs_wjets_->FindBin(qcd_score,w_score));
  double qcd_frac = ff_fracs_qcd_->GetBinContent(ff_fracs_qcd_->FindBin(qcd_score,w_score));
  //
  if(!isOS) {
    w_frac = ff_fracs_wjets_ss_->GetBinContent(ff_fracs_wjets_ss_->FindBin(qcd_score,w_score));
    qcd_frac = ff_fracs_qcd_ss_->GetBinContent(ff_fracs_qcd_ss_->FindBin(qcd_score,w_score));
    if(w_frac==0. && qcd_frac==0.) qcd_frac = 1.;
  }
  if(!isIso) {
    w_frac = ff_fracs_wjets_aiso_->GetBinContent(ff_fracs_wjets_aiso_->FindBin(qcd_score,w_score));
    qcd_frac = ff_fracs_qcd_aiso_->GetBinContent(ff_fracs_qcd_aiso_->FindBin(qcd_score,w_score));
    if(w_frac==0. && qcd_frac==0.) qcd_frac = 1.;
  }
  if(mt_1_>70) {
    w_frac = ff_fracs_wjets_highmt_->GetBinContent(ff_fracs_wjets_highmt_->FindBin(qcd_score,w_score));
    qcd_frac = ff_fracs_qcd_highmt_->GetBinContent(ff_fracs_qcd_highmt_->FindBin(qcd_score,w_score));
    if(w_frac==0. && qcd_frac==0.) w_frac = 1.;
  }
  double ttbar_frac = 1. - w_frac - qcd_frac;
  double os = 1.;
  if(!isOS) os = 0.;
  //
  args_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,pt_1_,os,met_var_qcd,met_var_w,mt_1_,iso_1_,pass_single,m_vis_,WpT,w_frac,qcd_frac,ttbar_frac};
  args_qcd_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,pt_1_,os,met_var_qcd,iso_1_,pass_single};
  args_w_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,pt_1_,met_var_w,mt_1_,pass_single,m_vis_,WpT};
  args_ttbar_ = {pt_2_,mva_dm_2_,ipsig,n_jets_,met_var_w};
}

std::map<std::string, double> FakeFactors::GetFakeFactors(TString sysType) {

  fake_factors_["ff_nominal"] = fns_["ff_lt_medium_mvadmbins"]->eval(args_.data());
  fake_factors_["ff_nominal_qcd"] = fns_["ff_lt_medium_mvadmbins_qcd"]->eval(args_qcd_.data()); 
  fake_factors_["ff_nominal_w"] = fns_["ff_lt_medium_mvadmbins_wjets"]->eval(args_w_.data());
  //
  if(sysType == "Nominal") {
    for(auto s : systs_mvadm_) {
      if(s == "") continue;
      fake_factors_["ff"+s] = fns_["ff_lt_medium_mvadmbins"+s]->eval(args_.data());
    }
  }
  return fake_factors_;
}

//////////////////////////////////////////
