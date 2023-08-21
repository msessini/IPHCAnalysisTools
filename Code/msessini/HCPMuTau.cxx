#include "HCPMuTau.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SkimConfig.h"
#include "Logger.h"


#include "TVector3.h"
#include "TMath.h"
//#include "Objects.h"
#include <algorithm>

HCPMuTau::HCPMuTau(TString Name_, TString id_, char* Channel_):
  Selection(Name_,id_)
  //DataMC_Corr(true,true,false),
  //tauTrgSF("tight")
{
  Channel = Channel_;
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;

  TFile *f_fracs=TFile::Open(((std::string)std::getenv("workdir")+"Code/CommonFiles/FakeFactors/mva_fract_mt_2018.root").c_str(), "READ");
  ff_fracs_qcd_ = (TH2D*)f_fracs->Get("QCD");
  ff_fracs_wjets_ = (TH2D*)f_fracs->Get("W");
  ff_fracs_qcd_->SetDirectory(0);
  ff_fracs_wjets_->SetDirectory(0);
  f_fracs->Close();

  TFile *f_fracs_ss=TFile::Open(((std::string)std::getenv("workdir")+"Code/CommonFiles/FakeFactors/mva_fract_mt_2018_ss.root").c_str(), "READ");
  ff_fracs_qcd_ss_ = (TH2D*)f_fracs_ss->Get("QCD");
  ff_fracs_wjets_ss_ = (TH2D*)f_fracs_ss->Get("W");
  ff_fracs_qcd_ss_->SetDirectory(0);
  ff_fracs_wjets_ss_->SetDirectory(0);
  f_fracs_ss->Close();

  TFile *f_fracs_aiso=TFile::Open(((std::string)std::getenv("workdir")+"Code/CommonFiles/FakeFactors/mva_fract_mt_2018_aiso.root").c_str(), "READ");
  ff_fracs_qcd_aiso_ = (TH2D*)f_fracs_aiso->Get("QCD");
  ff_fracs_wjets_aiso_ = (TH2D*)f_fracs_aiso->Get("W");
  ff_fracs_qcd_aiso_->SetDirectory(0);
  ff_fracs_wjets_aiso_->SetDirectory(0);
  f_fracs_aiso->Close();

  TFile *f_fracs_highmt=TFile::Open(((std::string)std::getenv("workdir")+"Code/CommonFiles/FakeFactors/mva_fract_mt_2018_highmt.root").c_str(), "READ");
  ff_fracs_qcd_highmt_ = (TH2D*)f_fracs_highmt->Get("QCD");
  ff_fracs_wjets_highmt_ = (TH2D*)f_fracs_highmt->Get("W");
  ff_fracs_qcd_highmt_->SetDirectory(0);
  ff_fracs_wjets_highmt_->SetDirectory(0);
  f_fracs_highmt->Close();

  TFile *f = new TFile();
  if(theYear == 2016) f=TFile::Open(((std::string)std::getenv("workdir")+"Code/CommonFiles/FakeFactors/fakefactors_ws_mt_lite_2016.root").c_str(), "READ");
  if(theYear == 2017) f=TFile::Open(((std::string)std::getenv("workdir")+"Code/CommonFiles/FakeFactors/fakefactors_ws_mt_lite_2017.root").c_str(), "READ");
  if(theYear == 2018) f=TFile::Open(((std::string)std::getenv("workdir")+"Code/CommonFiles/FakeFactors/fakefactors_ws_mt_lite_2018.root").c_str(), "READ");
  ff_ws_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  f->Close();

  _FF = new FakeFactors(theYear, ff_fracs_qcd_, ff_fracs_wjets_, ff_fracs_qcd_ss_, ff_fracs_wjets_ss_, ff_fracs_qcd_aiso_, ff_fracs_wjets_aiso_, ff_fracs_qcd_highmt_, ff_fracs_wjets_highmt_, ff_ws_);

  BDT=new BDTClassification();
  BDT->PreAnalysis();

}

HCPMuTau::~HCPMuTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  HCPMuTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    //if(i==Trigger)             cut.at(Trigger)=1;
    //if(i==Id_and_Kin)            cut.at(Id_and_Kin)=1;
    //if(i==NPairsFound)         cut.at(NPairsFound)=1;
    //if(i==GoodIndex)           cut.at(GoodIndex)=1.;
    //if(i==ZTTMC)                 cut.at(ZTTMC)=1.;
    //if(i==METFilters)            cut.at(METFilters)=1.;
    //if(i==genmatch)              cut.at(genmatch)=1;
    if(i==TausIsolation)         cut.at(TausIsolation)=1;
    //if(i==AgainstEleMu)          cut.at(AgainstEleMu)=1;
    //if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    //if(i==LeptonVeto)            cut.at(LeptonVeto)=0;
    //if(i==PairCharge)            cut.at(PairCharge)=1.;
    //if(i==PairMass)              cut.at(PairMass)=40.;
    //if(i==MTM)                 cut.at(MTM)=40;

  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // if(i==PrimeVtx){
    //   title.at(i)="Number of Prime Vertices $(N>$";
    //   title.at(i)+=cut.at(PrimeVtx);
    //   title.at(i)+=")";
    //   htitle=title.at(i);
    //   htitle.ReplaceAll("$","");
    //   htitle.ReplaceAll("\\","#");
    //   hlabel="Number of Prime Vertices";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    // }
    //if(i==Trigger){
    //  title.at(i)="Trigger Matching";
    //  hlabel="At least 1 good pair with Trig+Matching";
    //  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //}
    // else if(i==Id_and_Kin){
    //   title.at(i)="Id and Kinematic";
    //   hlabel="Number of Event with good particles";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==NPairsFound){
    //   title.at(i)="Pairs with good DeltaR";
    //   hlabel="Pairs with good DeltaR";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    //if(i==GoodIndex){
    //title.at(i)="Valid Index";
    //hlabel="Valid Index";
    //Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GoodIndex_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GoodIndex_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // if(i==ZTTMC){
    //   title.at(i)="ZTT MC";
    //   hlabel="ZTT MC";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZTTMC_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZTTMC_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==METFilters){
    //   title.at(i)="MET Filters";
    //   hlabel="MET Filters";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_METFilters_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_METFilters_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    //else if(i==Id_and_Kin){
    //  title.at(i)="Id and Kinematic";
    //  hlabel="Number of Event with good particles";
    //  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //}
    // else if(i==genmatch){
    //   title.at(i)="genmatch";
    //   hlabel="genmatch";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    if(i==TausIsolation){
      title.at(i)="Taus Isolation";
      hlabel="Isolation of Taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    //else if(i==AgainstEleMu){
    //  title.at(i)="Against Electrons and Muons";
    //  hlabel="Against Electrons and Muons";
    //  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //}
    // else if(i==Tau2Isolation){
    //   title.at(i)="Tau2 Isolation";
    //   hlabel="Isolation of Tau2";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    //else if(i==LeptonVeto){
    //  title.at(i)="Third Lepton Veto";
    //  hlabel="Third Lepton Veto";
    //  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //}
    //else if(i==PairCharge){
    //  title.at(i)="Pair Charge";
    //  hlabel="is pair OS";
    //  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //}
    //else if(i==PairMass){
    //  title.at(i)="Pair Visible Mass";
    //  hlabel="M(tau-tau)";
    //  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,30,0,150,hlabel,"Events"));
    //  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,30,0,150,hlabel,"Events"));
    //}

  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  
  //Control plots
  MuonpT=HConfig.GetTH1D(Name+"_MuonpT","Transverse momentum of selected muon",30,20,80,"p_{T}(#mu) (GeV)","Events");
  TaupT=HConfig.GetTH1D(Name+"_TaupT","Visible transverse momentum of selected tau",30,25,85,"p_{T}(#tau) (GeV)","Events");
  DitaupT=HConfig.GetTH1D(Name+"_DitaupT","Visible transverse momentum of selected tau pair",30,0,300,"p_{T}(#tau#tau) (GeV)","Events");
  Njets=HConfig.GetTH1D(Name+"_Njets","Number N of jets",6,0,6,"N_{jets}","Events");
  LeadingJetpT=HConfig.GetTH1D(Name+"_LeadingJetpT","Transverse momentum of leading jet",30,30,150,"p_{T}(lead.j) (GeV)","Events");
  SubleadingJetpT=HConfig.GetTH1D(Name+"_SubleadingJetpT","Transverse momentum of subleading jet",30,30,150,"p_{T}(sublead.j) (GeV)","Events");
  DijetpT=HConfig.GetTH1D(Name+"_DijetpT","Transverse momentum of the two leading jets",50,0,250,"p_{T}(jj) (GeV)","Events");
  DijetMass=HConfig.GetTH1D(Name+"_DijetMass","Invariant mass of the two leading jets",50,0,250,"m(jj) (GeV)","Events");
  DijetDeltaEta=HConfig.GetTH1D(Name+"_DijetDeltaEta","Eta separation between the two leading jets",12,0,6,"#Delta#eta_{jj}","Events");
  VisibleMass=HConfig.GetTH1D(Name+"_VisibleMass","Visible mass of selected tau pair",30,0,300,"m_{vis} (GeV)","Events");
  FastMTTditauMass=HConfig.GetTH1D(Name+"_fastMTTditauMass","Invariant mass of reconstructed tau pair (FastMTT)",35,0,350,"m_{#tau#tau} (GeV)","Events");
  PUPPImet=HConfig.GetTH1D(Name+"_PUPPImet","Missing transverse energy",30,0,150,"PUPPI-MET (GeV)","Events");
  MuMETmt=HConfig.GetTH1D(Name+"_MuMETmt","Transverse mass of #mu + MET system",30,0,150,"mt(#mu+MET) (GEV)","Events");
  BDTscoreHiggs=HConfig.GetTH1D(Name+"_BDTscoreHiggs","BDTscoreHiggs",20,0,1,"score","Events");
  BDTscoreJetFakes=HConfig.GetTH1D(Name+"_BDTscoreJetFakes","BDTscoreJetFakes",20,0,1,"score","Events");
  BDTscoreZTT=HConfig.GetTH1D(Name+"_BDTscoreZTT","BDTscoreZTT",20,0,1,"score","Events");
  BDTscoreA1MUHiggs=HConfig.GetTH1D(Name+"_BDTscoreA1MUHiggs","BDTscoreA1MUHiggs",20,0,1,"score","Events");
  BDTscoreA1MUJetFakes=HConfig.GetTH1D(Name+"_BDTscoreA1MUJetFakes","BDTscoreA1MUJetFakes",20,0,1,"score","Events");
  BDTscoreA1MUZTT=HConfig.GetTH1D(Name+"_BDTscoreA1MUZTT","BDTscoreA1MUZTT",20,0,1,"score","Events");
  PhiCPEvenPV=HConfig.GetTH1D(Name+"_PhiCPEvenPV","acop",10,0,2*TMath::Pi(),"#phi_{CP} (rad)","Events");
  PhiCPOddPV=HConfig.GetTH1D(Name+"_PhiCPOddPV","acop",10,0,2*TMath::Pi(),"#phi_{CP} (rad)","Events");
  PhiCPMMPV=HConfig.GetTH1D(Name+"_PhiCPMMPV","acop",10,0,2*TMath::Pi(),"#phi_{CP} (rad)","Events");
  PhiCPEvenDP=HConfig.GetTH1D(Name+"_PhiCPEvenDP","acop",10,0,2*TMath::Pi(),"#phi_{CP} (rad)","Events");
  PhiCPOddDP=HConfig.GetTH1D(Name+"_PhiCPOddDP","acop",10,0,2*TMath::Pi(),"#phi_{CP} (rad)","Events");
  PhiCPMMDP=HConfig.GetTH1D(Name+"_PhiCPMMDP","acop",10,0,2*TMath::Pi(),"#phi_{CP} (rad)","Events");
  ResMuonPt=HConfig.GetTH1D(Name+"_ResMuonPt","res",50,-0.2,0.2,"pull","Events");
  ResMuonEta=HConfig.GetTH1D(Name+"_ResMuonEta","res",50,-0.05,0.05,"pull","Events");
  ResTauPt=HConfig.GetTH1D(Name+"_ResTauPt","res",50,-1,1,"pull","Events");
  ResTauEta=HConfig.GetTH1D(Name+"_ResTauEta","res",50,-0.2,0.2,"pull","Events");
  ResPVx=HConfig.GetTH1D(Name+"_ResPVx","res",50,-0.2,0.2,"pull","Events");
  ResPVy=HConfig.GetTH1D(Name+"_ResPVy","res",50,-0.2,0.2,"pull","Events");
  ResPVz=HConfig.GetTH1D(Name+"_ResPVz","res",50,-0.2,0.2,"pull","Events");
  ResSVx=HConfig.GetTH1D(Name+"_ResSVx","res",50,-1,1,"pull","Events");
  ResSVy=HConfig.GetTH1D(Name+"_ResSVy","res",50,-1,1,"pull","Events");
  ResSVz=HConfig.GetTH1D(Name+"_ResSVz","res",50,-1,1,"pull","Events");
  //
  ttbar_contamination_hgs=HConfig.GetTH2D(Name+"_ttbar_contamination_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ttbar_contamination_ztt=HConfig.GetTH2D(Name+"_ttbar_contamination_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ttbar_contamination_fkj=HConfig.GetTH2D(Name+"_ttbar_contamination_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //PhiCP even
  phiCPeven_nominal_hgs=HConfig.GetTH2D(Name+"_phiCPeven_nominal_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_nominal_ztt=HConfig.GetTH2D(Name+"_phiCPeven_nominal_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_nominal_fkj=HConfig.GetTH2D(Name+"_phiCPeven_nominal_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_CMS_scale_t_3prong_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_t_3prong_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //
  phiCPeven_CMS_scale_t_3prong_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_t_3prong_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_CMS_scale_t_3prong_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_t_3prong_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //PhiCP odd
  phiCPodd_nominal_hgs=HConfig.GetTH2D(Name+"_phiCPodd_nominal_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_nominal_ztt=HConfig.GetTH2D(Name+"_phiCPodd_nominal_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_nominal_fkj=HConfig.GetTH2D(Name+"_phiCPodd_nominal_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //
  phiCPodd_CMS_scale_t_3prong_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_t_3prong_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_CMS_scale_t_3prong_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_t_3prong_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_CMS_scale_t_3prong_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_t_3prong_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //PhiCP MM
  phiCPMM_nominal_hgs=HConfig.GetTH2D(Name+"_phiCPMM_nominal_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Up_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Up_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Down_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Down_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_nominal_ztt=HConfig.GetTH2D(Name+"_phiCPMM_nominal_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Up_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Up_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Down_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Down_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_nominal_fkj=HConfig.GetTH2D(Name+"_phiCPMM_nominal_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Up_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Up_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Down_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Down_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //
  phiCPMM_CMS_scale_t_3prong_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_t_3prong_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgs=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgs","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_CMS_scale_t_3prong_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_t_3prong_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVUp_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVUp_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVDown_ztt=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVDown_ztt","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_CMS_scale_t_3prong_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_t_3prong_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkj=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkj","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //Application Region histograms for QCD substraction
  MuonpTAR=HConfig.GetTH1D(Name+"_MuonpTAR","Transverse momentum of selected muon",30,20,80,"p_{T}(#mu) (GeV)","Events");
  TaupTAR=HConfig.GetTH1D(Name+"_TaupTAR","Visible transverse momentum of selected tau",30,25,85,"p_{T}(#tau) (GeV)","Events");
  DitaupTAR=HConfig.GetTH1D(Name+"_DitaupTAR","Visible transverse momentum of selected tau pair",30,0,300,"p_{T}(#tau#tau) (GeV)","Events");
  NjetsAR=HConfig.GetTH1D(Name+"_NjetsAR","Number N of jets",6,0,6,"N_{jets}","Events");
  LeadingJetpTAR=HConfig.GetTH1D(Name+"_LeadingJetpTAR","Transverse momentum of leading jet",30,30,150,"p_{T}(lead.j) (GeV)","Events");
  SubleadingJetpTAR=HConfig.GetTH1D(Name+"_SubleadingJetpTAR","Transverse momentum of subleading jet",30,30,150,"p_{T}(sublead.j) (GeV)","Events");
  DijetpTAR=HConfig.GetTH1D(Name+"_DijetpTAR","Transverse momentum of the two leading jets",50,0,250,"p_{T}(jj) (GeV)","Events");
  DijetMassAR=HConfig.GetTH1D(Name+"_DijetMassAR","Invariant mass of the two leading jets",50,0,250,"m(jj) (GeV)","Events");
  DijetDeltaEtaAR=HConfig.GetTH1D(Name+"_DijetDeltaEtaAR","Eta separation between the two leading jets",12,0,6,"#Delta#eta_{jj}","Events");
  VisibleMassAR=HConfig.GetTH1D(Name+"_VisibleMassAR","Visible mass of selected tau pair",30,0,300,"m_{vis} (GeV)","Events");
  FastMTTditauMassAR=HConfig.GetTH1D(Name+"_fastMTTditauMassAR","Invariant mass of reconstructed tau pair (FastMTT)",35,0,350,"m_{#tau#tau} (GeV)","Events");
  PUPPImetAR=HConfig.GetTH1D(Name+"_PUPPImetAR","Missing transverse energy",30,0,150,"PUPPI-MET (GeV)","Events");
  MuMETmtAR=HConfig.GetTH1D(Name+"_MuMETmtAR","Transverse mass of #mu + MET system",30,0,150,"mt(#mu+MET) (GEV)","Events");
  BDTscoreHiggsAR=HConfig.GetTH1D(Name+"_BDTscoreHiggsAR","BDTscoreHiggsAR",20,0,1,"score","Events");
  BDTscoreJetFakesAR=HConfig.GetTH1D(Name+"_BDTscoreJetFakesAR","BDTscoreJetFakesAR",20,0,1,"score","Events");
  BDTscoreZTTAR=HConfig.GetTH1D(Name+"_BDTscoreZTTAR","BDTscoreZTTAR",20,0,1,"score","Events");
  BDTscoreA1MUHiggsAR=HConfig.GetTH1D(Name+"_BDTscoreA1MUHiggsAR","BDTscoreA1MUHiggsAR",20,0,1,"score","Events");
  BDTscoreA1MUJetFakesAR=HConfig.GetTH1D(Name+"_BDTscoreA1MUJetFakesAR","BDTscoreA1MUJetFakesAR",20,0,1,"score","Events");
  BDTscoreA1MUZTTAR=HConfig.GetTH1D(Name+"_BDTscoreA1MUZTTAR","BDTscoreA1MUZTTAR",20,0,1,"score","Events");
  //
  //PhiCP even
  phiCPeven_nominal_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_nominal_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_nominal_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_nominal_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_nominal_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_nominal_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_eff_b_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_eff_b_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_dyShape_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_dyShape_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_gg_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_gg_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_qcd_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_qcd_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_wjets_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_wjets_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_ff_mt_sub_systDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_ff_mt_sub_systDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_ttbar_embeded_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_ttbar_embeded_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_PreFire_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_PreFire_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //
  phiCPeven_CMS_scale_t_3prong_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_t_3prong_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_CMS_scale_t_3prong_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_t_3prong_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPeven_CMS_scale_t_3prong_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_t_3prong_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_t_3prong_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_mu_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_mu_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_res_j_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_res_j_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //PhiCP odd
  phiCPodd_nominal_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_nominal_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_nominal_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_nominal_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_nominal_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_nominal_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_eff_b_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_eff_b_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_dyShape_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_dyShape_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_gg_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_gg_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_qcd_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_qcd_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_wjets_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_wjets_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_ff_mt_sub_systDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_ff_mt_sub_systDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_ttbar_embeded_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_ttbar_embeded_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_PreFire_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_PreFire_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //
  phiCPodd_CMS_scale_t_3prong_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_t_3prong_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_CMS_scale_t_3prong_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_t_3prong_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPodd_CMS_scale_t_3prong_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_t_3prong_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_t_3prong_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_mu_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_mu_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_res_j_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_res_j_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //PhiCP MM
  phiCPMM_nominal_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_nominal_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Up_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Up_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Down_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Down_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_nominal_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_nominal_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Up_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Up_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Down_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Down_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_nominal_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_nominal_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_eff_b_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_eff_b_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_dyShape_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_dyShape_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_gg_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_gg_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_qcd_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_qcd_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets0Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets0Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets1Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets1Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Up_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Up_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_wjets_syst_njets2Down_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_wjets_syst_njets2Down_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_ff_mt_sub_systDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_ff_mt_sub_systDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_ttbar_embeded_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_ttbar_embeded_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_PreFire_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_PreFire_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  //
  phiCPMM_CMS_scale_t_3prong_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_t_3prong_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgsAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgsAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_CMS_scale_t_3prong_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_t_3prong_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVUp_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVUp_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVDown_zttAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVDown_zttAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  ////
  phiCPMM_CMS_scale_t_3prong_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_t_3prong_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_t_3prong_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_mu_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_mu_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_res_j_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_res_j_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");
  phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkjAR=HConfig.GetTH2D(Name+"_phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkjAR","#phi_{CP}",60,0,2*TMath::Pi(),20,0,1,"#phi_{CP} (rad)","BDT score");

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove

}

void  HCPMuTau::Store_ExtraDist(){

  Extradist1d.push_back(&MuonpT);
  Extradist1d.push_back(&TaupT);
  Extradist1d.push_back(&DitaupT);
  Extradist1d.push_back(&Njets);
  Extradist1d.push_back(&LeadingJetpT);
  Extradist1d.push_back(&SubleadingJetpT);
  Extradist1d.push_back(&DijetpT);
  Extradist1d.push_back(&DijetMass);
  Extradist1d.push_back(&DijetDeltaEta);
  Extradist1d.push_back(&VisibleMass);
  Extradist1d.push_back(&FastMTTditauMass);
  Extradist1d.push_back(&PUPPImet);
  Extradist1d.push_back(&MuMETmt);
  Extradist1d.push_back(&PhiCPEvenPV);
  Extradist1d.push_back(&PhiCPOddPV);
  Extradist1d.push_back(&PhiCPMMPV);
  Extradist1d.push_back(&PhiCPEvenDP);
  Extradist1d.push_back(&PhiCPOddDP);
  Extradist1d.push_back(&PhiCPMMDP);
  Extradist1d.push_back(&ResMuonPt);
  Extradist1d.push_back(&ResMuonEta);
  Extradist1d.push_back(&ResTauPt);
  Extradist1d.push_back(&ResTauEta);
  Extradist1d.push_back(&ResPVx);
  Extradist1d.push_back(&ResPVy);
  Extradist1d.push_back(&ResPVz);
  Extradist1d.push_back(&ResSVx);
  Extradist1d.push_back(&ResSVy);
  Extradist1d.push_back(&ResSVz);
  Extradist1d.push_back(&BDTscoreHiggs);
  Extradist1d.push_back(&BDTscoreJetFakes);
  Extradist1d.push_back(&BDTscoreZTT);
  Extradist1d.push_back(&BDTscoreA1MUHiggs);
  Extradist1d.push_back(&BDTscoreA1MUJetFakes);
  Extradist1d.push_back(&BDTscoreA1MUZTT);
  //
  if(Selection::Get_SysType() == "default") {
    Extradist2d.push_back(&phiCPeven_nominal_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_eff_b_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_eff_b_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_dyShape_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_dyShape_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_gg_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_gg_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_sub_systUp_hgs);
    Extradist2d.push_back(&phiCPeven_ff_mt_sub_systDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_ttbar_embeded_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_ttbar_embeded_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_PreFire_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_PreFire_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPeven_nominal_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_eff_b_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_eff_b_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_dyShape_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_dyShape_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_ttbarShape_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_ttbarShape_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_gg_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_gg_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_sub_systUp_ztt);
    Extradist2d.push_back(&phiCPeven_ff_mt_sub_systDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_ttbar_embeded_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_ttbar_embeded_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_PreFire_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_PreFire_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPeven_nominal_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_eff_b_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_eff_b_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_htt_dyShape_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_htt_dyShape_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_gg_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_gg_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_qcd_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_wjets_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_sub_systUp_fkj);
    Extradist2d.push_back(&phiCPeven_ff_mt_sub_systDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_ttbar_embeded_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_ttbar_embeded_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_PreFire_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_PreFire_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "TESUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_t_3prong_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_t_3prong_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_t_3prong_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "TESDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_t_3prong_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_t_3prong_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_t_3prong_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "MESUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_mu_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_mu_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_mu_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "MESDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_mu_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_mu_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_mu_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "FlavorQCDUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "FlavorQCDDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "RelativeBalUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "RelativeBalDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "HFUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "HFDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "HF_YEARUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "HF_YEARDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1Up") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1Down") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1_YEARUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1_YEARDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "EC2Up") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "EC2Down") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "EC2_YEARUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "EC2_YEARDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "AbsoluteUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "AbsoluteDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "Absolute_YEARUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "Absolute_YEARDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "JERUp") {
    Extradist2d.push_back(&phiCPeven_CMS_res_j_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_res_j_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_res_j_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "JERDown") {
    Extradist2d.push_back(&phiCPeven_CMS_res_j_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_res_j_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_res_j_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METResoUp") {
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METResoDown") {
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METScaleUp") {
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METScaleDown") {
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METUnclusteredUp") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METUnclusteredDown") {
    Extradist2d.push_back(&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj);
  }
  //
  if(Selection::Get_SysType() == "default") {
    Extradist2d.push_back(&phiCPodd_nominal_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_eff_b_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_eff_b_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_dyShape_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_dyShape_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_gg_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_gg_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_sub_systUp_hgs);
    Extradist2d.push_back(&phiCPodd_ff_mt_sub_systDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_ttbar_embeded_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_ttbar_embeded_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_PreFire_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_PreFire_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPodd_nominal_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_eff_b_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_eff_b_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_dyShape_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_dyShape_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_ttbarShape_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_ttbarShape_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_gg_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_gg_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_sub_systUp_ztt);
    Extradist2d.push_back(&phiCPodd_ff_mt_sub_systDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_ttbar_embeded_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_ttbar_embeded_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_PreFire_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_PreFire_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPodd_nominal_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_eff_b_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_eff_b_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_htt_dyShape_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_htt_dyShape_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_gg_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_gg_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_qcd_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_wjets_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_sub_systUp_fkj);
    Extradist2d.push_back(&phiCPodd_ff_mt_sub_systDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_ttbar_embeded_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_ttbar_embeded_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_PreFire_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_PreFire_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "TESUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_t_3prong_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_t_3prong_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_t_3prong_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "TESDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_t_3prong_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_t_3prong_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_t_3prong_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "MESUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_mu_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_mu_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_mu_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "MESDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_mu_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_mu_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_mu_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "FlavorQCDUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "FlavorQCDDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "RelativeBalUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "RelativeBalDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "HFUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "HFDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "HF_YEARUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "HF_YEARDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1Up") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1Down") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1_YEARUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1_YEARDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "EC2Up") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "EC2Down") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "EC2_YEARUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "EC2_YEARDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "AbsoluteUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "AbsoluteDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "Absolute_YEARUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "Absolute_YEARDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "JERUp") {
    Extradist2d.push_back(&phiCPodd_CMS_res_j_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_res_j_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_res_j_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "JERDown") {
    Extradist2d.push_back(&phiCPodd_CMS_res_j_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_res_j_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_res_j_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METResoUp") {
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METResoDown") {
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METScaleUp") {
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METScaleDown") {
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METUnclusteredUp") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_met_unclustered_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METUnclusteredDown") {
    Extradist2d.push_back(&phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPodd_CMS_scale_met_unclustered_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkj);
  }
  //
  if(Selection::Get_SysType() == "default") {
    Extradist2d.push_back(&phiCPMM_nominal_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_eff_b_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_eff_b_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_dyShape_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_dyShape_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_gg_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_gg_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets0Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets0Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets1Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets1Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets2Up_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets2Down_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_sub_systUp_hgs);
    Extradist2d.push_back(&phiCPMM_ff_mt_sub_systDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_ttbar_embeded_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_ttbar_embeded_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_PreFire_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_PreFire_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPMM_nominal_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_eff_b_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_eff_b_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_dyShape_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_dyShape_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_ttbarShape_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_ttbarShape_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_gg_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_gg_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets0Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets0Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets1Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets1Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets2Up_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets2Down_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_sub_systUp_ztt);
    Extradist2d.push_back(&phiCPMM_ff_mt_sub_systDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_ttbar_embeded_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_ttbar_embeded_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_PreFire_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_PreFire_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPMM_nominal_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_eff_b_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_eff_b_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_htt_dyShape_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_htt_dyShape_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_gg_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_gg_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_qcd_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets0Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets0Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets1Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets1Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets2Up_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_wjets_syst_njets2Down_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_sub_systUp_fkj);
    Extradist2d.push_back(&phiCPMM_ff_mt_sub_systDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_ttbar_embeded_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_ttbar_embeded_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_PreFire_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_PreFire_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "TESUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_t_3prong_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_t_3prong_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_t_3prong_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "TESDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_t_3prong_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_t_3prong_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_t_3prong_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "MESUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_mu_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_mu_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_mu_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "MESDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_mu_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_mu_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_mu_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "FlavorQCDUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "FlavorQCDDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "RelativeBalUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "RelativeBalDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "HFUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "HFDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "HF_YEARUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "HF_YEARDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1Up") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1Down") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1_YEARUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "BBEC1_YEARDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "EC2Up") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "EC2Down") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "EC2_YEARUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "EC2_YEARDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "AbsoluteUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "AbsoluteDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "Absolute_YEARUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "Absolute_YEARDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt);
    //
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj);
    Extradist2d.push_back(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "JERUp") {
    Extradist2d.push_back(&phiCPMM_CMS_res_j_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_res_j_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_res_j_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "JERDown") {
    Extradist2d.push_back(&phiCPMM_CMS_res_j_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_res_j_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_res_j_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METResoUp") {
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METResoDown") {
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METScaleUp") {
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METScaleDown") {
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkj);
  }
  if(Selection::Get_SysType() == "METUnclusteredUp") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_met_unclustered_13TeVUp_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkj);
  }
  if(Selection::Get_SysType() == "METUnclusteredDown") {
    Extradist2d.push_back(&phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgs);
    Extradist2d.push_back(&phiCPMM_CMS_scale_met_unclustered_13TeVDown_ztt);
    Extradist2d.push_back(&phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkj);
  }
}

void  HCPMuTau::doEvent()  { //  Method called on every event

  if(Ntp->pvx() == -99 && Ntp->pvy() == -99 && Ntp->pvz() == -99) return; 
  //
  unsigned int t;
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  //
  TLorentzVector taup4, mup4, metp4;
  taup4.SetPtEtaPhiE(Ntp->tauPt(), Ntp->tauEta(), Ntp->tauPhi(), Ntp->tauE());
  mup4.SetPtEtaPhiE(Ntp->muPt(), Ntp->muEta(), Ntp->muPhi(), Ntp->muE());
  metp4.SetPtEtaPhiE(Ntp->PUPPImet(), 0., Ntp->PUPPImetphi(), 0.);
  //
  TLorentzVector genmup4(Ntp->genMuonpx(), Ntp->genMuonpy(), Ntp->genMuonpz(), Ntp->genMuonE());
  TLorentzVector gentaup4(Ntp->genTaupx(), Ntp->genTaupy(), Ntp->genTaupz(), Ntp->genTauE());
  //
  //Compute FakeFactors here
  _FF->Initialize(taup4, mup4, metp4, Ntp->tauDM(), Ntp->Njets(), Ntp->dijetMass(), Ntp->muMETmt(), Ntp->muIso(), Ntp->tauIPsignificance(), Ntp->isOSpair(), Ntp->isIso());
  std::map<std::string, double> FFmap = _FF->GetFakeFactors(Selection::Get_SysType());
  //
  bool isZEWK = (id==203);
  bool isWEWK = (id==201 || id==202);
  bool isSignalRegion = Ntp->isOSpair() && Ntp->isIso() && Ntp->isTightJetID() && Ntp->Nbtags() == 0 && Ntp->isMediumID();
  bool isApplicationRegion = Ntp->isOSpair() && Ntp->isIso() && Ntp->isTightJetID() && Ntp->Nbtags() == 0 && !Ntp->isMediumID();
  bool isA1MU = (Ntp->tauDM() == 10 && Ntp->muIPsignificance() >= 1.5);
  //
  int NJets = Ntp->Njets();
  //
  bool genMatched = false;
  if(Ntp->isSignal()) genMatched = (Ntp->tauGenMatch() == 5 && Ntp->muGenMatch() == 4); //Real taus from signal
  else if(Ntp->isZ() || Ntp->isVV() || Ntp->isTTbar() || Ntp->isSingleTop() || isZEWK) genMatched = (!(Ntp->tauGenMatch() == 5 && Ntp->muGenMatch() == 4) && !(Ntp->tauGenMatch() == 6 || Ntp->muGenMatch() == 6)); //All fake taus except the ones from fake jets
  else if(Ntp->isW()) genMatched = false; //Already accounted
  else if(isWEWK) genMatched = !((Ntp->tauGenMatch() == 6 && Ntp->muGenMatch() == 2) || (Ntp->tauGenMatch() == 5 && Ntp->muGenMatch() == 6)); //You don't want real W+jets
  else if(id == 35) genMatched = (Ntp->tauGenMatch() == 5 && Ntp->muGenMatch() == 4); //Real taus from embedded samples
  else genMatched = true;
  //
  double w          = Ntp->wTot();
  double wTrg       = Ntp->wTrg();
  double wIDvsJet   = Ntp->wIDvsJet();
  double wIDvsEle   = Ntp->wIDvsEle();
  double wIDvsMu    = Ntp->wIDvsMu();
  double wIDMu      = Ntp->wIDMu();
  double wTrkMu     = Ntp->wTrkMu();
  double wMC        = Ntp->wMC();
  double wPU        = Ntp->wPU();
  double wPrefiring = Ntp->wPrefiring();
  double wBtag      = Ntp->wBtag();
  double wZpT       = Ntp->wZpT();
  double wToppT     = Ntp->wToppT();
  double wSignal    = Ntp->wSignal();
  //
  if(theYear == 2016) { //some bugs...
    if(id == 24) wMC = 1.;
    if(Ntp->isSignal()) w *= (1./wSignal);
    if(id == 35 && w>1000) w = 0.;
  }
  //
  double wTrgUp =1.;
  double wTrgDown = 1.;
  double wZpTUp = 1.;
  double wZpTDown = 1.;
  double wToppTUp = 1.;
  double wToppTDown = 1.;
  double wBtagUp = 1.;
  double wBtagDown = 1.;
  double wScaleUp = 1.;
  double wScaleDown = 1.;
  double wPSISRUp = 1.;
  double wPSISRDown = 1.;
  double wPSFSRUp = 1.;
  double wPSFSRDown = 1.;
  double wPrefiringUp = 1.;
  double wPrefiringDown = 1.;

  if(Selection::Get_SysType() == "default") {
    wTrgUp = Ntp->wTrgUp();
    wTrgDown = Ntp->wTrgDown();
    wZpTUp = Ntp->wZpTUp();
    wZpTDown = Ntp->wZpTDown();
    wToppTUp = Ntp->wToppTUp(); 
    wToppTDown = Ntp->wToppTDown();
    wBtagUp = Ntp->wBtagUp();
    wBtagDown = Ntp->wBtagDown();
    wScaleUp = Ntp->wScaleUp();
    wScaleDown = Ntp->wScaleDown();
    wPSISRUp = Ntp->wPSISRUp();
    wPSISRDown = Ntp->wPSISRDown();
    wPSFSRUp = Ntp->wPSFSRUp();
    wPSFSRDown = Ntp->wPSFSRDown();
    wPrefiringUp = Ntp->wPrefiringUp();
    wPrefiringDown = Ntp->wPrefiringDown();
  }

  //Stitching
  if(theYear == 2016) {
    if(Ntp->isZ()){
      if(isDY10to50 && id==30) w*=19.0080307;
      else if (isDY10to50 && id!=30) w*=-9999;
      if(Ntp->Npartons()==0||Ntp->Npartons() >= 5) w*=1.49005321266736;
      if(Ntp->Npartons()==1) w*=0.47521301562902;
      if(Ntp->Npartons()==2) w*=0.492313539105283;
      if(Ntp->Npartons()==3) w*=0.504730998787436;
      if(Ntp->Npartons()==4) w*=0.414018612660014;
    }
    if(Ntp->isW()){
      if(Ntp->Npartons()==0||Ntp->Npartons() >= 5) w*=25.3889395535057;
      if(Ntp->Npartons()==1) w*=6.82217925376281;
      if(Ntp->Npartons()==2) w*=2.09118705989359;
      if(Ntp->Npartons()==3) w*=0.686191177275517;
      if(Ntp->Npartons()==4) w*=0.691068406413718;
    }
  }
  if(theYear == 2017) {
    if(Ntp->isZ()){
      if(isDY10to50 && id==30) w*=19.5191962215717;
      else if (isDY10to50 && id!=30) w*=-9999;
      if(Ntp->Npartons()==0||Ntp->Npartons() >= 5) w*=2.5973962888025;
      if(Ntp->Npartons()==1) w*=0.453927720730712;
      if(Ntp->Npartons()==2) w*=0.922892174009934;
      if(Ntp->Npartons()==3) w*=0.5938168058852;
      if(Ntp->Npartons()==4) w*=0.408277745068798;
    }
    if(Ntp->isW()){
      if(Ntp->Npartons()==0||Ntp->Npartons() >= 5) w*=4.96803789167669;
      if(Ntp->Npartons()==1) w*=4.96803789167669;
      if(Ntp->Npartons()==2) w*=14.8987957676995;
      if(Ntp->Npartons()==3) w*=2.32379673645138;
      if(Ntp->Npartons()==4) w*=2.14597542644866;
    }
  }
  if(theYear == 2018) {
    if(Ntp->isZ()){
      if(isDY10to50 && id==30) w*=28.2040833505999;
      else if (isDY10to50 && id!=30) w*=-9999;
      if(Ntp->Npartons()==0||Ntp->Npartons() >= 5) w*=3.62105387778545;
      if(Ntp->Npartons()==1) w*=0.69924705836082;
      if(Ntp->Npartons()==2) w*=0.791046490589383;
      if(Ntp->Npartons()==3) w*=0.978992685464496;
      if(Ntp->Npartons()==4) w*=0.821659745356239;
    }
    if(Ntp->isW()){
      if(Ntp->Npartons()==0||Ntp->Npartons() >= 5) w*=51.714857425559;
      if(Ntp->Npartons()==1) w*=10.8715284468658;
      if(Ntp->Npartons()==2) w*=8.26944354355374;
      if(Ntp->Npartons()==3) w*=4.39072366992094;
      if(Ntp->Npartons()==4) w*=3.28821099006665;
    }
  }
  if(id==35) {
    if(theYear == 2018) w*=1.0415;
    if(theYear == 2017) w*=1.0396;
  }
  //
  double wEven = w;
  double wOdd = w;
  double wMM = w;
  //
  double FF = FFmap["ff_nominal"];
  double FFqcd1jet0Up = FFmap["ff_qcd_stat_unc1_njet0_mvadm10_up"];
  double FFqcd1jet0Down = FFmap["ff_qcd_stat_unc1_njet0_mvadm10_down"];
  double FFqcd1jet1Up = FFmap["ff_qcd_stat_unc1_njet1_mvadm10_up"];
  double FFqcd1jet1Down = FFmap["ff_qcd_stat_unc1_njet1_mvadm10_down"];
  double FFqcd1jet2Up = FFmap["ff_qcd_stat_unc1_njet2_mvadm10_up"];
  double FFqcd1jet2Down = FFmap["ff_qcd_stat_unc1_njet2_mvadm10_down"];
  double FFqcd2jet0Up = FFmap["ff_qcd_stat_unc2_njet0_mvadm10_up"];
  double FFqcd2jet0Down = FFmap["ff_qcd_stat_unc2_njet0_mvadm10_down"];
  double FFqcd2jet1Up = FFmap["ff_qcd_stat_unc2_njet1_mvadm10_up"];
  double FFqcd2jet1Down = FFmap["ff_qcd_stat_unc2_njet1_mvadm10_down"];
  double FFqcd2jet2Up = FFmap["ff_qcd_stat_unc2_njet2_mvadm10_up"];
  double FFqcd2jet2Down = FFmap["ff_qcd_stat_unc2_njet2_mvadm10_down"];
  double FFwjets1jet0Up = FFmap["ff_wjets_stat_unc1_njet0_mvadm10_up"];
  double FFwjets1jet0Down = FFmap["ff_wjets_stat_unc1_njet0_mvadm10_down"];
  double FFwjets1jet1Up = FFmap["ff_wjets_stat_unc1_njet1_mvadm10_up"];
  double FFwjets1jet1Down = FFmap["ff_wjets_stat_unc1_njet1_mvadm10_down"];
  double FFwjets1jet2Up = FFmap["ff_wjets_stat_unc1_njet2_mvadm10_up"];
  double FFwjets1jet2Down = FFmap["ff_wjets_stat_unc1_njet2_mvadm10_down"];
  double FFwjets2jet0Up = FFmap["ff_wjets_stat_unc2_njet0_mvadm10_up"];
  double FFwjets2jet0Down = FFmap["ff_wjets_stat_unc2_njet0_mvadm10_down"];
  double FFwjets2jet1Up = FFmap["ff_wjets_stat_unc2_njet1_mvadm10_up"];
  double FFwjets2jet1Down = FFmap["ff_wjets_stat_unc2_njet1_mvadm10_down"];
  double FFwjets2jet2Up = FFmap["ff_wjets_stat_unc2_njet2_mvadm10_up"];
  double FFwjets2jet2Down = FFmap["ff_wjets_stat_unc2_njet2_mvadm10_down"];
  double FFqcdmetUp = FFmap["ff_qcd_met_up"];
  double FFqcdmetDown = FFmap["ff_qcd_met_down"];
  double FFqcdlptUp = FFmap["ff_qcd_l_pt_up"];
  double FFqcdlptDown = FFmap["ff_qcd_l_pt_down"];
  double FFqcdsystUp = FFmap["ff_qcd_syst_up"];
  double FFqcdsystDown = FFmap["ff_qcd_syst_down"];
  double FFwjetsmetUp = FFmap["ff_wjets_met_up"];
  double FFwjetsmetDown = FFmap["ff_wjets_met_down"];
  double FFwjetslptUp = FFmap["ff_wjets_l_pt_up"];
  double FFwjetslptDown = FFmap["ff_wjets_l_pt_down"];
  double FFwjetssystUp = FFmap["ff_wjets_syst_up"];
  double FFwjetssystDown = FFmap["ff_wjets_syst_down"];
  //
  std::pair<float, int> max_pair;
  std::vector<float> scores = {};
  BDT->Execute(Ntp->muPt(),Ntp->tauPt(),Ntp->ditauPt(),Ntp->Njets(),Ntp->leadingjetPt(),Ntp->subleadingjetPt(),Ntp->dijetPt(),Ntp->dijetMass(),Ntp->dijetdeltaEta(),Ntp->pairvisMass(),Ntp->fastMTTmass(),Ntp->muMETmt(),Ntp->PUPPImet(),Ntp->EventNumber(),scores,max_pair);
  //
  double Angle = Ntp->dpPhiCP();
  //
  value.at(TausIsolation) = isSignalRegion;
  pass.at(TausIsolation) = value.at(TausIsolation);
  if(isApplicationRegion) {
    if(Ntp->isData()) {
      MuonpT.at(1).Fill(Ntp->muPt(),FF);
      TaupT.at(1).Fill(Ntp->tauPt(),FF);
      DitaupT.at(1).Fill(Ntp->ditauPt(),FF);
      Njets.at(1).Fill(Ntp->Njets(),FF);
      LeadingJetpT.at(1).Fill(Ntp->leadingjetPt(),FF);
      SubleadingJetpT.at(1).Fill(Ntp->subleadingjetPt(),FF);
      DijetpT.at(1).Fill(Ntp->dijetPt(),FF);
      DijetMass.at(1).Fill(Ntp->dijetMass(),FF);
      DijetDeltaEta.at(1).Fill(Ntp->dijetdeltaEta(),FF);
      VisibleMass.at(1).Fill(Ntp->pairvisMass(),FF);
      FastMTTditauMass.at(1).Fill(Ntp->fastMTTmass(),FF);
      PUPPImet.at(1).Fill(Ntp->PUPPImet(),FF);
      MuMETmt.at(1).Fill(Ntp->muMETmt(),FF);
      if(max_pair.second == 0) {
	BDTscoreHiggs.at(1).Fill(max_pair.first,FF);
	if(isA1MU) BDTscoreA1MUHiggs.at(1).Fill(max_pair.first,FF);
      }
      else if(max_pair.second == 1) {
	BDTscoreZTT.at(1).Fill(max_pair.first,FF);
	if(isA1MU) BDTscoreA1MUZTT.at(1).Fill(max_pair.first,FF);
      }
      else if(max_pair.second == 2) {
	BDTscoreJetFakes.at(1).Fill(max_pair.first,FF);
	if(isA1MU) BDTscoreA1MUJetFakes.at(1).Fill(max_pair.first,FF);
      }
      //
      if(isA1MU) {
	if(Selection::Get_SysType() == "default") {
	  Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_nominal_hgs,&phiCPeven_nominal_ztt,&phiCPeven_nominal_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_eff_b_13TeVUp_hgs,&phiCPeven_CMS_eff_b_13TeVUp_ztt,&phiCPeven_CMS_eff_b_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_eff_b_13TeVDown_hgs,&phiCPeven_CMS_eff_b_13TeVDown_ztt,&phiCPeven_CMS_eff_b_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_dyShape_13TeVUp_hgs,&phiCPeven_CMS_htt_dyShape_13TeVUp_ztt,&phiCPeven_CMS_htt_dyShape_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_dyShape_13TeVDown_hgs,&phiCPeven_CMS_htt_dyShape_13TeVDown_ztt,&phiCPeven_CMS_htt_dyShape_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgs,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_ztt,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgs,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_ztt,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_gg_13TeVUp_hgs,&phiCPeven_CMS_scale_gg_13TeVUp_ztt,&phiCPeven_CMS_scale_gg_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_gg_13TeVDown_hgs,&phiCPeven_CMS_scale_gg_13TeVDown_ztt,&phiCPeven_CMS_scale_gg_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgs,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_ztt,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgs,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_ztt,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgs,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_ztt,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgs,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_ztt,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_ff_mt_sub_systUp_hgs,&phiCPeven_ff_mt_sub_systUp_ztt,&phiCPeven_ff_mt_sub_systUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_ff_mt_sub_systDown_hgs,&phiCPeven_ff_mt_sub_systDown_ztt,&phiCPeven_ff_mt_sub_systDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_ttbar_embeded_13TeVUp_hgs,&phiCPeven_CMS_ttbar_embeded_13TeVUp_ztt,&phiCPeven_CMS_ttbar_embeded_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_ttbar_embeded_13TeVDown_hgs,&phiCPeven_CMS_ttbar_embeded_13TeVDown_ztt,&phiCPeven_CMS_ttbar_embeded_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_PreFire_13TeVUp_hgs,&phiCPeven_CMS_PreFire_13TeVUp_ztt,&phiCPeven_CMS_PreFire_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_PreFire_13TeVDown_hgs,&phiCPeven_CMS_PreFire_13TeVDown_ztt,&phiCPeven_CMS_PreFire_13TeVDown_fkj,true);
	  if(NJets == 0) {
	    Ntp->FillHist(t,Angle,max_pair,FFqcd1jet0Up,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcd1jet0Down,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjets1jet0Up,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjets1jet0Down,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcd2jet0Up,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcd2jet0Down,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjets2jet0Up,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjets2jet0Down,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcdmetUp,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcdmetDown,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcdlptUp,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcdlptDown,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcdsystUp,&phiCPeven_ff_mt_qcd_syst_njets0Up_hgs,&phiCPeven_ff_mt_qcd_syst_njets0Up_ztt,&phiCPeven_ff_mt_qcd_syst_njets0Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFqcdsystDown,&phiCPeven_ff_mt_qcd_syst_njets0Down_hgs,&phiCPeven_ff_mt_qcd_syst_njets0Down_ztt,&phiCPeven_ff_mt_qcd_syst_njets0Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPeven_ff_mt_wjets_syst_njets0Up_hgs,&phiCPeven_ff_mt_wjets_syst_njets0Up_ztt,&phiCPeven_ff_mt_wjets_syst_njets0Up_fkj,true);
	    Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPeven_ff_mt_wjets_syst_njets0Down_hgs,&phiCPeven_ff_mt_wjets_syst_njets0Down_ztt,&phiCPeven_ff_mt_wjets_syst_njets0Down_fkj,true);
	  }
	  if(NJets == 1) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet1Up,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet1Down,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet1Up,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet1Down,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet1Up,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet1Down,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet1Up,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet1Down,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetUp,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetDown,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptUp,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptDown,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystUp,&phiCPeven_ff_mt_qcd_syst_njets1Up_hgs,&phiCPeven_ff_mt_qcd_syst_njets1Up_ztt,&phiCPeven_ff_mt_qcd_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystDown,&phiCPeven_ff_mt_qcd_syst_njets1Down_hgs,&phiCPeven_ff_mt_qcd_syst_njets1Down_ztt,&phiCPeven_ff_mt_qcd_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPeven_ff_mt_wjets_syst_njets1Up_hgs,&phiCPeven_ff_mt_wjets_syst_njets1Up_ztt,&phiCPeven_ff_mt_wjets_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPeven_ff_mt_wjets_syst_njets1Down_hgs,&phiCPeven_ff_mt_wjets_syst_njets1Down_ztt,&phiCPeven_ff_mt_wjets_syst_njets1Down_fkj,true);
	  }
	  if(NJets == 2) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet2Up,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet2Down,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet2Up,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet2Down,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet2Up,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet2Down,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet2Up,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet2Down,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPeven_ff_mt_wjets_syst_njets2Up_hgs,&phiCPeven_ff_mt_wjets_syst_njets2Up_ztt,&phiCPeven_ff_mt_wjets_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPeven_ff_mt_wjets_syst_njets2Down_hgs,&phiCPeven_ff_mt_wjets_syst_njets2Down_ztt,&phiCPeven_ff_mt_wjets_syst_njets2Down_fkj,true);
 	  }
	}   
        if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_t_3prong_13TeVUp_hgs,&phiCPeven_CMS_scale_t_3prong_13TeVUp_ztt,&phiCPeven_CMS_scale_t_3prong_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_t_3prong_13TeVDown_hgs,&phiCPeven_CMS_scale_t_3prong_13TeVDown_ztt,&phiCPeven_CMS_scale_t_3prong_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_mu_13TeVUp_hgs,&phiCPeven_CMS_scale_mu_13TeVUp_ztt,&phiCPeven_CMS_scale_mu_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_mu_13TeVDown_hgs,&phiCPeven_CMS_scale_mu_13TeVDown_ztt,&phiCPeven_CMS_scale_mu_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgs,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_ztt,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgs,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_ztt,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkj,true);
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj,true);
        }
	if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_res_j_13TeVUp_hgs,&phiCPeven_CMS_res_j_13TeVUp_ztt,&phiCPeven_CMS_res_j_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_res_j_13TeVDown_hgs,&phiCPeven_CMS_res_j_13TeVDown_ztt,&phiCPeven_CMS_res_j_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgs,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_ztt,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgs,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_ztt,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgs,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_ztt,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgs,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_ztt,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj,true);
	//
	if(Selection::Get_SysType() == "default") {
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_nominal_hgs,&phiCPodd_nominal_ztt,&phiCPodd_nominal_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_eff_b_13TeVUp_hgs,&phiCPodd_CMS_eff_b_13TeVUp_ztt,&phiCPodd_CMS_eff_b_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_eff_b_13TeVDown_hgs,&phiCPodd_CMS_eff_b_13TeVDown_ztt,&phiCPodd_CMS_eff_b_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_dyShape_13TeVUp_hgs,&phiCPodd_CMS_htt_dyShape_13TeVUp_ztt,&phiCPodd_CMS_htt_dyShape_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_dyShape_13TeVDown_hgs,&phiCPodd_CMS_htt_dyShape_13TeVDown_ztt,&phiCPodd_CMS_htt_dyShape_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgs,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_ztt,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgs,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_ztt,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_gg_13TeVUp_hgs,&phiCPodd_CMS_scale_gg_13TeVUp_ztt,&phiCPodd_CMS_scale_gg_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_gg_13TeVDown_hgs,&phiCPodd_CMS_scale_gg_13TeVDown_ztt,&phiCPodd_CMS_scale_gg_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgs,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_ztt,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgs,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_ztt,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgs,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_ztt,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgs,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_ztt,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_ff_mt_sub_systUp_hgs,&phiCPodd_ff_mt_sub_systUp_ztt,&phiCPodd_ff_mt_sub_systUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_ff_mt_sub_systDown_hgs,&phiCPodd_ff_mt_sub_systDown_ztt,&phiCPodd_ff_mt_sub_systDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_ttbar_embeded_13TeVUp_hgs,&phiCPodd_CMS_ttbar_embeded_13TeVUp_ztt,&phiCPodd_CMS_ttbar_embeded_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_ttbar_embeded_13TeVDown_hgs,&phiCPodd_CMS_ttbar_embeded_13TeVDown_ztt,&phiCPodd_CMS_ttbar_embeded_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_PreFire_13TeVUp_hgs,&phiCPodd_CMS_PreFire_13TeVUp_ztt,&phiCPodd_CMS_PreFire_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_PreFire_13TeVDown_hgs,&phiCPodd_CMS_PreFire_13TeVDown_ztt,&phiCPodd_CMS_PreFire_13TeVDown_fkj,true);
          if(NJets == 0) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet0Up,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet0Down,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet0Up,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet0Down,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet0Up,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet0Down,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet0Up,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet0Down,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetUp,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetDown,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptUp,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptDown,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystUp,&phiCPodd_ff_mt_qcd_syst_njets0Up_hgs,&phiCPodd_ff_mt_qcd_syst_njets0Up_ztt,&phiCPodd_ff_mt_qcd_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystDown,&phiCPodd_ff_mt_qcd_syst_njets0Down_hgs,&phiCPodd_ff_mt_qcd_syst_njets0Down_ztt,&phiCPodd_ff_mt_qcd_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPodd_ff_mt_wjets_syst_njets0Up_hgs,&phiCPodd_ff_mt_wjets_syst_njets0Up_ztt,&phiCPodd_ff_mt_wjets_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPodd_ff_mt_wjets_syst_njets0Down_hgs,&phiCPodd_ff_mt_wjets_syst_njets0Down_ztt,&phiCPodd_ff_mt_wjets_syst_njets0Down_fkj,true);
          }
          if(NJets == 1) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet1Up,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet1Down,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet1Up,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet1Down,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet1Up,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet1Down,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet1Up,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet1Down,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetUp,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetDown,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptUp,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptDown,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystUp,&phiCPodd_ff_mt_qcd_syst_njets1Up_hgs,&phiCPodd_ff_mt_qcd_syst_njets1Up_ztt,&phiCPodd_ff_mt_qcd_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystDown,&phiCPodd_ff_mt_qcd_syst_njets1Down_hgs,&phiCPodd_ff_mt_qcd_syst_njets1Down_ztt,&phiCPodd_ff_mt_qcd_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPodd_ff_mt_wjets_syst_njets1Up_hgs,&phiCPodd_ff_mt_wjets_syst_njets1Up_ztt,&phiCPodd_ff_mt_wjets_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPodd_ff_mt_wjets_syst_njets1Down_hgs,&phiCPodd_ff_mt_wjets_syst_njets1Down_ztt,&phiCPodd_ff_mt_wjets_syst_njets1Down_fkj,true);
          }
          if(NJets == 2) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet2Up,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet2Down,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet2Up,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet2Down,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet2Up,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet2Down,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet2Up,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet2Down,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPodd_ff_mt_wjets_syst_njets2Up_hgs,&phiCPodd_ff_mt_wjets_syst_njets2Up_ztt,&phiCPodd_ff_mt_wjets_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPodd_ff_mt_wjets_syst_njets2Down_hgs,&phiCPodd_ff_mt_wjets_syst_njets2Down_ztt,&phiCPodd_ff_mt_wjets_syst_njets2Down_fkj,true);
          }
	}
        if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_t_3prong_13TeVUp_hgs,&phiCPodd_CMS_scale_t_3prong_13TeVUp_ztt,&phiCPodd_CMS_scale_t_3prong_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_t_3prong_13TeVDown_hgs,&phiCPodd_CMS_scale_t_3prong_13TeVDown_ztt,&phiCPodd_CMS_scale_t_3prong_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_mu_13TeVUp_hgs,&phiCPodd_CMS_scale_mu_13TeVUp_ztt,&phiCPodd_CMS_scale_mu_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_mu_13TeVDown_hgs,&phiCPodd_CMS_scale_mu_13TeVDown_ztt,&phiCPodd_CMS_scale_mu_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgs,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_ztt,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgs,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_ztt,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkj,true);
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_res_j_13TeVUp_hgs,&phiCPodd_CMS_res_j_13TeVUp_ztt,&phiCPodd_CMS_res_j_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_res_j_13TeVDown_hgs,&phiCPodd_CMS_res_j_13TeVDown_ztt,&phiCPodd_CMS_res_j_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgs,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_ztt,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgs,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_ztt,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgs,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_ztt,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgs,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_ztt,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgs,&phiCPodd_CMS_scale_met_unclustered_13TeVUp_ztt,&phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgs,&phiCPodd_CMS_scale_met_unclustered_13TeVDown_ztt,&phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkj,true);
	//
	if(Selection::Get_SysType() == "default") {
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_nominal_hgs,&phiCPMM_nominal_ztt,&phiCPMM_nominal_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_eff_b_13TeVUp_hgs,&phiCPMM_CMS_eff_b_13TeVUp_ztt,&phiCPMM_CMS_eff_b_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_eff_b_13TeVDown_hgs,&phiCPMM_CMS_eff_b_13TeVDown_ztt,&phiCPMM_CMS_eff_b_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_dyShape_13TeVUp_hgs,&phiCPMM_CMS_htt_dyShape_13TeVUp_ztt,&phiCPMM_CMS_htt_dyShape_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_dyShape_13TeVDown_hgs,&phiCPMM_CMS_htt_dyShape_13TeVDown_ztt,&phiCPMM_CMS_htt_dyShape_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgs,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_ztt,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgs,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_ztt,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_gg_13TeVUp_hgs,&phiCPMM_CMS_scale_gg_13TeVUp_ztt,&phiCPMM_CMS_scale_gg_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_gg_13TeVDown_hgs,&phiCPMM_CMS_scale_gg_13TeVDown_ztt,&phiCPMM_CMS_scale_gg_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgs,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_ztt,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgs,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_ztt,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgs,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_ztt,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgs,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_ztt,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_ff_mt_sub_systUp_hgs,&phiCPMM_ff_mt_sub_systUp_ztt,&phiCPMM_ff_mt_sub_systUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_ff_mt_sub_systDown_hgs,&phiCPMM_ff_mt_sub_systDown_ztt,&phiCPMM_ff_mt_sub_systDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_ttbar_embeded_13TeVUp_hgs,&phiCPMM_CMS_ttbar_embeded_13TeVUp_ztt,&phiCPMM_CMS_ttbar_embeded_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_ttbar_embeded_13TeVDown_hgs,&phiCPMM_CMS_ttbar_embeded_13TeVDown_ztt,&phiCPMM_CMS_ttbar_embeded_13TeVDown_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_PreFire_13TeVUp_hgs,&phiCPMM_CMS_PreFire_13TeVUp_ztt,&phiCPMM_CMS_PreFire_13TeVUp_fkj,true);
          Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_PreFire_13TeVDown_hgs,&phiCPMM_CMS_PreFire_13TeVDown_ztt,&phiCPMM_CMS_PreFire_13TeVDown_fkj,true);
          if(NJets == 0) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet0Up,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet0Down,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet0Up,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet0Down,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet0Up,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet0Down,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet0Up,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet0Down,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetUp,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetDown,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptUp,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptDown,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystUp,&phiCPMM_ff_mt_qcd_syst_njets0Up_hgs,&phiCPMM_ff_mt_qcd_syst_njets0Up_ztt,&phiCPMM_ff_mt_qcd_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystDown,&phiCPMM_ff_mt_qcd_syst_njets0Down_hgs,&phiCPMM_ff_mt_qcd_syst_njets0Down_ztt,&phiCPMM_ff_mt_qcd_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPMM_ff_mt_wjets_syst_njets0Up_hgs,&phiCPMM_ff_mt_wjets_syst_njets0Up_ztt,&phiCPMM_ff_mt_wjets_syst_njets0Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPMM_ff_mt_wjets_syst_njets0Down_hgs,&phiCPMM_ff_mt_wjets_syst_njets0Down_ztt,&phiCPMM_ff_mt_wjets_syst_njets0Down_fkj,true);
          }
          if(NJets == 1) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet1Up,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet1Down,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet1Up,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet1Down,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet1Up,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet1Down,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet1Up,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet1Down,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetUp,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdmetDown,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptUp,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdlptDown,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystUp,&phiCPMM_ff_mt_qcd_syst_njets1Up_hgs,&phiCPMM_ff_mt_qcd_syst_njets1Up_ztt,&phiCPMM_ff_mt_qcd_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcdsystDown,&phiCPMM_ff_mt_qcd_syst_njets1Down_hgs,&phiCPMM_ff_mt_qcd_syst_njets1Down_ztt,&phiCPMM_ff_mt_qcd_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPMM_ff_mt_wjets_syst_njets1Up_hgs,&phiCPMM_ff_mt_wjets_syst_njets1Up_ztt,&phiCPMM_ff_mt_wjets_syst_njets1Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPMM_ff_mt_wjets_syst_njets1Down_hgs,&phiCPMM_ff_mt_wjets_syst_njets1Down_ztt,&phiCPMM_ff_mt_wjets_syst_njets1Down_fkj,true);
          }
          if(NJets == 2) {
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet2Up,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd1jet2Down,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet2Up,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets1jet2Down,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet2Up,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFqcd2jet2Down,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet2Up,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjets2jet2Down,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetUp,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetsmetDown,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptUp,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetslptDown,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystUp,&phiCPMM_ff_mt_wjets_syst_njets2Up_hgs,&phiCPMM_ff_mt_wjets_syst_njets2Up_ztt,&phiCPMM_ff_mt_wjets_syst_njets2Up_fkj,true);
            Ntp->FillHist(t,Angle,max_pair,FFwjetssystDown,&phiCPMM_ff_mt_wjets_syst_njets2Down_hgs,&phiCPMM_ff_mt_wjets_syst_njets2Down_ztt,&phiCPMM_ff_mt_wjets_syst_njets2Down_fkj,true);
          }
	}
        if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_t_3prong_13TeVUp_hgs,&phiCPMM_CMS_scale_t_3prong_13TeVUp_ztt,&phiCPMM_CMS_scale_t_3prong_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_t_3prong_13TeVDown_hgs,&phiCPMM_CMS_scale_t_3prong_13TeVDown_ztt,&phiCPMM_CMS_scale_t_3prong_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_mu_13TeVUp_hgs,&phiCPMM_CMS_scale_mu_13TeVUp_ztt,&phiCPMM_CMS_scale_mu_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_mu_13TeVDown_hgs,&phiCPMM_CMS_scale_mu_13TeVDown_ztt,&phiCPMM_CMS_scale_mu_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgs,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_ztt,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgs,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_ztt,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkj,true);
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj,true);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj,true);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj,true);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj,true);
        }
        if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_res_j_13TeVUp_hgs,&phiCPMM_CMS_res_j_13TeVUp_ztt,&phiCPMM_CMS_res_j_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_res_j_13TeVDown_hgs,&phiCPMM_CMS_res_j_13TeVDown_ztt,&phiCPMM_CMS_res_j_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgs,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_ztt,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgs,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_ztt,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgs,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_ztt,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgs,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_ztt,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkj,true);
        if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgs,&phiCPMM_CMS_scale_met_unclustered_13TeVUp_ztt,&phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkj,true);
        if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,FF,&phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgs,&phiCPMM_CMS_scale_met_unclustered_13TeVDown_ztt,&phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkj,true);
      }
    }
    if(!Ntp->isData()) {
      MuonpTAR.at(t).Fill(Ntp->muPt(),w*FF);
      TaupTAR.at(t).Fill(Ntp->tauPt(),w*FF);
      DitaupTAR.at(t).Fill(Ntp->ditauPt(),w*FF);
      NjetsAR.at(t).Fill(Ntp->Njets(),w*FF);
      LeadingJetpTAR.at(t).Fill(Ntp->leadingjetPt(),w*FF);
      SubleadingJetpTAR.at(t).Fill(Ntp->subleadingjetPt(),w*FF);
      DijetpTAR.at(t).Fill(Ntp->dijetPt(),w*FF);
      DijetMassAR.at(t).Fill(Ntp->dijetMass(),w*FF);
      DijetDeltaEtaAR.at(t).Fill(Ntp->dijetdeltaEta(),w*FF);
      VisibleMassAR.at(t).Fill(Ntp->pairvisMass(),w*FF);
      FastMTTditauMassAR.at(t).Fill(Ntp->fastMTTmass(),w*FF);
      PUPPImetAR.at(t).Fill(Ntp->PUPPImet(),w*FF);
      MuMETmtAR.at(t).Fill(Ntp->muMETmt(),w*FF);
      if(max_pair.second == 0) {
        BDTscoreHiggsAR.at(t).Fill(max_pair.first,w*FF);
        if(isA1MU) BDTscoreA1MUHiggsAR.at(t).Fill(max_pair.first,w*FF);
      }
      else if(max_pair.second == 1) {
        BDTscoreZTTAR.at(t).Fill(max_pair.first,w*FF);
        if(isA1MU) BDTscoreA1MUZTTAR.at(t).Fill(max_pair.first,w*FF);
      }
      else if(max_pair.second == 2) {
        BDTscoreJetFakesAR.at(t).Fill(max_pair.first,w*FF);
        if(isA1MU) BDTscoreA1MUJetFakesAR.at(t).Fill(max_pair.first,w*FF);
      }
      //
      if(isA1MU) {
        if(Ntp->isSignal()) {
          wOdd *= Ntp->wOdd();
          wEven *= Ntp->wEven();
          wMM *= Ntp->wMM();
        }
	if(Selection::Get_SysType() == "default") {
          Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_nominal_hgsAR,&phiCPeven_nominal_zttAR,&phiCPeven_nominal_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wTrgUp)/wTrg)*FF,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wTrgDown)/wTrg)*FF,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wBtagUp)/wBtag)*FF,&phiCPeven_CMS_eff_b_13TeVUp_hgsAR,&phiCPeven_CMS_eff_b_13TeVUp_zttAR,&phiCPeven_CMS_eff_b_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wBtagDown)/wBtag)*FF,&phiCPeven_CMS_eff_b_13TeVDown_hgsAR,&phiCPeven_CMS_eff_b_13TeVDown_zttAR,&phiCPeven_CMS_eff_b_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wZpTUp)/wZpT)*FF,&phiCPeven_CMS_htt_dyShape_13TeVUp_hgsAR,&phiCPeven_CMS_htt_dyShape_13TeVUp_zttAR,&phiCPeven_CMS_htt_dyShape_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wZpTDown)/wZpT)*FF,&phiCPeven_CMS_htt_dyShape_13TeVDown_hgsAR,&phiCPeven_CMS_htt_dyShape_13TeVDown_zttAR,&phiCPeven_CMS_htt_dyShape_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wToppTUp)/wToppT)*FF,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgsAR,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_zttAR,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wToppTDown)/wToppT)*FF,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgsAR,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_zttAR,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*wScaleUp*FF,&phiCPeven_CMS_scale_gg_13TeVUp_hgsAR,&phiCPeven_CMS_scale_gg_13TeVUp_zttAR,&phiCPeven_CMS_scale_gg_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*wScaleDown*FF,&phiCPeven_CMS_scale_gg_13TeVDown_hgsAR,&phiCPeven_CMS_scale_gg_13TeVDown_zttAR,&phiCPeven_CMS_scale_gg_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*wPSISRUp*FF,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgsAR,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_zttAR,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*wPSISRDown*FF,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgsAR,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_zttAR,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*wPSFSRUp*FF,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgsAR,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_zttAR,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*wPSFSRDown*FF,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgsAR,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_zttAR,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_ff_mt_sub_systUp_hgsAR,&phiCPeven_ff_mt_sub_systUp_zttAR,&phiCPeven_ff_mt_sub_systUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_ff_mt_sub_systDown_hgsAR,&phiCPeven_ff_mt_sub_systDown_zttAR,&phiCPeven_ff_mt_sub_systDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_ttbar_embeded_13TeVUp_hgsAR,&phiCPeven_CMS_ttbar_embeded_13TeVUp_zttAR,&phiCPeven_CMS_ttbar_embeded_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_ttbar_embeded_13TeVDown_hgsAR,&phiCPeven_CMS_ttbar_embeded_13TeVDown_zttAR,&phiCPeven_CMS_ttbar_embeded_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wPrefiringUp)/wPrefiring)*FF,&phiCPeven_CMS_PreFire_13TeVUp_hgsAR,&phiCPeven_CMS_PreFire_13TeVUp_zttAR,&phiCPeven_CMS_PreFire_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wEven*wPrefiringDown)/wPrefiring)*FF,&phiCPeven_CMS_PreFire_13TeVDown_hgsAR,&phiCPeven_CMS_PreFire_13TeVDown_zttAR,&phiCPeven_CMS_PreFire_13TeVDown_fkjAR,false);
          if(NJets == 0) {
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd1jet0Up,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd1jet0Down,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets1jet0Up,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets1jet0Down,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd2jet0Up,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd2jet0Down,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets2jet0Up,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets2jet0Down,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdmetUp,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_zttAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdmetDown,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_zttAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdlptUp,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdlptDown,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdsystUp,&phiCPeven_ff_mt_qcd_syst_njets0Up_hgsAR,&phiCPeven_ff_mt_qcd_syst_njets0Up_zttAR,&phiCPeven_ff_mt_qcd_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdsystDown,&phiCPeven_ff_mt_qcd_syst_njets0Down_hgsAR,&phiCPeven_ff_mt_qcd_syst_njets0Down_zttAR,&phiCPeven_ff_mt_qcd_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetsmetUp,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_zttAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetsmetDown,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_zttAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetslptUp,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetslptDown,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetssystUp,&phiCPeven_ff_mt_wjets_syst_njets0Up_hgsAR,&phiCPeven_ff_mt_wjets_syst_njets0Up_zttAR,&phiCPeven_ff_mt_wjets_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetssystDown,&phiCPeven_ff_mt_wjets_syst_njets0Down_hgsAR,&phiCPeven_ff_mt_wjets_syst_njets0Down_zttAR,&phiCPeven_ff_mt_wjets_syst_njets0Down_fkjAR,false);
          }
          if(NJets == 1) {
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd1jet1Up,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd1jet1Down,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets1jet1Up,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets1jet1Down,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd2jet1Up,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd2jet1Down,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets2jet1Up,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets2jet1Down,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdmetUp,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_zttAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdmetDown,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_zttAR,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdlptUp,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdlptDown,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdsystUp,&phiCPeven_ff_mt_qcd_syst_njets1Up_hgsAR,&phiCPeven_ff_mt_qcd_syst_njets1Up_zttAR,&phiCPeven_ff_mt_qcd_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcdsystDown,&phiCPeven_ff_mt_qcd_syst_njets1Down_hgsAR,&phiCPeven_ff_mt_qcd_syst_njets1Down_zttAR,&phiCPeven_ff_mt_qcd_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetsmetUp,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_zttAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetsmetDown,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_zttAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetslptUp,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetslptDown,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetssystUp,&phiCPeven_ff_mt_wjets_syst_njets1Up_hgsAR,&phiCPeven_ff_mt_wjets_syst_njets1Up_zttAR,&phiCPeven_ff_mt_wjets_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetssystDown,&phiCPeven_ff_mt_wjets_syst_njets1Down_hgsAR,&phiCPeven_ff_mt_wjets_syst_njets1Down_zttAR,&phiCPeven_ff_mt_wjets_syst_njets1Down_fkjAR,false);
          }
          if(NJets == 2) {
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd1jet2Up,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd1jet2Down,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets1jet2Up,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets1jet2Down,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd2jet2Up,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFqcd2jet2Down,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets2jet2Up,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjets2jet2Down,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetsmetUp,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_zttAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetsmetDown,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_zttAR,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetslptUp,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetslptDown,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetssystUp,&phiCPeven_ff_mt_wjets_syst_njets2Up_hgsAR,&phiCPeven_ff_mt_wjets_syst_njets2Up_zttAR,&phiCPeven_ff_mt_wjets_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wEven*FFwjetssystDown,&phiCPeven_ff_mt_wjets_syst_njets2Down_hgsAR,&phiCPeven_ff_mt_wjets_syst_njets2Down_zttAR,&phiCPeven_ff_mt_wjets_syst_njets2Down_fkjAR,false);
          }
	}  
        if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_t_3prong_13TeVUp_hgsAR,&phiCPeven_CMS_scale_t_3prong_13TeVUp_zttAR,&phiCPeven_CMS_scale_t_3prong_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_t_3prong_13TeVDown_hgsAR,&phiCPeven_CMS_scale_t_3prong_13TeVDown_zttAR,&phiCPeven_CMS_scale_t_3prong_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_mu_13TeVUp_hgsAR,&phiCPeven_CMS_scale_mu_13TeVUp_zttAR,&phiCPeven_CMS_scale_mu_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_mu_13TeVDown_hgsAR,&phiCPeven_CMS_scale_mu_13TeVDown_zttAR,&phiCPeven_CMS_scale_mu_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgsAR,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_zttAR,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgsAR,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_zttAR,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_HF_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_HF_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_HF_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_HF_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkjAR,false);
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_EC2_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_EC2_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_EC2_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_EC2_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_res_j_13TeVUp_hgsAR,&phiCPeven_CMS_res_j_13TeVUp_zttAR,&phiCPeven_CMS_res_j_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_res_j_13TeVDown_hgsAR,&phiCPeven_CMS_res_j_13TeVDown_zttAR,&phiCPeven_CMS_res_j_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgsAR,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_zttAR,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgsAR,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_zttAR,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgsAR,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_zttAR,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgsAR,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_zttAR,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgsAR,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_zttAR,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,wEven*FF,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgsAR,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_zttAR,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkjAR,false);
	//
	if(Selection::Get_SysType() == "default") {
          Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_nominal_hgsAR,&phiCPodd_nominal_zttAR,&phiCPodd_nominal_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wTrgUp)/wTrg)*FF,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wTrgDown)/wTrg)*FF,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wBtagUp)/wBtag)*FF,&phiCPodd_CMS_eff_b_13TeVUp_hgsAR,&phiCPodd_CMS_eff_b_13TeVUp_zttAR,&phiCPodd_CMS_eff_b_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wBtagDown)/wBtag)*FF,&phiCPodd_CMS_eff_b_13TeVDown_hgsAR,&phiCPodd_CMS_eff_b_13TeVDown_zttAR,&phiCPodd_CMS_eff_b_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wZpTUp)/wZpT)*FF,&phiCPodd_CMS_htt_dyShape_13TeVUp_hgsAR,&phiCPodd_CMS_htt_dyShape_13TeVUp_zttAR,&phiCPodd_CMS_htt_dyShape_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wZpTDown)/wZpT)*FF,&phiCPodd_CMS_htt_dyShape_13TeVDown_hgsAR,&phiCPodd_CMS_htt_dyShape_13TeVDown_zttAR,&phiCPodd_CMS_htt_dyShape_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wToppTUp)/wToppT)*FF,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgsAR,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_zttAR,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wToppTDown)/wToppT)*FF,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgsAR,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_zttAR,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*wScaleUp*FF,&phiCPodd_CMS_scale_gg_13TeVUp_hgsAR,&phiCPodd_CMS_scale_gg_13TeVUp_zttAR,&phiCPodd_CMS_scale_gg_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*wScaleDown*FF,&phiCPodd_CMS_scale_gg_13TeVDown_hgsAR,&phiCPodd_CMS_scale_gg_13TeVDown_zttAR,&phiCPodd_CMS_scale_gg_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*wPSISRUp*FF,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgsAR,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_zttAR,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*wPSISRDown*FF,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgsAR,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_zttAR,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*wPSFSRUp*FF,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgsAR,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_zttAR,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*wPSFSRDown*FF,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgsAR,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_zttAR,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_ff_mt_sub_systUp_hgsAR,&phiCPodd_ff_mt_sub_systUp_zttAR,&phiCPodd_ff_mt_sub_systUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_ff_mt_sub_systDown_hgsAR,&phiCPodd_ff_mt_sub_systDown_zttAR,&phiCPodd_ff_mt_sub_systDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_ttbar_embeded_13TeVUp_hgsAR,&phiCPodd_CMS_ttbar_embeded_13TeVUp_zttAR,&phiCPodd_CMS_ttbar_embeded_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_ttbar_embeded_13TeVDown_hgsAR,&phiCPodd_CMS_ttbar_embeded_13TeVDown_zttAR,&phiCPodd_CMS_ttbar_embeded_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wPrefiringUp)/wPrefiring)*FF,&phiCPodd_CMS_PreFire_13TeVUp_hgsAR,&phiCPodd_CMS_PreFire_13TeVUp_zttAR,&phiCPodd_CMS_PreFire_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wOdd*wPrefiringDown)/wPrefiring)*FF,&phiCPodd_CMS_PreFire_13TeVDown_hgsAR,&phiCPodd_CMS_PreFire_13TeVDown_zttAR,&phiCPodd_CMS_PreFire_13TeVDown_fkjAR,false);
          if(NJets == 0) {
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd1jet0Up,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd1jet0Down,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets1jet0Up,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets1jet0Down,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd2jet0Up,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd2jet0Down,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets2jet0Up,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets2jet0Down,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdmetUp,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_zttAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdmetDown,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_zttAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdlptUp,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdlptDown,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdsystUp,&phiCPodd_ff_mt_qcd_syst_njets0Up_hgsAR,&phiCPodd_ff_mt_qcd_syst_njets0Up_zttAR,&phiCPodd_ff_mt_qcd_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdsystDown,&phiCPodd_ff_mt_qcd_syst_njets0Down_hgsAR,&phiCPodd_ff_mt_qcd_syst_njets0Down_zttAR,&phiCPodd_ff_mt_qcd_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetsmetUp,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_zttAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetsmetDown,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_zttAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetslptUp,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetslptDown,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetssystUp,&phiCPodd_ff_mt_wjets_syst_njets0Up_hgsAR,&phiCPodd_ff_mt_wjets_syst_njets0Up_zttAR,&phiCPodd_ff_mt_wjets_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetssystDown,&phiCPodd_ff_mt_wjets_syst_njets0Down_hgsAR,&phiCPodd_ff_mt_wjets_syst_njets0Down_zttAR,&phiCPodd_ff_mt_wjets_syst_njets0Down_fkjAR,false);
          }
          if(NJets == 1) {
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd1jet1Up,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd1jet1Down,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets1jet1Up,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets1jet1Down,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd2jet1Up,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd2jet1Down,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets2jet1Up,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets2jet1Down,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdmetUp,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_zttAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdmetDown,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_zttAR,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdlptUp,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdlptDown,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdsystUp,&phiCPodd_ff_mt_qcd_syst_njets1Up_hgsAR,&phiCPodd_ff_mt_qcd_syst_njets1Up_zttAR,&phiCPodd_ff_mt_qcd_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcdsystDown,&phiCPodd_ff_mt_qcd_syst_njets1Down_hgsAR,&phiCPodd_ff_mt_qcd_syst_njets1Down_zttAR,&phiCPodd_ff_mt_qcd_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetsmetUp,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_zttAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetsmetDown,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_zttAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetslptUp,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetslptDown,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetssystUp,&phiCPodd_ff_mt_wjets_syst_njets1Up_hgsAR,&phiCPodd_ff_mt_wjets_syst_njets1Up_zttAR,&phiCPodd_ff_mt_wjets_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetssystDown,&phiCPodd_ff_mt_wjets_syst_njets1Down_hgsAR,&phiCPodd_ff_mt_wjets_syst_njets1Down_zttAR,&phiCPodd_ff_mt_wjets_syst_njets1Down_fkjAR,false);
          }
          if(NJets == 2) {
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd1jet2Up,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd1jet2Down,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets1jet2Up,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets1jet2Down,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd2jet2Up,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFqcd2jet2Down,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets2jet2Up,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjets2jet2Down,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetsmetUp,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_zttAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetsmetDown,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_zttAR,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetslptUp,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetslptDown,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetssystUp,&phiCPodd_ff_mt_wjets_syst_njets2Up_hgsAR,&phiCPodd_ff_mt_wjets_syst_njets2Up_zttAR,&phiCPodd_ff_mt_wjets_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wOdd*FFwjetssystDown,&phiCPodd_ff_mt_wjets_syst_njets2Down_hgsAR,&phiCPodd_ff_mt_wjets_syst_njets2Down_zttAR,&phiCPodd_ff_mt_wjets_syst_njets2Down_fkjAR,false);
          }
	}  
        if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_t_3prong_13TeVUp_hgsAR,&phiCPodd_CMS_scale_t_3prong_13TeVUp_zttAR,&phiCPodd_CMS_scale_t_3prong_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_t_3prong_13TeVDown_hgsAR,&phiCPodd_CMS_scale_t_3prong_13TeVDown_zttAR,&phiCPodd_CMS_scale_t_3prong_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_mu_13TeVUp_hgsAR,&phiCPodd_CMS_scale_mu_13TeVUp_zttAR,&phiCPodd_CMS_scale_mu_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_mu_13TeVDown_hgsAR,&phiCPodd_CMS_scale_mu_13TeVDown_zttAR,&phiCPodd_CMS_scale_mu_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgsAR,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_zttAR,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgsAR,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_zttAR,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_HF_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_HF_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_HF_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_HF_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkjAR,false);
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_EC2_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_EC2_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_EC2_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_EC2_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_res_j_13TeVUp_hgsAR,&phiCPodd_CMS_res_j_13TeVUp_zttAR,&phiCPodd_CMS_res_j_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_res_j_13TeVDown_hgsAR,&phiCPodd_CMS_res_j_13TeVDown_zttAR,&phiCPodd_CMS_res_j_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgsAR,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_zttAR,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgsAR,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_zttAR,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgsAR,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_zttAR,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgsAR,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_zttAR,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgsAR,&phiCPodd_CMS_scale_met_unclustered_13TeVUp_zttAR,&phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,wOdd*FF,&phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgsAR,&phiCPodd_CMS_scale_met_unclustered_13TeVDown_zttAR,&phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkjAR,false);
	//
	if(Selection::Get_SysType() == "default") {
          Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_nominal_hgsAR,&phiCPMM_nominal_zttAR,&phiCPMM_nominal_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wTrgUp)/wTrg)*FF,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wTrgDown)/wTrg)*FF,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wBtagUp)/wBtag)*FF,&phiCPMM_CMS_eff_b_13TeVUp_hgsAR,&phiCPMM_CMS_eff_b_13TeVUp_zttAR,&phiCPMM_CMS_eff_b_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wBtagDown)/wBtag)*FF,&phiCPMM_CMS_eff_b_13TeVDown_hgsAR,&phiCPMM_CMS_eff_b_13TeVDown_zttAR,&phiCPMM_CMS_eff_b_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wZpTUp)/wZpT)*FF,&phiCPMM_CMS_htt_dyShape_13TeVUp_hgsAR,&phiCPMM_CMS_htt_dyShape_13TeVUp_zttAR,&phiCPMM_CMS_htt_dyShape_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wZpTDown)/wZpT)*FF,&phiCPMM_CMS_htt_dyShape_13TeVDown_hgsAR,&phiCPMM_CMS_htt_dyShape_13TeVDown_zttAR,&phiCPMM_CMS_htt_dyShape_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wToppTUp)/wToppT)*FF,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgsAR,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_zttAR,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wToppTDown)/wToppT)*FF,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgsAR,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_zttAR,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*wScaleUp*FF,&phiCPMM_CMS_scale_gg_13TeVUp_hgsAR,&phiCPMM_CMS_scale_gg_13TeVUp_zttAR,&phiCPMM_CMS_scale_gg_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*wScaleDown*FF,&phiCPMM_CMS_scale_gg_13TeVDown_hgsAR,&phiCPMM_CMS_scale_gg_13TeVDown_zttAR,&phiCPMM_CMS_scale_gg_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*wPSISRUp*FF,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgsAR,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_zttAR,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*wPSISRDown*FF,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgsAR,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_zttAR,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*wPSFSRUp*FF,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgsAR,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_zttAR,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*wPSFSRDown*FF,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgsAR,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_zttAR,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_ff_mt_sub_systUp_hgsAR,&phiCPMM_ff_mt_sub_systUp_zttAR,&phiCPMM_ff_mt_sub_systUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_ff_mt_sub_systDown_hgsAR,&phiCPMM_ff_mt_sub_systDown_zttAR,&phiCPMM_ff_mt_sub_systDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_ttbar_embeded_13TeVUp_hgsAR,&phiCPMM_CMS_ttbar_embeded_13TeVUp_zttAR,&phiCPMM_CMS_ttbar_embeded_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_ttbar_embeded_13TeVDown_hgsAR,&phiCPMM_CMS_ttbar_embeded_13TeVDown_zttAR,&phiCPMM_CMS_ttbar_embeded_13TeVDown_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wPrefiringUp)/wPrefiring)*FF,&phiCPMM_CMS_PreFire_13TeVUp_hgsAR,&phiCPMM_CMS_PreFire_13TeVUp_zttAR,&phiCPMM_CMS_PreFire_13TeVUp_fkjAR,false);
          Ntp->FillHist(t,Angle,max_pair,((wMM*wPrefiringDown)/wPrefiring)*FF,&phiCPMM_CMS_PreFire_13TeVDown_hgsAR,&phiCPMM_CMS_PreFire_13TeVDown_zttAR,&phiCPMM_CMS_PreFire_13TeVDown_fkjAR,false);
          if(NJets == 0) {
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd1jet0Up,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd1jet0Down,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets1jet0Up,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets1jet0Down,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd2jet0Up,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd2jet0Down,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR,false); //
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets2jet0Up,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets2jet0Down,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdmetUp,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_zttAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdmetDown,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_zttAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdlptUp,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdlptDown,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdsystUp,&phiCPMM_ff_mt_qcd_syst_njets0Up_hgsAR,&phiCPMM_ff_mt_qcd_syst_njets0Up_zttAR,&phiCPMM_ff_mt_qcd_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdsystDown,&phiCPMM_ff_mt_qcd_syst_njets0Down_hgsAR,&phiCPMM_ff_mt_qcd_syst_njets0Down_zttAR,&phiCPMM_ff_mt_qcd_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetsmetUp,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_zttAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetsmetDown,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_zttAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetslptUp,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetslptDown,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetssystUp,&phiCPMM_ff_mt_wjets_syst_njets0Up_hgsAR,&phiCPMM_ff_mt_wjets_syst_njets0Up_zttAR,&phiCPMM_ff_mt_wjets_syst_njets0Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetssystDown,&phiCPMM_ff_mt_wjets_syst_njets0Down_hgsAR,&phiCPMM_ff_mt_wjets_syst_njets0Down_zttAR,&phiCPMM_ff_mt_wjets_syst_njets0Down_fkjAR,false);
          }
          if(NJets == 1) {
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd1jet1Up,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd1jet1Down,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets1jet1Up,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets1jet1Down,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd2jet1Up,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd2jet1Down,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets2jet1Up,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets2jet1Down,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdmetUp,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_zttAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdmetDown,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_zttAR,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdlptUp,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdlptDown,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdsystUp,&phiCPMM_ff_mt_qcd_syst_njets1Up_hgsAR,&phiCPMM_ff_mt_qcd_syst_njets1Up_zttAR,&phiCPMM_ff_mt_qcd_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcdsystDown,&phiCPMM_ff_mt_qcd_syst_njets1Down_hgsAR,&phiCPMM_ff_mt_qcd_syst_njets1Down_zttAR,&phiCPMM_ff_mt_qcd_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetsmetUp,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_zttAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetsmetDown,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_zttAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetslptUp,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetslptDown,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetssystUp,&phiCPMM_ff_mt_wjets_syst_njets1Up_hgsAR,&phiCPMM_ff_mt_wjets_syst_njets1Up_zttAR,&phiCPMM_ff_mt_wjets_syst_njets1Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetssystDown,&phiCPMM_ff_mt_wjets_syst_njets1Down_hgsAR,&phiCPMM_ff_mt_wjets_syst_njets1Down_zttAR,&phiCPMM_ff_mt_wjets_syst_njets1Down_fkjAR,false);
          }
          if(NJets == 2) {
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd1jet2Up,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd1jet2Down,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets1jet2Up,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets1jet2Down,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd2jet2Up,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFqcd2jet2Down,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets2jet2Up,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjets2jet2Down,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetsmetUp,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_zttAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetsmetDown,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_zttAR,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetslptUp,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetslptDown,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetssystUp,&phiCPMM_ff_mt_wjets_syst_njets2Up_hgsAR,&phiCPMM_ff_mt_wjets_syst_njets2Up_zttAR,&phiCPMM_ff_mt_wjets_syst_njets2Up_fkjAR,false);
            Ntp->FillHist(t,Angle,max_pair,wMM*FFwjetssystDown,&phiCPMM_ff_mt_wjets_syst_njets2Down_hgsAR,&phiCPMM_ff_mt_wjets_syst_njets2Down_zttAR,&phiCPMM_ff_mt_wjets_syst_njets2Down_fkjAR,false);
          }
	}  
        if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_t_3prong_13TeVUp_hgsAR,&phiCPMM_CMS_scale_t_3prong_13TeVUp_zttAR,&phiCPMM_CMS_scale_t_3prong_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_t_3prong_13TeVDown_hgsAR,&phiCPMM_CMS_scale_t_3prong_13TeVDown_zttAR,&phiCPMM_CMS_scale_t_3prong_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_mu_13TeVUp_hgsAR,&phiCPMM_CMS_scale_mu_13TeVUp_zttAR,&phiCPMM_CMS_scale_mu_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_mu_13TeVDown_hgsAR,&phiCPMM_CMS_scale_mu_13TeVDown_zttAR,&phiCPMM_CMS_scale_mu_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgsAR,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_zttAR,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgsAR,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_zttAR,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_HF_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_HF_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_HF_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_HF_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkjAR,false);
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_EC2_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_EC2_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_EC2_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_EC2_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR,false);
        }
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR,false);
	  if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR,false);
	  if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR,false);
        }
        if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_res_j_13TeVUp_hgsAR,&phiCPMM_CMS_res_j_13TeVUp_zttAR,&phiCPMM_CMS_res_j_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_res_j_13TeVDown_hgsAR,&phiCPMM_CMS_res_j_13TeVDown_zttAR,&phiCPMM_CMS_res_j_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgsAR,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_zttAR,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgsAR,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_zttAR,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgsAR,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_zttAR,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgsAR,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_zttAR,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkjAR,false);
        if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgsAR,&phiCPMM_CMS_scale_met_unclustered_13TeVUp_zttAR,&phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkjAR,false);
        if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,wMM*FF,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgsAR,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_zttAR,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkjAR,false);
      }
    }
  }
  if(isSignalRegion && genMatched == true) {
    MuonpT.at(t).Fill(Ntp->muPt(),w);
    TaupT.at(t).Fill(Ntp->tauPt(),w);
    DitaupT.at(t).Fill(Ntp->ditauPt(),w);
    Njets.at(t).Fill(Ntp->Njets(),w);
    LeadingJetpT.at(t).Fill(Ntp->leadingjetPt(),w);
    SubleadingJetpT.at(t).Fill(Ntp->subleadingjetPt(),w);
    DijetpT.at(t).Fill(Ntp->dijetPt(),w);
    DijetMass.at(t).Fill(Ntp->dijetMass(),w);
    DijetDeltaEta.at(t).Fill(Ntp->dijetdeltaEta(),w);
    VisibleMass.at(t).Fill(Ntp->pairvisMass(),w);
    FastMTTditauMass.at(t).Fill(Ntp->fastMTTmass(),w);
    PUPPImet.at(t).Fill(Ntp->PUPPImet(),w);
    MuMETmt.at(t).Fill(Ntp->muMETmt(),w);
    ResMuonPt.at(t).Fill((genmup4.Pt() - mup4.Pt())/genmup4.Pt(),w);
    ResMuonEta.at(t).Fill((genmup4.Eta() - mup4.Eta())/genmup4.Eta(),w);
    ResTauPt.at(t).Fill((gentaup4.Pt() - taup4.Pt())/gentaup4.Pt(),w);
    ResTauEta.at(t).Fill((gentaup4.Eta() - taup4.Eta())/gentaup4.Eta(),w);
    ResPVx.at(t).Fill((Ntp->genPVx() - Ntp->pvx())/Ntp->genPVx(),w);
    ResPVy.at(t).Fill((Ntp->genPVy() - Ntp->pvy())/Ntp->genPVy(),w);
    ResPVz.at(t).Fill((Ntp->genPVz() - Ntp->pvz())/Ntp->genPVz(),w);
    ResSVx.at(t).Fill((Ntp->genTauSVx() - Ntp->tauSVx())/Ntp->genTauSVx(),w);
    ResSVy.at(t).Fill((Ntp->genTauSVy() - Ntp->tauSVy())/Ntp->genTauSVy(),w);
    ResSVz.at(t).Fill((Ntp->genTauSVz() - Ntp->tauSVz())/Ntp->genTauSVz(),w);
    if(max_pair.second == 0) {
      BDTscoreHiggs.at(t).Fill(max_pair.first,w);
      if(isA1MU) BDTscoreA1MUHiggs.at(t).Fill(max_pair.first,w);
    }
    else if(max_pair.second == 1) {
      BDTscoreZTT.at(t).Fill(max_pair.first,w);
      if(isA1MU) BDTscoreA1MUZTT.at(t).Fill(max_pair.first,w);
    }
    else if(max_pair.second == 2) {
      BDTscoreJetFakes.at(t).Fill(max_pair.first,w);
      if(isA1MU) BDTscoreA1MUJetFakes.at(t).Fill(max_pair.first,w);
    }
    //
    if(isA1MU) {
      if(Ntp->isSignal()) {
        wOdd *= Ntp->wOdd();
        wEven *= Ntp->wEven();
        wMM *= Ntp->wMM();
      }
      if(Selection::Get_SysType() == "default") {
        PhiCPEvenPV.at(t).Fill(Ntp->pvPhiCP(),wEven);
        PhiCPOddPV.at(t).Fill(Ntp->pvPhiCP(),wOdd);
        PhiCPMMPV.at(t).Fill(Ntp->pvPhiCP(),wMM);
        PhiCPEvenDP.at(t).Fill(Ntp->dpPhiCP(),wEven);
        PhiCPOddDP.at(t).Fill(Ntp->dpPhiCP(),wOdd);
        PhiCPMMDP.at(t).Fill(Ntp->dpPhiCP(),wMM);
	//
	if(Ntp->isVV() || Ntp->isTTbar()) {
	  Ntp->FillHist(t,Angle,max_pair,w,&ttbar_contamination_hgs,&ttbar_contamination_ztt,&ttbar_contamination_fkj,false);
	}
	//
        Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_nominal_hgs,&phiCPeven_nominal_ztt,&phiCPeven_nominal_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wTrgUp)/wTrg,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wTrgDown)/wTrg,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt,&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wBtagUp)/wBtag,&phiCPeven_CMS_eff_b_13TeVUp_hgs,&phiCPeven_CMS_eff_b_13TeVUp_ztt,&phiCPeven_CMS_eff_b_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wBtagDown)/wBtag,&phiCPeven_CMS_eff_b_13TeVDown_hgs,&phiCPeven_CMS_eff_b_13TeVDown_ztt,&phiCPeven_CMS_eff_b_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wZpTUp)/wZpT,&phiCPeven_CMS_htt_dyShape_13TeVUp_hgs,&phiCPeven_CMS_htt_dyShape_13TeVUp_ztt,&phiCPeven_CMS_htt_dyShape_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wZpTDown)/wZpT,&phiCPeven_CMS_htt_dyShape_13TeVDown_hgs,&phiCPeven_CMS_htt_dyShape_13TeVDown_ztt,&phiCPeven_CMS_htt_dyShape_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wToppTUp)/wToppT,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgs,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_ztt,&phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wToppTDown)/wToppT,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgs,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_ztt,&phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven*wScaleUp,&phiCPeven_CMS_scale_gg_13TeVUp_hgs,&phiCPeven_CMS_scale_gg_13TeVUp_ztt,&phiCPeven_CMS_scale_gg_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven*wScaleDown,&phiCPeven_CMS_scale_gg_13TeVDown_hgs,&phiCPeven_CMS_scale_gg_13TeVDown_ztt,&phiCPeven_CMS_scale_gg_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven*wPSISRUp,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgs,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_ztt,&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven*wPSISRDown,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgs,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_ztt,&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven*wPSFSRUp,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgs,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_ztt,&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven*wPSFSRDown,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgs,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_ztt,&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_sub_systUp_hgs,&phiCPeven_ff_mt_sub_systUp_ztt,&phiCPeven_ff_mt_sub_systUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_sub_systDown_hgs,&phiCPeven_ff_mt_sub_systDown_ztt,&phiCPeven_ff_mt_sub_systDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_ttbar_embeded_13TeVUp_hgs,&phiCPeven_CMS_ttbar_embeded_13TeVUp_ztt,&phiCPeven_CMS_ttbar_embeded_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_ttbar_embeded_13TeVDown_hgs,&phiCPeven_CMS_ttbar_embeded_13TeVDown_ztt,&phiCPeven_CMS_ttbar_embeded_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wPrefiringUp)/wPrefiring,&phiCPeven_CMS_PreFire_13TeVUp_hgs,&phiCPeven_CMS_PreFire_13TeVUp_ztt,&phiCPeven_CMS_PreFire_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wEven*wPrefiringDown)/wPrefiring,&phiCPeven_CMS_PreFire_13TeVDown_hgs,&phiCPeven_CMS_PreFire_13TeVDown_ztt,&phiCPeven_CMS_PreFire_13TeVDown_fkj,false);
	if(NJets == 0) {
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_syst_njets0Up_hgs,&phiCPeven_ff_mt_qcd_syst_njets0Up_ztt,&phiCPeven_ff_mt_qcd_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_syst_njets0Down_hgs,&phiCPeven_ff_mt_qcd_syst_njets0Down_ztt,&phiCPeven_ff_mt_qcd_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_syst_njets0Up_hgs,&phiCPeven_ff_mt_wjets_syst_njets0Up_ztt,&phiCPeven_ff_mt_wjets_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_syst_njets0Down_hgs,&phiCPeven_ff_mt_wjets_syst_njets0Down_ztt,&phiCPeven_ff_mt_wjets_syst_njets0Down_fkj,false);
	}
	if(NJets == 1) {
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_syst_njets1Up_hgs,&phiCPeven_ff_mt_qcd_syst_njets1Up_ztt,&phiCPeven_ff_mt_qcd_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_syst_njets1Down_hgs,&phiCPeven_ff_mt_qcd_syst_njets1Down_ztt,&phiCPeven_ff_mt_qcd_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_syst_njets1Up_hgs,&phiCPeven_ff_mt_wjets_syst_njets1Up_ztt,&phiCPeven_ff_mt_wjets_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_syst_njets1Down_hgs,&phiCPeven_ff_mt_wjets_syst_njets1Down_ztt,&phiCPeven_ff_mt_wjets_syst_njets1Down_fkj,false);
	}
	if(NJets == 2) {
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt,&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgs,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_ztt,&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt,&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_syst_njets2Up_hgs,&phiCPeven_ff_mt_wjets_syst_njets2Up_ztt,&phiCPeven_ff_mt_wjets_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_ff_mt_wjets_syst_njets2Down_hgs,&phiCPeven_ff_mt_wjets_syst_njets2Down_ztt,&phiCPeven_ff_mt_wjets_syst_njets2Down_fkj,false);
	}
      }
      if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_t_3prong_13TeVUp_hgs,&phiCPeven_CMS_scale_t_3prong_13TeVUp_ztt,&phiCPeven_CMS_scale_t_3prong_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_t_3prong_13TeVDown_hgs,&phiCPeven_CMS_scale_t_3prong_13TeVDown_ztt,&phiCPeven_CMS_scale_t_3prong_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_mu_13TeVUp_hgs,&phiCPeven_CMS_scale_mu_13TeVUp_ztt,&phiCPeven_CMS_scale_mu_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_mu_13TeVDown_hgs,&phiCPeven_CMS_scale_mu_13TeVDown_ztt,&phiCPeven_CMS_scale_mu_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgs,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_ztt,&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgs,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_ztt,&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "HF_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "HF_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "EC2_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "EC2_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "Absolute_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "Absolute_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt,&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_res_j_13TeVUp_hgs,&phiCPeven_CMS_res_j_13TeVUp_ztt,&phiCPeven_CMS_res_j_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_res_j_13TeVDown_hgs,&phiCPeven_CMS_res_j_13TeVDown_ztt,&phiCPeven_CMS_res_j_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgs,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_ztt,&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgs,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_ztt,&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgs,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_ztt,&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgs,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_ztt,&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,wEven,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj,false);
      //
      if(Selection::Get_SysType() == "default") {
        Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_nominal_hgs,&phiCPodd_nominal_ztt,&phiCPodd_nominal_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wTrgUp)/wTrg,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wTrgDown)/wTrg,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt,&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wBtagUp)/wBtag,&phiCPodd_CMS_eff_b_13TeVUp_hgs,&phiCPodd_CMS_eff_b_13TeVUp_ztt,&phiCPodd_CMS_eff_b_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wBtagDown)/wBtag,&phiCPodd_CMS_eff_b_13TeVDown_hgs,&phiCPodd_CMS_eff_b_13TeVDown_ztt,&phiCPodd_CMS_eff_b_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wZpTUp)/wZpT,&phiCPodd_CMS_htt_dyShape_13TeVUp_hgs,&phiCPodd_CMS_htt_dyShape_13TeVUp_ztt,&phiCPodd_CMS_htt_dyShape_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wZpTDown)/wZpT,&phiCPodd_CMS_htt_dyShape_13TeVDown_hgs,&phiCPodd_CMS_htt_dyShape_13TeVDown_ztt,&phiCPodd_CMS_htt_dyShape_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wToppTUp)/wToppT,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgs,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_ztt,&phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wToppTDown)/wToppT,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgs,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_ztt,&phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd*wScaleUp,&phiCPodd_CMS_scale_gg_13TeVUp_hgs,&phiCPodd_CMS_scale_gg_13TeVUp_ztt,&phiCPodd_CMS_scale_gg_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd*wScaleDown,&phiCPodd_CMS_scale_gg_13TeVDown_hgs,&phiCPodd_CMS_scale_gg_13TeVDown_ztt,&phiCPodd_CMS_scale_gg_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd*wPSISRUp,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgs,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_ztt,&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd*wPSISRDown,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgs,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_ztt,&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd*wPSFSRUp,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgs,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_ztt,&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd*wPSFSRDown,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgs,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_ztt,&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_sub_systUp_hgs,&phiCPodd_ff_mt_sub_systUp_ztt,&phiCPodd_ff_mt_sub_systUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_sub_systDown_hgs,&phiCPodd_ff_mt_sub_systDown_ztt,&phiCPodd_ff_mt_sub_systDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_ttbar_embeded_13TeVUp_hgs,&phiCPodd_CMS_ttbar_embeded_13TeVUp_ztt,&phiCPodd_CMS_ttbar_embeded_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_ttbar_embeded_13TeVDown_hgs,&phiCPodd_CMS_ttbar_embeded_13TeVDown_ztt,&phiCPodd_CMS_ttbar_embeded_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wPrefiringUp)/wPrefiring,&phiCPodd_CMS_PreFire_13TeVUp_hgs,&phiCPodd_CMS_PreFire_13TeVUp_ztt,&phiCPodd_CMS_PreFire_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wOdd*wPrefiringDown)/wPrefiring,&phiCPodd_CMS_PreFire_13TeVDown_hgs,&phiCPodd_CMS_PreFire_13TeVDown_ztt,&phiCPodd_CMS_PreFire_13TeVDown_fkj,false);
	if(NJets == 0) {
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_syst_njets0Up_hgs,&phiCPodd_ff_mt_qcd_syst_njets0Up_ztt,&phiCPodd_ff_mt_qcd_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_syst_njets0Down_hgs,&phiCPodd_ff_mt_qcd_syst_njets0Down_ztt,&phiCPodd_ff_mt_qcd_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_syst_njets0Up_hgs,&phiCPodd_ff_mt_wjets_syst_njets0Up_ztt,&phiCPodd_ff_mt_wjets_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_syst_njets0Down_hgs,&phiCPodd_ff_mt_wjets_syst_njets0Down_ztt,&phiCPodd_ff_mt_wjets_syst_njets0Down_fkj,false);
	}
	if(NJets == 1) {
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_syst_njets1Up_hgs,&phiCPodd_ff_mt_qcd_syst_njets1Up_ztt,&phiCPodd_ff_mt_qcd_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_syst_njets1Down_hgs,&phiCPodd_ff_mt_qcd_syst_njets1Down_ztt,&phiCPodd_ff_mt_qcd_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_syst_njets1Up_hgs,&phiCPodd_ff_mt_wjets_syst_njets1Up_ztt,&phiCPodd_ff_mt_wjets_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_syst_njets1Down_hgs,&phiCPodd_ff_mt_wjets_syst_njets1Down_ztt,&phiCPodd_ff_mt_wjets_syst_njets1Down_fkj,false);
	}
	if(NJets == 2) {
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt,&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgs,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_ztt,&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt,&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_syst_njets2Up_hgs,&phiCPodd_ff_mt_wjets_syst_njets2Up_ztt,&phiCPodd_ff_mt_wjets_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_ff_mt_wjets_syst_njets2Down_hgs,&phiCPodd_ff_mt_wjets_syst_njets2Down_ztt,&phiCPodd_ff_mt_wjets_syst_njets2Down_fkj,false);
	}
      }
      if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_t_3prong_13TeVUp_hgs,&phiCPodd_CMS_scale_t_3prong_13TeVUp_ztt,&phiCPodd_CMS_scale_t_3prong_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_t_3prong_13TeVDown_hgs,&phiCPodd_CMS_scale_t_3prong_13TeVDown_ztt,&phiCPodd_CMS_scale_t_3prong_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_mu_13TeVUp_hgs,&phiCPodd_CMS_scale_mu_13TeVUp_ztt,&phiCPodd_CMS_scale_mu_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_mu_13TeVDown_hgs,&phiCPodd_CMS_scale_mu_13TeVDown_ztt,&phiCPodd_CMS_scale_mu_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgs,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_ztt,&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgs,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_ztt,&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "HF_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "HF_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "EC2_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "EC2_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "Absolute_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "Absolute_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt,&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_res_j_13TeVUp_hgs,&phiCPodd_CMS_res_j_13TeVUp_ztt,&phiCPodd_CMS_res_j_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_res_j_13TeVDown_hgs,&phiCPodd_CMS_res_j_13TeVDown_ztt,&phiCPodd_CMS_res_j_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgs,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_ztt,&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgs,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_ztt,&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgs,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_ztt,&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgs,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_ztt,&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,wOdd,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj,false);
      //
      if(Selection::Get_SysType() == "default") {
        Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_nominal_hgs,&phiCPMM_nominal_ztt,&phiCPMM_nominal_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wTrgUp)/wTrg,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wTrgDown)/wTrg,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt,&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wBtagUp)/wBtag,&phiCPMM_CMS_eff_b_13TeVUp_hgs,&phiCPMM_CMS_eff_b_13TeVUp_ztt,&phiCPMM_CMS_eff_b_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wBtagDown)/wBtag,&phiCPMM_CMS_eff_b_13TeVDown_hgs,&phiCPMM_CMS_eff_b_13TeVDown_ztt,&phiCPMM_CMS_eff_b_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wZpTUp)/wZpT,&phiCPMM_CMS_htt_dyShape_13TeVUp_hgs,&phiCPMM_CMS_htt_dyShape_13TeVUp_ztt,&phiCPMM_CMS_htt_dyShape_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wZpTDown)/wZpT,&phiCPMM_CMS_htt_dyShape_13TeVDown_hgs,&phiCPMM_CMS_htt_dyShape_13TeVDown_ztt,&phiCPMM_CMS_htt_dyShape_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wToppTUp)/wToppT,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgs,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_ztt,&phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wToppTDown)/wToppT,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgs,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_ztt,&phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM*wScaleUp,&phiCPMM_CMS_scale_gg_13TeVUp_hgs,&phiCPMM_CMS_scale_gg_13TeVUp_ztt,&phiCPMM_CMS_scale_gg_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM*wScaleDown,&phiCPMM_CMS_scale_gg_13TeVDown_hgs,&phiCPMM_CMS_scale_gg_13TeVDown_ztt,&phiCPMM_CMS_scale_gg_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM*wPSISRUp,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgs,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_ztt,&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM*wPSISRDown,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgs,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_ztt,&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM*wPSFSRUp,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgs,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_ztt,&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM*wPSFSRDown,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgs,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_ztt,&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_sub_systUp_hgs,&phiCPMM_ff_mt_sub_systUp_ztt,&phiCPMM_ff_mt_sub_systUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_sub_systDown_hgs,&phiCPMM_ff_mt_sub_systDown_ztt,&phiCPMM_ff_mt_sub_systDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_ttbar_embeded_13TeVUp_hgs,&phiCPMM_CMS_ttbar_embeded_13TeVUp_ztt,&phiCPMM_CMS_ttbar_embeded_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_ttbar_embeded_13TeVDown_hgs,&phiCPMM_CMS_ttbar_embeded_13TeVDown_ztt,&phiCPMM_CMS_ttbar_embeded_13TeVDown_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wPrefiringUp)/wPrefiring,&phiCPMM_CMS_PreFire_13TeVUp_hgs,&phiCPMM_CMS_PreFire_13TeVUp_ztt,&phiCPMM_CMS_PreFire_13TeVUp_fkj,false);
        Ntp->FillHist(t,Angle,max_pair,(wMM*wPrefiringDown)/wPrefiring,&phiCPMM_CMS_PreFire_13TeVDown_hgs,&phiCPMM_CMS_PreFire_13TeVDown_ztt,&phiCPMM_CMS_PreFire_13TeVDown_fkj,false);
	if(NJets == 0) {
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_syst_njets0Up_hgs,&phiCPMM_ff_mt_qcd_syst_njets0Up_ztt,&phiCPMM_ff_mt_qcd_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_syst_njets0Down_hgs,&phiCPMM_ff_mt_qcd_syst_njets0Down_ztt,&phiCPMM_ff_mt_qcd_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_syst_njets0Up_hgs,&phiCPMM_ff_mt_wjets_syst_njets0Up_ztt,&phiCPMM_ff_mt_wjets_syst_njets0Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_syst_njets0Down_hgs,&phiCPMM_ff_mt_wjets_syst_njets0Down_ztt,&phiCPMM_ff_mt_wjets_syst_njets0Down_fkj,false);
	}
	if(NJets == 1) {
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_syst_njets1Up_hgs,&phiCPMM_ff_mt_qcd_syst_njets1Up_ztt,&phiCPMM_ff_mt_qcd_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_syst_njets1Down_hgs,&phiCPMM_ff_mt_qcd_syst_njets1Down_ztt,&phiCPMM_ff_mt_qcd_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_syst_njets1Up_hgs,&phiCPMM_ff_mt_wjets_syst_njets1Up_ztt,&phiCPMM_ff_mt_wjets_syst_njets1Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_syst_njets1Down_hgs,&phiCPMM_ff_mt_wjets_syst_njets1Down_ztt,&phiCPMM_ff_mt_wjets_syst_njets1Down_fkj,false);
	}
	if(NJets == 2) {
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt,&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgs,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_ztt,&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt,&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_syst_njets2Up_hgs,&phiCPMM_ff_mt_wjets_syst_njets2Up_ztt,&phiCPMM_ff_mt_wjets_syst_njets2Up_fkj,false);
	  Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_ff_mt_wjets_syst_njets2Down_hgs,&phiCPMM_ff_mt_wjets_syst_njets2Down_ztt,&phiCPMM_ff_mt_wjets_syst_njets2Down_fkj,false);
	}
      }
      if(Selection::Get_SysType() == "TESUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_t_3prong_13TeVUp_hgs,&phiCPMM_CMS_scale_t_3prong_13TeVUp_ztt,&phiCPMM_CMS_scale_t_3prong_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "TESDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_t_3prong_13TeVDown_hgs,&phiCPMM_CMS_scale_t_3prong_13TeVDown_ztt,&phiCPMM_CMS_scale_t_3prong_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "MESUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_mu_13TeVUp_hgs,&phiCPMM_CMS_scale_mu_13TeVUp_ztt,&phiCPMM_CMS_scale_mu_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "MESDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_mu_13TeVDown_hgs,&phiCPMM_CMS_scale_mu_13TeVDown_ztt,&phiCPMM_CMS_scale_mu_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "FlavorQCDUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgs,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_ztt,&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "FlavorQCDDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgs,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_ztt,&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "RelativeBalUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "RelativeBalDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "HFUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "HFDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "HF_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "HF_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "BBEC1Up") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "BBEC1Down") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "EC2Up") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "EC2Down") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "EC2_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "EC2_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "AbsoluteUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "AbsoluteDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "Absolute_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "Absolute_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj,false);
      }
      if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	if(theYear == 2016) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj,false);
	if(theYear == 2017) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj,false);
	if(theYear == 2018) Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt,&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj,false);
      }
      if(Selection::Get_SysType() == "JERUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_res_j_13TeVUp_hgs,&phiCPMM_CMS_res_j_13TeVUp_ztt,&phiCPMM_CMS_res_j_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "JERDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_res_j_13TeVDown_hgs,&phiCPMM_CMS_res_j_13TeVDown_ztt,&phiCPMM_CMS_res_j_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METResoUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgs,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_ztt,&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METResoDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgs,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_ztt,&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METScaleUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgs,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_ztt,&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METScaleDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgs,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_ztt,&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkj,false);
      if(Selection::Get_SysType() == "METUnclusteredUp") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj,false);
      if(Selection::Get_SysType() == "METUnclusteredDown") Ntp->FillHist(t,Angle,max_pair,wMM,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt,&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj,false);
    }
  }
} //do event

//  This is a function if you want to do something after the event loop
void HCPMuTau::Finish() {

  if(mode == RECONSTRUCT) {
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);

    double norm=1.;
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(CrossSectionandAcceptance.at(i)>0 || HConfig.GetID(i)==35 || HConfig.GetID(i)==20 || HConfig.GetID(i)==23 || HConfig.GetID(i)==30 || HConfig.GetID(i)==33){
        if(CrossSectionandAcceptance.at(i)>0)norm= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
        else norm=1.;

	// Add/Substract 10% of VV/TTbar contribution to embedded templates
	phiCPeven_CMS_ttbar_embeded_13TeVUp_hgs.at(HConfig.GetType(35)).Add(&ttbar_contamination_hgs.at(i),+0.1*norm);
        phiCPeven_CMS_ttbar_embeded_13TeVDown_hgs.at(HConfig.GetType(35)).Add(&ttbar_contamination_hgs.at(i),-0.1*norm);
        phiCPodd_CMS_ttbar_embeded_13TeVUp_hgs.at(HConfig.GetType(35)).Add(&ttbar_contamination_hgs.at(i),+0.1*norm);
        phiCPodd_CMS_ttbar_embeded_13TeVDown_hgs.at(HConfig.GetType(35)).Add(&ttbar_contamination_hgs.at(i),-0.1*norm);
        phiCPMM_CMS_ttbar_embeded_13TeVUp_hgs.at(HConfig.GetType(35)).Add(&ttbar_contamination_hgs.at(i),+0.1*norm);
        phiCPMM_CMS_ttbar_embeded_13TeVDown_hgs.at(HConfig.GetType(35)).Add(&ttbar_contamination_hgs.at(i),-0.1*norm);
        phiCPeven_CMS_ttbar_embeded_13TeVUp_ztt.at(HConfig.GetType(35)).Add(&ttbar_contamination_ztt.at(i),+0.1*norm);
        phiCPeven_CMS_ttbar_embeded_13TeVDown_ztt.at(HConfig.GetType(35)).Add(&ttbar_contamination_ztt.at(i),-0.1*norm);
        phiCPodd_CMS_ttbar_embeded_13TeVUp_ztt.at(HConfig.GetType(35)).Add(&ttbar_contamination_ztt.at(i),+0.1*norm);
        phiCPodd_CMS_ttbar_embeded_13TeVDown_ztt.at(HConfig.GetType(35)).Add(&ttbar_contamination_ztt.at(i),-0.1*norm);
        phiCPMM_CMS_ttbar_embeded_13TeVUp_ztt.at(HConfig.GetType(35)).Add(&ttbar_contamination_ztt.at(i),+0.1*norm);
        phiCPMM_CMS_ttbar_embeded_13TeVDown_ztt.at(HConfig.GetType(35)).Add(&ttbar_contamination_ztt.at(i),-0.1*norm);
        phiCPeven_CMS_ttbar_embeded_13TeVUp_fkj.at(HConfig.GetType(35)).Add(&ttbar_contamination_fkj.at(i),+0.1*norm);
        phiCPeven_CMS_ttbar_embeded_13TeVDown_fkj.at(HConfig.GetType(35)).Add(&ttbar_contamination_fkj.at(i),-0.1*norm);
        phiCPodd_CMS_ttbar_embeded_13TeVUp_fkj.at(HConfig.GetType(35)).Add(&ttbar_contamination_fkj.at(i),+0.1*norm);
        phiCPodd_CMS_ttbar_embeded_13TeVDown_fkj.at(HConfig.GetType(35)).Add(&ttbar_contamination_fkj.at(i),-0.1*norm);
        phiCPMM_CMS_ttbar_embeded_13TeVUp_fkj.at(HConfig.GetType(35)).Add(&ttbar_contamination_fkj.at(i),+0.1*norm);
        phiCPMM_CMS_ttbar_embeded_13TeVDown_fkj.at(HConfig.GetType(35)).Add(&ttbar_contamination_fkj.at(i),-0.1*norm);
	//
        cout<<"Soustraction du QCD: "<<endl;
        //
	MuonpT.at(1).Add(&MuonpTAR.at(i),-norm);
        TaupT.at(1).Add(&TaupTAR.at(i),-norm);
        DitaupT.at(1).Add(&DitaupTAR.at(i),-norm);
   	Njets.at(1).Add(&NjetsAR.at(i),-norm);
  	LeadingJetpT.at(1).Add(&LeadingJetpTAR.at(i),-norm);
   	SubleadingJetpT.at(1).Add(&SubleadingJetpTAR.at(i),-norm);
   	DijetpT.at(1).Add(&DijetpTAR.at(i),-norm);
   	DijetMass.at(1).Add(&DijetMassAR.at(i),-norm);
   	DijetDeltaEta.at(1).Add(&DijetDeltaEtaAR.at(i),-norm);
   	VisibleMass.at(1).Add(&VisibleMassAR.at(i),-norm);
   	FastMTTditauMass.at(1).Add(&FastMTTditauMassAR.at(i),-norm);
   	PUPPImet.at(1).Add(&PUPPImetAR.at(i),-norm);
   	MuMETmt.at(1).Add(&MuMETmtAR.at(i),-norm);
        BDTscoreHiggs.at(1).Add(&BDTscoreHiggsAR.at(i),-norm);
        BDTscoreA1MUHiggs.at(1).Add(&BDTscoreA1MUHiggsAR.at(i),-norm);
        BDTscoreZTT.at(1).Add(&BDTscoreZTTAR.at(i),-norm);
        BDTscoreA1MUZTT.at(1).Add(&BDTscoreA1MUZTTAR.at(i),-norm);
        BDTscoreJetFakes.at(1).Add(&BDTscoreJetFakesAR.at(i),-norm);
        BDTscoreA1MUJetFakes.at(1).Add(&BDTscoreA1MUJetFakesAR.at(i),-norm);
	//
	if(Selection::Get_SysType() == "default") {
          phiCPeven_nominal_hgs.at(1).Add(&phiCPeven_nominal_hgsAR.at(i),-norm);
	  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR.at(i),-norm);
	  phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_eff_b_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_eff_b_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_eff_b_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_eff_b_13TeVDown_hgsAR.at(i),-norm);
	  phiCPeven_CMS_htt_dyShape_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_htt_dyShape_13TeVUp_hgsAR.at(i),-norm);
	  phiCPeven_CMS_htt_dyShape_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_htt_dyShape_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_htt_ttbarShape_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_htt_ttbarShape_13TeVDown_hgsAR.at(i),-norm);
	  phiCPeven_CMS_scale_gg_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_gg_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_gg_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_gg_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_syst_njets0Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets0Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_syst_njets0Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets0Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_syst_njets1Up_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets1Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_qcd_syst_njets1Down_hgs.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets1Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_syst_njets0Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets0Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_syst_njets0Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets0Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_syst_njets1Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets1Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_syst_njets1Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets1Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_syst_njets2Up_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets2Up_hgsAR.at(i),-norm);
	  phiCPeven_ff_mt_wjets_syst_njets2Down_hgs.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets2Down_hgsAR.at(i),-norm);
          phiCPeven_ff_mt_sub_systUp_hgs.at(1).Add(&phiCPeven_ff_mt_sub_systUp_hgsAR.at(i),-norm*1.1);
          phiCPeven_ff_mt_sub_systDown_hgs.at(1).Add(&phiCPeven_ff_mt_sub_systDown_hgsAR.at(i),-norm*0.9);
          phiCPeven_CMS_ttbar_embeded_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_ttbar_embeded_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_ttbar_embeded_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_ttbar_embeded_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_PreFire_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_PreFire_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_PreFire_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_PreFire_13TeVDown_hgsAR.at(i),-norm);
          //
          phiCPeven_nominal_ztt.at(1).Add(&phiCPeven_nominal_zttAR.at(i),-norm);
          phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_eff_b_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_eff_b_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_eff_b_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_eff_b_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_dyShape_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_htt_dyShape_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_dyShape_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_htt_dyShape_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_ttbarShape_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_htt_ttbarShape_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_ttbarShape_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_htt_ttbarShape_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_gg_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_gg_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_gg_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_gg_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_PS_ISR_ggH_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_PS_ISR_ggH_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_PS_FSR_ggH_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_PS_FSR_ggH_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets0Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets0Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets0Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets0Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets1Up_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets1Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets1Down_ztt.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets1Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets0Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets0Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets0Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets0Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets1Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets1Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets1Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets1Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets2Up_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets2Up_zttAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets2Down_ztt.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets2Down_zttAR.at(i),-norm);
          phiCPeven_ff_mt_sub_systUp_ztt.at(1).Add(&phiCPeven_ff_mt_sub_systUp_zttAR.at(i),-norm*1.1);
          phiCPeven_ff_mt_sub_systDown_ztt.at(1).Add(&phiCPeven_ff_mt_sub_systDown_zttAR.at(i),-norm*0.9);
          phiCPeven_CMS_ttbar_embeded_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_ttbar_embeded_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_ttbar_embeded_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_ttbar_embeded_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_PreFire_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_PreFire_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_PreFire_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_PreFire_13TeVDown_zttAR.at(i),-norm);
	  //
	  phiCPeven_nominal_fkj.at(1).Add(&phiCPeven_nominal_fkjAR.at(i),-norm);
          phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_CMS_eff_b_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_eff_b_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_eff_b_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_eff_b_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_CMS_htt_dyShape_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_htt_dyShape_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_htt_dyShape_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_htt_dyShape_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_htt_ttbarShape_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_htt_ttbarShape_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_CMS_scale_gg_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_gg_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_scale_gg_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_gg_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_PS_ISR_ggH_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_PS_ISR_ggH_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_PS_FSR_ggH_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_PS_FSR_ggH_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets0Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets0Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets1Up_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_qcd_syst_njets1Down_fkj.at(1).Add(&phiCPeven_ff_mt_qcd_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets0Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets0Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets1Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets1Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets2Up_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_wjets_syst_njets2Down_fkj.at(1).Add(&phiCPeven_ff_mt_wjets_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPeven_ff_mt_sub_systUp_fkj.at(1).Add(&phiCPeven_ff_mt_sub_systUp_fkjAR.at(i),-norm*1.1);
          phiCPeven_ff_mt_sub_systDown_fkj.at(1).Add(&phiCPeven_ff_mt_sub_systDown_fkjAR.at(i),-norm*0.9);
          phiCPeven_CMS_ttbar_embeded_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_ttbar_embeded_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_ttbar_embeded_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_ttbar_embeded_13TeVDown_fkjAR.at(i),-norm);
          phiCPeven_CMS_PreFire_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_PreFire_13TeVUp_fkjAR.at(i),-norm);
          phiCPeven_CMS_PreFire_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_PreFire_13TeVDown_fkjAR.at(i),-norm);
	}
        if(Selection::Get_SysType() == "TESUp") {
          phiCPeven_CMS_scale_t_3prong_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_t_3prong_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_t_3prong_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_t_3prong_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_t_3prong_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_t_3prong_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "TESDown") {
          phiCPeven_CMS_scale_t_3prong_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_t_3prong_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_t_3prong_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_t_3prong_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_t_3prong_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_t_3prong_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "MESUp") {
          phiCPeven_CMS_scale_mu_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_mu_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_mu_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_mu_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_mu_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_mu_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "MESDown") {
          phiCPeven_CMS_scale_mu_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_mu_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_mu_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_mu_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_mu_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_mu_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "FlavorQCDUp") {
          phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_FlavorQCD_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_FlavorQCD_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "FlavorQCDDown") {
          phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_FlavorQCD_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_FlavorQCD_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "RelativeBalUp") {
          phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeBal_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "RelativeBalDown") {
          phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeBal_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HFUp") {
          phiCPeven_CMS_scale_j_HF_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_HF_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_HF_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HFDown") {
          phiCPeven_CMS_scale_j_HF_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_HF_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_HF_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2016_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2017_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2018_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2016_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2017_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2018_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_HF_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "BBEC1Up") {
          phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_BBEC1_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "BBEC1Down") {
          phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_BBEC1_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "EC2Up") {
          phiCPeven_CMS_scale_j_EC2_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_EC2_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_EC2_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "EC2Down") {
          phiCPeven_CMS_scale_j_EC2_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_EC2_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_EC2_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_EC2_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "AbsoluteUp") {
          phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_Absolute_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "AbsoluteDown") {
          phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_j_Absolute_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) {
	    phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "JERUp") {
          phiCPeven_CMS_res_j_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_res_j_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_res_j_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_res_j_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_res_j_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_res_j_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "JERDown") {
          phiCPeven_CMS_res_j_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_res_j_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_res_j_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_res_j_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_res_j_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_res_j_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METResoUp") {
          phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_reso_met_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_htt_boson_reso_met_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METResoDown") {
          phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_reso_met_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_htt_boson_reso_met_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METScaleUp") {
          phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_scale_met_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_htt_boson_scale_met_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METScaleDown") {
          phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_scale_met_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_htt_boson_scale_met_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METUnclusteredUp") {
          phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgs.at(1).Add(&phiCPeven_CMS_scale_met_unclustered_13TeVUp_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_met_unclustered_13TeVUp_ztt.at(1).Add(&phiCPeven_CMS_scale_met_unclustered_13TeVUp_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkj.at(1).Add(&phiCPeven_CMS_scale_met_unclustered_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METUnclusteredDown") {
          phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgs.at(1).Add(&phiCPeven_CMS_scale_met_unclustered_13TeVDown_hgsAR.at(i),-norm);
          phiCPeven_CMS_scale_met_unclustered_13TeVDown_ztt.at(1).Add(&phiCPeven_CMS_scale_met_unclustered_13TeVDown_zttAR.at(i),-norm);
          phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkj.at(1).Add(&phiCPeven_CMS_scale_met_unclustered_13TeVDown_fkjAR.at(i),-norm);
        }
	//
	if(Selection::Get_SysType() == "default") {
	  phiCPodd_nominal_hgs.at(1).Add(&phiCPodd_nominal_hgsAR.at(i),-norm);
	  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR.at(i),-norm);
	  phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR.at(i),-norm);
	  phiCPodd_CMS_eff_b_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_eff_b_13TeVUp_hgsAR.at(i),-norm);
	  phiCPodd_CMS_eff_b_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_eff_b_13TeVDown_hgsAR.at(i),-norm);
	  phiCPodd_CMS_htt_dyShape_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_htt_dyShape_13TeVUp_hgsAR.at(i),-norm);
	  phiCPodd_CMS_htt_dyShape_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_htt_dyShape_13TeVDown_hgsAR.at(i),-norm);
	  phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_htt_ttbarShape_13TeVUp_hgsAR.at(i),-norm);
	  phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_htt_ttbarShape_13TeVDown_hgsAR.at(i),-norm);
	  phiCPodd_CMS_scale_gg_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_gg_13TeVUp_hgsAR.at(i),-norm);
	  phiCPodd_CMS_scale_gg_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_gg_13TeVDown_hgsAR.at(i),-norm);
	  phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_hgsAR.at(i),-norm);
	  phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_hgsAR.at(i),-norm);
	  phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_hgsAR.at(i),-norm);
	  phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets0Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets0Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets1Up_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets1Down_hgs.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets0Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets0Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets1Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets1Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets2Up_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets2Up_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets2Down_hgs.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets2Down_hgsAR.at(i),-norm);
          phiCPodd_ff_mt_sub_systUp_hgs.at(1).Add(&phiCPodd_ff_mt_sub_systUp_hgsAR.at(i),-norm*1.1);
          phiCPodd_ff_mt_sub_systDown_hgs.at(1).Add(&phiCPodd_ff_mt_sub_systDown_hgsAR.at(i),-norm*0.9);
          phiCPodd_CMS_ttbar_embeded_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_ttbar_embeded_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_ttbar_embeded_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_ttbar_embeded_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_PreFire_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_PreFire_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_PreFire_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_PreFire_13TeVDown_hgsAR.at(i),-norm);
          //
          phiCPodd_nominal_ztt.at(1).Add(&phiCPodd_nominal_zttAR.at(i),-norm);
          phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_eff_b_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_eff_b_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_eff_b_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_eff_b_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_dyShape_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_htt_dyShape_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_dyShape_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_htt_dyShape_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_ttbarShape_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_htt_ttbarShape_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_ttbarShape_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_htt_ttbarShape_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_gg_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_gg_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_gg_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_gg_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_PS_ISR_ggH_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_PS_ISR_ggH_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_PS_FSR_ggH_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_PS_FSR_ggH_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets0Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets0Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets0Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets0Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets1Up_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets1Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets1Down_ztt.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets1Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets0Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets0Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets0Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets0Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets1Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets1Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets1Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets1Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets2Up_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets2Up_zttAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets2Down_ztt.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets2Down_zttAR.at(i),-norm);
          phiCPodd_ff_mt_sub_systUp_ztt.at(1).Add(&phiCPodd_ff_mt_sub_systUp_zttAR.at(i),-norm*1.1);
          phiCPodd_ff_mt_sub_systDown_ztt.at(1).Add(&phiCPodd_ff_mt_sub_systDown_zttAR.at(i),-norm*0.9);
          phiCPodd_CMS_ttbar_embeded_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_ttbar_embeded_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_ttbar_embeded_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_ttbar_embeded_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_PreFire_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_PreFire_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_PreFire_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_PreFire_13TeVDown_zttAR.at(i),-norm);
          //
          phiCPodd_nominal_fkj.at(1).Add(&phiCPodd_nominal_fkjAR.at(i),-norm);
          phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_CMS_eff_b_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_eff_b_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_eff_b_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_eff_b_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_CMS_htt_dyShape_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_htt_dyShape_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_htt_dyShape_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_htt_dyShape_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_htt_ttbarShape_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_htt_ttbarShape_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_CMS_scale_gg_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_gg_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_scale_gg_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_gg_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_PS_ISR_ggH_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_PS_ISR_ggH_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_PS_FSR_ggH_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_PS_FSR_ggH_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets0Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets0Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets1Up_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_qcd_syst_njets1Down_fkj.at(1).Add(&phiCPodd_ff_mt_qcd_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets0Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets0Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets1Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets1Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets2Up_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_wjets_syst_njets2Down_fkj.at(1).Add(&phiCPodd_ff_mt_wjets_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPodd_ff_mt_sub_systUp_fkj.at(1).Add(&phiCPodd_ff_mt_sub_systUp_fkjAR.at(i),-norm*1.1);
          phiCPodd_ff_mt_sub_systDown_fkj.at(1).Add(&phiCPodd_ff_mt_sub_systDown_fkjAR.at(i),-norm*0.9);
          phiCPodd_CMS_ttbar_embeded_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_ttbar_embeded_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_ttbar_embeded_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_ttbar_embeded_13TeVDown_fkjAR.at(i),-norm);
          phiCPodd_CMS_PreFire_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_PreFire_13TeVUp_fkjAR.at(i),-norm);
          phiCPodd_CMS_PreFire_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_PreFire_13TeVDown_fkjAR.at(i),-norm);
	}
        if(Selection::Get_SysType() == "TESUp") {
          phiCPodd_CMS_scale_t_3prong_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_t_3prong_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_t_3prong_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_t_3prong_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_t_3prong_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_t_3prong_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "TESDown") {
          phiCPodd_CMS_scale_t_3prong_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_t_3prong_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_t_3prong_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_t_3prong_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_t_3prong_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_t_3prong_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "MESUp") {
          phiCPodd_CMS_scale_mu_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_mu_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_mu_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_mu_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_mu_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_mu_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "MESDown") {
          phiCPodd_CMS_scale_mu_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_mu_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_mu_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_mu_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_mu_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_mu_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "FlavorQCDUp") {
          phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_FlavorQCD_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_FlavorQCD_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "FlavorQCDDown") {
          phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_FlavorQCD_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_FlavorQCD_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "RelativeBalUp") {
          phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeBal_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "RelativeBalDown") {
          phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeBal_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HFUp") {
          phiCPodd_CMS_scale_j_HF_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_HF_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_HF_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HFDown") {
          phiCPodd_CMS_scale_j_HF_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_HF_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_HF_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2016_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2017_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2018_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2016_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2017_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2018_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_HF_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "BBEC1Up") {
          phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_BBEC1_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "BBEC1Down") {
          phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_BBEC1_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "EC2Up") {
          phiCPodd_CMS_scale_j_EC2_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_EC2_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_EC2_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "EC2Down") {
          phiCPodd_CMS_scale_j_EC2_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_EC2_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_EC2_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_EC2_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "AbsoluteUp") {
          phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_Absolute_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "AbsoluteDown") {
          phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_j_Absolute_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) {
	    phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "JERUp") {
          phiCPodd_CMS_res_j_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_res_j_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_res_j_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_res_j_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_res_j_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_res_j_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "JERDown") {
          phiCPodd_CMS_res_j_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_res_j_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_res_j_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_res_j_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_res_j_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_res_j_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METResoUp") {
          phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_reso_met_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_htt_boson_reso_met_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METResoDown") {
          phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_reso_met_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_htt_boson_reso_met_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METScaleUp") {
          phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_scale_met_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_htt_boson_scale_met_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METScaleDown") {
          phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_scale_met_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_htt_boson_scale_met_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METUnclusteredUp") {
          phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgs.at(1).Add(&phiCPodd_CMS_scale_met_unclustered_13TeVUp_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_met_unclustered_13TeVUp_ztt.at(1).Add(&phiCPodd_CMS_scale_met_unclustered_13TeVUp_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkj.at(1).Add(&phiCPodd_CMS_scale_met_unclustered_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METUnclusteredDown") {
          phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgs.at(1).Add(&phiCPodd_CMS_scale_met_unclustered_13TeVDown_hgsAR.at(i),-norm);
          phiCPodd_CMS_scale_met_unclustered_13TeVDown_ztt.at(1).Add(&phiCPodd_CMS_scale_met_unclustered_13TeVDown_zttAR.at(i),-norm);
          phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkj.at(1).Add(&phiCPodd_CMS_scale_met_unclustered_13TeVDown_fkjAR.at(i),-norm);
        }
	//
	if(Selection::Get_SysType() == "default") {
	  phiCPMM_nominal_hgs.at(1).Add(&phiCPMM_nominal_hgsAR.at(i),-norm);
	  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_hgsAR.at(i),-norm);
	  phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_hgsAR.at(i),-norm);
	  phiCPMM_CMS_eff_b_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_eff_b_13TeVUp_hgsAR.at(i),-norm);
	  phiCPMM_CMS_eff_b_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_eff_b_13TeVDown_hgsAR.at(i),-norm);
	  phiCPMM_CMS_htt_dyShape_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_htt_dyShape_13TeVUp_hgsAR.at(i),-norm);
	  phiCPMM_CMS_htt_dyShape_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_htt_dyShape_13TeVDown_hgsAR.at(i),-norm);
	  phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_htt_ttbarShape_13TeVUp_hgsAR.at(i),-norm);
	  phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_htt_ttbarShape_13TeVDown_hgsAR.at(i),-norm);
	  phiCPMM_CMS_scale_gg_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_gg_13TeVUp_hgsAR.at(i),-norm);
	  phiCPMM_CMS_scale_gg_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_gg_13TeVDown_hgsAR.at(i),-norm);
	  phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_hgsAR.at(i),-norm);
	  phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_hgsAR.at(i),-norm);
	  phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_hgsAR.at(i),-norm);
	  phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets0Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets0Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets1Up_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets1Down_hgs.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets0Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets0Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets0Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets0Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets1Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets1Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets1Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets1Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets2Up_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets2Up_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets2Down_hgs.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets2Down_hgsAR.at(i),-norm);
          phiCPMM_ff_mt_sub_systUp_hgs.at(1).Add(&phiCPMM_ff_mt_sub_systUp_hgsAR.at(i),-norm*1.1);
          phiCPMM_ff_mt_sub_systDown_hgs.at(1).Add(&phiCPMM_ff_mt_sub_systDown_hgsAR.at(i),-norm*0.9);
          phiCPMM_CMS_ttbar_embeded_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_ttbar_embeded_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_ttbar_embeded_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_ttbar_embeded_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_PreFire_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_PreFire_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_PreFire_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_PreFire_13TeVDown_hgsAR.at(i),-norm);
          //
          phiCPMM_nominal_ztt.at(1).Add(&phiCPMM_nominal_zttAR.at(i),-norm);
          phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_eff_b_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_eff_b_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_eff_b_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_eff_b_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_dyShape_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_htt_dyShape_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_dyShape_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_htt_dyShape_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_ttbarShape_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_htt_ttbarShape_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_ttbarShape_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_htt_ttbarShape_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_gg_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_gg_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_gg_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_gg_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_PS_ISR_ggH_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_PS_ISR_ggH_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_PS_FSR_ggH_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_PS_FSR_ggH_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets0Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets0Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets0Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets0Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets1Up_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets1Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets1Down_ztt.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets1Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets0Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets0Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets0Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets0Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets1Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets1Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets1Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets1Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets2Up_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets2Up_zttAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets2Down_ztt.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets2Down_zttAR.at(i),-norm);
          phiCPMM_ff_mt_sub_systUp_ztt.at(1).Add(&phiCPMM_ff_mt_sub_systUp_zttAR.at(i),-norm*1.1);
          phiCPMM_ff_mt_sub_systDown_ztt.at(1).Add(&phiCPMM_ff_mt_sub_systDown_zttAR.at(i),-norm*0.9);
          phiCPMM_CMS_ttbar_embeded_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_ttbar_embeded_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_ttbar_embeded_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_ttbar_embeded_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_PreFire_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_PreFire_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_PreFire_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_PreFire_13TeVDown_zttAR.at(i),-norm);
          //
          phiCPMM_nominal_fkj.at(1).Add(&phiCPMM_nominal_fkjAR.at(i),-norm);
          phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_eff_Xtrigger_mt_MVADM10_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_CMS_eff_b_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_eff_b_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_eff_b_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_eff_b_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_CMS_htt_dyShape_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_htt_dyShape_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_htt_dyShape_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_htt_dyShape_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_htt_ttbarShape_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_htt_ttbarShape_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_CMS_scale_gg_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_gg_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_scale_gg_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_gg_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_PS_ISR_ggH_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_PS_ISR_ggH_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_PS_FSR_ggH_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_PS_FSR_ggH_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc1_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_stat_unc2_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc1_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets0_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets1_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_stat_unc2_njets2_mvadm10Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets0Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets0Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_met_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_l_pt_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets1Up_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_qcd_syst_njets1Down_fkj.at(1).Add(&phiCPMM_ff_mt_qcd_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets0Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets0Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets0Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets0Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets1Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets1Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets1Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets1Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_met_closure_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_l_pt_closure_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets2Up_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets2Up_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_wjets_syst_njets2Down_fkj.at(1).Add(&phiCPMM_ff_mt_wjets_syst_njets2Down_fkjAR.at(i),-norm);
          phiCPMM_ff_mt_sub_systUp_fkj.at(1).Add(&phiCPMM_ff_mt_sub_systUp_fkjAR.at(i),-norm*1.1);
          phiCPMM_ff_mt_sub_systDown_fkj.at(1).Add(&phiCPMM_ff_mt_sub_systDown_fkjAR.at(i),-norm*0.9);
          phiCPMM_CMS_ttbar_embeded_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_ttbar_embeded_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_ttbar_embeded_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_ttbar_embeded_13TeVDown_fkjAR.at(i),-norm);
          phiCPMM_CMS_PreFire_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_PreFire_13TeVUp_fkjAR.at(i),-norm);
          phiCPMM_CMS_PreFire_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_PreFire_13TeVDown_fkjAR.at(i),-norm);
	}
        if(Selection::Get_SysType() == "TESUp") {
          phiCPMM_CMS_scale_t_3prong_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_t_3prong_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_t_3prong_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_t_3prong_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_t_3prong_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_t_3prong_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "TESDown") {
          phiCPMM_CMS_scale_t_3prong_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_t_3prong_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_t_3prong_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_t_3prong_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_t_3prong_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_t_3prong_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "MESUp") {
          phiCPMM_CMS_scale_mu_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_mu_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_mu_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_mu_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_mu_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_mu_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "MESDown") {
          phiCPMM_CMS_scale_mu_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_mu_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_mu_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_mu_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_mu_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_mu_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "FlavorQCDUp") {
          phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_FlavorQCD_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_FlavorQCD_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "FlavorQCDDown") {
          phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_FlavorQCD_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_FlavorQCD_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "RelativeBalUp") {
          phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeBal_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "RelativeBalDown") {
          phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeBal_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HFUp") {
          phiCPMM_CMS_scale_j_HF_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_HF_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_HF_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HFDown") {
          phiCPMM_CMS_scale_j_HF_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_HF_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_HF_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "HF_YEARUp") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2016_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2017_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2018_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "HF_YEARDown") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2016_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2017_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2018_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_HF_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "BBEC1Up") {
          phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_BBEC1_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "BBEC1Down") {
          phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_BBEC1_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "BBEC1_YEARUp") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "BBEC1_YEARDown") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_BBEC1_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "EC2Up") {
          phiCPMM_CMS_scale_j_EC2_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_EC2_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_EC2_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "EC2Down") {
          phiCPMM_CMS_scale_j_EC2_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_EC2_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_EC2_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "EC2_YEARUp") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "EC2_YEARDown") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_EC2_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "AbsoluteUp") {
          phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_Absolute_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "AbsoluteDown") {
          phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_j_Absolute_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "Absolute_YEARUp") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "Absolute_YEARDown") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_Absolute_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "RelativeSample_YEARUp") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVUp_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVUp_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "RelativeSample_YEARDown") {
	  if(theYear == 2016) {
	    phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2016_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2017) {
	    phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2017_13TeVDown_fkjAR.at(i),-norm);
	  }
	  if(theYear == 2018) {
	    phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_hgsAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_zttAR.at(i),-norm);
	    phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_j_RelativeSample_2018_13TeVDown_fkjAR.at(i),-norm);
	  }
	}
        if(Selection::Get_SysType() == "JERUp") {
          phiCPMM_CMS_res_j_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_res_j_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_res_j_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_res_j_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_res_j_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_res_j_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "JERDown") {
          phiCPMM_CMS_res_j_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_res_j_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_res_j_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_res_j_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_res_j_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_res_j_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METResoUp") {
          phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_reso_met_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_htt_boson_reso_met_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METResoDown") {
          phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_reso_met_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_htt_boson_reso_met_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METScaleUp") {
          phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_scale_met_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_htt_boson_scale_met_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METScaleDown") {
          phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_scale_met_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_htt_boson_scale_met_13TeVDown_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METUnclusteredUp") {
          phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgs.at(1).Add(&phiCPMM_CMS_scale_met_unclustered_13TeVUp_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_met_unclustered_13TeVUp_ztt.at(1).Add(&phiCPMM_CMS_scale_met_unclustered_13TeVUp_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkj.at(1).Add(&phiCPMM_CMS_scale_met_unclustered_13TeVUp_fkjAR.at(i),-norm);
        }
        if(Selection::Get_SysType() == "METUnclusteredDown") {
          phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgs.at(1).Add(&phiCPMM_CMS_scale_met_unclustered_13TeVDown_hgsAR.at(i),-norm);
          phiCPMM_CMS_scale_met_unclustered_13TeVDown_ztt.at(1).Add(&phiCPMM_CMS_scale_met_unclustered_13TeVDown_zttAR.at(i),-norm);
          phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkj.at(1).Add(&phiCPMM_CMS_scale_met_unclustered_13TeVDown_fkjAR.at(i),-norm);
        }
        std::cout << "end " << std::endl;
      }
    }
  }
  Selection::Finish(Channel);
}
