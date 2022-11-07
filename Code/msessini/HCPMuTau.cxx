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

  /*TFile *f_fracs=TFile::Open("/opt/sbg/cms/safe1/cms/msessini/IPHCAnalysisTools/Code/CommonFiles/FakeFactors/mva_fract_mt_2018.root", "READ");
  ff_fracs_qcd_ = (TH2D*)f_fracs->Get("QCD");
  ff_fracs_wjets_ = (TH2D*)f_fracs->Get("W");
  ff_fracs_qcd_->SetDirectory(0);
  ff_fracs_wjets_->SetDirectory(0);
  f_fracs->Close();

  TFile *f_fracs_ss=TFile::Open("/opt/sbg/cms/safe1/cms/msessini/IPHCAnalysisTools/Code/CommonFiles/FakeFactors/mva_fract_mt_2018_ss.root", "READ");
  ff_fracs_qcd_ss_ = (TH2D*)f_fracs_ss->Get("QCD");
  ff_fracs_wjets_ss_ = (TH2D*)f_fracs_ss->Get("W");
  ff_fracs_qcd_ss_->SetDirectory(0);
  ff_fracs_wjets_ss_->SetDirectory(0);
  f_fracs_ss->Close();

  TFile *f_fracs_aiso=TFile::Open("/opt/sbg/cms/safe1/cms/msessini/IPHCAnalysisTools/Code/CommonFiles/FakeFactors/mva_fract_mt_2018_aiso.root", "READ");
  ff_fracs_qcd_aiso_ = (TH2D*)f_fracs_aiso->Get("QCD");
  ff_fracs_wjets_aiso_ = (TH2D*)f_fracs_aiso->Get("W");
  ff_fracs_qcd_aiso_->SetDirectory(0);
  ff_fracs_wjets_aiso_->SetDirectory(0);
  f_fracs_aiso->Close();

  TFile *f_fracs_highmt=TFile::Open("/opt/sbg/cms/safe1/cms/msessini/IPHCAnalysisTools/Code/CommonFiles/FakeFactors/mva_fract_mt_2018_highmt.root", "READ");
  ff_fracs_qcd_highmt_ = (TH2D*)f_fracs_highmt->Get("QCD");
  ff_fracs_wjets_highmt_ = (TH2D*)f_fracs_highmt->Get("W");
  ff_fracs_qcd_highmt_->SetDirectory(0);
  ff_fracs_wjets_highmt_->SetDirectory(0);
  f_fracs_highmt->Close();

  TFile *f=TFile::Open("/opt/sbg/cms/safe1/cms/msessini/IPHCAnalysisTools/Code/CommonFiles/FakeFactors/fakefactors_ws_mt_lite_2018.root", "READ");
  ff_ws_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  f->Close();

  _FF = new FakeFactors(2018, ff_fracs_qcd_, ff_fracs_wjets_, ff_fracs_qcd_ss_, ff_fracs_wjets_ss_, ff_fracs_qcd_aiso_, ff_fracs_wjets_aiso_, ff_fracs_qcd_highmt_, ff_fracs_wjets_highmt_, ff_ws_);*/

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
    if(i==Trigger)             cut.at(Trigger)=1;
    if(i==Id_and_Kin)            cut.at(Id_and_Kin)=1;
    //if(i==NPairsFound)         cut.at(NPairsFound)=1;
    //if(i==GoodIndex)           cut.at(GoodIndex)=1.;
    //if(i==ZTTMC)                 cut.at(ZTTMC)=1.;
    //if(i==METFilters)            cut.at(METFilters)=1.;
    //if(i==genmatch)              cut.at(genmatch)=1;
    if(i==TausIsolation)         cut.at(TausIsolation)=1;
    if(i==AgainstEleMu)          cut.at(AgainstEleMu)=1;
    //if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    if(i==LeptonVeto)            cut.at(LeptonVeto)=0;
    if(i==PairCharge)            cut.at(PairCharge)=1.;
    if(i==PairMass)              cut.at(PairMass)=40.;
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
    if(i==Trigger){
      title.at(i)="Trigger Matching";
      hlabel="At least 1 good pair with Trig+Matching";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
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
    else if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==genmatch){
    //   title.at(i)="genmatch";
    //   hlabel="genmatch";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==TausIsolation){
      title.at(i)="Taus Isolation";
      hlabel="Isolation of Taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==AgainstEleMu){
      title.at(i)="Against Electrons and Muons";
      hlabel="Against Electrons and Muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==Tau2Isolation){
    //   title.at(i)="Tau2 Isolation";
    //   hlabel="Isolation of Tau2";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==LeptonVeto){
      title.at(i)="Third Lepton Veto";
      hlabel="Third Lepton Veto";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PairMass){
      title.at(i)="Pair Visible Mass";
      hlabel="M(tau-tau)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,30,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,30,0,150,hlabel,"Events"));
    }

  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  polarimetricAcopAngleEven=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEven"," ",5,0.,2*TMath::Pi()," "," ");
  polarimetricAcopAngleOdd=HConfig.GetTH1D(Name+"_polarimetricAcopAngleOdd"," ",5,0.,2*TMath::Pi()," "," ");
  polarimetricAcopAngleMM=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMM"," ",5,0.,2*TMath::Pi()," "," ");
  AcopAngleEven=HConfig.GetTH1D(Name+"_AcopAngleEven"," ",5,0.,2*TMath::Pi()," "," ");
  AcopAngleOdd=HConfig.GetTH1D(Name+"_AcopAngleOdd"," ",5,0.,2*TMath::Pi()," "," ");
  AcopAngleMM=HConfig.GetTH1D(Name+"_AcopAngleMM"," ",5,0.,2*TMath::Pi()," "," ");
  genpolarimetricAcopAngleEven=HConfig.GetTH1D(Name+"_genpolarimetricAcopAngleEven"," ",5,0.,2*TMath::Pi()," "," ");
  genpolarimetricAcopAngleOdd=HConfig.GetTH1D(Name+"_genpolarimetricAcopAngleOdd"," ",5,0.,2*TMath::Pi()," "," ");
  genpolarimetricAcopAngleMM=HConfig.GetTH1D(Name+"_genpolarimetricAcopAngleMM"," ",5,0.,2*TMath::Pi()," "," ");
  genAcopAngleEven=HConfig.GetTH1D(Name+"_genAcopAngleEven"," ",5,0.,2*TMath::Pi()," "," ");
  genAcopAngleOdd=HConfig.GetTH1D(Name+"_genAcopAngleOdd"," ",5,0.,2*TMath::Pi()," "," ");
  genAcopAngleMM=HConfig.GetTH1D(Name+"_genAcopAngleMM"," ",5,0.,2*TMath::Pi()," "," ");
  pullPVx=HConfig.GetTH1D(Name+"_pullPVx"," ",50,-1,1," "," ");
  pullPVy=HConfig.GetTH1D(Name+"_pullPVy"," ",50,-1,1," "," ");
  pullPVz=HConfig.GetTH1D(Name+"_pullPVz"," ",50,-1,1," "," ");
  pullTauSVx=HConfig.GetTH1D(Name+"_pullTausVx"," ",50,-1,1," "," ");
  pullTauSVy=HConfig.GetTH1D(Name+"_pullTauSVy"," ",50,-1,1," "," ");
  pullTauSVz=HConfig.GetTH1D(Name+"_pullTauSVz"," ",50,-1,1," "," ");
  pullTauE=HConfig.GetTH1D(Name+"_pullTauE"," ",50,-1,1," "," ");
  pullTauPt=HConfig.GetTH1D(Name+"_pullTauPt"," ",50,-1,1," "," ");
  pullTauPhi=HConfig.GetTH1D(Name+"_pullTauPhi"," ",50,-1,1," "," ");
  pullTauEta=HConfig.GetTH1D(Name+"_pullTauEta"," ",50,-1,1," "," ");
  pullMuonE=HConfig.GetTH1D(Name+"_pullMuonE"," ",50,-1,1," "," ");
  pullMuonPt=HConfig.GetTH1D(Name+"_pullMuonPt"," ",50,-1,1," "," ");
  pullMuonPhi=HConfig.GetTH1D(Name+"_pullMuonPhi"," ",50,-1,1," "," ");
  pullMuonEta=HConfig.GetTH1D(Name+"_pullMuonEta"," ",50,-1,1," "," ");

  //AcopAngle uncategorized
  /*polarimetricAcopAngleEven=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEven","PhiCP",60,0.,2*TMath::Pi(),"PhiCP","Events");
  polarimetricAcopAngleMM=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMM","PhiCP",60,0.,2*TMath::Pi(),"PhiCP","Events");
  polarimetricAcopAngleOdd=HConfig.GetTH1D(Name+"_polarimetricAcopAngleOdd","PhiCP",60,0.,2*TMath::Pi(),"PhiCP","Events");
  //
  decayplaneAcopAngleEven=HConfig.GetTH1D(Name+"_decayplaneAcopAngleEven","PhiCP",60,0.,2*TMath::Pi(),"PhiCP","Events");
  decayplaneAcopAngleMM=HConfig.GetTH1D(Name+"_decayplaneAcopAngleMM","PhiCP",60,0.,2*TMath::Pi(),"PhiCP","Events");
  decayplaneAcopAngleOdd=HConfig.GetTH1D(Name+"_decayplaneAcopAngleOdd","PhiCP",60,0.,2*TMath::Pi(),"PhiCP","Events");
  //AcopAngle uncategorized QCD
  polarimetricAcopAngleEvenQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEvenQCDMC","QCDMC PhiCP",60,0.,2*TMath::Pi(),"PhiCP","Events");
  //AcopAngle categorized
  polarimetricAcopAngleEvenHiggs=HConfig.GetTH2D(Name+"_polarimetricAcopAngleEvenHiggs","Higgs BDT Score VS phiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","BDT score");
  polarimetricAcopAngleEvenJetFakes=HConfig.GetTH2D(Name+"_polarimetricAcopAngleEvenJetFakes","JetFakes BDT Score VS phiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","BDT score");
  polarimetricAcopAngleEvenZTT=HConfig.GetTH2D(Name+"_polarimetricAcopAngleEvenZTT","ZTT BDT Score VS phiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","BDT score");
  //AcopAngle categorized QCD
  polarimetricAcopAngleEvenHiggsQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAngleEvenHiggsQCDMC","Higgs BDT Score VS phiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","BDT score");
  polarimetricAcopAngleEvenJetFakesQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAngleEvenJetFakesQCDMC","JetFakes BDT Score VS phiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","BDT score");
  polarimetricAcopAngleEvenZTTQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAngleEvenZTTQCDMC","ZTT BDT Score VS phiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","BDT score");
  //AcopAngle categorized unrolled
  polarimetricAcopAngleEvenHiggsUnrolled=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEvenHiggsUnrolled","Higgs unrolled PhiCP Even",180,0.,3*2*TMath::Pi(),"PhiCP","Events");
  polarimetricAcopAngleEvenJetFakesUnrolled=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEvenJetFakesUnrolled","JetFakes unrolled PhiCP Even",180,0.,3*2*TMath::Pi(),"PhiCP","Events");
  polarimetricAcopAngleEvenZTTUnrolled=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEvenZTTUnrolled","ZTT unrolled PhiCP",180,0.,3*2*TMath::Pi(),"PhiCP","Events");
  //AcopAngle categorized unrolled QCD
  polarimetricAcopAngleEvenHiggsUnrolledQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEvenHiggsUnrolledQCDMC","Higgs unrolled PhiCP Even",180,0.,3*2*TMath::Pi(),"PhiCP","Events");
  polarimetricAcopAngleEvenJetFakesUnrolledQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEvenJetFakesUnrolledQCDMC","JetFakes unrolled PhiCP Even",180,0.,3*2*TMath::Pi(),"PhiCP","Events");
  polarimetricAcopAngleEvenZTTUnrolledQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleEvenZTTUnrolledQCDMC","ZTT unrolled PhiCP",180,0.,3*2*TMath::Pi(),"PhiCP","Events");
  //Wfakes
  polarimetricWfakesHiggs=HConfig.GetTH2D(Name+"_polarimetricWfakesHiggs","WfakesHiggs",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  polarimetricWfakesJetFakes=HConfig.GetTH2D(Name+"_polarimetricWfakesJetFakes","WfakesJetFakes",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  polarimetricWfakesZTT=HConfig.GetTH2D(Name+"_polarimetricWfakesZTT","WfakesZTT",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  //BDT scores
  HiggsBDTScore=HConfig.GetTH1D(Name+"_HiggsBDTScore","HiggsBDTScore",7,0.3,1,"BDT Score","Events");
  JetFakesBDTScore=HConfig.GetTH1D(Name+"_JetFakesBDTScore","JetFakesBDTScore",7,0.3,1,"BDT Score","Events");
  ZTTBDTScore=HConfig.GetTH1D(Name+"_ZTTBDTScore","ZTTBDTScore",7,0.3,1,"BDT Score","Events");
  //BDT scores QCD
  HiggsBDTScoreQCDMC=HConfig.GetTH1D(Name+"_HiggsBDTScoreQCDMC","HiggsBDTScore",7,0.3,1,"BDT Score","Events");
  JetFakesBDTScoreQCDMC=HConfig.GetTH1D(Name+"_JetFakesBDTScoreQCDMC","JetFakesBDTScore",7,0.3,1,"BDT Score","Events");
  ZTTBDTScoreQCDMC=HConfig.GetTH1D(Name+"_ZTTBDTScoreQCDMC","ZTTBDTScore",7,0.3,1,"BDT Score","Events");
  //
  ShapeSystHiggs=HConfig.GetTH2D(Name+"_ShapeSystHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  //trg eff high pT
  polarimetric_trgeff_highpT_mvadm10_UpEvenHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_highpT_mvadm10_UpEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_highpT_mvadm10_UpMMHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_highpT_mvadm10_UpMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_highpT_mvadm10_UpOddHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_highpT_mvadm10_UpOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  //
  decayplane_trgeff_highpT_mvadm10_UpEvenHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_highpT_mvadm10_UpEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_highpT_mvadm10_UpMMHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_highpT_mvadm10_UpMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_highpT_mvadm10_UpOddHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_highpT_mvadm10_UpOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  ////
  polarimetric_trgeff_highpT_mvadm10_DownEvenHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_highpT_mvadm10_DownEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_highpT_mvadm10_DownMMHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_highpT_mvadm10_DownMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_highpT_mvadm10_DownOddHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_highpT_mvadm10_DownOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  //
  decayplane_trgeff_highpT_mvadm10_DownEvenHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_highpT_mvadm10_DownEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_highpT_mvadm10_DownMMHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_highpT_mvadm10_DownMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_highpT_mvadm10_DownOddHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_highpT_mvadm10_DownOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  //trg eff
  polarimetric_trgeff_mvadm10_UpEvenHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_mvadm10_UpEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_mvadm10_UpMMHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_mvadm10_UpMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_mvadm10_UpOddHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_mvadm10_UpOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  //
  decayplane_trgeff_mvadm10_UpEvenHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_mvadm10_UpEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_mvadm10_UpMMHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_mvadm10_UpMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_mvadm10_UpOddHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_mvadm10_UpOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  ////
  polarimetric_trgeff_mvadm10_DownEvenHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_mvadm10_DownEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_mvadm10_DownMMHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_mvadm10_DownMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  polarimetric_trgeff_mvadm10_DownOddHiggs=HConfig.GetTH2D(Name+"_polarimetric_trgeff_mvadm10_DownOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  //
  decayplane_trgeff_mvadm10_DownEvenHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_mvadm10_DownEvenHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_mvadm10_DownMMHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_mvadm10_DownMMHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  decayplane_trgeff_mvadm10_DownOddHiggs=HConfig.GetTH2D(Name+"_decayplane_trgeff_mvadm10_DownOddHiggs","PhiCP",60,0.,2*TMath::Pi(),7,0.3,1,"PhiCP","Events");
  //
  */
  

  


  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove

}

void  HCPMuTau::Store_ExtraDist(){

  Extradist1d.push_back(&polarimetricAcopAngleEven);
  Extradist1d.push_back(&polarimetricAcopAngleOdd);
  Extradist1d.push_back(&polarimetricAcopAngleMM);
  Extradist1d.push_back(&AcopAngleEven);
  Extradist1d.push_back(&AcopAngleOdd);
  Extradist1d.push_back(&AcopAngleMM);
  Extradist1d.push_back(&genpolarimetricAcopAngleEven);
  Extradist1d.push_back(&genpolarimetricAcopAngleOdd);
  Extradist1d.push_back(&genpolarimetricAcopAngleMM);
  Extradist1d.push_back(&genAcopAngleEven);
  Extradist1d.push_back(&genAcopAngleOdd);
  Extradist1d.push_back(&genAcopAngleMM);
  Extradist1d.push_back(&pullPVx);
  Extradist1d.push_back(&pullPVy);
  Extradist1d.push_back(&pullPVz);
  Extradist1d.push_back(&pullTauSVx);
  Extradist1d.push_back(&pullTauSVy);
  Extradist1d.push_back(&pullTauSVz);
  Extradist1d.push_back(&pullTauE);
  Extradist1d.push_back(&pullTauPt);
  Extradist1d.push_back(&pullTauPhi);
  Extradist1d.push_back(&pullTauEta);
  Extradist1d.push_back(&pullMuonE);
  Extradist1d.push_back(&pullMuonPt);
  Extradist1d.push_back(&pullMuonPhi);
  Extradist1d.push_back(&pullMuonEta);

}

void  HCPMuTau::doEvent()  { //  Method called on every event
 unsigned int t;
 int id(Ntp->GetMCID());  //read event ID of a sample

 int TauGenMatch = Ntp->tauGenMatch(), MuGenMatch = Ntp->muGenMatch();
 bool GenMatchSelection = false;
 if(Ntp->Id() == 35) GenMatchSelection = (TauGenMatch == 5 && MuGenMatch == 4); //Embedded mu+tau
 else if(Ntp->isZ() || Ntp->Id() == 212) GenMatchSelection = true;

 

 if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide
 
 if(Ntp->isOSpair() && Ntp->isMediumID() && Ntp->isIso()) {
 polarimetricAcopAngleEven.at(t).Fill(Ntp->pvPhiCP(), Ntp->wEven());
 polarimetricAcopAngleOdd.at(t).Fill(Ntp->pvPhiCP(), Ntp->wOdd());
 polarimetricAcopAngleMM.at(t).Fill(Ntp->pvPhiCP(), Ntp->wMM());
 AcopAngleEven.at(t).Fill(Ntp->dpPhiCP(), Ntp->wEven());
 AcopAngleOdd.at(t).Fill(Ntp->dpPhiCP(), Ntp->wOdd());
 AcopAngleMM.at(t).Fill(Ntp->dpPhiCP(), Ntp->wMM());
 }
 TLorentzVector gentaup4(Ntp->genTaupx(), Ntp->genTaupy(), Ntp->genTaupz(), Ntp->genTauE());
 TLorentzVector genmuonp4(Ntp->genMuonpx(), Ntp->genMuonpy(), Ntp->genMuonpz(), Ntp->genMuonE());

 if(Selection::Get_SysType() == "default") {
   genpolarimetricAcopAngleEven.at(t).Fill(Ntp->genpvPhiCP(), Ntp->wEven());
   genpolarimetricAcopAngleOdd.at(t).Fill(Ntp->genpvPhiCP(), Ntp->wOdd());
   genpolarimetricAcopAngleMM.at(t).Fill(Ntp->genpvPhiCP(), Ntp->wMM());
   genAcopAngleEven.at(t).Fill(Ntp->gendpPhiCP(), Ntp->wEven());
   genAcopAngleOdd.at(t).Fill(Ntp->gendpPhiCP(), Ntp->wOdd());
   genAcopAngleMM.at(t).Fill(Ntp->gendpPhiCP(), Ntp->wMM());
   pullPVx.at(t).Fill((Ntp->genPVx() - Ntp->pvx())/Ntp->genPVx());
   pullPVy.at(t).Fill((Ntp->genPVy() - Ntp->pvy())/Ntp->genPVy());
   pullPVz.at(t).Fill((Ntp->genPVz() - Ntp->pvz())/Ntp->genPVz());
   pullTauSVx.at(t).Fill((Ntp->genTauSVx() - Ntp->tauSVx())/Ntp->genTauSVx());
   pullTauSVy.at(t).Fill((Ntp->genTauSVy() - Ntp->tauSVy())/Ntp->genTauSVy());
   pullTauSVz.at(t).Fill((Ntp->genTauSVz() - Ntp->tauSVz())/Ntp->genTauSVz());
   pullTauE.at(t).Fill((gentaup4.E() - Ntp->GEFtauE())/gentaup4.E());
   pullTauPt.at(t).Fill((gentaup4.Pt() - Ntp->GEFtauPt())/gentaup4.Pt());
   pullTauPhi.at(t).Fill((gentaup4.Phi() - Ntp->GEFtauPhi())/gentaup4.Phi());
   pullTauEta.at(t).Fill((gentaup4.Eta() - Ntp->GEFtauEta())/gentaup4.Eta());
   pullMuonE.at(t).Fill((genmuonp4.E() - Ntp->muE())/genmuonp4.E());
   pullMuonPt.at(t).Fill((genmuonp4.Pt() - Ntp->muPt())/genmuonp4.Pt());
   pullMuonPhi.at(t).Fill((genmuonp4.Phi() - Ntp->muPhi())/genmuonp4.Phi());
   pullMuonEta.at(t).Fill((genmuonp4.Eta() - Ntp->muEta())/genmuonp4.Eta());
 }

 TLorentzVector taup4, mup4, metp4;
 taup4.SetPtEtaPhiE(Ntp->tauPt(), Ntp->tauEta(), Ntp->tauPhi(), Ntp->tauE());
 mup4.SetPtEtaPhiE(Ntp->muPt(), Ntp->muEta(), Ntp->muPhi(), Ntp->muE());
 metp4.SetPtEtaPhiE(Ntp->PUPPImet(), 0., Ntp->PUPPImetphi(), 0.);

 //_FF->Initialize(taup4, mup4, metp4, Ntp->tauDM(), Ntp->Njets(), Ntp->dijetMass(), Ntp->muMETmt(), Ntp->muIso(), Ntp->tauIPsignificance(), Ntp->isOSpair(), Ntp->isIso());
 //std::map<std::string, double> FFmap = _FF->GetFakeFactors(Selection::Get_SysType());
 if(!Ntp->isData() && Ntp->isIso() && Ntp->isMediumID()) {
   

 }

} //do event

//  This is a function if you want to do something after the event loop
void HCPMuTau::Finish() {

  if(mode == RECONSTRUCT) {
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);

    double norm=1.;

    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(CrossSectionandAcceptance.at(i)>0 || HConfig.GetID(i)==36 || HConfig.GetID(i)==20 || HConfig.GetID(i)==23 || HConfig.GetID(i)==30 || HConfig.GetID(i)==33){
        if(CrossSectionandAcceptance.at(i)>0)norm= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
        else norm=1.;

        cout<<"Soustraction du QCD: "<<endl;
        std::cout << "end " << std::endl;
      }
    }
  }
  Selection::Finish(Channel);
}
