#include "HCPMuTau.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
//#include "SVfitProvider.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "SkimConfig.h"
#include "TauSpinerInterface.h"


#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TMath.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
//#include "Objects.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauPolSoftware/TauDecaysInterface/interface/fonction_a1.h"
#include "TauPolSoftware/TauDecaysInterface/interface/SCalculator.h"
#include <algorithm>

HCPMuTau::HCPMuTau(TString Name_, TString id_, char* Channel_, char* CPstate_):
  Selection(Name_,id_)
  //DataMC_Corr(true,true,false),
  //tauTrgSF("tight")
{
  Channel = Channel_;
  CPstate = CPstate_; 
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;

  TFile *f_fracs=TFile::Open("/opt/sbg/cms/safe1/cms/msessini/IPHCAnalysisTools/Code/CommonFiles/FakeFactors/mva_fract_mt_2018.root", "READ");
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

  _FF = new FakeFactors(2018, ff_fracs_qcd_, ff_fracs_wjets_, ff_fracs_qcd_ss_, ff_fracs_wjets_ss_, ff_fracs_qcd_aiso_, ff_fracs_wjets_aiso_, ff_fracs_qcd_highmt_, ff_fracs_wjets_highmt_, ff_ws_);

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

  polarimetricAcopAngle=HConfig.GetTH1D(Name+"_polarimetricAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove

}

void  HCPMuTau::Store_ExtraDist(){

  Extradist1d.push_back(&polarimetricAcopAngle);
}

void  HCPMuTau::doEvent()  { //  Method called on every event
 unsigned int t;
 int id(Ntp->GetMCID());  //read event ID of a sample

 int TauGenMatch = Ntp->tauGenMatch(), MuGenMatch = Ntp->muGenMatch();
 bool GenMatchSelection = false;
 if(Ntp->Id() == 35) GenMatchSelection = (TauGenMatch == 5 && MuGenMatch == 4); //Embedded mu+tau
 else if(Ntp->isZ() || Ntp->Id() == 212) GenMatchSelection = true;

 

 if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide

 polarimetricAcopAngle.at(t).Fill(Ntp->pvPhiCP(), Ntp->wEven());


  TLorentzVector taup4, mup4, metp4;
  taup4.SetPtEtaPhiE(Ntp->tauPt(), Ntp->tauEta(), Ntp->tauPhi(), Ntp->tauE());
  mup4.SetPtEtaPhiE(Ntp->muPt(), Ntp->muEta(), Ntp->muPhi(), Ntp->muE());
  metp4.SetPtEtaPhiE(Ntp->PUPPImet(), 0., Ntp->PUPPImetphi(), 0.);

 _FF->Initialize(taup4, mup4, metp4, Ntp->tauDM(), Ntp->Njets(), Ntp->dijetMass(), Ntp->muMETmt(), Ntp->muIso, Ntp->tauIPsignificance(), Ntp->isOSpair(), Ntp->isIso());
 std::map<std::string, double> FFmap = _FF->GetFakeFactors(Selection::Get_SysType());

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
  Selection::Finish(Channel,CPstate);
}
