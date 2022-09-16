#include "HCPTauTau.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
//#include "SVfitProvider.h"
#include "TLorentzVector.h"
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

HCPTauTau::HCPTauTau(TString Name_, TString id_, char* Channel_, char* CPstate_):
  Selection(Name_,id_)
  //DataMC_Corr(true,true,false),
  //tauTrgSF("tight")
{
  Channel = Channel_;
  CPstate = CPstate_; 
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
  WorkSpaceFF2016=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2016_dR_corr.root").c_str(), "READ");
  wFF2016= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2016->Close();
  WorkSpaceFF2017=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2017_dR_corr.root").c_str(), "READ");
  wFF2017= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2017->Close();
  WorkSpaceFF2018=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2018_dR_corr.root").c_str(), "READ");
  wFF2018= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2018->Close();
  BDT=new BDTClassification();
  BDT->PreAnalysis();
}

HCPTauTau::~HCPTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
  delete wFF2016;
  delete wFF2017;
  delete wFF2018;
}

void  HCPTauTau::Configure(){
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
  polarimetricGEFAcopAngle=HConfig.GetTH1D(Name+"_polarimetricGEFAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  decayplaneAcopAngle=HConfig.GetTH1D(Name+"_decayplaneAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  impactparameterAcopAngle=HConfig.GetTH1D(Name+"_impactparameterAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  DPIPAcopAngle=HConfig.GetTH1D(Name+"_DPIPAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  PVIPAcopAngle=HConfig.GetTH1D(Name+"_PVIPAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");

  polarimetricAcopAngleTruth=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  decayplaneAcopAngleTruth=HConfig.GetTH1D(Name+"_decayplaneAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  impactparameterAcopAngleTruth=HConfig.GetTH1D(Name+"_impactparameterAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  DPIPAcopAngleTruth=HConfig.GetTH1D(Name+"_DPIPAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  PVIPAcopAngleTruth=HConfig.GetTH1D(Name+"_PVIPAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");

  // mtt
  DeltaPhitau1MTT=HConfig.GetTH1D(Name+"_DeltaPhitau1MTT","DeltaPhitau1MTT",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1MTT=HConfig.GetTH1D(Name+"_DeltaEtatau1MTT","DeltaEtatau1MTT",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1MTT=HConfig.GetTH1D(Name+"_DeltaPtau1MTT","DeltaPtau1MTT",50,-1,1,"#Delta P","Events");
  DeltaEtau1MTT=HConfig.GetTH1D(Name+"_DeltaEtau1MTT","DeltaEtau1MTT",50,-1,1,"#Delta E","Events");
  DeltaPhitau2MTT=HConfig.GetTH1D(Name+"_DeltaPhitau2MTT","DeltaPhitau2MTT",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2MTT=HConfig.GetTH1D(Name+"_DeltaEtatau2MTT","DeltaEtatau2MTT",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2MTT=HConfig.GetTH1D(Name+"_DeltaPtau2MTT","DeltaPtau2MTT",50,-1,1,"#Delta P","Events");
  DeltaEtau2MTT=HConfig.GetTH1D(Name+"_DeltaEtau2MTT","DeltaEtau2MTT",50,-1,1,"#Delta E","Events");
  //svfit
  DeltaPhitau1SVFit=HConfig.GetTH1D(Name+"_DeltaPhitau1SVFit","DeltaPhitau1SVFit",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1SVFit=HConfig.GetTH1D(Name+"_DeltaEtatau1SVFit","DeltaEtatau1SVFit",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1SVFit=HConfig.GetTH1D(Name+"_DeltaPtau1SVFit","DeltaPtau1SVFit",50,-1,1,"#Delta P","Events");
  DeltaEtau1SVFit=HConfig.GetTH1D(Name+"_DeltaEtau1SVFit","DeltaEtau1SVFit",50,-1,1,"#Delta E","Events");
  DeltaPhitau2SVFit=HConfig.GetTH1D(Name+"_DeltaPhitau2SVFit","DeltaPhitau2SVFit",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2SVFit=HConfig.GetTH1D(Name+"_DeltaEtatau2SVFit","DeltaEtatau2SVFit",100,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2SVFit=HConfig.GetTH1D(Name+"_DeltaPtau2SVFit","DeltaPtau2SVFit",50,-1,1,"#Delta P","Events");
  DeltaEtau2SVFit=HConfig.GetTH1D(Name+"_DeltaEtau2SVFit","DeltaEtau2SVFit",50,-1,1,"#Delta E","Events");
  //mixed
  DeltaPhitau1Mixed=HConfig.GetTH1D(Name+"_DeltaPhitau1Mixed","DeltaPhitau1Mixed",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1Mixed=HConfig.GetTH1D(Name+"_DeltaEtatau1Mixed","DeltaEtatau1Mixed",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1Mixed=HConfig.GetTH1D(Name+"_DeltaPtau1Mixed","DeltaPtau1Mixed",50,-1,1,"#Delta P","Events");
  DeltaEtau1Mixed=HConfig.GetTH1D(Name+"_DeltaEtau1Mixed","DeltaEtau1Mixed",50,-1,1,"#Delta E","Events");
  DeltaPhitau2Mixed=HConfig.GetTH1D(Name+"_DeltaPhitau2Mixed","DeltaPhitau2Mixed",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2Mixed=HConfig.GetTH1D(Name+"_DeltaEtatau2Mixed","DeltaEtatau2Mixed",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2Mixed=HConfig.GetTH1D(Name+"_DeltaPtau2Mixed","DeltaPtau2Mixed",50,-1,1,"#Delta P","Events");
  DeltaEtau2Mixed=HConfig.GetTH1D(Name+"_DeltaEtau2Mixed","DeltaEtau2Mixed",50,-1,1,"#Delta E","Events");
  //GEF
  DeltaPhitau1GEF=HConfig.GetTH1D(Name+"_DeltaPhitau1GEF","DeltaPhitau1GEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1GEF=HConfig.GetTH1D(Name+"_DeltaEtatau1GEF","DeltaEtatau1GEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1GEF=HConfig.GetTH1D(Name+"_DeltaPtau1GEF","DeltaPtau1GEF",50,-1,1,"#Delta P","Events");
  DeltaEtau1GEF=HConfig.GetTH1D(Name+"_DeltaEtau1GEF","DeltaEtau1GEF",50,-1,1,"#Delta E","Events");
  DeltaPhitau2GEF=HConfig.GetTH1D(Name+"_DeltaPhitau2GEF","DeltaPhitau2GEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2GEF=HConfig.GetTH1D(Name+"_DeltaEtatau2GEF","DeltaEtatau2GEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2GEF=HConfig.GetTH1D(Name+"_DeltaPtau2GEF","DeltaPtau2GEF",50,-1,1,"#Delta P","Events");
  DeltaEtau2GEF=HConfig.GetTH1D(Name+"_DeltaEtau2GEF","DeltaEtau2GEF",50,-1,1,"#Delta E","Events");
  //
  DeltaPhitauPiGEF=HConfig.GetTH1D(Name+"_DeltaPhitauPiGEF","DeltaPhitauPiGEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauPiGEF=HConfig.GetTH1D(Name+"_DeltaEtatauPiGEF","DeltaEtatauPiGEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauPiGEF=HConfig.GetTH1D(Name+"_DeltaPtauPiGEF","DeltaPtauPiGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauPiGEF=HConfig.GetTH1D(Name+"_DeltaEtauPiGEF","DeltaEtauPiGEF",50,-1,1,"#Delta E","Events");
  DeltaPhitauHGEF=HConfig.GetTH1D(Name+"_DeltaPhitauHGEF","DeltaPhitauHGEF",50,-0.1,0.1,"#Delta#phi","Events");
  DeltaEtatauHGEF=HConfig.GetTH1D(Name+"_DeltaEtatauHGEF","DeltaEtatauHGEF",50,-0.1,0.1,"#Delta#eta","Events");
  DeltaPtauHGEF=HConfig.GetTH1D(Name+"_DeltaPtauHGEF","DeltaPtauHGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauHGEF=HConfig.GetTH1D(Name+"_DeltaEtauHGEF","DeltaEtauHGEF",50,-1,1,"#Delta E","Events");

  //REF
  RefX=HConfig.GetTH1D(Name+"_RefX","RefX",50,-1,1,"ref","Events");
  RefY=HConfig.GetTH1D(Name+"_RefY","RefY",50,-1,1,"ref","Events");
  RefZ=HConfig.GetTH1D(Name+"_RefZ","RefZ",50,-1,1,"ref","Events");

  SVfitMTTdR1=HConfig.GetTH1D(Name+"_SVfitMTTdR1","SVfitMTTdR1",50,0,0.06,"#Delta (R) #tau 1","a.u");
  SVfitMTTdR2=HConfig.GetTH1D(Name+"_SVfitMTTdR2","SVfitMTTdR2",50,0,0.06,"#Delta (R) #tau 2","a.u");

  PcorrEtaSVfitMTT1=HConfig.GetTH2D(Name+"_PcorrEtaSVfitMTT1","PcorrEtaSVfitMTT1",50,-1,1,50,-0.2,0.2,"P SVfit","#Eta MTT");
  PcorrPhiSVfitMTT1=HConfig.GetTH2D(Name+"_PcorrPhiSVfitMTT1","PcorrPhiSVfitMTT1",50,-1,1,50,-0.2,0.2,"P SVfit","#Phi MTT");

  PcorrEtaSVfitMTT2=HConfig.GetTH2D(Name+"_PcorrEtaSVfitMTT2","PcorrEtaSVfitMTT2",50,-1,1,50,-0.2,0.2,"P SVfit","#Eta MTT");
  PcorrPhiSVfitMTT2=HConfig.GetTH2D(Name+"_PcorrPhiSVfitMTT2","PcorrPhiSVfitMTT2",50,-1,1,50,-0.2,0.2,"P SVfit","#Phi MTT");

  dRandPcorrEta1=HConfig.GetTH3D(Name+"_dRandPcorrEta1","dRandPcorrEta1",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Eta MTT","#Delta R");
  dRandPcorrPhi1=HConfig.GetTH3D(Name+"_dRandPcorrPhi1","dRandPcorrPhi1",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Phi MTT","#Delta R");

  dRandPcorrEta2=HConfig.GetTH3D(Name+"_dRandPcorrEta2","dRandPcorrEta2",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Eta MTT","#Delta R");
  dRandPcorrPhi2=HConfig.GetTH3D(Name+"_dRandPcorrPhi2","dRandPcorrPhi2",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Phi MTT","#Delta R");

  PullAcopPV=HConfig.GetTH1D(Name+"_PullAcopPV","PullAcopPV",50,-1,1,"pull AcopPV","u.a");
  dR1vsAcopPV=HConfig.GetTH2D(Name+"_dR1vsAcopPV","dR1vsAcopPV",50,0,0.06,50,-1,1,"dR","acop");
  dR2vsAcopPV=HConfig.GetTH2D(Name+"_dR2vsAcopPV","dR2vsAcopPV",50,0,0.06,50,-1,1,"dR","acop");
  P1vsAcopPV=HConfig.GetTH2D(Name+"_P1vsAcopPV","P1vsAcopPV",50,-1,1,50,-1,1,"P svfit","acop");
  P2vsAcopPV=HConfig.GetTH2D(Name+"_P2vsAcopPV","P2vsAcopPV",50,-1,1,50,-1,1,"P svfit","acop");
  Phi1vsAcopPV=HConfig.GetTH2D(Name+"_Phi1vsAcopPV","Phi1vsAcopPV",50,-0.2,0.2,50,-1,1,"phi mtt","acop");
  Eta1vsAcopPV=HConfig.GetTH2D(Name+"_Eta1vsAcopPV","Eta1vsAcopPV",50,-0.2,0.2,50,-1,1,"eta mtt","acop");
  Phi2vsAcopPV=HConfig.GetTH2D(Name+"_Phi2vsAcopPV","Phi2vsAcopPV",50,-0.2,0.2,50,-1,1,"phi mtt","acop");
  Eta2vsAcopPV=HConfig.GetTH2D(Name+"_Eta2vsAcopPV","Eta2vsAcopPV",50,-0.2,0.2,50,-1,1,"eta mtt","acop");

  Fraction1=HConfig.GetTH1D(Name+"_Fraction1","Fraction1",2,0,2," "," ");
  Fraction2=HConfig.GetTH1D(Name+"_Fraction2","Fraction2",2,0,2," "," ");

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove

}

void  HCPTauTau::Store_ExtraDist(){

  Extradist1d.push_back(&polarimetricAcopAngle);
  Extradist1d.push_back(&polarimetricGEFAcopAngle);
  Extradist1d.push_back(&decayplaneAcopAngle);
  Extradist1d.push_back(&impactparameterAcopAngle);
  Extradist1d.push_back(&DPIPAcopAngle);
  Extradist1d.push_back(&PVIPAcopAngle);

  Extradist1d.push_back(&polarimetricAcopAngleTruth);
  Extradist1d.push_back(&decayplaneAcopAngleTruth);
  Extradist1d.push_back(&impactparameterAcopAngleTruth);
  Extradist1d.push_back(&DPIPAcopAngleTruth);
  Extradist1d.push_back(&PVIPAcopAngleTruth);

  //mtt
  Extradist1d.push_back(&DeltaPhitau1MTT);
  Extradist1d.push_back(&DeltaEtatau1MTT);
  Extradist1d.push_back(&DeltaPtau1MTT);
  Extradist1d.push_back(&DeltaEtau1MTT);
  Extradist1d.push_back(&DeltaPhitau2MTT);
  Extradist1d.push_back(&DeltaEtatau2MTT);
  Extradist1d.push_back(&DeltaPtau2MTT);
  Extradist1d.push_back(&DeltaEtau2MTT);
  //svfit
  Extradist1d.push_back(&DeltaPhitau1SVFit);
  Extradist1d.push_back(&DeltaEtatau1SVFit);
  Extradist1d.push_back(&DeltaPtau1SVFit);
  Extradist1d.push_back(&DeltaEtau1SVFit);
  Extradist1d.push_back(&DeltaPhitau2SVFit);
  Extradist1d.push_back(&DeltaEtatau2SVFit);
  Extradist1d.push_back(&DeltaPtau2SVFit);
  Extradist1d.push_back(&DeltaEtau2SVFit);
  //mixed
  Extradist1d.push_back(&DeltaPhitau1Mixed);
  Extradist1d.push_back(&DeltaEtatau1Mixed);
  Extradist1d.push_back(&DeltaPtau1Mixed);
  Extradist1d.push_back(&DeltaEtau1Mixed);
  Extradist1d.push_back(&DeltaPhitau2Mixed);
  Extradist1d.push_back(&DeltaEtatau2Mixed);
  Extradist1d.push_back(&DeltaPtau2Mixed);
  Extradist1d.push_back(&DeltaEtau2Mixed);
  //GEF
  Extradist1d.push_back(&DeltaPhitau1GEF);
  Extradist1d.push_back(&DeltaEtatau1GEF);
  Extradist1d.push_back(&DeltaPtau1GEF);
  Extradist1d.push_back(&DeltaEtau1GEF);
  Extradist1d.push_back(&DeltaPhitau2GEF);
  Extradist1d.push_back(&DeltaEtatau2GEF);
  Extradist1d.push_back(&DeltaPtau2GEF);
  Extradist1d.push_back(&DeltaEtau2GEF);
  //
  Extradist1d.push_back(&DeltaPhitauPiGEF);
  Extradist1d.push_back(&DeltaEtatauPiGEF);
  Extradist1d.push_back(&DeltaPtauPiGEF);
  Extradist1d.push_back(&DeltaEtauPiGEF);
  Extradist1d.push_back(&DeltaPhitauHGEF);
  Extradist1d.push_back(&DeltaEtatauHGEF);
  Extradist1d.push_back(&DeltaPtauHGEF);
  Extradist1d.push_back(&DeltaEtauHGEF);
  //ref
  Extradist1d.push_back(&RefX);
  Extradist1d.push_back(&RefY);
  Extradist1d.push_back(&RefZ);

  Extradist1d.push_back(&SVfitMTTdR1);
  Extradist1d.push_back(&SVfitMTTdR2);

  Extradist2d.push_back(&PcorrEtaSVfitMTT1);
  Extradist2d.push_back(&PcorrPhiSVfitMTT1);

  Extradist2d.push_back(&PcorrEtaSVfitMTT2);
  Extradist2d.push_back(&PcorrPhiSVfitMTT2);

  Extradist3d.push_back(&dRandPcorrEta1);
  Extradist3d.push_back(&dRandPcorrPhi1);

  Extradist3d.push_back(&dRandPcorrEta2);
  Extradist3d.push_back(&dRandPcorrPhi2);
 
  Extradist1d.push_back(&PullAcopPV);
  Extradist2d.push_back(&dR1vsAcopPV);
  Extradist2d.push_back(&dR2vsAcopPV);
  Extradist2d.push_back(&P1vsAcopPV);
  Extradist2d.push_back(&P2vsAcopPV);
  Extradist2d.push_back(&Phi1vsAcopPV);
  Extradist2d.push_back(&Eta1vsAcopPV);
  Extradist2d.push_back(&Phi2vsAcopPV);
  Extradist2d.push_back(&Eta2vsAcopPV);

  Extradist1d.push_back(&Fraction1);
  Extradist1d.push_back(&Fraction2);
}

void  HCPTauTau::doEvent()  { //  Method called on every event
} //do event

//  This is a function if you want to do something after the event loop
void HCPTauTau::Finish() {

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
