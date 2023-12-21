//Ntuple_Controller.cxx IMPLEMENTATION FILE
 

#include "Ntuple_Controller.h"
#include "Tools.h"
#include "TF1.h"
#include "Parameters.h"
#include "Logger.h"
#include <tuple>

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles, TString SysTree):
  copyTree(false)
  ,cannotObtainHiggsMass(false)
  ,ObjEvent(-1)
  ,isInit(false)
{  
  // TChains the ROOTuple file
  TChain *chain = NULL;
  if(GetInputNtuplePath().Contains("SingleMuon")) {
    chain = new TChain("HTauTauTree/Nominal");
  }
  else {
    if(SysTree == "TESUp") chain = new TChain("HTauTauTree/TESUp");
    else if(SysTree == "TESDown") chain = new TChain("HTauTauTree/TESDown");
    else if(SysTree == "MESUp") chain = new TChain("HTauTauTree/MESUp");
    else if(SysTree == "MESDown") chain = new TChain("HTauTauTree/MESDown");
    else if(SysTree == "JERUp") chain = new TChain("HTauTauTree/JERUp");
    else if(SysTree == "JERDown") chain = new TChain("HTauTauTree/JERDown");
    else if(SysTree == "METResoUp") chain = new TChain("HTauTauTree/METResoUp");
    else if(SysTree == "METResoDown") chain = new TChain("HTauTauTree/METResoDown");
    else if(SysTree == "METScaleUp") chain = new TChain("HTauTauTree/METScaleUp");
    else if(SysTree == "METScaleDown") chain = new TChain("HTauTauTree/METScaleDown");
    else if(SysTree == "METUnclusteredUp") chain = new TChain("HTauTauTree/METUnclusteredUp");
    else if(SysTree == "METUnclusteredDown") chain = new TChain("HTauTauTree/METUnclusteredDown");
    else if(SysTree == "FlavorQCDUp") chain = new TChain("HTauTauTree/FlavorQCDUp");
    else if(SysTree == "FlavorQCDDown") chain = new TChain("HTauTauTree/FlavorQCDDown");
    else if(SysTree == "RelativeBalUp") chain = new TChain("HTauTauTree/RelativeBalUp");
    else if(SysTree == "RelativeBalDown") chain = new TChain("HTauTauTree/RelativeBalDown");
    else if(SysTree == "HFUp") chain = new TChain("HTauTauTree/HFUp");
    else if(SysTree == "HFDown") chain = new TChain("HTauTauTree/HFDown");
    else if(SysTree == "BBEC1Up") chain = new TChain("HTauTauTree/BBEC1Up");
    else if(SysTree == "BBEC1Down") chain = new TChain("HTauTauTree/BBEC1Down");
    else if(SysTree == "EC2Up") chain = new TChain("HTauTauTree/EC2Up");
    else if(SysTree == "EC2Down") chain = new TChain("HTauTauTree/EC2Down");
    else if(SysTree == "AbsoluteUp") chain = new TChain("HTauTauTree/AbsoluteUp");
    else if(SysTree == "AbsoluteDown") chain = new TChain("HTauTauTree/AbsoluteDown");
    else if(SysTree == "BBEC1_YEARUp") chain = new TChain("HTauTauTree/BBEC1_YEARUp");
    else if(SysTree == "BBEC1_YEARDown") chain = new TChain("HTauTauTree/BBEC1_YEARDown");
    else if(SysTree == "EC2_YEARUp") chain = new TChain("HTauTauTree/EC2_YEARUp");
    else if(SysTree == "EC2_YEARDown") chain = new TChain("HTauTauTree/EC2_YEARDown");
    else if(SysTree == "Absolute_YEARUp") chain = new TChain("HTauTauTree/Absolute_YEARUp");
    else if(SysTree == "Absolute_YEARDown") chain = new TChain("HTauTauTree/Absolute_YEARDown");
    else if(SysTree == "HF_YEARUp") chain = new TChain("HTauTauTree/HF_YEARUp");
    else if(SysTree == "HF_YEARDown") chain = new TChain("HTauTauTree/HF_YEARDown");
    else if(SysTree == "RelativeSample_YEARUp") chain = new TChain("HTauTauTree/RelativeSample_YEARUp");
    else if(SysTree == "RelativeSample_YEARDown") chain = new TChain("HTauTauTree/RelativeSample_YEARDown");
    else chain = new TChain("HTauTauTree/Nominal");
  }
  Logger(Logger::Verbose) << "Loading " << RootFiles.size() << " files" << std::endl;
  int chainsize=0;
  bool failed=false;
  for(unsigned int i=0; i<RootFiles.size(); i++){
    chainsize=chain->Add(RootFiles[i]);
    if(chainsize!=1)
      {
	sleep(10);
	chainsize=chain->Add(RootFiles[i]);
	if(chainsize!=1)
	  {
	    sleep(10);
	    chainsize=chain->Add(RootFiles[i]);
	    if(chainsize!=1)
	      {
		sleep(10);
		chainsize=chain->Add(RootFiles[i]);
		if(chainsize!=1)failed=true; 
	      }
	  }
      }
  }
  if(failed==true)Logger(Logger::Error) << "Some root files are not added!!!" << std::endl;
  TTree *tree = (TTree*)chain;
  if(chain==0){
    Logger(Logger::Error) << "chain points to NULL" << std::endl;
  }
  //tree->SetBranchStatus("daughters_byIsolationMVA*old*",0);
  //tree->SetBranchStatus("trigger_name",0);
  Logger(Logger::Info) << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree, SysTree);
  nbytes=0; 
  nb=0;
  Logger(Logger::Info) << "Ntuple Configured" << std::endl;


  // TFile *currentFile = chain->GetCurrentFile();
  // hLLRCounters = (TH1F*)currentFile->Get("HTauTauTree/Counters");
  // if(!hLLRCounters) std::cout<<"Counters histogram not found!"<<std::endl;
 


  // Fit setup 

  // Resolution uncertainty setup

  gRandom->SetSeed(1234);

  // Rochester muon momentum corrections

  //  rmcor = new rochcor2012(); // For systematics use rmcor = new rochcor2012(seed!=1234);

  // Set object correction flags to default values
  tauCorrection = "";
  muonCorrection = "";
  elecCorrection = "";
  jetCorrection = "";

}

///////////////////////////////////////////////////////////////////////
//
// Function: void InitEvent()
//
// Purpose: Initialize variables etc on event base
//
///////////////////////////////////////////////////////////////////////

void Ntuple_Controller::InitEvent(){
  Muon_corrected_p4.clear();
  //	Muon_corrected_p4.resize(NMuons());
  Muon_isCorrected = false;

  // after everything is initialized
  isInit = true;
}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
  return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
  jentry=_jentry;
  Ntp->LoadTree(jentry);
  nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
  isInit = false;
  InitEvent();
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
  return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
  return Ntp->fChain->GetCurrentFile()->GetName();
}

TString Ntuple_Controller::GetInputNtuplePath(){
  Parameters Par; // assumes configured in Analysis.cxx
  TString dsPath;
  Par.GetString("InputNtuples:",dsPath);
  return dsPath;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
  Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
  Logger(Logger::Verbose) << "Cleaning up" << std::endl;
  delete Ntp;
  //  delete rmcor;
  Logger(Logger::Verbose) << "Complete." << std::endl;
}

void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
    Logger(Logger::Info) << "Starting D3PD cloning" << std::endl;
    
    newfile = TFile::Open(n+".root","recreate");
    if(newfile->IsZombie())Logger(Logger::Info) << "Zombie Output" << endl;
    if(!newfile->IsOpen())Logger(Logger::Info) << "Output not openned" << endl;
    Ntp->fChain->SetBranchStatus("daughters_byIsolationMVArun2v1DBoldDMwLTrawc",0);
    Ntp->fChain->SetBranchStatus("isOSCand",0);
    Ntp->fChain->SetBranchStatus("*_UP_*",0);
    Ntp->fChain->SetBranchStatus("*TauUp*",0);
    Ntp->fChain->SetBranchStatus("*EleUp*",0);
    Ntp->fChain->SetBranchStatus("*_DOWN_*",0);
    Ntp->fChain->SetBranchStatus("*Down*",0);
    Ntp->fChain->SetBranchStatus("*_cov*",0);
    Ntp->fChain->SetBranchStatus("MET_significance",0);
    Ntp->fChain->SetBranchStatus("mT_Dau*",0);
    Ntp->fChain->SetBranchStatus("aMCatNLOweight",0);
    Ntp->fChain->SetBranchStatus("PFTau_Track_*",0);
    Ntp->fChain->SetBranchStatus("susyModel",0);
    Ntp->fChain->SetBranchStatus("*jets_deepFlavor_*",0);
    Ntp->fChain->SetBranchStatus("ak8jets_PrunedMass",0);
    Ntp->fChain->SetBranchStatus("ak8jets_TrimmedMass",0);
    Ntp->fChain->SetBranchStatus("ak8jets_FilteredMass",0);
    Ntp->fChain->SetBranchStatus("VertexHash*TracksRemovedOld*",0);
    Ntp->fChain->SetBranchStatus("trg_*",0);
    Ntp->fChain->SetBranchStatus("trg_doubletau",1);
    Ntp->fChain->SetBranchStatus("byIsolationMVA3oldDMwLTraw_*",0);
    Ntp->fChain->SetBranchStatus("ptvis",0);
    Ntp->fChain->SetBranchStatus("*p*id_*",0);
    Ntp->fChain->SetBranchStatus("*_sv",0);
    Ntp->fChain->SetBranchStatus("deepTauVs*Raw_*",0);
    Ntp->fChain->SetBranchStatus("pvx",0);
    Ntp->fChain->SetBranchStatus("pvy",0);
    Ntp->fChain->SetBranchStatus("pvz",0);
    Ntp->fChain->SetBranchStatus("n*_*",0);
    Ntp->fChain->SetBranchStatus("sv*_*",0);
    Ntp->fChain->SetBranchStatus("pt_*",0);
    Ntp->fChain->SetBranchStatus("eta_*",0);
    Ntp->fChain->SetBranchStatus("phi_*",0);
    Ntp->fChain->SetBranchStatus("m_*",0);
    Ntp->fChain->SetBranchStatus("q_*",0);
    Ntp->fChain->SetBranchStatus("d0_*",0);
    Ntp->fChain->SetBranchStatus("dz_*",0);
    Ntp->fChain->SetBranchStatus("iso_*",0);
    
    SkimmedTree=Ntp->fChain->CloneTree(0);
    copyTree=true;
  }
}

void Ntuple_Controller::SaveCloneTree(){
  if(copyTree){
    SkimmedTree->AutoSave();
    newfile->Close();
  }
  Logger(Logger::Info) << "Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  Logger(Logger::Warning) << "ThinTree not implemented." << std::endl;
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
    doElectrons();
    doPhotons();
    doJets();
    doMuons();
    doTaus();
    doMET();
  }
}

void Ntuple_Controller::doElectrons(){
  electrons.clear();
  electrons_default.clear();
}

void Ntuple_Controller::doPhotons(){
  photons.clear();
  photons_default.clear();

}

void Ntuple_Controller::doJets(){
  jets.clear();
  jets_default.clear();
}

void Ntuple_Controller::doMuons(){
  muons.clear();
  muons_default.clear();
}

void Ntuple_Controller::doTaus(){
  taus.clear();
  taus_default.clear();
}

void Ntuple_Controller::doMET(){
}

//Physics get Functions
Long64_t Ntuple_Controller::GetMCID() {
	
  Long64_t DataMCTypeFromTupel = Ntp->_Id;
  int dmcType = -999;

  if (HistoC.hasID(DataMCTypeFromTupel % 100)) {
    dmcType = DataMCTypeFromTupel % 100;
  }
  else dmcType=DataMCTypeFromTupel % 100;
  if (DataMCTypeFromTupel==701)return 701;
  else if (DataMCTypeFromTupel==702)return 702;
  else if (DataMCTypeFromTupel==703)return 703;
  else if (dmcType==2)return 202;
  else if (dmcType==1 && DataMCTypeFromTupel!=1)return 201;
  else if (dmcType==3)return 203;
  else if (dmcType==61)return 461;
  else if (dmcType==60 && DataMCTypeFromTupel!=60)return 460;
  else return dmcType;
}


