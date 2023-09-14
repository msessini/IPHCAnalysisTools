//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 15 17:49:09 2017 by ROOT version 5.34/18
// from TTree HTauTauTree/HTauTauTree
// found on file: /home-pbs/vcherepa/cms_work/CMSSW_8_0_25/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis.root
//////////////////////////////////////////////////////////

#ifndef NtupleReader_h
#define NtupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupleReader {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Leaf type declarations
  int _tauIndex;
  int _tauGenMatch;
  float _tauDM;
  float _tauPt;
  float _tauEta;
  float _tauPhi;
  float _tauE;
  float _GEFtauPt;
  float _GEFtauEta;
  float _GEFtauPhi;
  float _GEFtauE;
  float _tauSVx;
  float _tauSVy;
  float _tauSVz;
  float _tauIPx;
  float _tauIPy;
  float _tauIPz;
  float _tauIPsignificance;
  int _muIndex;
  int _muGenMatch;
  float _muPt;
  float _muEta;
  float _muPhi;
  float _muE;
  double _muIso;
  float _muIPx;
  float _muIPy;
  float _muIPz;
  float _muIPsignificance;
  float _genTaupx;
  float _genTaupy;
  float _genTaupz;
  float _genTauE;
  float _genTauSVx;
  float _genTauSVy;
  float _genTauSVz;
  float _genMuonpx;
  float _genMuonpy;
  float _genMuonpz;
  float _genMuonE;
  float _genMuonIPx;
  float _genMuonIPy;
  float _genMuonIPz;
  float _genPVx;
  float _genPVy;
  float _genPVz;
  double _gendpPhiCP;
  double _genpvPhiCP;
  bool _isOSpair;
  bool _isIso;
  bool _isMediumID;
  bool _isTightJetID;
  bool _trgIsoMu;
  bool _trgXMuTau;
  double _pairvisMass;
  int _Njets;
  int _Nbtags;
  double _leadingjetPt;
  double _subleadingjetPt;
  double _dijetMass;
  double _dijetPt;
  double _dijetdeltaEta;
  double _ditauPt;
  float _muMETmt;
  double _fastMTTmass;
  float _PUPPImet;
  float _PUPPImetphi;
  float _PUPPIMETCov00;
  float _PUPPIMETCov10;
  float _PUPPIMETCov11;
  double _pvPhiCP;
  double _dpPhiCP;
  double _wEven;
  double _wOdd;
  double _wMM;
  double _wPrefiring;
  double _wPrefiringUp;
  double _wPrefiringDown;
  double _wIDvsJet;
  double _wIDvsJetUp;
  double _wIDvsJetDown;
  double _wIDvsEle;
  double _wIDvsEleUp;
  double _wIDvsEleDown;
  double _wIDvsMu;
  double _wIDvsMuUp;
  double _wIDvsMuDown;
  double _wTrg;
  double _wTrgUp;
  double _wTrgDown;
  double _wIDMu;
  double _wTrkMu;
  double _wPU;
  double _wZpT;
  double _wZpTUp;
  double _wZpTDown;
  double _wToppT;
  double _wToppTUp;
  double _wToppTDown;
  double _wBtag;
  double _wBtagUp;
  double _wBtagDown;
  double _wPSISRUp;
  double _wPSISRDown;
  double _wPSFSRUp;
  double _wPSFSRDown;
  double _wScaleUp;
  double _wScaleDown;
  double _wMC;
  double _wSignal;
  double _wTot;
  float _pvx;
  float _pvy;
  float _pvz;
  float _pvCov00;
  float _pvCov11;
  float _pvCov22;
  float _pvCov01;
  float _pvCov02;
  float _pvCov12;
  Int_t _Id;
  Int_t _trueId;
  int _Npartons;
  bool _isData;
  //bool _isEmbed;
  Int_t _theYear;
  ULong64_t _runNumber;
  ULong64_t _evt;
  Int_t _lumi;
  bool _isZ;
  bool _isW;
  bool _isSignal;
  bool _isQCD;
  bool _isVV;
  bool _isTTbar;
  bool _isSingleTop;

  // Branch declarations
  TBranch *b_tauIndex;
  TBranch *b_tauGenMatch;
  TBranch *b_tauDM;
  TBranch *b_tauPt;
  TBranch *b_tauEta;
  TBranch *b_tauPhi;
  TBranch *b_tauE;
  TBranch *b_GEFtauE;
  TBranch *b_GEFtauPt;
  TBranch *b_GEFtauPhi;
  TBranch *b_GEFtauEta;
  TBranch *b_tauSVx;
  TBranch *b_tauSVy;
  TBranch *b_tauSVz;
  TBranch *b_muIndex;
  TBranch *b_muGenMatch;
  TBranch *b_muPt;
  TBranch *b_muEta;
  TBranch *b_muPhi;
  TBranch *b_muE;
  TBranch *b_muIso;
  TBranch *b_muIPx;
  TBranch *b_muIPy;
  TBranch *b_muIPz;
  TBranch *b_muIPsignificance;
  TBranch *b_genTaupx;
  TBranch *b_genTaupy;
  TBranch *b_genTaupz;
  TBranch *b_genTauE;
  TBranch *b_genTauSVx;
  TBranch *b_genTauSVy;
  TBranch *b_genTauSVz;
  TBranch *b_genMuonpx;
  TBranch *b_genMuonpy;
  TBranch *b_genMuonpz;
  TBranch *b_genMuonE;
  TBranch *b_genMuonIPx;
  TBranch *b_genMuonIPy;
  TBranch *b_genMuonIPz;
  TBranch *b_genPVx;
  TBranch *b_genPVy;
  TBranch *b_genPVz;
  TBranch *b_gendpPhiCP;
  TBranch *b_genpvPhiCP;
  TBranch *b_isOSpair;
  TBranch *b_isIso;
  TBranch *b_isMediumID;
  TBranch *b_isTightJetID;
  TBranch *b_trgIsoMu;
  TBranch *b_trgXMuTau;
  TBranch *b_pairvisMass;
  TBranch *b_Njets;
  TBranch *b_Nbtags;
  TBranch *b_leadingjetPt;
  TBranch *b_subleadingjetPt;
  TBranch *b_dijetMass;
  TBranch *b_dijetPt;
  TBranch *b_dijetdeltaEta;
  TBranch *b_ditauPt;
  TBranch *b_muMETmt;
  TBranch *b_fastMTTmass;
  TBranch *b_PUPPImet;
  TBranch *b_PUPPImetphi;
  TBranch *b_PUPPIMETCov00;
  TBranch *b_PUPPIMETCov10;
  TBranch *b_PUPPIMETCov11;
  TBranch *b_pvPhiCP;
  TBranch *b_dpPhiCP;
  TBranch *b_wEven;
  TBranch *b_wOdd;
  TBranch *b_wMM;
  TBranch *b_wPrefiring;
  TBranch *b_wPrefiringUp;
  TBranch *b_wPrefiringDown;
  TBranch *b_wIDvsJet;
  TBranch *b_wIDvsJetUp;
  TBranch *b_wIDvsJetDown;
  TBranch *b_wIDvsEle;
  TBranch *b_wIDvsEleUp;
  TBranch *b_wIDvsEleDown;
  TBranch *b_wIDvsMu;
  TBranch *b_wIDvsMuUp;
  TBranch *b_wIDvsMuDown;
  TBranch *b_wTrg;
  TBranch *b_wTrgUp;
  TBranch *b_wTrgDown;
  TBranch *b_wIDMu;
  TBranch *b_wTrkMu;
  TBranch *b_wPU;
  TBranch *b_wZpT;
  TBranch *b_wZpTUp;
  TBranch *b_wZpTDown;
  TBranch *b_wToppT;
  TBranch *b_wToppTUp;
  TBranch *b_wToppTDown;
  TBranch *b_wBtag;
  TBranch *b_wBtagUp;
  TBranch *b_wBtagDown;
  TBranch *b_wPSISRUp;
  TBranch *b_wPSISRDown;
  TBranch *b_wPSFSRUp;
  TBranch *b_wPSFSRDown;
  TBranch *b_wScaleUp;
  TBranch *b_wScaleDown;
  TBranch *b_wMC;
  TBranch *b_wSignal;
  TBranch *b_wTot;
  TBranch *b_pvx;
  TBranch *b_pvy;
  TBranch *b_pvz;
  TBranch *b_pvCov00;
  TBranch *b_pvCov11;
  TBranch *b_pvCov22;
  TBranch *b_pvCov01;
  TBranch *b_pvCov02;
  TBranch *b_pvCov12;
  TBranch *b_Id;
  TBranch *b_trueId;
  TBranch *b_Npartons;
  TBranch *b_isData;
  //TBranch *b_isEmbed;
  TBranch *b_theYear;
  TBranch *b_runNumber;
  TBranch *b_evt;
  TBranch *b_lumi;
  TBranch *b_isZ;
  TBranch *b_isW;
  TBranch *b_isSignal;
  TBranch *b_isQCD;
  TBranch *b_isVV;
  TBranch *b_isTTbar;
  TBranch *b_isSingleTop;

  NtupleReader(TTree *tree=0, TString Sys=" ");
  virtual ~NtupleReader();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree, TString Sys);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupleReader_cxx
NtupleReader::NtupleReader(TTree *tree, TString Sys) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {


#ifdef SINGLE_TREE
    // The following code should be used if you want this class to access
    // a single tree instead of a chain
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
    if (!f || !f->IsOpen()) {
      f = new TFile("Memory Directory");
    }
    f->GetObject("HTauTauTree/HTauTauTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
    TChain * chain = new TChain("HTauTauTree/HTauTauTree","");
    tree = chain;
#endif // SINGLE_TREE

  }
  Init(tree, Sys);
}



NtupleReader::~NtupleReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t NtupleReader::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t NtupleReader::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void NtupleReader::Init(TTree *tree, TString Sys)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  //

  // Set object pointer
  _tauIndex = 0;
  _tauGenMatch = 0;
  _tauDM = 0;
  _tauPt = 0;
  _tauEta = 0;
  _tauPhi = 0;
  _tauE = 0;
  _GEFtauPt = 0;
  _GEFtauEta = 0;
  _GEFtauPhi = 0;
  _GEFtauE = 0;
  _tauSVx = 0;
  _tauSVy = 0;
  _tauSVz = 0;
  _muIndex = 0;
  _muGenMatch = 0;
  _muPt = 0;
  _muEta = 0;
  _muPhi = 0;
  _muE = 0;
  _muIso = 0;
  _muIPx = 0;
  _muIPy = 0;
  _muIPz = 0;
  _muIPsignificance = 0;
  _genTaupx = 0;
  _genTaupy = 0;
  _genTaupz = 0;
  _genTauE = 0;
  _genTauSVx = 0;
  _genTauSVy = 0;
  _genTauSVz = 0;
  _genMuonpx = 0;
  _genMuonpy = 0;
  _genMuonpz = 0;
  _genMuonE = 0;
  _genMuonIPx = 0;
  _genMuonIPy = 0;
  _genMuonIPz = 0;
  _genPVx = 0;
  _genPVy = 0;
  _genPVz = 0;
  _gendpPhiCP = 0;
  _genpvPhiCP = 0;
  _isOSpair = 0;
  _isIso = 0;
  _isMediumID = 0;
  _isTightJetID = 0;
  _trgIsoMu = 0;
  _trgXMuTau = 0;
  _pairvisMass = 0;
  _Njets = 0;
  _Nbtags = 0;
  _leadingjetPt = 0;
  _subleadingjetPt = 0;
  _dijetMass = 0;
  _dijetPt = 0;
  _dijetdeltaEta = 0;
  _ditauPt = 0;
  _muMETmt = 0;
  _fastMTTmass = 0;
  _PUPPImet = 0;
  _PUPPImetphi = 0;
  _PUPPIMETCov00 = 0;
  _PUPPIMETCov10 = 0;
  _PUPPIMETCov11 = 0;
  _pvPhiCP = 0;
  _dpPhiCP = 0;
  _wEven = 0;
  _wOdd = 0;
  _wMM = 0;
  _wPrefiring = 0;
  _wPrefiringUp = 0;
  _wPrefiringDown = 0;
  _wIDvsJet = 0;
  _wIDvsJetUp = 0;
  _wIDvsJetDown = 0;
  _wIDvsEle = 0;
  _wIDvsEleUp = 0;
  _wIDvsEleDown = 0;
  _wIDvsMu = 0;
  _wIDvsMuUp = 0;
  _wIDvsMuDown = 0;
  _wTrg = 0;
  _wTrgUp = 0;
  _wTrgDown = 0;
  _wIDMu = 0;
  _wTrkMu = 0;
  _wPU = 0;
  _wZpT = 0;
  _wZpTUp = 0;
  _wZpTDown = 0;
  _wToppT = 0;
  _wToppTUp = 0;
  _wToppTDown = 0;
  _wBtag = 0;
  _wBtagUp = 0;
  _wBtagDown = 0;
  _wPSISRUp = 0;
  _wPSISRDown = 0;
  _wPSFSRUp = 0;
  _wPSFSRDown = 0;
  _wScaleUp = 0;
  _wScaleDown = 0;
  _wMC = 0;
  _wSignal = 0;
  _wTot = 0;
  _pvx = 0;
  _pvy = 0;
  _pvz = 0;
  _pvCov00 = 0;
  _pvCov11 = 0;
  _pvCov22 = 0;
  _pvCov01 = 0;
  _pvCov02 = 0;
  _pvCov12 = 0;
  _Id = 0;
  _trueId = 0;
  _Npartons = 0;
  _isData = 0;
  //_isEmbed = 0;
  _theYear = 0;
  _runNumber = 0;
  _evt = 0;
  _lumi = 0;
  _isZ = 0;
  _isW = 0;
  _isSignal = 0;
  _isQCD = 0;
  _isVV = 0;
  _isTTbar = 0;
  _isSingleTop = 0;
   
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("run", &_runNumber, &b_runNumber);
  fChain->SetBranchAddress("evt", &_evt, &b_evt);
  fChain->SetBranchAddress("lumi", &_lumi, &b_lumi);
  fChain->SetBranchAddress("year", &_theYear, &b_theYear);
  fChain->SetBranchAddress("tauIndex", &_tauIndex, &b_tauIndex);
  fChain->SetBranchAddress("tauGenMatch", &_tauGenMatch, &b_tauGenMatch);
  fChain->SetBranchAddress("tauPt", &_tauPt, &b_tauPt);
  fChain->SetBranchAddress("tauEta", &_tauEta, &b_tauEta);
  fChain->SetBranchAddress("tauPhi", &_tauPhi, &b_tauPhi);
  fChain->SetBranchAddress("tauE", &_tauE, &b_tauE);
  fChain->SetBranchAddress("GEFtauE", &_GEFtauE, &b_GEFtauE);
  fChain->SetBranchAddress("GEFtauPt", &_GEFtauPt, &b_GEFtauPt);
  fChain->SetBranchAddress("GEFtauPhi", &_GEFtauPhi, &b_GEFtauPhi);
  fChain->SetBranchAddress("GEFtauEta", &_GEFtauEta, &b_GEFtauEta);
  fChain->SetBranchAddress("tauSVx", &_tauSVx, &b_tauSVx);
  fChain->SetBranchAddress("tauSVy", &_tauSVy, &b_tauSVy);
  fChain->SetBranchAddress("tauSVz", &_tauSVz, &b_tauSVz);
  fChain->SetBranchAddress("tauDM", &_tauDM, &b_tauDM);
  fChain->SetBranchAddress("muIndex", &_muIndex, &b_muIndex);
  fChain->SetBranchAddress("muGenMatch", &_muGenMatch, &b_muGenMatch);
  fChain->SetBranchAddress("muPt", &_muPt, &b_muPt);
  fChain->SetBranchAddress("muEta", &_muEta, &b_muEta);
  fChain->SetBranchAddress("muPhi", &_muPhi, &b_muPhi);
  fChain->SetBranchAddress("muE", &_muE, &b_muE);
  fChain->SetBranchAddress("muIso", &_muIso, &b_muIso);
  fChain->SetBranchAddress("isOSpair", &_isOSpair, &b_isOSpair);
  fChain->SetBranchAddress("isIso", &_isIso, &b_isIso);
  fChain->SetBranchAddress("isMediumID", &_isMediumID, &b_isMediumID);
  fChain->SetBranchAddress("isTightJetID", &_isTightJetID, &b_isTightJetID);
  fChain->SetBranchAddress("trgIsoMu", &_trgIsoMu, &b_trgIsoMu);
  fChain->SetBranchAddress("trgXMuTau", &_trgXMuTau, &b_trgXMuTau);
  fChain->SetBranchAddress("pairvisMass", &_pairvisMass, &b_pairvisMass);
  fChain->SetBranchAddress("Njets", &_Njets, &b_Njets);
  fChain->SetBranchAddress("Nbtags", &_Nbtags, &b_Nbtags);
  fChain->SetBranchAddress("leadingjetPt", &_leadingjetPt, &b_leadingjetPt);
  fChain->SetBranchAddress("subleadingjetPt", &_subleadingjetPt, &b_subleadingjetPt);
  fChain->SetBranchAddress("dijetPt", &_dijetPt, &b_dijetPt);
  fChain->SetBranchAddress("dijetMass", &_dijetMass, &b_dijetMass);
  fChain->SetBranchAddress("dijetdeltaEta", &_dijetdeltaEta, &b_dijetdeltaEta);
  fChain->SetBranchAddress("ditauPt", &_ditauPt, &b_ditauPt);
  fChain->SetBranchAddress("muMETmt", &_muMETmt, &b_muMETmt);
  fChain->SetBranchAddress("fastMTTmass", &_fastMTTmass, &b_fastMTTmass);
  fChain->SetBranchAddress("PUPPImet", &_PUPPImet, &b_PUPPImet);
  fChain->SetBranchAddress("PUPPImetphi", &_PUPPImetphi, &b_PUPPImetphi);
  fChain->SetBranchAddress("PUPPIMETCov00", &_PUPPIMETCov00, &b_PUPPIMETCov00);
  fChain->SetBranchAddress("PUPPIMETCov10", &_PUPPIMETCov10, &b_PUPPIMETCov10);
  fChain->SetBranchAddress("PUPPIMETCov11", &_PUPPIMETCov11, &b_PUPPIMETCov11);
  fChain->SetBranchAddress("pvPhiCP", &_pvPhiCP, &b_pvPhiCP);
  fChain->SetBranchAddress("dpPhiCP", &_dpPhiCP, &b_dpPhiCP);
  fChain->SetBranchAddress("MCId", &_Id, &b_Id);
  fChain->SetBranchAddress("Npartons", &_Npartons, &b_Npartons);
  fChain->SetBranchAddress("isData", &_isData, &b_isData);
  //fChain->SetBranchAddress("isEmbed", &_isEmbed, &b_isEmbed);
  fChain->SetBranchAddress("isZ", &_isZ, &b_isZ);
  fChain->SetBranchAddress("isW", &_isW, &b_isW);
  fChain->SetBranchAddress("isSignal", &_isSignal, &b_isSignal);
  fChain->SetBranchAddress("isQCD", &_isQCD, &b_isQCD);
  fChain->SetBranchAddress("isVV", &_isVV, &b_isVV);
  fChain->SetBranchAddress("isTTbar", &_isTTbar, &b_isTTbar);
  fChain->SetBranchAddress("isSingleTop", &_isSingleTop, &b_isSingleTop);
  fChain->SetBranchAddress("wEven", &_wEven, &b_wEven);
  fChain->SetBranchAddress("wOdd", &_wOdd, &b_wOdd);
  fChain->SetBranchAddress("wMM", &_wMM, &b_wMM);
  fChain->SetBranchAddress("wPrefiring", &_wPrefiring, &b_wPrefiring);
  fChain->SetBranchAddress("wIDvsJet", &_wIDvsJet, &b_wIDvsJet);
  fChain->SetBranchAddress("wIDvsEle", &_wIDvsEle, &b_wIDvsEle);
  fChain->SetBranchAddress("wIDvsMu", &_wIDvsMu, &b_wIDvsMu);
  fChain->SetBranchAddress("wTrg", &_wTrg, &b_wTrg);
  fChain->SetBranchAddress("wIDMu", &_wIDMu, &b_wIDMu);
  fChain->SetBranchAddress("wTrkMu", &_wTrkMu, &b_wTrkMu);
  fChain->SetBranchAddress("wPU", &_wPU, &b_wPU);
  fChain->SetBranchAddress("wZpT", &_wZpT, &b_wZpT);
  fChain->SetBranchAddress("wToppT", &_wToppT, &b_wToppT);
  fChain->SetBranchAddress("wBtag", &_wBtag, &b_wBtag);
  fChain->SetBranchAddress("wMC", &_wMC, &b_wMC);
  fChain->SetBranchAddress("wSignal", &_wSignal, &b_wSignal);
  fChain->SetBranchAddress("wTot", &_wTot, &b_wTot);
  fChain->SetBranchAddress("muIPx", &_muIPx, &b_muIPx);
  fChain->SetBranchAddress("muIPy", &_muIPy, &b_muIPy);
  fChain->SetBranchAddress("muIPz", &_muIPz, &b_muIPz);
  fChain->SetBranchAddress("muIPsignificance", &_muIPsignificance, &b_muIPsignificance);
  fChain->SetBranchAddress("pvx", &_pvx, &b_pvx);
  fChain->SetBranchAddress("pvy", &_pvy, &b_pvy);
  fChain->SetBranchAddress("pvz", &_pvz, &b_pvz);
  fChain->SetBranchAddress("pvCov00", &_pvCov00, &b_pvCov00);
  fChain->SetBranchAddress("pvCov11", &_pvCov11, &b_pvCov11);
  fChain->SetBranchAddress("pvCov22", &_pvCov22, &b_pvCov22);
  fChain->SetBranchAddress("pvCov01", &_pvCov01, &b_pvCov01);
  fChain->SetBranchAddress("pvCov02", &_pvCov02, &b_pvCov02);
  fChain->SetBranchAddress("pvCov12", &_pvCov12, &b_pvCov12);
  if(Sys == "default") {
    fChain->SetBranchAddress("wPrefiringUp", &_wPrefiringUp, &b_wPrefiringUp);
    fChain->SetBranchAddress("wPrefiringDown", &_wPrefiringDown, &b_wPrefiringDown);
    fChain->SetBranchAddress("wIDvsJetUp", &_wIDvsJetUp, &b_wIDvsJetUp);
    fChain->SetBranchAddress("wIDvsJetDown", &_wIDvsJetDown, &b_wIDvsJetDown);
    fChain->SetBranchAddress("wIDvsEleUp", &_wIDvsEleUp, &b_wIDvsEleUp);
    fChain->SetBranchAddress("wIDvsEleDown", &_wIDvsEleDown, &b_wIDvsEleDown);
    fChain->SetBranchAddress("wIDvsMuUp", &_wIDvsMuUp, &b_wIDvsMuUp);
    fChain->SetBranchAddress("wIDvsMuDown", &_wIDvsMuDown, &b_wIDvsMuDown);
    fChain->SetBranchAddress("wTrgUp", &_wTrgUp, &b_wTrgUp);
    fChain->SetBranchAddress("wTrgDown", &_wTrgDown, &b_wTrgDown);
    fChain->SetBranchAddress("wZpTUp", &_wZpTUp, &b_wZpTUp);
    fChain->SetBranchAddress("wZpTDown", &_wZpTDown, &b_wZpTDown);
    fChain->SetBranchAddress("wToppTUp", &_wToppTUp, &b_wToppTUp);
    fChain->SetBranchAddress("wToppTDown", &_wToppTDown, &b_wToppTDown);
    fChain->SetBranchAddress("wBtagUp", &_wBtagUp, &b_wBtagUp);
    fChain->SetBranchAddress("wBtagDown", &_wBtagDown, &b_wBtagDown);
    fChain->SetBranchAddress("wPSISRUp", &_wPSISRUp, &b_wPSISRUp);
    fChain->SetBranchAddress("wPSISRDown", &_wPSISRDown, &b_wPSISRDown);
    fChain->SetBranchAddress("wPSFSRUp", &_wPSFSRUp, &b_wPSFSRUp);
    fChain->SetBranchAddress("wPSFSRDown", &_wPSFSRDown, &b_wPSFSRDown);
    fChain->SetBranchAddress("wScaleUp", &_wScaleUp, &b_wScaleUp);
    fChain->SetBranchAddress("wScaleDown", &_wScaleDown, &b_wScaleDown);
    fChain->SetBranchAddress("genTaupx", &_genTaupx, &b_genTaupx);
    fChain->SetBranchAddress("genTaupy", &_genTaupy, &b_genTaupy);
    fChain->SetBranchAddress("genTaupz", &_genTaupz, &b_genTaupz);
    fChain->SetBranchAddress("genTauE", &_genTauE, &b_genTauE);
    fChain->SetBranchAddress("genTauSVx", &_genTauSVx, &b_genTauSVx);
    fChain->SetBranchAddress("genTauSVy", &_genTauSVy, &b_genTauSVy);
    fChain->SetBranchAddress("genTauSVz", &_genTauSVz, &b_genTauSVz);
    fChain->SetBranchAddress("genMuonpx", &_genMuonpx, &b_genMuonpx);
    fChain->SetBranchAddress("genMuonpy", &_genMuonpy, &b_genMuonpy);
    fChain->SetBranchAddress("genMuonpz", &_genMuonpz, &b_genMuonpz);
    fChain->SetBranchAddress("genMuonE", &_genMuonE, &b_genMuonE);
    fChain->SetBranchAddress("genPVx", &_genPVx, &b_genPVx);
    fChain->SetBranchAddress("genPVy", &_genPVy, &b_genPVy);
    fChain->SetBranchAddress("genPVz", &_genPVz, &b_genPVz);
    fChain->SetBranchAddress("gendpPhiCP", &_gendpPhiCP, &b_gendpPhiCP);
    fChain->SetBranchAddress("genpvPhiCP", &_genpvPhiCP, &b_genpvPhiCP);
  } 
  Notify();
}

Bool_t NtupleReader::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void NtupleReader::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t NtupleReader::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef NtupleReader_cxx
