//#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "ScaleFactor.h"

 void ScaleFactor::init_ScaleFactor(TString inputRootFile){
   //	std::cout<<"check 1  "<<std::endl;
 	TFile *fileIn = new TFile(inputRootFile, "READ");
 	// if root file not found
 	if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from NTupleMaker/src/ScaleFactor.cc : ‎File " <<inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); };
 	//	std::cout<<"check 2  "<<std::endl;
 	std::string HistoBaseName = "ZMass";
 	std::cout<<"check 21  "<<std::endl;
 	etaBinsH =(TH1D*)(fileIn->Get("etaBinsH"));
 	std::cout<<"check 215  "<<std::endl;
 	//	etaBinsH = etaBinsHTemp;
 	//	std::cout<<"check 22  "<<std::endl;

 	// ETrigIdEffFile = new TFile(basedir+"ElectronEfficiencies_Run2012ReReco_53X_Trig.root", "READ");
 	// ENonTrigIdEffFile = new TFile(basedir+"ElectronEfficiencies_Run2012ReReco_53X_NonTrig.root", "READ");
 	// ERecoEffFile = new TFile(basedir+"Electrons_ScaleFactors_Reco_8TeV.root", "READ");
 	// // load histograms
 	// ElectronTrigEff = (TH2D*)(ETrigIdEffFile->Get("electronsDATAMCratio_FO_ID_ISO"));
 	// ElectronNonTrigEff = (TH2D*)(ENonTrigIdEffFile->Get("h_electronScaleFactor_IdIsoSip"));


 	std::string etaLabel, GraphName;
 	//	std::cout<<"check 23  "<<std::endl;
 	int nEtaBins = etaBinsH->GetNbinsX();
 	std::cout<<"check 3  "<<std::endl;
  	for (int iBin=0; iBin<nEtaBins; iBin++){
	  std::cout<<"check 31  "<<std::endl;    
	  etaLabel = etaBinsH->GetXaxis()->GetBinLabel(iBin+1);
	  std::cout<<"check 32  "<<std::endl;
	  GraphName = HistoBaseName+etaLabel+"_Data";
	  std::cout<<"check 33  "<<std::endl;
	  eff_data[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName)); 
	  std::cout<<"check 34  "<<std::endl;
	  SetAxisBins(eff_data[etaLabel]);
	  std::cout<<"check 35  "<<std::endl;
	  GraphName = HistoBaseName+etaLabel+"_MC";
	  std::cout<<"check 36  "<<std::endl;
	  eff_mc[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
	  std::cout<<"check 4  "<<std::endl;
	  SetAxisBins(eff_mc[etaLabel]); 
	  bool sameBinning = check_SameBinning(eff_data[etaLabel], eff_mc[etaLabel]);
	  if (!sameBinning) {std::cout<< "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc . Can not proceed because ScaleFactor::check_SameBinning returned different pT binning for data and MC for eta label " << etaLabel << std::endl; exit(1); }; 
 	}
	
 	return;
 }





void ScaleFactor::init_ScaleFactor(TString inputRootFile, std::string HistoBaseName){

  TFile * fileIn = new TFile(inputRootFile, "read");
  // if root file not found                                                                                                                                                                          
  if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc : File " <<inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); };

  TH1D *etaBinsH = (TH1D*)fileIn->Get("etaBinsH");

  std::string etaLabel, GraphName;
  int nEtaBins = etaBinsH->GetNbinsX();
  for (int iBin=0; iBin<nEtaBins; iBin++){
    etaLabel = etaBinsH->GetXaxis()->GetBinLabel(iBin+1);
    GraphName = HistoBaseName+etaLabel+"_Data";
    eff_data[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
    SetAxisBins(eff_data[etaLabel]);
    GraphName = HistoBaseName+etaLabel+"_MC";
    eff_mc[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
    SetAxisBins(eff_mc[etaLabel]);
    bool sameBinning = check_SameBinning(eff_data[etaLabel], eff_mc[etaLabel]);
    if (!sameBinning) {std::cout<< "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc . Can not proceed because ScaleFactor::check_SameBinning returned different pT binning for data and MC for eta label " << etaLabel << std::endl; exit(1); };
  }

  return;
}



void ScaleFactor::SetAxisBins(TGraphAsymmErrors* graph) {

	int NPOINTS = graph->GetN(); 
	double AXISBINS[NPOINTS+1];
	for (int i=0; i<NPOINTS; i++) { AXISBINS[i] = (graph->GetX()[i] - graph->GetErrorXlow(i)); }
	AXISBINS[NPOINTS] = (graph->GetX()[NPOINTS-1] + graph->GetErrorXhigh(NPOINTS-1));
	graph->GetXaxis()->Set(NPOINTS, AXISBINS);
	return;
}

bool ScaleFactor::check_SameBinning(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2){
	bool haveSameBins = false;
	int n1 = graph1->GetXaxis()->GetNbins();
	int n2 = graph2->GetXaxis()->GetNbins();
	if (n1 != n2 ) {return false;}
	else {
		haveSameBins = true;
		const int nbins = n1;
		double x1, x2;
		for (int i=0; i<nbins; i++){ 
			x1 = (graph1->GetXaxis()->GetXbins())->GetArray()[i];
			x2 = (graph2->GetXaxis()->GetXbins())->GetArray()[i]; 
			haveSameBins = haveSameBins and (x1== x2) ;
		}
	}

	return haveSameBins;
}


std::string ScaleFactor::FindEtaLabel(double Eta){
	Eta = fabs(Eta);
	int binNumber = etaBinsH->GetXaxis()->FindFixBin(Eta);
	std::string EtaLabel = etaBinsH->GetXaxis()->GetBinLabel(binNumber);
	std::map<std::string, TGraphAsymmErrors*>::iterator it;
	it =  eff_data.find(EtaLabel);
	if ( it == eff_data.end()) { 
	std::cout << "ERROR in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : no object corresponding to eta label "<< EtaLabel << " for data " << std::endl; exit(1);
	}
	else return EtaLabel;
}


int ScaleFactor::FindPtBin( std::map<std::string, TGraphAsymmErrors *> eff_map, std::string EtaLabel, double Pt){

        int Npoints = eff_map[EtaLabel]->GetN();
	double ptMAX = (eff_map[EtaLabel]->GetX()[Npoints-1])+(eff_map[EtaLabel]->GetErrorXhigh(Npoints-1));
	double ptMIN = (eff_map[EtaLabel]->GetX()[0])-(eff_map[EtaLabel]->GetErrorXlow(0));
	// if pt is overflow, return last pt bin
 	if (Pt >= ptMAX ) return Npoints; 
	// if pt is underflow, return nonsense number and warning
	else if (Pt < ptMIN){ 
 	std::cout<< "WARNING in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: pT too low (pt = " << Pt << "), min value is " << ptMIN << ". Returned efficiency =1. Weight will be 1. " << std::endl;
	return -99;}
	// if pt is in range
	else {return eff_map[EtaLabel]->GetXaxis()->FindFixBin(Pt);} 
	}


double ScaleFactor::get_EfficiencyData(double pt, double eta){

        double eff;
	std::string label = FindEtaLabel(eta);

	int ptbin = FindPtBin(eff_data, label, pt); 
	if (ptbin == -99){eff =1;} // if pt is underflow 
	else eff = eff_data[label]->GetY()[ptbin-1];

	if (eff > 1.) {std::cout<< "WARNING in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: Returned efficiency in data > 1. " << std::endl;} 
	if (eff < 0 ) {std::cout<<"WARNING in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: Returned negative efficiency in data" <<std::endl;}

	return eff;
	
}


double ScaleFactor::get_EfficiencyMC(double pt, double eta) {

	double eff;		
	std::string label = FindEtaLabel(eta);

	int ptbin = FindPtBin(eff_mc, label, pt); 
	if (ptbin == -99){eff =1;} // if pt is underflow 
	else eff= eff_mc[label]->GetY()[ptbin-1];

	if (eff > 1. ) {std::cout << "WARNING in ScaleFactor::get_EfficiencyMC(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : Returned efficiency in MC > 1. " << std::endl;} 		
	if (eff < 0 ) {std::cout<<"WARNING in ScaleFactor::get_EfficiencyMC(double pt, double eta) from LepEffIntrface/src/ScaleFactor.cc : Returned negative efficiency in MC. " <<std::endl;}
	

	return eff;

}



double ScaleFactor::get_ScaleFactor(double pt, double eta){
	
	double efficiency_data = get_EfficiencyData(pt, eta);
	double efficiency_mc = get_EfficiencyMC(pt, eta);
	double SF;

	if ( efficiency_mc != 0) {SF = efficiency_data/efficiency_mc;}
	else {
	SF=0.; std::cout << "WARNING in ScaleFactor::get_ScaleFactor(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : MC efficiency = 0. Scale Factor set to 0. ";
	}

	return SF;	
	
}


double ScaleFactor::get_EfficiencyDataError(double pt, double eta){

	double eff_error;
	std::string label = FindEtaLabel(eta);
	int ptbin = FindPtBin(eff_data, label, pt); 
	if (ptbin == -99){eff_error =0.;} // if pt is underflow 
	else eff_error= eff_data[label]->GetErrorYhigh(ptbin-1); 
        // errors are supposed to be symmetric, can use GetErrorYhigh or GetErrorYlow

	double effData = get_EfficiencyData(pt,eta);
	if (eff_error > effData) eff_error = 0.5*effData;
	return eff_error;
}
	
	

double ScaleFactor::get_EfficiencyMCError(double pt, double eta){

	double eff_error;
	std::string label = FindEtaLabel(eta);
	int ptbin = FindPtBin(eff_mc, label, pt); 
	if (ptbin == -99){eff_error =0.;} // if pt is underflow 
	else eff_error= eff_mc[label]->GetErrorYhigh(ptbin-1); 
	// errors are supposed to be symmetric, can use GetErrorYhigh or GetErrorYlow

	double effMC = get_EfficiencyMC(pt,eta);
	if (eff_error > effMC ) eff_error = 0.5*effMC;
	return eff_error;
}

double ScaleFactor::get_ScaleFactorError(double pt, double eta){

	double SF_error = 0.;
	
	double effData = get_EfficiencyData(pt, eta);
	double effMC = get_EfficiencyMC(pt, eta);
	double errData = get_EfficiencyDataError(pt, eta);
	double errMC =  get_EfficiencyMCError(pt, eta);

	if (errData == 0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: uncertainty on data point = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (errMC ==0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: uncertainty on MC = 0, can not calculate uncerttainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (effData ==0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: efficiency in data = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (effMC ==0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: efficiency in MC = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	else {	
	SF_error = pow((errData/effData),2) + pow((errMC/effMC),2);
	SF_error = pow(SF_error, 0.5)*(effData/effMC);
	}
	return SF_error;
}
	

