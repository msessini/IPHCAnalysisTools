#ifndef BDTClassification_h
#define BDTClassification_h

#include <string>
#include "boost/format.hpp"
#include "boost/bind.hpp"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class BDTClassification{
 private:
  
  //int year_;
  float muPt_;
  float tauPt_;
  double ditauPt_;
  int Njets_;
  double leadingjetPt_;
  double subleadingjetPt_;
  double dijetPt_;
  double dijetMass_;
  double dijetdeltaEta_;
  double pairvisMass_;
  double fastMTTmass_;
  float muMETmt_;
  float PUPPImet_;

  TMVA::Reader *reader_even_;
  TMVA::Reader *reader_odd_;

 public:
  BDTClassification();
  virtual ~BDTClassification();
  virtual std::vector<float> read_mva_scores(unsigned isEven, std::vector<float> vars);
  virtual std::pair<float,int> getMaxScoreWithIndex(std::vector<float> vec);
  virtual int PreAnalysis();
  virtual int Execute(float muPt,float tauPt,double ditauPt,int Njets,double leadingjetPt,double subleadingjetPt,double dijetPt,double dijetMass,double dijetdeltaEta,double pairvisMass,double fastMTTmass,float muMETmt,float PUPPImet, unsigned long long evt_,std::vector<float> &score, std::pair<float, int> &max_pair);

  //int year_;
  unsigned isEven_;
  float event_;
  unsigned long long evt_;

  float var0_, var1_, var2_, var3_, var4_, var5_, var6_, var7_, var8_, var9_, var10_, var11_, var12_;

};

#endif
