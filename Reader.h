#ifndef __VRT_READER__
#define __VRT_READER__

#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "WaveForm.h"

const Int_t kMAXCHANNELS=100;

class Reader {
public:
  Reader(TString,Double_t min=+1, Double_t max=-1);
  virtual ~Reader();
  virtual void  ReadHeader()=0;
  virtual bool  ReadEvent()=0;
  void  ResetReading();
  inline WaveForm* GetTrace(Int_t i=0) {return i<kMAXCHANNELS?fTrace[i]:NULL;}
  inline TH2D* GetSummaryPlot(Int_t i=0) {return i<kMAXCHANNELS?fAll[i]:NULL;}

protected:
  TString  fFileName;
  std::ifstream fIFS;
  Int_t    fStart;
  Int_t    fSamples;
  Double_t fMinSummaryRange;
  Double_t fMaxSummaryRange;

  WaveForm *fTrace[kMAXCHANNELS];
  TH2D     *fAll[kMAXCHANNELS];

  ClassDef(Reader,0);
};
#endif
