#ifndef __VRT_READER__
#define __VRT_READER__

#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"


const Int_t kMAXCHANNELS=4;

class Reader {
public:
  Reader(TString);
  ~Reader();
  virtual void  ReadHeader()=0;
  virtual bool  ReadEvent()=0;
  void  ResetReading();
  TH1D* GetTrace(Int_t i=0) {return i<kMAXCHANNELS?fTrace[i]:NULL;}
  TH2D* GetSummaryPlot(Int_t i=0) {return i<kMAXCHANNELS?fAll[i]:NULL;}

protected:
  TString  fFileName;
  ifstream fIFS;
  Int_t    fStart;
  Int_t    fSamples;

  TH1D    *fTrace[kMAXCHANNELS];
  TH2D    *fAll[kMAXCHANNELS];
};
//=======
Reader::Reader(TString filename) {
  fFileName = filename;
  fIFS.close();
  fIFS.open( fFileName.Data() );
  fStart = 0;
  fSamples = 0;
  for(int i=0; i!=kMAXCHANNELS; ++i) {
    fTrace[i] = NULL;
    fAll[i] = NULL;
  }
}
//=======
Reader::~Reader() {
  fIFS.close();
  for(int i=0; i!=kMAXCHANNELS; ++i) {
    if(fAll[i]) delete fAll[i];
    if(fTrace[i]) delete fTrace[i];
  }
}
//=======
void Reader::ResetReading() {
  fIFS.seekg(fStart);
}

#endif
