#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Reader.h"

//=======
Reader::Reader(TString filename,Double_t minmV, Double_t maxmV) {
  fFileName = filename;
  fIFS.close();
  fIFS.open( fFileName.Data() );
  fStart = 0;
  fSamples = 0;
  fMinSummaryRange = minmV;
  fMaxSummaryRange = maxmV;
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
