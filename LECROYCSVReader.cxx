#ifndef __LECROYCSV_READER__
#define __LECROYCSV_READER__

#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Reader.cxx"
#include "WaveForm.cxx"

/*=====================
Reader for TRC: Lecroy binary files
keep in mind that Lecroy OSC stores one .trc
file per channel even if there were many channels
saved at the same time. In order to read/correlate
several channels: 1) create as many LECROYCSVReader objects
as channels there are. 2) invoke the ReadHeader for
all of them. 3) Create a loop in which you ReadEvent
for each object. 4) Voila
=====================*/

class LECROYCSVReader : public Reader {
public:
  LECROYCSVReader(TString);
  ~LECROYCSVReader();
  void  ReadHeader();
  bool  ReadEvent();

private:
  Int_t    fNumberOfEvents;
};

LECROYCSVReader::LECROYCSVReader(TString filename) : Reader(filename) {
  fNumberOfEvents = 0;
}
//=======
LECROYCSVReader::~LECROYCSVReader() {
}
//=======
void LECROYCSVReader::ReadHeader() {
  fIFS.seekg(0);
  std::string line;
  std::getline(fIFS,line);
  std::getline(fIFS,line);
  TString sLine = line;
  TString sMAX_WAVES = ((TObjString*) sLine.Tokenize(",")->At(1))->GetString();
  TString sSAMPLES   = ((TObjString*) sLine.Tokenize(",")->At(3))->GetString();
  std::getline(fIFS,line);
  fNumberOfEvents = sMAX_WAVES.Atoi();
  fSamples = sSAMPLES.Atoi();
  for(int n=0;n!=fNumberOfEvents;++n) {
    if(!std::getline(fIFS,line)) break;
  }
  for(int n=0;n!=1;++n) {
    if(!std::getline(fIFS,line)) break;
  }
  fStart = fIFS.tellg();
  TString lineFirst, lineLast;
  for(int y=0;y!=fSamples;++y) {
    std::getline(fIFS,line);
    if(y==0) lineFirst = line;
    if(y==fSamples-1) lineLast = line;
  }
  lineFirst = ((TObjString*) lineFirst.Tokenize("'")->At(0))->GetString();
  lineLast  = ((TObjString*) lineLast.Tokenize("'")->At(0))->GetString();
  double minX = lineFirst.Atof()*1e9;
  double maxX = lineLast.Atof()*1e9;
  if(fTrace[0]) delete fTrace[0];
  if(fAll[0]) delete fAll[0];
  fTrace[0] = new WaveForm("LECROYCH","LECROY;ns;mV", fSamples, minX, maxX );
  fAll[0] = new TH2D("LECROYSumCH","LECROY;ns;mV",100,minX,maxX, 100, -50, 50);
  ResetReading();
  cout << "NumberOfEvents " << fNumberOfEvents << endl;
}
//=======
bool LECROYCSVReader::ReadEvent() {
  if(!fTrace[0]) {
    cout << "call LECROYCSVReader::ReadHeader first " << endl;
    return false;
  }
  std::string line;
  Char_t  charsize1;   // byte
  Short_t shortsize2;  // int16 word
  Double_t V,s;
  for(int i=0; i!=fSamples; ++i) {
    std::getline(fIFS,line);
    TString rline = line;
    TObjArray *lst = rline.Tokenize(",");
    double ix = ((TObjString*) lst->At(0))->GetString().Atof()*1e9; // nS
    double iv = ((TObjString*) lst->At(1))->GetString().Atof()*1e3; // mV
    if ( (fIFS.rdstate() & std::ifstream::eofbit ) != 0 )
      return false;
    fTrace[0]->SetBinContent(i+1,iv);
    fAll[0]->Fill( ix, iv );
  }
  fTrace[0]->NonFitted();
  return true;
}

#endif
