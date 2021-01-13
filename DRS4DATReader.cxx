#ifndef __DRS4DAT_READER__
#define __DRS4DAT_READER__

#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Reader.cxx"

/*=====================
Reader for DAT: DRS4 binary files
keep in mind that DRS4 stores in one .dat file as
many channels as there were saved. So the unpacking
is a bit tricky. I addapted the Reader to accomodate
for buffering several channels at the same time.
=====================*/

class DRS4DATReader : public Reader {
public:
  DRS4DATReader(TString,Int_t nsamples=1024);
  ~DRS4DATReader();
  void  ReadHeader();
  bool  ReadEvent();

private:
  Int_t fStart;
  Int_t fNumberOfChannels;
};

DRS4DATReader::DRS4DATReader(TString filename, Int_t nsamples) : Reader(filename) {
  fSamples = nsamples;
  fStart = 0;
  fNumberOfChannels = 0;
  cout << "DRS configured for " << nsamples << " samples" << endl << endl;
}
//=======
DRS4DATReader::~DRS4DATReader() {
}
//=======
void DRS4DATReader::ReadHeader() {
  fIFS.seekg(0);
  Char_t  c4Tmp[4];
  Char_t  c2Tmp[2];
  Int_t   iTmp;
  //RUN HEADER
  fIFS.read((char*) &c4Tmp, 4);
  //std::cout << c4Tmp << " | "; // DRSN
  fIFS.read((char*) &c4Tmp, 4);
  //std::cout << c4Tmp << endl; // TIME

  fNumberOfChannels = 0;
  for(;;) {
    int prev = fIFS.tellg();
    fIFS.read((char*) &c2Tmp, 2);
    if(c2Tmp[0] != 'B' || c2Tmp[1] != '#') {
      fIFS.seekg(prev);
      break;
    }
    //std::cout << c2Tmp << " | "; // B#
    fIFS.read((char*) &iTmp, 2);
    //std::cout << iTmp << std::endl; // BoardNumber
    //discovering channels
    for(int ich=0; ich!=5; ++ich) {
      prev = fIFS.tellg();
      fIFS.read( (char*) &c4Tmp, 4);
      //cout << fc4Tmp << endl;
      if(c4Tmp[0] != 'C') {
	fIFS.seekg(prev);
	break;
      }
      cout << "Found channel " << c4Tmp[0] << c4Tmp[3] << endl;
      float tch[fSamples];
      fIFS.read((char*) &tch, 4*fSamples); // Effective time bin width in ns encoded in 4-byte floating point format
      //quick and dirty
      float endtime = tch[0];
      for(int isa=1; isa!=fSamples; ++isa) {
	endtime += tch[isa];
      }
      if(fTrace[fNumberOfChannels]) delete fTrace[fNumberOfChannels];
      if(fAll[fNumberOfChannels]) delete fAll[fNumberOfChannels];
      fTrace[fNumberOfChannels] = new WaveForm(Form("DRS4CH%d",fNumberOfChannels),"DRS4;s;V", fSamples, tch[0], endtime);
      fAll[fNumberOfChannels] = new TH2D(Form("DRS4SumCH%d",fNumberOfChannels),"DRS4;s;V", 100, tch[0], endtime, 100, -50, +50);
      fNumberOfChannels++;
    }
  }
  fStart = fIFS.tellg();
  ResetReading();
}
//=======
bool DRS4DATReader::ReadEvent() {
  if(!fTrace[0]) {
    cout << "call DRS4DATReader::ReadHeader first " << endl;
    return false;
  }
  typedef struct {
    Char_t   control[4];
    UInt_t   serial_number;
    UShort_t year;
    UShort_t month;
    UShort_t day;
    UShort_t hour;
    UShort_t minute;
    UShort_t second;
    UShort_t milisecond;
    UShort_t range;
  } HEADERTYPE;
  HEADERTYPE eheader;
  Double_t mV, ns;
  char  c4Tmp[4];
  char  c2Tmp[2];
  int   iTmp;
  UShort_t ch[1024];

  //EVENT HEADER
  fIFS.read((char*) &eheader, sizeof(eheader));
  if ( (fIFS.rdstate() & std::ifstream::eofbit ) != 0 )
    return false;
  //cout << "SN " << eheader.serial_number << endl;
  
  for(int ich=0; ich<fNumberOfChannels;) {
    fIFS.read((char*) &c2Tmp, 2);
    //std::cout <<c2Tmp << " "; // B#
    fIFS.read((char*) &iTmp, 2);
    //std::cout << iTmp << " | "; // BoardNumber
    fIFS.read((char*) &c2Tmp, 2);
    //std::cout << c2Tmp << " "; // T#
    fIFS.read((char*) &iTmp, 2);
    //std::cout << iTmp << std::endl; // TriggerCell
    int tcell = iTmp;
    //Reading the values of each sample in the buffer for an event
    for(int i=0; i!=5; ++i) {
      int prev = fIFS.tellg();
      fIFS.read((char*) &c4Tmp, 4); // C00i
      //std::cout <<c4Tmp << " "; // C00i
      if(c4Tmp[0] != 'C') {
        fIFS.seekg(prev);
        break;
      }
      fIFS.read((char*) &iTmp, 4); // Scaler in Hz
      fIFS.read((char*) &ch, fSamples*2); // SAMPLES*2 (2-byte integer. 0=RC-0.5V and 65535=RC+0.5V. RC see header.
      //Converting raw value to mV
      //cout << ich << endl;
      for (int isa=0; isa!=fSamples; ++isa) {
	mV = 1000*(ch[isa]/65535.-0.5)+eheader.range;
	/*
	float timeval = 0;
	for(int j=0; j<isa-1; ++j)
	  timeval += tch[ich][ (j+tcell)%1024 ];
	axis[ich]->Fill( timeval, chmV[ich][isa] );
	*/
	fTrace[ich]->SetBinContent(isa+1,mV);
	ns = fTrace[ich]->GetXaxis()->GetBinCenter(isa+1);
	fAll[ich]->Fill(ns,mV);
      }
      ich++;
      if(ich>=fNumberOfChannels) break; // needed for last event (?)
    }
  }
  
  return true;
}

#endif
