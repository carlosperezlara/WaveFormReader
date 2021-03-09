#ifndef __V1742DAT_READER__
#define __V1742DAT_READER__

#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Reader.cxx"

/*=====================
Reader for DAT: DRS4 in CAEN V1742 binary files
keep in mind that DRS4 stores in one .dat file as
many channels as there were saved. So the unpacking
is a bit tricky. I addapted the Reader to accomodate
for buffering several channels at the same time.
=====================*/

class V1742DATReader : public Reader {
public:
  V1742DATReader(TString,Double_t min=-500,Double_t max=+500);
  ~V1742DATReader();
  void  ReadHeader();
  bool  ReadEvent();

private:
  bool ReadGroup(int group);
  Int_t fNumberOfChannels;
  UInt_t fGroupMask;
};

V1742DATReader::V1742DATReader(TString filename, Double_t min, Double_t max) : Reader(filename,min,max) {
  fNumberOfChannels = 0;
  fGroupMask = 0b1111;
}
//=======
V1742DATReader::~V1742DATReader() {
}
//=======
void V1742DATReader::ReadHeader() {
  fIFS.seekg(0);
  fStart = fIFS.tellg();

  fNumberOfChannels = 32;
  fSamples = 1024;
  Double_t minX = -0.5;
  Double_t maxX = 1023.5;
  Double_t minY = fMinSummaryRange;
  Double_t maxY = fMaxSummaryRange;
  if(fMaxSummaryRange<fMinSummaryRange) {
    minY = -500;
    maxY = +500;
  }
  for(int ich=0;ich!=32;++ich) {
    if(fTrace[ich]) delete fTrace[ich];
    if(fAll[ich]) delete fAll[ich];
    fTrace[ich] = new WaveForm(Form("Trace_CH%d_%s",ich,fFileName.Data()),
			       Form("Trace  CH%d  %s;s;V",ich,fFileName.Data()),
			       fSamples, minX, maxX);
    fAll[ich] = new TH2D(Form("Summary_CH%d_%s",ich,fFileName.Data()),
			 Form("Summary  CH%d  %s;s;V",ich,fFileName.Data()),
			 fSamples, minX, maxX, 100, minY, maxY);
  }
  ResetReading();
}
//=======
bool V1742DATReader::ReadEvent() {
  UInt_t  uint;
  //RUN HEADER
  fIFS.read((char*) &uint, 4);
  //std::cout << "First word is " << uint << std::endl;
  UInt_t INIT = uint>>28;
  UInt_t TOTAL_EVENT_SIZE = (uint << 4)>>4;
  //std::cout << " INIT: " << INIT << (INIT==0b1010?" GOOD":" OOPPS") << std::endl;
  //std::cout << " TOTAL EVENT SIZE: " << TOTAL_EVENT_SIZE << std::endl;

  fIFS.read((char*) &uint, 4);
  //std::cout << "Second word is " << uint << std::endl;
  //UInt_t BOARDID = uint>>27;
  UInt_t BF = (uint>>26) & 0b1;
  //UInt_t RES1 = (uint>>24) & 0b11;
  //UInt_t PATTERN = ((uint<<8)>>16);
  //UInt_t RES2 = ((uint<<24)>>28);
  UInt_t GROUPMASK = uint & 0b1111;
  //std::cout << " BOARD ID: " << BOARDID << std::endl;
  //std::cout << " BF: " << BF << (BF!=0?"    W A R N I N G ! ! !    BOARD FAILED FLAG FIRED":"") << std::endl;
  //std::cout << " RESERVED: " << RES1 << std::endl;
  //std::cout << " PATTERN: " << PATTERN << std::endl;
  //std::cout << " RESERVED: " << RES2 << std::endl;
  //std::cout << " GROUP MASK: " << GROUPMASK << std::endl;
  fGroupMask = GROUPMASK;
  
  fIFS.read((char*) &uint, 4);
  //std::cout << "Third word is " << uint << std::endl;
  //UInt_t RES3 = uint>>24;
  UInt_t EVENTCOUNTER = (uint<<8)>>8;
  //std::cout << " RESERVED: " << RES3 << std::endl;
  //std::cout << " EVENT COUNTER: " << EVENTCOUNTER << std::endl;

  fIFS.read((char*) &uint, 4);
  //std::cout << "Fourth word is " << uint << std::endl;

  bool allGood = true;
  if(fGroupMask & 0b0001)  allGood = allGood && ReadGroup(0);
  if(fGroupMask & 0b0010)  allGood = allGood && ReadGroup(1);
  if(fGroupMask & 0b0100)  allGood = allGood && ReadGroup(2);
  if(fGroupMask & 0b1000)  allGood = allGood && ReadGroup(3);
  return allGood;
}
//=======
bool V1742DATReader::ReadGroup(int iGroup) {
  UInt_t  uint;
  fIFS.read((char*) &uint, 4);
  UInt_t CONTROL1 = uint>>30;
  //UInt_t STARTINDEXCELL = (uint<<2)>>22;
  UInt_t CONTROL2 = (uint>>18) & 0b11;
  UInt_t CONTROL3 = (uint>>13) & 0b111;
  UInt_t FREQ = (uint>>16) & 0b11;
  UInt_t TR = (uint>>12) & 0b1;
  //UInt_t SIZE = (uint<<20)>>20;
  //std::cout << std::endl;
  //cout << " READING GROUP " << iGroup << endl;
  //std::cout << " START INDEX CELL: " << STARTINDEXCELL << std::endl;
  //std::cout << " FREQ: " << FREQ << " { 5 GS/s, 2.5 GS/s, 1 GS/s, 750 MS/s } " << std::endl;
  //std::cout << " TR: " << TR << (TR==1?" (PRESENT)":"") << std::endl;
  //std::cout << " SIZE: " << SIZE << std::endl;
  //std::cout << " CONTROL1: " << CONTROL1 << std::endl;
  //std::cout << " CONTROL2: " << CONTROL2 << std::endl;
  //std::cout << " CONTROL3: " << CONTROL3 << std::endl;
  if((CONTROL1+CONTROL2+CONTROL3)!=0) {
    return false;
  }
  UInt_t buffer[3];
  for(int isa=0; isa!=1024; ++isa) {
    fIFS.read((char*) &buffer[2], 4);
    fIFS.read((char*) &buffer[1], 4);
    fIFS.read((char*) &buffer[0], 4);
    //cout << "ISA " << isa << endl;
    for(int ich=0; ich!=8; ++ich) {
      UInt_t adc = buffer[0] & 0b111111111111;
      Double_t mV = 1000*(adc/4095.-0.5)+0;
      fTrace[ich+(iGroup*8)]->SetBinContent(isa+1,mV);
      fAll[ich+(iGroup*8)]->Fill(isa,mV);
      //cout << " " << ich << ": " << adc << " " << mV << endl;
      buffer[0] = buffer[0] >> 12;
    }
  }
  if(TR) {
    for(int ichunk=0; ichunk!=128; ++ichunk) { //128*8 = 1024
      fIFS.read((char*) &buffer[2], 4);
      fIFS.read((char*) &buffer[1], 4);
      fIFS.read((char*) &buffer[0], 4);
      //cout << "CHUNK " << ichunk << endl;
      for(int isa=0; isa!=8; ++isa) {
	UInt_t adc = buffer[0] & 0b111111111111;
	//cout << " " << isa << ": " << adc << endl;
	buffer[0] = buffer[0] >> 12;
      }
    }
  }
  
  fIFS.read((char*) &uint, 4);
  //std::cout << " GROUP TRIGGER TIME TAG: " << uint << std::endl;

  return true;
}

  /*  
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
	float timeval = 0;
	for(int j=0; j<isa-1; ++j)
	  timeval += tch[ich][ (j+tcell)%1024 ];
	axis[ich]->Fill( timeval, chmV[ich][isa] );
	fTrace[ich]->SetBinContent(isa+1,mV);
	ns = fTrace[ich]->GetXaxis()->GetBinCenter(isa+1);
	fAll[ich]->Fill(ns,mV);
      }
      fTrace[ich]->NonFitted();
      ich++;
      if(ich>=fNumberOfChannels) break; // needed for last event (?)
    }
  }
  
  return true;
}
*/
#endif
