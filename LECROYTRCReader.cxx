#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "LECROYTRCReader.h"
#include "WaveForm.h"

/*=====================
Reader for TRC: Lecroy binary files
keep in mind that Lecroy OSC stores one .trc
file per channel even if there were many channels
saved at the same time. In order to read/correlate
several channels: 1) create as many LECROYTRCReader objects
as channels there are. 2) invoke the ReadHeader for
all of them. 3) Create a loop in which you ReadEvent
for each object. 4) Voila
=====================*/

LECROYTRCReader::LECROYTRCReader(TString filename,Double_t min, Double_t max) : Reader(filename,min,max) {
  fNumberOfEvents = 0;
  fSize = 0;
  fGain = 0;
  fOffset = 0;
}
//=======
LECROYTRCReader::~LECROYTRCReader() {
}
//=======
void LECROYTRCReader::ReadHeader() {
  fIFS.seekg(0);
  char  *c50Tmp = new char [50];
  char  *c16Tmp = new char [16];
  char  *c12Tmp = new char [12];
  UChar_t   charsize1;   // byte
  UShort_t  ushortsize2; // int16 word
  Short_t   shortsize2;  // int16 word
  UInt_t    uintsize4;   // int32
  Int_t     intsize4;    // int32
  Float_t   floatsize4;  // float
  Double_t  doublesize8; // double
  //fIFS.read(c50Tmp, 50);
  fStart = 11;
  Int_t start=11;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+16); //27
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read(c16Tmp, 16); // template name
  //std::cout << "template name " << c16Tmp << std::endl;
  fIFS.read((char*) &ushortsize2, 2); // 8 or 16 bit
  fSize = ushortsize2; // 0:byte 1:word(2bytes)
  //std::cout << "size {8,16}" << ushortsize2 << std::endl;
  fIFS.read((char*) &ushortsize2, 2); // endianness
  //std::cout << "endianness {>,<} "<< ushortsize2 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "wave descriptor " << uintsize4 << std::endl;
  fStart += uintsize4;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "usertext " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "resdesc1 " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); // 
  //std::cout << "trigtimearray " << uintsize4 << std::endl;
  fStart += uintsize4;
  fIFS.read((char*) &uintsize4, 4); // 
  //std::cout << "ristimearray " << uintsize4 << std::endl;
  fStart += uintsize4;
  fIFS.read((char*) &uintsize4, 4); // 
  //std::cout << "resarray1 " << uintsize4 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+60); //12
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &uintsize4, 4);
  int samplesize = uintsize4;
  //std::cout << "wavearray01 " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4);
  //std::cout << "wavearray02 " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4);
  //std::cout << "resarray2 " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4);
  //std::cout << "resarray3 " << uintsize4 << std::endl;
  fIFS.read(c16Tmp, 16); // instrument name
  std::cout << "instrument name " << c16Tmp << std::endl;
  fIFS.read((char*) &uintsize4, 4); // instrument number
  std::cout << "instrument number " << uintsize4 << std::endl;
  fIFS.read(c16Tmp, 16); //
  //std::cout << "trace label " << c16Tmp << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "4bytes reserved " << uintsize4 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+116);//20
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &uintsize4, 4); // wave array count
  //std::cout << "wavearraycount " << uintsize4 << std::endl;
  int wavearraycount = uintsize4;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "pntsPerScreen " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "firstvalidpnt " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "lastvalidpnt " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "firstpnt " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "sparsingfactor " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "segmentindex " << uintsize4 << std::endl;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "subarraycount " << uintsize4 << std::endl;
  fNumberOfEvents = uintsize4;
  fIFS.read((char*) &uintsize4, 4); //
  //std::cout << "sweepsperacq " << uintsize4 << std::endl;
  fIFS.read((char*) &ushortsize2, 2); //
  //std::cout << "pointsperpair " << ushortsize2 << std::endl;
  fIFS.read((char*) &ushortsize2, 2); //
  //std::cout << "pairoffset " << ushortsize2 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+156);//36
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &floatsize4, 4); //
  fGain = floatsize4;
  //std::cout << "vertical gain " << floatsize4 << std::endl;
  fIFS.read((char*) &floatsize4, 4); //
  fOffset = floatsize4;
  //std::cout << "vertical offset " << floatsize4 << std::endl;
  fIFS.read((char*) &floatsize4, 4); //
  double maxV = floatsize4*fGain;
  //std::cout << "maxvalue " << floatsize4 << std::endl;
  fIFS.read((char*) &floatsize4, 4); //
  double minV = floatsize4*fGain;
  //std::cout << "minvalue " << floatsize4 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+172);//8
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &ushortsize2, 2); // nominalBits
  //std::cout << "nominalBits " << ushortsize2 << std::endl;
  fIFS.read((char*) &ushortsize2, 2); //
  //std::cout << "nomsubarraycount " << ushortsize2 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+176);//2
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &floatsize4, 4); // horizInterval
  float horizinterval = floatsize4;
  //std::cout << "horizinterval " << floatsize4 << std::endl;
  fIFS.read((char*) &doublesize8, 8); // horizOffset
  float horizoffset = floatsize4;
  //std::cout << "horizoffset " << doublesize8 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+316);//128
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &ushortsize2, 2); // 
  //std::cout << "recordtypelist  {single_sweep,interleaved,histogram,graph,filter_coefficient,coplex,extrema,sequence_obsolote} " << ushortsize2 << std::endl;
  fIFS.read((char*) &ushortsize2, 2); // 
  //std::cout << "processingdone  {off on} " << ushortsize2 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+326);//6
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &ushortsize2, 2); // 
  std::cout << "vertical coupling {dc50, gnd, dc1m, gnd, ac1m} " << ushortsize2 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.seekg(start+344);//16
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fIFS.read((char*) &ushortsize2, 2); // source {ch1,ch2,ch3,ch4,unk}
  std::cout << "source {ch1, ch2, ch3, ch4, unk} " << ushortsize2 << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  //std::cout << "pos: " << fIFS.tellg() << std::endl;
  fSamples = wavearraycount/fNumberOfEvents;
  std::cout << "samples: " << fSamples << " || waves: " << fNumberOfEvents << std::endl;
  if(fTrace[0]) delete fTrace[0];
  if(fAll[0]) delete fAll[0];
  horizoffset *=1e9;
  horizinterval *=1e9;
  minV *= 1e3;
  maxV *= 1e3;
  if(fMinSummaryRange<fMaxSummaryRange) {
    minV = fMinSummaryRange;
    maxV = fMaxSummaryRange;
  }
  fTrace[0] = new WaveForm(Form("fTrace0_%s",fFileName.Data()),
			   Form("Trace0  %s;ns;mV",fFileName.Data()),fSamples,horizoffset,fSamples*horizinterval+horizoffset);
  fAll[0] = new TH2D(Form("fSummary0_%s",fFileName.Data()),
		     Form("Summary0  %s;ns;mV",fFileName.Data()),
		     100,horizoffset,fSamples*horizinterval + horizoffset, 100, minV, maxV);
  ResetReading();
}
//=======
bool LECROYTRCReader::ReadEvent() {
  if(!fTrace[0]) {
    std::cout << "call LECROYTRCReader::ReadHeader first " << std::endl;
    return false;
  }
  Char_t  charsize1;   // byte
  Short_t shortsize2;  // int16 word
  Double_t mV,ns;
  for(int i=0; i!=fSamples; ++i) {
    if( fSize==1 ) { //16bit packet
      fIFS.read((char*) &shortsize2, 2); // pack
      mV = shortsize2*fGain + fOffset;
    } else { // 8bit packet
      fIFS.read((char*) &charsize1, 1); // pack
      mV = charsize1*fGain + fOffset;
    }
    mV *= 1e3;
    if ( (fIFS.rdstate() & std::ifstream::eofbit ) != 0 )
      return false;
    //std::cout << i << " "<< shortsize2 << " " << mV << std::endl;
    fTrace[0]->SetBinContent(i+1,mV);
    ns = fTrace[0]->GetXaxis()->GetBinCenter( i+1 );
    fAll[0]->Fill(ns,mV);
  }
  fTrace[0]->NonFitted();
  return true;
}
