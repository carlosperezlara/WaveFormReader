#ifndef __V1742DAT_READER__
#define __V1742DAT_READER__

#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Reader.h"

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
  ClassDef(V1742DATReader,0);
};

#endif
