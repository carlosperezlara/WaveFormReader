#ifndef __LECROYTRC_READER__
#define __LECROYTRC_READER__

#include "TString.h"
#include "Reader.h"

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

class LECROYTRCReader : public Reader {
public:
  LECROYTRCReader(TString,Double_t min=+1, Double_t max=-1);
  ~LECROYTRCReader();
  void  ReadHeader();
  bool  ReadEvent();

private:
  Int_t    fNumberOfEvents;
  Int_t    fSize;
  Float_t  fGain;
  Float_t  fOffset;
  ClassDef(LECROYTRCReader,0);
};

#endif
