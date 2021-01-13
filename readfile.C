#include "LECROYTRCReader.cxx"

int readfile() {
  LECROYTRCReader file("C2--JanBoard-HG-44V--00000.trc");
  file.ReadHeader();

  TH1D *trace = file.GetTrace();
  TH2D *summary = file.GetSummaryPlot();
  for(int nev=0;nev!=6;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    // i have access through trace
  }
  new TCanvas();
  trace->DrawCopy();
  new TCanvas();
  summary->DrawCopy("colz");
  return 0;
}
