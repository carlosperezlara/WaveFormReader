#include "LECROYTRCReader.cxx"
#include "LECROYCSVReader.cxx"
#include "DRS4DATReader.cxx"

int reader() {
  LECROYTRCReader file("C2--2e14-40.15V-n30-Green--00000.trc");
  //LECROYCSVReader file("C3--scope_45VHighGain.txt");
  //DRS4DATReader file("same.dat");
  file.ReadHeader();

  TH1D *trace = file.GetTrace();
  TH2D *summary = file.GetSummaryPlot();
  //for(int nev=0;nev!=2000;++nev) {
  for(int nev=0;;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    // i have access through trace
  }
  //trace->DrawCopy();
  summary->DrawCopy("colz");
  return 0;
}
