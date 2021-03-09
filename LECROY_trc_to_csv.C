#include "WaveForm.cc"
#include "Reader.cc"
#include "LECROYTRCReader.cc"

int LECROY_trc_to_csv(TString filesource="C1--47V-8.0--00000.trc") {
  TString filetarget = Form("%s.csv",filesource.Data());
  LECROYTRCReader file(filesource.Data());
  ofstream fout( filetarget.Data() );
  file.ReadHeader();
  TH1D *trace = file.GetTrace();
  TH2D *summary = file.GetSummaryPlot();
  int bins = trace->GetNbinsX();
  int nwaves = 100;
  fout << Form("LECROYWR640Zi,63891,Waveform") << endl;
  fout << Form("Segments,%d,SegmentSize,%d",nwaves,bins) << endl;
  fout << Form("Segment,TrigTime,TimeSinceSegment1") << endl;
  for(int i=0; i!=nwaves; ++i)
    fout << Form("#%d,25-Ene-2012 01:30:00,0",i+1) << endl; // per wave
  fout << "Time,Ampl" << endl;
  for(int nev=0;nev!=6;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      fout << "Events read so far: " << nev << endl;
    // i have access through trace
    for(int b=0; b!=bins; ++b) {
      double time = trace->GetBinCenter(b+1) * 1e-9; // ns -> s
      double ampl = trace->GetBinContent(b+1) * 1e-3; // mV -> V
      fout << Form("%e, %e",time,ampl) << endl;
    }
  }
  fout.close();
  new TCanvas();
  summary->DrawCopy("colz");
  cout << "File " << filetarget.Data() << " was created " << endl;
  return 0;
}
