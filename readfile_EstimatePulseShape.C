#include "LECROYTRCReader.cxx"
#include "WaveForm.cxx"

int readfile_EstimatePulseShape() {
  LECROYTRCReader file("C2--JanBoard-HG-44V--00000.trc");
  file.ReadHeader();

  WaveForm *trace = file.GetTrace();
  TH2D *summary = file.GetSummaryPlot();
  Int_t bins = trace->GetXaxis()->GetNbins();
  Double_t bins_per_ns = bins / ( trace->GetXaxis()->GetBinLowEdge(bins+1) - trace->GetXaxis()->GetBinLowEdge(1) );
  cout << "Bins per ns " << bins_per_ns << endl;
  
  Double_t meanPed, rmsPed;
  TH1D *hMeanPed = new TH1D("hMeanPed","MeanPed;mV",100,-50,+50);
  TH1D *hRMSPed = new TH1D("hRMSPed","RMSPed;mV",100,0,+15);

  Double_t min_mV;
  Int_t min_ns;
  TH1D *hMinmV = new TH1D("hMinmV","MinmV;mV",100,-100,+100);
  TH1D *hMinns = new TH1D("hMinns","MinnS;ns",100,-20,+200);

  for(int nev=0;;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;

    // second iteration: remove baseline
    trace->Subtract(  1.47257e+01 );

    trace->ComputePedestal(1,20*bins_per_ns,meanPed,rmsPed);
    hMeanPed->Fill( meanPed );
    hRMSPed->Fill( rmsPed );

    trace->ComputeMin(30*bins_per_ns,60*bins_per_ns,min_mV,min_ns);
    //trace->ComputeMin(1,bins,min_mV,min_ns);
    hMinmV->Fill( min_mV );
    hMinns->Fill( trace->GetXaxis()->GetBinCenter(min_ns) );
  }
  new TCanvas();
  trace->DrawCopy();
  new TCanvas();
  summary->DrawCopy("colz");

  new TCanvas();
  hMeanPed->Draw();
  hMeanPed->Fit("gaus","I");
  new TCanvas();
  hRMSPed->Draw();

  new TCanvas();
  hMinmV->Draw();
  new TCanvas();
  hMinns->Draw();
  
  return 0;
}
