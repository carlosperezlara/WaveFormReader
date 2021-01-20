#include "LECROYTRCReader.cxx"

int readfile_OverlayEvents(TString filename="C2--Jan21Board2-45V-0.45mA-00000--00000.trc",Int_t ev0=0, Int_t ev1=140, Int_t ev2=127) {
  gStyle->SetOptStat(0);
  LECROYTRCReader file("C2--Jan21Board2-45V-0.45mA-00000--00000.trc");
  file.ReadHeader();

  TH1D *trace = file.GetTrace();
  TH1D *evShow[3] = {NULL,NULL,NULL};
  int color[3] = { kBlue-3, kRed-3, kGreen-3 };
  TH2D *summary = file.GetSummaryPlot();
  for(int nev=0;;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    // i have access through trace
    if(nev==ev0) evShow[0] = (TH1D*) trace->Clone( Form("Event_%05d",nev) );
    if(nev==ev1) evShow[1] = (TH1D*) trace->Clone( Form("Event_%05d",nev) );
    if(nev==ev2) evShow[2] = (TH1D*) trace->Clone( Form("Event_%05d",nev) );
    if(evShow[0]&&evShow[1]&&evShow[2]) break;
  }
  new TCanvas();
  summary->DrawCopy("colz");
  summary->Reset();
  
  new TCanvas();
  summary->DrawCopy();
  TLegend *leg = new TLegend(0.1,0.7,0.35,0.9);
  for(int i=0; i!=3; ++i) {
    if(evShow[i]) {
      evShow[i]->SetLineColor( color[i] );
      evShow[i]->Draw("hsame");
      leg->AddEntry( evShow[i] , evShow[i]->GetName() );
    }
  }
  leg->Draw();
  return 0;
}
