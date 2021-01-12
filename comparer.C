#include "LECROYTRCReader.cxx"
#include "LECROYCSVReader.cxx"
#include "DRS4DATReader.cxx"

int comparer() {
  //LECROYTRCReader file("C2--2e14-40.15V-n30-Green--00000.trc");
  //LECROYCSVReader file("C3--scope_45VHighGain.txt");
  //DRS4DATReader file("same.dat");

  LECROYTRCReader file1("12-18-2020/C2--2e14-40.65V-n30-Green--00000.trc");
  LECROYTRCReader file2("12-18-2020/C3--2e14-40.65V-n30-Green--00000.trc");
  file1.ReadHeader();
  file2.ReadHeader();

  TH1D *trace1 = file1.GetTrace();
  TH1D *trace2 = file2.GetTrace();
  TH2D *summary1 = file1.GetSummaryPlot();
  TH2D *summary2 = file2.GetSummaryPlot();
  TH2D *histcorr = new TH2D("corr","Minimum;C2  (mV);C3  (mV)",100,-10,-600,100,-10,-600);
  for(int nev=0;;++nev) {
    if(!file1.ReadEvent()) break;
    if(!file2.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    // i have access through trace
    histcorr->Fill( trace1->GetMinimum(), trace2->GetMinimum() );
  }
  //trace->DrawCopy();
  TCanvas *csummary = new TCanvas("summary");
  csummary->Divide(2,1);
  csummary->cd(1); summary1->DrawCopy("colz");
  csummary->cd(2); summary2->DrawCopy("colz");
  TCanvas *ccorr = new TCanvas("correlation");
  histcorr->Draw("colz");
  TLatex *tex = new TLatex();
  tex->DrawLatexNDC(0.12,0.82,Form("C3/C2 = %.1f",histcorr->GetMean(2)/histcorr->GetMean(1)));
  return 0;
}
