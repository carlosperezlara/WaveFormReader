#include "LECROYTRCReader.cxx"
#include "LECROYCSVReader.cxx"
#include "DRS4DATReader.cxx"
#include "WaveForm.cxx"

int comparer() {
  LECROYTRCReader file1("C2--41V--00000.trc");
  LECROYTRCReader file2("C3--41V--00000.trc");
  file1.ReadHeader();
  file2.ReadHeader();

  WaveForm *trace1 = file1.GetTrace();
  WaveForm *trace2 = file2.GetTrace();
  TH2D *summary1 = file1.GetSummaryPlot();
  TH2D *summary2 = file2.GetSummaryPlot();
  TH2D *histcorr = new TH2D("corr","Minimum;C2  (mV);C3  (mV)",100,-10,-600,100,-10,-600);
  double meanC1, meanC2, rmsC1, rmsC2;
  for(int nev=0;;++nev) {
    if(!file1.ReadEvent()) break;
    if(!file2.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;
    // i have access through trace
    trace1->ComputePedestal(1,200,meanC1,rmsC1);
    trace2->ComputePedestal(1,200,meanC2,rmsC2);
    histcorr->Fill( trace1->GetMinimum()-meanC1, trace2->GetMinimum()-meanC2 );
  }
  //trace->DrawCopy();
  TCanvas *csummary = new TCanvas("summary");
  csummary->Divide(2,1);
  csummary->cd(1); summary1->DrawCopy("colz");
  csummary->cd(2); summary2->DrawCopy("colz");
  TCanvas *ccorr = new TCanvas("correlation");
  histcorr->Draw("colz");
  TLatex *tex = new TLatex();
  tex->DrawLatexNDC(0.12,0.82,Form("C2/C3 = %.3f",histcorr->GetMean(2)/histcorr->GetMean(1)));
  return 0;
}
