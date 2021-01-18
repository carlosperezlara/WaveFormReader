int drawPeaks( TString filename="C2--4Layer-45--00000.trc" ) {

  const int pulses=12;
  int color[7] = {
    kYellow-3,kCyan-3,kGreen-3,
    kMagenta-3,kRed-3,kBlue-3,
    kBlack };

  TH1D *gr[pulses];
  for(int i=2; i!=pulses; ++i) {
    TFile *file = new TFile( Form("%s_Step%02d.root",filename.Data(),i) );
    gr[i] = (TH1D*) ( (TH1D*) file->Get("hFitAmplitude_all") )->Clone( Form("D%d",i)  );
    gr[i]->Rebin(4);
    gr[i]->SetMarkerStyle(24);
    gr[i]->SetMarkerColor( color[i%7] );
    gr[i]->SetLineColor( color[i%7] );
    gr[i]->SetLineWidth(2);
  }

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  gr[2]->Draw("h");
  for(int i=2; i<pulses; ++i) {
    gr[i]->Draw("hsame");
    leg->AddEntry( gr[i], Form("Step %d",i) );
  }
  leg->Draw();
  
  return 0;
}
