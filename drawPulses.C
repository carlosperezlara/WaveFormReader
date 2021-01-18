int drawPulses( TString filename="C2--4Layer-45--00000.trc" ) {

  const int pulses=12;
  int color[7] = {
    kYellow-3,kCyan-3,kGreen-3,
    kMagenta-3,kRed-3,kBlue-3,
    kBlack };

  TGraph *gr[pulses];
  for(int i=0; i!=pulses; ++i) {
    gr[i] = new TGraph( Form("%s_Step%02d.dat",filename.Data(),i+1) );
    /*
    int n = gr[i]->GetN();
    double ymin = 999;
    double xmin = 0;
    for(int b=0; b!=n; ++b) {
      double val = gr[i]->GetPointY(b);
      if(val<ymin) {
	xmin = gr[i]->GetPointX(b);
	ymin = val;
      }
    }
    Double_t *equis = gr[i]->GetX();
    for(int b=0; b!=n; ++b) {
      equis[b] += 60-xmin;
    }
    */
    gr[i]->SetMarkerStyle(24);
    gr[i]->SetMarkerColor( color[i%7] );
    gr[i]->SetLineColor( color[i%7] );
    gr[i]->SetLineWidth(2);
  }

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  gr[0]->Draw("Al");
  for(int i=0; i<pulses; ++i) {
    gr[i]->Draw("lsame");
    leg->AddEntry( gr[i], Form("Step %d",i+1) );
  }
  leg->Draw();
  
  return 0;
}
