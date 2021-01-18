Double_t first_min[7]  = { 1, 1, 1, 0.5, 1.0, 2.0, 3.0 };
Double_t first_max[7]  = { 2, 2, 2, 1.1, 1.8, 3.0, 4.0 };
Double_t second_min[7] = { 2, 2, 2, 1.4, 2.2, 4.0, 6.2 };
Double_t second_max[7] = { 3, 3, 3, 1.8, 3.0, 5.6, 7.6 };
Double_t third_min[7]  = { 3, 3, 3, 2.2, 3.5, 6.0, 9.5 };
Double_t third_max[7]  = { 4, 4, 4, 2.8, 4.5, 8.0,11.5 };

double xmin[7] = {1,1,1, 0.3, 0.5, 1.5, 2.0};
double xmax[7] = {4,4,4, 3.3, 5.0, 8.2, 11.5};

Double_t TMK(Double_t *x, Double_t *p, Int_t fRejectFit) {

  bool selectRegion = false;
  if( x[0] > first_min[fRejectFit]  && x[0] < first_max[fRejectFit]  ) selectRegion = true;
  if( x[0] > second_min[fRejectFit] && x[0] < second_max[fRejectFit] ) selectRegion = true;
  if( x[0] > third_min[fRejectFit]  && x[0] < third_max[fRejectFit]  ) selectRegion = true;

  if(!selectRegion) {
    TF1::RejectPoint();
    return 0;
  }

  Double_t first  = p[0] * TMath::Gaus( x[0] - p[3], p[6], 1 );
  Double_t second = p[1] * TMath::Gaus( x[0] - p[4], p[7], 1 );
  Double_t third  = p[2] * TMath::Gaus( x[0] - p[5], p[8], 1 );
  return first + second + third;
}

Double_t TMK0(Double_t *x, Double_t *p) {return TMK(x,p,0);}
Double_t TMK1(Double_t *x, Double_t *p) {return TMK(x,p,1);}
Double_t TMK2(Double_t *x, Double_t *p) {return TMK(x,p,2);}
Double_t TMK3(Double_t *x, Double_t *p) {return TMK(x,p,3);}
Double_t TMK4(Double_t *x, Double_t *p) {return TMK(x,p,4);}
Double_t TMK5(Double_t *x, Double_t *p) {return TMK(x,p,5);}
Double_t TMK6(Double_t *x, Double_t *p) {return TMK(x,p,6);}

int drawPeakTrend() {
  gStyle->SetOptStat(0);
  const int files=7;
  TString filename[files];
  filename[0] = "C2--4Layer-38.5--00000.trc";
  filename[1] = "C2--4Layer-39.0--00000.trc";
  filename[2] = "C2--4Layer-39.5--00000.trc";
  filename[3] = "C2--4Layer-40--00000.trc";
  filename[4] = "C2--4Layer-41--00000.trc";
  filename[5] = "C2--4Layer-43--00000.trc";
  filename[6] = "C2--4Layer-45--00000.trc";

  double xmin[7] = {1,1,1, 0.3, 0.5, 1.5, 2.0};
  double xmax[7] = {4,4,4, 3.3, 5.0, 8.2, 11.5};
  
  int color[7] = {
    kYellow-3,kCyan-3,kGreen-3,
    kMagenta-3,kRed-3,kBlue-3,
    kBlack };

  double f_est[files] = { 0.0, 0.0, 0.0, 1.00, 1.50, 2.50, 3.50 };
  double f_gat[files] = { 0.5, 0.5, 0.5, 0.30, 0.37, 0.44, 0.51 };
  
  TH1D *gr[files];
  TF1  *ft[files];
  ft[0] = new TF1( Form("fit1_%d",0), TMK0, 0, 13, 9 );
  ft[1] = new TF1( Form("fit1_%d",1), TMK1, 0, 13, 9 );
  ft[2] = new TF1( Form("fit1_%d",2), TMK2, 0, 13, 9 );
  ft[3] = new TF1( Form("fit1_%d",3), TMK3, 0, 13, 9 );
  ft[4] = new TF1( Form("fit1_%d",4), TMK4, 0, 13, 9 );
  ft[5] = new TF1( Form("fit1_%d",5), TMK5, 0, 13, 9 );
  ft[6] = new TF1( Form("fit1_%d",6), TMK6, 0, 13, 9 );

  TF1  *ft1[files];
  TF1  *ft2[files];
  TF1  *ft3[files];
  
  for(int i=4; i!=files; ++i) {
    TFile *file = new TFile( Form("%s_PEAKS.root",filename[i].Data()) );
    gr[i] = (TH1D*) ( (TH1D*) file->Get("hFitAmplitude_all") )->Clone( Form("D%d",i)  );

    cout << endl << endl << endl;
    
    gr[i]->Rebin(2);

    if(0) {
      ft[i]->SetLineColor( color[i%7] );

      ft[i]->SetParameter(0,100);  ft[i]->SetParLimits(0,1,5e3);  ft[i]->SetParName(0,"1PE_yield");
      ft[i]->SetParameter(1,100);  ft[i]->SetParLimits(1,1,5e3);  ft[i]->SetParName(1,"2PE_yield");
      ft[i]->SetParameter(2,100);  ft[i]->SetParLimits(2,1,5e3);  ft[i]->SetParName(2,"3PE_yield");

      ft[i]->SetParameter(3,0.5*(first_min[i]+first_max[i]));   ft[i]->SetParLimits(3,first_min[i],first_max[i]);   ft[i]->SetParName(3,"1PE_mean");
      ft[i]->SetParameter(4,0.5*(second_min[i]+second_max[i])); ft[i]->SetParLimits(4,second_min[i],second_max[i]); ft[i]->SetParName(4,"2PE_mean");
      ft[i]->SetParameter(5,0.5*(third_min[i]+third_max[i]));   ft[i]->SetParLimits(5,third_min[i],third_max[i]);   ft[i]->SetParName(5,"3PE_mean");

      ft[i]->SetParameter(6,0.01); ft[i]->SetParLimits(6,0.00001,0.1);  ft[i]->SetParName(6,"1PE_sigma");
      ft[i]->SetParameter(7,0.01); ft[i]->SetParLimits(7,0.00001,0.2);  ft[i]->SetParName(7,"2PE_sigma");
      ft[i]->SetParameter(8,0.01); ft[i]->SetParLimits(8,0.00001,0.3);  ft[i]->SetParName(8,"3PE_sigma");
      gr[i]->Fit( ft[i], "IEMR", "", xmin[i], xmax[i] );
      gr[i]->Fit( ft[i], "IEMR", "", xmin[i], xmax[i] );
    } else {
      ft1[i] = new TF1( Form("fit1_%d",i), "gaus", first_min[i],  first_max[i] );
      ft2[i] = new TF1( Form("fit2_%d",i), "gaus", second_min[i], second_max[i] );
      ft3[i] = new TF1( Form("fit3_%d",i), "gaus", third_min[i],  third_max[i] );

      ft1[i]->SetLineColor( color[i%7] );
      ft2[i]->SetLineColor( color[i%7] );
      ft3[i]->SetLineColor( color[i%7] );

      ft1[i]->SetParameter(0,100);  ft1[i]->SetParLimits(0,1,5e3);  ft1[i]->SetParName(0,"1PE_yield");
      ft2[i]->SetParameter(1,100);  ft2[i]->SetParLimits(1,1,5e3);  ft2[i]->SetParName(1,"2PE_yield");
      ft3[i]->SetParameter(2,100);  ft3[i]->SetParLimits(2,1,5e3);  ft3[i]->SetParName(2,"3PE_yield");

      ft1[i]->SetParameter(3,0.5*(first_min[i]+first_max[i]));   ft1[i]->SetParLimits(3,first_min[i],first_max[i]);   ft1[i]->SetParName(3,"1PE_mean");
      ft2[i]->SetParameter(4,0.5*(second_min[i]+second_max[i])); ft2[i]->SetParLimits(4,second_min[i],second_max[i]); ft2[i]->SetParName(4,"2PE_mean");
      ft3[i]->SetParameter(5,0.5*(third_min[i]+third_max[i]));   ft3[i]->SetParLimits(5,third_min[i],third_max[i]);   ft3[i]->SetParName(5,"3PE_mean");

      ft1[i]->SetParameter(6,0.01); ft1[i]->SetParLimits(6,0.00001,0.1);  ft1[i]->SetParName(6,"1PE_sigma");
      ft2[i]->SetParameter(7,0.01); ft2[i]->SetParLimits(7,0.00001,0.2);  ft2[i]->SetParName(7,"2PE_sigma");
      ft3[i]->SetParameter(8,0.01); ft3[i]->SetParLimits(8,0.00001,0.3);  ft3[i]->SetParName(8,"3PE_sigma");

      gr[i]->Fit( ft1[i], "IEMR", "", first_min[i],  first_max[i] );
      gr[i]->Fit( ft2[i], "IEMR", "", second_min[i], second_max[i] );
      gr[i]->Fit( ft3[i], "IEMR", "", third_min[i],  third_max[i] );
      //ft1[i]->Draw("same");
      //ft2[i]->Draw("same");
      //ft3[i]->Draw("same");
    }
    
    //gr[i]->Sumw2();
    gr[i]->SetMarkerStyle( 1 );
    //gr[i]->SetMarkerStyle(24);
    gr[i]->SetMarkerColor( color[i%7] );
    gr[i]->SetLineColor( color[i%7] );
    //gr[i]->SetLineWidth(2);
    //gr[i]->Fit( ft[i], "IPEM" );
    //gr[i]->Fit( ft[i], "IPEM" );
    //gr[i]->Fit( ft[i], "IPEM" );
    //gr[i]->Fit( ft[i], "IPEM" );
    //gr[i]->Fit( ft[i], "IPEM" );
  }

  TLegend *leg = new TLegend(0.45,0.6,0.9,0.9);
  gr[6]->Draw("E");
  for(int i=4; i!=files; ++i) {
    gr[i]->Draw("Esame");
    ft1[i]->Draw("same");
    ft2[i]->Draw("same");
    ft3[i]->Draw("same");
    leg->AddEntry( gr[i], Form("%s",filename[i].Data()) );
  }
  leg->Draw();

  TCanvas *main = new TCanvas();
  TF1 *mVperPE[files];
  Double_t x[3];
  Double_t y[files][3];
  TGraph *grL[files];
  for(int i=0; i!=3; ++i) {
    x[i] = i+1;
  }
  TLegend *leg2 = new TLegend(0.1,0.5,0.5,0.9);
  TH2D *axis = new TH2D("axis",";PE;mV  per  PE",100,0.5,3.5,100,0,15);
  for(int file=4; file!=files; ++file) {
    y[file][0] = ft1[file]->GetParameter(1);
    y[file][1] = ft2[file]->GetParameter(1);
    y[file][2] = ft3[file]->GetParameter(1);
    grL[file] = new TGraph(3,x,y[file]);
    grL[file]->SetMarkerStyle(20);
    grL[file]->SetMarkerColor( color[file%7] );
    mVperPE[file] = new TF1(Form("mVperPE_%d",file),"[0]+[1]*x",0,4);
    mVperPE[file]->SetLineColor( color[file%7] );
    grL[file]->Fit(mVperPE[file]);
    leg2->AddEntry( grL[file], Form("%s",filename[file].Data()) );
  }
  axis->Draw();
  grL[4]->Draw("PSAME");
  grL[5]->Draw("PSAME");
  grL[6]->Draw("PSAME");
  leg2->Draw();

  TLatex *tex = new TLatex();
  TH2D *axis2 = new TH2D("axis2",";Bias Voltage;mV per PE",100,38,46,100,0,4);
  TF1 *fitfit = new TF1("fitfit","[1]*(x-[0])");
  Double_t xxx[3] = { 41.0, 43.0, 45.0 };
  Double_t yyy[3];
  yyy[0] = mVperPE[4]->GetParameter( 1 );
  yyy[1] = mVperPE[5]->GetParameter( 1 );
  yyy[2] = mVperPE[6]->GetParameter( 1 );
  TCanvas *main2 = new TCanvas();
  TGraph *grF = new TGraph(3,xxx,yyy);
  axis2->Draw();
  grF->Draw("*SAME");
  grF->Fit(fitfit);
  tex->DrawLatexNDC( 0.15, 0.8, Form("%.2f (x - %.2f)",fitfit->GetParameter(1),fitfit->GetParameter(0)) );
  
  return 0;
}
