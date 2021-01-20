#include "LECROYTRCReader.cxx"
#include "WaveForm.cxx"

int readfile_ComputeGain(TString filename="C2--4Layer-43--00000.trc",
			 TString pulsefile="C2--4Layer-45--00000.trc_step05.dat") {
  cout << "InputFile: " << filename.Data() << endl;
  LECROYTRCReader file(filename.Data());
  file.ReadHeader();
  WaveForm *trace = file.GetTrace();
  TH2D *summary = file.GetSummaryPlot();

  cout << " RANGES from " << "ranges.dat" << endl;
  ifstream fin;
  fin.open( "ranges.dat" );
  double ranges_baseline=20;
  double ranges_gate_low=40;
  double ranges_gate_high=60;
  double ranges_extreme_low=-14.0;
  double ranges_extreme_high=-8.0;
  double ranges_amplitude_low=2.0;
  double ranges_amplitude_high=4.0;
  double ranges_fit_maxchi2=2.0;
  fin >> ranges_baseline >> ranges_gate_low >> ranges_gate_high >> ranges_extreme_low >> ranges_extreme_high;
  fin >> ranges_amplitude_low >> ranges_amplitude_high >> ranges_fit_maxchi2;
  fin.close();
  cout << "   => RANGES:BASELINE " << ranges_baseline << endl;
  cout << "   => RANGES:GATE_LOW " << ranges_gate_low << endl;
  cout << "   => RANGES:GATE_HIGH " << ranges_gate_high << endl;
  cout << "   => RANGES:EXTREME_LOW " << ranges_extreme_low << endl;
  cout << "   => RANGES:EXTREME_HIGH " << ranges_extreme_high << endl;
  cout << "   => RANGES:AMPLITUDE_LOW " << ranges_amplitude_low << endl;
  cout << "   => RANGES:AMPLITUDE_HIGH " << ranges_amplitude_high << endl;
  cout << "   => RANGES:FIT_MAXCHI2 " << ranges_fit_maxchi2 << endl;

  double baseline=0;
  cout << " BASELINE from " << filename.Data() << Form("_step%02d.dat",0) << endl;
  fin.open( Form("%s_step00.dat",filename.Data()) );
  fin >> baseline;
  fin.close();
  cout << " PULSE from " << pulsefile.Data() << endl;
  trace->LoadTemplate( Form("%s",pulsefile.Data()) );

  Int_t bins = trace->GetXaxis()->GetNbins();
  Double_t bins_per_ns = bins / ( trace->GetXaxis()->GetBinLowEdge(bins+1) - trace->GetXaxis()->GetBinLowEdge(1) );
  cout << "Bins per ns " << bins_per_ns << endl;
  Double_t max_mV = summary->GetYaxis()->GetBinLowEdge( summary->GetYaxis()->GetNbins()+1 );
  Double_t min_mV = summary->GetYaxis()->GetBinLowEdge( 1 );
  Double_t max_ns = summary->GetXaxis()->GetBinLowEdge( summary->GetXaxis()->GetNbins()+1 );
  Double_t min_ns = summary->GetXaxis()->GetBinLowEdge( 1 );
  
  Double_t meanPed, rmsPed;
  TH1D *hMeanPed = new TH1D("hMeanPed","MeanPed;mV",100,min_mV-baseline,max_mV-baseline);
  TH1D *hRMSPed = new TH1D("hRMSPed","RMSPed;mV",100,0,+15);

  Double_t extreme_mV;
  Int_t extreme_timebin;
  //TH1D *hMinmV = new TH1D("hMinmV","MinmV;mV",1000,min_mV-baseline,max_mV-baseline);
  TH1D *hMinmV = new TH1D("hMinmV","MinmV;mV",1000,min_mV-baseline,0);
  TH1D *hMinns = new TH1D("hMinns","MinnS;ns",1000,min_ns,max_ns);

  TH1D *hFitChi2_all      = new TH1D("hFitChi2_all", "Fit_Chi2_all",         1000,0,100);
  TH1D *hFitBaseline_all  = new TH1D("hFitBaseline_all", "Fit_Baseline;mV",  100,-5,+5);
  TH1D *hFitAmplitude_all = new TH1D("hFitAmplitude_all","Fit_Amplitude;mV", 1000,-5,+20);
  TH1D *hFitWalk_all      = new TH1D("hFitWalk_all",     "Fit_Walk;ns",      1000,-10,+10);

  TH1D *hFitChi2_sel      = new TH1D("hFitChi2_sel", "Fit_Chi2_sel",         1000,0,100);
  TH1D *hFitBaseline_sel  = new TH1D("hFitBaseline_sel", "Fit_Baseline;mV",  100,-5,+5);
  TH1D *hFitAmplitude_sel = new TH1D("hFitAmplitude_sel","Fit_Amplitude;mV", 1000,-5,+20);
  TH1D *hFitWalk_sel      = new TH1D("hFitWalk_sel",     "Fit_Walk;ns",      1000,-10,+10);

  TH2D *hFitBaseline2D_all  = new TH2D("hFitBaseline2D_all", "Fit_Baseline;fit mV;trace mV",       100,-5,+5,   100,-5,+5);
  TH2D *hFitAmplitude2D_all = new TH2D("hFitAmplitude2D_all","Fit_Amplitude;fit mV;trace_minmV mV",100,-5,+20,  100,min_mV-baseline,0);
  TH2D *hFitAmplChi2_all    = new TH2D("hFitAmplChi2_all",   "Fit_AmplChi2;mV",                    1000,-5,+20, 100,0,100);
  TH2D *hFitWalk2D_all      = new TH2D("hFitWalk2D_all",     "Fit_Walk;fit_walk ns;trace_tmax ns", 100,-10,+10, 100,min_ns,max_ns);

  TH2D *hFitBaseline2D_sel  = new TH2D("hFitBaseline2D_sel", "Fit_Baseline;fit mV;trace mV",       100,-5,+5,   100,-5,+5);
  TH2D *hFitAmplitude2D_sel = new TH2D("hFitAmplitude2D_sel","Fit_Amplitude;fit mV;trace_minmV mV",100,-5,+20,  100,min_mV-baseline,0);
  TH2D *hFitAmplChi2_sel    = new TH2D("hFitAmplChi2_sel",   "Fit_AmplChi2;mV",                    1000,-5,+20, 100,0,ranges_fit_maxchi2);
  TH2D *hFitWalk2D_sel      = new TH2D("hFitWalk2D_sel",     "Fit_Walk;fit_walk ns;trace_tmax ns", 100,-10,+10, 100,min_ns,max_ns);

  TH1D *hExample[100];
  for(int i=0; i!=100; ++i) hExample[i] = NULL;

  TList *list = new TList();
  list->SetName("AcceptedFitExamples");
  list->SetOwner();
  
  Int_t nSinglesTraces = 0;
  int nExample = 0;
  for(int nev=0;;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;

    trace->ComputePedestal(1,ranges_baseline*bins_per_ns,meanPed,rmsPed);
    hMeanPed->Fill( meanPed );
    hRMSPed->Fill( rmsPed );

    trace->ComputeMin(ranges_gate_low*bins_per_ns,ranges_gate_high*bins_per_ns,extreme_mV,extreme_timebin);
    hMinmV->Fill( extreme_mV );
    Double_t extreme_ns = trace->GetXaxis()->GetBinCenter(extreme_timebin);
    hMinns->Fill( extreme_ns );

    Double_t ampl = abs(extreme_mV - baseline);
    Double_t gate = 2*ampl;
    if(ampl<3) gate=5;
    //cout << "  AMPL [" << ampl-gate << " " << ampl << " " << ampl+gate << "]" << endl;
    
    //int res = trace->FitTemplate(meanPed,ampl,   0,
    int res = trace->FitTemplate(meanPed,ampl,  -5, //trick: baseline set to -5 (instead of 0) to force randomization
				 3,      gate,  10,
				 6, "WWRMQ");
    //cout << "*********** " << res << endl;
    Double_t chi2 = trace->GetLastReducedChiSquared();
    hFitChi2_all->Fill( chi2 );
    hFitBaseline_all   ->Fill( trace->EstimateBaseline() );
    hFitAmplitude_all  ->Fill( trace->EstimateAmplitude() );
    hFitWalk_all       ->Fill( trace->EstimateWalk() );
    hFitBaseline2D_all ->Fill( trace->EstimateBaseline(), meanPed );
    hFitAmplitude2D_all->Fill( trace->EstimateAmplitude(), extreme_mV );
    hFitAmplChi2_all   ->Fill( trace->EstimateAmplitude(), chi2 );
    hFitWalk2D_all     ->Fill( trace->EstimateWalk(), extreme_ns );

    ampl = trace->EstimateAmplitude();
    if(chi2<ranges_fit_maxchi2 && res==0) {
      hFitChi2_sel->Fill( chi2 );
      hFitBaseline_sel   ->Fill( trace->EstimateBaseline() );
      hFitAmplitude_sel  ->Fill( trace->EstimateAmplitude() );
      hFitWalk_sel       ->Fill( trace->EstimateWalk() );
      hFitBaseline2D_sel ->Fill( trace->EstimateBaseline(), meanPed );
      hFitAmplitude2D_sel->Fill( trace->EstimateAmplitude(), extreme_mV );
      hFitAmplitude2D_sel->Fill( trace->EstimateAmplitude(), chi2 );
      hFitWalk2D_sel     ->Fill( trace->EstimateWalk(), extreme_ns );
      nSinglesTraces++;
      if(nExample<100) {
	hExample[nExample] = (TH1D*) trace->Clone( Form("SelectedFitExample_%d_%d",nExample,nev) );
	nExample++;
      }
    }
  }
  cout << " Number of good singles used in profile " << nSinglesTraces << endl;
  new TCanvas();
  trace->DrawCopy();

  hMeanPed->Fit("gaus","I");

  TFile *foutroot = new TFile( Form("%s_PEAKS.root",filename.Data()), "RECREATE" );
  summary->Write();
  hMeanPed->Write();
  hRMSPed->Write();
  hMinmV->Write();
  hMinns->Write();

  hFitChi2_all->Write();
  hFitBaseline_all->Write();
  hFitAmplitude_all->Write();
  hFitWalk_all->Write();

  hFitBaseline2D_all->Write();
  hFitAmplitude2D_all->Write();
  hFitAmplChi2_all->Write();
  hFitWalk2D_all->Write();

  hFitChi2_sel->Write();
  hFitBaseline_sel->Write();
  hFitAmplitude_sel->Write();
  hFitWalk_sel->Write();

  hFitBaseline2D_sel->Write();
  hFitAmplitude2D_sel->Write();
  hFitAmplChi2_sel->Write();
  hFitWalk2D_sel->Write();

  for(int i=0; i!=nExample; ++i) {
    if(hExample[i])
      list->Add( hExample[i] );
  }
  cout << "  Number of examples: " << nExample << endl;
  list->Write("FitExamples",TObject::kSingleKey);

  foutroot->Close();

  return 0;
}
