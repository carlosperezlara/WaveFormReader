#include "LECROYTRCReader.cxx"
#include "WaveForm.cxx"

int readfile_EstimatePulseShape(TString filename="C2--4Layer-45--00000.trc",Int_t nstep=0) {
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

  if(nstep>0) {    // STEP1++: FIT PULSE
    cout << " PULSE from " << filename.Data() << Form("_step%02d.dat",nstep-1) << endl;
    trace->LoadTemplate( Form("%s_step%02d.dat",filename.Data(), nstep-1) );
  }

  Int_t bins = trace->GetXaxis()->GetNbins();
  Double_t bins_per_ns = bins / ( trace->GetXaxis()->GetBinLowEdge(bins+1) - trace->GetXaxis()->GetBinLowEdge(1) );
  cout << "Bins per ns " << bins_per_ns << endl;
  Double_t max_mV = summary->GetYaxis()->GetBinLowEdge( summary->GetYaxis()->GetNbins()+1 );
  Double_t min_mV = summary->GetYaxis()->GetBinLowEdge( 1 );
  Double_t max_ns = summary->GetXaxis()->GetBinLowEdge( summary->GetXaxis()->GetNbins()+1 );
  Double_t min_ns = summary->GetXaxis()->GetBinLowEdge( 1 );
  
  Double_t meanPed, rmsPed;
  Double_t extreme_mV;
  Int_t extreme_timebin;

  TH1D *hMeanPed   = new TH1D("hMeanPed_RAW", "MeanPed_RAW;mV",100,min_mV,max_mV);
  TH1D *hRMSPed    = new TH1D("hRMSPed_RAW",  "RMSPed_RAW;mV",100,0,+8);
  
  TList *lstNoCuts = new TList();
  TH1D *hMeanCent0 = new TH1D("hMeanCent_SUB","MeanCent_SUB;ns",   100,ranges_gate_low,ranges_gate_high);
  TH1D *hRMSCent0  = new TH1D("hRMSCent_SUB", "RMSCent_SUB;ns",    100,-18,+18);
  TH2D *hCent2D0   = new TH2D("hCent2D_SUB",  "MeanCent_SUB;ns;ns",100,ranges_gate_low,ranges_gate_high,100,-18,+18);
  TH1D *hMinmV0    = new TH1D("hMinmV_SUB",   "MinmV_SUB;mV",      1000,min_mV,0);
  TH1D *hMinns0    = new TH1D("hMinns_SUB",   "MinnS_RAW;ns",      1000,ranges_gate_low,ranges_gate_high);
  lstNoCuts->Add( hMinmV0 );
  lstNoCuts->Add( hMinns0 );
  lstNoCuts->Add( hMeanCent0 );
  lstNoCuts->Add( hRMSCent0 );
  lstNoCuts->Add( hCent2D0 );
  lstNoCuts->SetOwner();

  TList *lstCuts = new TList();
  TH1D *hMeanCent1 = new TH1D("hMeanCent_SUB_CUT","MeanCent_SUB;ns",   100,ranges_gate_low,ranges_gate_high);
  TH1D *hRMSCent1  = new TH1D("hRMSCent_SUB_CUT", "RMSCent_SUB;ns",    100,-18,+18);
  TH2D *hCent2D1   = new TH2D("hCent2D_SUB_CUT",  "MeanCent_SUB;ns;ns",100,ranges_gate_low,ranges_gate_high,100,-18,+18);
  TH1D *hMinmV1    = new TH1D("hMinmV_SUB_CUT",   "MinmV_SUB;mV",      1000,min_mV,0);
  TH1D *hMinns1    = new TH1D("hMinns_SUB_CUT",   "MinnS_RAW;ns",      1000,ranges_gate_low,ranges_gate_high);
  lstCuts->Add( hMinmV1 );
  lstCuts->Add( hMinns1 );
  lstCuts->Add( hMeanCent1 );
  lstCuts->Add( hRMSCent1 );
  lstCuts->Add( hCent2D1 );
  lstCuts->SetOwner();
  
  TList *lstFitNoCuts = new TList();
  TH1D *hFitChi2_all        = new TH1D("hFitChi2_all",       "Fit_Chi2_all",     1000,0,+30);
  TH1D *hFitBaseline_all    = new TH1D("hFitBaseline_all",   "Fit_Baseline;mV",  100,-5,+5);
  TH1D *hFitAmplitude_all   = new TH1D("hFitAmplitude_all",  "Fit_Amplitude;mV", 1000,-5,+20);
  TH1D *hFitWalk_all        = new TH1D("hFitWalk_all",       "Fit_Walk;ns",      1000,-10,+10);
  TH2D *hFitBaseline2D_all  = new TH2D("hFitBaseline2D_all", "Fit_Baseline;fit mV;trace_baseline mV", 100,-5,+5,   100,-500,+500);
  TH2D *hFitAmplitude2D_all = new TH2D("hFitAmplitude2D_all","Fit_Amplitude;fit mV;trace_ampl mV",    100,-5,+20,  100,min_mV,3);
  TH2D *hFitWalk2D_all      = new TH2D("hFitWalk2D_all",     "Fit_Walk;fit_walk ns;trace_tmax ns",    100,-10,+10, 100,min_ns,max_ns);
  TH2D *hFitWalk2DC_all     = new TH2D("hFitWalk2DC_all",    "Fit_Walk;fit_walk ns;Fit_Amplitude mV", 100,-10,+10, 100,-5,+20);
  lstFitNoCuts->Add( hFitChi2_all );
  lstFitNoCuts->Add( hFitBaseline_all );
  lstFitNoCuts->Add( hFitAmplitude_all );
  lstFitNoCuts->Add( hFitWalk_all );
  lstFitNoCuts->Add( hFitBaseline2D_all );
  lstFitNoCuts->Add( hFitAmplitude2D_all );
  lstFitNoCuts->Add( hFitWalk2D_all );
  lstFitNoCuts->Add( hFitWalk2DC_all );
  lstFitNoCuts->SetOwner();
  
  TList *lstFitCuts = new TList();
  TH1D *hFitChi2_sel        = new TH1D("hFitChi2_sel",       "Fit_Chi2_sel",     1000,0,+30);
  TH1D *hFitBaseline_sel    = new TH1D("hFitBaseline_sel",   "Fit_Baseline;mV",  100,-5,+5);
  TH1D *hFitAmplitude_sel   = new TH1D("hFitAmplitude_sel",  "Fit_Amplitude;mV", 1000,-5,+20);
  TH1D *hFitWalk_sel        = new TH1D("hFitWalk_sel",       "Fit_Walk;ns",      1000,-10,+10);
  TH2D *hFitBaseline2D_sel  = new TH2D("hFitBaseline2D_sel", "Fit_Baseline;fit mV;trace_baseline mV", 100,-5,+5,   100,-500,+500);
  TH2D *hFitAmplitude2D_sel = new TH2D("hFitAmplitude2D_sel","Fit_Amplitude;fit mV;trace_ampl mV",    100,-5,+20,  100,min_mV,0);
  TH2D *hFitWalk2D_sel      = new TH2D("hFitWalk2D_sel",     "Fit_Walk;fit_walk ns;trace_tmax ns",    100,-10,+10, 100,min_ns,max_ns);
  TH2D *hFitWalk2DC_sel     = new TH2D("hFitWalk2DC_sel",    "Fit_Walk;fit_walk ns;Fit_Amplitude mV", 100,-10,+10, 100,-5,+20);
  lstFitCuts->Add( hFitChi2_sel );
  lstFitCuts->Add( hFitBaseline_sel );
  lstFitCuts->Add( hFitAmplitude_sel );
  lstFitCuts->Add( hFitWalk_sel );
  lstFitCuts->Add( hFitBaseline2D_sel );
  lstFitCuts->Add( hFitAmplitude2D_sel );
  lstFitCuts->Add( hFitWalk2D_sel );
  lstFitCuts->Add( hFitWalk2DC_sel );
  lstFitCuts->SetOwner();

  TH1D *hExample[100];
  for(int i=0; i!=100; ++i) hExample[i] = NULL;
  TList *list = new TList();
  list->SetOwner();
  int nExample = 0;
  
  trace->CreateProfile();
  Int_t nSinglesTraces = 0;
  for(int nev=0;;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;

    trace->ComputePedestal(1,ranges_baseline*bins_per_ns,meanPed,rmsPed);
    hMeanPed->Fill( meanPed );
    hRMSPed->Fill( rmsPed );
    trace->Subtract( meanPed );
    trace->ComputeMin(ranges_gate_low*bins_per_ns,ranges_gate_high*bins_per_ns,extreme_mV,extreme_timebin);
    Double_t extreme_ns = trace->GetXaxis()->GetBinCenter(extreme_timebin);
    Double_t centroid, rmscentroid;
    trace->ComputeCentroid(ranges_gate_low*bins_per_ns,ranges_gate_high*bins_per_ns,centroid,rmscentroid);
    hMinmV0->Fill( extreme_mV );
    hMinns0->Fill( extreme_ns );
    hMeanCent0->Fill(centroid);
    hRMSCent0->Fill(rmscentroid);
    hCent2D0->Fill(centroid,rmscentroid);

    bool passed = true;
    if(extreme_mV<ranges_extreme_low) passed = false;
    if(rmscentroid<0.5) passed = false;
    if(rmscentroid>2.5) passed = false;
    if(extreme_mV>ranges_extreme_high) passed = false;
    if(!passed) continue;

    hMinmV1->Fill( extreme_mV );
    hMinns1->Fill( extreme_ns );
    hMeanCent1->Fill(centroid);
    hRMSCent1->Fill(rmscentroid);
    hCent2D1->Fill(centroid,rmscentroid);
    
    if(nstep>0) {    // STEP1++: FIT PULSE
      Double_t ampl = 0.5*(ranges_amplitude_high + ranges_amplitude_low);
      Double_t gate = 5*0.5*(ranges_amplitude_high - ranges_amplitude_low); // fit up to five peaks
      int res = trace->FitTemplate(0,      ampl,   0,
				   3,      gate,   5,
				   20.0,  180.0, "WWRMEQ");
      //cout << "*********** " << res << endl;
      Double_t chi2 = trace->GetLastReducedChiSquared();
      hFitChi2_all->Fill( chi2 );
      hFitBaseline_all  ->Fill( trace->EstimateBaseline() );
      hFitAmplitude_all ->Fill( trace->EstimateAmplitude() );
      hFitWalk_all      ->Fill( trace->EstimateWalk() );
      hFitBaseline2D_all ->Fill( trace->EstimateBaseline(), meanPed );
      hFitAmplitude2D_all->Fill( trace->EstimateAmplitude(), extreme_mV );
      hFitWalk2D_all     ->Fill( trace->EstimateWalk(), extreme_ns );
      hFitWalk2DC_all    ->Fill( trace->EstimateWalk(), trace->EstimateAmplitude() );
	
      ampl = trace->EstimateAmplitude();
      bool pass = true;
      if(ampl<ranges_amplitude_low) pass = false;
      if(ampl>ranges_amplitude_high) pass = false;
      if(fabs(trace->EstimateBaseline()) > 0.5) pass = false;
      if(chi2>ranges_fit_maxchi2) pass = false;
      if(fabs(trace->EstimateWalk()) > 0.5 ) pass = false;
      if(nstep>=2) {
	if(chi2>ranges_fit_maxchi2/2.0) pass = false;
      }
      if(res!=0) pass = false;
      if(pass) {
	hFitChi2_sel->Fill( chi2 );
	hFitBaseline_sel  ->Fill( trace->EstimateBaseline() );
	hFitAmplitude_sel ->Fill( trace->EstimateAmplitude() );
	hFitWalk_sel      ->Fill( trace->EstimateWalk() );
	hFitBaseline2D_sel ->Fill( trace->EstimateBaseline(), meanPed );
	hFitAmplitude2D_sel->Fill( trace->EstimateAmplitude(), extreme_mV );
	hFitWalk2D_sel     ->Fill( trace->EstimateWalk(), extreme_ns );
	hFitWalk2DC_sel    ->Fill( trace->EstimateWalk(), trace->EstimateAmplitude() );
	trace->FillProfile(); // fill selected
	nSinglesTraces++;
	if(nExample<100) {
	  hExample[nExample] = (TH1D*) trace->Clone( Form("SelectedFitExample_%d_%d",nExample,nev) );
	  nExample++;
	}
      }
    } else {
      trace->FillProfile(); // fill selected
      nSinglesTraces++;
    }
  }
  cout << " Number of good singles used in profile " << nSinglesTraces << endl;
  new TCanvas();
  trace->DrawCopy();

  for(int i=0; i!=nExample; ++i) {
    if(hExample[i])
      list->Add( hExample[i] );
  }
  cout << "  Number of examples: " << nExample << endl;

  cout << " PULSE to " << filename.Data() << Form("_step%02d.dat",nstep) << endl;
  trace->SaveProfileToTemplate( Form("%s_step%02d.dat",filename.Data(),nstep), int(bins_per_ns*nstep) );

  TFile *foutroot = new TFile( Form("%s_step%02d.root",filename.Data(),nstep), "RECREATE" );
  summary->Write();
  hMeanPed->Write();
  hRMSPed->Write();
  lstNoCuts->Write("BeforeCuts",TObject::kSingleKey);
  lstCuts->Write("AfterCuts",TObject::kSingleKey);
  if(nstep>0) {
    lstFitNoCuts->Write("BeforeFitCuts",TObject::kSingleKey);
    lstFitCuts->Write("AfterFitCuts",TObject::kSingleKey);
    list->Write("FitExamples",TObject::kSingleKey);
  }  
  trace->GetProfile()->Write();
  foutroot->Close();

  return 0;
}
