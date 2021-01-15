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
  fin >> ranges_baseline >> ranges_gate_low >> ranges_gate_high >> ranges_extreme_low >> ranges_extreme_high;
  fin.close();
  cout << "   => RANGES:BASELINE " << ranges_baseline << endl;
  cout << "   => RANGES:GATE_LOW " << ranges_gate_low << endl;
  cout << "   => RANGES:GATE_HIGH " << ranges_gate_high << endl;
  cout << "   => RANGES:EXTREME_LOW " << ranges_extreme_low << endl;
  cout << "   => RANGES:EXTREME_HIGH " << ranges_extreme_high << endl;

  double baseline=0;
  if(nstep>0) {    // STEP 1++: REMOVE BASELINE
    cout << " BASELINE from " << filename.Data() << Form("_step%02d.dat",0) << endl;
    fin.open( Form("%s_step00.dat",filename.Data()) );
    fin >> baseline;
    fin.close();
  }

  if(nstep>1) {    // STEP2++: FIT PULSE
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
  TH1D *hMeanPed = new TH1D("hMeanPed","MeanPed;mV",100,min_mV-baseline,max_mV-baseline);
  TH1D *hRMSPed = new TH1D("hRMSPed","RMSPed;mV",100,0,+15);

  Double_t extreme_mV;
  Int_t extreme_timebin;
  //TH1D *hMinmV = new TH1D("hMinmV","MinmV;mV",1000,min_mV-baseline,max_mV-baseline);
  TH1D *hMinmV = new TH1D("hMinmV","MinmV;mV",1000,min_mV-baseline,0);
  TH1D *hMinns = new TH1D("hMinns","MinnS;ns",1000,min_ns,max_ns);

  TH1D *hFitBaseline  = new TH1D("hFitBaseline", "Fit_Baseline;mV",  100,min_mV-baseline,max_mV-baseline);
  TH1D *hFitAmplitude = new TH1D("hFitAmplitude","Fit_Amplitude;ns", 1000,min_mV-baseline,0);
  TH1D *hFitWalk      = new TH1D("hFitWalk",     "Fit_Walk;mV",      1000,min_ns,max_ns);
  
  trace->CreateProfile();
  Int_t nSinglesTraces = 0;
  for(int nev=0;;++nev) {
    if(!file.ReadEvent()) break;
    if(nev%500==0)
      cout << "Events read so far: " << nev << endl;

    trace->Subtract( baseline );

    trace->ComputePedestal(1,ranges_baseline*bins_per_ns,meanPed,rmsPed);
    hMeanPed->Fill( meanPed );
    hRMSPed->Fill( rmsPed );

    trace->ComputeMin(ranges_gate_low*bins_per_ns,ranges_gate_high*bins_per_ns,extreme_mV,extreme_timebin);
    hMinmV->Fill( extreme_mV );
    Double_t extreme_ns = trace->GetXaxis()->GetBinCenter(extreme_timebin);
    hMinns->Fill( extreme_ns );

    if( (extreme_mV>ranges_extreme_low) && (extreme_mV<ranges_extreme_high) ) {

      if(nstep>1) {    // STEP2++: FIT PULSE
	trace->FitTemplate(0,extreme_mV,extreme_ns);
	hFitBaseline  ->Fill( trace->EstimateBaseline());
	hFitAmplitude ->Fill( trace->EstimateAmplitude());
	hFitWalk      ->Fill( trace->EstimateWalk());
      }

      trace->FillProfile();
      nSinglesTraces++;
    }
  }
  cout << " Number of good singles used in profile " << nSinglesTraces << endl;
  new TCanvas();
  trace->DrawCopy();

  hMeanPed->Fit("gaus","I");
  if(nstep==0) { // STEP 0: SAVING BASELINE
    ofstream fout( Form("%s_step%02d.dat",filename.Data(),nstep) );
    fout << Form("%e",((TF1*) hMeanPed->GetListOfFunctions()->At(0))->GetParameter(1)) << endl;
    fout.close();
  }

  TFile *foutroot = new TFile( Form("%s_step%02d.root",filename.Data(),nstep), "RECREATE" );
  summary->Write();
  hMeanPed->Write();
  hRMSPed->Write();
  hMinmV->Write();
  hMinns->Write();
  trace->GetProfile()->Write();
  foutroot->Close();
  if(nstep>0) {
    cout << " PULSE to " << filename.Data() << Form("_step%02d.dat",nstep) << endl;
    trace->SaveProfileToTemplate( Form("%s_step%02d.dat",filename.Data(),nstep));
  }
  if(nstep>1) {
    hFitBaseline->Write();
    hFitAmplitude->Write();
    hFitWalk->Write();
  }  
  return 0;
}
