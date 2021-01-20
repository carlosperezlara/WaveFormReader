#ifndef __WAVEFORM__
#define __WAVEFORM__
#include "TH1D.h"
#include "TProfile.h"
#include "TString.h"

class WaveForm : public TH1D {
 public:
  WaveForm(TString,TString,Int_t,Double_t,Double_t);
  ~WaveForm();

  // pasive methods: get info from data but does not modify them
  void ComputePedestal(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  void ComputeMax(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time);
  void ComputeMin(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time);
  void ComputeCentroid(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  void NonFitted() {fFitted=kFALSE;}
  
  // based on template
  void LoadTemplate(TString filename);
  Int_t FitTemplate(Double_t p0, Double_t p1, Double_t p2, Double_t r0, Double_t r1, Double_t r2, Double_t clip=0, TString option="0WWRMQ");
  Int_t FitTemplate(Double_t p0, Double_t p1, Double_t p2, Double_t r0, Double_t r1, Double_t r2, Double_t fit1=20, Double_t range2=100, TString option="0WWRMQ");
  Double_t EstimateBaseline();
  Double_t EstimateAmplitude();
  Double_t EstimateCharge();
  Double_t EstimateWalk();
  Double_t GetLastReducedChiSquared();
  
  // active methods: modifies data in some way
  void Subtract(Double_t ped);

  // profiling
  TProfile* CreateProfile();
  TProfile* GetProfile() {return fProfile;}
  void FillProfile();
  void SaveProfileToTemplate(TString filename, Int_t clip=0);
  
 private:
  // helpers
  void StatY(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  void StatX(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  void LocalExtreme(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time, bool FindMax);

  // datamembers
  TSpline3 *fPulse;
  TF1 *fFit;
  Bool_t fFitted;
  Double_t fTemplate_HeightToCharge;
  TProfile *fProfile;
};

//====================================
WaveForm::WaveForm(TString name, TString title, Int_t bins , Double_t minX, Double_t maxX) :
  TH1D(name,title,bins,minX,maxX),
  fPulse(NULL),
  fFit(NULL),
  fFitted(kFALSE),
  fTemplate_HeightToCharge(1),
  fProfile(NULL) {
}
//====================================
WaveForm::~WaveForm() {
  if(fPulse) delete fPulse;
  if(fFit) delete fFit;
  if(fProfile) delete fProfile;
}
//====================================
TProfile* WaveForm::CreateProfile() {
  if(fProfile) return fProfile;
  fProfile = new TProfile( Form("%s_pX",GetName()), Form("%s (PX)",GetTitle()), GetNbinsX(), GetBinLowEdge(1), GetBinLowEdge(GetNbinsX()+1) );
  return fProfile;
}
//====================================
void WaveForm::FillProfile() {
  if(!fProfile) CreateProfile();
  // FillProfile with the trace raw data
  // if available, use the found walk to
  // correct for jitter
  double jitter = 0;
  if(fFitted) {
    jitter = EstimateWalk();
  }
  for(Int_t bin = 0; bin!=GetNbinsX(); ++bin) {
    fProfile->Fill( GetBinCenter(bin+1)-jitter, GetBinContent(bin+1) );
  }
}
//====================================
void WaveForm::SaveProfileToTemplate(TString filename, Int_t clip) {
  if(!fProfile) return;
  fProfile->GetXaxis()->SetRange(clip,GetNbinsX()-clip);
  Double_t scale = 1.0/abs(fProfile->GetMinimum());
  std::cout << "WAVEFORM :: SCALE " << scale << std::endl;
  ofstream fout( filename.Data() );
  for(Int_t bin = 0+clip; bin<GetNbinsX()-clip; ++bin) {
    
    fout <<  Form("%e %e", fProfile->GetBinCenter(bin+1), fProfile->GetBinContent(bin+1)*scale ) << endl;
  }
  fout.close();
}
//====================================
void WaveForm::ComputeCentroid(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms) {
  StatX(binMin, binMax, mean, rms);
}
//====================================
void WaveForm::ComputeMax(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time) {
  LocalExtreme(binMin, binMax, ampl, time, true);
}
//====================================
void WaveForm::ComputeMin(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time) {
  LocalExtreme(binMin, binMax, ampl, time, false);
}
//====================================
void WaveForm::ComputePedestal(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms) {
  StatY(binMin, binMax, mean, rms);
}
//====================================
void WaveForm::Subtract(Double_t ped) {
  for(int i=0; i!=this->GetNbinsX(); ++i) {
    this->SetBinContent(i+1, this->GetBinContent(i+1) - ped );
  }
}
//====================================
void WaveForm::LoadTemplate(TString filename) {
  TGraph *gr = new TGraph(filename.Data());
  fPulse = new TSpline3("pulse",gr);
  delete gr;
  fFit = new TF1("fFit", [&](double*x, double *p){ return p[0] + p[1] * fPulse->Eval(x[0] + p[2]); },
		 this->GetBinLowEdge(1),
		 this->GetBinLowEdge( this->GetNbinsX()+1 ),
		 3);
  fTemplate_HeightToCharge = 1/50; // mV*ns / 50 Ohms ==> pC
}
//====================================
Int_t WaveForm::FitTemplate(Double_t p0, Double_t p1, Double_t p2, // base ampl walk
			    Double_t r0, Double_t r1, Double_t r2, // windows
			    Double_t minR, Double_t maxR, TString option) {
  fFitted = kTRUE;
  fFit->SetParameters(p0,p1,p2);
  fFit->SetParLimits(0,p0-r0,p0+r0);
  fFit->SetParLimits(1,p1-r1,p1+r1);
  fFit->SetParLimits(2,p2-r2,p2+r2);
  this->Fit(fFit,option.Data(),"",minR,maxR);
  if( ( fFit->GetParameter(0)-1e-6 < p0-r0 ) ||
      ( fFit->GetParameter(0)+1e-6 > p0+r0 ) ||
      ( fFit->GetParameter(1)-1e-6 < p1-r1 ) ||
      ( fFit->GetParameter(1)+1e-6 > p1+r1 ) ||
      ( fFit->GetParameter(2)-1e-6 < p2-r2 ) ||
      ( fFit->GetParameter(2)+1e-6 > p2+r2 )  )
    return 1;
  return 0;
}
//====================================
Int_t WaveForm::FitTemplate(Double_t p0, Double_t p1, Double_t p2, // base ampl walk
			    Double_t r0, Double_t r1, Double_t r2, // windows
			    Double_t clip, TString option) {
  Double_t minR = GetBinLowEdge( 1 ) + clip;
  Double_t maxR = GetBinLowEdge( GetNbinsX()+1 ) - clip;
  return FitTemplate(p0, p1, p2, // base ampl walk                                                                          
		     r0, r1, r2, // windows                                                                                 
		     minR, maxR, option);
}
//====================================
Double_t WaveForm::GetLastReducedChiSquared() {
  if(fFit->GetNDF()!=0)
    return fFit->GetChisquare()/fFit->GetNDF();
  else
    return -9999;
}
//====================================
Double_t WaveForm::EstimateBaseline() {
  return fFit->GetParameter(0);
}
//====================================
Double_t WaveForm::EstimateAmplitude() {
  return fFit->GetParameter(1);
}
//====================================
Double_t WaveForm::EstimateWalk() {
  return fFit->GetParameter(2);
}
//====================================
Double_t WaveForm::EstimateCharge() {
  return EstimateAmplitude() * fTemplate_HeightToCharge;
}

//====================================
void WaveForm::StatY(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms) {
  Int_t range = binMax - binMin;
  if(range<0) {
    std::cout << "WaveForm::StatY( bad range )!" << std::endl;
    return;
  }
  double sx=0;
  double sxx=0;
  for(int i=binMin; i!=binMax; ++i) {
    double val = this->GetBinContent(i);
    sx += val;
    sxx += val*val;
  }
  mean = sx / range;
  rms = TMath::Sqrt( sxx/range - mean*mean );
  return;
}
//====================================
void WaveForm::StatX(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms) {
  Int_t range = binMax - binMin;
  if(range<0) {
    std::cout << "WaveForm::StatX( bad range )!" << std::endl;
    return;
  }
  double s=0;
  double sx=0;
  double sxx=0;
  for(int i=binMin; i!=binMax; ++i) {
    double yval = this->GetBinContent( i );
    double xval = this->GetBinCenter( i );
    s += yval;
    sx += yval * xval;
    sxx += yval * xval*xval;
  }
  mean = sx / s;
  rms = sxx/s - mean*mean;
  if(rms>0) rms = TMath::Sqrt(rms);
  else rms = -TMath::Sqrt(-rms);
  return;
}
//====================================
void WaveForm::LocalExtreme(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time, bool Positive) {
  Int_t range =	binMax - binMin;
  if(range<0) {
    std::cout << "WaveForm::LocalExtreme( bad range )!" << std::endl;
    return;
  }
  time=0;
  ampl=-1e16;
  if(!Positive) ampl = 1e16;
  for(int i=binMin; i!=binMax; ++i) {
    double val = this->GetBinContent(i);
    bool accept = false;
    if(Positive) {
      accept = (val>ampl);
    } else {
      accept = (val<ampl);
    }
    if(accept) {
      ampl = val;
      time = i;
    }
  }
  return;
}

#endif
