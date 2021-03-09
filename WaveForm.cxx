#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TProfile.h"
#include "TString.h"
#include "TSpline.h"
#include "TF1.h"
#include "TMath.h"
#include "WaveForm.h"

//====================================
TH1D* WaveForm::DelayInvertAndSubtract(Int_t delta) {
  if(!fTransformed) fTransformed = (TH1D*) this->Clone( Form("%s_DIS",GetName()) );
  fTransformed->Reset();
  for(int i=delta; i!=GetNbinsX(); ++i) {
    fTransformed->SetBinContent(i+1, GetBinContent(i+1)-GetBinContent(i+1-delta) );
  }
  return fTransformed;
}
//====================================
WaveForm::WaveForm(TString name, TString title, Int_t bins , Double_t minX, Double_t maxX) :
  TH1D(name,title,bins,minX,maxX),
  fPulse(NULL),
  fTransformed(NULL),
  fFit(NULL),
  fFitted(kFALSE),
  fTemplate_HeightToCharge(1),
  fProfile(NULL),
  fProfileDFT(NULL) {
}
//====================================
WaveForm::~WaveForm() {
  if(fPulse) delete fPulse;
  if(fTransformed) delete fTransformed;
  if(fFit) delete fFit;
  if(fProfile) delete fProfile;
  if(fProfileDFT) delete fProfileDFT;
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
TProfile* WaveForm::CreateProfileDFT() {
  if(fProfileDFT) return fProfileDFT;
  Int_t N = GetNbinsX();
  // turns out that computing DFT is computationaly expensive (who knows!)
  // and looking into a full map between nsample -> nkths may not be needed
  // specially when samples are higher than 100 (which is most of the case)
  // so we reduce the powers we look at by rebining to match 100 or N/2
  // whichever is smallest.
  Int_t kbins = TMath::Min( 100, N/2 );
  Int_t kmax = TMath::Floor(N/kbins/2) * kbins; // making sure each bincenter is INT
  /*
  Double_t min = -N/2-0.5;
  Double_t max = +N/2-1.5;
  fProfileDFT = new TProfile( Form("%s_pXDFT",GetName()), Form("%s DFT (PX)",GetTitle()), N, min, max );
  */
  fProfileDFT = new TProfile( Form("%s_pXDFT",GetName()), Form("%s DFT (PX)",GetTitle()),
			      kbins, -0.5, kmax-0.5 );
  Double_t timeWindow = GetBinLowEdge( N+1 ) - GetBinLowEdge(1);
  for(int i=0; i!=kbins; ++i) {
   if(i%10==0) fProfileDFT->GetXaxis()->SetBinLabel( i+1, Form("%.2f GHz", fProfileDFT->GetBinCenter(i+1)/timeWindow) );
  }
  return fProfileDFT;
}
//====================================
Double_t WaveForm::GetKthSquared(Int_t k) {
  // Computes the k-th power of the
  // Discrete Fourier Transform of
  // the current trace
  // WARNING: returns the Power Squared
  // (sqrt is expensive)
  Int_t N = GetNbinsX();
  static Double_t TwoPiOverN = TMath::TwoPi()/N;
  double re = 0;
  double im = 0;
  for(int n=0; n!=N; ++n) {
    double x = GetBinContent( n+1 );
    re += x*TMath::Cos( TwoPiOverN*k*n );
    im += x*TMath::Sin( TwoPiOverN*k*n );
  }
  return re*re + im*im;
}

//====================================
void WaveForm::FillProfileDFT() {
  if(!fProfileDFT) CreateProfileDFT();
  // FillProfileDFT with DFT of trace raw data
  Int_t N = fProfileDFT->GetNbinsX();
  for(int kb=0; kb!=N; ++kb) {
    Int_t k = fProfileDFT->GetBinCenter( kb+1 );
    double pow2 = GetKthSquared( k );
    //fProfileDFT->Fill( k>N/2?k-N:k, pow2 );
    fProfileDFT->Fill( k, pow2 );
  }
}
//====================================
void WaveForm::SaveProfileToTemplate(TString filename, Int_t clip) {
  if(!fProfile) return;
  fProfile->GetXaxis()->SetRange(clip,GetNbinsX()-clip);
  Double_t scale = 1.0/abs(fProfile->GetMinimum());
  std::cout << "WAVEFORM :: SCALE " << scale << std::endl;
  std::ofstream fout( filename.Data() );
  for(Int_t bin = 0+clip; bin<GetNbinsX()-clip; ++bin) {
    
    fout <<  Form("%e %e", fProfile->GetBinCenter(bin+1), fProfile->GetBinContent(bin+1)*scale ) << std::endl;
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
//====================================
Double_t WaveForm::FindCrossing(Int_t binMin, Int_t binMax, Double_t thr, bool Positive) {
  Int_t range =	binMax - binMin;
  if(range<0) {
    std::cout << "WaveForm::LocalExtreme( bad range )!" << std::endl;
    return -99998;
  }
  int time=0;
  // first iteration
  for(int i=binMin; i!=binMax; ++i) {
    double val = GetBinContent(i);
    //std::cout << ">> " << i << " << " << val << " " << thr << std::endl;
    bool accept = false;
    if(Positive) {
      accept = (val>thr);
    } else {
      accept = (val<thr);
    }
    if(accept) {
      time = i;
      //std::cout << "FOUND!" << std::endl;
      break;
    }
  }
  if(time==0) {
    //std::cout << "****** NOT FOUND!!!!!!!" << Form("[%d,%d] thr %f",binMin,binMax,thr) << std::endl;
    return -99997; // failed
  }
  // second iteration
  // interpolate
  double y0 = GetBinContent(time-1);
  double y1 = GetBinContent(time);
  double bw = GetBinWidth(time);
  double dx = (thr-y0)*bw/(y1-y0);
  Double_t ret = GetBinCenter(time)-0.5*bw+dx;
  return ret;
}
