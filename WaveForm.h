#ifndef __WAVEFORM__
#define __WAVEFORM__

#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TProfile.h"
#include "TString.h"
#include "TSpline.h"
#include "TF1.h"
#include "TMath.h"

class WaveForm : public TH1D {
 public:
  WaveForm(TString,TString,Int_t,Double_t,Double_t);
  ~WaveForm();

  // pasive methods: get info from data but does not modify them
  void ComputePedestal(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  void ComputeMax(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time);
  void ComputeMin(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time);
  void ComputeCentroid(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  inline void NonFitted() {fFitted=kFALSE;}
  Double_t GetKthSquared(Int_t k);
  inline TH1D* GetDISTrace() {return fTransformed;}
  TH1D* DelayInvertAndSubtract(Int_t delta);
  Double_t FindCrossing(Int_t binMin, Int_t binMax, Double_t thr, bool Positive);
  
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
  //
  TProfile* CreateProfileDFT();
  inline TProfile* GetProfileDFT() {return fProfileDFT;}
  void FillProfileDFT(); // extremely expensive
  
 private:
  // helpers
  void StatY(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  void StatX(Int_t binMin, Int_t binMax, Double_t &mean, Double_t &rms);
  void LocalExtreme(Int_t binMin, Int_t binMax, Double_t &ampl, Int_t &time, bool FindMax);

  // datamembers
  TSpline3 *fPulse;
  TH1D *fTransformed;
  TF1 *fFit;
  Bool_t fFitted;
  Double_t fTemplate_HeightToCharge;
  TProfile *fProfile;
  TProfile *fProfileDFT;

  ClassDef(WaveForm, 0);
};

#endif
