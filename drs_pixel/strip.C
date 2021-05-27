struct pulse_t {
  Int_t event;
  Short_t tc[2];
  Float_t channel[18][1024];
  Float_t times[2][1024];
  Float_t xIntercept;
  Float_t yIntercept;
  Float_t xSlope;
  Float_t ySlope;
  Float_t chi2;
  Int_t   ntracks;
} pulse;
  
int strip(int run) {
  TFile *file = new TFile(Form("merged/Run_%d.root",run));
  TTree *tree = (TTree*) file->Get("pulse");
  tree->SetBranchAddress("channel",&pulse.channel);
  tree->SetBranchAddress("times",&pulse.times);
  tree->SetBranchAddress("xIntercept",&pulse.xIntercept);
  tree->SetBranchAddress("yIntercept",&pulse.yIntercept);
  tree->SetBranchAddress("xSlope",&pulse.xSlope);
  tree->SetBranchAddress("ySlope",&pulse.ySlope);
  tree->SetBranchAddress("chi2",&pulse.chi2);
  tree->SetBranchAddress("ntracks",&pulse.ntracks);

  TCanvas *main = new TCanvas();
  main->Divide(7,7);
  TGraph *gr[49];
  Long_t nevts = tree->GetEntries();
  ofstream fout[17];
  for(int i=0; i!=16; ++i) {
    fout[i].open( Form("inCSV/RUN%d_CH%d.csv",run,i) );
  }
  fout[16].open( Form("inCSV/RUN%d_PIXEL.csv",run) );
  for(int i=0; i!=nevts; ++i) {
    tree->GetEntry(i);
    if(i<49) {
      gr[i] = new TGraph( 1024, pulse.times[0], pulse.channel[2] );
      main->cd(1+i);
      gr[i]->Draw("AL");
    }
    if((i%500)==0)
      cout << "Events processed: " << i << endl;
    for(int ich=0; ich!=18; ++ich) {
      int xch = ich;
      if(ich==8 || ich==17) continue; // skipping TR0 (8, 17)
      int group = ich/9;
      if(group==1) xch = ich - 1; // 9 => 8, 16 => 15
      for(int isa=0; isa!=1024; ++isa) {
	fout[xch] << pulse.times[group][isa] << ", ";
	fout[xch] << pulse.channel[ich][isa]; //fixed typo (may 4th)
	if(isa<1023) fout[xch] << ", ";
      }
      fout[xch] << endl;
    }
    fout[16] << pulse.xIntercept << ", ";
    fout[16] << pulse.yIntercept << ", ";
    fout[16] << pulse.xSlope << ", ";
    fout[16] << pulse.ySlope << ", ";
    fout[16] << pulse.chi2 << ", ";
    fout[16] << pulse.ntracks << endl;
  }
  
  
  return 0;
}
