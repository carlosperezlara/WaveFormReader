const double firstPeakG  = 0.5110;
const double secondPeakG = 1.2745;
const int CHANNELS = 4;
const int SAMPLES = 1024;
bool debug = false;

int waveTrainDRS(TString inputfile="20200930_DRS/9-30-2020-4ch-SiPM-LYSO-69V.dat", int dofit=0, float xc0=-5, float xc1=-5, float xc2=-3, float xc3=-3) {
  float xintmax[4] = { xc0, xc1, xc2, xc3 };
  std::cout << "=====================" << std::endl;
  std::cout << "= Source file: " << inputfile.Data() << std::endl;
  std::cout << "=====================" << std::endl;

  std::ifstream inf( inputfile.Data() , std::ifstream::binary);
  
  gStyle->SetOptStat(0);
  int col = 0;
  int color[4] = { kOrange-3, kCyan-3, kMagenta-3, kGreen-3 };
  int nevt = 0;

  char  *c4Tmp = new char [4];
  char  *c2Tmp = new char [2];
  int   iTmp;
  //RUN HEADER
  inf.read(c4Tmp, 4);
  if(inf.rdstate()!=0) {
    std::cout << " Error while reading file. Bye bye" << std::endl;
    return 1;
  }
  std::cout << c4Tmp << " | "; // DRSN
  inf.read(c4Tmp, 4);
  std::cout << c4Tmp << " | "; // TIME
  inf.read(c2Tmp, 2);
  std::cout << c2Tmp << " | "; // B#
  inf.read((char*) &iTmp, 2);
  std::cout << iTmp << std::endl; // BoardNumber
  //looping over each channel getting time
  float tch[CHANNELS][SAMPLES];  // Timing calibration
  float time[CHANNELS][SAMPLES];
  for(int ich=0; ich!=CHANNELS; ++ich) {
    inf.read(c4Tmp, 4);
    std::cout << "CHN" << ich << ": " << c4Tmp << std::endl; // C00i
    inf.read((char*) &tch[ich], 4*SAMPLES); // Effective time bin width in ns encoded in 4-byte floating point format
    // moving to axis
    for(int isa=0; isa!=SAMPLES; ++isa) {
      time[ich][isa] = tch[ich][isa];
    }
    for(int isa=1; isa!=SAMPLES; ++isa) {
      time[ich][isa] += time[ich][isa-1];
    }
    std::cout << "   TIME [ " << time[ich][0] << ", " << time[ich][SAMPLES-1] << " ]" << std::endl;
  }
  
  double minX = time[0][0];
  double maxX = time[0][1023];
  TH1D* fTmp = new TH1D( "Wave", "Wave;Time  (ns);Signal  (mV)", SAMPLES, minX, maxX );
  TH2D *axis[CHANNELS];
  TH1D *integral[CHANNELS];
  TH1D *amplitude[CHANNELS];
  TH1D *integralControl[CHANNELS];
  TH1D *amplitudeControl[CHANNELS];
  for(int ich=0; ich!=CHANNELS; ++ich) {
    axis[ich] = new TH2D(Form("axis_%d",ich),Form("File: %s  Channel %d;Time  (ns);Signal  (mV)",inputfile.Data(),ich),512,minX,maxX,150,-550,+550);
    integral[ich] = new TH1D(Form("charge_%d",ich),Form("File: %s  Channel %d;charge  (nC)",inputfile.Data(),ich),100, xintmax[ich] ,0);
    amplitude[ich]= new TH1D(Form("amplitude_%d",ich),Form("File: %s  Channel %d;charge  (nC)",inputfile.Data(),ich),100, xintmax[ich]/300 ,0);
    integralControl[ich] = new TH1D(Form("chargeControl_%d",ich),Form("File: %s  Channel %d;fraction",inputfile.Data(),ich),100,0.6,1.3);
    amplitudeControl[ich]= new TH1D(Form("amplitudeControl_%d",ich),Form("File: %s  Channel %d;fraction",inputfile.Data(),ich),100,0.0,0.01);
    integral[ich]->SetLineColor(color[ich]);
    integral[ich]->SetMarkerColor(color[ich]);
    integral[ich]->SetMarkerStyle(24);
    integral[ich]->Sumw2();
  }

  int corrIdx[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
  TH2D *corrHis[6];
  for(int ico=0; ico!=6; ++ico) {
    corrHis[ico] = new TH2D( Form("corr_%d_%d", corrIdx[ico][0], corrIdx[ico][1]),
			     Form("File: %s;charge_ch%d   (pC);charge_ch%d   (pC)",inputfile.Data(), corrIdx[ico][0], corrIdx[ico][1]),
			     100,xintmax[corrIdx[ico][0]],0,100,xintmax[corrIdx[ico][1]],0);
  }

  std::cout << "*****" << std::endl;
  //TRAILER
  char *c5Tmp = new char[5];
  short s8Tmp[8];
  ushort ch[CHANNELS][SAMPLES];
  double chmV[CHANNELS][SAMPLES]; //This is the voltage recorded on channel 1 in mV
  int nsamp = 50;
  double meanrms[4] = {0,0,0,0};
  while(!inf.eof()) {
    if (nevt%10000==0) {
      std::cout << "Processing event #" << nevt << " ";
      std::cout << inf.rdstate() << std::endl;
    }
    //EVENT HEADER
    inf.read(c4Tmp, 4);
    if(debug) std::cout << c4Tmp << " | "; //EHDR
    inf.read((char *) &iTmp, 4);
    if(debug) std::cout << iTmp << " | ["; //EventNumber
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << "/"; //Year
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << "/"; //Month
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << " "; //Day
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << ":"; //Hour
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << ":"; //Minute
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << "::"; //Second
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << "] | "; //Milisecond
    inf.read((char *) &iTmp, 2);
    if(debug) std::cout << iTmp << " | "; //RangeCenter in mV
    int rc = iTmp;
    inf.read(c2Tmp, 2);
    if(debug) std::cout << c2Tmp << " "; // B#
    inf.read((char*) &iTmp, 2);
    if(debug) std::cout << iTmp << " | "; // BoardNumber
    inf.read(c2Tmp, 2);
    if(debug) std::cout << c2Tmp << " "; // T#
    inf.read((char*) &iTmp, 2);
    if(debug) std::cout << iTmp << std::endl; // TriggerCell
    int tcell = iTmp;
    //Reading the values of each sample in the buffer for an event
    double charge[CHANNELS];
    for(int ich=0; ich!=CHANNELS; ++ich) {
      inf.read(c4Tmp, 4); // C00i
      inf.read((char*) &iTmp, 4); // Scaler in Hz
      inf.read((char*) &ch[ich], SAMPLES*2); // SAMPLES*2 (2-byte integer. 0=RC-0.5V and 65535=RC+0.5V. RC see header.
      //Converting raw value to mV
      double sxx = 0;
      for (int isa=0; isa!=SAMPLES; ++isa) {
	chmV[ich][isa] = 1000*(ch[ich][isa]/65535.-0.5)+rc;
	fTmp->SetBinContent(isa+1,chmV[ich][isa]);
	axis[ich]->Fill( time[ich][isa], chmV[ich][isa] );
	/*
	float timeval = 0;
	for(int j=0; j<isa-1; ++j)
	  timeval += tch[ich][ (j+tcell)%1024 ];
	axis[ich]->Fill( timeval, chmV[ich][isa] );
	*/
	if(isa<nsamp) {
	  sxx += chmV[ich][isa]*chmV[ich][isa];
	}
	//std::cout << "******> " << chmV[ich][isa] << std::endl;
      }
      meanrms[ich] += sqrt( (sxx/nsamp) ); // adding rms of wave x
      int binCenter = fTmp->GetMinimumBin();
      double valA = fTmp->GetBinContent( binCenter )*fTmp->GetBinWidth(1) * 1e-3 / 50;
      double FullInt = fTmp->Integral(1,fTmp->GetNbinsX()) * 1e-3 / 50;
      amplitude[ich]->Fill( valA );
      amplitudeControl[ich]->Fill( valA/FullInt );
      if(binCenter<200) continue;
      if(binCenter>1024-300) continue;
      int binMin = binCenter -200;// approx 40ns;   // was -200
      int binMax = binCenter +300;// approx 60ns;   // was +300
      double valI = fTmp->Integral( binMin, binMax ) * 1e-3 / 50; // nS*mV/Ohm*1e-3 -> nC
      charge[ich] = valI;
      integral[ich]->Fill( valI );
      integralControl[ich]->Fill( valI/FullInt );
    }
    if(CHANNELS==4) {
      for(int ico=0; ico!=6; ++ico) {
	corrHis[ico]->Fill( charge[ corrIdx[ico][0] ], charge[ corrIdx[ico][1] ] );
      }
    }

    nevt++;
  }
  inf.close();
  std::cout << "FOUND " << nevt << " EVENTS." << std::endl;
  
  for(int i=0; i!=CHANNELS; ++i)
    meanrms[i] /= nevt;

  TLatex *tex = new TLatex();
  TBox *lin = new TBox();

  TCanvas *main = new TCanvas();
  main->Divide(2,2);
  for(int i=0; i!=CHANNELS; ++i) {
    main->cd(i+1)->SetLogz(1);
    axis[i]->Draw("colz");
    tex->SetTextSize(0.04);
    tex->SetTextColor( kBlack );
    tex->DrawLatexNDC(0.13,0.18, Form("< RMS_{%d} > = %.2f mV",nsamp,meanrms[i]) );
    double nsamp_x = fTmp->GetXaxis()->GetBinLowEdge(nsamp+1);
    cout << " ***** " << minX << " |||| " << nsamp_x << endl;
    lin->SetFillColor(kBlack);
    lin->DrawBox(minX,-460,nsamp_x,-470);
  }
  //main->SaveAs( Form("%s_OSC.png",file.Data()), "PNG" );

  TCanvas *main2 = new TCanvas();
  main2->Divide(2,2);
  TF1 *fit[CHANNELS];
  TF1 *fit1[CHANNELS];
  TF1 *fit2[CHANNELS];
  TF1 *fit3[CHANNELS];
  TF1 *fit4[CHANNELS];
  TF1 *fitA[CHANNELS];
  for(int i=0; i!=CHANNELS; ++i) {
    main2->cd(i+1);//->SetLogy(1);
    integral[i]->GetYaxis()->SetRangeUser(0,integral[i]->GetMaximum()*1.5);
    integral[i]->Draw("E");
    if(dofit>0) {
      //double scale = nevt*integral[i]->GetBinWidth(1);
      double scale = integral[i]->GetEntries()*integral[i]->GetBinWidth(1);
      std::cout << " *****> " << scale << std::endl;
      if(dofit>=2) {
	fit[i] = new TF1(Form("fit_%d",i),Form("[0]*([4]*(100/280.*TMath::Gaus(x,[1],[2],1)+180/280.*TMath::Gaus(x,[1]*%f,[3],1))+(1-[4])*(100/280.*TMath::Gaus(x,[5],[2]*[5]/[1],1)+180/280.*TMath::Gaus(x,[5]*%f,[3]*[5]/[1],1)))",
					       firstPeakG/secondPeakG,firstPeakG/secondPeakG));//,xintmax[i],0);
	fit[i]->SetParNames("evts","mu1","sgG1","sgG2","na","muB1");
	fit[i]->SetParameter(0,scale);                fit[i]->SetParLimits(0,scale*0.6,scale*1.2);                        //eN
	fit[i]->SetParameter(1,0.80*xintmax[i]);      fit[i]->SetParLimits(1,0.90*xintmax[i],0.60*xintmax[i]);            //mu1
	fit[i]->SetParameter(2,0.05*abs(xintmax[i])); fit[i]->SetParLimits(2,0.04*abs(xintmax[i]),0.07*abs(xintmax[i]));  //sigmaG1
	fit[i]->SetParameter(3,0.05*abs(xintmax[i])); fit[i]->SetParLimits(3,0.04*abs(xintmax[i]),0.18*abs(xintmax[i]));  //sigmaG2
	fit[i]->SetParameter(4,0.8);                  fit[i]->SetParLimits(4,0.6,0.99);                                    //na
	fit[i]->SetParameter(5,0.40*xintmax[i]);      fit[i]->SetParLimits(5,0.55*xintmax[i],0.30*xintmax[i]);            //muB1
      } else if(dofit==1) {
	fit[i] = new TF1(Form("fit_%d",i),"[0]*180./280.*TMath::Gaus(x,[1],[2],1)", xintmax[i],0);
	//fit[i] = new TF1(Form("fit_%d",i),"[0]*180./280.*TMath::Gaus(x,[1],[2],1)");
	fit[i]->SetParNames("evts","muG1","sgG2");
	fit[i]->SetParameter(0,0.5*scale);            fit[i]->SetParLimits(0,scale*0.2,1.2*scale);                        //eN
	fit[i]->SetParameter(1,0.80*xintmax[i]);      fit[i]->SetParLimits(1,0.90*xintmax[i],0.60*xintmax[i]);            //muG2
	fit[i]->SetParameter(2,0.05*abs(xintmax[i])); fit[i]->SetParLimits(2,0.03*abs(xintmax[i]),0.09*abs(xintmax[i]));  //sigmaG2
        fitA[i] = new TF1(Form("fitA_%d",i),"[0]*TMath::Gaus(x,[1],[2],1)", xintmax[i]/100,0);
        fitA[i]->SetParNames("evts","muG1","sgG1");
      }
      fit[i]->SetLineStyle(9);
      if(dofit==1) integral[i]->Fit(fit[i],"lqr","",xintmax[i],0.50*xintmax[i]);
      else         integral[i]->Fit(fit[i],"lqr","",xintmax[i],0.12*xintmax[i]);
      double chi2 =  fit[i]->GetChisquare() / fit[i]->GetNDF();
      for(int k=0;k!=1;++k) {
	chi2 =  fit[i]->GetChisquare() / fit[i]->GetNDF();
	cout << "last chi2 " << chi2 << endl;
	if(chi2<2) break;
	if(dofit==1) integral[i]->Fit(fit[i],"emiqr","",
				      fit[i]->GetParameter(1)-fit[i]->GetParameter(2),
				      fit[i]->GetParameter(1)+fit[i]->GetParameter(2));
	else         integral[i]->Fit(fit[i],"emiqr","",xintmax[i],0.12*xintmax[i]);
      }
      if(dofit==1) integral[i]->Fit(fit[i],"emir","",
				    fit[i]->GetParameter(1)-fit[i]->GetParameter(2),
				    fit[i]->GetParameter(1)+fit[i]->GetParameter(2));
      else         integral[i]->Fit(fit[i],"emir","",xintmax[i],0.12*xintmax[i]);
      if(dofit>=2) {
	fit1[i] = new TF1(Form("firstCompG_%d",i), Form("%f*TMath::Gaus(x,%f,%f,1)",
							fit[i]->GetParameter(0)*fit[i]->GetParameter(4)*100./280.,
							fit[i]->GetParameter(1),
							fit[i]->GetParameter(2)),xintmax[i],0);
	fit2[i] = new TF1(Form("seconCompG_%d",i), Form("%f*TMath::Gaus(x,%f,%f,1)",
							fit[i]->GetParameter(0)*fit[i]->GetParameter(4)*180./280.,
							fit[i]->GetParameter(1)*firstPeakG/secondPeakG,
							fit[i]->GetParameter(3)),xintmax[i],0);
	fit3[i] = new TF1(Form("firstCompB_%d",i), Form("%f*TMath::Gaus(x,%f,%f,1)",
							fit[i]->GetParameter(0)*(1-fit[i]->GetParameter(4))*100./280.,
							fit[i]->GetParameter(5),
							fit[i]->GetParameter(2)*fit[i]->GetParameter(5)/fit[i]->GetParameter(1)),xintmax[i],0);
	fit4[i] = new TF1(Form("secondCompB_%d",i),Form("%f*TMath::Gaus(x,%f,%f,1)",
							fit[i]->GetParameter(0)*(1-fit[i]->GetParameter(4))*180./280.,
							fit[i]->GetParameter(5)*firstPeakG/secondPeakG,
							fit[i]->GetParameter(3)*fit[i]->GetParameter(5)/fit[i]->GetParameter(1)),xintmax[i],0);
	fit1[i]->SetLineColor( kBlack );
	fit2[i]->SetLineColor( kBlack );
	fit3[i]->SetLineColor( kGray );
	fit4[i]->SetLineColor( kGray );
	fit4[i]->Draw("same");
	fit3[i]->Draw("same");
	fit2[i]->Draw("same");
	fit1[i]->Draw("same");
      }
      fit[i]->Draw("same");
      tex->SetTextColor( kRed-3 );
      tex->SetTextSize(0.045);
      if(dofit==2) {
	tex->DrawLatexNDC(0.12,0.82, "eN ( ");
	tex->DrawLatexNDC(0.17,0.82, Form("Na22 ( #frac{100}{280} gauss( #mu_{1.3}, #sigma_{1.3} ) + #frac{180}{280} gaus( #mu_{1.3} #frac{%.4f}{%.3f}, #sigma_{0.5}) ) + ",firstPeakG,secondPeakG));
	tex->DrawLatexNDC(0.15,0.72, Form("(1-Na22) ( #frac{100}{280} gaus( #mu_{B1.3}, #sigma_{1.3} #frac{#mu_{1.3}}{#mu_{B1.3}} ) + #frac{180}{280} gaus( #mu_{B1.3} #frac{%.4f}{%.3f}, #sigma_{0.5} #frac{#mu_{1.3}}{#mu_{B1.3}} ))",firstPeakG,secondPeakG));

      } if (dofit==3) {
	tex->SetTextSize(0.04);
	tex->SetTextColor( kBlack );
	tex->DrawLatexNDC(0.13,0.82, Form( "#mu_{1.3} = %.0f nC  --> %.0f keV",   fit[i]->GetParameter(1)*1e3, fit[i]->GetParameter(1)*secondPeakG/fit[i]->GetParameter(1)*1000));
	tex->DrawLatexNDC(0.13,0.75, Form( "#sigma_{1.3} = %.0f nC  --> %.0f keV",fit[i]->GetParameter(2)*1e3, fit[i]->GetParameter(2)*secondPeakG/abs(fit[i]->GetParameter(1))*1000));
	if(dofit>=2) {
	  tex->DrawLatexNDC(0.13,0.68, Form( "#sigma_{0.5} = %.0f nC  --> %.0f keV",fit[i]->GetParameter(3)*1e3, fit[i]->GetParameter(3)*secondPeakG/abs(fit[i]->GetParameter(1))*1000));
	  tex->DrawLatexNDC(0.53,0.82, Form( "#mu_{B1.3} = %.0f nC  --> %.0f keV",   fit[i]->GetParameter(5)*1e3, fit[i]->GetParameter(5)*secondPeakG/fit[i]->GetParameter(1)*1000));
	  tex->DrawLatexNDC(0.13,0.61, Form( "Na22 = %.0f%%",fit[i]->GetParameter(4)*100));
	}
	tex->DrawLatexNDC(0.13,0.18, Form( "e = %.0f%%",fit[i]->GetParameter(0)*100/scale));
	//tex->DrawLatexNDC(0.13,0.28, Form( "#chi^{2}/NDF = %.2f",chi2));
      }
    }
  }
  for(int i=0; i!=CHANNELS; ++i)
    cout << " mu_" << i << " " << fit[i]->GetParameter(1) << endl;

  main2->SaveAs( Form("%s_CHARGE.pdf",inputfile.Data()), "PDF" );

  TCanvas *main3 = new TCanvas();
  main3->Divide(2,2);
  for(int i=0; i!=CHANNELS; ++i) {
    main3->cd(i+1);//->SetLogy(1);
    amplitude[i]->GetYaxis()->SetRangeUser(0,amplitude[i]->GetMaximum()*1.5);
    amplitude[i]->Draw("E");
    if(dofit==1) {
      double scale = amplitude[i]->GetEntries()* amplitude[i]->GetBinWidth(1);
      double xsc = xintmax[i]/300;
      int binM = amplitude[i]->GetMaximumBin();
      //double xm = amplitude[i]->GetBinCenter( binM-8 );
      //double xc = amplitude[i]->GetBinCenter( binM );
      //double xM = amplitude[i]->GetBinCenter( binM+8 );
      double xm = xsc*0.9;
      double xc = xsc*0.5;
      double xM = xsc*0.4;
      fitA[i]->SetParameter(0,0.5*scale);     fitA[i]->SetParLimits(0,scale*0.2,1.2*scale);          //eN
      fitA[i]->SetParameter(1,xc);            fitA[i]->SetParLimits(1,xm,xM);            //muG2
      fitA[i]->SetParameter(2,0.05*abs(xsc)); fitA[i]->SetParLimits(2,0.03*abs(xsc),0.09*abs(xsc));  //sigmaG2
      amplitude[i]->Fit( fitA[i], "RQ", "", xm, xM );
      for(int jjj=0; jjj!=3; ++jjj)
	amplitude[i]->Fit( fitA[i], "RQ", "",
			   fitA[i]->GetParameter(1)-fitA[i]->GetParameter(2),
			   fitA[i]->GetParameter(1)+fitA[i]->GetParameter(2) );
      amplitude[i]->Fit( fitA[i], "R", "",
			 fitA[i]->GetParameter(1)-fitA[i]->GetParameter(2),
			 fitA[i]->GetParameter(1)+fitA[i]->GetParameter(2) );
      cout << "********" << endl;
    }
  }
  
  TCanvas *main5 = new TCanvas();
  main5->Divide(2,2);
  TCanvas *main6 = new TCanvas();
  main6->Divide(2,2);
  for(int i=0; i!=CHANNELS; ++i) {
    main5->cd(i+1);
    amplitudeControl[i]->Draw();
    tex->DrawLatexNDC(0.12,0.8, Form("%f",amplitudeControl[i]->GetMean()) );
    main6->cd(i+1);
    integralControl[i]->Draw();
    tex->DrawLatexNDC(0.12,0.8, Form("%f",integralControl[i]->GetMean()) );
  }

  if(CHANNELS!=4) return 0;
  TCanvas *main4 = new TCanvas();
  main4->Divide(4,4);
  main4->cd(2)->SetLogz(1);
  corrHis[0]->Draw("colz");
  main4->cd(3)->SetLogz(1);
  corrHis[1]->Draw("colz");
  main4->cd(4)->SetLogz(1);
  corrHis[2]->Draw("colz");
  main4->cd(7)->SetLogz(1);
  corrHis[3]->Draw("colz");
  main4->cd(8)->SetLogz(1);
  corrHis[4]->Draw("colz");
  main4->cd(12)->SetLogz(1);
  corrHis[5]->Draw("colz");
  
  return 0;

}
