// C++ includes
#include <fstream>
#include <string>
#include <iostream>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

//LOCAL INCLUDES
#include "Aux.hh"
//#include "Config.hh"

using namespace std;

struct FTBFPixelEvent {
    double xSlope;
    double ySlope;
    double xIntercept;
    double yIntercept;
    double chi2;
    double xResidBack;
    double yResidBack;
    int trigger;
    int runNumber;
    int nPlanes;
    int numPixels;
    int numBackPlanes;  
    Long64_t timestamp;
    Long64_t bco;    
};


TStyle* style;



int main(int argc, char **argv) {
  gROOT->SetBatch();
  
  FILE* fp1;
  char stitle[200];
  int dummy;

  //**************************************
  // Parse command line arguments
  //**************************************

  int numRequiredArgs = 4;
  if (argc - 1 < numRequiredArgs) {
    std::cerr <<  "[ERROR]: Usage --> ./maketree --inputFileName=<your_caen_drs4_file.dat> --pixelInputFileName=<your_pixel_file.root> --outputFileName=<your_output_file.root> --nEvents=<number_of_events_to_be_processed>" << std::endl;
      return -1;
  }
  
  std::cout << "\n=== Beginning program ===\n" << std::endl;

  //****************************************
  // Getting Input FileName (CAEN DRS4 file)
  //****************************************
  std::string inputFileName = ParseCommandLine( argc, argv, "--inputFileName=" );
  if ( inputFileName == "" )
    {
      std::cerr << "[ERROR]: please provide a valid inputFileName. Use: --inputFileName=<your_input_file_name> " << std::endl;
      exit(0);
    }
  
  //*****************************
  // Getting Pixel Input FileName
  //*****************************
  std::string pixelInputFileName = ParseCommandLine( argc, argv, "--pixelInputFileName=" );
  if ( pixelInputFileName == "" )
  {
      std::cerr << "[ERROR]: please provide a valid pixelInputFileName. Use: --pixelInputFileName=<your_pixel_file_name> " << std::endl;
      exit(0);
  }
  
  //*******************************
  // Getting Output FileName (ROOT)
  //*******************************
  std::string outputFileName = ParseCommandLine( argc, argv, "--outputFileName=" );
  if ( outputFileName == "" )
  {
      std::cerr << "[ERROR]: please provide a valid outputFileName. Use: --outputFileName=<your_output_file_name> " << std::endl;
      exit(0);
  }

  //*****************************************
  // Getting number of events to be processed
  //*****************************************
  std::string nEvents = ParseCommandLine( argc, argv, "--nEvents=" );
  if ( nEvents == "" )
    {
      std::cerr << "[ERROR]: please provide a valid nEvents. Use: --nEvents=<events_to_process> " << std::endl;
      exit(0);
    }
  int nevents = atoi(nEvents.c_str());
  
  std::cout << "[INFO]: input file --> " << inputFileName << std::endl;
  std::cout << "[INFO]: output file --> " << outputFileName << std::endl;
  std::cout << "[INFO]: processing --> " << nEvents << " events" << std::endl;
  
  //****************************
  // Getting channel config file
  //****************************
  std::string configName = ParseCommandLine( argc, argv, "--config=" );
  if ( configName == "" )
    {
      configName = "config/15may2017.config";
      std::cout << "[INFO]: using default config file: " << configName << ". Use: --config=<config_file> to set your config" << std::endl;
    }

  //*********************************************************
  // Check if has valid input file, otherwise exit with error
  //*********************************************************
  ifstream ifile(inputFileName);
  if (!ifile)
    {
      std::cerr << "[ERROR]: !USAGE! Input file does not exist. Please enter valid file name" << std::endl;
      exit(0);
    }
  

  
  //*********************************************************
  // Check if has valid input file, otherwise exit with error
  //*********************************************************
  bool saveRaw = false;
  std::string _saveRaw = ParseCommandLine( argc, argv, "--saveRaw" );
  if ( _saveRaw == "yes" )
    {
      saveRaw = true;
      std::cout << "[INFO]: Will save raw pulses" << std::endl;
    }

  bool drawDebugPulses = false;
  std::string _drawDebugPulses = ParseCommandLine( argc, argv, "--debug" );
  if ( _drawDebugPulses == "yes" )
    {
      drawDebugPulses = true;
      std::cout << "[INFO]: Will draw pulse to check procedure" << std::endl;
    }

    /*
  std::cout << "\n=== Parsing configuration file " << configName << " ===\n" << std::endl;
  Config config(configName);
  if ( !config.hasChannels() || !config.isValid() ) {
    std::cerr << "\nFailed to load channel information from config " << configName << std::endl;
    return -1;
  }
    */
  
  //**************************************
  // Load Voltage Calibration
  //**************************************

  std::cout << "\n=== Loading voltage calibration ===\n" << std::endl;
  double off_mean[2][9][1024];
  for( int i = 0; i < 2; i++ ){
      sprintf( stitle, "dt5742_offset_%d.txt", i );
      fp1 = fopen( stitle, "r" );
      printf("Loading offset data from %s\n", stitle);

      for( int k = 0; k < 1024; k++ ) {     
          for( int j = 0; j < 9; j++ ){      
              dummy = fscanf( fp1, "%lf ", &off_mean[i][j][k] ); 
          }
      }
      fclose(fp1);
  }

  //**************************************
  // Load Time Calibration
  //**************************************

  std::cout << "\n=== Loading time calibration ===\n" << std::endl;
  double fdummy;
  double tcal_dV[2][1024];
  for( int i = 0; i < 2; i++ ) {
      sprintf( stitle, "dt5742_axis_%d.txt", i );
      fp1 = fopen( stitle, "r" );
      printf("Loading dV data from %s\n", stitle);

      for( int k = 0; k < 1024; k++)      
          dummy = fscanf( fp1, "%lf", &tcal_dV[i][k] ); 
      fclose(fp1);
  }
  double dV_sum[2] = {0, 0};
  for( int i = 0; i < 2; i++ ) {
      for( int j = 0; j < 1024; j++ )
          dV_sum[i] += tcal_dV[i][j];
  }

  double tcal[2][1024];
  for( int i = 0; i < 2; i++) {
      for( int j = 0; j < 1024; j++) {
          tcal[i][j] = tcal_dV[i][j] / dV_sum[i] * 200.0;
      }
  }
  
  //**************************************
  // Define output
  //**************************************

  TFile* file = new TFile( outputFileName.c_str(), "RECREATE", "CAEN V1742");
  TTree* tree = new TTree("pulse", "Digitized waveforms");

  int event;
  short tc[2]; // trigger counter bin
  float times[2][1024]; // calibrated time
  short raw[18][1024]; // ADC counts
  float channel[18][1024]; // calibrated input (in V)
  float xIntercept;
  float yIntercept;
  float xSlope;
  float ySlope;
  float chi2;
  int ntracks;
 
  tree->Branch("event", &event, "event/I");
  tree->Branch("tc", tc, "tc[2]/s");
  if (saveRaw) {
    tree->Branch("raw", raw, "raw[18][1024]/S");   
  }
  tree->Branch("channel", channel, "channel[18][1024]/F");
  tree->Branch("times", times, "times[2][1024]/F");
  tree->Branch("xIntercept", &xIntercept, "xIntercept/F");
  tree->Branch("yIntercept", &yIntercept, "yIntercept/F");
  tree->Branch("xSlope", &xSlope, "xSlope/F");
  tree->Branch("ySlope", &ySlope, "ySlope/F");
  tree->Branch("chi2", &chi2, "chi2/F");
  tree->Branch("ntracks", &ntracks, "ntracks/I");

  //*************************
  // Open Pixel Tree
  //*************************
  TFile *pixelDataFile = TFile::Open(pixelInputFileName.c_str(),"READ");
  if (!pixelDataFile)
    {
      std::cout << "[ERROR]: Pixel file not found" << std::endl;
      return -1;
    }
  TTree *pixelTree = (TTree*)pixelDataFile->Get("CMSTiming");
  if (!pixelTree) {
    cout << "Error: Pixel Tree not found\n";
    return 0;
  }
 
  FTBFPixelEvent pixelEvent;
  pixelTree->SetBranchAddress("event",&pixelEvent);


  // temp variables for data input
  uint   event_header;
  uint   temp[3];
  ushort samples[9][1024];

  //*************************
  // Open Input File
  //*************************

  FILE* fpin = fopen( inputFileName.c_str(), "r" );

  //*************************
  //Event Loop
  //*************************

  std::cout << "\n=== Processing input data ===\n" << std::endl;
  int nGoodEvents = 0;
  int maxEvents = nevents;
  if (nevents < 0) maxEvents = 999999;

  int pixLastEvent=0;

  for( int iEvent = 0; iEvent < maxEvents; iEvent++){ 
      
      if ( iEvent % 100 == 0 ) {
          if (nevents >= 0) {
              std::cout << "Event " << iEvent << " of " << nevents << std::endl;
          } else {
              std::cout << "Event " << iEvent << "\n";
          }
      }
      event = nGoodEvents; // for output tree
      
      // first header word
      dummy = fread( &event_header, sizeof(uint), 1, fpin);

      // second header word
      dummy = fread( &event_header, sizeof(uint), 1, fpin);  
      uint grM     = event_header & 0b11; // 2-bit channel group mask
      // std::cout << " Group mask = " << grM << std::endl;

      // third and fourth header words
      dummy = fread( &event_header, sizeof(uint), 1, fpin);  
      dummy = fread( &event_header, sizeof(uint), 1, fpin);  
      
      // check for end of file
      if (feof(fpin)) break;
      
      //*************************
      // Parse group mask into channels
      //*************************
      
      bool _isGR_On[2];
      _isGR_On[0] = (grM & 0x01);
      _isGR_On[1] = (grM & 0x02);
      
      int activeGroupsN = 0;
      int realGroup[2] = {-1, -1 };
      for ( int l = 0; l < 2; l++ ) {
          if ( _isGR_On[l] ) 
          {
              realGroup[activeGroupsN] = l; 
              activeGroupsN++;
          }
      }
      
      //************************************
      // Loop over channel groups
      //************************************
      
      for ( int group = 0; group < activeGroupsN; group++ ) {
          // Read group header
          dummy = fread( &event_header, sizeof(uint), 1, fpin);  
          // ushort tcn = (event_header >> 20) & 0xfff; // trigger counter bin
          ushort tcn = (event_header << 2) >> 22; // trigger counter bin
          tc[realGroup[group]] = tcn;
          uint TR = (event_header >> 12) & 0b1;
          // std::cout << " Group =  " << group << "   TR0 = " << TR << std::endl;
          
          // Check if all channels were active (if 8 channels active return 3072)
          int nsample = (event_header & 0xfff) / 3;
          // std::cout << " Group =  " << group << "   samples = " << nsample << std::endl;
          
          // Define time coordinate
          times[realGroup[group]][0] = 0.0;
          for( int i = 1; i < 1024; i++ ){
          //    times[realGroup[group]][i] = float(i);
              times[realGroup[group]][i] = float(tcal[realGroup[group]][(i-1+tcn)%1024] + times[realGroup[group]][i-1]);
          }      
          
          //************************************
          // Read sample info for group
          //************************************      
          
          for ( int i = 0; i < nsample; i++ ) {
              dummy = fread( &temp, sizeof(uint), 3, fpin );  
              samples[0][i] =  temp[0] & 0xfff;
              samples[1][i] = (temp[0] >> 12) & 0xfff;
              samples[2][i] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
              samples[3][i] = (temp[1] >>  4) & 0xfff;
              samples[4][i] = (temp[1] >> 16) & 0xfff;
              samples[5][i] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
              samples[6][i] = (temp[2] >>  8) & 0xfff;
              samples[7][i] =  temp[2] >> 20;	
          }
          
          // Trigger channel
          if(TR){
              for(int j = 0; j < nsample/8; j++){
                  fread( &temp, sizeof(uint), 3, fpin);  
                  samples[8][j*8+0] =  temp[0] & 0xfff;
                  samples[8][j*8+1] = (temp[0] >> 12) & 0xfff;
                  samples[8][j*8+2] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
                  samples[8][j*8+3] = (temp[1] >>  4) & 0xfff;
                  samples[8][j*8+4] = (temp[1] >> 16) & 0xfff;
                  samples[8][j*8+5] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
                  samples[8][j*8+6] = (temp[2] >>  8) & 0xfff;
                  samples[8][j*8+7] =  temp[2] >> 20;
              }
          }else{
              for(int j = 0; j < nsample/8; j++){
                  samples[8][j*8+0] =  0;
                  samples[8][j*8+1] =  0;
                  samples[8][j*8+2] =  0;
                  samples[8][j*8+3] =  0;
                  samples[8][j*8+4] =  0;
                  samples[8][j*8+5] =  0;
                  samples[8][j*8+6] =  0;
                  samples[8][j*8+7] =  0;
              }
          }
          
          //************************************
          // Loop over channels 0-8
          //************************************      
          
          for(int i = 0; i < 9; i++) {
              
              int totalIndex = realGroup[group]*9 + i;
              
              for ( int j = 0; j < 1024; j++ ) {
                  raw[totalIndex][j] = 0;
                  channel[totalIndex][j] = 0;
              }
              
              // Fill pulses
              for ( int j = 0; j < 1024; j++ ) {
                  raw[totalIndex][j] = (short)(samples[i][j]);
                  channel[totalIndex][j] = (double)(samples[i][j]) - (double)(off_mean[realGroup[group]][i][(j+tcn)%1024]);
                  channel[totalIndex][j] = 1000.*((channel[totalIndex][j] - 0)/4095. - 0.5) + 0;
              }
          }
          dummy = fread( &event_header, sizeof(uint), 1, fpin);
      }

      //find corresponding pixel event    
      xIntercept = -1e+9;
      yIntercept = -1e+9;
      xSlope = -999;
      ySlope = -999;
      chi2 = -999.;
      ntracks = 0;
      
      for( int iPixelEvent = pixLastEvent; iPixelEvent < pixelTree->GetEntries(); iPixelEvent++){ 
          pixelTree->GetEntry(iPixelEvent);
          if (pixelEvent.trigger == iEvent) {
              if(pixelEvent.chi2>1e-6){
                  xIntercept = pixelEvent.xIntercept;
                  yIntercept = pixelEvent.yIntercept;
                  xSlope = pixelEvent.xSlope;
                  ySlope = pixelEvent.ySlope;
                  chi2 = pixelEvent.chi2;
                  ntracks++;
              }
              pixLastEvent = iPixelEvent;
          }else if(pixelEvent.trigger > iEvent){
              break;
          }
      }        
      
      tree->Fill();
      nGoodEvents++;
  }
  
  fclose(fpin);
  cout << "\nProcessed total of " << nGoodEvents << " events\n";
  
  file->Write();
  file->Close();
  
  return 0;
}



