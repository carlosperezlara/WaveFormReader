char   charsize1;   // byte
ushort ushortsize2; // int16 word
short  shortsize2;  // int16 word
uint   uintsize4;   // int32
int    intsize4;    // int32
float  floatsize4;  // float
double doublesize8; // double
struct timestamp_t {
  double second;
  char  minute;
  char  hour;
  char  day;
  char  month;
  ushort  year;
} timestampsize14;

int readLecroy() {
  //ifstream fin("../aaa/C3--40.31V-GreenLaser-Mid.trc");
  //ifstream fin("../aaa/C2--40.31V-GreenLaser-Mid.trc");
  //ifstream fin("../Dec16/C4--GenPulse_1--00000.trc");
  //ifstream fin("../Dec16/C4--GenPulse_1--00001.trc");
  //ifstream fin("../Dec16/C4--GenPulse_1--00002.trc");
  //ifstream fin("C3--2e14-40.65V-n30-Green--00001.trc");
  ifstream fin("C2--2e14-40.65V-n30-Green--00001.trc");
  char  *c50Tmp = new char [50];
  fin.read(c50Tmp, 50);
  std::cout << c50Tmp << endl;
  int start = 11;
  char  *c16Tmp = new char [16];
  char  *c12Tmp = new char [12];

  //std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+16); //27
  //std::cout << "pos: " << fin.tellg() << endl;

  fin.read(c16Tmp, 16); // template name
  std::cout << "template name " << c16Tmp << endl;
  fin.read((char*) &ushortsize2, 2); // 8 or 16 bit
  int size = ushortsize2; // 0:byte 1:word(2bytes)
  std::cout << "size {8,16}" << ushortsize2 << endl;
  fin.read((char*) &ushortsize2, 2); // endianness
  std::cout << "endianness {>,<} "<< ushortsize2 << endl;

  fin.read((char*) &uintsize4, 4); //
  std::cout << "wave descriptor " << uintsize4 << endl;
  int skip = uintsize4;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "usertext " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "resdesc1 " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); // 
  std::cout << "trigtimearray " << uintsize4 << endl;
  skip += uintsize4;
  fin.read((char*) &uintsize4, 4); // 
  std::cout << "ristimearray " << uintsize4 << endl;
  skip += uintsize4;
  fin.read((char*) &uintsize4, 4); // 
  std::cout << "resarray1 " << uintsize4 << endl;

  std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+60); //12
  std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &uintsize4, 4);
  int samplesize = uintsize4;
  std::cout << "wavearray01 " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4);
  std::cout << "wavearray02 " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4);
  std::cout << "resarray2 " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4);
  std::cout << "resarray3 " << uintsize4 << endl;

  fin.read(c16Tmp, 16); // instrument name
  std::cout << "instrument name " << c16Tmp << endl;
  fin.read((char*) &uintsize4, 4); // instrument number
  std::cout << "instrument number " << uintsize4 << endl;

  fin.read(c16Tmp, 16); //
  std::cout << "trace label " << c16Tmp << endl;

  fin.read((char*) &uintsize4, 4); //
  std::cout << "4bytes reserved " << uintsize4 << endl;

  std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+116);//20
  std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &uintsize4, 4); // wave array count
  std::cout << "wavearraycount " << uintsize4 << endl;
  int wavearraycount = uintsize4;

  fin.read((char*) &uintsize4, 4); //
  std::cout << "pntsPerScreen " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "firstvalidpnt " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "lastvalidpnt " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "firstpnt " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "sparsingfactor " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "segmentindex " << uintsize4 << endl;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "subarraycount " << uintsize4 << endl;
  int waves = uintsize4;
  fin.read((char*) &uintsize4, 4); //
  std::cout << "sweepsperacq " << uintsize4 << endl;
  fin.read((char*) &ushortsize2, 2); //
  std::cout << "pointsperpair " << ushortsize2 << endl;
  fin.read((char*) &ushortsize2, 2); //
  std::cout << "pairoffset " << ushortsize2 << endl;

  std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+156);//36
  std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &floatsize4, 4); //
  float gain = floatsize4;
  std::cout << "vertical gain " << floatsize4 << endl;
  fin.read((char*) &floatsize4, 4); //
  float offset = floatsize4;
  std::cout << "vertical offset " << floatsize4 << endl;
  fin.read((char*) &floatsize4, 4); //
  double maxV = floatsize4*gain;
  std::cout << "maxvalue " << floatsize4 << endl;
  fin.read((char*) &floatsize4, 4); //
  double minV = floatsize4*gain;
  std::cout << "minvalue " << floatsize4 << endl;

  std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+172);//8
  std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &ushortsize2, 2); // nominalBits
  std::cout << "nominalBits " << ushortsize2 << endl;
  fin.read((char*) &ushortsize2, 2); //
  std::cout << "nomsubarraycount " << ushortsize2 << endl;

  std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+176);//2
  std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &floatsize4, 4); // horizInterval
  float horizinterval = floatsize4;
  std::cout << "horizinterval " << floatsize4 << endl;
  fin.read((char*) &doublesize8, 8); // horizOffset
  float horizoffset = floatsize4;
  std::cout << "horizoffset " << doublesize8 << endl;
  
  //std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+316);//128
  //std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &ushortsize2, 2); // 
  std::cout << "recordtypelist  {single_sweep,interleaved,histogram,graph,filter_coefficient,coplex,extrema,sequence_obsolote} " << ushortsize2 << endl;
  fin.read((char*) &ushortsize2, 2); // 
  std::cout << "processingdone  {off on} " << ushortsize2 << endl;

  //std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+326);//6
  //std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &ushortsize2, 2); // 
  std::cout << "vertical coupling {dc50, gnd, dc1m, gnd, ac1m} " << ushortsize2 << endl;
  
  //std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+344);//16
  //std::cout << "pos: " << fin.tellg() << endl;

  fin.read((char*) &ushortsize2, 2); // source {ch1,ch2,ch3,ch4,unk}
  std::cout << "source {ch1, ch2, ch3, ch4, unk} " << ushortsize2 << endl;

  //std::cout << "pos: " << fin.tellg() << endl;
  fin.seekg(start+skip);//0
  //std::cout << "pos: " << fin.tellg() << endl;

  int samples = wavearraycount/waves;
  cout << "samples: " << samples << " || waves: " << waves << endl;
  TH1D *wave = new TH1D("wave","wave;s;V",samples,horizoffset,samples*horizinterval + horizoffset);
  TH2D *allwaves = new TH2D("allwaves","wave;s;V",samples,horizoffset,samples*horizinterval + horizoffset, 100, minV, maxV);
  for(int xx=0; xx!=waves; ++xx) {
    for(int i=0; i!=samples; ++i) {
      float mV;
      if( size==1 ) { //16bit packet
	fin.read((char*) &shortsize2, 2); // pack
	mV = shortsize2*gain + offset;
      } else { // 8bit packet
	fin.read((char*) &charsize1, 1); // pack
	mV = charsize1*gain + offset;
      }
      //std::cout << i << " "<< shortsize2 << " " << mV << endl;
      wave->SetBinContent(i+1,mV);
      float x = wave->GetBinCenter(i+1);
      allwaves->Fill( x, mV );
    }
    //if(xx==3) return 1;
  }
  new TCanvas();
  wave->Draw();
  new TCanvas();
  allwaves->Draw("colz");
  
  fin.close();
  return 0;
}
