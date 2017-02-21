#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <climits>
#include <TMath.h>
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalDeterministicFit.h"

constexpr float HcalDeterministicFit::invGpar[3];
constexpr float HcalDeterministicFit::negThresh[2];
constexpr int HcalDeterministicFit::HcalRegion[2];
constexpr float HcalDeterministicFit::rCorr[2];
constexpr float HcalDeterministicFit::rCorrSiPM[2];

using namespace std;

HcalDeterministicFit::HcalDeterministicFit() {
}

HcalDeterministicFit::~HcalDeterministicFit() { 
}

void HcalDeterministicFit::init(HcalTimeSlew::ParaSource tsParam, HcalTimeSlew::BiasSetting bias, bool iApplyTimeSlew, PedestalSub pedSubFxn_, std::vector<double> pars, double respCorr) {
  for(int fi=0; fi<9; fi++){
	fpars[fi] = pars.at(fi);
  }

  applyTimeSlew_=iApplyTimeSlew;
  useExtPulse_=false;
  fTimeSlew=tsParam;
  fTimeSlewBias=bias;
  fPedestalSubFxn_=pedSubFxn_;
  frespCorr=respCorr;
  shape_ = 0;
  std::fill(&pulseFrac1_[0][0], &pulseFrac1_[0][0] + 58*10, 0.);
  std::fill(&pulseFrac2_[0][0], &pulseFrac2_[0][0] + 58*10, 0.);
  std::fill(&pulseFrac_[0][0], &pulseFrac_[0][0] + 58*10, 0.);
  std::fill(&pulseFracDeriv1_[0][0], &pulseFracDeriv1_[0][0] + 58*10, 0.);
  std::fill(&pulseFracDeriv2_[0][0], &pulseFracDeriv2_[0][0] + 58*10, 0.);
  std::fill(&pulseFracDeriv_[0][0], &pulseFracDeriv_[0][0] + 58*10, 0.);
  std::fill(minCharge1_, minCharge1_ + 58, 0.);
  std::fill(minCharge2_, minCharge2_ + 58, 0.);
  std::fill(minCharge_, minCharge_ + 58, 0.);
  std::fill(maxCharge1_, maxCharge1_ + 58, 0.);
  std::fill(maxCharge2_, maxCharge2_ + 58, 0.);
  std::fill(maxCharge_, maxCharge_ + 58, 0.);
  std::fill(timeSlew1_, timeSlew1_ + 58, 0.);
  std::fill(timeSlew2_, timeSlew2_ + 58, 0.);
  std::fill(timeSlew_, timeSlew_ + 58, 0.);

  // ARRAY FOR LANDAU
  // Combine [0,10]+[10,20] into [0,20] and [580,590]+[590,600] into [580,600] because there's only 58 rows in pulseFrac
  minCharge1_[0] = 0.;
  maxCharge1_[0] = 20.;
  float tstart0 = TMath::Min(6.0, 12.2999-2.19142*log(5));
  for (int j = 0; j < 4; j++) pulseFrac1_[0][j] = 0;
  pulseFrac1_[0][4] = landauFrac[ int( ceil( -tstart0 + 25 ) ) ];
  pulseFrac1_[0][5] = landauFrac[ int( ceil( -tstart0 + 50 ) ) ];
  pulseFrac1_[0][6] = landauFrac[ int( ceil( -tstart0 + 75 ) ) ];
  for (int j = 7; j < 10; j++) pulseFrac1_[0][j] = 0;
  for (int j = 0; j < 10; j++) pulseFracDeriv1_[0][j] = 0;
  if (pulseFrac1_[0][4] != 0) 
      timeSlew1_[0] = pulseFrac1_[0][5]/pulseFrac1_[0][4];
  else timeSlew1_[0] = 0;

  minCharge1_[57] = 580.;
  maxCharge1_[57] = 600.;
  float tstart1 = TMath::Min(6.0, 12.2999-2.19142*log(590));
  for (int j = 0; j < 4; j++) pulseFrac1_[57][j] = 0;
  pulseFrac1_[57][4] = landauFrac[ int( ceil( -tstart1 + 25 ) ) ];
  pulseFrac1_[57][5] = landauFrac[ int( ceil( -tstart1 + 50 ) ) ];
  pulseFrac1_[57][6] = landauFrac[ int( ceil( -tstart1 + 75 ) ) ];
  for (int j = 7; j < 10; j++) pulseFrac1_[57][j] = 0;
  for (int j = 0; j < 10; j++) pulseFracDeriv1_[57][j] = 0;
  if (pulseFrac1_[57][4] != 0) 
      timeSlew1_[57] = pulseFrac1_[57][5]/pulseFrac1_[57][4];
  else timeSlew1_[57] = 0;

  int kay = 1; 
  for (int i=25; i<580; i+=10) 
  {
    float tstart=TMath::Min(6.0, 12.2999-2.19142*log(i));
    minCharge1_[kay] = i-5.;
    maxCharge1_[kay] = i+5.;
    for (int j = 0; j < 4; j++) pulseFrac1_[kay][j] = 0;
    pulseFrac1_[kay][4] = landauFrac[ int( ceil( -tstart + 25 ) ) ];
    pulseFrac1_[kay][5] = landauFrac[ int( ceil( -tstart + 50 ) ) ];
    pulseFrac1_[kay][6] = landauFrac[ int( ceil( -tstart + 75 ) ) ];
    for (int j = 7; j < 10; j++) pulseFrac1_[kay][j] = 0;
    for (int j = 0; j < 10; j++) pulseFracDeriv1_[kay][j] = 0;
    if (pulseFrac1_[kay][4] != 0) 
        timeSlew1_[kay] = pulseFrac1_[kay][5]/pulseFrac1_[kay][4];
    else timeSlew1_[kay] = 0;
    kay++;
  }

  // ARRAY FOR 105
  ts1_ = 8;
  ts2_ = 19;
  ts3_ = 29.3;
  thpd_ = 4.0;
  tpre_ = 9.0;
  wd1_ = 2.0;
  wd2_ = 0.7;
  wd3_ = 0.32;
  makeShape105();
  for (unsigned int i = 0; i < 58; i++)
  {
      minCharge2_[i] = 20.+10.*i;
      maxCharge2_[i] = 30.+10.*i;
      
      double rawDelay=tzero[1]+slope[1]*log(25+10*i);
      float tstart = (rawDelay<0)?(0):((rawDelay>tmax[1])?(tmax[1]):(rawDelay));
       
      std::array<double,HcalConst::maxSamples> tmpPulseFrac, tmpIni, tmpFin;
      compute105(tmpPulseFrac, tstart);
      compute105(tmpIni, tstart-1);
      compute105(tmpFin, tstart+1);
      for (int j = 0; j<10; j++)
      {
          pulseFrac2_[i][j] = tmpPulseFrac[j];
          pulseFracDeriv2_[i][j] = 0.5*(tmpFin[j] - tmpIni[j]);
      }
      if (pulseFrac2_[i][4] != 0) timeSlew2_[i] = pulseFrac2_[i][5]/pulseFrac2_[i][4];
      else timeSlew2_[i] = 0;
  }
  
  

}

void HcalDeterministicFit::setExternalPulseShape(int shape) {
    std::cout << "deterministic " << shape << std::endl;
  shape_ = shape;
  if (shape > 0) useExtPulse_=true;
  if (shape == 0) useExtPulse_=false; // just to be sure
  if (shape == 1) 
  {
     std::copy(&pulseFrac1_[0][0], &pulseFrac1_[0][0]+58*10, &pulseFrac_[0][0]);
     std::copy(&pulseFracDeriv1_[0][0], &pulseFracDeriv1_[0][0]+58*10, &pulseFracDeriv_[0][0]);
     std::copy(minCharge1_, minCharge1_+58, minCharge_);
     std::copy(maxCharge1_, maxCharge1_+58, maxCharge_);
     std::copy(timeSlew1_, timeSlew1_+58, timeSlew_);
  }
  else if (shape == 2) 
  {
     std::copy(&pulseFrac2_[0][0], &pulseFrac2_[0][0]+58*10, &pulseFrac_[0][0]);
     std::copy(&pulseFracDeriv2_[0][0], &pulseFracDeriv2_[0][0]+58*10, &pulseFracDeriv_[0][0]);
     std::copy(minCharge2_, minCharge2_+58, minCharge_);
     std::copy(maxCharge2_, maxCharge2_+58, maxCharge_);
     std::copy(timeSlew2_, timeSlew2_+58, timeSlew_);
  }
}



constexpr float HcalDeterministicFit::landauFrac[];
// Landau function integrated in 1 ns intervals
//Landau pulse shape from https://indico.cern.ch/event/345283/contribution/3/material/slides/0.pdf
//Landau turn on by default at left edge of time slice 
// normalized to 1 on [0,10000]
void HcalDeterministicFit::getLandauFrac(float tStart, float tEnd, float &sum) const{
  if (std::abs(tStart-tEnd-tsWidth)<0.1) {
    sum=0;
    return;
  }
  sum= landauFrac[int(ceil(tStart+tsWidth))];
  return;
}


void HcalDeterministicFit::getLandauFrac(float fC, int offset, double fpar0, double fpar1, double fpar2, float &sum) const{

  if (useExtPulse_) 
  {
    // hardcoded array :(
    int chargeBin = -1;
    for (int i=0; i<58; i++) {
      if (fC>minCharge_[i] && fC<maxCharge_[i]) chargeBin=i;
    }
    if (fC>maxCharge_[57]) chargeBin=57;

    if (chargeBin==-1) chargeBin=0;

    sum=pulseFrac_[chargeBin][offset+3];
    return;

  }
}


void HcalDeterministicFit::phase1Apply(const HBHEChannelInfo& channelData,
				       float& reconstructedEnergy,
				       float& reconstructedTime,
                       float& previousEnergy,
                       float& nextEnergy) const
{

  std::vector<double> corrCharge;
  std::vector<double> inputCharge;
  std::vector<double> inputPedestal;
  double gainCorr = 0;
  double respCorr = 0;

  for(unsigned int ip=0; ip<channelData.nSamples(); ip++){

    double charge = channelData.tsRawCharge(ip);
    double ped = channelData.tsPedestal(ip); 
    double gain = channelData.tsGain(ip);

    gainCorr = gain;
    inputCharge.push_back(charge);
    inputPedestal.push_back(ped);

  }

  fPedestalSubFxn_.calculate(inputCharge, inputPedestal, corrCharge);

  const HcalDetId& cell = channelData.id();

  double fpar0, fpar1, fpar2;
  if(std::abs(cell.ieta())<HcalRegion[0]){
    fpar0 = fpars[0];
    fpar1 = fpars[1];
    fpar2 = fpars[2];
  }else if(std::abs(cell.ieta())==HcalRegion[0]||std::abs(cell.ieta())==HcalRegion[1]){
    fpar0 = fpars[3];
    fpar1 = fpars[4];
    fpar2 = fpars[5];
  }else{
    fpar0 = fpars[6];
    fpar1 = fpars[7];
    fpar2 = fpars[8];
  }

  if (fTimeSlew==0) respCorr=1.0;
  else if (fTimeSlew==1) channelData.hasTimeInfo()?respCorr=rCorrSiPM[0]:respCorr=rCorr[0];
  else if (fTimeSlew==2) channelData.hasTimeInfo()?respCorr=rCorrSiPM[1]:respCorr=rCorr[1];
  else if (fTimeSlew==3)respCorr=frespCorr;

  float ch3=0;
  float ch4=0;
  float ch5=0;

  float tsShift3=0;
  float tsShift4=0;
  float tsShift5=0;

  if (useExtPulse_) 
  {
    float i3=0;
    getLandauFrac(inputCharge[3], 1, fpar0, fpar1, fpar2, i3);
    float n3=0;
    getLandauFrac(inputCharge[3], 2, fpar0, fpar1, fpar2, n3);
    float nn3=0;
    getLandauFrac(inputCharge[3], 3, fpar0, fpar1, fpar2, nn3);

    float i4=0;
    getLandauFrac(inputCharge[4], 1, fpar0, fpar1, fpar2, i4);
    float n4=0;
    getLandauFrac(inputCharge[4], 2, fpar0, fpar1, fpar2, n4);

    float i5=0;
    getLandauFrac(inputCharge[5], 1, fpar0, fpar1, fpar2, i5);
    float n5=0;
    getLandauFrac(inputCharge[5], 2, fpar0, fpar1, fpar2, n5);

    if (i3 != 0 && i4 != 0 && i5 != 0) 
    {
      ch3=corrCharge[3]/i3;
      ch4=(i3*corrCharge[4]-n3*corrCharge[3])/(i3*i4);
      ch5=(n3*n4*corrCharge[3]-i4*nn3*corrCharge[3]-i3*n4*corrCharge[4]+i3*i4*corrCharge[5])/(i3*i4*i5);

      if (ch3<negThresh[0]) 
      {
            ch3=negThresh[0];
            ch4=corrCharge[4]/i4;
            ch5=(i4*corrCharge[5]-n4*corrCharge[4])/(i4*i5);
      }
      
      if (ch5<negThresh[0] && ch4>negThresh[1]) 
      {
            float newTS = (corrCharge[5]-negThresh[0]*i5)/(corrCharge[4]-ch3*i3);
            int newBin=0;
            for (int k=0; k<58; k++) 
            {
              if (newTS < timeSlew_[k]) newBin=k;
            }
            float i4_new = pulseFrac_[newBin][4];
            
            if (i4_new!=0)
              {
                ch5=negThresh[0];
                ch4=(corrCharge[4]-ch3*n3)/(i4_new);
              }
      }
    }

    if (ch4<1) 
    {
      ch4=0;
    }

  } 
  else if (shape_ == 0)  // just to make sure
  {
    //default Run2
    if (applyTimeSlew_) 
    {
      tsShift3=HcalTimeSlew::delay(inputCharge[3], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());
      tsShift4=HcalTimeSlew::delay(inputCharge[4], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());
      tsShift5=HcalTimeSlew::delay(inputCharge[5], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());
    }

    float i3=0;
    getLandauFrac(-tsShift3,-tsShift3+tsWidth,i3);
    float n3=0;
    getLandauFrac(-tsShift3+tsWidth,-tsShift3+tsWidth*2,n3);
    float nn3=0;
    getLandauFrac(-tsShift3+tsWidth*2,-tsShift3+tsWidth*3,nn3);

    float i4=0;
    getLandauFrac(-tsShift4,-tsShift4+tsWidth,i4);
    float n4=0;
    getLandauFrac(-tsShift4+tsWidth,-tsShift4+tsWidth*2,n4);

    float i5=0;
    getLandauFrac(-tsShift5,-tsShift5+tsWidth,i5);
    float n5=0;
    getLandauFrac(-tsShift5+tsWidth,-tsShift5+tsWidth*2,n5);

    if (i3 != 0 && i4 != 0 && i5 != 0) 
    {
      ch3=corrCharge[3]/i3;
      ch4=(i3*corrCharge[4]-n3*corrCharge[3])/(i3*i4);
      ch5=(n3*n4*corrCharge[3]-i4*nn3*corrCharge[3]-i3*n4*corrCharge[4]+i3*i4*corrCharge[5])/(i3*i4*i5);

      if (ch3<negThresh[0]) 
      {
            ch3=negThresh[0];
            ch4=corrCharge[4]/i4;
            ch5=(i4*corrCharge[5]-n4*corrCharge[4])/(i4*i5);
      }
      if (ch5<negThresh[0] && ch4>negThresh[1]) 
      {
        float newTS = (corrCharge[5]-negThresh[0]*i5)/(corrCharge[4]-ch3*i3);
        int newBin=0;
        for (int k=0; k<58; k++) 
        {
          if (newTS < timeSlew_[k]) newBin=k;
        }
        float i4_new = pulseFrac_[newBin][4];
        
        if (i4_new!=0)
          {
            ch5=negThresh[0];
            ch4=(corrCharge[4]-ch3*n3)/(i4_new);
          }
      }
    }

    if (ch4<1) ch4=0;

  } ////

  reconstructedEnergy=ch4*respCorr*gainCorr;
  previousEnergy=ch3*respCorr*gainCorr;
  nextEnergy=ch5*respCorr*gainCorr;
  reconstructedTime=tsShift4;
}

void HcalDeterministicFit::compute105(std::array<double,HcalConst::maxSamples> &ntmpbin, float pulseTime)  
{
    float slew=0;
    constexpr int ns_per_bx = HcalConst::nsPerBX;
    constexpr int num_ns = HcalConst::nsPerBX*HcalConst::maxSamples;
    constexpr int num_bx = num_ns/ns_per_bx;
    //Get the starting time
    int i_start         = ( -HcalConst::iniTimeShift - pulseTime - slew >0 ? 0 : (int)std::abs(-HcalConst::iniTimeShift-pulseTime-slew) + 1);
    double offset_start = i_start - HcalConst::iniTimeShift - pulseTime - slew; //-199-2*pars[0]-2.*slew (for pars[0] > 98.5) or just -98.5-pars[0]-slew;
    if( offset_start == 1.0 ){ offset_start = 0.; i_start-=1; } //Deal with boundary
    const int bin_start        = (int) offset_start; //bin off to integer
    const int bin_0_start      = ( offset_start < bin_start + 0.5 ? bin_start -1 : bin_start ); //Round it
    const int iTS_start        = i_start/ns_per_bx;         //Time Slice for time shift
    const int distTo25ns_start = HcalConst::nsPerBX - 1 - i_start%ns_per_bx;    //Delta ns 
    const double factor = offset_start - bin_0_start - 0.5; //Small correction?
    
    //Build the new pulse
    ntmpbin[iTS_start] = (bin_0_start == -1 ? // Initial bin (I'm assuming this is ok)
			  accVarLenIdxMinusOneVec[distTo25ns_start] + factor * diffVarItvlIdxMinusOneVec[distTo25ns_start]
			  : accVarLenIdxZEROVec    [distTo25ns_start] + factor * diffVarItvlIdxZEROVec    [distTo25ns_start]);
    //Fill the rest of the bins
    for(int iTS = iTS_start+1; iTS < num_bx; ++iTS){
	int bin_idx = distTo25ns_start + 1 + (iTS-iTS_start-1)*ns_per_bx + bin_0_start;
	ntmpbin[iTS] = acc25nsVec[bin_idx] + factor * diff25nsItvlVec[bin_idx];
    }
}

void HcalDeterministicFit::makeShape105()
{
    // pulse shape componnts over a range of time 0 ns to 255 ns in 1 ns steps
    unsigned int nbin = 256;

    for (unsigned int i=0; i<nbin; i++) {Shape_.push_back(0.0);}
    //Shape_(nbin,0.0);  // zeroing output pulse shape
    std::vector<float> nth(nbin,0.0);   // zeroing HPD drift shape
    std::vector<float> ntp(nbin,0.0);   // zeroing Binkley preamp shape
    std::vector<float> ntd(nbin,0.0);   // zeroing Scintillator decay shape

    unsigned int i,j,k;
    float norm;

    // HPD starts at I and rises to 2I in thpd of time
    norm=0.0;
    for(j=0;j<thpd_ && j<nbin;j++){
      nth[j] = 1.0 + ((float)j)/thpd_;
      norm += nth[j];
    }
    // normalize integrated current to 1.0
    for(j=0;j<thpd_ && j<nbin;j++){
      nth[j] /= norm;
    }
  
    // Binkley shape over 6 time constants
    norm=0.0;
    for(j=0;j<6*tpre_ && j<nbin;j++){
      ntp[j] = ((float)j)*exp(-((float)(j*j))/(tpre_*tpre_));
      norm += ntp[j];
    }
    // normalize pulse area to 1.0
    for(j=0;j<6*tpre_ && j<nbin;j++){
      ntp[j] /= norm;
    }

    // effective tile plus wave-length shifter decay time over 4 time constants
    unsigned int tmax_ = 6 * (int)ts3_;
 
    norm=0.0;
    for(j=0;j<tmax_ && j<nbin;j++){
      ntd[j] = wd1_ * exp(-((float)j)/ts1_) + 
	wd2_ * exp(-((float)j)/ts2_) + 
	wd3_ * exp(-((float)j)/ts3_) ; 
      norm += ntd[j];
    }
    // normalize pulse area to 1.0
    for(j=0;j<tmax_ && j<nbin;j++){
      ntd[j] /= norm;
    }
  
    unsigned int t1,t2,t3,t4;
    for(i=0;i<tmax_ && i<nbin;i++){
      t1 = i;
      //    t2 = t1 + top*rand;
      // ignoring jitter from optical path length
      t2 = t1;
      for(j=0;j<thpd_ && j<nbin;j++){
	t3 = t2 + j;
	for(k=0;k<4*tpre_ && k<nbin;k++){       // here "4" is set deliberately,
	  t4 = t3 + k;                         // as in test fortran toy MC ...
	  if(t4<nbin){                         
	    unsigned int ntb=t4;                        
	    Shape_[ntb] += ntd[i]*nth[j]*ntp[k];
	  }
	}
      }
    }
  
    // normalize for 1 GeV pulse height
    norm = 0.;
    for(i=0;i<nbin;i++){
      norm += Shape_[i];
    }

    //cout << " Convoluted SHAPE ==============  " << endl;
    for(i=0; i<nbin; i++){
      Shape_[i] /= norm;
      //std::cout << i << ",  " << ntmp[i] << std::endl;   
    }

    std::array<float,256> pulse_hist;
    for (int i=0; i<HcalConst::maxPSshapeBin; i++) {
      acc25nsVec.push_back(0);
      diff25nsItvlVec.push_back(0);
    }
    for (int i=0; i<HcalConst::nsPerBX; i++) {
      accVarLenIdxZEROVec.push_back(0);
      diffVarItvlIdxZEROVec.push_back(0);
      accVarLenIdxMinusOneVec.push_back(0);
      diffVarItvlIdxMinusOneVec.push_back(0);
    }

    for(int i=0;i<HcalConst::maxPSshapeBin;++i) pulse_hist[i] = Shape_[i];
    // Accumulate 25ns for each starting point of 0, 1, 2, 3...
    for(int i=0; i<HcalConst::maxPSshapeBin; ++i){
      for(int j=i; j<i+HcalConst::nsPerBX; ++j){  //sum over HcalConst::nsPerBXns from point i
	acc25nsVec[i] += ( j < HcalConst::maxPSshapeBin? pulse_hist[j] : pulse_hist[HcalConst::maxPSshapeBin-1]);
      }
      diff25nsItvlVec[i] = ( i+HcalConst::nsPerBX < HcalConst::maxPSshapeBin? pulse_hist[i+HcalConst::nsPerBX] - pulse_hist[i] : pulse_hist[HcalConst::maxPSshapeBin-1] - pulse_hist[i]);
    }
    // Accumulate different ns for starting point of index either 0 or -1
    for(int i=0; i<HcalConst::nsPerBX; ++i){
      if( i==0 ){
	accVarLenIdxZEROVec[0] = pulse_hist[0];
	accVarLenIdxMinusOneVec[i] = pulse_hist[0];
      } else{
	accVarLenIdxZEROVec[i] = accVarLenIdxZEROVec[i-1] + pulse_hist[i];
	accVarLenIdxMinusOneVec[i] = accVarLenIdxMinusOneVec[i-1] + pulse_hist[i-1];
      }
      diffVarItvlIdxZEROVec[i] = pulse_hist[i+1] - pulse_hist[0];
      diffVarItvlIdxMinusOneVec[i] = pulse_hist[i] - pulse_hist[0];
    }

}

