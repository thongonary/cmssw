#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <climits>
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalDeterministicFit.h"

constexpr float HcalDeterministicFit::invGpar[3];
constexpr float HcalDeterministicFit::negThresh[2];
constexpr int HcalDeterministicFit::HcalRegion[2];
constexpr float HcalDeterministicFit::rCorr[2];

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

}

void HcalDeterministicFit::setExternalPulseShape(std::string filename) {

  if (useExtPulse_) return;
  std::ifstream ifs;
  ifs.open(filename.c_str());
  assert(ifs.is_open());
  std::string line;

  int i = 0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;

    std::string tmpStr;
    std::stringstream ss(line);
    ss >> tmpStr;
    minCharge_[i] = std::atoi(tmpStr.c_str());
    ss >> tmpStr;
    maxCharge_[i] = std::atoi(tmpStr.c_str());
    for (int k=0; k<10; k++) { ss >> tmpStr; pulseFrac_[i][k] = std::atof(tmpStr.c_str()); }
    for (int k=0; k<10; k++) { ss >> tmpStr; pulseFracDeriv_[i][k] = std::atof(tmpStr.c_str()); }

    i++;
  }
  useExtPulse_=true;

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

  if (useExtPulse_ == true) {
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

  if (fTimeSlew==0)respCorr=1.0;
  else if (fTimeSlew==1) channelData.hasTimeInfo()?respCorr=rCorrSiPM[0]:respCorr=rCorr[0];
  else if (fTimeSlew==2) channelData.hasTimeInfo()?respCorr=rCorrSiPM[1]:respCorr=rCorr[1];
  else if (fTimeSlew==3)respCorr=frespCorr;

  float ch3=0;
  float ch4=0;
  float ch5=0;

  float tsShift3=0;
  float tsShift4=0;
  float tsShift5=0;

  if(useExtPulse_) {

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

    if (i3 != 0 && i4 != 0 && i5 != 0) {

      ch3=corrCharge[3]/i3;
      ch4=(i3*corrCharge[4]-n3*corrCharge[3])/(i3*i4);
      ch5=(n3*n4*corrCharge[3]-i4*nn3*corrCharge[3]-i3*n4*corrCharge[4]+i3*i4*corrCharge[5])/(i3*i4*i5);

      if (ch3<negThresh[0]) {
	ch3=negThresh[0];
	ch4=corrCharge[4]/i4;
	ch5=(i4*corrCharge[5]-n4*corrCharge[4])/(i4*i5);
      }
      if (ch5<negThresh[0] && ch4>negThresh[1]) {
	double ratio = (corrCharge[4]-ch3*i3)/(corrCharge[5]-negThresh[0]*i5);
	if (ratio < 5 && ratio > 0.5) {
	  //double invG = invGpar[0]+invGpar[1]*std::sqrt(2*std::log(invGpar[2]/ratio));
	  float iG=0;
	  //getLandauFrac(-invG,-invG+tsWidth,iG);
	  if (iG != 0 ) {
	    ch4=(corrCharge[4]-ch3*n3)/(iG);
	    //tsShift4=invG;
	  }
	}
      }
    }

    if (ch4<1) {
      ch4=0;
    }

  } else {
    //default Run2

    if(applyTimeSlew_) {

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

    if (i3 != 0 && i4 != 0 && i5 != 0) {

      ch3=corrCharge[3]/i3;
      ch4=(i3*corrCharge[4]-n3*corrCharge[3])/(i3*i4);
      ch5=(n3*n4*corrCharge[3]-i4*nn3*corrCharge[3]-i3*n4*corrCharge[4]+i3*i4*corrCharge[5])/(i3*i4*i5);

      if (ch3<negThresh[0]) {
	ch3=negThresh[0];
	ch4=corrCharge[4]/i4;
	ch5=(i4*corrCharge[5]-n4*corrCharge[4])/(i4*i5);
      }
      if (ch5<negThresh[0] && ch4>negThresh[1]) {
	double ratio = (corrCharge[4]-ch3*i3)/(corrCharge[5]-negThresh[0]*i5);
	if (ratio < 5 && ratio > 0.5) {
	  double invG = invGpar[0]+invGpar[1]*std::sqrt(2*std::log(invGpar[2]/ratio));
	  float iG=0;
	  getLandauFrac(-invG,-invG+tsWidth,iG);
	  if (iG != 0 ) {
	    ch4=(corrCharge[4]-ch3*n3)/(iG);
	    tsShift4=invG;
	  }
	}
      }
    }

    if (ch4<1) {
      ch4=0;
    }

  } ////

  reconstructedEnergy=ch4*respCorr*gainCorr;
  previousEnergy=ch3*respCorr*gainCorr;
  nextEnergy=ch5*respCorr*gainCorr;
  reconstructedTime=tsShift4;

}
