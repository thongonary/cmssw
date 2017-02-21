#include <algorithm>

#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/SimpleHBHEPhase1Algo.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCorrectionFunctions.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HBHERecHitAuxSetter.h"

#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/METReco/interface/HcalPhase1FlagLabels.h"


// Maximum fractional error for calculating Method 0
// pulse containment correction
constexpr float PulseContainmentFractionalError = 0.002f;

bool doCout=false;

SimpleHBHEPhase1Algo::SimpleHBHEPhase1Algo(
    const int firstSampleShift,
    const int samplesToAdd,
    const float phaseNS,
    const float timeShift,
    const bool correctForPhaseContainment,
    const int pulseShapeType,
    std::unique_ptr<PulseShapeFitOOTPileupCorrection> m2,
    std::unique_ptr<HcalDeterministicFit> detFit,
    std::unique_ptr<DoMahiAlgo> mahi)
    : pulseCorr_(PulseContainmentFractionalError),
      firstSampleShift_(firstSampleShift),
      samplesToAdd_(samplesToAdd),
      phaseNS_(phaseNS),
      timeShift_(timeShift),
      runnum_(0),
      corrFPC_(correctForPhaseContainment),
      pulseShapeType_(pulseShapeType),
      psFitOOTpuCorr_(std::move(m2)),
      hltOOTpuCorr_(std::move(detFit)),
      psFitMAHIOOTpuCorr_(std::move(mahi))
{
}

void SimpleHBHEPhase1Algo::beginRun(const edm::Run& r,
                                    const edm::EventSetup& es)
{
    runnum_ = r.run();
    pulseCorr_.beginRun(es);
}

void SimpleHBHEPhase1Algo::endRun()
{
    runnum_ = 0;
    pulseCorr_.endRun();
}

HBHERecHit SimpleHBHEPhase1Algo::reconstruct(const HBHEChannelInfo& info,
                                             const HcalRecoParam* params,
                                             const HcalCalibrations& calibs,
                                             const bool isData)
{
    HBHERecHit rh;

    const HcalDetId channelId(info.id());
    char *cmssw = getenv("CMSSW_BASE");

    const unsigned cssize = info.nSamples();
    double tsTOTen = 0;
    for(unsigned int ip=0; ip<cssize; ++ip){
      if( ip >= (unsigned) 10 ) continue;

      double charge = info.tsRawCharge(ip);
      double ped = info.tsPedestal(ip);
      double gain = info.tsGain(ip);

      double energy = charge*gain;
      double peden = ped*gain;

      tsTOTen += energy - peden;

    }


    // Calculate "Method 0" quantities
    float m0t = 0.f, m0E = 0.f;
    {
        int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int nSamplesToAdd = params ? params->samplesToAdd() : samplesToAdd_;
        const double fc_ampl = info.chargeInWindow(ibeg, ibeg + nSamplesToAdd);
        const bool applyContainment = params ? params->correctForPhaseContainment() : corrFPC_;
        const float phasens = params ? params->correctionPhaseNS() : phaseNS_;
        m0E = m0Energy(info, fc_ampl, applyContainment, phasens, nSamplesToAdd);
        m0E *= hbminusCorrectionFactor(channelId, m0E, isData);
        m0t = m0Time(info, fc_ampl, calibs, nSamplesToAdd);
    }

    if(doCout && tsTOTen>20 && (!info.hasTimeInfo())) std::cout << " ============================================================" << std::endl;

    float m2t = 0.f, m2E = 0.f, chi2 = -1.f, m3Ets3 = 0.f, m3Ets5 = 0.f;
    float m10E = 0.f, chi2_mahi = -1.f;
    float m3t = 0.f, m3E = 0.f;
    const PulseShapeFitOOTPileupCorrection* method2 = psFitOOTpuCorr_.get();
    const HcalDeterministicFit* method3 = hltOOTpuCorr_.get();
    DoMahiAlgo* mahi = psFitMAHIOOTpuCorr_.get();
    
    bool useTriple = false;
    
    // for now run only on the barrel
    if(!info.hasTimeInfo()) {

    // Run "Method 2"
    if (method2)
    {
    std::cout << "m2 enabled\n";
	if(pulseShapeType_==1) {
	  if(doCout && tsTOTen>20) std::cout << "METHOD2 = setting up the default shape=" << pulseShapeType_ << std::endl;
	  psFitOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(info.recoShape()),!info.hasTimeInfo()); // this is the standard 105
	}

	if(pulseShapeType_==2) {
	  if(doCout && tsTOTen>20) std::cout << "METHOD2 = setting up the CSV=" << pulseShapeType_ << std::endl;
	  if(!info.hasTimeInfo()) psFitOOTpuCorr_->newSetPulseShapeTemplate(((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HBHE.csv").c_str(),!info.hasTimeInfo()); // this is the CSV 105
	  if(info.hasTimeInfo()) psFitOOTpuCorr_->newSetPulseShapeTemplate(((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HE_SIPM.csv").c_str(),!info.hasTimeInfo()); // here is the CSV 203
	}
	
	if(pulseShapeType_==3) 
    {
	  if(doCout && tsTOTen>20) std::cout << "METHOD3 = setting up the LAG pulse type=" << pulseShapeType_ << std::endl;
        
      if(!isData )
      {

	    // this means MC
	    if(!info.hasTimeInfo()) psFitOOTpuCorr_->newSetPulseShapeTemplate(((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HB_MC.csv").c_str(),!info.hasTimeInfo()); // this is the LAG, MC
	    if(info.hasTimeInfo()) psFitOOTpuCorr_->newSetPulseShapeTemplate(((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HE_MC_HPD.csv").c_str(),!info.hasTimeInfo()); // this is the LAG, MC // this is a placeholder
	  }
	  if(isData)
      {
	    // this means data
	    if(!info.hasTimeInfo()) psFitOOTpuCorr_->newSetPulseShapeTemplate(((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HB_Dat.csv").c_str(),!info.hasTimeInfo()); // this is the LAG, Data
	    if(info.hasTimeInfo()) psFitOOTpuCorr_->newSetPulseShapeTemplate(((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HE_Dat_HPD.csv").c_str(),!info.hasTimeInfo()); // this is the LAG, Data // this is a placeholder
	  }
	}

        // "phase1Apply" call below sets m2E, m2t, useTriple, and chi2.
        // These parameters are pased by non-const reference.
        method2->phase1Apply(info, m2E, m2t, useTriple, chi2);
        m2E *= hbminusCorrectionFactor(channelId, m2E, isData);
    }

    // Run "Method 3"
    if (method3)
    {
     std::cout << "simple " << pulseShapeType_ << std::endl;
      if(pulseShapeType_==0) {
	if(doCout && tsTOTen>20) std::cout << "METHOD3 = setting up the default M3 landau shape =" << pulseShapeType_ << std::endl;
	hltOOTpuCorr_->setExternalPulseShape(0);
      } else if (pulseShapeType_==1) {
	if(doCout && tsTOTen>20) std::cout << "METHOD3 = setting up the csv M3 landau =" << pulseShapeType_ << std::endl;
//	hltOOTpuCorr_->setExternalPulseShape(((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_M3_HPD.csv").c_str());
	hltOOTpuCorr_->setExternalPulseShape(1);
      } else if (pulseShapeType_==2) {
	if(doCout && tsTOTen>20) std::cout << "METHOD3 = setting up the M2 105 CSV =" << pulseShapeType_ << std::endl;
	if (!info.hasTimeInfo()) hltOOTpuCorr_->setExternalPulseShape(2);  // this is the CSV 105  HB HPD.
	if(info.hasTimeInfo()) hltOOTpuCorr_->setExternalPulseShape(4); // this is the CSV 203  HE SiPM. 
      }

    //   "phase1Apply" sets m3E and m3t (pased by non-const reference)

      method3->phase1Apply(info, m3E, m3t,m3Ets3,m3Ets5);
      m3E *= hbminusCorrectionFactor(channelId, m3E, isData);
      std::cout << "shape = " << pulseShapeType_ << "\t m3E = " << m3E << std::endl;
    }

    // Run "Mahi"
    /*float m10t = 0.f, */
    //    bool useTriple_mahi = false;
    //    const DoMahiAlgo* mahi = psFitMAHIOOTpuCorr_.get();
    
    
    if (mahi)
    {
    std::cout << "mahi enabled\n";

      if(!info.hasTimeInfo()) {
	if(doCout && tsTOTen>20) std::cout << "MAHI = setting up the LAG pulse type=" << pulseShapeType_ << std::endl;

	if(!isData ){
	  // this means MC
	  if(!info.hasTimeInfo()) psFitMAHIOOTpuCorr_->setPulseShapeTemplate(true,((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HB_MC.csv").c_str()); // this is the LAG, MC
	  if(info.hasTimeInfo()) psFitMAHIOOTpuCorr_->setPulseShapeTemplate(true,((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HE_MC_HPD.csv").c_str()); // this is the LAG, MC // this is a placeholder
	}
	if(isData){
	  // this means data
	  if(!info.hasTimeInfo()) psFitMAHIOOTpuCorr_->setPulseShapeTemplate(true,((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HB_Dat.csv").c_str()); // this is the LAG, Data
	  if(info.hasTimeInfo()) psFitMAHIOOTpuCorr_->setPulseShapeTemplate(true,((std::string)cmssw+"/src/CalibCalorimetry/HcalAlgos/data/pulse_shape_HE_Dat_HPD.csv").c_str()); // this is the LAG, Data // this is a placeholder
	}
	mahi->phase1Apply(info,m10E,chi2_mahi);
	m10E *= hbminusCorrectionFactor(channelId, m10E, isData);
      }
    }
    } // end if HB

    // Finally, construct the rechit
    float rhE = m0E;
    float rht = m0t;
    if (method2)
    {
        rhE = m2E;
        rht = m2t;
    }
    else if (method3)
    {
        rhE = m3E;
        rht = m3t;
    }

    else if (mahi)
    {
        rhE = m10E;
    }

    float tdcTime = info.soiRiseTime();
    if (!HcalSpecialTimes::isSpecial(tdcTime))
        tdcTime += timeShift_;
    rh = HBHERecHit(channelId, rhE, rht, tdcTime);
    rh.setRawEnergy(m0E);
    rh.setAuxEnergy(m3E);
    rh.setChiSquared(chi2);
//
//    rh.setRawEnergy(m3Ets3);
//    rh.setAuxEnergy(m3E);
//    rh.setChiSquared(m3Ets5);

    // Set rechit aux words
    HBHERecHitAuxSetter::setAux(info, &rh);

    // Set some rechit flags (here, for Method 2)
    if (useTriple)
       rh.setFlagField(1, HcalPhase1FlagLabels::HBHEPulseFitBit);

    return rh;
}

float SimpleHBHEPhase1Algo::hbminusCorrectionFactor(const HcalDetId& cell,
                                                    const float energy,
                                                    const bool isRealData) const
{
    float corr = 1.f;
    if (isRealData && runnum_ > 0)
        if (cell.subdet() == HcalBarrel)
        {
            const int ieta = cell.ieta();
            const int iphi = cell.iphi();
            corr = hbminus_special_ecorr(ieta, iphi, energy, runnum_);
        }
    return corr;
}

float SimpleHBHEPhase1Algo::m0Energy(const HBHEChannelInfo& info,
                                     const double fc_ampl,
                                     const bool applyContainmentCorrection,
                                     const double phaseNs,
                                     const int nSamplesToAdd)
{
    int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
    if (ibeg < 0)
        ibeg = 0;
    double e = info.energyInWindow(ibeg, ibeg + nSamplesToAdd);

    // Pulse containment correction
    {    
        double corrFactor = 1.0;
        if (applyContainmentCorrection)
            corrFactor = pulseCorr_.get(info.id(), nSamplesToAdd, phaseNs)->getCorrection(fc_ampl);
        e *= corrFactor;
    }

    return e;
}

float SimpleHBHEPhase1Algo::m0Time(const HBHEChannelInfo& info,
                                   const double fc_ampl,
                                   const HcalCalibrations& calibs,
                                   const int nSamplesToExamine) const
{
    float time = -9999.f; // historic value

    const unsigned nSamples = info.nSamples();
    if (nSamples > 2U)
    {
        const int soi = info.soi();
        int ibeg = soi + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int iend = ibeg + nSamplesToExamine;
        unsigned maxI = info.peakEnergyTS(ibeg, iend);
        if (maxI < HBHEChannelInfo::MAXSAMPLES)
        {
            if (!maxI)
                maxI = 1U;
            else if (maxI >= nSamples - 1U)
                maxI = nSamples - 2U;

            // The remaining code in this scope emulates
            // the historic algorithm
            float t0 = info.tsEnergy(maxI - 1U);
            float maxA = info.tsEnergy(maxI);
            float t2 = info.tsEnergy(maxI + 1U);

            // Handle negative excursions by moving "zero"
            float minA = t0;
            if (maxA < minA) minA = maxA;
            if (t2 < minA)   minA=t2;
            if (minA < 0.f) { maxA-=minA; t0-=minA; t2-=minA; }
            float wpksamp = (t0 + maxA + t2);
            if (wpksamp) wpksamp = (maxA + 2.f*t2) / wpksamp;
            time = (maxI - soi)*25.f + timeshift_ns_hbheho(wpksamp);

            // Legacy QIE8 timing correction
            time -= HcalTimeSlew::delay(std::max(1.0, fc_ampl),
                                        HcalTimeSlew::Medium);
            // Time calibration
            time -= calibs.timecorr();
        }
    }
    return time;
}
