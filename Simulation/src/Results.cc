#include "ad_type.h"


#include "Results.hh"

#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>

void WriteResults(struct Results& res, int numEvents, int seed) {
  // for the histograms, bring them to be mean per event and write
  const G4double norm = numEvents > 0 ? 1.0/numEvents : 1.0;
  res.fEdepPerLayer.Scale(norm);
  res.fEdepGapPerLayer.Scale(norm);
  res.fGammaTrackLenghtPerLayer.Scale(norm);
  res.fElPosTrackLenghtPerLayer.Scale(norm);

  res.fEdepPerLayer.WriteToFile(false);
  std::ofstream edeps("edeps_" + std::to_string(seed));
  //std::ofstream edeps("edeps");
  const int nLayers = (int)res.fEdepPerLayer_Acc.size();
  for(int i=0; i<nLayers; i++){
     edeps << std::setprecision(14) << res.fEdepPerLayer_Acc[i].getMean() << " " << res.fEdepPerLayer_Acc[i].getMeanSq();
     #if CODI_FORWARD
        edeps << " " << res.fEdepPerLayer_AccD[i].getMean() << " " << res.fEdepPerLayer_AccD[i].getMeanSq();
     #endif
     edeps << "\n";
  }
  edeps.close();

  // per-layer GAP energy (sampled signal), same column layout as edeps_<seed>:
  // "mean meanSq [mean_dot meanSq_dot]" per layer. Absorber energy = combined - gap.
  std::ofstream edepsGap("edeps_gap_" + std::to_string(seed));
  const int nLayersGap = (int)res.fEdepGapPerLayer_Acc.size();
  for(int i=0; i<nLayersGap; i++){
     edepsGap << std::setprecision(14) << res.fEdepGapPerLayer_Acc[i].getMean() << " " << res.fEdepGapPerLayer_Acc[i].getMeanSq();
     #if CODI_FORWARD
        edepsGap << " " << res.fEdepGapPerLayer_AccD[i].getMean() << " " << res.fEdepGapPerLayer_AccD[i].getMeanSq();
     #endif
     edepsGap << "\n";
  }
  edepsGap.close();

  #ifdef CODI_REVERSE
     // legacy aggregate format (kept for backward compatibility / regression checks):
     // sum over layers == derivative w.r.t. the shared uniform thickness.
     std::ofstream barInputs("barInputs");
     barInputs << std::setprecision(14);
     barInputs << res.barThicknessAbsorber.getMean() << " " << res.barThicknessAbsorber.getVar() << "\n";
     barInputs << res.barThicknessGap.getMean() << " " << res.barThicknessGap.getVar() << "\n";
     barInputs << res.barParticleEnergy.getMean() << " " << res.barParticleEnergy.getVar() << "\n";
     barInputs.close();

     // per-layer gradients: rows 0..N-1 -> d/d_absThick[i], rows N..2N-1 -> d/d_gapThick[i],
     // final row -> d/d_energy. Each row is "mean var".
     const int NL = (int)res.barAbsThick.size();
     std::ofstream barPL("barInputsPerLayer");
     barPL << std::setprecision(14);
     barPL << "# N=" << NL << " ; rows 0..N-1: d/d_absThick[i] (mean var); rows N..2N-1: d/d_gapThick[i]; last row: d/d_energy\n";
     for(int i=0; i<NL; i++){
        barPL << res.barAbsThick[i].getMean() << " " << res.barAbsThick[i].getVar() << "\n";
     }
     for(int i=0; i<NL; i++){
        barPL << res.barGapThick[i].getMean() << " " << res.barGapThick[i].getVar() << "\n";
     }
     barPL << res.barParticleEnergy.getMean() << " " << res.barParticleEnergy.getVar() << "\n";
     barPL.close();
  #endif


  res.fGammaTrackLenghtPerLayer.WriteToFile(false);
  res.fElPosTrackLenghtPerLayer.WriteToFile(false);

  //
  res.fEdepAbs  = res.fEdepAbs*norm;
  res.fEdepAbs2 = res.fEdepAbs2*norm;
  const G4double rmsEAbs = std::sqrt(std::abs(res.fEdepAbs2 - res.fEdepAbs*res.fEdepAbs));

  res.fEdepGap  = res.fEdepGap*norm;
  res.fEdepGap2 = res.fEdepGap2*norm;
  const G4double rmsEGap = std::sqrt(std::abs(res.fEdepGap2 - res.fEdepGap*res.fEdepGap));


  // the secondary type and step number statistics
  std::cout << std::endl;
  std::cout << " --- Results::WriteResults ---------------------------------- " << std::endl;
  std::cout << std::setprecision(6);
  std::cout << " Absorber: mean Edep = " << res.fEdepAbs << " [MeV] and  Std-dev = " << rmsEAbs << " [MeV]"<< std::endl;
  std::cout << " Gap     : mean Edep = " << res.fEdepGap << " [MeV] and  Std-dev = " << rmsEGap << " [MeV]"<< std::endl;

  std::cout << std::endl;
  std::cout << std::setprecision(14);
  std::cout << " Mean number of gamma       " << res.fNumSecGamma*norm    << std::endl;
  std::cout << " Mean number of e-          " << res.fNumSecElectron*norm << std::endl;
  std::cout << " Mean number of e+          " << res.fNumSecPositron*norm << std::endl;

  std::cout << std::endl;
  std::cout << std::setprecision(6)
            << " Mean number of e-/e+ steps " << res.fNumStepsElPos*norm  << std::endl;
  std::cout << " Mean number of gamma steps " << res.fNumStepsGamma*norm  << std::endl;
  std::cout << " ------------------------------------------------------------\n";

  #ifdef CODI_REVERSE
    G4double::getTape().printStatistics(std::cout);
  #endif

}
