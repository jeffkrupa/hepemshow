#include "ad_type.h"

#ifndef VX_EDEP_HIST_HH
#define VX_EDEP_HIST_HH

// Simple per-event histogram of edep vs vx per layer.
namespace VxEdepHist {
  void ResetRun(int nLayers=50, int nBins=40, double vxMin=-1.0, double vxMax=1.0);
  void Fill(int layer, double vx, const G4double& edep);
  void EndEvent();
  void FlushRun();
}

#endif // VX_EDEP_HIST_HH
