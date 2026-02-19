#include "VxEdepHist.hh"

#include <algorithm>
#include <fstream>
#include <vector>

namespace {
  struct VxHistState {
    int nLayers = 50;
    int nBins = 20;
    double vxMin = -1.0;
    double vxMax = 1.0;
    double invBinWidth = 1.0;
    std::vector<double> bins;
    std::vector<double> layerSum;
    std::vector<double> meanEdep;
    std::vector<double> m2Edep;
    std::vector<double> meanFrac;
    std::vector<double> m2Frac;
    int nEvents = 0;
    bool inited = false;
    std::ofstream stream;

    void ensureStream() {
      if (inited) return;
      stream.open("vx_edep_hist.csv", std::ios::out);
      stream << "layer,bin,vx_low,vx_high,mean_edep,var_edep,mean_frac,var_frac,n_events\n";
      inited = true;
    }

    void reset(int layers, int binsCount, double vmin, double vmax) {
      nLayers = layers;
      nBins = binsCount;
      vxMin = vmin;
      vxMax = vmax;
      invBinWidth = (vxMax > vxMin) ? (nBins / (vxMax - vxMin)) : 1.0;
      bins.assign(nLayers * nBins, 0.0);
      layerSum.assign(nLayers, 0.0);
      meanEdep.assign(nLayers * nBins, 0.0);
      m2Edep.assign(nLayers * nBins, 0.0);
      meanFrac.assign(nLayers * nBins, 0.0);
      m2Frac.assign(nLayers * nBins, 0.0);
      nEvents = 0;
    }

    int binIndex(double vx) const {
      if (vx <= vxMin) return 0;
      if (vx >= vxMax) return nBins - 1;
      int idx = static_cast<int>((vx - vxMin) * invBinWidth);
      if (idx < 0) idx = 0;
      if (idx >= nBins) idx = nBins - 1;
      return idx;
    }

    double binLow(int idx) const {
      const double width = (vxMax - vxMin) / nBins;
      return vxMin + idx * width;
    }

    double binHigh(int idx) const {
      const double width = (vxMax - vxMin) / nBins;
      return vxMin + (idx + 1) * width;
    }
  };

  VxHistState gState;
}

namespace VxEdepHist {

void ResetRun(int nLayers, int nBins, double vxMin, double vxMax) {
  gState.reset(nLayers, nBins, vxMin, vxMax);
}

void Fill(int layer, double vx, const G4double& edep) {
  if (layer < 0 || layer >= gState.nLayers) return;
  const double e = GET_VALUE(edep);
  if (e <= 0.0) return;
  const int bin = gState.binIndex(vx);
  gState.bins[layer * gState.nBins + bin] += e;
  gState.layerSum[layer] += e;
}

void EndEvent() {
  gState.nEvents += 1;
  for (int layer = 0; layer < gState.nLayers; ++layer) {
    const double denom = gState.layerSum[layer];
    for (int bin = 0; bin < gState.nBins; ++bin) {
      const int idx = layer * gState.nBins + bin;
      const double edep = gState.bins[idx];
      const double frac = (denom > 0.0) ? (edep / denom) : 0.0;

      // Welford update for edep
      const double deltaE = edep - gState.meanEdep[idx];
      gState.meanEdep[idx] += deltaE / gState.nEvents;
      gState.m2Edep[idx] += deltaE * (edep - gState.meanEdep[idx]);

      // Welford update for frac
      const double deltaF = frac - gState.meanFrac[idx];
      gState.meanFrac[idx] += deltaF / gState.nEvents;
      gState.m2Frac[idx] += deltaF * (frac - gState.meanFrac[idx]);
    }
  }

  std::fill(gState.bins.begin(), gState.bins.end(), 0.0);
  std::fill(gState.layerSum.begin(), gState.layerSum.end(), 0.0);
}

void FlushRun() {
  gState.ensureStream();
  for (int layer = 0; layer < gState.nLayers; ++layer) {
    for (int bin = 0; bin < gState.nBins; ++bin) {
      const int idx = layer * gState.nBins + bin;
      const double meanE = gState.meanEdep[idx];
      const double varE = (gState.nEvents > 1) ? (gState.m2Edep[idx] / (gState.nEvents - 1)) : 0.0;
      const double meanF = gState.meanFrac[idx];
      const double varF = (gState.nEvents > 1) ? (gState.m2Frac[idx] / (gState.nEvents - 1)) : 0.0;
      gState.stream << layer << ',' << bin << ','
                    << gState.binLow(bin) << ',' << gState.binHigh(bin) << ','
                    << meanE << ',' << varE << ','
                    << meanF << ',' << varF << ','
                    << gState.nEvents << '\n';
    }
  }
}

}  // namespace VxEdepHist
