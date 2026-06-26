#include "ad_type.h"


#include "SteppingLoop.hh"

// G4HepEm includes
#include "G4HepEmTLData.hh"
#include "G4HepEmState.hh"
#include "G4HepEmTrack.hh"


#include "G4HepEmData.hh"
#include "G4HepEmMatCutData.hh"


// application local includes
#include "TrackStack.hh"
#include "Physics.hh"
#include "Geometry.hh"
#include "Box.hh"
#include "Results.hh"
#include "VxEdepHist.hh"

#include <fstream>
#include <atomic>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <cstdlib>

#ifndef MICRO_AUDIT_MAX
#define MICRO_AUDIT_MAX 2e6   // hard cap on lines we'll write
#endif

const double Ldot_thr   = 5e999;         // tweak to match your scale
const double Edot_thr   = 5e999;         // idem

namespace {
  G4double kDefaultKECut = 0.5;
  G4double gKECutValue = kDefaultKECut;
  bool gEnableKECut = false;
}

bool outputall = (std::getenv("HEPEMSHOW_OUTPUT_ALL") && std::string(std::getenv("HEPEMSHOW_OUTPUT_ALL")) != "0");
bool outputboundarylayers = (std::getenv("HEPEMSHOW_OUTPUT_ALL") && std::string(std::getenv("HEPEMSHOW_OUTPUT_ALL")) != "0");
bool outputpingpongtracks = false;
bool outputvx = false;
;
struct MicroAudit {
  static std::ofstream& stream() {
    static std::ofstream s("micro_audit.csv", std::ios::out);
    static bool inited = false;
    if (!inited) {
      s << "event,charge,trackID,parentID,parentStep,step,layer,winnerIdx,is_stopgrad,"
        << "onBoundary,blAxis,"
        << "gX_pre,gY_pre,gZ_pre,"
        << "gX_pre_dot,gY_pre_dot,gZ_pre_dot,"
        << "gX,gY,gZ,"
        << "gX_postMSC,gY_postMSC,gZ_postMSC,"
        << "vx_prestep,vy_prestep,vz_prestep,"
        << "vx,vy,vz,preStepSafety,postStepSafety,"
        << "numIAleft_0_prePerform,numIAleft_1_prePerform,numIAleft_2_prePerform,"
        << "numIAleft_0_prePerform_dot,numIAleft_1_prePerform_dot,numIAleft_2_prePerform_dot,"
        << "numIAleft_0,numIAleft_1,numIAleft_2,"
        << "numIAleft_0_dot,numIAleft_1_dot,numIAleft_2_dot,"
        << "mfp_0,mfp_1,mfp_2,"
        << "mfp_0_dot,mfp_1_dot,mfp_2_dot,"
        << "stepLength,stepLength_dot,"
        << "pStepLength,pStepLength_dot,"
        << "edep,edep_dot,EKin,EKin_dot,"
        << "distB,distB_dot,distP,distP_dot\n";
      inited = true;
    }
    return s;
  }
  static std::atomic<int>& counter() {
    static std::atomic<int> c{0};
    return c;
  }
  static void logLine(const std::string& line) {
    int k = counter().fetch_add(1);
    if (k < MICRO_AUDIT_MAX) stream() << line << '\n';
  }
};



//template<typename Expr>
//inline G4double stop_grad(const Expr& x) {
  ////return G4double(x);
//  return G4double(GET_VALUE(x));
//}

namespace {
  int ReadEnvInt(const char* key, int defaultValue) {
    const char* raw = std::getenv(key);
    if (raw == nullptr || *raw == '\0') {
      return defaultValue;
    }
    char* endPtr = nullptr;
    const long parsed = std::strtol(raw, &endPtr, 10);
    if (endPtr == raw) {
      return defaultValue;
    }
    return static_cast<int>(parsed);
  }

  bool IsDebugLoggingEnabled() {
    static const int kEnabled = ReadEnvInt("HEPEMSHOW_ENABLE_DEBUG_LOGS", 0);
    return kEnabled != 0;
  }

  inline double DotValueOrZero(const G4double& v) {
#ifdef CODI_FORWARD
    return static_cast<double>(GET_DOTVALUE(v));
#else
    (void)v;
    return 0.0;
#endif
  }

  struct StepChoiceAudit {
    static int TargetTrackId() {
      static const int kTargetTrack = ReadEnvInt("HEPEMSHOW_STEP_CHOICE_TRACK_ID", -1);
      return kTargetTrack;
    }

    static bool ShouldLog(const G4HepEmTrack& track) {
      if (!IsDebugLoggingEnabled()) {
        return false;
      }
      const int target = TargetTrackId();
      return target >= 0 && track.GetID() == target;
    }

    static std::ofstream& stream() {
      static std::ofstream s("step_choice_debug.csv", std::ios::out | std::ios::trunc);
      static bool inited = false;
      if (!inited) {
        s << "event,trackID,parentID,step,stage,physicsWinnerIdx,"
          << "branchWinner,isBoundaryBranch,"
          << "onBoundaryBefore,onBoundaryAfter,wasOnBoundary,"
          << "distToBoundary,distToBoundary_dot,"
          << "distToPhysics,distToPhysics_dot,"
          << "stepLength,stepLength_dot,"
          << "preStepSafety,preStepSafety_dot,"
          << "vx,vx_dot,vy,vy_dot,vz,vz_dot,"
          << "localX,localX_dot,localY,localY_dot,localZ,localZ_dot,"
          << "tx,tx_dot,ty,ty_dot,tz,tz_dot,tmin,tmin_dot,"
          << "winnerAxis,formulaMatchesDistB\n";
        inited = true;
      }
      return s;
    }

    static void LogDecision(int eventID, G4HepEmTrack& track, int step, const char* stage,
                            bool onBoundaryBefore, bool onBoundaryAfter, bool wasOnBoundary,
                            const G4double& distToBoundary, const G4double& distToPhysics,
                            const G4double& stepLength, const G4double& preStepSafety,
                            const Box* currentVolume, const G4double* localPosition, const G4double* curDirection) {
      if (!ShouldLog(track)) {
        return;
      }
      G4double tx = 1.0E+20;
      G4double ty = 1.0E+20;
      G4double tz = 1.0E+20;
      G4double tmin = 1.0E+20;
      char winnerAxis = '?';
      int formulaMatchesDistB = 0;
      if (currentVolume != nullptr && localPosition != nullptr && curDirection != nullptr) {
        const G4double hx = currentVolume->GetHalfLength(0);
        const G4double hy = currentVolume->GetHalfLength(1);
        const G4double hz = currentVolume->GetHalfLength(2);
        const G4double vxLoc = curDirection[0];
        const G4double vyLoc = curDirection[1];
        const G4double vzLoc = curDirection[2];
        tx = (vxLoc == 0) ? 1.0E+20 : (G4double)((std::copysign(hx, vxLoc) - localPosition[0]) / vxLoc);
        ty = (vyLoc == 0) ? tx : (G4double)((std::copysign(hy, vyLoc) - localPosition[1]) / vyLoc);
        const G4double txy = std::min(tx, ty);
        tz = (vzLoc == 0) ? txy : (G4double)((std::copysign(hz, vzLoc) - localPosition[2]) / vzLoc);
        tmin = std::min(txy, tz);
        const double txVal = GET_VALUE(tx);
        const double tyVal = GET_VALUE(ty);
        const double tzVal = GET_VALUE(tz);
        if (txVal <= tyVal && txVal <= tzVal) {
          winnerAxis = 'x';
        } else if (tyVal <= txVal && tyVal <= tzVal) {
          winnerAxis = 'y';
        } else {
          winnerAxis = 'z';
        }
        const double absDiff = std::abs(GET_VALUE(tmin - distToBoundary));
        formulaMatchesDistB = (absDiff < 1.0E-10) ? 1 : 0;
      }
      const bool isBoundaryBranch = GET_VALUE(distToPhysics) >= GET_VALUE(distToBoundary);
      std::ofstream& s = stream();
      const G4double vx = track.GetDirection()[0];
      const G4double vy = track.GetDirection()[1];
      const G4double vz = track.GetDirection()[2];
      s << eventID << ','
        << track.GetID() << ','
        << track.GetParentID() << ','
        << step << ','
        << stage << ','
        << track.GetWinnerProcessIndex() << ','
        << (isBoundaryBranch ? "boundary" : "physics") << ','
        << (isBoundaryBranch ? 1 : 0) << ','
        << (onBoundaryBefore ? 1 : 0) << ','
        << (onBoundaryAfter ? 1 : 0) << ','
        << (wasOnBoundary ? 1 : 0) << ','
        << GET_VALUE(distToBoundary) << ',' << DotValueOrZero(distToBoundary) << ','
        << GET_VALUE(distToPhysics) << ',' << DotValueOrZero(distToPhysics) << ','
        << GET_VALUE(stepLength) << ',' << DotValueOrZero(stepLength) << ','
        << GET_VALUE(preStepSafety) << ',' << DotValueOrZero(preStepSafety) << ','
        << GET_VALUE(vx) << ',' << DotValueOrZero(vx) << ','
        << GET_VALUE(vy) << ',' << DotValueOrZero(vy) << ','
        << GET_VALUE(vz) << ',' << DotValueOrZero(vz) << ','
        << GET_VALUE(localPosition[0]) << ',' << DotValueOrZero(localPosition[0]) << ','
        << GET_VALUE(localPosition[1]) << ',' << DotValueOrZero(localPosition[1]) << ','
        << GET_VALUE(localPosition[2]) << ',' << DotValueOrZero(localPosition[2]) << ','
        << GET_VALUE(tx) << ',' << DotValueOrZero(tx) << ','
        << GET_VALUE(ty) << ',' << DotValueOrZero(ty) << ','
        << GET_VALUE(tz) << ',' << DotValueOrZero(tz) << ','
        << GET_VALUE(tmin) << ',' << DotValueOrZero(tmin) << ','
        << winnerAxis << ','
        << formulaMatchesDistB
        << '\n';
    }
  };

  std::unordered_set<int> gDisabledTrackGradients;
  int gGradientStopMode = 2;
  bool gGrazingStopsFullTrack = true;
  bool gEnableBackwardBoundaryStop = false;
  bool gEnableMscDisplacement = true;
  G4double gMscDisplacementSafetyFloor = 0.0;
  int gSameBoundaryStopThreshold = 0;
  G4double gSameBoundaryPositionTolerance = 1.0E-6;
  int gSameBoundaryMinFlips = 1;
  bool gSameBoundaryFullTrackStop = false;
  int gSameBoundaryHardStopThreshold = 0;
  const G4double kDefaultNearBoundarySafety = 5.0E-6;

  struct VxStats {
    static std::ofstream& stream() {
      static std::ofstream s("vx_edep_KE.csv", std::ios::out);
      static bool inited = false;
      if (!inited) {
        s << "vx,edep,KE\n";
        inited = true;
      }
      return s;
    }
    static void logLine(const std::string& line) {
      stream() << line << '\n';
    }
  };

  struct BoundaryStats {
    static std::ofstream& stream() {
      static std::ofstream s("boundary_stats.csv", std::ios::out);
      static bool inited = false;
      if (!inited) {
        s << "event,numTracks,totalSteps,totalStepsBack,totalStepsForward,boundarySteps,boundaryStepsBack,boundaryStepsForward,edep,steplength_dot";
        for (int i = 0; i < 50; ++i) {
          s << ",layer_edep_" << i;
        }
        s << ",maxConsecBoundary,maxConsecBoundaryBack,maxConsecBoundaryForward\n";
        inited = true;
      }
      return s;
    }
    static void logLine(const std::string& line) {
      stream() << line << '\n';
    }
  };

  struct PingPongTrackStats {
    static std::ofstream& stream() {
      static std::ofstream s("track_pingpong_stats.csv", std::ios::out);
      static bool inited = false;
      if (!inited) {
        s << "event,trackID,parentID,creationStep,charge,"
          << "totalSteps,totalStepsBack,totalStepsForward,"
          << "boundarySteps,boundaryStepsBack,boundaryStepsForward,"
          << "nFlipVx,maxConsecBL,maxConsecBLBack,maxConsecBLForward,"
          << "initVx,finalVx,isInitBackward,isFinalBackward,"
          << "blEkinCount,blEkinMean,blEkinStd,blEkinMin,blEkinMax\n";
        inited = true;
      }
      return s;
    }
    static void logLine(const std::string& line) {
      stream() << line << '\n';
    }
  };

  struct BoundaryEventAgg {
    int numTracks = 0;
    long long totalSteps = 0;
    long long totalStepsBack = 0;
    long long totalStepsForward = 0;
    long long boundarySteps = 0;
    long long boundaryStepsBack = 0;
    long long boundaryStepsForward = 0;
    double stepLengthDot = 0.0;
    int maxConsecBoundary = 0;
    int maxConsecBoundaryBack = 0;
    int maxConsecBoundaryForward = 0;
  };

  std::unordered_map<int, BoundaryEventAgg> gBoundaryStatsByEvent;

  inline bool IsFiniteNumber(double x) {
    return std::isfinite(x);
  }

  inline double ScalarDotValue(const G4double& x) {
    #ifdef CODI_FORWARD
      return GET_DOTVALUE(x);
    #else
      (void)x;
      return 0.0;
    #endif
  }

  inline void AccumulateBoundaryTrackStats(
      int eventID,
      int totalSteps,
      int totalStepsBack,
      int totalStepsForward,
      int boundarySteps,
      int boundaryStepsBack,
      int boundaryStepsForward,
      double stepLengthDot,
      int maxConsecBoundary,
      int maxConsecBoundaryBack,
      int maxConsecBoundaryForward) {
    auto& agg = gBoundaryStatsByEvent[eventID];
    ++agg.numTracks;
    agg.totalSteps += totalSteps;
    agg.totalStepsBack += totalStepsBack;
    agg.totalStepsForward += totalStepsForward;
    agg.boundarySteps += boundarySteps;
    agg.boundaryStepsBack += boundaryStepsBack;
    agg.boundaryStepsForward += boundaryStepsForward;
    if (IsFiniteNumber(stepLengthDot)) {
      agg.stepLengthDot += stepLengthDot;
    }
    if (maxConsecBoundary > agg.maxConsecBoundary) {
      agg.maxConsecBoundary = maxConsecBoundary;
    }
    if (maxConsecBoundaryBack > agg.maxConsecBoundaryBack) {
      agg.maxConsecBoundaryBack = maxConsecBoundaryBack;
    }
    if (maxConsecBoundaryForward > agg.maxConsecBoundaryForward) {
      agg.maxConsecBoundaryForward = maxConsecBoundaryForward;
    }
  }

  inline void SanitizeTrackState(G4HepEmTrack& track) {
    G4double* position = track.GetPosition();
    track.SetPosition(stop_grad(position[0]), stop_grad(position[1]), stop_grad(position[2]));
    G4double* direction = track.GetDirection();
    track.SetDirection(stop_grad(direction[0]), stop_grad(direction[1]), stop_grad(direction[2]));
  }

  inline void SanitizeTrackStateForFullStop(G4HepEmTrack& track) {
    SanitizeTrackState(track);
    track.SetGStepLength(stop_grad(track.GetGStepLength()));
    track.SetSafety(stop_grad(track.GetSafety()));
  }

  inline G4double AbsValue(G4double v) {
    return (v < 0.0) ? -v : v;
  }

  inline int VxSign(double vx) {
    if (vx < 0.0) return -1;
    if (vx > 0.0) return 1;
    return 0;
  }

  inline double BoundaryXValue(const G4double* globalPosition, const G4double* direction, const G4double& stepLength) {
    return GET_VALUE(globalPosition[0] + stepLength * direction[0]);
  }

  inline long long BoundaryKeyFromX(double boundaryX) {
    const double tolValue = GET_VALUE(gSameBoundaryPositionTolerance);
    const double tol = tolValue > 0.0 ? tolValue : 1.0E-6;
    return static_cast<long long>(std::llround(boundaryX / tol));
  }

  inline void LogTrackPingPongStats(
      int eventID,
      G4HepEmTrack& track,
      int creationStep,
      int totalSteps,
      int totalStepsBack,
      int boundarySteps,
      int boundaryStepsBack,
      int nFlipVx,
      int maxConsecBoundary,
      int maxConsecBoundaryBack,
      int maxConsecBoundaryForward,
      double initVx,
      int blEkinCount,
      double blEkinSum,
      double blEkinSqSum,
      double blEkinMin,
      double blEkinMax) {
    if (!outputpingpongtracks) {
      return;
    }
    const int totalStepsForward = totalSteps - totalStepsBack;
    const int boundaryStepsForward = boundarySteps - boundaryStepsBack;
    const double finalVx = GET_VALUE(track.GetDirection()[0]);
    const int isInitBackward = initVx < 0.0 ? 1 : 0;
    const int isFinalBackward = finalVx < 0.0 ? 1 : 0;

    double blEkinMean = NAN;
    double blEkinStd = NAN;
    double blEkinMinOut = NAN;
    double blEkinMaxOut = NAN;
    if (blEkinCount > 0) {
      blEkinMean = blEkinSum / static_cast<double>(blEkinCount);
      const double meanSq = blEkinSqSum / static_cast<double>(blEkinCount);
      const double var = meanSq - blEkinMean * blEkinMean;
      blEkinStd = var > 0.0 ? std::sqrt(var) : 0.0;
      blEkinMinOut = blEkinMin;
      blEkinMaxOut = blEkinMax;
    }

    std::ostringstream s;
    s << eventID << ','
      << track.GetID() << ','
      << track.GetParentID() << ','
      << creationStep << ','
      << GET_VALUE(track.GetCharge()) << ','
      << totalSteps << ','
      << totalStepsBack << ','
      << totalStepsForward << ','
      << boundarySteps << ','
      << boundaryStepsBack << ','
      << boundaryStepsForward << ','
      << nFlipVx << ','
      << maxConsecBoundary << ','
      << maxConsecBoundaryBack << ','
      << maxConsecBoundaryForward << ','
      << initVx << ','
      << finalVx << ','
      << isInitBackward << ','
      << isFinalBackward << ','
      << blEkinCount << ','
      << blEkinMean << ','
      << blEkinStd << ','
      << blEkinMinOut << ','
      << blEkinMaxOut;
    PingPongTrackStats::logLine(s.str());
  }
}

void SteppingLoop::ResetDisabledGradients() {
  gDisabledTrackGradients.clear();
}

bool SteppingLoop::IsTrackGradientDisabled(const G4HepEmTrack& track) {
  if (gGradientStopMode == 0 || track.GetID() < 0) {
    return false;
  }
  return gDisabledTrackGradients.count(track.GetID()) > 0;
}

void SteppingLoop::DisableTrackGradient(G4HepEmTrack& track) {
  if (gGradientStopMode == 0) {
    return;
  }
  SanitizeTrackStateForFullStop(track);
  if (track.GetID() >= 0) {
    gDisabledTrackGradients.insert(track.GetID());
  }
}

void SteppingLoop::SetGradientStopMode(int mode) {
  if (mode < 0) mode = 0;
  if (mode > 2) mode = 2;
  gGradientStopMode = mode;
  if (gGradientStopMode == 0) {
    gDisabledTrackGradients.clear();
  }
}

void SteppingLoop::ConfigureKECut(bool enable, G4double threshold) {
  gEnableKECut = enable;
  if (gEnableKECut) {
    if (threshold < 0.0) {
      threshold = 0.0;
    }
    gKECutValue = threshold;
  } else {
    gKECutValue = kDefaultKECut;
  }
}

void SteppingLoop::ConfigureGrazingStopPolicy(bool disableFullTrack) {
  gGrazingStopsFullTrack = disableFullTrack;
}

void SteppingLoop::ConfigureBackwardBoundaryStop(bool enable) {
  gEnableBackwardBoundaryStop = enable;
}

void SteppingLoop::ConfigureMscDisplacement(bool enable, G4double postSafetyFloor) {
  gEnableMscDisplacement = enable;
  gMscDisplacementSafetyFloor = postSafetyFloor < 0.0 ? 0.0 : postSafetyFloor;
}

void SteppingLoop::ConfigureSameBoundaryStop(int repeatThreshold, G4double positionTolerance, int minFlips, bool fullTrackStop, int hardStopThreshold) {
  gSameBoundaryStopThreshold = repeatThreshold > 0 ? repeatThreshold : 0;
  gSameBoundaryPositionTolerance = positionTolerance > 0.0 ? positionTolerance : 1.0E-6;
  gSameBoundaryMinFlips = minFlips > 0 ? minFlips : 0;
  gSameBoundaryFullTrackStop = fullTrackStop;
  gSameBoundaryHardStopThreshold = hardStopThreshold > 0 ? hardStopThreshold : 0;
}

void SteppingLoop::FlushBoundaryStatsForEvent(int eventID, double eventEdep, const G4double* layerEdep, int numLayers) {
  if (!outputboundarylayers) {
    return;
  }
  auto it = gBoundaryStatsByEvent.find(eventID);
  if (it == gBoundaryStatsByEvent.end()) {
    return;
  }
  const BoundaryEventAgg& agg = it->second;
  std::ostringstream statss;
  statss << eventID << ','
         << agg.numTracks << ','
         << agg.totalSteps << ','
         << agg.totalStepsBack << ','
         << agg.totalStepsForward << ','
         << agg.boundarySteps << ','
         << agg.boundaryStepsBack << ','
         << agg.boundaryStepsForward << ','
         << eventEdep << ','
         << agg.stepLengthDot;
  for (int i = 0; i < 50; ++i) {
    if (layerEdep != nullptr && i < numLayers) {
      statss << ',' << GET_VALUE(layerEdep[i]);
    } else {
      statss << ',' << NAN;
    }
  }
  statss << ','
         << agg.maxConsecBoundary << ','
         << agg.maxConsecBoundaryBack << ','
         << agg.maxConsecBoundaryForward;
  BoundaryStats::logLine(statss.str());
  gBoundaryStatsByEvent.erase(it);
}


//
// NOTE: we always calculate the distance to boundary and the pre-step point safety
//       that is very far from being optimal. In real g4 tracking, the safety is
//       updated after each step (post-stepSafety = pre-stepSafety - "stepLength")
//       So as long as we the step-length is within the up-to-date safety we do not
//       need to re-calculate the safety and we do not need to calculate the distance
//       to boundary as for sure the step will end up far from the boundaries.
//       But here we have a simplified gometry and navigation....

void SteppingLoop::GammaStepper(G4HepEmTLData& theTLData, G4HepEmState& theState, TrackStack& theTrackStack, Geometry& theGeometry, Results& theResult, int eventID, G4double threshold, G4double threshold2) {
  // NOTE: the start tracking procedure (reset the track and the rng) was done
  G4HepEmTrack* theTrack = theTLData.GetPrimaryGammaTrack()->GetTrack();

  //
  // if this is a real primary track then I need to locate it
  // if this is a secondary then I could already know, but
  // anyway: locate in all cases to keep it simply (but slower anyway)
  //
  int  numStep       = 0;
  Box* currentVolume = nullptr;
  bool onBoundary    = false;
  int  indxLayer     = -1;
  int  indxAbs       = -1;
  G4double  localPosition[3];
  const G4double vxThreshold = AbsValue(threshold);
  const G4double nearBoundarySafety = threshold2 > 0.0 ? threshold2 : kDefaultNearBoundarySafety;
  bool stop_tracking = IsTrackGradientDisabled(*theTrack);
  if (stop_tracking) {
    DisableTrackGradient(*theTrack);
  }
  const int creationStep = theTrackStack.GetTrackCreationStep(theTrack->GetID());
  const double initVx = GET_VALUE(theTrack->GetDirection()[0]);
  int totalSteps = 0;
  int totalStepsBack = 0;
  int totalStepsForward = 0;
  int boundarySteps = 0;
  int boundaryStepsBack = 0;
  int boundaryStepsForward = 0;
  int consecBoundary = 0;
  int consecBoundaryBack = 0;
  int consecBoundaryForward = 0;
  int maxConsecBoundary = 0;
  int maxConsecBoundaryBack = 0;
  int maxConsecBoundaryForward = 0;
  int nFlipVx = 0;
  int prevVxSign = 0;
  bool hasPrevVxSign = false;
  int blEkinCount = 0;
  double blEkinSum = 0.0;
  double blEkinSqSum = 0.0;
  double blEkinMin = std::numeric_limits<double>::infinity();
  double blEkinMax = -std::numeric_limits<double>::infinity();
  std::unordered_map<long long, int> sameBoundaryHitCount;
  std::unordered_map<long long, int> sameBoundaryFlipCount;
  std::unordered_map<long long, int> sameBoundaryLastVxSign;

  while (theTrack->GetEKin() > 0.0) {
    if (stop_tracking) {
      DisableTrackGradient(*theTrack);
    }

    // calculate distance to boundary from the pre-step point: will locate the pont
    // NOTE: this should never be zero as zero means that the point is outside of the volume
    //       (taking into account the direction and tolerance)
    // NOTE: the given position will be in local coordiantes at return
    G4double* globalPosition = theTrack->GetPosition();
    G4double* curDirection   = theTrack->GetDirection();
    G4double pre_step_gX = globalPosition[0];
    G4double pre_step_gY = globalPosition[1];
    G4double pre_step_gZ = globalPosition[2];
    #ifdef CODI_FORWARD
      G4double pre_step_gX_dot = GET_DOTVALUE(globalPosition[0]);
      G4double pre_step_gY_dot = GET_DOTVALUE(globalPosition[1]);
      G4double pre_step_gZ_dot = GET_DOTVALUE(globalPosition[2]);
    #endif

    G4double pre_step_vx = curDirection[0];
    G4double pre_step_vy = curDirection[1];
    G4double pre_step_vz = curDirection[2];
    // set the local position = global position (will be local after CalculateDistanceToOut)
    Set3Vect(localPosition, globalPosition);
    Geometry::SetDistanceDebugContext(eventID, theTrack->GetID(), theTrack->GetParentID(), numStep, "gamma_pre");
    G4double distToBoundary = theGeometry.CalculateDistanceToOut(localPosition, curDirection, &currentVolume, &indxLayer, &indxAbs);
    Geometry::ClearDistanceDebugContext();
    G4double localPositionAtDistCalc[3];
    Set3Vect(localPositionAtDistCalc, localPosition);
    G4double directionAtDistCalc[3];
    Set3Vect(directionAtDistCalc, curDirection);
    // STOP HERE IF `distToBoundary = 1.0E+20` i.e. we are going out from the Calorimeter
    if (distToBoundary > 1.0E+10) {
      LogTrackPingPongStats(
          eventID,
          *theTrack,
          creationStep,
          totalSteps,
          totalStepsBack,
          boundarySteps,
          boundaryStepsBack,
          nFlipVx,
          maxConsecBoundary,
          maxConsecBoundaryBack,
          maxConsecBoundaryForward,
          initVx,
          blEkinCount,
          blEkinSum,
          blEkinSqSum,
          blEkinMin,
          blEkinMax);
      if (outputboundarylayers) {
        AccumulateBoundaryTrackStats(
            eventID,
            totalSteps,
            totalStepsBack,
            totalStepsForward,
            boundarySteps,
            boundaryStepsBack,
            boundaryStepsForward,
            NAN,
            maxConsecBoundary,
            maxConsecBoundaryBack,
            maxConsecBoundaryForward);
      }
      return;
    }
    // calculate pre-step point safety
    localPosition[0] = stop_grad(localPosition[0]);  //FIX
    localPosition[1] = stop_grad(localPosition[1]);  //FIX
    localPosition[2] = stop_grad(localPosition[2]);  //FIX
    const G4double preStepSafety  = currentVolume->DistanceToOut(localPosition);
    const G4double boundaryTol = Geometry::GetBoundaryTolerance();
    bool onBoundary = (preStepSafety <= boundaryTol);
    // get the material-cuts couple index from the volume
    const int indxMaterial = currentVolume->GetMaterialIndx();
    // set the fields needed for computing the physics step limit:
    // - material-cuts couple index and onBoundary falg
    const int hepEmIMC = theState.fData->fTheMatCutData->fG4MCIndexToHepEmMCIndex[indxMaterial];
    theTrack->SetMCIndex(hepEmIMC);
    theTrack->SetOnBoundary(onBoundary);
    //
    // Invoke the G4HepEmGammaManager to compute how far this photon goes till the next interaction
    // NOTE: 1. result of step limit will be written into `theTLData` PrimaryTrack HepEmTrack object
    //       2. the result is the straight line distance that the photon needs to travel along the current
    //          direction till the next physics interaction (assuming the same material along)
    G4HepEmGammaManager::HowFar(theState.fData, theState.fParameters, &theTLData, numStep);
    const G4double distToPhysics = theTrack->GetGStepLength();
    //
    // take the shortest from the geometry and the physics step limits as the current (straight line) step length
    G4double stepLength = distToBoundary;
    onBoundary        = true;
    if (distToPhysics < distToBoundary) {
      stepLength = distToPhysics;
      onBoundary = false;
    }

    bool hasBoundaryKey = false;
    long long boundaryKey = 0;
    int sameBoundaryHits = 0;
    int sameBoundaryFlips = 0;
    if (onBoundary) {
      const double boundaryX = BoundaryXValue(globalPosition, curDirection, stepLength);
      boundaryKey = BoundaryKeyFromX(boundaryX);
      hasBoundaryKey = true;
      int& hitCount = sameBoundaryHitCount[boundaryKey];
      ++hitCount;
      sameBoundaryHits = hitCount;
    }

    const G4double stepVx = theTrack->GetDirection()[0];
    const int vxSign = VxSign(GET_VALUE(stepVx));
    if (vxSign != 0) {
      if (hasPrevVxSign && vxSign != prevVxSign) {
        ++nFlipVx;
      }
      prevVxSign = vxSign;
      hasPrevVxSign = true;
    }
    if (onBoundary && hasBoundaryKey) {
      sameBoundaryFlips = sameBoundaryFlipCount[boundaryKey];
      if (vxSign != 0) {
        int& lastSign = sameBoundaryLastVxSign[boundaryKey];
        if (lastSign != 0 && lastSign != vxSign) {
          int& flipCount = sameBoundaryFlipCount[boundaryKey];
          ++flipCount;
          sameBoundaryFlips = flipCount;
        }
        lastSign = vxSign;
      }
    }
    const bool isBackward = stepVx < 0.0;
    const bool isGrazing = AbsValue(stepVx) < vxThreshold;
    const bool isUnsafeBackward = gEnableBackwardBoundaryStop && isBackward && onBoundary;
    const bool isUnsafeGrazing = (onBoundary || (preStepSafety < nearBoundarySafety)) && isGrazing;
    const bool reachedRepeatThreshold = onBoundary && (gSameBoundaryStopThreshold > 0) && (sameBoundaryHits >= gSameBoundaryStopThreshold);
    const bool reachedFlipThreshold = (gSameBoundaryMinFlips <= 0) || (sameBoundaryFlips >= gSameBoundaryMinFlips);
    const bool isUnsafeRepeatBoundary = reachedRepeatThreshold && reachedFlipThreshold;
    const bool isUnsafeStep = isUnsafeBackward || isUnsafeGrazing || isUnsafeRepeatBoundary;
    if (isUnsafeStep) {
      const bool fullTrackDueToBackward = isUnsafeBackward;
      const bool fullTrackDueToGrazing = isUnsafeGrazing && gGrazingStopsFullTrack;
      const bool fullTrackDueToRepeat = isUnsafeRepeatBoundary &&
        (gSameBoundaryFullTrackStop ||
         (gSameBoundaryHardStopThreshold > 0 && sameBoundaryHits >= gSameBoundaryHardStopThreshold));
      if (fullTrackDueToBackward || fullTrackDueToGrazing || fullTrackDueToRepeat) {
        DisableTrackGradient(*theTrack);
        stop_tracking = true;
      } else {
        SanitizeTrackState(*theTrack);
      }
      globalPosition = theTrack->GetPosition();
      curDirection   = theTrack->GetDirection();
      if (onBoundary) {
        // Keep boundary-limited steps boundary-limited after sanitizing/stopgrad.
        Set3Vect(localPosition, globalPosition);
        Geometry::SetDistanceDebugContext(eventID, theTrack->GetID(), theTrack->GetParentID(), numStep, "gamma_post_unsafe");
        distToBoundary = theGeometry.CalculateDistanceToOut(localPosition, curDirection, &currentVolume, &indxLayer, &indxAbs);
        Geometry::ClearDistanceDebugContext();
        Set3Vect(localPositionAtDistCalc, localPosition);
        Set3Vect(directionAtDistCalc, curDirection);
        stepLength = distToBoundary;
      }
    }

    std::ostringstream oss;
    oss.setf(std::ios::scientific);
    oss.precision(9);

    // Apply a small push if the step length is zero.
    // NOTE: it can happen that we are actually (logically) out of the volume
    //       where we located to be (due to this simplified "navigaton"). So
    //       just apply a small push to the current direction and relocate.
    if (stepLength==0.0) {
      stepLength = 1.0E-6;
      AddTo3Vect(globalPosition, curDirection, stepLength);
      continue;
    }
    // move the track to the corresponding post-step point
    AddTo3Vect(globalPosition, curDirection, stepLength);
    // update the geometrical step length (taking the selected)
    theTrack->SetGStepLength(stepLength);
    // update the `onBoundary` falg
    theTrack->SetOnBoundary(onBoundary);
    ++totalSteps;
    if (stepVx < 0.0) {
      ++totalStepsBack;
    } else {
      ++totalStepsForward;
    }
    if (onBoundary) {
      ++boundarySteps;
      ++consecBoundary;
      if (consecBoundary > maxConsecBoundary) maxConsecBoundary = consecBoundary;
      const double stepEkin = GET_VALUE(theTrack->GetEKin());
      ++blEkinCount;
      blEkinSum += stepEkin;
      blEkinSqSum += stepEkin * stepEkin;
      if (stepEkin < blEkinMin) blEkinMin = stepEkin;
      if (stepEkin > blEkinMax) blEkinMax = stepEkin;
      if (theTrack->GetDirection()[0] < 0.0) {
        ++boundaryStepsBack;
        ++consecBoundaryBack;
        consecBoundaryForward = 0;
        if (consecBoundaryBack > maxConsecBoundaryBack) maxConsecBoundaryBack = consecBoundaryBack;
      } else {
        ++boundaryStepsForward;
        ++consecBoundaryForward;
        consecBoundaryBack = 0;
        if (consecBoundaryForward > maxConsecBoundaryForward) maxConsecBoundaryForward = consecBoundaryForward;
      }
    } else {
      consecBoundary = 0;
      consecBoundaryBack = 0;
      consecBoundaryForward = 0;
    }
    // Then call `Perform` to do evything needs to be done with the track regarding physics
    // NOTE:
    //  - in case of boundary limited steps: no physics interaction just update
    //       of the `number of interaction left` based on the current step length
    //  - in case of physics limited step: interaction happens additionaly
    auto numIAleft_0_prePerform = GET_VALUE(theTrack->GetNumIALeft()[0]);
    auto numIAleft_1_prePerform = GET_VALUE(theTrack->GetNumIALeft()[1]);
    auto numIAleft_2_prePerform = GET_VALUE(theTrack->GetNumIALeft()[2]);
    #ifdef CODI_FORWARD
      auto numIAleft_0_prePerform_dot = GET_DOTVALUE(theTrack->GetNumIALeft()[0]);
      auto numIAleft_1_prePerform_dot = GET_DOTVALUE(theTrack->GetNumIALeft()[1]);
      auto numIAleft_2_prePerform_dot = GET_DOTVALUE(theTrack->GetNumIALeft()[2]);
    #endif 
    G4HepEmGammaManager::Perform(theState.fData, theState.fParameters, &theTLData);
    //if (stop_tracking) {
    //  DisableTrackGradient(*theTrack);
    // }



    #ifdef CODI_FORWARD
      const bool big_grad = IsDebugLoggingEnabled() &&
        (outputall || (std::abs(GET_DOTVALUE(stepLength)) > Ldot_thr) || (std::abs(GET_DOTVALUE(theTrack->GetEnergyDeposit())) > Edot_thr));

      if(big_grad){
        oss
          << eventID << ',' << theTrack->GetCharge() << ',' << theTrack->GetID() << ',' << theTrack->GetParentID() << ',' << creationStep << ','
          << numStep << ',' << indxLayer << ',' << theTrack->GetWinnerProcessIndex() << ',' << (stop_tracking ?1:0) << ','
          << (onBoundary ?1:0) << ',' << -1 << ','
          << pre_step_gX << ',' << pre_step_gY << ',' << pre_step_gZ << ','
          << pre_step_gX_dot << ',' << pre_step_gY_dot << ',' << pre_step_gZ_dot << ','
          << theTrack->GetPosition()[0] << ',' << theTrack->GetPosition()[1] << ',' << theTrack->GetPosition()[2] << ','
          << '-999' << ',' << '-999' << ',' << '-999' << ','
          << pre_step_vx << ',' << pre_step_vy << ',' << pre_step_vz << ','
          << theTrack->GetDirection()[0] << ',' << theTrack->GetDirection()[1] << ',' << theTrack->GetDirection()[2] << ',' << preStepSafety << ',' << '-999' << ','
          << numIAleft_0_prePerform << ',' << numIAleft_1_prePerform << ',' << numIAleft_2_prePerform << ','
          << numIAleft_0_prePerform_dot << ',' << numIAleft_1_prePerform_dot << ',' << numIAleft_2_prePerform_dot << ','
          << theTrack->GetNumIALeft()[0] << ',' << theTrack->GetNumIALeft()[1] << ',' << theTrack->GetNumIALeft()[2] << ','
          << GET_DOTVALUE(theTrack->GetNumIALeft()[0]) << ',' << GET_DOTVALUE(theTrack->GetNumIALeft()[1]) << ',' << GET_DOTVALUE(theTrack->GetNumIALeft()[2]) << ','
          << theTrack->GetMFP()[0] << ',' << theTrack->GetMFP()[1] << ',' << theTrack->GetMFP()[2] << ',' 
          << GET_DOTVALUE(theTrack->GetMFP()[0]) << ',' << GET_DOTVALUE(theTrack->GetMFP()[1]) << ',' << GET_DOTVALUE(theTrack->GetMFP()[2]) << ',' 
          << stepLength << ',' << GET_DOTVALUE(stepLength) << ','
          << '-999' << ',' << '-999' << ',' 
          << theTrack->GetEnergyDeposit() << ',' << GET_DOTVALUE(theTrack->GetEnergyDeposit()) << ','
          << theTrack->GetEKin() << ',' << GET_DOTVALUE(theTrack->GetEKin())<< ','
          << distToBoundary << ',' << GET_DOTVALUE(distToBoundary) << ',' << distToPhysics << ',' << GET_DOTVALUE(distToPhysics);
        MicroAudit::logLine(oss.str());

      }
    #endif
    //
    // Take and stack all secondaries (if any) that has been produced.
    if (theTLData.GetNumSecondaryElectronTrack() + theTLData.GetNumSecondaryGammaTrack() > 0 ) {
      StackSecondaries(theTLData, theTrackStack, *theTrack, numStep);
    }

   

    // call the SteppingAction (whenever a step was done in the calorimeter)
    SteppingAction(theResult, *theTrack, currentVolume, stepLength, indxLayer, indxAbs, eventID, numStep);
    ++numStep;
  }
  LogTrackPingPongStats(
      eventID,
      *theTrack,
      creationStep,
      totalSteps,
      totalStepsBack,
      boundarySteps,
      boundaryStepsBack,
      nFlipVx,
      maxConsecBoundary,
      maxConsecBoundaryBack,
      maxConsecBoundaryForward,
      initVx,
      blEkinCount,
      blEkinSum,
      blEkinSqSum,
      blEkinMin,
      blEkinMax);
  if (outputboundarylayers) {
    AccumulateBoundaryTrackStats(
        eventID,
        totalSteps,
        totalStepsBack,
        totalStepsForward,
        boundarySteps,
        boundaryStepsBack,
        boundaryStepsForward,
        ScalarDotValue(theTrack->GetGStepLength()),
        maxConsecBoundary,
        maxConsecBoundaryBack,
        maxConsecBoundaryForward);
  }
}


void SteppingLoop::ElectronStepper(G4HepEmTLData& theTLData, G4HepEmState& theState, TrackStack& theTrackStack, Geometry& theGeometry, Results& theResult, int eventID, G4double threshold, G4double threshold2) {
  // NOTE: the start tracking procedure (reset the track and the rng) was already done in the EventLoop
  G4HepEmTrack*           theTrack = theTLData.GetPrimaryElectronTrack()->GetTrack();
  G4HepEmMSCTrackData*  theMSCData = theTLData.GetPrimaryElectronTrack()->GetMSCTrackData();
  //
  // if this is a real primary track then I need to locate it
  // if this is a secondary then I could already know, but
  // anyway: locate in all cases to keep it simply (but slower anyway)
  //
  int  numStep       = 0;
  Box* currentVolume = nullptr;
  bool onBoundary    = false;
  int  indxLayer     = -1;
  int  indxAbs       = -1;
  G4double  localPosition[3];
  const G4double vxThreshold = AbsValue(threshold);
  const G4double nearBoundarySafety = threshold2 > 0.0 ? threshold2 : kDefaultNearBoundarySafety;
  bool wasOnBoundary = false;
//  bool wasPushed     = false;
  const int creationStep = theTrackStack.GetTrackCreationStep(theTrack->GetID());
  const double initVx = GET_VALUE(theTrack->GetDirection()[0]);
  int totalSteps = 0;
  int totalStepsBack = 0;
  int totalStepsForward = 0;
  int boundarySteps = 0;
  int boundaryStepsBack = 0;
  int boundaryStepsForward = 0;
  int consecBoundary = 0;
  int consecBoundaryBack = 0;
  int consecBoundaryForward = 0;
  int maxConsecBoundary = 0;
  int maxConsecBoundaryBack = 0;
  int maxConsecBoundaryForward = 0;
  int nFlipVx = 0;
  int prevVxSign = 0;
  bool hasPrevVxSign = false;
  int blEkinCount = 0;
  double blEkinSum = 0.0;
  double blEkinSqSum = 0.0;
  double blEkinMin = std::numeric_limits<double>::infinity();
  double blEkinMax = -std::numeric_limits<double>::infinity();
  std::unordered_map<long long, int> sameBoundaryHitCount;
  std::unordered_map<long long, int> sameBoundaryFlipCount;
  std::unordered_map<long long, int> sameBoundaryLastVxSign;

  // keep tracking while the kinetic energy drops to zero (i.e. e-/e+ lose all its energy; e+ annihilates)
  // unless the track is going out of the Calorimeter
  bool stop_tracking = IsTrackGradientDisabled(*theTrack);
  if (stop_tracking) {
    DisableTrackGradient(*theTrack);
  }
  while (theTrack->GetEKin() > 0.0) {
    if (stop_tracking) {
      DisableTrackGradient(*theTrack);
    }
    // calculate distance to boundary from the pre-step point: will locate the pont
    // NOTE: this should never be zero as zero means that the point is outside of the volume
    //       (taking into account the direction and tolerance)
    // NOTE: the given position will be in local coordiantes at return
    G4double* globalPosition = theTrack->GetPosition();
    G4double* curDirection   = theTrack->GetDirection();
    // set the local position = global position (will be local after CalculateDistanceToOut)
    Set3Vect(localPosition, globalPosition);

    G4double pre_step_gX = globalPosition[0];
    G4double pre_step_gY = globalPosition[1];
    G4double pre_step_gZ = globalPosition[2];

    #ifdef CODI_FORWARD
      G4double pre_step_gX_dot = GET_DOTVALUE(globalPosition[0]);
      G4double pre_step_gY_dot = GET_DOTVALUE(globalPosition[1]);
      G4double pre_step_gZ_dot = GET_DOTVALUE(globalPosition[2]);
    #endif
    G4double pre_step_vx = curDirection[0];
    G4double pre_step_vy = curDirection[1];
    G4double pre_step_vz = curDirection[2];

    Geometry::SetDistanceDebugContext(eventID, theTrack->GetID(), theTrack->GetParentID(), numStep, "electron_pre");
    G4double distToBoundary = theGeometry.CalculateDistanceToOut(localPosition, curDirection, &currentVolume, &indxLayer, &indxAbs);
    Geometry::ClearDistanceDebugContext();
    G4double localPositionAtDistCalc[3];
    Set3Vect(localPositionAtDistCalc, localPosition);
    G4double directionAtDistCalc[3];
    Set3Vect(directionAtDistCalc, curDirection);
    // STOP HERE IF `distToBoundary = 1.0E+20` i.e. we are going out from the Calorimeter
    if (distToBoundary > 1.0E+10) {
      LogTrackPingPongStats(
          eventID,
          *theTrack,
          creationStep,
          totalSteps,
          totalStepsBack,
          boundarySteps,
          boundaryStepsBack,
          nFlipVx,
          maxConsecBoundary,
          maxConsecBoundaryBack,
          maxConsecBoundaryForward,
          initVx,
          blEkinCount,
          blEkinSum,
          blEkinSqSum,
          blEkinMin,
          blEkinMax);
      if (outputboundarylayers){
        AccumulateBoundaryTrackStats(
            eventID,
            totalSteps,
            totalStepsBack,
            totalStepsForward,
            boundarySteps,
            boundaryStepsBack,
            boundaryStepsForward,
            NAN,
            maxConsecBoundary,
            maxConsecBoundaryBack,
            maxConsecBoundaryForward);
      } 
      return;
    }
    // at the pre-step point: calculate safety and check if on-boundary (use only if we do not know that the
    // previous step ended up on boundary i.e. use only in the very first or pushed steps)
    localPosition[0] = stop_grad(localPosition[0]);  //FIX
    localPosition[1] = stop_grad(localPosition[1]);  //FIX
    localPosition[2] = stop_grad(localPosition[2]);  //FIX
    G4double safety   = currentVolume->DistanceToOut(localPosition);
    const G4double boundaryTol = Geometry::GetBoundaryTolerance();
    const G4double onBoundaryTol = boundaryTol > 0.0 ? boundaryTol : 5.0E-10;
    bool onBoundary = numStep == 0 ? (safety <= onBoundaryTol) : wasOnBoundary;
    const G4double preStepSafety = onBoundary ? 0.0 : safety;

    // get the material-cuts couple index from the volume
    const int indxMaterial = currentVolume->GetMaterialIndx();
    // set the fields needed for computing the physics step limit:
    // - material-cuts couple index and onBoundary falg and the additional Safety for e-/e+
    const int hepEmIMC = theState.fData->fTheMatCutData->fG4MCIndexToHepEmMCIndex[indxMaterial];
    theTrack->SetMCIndex(hepEmIMC);
    theTrack->SetOnBoundary(onBoundary);
    // the additional pre-step-point safety that is used in the MSC
    theTrack->SetSafety(preStepSafety);

    //
    // Invoke the G4HepEmElectronManager to compute how far this e-/e+ goes till the next interaction
    // (that might be simply continuous step limit due to energy loss or MSC that do not produce seconday)
    // NOTE: 1. result of step limit will be written into `theTLData` PrimaryTrack HepEmTrack object
    //       2. the result is the straight line distance that the e-/e+ needs to travel along the current
    //          direction
    //       3. at the end, an additional lateral displacement might be applied (along the perpendicular plane)
    //          due to MSC
    //       4. also note, that the real length (physical) of the step is longer than the straight light along the
    //          original direction (geometrical) step length due to MSC
    G4HepEmElectronManager::HowFar(theState.fData, theState.fParameters, &theTLData);
    const G4double distToPhysics = theTrack->GetGStepLength();
    //
    // take the shortest from the geometry and physics step limits as current (straight line) step length
    // along the original direction and see if the post-step point is on-boundary
    G4double stepLength = distToBoundary;
    const bool onBoundaryBeforeSelect = onBoundary;
    onBoundary        = true;
    if (distToPhysics < distToBoundary) {
      stepLength = distToPhysics;
      onBoundary = false;
    }
    StepChoiceAudit::LogDecision(
        eventID,
        *theTrack,
        numStep,
        "pre_unsafe",
        onBoundaryBeforeSelect,
        onBoundary,
        wasOnBoundary,
        distToBoundary,
        distToPhysics,
        stepLength,
        preStepSafety,
        currentVolume,
        localPositionAtDistCalc,
        directionAtDistCalc);

    bool hasBoundaryKey = false;
    long long boundaryKey = 0;
    int sameBoundaryHits = 0;
    int sameBoundaryFlips = 0;
    if (onBoundary) {
      const double boundaryX = BoundaryXValue(globalPosition, curDirection, stepLength);
      boundaryKey = BoundaryKeyFromX(boundaryX);
      hasBoundaryKey = true;
      int& hitCount = sameBoundaryHitCount[boundaryKey];
      ++hitCount;
      sameBoundaryHits = hitCount;
    }
    
    const G4double stepVx = theTrack->GetDirection()[0];
    const int vxSign = VxSign(GET_VALUE(stepVx));
    if (vxSign != 0) {
      if (hasPrevVxSign && vxSign != prevVxSign) {
        ++nFlipVx;
      }
      prevVxSign = vxSign;
      hasPrevVxSign = true;
    }
    if (onBoundary && hasBoundaryKey) {
      sameBoundaryFlips = sameBoundaryFlipCount[boundaryKey];
      if (vxSign != 0) {
        int& lastSign = sameBoundaryLastVxSign[boundaryKey];
        if (lastSign != 0 && lastSign != vxSign) {
          int& flipCount = sameBoundaryFlipCount[boundaryKey];
          ++flipCount;
          sameBoundaryFlips = flipCount;
        }
        lastSign = vxSign;
      }
    }
    const bool isBackward = stepVx < 0.0;
    const bool isGrazing = AbsValue(stepVx) < vxThreshold;
    const bool isUnsafeBackward = gEnableBackwardBoundaryStop && isBackward && onBoundary;
    
    const bool isUnsafeGrazing = (onBoundary || (preStepSafety < nearBoundarySafety)) && isGrazing;
    const bool reachedRepeatThreshold = onBoundary && (gSameBoundaryStopThreshold > 0) && (sameBoundaryHits >= gSameBoundaryStopThreshold);
    const bool reachedFlipThreshold = (gSameBoundaryMinFlips <= 0) || (sameBoundaryFlips >= gSameBoundaryMinFlips);
    const bool isUnsafeRepeatBoundary = reachedRepeatThreshold && reachedFlipThreshold;
    const bool isUnsafeStep = isUnsafeBackward || isUnsafeGrazing || isUnsafeRepeatBoundary;
    if (isUnsafeStep) {
      const bool fullTrackDueToBackward = isUnsafeBackward;
      const bool fullTrackDueToGrazing = isUnsafeGrazing && gGrazingStopsFullTrack;
      const bool fullTrackDueToRepeat = isUnsafeRepeatBoundary &&
        (gSameBoundaryFullTrackStop ||
         (gSameBoundaryHardStopThreshold > 0 && sameBoundaryHits >= gSameBoundaryHardStopThreshold));
      if (fullTrackDueToBackward || fullTrackDueToGrazing || fullTrackDueToRepeat) {
        DisableTrackGradient(*theTrack);
        stop_tracking = true;
      } else {
        SanitizeTrackState(*theTrack);
      }
      globalPosition = theTrack->GetPosition();
      curDirection   = theTrack->GetDirection();
      if (onBoundary) {
        // Keep boundary-limited steps boundary-limited after sanitizing/stopgrad.
        Set3Vect(localPosition, globalPosition);
        Geometry::SetDistanceDebugContext(eventID, theTrack->GetID(), theTrack->GetParentID(), numStep, "electron_post_unsafe");
        distToBoundary = theGeometry.CalculateDistanceToOut(localPosition, curDirection, &currentVolume, &indxLayer, &indxAbs);
        Geometry::ClearDistanceDebugContext();
        Set3Vect(localPositionAtDistCalc, localPosition);
        Set3Vect(directionAtDistCalc, curDirection);
        stepLength = distToBoundary;
      }
      StepChoiceAudit::LogDecision(
          eventID,
          *theTrack,
          numStep,
          "post_unsafe",
          onBoundaryBeforeSelect,
          onBoundary,
          wasOnBoundary,
          distToBoundary,
          distToPhysics,
          stepLength,
          preStepSafety,
          currentVolume,
          localPositionAtDistCalc,
          directionAtDistCalc);
    }


    // Apply a small push if the step length is zero.
    // NOTE: it can happen that we are actually (logically) out of the volume
    //       where we located to be (due to this simplified "navigaton"). So
    //       just apply a small push to the current direction and relocate.
//    wasPushed = false;
    std::ostringstream oss;
    oss.setf(std::ios::scientific);
    oss.precision(9);

    if (stepLength==0.0) {
//      wasPushed  = true;
      stepLength = 1.0E-6;
      AddTo3Vect(globalPosition, curDirection, stepLength);
      continue;
    }
    // move the track to the corresponding post-step point
    AddTo3Vect(globalPosition, curDirection, stepLength);
    // update the geometrical step length (taking the selected)
    theTrack->SetGStepLength(stepLength);
    // update the `onBoundary` falg
    theTrack->SetOnBoundary(onBoundary);
    ++totalSteps;
    if (stepVx < 0.0) {
      ++totalStepsBack;
    } else {
      ++totalStepsForward;
    }
    if (onBoundary) {
      ++boundarySteps;
      ++consecBoundary;
      if (consecBoundary > maxConsecBoundary) maxConsecBoundary = consecBoundary;
      const double stepEkin = GET_VALUE(theTrack->GetEKin());
      ++blEkinCount;
      blEkinSum += stepEkin;
      blEkinSqSum += stepEkin * stepEkin;
      if (stepEkin < blEkinMin) blEkinMin = stepEkin;
      if (stepEkin > blEkinMax) blEkinMax = stepEkin;
      if (theTrack->GetDirection()[0] < 0.0) {
        ++boundaryStepsBack;
        ++consecBoundaryBack;
        consecBoundaryForward = 0;
        if (consecBoundaryBack > maxConsecBoundaryBack) maxConsecBoundaryBack = consecBoundaryBack;
      } else {
        ++boundaryStepsForward;
        ++consecBoundaryForward;
        consecBoundaryBack = 0;
        if (consecBoundaryForward > maxConsecBoundaryForward) maxConsecBoundaryForward = consecBoundaryForward;
      }
    } else {
      consecBoundary = 0;
      consecBoundaryBack = 0;
      consecBoundaryForward = 0;
    }

    // Then call `Perform` to do evything needs to be done with the track regarding physics
    //  - the continuous interactions will be performed in all cases (i.e. independently
    //    if geometry or physics limited the step):
    //    = these continuous interactions are:
    //       a. first the geometrical step is converted to physical by accounting the effects of MSC
    //       b. this real physical step length is used to compute the energy loss due to sub-threshold
    //          interactions (the mean energy loss is comuted then fluctuation is added) `
    //  - in case of continuous physics or boundary limited the step:
    //    = no further physics interaction just update of the `number of interaction left`
    //      based on the current real (i.e. physical) step length
    //  - in case of physics limited step: discrete interaction, producing seondary particle(s), happens additionaly
    // keep the original direction as it will be changed during the physics (even without discrete interaction due to MSC)
    G4double orgDirection[3];
    Set3Vect(orgDirection, curDirection);
    auto numIAleft_0_prePerform = GET_VALUE(theTrack->GetNumIALeft()[0]);
    auto numIAleft_1_prePerform = GET_VALUE(theTrack->GetNumIALeft()[1]);
    auto numIAleft_2_prePerform = GET_VALUE(theTrack->GetNumIALeft()[2]);
    #ifdef CODI_FORWARD
      auto numIAleft_0_prePerform_dot = GET_DOTVALUE(theTrack->GetNumIALeft()[0]);
      auto numIAleft_1_prePerform_dot = GET_DOTVALUE(theTrack->GetNumIALeft()[1]);
      auto numIAleft_2_prePerform_dot = GET_DOTVALUE(theTrack->GetNumIALeft()[2]);
    #endif 

    G4HepEmElectronManager::Perform(theState.fData, theState.fParameters, &theTLData);
    // take the real, i.e. physical step length (only if MSC is active in G4HepEmElectronManager because the
    // physical step length stays zero when MSC is not active as physical = geometrical in that case)
    const G4double pStepLength = theMSCData->fTrueStepLength > 0.0 ? theMSCData->fTrueStepLength : stepLength;

    // get the displacement and check if we need to apply (should not if the energy is zero but ok keep its simply)
    // we apply it if its length is lonegr than a minimum and we are not on boudnry (i.e. the current post-step point)

    G4double pre_MSC_gX = globalPosition[0];
    G4double pre_MSC_gY = globalPosition[1];
    G4double pre_MSC_gZ = globalPosition[2];
    G4double postSafety = -999;

    if (gEnableMscDisplacement && !onBoundary) {
      const G4double* displacement    = theMSCData->GetDisplacement();
      const G4double  dLength2        = displacement[0]*displacement[0] + displacement[1]*displacement[1] + displacement[2]*displacement[2];
      const G4double  kGeomMinLength  = 5.0e-8;  // 0.05 [nm]
      const G4double  kGeomMinLength2 = kGeomMinLength*kGeomMinLength; // (0.05 [nm])^2
      if (dLength2 > kGeomMinLength2) {
        // apply displacement
        // bool isPositionChanged  = true;
        const G4double dispR = std::sqrt(dLength2);
        // update local position by moving to the local longitudinal (i.e. along the original direction) post step-point
        // just to be able to compute the safety at that point
        AddTo3Vect(localPosition, orgDirection, stepLength);
        // compute the current post-step point safety and reduce a bit
        localPosition[0] = stop_grad(localPosition[0]);  //FIX
        localPosition[1] = stop_grad(localPosition[1]);  //FIX
        localPosition[2] = stop_grad(localPosition[2]);  //FIX
        postSafety = 0.99*currentVolume->DistanceToOut(localPosition);
        G4double clipSafety = postSafety;
        if (clipSafety > 0.0 && gMscDisplacementSafetyFloor > 0.0 && clipSafety < gMscDisplacementSafetyFloor) {
          clipSafety = gMscDisplacementSafetyFloor;
        }
        if (clipSafety > 0.0 && dispR < clipSafety) {
          // far away from boundary: can be applied safely i.e. we won't get to boundary
          AddTo3Vect(globalPosition, displacement);
          //near the boundary
        } else {
          // displaced point is definitely within the volume
          if (dispR < clipSafety) {
            AddTo3Vect(globalPosition, displacement);
          } else if (clipSafety > kGeomMinLength) {
            // reduced displacement
            const G4double scale = (clipSafety/dispR);
            AddTo3Vect(globalPosition, displacement, scale);
          } // else {
            // very small postSafety
            // isPositionChanged = false;
          // }
        }
      }
    }


    #ifdef CODI_FORWARD
      const bool big_grad = IsDebugLoggingEnabled() &&
        (outputall || (std::abs(GET_DOTVALUE(stepLength)) > Ldot_thr) || (std::abs(GET_DOTVALUE(theTrack->GetEnergyDeposit())) > Edot_thr));

      if(big_grad){
        oss
          << eventID << ',' << theTrack->GetCharge() << ',' << theTrack->GetID() << ',' << theTrack->GetParentID() << ',' << creationStep << ','
          << numStep << ',' << indxLayer << ',' << theTrack->GetWinnerProcessIndex() << ',' << (stop_tracking ?1:0) << ','
          << (onBoundary ?1:0) << ',' << -1 << ','
          << pre_step_gX << ',' << pre_step_gY << ',' << pre_step_gZ << ','
          << pre_step_gX_dot << ',' << pre_step_gY_dot << ',' << pre_step_gZ_dot << ','
          << pre_MSC_gX << ',' << pre_MSC_gY << ',' << pre_MSC_gZ << ','
          << theTrack->GetPosition()[0] << ',' << theTrack->GetPosition()[1] << ',' << theTrack->GetPosition()[2] << ','
          << pre_step_vx << ',' << pre_step_vy << ',' << pre_step_vz << ','
          << theTrack->GetDirection()[0] << ',' << theTrack->GetDirection()[1] << ',' << theTrack->GetDirection()[2] << ',' << preStepSafety << ',' << postSafety << ','
          << numIAleft_0_prePerform << ',' << numIAleft_1_prePerform << ',' << numIAleft_2_prePerform << ','
          << numIAleft_0_prePerform_dot << ',' << numIAleft_1_prePerform_dot << ',' << numIAleft_2_prePerform_dot << ','
          << theTrack->GetNumIALeft()[0] << ',' << theTrack->GetNumIALeft()[1] << ',' << theTrack->GetNumIALeft()[2] << ','
          << GET_DOTVALUE(theTrack->GetNumIALeft()[0]) << ',' << GET_DOTVALUE(theTrack->GetNumIALeft()[1]) << ',' << GET_DOTVALUE(theTrack->GetNumIALeft()[2]) << ','
          << theTrack->GetMFP()[0] << ',' << theTrack->GetMFP()[1] << ',' << theTrack->GetMFP()[2] << ',' 
          << GET_DOTVALUE(theTrack->GetMFP()[0]) << ',' << GET_DOTVALUE(theTrack->GetMFP()[1]) << ',' << GET_DOTVALUE(theTrack->GetMFP()[2]) << ',' 
          << stepLength << ',' << GET_DOTVALUE(stepLength) << ','
          << pStepLength << ',' << GET_DOTVALUE(pStepLength) << ','
          << theTrack->GetEnergyDeposit() << ',' << GET_DOTVALUE(theTrack->GetEnergyDeposit()) << ','
          << theTrack->GetEKin() << ',' << GET_DOTVALUE(theTrack->GetEKin())<< ','
          << distToBoundary << ',' << GET_DOTVALUE(distToBoundary) << ',' << distToPhysics << ',' << GET_DOTVALUE(distToPhysics);
        MicroAudit::logLine(oss.str());

      }
    #endif
    //
    // stack all secondaries (if any) that has been produced in this step
    if (theTLData.GetNumSecondaryElectronTrack() + theTLData.GetNumSecondaryGammaTrack() > 0 ) {
      StackSecondaries(theTLData, theTrackStack, *theTrack, numStep);
    }

    SteppingAction(theResult, *theTrack, currentVolume, pStepLength, indxLayer, indxAbs, eventID, numStep);
    wasOnBoundary = onBoundary;
    ++numStep;
  }
  LogTrackPingPongStats(
      eventID,
      *theTrack,
      creationStep,
      totalSteps,
      totalStepsBack,
      boundarySteps,
      boundaryStepsBack,
      nFlipVx,
      maxConsecBoundary,
      maxConsecBoundaryBack,
      maxConsecBoundaryForward,
      initVx,
      blEkinCount,
      blEkinSum,
      blEkinSqSum,
      blEkinMin,
      blEkinMax);
  if (outputboundarylayers) {
    AccumulateBoundaryTrackStats(
        eventID,
        totalSteps,
        totalStepsBack,
        totalStepsForward,
        boundarySteps,
        boundaryStepsBack,
        boundaryStepsForward,
        ScalarDotValue(theTrack->GetGStepLength()),
        maxConsecBoundary,
        maxConsecBoundaryBack,
        maxConsecBoundaryForward);
  }
}


void SteppingLoop::StackSecondaries(G4HepEmTLData& theTLData, TrackStack& theTrackStack, G4HepEmTrack& thePrimary, int parentStep) {
  // secondary: only possible is e-/e+ or gamma at the moemnt
  const int numSecElectron = theTLData.GetNumSecondaryElectronTrack();
  const int numSecGamma    = theTLData.GetNumSecondaryGammaTrack();
  const int numSecondaries = numSecElectron+numSecGamma;
  const bool disableSecondaries = (gGradientStopMode == 2) && IsTrackGradientDisabled(thePrimary);
  if (numSecondaries>0) {
    for (int is=0; is<numSecElectron; ++is) {
      G4HepEmTrack* secTrack = theTLData.GetSecondaryElectronTrack(is)->GetTrack();
      secTrack->SetID(theTrackStack.GetNextTrackID());
      secTrack->SetParentID(thePrimary.GetID());
      secTrack->SetPosition(thePrimary.GetPosition());
      secTrack->SetMCIndex(thePrimary.GetMCIndex());
      theTrackStack.SetTrackCreationStep(secTrack->GetID(), parentStep);
      if (disableSecondaries) {
        DisableTrackGradient(*secTrack);
      }
      theTrackStack.Copy(*secTrack, theTrackStack.Insert());
    }
    theTLData.ResetNumSecondaryElectronTrack();

    for (int is=0; is<numSecGamma; ++is) {
      G4HepEmTrack* secTrack = theTLData.GetSecondaryGammaTrack(is)->GetTrack();
      secTrack->SetID(theTrackStack.GetNextTrackID());
      secTrack->SetParentID(thePrimary.GetID());
      secTrack->SetPosition(thePrimary.GetPosition());
      secTrack->SetMCIndex(thePrimary.GetMCIndex());
      theTrackStack.SetTrackCreationStep(secTrack->GetID(), parentStep);
      if (disableSecondaries) {
        DisableTrackGradient(*secTrack);
      }
      theTrackStack.Copy(*secTrack, theTrackStack.Insert());
    }
    theTLData.ResetNumSecondaryGammaTrack();
  }
}


void SteppingLoop::SteppingAction(Results& theResult, const G4HepEmTrack& theTrack, const Box* /*currentVolume*/, G4double currentPhysStepLength, int indxLayer, int indxAbsorber, int /*eventID*/, int /*stepID*/) {
  if (indxLayer < 0) return;
  //
  G4double edep = theTrack.GetEnergyDeposit();
  G4double ke = theTrack.GetEKin();

  
  if (gEnableKECut && (ke < gKECutValue)){
    edep.setGradient(0.);
    const_cast<G4HepEmTrack&>(theTrack).SetEnergyDeposit(edep); // write-back

  }
  if (edep > 0.0) {
    theResult.fEdepPerLayer_CurrentEvent.Fill(indxLayer, edep);
    auto* dir = const_cast<G4HepEmTrack&>(theTrack).GetDirection();
    VxEdepHist::Fill(indxLayer, GET_VALUE(dir[0]), edep);
    switch (indxAbsorber) {
      case 0: theResult.fPerEventRes.fEdepAbs += edep;
              break;
      case 1: theResult.fPerEventRes.fEdepGap += edep;
              // per-layer GAP energy: accumulate only deposits in the gap region
              theResult.fEdepGapPerLayer_CurrentEvent.Fill(indxLayer, edep);
              break;
      default: //
              break;
    }
  }

  //
  if (currentPhysStepLength <= 0.0) return;
  if (theTrack.GetCharge() == 0.0) {
    theResult.fGammaTrackLenghtPerLayer.Fill(indxLayer, currentPhysStepLength);
    theResult.fPerEventRes.fNumStepsGamma += 1.0;
  } else {
    theResult.fElPosTrackLenghtPerLayer.Fill(indxLayer, currentPhysStepLength);
    theResult.fPerEventRes.fNumStepsElPos += 1.0;
  }
}


// some utilities to modify 3vectors
void SteppingLoop::Set3Vect(G4double* v, G4double to) {
  v[0] = to;
  v[1] = to;
  v[2] = to;
}

void SteppingLoop::Set3Vect(G4double* v, const G4double* to) {
  v[0] = to[0];
  v[1] = to[1];
  v[2] = to[2];
}

void SteppingLoop::AddTo3Vect(G4double* v, const G4double* u, G4double scale) {
  v[0] += scale*u[0];
  v[1] += scale*u[1];
  v[2] += scale*u[2];
}
