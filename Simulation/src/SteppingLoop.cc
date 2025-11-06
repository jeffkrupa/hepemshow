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

#include <fstream>
#include <atomic>
#include <unordered_set>

#ifndef MICRO_AUDIT_MAX
#define MICRO_AUDIT_MAX 2e6   // hard cap on lines we'll write
#endif

const double Ldot_thr   = 5e6;         // tweak to match your scale
const double Edot_thr   = 5e6;         // idem

namespace {
  G4double kDefaultKECut = 0.5;
  G4double gKECutValue = kDefaultKECut;
  bool gEnableKECut = false;
}

bool outputall = false;

struct MicroAudit {
  static std::ofstream& stream() {
    static std::ofstream s("micro_audit.csv", std::ios::out);
    static bool inited = false;
    if (!inited) {
      s << "event,charge,trackID,parentID,step,layer,winnerIdx,"
        << "onBoundary,blAxis,"
        << "gX,gY,gZ,"
        << "vx,vy,vz,"
        << "stepLength,stepLength_dot,"
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



template<typename Expr>
inline G4double stop_grad(const Expr& x) {
  //return G4double(x);
  return G4double(GET_VALUE(x));
}

static std::ofstream debugFile("electron_debug.csv");
static bool debugFileInitialized = [](){
  debugFile << "step,charge,trackID,parentID,localX,localX_dot,localY,localY_dot,localZ,localZ_dot,globalX,globalX_dot,globalY,globalY_dot,globalZ,globalZ_dot,distToBoundary,distToBoundary_dot,safety,safety_dot,preStepSafety,preStepSafety_dot,distToPhysics,distToPhysics_dot,stepLength,stepLength_dot,pStepLength,pStepLength_dot,KE,KE_dot,winnerIdx,onBoundary,wasOnBoundary,directionX,directionX_dot,directionY,directionY_dot,directionZ,directionZ_dot\n";
  return true;
}();

namespace {
  std::unordered_set<int> gDisabledTrackGradients;
  int gGradientStopMode = 2;

  inline void SanitizeTrackState(G4HepEmTrack& track) {
    G4double* position = track.GetPosition();
    track.SetPosition(stop_grad(position[0]), stop_grad(position[1]), stop_grad(position[2]));
    G4double* direction = track.GetDirection();
    track.SetDirection(stop_grad(direction[0]), stop_grad(direction[1]), stop_grad(direction[2]));
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
  SanitizeTrackState(track);
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
  bool stop_tracking = IsTrackGradientDisabled(*theTrack);
  if (stop_tracking) {
    DisableTrackGradient(*theTrack);
  }
  int  nBackScatter  = 0;  //FIX
  G4double lastDirection = theTrack->GetDirection()[0];  //FIX

  while (theTrack->GetEKin() > 0.0) {
    if (stop_tracking) {
      DisableTrackGradient(*theTrack);
    }
    if (lastDirection * theTrack->GetDirection()[0] < -1e-8) nBackScatter++;  //FIX
    if ((nBackScatter>0) || (theTrack->GetDirection()[0] < threshold && theTrack->GetDirection()[0] > threshold2)) //FIX
    {
      DisableTrackGradient(*theTrack);
      stop_tracking = true;
    }

    // calculate distance to boundary from the pre-step point: will locate the pont
    // NOTE: this should never be zero as zero means that the point is outside of the volume
    //       (taking into account the direction and tolerance)
    // NOTE: the given position will be in local coordiantes at return
    G4double* globalPosition = theTrack->GetPosition();
    G4double* curDirection   = theTrack->GetDirection();
    // set the local position = global position (will be local after CalculateDistanceToOut)
    Set3Vect(localPosition, globalPosition);
    const G4double distToBoundary = theGeometry.CalculateDistanceToOut(localPosition, curDirection, &currentVolume, &indxLayer, &indxAbs);
    // STOP HERE IF `distToBoundary = 1.0E+20` i.e. we are going out from the Calorimeter
    if (distToBoundary > 1.0E+10) {
      return;
    }
    // calculate pre-step point safety
    const G4double preStepSafety  = currentVolume->DistanceToOut(localPosition);
    bool onBoundary = (preStepSafety == 0.0);
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
    G4HepEmGammaManager::HowFar(theState.fData, theState.fParameters, &theTLData);
    const G4double distToPhysics = theTrack->GetGStepLength();
    //
    // take the shortest from the geometry and the physics step limits as the current (straight line) step length
    G4double stepLength = distToBoundary;
    onBoundary        = true;
    if (distToPhysics < distToBoundary) {
      stepLength = distToPhysics;
      onBoundary = false;
    }
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
    // Then call `Perform` to do evything needs to be done with the track regarding physics
    // NOTE:
    //  - in case of boundary limited steps: no physics interaction just update
    //       of the `number of interaction left` based on the current step length
    //  - in case of physics limited step: interaction happens additionaly
    G4HepEmGammaManager::Perform(theState.fData, theState.fParameters, &theTLData);
    //if (stop_tracking) {
    //  DisableTrackGradient(*theTrack);
    // }



    bool big_grad = outputall || (std::abs(GET_DOTVALUE(stepLength)) > Ldot_thr) || (std::abs(GET_DOTVALUE(theTrack->GetEnergyDeposit())) > Edot_thr);

    if(big_grad){
      std::ostringstream oss;
      oss.setf(std::ios::scientific);
      oss.precision(9);
      oss
        << eventID << ',' << theTrack->GetCharge() << ',' << theTrack->GetID() << ',' << theTrack->GetParentID() << ','
        << numStep << ',' << indxLayer << ',' << ((distToPhysics < distToBoundary) ? 1 : 0) << ','
        << (onBoundary ?1:0) << ',' << -1 << ','
        << theTrack->GetPosition()[0] << ',' << theTrack->GetPosition()[1] << ',' << theTrack->GetPosition()[2] << ','
        << theTrack->GetDirection()[0] << ',' << theTrack->GetDirection()[1] << ',' << theTrack->GetDirection()[2] << ','
        << stepLength << ',' << GET_DOTVALUE(stepLength) << ','
        << theTrack->GetEnergyDeposit() << ',' << GET_DOTVALUE(theTrack->GetEnergyDeposit()) << ','
        << theTrack->GetEKin() << ',' << GET_DOTVALUE(theTrack->GetEKin())<< ','
        << distToBoundary << ',' << GET_DOTVALUE(distToBoundary) << ',' << distToPhysics << ',' << GET_DOTVALUE(distToPhysics);
      MicroAudit::logLine(oss.str());

    }
    //
    // Take and stack all secondaries (if any) that has been produced.
    if (theTLData.GetNumSecondaryElectronTrack() + theTLData.GetNumSecondaryGammaTrack() > 0 ) {
      StackSecondaries(theTLData, theTrackStack, *theTrack);
    }
    // call the SteppingAction (whenever a step was done in the calorimeter)
    SteppingAction(theResult, *theTrack, currentVolume, stepLength, indxLayer, indxAbs, eventID, numStep);
    lastDirection = theTrack->GetDirection()[0];  //FIX

    ++numStep;
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
  bool wasOnBoundary = false;
  int  nBackScatter  = 0;  //FIX
  G4double lastDirection = theTrack->GetDirection()[0];  //FIX
//  bool wasPushed     = false;

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
    if (lastDirection * theTrack->GetDirection()[0] < -1e-8) nBackScatter++;  //FIX
    lastDirection = theTrack->GetDirection()[0];  //FIX
    if ((nBackScatter>0) || (theTrack->GetDirection()[0] < threshold && theTrack->GetDirection()[0] > threshold2)) //FIX
    {
      DisableTrackGradient(*theTrack);
      stop_tracking = true;
    }
    // calculate distance to boundary from the pre-step point: will locate the pont
    // NOTE: this should never be zero as zero means that the point is outside of the volume
    //       (taking into account the direction and tolerance)
    // NOTE: the given position will be in local coordiantes at return
    G4double* globalPosition = theTrack->GetPosition();
    G4double* curDirection   = theTrack->GetDirection();
    // set the local position = global position (will be local after CalculateDistanceToOut)
    Set3Vect(localPosition, globalPosition);

    const G4double distToBoundary = theGeometry.CalculateDistanceToOut(localPosition, curDirection, &currentVolume, &indxLayer, &indxAbs);
    // STOP HERE IF `distToBoundary = 1.0E+20` i.e. we are going out from the Calorimeter
    if (distToBoundary > 1.0E+10) {
      return;
    }
    // at the pre-step point: calculate safety and check if on-boundary (use only if we do not know that the
    // previous step ended up on boundary i.e. use only in the very first or pushed steps)
    localPosition[0] = stop_grad(localPosition[0]);  //FIX
    localPosition[1] = stop_grad(localPosition[1]);  //FIX
    localPosition[2] = stop_grad(localPosition[2]);  //FIX
    G4double safety   = currentVolume->DistanceToOut(localPosition);
    bool onBoundary = numStep == 0 ? (safety<5.0E-10) : wasOnBoundary;
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

    onBoundary        = true;
    if (distToPhysics < distToBoundary) {
      stepLength = distToPhysics;
      onBoundary = false;
    }
    // Apply a small push if the step length is zero.
    // NOTE: it can happen that we are actually (logically) out of the volume
    //       where we located to be (due to this simplified "navigaton"). So
    //       just apply a small push to the current direction and relocate.
//    wasPushed = false;
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
    G4HepEmElectronManager::Perform(theState.fData, theState.fParameters, &theTLData);
    // take the real, i.e. physical step length (only if MSC is active in G4HepEmElectronManager because the
    // physical step length stays zero when MSC is not active as physical = geometrical in that case)
    const G4double pStepLength = theMSCData->fTrueStepLength > 0.0 ? theMSCData->fTrueStepLength : stepLength;

    // get the displacement and check if we need to apply (should not if the energy is zero but ok keep its simply)
    // we apply it if its length is lonegr than a minimum and we are not on boudnry (i.e. the current post-step point)


    bool big_grad = outputall || (std::abs(GET_DOTVALUE(stepLength)) > Ldot_thr) || (std::abs(GET_DOTVALUE(theTrack->GetEnergyDeposit())) > Edot_thr);

    if(big_grad){
      std::ostringstream oss;
      oss.setf(std::ios::scientific);
      oss.precision(9);
      oss
        << eventID << ',' << theTrack->GetCharge() << ',' << theTrack->GetID() << ',' << theTrack->GetParentID() << ','
        << numStep << ',' << indxLayer << ',' << ((distToPhysics < distToBoundary) ? 1 : 0) << ','
        << (onBoundary ?1:0) << ',' << -1 << ','
        << theTrack->GetPosition()[0] << ',' << theTrack->GetPosition()[1] << ',' << theTrack->GetPosition()[2] << ','
        << theTrack->GetDirection()[0] << ',' << theTrack->GetDirection()[1] << ',' << theTrack->GetDirection()[2] << ','
        << stepLength << ',' << GET_DOTVALUE(stepLength) << ','
        << theTrack->GetEnergyDeposit() << ',' << GET_DOTVALUE(theTrack->GetEnergyDeposit()) << ','
        << theTrack->GetEKin() << ',' << GET_DOTVALUE(theTrack->GetEKin())<< ','
        << distToBoundary << ',' << GET_DOTVALUE(distToBoundary) << ',' << distToPhysics << ',' << GET_DOTVALUE(distToPhysics);
      MicroAudit::logLine(oss.str());

    }

    if (!onBoundary) {
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
        const G4double postSafety = 0.99*currentVolume->DistanceToOut(localPosition);
        if (postSafety > 0.0 && dispR < postSafety) {
          // far away from boundary: can be applied safely i.e. we won't get to boundary
          AddTo3Vect(globalPosition, displacement);
          //near the boundary
        } else {
          // displaced point is definitely within the volume
          if (dispR < postSafety) {
            AddTo3Vect(globalPosition, displacement);
          } else if(postSafety > kGeomMinLength) {
            // reduced displacement
            const G4double scale = (postSafety/dispR);
            AddTo3Vect(globalPosition, displacement, scale);
          } // else {
            // very small postSafety
            // isPositionChanged = false;
          // }
        }
      }
    }
    //
    // stack all secondaries (if any) that has been produced in this step
    if (theTLData.GetNumSecondaryElectronTrack() + theTLData.GetNumSecondaryGammaTrack() > 0 ) {
      StackSecondaries(theTLData, theTrackStack, *theTrack);
    }

    SteppingAction(theResult, *theTrack, currentVolume, pStepLength, indxLayer, indxAbs, eventID, numStep);
    wasOnBoundary = onBoundary;
    lastDirection = theTrack->GetDirection()[0];  //FIX
    ++numStep;
  }
}


void SteppingLoop::StackSecondaries(G4HepEmTLData& theTLData, TrackStack& theTrackStack, G4HepEmTrack& thePrimary) {
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
    switch (indxAbsorber) {
      case 0: theResult.fPerEventRes.fEdepAbs += edep;
              break;
      case 1: theResult.fPerEventRes.fEdepGap += edep;
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
