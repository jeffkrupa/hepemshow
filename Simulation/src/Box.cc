#include "ad_type.h"


#include "Box.hh"


#include <iostream>
#include <sstream>
#include <cmath>

namespace {

G4double& BoxDistanceDirDenFloor() {
  static G4double floor = static_cast<G4double>(0.0);
  return floor;
}

inline G4double BoxStopGrad(const G4double& x) {
  return static_cast<G4double>(GET_VALUE(x));
}

inline G4double BoxSignedFloorFromValue(const G4double& val, const G4double& floor) {
  if (!(GET_VALUE(floor) > 0.0)) {
    return val;
  }
  const double absVal = std::abs(GET_VALUE(val));
  if (absVal >= GET_VALUE(floor)) {
    return val;
  }
  return GET_VALUE(val) < 0.0 ? -floor : floor;
}

inline G4double BoxRegularizedRatioKeepPrimal(const G4double& num, const G4double& den, const G4double& denFloor) {
  if (!(GET_VALUE(denFloor) > 0.0)) {
    return num / den;
  }
  const G4double numVal = BoxStopGrad(num);
  const G4double denVal = BoxStopGrad(den);
  const G4double denReg = BoxSignedFloorFromValue(denVal, denFloor);
  const G4double primal = BoxStopGrad(num / den);
  const G4double dfdn = BoxStopGrad(1.0 / denReg);
  const G4double dfdd = BoxStopGrad(-numVal / (denReg * denReg));
  const G4double corr = (num - numVal) * dfdn + (den - denVal) * dfdd;
  return primal + (corr - BoxStopGrad(corr));
}

}  // namespace


Box::Box (const std::string& name, int indxMat, G4double pX, G4double pY, G4double pZ)
: fName(name),
  fMaterialIndx(indxMat),
  fDx(pX),
  fDy(pY),
  fDz(pZ) {
  fDelta = 0.5*kCarTolerance;
  // check minimum size
//  if (pX < 2*kCarTolerance ||
//      pY < 2*kCarTolerance ||
//      pZ < 2*kCarTolerance)  {
//    std::ostringstream message;
//    message << "Dimensions too small for Solid: " << GetName() << "!" << std::endl
//            << "     hX, hY, hZ = " << pX << ", " << pY << ", " << pZ;
//    std::cout << message.str();
//  }
}


void Box::SetHalfLength(G4double val, int idx) {
   // limit to thickness of surfaces
  if (val > 2*kCarTolerance) {
    switch (idx) {
      case 0: fDx = val;
            break;
      case 1: fDy = val;
              break;
      case 2: fDz = val;
              break;
    };
  } else {
//    std::ostringstream message;
//    message << "Dimension too small for solid: " << GetName() << "!"
//            << std::endl
//            << "      val = " << val << std::endl
//            << "      idx = " << idx;
//    std::cout << message.str();
  }
}

G4double Box::GetHalfLength(int idx) const {
  switch (idx) {
    case 0: return fDx;
    case 1: return fDy;
    case 2: return fDz;
  };
  return 0;
}

void Box::ConfigureDistanceToOutDerivativeRegularization(G4double dirDenFloor) {
  BoxDistanceDirDenFloor() = GET_VALUE(dirDenFloor) > 0.0 ? dirDenFloor : static_cast<G4double>(0.0);
}

G4double Box::GetDistanceToOutDerivativeRegularization() {
  return BoxDistanceDirDenFloor();
}


// p should be in local coordinates
// returns zero if p is outside of the box or within tolerance
G4double Box::DistanceToOut(G4double* p, G4double *v) const {
  // Check if point is not inside and traveling away: zero
  // Note: eitehr in surafece or outside
  if ((std::abs(p[0]) - fDx) >= -fDelta && p[0]*v[0] > 0) {
    return 0.0;
  }
  if ((std::abs(p[1]) - fDy) >= -fDelta && p[1]*v[1] > 0) {
    return 0.0;
  }
  if ((std::abs(p[2]) - fDz) >= -fDelta && p[2]*v[2] > 0) {
    return 0.0;
  }
  // Find intersection
  //
  const G4double vx = v[0];
  const G4double dirDenFloor = BoxDistanceDirDenFloor();
  const G4double txNum = static_cast<G4double>(std::copysign(GET_VALUE(fDx), GET_VALUE(vx))) - p[0];
  const G4double tx = (GET_VALUE(vx) == 0.0) ? 1.0E+20 : BoxRegularizedRatioKeepPrimal(txNum, vx, dirDenFloor);
  //
  const G4double vy = v[1];
  const G4double tyNum = static_cast<G4double>(std::copysign(GET_VALUE(fDy), GET_VALUE(vy))) - p[1];
  const G4double ty = (GET_VALUE(vy) == 0.0) ? tx : BoxRegularizedRatioKeepPrimal(tyNum, vy, dirDenFloor);
  const G4double txy = std::min(tx,ty);
  //
  const G4double vz = v[2];
  const G4double tzNum = static_cast<G4double>(std::copysign(GET_VALUE(fDz), GET_VALUE(vz))) - p[2];
  const G4double tz = (GET_VALUE(vz) == 0.0) ? txy : BoxRegularizedRatioKeepPrimal(tzNum, vz, dirDenFloor);
  const G4double tmax = std::min(txy,tz);
  //
  return tmax;
}


G4double Box::DistanceToOut(G4double* p) const {
  G4double dist = std::min( std::min(
                   fDx-std::abs(p[0]),
                   fDy-std::abs(p[1]) ),
                   fDz-std::abs(p[2]) );
  return (dist > 0) ? dist : 0.0;
}


/*
EInside Box::Inside(G4double rx, G4double ry, G4double rz) const {
  G4double dist = std::max ( std::max (
                  std::abs(rx)-fDx,
                  std::abs(ry)-fDy),
                  std::abs(rz)-fDz);
  return (dist > fDelta) ? kOutside : ((dist > -fDelta) ? kSurface : kInside);
}
EInside Box::Inside(G4double* r) const {
  G4double dist = std::max ( std::max (
                  std::abs(r[0])-fDx,
                  std::abs(r[1])-fDy),
                  std::abs(r[2])-fDz);
  return (dist > fDelta) ? kOutside : ((dist > -fDelta) ? kSurface : kInside);
}
*/
