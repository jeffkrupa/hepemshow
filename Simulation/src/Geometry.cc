#include "ad_type.h"


#include "Geometry.hh"

#include "Box.hh"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <string>
#include <limits>

template<typename Expr>
inline G4double stop_grad(const Expr& x) {
  //return G4double(x);
  return G4double(GET_VALUE(x));
}
G4double fCaloOffsetX;
namespace {
  G4double gBoundaryTolerance = 0.0;

  struct GeometryDebugContext {
    bool active = false;
    int eventID = -1;
    int trackID = -1;
    int parentID = -1;
    int stepID = -1;
    std::string stage;
  };

  GeometryDebugContext& GetGeometryDebugContext() {
    static thread_local GeometryDebugContext ctx;
    return ctx;
  }

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

  const char* ReadEnvStr(const char* key, const char* defaultValue) {
    const char* raw = std::getenv(key);
    if (raw == nullptr || *raw == '\0') {
      return defaultValue;
    }
    return raw;
  }

  int GeomDebugTrackId() {
    static const int kTrackId = ReadEnvInt("HEPEMSHOW_GEOM_DEBUG_TRACK_ID", -1);
    return kTrackId;
  }

  bool IsGeometryDebugEnabled() {
    static const int kEnabled = ReadEnvInt("HEPEMSHOW_ENABLE_DEBUG_LOGS", 0);
    return kEnabled != 0;
  }

  bool ShouldLogGeometry() {
    if (!IsGeometryDebugEnabled()) {
      return false;
    }
    const GeometryDebugContext& ctx = GetGeometryDebugContext();
    const int target = GeomDebugTrackId();
    return ctx.active && target >= 0 && ctx.trackID == target;
  }

  inline double ValueOf(const G4double& v) {
    return static_cast<double>(GET_VALUE(v));
  }

#ifdef CODI_FORWARD
  inline double DotOf(const G4double& v) {
    return static_cast<double>(GET_DOTVALUE(v));
  }
#else
  inline double DotOf(const G4double&) {
    return 0.0;
  }
#endif

  G4double DebugNaN() {
    return static_cast<G4double>(std::numeric_limits<double>::quiet_NaN());
  }

  std::ofstream& GeomDebugStream() {
    static std::ofstream out;
    static bool initialised = false;
    if (!initialised) {
      const char* filePath = ReadEnvStr("HEPEMSHOW_GEOM_DEBUG_FILE", "geom_distance_debug.csv");
      out.open(filePath, std::ios::out | std::ios::trunc);
      if (out.good()) {
        out << "event,trackID,parentID,step,ctxStage,geomStage,reason,volume,layer,abs,"
            << "rX,rX_dot,rY,rY_dot,rZ,rZ_dot,"
            << "vX,vX_dot,vY,vY_dot,vZ,vZ_dot,"
            << "rxCalo,rxCalo_dot,trLayer,trLayer_dot,rxLayer,rxLayer_dot,"
            << "dToCalo,dToCalo_dot,distOut,distOut_dot\n";
      }
      initialised = true;
    }
    return out;
  }

  void LogGeomDistanceStage(const char* geomStage, const char* reason,
                            const Box* volume, int indxLayer, int indxAbs,
                            const G4double* r, const G4double* v,
                            const G4double& rxCalo, const G4double& trLayeri, const G4double& rxLayer,
                            const G4double& dToCalo, const G4double& distOut) {
    if (!ShouldLogGeometry()) {
      return;
    }
    std::ofstream& out = GeomDebugStream();
    if (!out.good()) {
      return;
    }
    const GeometryDebugContext& ctx = GetGeometryDebugContext();
    const std::string volumeName = volume == nullptr ? "null" : volume->GetName();
    out
      << ctx.eventID << ','
      << ctx.trackID << ','
      << ctx.parentID << ','
      << ctx.stepID << ','
      << ctx.stage << ','
      << geomStage << ','
      << reason << ','
      << volumeName << ','
      << indxLayer << ','
      << indxAbs << ','
      << ValueOf(r[0]) << ',' << DotOf(r[0]) << ','
      << ValueOf(r[1]) << ',' << DotOf(r[1]) << ','
      << ValueOf(r[2]) << ',' << DotOf(r[2]) << ','
      << ValueOf(v[0]) << ',' << DotOf(v[0]) << ','
      << ValueOf(v[1]) << ',' << DotOf(v[1]) << ','
      << ValueOf(v[2]) << ',' << DotOf(v[2]) << ','
      << ValueOf(rxCalo) << ',' << DotOf(rxCalo) << ','
      << ValueOf(trLayeri) << ',' << DotOf(trLayeri) << ','
      << ValueOf(rxLayer) << ',' << DotOf(rxLayer) << ','
      << ValueOf(dToCalo) << ',' << DotOf(dToCalo) << ','
      << ValueOf(distOut) << ',' << DotOf(distOut)
      << '\n';
  }
}

void Geometry::SetBoundaryTolerance(G4double tol) {
  gBoundaryTolerance = std::max<G4double>(0.0, tol);
}

G4double Geometry::GetBoundaryTolerance() {
  return gBoundaryTolerance;
}

void Geometry::SetDistanceDebugContext(int eventID, int trackID, int parentID, int stepID, const char* stage) {
  GeometryDebugContext& ctx = GetGeometryDebugContext();
  ctx.active = true;
  ctx.eventID = eventID;
  ctx.trackID = trackID;
  ctx.parentID = parentID;
  ctx.stepID = stepID;
  ctx.stage = stage == nullptr ? "" : stage;
}

void Geometry::ClearDistanceDebugContext() {
  GetGeometryDebugContext() = GeometryDebugContext{};
}

Geometry::Geometry() {
  // default values: 50 layers of 2.3 [mm] absorber (PbWO4) and 5.7 [mm] gap (lAr)
  fNumLayers  =  50;
  fAbsThick   = 2.3; // defult value [mm]
  fGapThick   = 5.7; // defult value [mm]
  fCaloSizeYZ = 400; // defult value [mm]

  // these will be computed automatically in the `UpdateParameters`
  fCaloStartX       = 0.0;
  fPrimaryXPosition = 0.0;

  // crate shapes here for all objects:
  // - their proper size is set when calling `UpdateParameters` below
  // - material index is set to 0, 1 or 2 that corresponds to (using the default
  //   material list) {0 - G4_Galactic; 1 - G4_PbWO4, 2 - G4_lAr}
  fBoxWorld  = new Box("World", 0, 1,1,1); // material index = 0 G4_Galactic
  fBoxCalo   = new Box("Calo" , 0, 1,1,1); // material index = 0 G4_Galactic
  fBoxLayer  = new Box("Layer", 0, 1,1,1); // material index = 0 G4_Galactic
  fBoxAbs    = new Box("Abs"  , 1, 1,1,1); // material index = 1 G4_PbWO4
  fBoxGap    = new Box("Gap"  , 2, 1,1,1); // material index = 2 G4_lAr

  UpdateParameters();
}


Geometry::~Geometry() {
  delete fBoxWorld;
  delete fBoxCalo;
  delete fBoxLayer;
  delete fBoxGap;
  delete fBoxAbs;
}


void Geometry::UpdateParameters() {
  // calculate the layer and calorimeter thicknesses based on the `absorber`,
  // `gap` thinkesses and the number of layers
  fLayerThick = fAbsThick + fGapThick;
  fCaloThick  = fNumLayers*fLayerThick;

  // set/calculate the left hand side x point where the calorimeter starts
  //fCaloStartX = -0.5*fCaloThick;
  fCaloOffsetX = 0.5 * fCaloThick;  //FIX
  fCaloStartX = 0.0;  //FIX
  // set/calculate a world size such that everything fits inside
  const G4double worldThick = 1.1*fCaloThick;
  // set/calculate the mid-point between the world and calorimeter on the left
  //fPrimaryXPosition = -0.25*(worldThick + fCaloThick);
  fPrimaryXPosition = GET_VALUE(-0.25*(worldThick + fCaloThick) + fCaloOffsetX);  //FIX

  // half size of all (but the world) along the YZ plane
  const G4double halfCaloYZ = 0.5*fCaloSizeYZ;

  fBoxWorld->SetHalfLength(0.5*worldThick, 0);
  fBoxWorld->SetHalfLength(1.1*halfCaloYZ, 1);
  fBoxWorld->SetHalfLength(1.1*halfCaloYZ, 2);

  int iMatCalo = 0;
  fBoxCalo->SetMaterialIndx(iMatCalo);
  fBoxCalo->SetHalfLength(0.5*fCaloThick, 0);
  fBoxCalo->SetHalfLength(halfCaloYZ, 1);
  fBoxCalo->SetHalfLength(halfCaloYZ, 2);


  fBoxLayer->SetHalfLength(0.5*fLayerThick, 0);
  fBoxLayer->SetHalfLength(halfCaloYZ, 1);
  fBoxLayer->SetHalfLength(halfCaloYZ, 2);

  fBoxAbs->SetHalfLength(0.5*fAbsThick, 0);
  fBoxAbs->SetHalfLength(halfCaloYZ, 1);
  fBoxAbs->SetHalfLength(halfCaloYZ, 2);

  fBoxGap->SetHalfLength(0.5*fGapThick, 0);
  fBoxGap->SetHalfLength(halfCaloYZ, 1);
  fBoxGap->SetHalfLength(halfCaloYZ, 2);
}


// note: try to keep this more verbose than fast to keep it clear
G4double Geometry::CalculateDistanceToOut(G4double* r, G4double *v, Box** currentVolume, int* indxLayer, int* indxAbs) {
  const G4double boundaryTol = gBoundaryTolerance;
  const G4double dbgNaN      = DebugNaN();
  // init everything to a step in the `world` case
  *currentVolume = fBoxWorld;
  *indxLayer     = -1;
  *indxAbs       = -1;
  LogGeomDistanceStage("enter_global", "start", *currentVolume, *indxLayer, *indxAbs, r, v, dbgNaN, dbgNaN, dbgNaN, dbgNaN, dbgNaN);

  // calculate position in the `calorimeter` system:
  // - only x-coordinate is need as everything is centered along the yz
  // - actually its the same as the global: the calorimeter is not translated nor rotated
  const G4double rx_Calo = r[0];
  r[0] = r[0] - 0.5 * fCaloThick; //FIX
  const G4double dToCalo = fBoxCalo->DistanceToOut(r, v);
  LogGeomDistanceStage("after_calo_shift", "computed_dToCalo", fBoxCalo, *indxLayer, *indxAbs, r, v, rx_Calo, dbgNaN, dbgNaN, dToCalo, dbgNaN);
  // check if about leaving the calorimeter volume: distance to out is zero
  if (dToCalo <= boundaryTol) {
    // currentVolume is already set to `world`
    const G4double distOut = 1.0E+20;
    LogGeomDistanceStage("return_world", "dToCalo_le_boundaryTol", *currentVolume, *indxLayer, *indxAbs, r, v, rx_Calo, dbgNaN, dbgNaN, dToCalo, distOut);
    return 1.0E+20;
  }

  // calculate the position in the `layer` system:
  // - first calculate the index of the `layer` in which the point is located
  const int iLayer = int( GET_VALUE(((rx_Calo)/fLayerThick)) );  //FIX +0.5*fCaloThick
  *indxLayer = iLayer;
  // - then the corresponding translation vector and transform the point
  const G4double trLayeri = (iLayer)*fLayerThick; //FIX -0.5*fCaloThick +  +0.5
  const G4double rx_Layer = rx_Calo - trLayeri;
  r[0] = rx_Layer - 0.5*fLayerThick; //FIX

  // calculate the distance to the `layer` boundary along the given direction
  // why: tolerance and direction was not considered! So to detect here that
  //      the point is actually miss-located (distance is zero in that case.)
  const G4double dToLayer = fBoxLayer->DistanceToOut(r, v);
  LogGeomDistanceStage("after_layer_map", "computed_dToLayer", fBoxLayer, *indxLayer, *indxAbs, r, v, rx_Calo, trLayeri, rx_Layer, dToCalo, dToLayer);
  if (dToLayer <= boundaryTol) {
    const G4double distOut = 0.0;
    LogGeomDistanceStage("return_layer_boundary", "dToLayer_le_boundaryTol", fBoxLayer, *indxLayer, *indxAbs, r, v, rx_Calo, trLayeri, rx_Layer, dToCalo, distOut);
    return 0.0;
    // NOTE: I could also push here and do recursion but keep it clear and push only in the steppers
  }

  // calculate if the point is in the `absorber` or the `gap` part of the `layer`
  if (rx_Layer < fAbsThick || fGapThick == 0) { // in the `absorber`
    // calculate the position in the `absorber` system:
    // - the translation vector and transform the point
    const G4double trAbs = 0.5 * fAbsThick;  //FIX -0.5*(fLayerThick - fAbsThick)
    r[0] = rx_Layer - trAbs;
    // set what is left and calculate the distance to the `absorber` boundary along
    // the given direction (again, I could push here and do recursion whenever it's zero)
    *currentVolume = fBoxAbs;
    *indxAbs       = 0;
    const G4double distOut = fBoxAbs->DistanceToOut(r, v);
    LogGeomDistanceStage("return_abs", "inside_absorber", *currentVolume, *indxLayer, *indxAbs, r, v, rx_Calo, trLayeri, rx_Layer, dToCalo, distOut);
    return distOut;
  } else { // in the `gap`
    // calculate the position in the `gap` system:
    // - the translation vector and transform the point
    const G4double  trGap = fAbsThick + 0.5 * fGapThick; //FIX -0.5*(fLayerThick -fGapThick) + 
    r[0] = rx_Layer - trGap;
    // set what is left and calculate the distance to the `gap` boundary along
    // the given direction (again, I could push here and do recursion whenever it's zero)
    *currentVolume = fBoxGap;
    *indxAbs       = 1;
    const G4double distOut = fBoxGap->DistanceToOut(r, v);
    LogGeomDistanceStage("return_gap", "inside_gap", *currentVolume, *indxLayer, *indxAbs, r, v, rx_Calo, trLayeri, rx_Layer, dToCalo, distOut);
    return distOut;
  }
}
