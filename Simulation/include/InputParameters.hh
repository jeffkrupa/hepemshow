#include "ad_type.h"


#ifndef INPUTPARAMETERS_HH
#define INPUTPARAMETERS_HH

/**
 * @file    InputParameters.hh
 * @struct  InputParameters
 * @author  M. Novak
 * @date    Aug 2023
 *
 * @brief A data structure that encapsulates all the possible input arguments of the `HepEmShow` application.
 */

#include <iostream>
#include <string>
#include <vector>

// NOTE: this is Unix specific!
#include <getopt.h>


struct InputParameters {

  /** CTR with default values: default geometry, primary and event configuirations (see below) with
    * pre-generated data files expected at `../data/hepem_data` relative to the `HepEmShow` executable.*/
  InputParameters()
  : fG4HepEmDataFile("../data/hepem_data"),
    fRunVerbosity(1),
    fThreshold(0.1),
    fThreshold2(-1.),
    fGrazingStopsTrack(true),
    fBackwardBoundaryStop(false),
    fGradientStopMode(2),
    fEnableKECut(false),
    fKECut(0.5),
    fEnableMscDisplacement(true),
    fEnableMscStepRandomization(true),
    fBoundaryTolerance(0.0),
    fMscDisplacementSafetyFloor(0.0),
    fSameBoundaryStop(0),
    fSameBoundaryPosTolerance(1.0E-6),
    fSameBoundaryMinFlips(1),
    fSameBoundaryFullTrackStop(false),
    fSameBoundaryHardStop(0),
    fNumIALeftMfpFloor(0.0),
    fGammaNumIALeftMfpFloor(0.0),
    fGammaPhotoelectricEkinFloor(0.0),
    fBoxDirDenFloor(0.0),
    fRotateUpDerivativeFloor(0.0),
    fConversionDerivativeEpsilon(0.0),
    fUmscCosThetaDenFloor(0.0),
    fUmscTauBlendEpsilon(0.0),
    fUmscSimpleDenFloor(0.0),
    fUmscDispRadFloor(0.0) {}


  /** The geometry related input arguments.*/
  struct Geometry {
    /** CTR with default values: 50 layers of 2.3 mm absorber and 5.7 mm gap with 400 mm transvers size.*/
    Geometry()
    : fNumLayers(50),
      fThicknessAbsorber(2.3),
      fThicknessGap(5.7),
      fThicknessCalo(0),
      fSizeTransverse(400.0) {}

    int    fNumLayers;         ///< number of layers in the calorimeter
    G4double fThicknessAbsorber; ///< absorber thickness along X in [mm]
    G4double fThicknessGap;      ///< gap thickness along X in [mm]
    G4double fThicknessCalo;     ///< calorimeter thickness along X [mm] ONLY if number of layers is zero
    G4double fSizeTransverse;    ///< calorimeter full size along YZ in [mm]
  };


  /** The primary partcile and events related input arguments. */
  struct PrimaryAndEvents {
    /**CTR with default values: simulate 1000 events, starting with an electron of 10 GeV each (do not report progress).*/
    PrimaryAndEvents( )
    : fParticleName("e-"),
      fParticleEnergy(10000.0),
      fNumEvents(1000),
      fRandomSeed(1234) {}

    std::string  fParticleName;   ///< primary particle name: {"e-", "e+" or "gamma"}
    G4double       fParticleEnergy; ///< primary particle energy in [MeV]
    int          fNumEvents;      ///< number of events to simulate (each will start with a single primary)
    G4double       fRandomSeed;     ///< seed for the random number generator
  };

  // all members
  Geometry         fGeometry;         ///< the geometry related configuration
  PrimaryAndEvents fPrimaryAndEvents; ///< the primary partcile and events related configuration
  std::string      fG4HepEmDataFile;  ///< the pre-generated data file (with path)
  int              fRunVerbosity;     ///< level of printout verbosity duing setting up: nothing when < 1.
  G4double fThreshold;        // FIX
  G4double fThreshold2;       // FIX
  bool fGrazingStopsTrack; ///< when true, unsafe grazing steps disable gradients for the rest of the track
  bool fBackwardBoundaryStop; ///< when true, boundary-limited backward (vx<0) steps are treated as unsafe for stop-grad logic
  int fGradientStopMode;   ///< 0: keep gradients, 1: stop on current track, 2: stop on track and descendants
  bool fEnableKECut;       ///< apply kinetic-energy cut when true
  G4double fKECut;         ///< threshold value for kinetic-energy cut in [MeV]
  bool fEnableMscDisplacement;     ///< apply MSC lateral displacement in transport when true
  bool fEnableMscStepRandomization; ///< randomize UMSC step limit with Gaussian smearing when true
  G4double fBoundaryTolerance;      ///< boundary tolerance [mm] for zero-distance checks
  G4double fMscDisplacementSafetyFloor; ///< floor [mm] applied to displacement clipping safety
  int fSameBoundaryStop;            ///< total same-boundary hits over track history that trigger stop-grad (0 disables)
  G4double fSameBoundaryPosTolerance; ///< tolerance [mm] for matching repeated boundary x-planes
  int fSameBoundaryMinFlips;        ///< minimum number of vx sign flips observed on the same boundary over track history
  bool fSameBoundaryFullTrackStop;  ///< if true, same-boundary trigger disables full track instead of step-local sanitize
  int fSameBoundaryHardStop;        ///< hard full-track stop threshold on total hits of the same boundary (0 disables)
  G4double fNumIALeftMfpFloor;      ///< derivative-only mfp floor [mm] in UpdateNumIALeft (0 disables)
  G4double fGammaNumIALeftMfpFloor; ///< derivative-only mfp floor [mm] in Gamma UpdateNumIALeft (0 disables)
  G4double fGammaPhotoelectricEkinFloor; ///< derivative-only ekin floor [MeV] for 1/ekin in gamma photoelectric xsec (0 disables)
  G4double fBoxDirDenFloor;         ///< derivative-only signed floor for Box::DistanceToOut direction denominators vx/vy/vz
  G4double fRotateUpDerivativeFloor; ///< derivative-only floor for RotateToReferenceFrame denominator sqrt(refDir_x^2+refDir_y^2)
  G4double fConversionDerivativeEpsilon; ///< derivative-only epsilon for near-singular MSC true/geom conversion formulas
  G4double fUmscCosThetaDenFloor;   ///< derivative-only denominator floor in UMSC SampleCosineTheta ratio chain
  G4double fUmscTauBlendEpsilon;    ///< derivative-only smoothing width around UMSC tau branch threshold
  G4double fUmscSimpleDenFloor;     ///< derivative-only denominator floor in UMSC SimpleScattering
  G4double fUmscDispRadFloor;       ///< derivative-only floor for UMSC displacement sqrt radicand derivative
  #ifdef CODI_REVERSE
    std::vector<double> barEdep;     ///< Bar values of the energy depositions
  #endif
};


void PrintParameters (const struct InputParameters& theParam) {

  std::cout << " \n === HepEmShow input parameters: "    << std::endl;
  std::cout << "     --- Geometry configuration: " << std::endl;
  std::cout << "         - number-of-layers      : "     << theParam.fGeometry.fNumLayers         << std::endl;
  std::cout << "         - absorber-thickness    : "     << theParam.fGeometry.fThicknessAbsorber << " [mm]" << std::endl;
  std::cout << "         - gap-thickness         : "     << theParam.fGeometry.fThicknessGap      << " [mm]" << std::endl;
  std::cout << "         - transverse-size       : "     << theParam.fGeometry.fSizeTransverse    << " [mm]" << std::endl;

  std::cout << "     --- Primary and Event configuration: " << std::endl;
  std::cout << "         - primary-particle      : "     << theParam.fPrimaryAndEvents.fParticleName   << std::endl;
  std::cout << "         - primary-energy        : "     << theParam.fPrimaryAndEvents.fParticleEnergy << " [MeV]" << std::endl;
  std::cout << "         - number-of-events      : "     << theParam.fPrimaryAndEvents.fNumEvents      <<  std::endl;
  std::cout << "         - random-seed           : "     << theParam.fPrimaryAndEvents.fRandomSeed     <<  std::endl;

  std::cout << "     --- Additional configuration: " << std::endl;
  std::cout << "         - g4hepem-data-file    : "     << theParam.fG4HepEmDataFile  << std::endl;
  std::cout << "         - run-verbosity        : "     << theParam.fRunVerbosity     << std::endl;
  std::cout << "         - threshold            : "     << theParam.fThreshold        << " (|vx| grazing threshold)" << std::endl;
  std::cout << "         - threshold2           : "     << theParam.fThreshold2       << " (near-boundary safety [mm], <=0 uses default)" << std::endl;
  std::cout << "         - grazing-stop-track   : "     << (theParam.fGrazingStopsTrack ? 1 : 0) << std::endl;
  std::cout << "         - backward-boundary-stop: "    << (theParam.fBackwardBoundaryStop ? 1 : 0) << std::endl;
  std::cout << "         - stop-grad-mode       : "     << theParam.fGradientStopMode << std::endl;
  std::cout << "         - msc-displacement     : "     << (theParam.fEnableMscDisplacement ? 1 : 0) << std::endl;
  std::cout << "         - msc-step-random      : "     << (theParam.fEnableMscStepRandomization ? 1 : 0) << std::endl;
  std::cout << "         - boundary-tolerance   : "     << theParam.fBoundaryTolerance << " [mm]" << std::endl;
  std::cout << "         - msc-disp-safe-floor  : "     << theParam.fMscDisplacementSafetyFloor << " [mm]" << std::endl;
  std::cout << "         - same-boundary-stop   : "     << theParam.fSameBoundaryStop << " (total hits over track history, 0=off)" << std::endl;
  std::cout << "         - same-boundary-pos-tol: "     << theParam.fSameBoundaryPosTolerance << " [mm]" << std::endl;
  std::cout << "         - same-boundary-min-flips: "   << theParam.fSameBoundaryMinFlips << std::endl;
  std::cout << "         - same-boundary-full-track: "  << (theParam.fSameBoundaryFullTrackStop ? 1 : 0) << std::endl;
  std::cout << "         - same-boundary-hard-stop: "   << theParam.fSameBoundaryHardStop << " (0=off)" << std::endl;
  std::cout << "         - numia-mfp-floor      : "     << theParam.fNumIALeftMfpFloor << " [mm] (derivative-only floor, 0=off)" << std::endl;
  std::cout << "         - gamma-numia-mfp-floor: "     << theParam.fGammaNumIALeftMfpFloor << " [mm] (derivative-only floor, 0=off)" << std::endl;
  std::cout << "         - gamma-pe-ekin-floor  : "     << theParam.fGammaPhotoelectricEkinFloor << " [MeV] (derivative-only floor, 0=off)" << std::endl;
  std::cout << "         - box-dir-den-floor    : "     << theParam.fBoxDirDenFloor << " (derivative-only, 0=off)" << std::endl;
  std::cout << "         - rotate-up-floor      : "     << theParam.fRotateUpDerivativeFloor << " (derivative-only, 0=off)" << std::endl;
  std::cout << "         - conversion-reg-eps   : "     << theParam.fConversionDerivativeEpsilon << " (derivative-only, 0=off)" << std::endl;
  std::cout << "         - umsc-cos-den-floor   : "     << theParam.fUmscCosThetaDenFloor << " (derivative-only, 0=off)" << std::endl;
  std::cout << "         - umsc-tau-blend-eps   : "     << theParam.fUmscTauBlendEpsilon << " (derivative-only, 0=off)" << std::endl;
  std::cout << "         - umsc-simple-den-floor: "     << theParam.fUmscSimpleDenFloor << " (derivative-only, 0=off)" << std::endl;
  std::cout << "         - umsc-disp-rad-floor  : "     << theParam.fUmscDispRadFloor << " (derivative-only, 0=off)" << std::endl;
  std::cout << "         - ke-cut-threshold    : ";
  if (theParam.fEnableKECut) {
    std::cout << theParam.fKECut << " [MeV]" << std::endl;
  } else {
    std::cout << "disabled" << std::endl;
  }

}


//
// options for providign input arguments to the `HepEmShow` application
static struct option options[] = {
  {"number-of-layers      (number of layers in the calorimeter)           - default: 50"     , required_argument, 0, 'l'},
  {"absorber-thickness    (in [mm] units)                                 - default: 2.3"    , required_argument, 0, 'a'},
  {"gap-thickness         (in [mm] units)                                 - default: 5.7"    , required_argument, 0, 'g'},
  {"transverse-size       (of the calorimeter in [mm] units)              - default: 400"    , required_argument, 0, 't'},

  {"primary-particle      (possible particle names: e-, e+ and gamma)     - default: e-"     , required_argument, 0, 'p'},
  {"primary-energy        (in [MeV] units)                                - default: 10 000" , required_argument, 0, 'e'},
  {"number-of-events      (number of primary events to simulate)          - default: 1000"   , required_argument, 0, 'n'},
  {"random-seed                                                           - default: 1234"   , required_argument, 0, 's'},

  {"g4hepem-data-file     (the pre-generated data file with its path)     - default: ../data/hepem_data" , required_argument, 0, 'd'},
  #ifdef CODI_REVERSE
    {"edep-bars             (bar values of edeps, in [MeV] units)           - default:: 0:0:...:0", required_argument, 0, 'b'},
  #endif
  {"run-verbosity         (verbosity of run information: nothing when 0)  - default: 1"      , required_argument, 0, 'v'},
  {"threshold             (|vx| threshold for grazing condition) - default: 0.1"              , required_argument, 0, 'f'},
  {"threshold2            (near-boundary safety threshold [mm], <=0 uses internal default) - default:-1.0", required_argument, 0, 'k'},
  {"grazing-stop-track    (1: disable gradient for full track, 0: current-step sanitize only) - default: 1", required_argument, 0, 'y'},
  {"backward-boundary-stop (1: treat boundary-limited backward vx<0 as unsafe, 0: disable this trigger) - default: 0", required_argument, 0, 'B'},
  {"msc-displacement      (1: enable MSC lateral displacement, 0: disable) - default: 1"    , required_argument, 0, 'm'},
  {"msc-step-random       (1: enable UMSC step-limit randomization, 0: deterministic) - default: 1", required_argument, 0, 'r'},
  {"boundary-tolerance    (distance-to-boundary tolerance in [mm]) - default: 0.0"           , required_argument, 0, 'u'},
  {"msc-disp-safe-floor   (post-step safety floor [mm] for displacement clipping) - default: 0.0", required_argument, 0, 'w'},
  {"same-boundary-stop    (stop-grad after N total hits of same x-boundary plane over track history, 0 disables) - default: 0", required_argument, 0, 'q'},
  {"same-boundary-pos-tol (tolerance [mm] to identify same x-boundary plane) - default: 1e-6", required_argument, 0, 'z'},
  {"same-boundary-min-flips (minimum vx sign flips on same boundary over track history before triggering) - default: 1", required_argument, 0, 'j'},
  {"same-boundary-full-track (1: full-track stopgrad on same-boundary trigger, 0: step-local sanitize) - default: 0", required_argument, 0, 'o'},
  {"same-boundary-hard-stop (full-track stopgrad when total hits on same boundary reach this count; 0 disables) - default: 0", required_argument, 0, 'i'},
  {"ke-cut-threshold      (set kinetic energy cut in [MeV], disabled when absent)"            , required_argument, 0, 'c'},
  {"numia-mfp-floor       (derivative-only mfp floor [mm] for UpdateNumIALeft, 0 disables) - default: 0.0", required_argument, 0, 'A'},
  {"gamma-numia-mfp-floor (derivative-only mfp floor [mm] for Gamma UpdateNumIALeft, 0 disables) - default: 0.0", required_argument, 0, 'T'},
  {"gamma-pe-ekin-floor   (derivative-only ekin floor [MeV] for 1/ekin in gamma photoelectric xsec, 0 disables) - default: 0.0", required_argument, 0, 'U'},
  {"box-dir-den-floor     (derivative-only signed floor for Box::DistanceToOut direction denominators vx/vy/vz, 0 disables) - default: 0.0", required_argument, 0, 'V'},
  {"rotate-up-floor       (derivative-only floor for RotateToReferenceFrame denominator, 0 disables) - default: 0.0", required_argument, 0, 'F'},
  {"conversion-reg-eps    (derivative-only epsilon for MSC true/geom conversion regularization, 0 disables) - default: 0.0", required_argument, 0, 'N'},
  {"umsc-cos-den-floor    (derivative-only denominator floor in UMSC SampleCosineTheta ratios, 0 disables) - default: 0.0", required_argument, 0, 'P'},
  {"umsc-tau-blend-eps    (derivative-only smoothing width around UMSC tau branch threshold, 0 disables) - default: 0.0", required_argument, 0, 'Q'},
  {"umsc-simple-den-floor (derivative-only denominator floor in UMSC SimpleScattering, 0 disables) - default: 0.0", required_argument, 0, 'R'},
  {"umsc-disp-rad-floor   (derivative-only floor for UMSC displacement sqrt radicand derivative, 0 disables) - default: 0.0", required_argument, 0, 'S'},
  {"stop-grad-mode (0:none,1:track,2:track+desc) - default: 2"         , required_argument, 0, 'x'},
  {"help"                                                                                    , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};


// auxiliary functions for obtaining input arguments
void Help() {
  std::cout<<"\n === Usage: HepEmShow [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
}

// Parses string 'number1:number2:...:numberN' into vector {number1,...,numberN}.
static inline std::vector<double> stod_array(const char* arg){
  std::vector<double> ret;
  std::string arg_s(arg);
  arg_s = arg_s + ":";
  size_t pos = 0;
  do {
     int sep = arg_s.find(":",pos);
     std::string s = arg_s.substr(pos,sep-pos);
     ret.push_back(std::stod(s.c_str()));
     pos = sep+1;
  } while(pos<arg_s.size());
  return ret;
}

// In forward-mode AD, allow real-number arguments to consist of two numbers separated by ':'.
// If existent, the second number specifies the dot value of the input argument.
static inline G4double parseRealInput(const char* arg){
  std::vector<double> arg_d = stod_array(arg);
  if(arg_d.size()==1){
     return arg_d[0];
  } else {
     if(arg_d.size()>2){
        std::cerr << "Specify 'number' or 'number:number', not more than two elements." << std::endl;
     }
     G4double ret = arg_d[0];
     #ifdef CODI_FORWARD
        SET_DOTVALUE(ret, arg_d[1]);
     #else
        std::cerr << "Ignoring specification of dot value in argument, as this is not a forward-mode AD build." << std::endl;
     #endif
     return ret;
  }
}


void GetOpt(int argc, char *argv[], InputParameters& param) {
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "hl:a:g:t:p:e:n:s:d:v:b:f:k:y:B:m:r:u:w:q:z:j:o:i:c:A:T:U:V:F:N:P:Q:R:S:x:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
       c = options[optidx].val;
       /* fall through */

    case 'l':
       param.fGeometry.fNumLayers = std::stoi(optarg);
       break;
    case 'a':
       param.fGeometry.fThicknessAbsorber = parseRealInput(optarg);
       break;
    case 'g':
       param.fGeometry.fThicknessGap = parseRealInput(optarg);
       break;
    case 't':
       param.fGeometry.fSizeTransverse = std::stod(optarg);
       break;

    case 'x': {
       const int mode = std::stoi(optarg);
       if (mode < 0 || mode > 2) {
         std::cerr << "Unsupported stop-grad-mode value: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fGradientStopMode = mode;
       break;
     }

    case 'p':
       param.fPrimaryAndEvents.fParticleName = optarg;
       if ( !(param.fPrimaryAndEvents.fParticleName=="e-" || param.fPrimaryAndEvents.fParticleName=="e+" || param.fPrimaryAndEvents.fParticleName=="gamma") ) {
         std::cout << "\n *** Unknown primary particle name -p: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       break;
    case 'e':
       param.fPrimaryAndEvents.fParticleEnergy = parseRealInput(optarg);
       break;
    case 'n':
       param.fPrimaryAndEvents.fNumEvents = std::stoi(optarg);
       break;
    case 's':
       param.fPrimaryAndEvents.fRandomSeed = std::stod(optarg);
       break;

    case 'd':
       param.fG4HepEmDataFile = optarg;
       break;
    case 'v':
       param.fRunVerbosity = std::stoi(optarg);
       break;

    case 'b':
       #ifdef CODI_REVERSE
          param.barEdep = stod_array(optarg);
       #else
          std::cerr << "Ignoring -b argument, as this is a not a reverse-AD build." << std::endl;
       #endif
       break;
    
    case 'f':
       param.fThreshold = parseRealInput(optarg); // FIX
       break;
    case 'k':
       param.fThreshold2 = parseRealInput(optarg); // FIX
       break;
    case 'y': {
       const int flag = std::stoi(optarg);
       if (flag != 0 && flag != 1) {
         std::cerr << "grazing-stop-track must be 0 or 1: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fGrazingStopsTrack = (flag == 1);
       break;
    }
    case 'B': {
       const int flag = std::stoi(optarg);
       if (flag != 0 && flag != 1) {
         std::cerr << "backward-boundary-stop must be 0 or 1: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fBackwardBoundaryStop = (flag == 1);
       break;
    }
    case 'm': {
       const int flag = std::stoi(optarg);
       if (flag != 0 && flag != 1) {
         std::cerr << "msc-displacement must be 0 or 1: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fEnableMscDisplacement = (flag == 1);
       break;
    }
    case 'r': {
       const int flag = std::stoi(optarg);
       if (flag != 0 && flag != 1) {
         std::cerr << "msc-step-random must be 0 or 1: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fEnableMscStepRandomization = (flag == 1);
       break;
    }
    case 'u': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "boundary-tolerance must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fBoundaryTolerance = value;
       break;
    }
    case 'w': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "msc-disp-safe-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fMscDisplacementSafetyFloor = value;
       break;
    }
    case 'q': {
       const int value = std::stoi(optarg);
       if (value < 0) {
         std::cerr << "same-boundary-stop must be >= 0: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fSameBoundaryStop = value;
       break;
    }
    case 'z': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "same-boundary-pos-tol must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fSameBoundaryPosTolerance = value;
       break;
    }
    case 'j': {
       const int value = std::stoi(optarg);
       if (value < 0) {
         std::cerr << "same-boundary-min-flips must be >= 0: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fSameBoundaryMinFlips = value;
       break;
    }
    case 'o': {
       const int flag = std::stoi(optarg);
       if (flag != 0 && flag != 1) {
         std::cerr << "same-boundary-full-track must be 0 or 1: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fSameBoundaryFullTrackStop = (flag == 1);
       break;
    }
    case 'i': {
       const int value = std::stoi(optarg);
       if (value < 0) {
         std::cerr << "same-boundary-hard-stop must be >= 0: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fSameBoundaryHardStop = value;
       break;
    }
    case 'c': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "ke-cut-threshold must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fEnableKECut = true;
       param.fKECut = value;
       break;
    }
    case 'A': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "numia-mfp-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fNumIALeftMfpFloor = value;
       break;
    }
    case 'T': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "gamma-numia-mfp-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fGammaNumIALeftMfpFloor = value;
       break;
    }
    case 'U': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "gamma-pe-ekin-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fGammaPhotoelectricEkinFloor = value;
       break;
    }
    case 'V': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "box-dir-den-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fBoxDirDenFloor = value;
       break;
    }
    case 'F': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "rotate-up-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fRotateUpDerivativeFloor = value;
       break;
    }
    case 'N': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "conversion-reg-eps must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fConversionDerivativeEpsilon = value;
       break;
    }
    case 'P': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "umsc-cos-den-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fUmscCosThetaDenFloor = value;
       break;
    }
    case 'Q': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "umsc-tau-blend-eps must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fUmscTauBlendEpsilon = value;
       break;
    }
    case 'R': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "umsc-simple-den-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fUmscSimpleDenFloor = value;
       break;
    }
    case 'S': {
       G4double value = parseRealInput(optarg);
       if (value < 0.0) {
         std::cerr << "umsc-disp-rad-floor must be non-negative: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       param.fUmscDispRadFloor = value;
       break;
    }

    case 'h':
       Help();
       exit(-1);
       break;

    default:
      printf("\n *** Unknown input argument: %c\n",c);
      Help();
      exit(-1);
    }
   }
   // number of layers must be >= 1
   if (param.fGeometry.fNumLayers < 1 ) {
     printf("\n *** Calorimeter number of layers must be >= 1! \n");
     Help();
     exit(-1);
   }
   // check if the data file was given with/without extension
   if (param.fG4HepEmDataFile.find(".json")==std::string::npos) {
     param.fG4HepEmDataFile += ".json";
   }
   // print parameters if the verbosity > 0
   if (param.fRunVerbosity > 0) {
     PrintParameters(param);
   }
}



#endif // INPUTPARAMETERS_HH
