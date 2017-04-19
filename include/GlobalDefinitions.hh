#ifndef global_definitions_hh
#define global_definitions_hh

//#define DEBUG_MC_NODES

//#define TOP_MESH_TEST
//^if defined, then additional detectors boxes are created:
// above the topmost pseudo GEM (but below cell volume) and LAr layer also becomes a detector
//Photon posistions upon the hit the are written in bmp (one file for LAr, second for additional top box). 
//Also photons are generated orthogonally to the plate in order to 'scan' it. 
//Position of initial photons follows the square pattern with small fluctuations. (initial distributions are overridden)

//#define TEST_MESH_SIDEWAYS
//^for ckecking that moving from cell to cell is right, check are solely by debug points

#ifdef TEST_MESH_SIDEWAYS
#undef TOP_MESH_TEST
#endif

#define FIX_TRANSMITTANCE_
//^if defined, then transmittance spectras are normalized to their top plateau
//this is correct, because the deviation of the plateau from unity is caused by refractive index,
//not absorbtion. Since n is used separeately transmittance spectras must be adjusted. They were
//obtained experimentally, with several boundary processes. Adding the calculated from n reflectance
//wouldn't be right, because if there is no transmittance there would be several boundary proceses and 
//hence higher R. The validity of n can be checked from R at zero transmittance (which means light comes throgh one boundary
//to the WLS and gets absorbed there), when it can be calculated by fresnel law.
//\/DONE: T in such configuration (wls+acrylic, ideal PMT, 2.9eV gamma, only Fresnel) is 0.91132 (direct 
//\/MC simulation), which is quite close to 0.8967, so the precedure above is valid
//DONE: TODO: check renormalization by plateau with R calclated for wls-acrilyc plate (as in source of data)
//TODO: ideally transmittance adjustment should be done by simulating T without absorbtion in the same configuration
//by refractive indices, but I think that would be an overkill

#define REFLECTIVITY_G10_MODEL
//^if defined then surface with reflectance is created for g10 (referred as fr4 sometimes, but g10 is used)
#if defined(TOP_MESH_TEST)||defined(TEST_MESH_SIDEWAYS)
#undef REFLECTIVITY_G10_MODEL
#endif
#define REFLECTIVITY_COPPER_MODEL
//^if defined then copper otical properties are set by reflectivity, not complex refractive index
#ifdef REFLECTIVITY_COPPER_MODEL
#define REFINED_CU_MODEL
//^if defined, then zero reflectivity is taken for VUV photons.
#if defined(TOP_MESH_TEST)||defined(TEST_MESH_SIDEWAYS)
#undef REFINED_CU_MODEL
#endif
#endif

//#define NO_QE_
//#define NO_WLS_PROC
//#define ZERO_CE_
//#define UNITY_CE_
//#define NO_BACK_SCATTERING
//#define NO_VUV_WLS_FRESNEL

//#define AR_EMISSION_NITRO
//^if defined then initial photon's energy is defined by a continoius spectrum of N2 admixture
//otherwise 128nm line is taken
#ifdef AR_EMISSION_NITRO
//#define AR_SPEC_TEST
//^if defiened, then no memory decrease, top absorbed behind mesh is detector, all photons go through GEM to it;
//^and get_detected_spectrum called after simulation. All in order to check Ar emission spectrum generation
#endif

//#define SPATIAL_ANGLE_
#ifdef SPATIAL_ANGLE_
#define NO_QE_
#define UNITY_CE_
#endif

//\/in mm, diameter of cylindrical area where photons are generated
#define CONVERSION_DIAM 52

#define TEMP_CODE_
#define WLS_FILM_WIDTH 100*micrometer
#define PMMA_WIDTH 1.5*mm
#define PMT_DIAMETER 51*mm
#define plate_W	0.6
#define plate_real_W 0.5
//\/from top GEM's container (plate_W) to the center
#define CELL_OFFSET_TOP 10
#define ABSORBER_OFFSET_TOP 5
//\/from bottom GEM's container (plate_W) to the center
#define CELL_OFFSET_BOT 10
#define ABSORBER_OFFSET_BOT 5
//\/in mm
#define CONVERSION_AREA_DIAM 52
//\/in mm
#define GEM_SIZE_TOTAL 100
//\/if defined then proper holes are used - withoud slope, but drilled
#define FIX_GEM_HOLE_
#define MIN_ALLOWED_PROBABILITY 1e-6
//^was 4e-6

#define MIN_ALLOWED_STEPPING 4e-6
#define SPEC_INTEGRATION_STEPS 1000
#define RM_PHOTON_UNDEFINED	-2
#define RM_CHOOSE_KILL	0
#define RM_CHOOSE_DEFL	1
#define RM_CHOOSE_REFL	2	
#define RM_CHOOSE_BOTH	3

#define WLS_SPECTRUM_FILE "WLS_spec.txt"
#define PMT_QE_FILE "PMT_QE.txt"
#define ARGON_SPECTRUM_FILE "Ar_spec.txt"
#define TEST_WLS_OUT_SPEC "WLS_spec_out.txt"
#define TEST_WLS_OUT_I_SPEC "WLS_integral_spec_out.txt"
#define TEST_OUT_I_SPEC1 "WLS_spec_out_reverse.txt"
#ifdef REFINED_CU_MODEL
#define COPPER_REFLECTIVITY "CU_reflectivity_refined.txt"
#else
#define COPPER_REFLECTIVITY "CU_reflectivity.txt"
#endif
#define G10_REFLECTIVITY "G10_reflectivity.txt"

#define RM_PHOTON_NO_TRACK	-1
#define RM_PHOTON_PROPOG	0
#define RM_PHOTON_BOUND		1
#define RM_PHOTON_TOT_REFL	2
#define RM_PHOTON_REFL		3
#define RM_PHOTON_DEFL		4
#define RM_PHOTON_ABSORB	5
#define RM_PHOTON_DETECTION	6

#endif