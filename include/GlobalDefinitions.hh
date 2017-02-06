#ifndef global_definitions_hh
#define global_definitions_hh

//#define DEBUG_MC_NODES
//#define TOP_MESH_TEST
//^if defined, then additional detectors boxes are created:
// above the topmost pseudo GEM (but below cell volume) and LAr layer which also becomes a detector
//Photon posistions upon the hit the are written in bmp (one file for LAr, second for additional top box). 
//Also photons are generated orthogonally to the plate in order to 'scan' it. 
//Position of initial photons follows the square pattern with small fluctuations. (initial distributions are overridden)
//#define TEST_MESH_SIDEWAYS

#ifdef TEST_MESH_SIDEWAYS
#undef TOP_MESH_TEST
#endif

#define TEMP_CODE_
//^marks everything that is temporary so I don't forget
#define WLS_FILM_WIDTH 100*micrometer
#define PMMA_WIDTH 1.5*mm
#define PMT_DIAMETER 51*mm
#define MIN_ALLOWED_PROBABILITY 4e-6
#define MIN_ALLOWED_STEPPING 4e-6
#define SPEC_INTEGRATION_STEPS 500
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

#define RM_PHOTON_NO_TRACK	-1
#define RM_PHOTON_PROPOG	0
#define RM_PHOTON_BOUND		1
#define RM_PHOTON_TOT_REFL	2
#define RM_PHOTON_REFL		3
#define RM_PHOTON_DEFL		4
#define RM_PHOTON_ABSORB	5
#define RM_PHOTON_DETECTION	6

#endif