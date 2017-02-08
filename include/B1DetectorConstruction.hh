#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "PseudoMesh.hh"
#include "CustomRunManager.hh"
#include <list>

class G4VPhysicalVolume;
class G4LogicalVolume;
class PseudoMesh;
class PseudoMeshData;
class CustomRunManager;

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

	G4double GetHitProbability(G4StepPoint* post_point);
	//^returns -1 if there is no hit (continue run) or probability of photon detection [0;1]
	std::list<G4LogicalVolume*> detector_volumes;
	std::list<G4LogicalVolume*> absorbtion_volumes;

private:
	G4RotationMatrix *x_rot_m;
	G4RotationMatrix *y_rot_m;
	G4Material* _Argon_mat(void);
	G4Material* _WLS_mat(void);
	G4Material* _Copper_mat(void);
	G4Material* _Acrylic_mat(void);
	G4Material* _LArgon_mat(void);
	G4Material* _FusedSilica_mat(void); //PMT
	G4Material* _fr4_mat(void);
	void _SetCopperSurface(G4OpticalSurface* surface);
	void _SetVisibilityParameters(void);
public:
	G4double GetSafeOffset(void);
	//returns small length which can be taken from volume boundary and be guaranteed to lead to a neighbour one.
	//e.g. auxilary volumes are always have offset of 0.5-1 mm from thier content 
	void OnEventStartProc(CustomRunManager* manman); //not used now
	G4ThreeVector MappingProc(PseudoMeshData *psm_data, const G4Track& track, 
		const G4Step& aStep, G4TouchableHandle &fCurrentTouchableHandle);
	void GetPseudoMeshByPoint(PseudoMeshData* data, G4ThreeVector pos, G4ThreeVector momDir); //changes only data->curr_mesh
	G4VPhysicalVolume* GetPVolumeByPoint(G4ThreeVector pos, G4ThreeVector momDir);
	G4bool FindDaughterPhysicalInL(G4VPhysicalVolume* to_find, G4LogicalVolume *origin); //recursive

	G4LogicalVolume* top_ps_plate; //just a solid plate
	G4LogicalVolume* bot_ps_plate;
	G4LogicalVolume* top_cu_plate; //just a solid plate
	G4LogicalVolume* bot_cu_plate;
	G4LogicalVolume* LAr_layer;

	G4LogicalVolume* top_cell_container;
	G4LogicalVolume* top_cell;
	G4LogicalVolume* top_cell_hole;
	G4LogicalVolume* top_cell_hole_dielectric;
	
	G4LogicalVolume* bot_cell_container;
	G4LogicalVolume* bot_cell;
	G4LogicalVolume* bot_cell_hole;
	G4LogicalVolume* bot_cell_hole_dielectric;
#ifdef TOP_MESH_TEST
	G4LogicalVolume* top_mesh_test_detector;
#else
	G4LogicalVolume* top_mesh_absorber; //between cell and GEM, absorbs light passed through GEM
#endif
	G4LogicalVolume* bot_mesh_absorber;
	PseudoMesh* top_GEM;
	PseudoMesh* bot_GEM;

	G4LogicalVolume* Xp_wls; //p-plus, m-minus
	G4LogicalVolume* Xn_wls;
	G4LogicalVolume* Yp_wls;
	G4LogicalVolume* Yn_wls;

	G4LogicalVolume* Xp_acrylic;
	G4LogicalVolume* Xn_acrylic;
	G4LogicalVolume* Yp_acrylic;
	G4LogicalVolume* Yn_acrylic;

	G4LogicalVolume* Xp_LAr_gap;
	G4LogicalVolume* Xn_LAr_gap;
	G4LogicalVolume* Yp_LAr_gap;
	G4LogicalVolume* Yn_LAr_gap;

	//deprecated:
	//G4LogicalVolume* Xp_PMT_window;
	//G4LogicalVolume* Xn_PMT_window;
	//G4LogicalVolume* Yp_PMT_window;
	//G4LogicalVolume* Yn_PMT_window;
	//I needed this to be placed inside windows and used as detectors 
	//(because the code requires the detector volumes to be slightly inside its physical representation so that optics works as required)
	//because in case of reflection end step point is placed inside volume photon is reflected form and this leads to Hit if that volume is detector 
	//even though reflection occured. DONE TODO: This can also be eluminated via this->GetHitProbability() by considering photon direction

	G4LogicalVolume* Xp_PMT;
	G4LogicalVolume* Xn_PMT;
	G4LogicalVolume* Yp_PMT;
	G4LogicalVolume* Yn_PMT;

	G4LogicalVolume* envelope;
	G4LogicalVolume* world;
	G4LogicalVolume* box_interior;
	G4LogicalVolume* box;

private:
	void wavelen_transmittance_to_table(G4MaterialPropertiesTable* table, const G4double* wl, const G4double* tr, G4int size, G4double width);
	void gen_integral_en_spec (G4String file, G4MaterialPropertiesTable* table);
	void gen_PMT_QE(G4String file, G4MaterialPropertiesTable* table);
	void read_table_En(G4String file, G4double* &Ens, G4double* &ys, G4int &size);
	G4ThreeVector get_global_normal(G4StepPoint* point, G4int *validity);
};

#endif

