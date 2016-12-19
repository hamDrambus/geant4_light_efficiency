#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "CustomRunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "globals.hh"
#include <list>

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

	G4double GetHitProbability(G4LogicalVolume* in_volume, G4double particle_energy, G4ThreeVector Momentum_dir, G4ThreeVector position);
	//^returns -1 if there is no hit (continue run) or probability of photon detection [0;1]
	std::list<G4LogicalVolume*> detector_volumes;
	std::list<G4LogicalVolume*> absorbtion_volumes;

private:
	G4Material* _Argon_mat(void);
	G4Material* _WLS_mat(void);
	G4Material* _Copper_mat(void);
	G4Material* _Acrylic_mat(void);
	G4Material* _LArgon_mat(void);
	G4Material* _FusedSilica_mat(void); //PMT
	void _SetCopperSurface(G4OpticalSurface* surface);
	void _SetVisibilityParameters(void);
public:
	G4LogicalVolume* top_plate; //gem, but at the moment is just a solid plate
	G4LogicalVolume* bot_plate;
	G4LogicalVolume* LAr_layer;

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

	G4LogicalVolume* Xp_PMT_window;
	G4LogicalVolume* Xn_PMT_window;
	G4LogicalVolume* Yp_PMT_window;
	G4LogicalVolume* Yn_PMT_window;

	//I need this to be placed inside windows and used as detectors 
	//(because the code requires the detector volumes to be slightly inside its physical representation so that optics works as required)
	//because in case of reflection end step point is placed inside volume photon is reflected form and this leads to Hit if that volume is detector 
	//even though reflection occured. TODO: This can also be eluminated via this->GetHitProbability() by considering photon direction
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
};

#endif

