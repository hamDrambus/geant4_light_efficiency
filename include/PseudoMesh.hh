#ifndef PSEUDO_MESH_HH
#define PSEUDO_MESH_HH
//TODO: should I make this via inhereting Navigator?
#include "CustomRunManager.hh"
#include "B1DetectorConstruction.hh"


//TODO: account for rotation (add transformation of local to global)? 
//meaning of cell_size may vary depending on its shape
//returns new position
typedef G4ThreeVector(*box_map_function)(G4ThreeVector GlobalParentPos /*center*/, G4double x_par, G4double y_par, G4double z_par,
	G4ThreeVector GlobalCellPos, G4double cell_size, G4int *curr_x_id, G4int *curr_y_id,
	const G4Track& aTrack, const G4Step& aStep);

//    .*.
//  *     *
// |    ___| cell_size== ___
// |       |
//  *.   .*
//     *
//TODO:
//1) account for entering the grid via sides
//2) acoount for transformations (rotations) (and make them via separate methods)
G4ThreeVector hexagonal_mapping(G4ThreeVector GlobalParentPos /*center*/, G4double x_par, G4double y_par, G4double z_par,
	G4ThreeVector GlobalCellPos, G4double cell_size, G4int *curr_x_id, G4int *curr_y_id,
	const G4Track& aTrack, const G4Step& aStep);
//step contains position and momentum, then changed by this function
//NOTE: this class and box_map_function is dependendant on geometry and the function must be rewritten for different shapes
class PseudoMesh
{
public:
	G4LogicalVolume* parent;
	G4LogicalVolume* cell;
	G4int *curr_x_id, *curr_y_id; //contain only single int value
	
	G4ThreeVector g_parent_pos;
	G4double par_x_size;
	G4double par_y_size;
	G4double par_z_size;
	
	G4ThreeVector g_cell_posit;
	G4double cell_size;
	
	box_map_function map_function;
	PseudoMesh(G4LogicalVolume* parent, G4ThreeVector g_par_pos, G4double par_x, G4double par_y, G4double par_z,
		G4LogicalVolume* cell, G4ThreeVector g_cell_pos, G4double cell_parameter, box_map_function function);
	G4ThreeVector PostSteppingAction(const G4Track& aTrack, const G4Step& aStep, G4TouchableHandle &post_step_handle);
	//^checks entering/leaving, called from UserSteppingAction, manages everything
	~PseudoMesh();
};

#endif
