#ifndef PSEUDO_MESH_HH
#define PSEUDO_MESH_HH
//depr. TODO: should I make this via inhereting Navigator?
#include "G4VUserDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4GeometryTolerance.hh"

class PseudoMesh;

class PseudoMeshData
{
public:
	PseudoMesh* curr_mesh;
	G4int curr_x_id, curr_y_id;
	PseudoMeshData(void);
	void SetDefauldInd(void);
};

//    .*.
//  *     *
// |    ___| cell_size== ___
// |       |
//  *.   .*
//     *
//TODO:
//1) account for entering the grid via sides
//2) acoount for transformations (rotations) (and make them via separate methods)
//G4ThreeVector PseudoMesh::hexagonal_mapping(G4ThreeVector GlobalParentPos /*center*/, G4double x_par, G4double y_par, G4double z_par,
//	G4ThreeVector GlobalCellPos, G4double cell_size, G4int *curr_x_id, G4int *curr_y_id,
//	const G4Track& aTrack, const G4Step& aStep);
//step contains position and momentum, then changed by this function
//NOTE: this class and box_map_function is dependendant on geometry and the function must be rewritten for different shapes
class PseudoMesh
{
public:
	G4LogicalVolume* parent;
	G4LogicalVolume* cell;
	//G4int *curr_x_id, *curr_y_id; //contain only single int value
	
#ifdef TEMP_CODE_
	G4ThreeVector last_ret;
#endif

	G4ThreeVector g_parent_pos;/*center*/
	G4double par_x_size;
	G4double par_y_size;
	G4double par_z_size;
	
	G4ThreeVector g_cell_posit;
	G4double cell_size;
	
	//box_map_function map_function;
	PseudoMesh(G4LogicalVolume* parent, G4ThreeVector g_par_pos, G4double par_x, G4double par_y, G4double par_z,
		G4LogicalVolume* cell, G4ThreeVector g_cell_pos, G4double cell_parameter/*, box_map_function function*/);
	G4ThreeVector PostSteppingAction(PseudoMeshData *psmData, const G4Track& aTrack, const G4Step& aStep, G4TouchableHandle &post_step_handle);
	//^checks entering/leaving, called from UserSteppingAction, manages everything
	G4ThreeVector OnStartOfEventProc(G4ThreeVector position, G4ThreeVector mom_direction,G4LogicalVolume* start_volume);
	//^ sets initial state after each start of ProcessOneEvent (current indices)
	//^ DONE !!! TODO: in case start_volume is inside the cell, there is no reset of indeces - for a secondary process occur normally,
	//^ when originating as secondary from previous mapping. This is not correct behaviour. State of mapping must be written as initial
	//^ parameter for secondary processes in order to allow for a multi-nested mapping or several secondaries from one mapping (i.g. from cell)
	//^ This is done quite easily because only two indices per PseudoMesh define the state, but the case of nested mappings
	//^ (e.g. photon mapped from top GEM and hitting copper surface and further adsorbed in the main event (Primary_MC_Node) will give rise to
	//^ a secondary, reflected photon, which in turn also may be mapped and so on) implies that the present class should be split in two: one 
	//^ carriying the state of mapping (just like probability(weight) of photon) and other modifying the state and track of photon.
	//^ Also any photon may be mapped just once at the same time (there is no pressing need for more because there is no nested periodic 
	//^ constructions like GEM inside every cell of a larger one), so the mapping data may contain
	//^ only pointer to PseudoMesh class and the two indices. Positions and touchables are already assigned.
	void GetDefaultMappingData(PseudoMeshData* data); //set some default indices, required 
	//^ when inital position is in the cell and there is no mapping data, which is actually kind of error
	G4ThreeVector PseudoMesh::hexagonal_mapping(PseudoMeshData *data, const G4Track& aTrack, const G4Step& aStep);
	~PseudoMesh();
	//TODO: account for rotation (add transformation of local to global)? 
	//^ meaning of cell_size may vary depending on its shape
	//^ returns new position

protected:
	G4ThreeVector PseudoMesh::hex_on_leave(PseudoMeshData* data, const G4Step& step,
		G4double active_x_par, G4double active_y_par);
	//^processes leaving the cell and transporting back to parent
	void hex_indices_by_pos(PseudoMeshData *data, G4double x, G4double y, G4int iX_max, G4int iY_max);
	//^x,y from [0,0] to [Xmax, Ymax]
};

#endif
