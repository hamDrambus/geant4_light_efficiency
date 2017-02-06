#include "PseudoMesh.hh"

PseudoMeshData::PseudoMeshData(void)
{
	curr_mesh = NULL;
}

void PseudoMeshData::SetDefauldInd(void)
{
	curr_x_id = -1;
	curr_y_id = -1;
}

void PseudoMesh::hex_indices_by_pos(PseudoMeshData *data, G4double x, G4double y, G4int iX_max, G4int iY_max)
//^x,y from [0,0] to [Xmax, Ymax]
{
	data->curr_x_id = -1;
	data->curr_y_id= -1;
	G4int approx_ix = (int)(x / (2 * cell_size)) + 1;
	G4int approx_iy = (int)((y - cell_size / sqrt(3)) / (sqrt(3)* cell_size)) + 1;
	//error in the above is +-1
	G4int found = 0;
	G4int index_x, index_y;
	for (int iix = -1; iix < 2; iix++) //finding the exact center (x,y) in in the hexagonal
	{
		for (int iiy = -1; iiy < 2; iiy++)
		{
			index_x = approx_ix + iix;
			index_y = approx_iy + iiy;
			G4double x_app_center = ((index_x - 1) * 2 * cell_size) + cell_size + ((index_y % 2) ? 0 : cell_size); //x center depends on index of y
			G4double y_app_center = ((index_y - 1)*cell_size*sqrt(3)) + (2 * cell_size / sqrt(3));
			G4double dx = x - x_app_center;
			G4double dy = y - y_app_center;
			if ((abs(dx) <= cell_size) && (abs(dy - (dx / sqrt(3))) < (2 * cell_size / sqrt(3))) 
				&& (abs(dy + (dx / sqrt(3))) < (2 * cell_size / sqrt(3))))
				//^condition of x,y being in the hexagonal with the center at (x_app_center y_app_center)
			{
				found = 1;
				goto for_out;
			}
		}
	}
	for_out:
	if (found)
	{
		if ((index_x <= 0) || (index_y <= 0) || (index_x>iX_max) || (index_y > iY_max))
			return;
		data->curr_x_id = index_x;
		data->curr_y_id = index_y;
	}
}

//processes leaving the cell and transporting back to parent
G4ThreeVector PseudoMesh::hex_on_leave(PseudoMeshData *data, const G4Step& step,
	G4double active_x_par, G4double active_y_par)
{
	G4ThreeVector pos= step.GetPostStepPoint()->GetPosition();
	G4double dz = (pos - g_cell_posit).z();
	G4double x_center = (data->curr_x_id - 1) * 2 * cell_size + cell_size + ((data->curr_y_id % 2) ? 0 : cell_size); //x center depends on index of y
	G4double y_center = ((data->curr_y_id - 1)*cell_size*sqrt(3)) + (2 * cell_size / sqrt(3)); //center of a cell, measured from edge of mapped area
	G4double x_rel = pos.x() - g_cell_posit.x() + x_center;
	G4double y_rel = pos.y() - g_cell_posit.y() + y_center; //relative to edge of mapped area
	x_rel = x_rel - (active_x_par/ 2);
	y_rel = y_rel - (active_y_par/ 2);
	G4ThreeVector new_pos(x_rel, y_rel, dz);
	new_pos = new_pos + g_parent_pos;
	data->curr_x_id = -1;
	data->curr_y_id = -1;
#ifdef DEBUG_MC_NODES
	G4cout << "MAPPING: Cell mapping to grid. Proposed position: " << new_pos << G4endl;
#endif
	return new_pos;
}

//    .*.
//  *     *
// |    ___| cell_size== ___
// |       |
//  *.   .*
//     *
//TODO:
//1) account for entering the grid via sides
//2) account for transformations (rotations) (and make them via separate methods)
G4ThreeVector PseudoMesh::hexagonal_mapping(PseudoMeshData *data,
	const G4Track& aTrack, const G4Step& aStep)
{
	G4ThreeVector ret = aStep.GetPostStepPoint()->GetPosition();
	G4int X_max_index = int (par_x_size/(2*cell_size));
	G4int Y_max_index = int((par_y_size - cell_size / (sqrt(3))) / (cell_size*sqrt(3)));
	if (X_max_index <= 0 || Y_max_index <= 0) return ret; //no grid - no changes
	//I have all indices starting from 1
	G4double active_x_par = X_max_index * 2 * cell_size;  //only this length is mapped
	G4double active_y_par = (Y_max_index * cell_size*sqrt(3))+(cell_size/sqrt(3));
	if ((data->curr_x_id < 0)||(data->curr_y_id<0))//photon from outside
	{
		G4ThreeVector pos = aStep.GetPostStepPoint()->GetPosition();
		G4ThreeVector dir = aStep.GetPostStepPoint()->GetMomentumDirection();
		G4double dz = (pos - g_parent_pos).z();
		G4double tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
		if ((abs((abs(dz) - par_z_size / 2)) < tolerance) && (G4ThreeVector(0, 0, dz)*dir < 0)) //on the surface and moves inside
		{
			G4double x_rel = pos.x() - g_parent_pos.x() + par_x_size / 2 - ((par_x_size - active_x_par) / 2); //x_rel starts form the edge of mapped plate 
			//(mapped not entire one but the part which contains integer amount of cells (hence offset (x_par - active_x_par) / 2))
			G4double y_rel = pos.y() - g_parent_pos.y() + par_y_size / 2 - ((par_y_size - active_y_par) / 2); //should be granted to be positive
			if ((x_rel < 0) || (y_rel < 0))
				return ret;
			if ((x_rel>active_x_par) || (y_rel > active_y_par))
				return ret;
			G4int x_i, y_i;
			hex_indices_by_pos(data,x_rel, y_rel, X_max_index, Y_max_index);
			if ((data->curr_x_id< 0) || (data->curr_y_id< 0))
				return ret;
			x_i = data->curr_x_id;
			y_i = data->curr_y_id;
			G4double a = cell_size;
			G4double x_center = (x_i - 1) * 2 * a + a + ((y_i % 2) ? 0 : a); //x center depends on index of y
			G4double y_center = ((y_i - 1)*a*sqrt(3))+ (2 * a / sqrt(3));
			G4ThreeVector new_pos(x_rel - x_center, y_rel - y_center, dz);
			new_pos = new_pos + g_cell_posit;
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Grid mapping to cell. Proposed position: " << new_pos<<G4endl;
#endif
			return new_pos;
		}
		else //not on the surface or moves outside
			//TODO: acoount for photons entering from sides
			//THIS IS QUITE IMPORTANT
			return ret; //no changes are proposed
	}
	else //photon is already in the cell
	{
		G4ThreeVector pos = aStep.GetPostStepPoint()->GetPosition();
		G4ThreeVector dir = aStep.GetPostStepPoint()->GetMomentumDirection();
		G4double dz = (pos - g_cell_posit).z();
		G4double tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
		if ((abs((abs(dz) - par_z_size / 2)) < tolerance) && (G4ThreeVector(0, 0, dz)*dir > 0)) //on the surface and moves outside
			return hex_on_leave(data, aStep,active_x_par, active_y_par);
#ifdef DEBUG_MC_NODES
		G4cout << "MAPPING: Cell to cell map. From side."<<G4endl;
#endif
		//process moving from cell's sides
		//x-y+ .*.  x+y+
		//   *     *
		//x-|       | x+
		//  |       |
		//   *.   .*
		//x-y-  *  x+y-
		G4double x_rel = pos.x() - g_cell_posit.x();
		G4double y_rel = pos.y() - g_cell_posit.y();
		G4int y_index_odd = data->curr_y_id - 2 * ((int)(data->curr_y_id / 2.0));
		if ((abs(cell_size - x_rel)<tolerance)&&(dir.x()>0)) //leaving x+ edge
		{
			data->curr_x_id++;
			if (data->curr_x_id > X_max_index)
			{//leave cell and go to parent
				data->curr_x_id--;
				return hex_on_leave(data, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() - 2 * cell_size, pos.y(), pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x+. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(cell_size + x_rel)<tolerance)&&(dir.x()<0)) //x- edge
		{
			data->curr_x_id--;
			if (data->curr_x_id < 1)
			{//leave cell and go to parent
				data->curr_x_id++;
				return hex_on_leave(data, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() + 2 * cell_size, pos.y(), pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x-. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel - 2*cell_size/sqrt(3) +x_rel / sqrt(3))<tolerance)&&(dir*G4ThreeVector(0.5,sqrt(3)/2,0)>0)) //x+,y+
		{
			data->curr_x_id+=(y_index_odd?0:1);
			data->curr_y_id++;
			if ((data->curr_y_id>Y_max_index) || (data->curr_x_id >X_max_index))
			{//leave cell and go to parent
				data->curr_x_id -= (y_index_odd ? 0 : 1);
				data->curr_y_id--;
				return hex_on_leave(data, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() - 2 * cell_size*0.5, pos.y() - 2 * cell_size*sqrt(3) / 2, pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x+y+. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel - 2*cell_size/sqrt(3) - x_rel / sqrt(3))<tolerance) && (dir*G4ThreeVector(-0.5, sqrt(3) / 2, 0)>0)) //x-,y+
		{
			data->curr_x_id -= (y_index_odd ? 1 : 0); //- picture required - this is the way numeration works
			data->curr_y_id++;
			if ((data->curr_y_id>Y_max_index) || (data->curr_x_id <1))
			{//leave cell and go to parent
				data->curr_x_id += (y_index_odd ? 1 : 0);
				data->curr_y_id--;
				return hex_on_leave(data, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() + 2 * cell_size*0.5, pos.y() - 2 * cell_size*sqrt(3) / 2, pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x-y+. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel + 2*cell_size/sqrt(3) + x_rel / sqrt(3))<tolerance) && (dir*G4ThreeVector(-0.5, -sqrt(3) / 2, 0)>0)) //x-,y-
		{
			data->curr_x_id -= (y_index_odd ? 1 : 0); //- picture required - this is the way numeration works
			data->curr_y_id--;
			if ((data->curr_y_id<1) || (data->curr_x_id <1))
			{//leave cell and go to parent
				data->curr_x_id += (y_index_odd ? 1 : 0);
				data->curr_y_id++;
				return hex_on_leave(data,aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() + 2 * cell_size*0.5, pos.y() + 2 * cell_size*sqrt(3) / 2, pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x-y-. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel + 2*cell_size/sqrt(3) - x_rel / sqrt(3))<tolerance) && (dir*G4ThreeVector(0.5, -sqrt(3) / 2, 0)>0)) //x+,y-
		{
			data->curr_x_id += (y_index_odd ? 0 : 1); //- picture required - this is the way numeration works
			data->curr_y_id--;
			if ((data->curr_y_id<1) || (data->curr_x_id <1))
			{//leave cell and go to parent
				data->curr_x_id -= (y_index_odd ? 0 : 1);
				data->curr_y_id++;
				return hex_on_leave(data, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() - 2 * cell_size*0.5, pos.y() + 2 * cell_size*sqrt(3) / 2, pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x+y-. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		//not at border - no changes
		return ret;
	}
}

void PseudoMesh::GetDefaultMappingData(PseudoMeshData* data)
{
	G4int X_max_index = int(par_x_size / (2 * cell_size));
	G4int Y_max_index = int((par_y_size - cell_size / (sqrt(3))) / (cell_size*sqrt(3)));
	data->curr_x_id = X_max_index / 2;
	data->curr_y_id = Y_max_index / 2; //~center
}

PseudoMesh::PseudoMesh(G4LogicalVolume* parent_vol, G4ThreeVector g_par_pos, G4double par_x, G4double par_y, G4double par_z,
	G4LogicalVolume* cell_vol, G4ThreeVector g_cell_pos, G4double cell_parameter/*, box_map_function function*/)
{
	//map_function = function;
	parent = parent_vol;
	g_parent_pos = g_par_pos;
	par_x_size = par_x;
	par_y_size = par_y;
	par_z_size = par_z;
	cell = cell_vol;
	g_cell_posit = g_cell_pos;
	cell_size = cell_parameter;
}

PseudoMesh::~PseudoMesh()
{
}

//returns global position
//detector_construction makes checks of entering/leaveing
G4ThreeVector PseudoMesh::PostSteppingAction(PseudoMeshData* psm_data, const G4Track& aTrack, const G4Step& aStep, G4TouchableHandle &post_step_handle) //checks entering/leaving, called from DetectorConstruction, manages everything
{
	G4ThreeVector ret = hexagonal_mapping(psm_data, aTrack, aStep);
	if ((psm_data->curr_x_id < 0) || (psm_data->curr_y_id < 0))
	{
		if (psm_data->curr_mesh == this)
			psm_data->curr_mesh = NULL;
	}
	else
		psm_data->curr_mesh = this;
	return ret;
}
