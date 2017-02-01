#include "PseudoMesh.hh"

void hex_indices_by_pos(G4double x, G4double y, G4double a, G4int iX_max, G4int iY_max, G4int *x_i, G4int *y_i)//x,y from [0,0] to [Xmax, Ymax]
{
	*x_i = -1;
	*y_i = -1;
	G4int approx_ix = (int)(x / (2 * a)) + 1;
	G4int approx_iy = (int)((y - a / sqrt(3)) / (sqrt(3)* a)) + 1;
	//error in the above is +-1
	G4int found = 0;
	G4int index_x, index_y;
	for (int iix = -1; iix < 2; iix++) //finding the exact center (x,y) in in the hexagonal
	{
		for (int iiy = -1; iiy < 2; iiy++)
		{
			index_x = approx_ix + iix;
			index_y = approx_iy + iiy;
			G4double x_app_center = (index_x - 1) * 2 * a + a + ((index_y - 2 * ((int)(index_y / 2.0)) == 0) ? a : 0); //x center depends on index of y
			G4double y_app_center = (index_y - 1)*a*sqrt(3) + 2 * a / sqrt(3);
			G4double dx = x - x_app_center;
			G4double dy = y - y_app_center;
			if ((abs(dx) <= a) && (dy < (a*sqrt(3) - dx / sqrt(3))) && (dy<(a*sqrt(3) + dx / sqrt(3)))
				&& (dy>(-a*sqrt(3) - dx / sqrt(3))) && (dy > (-a*sqrt(3) + dx / sqrt(3))))
				//condition of x,y being in the hexagonal with the center at (x_app_center y_app_center)
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
		*x_i = index_x;
		*y_i = index_y;
	}
}

//processes leaving the cell and transporting back to parent
G4ThreeVector hex_on_leave (G4ThreeVector GlobalParentPos /*center*/, G4double x_par, G4double y_par,
	G4ThreeVector GlobalCellPos, G4double a, G4int *curr_x_id, G4int *curr_y_id, const G4Step& step,
	G4double active_x_par, G4double active_y_par)
{
	G4ThreeVector pos= step.GetPostStepPoint()->GetPosition();
	G4double dz = (pos - GlobalCellPos).z();
	G4double x_center = (*curr_x_id - 1) * 2 * a + a + (*curr_y_id - 2 * ((int)(*curr_y_id / 2.0) == 0)) ? a : 0; //x center depends on index of y
	G4double y_center = (*curr_y_id - 1)*a*sqrt(3) + 2 * a / sqrt(3);
	G4double x_rel = pos.x() - GlobalCellPos.x() + x_center;
	G4double y_rel = pos.y() - GlobalCellPos.y() + y_center;
	x_rel = x_rel -active_x_par/ 2;
	y_rel = y_rel - active_y_par/ 2;
	G4ThreeVector new_pos(x_rel, y_rel, dz);
	new_pos = new_pos + GlobalParentPos;
	*curr_x_id = -1;
	*curr_y_id = -1;
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
//2) acoount for transformations (rotations) (and make them via separate methods)
G4ThreeVector hexagonal_mapping(G4ThreeVector GlobalParentPos /*center*/, G4double x_par, G4double y_par, G4double z_par,
	G4ThreeVector GlobalCellPos, G4double cell_size, G4int *curr_x_id, G4int *curr_y_id,
	const G4Track& aTrack, const G4Step& aStep)
{
	G4ThreeVector ret = aStep.GetPostStepPoint()->GetPosition();
	G4int X_max_index = int (x_par/(2*cell_size));
	G4int Y_max_index = int((y_par - cell_size / (sqrt(3))) / (cell_size*sqrt(3)));
	if (X_max_index <= 0 || Y_max_index <= 0) return ret; //no grid - no changes
	//I have all indices starting from 1
	G4double active_x_par = X_max_index * 2 * cell_size;  //only this length is mapped
	G4double active_y_par = Y_max_index * cell_size*sqrt(3)+cell_size/sqrt(3);
	if ((*curr_x_id < 0)||(*curr_y_id<0))//photon from outside
	{
		G4ThreeVector pos = aStep.GetPostStepPoint()->GetPosition();
		G4ThreeVector dir = aStep.GetPostStepPoint()->GetMomentumDirection();
		G4double dz = (pos - GlobalParentPos).z();
		G4double tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
		if ((abs((abs(dz) - z_par / 2)) < tolerance) && (G4ThreeVector(0, 0, dz)*dir < 0)) //on the surface and moves inside
		{
			G4double x_rel = pos.x() - GlobalParentPos.x() + x_par / 2 - (x_par - active_x_par) / 2; //x_rel starts form the edge of mapped plate 
			//(mapped not entire one but the part which contains integer amount of cells (hence offset (x_par - active_x_par) / 2))
			G4double y_rel = pos.y() - GlobalParentPos.y() + y_par / 2 - (y_par - active_x_par) / 2; //should be granted to be positive
			if ((x_rel < 0) || (y_rel < 0))
				return ret;
			if ((x_rel>active_x_par) || (y_rel > active_y_par))
				return ret;
			G4int x_i, y_i;
			hex_indices_by_pos(x_rel, y_rel, cell_size, X_max_index, Y_max_index, &x_i, &y_i);
			if ((x_i < 0) || (y_i < 0))
				return ret;
			*curr_x_id = x_i;
			*curr_y_id = y_i;

			G4double a = cell_size;
			G4double x_center = (x_i - 1) * 2 * a + a + ((y_i - 2 * ((int)(y_i / 2.0)) == 0) ? a : 0); //x center depends on index of y
			G4double y_center = (y_i - 1)*a*sqrt(3) + 2 * a / sqrt(3);

			G4ThreeVector new_pos(x_rel - x_center, y_rel - y_center, dz);
			new_pos = new_pos + GlobalCellPos;
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
		G4double dz = (pos - GlobalCellPos).z();
		G4double tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
		if ((abs((abs(dz) - z_par / 2)) < tolerance) && (G4ThreeVector(0, 0, dz)*dir > 0)) //on the surface and moves outside
			return hex_on_leave(GlobalParentPos, x_par, y_par,
				GlobalCellPos, cell_size, curr_x_id, curr_y_id,aStep,active_x_par, active_y_par);
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
		G4double x_rel = pos.x() - GlobalCellPos.x();
		G4double y_rel = pos.y() - GlobalCellPos.y();
		G4int y_index_odd = *curr_y_id - 2 * ((int)(*curr_y_id / 2.0));
		if ((abs(cell_size - x_rel)<tolerance)&&(dir.x()>0)) //leaving x+ edge
		{
			*curr_x_id++;
			if (*curr_x_id > X_max_index)
			{//leave cell and go to parent
				return hex_on_leave(GlobalParentPos, x_par, y_par,
					GlobalCellPos, cell_size, curr_x_id, curr_y_id, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() - 2 * cell_size, pos.y(), pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x+. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(cell_size + x_rel)<tolerance)&&(dir.x()<0)) //x- edge
		{
			*curr_x_id--;
			if (*curr_x_id < 1)
			{//leave cell and go to parent
				return hex_on_leave(GlobalParentPos, x_par, y_par,
					GlobalCellPos, cell_size, curr_x_id, curr_y_id, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() + 2 * cell_size, pos.y(), pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x-. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel - cell_size*sqrt(3) +x_rel / sqrt(3))<tolerance)&&(dir*G4ThreeVector(0.5,sqrt(3)/2,0)>0)) //x+,y+
		{
			*curr_x_id+=(y_index_odd?0:1);
			*curr_y_id++;
			if ((*curr_y_id>Y_max_index) || (*curr_x_id >X_max_index))
			{//leave cell and go to parent
				return hex_on_leave(GlobalParentPos, x_par, y_par,
					GlobalCellPos, cell_size, curr_x_id, curr_y_id, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() - 2 * cell_size*0.5, pos.y() - 2 * cell_size*sqrt(3) / 2, pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x+y+. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel - cell_size*sqrt(3) - x_rel / sqrt(3))<tolerance) && (dir*G4ThreeVector(-0.5, sqrt(3) / 2, 0)>0)) //x-,y+
		{
			*curr_x_id -= (y_index_odd ? 1 : 0); //- picture required - this is the way numeration works
			*curr_y_id++;
			if ((*curr_y_id>Y_max_index) || (*curr_x_id <1))
			{//leave cell and go to parent
				return hex_on_leave(GlobalParentPos, x_par, y_par,
					GlobalCellPos, cell_size, curr_x_id, curr_y_id, aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() + 2 * cell_size*0.5, pos.y() - 2 * cell_size*sqrt(3) / 2, pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x-y+. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel + cell_size*sqrt(3) + x_rel / sqrt(3))<tolerance) && (dir*G4ThreeVector(-0.5, -sqrt(3) / 2, 0)>0)) //x-,y-
		{
			*curr_x_id -= (y_index_odd ? 1 : 0); //- picture required - this is the way numeration works
			*curr_y_id--;
			if ((*curr_y_id<1) || (*curr_x_id <1))
			{//leave cell and go to parent
				return hex_on_leave(GlobalParentPos, x_par, y_par,
					GlobalCellPos, cell_size, curr_x_id, curr_y_id,  aStep, active_x_par, active_y_par);
			}
			pos = G4ThreeVector(pos.x() + 2 * cell_size*0.5, pos.y() + 2 * cell_size*sqrt(3) / 2, pos.z());
#ifdef DEBUG_MC_NODES
			G4cout << "MAPPING: Cell to cell map x-y-. Proposed position: " << pos << G4endl;
#endif
			return pos;
		}
		if ((abs(y_rel + cell_size*sqrt(3) - x_rel / sqrt(3))<tolerance) && (dir*G4ThreeVector(0.5, -sqrt(3) / 2, 0)>0)) //x+,y-
		{
			*curr_x_id += (y_index_odd ? 0 : 1); //- picture required - this is the way numeration works
			*curr_y_id--;
			if ((*curr_y_id<1) || (*curr_x_id <1))
			{//leave cell and go to parent
				return hex_on_leave(GlobalParentPos, x_par, y_par,
					GlobalCellPos, cell_size, curr_x_id, curr_y_id, aStep, active_x_par, active_y_par);
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

PseudoMesh::PseudoMesh(G4LogicalVolume* parent_vol, G4ThreeVector g_par_pos, G4double par_x, G4double par_y, G4double par_z,
	G4LogicalVolume* cell_vol, G4ThreeVector g_cell_pos, G4double cell_parameter, box_map_function function)
{
	map_function = function;
	parent = parent_vol;
	g_parent_pos = g_par_pos;
	par_x_size = par_x;
	par_y_size = par_y;
	par_z_size = par_z;
	cell = cell_vol;
	cell_size = cell_parameter;
	g_cell_posit = g_cell_pos;

	curr_x_id = new G4int(-1);
	curr_y_id = new G4int (-1);
}

PseudoMesh::~PseudoMesh()
{
	delete curr_x_id;
	delete curr_y_id;
}

//returns global position
G4ThreeVector PseudoMesh::PostSteppingAction(const G4Track& aTrack, const G4Step& aStep, G4TouchableHandle &post_step_handle) //checks entering/leaving, called from DetectorConstruction, manages everything
{
	G4LogicalVolume* postV = post_step_handle->GetVolume()->GetLogicalVolume();
	G4LogicalVolume* preV = aStep.GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	// TODO? change to (postV == parent) || (postV == cell) || (preV == cell) || (preV == parent) ?
	if ((aStep.GetPostStepPoint()->GetStepStatus()==fGeomBoundary)&&((postV == parent) || (preV == cell)))
		return (*map_function)(g_parent_pos, par_x_size, par_y_size, par_z_size, g_cell_posit, cell_size, curr_x_id, curr_y_id,
		aTrack, aStep);
	return aStep.GetPostStepPoint()->GetPosition();
}
