#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4VisExecutive.hh"
#include "G4VisExtent.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4TransportationManager.hh"

#define _light_blue_RGB 0.5, 0.5, 0.9
#define _red_RGB 0.9, 0.4, 0.4

B1DetectorConstruction::B1DetectorConstruction()
	: G4VUserDetectorConstruction()
{
	x_rot_m = NULL;
	y_rot_m = NULL;
	top_GEM = NULL;
}

B1DetectorConstruction::~B1DetectorConstruction()
{ 
	if (x_rot_m) delete x_rot_m;
	if (y_rot_m) delete y_rot_m;
	if (top_GEM) delete top_GEM;
}

G4ThreeVector B1DetectorConstruction::get_global_normal(G4StepPoint* point, G4int *validity)
{
	G4bool valid = 1;
	G4ThreeVector position = point->GetPosition();
	//obtaining of the global normal to entered volume (copied from boundary process)
	G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
	G4Navigator* Nav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	G4ThreeVector theGlobalNormal = Nav->GetLocalExitNormal(&valid);
	*validity = valid;
	if (valid)
	{
		theGlobalNormal = Nav->GetLocalToGlobalTransform().TransformAxis(theGlobalNormal);
	}
	else
	{
		G4ExceptionDescription ed;
		ed << " G4DetectorConstruction->GetHitProbability(): "
			<< " The Navigator reports that it returned an invalid normal"
			<< G4endl;
		G4Exception("G4DetectorConstruction->GetHitProbability()", "Det_contsr",
			EventMustBeAborted, ed,
			"Invalid Surface Normal - Geometry must return valid surface normal");
		return theGlobalNormal;
	}
	if (theGlobalNormal*theGlobalNormal != 0)
	{
		return theGlobalNormal;
	}
	*validity = 0;
	return theGlobalNormal;
}

G4double B1DetectorConstruction::GetHitProbability(G4StepPoint* post_point)
{
	G4VPhysicalVolume* tremor = post_point->GetTouchableHandle()->GetVolume();
	if (NULL == tremor)
	{// kind of error
		return -1;
	}
	G4LogicalVolume* in_volume = tremor->GetLogicalVolume();
	G4ThreeVector position = post_point->GetPosition();
	G4ThreeVector Momentum_dir = post_point->GetMomentumDirection();
	G4double particle_energy = post_point->GetTotalEnergy();
	const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();
	G4int is_on_boundary = 0;
	if (hStep)
		is_on_boundary = (hStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
	else
		is_on_boundary = (post_point->GetStepStatus() == fGeomBoundary);
	for (auto j = detector_volumes.begin(); j != detector_volumes.end(); j++)
	{
		if (*j == in_volume)
		{
//#ifdef TEMP_CODE_
//			return 1;
//#endif
			if (is_on_boundary)
			{
				G4int is_valid = 1;
				G4ThreeVector theGlobalNormal = get_global_normal(post_point, &is_valid);
				if (is_valid)
					if (theGlobalNormal*Momentum_dir < 0)
						return -1;
			}
			G4Material* mat_= in_volume->GetMaterial();
			if (mat_)
			{
				G4MaterialPropertiesTable* mat_t = mat_->GetMaterialPropertiesTable();
				if (mat_t)
				{
					G4MaterialPropertyVector* qe = mat_t->GetProperty("PMT_QE");
					if (qe)
					{
						G4double prob = qe->Value(particle_energy);
						if (prob < 0) return 0;
						else
							return prob;
					}
				}
			}
			return 1; //^DONE? TODO: modify considering momentum and position to avoid nesting detector volume in the larger volume of the same properties? 
		}
	}
	for (auto j = absorbtion_volumes.begin(); j != absorbtion_volumes.end(); j++)
	{
		if (*j == in_volume)
		{
			if (is_on_boundary)
			{
				G4int is_valid = 1;
				G4ThreeVector theGlobalNormal = get_global_normal(post_point, &is_valid);
				if (is_valid)
					if (theGlobalNormal*Momentum_dir < 0)
						return -1;
			}
			return 0; //kills event
		}
	}
	return -1;
}
void B1DetectorConstruction::wavelen_transmittance_to_table
(G4MaterialPropertiesTable* table, const G4double* wl, const G4double* tr, G4int size, G4double width)
{
	G4double* energies = new G4double [size];
	G4double* abs_lengths = new G4double[size];
	for (G4int counter = 0; counter < size; counter++)
	{
		energies[size - 1 - counter] = 1 * eV*(1.2398e3*nm / wl[counter]);//maintains same order (i.e. increasing)
		if (0 == tr[counter])
		{
			abs_lengths[size - 1 - counter] = 0;
		}
		else
		{
			if (1 != tr[counter])
				abs_lengths[size - 1 - counter] = -width / (log(tr[counter]));
			else
				abs_lengths[size - 1 - counter] = DBL_MAX; //no absorbtion
		}
#ifdef TEMP_CODE_
		int b = abs_lengths[size - 1 - counter];
		int a = b;
#endif
	}
	table->AddProperty("ABSORBTION_LENGTH", energies, abs_lengths, size);
	delete[] energies;
	delete[] abs_lengths;
}

//TODO: add reading ready integral spec
void B1DetectorConstruction::gen_integral_en_spec
(G4String file, G4MaterialPropertiesTable* table)
{
	std::list<G4double> wavelengths;
	std::list<G4double> probs;
	G4double min_x,max_x;
	std::ifstream fl(file);
	while (!fl.eof())
	{
		G4double x, y;
		fl >> x;
		if (fl.eof())
			break;
		fl >> y;
		wavelengths.push_back(x);
		probs.push_back(y);
	}
	fl.close();
	if (wavelengths.empty())
		return;
	min_x = *wavelengths.begin(); //sorted input assumed
	max_x = wavelengths.back();
	G4MaterialPropertiesTable* temp_prop_table = new G4MaterialPropertiesTable();
	G4double* temp_wavelengths = new G4double[wavelengths.size()];
	G4double* temp_weights = new G4double[wavelengths.size()];

	int fgg = 0;
	auto j = probs.begin();
	for (auto u = wavelengths.begin(); u != wavelengths.end(); fgg++,j++,u++)
	{
		temp_wavelengths[fgg] = *u;
		temp_weights[fgg] = *j;
	}
	G4MaterialPropertyVector file_vals = G4MaterialPropertyVector(temp_wavelengths, temp_weights, wavelengths.size());
	G4double* i_wavelengths = new G4double[SPEC_INTEGRATION_STEPS+1];
	G4double* i_probab = new G4double[SPEC_INTEGRATION_STEPS+1];
	wavelengths.erase(wavelengths.begin(),wavelengths.end());
	probs.erase(probs.begin(), probs.end());
//#ifdef TEMP_CODE_
//	std::ofstream diff_o(TEST_WLS_OUT_SPEC);
//	for (G4int counter = 0; counter <= SPEC_INTEGRATION_STEPS; counter++)
//	{
//		diff_o << min_x + counter*(max_x - min_x) / SPEC_INTEGRATION_STEPS
//			<< "\t" << file_vals.Value(min_x + counter*(max_x - min_x) / SPEC_INTEGRATION_STEPS) << G4endl;
//	}
//	diff_o.close();
//#endif //TEMP_CODE_
	G4double int_value=0;
	G4double delta_x = (max_x - min_x) / SPEC_INTEGRATION_STEPS;
	G4double point, point_left, point_right;
	G4double val_left, val_center, val_right;
	for (G4int counter = 0; counter <= SPEC_INTEGRATION_STEPS; counter++)
	{
		point = min_x + counter*delta_x;
		point_left = min_x + counter*delta_x - delta_x/2;
		point_right = min_x + counter*delta_x + delta_x / 2;
		point_left = point_left<min_x ? min_x : point_left;
		point_left = point_left>max_x ? max_x : point_left;
		
		val_left = file_vals.Value(point_left);
		val_center = file_vals.Value(point);
		val_right =file_vals.Value(point_right);
		if (val_left < 0) val_left = 0;
		if (val_center < 0) val_center = 0;
		if (val_right < 0) val_right = 0;
		int_value += delta_x*(val_center+(val_right-val_left)/8+(val_left+val_right-2*val_center)/12); //integration
		//for (G4int tt = 0; tt <= counter; tt++)
		//{
		//	G4double val = file_vals.Value(min_x + tt*(max_x - min_x) / SPEC_INTEGRATION_STEPS);
		//	if (val < 0) val = 0;
		//	int_value +=val* ((max_x - min_x) / SPEC_INTEGRATION_STEPS);
		//}
		i_wavelengths[counter]=(min_x+counter*delta_x);
		i_probab[counter]=int_value;
	}
	G4double integral = i_probab[SPEC_INTEGRATION_STEPS];
	for (G4int counter = 0; counter <= SPEC_INTEGRATION_STEPS; counter++)
	{
		i_probab[counter] = i_probab[counter] / integral;
		i_wavelengths[counter] = 1 * eV*(1.2398e3/**nm - file in nm already*/ / i_wavelengths[counter]);//become of reverse order (it IS OK)
	}
	G4MaterialPropertyVector *integral_spec =new G4MaterialPropertyVector(i_probab, i_wavelengths, SPEC_INTEGRATION_STEPS+1);
	table->AddProperty("WLS_ENERGY_SPECTRUM", integral_spec);
//#ifdef TEMP_CODE_
//	diff_o.open(TEST_WLS_OUT_I_SPEC,std::ios_base::trunc);
//	G4double probab;
//	for (G4int counter = 0; counter <= SPEC_INTEGRATION_STEPS; counter++)
//	{
//		probab = 1.0*counter / SPEC_INTEGRATION_STEPS;
//		diff_o << 1.2398e3*eV/integral_spec->Value(probab)/*/(1*eV)*/<< "\t" << probab<< G4endl;
//	}
//	diff_o.close();
//
//	diff_o.open(TEST_OUT_I_SPEC1, std::ios_base::trunc);
//	G4double diff;
//	for (G4int counter = 0; counter <= SPEC_INTEGRATION_STEPS; counter++)
//	{
//		//probab = (double)counter / SPEC_INTEGRATION_STEPS;
//		diff = (counter == SPEC_INTEGRATION_STEPS ? 0 : -(i_probab[counter] - i_probab[counter + 1]) /
//			(integral_spec->Value(i_probab[counter]) - integral_spec->Value(i_probab[counter + 1])));
//		diff_o << 1.2398e3*eV / integral_spec->Value(i_probab[counter])/*/(1*eV)*/ << "\t" << diff << G4endl;
//	}
//	diff_o.close();
//#endif //TEMP_CODE_
	delete[] temp_wavelengths;
	delete[] temp_weights;
	delete[] i_wavelengths;
	delete[] i_probab;
}

void B1DetectorConstruction::gen_PMT_QE(G4String file, G4MaterialPropertiesTable* table)
{
	std::list<G4double> wavelengths;
	std::list<G4double> probs;
	std::ifstream fl(file);
	while (!fl.eof())
	{
		G4double x, y;
		fl >> x;
		if (fl.eof())
			break;
		fl >> y;
		wavelengths.push_back(x);
		probs.push_back(y);
	}
	fl.close();
	if (wavelengths.empty())
		return;
	G4MaterialPropertiesTable* temp_prop_table = new G4MaterialPropertiesTable();
	G4double* _energies = new G4double[wavelengths.size()];
	G4double* _probs = new G4double[wavelengths.size()];

	int fgg = 0;
	auto j = probs.rbegin();
	for (auto u = wavelengths.rbegin(); u != wavelengths.rend(); fgg++, ++j, ++u)
	{
		_energies[fgg] = 1 * eV*(1.2398e3/**nm - file in nm already*/ / *u);
		_probs[fgg] = *j;
	}
	G4MaterialPropertyVector *qe_spec = new G4MaterialPropertyVector(_energies, _probs, wavelengths.size());
	table->AddProperty("PMT_QE", qe_spec);
	delete[] _energies;
	delete[] _probs;
}

G4Material* B1DetectorConstruction::_Argon_mat(void)
{
	G4Material* _Argon = new G4Material("Argon",18,39.95*g/mole, 1.784e-3*273/87*g/cm3,kStateGas,87*kelvin,1*atmosphere);
	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();
	G4double energies[8] = { 2.72 * eV, 2.95 * eV, 3.05*eV, 3.26*eV, 3.46*eV, 3.68*eV, 3.92*eV,9.68*eV};
	G4double rindexes[8] = { 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01};
	prop_table->AddProperty("RINDEX", energies, rindexes, 8);
	//TODO: gen_integral_en_spec(ARGON_SPECTRUM_FILE, prop_table); //emission parameters, used in CustomRunManger
	_Argon->SetMaterialPropertiesTable(prop_table);
	return _Argon;
}

G4Material* B1DetectorConstruction::_WLS_mat(void)
{
	G4Element* el_C = new G4Element("Carbon", "C", 6, 12.0*g/mole);
	G4Element* el_H = new G4Element("Hydrogen", "H", 1, 1*g/mole);
	G4Material* _WLS = new G4Material("WLS",1* g / cm3,2,kStateSolid,87*kelvin,1*atmosphere);
	_WLS->AddElement(el_C, 8);
	_WLS->AddElement(el_H, 8);
	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();

	G4double energies[4] = {2.72*eV, 2.95 * eV, 3.68*eV, 9.68 * eV }; //TODO: more data required
	G4double rindexes[4] = { 1.62, 1.62, 1.62, 1.60}; //TODO: figure out n complexity here

	prop_table->AddProperty("RINDEX", energies, rindexes, 4);

	G4double wavelen[22] = {340*nm, 350*nm,    360*nm, 362*nm, 366*nm, 369*nm, 372*nm, 375*nm, 378*nm, 381*nm, 384*nm, 387*nm, 390*nm, 
		393*nm, 396*nm, 400*nm, 403*nm, 408*nm, 412*nm, 415.8*nm, 420*nm, 424*nm};
	G4double transmittance[22] = {1e-6, 1e-4, 0.0048, 0.0096, 0.0120, 0.0168, 0.0229, 0.0361, 0.0614, 0.1024, 0.1697, 0.2557, 0.3600,
		0.5633, 0.6962, 0.7872, 0.8438, 0.8708, 0.8844, 0.8893, 0.8967, 0.8979};

	wavelen_transmittance_to_table(prop_table, wavelen, transmittance, 22,WLS_FILM_WIDTH);

	G4double emiss_en[8] = { 2.72 * eV, 2.95 * eV, 3.05*eV, 3.26*eV, 3.46*eV, 3.68*eV, 3.92*eV, 9.68*eV };
	G4double emiss_eff[8] = { 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.58};
	prop_table->AddProperty("WLS_EFFICIENCY", emiss_en, emiss_eff, 8);

	gen_integral_en_spec(WLS_SPECTRUM_FILE,prop_table);

	_WLS->SetMaterialPropertiesTable(prop_table);
	return _WLS;
}

G4Material*  B1DetectorConstruction::_Copper_mat(void)
{
	G4Material* _Copper = new G4Material("Copper", 29, 63.55*g / mole, 8.96* g / cm3, kStateSolid, 87 * kelvin, 1 * atmosphere);
	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();
	G4double energies[8] = { 2.72 * eV, 2.95 * eV, 3.05*eV, 3.26*eV, 3.46*eV, 3.68*eV, 3.92*eV, 9.68*eV };
	G4double rindexes[8] = { 1, 1, 1, 1, 1, 1, 1, 1}; //its optical properties are set via OpticalSurface
	prop_table->AddProperty("RINDEX", energies, rindexes, 8);
	_Copper->SetMaterialPropertiesTable(prop_table);
	return _Copper;
}

void B1DetectorConstruction::_SetCopperSurface(G4OpticalSurface* surface)
{
	surface->SetType(dielectric_metal);
	surface->SetFinish(polished);
	surface->SetModel(glisur);
	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();
	G4double en_index[17]={2.63*eV,2.75*eV,2.88*eV,3.37*eV,3.87*eV,4.12*eV,4.24*eV,4.36*eV,4.49*eV,4.61*eV,4.74*eV,4.86*eV,5.11*eV,5.60*eV,5.98*eV,6.60*eV,9.68*eV};
	G4double r_rindex[17] ={1.25,	1.24,	1.25,	1.36,	1.38,	1.40,	1.42,	1.45,	1.46,	1.45,	1.41,	1.41,	1.34,	1.13,	1.01,	0.94,	1.12};
	G4double i_rindex[17] ={2.483,	2.397,	2.305,	1.975,	1.783,	1.679,	1.633,	1.633,	1.646,	1.668,	1.691,	1.741,	1.799,	1.737,	1.599,	1.337,	0.7284};
	//^Phys. Rev. Vol.6, Number 12, P.B. Jonson and R. W. Christy "Optical Constants of the Noble Metals", table 1 p.4347
	prop_table->AddProperty("REALRINDEX", en_index, r_rindex, 17);
	prop_table->AddProperty("IMAGINARYRINDEX", en_index, i_rindex, 17);
	surface->SetMaterialPropertiesTable(prop_table);
}

G4Material* B1DetectorConstruction::_Acrylic_mat(void)
{
	G4Element* el_C = new G4Element("Carbon", "C", 6, 12*g/mole);
	G4Element* el_H = new G4Element("Hydrogen", "H", 1, 1 * g / mole);
	G4Element* el_O = new G4Element("Oxygen", "O", 8, 16 * g / mole);
	G4Material* _Acrylic = new G4Material("PMMA (Acrylic)", 1.18 * g / cm3, 3, kStateSolid, 87 * kelvin, 1 * atmosphere);
	_Acrylic->AddElement(el_C, 5);
	_Acrylic->AddElement(el_H, 8);
	_Acrylic->AddElement(el_O, 2);
	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();
	G4double energies[6] = { 1.55*eV, 2.817*eV, 2.95*eV, 3.10*eV, 3.54*eV, 4.13*eV};
	G4double rindexes[6] = { 1.485,   1.501,    1.503,   1.47,    1.46,    1.46}; //TODO: confirm the data

	prop_table->AddProperty("RINDEX", energies, rindexes, 6);
	G4double wavelen[26] = { 320*nm, 336 * nm, 343.8 * nm, 350 * nm, 360 * nm, 363.5 * nm, 366.4 * nm, 367 * nm,
		371.9 * nm, 373.5 * nm, 376 * nm, 378 * nm, 380 * nm, 382.2*nm, 385.7*nm, 387.3*nm, 390.6*nm, 
		394.4*nm, 398.6*nm, 403.4*nm, 414.1*nm, 427.6*nm, 444 * nm, 475.3*nm, 501.1*nm, 229.8*nm};
	G4double transmittance[26] = {5e-6, 1e-4, 0.015,0.023,0.067, 0.115,0.182,0.270, 0.356, 0.417, 0.5,0.586, 0.652,
		0.717, 0.782, 0.829, 0.857, 0.874, 0.886,0.892,0.894,0.897,0.9,0.902,0.905,0.908};
	wavelen_transmittance_to_table(prop_table, wavelen, transmittance, 26, PMMA_WIDTH);
	G4double emiss_en[8] = { 2.72 * eV, 2.95 * eV, 3.05*eV, 3.26*eV, 3.46*eV, 3.68*eV, 3.92*eV, 9.68*eV };
	G4double emiss_eff[8] = {0,0,0,0,0,0,0,0}; //just absortion, no reemission
	prop_table->AddProperty("WLS_EFFICIENCY", emiss_en, emiss_eff, 8);
	_Acrylic->SetMaterialPropertiesTable(prop_table);
	return _Acrylic;
}

G4Material* B1DetectorConstruction::_LArgon_mat(void)
{
	G4Material* _Argon = new G4Material("Liquid Argon", 18, 39.95*g / mole, 1.4* g / cm3, kStateLiquid, 87 * kelvin, 1 * atmosphere);
	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();
	G4double energies[10] = { 2.95 * eV, 3.2*eV, 3.4*eV, 3.8*eV, 4*eV, 5.63*eV, 6.89*eV, 7.75*eV, 8.86*eV, 9.68*eV };
	G4double rindexes[10] = { 1.23,		 1.23,	 1.23,	 1.23,	 1.23, 1.26,	1.29,	 1.31,	  1.34,    1.38 };
	/*n===+0.03 of Neumeier, "Optical Properties of Liquid Noble
	Gas Scintillators" - corrected for logbook data*/
	prop_table->AddProperty("RINDEX", energies, rindexes, 10);
	_Argon->SetMaterialPropertiesTable(prop_table);
	return _Argon;
}

G4Material* B1DetectorConstruction::_FusedSilica_mat(void)
{
	G4Element* el_Si = new G4Element("Silicon", "Si", 14, 28 * g / mole);
	G4Element* el_O = new G4Element("Oxygen", "O", 8, 16 * g / mole);
	G4Material* _Window = new G4Material("PMMA (Acrylic)", 2.203* g / cm3, 2, kStateSolid, 87 * kelvin, 1 * atmosphere);
	_Window->AddElement(el_Si, 1);
	_Window->AddElement(el_O, 2);
	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();
	G4double energies[10] = { 2.72 * eV, 2.95*eV, 3.06*eV, 3.4*eV, 4.11*eV, 4.43*eV, 5.79*eV, 6.2*eV, 6.7*eV, 7.29*eV };
	G4double rindexes[10] = { 1.465,	 1.468,	  1.469,   1.474,  1.487,	1.494,	 1.53,	  1.55,	  1.575,  1.615}; 
	//^most taken at room temperature: error is about 1e-3
	prop_table->AddProperty("RINDEX", energies, rindexes, 10);
	gen_PMT_QE(PMT_QE_FILE, prop_table);
	_Window->SetMaterialPropertiesTable(prop_table);
	return _Window;
}

G4Material* B1DetectorConstruction::_fr4_mat(void)
{
	G4Element* el_C = new G4Element("Carbon", "C", 6, 12 * g / mole);
	G4Element* el_O = new G4Element("Oxygen", "O", 8, 16 * g / mole);
	G4Element* el_H = new G4Element("Hydrogen", "H", 1, 1 * g / mole);
	G4Material* _fr4 = new G4Material("PMMA (Acrylic)", 1.85 * g / cm3, 3, kStateSolid, 87 * kelvin, 1 * atmosphere);
	_fr4->AddElement(el_C, 5);
	_fr4->AddElement(el_H, 8);
	_fr4->AddElement(el_O, 2);

	G4MaterialPropertiesTable* prop_table = new G4MaterialPropertiesTable();
	G4double energies[10] = { 2.72 * eV, 2.95*eV, 3.06*eV, 3.4*eV, 4.11*eV, 4.43*eV, 5.79*eV, 6.2*eV, 6.7*eV, 7.29*eV };
	G4double rindexes[10] = { 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
	//^most taken at room temperature: error is about 1e-3
	prop_table->AddProperty("RINDEX", energies, rindexes, 10);
	return _fr4;
}

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* WLS_mat = _WLS_mat();
  G4Material* Acrylic_mat= _Acrylic_mat();
  G4Material* Copper_mat = _Copper_mat();
  G4Material* LAr_mat = _LArgon_mat();
  G4Material* Ar_mat = _Argon_mat();
  G4Material* FusedSilica_mat = _FusedSilica_mat();
  G4Material* FR4_mat = _fr4_mat();

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;
  //     
  // World
  //
#define plate_W 0.6
#define plate_real_W 0.5
  G4double world_sizeXY = (141+2*(4+2+1+1)) * mm+2*(WLS_FILM_WIDTH+PMMA_WIDTH);
  G4double world_sizeZ  = (20+plate_W+18+4+plate_W+20+2*(1+1))*mm;
  G4double envelope_sizeXY = (141 + 2 * (4 + 2 + 1)) * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH);
  G4double envelope_sizeZ = (20 + plate_W + 18 + 4 + plate_W + 20 + 2 * (1))*mm;

  G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
  world = new G4LogicalVolume(solidWorld, Ar_mat,"World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), world, "World", 0, false, 0, checkOverlaps);
  
  //Envelope - for photons killing upon entering to world
  G4Box* solidEnv = new G4Box("Envelope", 0.5*envelope_sizeXY, 0.5*envelope_sizeXY, 0.5*envelope_sizeZ);
  envelope = new G4LogicalVolume(solidEnv, Ar_mat, "Envelope");
  G4VPhysicalVolume* physEnv = new G4PVPlacement(0, G4ThreeVector(), envelope, "Envelope", world, false, 0, checkOverlaps);

  G4double box_sizeXY = (141 + 2 * (2)) * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH); 
  //^includes gap between PMT (window) and acrylic. Liquid in that gap is placed here 
  G4double box_sizeZ = (20 + plate_W + 18 + 4 + plate_W + 20)*mm;
  G4Box* solid_box = new G4Box("MainBox", 0.5*box_sizeXY, 0.5*box_sizeXY, 0.5 *box_sizeZ); //volume exluding detectors
  box = new G4LogicalVolume(solid_box, Ar_mat, "MainBox");
  G4VPhysicalVolume* phys_box = new G4PVPlacement(0, G4ThreeVector(0,0,0), box, "MainBox", envelope, false, 0, checkOverlaps);

  G4double box_interior_sizeXY = (141) * mm;
  G4double box_interior_sizeZ = (plate_W + 18 + 4 + plate_W)*mm;
  G4Box* solid_box_interior = new G4Box("MainBoxInterior", 0.5*box_interior_sizeXY, 0.5*box_interior_sizeXY, 0.5 *box_interior_sizeZ);
  //volume contaning GEMs and a bit of liquid argon
  box_interior = new G4LogicalVolume(solid_box_interior, Ar_mat, "MainBoxInterior");
  G4VPhysicalVolume* phys_box_interior = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), box_interior, 
	  "MainBoxInterior", box, false, 0, checkOverlaps);

  G4Box* solid_top_pseudo_GEM = new G4Box("TopPlate", 0.5*(141 * mm), 0.5*(141 * mm), 0.5 *(plate_W*mm)); //mapped to cells
  G4Box* solid_bot_pseudo_GEM = new G4Box("BotPlate", 0.5*(141 * mm), 0.5*(141 * mm), 0.5 *(plate_W*mm)); //mapped to cells

  G4Box* solid_top_GEM = new G4Box("TopGEM", 0.5*(141*mm), 0.5*(141*mm), 0.5 *(plate_real_W*mm)); //placed in above pseudo volumes
  G4Box* solid_bot_GEM = new G4Box("BotGEM", 0.5*(141 * mm), 0.5*(141 * mm), 0.5 *(plate_real_W*mm));
  G4Box* solid_LAr_layer = new G4Box("LArLayer", 0.5*(141 * mm), 0.5*(141 * mm), 0.5 *(4*mm));
  top_ps_plate= new G4LogicalVolume(solid_top_pseudo_GEM, Ar_mat, "TopPseudoGEM");
  bot_ps_plate = new G4LogicalVolume(solid_bot_pseudo_GEM, LAr_mat, "BotPseudoGEM");
  top_cu_plate = new G4LogicalVolume(solid_top_GEM, Copper_mat, "TopGEM");
  bot_cu_plate = new G4LogicalVolume(solid_bot_GEM, Copper_mat, "BotGEM");
  LAr_layer = new G4LogicalVolume(solid_LAr_layer, LAr_mat, "LArLayer");
  G4VPhysicalVolume* phys_top_pseudo_GEM = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(18 + 4 + plate_W)*mm), top_ps_plate,
	  "TopPseudoGEM", box_interior, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_bot_pseudo_GEM = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*(18 + 4 + plate_W)*mm), bot_ps_plate,
	  "BotPseudoGEM", box_interior, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_LAr_layer = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*(18)*mm), LAr_layer, "LArLayer", box_interior,
	  false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_top_CU = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), top_cu_plate, //placed in above pseudo volumes
	  "TopCU", top_ps_plate, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_bot_CU = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), bot_ps_plate,
	  "BotCU", bot_ps_plate, false, 0, checkOverlaps);
	//================================ CELLs
#define CELL_SIZE 0.45
#define CELL_HOLE_DIAM 0.5
#define CELL_RIM_W 0.1
  G4double r_[2] = {CELL_SIZE+0.1,CELL_SIZE+0.1};
  G4double z_[2] = { -plate_W / 2, plate_W/2 };
  G4Polyhedra* solid_top_cell_container = new G4Polyhedra("top_cell_container", pi/6.0, twopi, 6, 2,r_,z_);
  top_cell_container = new G4LogicalVolume (solid_top_cell_container,Ar_mat,"top_cell_container");
  G4VPhysicalVolume* phys_top_cell_container = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(18 + 4 + plate_W+12)*mm), top_cell_container,
	  "TopCellContainer", box, false, 0, checkOverlaps);
  r_[0] = CELL_SIZE; r_[1] = CELL_SIZE;
  z_[0] = -plate_real_W / 2; z_[0] = plate_real_W / 2;
  G4Polyhedra* solid_top_cell= new G4Polyhedra("top_cell", pi / 6.0, twopi, 6, 2, r_, z_);
  top_cell = new G4LogicalVolume(solid_top_cell, FR4_mat, "top_cell");
  G4VPhysicalVolume* phys_top_cell = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), top_cell,
	  "TopCellContainer", top_cell_container, false, 0, checkOverlaps);
  G4double hole_rs[3] = { CELL_HOLE_DIAM / 2 + CELL_RIM_W, CELL_HOLE_DIAM / 2, CELL_HOLE_DIAM / 2 + CELL_RIM_W };
  G4double hole_zs[3] = { -plate_real_W / 2, 0, plate_real_W / 2 };
  G4Polycone* solid_top_cell_hole = new G4Polycone("cell_hole",0,twopi,3,hole_rs,hole_zs);
  top_cell_hole = new G4LogicalVolume(solid_top_cell_hole, Ar_mat, "top_cell_hole");
  G4VPhysicalVolume* phys_top_cell_hole = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), top_cell_hole,
	  "TopCellContainer", top_cell, false, 0, checkOverlaps);
  top_GEM = new PseudoMesh(top_ps_plate, G4ThreeVector(0, 0, 0.5*(18 + 4 + plate_W)*mm), 141 * mm, 141 * mm, plate_W,
	  top_cell_container, G4ThreeVector(0, 0, 0.5*(18 + 4 + plate_W + 12)*mm), CELL_SIZE / 2, hexagonal_mapping);
	//================================
#ifdef TOP_MESH_TEST
  G4Box* solid_top_mesh_test_detector = new G4Box("Top_mesh_test_detector", 0.5*(141 * mm), 0.5*(141 * mm), 0.5 *(plate_real_W*mm));
  top_mesh_test_detector = new G4LogicalVolume(solid_top_mesh_test_detector,Ar_mat,"top_mesh_test_detector");
  G4VPhysicalVolume* phys_top_mesh_test_detector = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(18 + 4 + plate_W + 12)*mm), top_cell_hole,
	  "TopCellContainer", box, false, 0, checkOverlaps);
#endif
  G4Box* solid_x_WLS = new G4Box("x_WLS", 0.5*(WLS_FILM_WIDTH), 0.5*(141 * mm + WLS_FILM_WIDTH), 0.5 *(box_sizeZ));
  G4Box* solid_y_WLS = new G4Box("y_WLS", 0.5*(141 * mm + WLS_FILM_WIDTH), 0.5*(WLS_FILM_WIDTH), 0.5 *(box_sizeZ));
  G4Box* solid_x_Acrylic = new G4Box("x_Acrylic", 0.5*(PMMA_WIDTH), 0.5*(141 * mm + PMMA_WIDTH + 2*WLS_FILM_WIDTH), 0.5 *(box_sizeZ));
  G4Box* solid_y_Acrylic = new G4Box("y_Acrylic", 0.5*(141 * mm + PMMA_WIDTH + 2*WLS_FILM_WIDTH), 0.5*(PMMA_WIDTH), 0.5 *(box_sizeZ));
  G4Box* solid_y_LArGap = new G4Box("y_LarGap", 0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH) + 2 * mm), 0.5*(2 * mm), 0.5 *(4 + 0.5 + 20)*mm);
  G4Box* solid_x_LArGap = new G4Box("x_LarGap", 0.5*(2 * mm), 0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH) + 2 * mm), 0.5 *(4 + 0.5 + 20)*mm);

  Xp_wls = new G4LogicalVolume(solid_x_WLS, WLS_mat, "X_PositiveWLS");
  Xn_wls = new G4LogicalVolume(solid_x_WLS, WLS_mat, "X_NegativeWLS");
  Yp_wls = new G4LogicalVolume(solid_y_WLS, WLS_mat, "Y_PositiveWLS");
  Yn_wls = new G4LogicalVolume(solid_y_WLS, WLS_mat, "Y_NegativeWLS");
  Xp_acrylic = new G4LogicalVolume(solid_x_Acrylic, Acrylic_mat, "X_PositiveAcrylic");
  Xn_acrylic = new G4LogicalVolume(solid_x_Acrylic, Acrylic_mat, "X_NegativeAcrylic");
  Yp_acrylic = new G4LogicalVolume(solid_y_Acrylic, Acrylic_mat, "Y_PositiveAcrylic");
  Yn_acrylic = new G4LogicalVolume(solid_y_Acrylic, Acrylic_mat, "Y_NegativeAcrylic");
  Xp_LAr_gap = new G4LogicalVolume(solid_x_LArGap, LAr_mat, "X_PositiveLArGap");
  Xn_LAr_gap = new G4LogicalVolume(solid_x_LArGap, LAr_mat, "X_NegativeLArGap");
  Yp_LAr_gap = new G4LogicalVolume(solid_y_LArGap, LAr_mat, "Y_PositiveLArGap");
  Yn_LAr_gap = new G4LogicalVolume(solid_y_LArGap, LAr_mat, "Y_NegativeLArGap");

  G4VPhysicalVolume* phys_Xp_WLS = new G4PVPlacement(0, G4ThreeVector(+0.5*(141*mm+WLS_FILM_WIDTH), -0.5*(WLS_FILM_WIDTH), 0),
	  Xp_wls, "X_PositiveWLS", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Xn_WLS = new G4PVPlacement(0, G4ThreeVector(-0.5*(141 * mm + WLS_FILM_WIDTH), +0.5*(WLS_FILM_WIDTH), 0),
	  Xn_wls, "X_NegativeWLS", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yp_WLS = new G4PVPlacement(0, G4ThreeVector(+0.5*(WLS_FILM_WIDTH), +0.5*(141 * mm + WLS_FILM_WIDTH), 0),
	  Yp_wls, "Y_PositiveWLS", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yn_WLS = new G4PVPlacement(0, G4ThreeVector(-0.5*(WLS_FILM_WIDTH), -0.5*(141 * mm + WLS_FILM_WIDTH), 0),
	  Yn_wls, "Y_NegativeWLS", box, false, 0, checkOverlaps);

  G4VPhysicalVolume* phys_Xp_Acrylic = new G4PVPlacement(0, G4ThreeVector(0.5*(141 * mm +2*WLS_FILM_WIDTH+PMMA_WIDTH), -0.5*(PMMA_WIDTH), 0),
	  Xp_acrylic, "X_PositiveAcrylic", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Xn_Acrylic = new G4PVPlacement(0, G4ThreeVector(-0.5*(141 * mm +2*WLS_FILM_WIDTH + PMMA_WIDTH), +0.5*(PMMA_WIDTH), 0),
	  Xn_acrylic, "X_NegativeAcrylic", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yp_Acrylic = new G4PVPlacement(0, G4ThreeVector(+0.5*(PMMA_WIDTH), +0.5*(141 * mm + 2 * WLS_FILM_WIDTH + PMMA_WIDTH), 0),
	  Yp_acrylic, "Y_PositiveAcrylic", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yn_Acrylic = new G4PVPlacement(0, G4ThreeVector(-0.5*(PMMA_WIDTH), -0.5*(141 * mm + 2 * WLS_FILM_WIDTH + PMMA_WIDTH), 0),
	  Yn_acrylic, "Y_NegativeAcrylic", box, false, 0, checkOverlaps);

  G4VPhysicalVolume* phys_Xp_LArGap = new G4PVPlacement(0, G4ThreeVector(+0.5*(141 * mm + 2*(WLS_FILM_WIDTH+PMMA_WIDTH)+2*mm), -0.5*(2*mm), -0.5*(18+0.5+20)*mm),
	  Xp_LAr_gap, "X_PositiveLArGap", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Xn_LArGap = new G4PVPlacement(0, G4ThreeVector(-0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH) + 2 * mm), +0.5*(2 * mm), -0.5*(18 + 0.5 + 20)*mm),
	  Xn_LAr_gap, "X_NegativeLArGap", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yp_LArGap = new G4PVPlacement(0, G4ThreeVector(+0.5*(2 * mm), +0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH) + 2 * mm), -0.5*(18 + 0.5 + 20)*mm),
	  Yp_LAr_gap, "Y_PositiveLArGap", box, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yn_LArGap = new G4PVPlacement(0, G4ThreeVector(-0.5*(2 * mm), -0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH) + 2 * mm), -0.5*(18 + 0.5 + 20)*mm),
	  Yn_LAr_gap, "Y_NegativeLArGap", box, false, 0, checkOverlaps);


  G4Tubs* solid_PMT = new G4Tubs("window", 0, (PMT_DIAMETER/*-DET_OFFSET*/) / 2, (4 * mm/*-DET_OFFSET*/)/2, 0 * deg, 360 * deg);
  Xp_PMT = new G4LogicalVolume(solid_PMT, FusedSilica_mat, "X_PositivePMT_Window");
  Xn_PMT = new G4LogicalVolume(solid_PMT, FusedSilica_mat, "X_NegativePMT_Window");
  Yp_PMT = new G4LogicalVolume(solid_PMT, FusedSilica_mat, "Y_PositivePMT_Window");
  Yn_PMT = new G4LogicalVolume(solid_PMT, FusedSilica_mat, "Y_NegativePMT_Window");
  G4RotationMatrix *x_rot_m = new G4RotationMatrix();
  G4RotationMatrix *y_rot_m = new G4RotationMatrix();
  x_rot_m->rotateY(90 * deg);
  y_rot_m->rotateX(90 * deg);

  G4VPhysicalVolume* phys_Xp_PMT = new G4PVPlacement(x_rot_m, G4ThreeVector(+0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH + 2 * mm) + 4 * mm), 0, 0),
	  Xp_PMT, "X_PositivePMT", envelope, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Xn_PMT = new G4PVPlacement(x_rot_m, G4ThreeVector(-0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH + 2 * mm) + 4 * mm), 0, 0),
	  Xn_PMT, "X_NegativePMT", envelope, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yp_PMT = new G4PVPlacement(y_rot_m, G4ThreeVector(0, +0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH + 2 * mm) + 4 * mm), 0),
	  Yp_PMT, "Y_PositivePMT", envelope, false, 0, checkOverlaps);
  G4VPhysicalVolume* phys_Yn_PMT = new G4PVPlacement(y_rot_m, G4ThreeVector(0, -0.5*(141 * mm + 2 * (WLS_FILM_WIDTH + PMMA_WIDTH + 2 * mm) + 4 * mm), 0),
	  Yn_PMT, "Y_NegativePMT", envelope, false, 0, checkOverlaps);
            
  G4OpticalSurface* top_plate_surface = new G4OpticalSurface("top_plate_surface");
  G4LogicalBorderSurface* logical_top_plate_surface = new G4LogicalBorderSurface("top_plate_surface", phys_top_pseudo_GEM, phys_top_CU, top_plate_surface);
  _SetCopperSurface(top_plate_surface);
  G4OpticalSurface* bot_plate_surface0 = new G4OpticalSurface("bot_plate_surface0");
  G4LogicalBorderSurface* logical_bot_plate_surface0 = new G4LogicalBorderSurface("bot_plate_surface0", phys_bot_pseudo_GEM, phys_bot_CU, bot_plate_surface0);
  _SetCopperSurface(bot_plate_surface0);
  G4OpticalSurface* top_plate_surface_r = new G4OpticalSurface("top_plate_real_surface");
  G4LogicalBorderSurface* logical_top_real_plate_surface = new G4LogicalBorderSurface("top_plate_real_surface", phys_top_cell,
	  phys_top_cell_container, top_plate_surface_r);
  _SetCopperSurface(top_plate_surface_r);

  _SetVisibilityParameters();

  detector_volumes.push_back(Xp_PMT);
  detector_volumes.push_back(Xn_PMT);
  detector_volumes.push_back(Yp_PMT);
  detector_volumes.push_back(Yn_PMT);
#ifdef TOP_MESH_TEST
  detector_volumes.push_back(top_mesh_test_detector);
#endif
  
  absorbtion_volumes.push_back(world);
  absorbtion_volumes.push_back(top_cu_plate);
  absorbtion_volumes.push_back(bot_cu_plate);
  absorbtion_volumes.push_back(top_cell);

  //always return the physical World
#undef plate_W
#undef plate_real_W
  return physWorld;
}

void B1DetectorConstruction :: _SetVisibilityParameters(void)
{
	//=====PSUEDO GEMS
	G4VisAttributes* temp = new G4VisAttributes(G4Color(0.5, 0.15, 0.15));
	temp->SetVisibility(false);//TODO: true
	top_cu_plate->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.5, 0.15, 0.15));
	temp->SetVisibility(true);
	bot_cu_plate->SetVisAttributes(temp);
	//=====PSEUDO GEMS
	//=====CELLS
	temp = new G4VisAttributes(G4Color(0.5, 0.15, 0.15));
	temp->SetVisibility(true);
	top_cell->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(_light_blue_RGB,0.1));
	temp->SetVisibility(false);
	top_cell_hole->SetVisAttributes(temp);
	//=====CELLS
	//=====LAr
	temp = new G4VisAttributes(G4Color(_light_blue_RGB,0.3));
	temp->SetVisibility(true);
	LAr_layer->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(_light_blue_RGB, 0.3));
	temp->SetVisibility(true);
	Xp_LAr_gap->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(_light_blue_RGB, 0.3));
	temp->SetVisibility(true);
	Xn_LAr_gap->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(_light_blue_RGB, 0.3));
	temp->SetVisibility(true);
	Yp_LAr_gap->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(_light_blue_RGB, 0.3));
	temp->SetVisibility(true);
	Yn_LAr_gap->SetVisAttributes(temp);
	//======LAr
	//======WSL
	temp = new G4VisAttributes(G4Color(0.15, 0.6, 0.15));
	temp->SetVisibility(true);
	Xp_wls->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.15, 0.6, 0.15));
	temp->SetVisibility(true);
	Xn_wls->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.15, 0.6, 0.15));
	temp->SetVisibility(true);
	Yp_wls->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.15, 0.6, 0.15));
	temp->SetVisibility(true);
	Yn_wls->SetVisAttributes(temp);
	//======WSL
	//======ACRYLIC
	temp = new G4VisAttributes(G4Color(0.7, 0.2, 0.7));
	temp->SetVisibility(true);
	Xp_acrylic->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.7, 0.2, 0.7));
	temp->SetVisibility(true);
	Xn_acrylic->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.7, 0.2, 0.7));
	temp->SetVisibility(true);
	Yp_acrylic->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.7, 0.2, 0.7));
	temp->SetVisibility(true);
	Yn_acrylic->SetVisAttributes(temp);
	//======ACRYLIC
	//======PMTS
	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(true);
	temp->SetForceWireframe(true);
	Xp_PMT->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(true);
	temp->SetForceWireframe(true);
	Xn_PMT->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(true);
	temp->SetForceWireframe(true);
	Yp_PMT->SetVisAttributes(temp);
	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(true);
	temp->SetForceWireframe(true);
	Yn_PMT->SetVisAttributes(temp);

#ifdef TOP_MESH_TEST
	temp = new G4VisAttributes(G4Color(0.8, 0.6, 0.7));
	temp->SetVisibility(true);
	temp->SetForceWireframe(true);
	top_mesh_test_detector->SetVisAttributes(temp);
#endif
	//======PMTS
	//======AUXILARY VOLUMES
	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(false);
	envelope->SetVisAttributes(temp);

	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(false);
	box_interior->SetVisAttributes(temp);

	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(false);
	box->SetVisAttributes(temp);

	temp = new G4VisAttributes(G4Color(0.7, 0.7, 0.7));
	temp->SetVisibility(false);
	world->SetVisAttributes(temp);
	//=======AUXILARY VOLUMES
}
