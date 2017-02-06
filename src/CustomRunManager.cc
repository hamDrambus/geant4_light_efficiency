#include "CustomRunManager.hh"
#include "G4UImanager.hh"

void CustomRunManager::ProcessOneEvent(G4int i_event)
{
	OnEventStartProc();
	G4RunManager::ProcessOneEvent(i_event);
}

void CustomRunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  InitializeEventLoop(n_event,macroFile,n_select);
  this->init_event();
// Event loop
  G4int interval = n_event / 100.0;
  G4int counter = 0;
  for(G4int i_event=0; i_event<n_event; i_event++ )
  {
	G4int were_secondaries = 0;
    ProcessOneEvent(i_event);
	while (is_secondary_MC()) //secondary events
	{
		if (!were_secondaries)
		{
			were_secondaries = 1;
			TerminateOneEvent();
			if (runAborted) goto g_out;
		}
		ProcessOneEvent(extra_run_id + i_event);
		TerminateOneEvent();
		if (runAborted) goto g_out;
	}
    if (!were_secondaries)
		TerminateOneEvent();
	if (interval*counter <= i_event)
	{
		counter++;
		G4cout << "simulated " << i_event << "/" << n_event << " ("<< G4int(100.0*i_event/n_event)<< "%)" << G4endl;
	}
    if(runAborted) break;
  }
  g_out:
  TerminateEventLoop();
}

G4int CustomRunManager::spawn_new_MC_node(const G4Step* step, G4double prob,
	G4ThreeVector momentum, G4ThreeVector polarization, G4int num_of_sims)
{
#ifdef DEBUG_MC_NODES
	G4cout << "RM: MC_spawn BOUNDARY "<< prob<< G4endl;
#endif
	if (last_event == NULL) return -1;
	G4double pppp = prob*get_new_spawn_prob();
#ifdef DEBUG_MC_NODES
	G4cout << "RM: spawn_probability= " << pppp << G4endl;
#endif
	if (pppp < MIN_ALLOWED_PROBABILITY)
	{
#ifdef DEBUG_MC_NODES
		G4cout << "RM: MC_spawn quenched" << G4endl;
#endif
		return 0;
	}
	G4int depth = get_sim_depth();
#if defined(TOP_MESH_TEST)||defined(TEST_MESH_SIDEWAYS)
	if (depth >= 3)
#endif
#if !(defined(TOP_MESH_TEST)||defined(TEST_MESH_SIDEWAYS))
	if (depth >= 8)
#endif

	{
		num_of_sims = 0;
		goto endf;
	}
	//if (depth >= 2)
	//	num_of_sims = (num_of_sims<8)? 1: num_of_sims / 8.0;
endf:
	if (num_of_sims == 0)
	{
#ifdef DEBUG_MC_NODES
		G4cout << "RM: MC_spawn quenched none to simulate" << G4endl;
#endif
		return 0;
	}
	has_finished_secondaries = 0;
	MC_node* node = new MC_node(last_event);
	node->set_MC_node(prob, curr_mapping_state, step->GetPostStepPoint()->GetPosition(), momentum, polarization, step->GetPreStepPoint()->GetTotalEnergy(), num_of_sims);
	return 0;
}

G4int CustomRunManager::spawn_new_MC_node(const G4Step* step, G4double prob, G4Material *WLS_pars, G4int num_of_sims)
{
#ifdef DEBUG_MC_NODES
	G4cout << "RM: MC_spawn WLS "<<prob << G4endl;
#endif
	if (last_event == NULL) return -1;
	G4double pppp = prob*get_new_spawn_prob();
#ifdef DEBUG_MC_NODES
	G4cout << "RM: spawn_probability= " << pppp << G4endl;
#endif
	if (pppp < MIN_ALLOWED_PROBABILITY)
	{
#ifdef DEBUG_MC_NODES
		G4cout << "RM: MC_spawn quenched" << G4endl;
#endif
		return 0;
	}
	G4int depth = get_sim_depth(MC_NODE_CONTINIOUS);
	//if (depth >= 3)
	//{
	//	num_of_sims = 0;
	//	goto endf;
	//}
	//if (depth >= 2)
	//{
	//	num_of_sims = 0;
	//	goto endf;
	//}
	if (depth >= 1)
		num_of_sims = 0;// num_of_sims / 10.0;
endf:
	if (num_of_sims == 0)
	{
#ifdef DEBUG_MC_NODES
		G4cout << "RM: MC_spawn quenched none to simulate" << G4endl;
#endif
		return 0;
	}
	MC_node* node = new MC_node(last_event);
	G4MaterialPropertiesTable* aMaterialPropertiesTable;
	G4MaterialPropertyVector* abs_len;
	G4MaterialPropertyVector* energy_spectrum;

	aMaterialPropertiesTable = WLS_pars->GetMaterialPropertiesTable();
	if (aMaterialPropertiesTable) {
		abs_len = aMaterialPropertiesTable->GetProperty("ABSORBTION_LENGTH");
		energy_spectrum = aMaterialPropertiesTable->GetProperty("WLS_ENERGY_SPECTRUM"); //at the moment contains single out par, but shall spectrum
	}
	else return -1; //error
	G4double ph_E = step->GetPreStepPoint()->GetTotalEnergy();
	node->set_MC_node(prob, curr_mapping_state, step->GetPreStepPoint()->GetPosition(), step->GetPostStepPoint()->GetPosition(),
		abs_len->Value(ph_E),energy_spectrum,num_of_sims);
	has_finished_secondaries = 0;
	return 0;
}

//below manage memory allocation for photon_events* chains (list that is)
//so that every photon process is written at the right place 
G4int CustomRunManager::init_event() //called at the very start af event (=simulation)
{
#ifdef DEBUG_MC_NODES
	G4cout << "RM: init_event()" << G4endl;
#endif
	if (!primary_Monte_Carlo)
		primary_Monte_Carlo = new net_sim_data(NULL);
	while (primary_Monte_Carlo->events.begin() != primary_Monte_Carlo->events.end()) //there could be left over data from earlier runs
	{
		primary_Monte_Carlo->events.pop_back();
	}
	extra_run_id = 0;
	primary_Monte_Carlo->num_of_succ_sims = -1;
	primary_Monte_Carlo->new_sequence();
	last_event=primary_Monte_Carlo->new_event();
	current_working_node = NULL;
	if (last_event == NULL) //should not happen
		return -1;
	return 0;
}
G4int CustomRunManager::next_event(const G4Step* step) //called in UserSteppingAction ()
{
	is_first_process = 1;
#ifdef DEBUG_MC_NODES
	G4cout << "RM: next_event()" << G4endl;
#endif
	MC_node** pp_MC;
	pp_MC = new MC_node*[1];
	if (current_working_node)
	{
		current_working_node->simulation_data.set_event(step, curr_ph_ev_prob, curr_ph_ev_type);
		last_event = current_working_node->new_event(pp_MC); //calls parental new_events() too
		current_working_node = *pp_MC;
		if (NULL == last_event) //not over, secondaries should be searched in primary_MC depr:primary simulation is over
		{
			current_working_node = primary_Monte_Carlo->find_next_to_simulate();
			if (NULL == current_working_node) //entire simulation is over
			{
				has_finished_secondaries = 1;
				curr_ph_ev_prob = 1;
				curr_ph_ev_type = RM_PHOTON_UNDEFINED;
				OnNewSimulationProc();
				primary_Monte_Carlo->new_sequence();
				last_event = primary_Monte_Carlo->new_event();
				delete pp_MC;
				return 0;
			}
			//
			last_event = current_working_node->new_event(pp_MC);
			current_working_node = *pp_MC;
			curr_ph_ev_prob = 1;
			curr_ph_ev_type = RM_PHOTON_UNDEFINED;
			//primary_Monte_Carlo->new_sequence(); //go to next primary sequence then
			//last_event = primary_Monte_Carlo->new_event();
			delete pp_MC;
			return 1;
		}
	}
	else
	{
		primary_Monte_Carlo->set_event(step, curr_ph_ev_prob, curr_ph_ev_type);
		last_event=primary_Monte_Carlo->new_event();
		if (NULL == last_event) //primary fork is over
		{
			current_working_node = primary_Monte_Carlo->find_next_to_simulate();
			if (NULL == current_working_node) //entire simulation is over
			{
				has_finished_secondaries = 1;
				curr_ph_ev_prob = 1;
				curr_ph_ev_type = RM_PHOTON_UNDEFINED;
				OnNewSimulationProc(); //TODO ?: transfer few lines above inside?
				primary_Monte_Carlo->new_sequence();
				last_event = primary_Monte_Carlo->new_event();
				delete pp_MC;
				return 0;
			}
			last_event = current_working_node->new_event(pp_MC);
			current_working_node = *pp_MC;
			curr_ph_ev_prob = 1;
			curr_ph_ev_type = RM_PHOTON_UNDEFINED;
			delete pp_MC;
			return 1;
		}
	}
	curr_ph_ev_prob=1;
	curr_ph_ev_type=RM_PHOTON_UNDEFINED;
	delete pp_MC;
	return 0;
}

G4int CustomRunManager::process_end(const G4Step* step)
{
#ifdef DEBUG_MC_NODES
		G4cout << "RM: process_end()" << G4endl;
#endif
	if (is_first_process)
	{
		is_first_process = 0;
		prev_ph_ev_prob = curr_ph_ev_prob;
		prev_ph_ev_type = curr_ph_ev_type;
		return 0;
	}
	G4double temp_prob=curr_ph_ev_prob;
	G4double temp_type=curr_ph_ev_type;
	curr_ph_ev_prob = prev_ph_ev_prob; //curr ph_ev pars set ot the previous ones, and next_event is called using them
	curr_ph_ev_type = prev_ph_ev_type; //so for the last process in sequence next_event is not called
	G4int res=next_event(step);
	curr_ph_ev_prob = temp_prob;
	curr_ph_ev_type = temp_type;
	prev_ph_ev_prob = curr_ph_ev_prob;
	prev_ph_ev_type = curr_ph_ev_type;
	return res;
}

//depr: SetHit is used for this purpose G4int close_event(); //called at EndOfEvent (end of simulation)

G4int CustomRunManager::is_secondary_MC(void) //true if there are secondary MC (MC_nodes) to be simulated
{
	return !has_finished_secondaries;
}

void CustomRunManager::SetPhEvType(G4int evType)
{
#ifdef DEBUG_MC_NODES
	G4cout << "RM: SetPhEvType "<<evType << G4endl;
#endif
	curr_ph_ev_type = evType;
}

void CustomRunManager::SetPhEvProb(G4double evProb)
{
#ifdef DEBUG_MC_NODES
	G4cout << "RM: SetPhEvProb "<<evProb << G4endl;
#endif
	curr_ph_ev_prob = evProb;
}

void CustomRunManager::SetHit(G4bool is_hit)
{
#ifdef DEBUG_MC_NODES
	G4cout << "RM: SetHit "<<is_hit << G4endl;
#endif
	if (current_working_node)
	{
		extra_run_id++; //probably not required
		current_working_node->SetHit(is_hit);
	}
	else
		primary_Monte_Carlo->SetHit(is_hit);
}

//!WARNING this is recursive function may fail at strongly nested event trees
G4double CustomRunManager::get_total_detetion_eff(void) //not checked whether simulation is completed, should be called after simulation
{
#ifdef DEBUG_MC_NODES
	G4cout << "RM: get_total_detecton_eff" << G4endl;
#endif
	return primary_Monte_Carlo->hit_prob();
}

void CustomRunManager::get_total_detetion_eff(G4double* no_reemiss, G4double *reemissed, G4double* total)
{
	*total = 0;
	*no_reemiss = 0;
	*reemissed = 0;
	primary_Monte_Carlo->hit_prob(no_reemiss, reemissed, total);
}

//0 - kill process, 1 - select reflection, 2 - select defraction, 3 - select both
G4int CustomRunManager::select_photon_BP(const G4Step* step, G4ThreeVector defl_momentum, G4ThreeVector refl_momentum)
{
	//TODO: major rework required (main point - cut off total reflection? depending on absortion length for a given energy?)
#ifdef DEBUG_MC_NODES
	G4cout << "RM: select_photon_BP()" << G4endl;
#endif
	const B1DetectorConstruction* detC
		= static_cast<const B1DetectorConstruction*> (this->GetUserDetectorConstruction());
	G4LogicalVolume* pre_volume=step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	G4LogicalVolume* post_volume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	if ((pre_volume == detC->bot_ps_plate) || (pre_volume == detC->top_ps_plate))
	{
		if ((post_volume != detC->bot_cu_plate) && (post_volume != detC->top_cu_plate))
			return RM_CHOOSE_DEFL;
		else
			return RM_CHOOSE_REFL;
	}
	if (((pre_volume == detC->top_cell_hole) || (pre_volume == detC->top_cell_container) || (pre_volume == detC->top_cell_hole_dielectric)))
	{
#if defined(TOP_MESH_TEST)||defined(TEST_MESH_SIDEWAYS)
		return RM_CHOOSE_BOTH;
#endif
		if (post_volume==detC->top_cell)

			return RM_CHOOSE_REFL;
	}
	if (post_volume == detC->bot_cell)
		return RM_CHOOSE_REFL;
	if (post_volume == detC->top_cell_hole_dielectric)
		return RM_CHOOSE_BOTH;
#ifdef TEMP_CODE_
	return RM_CHOOSE_DEFL;
#endif
	return RM_CHOOSE_BOTH; //TODO: W! too much treeing
}

G4ThreeVector CustomRunManager::FetchPosition()
{return initial_position;}

G4ThreeVector CustomRunManager::GenPosition()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (current_working_node != NULL)
		return initial_position = current_working_node->GenPostion();
#ifdef TOP_MESH_TEST
	if (t_counter < x_num*y_num)
	{
		G4double x = x_start+t_step*((int)(t_counter / y_num));
		G4double y = y_start + t_step*((int)(t_counter % y_num));
		x += G4RandGauss::shoot(G4Random::getTheEngine(), 0.0, t_uncert);
		y += G4RandGauss::shoot(G4Random::getTheEngine(),0.0,t_uncert);
		t_counter++;	
		return initial_position = G4ThreeVector(x, y, 0);
	}
	else
		return initial_position = G4ThreeVector(0, 0, 0);
#endif
#ifdef TEST_MESH_SIDEWAYS
	return initial_position = G4ThreeVector(-9, 0, 10.9); //CELL SIZE = 4.5
#endif
//#ifdef TEMP_CODE_
//	return initial_position = G4ThreeVector(3.6+0.42, 0, 0);
//#endif
	return initial_position = G4ThreeVector(0, 0, 0); //primary event parameters are such for a while
}

G4ThreeVector CustomRunManager::FetchMomentum()
{return initial_momentum_direction;}

G4ThreeVector CustomRunManager::GenMomentum()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (current_working_node != NULL)
		return initial_momentum_direction = current_working_node->GenMomentum();
#ifdef TOP_MESH_TEST
	return initial_momentum_direction = G4ThreeVector(0, 0, 1);
#endif
#ifdef TEST_MESH_SIDEWAYS
	//1*cell_size/(2*h+plate_width)<z<3*cell_size/(2*h+plate_width), h==z_ps_gem_boundary - z_initial 
	G4double dx = 3 * 4.5 / (0.2 + 0.6);
	return initial_momentum_direction = G4ThreeVector(dx*cos(CLHEP::pi / 3), dx*sin(CLHEP::pi / 3), 1).unit(); //CELL SIZE = 4.5
#endif
//#ifdef TEMP_CODE_
//	return initial_momentum_direction = G4ThreeVector(0, 0, 1);
//#endif
	G4double phi = CLHEP::twopi*G4UniformRand();
	G4double cos_theta = 2*(G4UniformRand())-1;
	return initial_momentum_direction = G4ThreeVector(sin(phi)*sqrt(1 - cos_theta*cos_theta), cos(phi)*sqrt(1 - cos_theta*cos_theta), cos_theta);
}

G4ThreeVector CustomRunManager::FetchPolarization()
{return initial_polarization;}

G4ThreeVector CustomRunManager::GenPolarization()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (current_working_node != NULL)
		return initial_polarization = current_working_node->GenPolarization();
	G4double phi = CLHEP::twopi*G4UniformRand();
	G4double cos_theta = 2 * (G4UniformRand()) - 1;
	return initial_polarization = G4ThreeVector(sin(phi)*sqrt(1 - cos_theta*cos_theta), cos(phi)*sqrt(1 - cos_theta*cos_theta), cos_theta);
}

G4double CustomRunManager::FetchEnergy()
{return initial_energy;}

G4double CustomRunManager::GenEnergy()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (current_working_node != NULL)
		return initial_energy = current_working_node->GenEnergy();
// 9.65*eV ==128nm
// 3.9236*eV == 316nm
// 3.6791*eV == 337nm
// 3.4632*eV == 358nm
// 3.2585*eV == 380.5nm
// 3.0538*eV == 406nm
	return initial_energy = 3.9236*eV; //primary event parameters are such for a while
}

PseudoMeshData*	CustomRunManager::GenMappingState() //Position is assumed to be already updated before the call
{
	PseudoMeshData* temp=NULL;
	if (current_working_node != NULL)
		temp = current_working_node->GenMappingData();
	if (temp == NULL) //in case current_working_node!=NULL this is erroneous situation but not fatal
	{
		B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(GetUserDetectorConstruction());
		detectorConstruction->GetPseudoMeshByPoint(curr_mapping_state, FetchPosition(),FetchMomentum());
		if (curr_mapping_state->curr_mesh != NULL)
			curr_mapping_state->curr_mesh->GetDefaultMappingData(curr_mapping_state);
		else
			curr_mapping_state->SetDefauldInd();
	}
	else
		*curr_mapping_state = *temp;
	return curr_mapping_state;
}

G4double CustomRunManager::get_new_spawn_prob()
{
	//last_event!=NULL
	G4int is_over = 0;
	G4double prob=-1;
	for (auto j = last_event->container->events.rbegin(); j != last_event->container->events.rend(); ++j)
	{
		if (&(*j) == last_event)
		{
			if ((++j) != last_event->container->events.rend())
			{
				prob = (j->net_prob);
				//j++;//forgot at first
				break;
			}
			else
			{
				prob = 1;//last_event is the first in the list
				break;
			}
		}
	}
	if (-1==prob ) return -1;//error
	if (last_event->container->container->parent == NULL)
	{
#ifdef DEBUG_MC_NODES
		//G4cout << "RM: return get_new_spawn_prob() " << prob << G4endl;
#endif
		return prob; //last_event belongs to primary simulation
	}
	MC_node* currP = last_event->container->container->parent;
	prob = prob*currP->node_probability;
	while (!is_over)
	{
		photon_event* parentt = currP->ev_parent;
		if (parentt == 0)
			return -1;//error
		for (auto j = parentt->container->events.rbegin(); j != parentt->container->events.rend(); ++j)
		{
			if (&(*j) == parentt)
			{
				if ((++j) != parentt->container->events.rend())//reverse iterators must be used here depr:WARNING! is it correct?
				{
					prob = prob*(j->net_prob);
					break; //'for' break - found last_photon_event probability
				}//if ==rend () then prob*=1;
				break;
			}
		}
		if (parentt->container->container->parent == NULL)
			break;
		currP = parentt->container->container->parent;
		prob *= currP->node_probability;
	}
#ifdef DEBUG_MC_NODES
	//G4cout << "RM: return get_new_spawn_prob() "<<prob << G4endl;
#endif
	return prob;
}

G4int CustomRunManager::get_sim_depth()
{
	G4int out = 0;
	if (current_working_node == NULL)
		return 0;
	MC_node* p=current_working_node;
	photon_event* par;
	while (NULL != p)
	{
		out++;
		par = p->ev_parent;
		p = par->container->container->parent;
	}
	return out;
}

G4int CustomRunManager::get_sim_depth(G4int node_type)
{
	G4int out = 0;
	if (current_working_node == NULL)
		return 0;
	MC_node* p = current_working_node;
	photon_event* par;
	while (NULL != p)
	{
		par = p->ev_parent;
		p = par->container->container->parent;
		if (p)
			if (p->chosen_type == node_type)
				out++;
	}
	return out;
}

G4double CustomRunManager::get_curr_event_probab()
{
	G4double prob = 1;
	prob *= last_event->net_prob;
	MC_node* p = current_working_node;
	photon_event* par;
	while (NULL != p)
	{
		par = p->ev_parent;
	//TODO: make function in photon_event to obtaint current probability
		for (auto j = par->container->events.rbegin(); j != par->container->events.rend(); ++j)
		{
			if (&(*j) == par)
			{
				if (par->container->events.rend() != (++j))
					prob *= j->net_prob;
				break;
			}
		}
		p = par->container->container->parent;
		if (p)
			prob *= p->node_probability;
	}
	return prob;
}

G4int CustomRunManager::get_detected_spectrum(std::list<G4double> *energies, std::list < G4double> *probabilities)
{
	G4double probab = 1;
	for (auto h = primary_Monte_Carlo->events.begin(); h != primary_Monte_Carlo->events.end(); h++)
	{
		G4double local_prob = probab;
		for (auto j = h->events.begin(); j != h->events.end(); j++)
		{
			if (NULL != j->daughter_node)
				get_detected_spectrum(local_prob, j->daughter_node, energies, probabilities);
			local_prob *= j->continue_prob;
		}
		if (h->events.size() != 0)
		{
			if ((1 == h->events.back().has_hit) && (0 != local_prob))
			{
				energies->push_back(h->events.back().energy);
				probabilities->push_back(local_prob);
			}
		}
	}
	return 0;
}

void CustomRunManager::get_detected_spectrum(G4double reach_node_prob, MC_node* node, std::list<G4double> *energies, std::list < G4double> *probabilities)
{
	G4double probab = reach_node_prob*node->node_probability;
	G4int fgfg = node->simulation_data.events.size();
	if (0 == fgfg)
		return;
	probab = probab / fgfg;
	for (auto h = node->simulation_data.events.begin(); h != node->simulation_data.events.end(); h++)
	{
		G4double local_prob = probab;
		for (auto j = h->events.begin(); j != h->events.end(); j++)
		{
			if (NULL != j->daughter_node)
				get_detected_spectrum(local_prob, j->daughter_node, energies, probabilities);
			local_prob *= j->continue_prob;
		}
		if (0!=h->events.size())
		{
			if ((1 == h->events.back().has_hit) && (0!=local_prob))
			{
				energies->push_back(h->events.back().energy);
				probabilities->push_back(local_prob);
			}
		}
	}
}

void CustomRunManager::get_detected_spectrum(std::string out_filename)
{
#define NUM_OF_BINS 150
	std::list<G4double> energies;
	std::list<G4double> probs;
	get_detected_spectrum(&energies, &probs);
	if (NUM_OF_BINS >= energies.size())
		return;
	G4double min_e=-1, max_e=-1, delta_e=-1;
	std::list<G4double> en_bins;
	std::list<G4double> sum_probs;
	for (auto i = energies.begin(); i != energies.end(); i++)
	{
		min_e = min_e < 0 ? *i : (*i < min_e ? *i : min_e);
		max_e = max_e < 0 ? *i : (*i > max_e ? *i : max_e);
	}
	delta_e = (max_e - min_e) / NUM_OF_BINS;
	for (int counter = 0; counter < NUM_OF_BINS; counter++)
	{
		en_bins.push_back(min_e + counter*(max_e - min_e) / NUM_OF_BINS);
		sum_probs.push_back(0);
	}
	for (auto i = energies.begin(), j = probs.begin(); (i != energies.end()) && (j != probs.end()); i++, j++)
		for (auto e = en_bins.begin(), sp=sum_probs.begin(); (e!= en_bins.end())&&(sp!=sum_probs.end()); sp++,e++)
			if (((*e-delta_e/2) < *i) && (*i < (*e + delta_e/2)))
			{
				*sp += *j;
				break;
			}
	std::ofstream file(out_filename);
	for (auto e = en_bins.rbegin(), sp = sum_probs.rbegin(); (e != en_bins.rend()) && (sp != sum_probs.rend()); ++sp, ++e)
	{
		file << 1.2398e3*(1 * eV / (*e))/*nm*/ << "\t" << *sp << std::endl;
	}
	file.close();
#undef NUM_OF_BINS
}

#ifdef TOP_MESH_TEST
//!!WARNING - horizontal positioni of test detector is assumed
void CustomRunManager::on_hit_proc(G4ThreeVector point, G4double prob, G4int bot_top) //called when test detector is hit with photon
{
	G4double net_prob = prob*get_curr_event_probab();
	if (net_prob < 0) return;
	if (bot_top == 0)
	{
		bot_hits_xs.push_back(point.x());
		bot_hits_ys.push_back(point.y());
		bot_hits_probs.push_back(net_prob);
	}
	if (bot_top == 1)
	{
		top_hits_xs.push_back(point.x());
		top_hits_ys.push_back(point.y());
		top_hits_probs.push_back(net_prob);
	}
}

void CustomRunManager::export_to_bmp(std::list<G4double>* hits_xs, std::list<G4double> *hits_ys, std::list<G4double> *hits_probs, G4String filename) 
//^first step is to map position of hits with arbitrary coordinates to bitmap. Then write to file 
{
	t_counter = 0; //nullifies history so "scan" can be done again
	G4double * bits = new G4double[x_num*y_num];
	for (int g = 0; g < x_num*y_num; g++)
		bits[g] = 0;
	G4double lx = t_step*x_num;
	G4double ly = t_step*y_num;
	while ((!hits_xs->empty()) && (!hits_ys->empty()) && (!hits_probs->empty()))
	{
		G4double x = hits_xs->back()-x_start;
		G4double y = hits_ys->back()-y_start;
		G4double prob = hits_probs->back();
		hits_xs->pop_back();
		hits_ys->pop_back();
		hits_probs->pop_back();
		if ((x<0) || (x>t_step*x_num))
			continue;
		if ((y<0) || (y>t_step*y_num))
			continue;
		G4int x_ind = (int)(x /t_step);
		G4int y_ind = (int)(y /t_step);
#define NUM_TO_CONSIDER 3
		//^should be odd
		G4double neighbours[NUM_TO_CONSIDER*NUM_TO_CONSIDER];
		//firstly - obtain squared lengths to the centers of bits in bmp, then exp distribution to the neighbours
		//then - renormalization and writing to global double array
		G4double sum=0;
		for (G4int h = 0; h < NUM_TO_CONSIDER*NUM_TO_CONSIDER; h++)
		{
			G4int te_x_ind = (h / NUM_TO_CONSIDER) - NUM_TO_CONSIDER/2;
			G4int te_y_ind = h % NUM_TO_CONSIDER - NUM_TO_CONSIDER / 2;
			neighbours[h] = (x - t_step*(x_ind+te_x_ind))*(x - t_step*(x_ind + te_x_ind)) 
				+ (y - t_step*(y_ind + te_y_ind))*(y - t_step*(y_ind + te_y_ind));
			//^lengths to centers of neighbour cells
			neighbours[h] = exp(-neighbours[h] / (t_step*t_step));
			if (((te_x_ind + x_ind) > 0) && ((te_x_ind + x_ind) < x_num) && ((te_y_ind + y_ind) > 0) && ((te_y_ind + y_ind) < y_num))
				sum += neighbours[h];
		}
		if (0 != sum)
		{
			sum = prob / sum;
			for (G4int h = 0; h < NUM_TO_CONSIDER*NUM_TO_CONSIDER; h++)
			{
				G4int te_x_ind = (h / NUM_TO_CONSIDER) - NUM_TO_CONSIDER / 2;
				G4int te_y_ind = h % NUM_TO_CONSIDER - NUM_TO_CONSIDER / 2;
				neighbours[h] *= sum;
				if (((te_x_ind + x_ind) >=0 ) && ((te_x_ind + x_ind) < x_num) && ((te_y_ind + y_ind) >= 0) && ((te_y_ind + y_ind) < y_num))
					bits[(te_x_ind + x_ind)*y_num + (te_y_ind + y_ind)] += neighbours[h];
			}
		}
#undef NUM_TO_CONSIDER
	}
	G4double max = -1;
	for (int g = 0; g < x_num*y_num; g++)
	{
		if (max < bits[g]) max = bits[g];
	}
	if (max <= 0)
	{
		delete[] bits;
		return;
	}
	unsigned char *_bits = new unsigned char[x_num*y_num];
	for (int g = 0; g < x_num*y_num; g++)
		_bits[g] = (bits[g] < 0 ? 0 : (int)(255.0*bits[g] / max));
	delete[] bits;
	//WORK WITH BMP
	unsigned int headers[13];
	G4int extra_bytes = (4 - ((x_num * 3 )% 4))%4;
	int paddle_size = (x_num * 3 + extra_bytes)*y_num;
	headers[0] = paddle_size + 54;
	headers[1] = 0;
	headers[2] = 54;
	headers[3] = 40;
	headers[4] = x_num;  
	headers[5] = y_num;
	headers[7] = 0;
	headers[8] = paddle_size;
	headers[9] = 0;
	headers[10] = 0;
	headers[11] = 0;
	headers[12] = 0;
	std::ofstream bmp;
	bmp.open(filename,std::ios_base::trunc|std::ios_base::binary);
	bmp << "BM";
	for (G4int h = 0; h <= 5; h++)
	{
		bmp << (char)(headers[h] & 0x000000ff);
		bmp << (char)((headers[h] & 0x0000ff00)>>8);
		bmp << (char)((headers[h] & 0x00ff0000)>>16);
		bmp << (char)((headers[h] & (unsigned int) 0xff000000)>>24);
	}
	bmp << (char)1 << (char)0 << (char)24 << (char)0;
	for (G4int h = 7; h <= 12; h++)
	{
		bmp << (char)(headers[h] & 0x000000ff);
		bmp << (char)((headers[h] & 0x0000ff00) >> 8);
		bmp << (char)((headers[h] & 0x00ff0000) >> 16);
		bmp << (char)((headers[h] & (unsigned int)0xff000000) >> 24);
	}
	for (G4int h = 0; h < y_num;h++)
	{
		for (G4int g = 0; g < x_num; g++)
			bmp << _bits[g*y_num + h] << _bits[g*y_num + h] << _bits[g*y_num + h];
		for (int g = 0; g < extra_bytes; g++)
			bmp << (char)0;
	}
	bmp.close();
	//END BMP OUTPUT
	delete[] _bits;
}S
#endif

void CustomRunManager::OnEventStartProc()
{
	B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(GetUserDetectorConstruction());
	GenPosition();
	GenPolarization();
	GenMomentum();
	GenEnergy();
	GenMappingState();
	detectorConstruction->OnEventStartProc(this);
}

G4ThreeVector CustomRunManager::MappingProc(const G4Track& track, const G4Step& aStep, G4TouchableHandle &fCurrentTouchableHandle)
{
	B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)GetUserDetectorConstruction();
	return detectorConstruction->MappingProc(curr_mapping_state, track, aStep, fCurrentTouchableHandle);
}

void CustomRunManager::OnNewSimulationProc(void)
{

}
