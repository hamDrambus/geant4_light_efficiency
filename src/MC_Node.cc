#include "MC_Node.hh"

G4double one_sim_data::hit_prob() //recursive function
{
#ifdef DEBUG_MC_NODES
	G4cout << "OSM: hit_prob()" << G4endl;
#endif
	G4int no_sim=0;
	if (events.begin() == events.end()) return -1;
	G4double pb = events.back().net_prob*events.back().has_hit;
	if (pb < 0) {pb = 0; no_sim = 1;}
	for (auto i = events.rbegin(); i != events.rend(); ++i)
	{
		if (i->daughter_node)
		{
			if (!(i->daughter_node->is_accounted))
			{
				G4double reach_prob;
				auto j = i;
				if ((++j) != events.rend())
					reach_prob = (j->net_prob);
				else reach_prob=1;
				G4double node_probab = i->daughter_node->hit_prob();
				if (node_probab >= 0) no_sim = 0;
				node_probab = node_probab < 0 ? 0 : node_probab;
				pb += reach_prob*node_probab;
			}
		}
	}
	return no_sim?-1:pb;
}

void one_sim_data::clear_prob_calc()
{
#ifdef DEBUG_MC_NODES
G4cout << "OSM: clear_prob_calc()" << G4endl;
#endif
	if (events.begin() == events.end()) return;
	for (auto i = events.rbegin(); i != events.rend(); ++i)
	{
		if (i->daughter_node)
		{
			i->daughter_node->clear_prob_calc();
		}
	}
}

void one_sim_data::set_event(const G4Step* step, G4double prob, G4int ev_type, G4int is_hit)
{
#ifdef DEBUG_MC_NODES
	G4cout << "OSM: set_event()" << G4endl;
#endif
	if (events.end()!=events.begin())
	{
		events.back().has_hit = (events.back().has_hit == -1) ? is_hit : events.back().has_hit;
		events.back().chosen_event_type = ev_type;
		events.back().post_step_pos = step->GetPostStepPoint()->GetPosition();
		events.back().pre_step_pos = step->GetPreStepPoint()->GetPosition();
		events.back().post_volume = step->GetPostStepPoint()->GetPhysicalVolume();
		events.back().pre_volume = step->GetPreStepPoint()->GetPhysicalVolume();
		events.back().continue_prob = prob;
		events.back().net_prob = prob*events.back().net_prob;
		events.back().energy = step->GetPreStepPoint()->GetTotalEnergy();
	}
}

void one_sim_data::SetHit(G4bool is_hit)
{
#ifdef DEBUG_MC_NODES
	//G4cout << "OSM: SetHit" << G4endl;
#endif
	if (events.end() != events.begin())
		events.back().has_hit = is_hit;
}

photon_event* one_sim_data::new_event()
{
#ifdef DEBUG_MC_NODES
	G4cout << "OSM: new_event()" << G4endl;
#endif
	G4bool was_empty = false;
	if (events.end() == events.begin()) was_empty = true;
	events.push_back(photon_event(this));
	if (!was_empty)
		events.back().net_prob = (++events.rbegin())->net_prob;
#ifdef DEBUG_MC_NODES
	G4cout << "OSM: return new_event()" << G4endl;
#endif
	return &(events.back());
}

G4int one_sim_data::is_hit()
{
#ifdef DEBUG_MC_NODES
	G4cout << "OSM: is_hit()" << G4endl;
#endif
	if (events.end() != events.begin())
		this->has_hit = events.back().has_hit;
	else return -1; //no events - no hit
	return this->has_hit;
}

net_sim_data::net_sim_data(MC_node* paren) :probability(-1), num_of_succ_sims(-1),parent(paren) {};

//recursive function, average here
G4double net_sim_data::hit_prob()
{
#ifdef DEBUG_MC_NODES
	G4cout << "NSM: hit_prob()" << G4endl;
#endif
	if (events.empty()) return -1;
	probability = 0;
	num_of_succ_sims = 0;
	G4double ppp;
	for (auto i = events.begin(); i != events.end(); i++)
	{
		ppp = i->hit_prob();
		if (ppp >= 0)
		{
			num_of_succ_sims++;
			probability += ppp;
		}
	}
	if (num_of_succ_sims == 0) return -1;
	return probability/num_of_succ_sims;
}

void net_sim_data::clear_prob_calc()
{
#ifdef DEBUG_MC_NODES
	G4cout << "NSM: clear_prob_calc()" << G4endl;
#endif
	probability = -1;
	num_of_succ_sims = -1;
	if (events.empty()) return;
	for (auto i = events.begin(); i != events.end(); i++)
	{
		i->clear_prob_calc();
	}
}

void net_sim_data::new_sequence(void)
{
#ifdef DEBUG_MC_NODES
	G4cout << "NSM: new_sequence" << G4endl;
#endif
	if (parent)
	{
		parent->simulated++;
		if (!parent->is_over()) //there is need to create a new sequence
			events.push_back(one_sim_data(this));
	}
	else //primary simulation
	{
		events.push_back(one_sim_data(this));
	}
}

void net_sim_data::SetHit(G4bool is_hit)
{
#ifdef DEBUG_MC_NODES
	//G4cout << "NSM: SetHit" << G4endl;
#endif
	if (events.end() != events.begin())
	{
		events.back().SetHit(is_hit);
	}
}

photon_event* net_sim_data::new_event()
{
#ifdef DEBUG_MC_NODES
	G4cout << "NSM: new_event()" << G4endl;
#endif
	if (events.end() == events.begin())
	{
		return NULL;//events.push_back(one_sim_data(this)); //depr: this must be called in new_sequence
	}
	if (events.back().is_hit()!=-1) return NULL; //0 or 1 - kill or hit
	return events.back().new_event();
}

void net_sim_data::set_event(const G4Step* step, G4double prob, G4int ev_type, G4int is_hit)
{
#ifdef DEBUG_MC_NODES
	G4cout << "NSM: set_event()" << G4endl;
#endif
	if (events.end() != events.begin())
		events.back().set_event(step, prob, ev_type, is_hit);
}

MC_node* net_sim_data::find_next_to_simulate(void)
{
#ifdef DEBUG_MC_NODES
	G4cout << "NSM: find_next_to_simulate()" << G4endl;
#endif
	if (events.begin() == events.end()) return NULL;
	for (auto g = events.back().events.begin(); g != events.back().events.end(); g++)
		if (NULL!=g->daughter_node)
			if (!g->daughter_node->is_over()) return g->daughter_node;
	return NULL;
}

MC_node::MC_node(photon_event* parent) :simulation_data(this), ev_parent(parent)
{
#ifdef DEBUG_MC_NODES
	//G4cout << "MCN: constructor" << G4endl;
#endif
	node_probability = 1;
	to_simulate=-1;
	simulated=-1;
	spacial_disrt_par=-1;
	PhotonEnergy=-1;
	parent_container=parent->container;
	parent_cont_cont = parent_container?parent_container->container:NULL;
	gl_parent = parent_cont_cont?parent_cont_cont->parent:NULL;
	chosen_type=MC_NODE_UNDEFINED;
	is_accounted = 0;
	parent->daughter_node = this;
	EnergySpectrum = NULL;
}

photon_event* MC_node::new_event(MC_node** pp_new_MC_node)//if pp_new_MC_node set to MC_node where new event is created
{
#ifdef DEBUG_MC_NODES
	G4cout << "MCN" << chosen_type << ": new_event()" << G4endl;
#endif
	MC_node* next_MC=NULL;
	if (is_over()) //happens when node is over. depr: WARNING! this should never happen, so this case must be reworked to simply return NULL
		//but then all functions calling this method must be changed to account for a possible returned NULL
	{
		if (NULL != this->gl_parent) //there is parent MC_node
			return this->gl_parent->new_event(pp_new_MC_node);
		else //no parent MC_node hence call new event from container (i.e. primary simulation)
			//this is managed in CustomRunManager->next-event, so return NULL
		{
			*pp_new_MC_node = NULL;
		return NULL;
//			*pp_new_MC_node = NULL;
//#ifdef DEBUG_MC_NODES
//			G4cout << "MCN" << chosen_type << ": return primary new_event()" << G4endl;
//#endif	
//			return (this->ev_parent ? (this->ev_parent->container ? (this->ev_parent->container->container ? 
//				this->ev_parent->container->container->new_event():NULL):NULL):NULL);
		}
	}
	else
	{
		if (simulation_data.events.empty())
		{
			simulation_data.new_sequence(); //just started 'this' node - new sequence creating required
#ifdef DEBUG_MC_NODES
			G4cout << "MCN" << chosen_type << ": first sequece new_event()" << G4endl;
#endif
		}
		photon_event* pp = simulation_data.new_event(); //first major simulation checked, then secondary MC's
		//^returns NULL if there was hit/kill or empty
		if (pp != NULL)
		{
			*pp_new_MC_node = this;
#ifdef DEBUG_MC_NODES
			G4cout << "MCN" << chosen_type << ": return secondary_same_sequence new_event()" << G4endl;
#endif
			return pp;
		}
		else //secondary MC's
		{
			next_MC = simulation_data.find_next_to_simulate();
			if (next_MC == NULL) //all daughter are simulated, go to next sequence in this node
			{
				//moved to new_sequence()// simulated++;
				simulation_data.new_sequence();
				if (is_over())
				{
					if (NULL != this->gl_parent) //there is parent MC_node
					{
#ifdef DEBUG_MC_NODES
						G4cout << "MCN" << chosen_type << ": return over_goto_parent new_event()" << G4endl;
#endif
						return this->gl_parent->new_event(pp_new_MC_node);
					}
					else
					{
#ifdef DEBUG_MC_NODES
						G4cout << "MCN" << chosen_type << ": return over_goto_primary new_event()" << G4endl;
#endif
						*pp_new_MC_node = NULL; //there is no parent MC_node (main simulation)
						return (this->ev_parent ? (this->ev_parent->container ? (this->ev_parent->container->container ?
							this->ev_parent->container->container->new_event() : NULL) : NULL) : NULL);
					}
				}
				else
				{
#ifdef DEBUG_MC_NODES
					G4cout << "MCN" << chosen_type << ": return secondary_new_sequence new_event()" << G4endl;
#endif
					*pp_new_MC_node = this;
					return simulation_data.new_event();
				}
			}
#ifdef DEBUG_MC_NODES
			G4cout << "MCN" << chosen_type << ": return secondary_daughter_MC new_event()" << G4endl;
#endif
			return next_MC->new_event(pp_new_MC_node);
		}
	}
	*pp_new_MC_node = this;
	return simulation_data.new_event(); //never called
}

void MC_node::SetHit(G4bool is_hit)
{
#ifdef DEBUG_MC_NODES
	G4cout << "MCN" << chosen_type << ": SetHit" << G4endl;
#endif
	if ((simulated == -1) || (chosen_type == MC_NODE_UNDEFINED)) return;
	simulation_data.SetHit(is_hit);
}

G4int MC_node::is_over(void)
{
#ifdef DEBUG_MC_NODES
	G4int ff=(simulated >= to_simulate);
	if (to_simulate == 0) ff = 1;
	G4cout << "MCN"<<chosen_type<<": is_over() "<<ff << G4endl;
	return ff;
#endif
	return (simulated >= to_simulate); //simulated ++ only when all daughter simulations are over -in call of new_sequence
	//&&(simulation_data.find_next_to_simulate()==NULL); 
}

G4ThreeVector MC_node::GenPostion()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (chosen_type == MC_NODE_UNDEFINED) return G4ThreeVector();
	if (chosen_type == MC_NODE_DESCRETE) return start_point1;
	G4double R = (start_point1-start_point2).getR();
	if (R == 0) return start_point1;
	else
	{
		G4double rand_R = -spacial_disrt_par*log(1.0 - G4UniformRand()*(1 - exp(-R / spacial_disrt_par)));
		rand_R = rand_R / R;
		//modification because of tolerances
		//(one must not simulate emission of photons near the boundaries (that is start_point1 and start_point2)
		//because boundary processes won't trigger at very small (<G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / 2) distances
		//and particle will be trasferred directly into a new volume)
		G4double Toler = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
		if (Toler >= R)
			rand_R = 0.5;//just in the middle then
		else
		{
			Toler = Toler / 2;
			rand_R = Toler / R + rand_R*(R - 2 * Toler);
		}
		//end modification
		return start_point1 + rand_R*(start_point2 - start_point1);
	}
}
G4ThreeVector MC_node::GenMomentum()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (chosen_type == MC_NODE_UNDEFINED) return G4ThreeVector();
	if (chosen_type == MC_NODE_DESCRETE)
		return PhotonMomentum;
	G4double phi = CLHEP::twopi*G4UniformRand();
	G4double cos_theta = 2 * (G4UniformRand()) - 1;
	return G4ThreeVector(sin(phi)*sqrt(1 - cos_theta*cos_theta), cos(phi)*sqrt(1 - cos_theta*cos_theta), cos_theta);
}
G4ThreeVector MC_node::GenPolarization()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (chosen_type == MC_NODE_DESCRETE)
		return PhotonPolarization;
	G4double phi = CLHEP::twopi*G4UniformRand();
	G4double cos_theta = 2 * (G4UniformRand()) - 1;
	return G4ThreeVector(sin(phi)*sqrt(1 - cos_theta*cos_theta), cos(phi)*sqrt(1 - cos_theta*cos_theta), cos_theta);
}

G4double MC_node::GenEnergy()
{
#ifdef DEBUG_MC_NODES
	//G4cout << "" << G4endl;
#endif
	if (chosen_type == MC_NODE_DESCRETE)
		return PhotonEnergy;
	if (chosen_type == MC_NODE_CONTINIOUS)
		return EnergySpectrum->Value(G4UniformRand()); //P.S. neat
	return PhotonEnergy;
}

//better not be called after the first sumulation is started
void MC_node::set_MC_node(G4double prob, G4ThreeVector _start_point1, G4ThreeVector momentum, G4ThreeVector polarization, G4double energy, G4int num_of_sims)
{
	if (num_of_sims == 0)
	{
#ifdef DEBUG_MC_NODES
		G4cout << "MCN: set_MC_node() DESCRETE none to simulate" << G4endl;
#endif
		return;
	}
#ifdef DEBUG_MC_NODES
	G4cout << "MCN: set_MC_node() DESCRETE" << G4endl;
#endif
	chosen_type = MC_NODE_DESCRETE;
	node_probability = prob;
	to_simulate = num_of_sims;
	simulated = simulated == -1 ? -1 : simulated; //means override when setting again, but without clearing sumulated data
	start_point1 = _start_point1;
	PhotonMomentum = momentum;
	PhotonPolarization = polarization;
	PhotonEnergy = energy;
	EnergySpectrum = NULL;
}

//better not be called after the first sumulation is started
void MC_node::set_MC_node(G4double prob, G4ThreeVector _start_point1, G4ThreeVector _start_point2, G4double abs_len,
	G4MaterialPropertyVector* energy_spec, G4int num_of_sims)
{
	if (num_of_sims == 0)
	{
#ifdef DEBUG_MC_NODES
		G4cout << "MCN: set_MC_node() CONTINUOUS none to simulate" << G4endl;
#endif
		return;
	}
#ifdef DEBUG_MC_NODES
	G4cout << "MCN: set_MC_node() CONTINUOUS" << G4endl;
#endif
	chosen_type = MC_NODE_CONTINIOUS;
	node_probability = prob;
	to_simulate = num_of_sims;
	simulated = simulated == -1 ? -1 : simulated;
	start_point1 = _start_point1;
	start_point2 = _start_point2;
	PhotonEnergy = -1;
	EnergySpectrum = energy_spec;
}

void MC_node::clear_prob_calc()
{
#ifdef DEBUG_MC_NODES
	G4cout << "MCN: clear_prob_calc()" << G4endl;
#endif
	is_accounted = 0;
	simulation_data.clear_prob_calc();
}

G4double MC_node::hit_prob() //net hit probability of node
{
	is_accounted = 1;
	G4double g = node_probability*simulation_data.hit_prob();
#ifdef DEBUG_MC_NODES
	G4cout << "MCN" << chosen_type << ": hit_prob() "<<g<< G4endl;
#endif
	return g;
}
