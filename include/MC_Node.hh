#ifndef MC_NODE_HH
#define MC_NODE_HH

//TODO: optimize get_hit_probb(double*,double*,double*) (excessive checks at the moment)

#include "PseudoMesh.hh"
#include "GlobalDefinitions.hh"
#include "G4GeometryTolerance.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <list>

#define MC_NODE_UNDEFINED 0
#define MC_NODE_DESCRETE 1
#define MC_NODE_CONTINIOUS 2

class MC_node;
class one_sim_data;
class net_sim_data;
class PseudoMeshData;

class photon_event
{
public:
	photon_event(one_sim_data* p_container) :/* pre_step_pos(), post_step_pos(), pre_volume(NULL), post_volume(NULL),*/
		continue_prob(1), net_prob(1), daughter_node(0), has_hit(-1), chosen_event_type(RM_PHOTON_UNDEFINED), container(p_container), energy(-1*eV){};
	~photon_event();
	//had to decrease memory consumption
	//G4ThreeVector pre_step_pos;
	//G4ThreeVector post_step_pos;
	//G4VPhysicalVolume* pre_volume;
	//G4VPhysicalVolume* post_volume;
	G4double energy;

	G4double continue_prob;
	G4double net_prob;
	MC_node* daughter_node; //this poiters are just a copy, mustn't be deleted
	one_sim_data* container;
	G4int has_hit; //-1, 0, 1
	G4int chosen_event_type;
};

class one_sim_data
{
public:
	net_sim_data* container;
	std::list<photon_event> events; //for debugging actually, and stores end probability of event sequence
	
	one_sim_data(net_sim_data* p_contaier) :container(p_contaier){};
	~one_sim_data();
	photon_event* new_event(); //returns NULL if last event has hit (0 or 1), managed in net_sim_data
	void set_event(const G4Step* step, G4double prob, G4int ev_type, G4int is_hit);
	G4int is_hit();
	G4int has_hit;
	void SetHit(G4bool is_hit); //set hit for last event, no extra logic here

	//G4double probability; //if -1 - not calculated yet
	G4double hit_prob(); //recursive function
	void one_sim_data::hit_prob(G4double* no_reemiss, G4double *reemissed, G4double* total);
	//void clear_prob_calc();
};

class net_sim_data
{
public:
	net_sim_data(MC_node* parent);
	~net_sim_data();
	std::list<one_sim_data> events;
	MC_node* parent;
	void new_sequence(void); //managed by MC_node
	photon_event* new_event();
	void set_event(const G4Step* step, G4double prob, G4int ev_type, G4int is_hit=-1); //do not generate new event, thus new_event must be called before
	void SetHit(G4bool is_hit); //does not call new_sequence
	MC_node* find_next_to_simulate(void);

	G4int num_of_succ_sims;
	G4double probability; //if -1 - not calculated yet
	G4double hit_prob(); //recursive function, average here
	void hit_prob(G4double* no_reemiss, G4double *reemissed, G4double* total); //recursive function, average here
	void clear_prob_calc();
};

class MC_node
{
public:
	//generating conditions:
	G4int to_simulate;
	G4int simulated; //is -1 at the start 'cause it is ++ in new_seqeunce which is called obviuosly _before_ full sequence is proceded
	G4double node_probability; //cosnt
	G4ThreeVector start_point1;
	G4ThreeVector start_point2; //if ==start_point1, from fixed point
	G4double spacial_disrt_par; //==absorption length //units?
	G4ThreeVector PhotonMomentum; //if =0 - isotropic distr
	G4ThreeVector PhotonPolarization;
	G4double PhotonEnergy;
	G4MaterialPropertyVector* EnergySpectrum;

	PseudoMeshData* n_mapping_data; //created in set node, delteted in destructor

	photon_event* ev_parent;
	one_sim_data* parent_container;
	net_sim_data* parent_cont_cont;
	MC_node* gl_parent; //contains photon_event;
	//wrong: one node may contain several MC_node* daughter_node;//deleted from its arent events, this just for a convinience
	G4int chosen_type;
	net_sim_data simulation_data;

	G4ThreeVector	GenPostion();
	G4ThreeVector	GenMomentum();
	G4ThreeVector	GenPolarization();
	G4double		GenEnergy();
	PseudoMeshData*	GenMappingData();

	MC_node(photon_event* parent);
	~MC_node();

	//better not be called after the first sumulation is started
	void set_MC_node(G4double prob,PseudoMeshData* mapping_state, G4ThreeVector _start_point1, G4ThreeVector _start_point2, G4double abs_len, 
		G4MaterialPropertyVector* energy_spec, G4int num_of_sims = 2);
	//better not be called after the first sumulation is started
	void set_MC_node(G4double prob, PseudoMeshData* mapping_state, G4ThreeVector _start_point1, G4ThreeVector momentum, 
		G4ThreeVector polariztion, G4double energy, G4int num_of_sims = 1);
	void SetHit(G4bool is_hit);
	photon_event* new_event(MC_node **pp_new_MC_node); //if pp_new_MC_node set no *NULL, then new event is crated in current MC_node 
	G4int is_over(void); //all simulations are over

	//G4int is_accounted; //0 - it will be used in hit_prob() of 
	//void clear_prob_calc();
	G4double hit_prob(); //net hit probability of node
	void hit_prob(G4double* no_reemiss, G4double *reemissed, G4double* total);

	int is_reemissed(void);
};

#endif
