#ifndef CustomRunManager_h
#define CustomRunManager_h 1

#include "G4SystemOfUnits.hh"
#include "PseudoMesh.hh"
#include "SimulationSummary.hh"
#include "MC_Node.hh"
#include "B1DetectorConstruction.hh"
#include "GlobalDefinitions.hh"
#include "G4RunManager.hh"
#include <list>

class MC_node;
class photon_event;
class one_sim_data;
class net_sim_data;
class PseudoMeshData;
class SimulationSummary;

//TODO: this class is becoming too complex and complicated thus is should be split
//Current mandates:
//	0) photon with statistical weights (probabilities) managment
//	1) Setting particle gun parameter accordingly to simulatoion state (MC's nodes managment)
//	2) Storing current simulation state (curr_photon_event, curr probability and so on)
//	3) Interface for modifying sim. state and updating current state accordingly
//	4) Storing mapping state and transitional functiion for work with it (being called from process,
//	call detector construction's methods)
//	5) Analisis of simulation
//	6) Supress growths of secondary photons (number of combinations) - this option is shared with detector construction
// issue: fuctions for primary simulation MC node are somewhat duplicated from MC_node class
// Also managing of MC node be better be reworked with proper interface

class CustomRunManager:public G4RunManager
{
protected:
	G4ThreeVector initial_position;
	G4ThreeVector initial_momentum_direction;
	G4ThreeVector initial_polarization;
	G4double initial_energy;

	G4double curr_ph_ev_prob;
	G4double curr_ph_ev_type;
	G4int has_finished_secondaries;

	G4int is_first_process; //for end process
	G4double prev_ph_ev_prob; //cause of sequence of processes, see process_end()
	G4double prev_ph_ev_type; //^
public: 
	G4int extra_run_id;
	//G4int current_sim_status;
	
	MC_node* current_working_node; //there secondary generation info is stored too
	photon_event* last_event; //current working event - for convinience
	PseudoMeshData* curr_mapping_state; //set from PreEventProcedure either based on primary photon posions or
	//mapping state of MC_Node
	
	SimulationSummary* sim_results;

	G4int spawn_new_MC_node(const G4Step* step, G4double prob, G4ThreeVector momentum, G4ThreeVector polarization, G4int num_of_sims = 1);
#if defined(TOP_MESH_TEST)||defined(TEST_MESH_SIDEWAYS)
	G4int spawn_new_MC_node(const G4Step* step, G4double prob, G4Material *WSL_pars, G4int num_of_sims =0);
#else
	G4int spawn_new_MC_node(const G4Step* step, G4double prob, G4Material *WSL_pars, G4int num_of_sims = 150);//150
#endif
#ifdef TOP_MESH_TEST
	std::list<G4double> top_hits_xs, top_hits_ys, top_hits_probs;
	std::list<G4double> bot_hits_xs, bot_hits_ys, bot_hits_probs;
	G4int x_num, y_num; //discretisation parameters (same as size of bmp)
	G4double x_start, y_start;
	G4double t_step, t_uncert;
	G4int t_counter;
	void on_hit_proc(G4ThreeVector point,G4double prob,G4int top_bot); //called when test detector is hit with photon
	void export_to_bmp(std::list<G4double>* hits_x, std::list<G4double> *hits_y, std::list<G4double> *probs, G4String filename);
#endif
	//below manage memory allocation for photon_events* chains (list that is)
	//so that every photon process is written at the right place 
	G4int init_event(); //called at a start of event (=simulation)
	G4int next_event(const G4Step* step); //called in UserSteppingAction (), writing data is here too
	G4int process_end(const G4Step* step); //called at the end of the process when something does happen to a photon
	//its purpose is to manage next_event generation in case of several sequent processes (because StepUserAction is not called between)
	//It generarates next_events for all but last process in sequence (for the last UserSteppingAction is used)
	//So it effectively calls next_event for the previous process to the one it is entered from 
	//('step' is the same for all process in the sequence)
	//also new_event() in next_event() must not return NULL in this function (because there is no hit inside processes, only in UserSteppingAction)
	//WARNING! must be called BEFORE spawn_new_MC_node
	//P.S. ugly

	virtual void OnNewSimulationProc(void); //called right before calling PrimaryMCnode->new_sequence
	//its present purpose is to decrease memory consumption by writing sim data and removing simulated sequence

	G4int is_secondary_MC(void); //true if there are secondary MC (MC_nodes) to be simulated
	net_sim_data* primary_Monte_Carlo; //has no parent MC - stores data (sequences) for initial simulation
	//number of simulations there^ is number of events in beamOn()

	void SetPhEvType(G4int evType);
	void SetPhEvProb(G4double evProb);
	void SetHit(G4bool is_hit=true);
	G4double get_total_detetion_eff(void);
	void get_total_detetion_eff(G4double* no_reemiss, G4double *reemissed, G4double* total);
	CustomRunManager() :G4RunManager(), extra_run_id(1), curr_ph_ev_prob(1), curr_ph_ev_type(RM_PHOTON_UNDEFINED),
		has_finished_secondaries(1), current_working_node(NULL), last_event(NULL), primary_Monte_Carlo(NULL) 
	{
		curr_mapping_state = new PseudoMeshData;
		sim_results = new SimulationSummary;
#ifdef TOP_MESH_TEST
		x_num = 600;
		y_num = 450;
		x_start = 8 * mm / 2;
		y_start = 0 * mm / 2;
		t_step = 0.1*mm;
		t_uncert = 0.05*mm;
		t_counter = 0;
#endif
	};
	~CustomRunManager()
	{
		if (primary_Monte_Carlo) delete primary_Monte_Carlo;
		if (curr_mapping_state) delete curr_mapping_state;
		if (sim_results) delete sim_results;
	};
	//depr - went to init_event //virtual void Initialize(); //clear MC_nodes and create parent one, clear all previous sim data
	//deleting last MC_nodes and sim data - in destructor!

    virtual void DoEventLoop(G4int n_event,const char* macroFile=0,G4int n_select=-1);
	//selects what process should be tracked further (called in OpticsBoundaryProcess)
	virtual G4int select_photon_BP(const G4Step* step, G4ThreeVector defl_momentum, G4ThreeVector refl_momentum); 
	//0 - kill process, 1 - select reflection, 2 - select difraction, 3 - select both
	G4ThreeVector	GenPosition(); //Gen called only at the start of event
	G4ThreeVector	GenMomentum();
	G4ThreeVector	GenPolarization();
	G4double		GenEnergy();
	PseudoMeshData*	GenMappingState();

	G4ThreeVector	FetchPosition();
	G4ThreeVector	FetchMomentum();
	G4ThreeVector	FetchPolarization();
	G4double		FetchEnergy();

	void OnEventStartProc();
	virtual void ProcessOneEvent(G4int i_event); //added OnEventStartProc call here
	//^ for setting parameters before call of GeneratePrimaries (EventAction doesn't fit this purpose)

	G4double get_new_spawn_prob();//called in spawn_new_MC_node and prevents the tree from being too deep
	//aborts spawnig if probability of a new node is less than MIN_ALLOWED_PROBABILITY
	G4int get_sim_depth(); //affects number of simulation for secondary MC's - the deeper the less
	G4int get_sim_depth(G4int node_type); //same as above but seaches just particular type of nodes (in order to prevent double reemission)
	G4double get_curr_event_probab();
	G4int get_detected_spectrum(std::list<G4double> *energies, std::list < G4double> *probabilities);
	void get_detected_spectrum(std::string out_filename = "WLS_spectrum_test.txt");
	//recursive function below
	void get_detected_spectrum(G4double reach_node_prob, MC_node* node, std::list<G4double> *energies, std::list < G4double> *probabilities);

	G4ThreeVector MappingProc(const G4Track& track, const G4Step& aStep, G4TouchableHandle &fCurrentTouchableHandle);
};

#endif

