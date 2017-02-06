#ifndef Simulation_Summary_hh
#define Simulation_Summary_hh
#include "MC_Node.hh"
#include "GlobalDefinitions.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

class net_sim_data;
//stores all simulation processed results (specta, probabilities, etc.)
class SimulationSummary
{
public:
	G4double total_probability;
	G4double no_reemission_probability;
	G4double reemission_probability;

	G4double En_min;
	G4double En_max;
	G4int num_of_En_bins; 
	
	G4double *tot_spec;
	G4double *no_reemiss_spec;
	G4double *with_reemiss_spec;
	//^differential spectra of probabilities over energies
	G4int curr_weight_pbs; //required for updating all probs and spectra by a new simulation sequence
	G4int curr_weight_spec;//different just in case
	//^same as number of already cosidered simulations.
	virtual void UpdateProbabilities(net_sim_data* primary_node); //called from Manager in OnNewSequenceProc
	virtual void UpdateSpectra	(net_sim_data* primary_node); //TODO: write this function (mc_node procs (hooks) required) 

	virtual G4int Num_of_events();

	SimulationSummary(void);
	~SimulationSummary();
	virtual void ExportEnSpectraToFile(G4String str1 = "total_spectrum.txt", G4String str2 = "no_reemiss_spectrum.txt",
		G4String str3 = "reemissed_spectrum.txt");
	virtual void ExportLSpectraToFile(G4String str1 = "l_total_spectrum.txt", G4String str2 = "l_no_reemiss_spectrum.txt",
		G4String str3 = "l_reemissed_spectrum.txt");
	virtual void ClearAllData(void);
	friend  std::ostream& operator<<(std::ostream& str, const SimulationSummary& data);
protected:
	//return number of succesful simulations
	virtual G4int GetProbabilities(G4double *no_reemiss_pb, G4double *reemiss_pb, G4double* tot_pb, net_sim_data *node);
	virtual G4int GetProbSpectra(G4double *no_re_spec, G4double *reem_spec, G4double* tot_spec, net_sim_data *node,
		G4double En_min,G4double En_max, G4int Nbins); //TODO: write this function (mc_node procs (hooks) required) 
};

#endif