#include "SimulationSummary.hh"

SimulationSummary::SimulationSummary(void)
{
	num_of_En_bins = 300;
	En_max = 9.69*eV;
	En_min = 3 * eV;
	no_reemiss_spec = new G4double [num_of_En_bins];
	with_reemiss_spec = new G4double[num_of_En_bins];
	tot_spec = new G4double[num_of_En_bins];
	ClearAllData();
}

SimulationSummary::~SimulationSummary()
{
	if (no_reemiss_spec) delete [] no_reemiss_spec;
	if (with_reemiss_spec) delete [] with_reemiss_spec;
	if (tot_spec) delete [] tot_spec;
}

void SimulationSummary::ExportEnSpectraToFile(G4String str1, G4String str2,G4String str3)
{
	std::ofstream file1, file2, file3;
	file1.open(str1,std::ios_base::trunc);
	file2.open(str2, std::ios_base::trunc);
	file3.open(str3, std::ios_base::trunc);
	G4double dEn = (En_max - En_min)/num_of_En_bins;
	for (G4int g = 0; g < num_of_En_bins; g++)
	{
		G4double En = dEn/2 + g*dEn; //centered bins
		file1 << En/eV << "\t" << no_reemiss_spec[g] << G4endl;
		file2 << En/eV << "\t" << with_reemiss_spec[g] << G4endl;
		file3 << En/eV << "\t" << tot_spec[g] << G4endl;
	}
	file1.close();
	file2.close();
	file3.close();
}

void SimulationSummary::ExportLSpectraToFile(G4String str1, G4String str2, G4String str3)
{
	std::ofstream file1, file2, file3;
	file1.open(str1, std::ios_base::trunc);
	file2.open(str2, std::ios_base::trunc);
	file3.open(str3, std::ios_base::trunc);
	G4double dEn = (En_max - En_min) / num_of_En_bins;
	for (G4int g = num_of_En_bins-1; g >= 0; g--)
	{
		G4double En = dEn / 2 + g*dEn; //centered bins
		En = 1.2398e3 / (En / eV); //[nm]
		file1 << En << "\t" << no_reemiss_spec[g] << G4endl;
		file2 << En << "\t" << with_reemiss_spec[g] << G4endl;
		file3 << En << "\t" << tot_spec[g] << G4endl;
	}
	file1.close();
	file2.close();
	file3.close();
}

void SimulationSummary::ClearAllData(void)
{
	no_reemission_probability = -1;
	reemission_probability = -1;
	total_probability = -1;
	curr_weight_pbs = 0;
	curr_weight_spec = 0;
	for (G4int h = 0; h < num_of_En_bins; h++)
	{
		tot_spec[h] = 0;
		no_reemiss_spec[h] = 0;
		with_reemiss_spec[h] = 0;
	}
}

G4int SimulationSummary::Num_of_events()
{
	return curr_weight_pbs;
}
//\/returns number of succesful simulations
G4int SimulationSummary::GetProbabilities(G4double *no_reemiss_pb, G4double *reemiss_pb, G4double* tot_pb, net_sim_data *node)
{
	if (NULL == node)
		return -1;
	*no_reemiss_pb = 0;
	*reemiss_pb = 0;
	*tot_pb = 0;
	node->hit_prob(no_reemiss_pb, reemiss_pb, tot_pb);
	return node->num_of_succ_sims; 
}

G4int SimulationSummary::GetProbSpectra(G4double *no_re_spec, G4double *reem_spec, G4double* tot_spec, net_sim_data *node,
	G4double En_min, G4double En_max, G4int Nbins) //TODO: write this function (mc_node procs (hooks) required)
{
	if (NULL == node)
		return -1;
	return 0;
}

void SimulationSummary::UpdateProbabilities(net_sim_data* primary_node) //called from Manager in OnNewSequenceProc
{
	G4double tot_pb;
	G4double no_reemiss_pb;
	G4double reemiss_pb;
	G4int weight = GetProbabilities(&no_reemiss_pb, &reemiss_pb, &tot_pb, primary_node);
	if (weight <= 0)
		return;
	if (0>=curr_weight_pbs)
	{
		curr_weight_pbs = weight;
		total_probability = tot_pb;
		no_reemission_probability = no_reemiss_pb;
		reemission_probability = reemiss_pb;
		return;
	}
	total_probability = (total_probability*curr_weight_pbs + tot_pb) / (curr_weight_pbs + weight);
	no_reemission_probability = (no_reemission_probability*curr_weight_pbs + no_reemiss_pb) / (curr_weight_pbs + weight);
	reemission_probability = (reemission_probability*curr_weight_pbs + reemiss_pb) / (curr_weight_pbs + weight);
	curr_weight_pbs += weight;
	return;
}

void SimulationSummary::UpdateSpectra(net_sim_data* primary_node)
{
}

std::ostream& operator<<(std::ostream& str, const SimulationSummary& data)
{
	str << "total efficieency: " << data.total_probability<< G4endl;
	str << "not reemissed efficieency: " << data.no_reemission_probability<< G4endl;
	str << "remissed efficieency: " << data.reemission_probability<< G4endl;
	return str;
}
