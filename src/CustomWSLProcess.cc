
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

#include "CustomWLSProcess.hh"
#include "G4GeometryTolerance.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ParallelWorldProcess.hh"

#include "G4SystemOfUnits.hh"

CustomWLSProcess::CustomWLSProcess(const G4String& processName,
	G4ProcessType type)
	: G4VDiscreteProcess(processName, type)
{
	SetVerboseLevel(0);
	if (verboseLevel > 0) {
		G4cout << GetProcessName() << " is created " << G4endl;
	}

	SetProcessSubType(fOpAbsorption);
	theStatus = Defalt;
}

CustomWLSProcess::~CustomWLSProcess(){}


// PostStepDoIt
// ------------
//

G4VParticleChange*
CustomWLSProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
#ifdef DEBUG_MC_NODES
	G4cout << "CWLSP: wls PostStepDoIt()" << G4endl;
#endif
	theStatus = Defalt;

	aParticleChange.Initialize(aTrack);
	aParticleChange.ProposeVelocity(aTrack.GetVelocity());
	const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
	aParticleChange.ProposeMomentumDirection(aParticle->GetMomentumDirection());
	aParticleChange.ProposePolarization(aParticle->GetPolarization());

	// Get hyperStep from  G4ParallelWorldProcess
	//  NOTE: PostSetpDoIt of this process should be
	//        invoked after G4ParallelWorldProcess!

	const G4Step* pStep = &aStep;
	const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();

	if (hStep) pStep = hStep;

	//G4bool isOnBoundary =
	//	(pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
	G4Material* WLS_mat = pStep->GetPreStepPoint()->GetMaterial();
	
	if (verboseLevel > 0) {
		G4cout << "WSL absorbtion simulation" << G4endl;
	}

	if (aTrack.GetStepLength() <= G4GeometryTolerance::GetInstance()
		->GetSurfaceTolerance()/ 2){
		theStatus = StepTooSmallet;
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	G4double photon_energy = aParticle->GetTotalEnergy();

	G4ThreeVector EndPoint = pStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector StartPoint = pStep->GetPreStepPoint()->GetPosition();
#ifdef DEBUG_MC_NODES
	G4cout << "CWLSP: start point: "<<StartPoint << G4endl;
	G4cout << "CWLSP: end point: " << EndPoint << G4endl;
#endif

	G4MaterialPropertiesTable* aMaterialPropertiesTable;
	G4MaterialPropertyVector* abs_len;
	G4MaterialPropertyVector* wls_efficiency;
	G4double absorb_l, wls_eff;

	aMaterialPropertiesTable = WLS_mat->GetMaterialPropertiesTable();
	if (aMaterialPropertiesTable) {
		abs_len = aMaterialPropertiesTable->GetProperty("ABSORBTION_LENGTH");
		wls_efficiency = aMaterialPropertiesTable->GetProperty("WLS_EFFICIENCY");
	}
	else {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	if ((abs_len)&&(wls_efficiency)) {
		absorb_l = GetMaterialAbsLength(abs_len,photon_energy);
		if (absorb_l==DBL_MAX)
			return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep); //TODO: error?
		wls_eff= wls_efficiency->Value(photon_energy);
	}
	else {
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	G4double path_l=(EndPoint- StartPoint).getR();
	G4double reemiss_prob = wls_eff*(1.0-exp(-path_l / absorb_l));
	G4double continue_prob = exp(-path_l / absorb_l);
	CustomRunManager* manman;
	manman = (CustomRunManager*)(G4RunManager::GetRunManager());
	manman->SetPhEvType(RM_PHOTON_PROPOG);
	manman->SetPhEvProb(continue_prob);
	manman->process_end(&aStep);
	if (reemiss_prob!=0) manman->spawn_new_MC_node(pStep, reemiss_prob, WLS_mat);

	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
// ---------------
//
G4double CustomWLSProcess::GetMeanFreePath(const G4Track&,
	G4double,
	G4ForceCondition* condition)
{
	*condition = Forced;
	return DBL_MAX;
}

G4bool CustomWLSProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
	return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

//abs_len_property is implied tobe sorted by energy;
G4double CustomWLSProcess::GetMaterialAbsLength(G4MaterialPropertyVector* abs_len_prop, G4double energy)
{
	G4int N = abs_len_prop->GetVectorLength();
	if (N <= 0)
		return DBL_MAX;
	if (N == 1)
		return (*abs_len_prop)(0);//TODO: check 
	//G4double min_en=(*abs_len_prop)[0];
	G4double min_en = (*abs_len_prop).Energy(0);
	G4double max_en = (*abs_len_prop).Energy(N-1);
	if (energy < min_en)
		return (*abs_len_prop)(0);
	if (energy>max_en)
		return (*abs_len_prop)[N - 1];
	return (abs_len_prop->Value(energy));
}
