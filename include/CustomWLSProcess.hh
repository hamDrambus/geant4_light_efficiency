
#ifndef CustomWSLProcess_h
#define CustomWSLProcess_h 1

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"

#include "G4RandomTools.hh"
#include "G4RandomDirection.hh"

#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"
#include "CustomRunManager.hh"

// Class Description:
// Discrete Process -- 
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:


enum CustomWLSProcessStatus {
	Defalt, StepTooSmallet
};

class CustomWLSProcess : public G4VDiscreteProcess
{

public:
	CustomWLSProcess(const G4String& processName = "OpWSL",
		G4ProcessType type = fOptical);
	~CustomWLSProcess();

private:

	CustomWLSProcess(const CustomWLSProcess &right);
	CustomWLSProcess& operator=(const CustomWLSProcess &right);

public:

	G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
	// Returns true -> 'is applicable' only for an optical photon.

	G4double GetMeanFreePath(const G4Track&,
		G4double,
		G4ForceCondition* condition);
	// Returns infinity; i. e. the process does not limit the step,
	// but sets the 'Forced' condition for the DoIt to be invoked at
	// every step. However, only at a boundary will any action be
	// taken.

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
		const G4Step&  aStep);
	// This is the method implementing boundary processes.

	CustomWLSProcessStatus GetStatus() const;
	// Returns the current status.
private:
	G4double GetMaterialAbsLength(G4MaterialPropertyVector* abs_len_prop, G4double energy);
	CustomWLSProcessStatus theStatus;
};

inline
CustomWLSProcessStatus CustomWLSProcess::GetStatus() const
{
	return theStatus;
}

#endif /* G4OpBoundaryProcess_h */
