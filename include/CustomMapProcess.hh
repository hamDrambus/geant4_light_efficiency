#ifndef CustomMapProcess_hh
#define CustomMapProcess_hh 1

#include "G4VProcess.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "CustomParticleChangeForTransport.hh"
#include "CustomRunManager.hh"

class G4SafetyHelper;

class CustomMapProcess : public G4VProcess
{
	// Concrete class that does the geometrical transport 

public:  // with description

	CustomMapProcess(G4int verbosityLevel = 1);
	~CustomMapProcess();

	G4double      AlongStepGetPhysicalInteractionLength(
		const G4Track& track,
		G4double  previousStepSize,
		G4double  currentMinimumStep,
		G4double& currentSafety,
		G4GPILSelection* selection
		);

	G4VParticleChange* AlongStepDoIt(
		const G4Track& track,
		const G4Step& stepData
		);

	G4VParticleChange* PostStepDoIt(
		const G4Track& track,
		const G4Step&  stepData
		);
	// Responsible for the relocation.

	G4double PostStepGetPhysicalInteractionLength(
		const G4Track&,
		G4double   previousStepSize,
		G4ForceCondition* pForceCond
		);
	// Forces the PostStepDoIt action to be called, 
	// but does not limit the step.

	inline void   SetVerboseLevel(G4int verboseLevel);
	inline G4int  GetVerboseLevel() const;
	// Level of warnings regarding eg energy conservation
	// in field integration.

public:  // without description

	G4double AtRestGetPhysicalInteractionLength(
		const G4Track&,
		G4ForceCondition*
		) {
		return -1.0;
	};
	// No operation in  AtRestDoIt.

	G4VParticleChange* AtRestDoIt(
		const G4Track&,
		const G4Step&
		) {
		return 0;
	};
	// No operation in  AtRestDoIt.

private:

	G4Navigator*         fLinearNavigator;

	G4ThreeVector        fTransportEndPosition;
	G4ThreeVector        fTransportEndMomentumDir;
	G4bool               fMomentumChanged;
	// The particle's state after this Step, Store for DoIt

	G4bool               isMapped;  // Flag to determine whether a boundary was reached.

	G4TouchableHandle    fCurrentTouchableHandle;

	G4ThreeVector  fPreviousSftOrigin;
	G4double       fPreviousSafety;
	// Remember last safety origin & value.

	CustomParticleChangeForTransport fParticleChange;
	// New ParticleChange

	G4double fEndPointDistance;

	G4SafetyHelper* fpSafetyHelper;  // To pass it the safety value obtained

	// Verbosity 
	G4int    fVerboseLevel;
	// Verbosity level for warnings
	// eg about energy non-conservation in magnetic field.

};

inline void CustomMapProcess::SetVerboseLevel(G4int verboseLev)
{
	fVerboseLevel = verboseLev;
}

inline G4int CustomMapProcess::GetVerboseLevel() const
{
	return fVerboseLevel;
}

#endif  
