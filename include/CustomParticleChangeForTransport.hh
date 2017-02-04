//the only difference for standart G4ParticleChangeTransport is that it changes
//position in UpdateForPostStep
//Could not inherit from G4ParticleChangeTransport because it contains important private members
//used in virtual UpdateFor...'s

#ifndef CustomParticleChangeForTransport_h
#define CustomParticleChangeForTransport_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4TouchableHandle.hh"
#include "G4ParticleChange.hh"

class G4MaterialCutsCouple;
class G4VSensitiveDetector;

class CustomParticleChangeForTransport : public G4ParticleChange
{
public:
	// default constructor
	CustomParticleChangeForTransport();

	// destructor
	virtual ~CustomParticleChangeForTransport();

protected:
	// hide copy constructor and assignment operator as protected
	CustomParticleChangeForTransport(const CustomParticleChangeForTransport &right);
	CustomParticleChangeForTransport & operator=(const CustomParticleChangeForTransport &right);

public: // with description
	// ----------------------------------------------------
	// --- the following methods are for updating G4Step -----   
	// Return the pointer to the G4Step after updating the Step information
	// by using final state information of the track given by a physics
	// process    
	virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
	virtual G4Step* UpdateStepForAtRest(G4Step* Step);
	virtual G4Step* UpdateStepForPostStep(G4Step* Step);
	// A physics process gives the final state of the particle 
	// based on information of G4Track (or equivalently the PreStepPoint)

	virtual void Initialize(const G4Track&);
	// Initialize all propoerties by using G4Track information

	// ----------------------------------------------------
	//--- methods to keep information of the final state--
	//  IMPORTANT NOTE: Although the name of the class and methods are
	//   "Change", what it stores (and returns in get) are the "FINAL" 
	//   values of the Position, Momentum, etc.

	const G4TouchableHandle& GetTouchableHandle() const;
	void  SetTouchableHandle(const G4TouchableHandle& fTouchable);
	//  Get/Set the touchable of the current particle.
	//  Note: Touchable in PostStepPoint will be updated only after PostStepDoIt

	G4Material* GetMaterialInTouchable() const;
	void SetMaterialInTouchable(G4Material* fMaterial);
	//  Get/Propose the material in the touchable of the current particle.

	const G4MaterialCutsCouple* GetMaterialCutsCoupleInTouchable() const;
	void SetMaterialCutsCoupleInTouchable(const G4MaterialCutsCouple* fMaterialCutsCouple);
	//  Get/Set the materialCutsCouple in the touchable of the current particle.

	G4VSensitiveDetector* GetSensitiveDetectorInTouchable() const;
	void SetSensitiveDetectorInTouchable(G4VSensitiveDetector* fSensitiveDetector);
	//  Get/Set the sensitive detector in the touchable of the current particle.

	G4bool GetMomentumChanged() const;
	void SetMomentumChanged(G4bool b);

public:
	virtual void DumpInfo() const;

protected:
	G4TouchableHandle theTouchableHandle;
	//  The changed touchable of a given particle.

public:

	// Prototype implementation of smooth representation of curved trajectories.
	// Auxiliary points are ThreeVectors for now; change to G4AuxiliaryPoints.

	inline void SetPointerToVectorOfAuxiliaryPoints(std::vector<G4ThreeVector>* theNewVectorPointer);
	inline std::vector<G4ThreeVector>* GetPointerToVectorOfAuxiliaryPoints() const;

private:
	G4bool     isMomentumChanged;
	//  The flag which is set if momentum is changed in current step
	G4Material* theMaterialChange;
	const G4MaterialCutsCouple* theMaterialCutsCoupleChange;
	G4VSensitiveDetector* theSensitiveDetectorChange;
	// The material (and MaterialCutsCouple) where given track
	// currently locates

private:
	std::vector<G4ThreeVector>* fpVectorOfAuxiliaryPointsPointer;
};

inline
void CustomParticleChangeForTransport::SetTouchableHandle(const G4TouchableHandle&
fTouchable)
{
	theTouchableHandle = fTouchable;
}

inline
const G4TouchableHandle& CustomParticleChangeForTransport::GetTouchableHandle() const
{
	return theTouchableHandle;
}


inline
void CustomParticleChangeForTransport::SetMaterialInTouchable(G4Material* fMaterial)
{
	theMaterialChange = fMaterial;
}

inline
G4Material* CustomParticleChangeForTransport::GetMaterialInTouchable() const
{
	return theMaterialChange;
}

inline
void CustomParticleChangeForTransport::SetMaterialCutsCoupleInTouchable(const G4MaterialCutsCouple* fMaterialCutsCouple)
{
	theMaterialCutsCoupleChange = fMaterialCutsCouple;
}

inline
const G4MaterialCutsCouple* CustomParticleChangeForTransport::GetMaterialCutsCoupleInTouchable() const
{
	return theMaterialCutsCoupleChange;
}

inline
void CustomParticleChangeForTransport::SetSensitiveDetectorInTouchable(G4VSensitiveDetector* fSensitiveDetector)
{
	theSensitiveDetectorChange = fSensitiveDetector;
}

inline
G4VSensitiveDetector* CustomParticleChangeForTransport::GetSensitiveDetectorInTouchable() const
{
	return theSensitiveDetectorChange;
}

inline
G4bool CustomParticleChangeForTransport::GetMomentumChanged() const
{
	return isMomentumChanged;
}

inline
void CustomParticleChangeForTransport::SetMomentumChanged(G4bool b)
{
	isMomentumChanged = b;
}

//----------------------------------------------------------------
// functions for Initialization
//

inline void CustomParticleChangeForTransport::Initialize(const G4Track& track)
{
	// use base class's method at first
	InitializeStatusChange(track);
	//  InitializeLocalEnergyDeposit(track);
	InitializeSteppingControl(track);
	//  InitializeTrueStepLength(track);
	//  InitializeSecondaries(track);

	// set Energy/Momentum etc. equal to those of the parent particle
	const G4DynamicParticle*  pParticle = track.GetDynamicParticle();
	//  theEnergyChange          = pParticle->GetKineticEnergy();
	//  theMomentumChange        = pParticle->GetMomentumDirection();
	theVelocityChange = track.GetVelocity();
	isVelocityChanged = false;
	thePolarizationChange = pParticle->GetPolarization();
	//  theProperTimeChange      = pParticle->GetProperTime();

	// set Position/Time etc. equal to those of the parent track
	//  thePositionChange      = track.GetPosition();
	// set TimeChange equal to local time of the parent track
	theTimeChange = track.GetLocalTime();
	// set initial Local/Global time of the parent track
	theLocalTime0 = track.GetLocalTime();
	theGlobalTime0 = track.GetGlobalTime();

	// set touchable equal to the next touchable of the parent track
	// not set as for now
	//theTouchableChange     = track.GetNextTouchable();

	// So almost nothing is initialized here.
	// theMomentumChange, theProperTimeChange, thePositionChange and theTimeChange
	// are set by G4Transportation::AlongStepDoIt;
	// the others are not needed.
	// Take care when implementing the PostStep related things!
	// (P. Urban)
}

// Prototype implementation of smooth representation of curved trajectories.

inline void
CustomParticleChangeForTransport::
SetPointerToVectorOfAuxiliaryPoints(std::vector<G4ThreeVector>*
theNewVectorPointer)
{
	fpVectorOfAuxiliaryPointsPointer = theNewVectorPointer;
}

inline std::vector<G4ThreeVector>*
CustomParticleChangeForTransport::GetPointerToVectorOfAuxiliaryPoints() const
{
	return fpVectorOfAuxiliaryPointsPointer;
}


#endif
