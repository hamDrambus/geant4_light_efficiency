#include "CustomMapProcess.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ProductionCutsTable.hh"

CustomMapProcess::CustomMapProcess(G4int verbosity)
	: G4VProcess(G4String("Mapping"), fGeneral),
	fTransportEndPosition(0.0, 0.0, 0.0),
	fTransportEndMomentumDir(0.0, 0.0, 0.0),
	fMomentumChanged(true),
	isMapped(false),
	fPreviousSftOrigin(0., 0., 0.),
	fPreviousSafety(0.0),
	// fParticleChange(),
	fEndPointDistance(-1.0),
	fVerboseLevel(verbosity)
{
	// set Process Sub Type
	pParticleChange = &fParticleChange;   // Required to conform to G4VProcess 

	G4TransportationManager* transportMgr;
	transportMgr = G4TransportationManager::GetTransportationManager();
	fLinearNavigator = transportMgr->GetNavigatorForTracking();
	fpSafetyHelper = transportMgr->GetSafetyHelper();  // New 

	// Cannot determine whether a field exists here, as it would 
	//  depend on the relative order of creating the detector's 
	//  field and this process. That order is not guaranted.
	// Instead later the method DoesGlobalFieldExist() is called

	static G4ThreadLocal G4TouchableHandle* pNullTouchableHandle = 0;
	if (!pNullTouchableHandle)  { pNullTouchableHandle = new G4TouchableHandle; }
	fCurrentTouchableHandle = *pNullTouchableHandle;
	// Points to (G4VTouchable*) 0
}

CustomMapProcess::~CustomMapProcess()
{}

// all deprecated - Touchable obtained via Navigator: Responsibilities:
//    sets fCurrentTouchableHandle to that of Transport process. Because that is not updated yet 
//	  in AlongStepDoIt PostStepPoint to contain real post step volume, but updated only after call of transport's PostStep
//    So "prediction" required, may be done due to track's correct status (boundary or not).
//	  TODO: WARNING!: does not account for field manager (==non-linear trajectories)
G4double CustomMapProcess:: AlongStepGetPhysicalInteractionLength(const G4Track&  track, G4double, //  previousStepSize
	G4double  currentMinimumStep, G4double& currentSafety,G4GPILSelection* selection)
{
	G4double safety = -1;
	G4double linearStepLength = fLinearNavigator->ComputeStep(track.GetPosition(),
		track.GetMomentumDirection(),currentMinimumStep, safety);
	*selection = CandidateForSelection;
	fTransportEndPosition = track.GetPosition() + linearStepLength*track.GetMomentumDirection();
	return DBL_MAX;
}

//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)
G4VParticleChange* CustomMapProcess::AlongStepDoIt(const G4Track& track,
	const G4Step&  stepData)
{
	fParticleChange.Initialize(track);
	fParticleChange.ProposeMomentumDirection(track.GetMomentumDirection());
	fParticleChange.ProposeEnergy(track.GetTotalEnergy());
	fParticleChange.SetMomentumChanged(false);
	fParticleChange.ProposePolarization(track.GetPolarization());
	fParticleChange.ProposeLocalTime(track.GetLocalTime());
	fParticleChange.ProposeProperTime(track.GetProperTime());
	fParticleChange.ProposePosition(stepData.GetPreStepPoint()->GetPosition());//not PostStep, because changes in AliongSteps are 
	//additive (+ with transport)
	
	//?? get Touchable by fTransportEndPosition but without harm to Navigator?
	//fLinearNavigator->CreateTouchableHistory();//may be done only in poststep
	//but I need new new position now, otherwise safety will be set wrongly after tranport along step?
	//fCurrentTouchableHandle = stepData.GetPostStepPoint()->GetTouchableHandle();
	//CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
	//B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(manman->GetUserDetectorConstruction());
	//if (fCurrentTouchableHandle!=0)
	//	fTransportEndPosition = detectorConstruction->top_GEM->PostSteppingAction(track, stepData, fCurrentTouchableHandle); //Here is the core of mapping
	//else fTransportEndPosition = stepData.GetPostStepPoint()->GetPosition();
	//fParticleChange.ProposePosition(fTransportEndPosition  - 
	//	(stepData.GetPostStepPoint()->GetPosition()-stepData.GetPreStepPoint()->GetPosition())); 
	////^may be the same. The subtraction is required in order to compesate for change in transport. (otherwise transport shift is doubling
	////because position changes in AlongStepDoIt are additive (see SteppingManager entrails))
	//G4VPhysicalVolume* temp_p = fLinearNavigator->LocateGlobalPointAndSetup(stepData.GetPreStepPoint()->GetPosition(), &track.GetMomentumDirection());
	//fCurrentTouchableHandle = fLinearNavigator->CreateTouchableHistory();
	//if (0 == temp_p) //something is not ok. TODO: this scenario should be tested.
	//	fCurrentTouchableHandle->UpdateYourself(temp_p);
	/*if ((fTransportEndPosition - stepData.GetPostStepPoint()->GetPosition()).r() > 0.0)
		isMapped = 1;
	else isMapped = 0;*/
	return &fParticleChange;
}

//  This ensures that the PostStep action is always called,
//  so that it can do the relocation if it is needed.
G4double CustomMapProcess::PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition* pForceCond)
{
	*pForceCond = Forced;
	return DBL_MAX;  // was kInfinity ; but convention now is DBL_MAX
}

G4VParticleChange* CustomMapProcess::PostStepDoIt(const G4Track& track, const G4Step& aStep)
{
	G4TouchableHandle retCurrentTouchable = aStep.GetPostStepPoint()->GetTouchableHandle();
	fParticleChange.ProposeTrackStatus(track.GetTrackStatus());
	fTransportEndPosition = aStep.GetPostStepPoint()->GetPosition();
	if (aStep.GetPostStepPoint()->GetStepStatus()!= fGeomBoundary)
		goto out_;
	fCurrentTouchableHandle = aStep.GetPostStepPoint()->GetTouchableHandle();
	CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
	//B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(manman->GetUserDetectorConstruction());
	if (fCurrentTouchableHandle!=0)
		fTransportEndPosition = manman->MappingProc(track, aStep, fCurrentTouchableHandle); //Here is the core of mapping
	else fTransportEndPosition = aStep.GetPostStepPoint()->GetPosition();
	if ((fTransportEndPosition - aStep.GetPostStepPoint()->GetPosition()).r() > 0.0)
		isMapped = 1;
	else isMapped = 0;
	if (isMapped) // If mapping has ocurred logically relocate the particle
	{
		//depr: fLinearNavigator->LocateGlobalPointAndUpdateTouchableHandle(fTransportEndPosition,track.GetMomentumDirection(),
		//	fCurrentTouchableHandle, true); //does not work for 'teleportation'

		//The following voodoo magic is manipulations with navigator to trick it into thinking that the new point to which 
		//mapping ocurred is on the boundary. This requires going from one touchable to neighbour one.
		//becuse new point is on the boundary, methods from navigator can get both touchables by taking two opposite directions
		//calls are taken from Navigator->LocateGlobalPointAndUpdateTouchableHandle
		G4double temp_s = fLinearNavigator->ComputeSafety(fTransportEndPosition);
		CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
		B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*)(manman->GetUserDetectorConstruction());
		G4double shift = detectorConstruction->GetSafeOffset();
		fLinearNavigator->SetGeometricallyLimitedStep();
		G4VPhysicalVolume* temp_p = fLinearNavigator->LocateGlobalPointAndSetup(fTransportEndPosition - shift*(track.GetMomentumDirection()),
			&-(track.GetMomentumDirection()), true);
		//^first geeting 'from' touchable hence -track...
		fCurrentTouchableHandle = fLinearNavigator->CreateTouchableHistory();
#ifdef TEMP_CODE_
		G4VPhysicalVolume* temp_pp_ =fCurrentTouchableHandle->GetVolume();
#endif
		if (0 == temp_p) //something is not ok. TODO: this scenario should be tested.
		{
			fCurrentTouchableHandle->UpdateYourself(temp_p);
			fLinearNavigator->SetGeometricallyLimitedStep();
			G4VPhysicalVolume* temp_p = fLinearNavigator->LocateGlobalPointAndSetup(fTransportEndPosition, &(track.GetMomentumDirection()), true);
			fCurrentTouchableHandle = fLinearNavigator->CreateTouchableHistory();
			if (0 == temp_p) //something is not ok. TODO: this scenario should be tested.
				fCurrentTouchableHandle->UpdateYourself(temp_p);
		}
		else
			fLinearNavigator->LocateGlobalPointAndUpdateTouchableHandle(fTransportEndPosition, track.GetMomentumDirection(),
			fCurrentTouchableHandle, true);
#ifdef TEMP_CODE_
		temp_pp_ = fCurrentTouchableHandle->GetVolume();
#endif
		if (fCurrentTouchableHandle == 0) // Check whether the particle is out of the world volume
		{
			fParticleChange.ProposeTrackStatus(fStopAndKill);
			fCurrentTouchableHandle = retCurrentTouchable;
		}
		else
			if (fCurrentTouchableHandle->GetVolume() == 0) // If so it has exited and must be killed.
				fParticleChange.ProposeTrackStatus(fStopAndKill);
		retCurrentTouchable = fCurrentTouchableHandle;
		//end of voodoo
	}
	out_:

	fParticleChange.SetTouchableHandle(retCurrentTouchable);
	fParticleChange.ProposePosition(fTransportEndPosition);
	//following is copied from transport process, seemingly sets materials
	const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume();
	const G4Material* pNewMaterial = 0;
	const G4VSensitiveDetector* pNewSensitiveDetector = 0;

	if (pNewVol != 0)
	{
		pNewMaterial = pNewVol->GetLogicalVolume()->GetMaterial();
		pNewSensitiveDetector = pNewVol->GetLogicalVolume()->GetSensitiveDetector();
	}
	// ( <const_cast> pNewMaterial ) ;
	// ( <const_cast> pNewSensitiveDetector) ;
	fParticleChange.SetMaterialInTouchable((G4Material *)pNewMaterial);
	fParticleChange.SetSensitiveDetectorInTouchable((G4VSensitiveDetector *)pNewSensitiveDetector);

	const G4MaterialCutsCouple* pNewMaterialCutsCouple = 0;
	if (pNewVol != 0)
		pNewMaterialCutsCouple = pNewVol->GetLogicalVolume()->GetMaterialCutsCouple();

	if (pNewVol != 0 && pNewMaterialCutsCouple != 0 && pNewMaterialCutsCouple->GetMaterial() != pNewMaterial)
	{
		// for parametrized volume
		pNewMaterialCutsCouple =
			G4ProductionCutsTable::GetProductionCutsTable()
			->GetMaterialCutsCouple(pNewMaterial,
			pNewMaterialCutsCouple->GetProductionCuts());
	}
	fParticleChange.SetMaterialCutsCoupleInTouchable(pNewMaterialCutsCouple);

	// Set the touchable in ParticleChange
	// this must always be done because the particle change always
	// uses this value to overwrite the current touchable pointer.
	fParticleChange.SetTouchableHandle(retCurrentTouchable);
	return &fParticleChange;
}
