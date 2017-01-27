#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "CustomRunManager.hh"
#include "G4LogicalVolume.hh"


B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  stopPhotonVolume(0)
{}


B1SteppingAction::~B1SteppingAction()
{}

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{

  // get volume of the current step
	G4VPhysicalVolume* tremor = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	if (NULL==tremor)
	{// kind of error
		CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
		manman->SetHit(0);
		step->GetTrack()->SetTrackStatus(fStopAndKill);
		manman->next_event(step);
		return;
	}
  G4LogicalVolume* post_log_volume = tremor->GetLogicalVolume();
  CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
  if (step->GetTrack()->GetTrackStatus() == fStopAndKill) //in case some process killed the track, I need to handle events properly
  {
	  manman->SetHit(0);
	  manman->next_event(step);
	  return;
  }
  if (manman->get_curr_event_probab()< MIN_ALLOWED_STEPPING)
  {
	  manman->SetHit(0);
	  step->GetTrack()->SetTrackStatus(fStopAndKill);
	  manman->next_event(step);
	  return;
  }
  //WARNING! changed casting here
  B1DetectorConstruction* detectorConstruction = (B1DetectorConstruction*) (manman->GetUserDetectorConstruction());
  detectorConstruction->top_GEM->PostSpeppingAction(step); //changes the positoin in a manner of "teleportation"
  G4double detect_prob = detectorConstruction->GetHitProbability(step->GetPostStepPoint());
  if (0 == detect_prob)
  {
	  manman->SetHit(0);
	  step->GetTrack()->SetTrackStatus(fStopAndKill);
	  manman->next_event(step);
	  return;
  }
  if (1 == detect_prob)
  {
	  manman->SetHit(1);
	  step->GetTrack()->SetTrackStatus(fStopAndKill);
#ifdef TOP_MESH_TEST
	  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == detectorConstruction->top_mesh_test_detector)
		  manman->on_hit_proc(step->GetPostStepPoint()->GetPosition(), detect_prob);
#endif
	  manman->next_event(step);
	  return;
  }
  if (-1 == detect_prob)
  {
	  manman->next_event(step);
	  return;
  }
  //if detection probability is not a 0 or 1, then it's treated as process
#ifdef TOP_MESH_TEST
  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == detectorConstruction->top_mesh_test_detector)
	  manman->on_hit_proc(step->GetPostStepPoint()->GetPosition(), detect_prob);
#endif
  manman->SetPhEvProb(detect_prob);
  manman->SetPhEvType(RM_PHOTON_DETECTION);
  manman->process_end(step);
  //after "detection process" hit as usual
  manman->SetHit(1);
  step->GetTrack()->SetTrackStatus(fStopAndKill);
  manman->next_event(step);
}

