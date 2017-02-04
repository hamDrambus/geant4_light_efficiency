#include "B1EventAction.hh"
#include "B1Run.hh"

#include "G4Event.hh"
#include "CustomRunManager.hh"


B1EventAction::B1EventAction()
: G4UserEventAction(),
  fEdep(0.)
{} 

B1EventAction::~B1EventAction()
{}

void B1EventAction::BeginOfEventAction(const G4Event* event)
{   
	/* depr - called after generate primaries: CustomRunManager* manman = (CustomRunManager*)(G4RunManager::GetRunManager());
	manman->OnEventStartProc(event);*/ //sets pesudoMeshes properly (depending on starting position of photon)
	//and generates new Gun's parameters
	fEdep = 0.;
}


void B1EventAction::EndOfEventAction(const G4Event* event)
{   
  //Draw trajectory
  G4VVisManager* pVis = G4VVisManager::GetConcreteInstance();
  if (pVis)
  {
	  G4TrajectoryContainer* pCont = event->GetTrajectoryContainer();
	  G4int n_of_traject = pCont?pCont->entries():0;
	  for (int ind = 0; ind < n_of_traject; ind++)
	  {
		  G4Trajectory* trj = (G4Trajectory*)((*pCont)[ind]);
		  trj->DrawTrajectory();
	  }
  }
  //CustomRunManager* manman = (CustomRunManager*)G4RunManager::GetRunManager();
  //manman->close_event();
}
