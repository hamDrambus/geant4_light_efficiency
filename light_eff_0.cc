#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "CustomRunManager.hh"
#endif

#include "G4UImanager.hh"
//#include "QBBC.hh"
#include "OpticsPhysicsList.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"
#include "G4VSteppingVerbose.hh"

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine(40));
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  CustomRunManager* runManager = new CustomRunManager;
#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetVerboseLevel(0);
  runManager->SetUserInitialization(new B1DetectorConstruction());

  // Physics list
  G4VUserPhysicsList* physicsList = new OpticPhysicsList();
  //physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new B1ActionInitialization());

  // Initialize G4 kernel
  //
  runManager->Initialize();
  runManager->BeamOn(300);
  G4cout<<"total efficieency: "<<runManager->get_total_detetion_eff()<<G4endl;
  //primary_Monte_Carlo->events.size() is always larger by 1 then actual number of primary runs, because 
  //primary_Monte_Carlo->new_sequence() and ->new_event() are called at the end of every sim. (see CustomRunManager->next_event())
  //(so that before every global run last_photon_event is allocated)
  G4cout << "total events proceded: " << (runManager->primary_Monte_Carlo->events.size()-1)+(runManager->extra_run_id) << G4endl;
  runManager->get_detected_spectrum();
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv,"qt");
	//G4UIExecutive* ui = new G4UIExecutive(0, 0);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
