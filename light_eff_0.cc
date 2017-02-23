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

#include <time.h>
#include <windows.h>

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

  // Set mandatory initialization classes:
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
  //runManager->BeamOn(runManager->x_num*runManager->y_num);
  
  time_t timer_start, timer_end;
  time(&timer_start);

#ifdef AR_SPEC_TEST
	  runManager->BeamOn(700000);
#else
#ifdef TOP_MESH_TEST
  runManager->BeamOn(runManager->x_num*runManager->y_num);
#else
  runManager->BeamOn(5000);
#endif
#endif
  //DONE: TODO: check reemiss/noreemiss code
  G4double pb_total, pb_no_reemiss, pb_reemissed;
  //runManager->get_total_detetion_eff(&pb_no_reemiss, &pb_reemissed, &pb_total);
  time(&timer_end);
  G4double seconds = difftime(timer_end, timer_start);
  G4cout<<*(runManager->sim_results);
  //DONE - no, no problem: TODO figure out: messed up with the names for some reason
  G4cout << "total events proceded: " << (runManager->sim_results->Num_of_events())+(runManager->extra_run_id) << G4endl;
  G4cout << "Time elapsed: "<< ((int)seconds/3600)<<"h "<<((int)seconds % 3600)/60<<"m "<<(int)seconds % 60<<"s"<<G4endl;
#ifdef AR_SPEC_TEST
  runManager->get_detected_spectrum();
#endif
#ifdef TOP_MESH_TEST
  runManager->export_to_bmp(&runManager->top_hits_xs, &runManager->top_hits_ys, &runManager->top_hits_probs, "top_mesh_test.bmp");
  runManager->export_to_bmp(&runManager->bot_hits_xs, &runManager->bot_hits_ys, &runManager->bot_hits_probs, "bot_mesh_test.bmp");
#endif
  Beep(1500, 300);
  Sleep(150);
  Beep(1500, 300);
  Sleep(150);
  Beep(1500, 300);
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

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
