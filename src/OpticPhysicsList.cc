#include "OpticsPhysicsList.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4OpticalPhoton.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "CustomTransportation.hh"
#include "CustomWLSProcess.hh"
#include "CustomOpBoundaryProcess.hh"
#include "CustomMapProcess.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

OpticPhysicsList::OpticPhysicsList() : G4VUserPhysicsList()
{}

OpticPhysicsList::~OpticPhysicsList()
{}

void OpticPhysicsList::ConstructParticle()
{
	G4OpticalPhoton::OpticalPhotonDefinition();
	G4Geantino::GeantinoDefinition();
}

void OpticPhysicsList::ConstructProcess()
{
	AddCustomTransportation();
	CustomWLSProcess* wls_abs_process = new CustomWLSProcess();
	CustomOpBoundaryProcess* boundaryProcess = new CustomOpBoundaryProcess();
	CustomMapProcess* mapProcess = new CustomMapProcess();
#ifdef DEBUG_MC_NODES
	boundaryProcess->SetVerboseLevel(2);
#else
	boundaryProcess->SetVerboseLevel(0);
#endif
	theParticleIterator->reset();
	while ((*theParticleIterator)())
	{
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		pmanager->AddProcess(mapProcess,InActivated,ordDefault,ordDefault); //added for every particle
		if (particleName == "opticalphoton") 
		{
			G4cout << " AddDiscreteProcesses to OpticalPhoton " << G4endl;
#if !defined(NO_WLS_PROC)
			pmanager->AddDiscreteProcess(wls_abs_process);
#endif
			pmanager->AddDiscreteProcess(boundaryProcess);
		}
	}
}

void OpticPhysicsList::SetCuts()
{
	SetCutsWithDefault();
}

//TODO: Does no account for parallel world
void OpticPhysicsList::AddCustomTransportation(void)
{
	//((this->subInstanceManager.offset[this->g4vuplInstanceID])._thePLHelper)->;
	G4int verboseLevelTransport = 0;

#ifdef G4VERBOSE
	if (verboseLevel >2){
		G4cout << "G4PhysicsListHelper::AddCustomTransportation()  " << G4endl;
	}
#endif
	G4VProcess *theTransportationProcess = new CustomTransportation(verboseLevelTransport);
	// loop over all particles in G4ParticleTable
	theParticleIterator->reset();
	while ((*theParticleIterator)()){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		// Add transportation process for all particles 
		if (pmanager == 0) {
			// Error !! no process manager
#ifdef G4VERBOSE
			if (verboseLevel>0){
				G4cout << "G4PhysicsListHelper::AddTransportation  "
					<< " : No Process Manager for "
					<< particle->GetParticleName() << G4endl;
			}
#endif
			G4Exception("G4PhysicsListHelper::AddTransportation",
				"Run0104", FatalException,
				"No process manager");
			continue;
		}
		// Molecule use different type transportation
		if (particle->GetParticleType() == "Molecule") continue;

		// add transportation with ordering = ( -1, "first", "first" )
		pmanager->AddProcess(theTransportationProcess);
		pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
		pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
	}
}