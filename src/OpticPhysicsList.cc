#include "OpticsPhysicsList.hh"

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
	AddTransportation();
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
			pmanager->AddDiscreteProcess(wls_abs_process);
			pmanager->AddDiscreteProcess(boundaryProcess);
		}
	}
}

void OpticPhysicsList::SetCuts()
{
	SetCutsWithDefault();
}