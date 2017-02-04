#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "CustomRunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="opticalphoton");
  fParticleGun->SetParticleDefinition(particle);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.1,1.));
  //fParticleGun->SetParticleEnergy(1.*eV);
}

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  CustomRunManager* manman = (CustomRunManager*)G4RunManager::GetRunManager();
  fParticleGun->SetParticlePosition(manman->FetchPosition());
  fParticleGun->SetParticleMomentumDirection(manman->FetchMomentum());
  fParticleGun->SetParticlePolarization(manman->FetchPolarization());
  fParticleGun->SetParticleEnergy(manman->FetchEnergy());
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
