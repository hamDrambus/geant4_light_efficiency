#ifndef OpticPhysicsList_h
#define OpticPhysicsList_h
#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class OpticPhysicsList : public G4VUserPhysicsList
{
public:
	OpticPhysicsList();
	virtual ~OpticPhysicsList();
	virtual void ConstructParticle();
	virtual void ConstructProcess();
	virtual void SetCuts();
};

#endif