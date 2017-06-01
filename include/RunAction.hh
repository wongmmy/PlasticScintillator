#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class PrimaryGeneratorAction;

class G4Run;
class G4Event;

class RunAction : public G4UserRunAction

{
public:
	RunAction(PrimaryGeneratorAction*);
	virtual ~RunAction();

	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);

	void RecordData(G4int evtNb, const G4Event*, G4int nElectron, G4double eDep, G4double doseDep, G4double trackLength, double qEdep); 

	
private:
	PrimaryGeneratorAction* particleGun;

	G4String outputFile;
	G4String outputFile_INFO;
	FILE* pFile;
	FILE* pFile_INFO;

};

#endif