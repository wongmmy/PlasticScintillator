#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4THitsMap.hh"

class RunAction;
class SteppingAction;
class EventAction : public G4UserEventAction
{
	public:
		//constructor/destruction
		EventAction(RunAction* runAction);
		virtual ~EventAction();

		virtual void BeginOfEventAction (const G4Event*);
		virtual void EndOfEventAction (const G4Event*);
		void SetSteppingAction(SteppingAction* sa) { steppingAction = sa; }

	private:

		RunAction* fRunAct;
		SteppingAction* steppingAction;
		G4int nElectronID;
		G4int eDepID;
		G4int doseDepID;
		G4int trackLengthID;
		G4int eDepID_ACD;

		G4THitsMap<G4double>* eventNElectron;
		G4THitsMap<G4double>* eventEDep;
		G4THitsMap<G4double>* eventDoseDep;
		G4THitsMap<G4double>* eventTrackLength;
		G4THitsMap<G4double>* eventEDep_ACD;

		G4int fPrintModulo; //print to screen progress every n events run
};

#endif



