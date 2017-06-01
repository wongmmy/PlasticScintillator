#include "EventAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"

#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction):G4UserEventAction(),
	fRunAct(runAction),fPrintModulo(10000)
{
	//steppingAction = static_cast<const SteppingAction* > (G4RunManager::GetRunManager()->Get);
}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* evt)
{	//returns event number
	G4int evtNb = evt->GetEventID();
	steppingAction->ResetqEdep();
	if (evtNb==0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();

		
		//ID is used to find the value within G4HCofThisEvent (both share same ID number for a value)
		nElectronID = SDMan->GetCollectionID("PlasScin/nElectron");
		eDepID = SDMan->GetCollectionID("PlasScin/eDep");
		doseDepID = SDMan->GetCollectionID("PlasScin/doseDep");
		trackLengthID = SDMan->GetCollectionID("PlasScin/trackLength");
		//add particle name call per event maybe
		eDepID_ACD = SDMan->GetCollectionID("ACD/eDep_ACD");
	}
	if(evtNb%fPrintModulo == 0) {
		//prints to screen on modulo number for progress of run
		G4cout << "\n------> Begin of event number: " << evtNb << G4endl;
	}

}

void EventAction::EndOfEventAction(const G4Event* evt)
{
	G4int evtNb = evt->GetEventID();

	//Hits collection
	//HCofThisEvent automatically stores information that is recorded 
	//in the sensitive detector definitions inside DetectorConstruction



	G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
	if (!HCE) return;

	eventNElectron = (G4THitsMap<G4double>*)(HCE->GetHC(nElectronID));
	eventEDep = (G4THitsMap<G4double>*)(HCE->GetHC(eDepID));
	eventDoseDep = (G4THitsMap<G4double>*)(HCE->GetHC(doseDepID));
	eventTrackLength = (G4THitsMap<G4double>*)(HCE->GetHC(trackLengthID));
	eventEDep_ACD = (G4THitsMap<G4double>*)(HCE->GetHC(eDepID_ACD));

	//zero out variables after previous run
	G4int nElectrons = 0;
	G4double eDep = 0.;
	G4double doseDep = 0.;
	G4double trackLength = 0.;
	G4double eDep_ACD = 0.;


	//generates map, itr, that is an iterator filled column 1=int (marker)
	//and column two is the double values used by simulation vars (nElectron, etc) all values are pointer
	std::map<G4int, G4double*>::iterator itr;

	//loops through itr=beginning of eventNElectron THitsMap to itr=end of eventNElectron THitsMap
	for (itr = eventNElectron->GetMap()->begin(); itr != eventNElectron->GetMap()->end(); itr++) {
		nElectrons = *(itr->second);
	}

	// Get the total energy deposited in this event
	for (itr = eventEDep->GetMap()->begin(); itr != eventEDep->GetMap()->end(); itr++) {
		eDep = *(itr->second);
	}

	// Get the dose deposited in this event
	for (itr = eventDoseDep->GetMap()->begin(); itr != eventDoseDep->GetMap()->end(); itr++) {
		doseDep = *(itr->second);
	}

	// Get the track length for this event
	for (itr = eventTrackLength->GetMap()->begin(); itr != eventTrackLength->GetMap()->end(); itr++) {
		trackLength = *(itr->second);
	}

	// Get the total energy deposited in this event in the ACD
	for (itr = eventEDep_ACD->GetMap()->begin(); itr != eventEDep_ACD->GetMap()->end(); itr++) {
		eDep_ACD = *(itr->second);
	}
	
	double qEdep = steppingAction->GetqEdep();

	if (eDep > 0){
		fRunAct->RecordData(evtNb, evt, nElectrons, eDep, doseDep, trackLength,qEdep);
	}
}
