#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

SteppingAction::SteppingAction(EventAction* eventAction) : G4UserSteppingAction(),		//constructor
sensV(0)
//in the name space of SteppingAction, the constructor SteppingAction takes an EventAction pointer names eventAction, extends G4UserSteppingAction

{
	scEdep = 0.;	//scoring variables
	scStepLength = 0.;
	correction = 0.;
	qEdep = 0.;
	kb = 1;
}

SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
	if (!sensV) //just to fill the sensitive volume with the logical volume
	{
		const DetectorConstruction* DC = static_cast<const DetectorConstruction* > (G4RunManager::GetRunManager()->GetUserDetectorConstruction());		
		//to get current run manager, saving DC to access pointer of DetectorConstruction to get to the function that gets LogicalVolume
		sensV = DC->GetScoringVolume();
	}
	if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == sensV)
	{
		scEdep = step->GetDeltaEnergy();
		scStepLength = step->GetStepLength();
		correction = (scEdep / scStepLength) / (1 + kb*(scEdep / scStepLength));
		qEdep += scEdep*correction;
	}
}