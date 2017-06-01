//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//// BetaDetector is a plastic scintillator run in coincidence for gathering beta spectra produced 
// by geant4's G4RadioactiveDecay process
// Richard Garnett



#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4PhysListFactory.hh"

// Simulation Files
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "PhysicsList.hh"
#include "SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  	// Choose the Random engine
  	G4Random::setTheEngine(new CLHEP::RanecuEngine);
     
  	// Construct the default run manager
	#ifdef G4MULTITHREADED
  		G4MTRunManager* runManager = new G4MTRunManager;
	#else
  		G4RunManager* runManager = new G4RunManager;
	#endif  

  	runManager->SetUserInitialization(new PhysicsList);

  	// Set mandatory initialization classes
  	DetectorConstruction* detector = new DetectorConstruction();
  	runManager->SetUserInitialization(detector);
  	
  	// Primary Generator Action
	PrimaryGeneratorAction* primary = new PrimaryGeneratorAction();
	runManager->SetUserAction(primary);
	
	// Run Action
	RunAction* runAction = new RunAction(primary);
	runManager->SetUserAction(runAction);
	
	// Event Action
	EventAction* eventAction = new EventAction(runAction);  //giving EventAction a name because we are using it twice, SteppingAction needs to know of EventAction
	runManager->SetUserAction(eventAction);

	SteppingAction* steppingAction = new SteppingAction(eventAction); //steppingAction is not a const here and can be changed
	runManager->SetUserAction(steppingAction);
	eventAction->SetSteppingAction(steppingAction); //accessing steppingAction through eventAction now to bypass the const problem
	//runManager->GetCurrentEvent()->GetEventID();
  	// Initialize G4 kernel
  	runManager->Initialize();

	
  
	#ifdef G4VIS_USE
  	// Initialize visualization
  	G4VisManager* visManager = new G4VisExecutive;
  	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  	// G4VisManager* visManager = new G4VisExecutive("Quiet");
  	visManager->Initialize();
	#endif

  	// Get the pointer to the User Interface manager
  	G4UImanager* UImanager = G4UImanager::GetUIpointer();

  	if (argc!=1)   // batch mode
    {
    	G4String command = "/control/execute ";
      	G4String fileName = argv[1];
      	UImanager->ApplyCommand(command+fileName);
    	}
  	else
    	{  // interactive mode : define UI session
		#ifdef G4UI_USE
      		G4UIExecutive* ui = new G4UIExecutive(argc, argv, "csh");
		#ifdef G4VIS_USE
      		UImanager->ApplyCommand("/control/execute init_vis.mac"); 
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


