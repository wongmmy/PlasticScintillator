#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

	public:
  		// Constructor
  		PrimaryGeneratorAction();    

  		// Destructor
  		virtual ~PrimaryGeneratorAction();
  
  		// Method
  		void GeneratePrimaries(G4Event*);
  
	public:
    	G4GeneralParticleSource* GetGPS() {return particleGun;};

	private:
		// Data member
 		G4GeneralParticleSource* particleGun;	 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif