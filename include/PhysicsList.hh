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
/// \file electromagnetic/TestEm0/include/PhysicsList.hh
/// \brief Definition of the PhysicsList class
//
//
// $Id: PhysicsList.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// 14.10.02 (V.Ivanchenko) provide modular list on base of old PhysicsList
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PhysicsListMessenger;
class G4VPhysicsConstructor;
class G4Scintillation;				//letting header know these classes exist

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
    virtual ~PhysicsList();
   
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    
    void AddPhysicsList(const G4String& name);
    
    virtual void SetCuts();
    
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);
	void SetCutForProton(G4double);
      
  private:
    G4double fCutForGamma;
    G4double fCutForElectron;
    G4double fCutForPositron;
	G4double fCutForProton;
    G4double fCurrentDefaultCut;
	
	G4VPhysicsConstructor* fRaddecayList;
    G4VPhysicsConstructor* fEmPhysicsList;
	G4VPhysicsConstructor* fParticleList;
	G4VPhysicsConstructor* fDecayPhysicsList;
	G4Scintillation* fScintillationProcess;		// make a real object for physics list, making a pointer!
    G4String                fEmName;

    
    PhysicsListMessenger*   fMessenger;         
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif