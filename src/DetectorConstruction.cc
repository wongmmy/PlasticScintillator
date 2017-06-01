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
// $Id: B1DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProductionCuts.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

//includes for Primitive Scorer and associated values to score
#include "G4PSNofSecondary.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSPassageTrackLength.hh"
#include "G4PSTrackLength.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fRegDetector(0)
{
  G4double cut = 0.5*um;
  fDetectorCuts = new G4ProductionCuts();
  fDetectorCuts->SetProductionCut(cut,"gamma");
  fDetectorCuts->SetProductionCut(cut,"e-");
  fDetectorCuts->SetProductionCut(cut,"e+");
  fDetectorCuts->SetProductionCut(cut,"proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
	delete fDetectorCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	G4GeometryManager::GetInstance()->OpenGeometry();

	//clear old G4Region data (if it exists) and creates new one for new cut values
	if(fRegDetector){delete fRegDetector;}
	fRegDetector = new G4Region("Plastic_Detector_Region");
	fRegDetector->SetProductionCuts(fDetectorCuts);

	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	// Option to switch on/off checking of volumes overlaps
	//
	G4bool checkOverlaps = true;

	//     
	// World
	//
	G4double world_sizeXY = 1.2*m;
	G4double world_sizeZ  = 1.2*m;
	G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

	G4Box* solidWorld =    
		new G4Box("World",                       //its name
		0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

	G4LogicalVolume* logicWorld =                         
		new G4LogicalVolume(solidWorld,          //its solid
		world_mat,           //its material
		"World");            //its name

	G4VPhysicalVolume* physWorld = 
		new G4PVPlacement(0,                     //no rotation
		G4ThreeVector(),       //at (0,0,0)
		logicWorld,            //its logical volume
		"World",               //its name
		0,                     //its mother  volume
		false,                 //no boolean operation
		0,                     //copy number
		checkOverlaps);        //overlaps checking

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// Dimension Definitions ///////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Central Aluminum Casing **ID-Inner Diameter, OD-Outer Diameter**
	G4double AlCasingThick = 1*mm;
	G4double AlCasingOD = 7.6*cm;
	G4double AlCasingID = AlCasingOD-AlCasingThick;

	// Left Lid, seals Mylar to main casing
	G4double LeftLidLeftOD = 7.6*cm;
	G4double LeftLidLeftID = LeftLidLeftOD-AlCasingThick;
	G4double LeftLidRightOD = 7.6*cm;
	G4double LeftLidRightID = 5.535*cm;
	G4double LeftLidThick = 4.17*mm;

	// Right lid for Aluminum housing
	G4double RightLidDiam = AlCasingOD;
	G4double RightLidThick = AlCasingThick;

	// Central Aluminum Casing length subtracting left and right ends
	G4double AlCasingLength = 24*cm - LeftLidThick - RightLidThick;

	// PMT
	G4double PMTWindowDiam = 51*mm;
	G4double PMTWindowThickness = 3*mm;
	G4double PMTWallThickness = 1.7*mm;
	G4double PMTLength = 15*cm;
	G4double PMTDiam = 51*mm;
	G4double PMTRightCapDiam = 51*mm;
	G4double PMTInnerDiam = PMTDiam-PMTWallThickness;

	// Mylar Window
	G4double MylarWindowDiam = 55.3*mm;
	G4double MylarThickness = 8.69*um;
	
	// silicon detector
	G4double SiDiam = 70*mm;
	G4double SiThickness = 100*um;

	// Scintillator, which will be a sensitive volume
	G4double ScintillatorDiam = 5*cm;
	G4double ScintillatorThick = 2*cm;
	
	// Light Guide
	G4double LightGuideDiam = 5*cm;
	G4double LightGuideThick = 2*cm;

	// Aluminum shield for betas
	G4double AluminumShieldDiam = 5 *cm;
	G4double AluminumShieldThick = 0.0001*cm;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// Materials definitions ///////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	// PMT-Borosilicate Glass
	G4Material* Si=nist->FindOrBuildMaterial("G4_Si");
	G4Material* B=nist->FindOrBuildMaterial("G4_B");
	G4Material* Al=nist->FindOrBuildMaterial("G4_Al");
	G4Material* Na=nist->FindOrBuildMaterial("G4_Na");
	G4Material* Mg=nist->FindOrBuildMaterial("G4_Mg");
	G4Material* Ca=nist->FindOrBuildMaterial("G4_Ca");
	G4Material* F=nist->FindOrBuildMaterial("G4_F");
	G4Material* Ti=nist->FindOrBuildMaterial("G4_Ti");
	G4Material* Fe=nist->FindOrBuildMaterial("G4_Fe");
	G4Material* O=nist->FindOrBuildMaterial("G4_O");
	G4String name2 = "Borosilicate Glass";
	G4int ncomponents2 = 10;
	G4double density2=2.60*g/cm3;
	G4Material* Glass_mat= new G4Material(name2,density2,ncomponents2);
	Glass_mat->AddMaterial(Si,0.2577);
	Glass_mat->AddMaterial(B,0.0213);
	Glass_mat->AddMaterial(Al,0.0762);
	Glass_mat->AddMaterial(Na,0.0045);
	Glass_mat->AddMaterial(Mg,0.0254);
	Glass_mat->AddMaterial(Ca,0.1265);
	Glass_mat->AddMaterial(F,0.0021);
	Glass_mat->AddMaterial(Ti,0.0034);
	Glass_mat->AddMaterial(Fe,0.0023);
	Glass_mat->AddMaterial(O,0.4806);

	// Light Guide- Poly(methyl methacrylate), Acrylic
	G4Material* C = nist->FindOrBuildMaterial("G4_C");
	G4Material* H = nist->FindOrBuildMaterial("G4_H");
	G4double AcrylicDensity = 1.18*g/cm3;
	G4int AcrylicNumComponents = 3;
	G4Material* Acrylic_mat = new G4Material("Acrylic",AcrylicDensity,AcrylicNumComponents);
	Acrylic_mat->AddMaterial(C,0.5998);
	Acrylic_mat->AddMaterial(H,0.0805);
	Acrylic_mat->AddMaterial(O,0.3196);

	// Scintillator material, 5.15e22 H per cm3 and 4.68e22 C per cm3, find # mol each, mult. by molecular mass of each element,
	// gives g/cm3 of each element works out to very similar mass fractions of Polyvinyl Toluene (base for plastic scintillator EJ-204)
	G4double ScintillatorDensity = 1.023*g/cm3;
	G4int ScintNumComponents = 2;
	G4Material* Scintillator_mat = new G4Material("PlasticScintillator",ScintillatorDensity,ScintNumComponents);
	Scintillator_mat->AddMaterial(H,0.0843/(0.0843+0.9123));
	Scintillator_mat->AddMaterial(C,0.9123/(0.0843+0.9123));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Optical Properties
	Scintillator_mat->GetIonisation()->SetBirksConstant(5*mm / MeV);
	/*
	const G4int NUMENTRIES = 14;
	G4double Scnt_PP[NUMENTRIES] = { 3.26*eV, 3.22*eV, 3.18*eV, 3.16*eV,
		3.12*eV, 3.10*eV, 3.03*eV, 2.97*eV, 2.93*eV, 2.91*eV, 2.82*eV, 2.74*eV, 2.58*eV, 2.48*eV};

	G4double Scnt_FAST[NUMENTRIES] = { 0.03, 0.2, 0.4, 0.6, 0.8, 0.87, 1.0, 0.8, 0.6, 0.52, 0.4, 0.2, 0.05, 0.02 };



	G4MaterialPropertiesTable* Scnt_MPT = new G4MaterialPropertiesTable();

	Scnt_MPT->AddProperty("FASTCOMPONENT", Scnt_PP, Scnt_FAST, NUMENTRIES);

	Scnt_MPT->AddConstProperty("SCINTILLATIONYIELD", 10400. / MeV);
	Scnt_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);		// 1.0 is a Gaussian
	Scnt_MPT->AddConstProperty("FASTTIMECONSTANT", 0.7*ns);		// using Rise Time for now
	Scnt_MPT->AddConstProperty("YIELDRATIO", 1.0);
	//Scnt_MPT->AddProperty()
	
	

	Scintillator_mat->SetMaterialPropertiesTable(Scnt_MPT);
	*/

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	G4Material* Aluminum = nist->FindOrBuildMaterial("G4_Al");

	G4Material* Mylar = nist->FindOrBuildMaterial("G4_MYLAR");

	G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////// Instantiation of logical and physical volumes ////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Aluminum Casing including left and right lid
	// define overlap thickness
	G4double OLT = 1.0*um;
	G4Cons* LeftLid = new G4Cons("LeftLid",0.5*LeftLidLeftID,0.5*LeftLidLeftOD,0.5*LeftLidRightID,0.5*LeftLidRightOD,
		0.5*LeftLidThick,0.,360*degree);
	
	G4Tubs* CentralAlCasing = new G4Tubs("CentralAlCasing",0.5*AlCasingID,0.5*AlCasingOD,0.5*AlCasingLength,0,360*degree);

	G4UnionSolid* intermediateCasing1 = new G4UnionSolid("intermediateCasing1",CentralAlCasing,LeftLid,0,
		G4ThreeVector(0,0,-0.5*(AlCasingLength+LeftLidThick-OLT)));

	G4Tubs* RightLid = new G4Tubs("RightLid",0.,0.5*RightLidDiam,0.5*RightLidThick,0,360*degree);

	G4UnionSolid* AlCasing = new G4UnionSolid("AlCasing",intermediateCasing1,RightLid,0,G4ThreeVector(0,0,0.5*(RightLidThick+AlCasingLength-OLT)));

	G4LogicalVolume* logicAlCasing = new G4LogicalVolume(AlCasing,Aluminum,"AlCasing");

	G4VPhysicalVolume* physAlCasing = new G4PVPlacement(0,					// Rotation vector
														G4ThreeVector(),					// Translation vector
														logicAlCasing,		// logic volume to be placed
														"physAlCasing",		// name for physical volume
														logicWorld,			// mother volume
														false,				// boolean operations
														0,					// copy number
														checkOverlaps);		// checks for overlaps when run

	// Plastic Scintillator
	G4Tubs* Scintillator = new G4Tubs("Scintillator", 0, 0.5*ScintillatorDiam, 0.5*ScintillatorThick, 0, 360*degree);

	G4LogicalVolume* logicScintillator = new G4LogicalVolume(Scintillator,Scintillator_mat,"Scintillator");
	fScoringVolume = logicScintillator;		//to return the value of the logical volume

	G4VPhysicalVolume* physScintillator = new G4PVPlacement(0,G4ThreeVector(0,0,-0.5*(AlCasingLength-ScintillatorThick)),
															logicScintillator,"physScintillator",logicWorld,false,0,checkOverlaps);

	// Light Guide
	G4Tubs* LightGuide = new G4Tubs("LightGuide",0,0.5*LightGuideDiam,0.5*LightGuideThick,0,360*degree);

	G4LogicalVolume* logicLightGuide = new G4LogicalVolume(LightGuide,Acrylic_mat,"logicLightGuide");

	G4VPhysicalVolume* physLightGuide = new G4PVPlacement(0,G4ThreeVector(0,0,-0.5*(AlCasingLength-LightGuideThick)+ScintillatorThick),
														logicLightGuide,"physLightGuide",logicWorld,false,0,checkOverlaps);

	// PMT
	G4Tubs* PMTOutside = new G4Tubs("PMTOutside",0,0.5*PMTDiam,0.5*PMTLength,0,360*degree);

	G4Tubs* PMTInside = new G4Tubs("PMTInside",0,0.5*PMTInnerDiam,0.5*(PMTLength-2*PMTWindowThickness),0,360*degree);

	G4SubtractionSolid* PMT = new G4SubtractionSolid("PMT",PMTOutside,PMTInside,0,G4ThreeVector());

	G4LogicalVolume* logicPMT = new G4LogicalVolume(PMT,Glass_mat,"logicPMT");

	G4VPhysicalVolume* physPMT = new G4PVPlacement(0,G4ThreeVector(0,0,(-0.5*(AlCasingLength-PMTLength)+ScintillatorThick+LightGuideThick)),
													logicPMT,"physPMT",logicWorld,false,0,checkOverlaps);

	// Mylar Window
	G4Tubs* MylarWindow = new G4Tubs("MylarWindow",0,0.5*MylarWindowDiam,0.5*MylarThickness,0,360*degree);

	G4LogicalVolume* logicMylarWindow = new G4LogicalVolume(MylarWindow,Mylar,"logicMylarWindow");

	G4VPhysicalVolume* physMylarWindow = new G4PVPlacement(0,G4ThreeVector(0,0,-0.5*(AlCasingLength+MylarThickness)),logicMylarWindow,
															"physMylarWindow",logicWorld,false,0,checkOverlaps);

	//ACD silicon detector, temporarily made of air until we want to apply ACD capabilities

	G4Tubs* SiliconDetector = new G4Tubs("SiliconDetector",0,0.5*SiDiam,0.5*SiThickness,0,360*degree);

	G4LogicalVolume* logicSiliconDetector = new G4LogicalVolume(SiliconDetector,world_mat,"logicSiDetector");

	G4VPhysicalVolume* physSiliconDetector = new G4PVPlacement(0,G4ThreeVector(0,0,-0.5*(AlCasingLength+2*LeftLidThick+5*mm)),logicSiliconDetector,
																"physSiliconDetector",logicWorld,false,0,checkOverlaps);

	// Aluminum Shield for the source
	//G4Tubs* AluminumShield = new G4Tubs("AluminumShield", 0, 0.5*AluminumShieldDiam, 0.5*AluminumShieldThick, 0, 360 * degree);

	//G4LogicalVolume* logicAluminumShield = new G4LogicalVolume(AluminumShield, Aluminum, "logicAluminumShield");

	//G4VPhysicalVolume* physAluminumShield = new G4PVPlacement(0, G4ThreeVector(0, 0, -12.9*cm), logicAluminumShield, "physAluminumShield", logicWorld, false, 0, checkOverlaps);

	//Register the sensitive volume, and apply the cuts to it, default cuts(physics list) used everywhere else
	fRegDetector->AddRootLogicalVolume(logicScintillator);
	G4String filterName, particleName;
	G4SDParticleFilter* electronFilter = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
	//////////////////////////////////////////////////////////////////////////////////////////
	//Generate G4MultiFunctionalDetector for plastic scintillator, insert into G4SDManager

	//instantiate G4MultiFunctionalDetector and a string for it's name for later referencing
	G4MultiFunctionalDetector* PlasticScintillator = new G4MultiFunctionalDetector("PlasScin");
	G4SDManager::GetSDMpointer()->AddNewDetector(PlasticScintillator);
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	logicScintillator->SetSensitiveDetector(PlasticScintillator);
	
	//Instantiate a specific G4PrimitiveScorer and give name for later reference
	G4PSNofSecondary* nElectron = new G4PSNofSecondary("nElectron");
	nElectron->SetFilter(electronFilter);
	PlasticScintillator->RegisterPrimitive(nElectron);

	G4PSEnergyDeposit* eDep = new G4PSEnergyDeposit("eDep");
	PlasticScintillator->RegisterPrimitive(eDep);

	G4PSDoseDeposit* doseDep = new G4PSDoseDeposit("doseDep");
	PlasticScintillator->RegisterPrimitive(doseDep);

	G4PSPassageTrackLength* trackLength = new G4PSPassageTrackLength("trackLength");
	PlasticScintillator->RegisterPrimitive(trackLength);

	G4MultiFunctionalDetector* ACDScorer = new G4MultiFunctionalDetector("ACD");
	G4SDManager::GetSDMpointer()->AddNewDetector(ACDScorer);
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	logicSiliconDetector->SetSensitiveDetector(ACDScorer);

	G4PSEnergyDeposit* eDep_ACD = new G4PSEnergyDeposit("eDep_ACD");
	ACDScorer->RegisterPrimitive(eDep_ACD);


	//always return the physical World
	//
	return physWorld;
}
