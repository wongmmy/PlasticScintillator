#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

#include "PrimaryGeneratorAction.hh"

#include <stdio.h>
#include <time.h>

#include "G4SystemOfUnits.hh"

#define FILENAME_FORMAT "output_%Y%m%d_%H%M%S.csv"
#define FILENAME_FORMAT_INFO "output_%Y%m%d_%H%M%S_INFO.info"
#define FILENAME_SIZE 60

RunAction::RunAction(PrimaryGeneratorAction* primary): particleGun(primary) {}

RunAction::~RunAction(){}

void RunAction::BeginOfRunAction(const G4Run* run)
{
	//set seed
	long seeds[2];
	time_t systime = time(NULL);
	seeds[0] = (long) systime;
	seeds[1] = (long) (systime*G4UniformRand());
	//CLHEP::HepRandom::setTheSeeds(seeds);
	seeds[0] = 1;
	seeds[1] = 2;
	CLHEP::HepRandom::setTheSeeds(seeds);
	
	// get the time for the file names
	time_t now = time(0);

	static char outputFileTemplate[FILENAME_SIZE];
	static char outputFileTemplate_INFO[FILENAME_SIZE];
    
	strftime(outputFileTemplate, sizeof(outputFileTemplate),FILENAME_FORMAT, localtime(&now));
	strftime(outputFileTemplate_INFO, sizeof(outputFileTemplate_INFO), FILENAME_FORMAT_INFO, localtime(&now));
 
	outputFile = outputFileTemplate;
	outputFile_INFO = outputFileTemplate_INFO;

	// open anti-coinc file and coinc file
	pFile = fopen (outputFile,"w");
	//outFile_ACOINC is set by default to use std::ios::out and write in location outpuFile_ACOINC)
	std::ofstream outFile(outputFile);

	outFile << "Event ID" << ","
			<< "Particle Type" << ","
			<< "TEPC Primary Electrons" << ","
			<< "TEPC Energy Deposited [keV]" << ","
			<< "TEPC Dose Deposited [pGy]" << ","
			<< "TEPC Track Length [mm]" << ","
		    << "TEPC Quenched Energy Deposited [keV]"
		    
			<< G4endl;

}

void RunAction::RecordData(G4int evtNb, const G4Event* evt, G4int nElectron, G4double eDep, G4double doseDep, G4double trackLength,double qEdep)
{	//takes event and returns the primary particle of that event (hopefully what ever was emitted in RDecay)
	const G4ParticleDefinition* part=(evt->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition());
	//generates new line in outputfile, either as beginning of file(default) or appended (ios::app)
	std::ofstream outFile(outputFile, std::ios::out|std::ios::app);
	outFile << evtNb << "," << part->GetParticleName() << "," << nElectron << "," << eDep/keV << "," << doseDep/(1.e-12*gray) << "," << trackLength/mm << "," << qEdep/keV << G4endl;
} 



void RunAction::EndOfRunAction(const G4Run* run)
{
	//closes acoinc and coinc files, done writing into them
	fclose(pFile);

	//creates the info file, using date/time information
	pFile_INFO = fopen (outputFile_INFO,"w");
	std::ofstream outFile_INFO(outputFile_INFO);

	// prints source info into the info file using GPS
	outFile_INFO << "========================	Source Information ===================" << G4endl;
	outFile_INFO << "Type:" << "\t" << particleGun->GetGPS()->GetParticleDefinition() << G4endl;
	outFile_INFO << "Energy: " << "\t" << G4BestUnit(particleGun->GetGPS()->GetParticleEnergy(), "Energy") << G4endl;
	outFile_INFO << "Number: " << "\t" << run->GetNumberOfEvent() << G4endl;
	outFile_INFO << "=================================================================" << G4endl;

	//lets me know run has ended
	G4cout << "\n--------------------End of run--------------------------------------\n" << G4endl;
	//closes info file
	fclose(pFile_INFO);
}
