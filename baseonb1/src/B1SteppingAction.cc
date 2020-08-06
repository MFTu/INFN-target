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
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4Track.hh"

#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"


#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
// G4cout << G4BestUnit(StepSize, "Length");



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }
 

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // get kinetic energy
  G4double kineticEfromStepping 
    = step->GetPreStepPoint()->GetKineticEnergy();	


  // get particle information
  G4Track* track = step->GetTrack();
  
  // From the track you can obtain the pointer to the dynamic particle:
  const G4DynamicParticle* dynParticle = track->GetDynamicParticle();

  // From the dynamic particle, retrieve the particle definition:
  G4ParticleDefinition* particle = dynParticle->GetDefinition();

  // The dynamic particle class contains e.g. the kinetic energy after the step:
  G4double kinEnergyfromTrack = dynParticle->GetKineticEnergy();

  G4ThreeVector momentumfromTrack = dynParticle->GetMomentum();
  G4double  momentumX = momentumfromTrack.getX();
  G4double  momentumY = momentumfromTrack.getY();
  G4double  momentumZ = momentumfromTrack.getZ();



  const G4ThreeVector& momentumdirection = dynParticle->GetMomentumDirection(); 
  G4double  momentumdirectionX = momentumdirection.getX();
  G4double  momentumdirectionY = momentumdirection.getY();
  G4double  momentumdirectionZ = momentumdirection.getZ();



  // Kinetic energy at the start point (vertex position) of the track
  const G4double vertexkinEnergy = track->GetVertexKineticEnergy();   

  // From the particle definition class you can retrieve static information, like the particle name:
  G4String particleName = particle->GetParticleName();


  G4int parentID = track->GetParentID();
  G4int trackID = track->GetTrackID();

  const G4ThreeVector& VertexPosition = track->GetVertexPosition();
  G4double VertexPositionZ = VertexPosition.getZ();

  const G4VProcess* creatorProcess = track->GetCreatorProcess(); 
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);  


  //G4cout << "From_steppingaction" << ' ' << particleName << ' ' << kineticEfromStepping/(CLHEP::MeV) << ' ' << "MeV" << ' ' << kinEnergyfromTrack/(CLHEP::MeV) << ' ' << "MeV" << G4endl;


  G4cout << "From_steppingaction" << ' ' << particleName << ' ' << kineticEfromStepping/(MeV) << ' ' << "MeV" << ' ' << momentumdirectionX << ' ' << momentumdirectionY << ' ' << momentumdirectionZ << ' ' << VertexPositionZ/(cm) << ' ' << "cm" << ' ' << parentID << ' ' << G4endl;


  //G4cout << " This came from steppingaction" << G4endl;


/*
  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();   
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;


  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);  
*/


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

