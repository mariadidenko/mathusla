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
/// \file B2bSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "B2bSteppingAction.hh"
#include "B2EventAction.hh"
#include "B2bDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2bSteppingAction::B2bSteppingAction(
                      const B2bDetectorConstruction* detectorConstruction,
                      B2EventAction* eventAction
)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2bSteppingAction::~B2bSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bSteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
    
//    if (volume == fDetConstruction->GetAnglePV()){
    if(step->GetTrack()->GetDefinition()->GetParticleName() == "mu-" ){
    if (volume == fDetConstruction->GetAnglePV()){
        G4ThreeVector ang = step->GetPreStepPoint()->GetMomentumDirection();
        G4ThreeVector *centerVector = new G4ThreeVector(0, 0, 1);
        G4double angle=ang.angle(*centerVector);
        
        angle = angle*180/3.14;
        
        G4cout << ">>> !!!!!!!!!!!!!!!!!!!!! ANGLE: " << angle  << G4endl;
        
        fEventAction->StoreAngle(angle);
        step->GetTrack()->SetTrackStatus(fStopAndKill);
      }
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
