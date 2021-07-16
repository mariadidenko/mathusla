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
/// \file B2ActionInitialization.cc
/// \brief Implementation of the B2ActionInitialization class

#include "B2ActionInitialization.hh"
#include "B2PrimaryGeneratorAction.hh"
#include "B2RunAction.hh"
#include "B2bSteppingAction.hh"
#include "B2EventAction.hh"
#include "B2bDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


B2bActionInitialization::B2bActionInitialization
                            (B2bDetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
  fDetConstruction(detConstruction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2bActionInitialization::~B2bActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bActionInitialization::BuildForMaster() const
{
    //SetUserAction(new B4RunAction);
    B2EventAction *eventAction = new B2EventAction();
    SetUserAction(new B2RunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bActionInitialization::Build() const
{
    auto eventAction = new B2EventAction;
    B2RunAction *runAction = new B2RunAction(eventAction);
  
    SetUserAction(new B2PrimaryGeneratorAction);
    SetUserAction(runAction);
  
  SetUserAction(new B2EventAction);
  SetUserAction(new B2bSteppingAction(fDetConstruction,eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
