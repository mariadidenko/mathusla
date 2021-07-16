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
/// \file B2bDetectorConstruction.cc
/// \brief Implementation of the B2bDetectorConstruction class
 
#include "B2bDetectorConstruction.hh"
#include "B2bDetectorMessenger.hh"
#include "B2bChamberParameterisation.hh"
#include "B2TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B2bDetectorConstruction::fMagFieldMessenger = 0;
 
B2bDetectorConstruction::B2bDetectorConstruction()
:G4VUserDetectorConstruction(),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fanglePV(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new B2bDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2bDetectorConstruction::~B2bDetectorConstruction()
{
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* B2bDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_Galactic");
  //  defined using NIST Manager
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
   // fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_W");
  // fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_AIR");

  // Wolframm defined using NIST Manager
//  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_W");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2bDetectorConstruction::DefineVolumes()
{
 
  auto env_mat = G4Material::GetMaterial("G4_Galactic");

  // Sizes of the principal geometrical components (solids)
  
  G4int NbOfChambers = 1;
  
    // Envelope parameters
    //
    G4double env_sizeXY = 10*m, env_sizeZ = 10*m;
    
    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;
    
  // World
    G4double world_sizeXY = 1.2*env_sizeXY;
    G4double world_sizeZ  = 1.2*env_sizeZ;
    
    G4double detSizeXY =  9.*m;
    G4double detSizeZ =  2.*cm;

    G4Material* world_mat  = G4Material::GetMaterial("G4_AIR");
    
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
   
    //
    // Envelope / traker
    //
    G4Box* solidEnv =
      new G4Box("Envelope",                    //its name
          0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
        
    G4LogicalVolume* logicEnv =
      new G4LogicalVolume(solidEnv,            //its solid
                          world_mat,             //its material
                          "Envelope");         //its name
                 
    new G4PVPlacement(0,                       //no rotation
                      G4ThreeVector(),         //at (0,0,0)
                      logicEnv,                //its logical volume
                      "Envelope",              //its name
                      logicWorld,              //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    G4ThreeVector pos2 = G4ThreeVector(0, 0, - 1*m);
    G4ThreeVector pos1 = G4ThreeVector(0, 0, 0);
    
    G4double targetSizeXY =  9.*m;
    G4double targetSizeZ =  2.*cm;
    
    // Module:
    //
    auto moduleS
      = new G4Box("Module",     // its name
                  0.5*targetSizeXY, 0.5*targetSizeXY,0.5*targetSizeZ); // its size
                           
    moduleLV
      = new G4LogicalVolume(
                   moduleS,     // its solid
                   fTargetMaterial,  // its material
                   "Module");   // its name
    
    fanglePV =
        new G4PVPlacement(
                     0,                // no rotation
                     pos1,//G4ThreeVector(),  // at (0,0,0)
                     moduleLV,          // its logical volume
                     "Module",    // its name
                     logicEnv,          // its mother  volume
                     false,            // no boolean operation
                     0,                // copy number
                     fCheckOverlaps);  // checking overlaps
    
    //
    // Chamber detector
    //
    //G4Material* det_mat = nist->FindOrBuildMaterial("G4_W");
   
    G4Box* detS =
      new G4Box("Detector",0.5*detSizeXY, 0.5*detSizeXY, 0.5*detSizeZ);
                        
    fLogicChamber =
      new G4LogicalVolume(detS,         //its solid
                          fChamberMaterial,          //its material
                          "Detector");           //its name
                 
//    new G4PVPlacement(0,                       //no rotation
//                      pos1,                    //at position
//                      fLogicChamber,             //its logical volume
//                      "Detector",                //its name
//                      logicEnv,                //its mother  volume
//                      false,                   //no boolean operation
//                      0,                       //copy number
//                      checkOverlaps);          //overlaps checking

    //
    // Target
    //
      
      G4Box* targetS =
        new G4Box("Target",0.5*targetSizeXY, 0.5*targetSizeXY, 0.5*targetSizeZ);
                  
      fLogicTarget =
       new G4LogicalVolume(targetS,         //its solid
                          fTargetMaterial,          //its material
                          "Target");           //its name
                 
     new G4PVPlacement(0,                       //no rotation
                      pos2,                    //at position
                      fLogicTarget,             //its logical volume
                      "Target",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld  ->SetVisAttributes(boxVisAtt);
  fLogicTarget ->SetVisAttributes(boxVisAtt);
  logicEnv ->SetVisAttributes(boxVisAtt);
  G4VisAttributes* mod= new G4VisAttributes(G4Colour(1.0,0,0));
  moduleLV -> SetVisAttributes(mod);

  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  fLogicChamber->SetVisAttributes(chamberVisAtt);
  
  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

//  G4double maxStep = 0.5*detSizeZ ;
//  fStepLimit = new G4UserLimits(maxStep);
//  logicEnv->SetUserLimits(fStepLimit);
  
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2bDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector(fLogicChamber,  aTrackerSD );//important !fLogicChamber!

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2bDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        if (fLogicChamber) fLogicChamber->SetMaterial(fChamberMaterial);
        G4cout
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout
          << G4endl
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2bDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
