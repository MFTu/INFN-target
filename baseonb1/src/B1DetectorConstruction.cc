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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4Element.hh"

//#include "globals.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include <cmath>

//#include "globals.hh"
//#include "G4ThreeVector.hh"
//#include <CLHEP/Vector/Rotation.h>

// Eye candy headers
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  // G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  // G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     /run/beamOn 10
  // World
  //
  G4double world_sizeXY = 80.*cm;
  G4double world_sizeZ  = 140.*cm;

  

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        vacuum,           //its material
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
  // Envelope
  //  
  // G4Box* solidEnv =    
  //   new G4Box("Envelope",                    //its name
  //       0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  // G4LogicalVolume* logicEnv =                         
  //   new G4LogicalVolume(solidEnv,            //its solid
  //                       env_mat,             //its material
  //                       "Envelope");         //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   G4ThreeVector(),         //at (0,0,0)
  //                   logicEnv,                //its logical volume
  //                   "Envelope",              //its name
  //                   logicWorld,              //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking



  //     
  // First LiF (0.5 cm)
  //  

  // build material LiF

  G4int LiFncomponents, LiFnatoms;
  G4double LiFdensity = 2.64*g/cm3;
  G4Element* elemLi = nist->FindOrBuildElement("Li");
  G4Element* elemF = nist->FindOrBuildElement("F");

  G4Material* LiF_material = new G4Material("materialLiF", LiFdensity, LiFncomponents = 2);
  LiF_material -> AddElement(elemLi, LiFnatoms = 1);
  LiF_material -> AddElement(elemF, LiFnatoms = 1);

  G4ThreeVector pos_LiF1st = G4ThreeVector(0, 0, (0.5/2)*cm);

  // Conical section for first LiF       
  G4double LiFfirst_rminb =  0.*cm, LiFfirst_rmaxb = (6+(38.5/3))*cm;
  G4double LiFfirst_rmina =  0.*cm, LiFfirst_rmaxa = (6+13.0)*cm;
  G4double LiFfirst_hz = (0.5/2)*cm;
  G4double LiFfirst_phimin = 0.*deg, LiFfirst_phimax = 360.*deg;
  G4Cons* solidLiFfirst =    
    new G4Cons("LiFfirst", 
    LiFfirst_rmina, LiFfirst_rmaxa, LiFfirst_rminb, LiFfirst_rmaxb, LiFfirst_hz,
    LiFfirst_phimin, LiFfirst_phimax);
                      
  G4LogicalVolume* logicLiFfirst =                         
    new G4LogicalVolume(solidLiFfirst,         //its solid
                        LiF_material,          //its material
                        "LiFfirst");           //its name
  logicLiFfirst->SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
               
//  new G4PVPlacement(0,                       //no rotation
//                    pos_LiF1st,                    //at position
//                    logicLiFfirst,             //its logical volume
//                    "LiFfirst",                //its name
//                    logicWorld,                //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking





  //     
  // AlF3 (36.5 cm)
  //  

  // build material AlF3

  G4int AlF3ncomponents, AlF3natoms;
  G4double AlF3density = 2.88*g/cm3;
  G4Element* elemAl = nist->FindOrBuildElement("Al");
  //G4Element* elemF = nist->FindOrBuildElement("F");

  G4Material* AlF3_material = new G4Material("materialAlF3", AlF3density, AlF3ncomponents = 2);
  AlF3_material -> AddElement(elemAl, AlF3natoms = 1);
  AlF3_material -> AddElement(elemF, AlF3natoms = 3);

  G4ThreeVector pos_AlF3 = G4ThreeVector(0, 0, (0.5+(36.5/2))*cm);

  // Conical section for AlF3       
  G4double AlF3_rmina =  0.*cm, AlF3_rmaxa = (6+(38.5/3))*cm;
  G4double AlF3_rminb =  0.*cm, AlF3_rmaxb = (6+(2.0/3))*cm;
  G4double AlF3_hz = (36.5/2)*cm;
  G4double AlF3_phimin = 0.*deg, AlF3_phimax = 360.*deg;
  G4Cons* solidAlF3 =    
    new G4Cons("AlF3install", 
    AlF3_rmina, AlF3_rmaxa, AlF3_rminb, AlF3_rmaxb, AlF3_hz,
    AlF3_phimin, AlF3_phimax);
                      
  G4LogicalVolume* logicAlF3 =                         
    new G4LogicalVolume(solidAlF3,         //its solid
                        AlF3_material,          //its material
                        "AlF3install");           //its name
  logicAlF3->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
               
//  new G4PVPlacement(0,                       //no rotation
//                    pos_AlF3,                    //at position
//                    logicAlF3,             //its logical volume
//                    "AlF3install",                //its name
//                    logicWorld,                //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking






  //     
  // Second LiF (0.5 cm)
  //  

  G4ThreeVector pos_LiF2nd = G4ThreeVector(0, 0, (0.5+36.5+(0.5/2))*cm);

  // Conical section for second LiF      
  G4double LiFsecond_rminb =  0.*cm, LiFsecond_rmaxb = (6+(1.5/3))*cm;
  G4double LiFsecond_rmina =  0.*cm, LiFsecond_rmaxa = (6+(2.0/3))*cm;
  G4double LiFsecond_hz = (0.5/2)*cm;
  G4double LiFsecond_phimin = 0.*deg, LiFsecond_phimax = 360.*deg;
  G4Cons* solidLiFsecond =    
    new G4Cons("LiFsecond", 
    LiFsecond_rmina, LiFsecond_rmaxa, LiFsecond_rminb, LiFsecond_rmaxb, LiFsecond_hz,
    LiFsecond_phimin, LiFsecond_phimax);
                      
  G4LogicalVolume* logicLiFsecond =                         
    new G4LogicalVolume(solidLiFsecond,         //its solid
                        LiF_material,          //its material
                        "LiFsecond");           //its name
  logicLiFsecond->SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
               
//  new G4PVPlacement(0,                       //no rotation
//                    pos_LiF2nd,                    //at position
//                    logicLiFsecond,             //its logical volume
//                    "LiFsecond",                //its name
//                    logicWorld,                //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking





  //
  // Ti install (1 cm)
  //  
  G4Material* Ti_mat = nist->FindOrBuildMaterial("G4_Ti"); //Ti
  G4ThreeVector pos_Ti = G4ThreeVector(0, 0, (0.5+36.5+0.5+(1.0/2))*cm);
        
  // Conical section for Ti      
  G4double Ti_rmina =  0.*cm, Ti_rmaxa = (6+(1.5/3))*cm;
  G4double Ti_rminb =  0.*cm, Ti_rmaxb = (6+(0.5/3))*cm;
  G4double Ti_hz = (1.0/2)*cm;
  G4double Ti_phimin = 0.*deg, Ti_phimax = 360.*deg;
  G4Cons* solidTi =    
    new G4Cons("TiInstall", 
    Ti_rmina, Ti_rmaxa, Ti_rminb, Ti_rmaxb, Ti_hz,
    Ti_phimin, Ti_phimax);
                      
  G4LogicalVolume* logicTi =                         
    new G4LogicalVolume(solidTi,         //its solid
                        Ti_mat,          //its material
                        "TiInstall");           //its name
  //logicTi->SetVisAttributes(new G4VisAttributes(G4Colour(255/255.,183/255.,169/255.)));
  logicTi->SetVisAttributes(new G4VisAttributes(G4Colour::Grey()));

//  new G4PVPlacement(0,                       //no rotation
//                    pos_Ti,                    //at position
//                    logicTi,             //its logical volume
//                    "TiInstall",                //its name
//                    logicWorld,                //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking




  //
  // Gamma_Ray_Sheilding (0.5cm)
  //  
  G4Material* grshielding_mat = nist->FindOrBuildMaterial("G4_Bi"); //Bi
  G4ThreeVector pos_grshielding = G4ThreeVector(0, 0, (0.5+36.5+0.5+1.0+(0.5/2))*cm);
        
  // Conical section for Gamma_Ray_Sheilding       
  G4double grshielding_rmina =  0.*cm, grshielding_rmaxa = (6+(0.5/3))*cm;
  G4double grshielding_rminb =  0.*cm, grshielding_rmaxb = 6.0*cm;
  G4double grshielding_hz = (0.5/2)*cm;
  G4double grshielding_phimin = 0.*deg, grshielding_phimax = 360.*deg;
  G4Cons* solidGrshielding =    
    new G4Cons("Grshielding", 
    grshielding_rmina, grshielding_rmaxa, grshielding_rminb, grshielding_rmaxb, grshielding_hz,
    grshielding_phimin, grshielding_phimax);
                      
  G4LogicalVolume* logicGrshielding =                         
    new G4LogicalVolume(solidGrshielding,         //its solid
                        grshielding_mat,          //its material
                        "Grshielding");           //its name
  logicGrshielding->SetVisAttributes(new G4VisAttributes(G4Colour(255/255.,183/255.,169/255.)));

//  new G4PVPlacement(0,                       //no rotation
//                    pos_grshielding,                    //at position
//                    logicGrshielding,             //its logical volume
//                    "Grshielding",                //its name
//                    logicWorld,                //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking



  //
  // full cone of assembly
  //  

  //G4ThreeVector pos_fullcone_assembly = G4ThreeVector(0, 0, ((0.5+36.5+0.5+1.0+0.5)/2)*cm);
             
  G4double fullconemp_rmina =  0.*cm, fullconemp_rmaxa = (6.0+13.0)*cm;
  G4double fullconemp_rminb =  0.*cm, fullconemp_rmaxb = (6.0)*cm;
  G4double fullconemp_hz = ((0.5+36.5+0.5+1.0+0.5)/2)*cm;
  G4double fullconemp_phimin = 0.*deg, fullconemp_phimax = 360.*deg;
  G4Cons* solidFullconemp =    
    new G4Cons("Fullconemp", 
    fullconemp_rmina, fullconemp_rmaxa, fullconemp_rminb, fullconemp_rmaxb, fullconemp_hz,
    fullconemp_phimin, fullconemp_phimax);


  //
  // right side Pd (reflector)  
  //  
  G4Material* rside_reflector = nist->FindOrBuildMaterial("G4_Pb"); //Pb
  G4ThreeVector pos_rside_reflector = G4ThreeVector(0, 0, ((0.5+36.5+0.5+1.0+0.5)/2)*cm);
        
  G4double rsbmp_sizeXY = 70.*cm;
  G4double rsbmp_sizeZ  = (0.5+36.5+0.5+1.0+0.5)*cm;     

  G4Box* solidRsbmp =    
   new G4Box("RsideEnvelopemp",                    //its name
       0.5*rsbmp_sizeXY, 0.5*rsbmp_sizeXY, 0.5*rsbmp_sizeZ); //its size

  G4SubtractionSolid* RsbsubtconeminusPL =
   new G4SubtractionSolid("RsbsubtconeminusPL", solidRsbmp, solidFullconemp);


  G4LogicalVolume* logicRsbsubtconeminusPL =                         
    new G4LogicalVolume(RsbsubtconeminusPL,         //its solid
                        rside_reflector,          //its material
                        "RsbsubtconeminusPL");               //its name
  logicRsbsubtconeminusPL->SetVisAttributes(new G4VisAttributes(G4Colour(204/255.,0/255.,204/255.)));
              
  new G4PVPlacement(0,                       //no rotation
                    pos_rside_reflector,                    //at position
                    logicRsbsubtconeminusPL,             //its logical volume
                    "RsbsubtconeminusPL",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  //
  // poly-Li(delimiter) (polyethylene)
  //  
  G4Material* delimiter_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE"); //POLYETHYLENE
  G4ThreeVector pos_delimiter = G4ThreeVector(0, 0, (0.5+36.5+0.5+1.0+0.5+(2.0/2))*cm);
        
  // poly-Li(delimiter) length
  G4double delimiter_sizeXY = 70.*cm;
  G4double delimiter_sizeZ  = 2.*cm;   

  // last section cone 
      
  G4double lasectioncone_rmina =  0.*cm, lasectioncone_rmaxb = (12.0/2)*cm;
  G4double lasectioncone_rminb =  0.*cm, lasectioncone_rmaxa = (12.0/2)*cm;
  G4double lasectioncone_hz = (2.0/2)*cm;
  G4double lasectioncone_phimin = 0.*deg, lasectioncone_phimax = 360.*deg;
  G4Cons* solidLasectioncone =    
    new G4Cons("Lasectioncone", 
    lasectioncone_rmina, lasectioncone_rmaxa, lasectioncone_rminb, lasectioncone_rmaxb, lasectioncone_hz,
    lasectioncone_phimin, lasectioncone_phimax);  


  G4Box* solidlastlid =    
   new G4Box("lastlid",                    //its name
       0.5*delimiter_sizeXY, 0.5*delimiter_sizeXY, 0.5*delimiter_sizeZ); //its size

  G4SubtractionSolid* DelimiterLid =
   new G4SubtractionSolid("DelimiterLid", solidlastlid, solidLasectioncone);


  G4LogicalVolume* logicDelimiterLid =                         
    new G4LogicalVolume(DelimiterLid,         //its solid
                        delimiter_mat,          //its material
                        "DelimiterLid");               //its name
  logicDelimiterLid->SetVisAttributes(new G4VisAttributes(G4Colour(255/255.,0/255.,127/255.)));
              
//  new G4PVPlacement(0,                       //no rotation
//                    pos_delimiter,                    //at position
//                    logicDelimiterLid,             //its logical volume
//                    "DelimiterLid",                //its name
//                    logicWorld,                //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking




  //
  //BSA exit detector
  //

  G4double BSAeDetector_innerRadius = 0.*cm; G4double BSAeDetector_outerRadius = (12.0/2)*cm; 
//  G4double BSAeDetector_hz = (1.0/2)*cm;
  G4double BSAeDetector_hz = (0.1/2)*cm;

  G4double BSAeDetector_startAngle = 0.*deg; G4double BSAeDetector_spanningAngle = 360.*deg;

  G4Material* BSA_exit_detector_mat = nist->FindOrBuildMaterial("G4_POLYSTYRENE"); //POLYSTYRENE

//  G4ThreeVector pos_BSA_exit_detector = G4ThreeVector(0, 0, (0.5+36.5+0.5+1.0+0.5+2.0+(1.0/2))*cm);

  G4ThreeVector pos_BSA_exit_detector = G4ThreeVector(0, 0, (0.0)*cm);  


  G4Tubs* solidBSAeDetector
     = new G4Tubs("BSAeDetector",
               BSAeDetector_innerRadius,
               BSAeDetector_outerRadius,
               BSAeDetector_hz,
               BSAeDetector_startAngle,
               BSAeDetector_spanningAngle);   


                      
  G4LogicalVolume* logicBSAeDetector =                         
    new G4LogicalVolume(solidBSAeDetector,         //its solid 
                        BSA_exit_detector_mat,          //its material
                    //    vacuum,          //material: vacuum
                        "BSAeDetector");           //its name
  logicBSAeDetector->SetVisAttributes(new G4VisAttributes(G4Colour:: Cyan()));

               
  new G4PVPlacement(0,                       //no rotation
                    pos_BSA_exit_detector,                    //at position
                    logicBSAeDetector,             //its logical volume
                    "BSAeDetector",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicBSAeDetector;



  
  //
  // left side reflector
  //

  //BeamPath box (27cm * 27cm, 46cm long)  
  G4double BeamPath_sizeYZ = 20.*cm;
  //G4double BeamPath_sizeYZ = 27.*cm;
  G4double BeamPath_sizeX  = 46.*cm;  

  G4Box* BeamPath =    
   new G4Box("BeamPathbox",                    //its name
       0.5*BeamPath_sizeX, 0.5*BeamPath_sizeYZ, 0.5*BeamPath_sizeYZ); //its size


  // left side box 
  G4double leftsb_sizeXY = 70.*cm;
  G4double leftsb_sizeZ  = (0.5+36.5+0.5+1.0+0.5)*cm; 
   
  G4Box* solidLeftsb =    
   new G4Box("Leftsb",                    //its name
       0.5*leftsb_sizeXY, 0.5*leftsb_sizeXY, 0.5*leftsb_sizeZ); //its size


  //positon and orientation of BeamTarget relative to left box 
  G4RotationMatrix* LsideyBeamRot = new G4RotationMatrix; // Rotates X and Z axes only 
  //LsideyBeamRot->rotateY(-(M_PI/4.)*rad); // Rotates 45 degrees 
  LsideyBeamRot->rotateY(0.0*rad);  

    // yBeamRot->rotateY(0.*rad); // Rotates 45 degrees 
    // G4ThreeVector zBeamTrans(0, 0, 0);   
    //G4ThreeVector zBeamTrans(((((-84-2)/2.)+1.)/sqrt(2))*cm, 0, (-((14+82+2)/2.)+1.)*cm); 
    //G4ThreeVector LsidezxBeamTrans(((-((84+target_thick)/2.)+(target_thick)/2.)/(sqrt(2)))*cm, 0, 0.*cm);


  //G4ThreeVector LsidezxBeamTrans(-(42./(sqrt(2)))*cm, 0, -1.5*cm);    
  G4ThreeVector LsidezxBeamTrans(-(leftsb_sizeXY/2 - BeamPath_sizeX/2), 0, (leftsb_sizeZ/2 - BeamPath_sizeYZ/2)); 


  //Left side box subtract beamplustatget
  G4SubtractionSolid* LeftsbminusBeam =
    new G4SubtractionSolid("LeftsbminusBeam", solidLeftsb, BeamPath, LsideyBeamRot, LsidezxBeamTrans);
     

  G4ThreeVector lside_refle_minus_beam = G4ThreeVector(0, 0, -(leftsb_sizeZ/2));

  G4LogicalVolume* logicLeftsbminusBeam =                         
    new G4LogicalVolume(LeftsbminusBeam,         //its solid
                        rside_reflector,          //its material
                        "LeftsbminusBeam");               //its name
  logicLeftsbminusBeam->SetVisAttributes(new G4VisAttributes(G4Colour(204/255.,0/255.,204/255.)));
              
  new G4PVPlacement(0,                       //no rotation
                    lside_refle_minus_beam,                    //at position
                    logicLeftsbminusBeam,             //its logical volume
                    "LeftsbminusBeam",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking




  //
  //proton to neutron conversion target Beryllium
  //

  //G4Material* tar_mat = nist->FindOrBuildMaterial("G4_Li"); //Li
  G4Material* tar_matBe = nist->FindOrBuildMaterial("G4_Be"); //Beryllium


  G4double Target_sizeXY = 16.*cm;
  //G4double Target_sizeZ  = (0.01)*cm; // 100 microns    
  G4double Target_sizeZ  = (0.1)*cm; // 1 mm 

  G4Box* solidPNConversionTarget =    
   new G4Box("PNConversionTarget",                    //its name
       0.5*Target_sizeXY, 0.5*Target_sizeXY, 0.5*Target_sizeZ); //its size

                      
  G4LogicalVolume* logicPNConversionTarget =                         
    new G4LogicalVolume(solidPNConversionTarget,         //its solid
                        tar_matBe,          //its material Beryllium                         
                        "PNConversionTarget");           //its name
  logicPNConversionTarget->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));


  G4RotationMatrix* yTarRot_1 = new G4RotationMatrix; // Rotates Y axe only 
   //yTarRot->rotateY(-M_PI/4.*rad); // Rotates 45 degrees 
   yTarRot_1->rotateY(-M_PI/6.*rad); // Rotates 30 degrees

  G4RotationMatrix* yTarRot_2 = new G4RotationMatrix; // Rotates Y axe only 
   yTarRot_2->rotateY(M_PI/6.*rad); // Rotates 30 degrees   


  G4ThreeVector pos_tar_1 = G4ThreeVector(0, 0, -(BeamPath_sizeYZ/2 - Target_sizeXY/4));
  G4ThreeVector pos_tar_2 = G4ThreeVector(0, 0, -(BeamPath_sizeYZ/2 + Target_sizeXY/4));
               
  new G4PVPlacement(yTarRot_1,                       //rotation
                    pos_tar_1,                    //at position
                    logicPNConversionTarget,             //its logical volume
                    "PNConversionTarget_1",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  new G4PVPlacement(yTarRot_2,                       //rotation
                    pos_tar_2,                    //at position
                    logicPNConversionTarget,             //its logical volume
                    "PNConversionTarget_1",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking             
 







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //     
  // Shifter
  //  
  // build spectrum shifter material TiF3

  G4int ncomponents, natoms;
  G4double TiF3density = 3.4*g/cm3;
  G4Element* elemTi = nist->FindOrBuildElement("Ti");
  //G4Element* elemF = nist->FindOrBuildElement("F");

  G4Material* shifter_mat = new G4Material("MatConTiF3", TiF3density, ncomponents = 2);
  shifter_mat -> AddElement(elemTi, natoms = 1);
  shifter_mat -> AddElement(elemF, natoms = 3);

  //G4double target_thick  = (0.1); //cm

  //G4ThreeVector pos_shifter = G4ThreeVector(0, 0, (14/2)*cm);

  ////beam(10cm diameter, 82cm long) plus target(2cm lithium)Tub
  // G4double BeaminnerRadius = 0.*cm; G4double BeamouterRadius = 10/2.*cm; G4double Beamhz = (84+target_thick)/2.*cm;
  // G4double BeamstartAngle = 0.*deg; G4double BeamspanningAngle = 360.*deg;
  // G4Tubs* BeamPlusTar
  //   = new G4Tubs("BeamPlusTar",
  //             BeaminnerRadius,
  //             BeamouterRadius,
  //             Beamhz,
  //             BeamstartAngle,   
  //             BeamspanningAngle);


  ////positon and orientation of BeamTarget "relative to shiter cone" 
  // G4RotationMatrix* yBeamRot = new G4RotationMatrix; // Rotates X and Z axes only 
  //  yBeamRot->rotateY(-M_PI/4.*rad); // Rotates 45 degrees 

   
  // G4ThreeVector zxBeamTrans(((-((14+84+target_thick)/2.)+(target_thick)/2.)/(sqrt(2)))*cm, 0, (-((84+target_thick)/2.)+(target_thick)/2.)*cm);
        
  // Conical section for shifter       
  //G4double shifterwosubeam_rmina =  0.*cm, shifterwosubeam_rmaxb = (6+43/3.)*cm;
  //G4double shifterwosubeam_rminb =  0.*cm, shifterwosubeam_rmaxa = 25.*cm;
  //G4double shifterwosubeam_hz = (14/2)*cm;
  //G4double shifterwosubeam_phimin = 0.*deg, shifterwosubeam_phimax = 360.*deg;
  //G4Cons* solidShifterwosubeam =    
  //  new G4Cons("Shifterwosubeam", 
  //  shifterwosubeam_rmina, shifterwosubeam_rmaxa, shifterwosubeam_rminb, shifterwosubeam_rmaxb, 
  //  shifterwosubeam_hz, shifterwosubeam_phimin, shifterwosubeam_phimax);

  //Conical section for shifter subtract beamplustatget
  //G4SubtractionSolid* solidShifter =
  //  new G4SubtractionSolid("solidShifter", solidShifterwosubeam, BeamPlusTar, yBeamRot, zxBeamTrans);

                  
  //G4LogicalVolume* logicShifter =                         
  //  new G4LogicalVolume(solidShifter,         //its solid
  //                      shifter_mat,          //its material
  //                      "Shifter");           //its name
  //logicShifter->SetVisAttributes(new G4VisAttributes(G4Colour::Grey()));
               
  //new G4PVPlacement(0,                       //no rotation
  //                  pos_shifter,                    //at position
  //                  logicShifter,             //its logical volume
  //                  "Shifter",                //its name
  //                  logicWorld,                //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking




  ////     
  //// Fast_Neutron_Filter
  ////  
  //G4Material* fnfilter_mat = nist->FindOrBuildMaterial("G4_Ni"); //Ni
  //G4ThreeVector pos_fnfilter = G4ThreeVector(0, 0, (14+30/2)*cm);
        
  //// Conical section for fast neutron filter       
  //G4double fnfilter_rmina =  0.*cm, fnfilter_rmaxb = (6+13/3.)*cm;
  //G4double fnfilter_rminb =  0.*cm, fnfilter_rmaxa = (6+43/3.)*cm;
  //G4double fnfilter_hz = (30/2)*cm;
  //G4double fnfilter_phimin = 0.*deg, fnfilter_phimax = 360.*deg;
  //G4Cons* solidFnfilter =    
    //new G4Cons("Fnfilter", 
    //fnfilter_rmina, fnfilter_rmaxa, fnfilter_rminb, fnfilter_rmaxb, fnfilter_hz,
    //fnfilter_phimin, fnfilter_phimax);
                      
  //G4LogicalVolume* logicFnfilter =                         
    //new G4LogicalVolume(solidFnfilter,         //its solid
    //                    fnfilter_mat,          //its material
     //                   "Fnfilter");           //its name
  //logicFnfilter->SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
               
  //new G4PVPlacement(0,                       //no rotation
  //                  pos_fnfilter,                    //at position
  //                  logicFnfilter,             //its logical volume
  //                  "Fnfilter",                //its name
  //                  logicWorld,                //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking



//// thermal_neutron_absorber
  //  
  //G4Material* thneuabs_mat = nist->FindOrBuildMaterial("G4_Cd"); //Cd
  //G4ThreeVector pos_thneuabs = G4ThreeVector(0, 0, (14+30+3.5+0.1/2)*cm);
        
  //// Conical section for thermal_neutron_absorber      
  //G4double thneuabs_rmina =  0.*cm, thneuabs_rmaxb = (6+9.4/3.)*cm;
  //G4double thneuabs_rminb =  0.*cm, thneuabs_rmaxa = (6+9.5/3.)*cm;
  //G4double thneuabs_hz = (0.1/2)*cm;
  //G4double thneuabs_phimin = 0.*deg, thneuabs_phimax = 360.*deg;
  //G4Cons* solidThneuabs =    
  //  new G4Cons("Thneuabs", 
  //  thneuabs_rmina, thneuabs_rmaxa, thneuabs_rminb, thneuabs_rmaxb, thneuabs_hz,
  //  thneuabs_phimin, thneuabs_phimax);
                      
  //G4LogicalVolume* logicThneuabs =                         
  //  new G4LogicalVolume(solidThneuabs,         //its solid
  //                      thneuabs_mat,          //its material
  //                      "Thneuabs");           //its name
  //logicThneuabs->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));

               
  //new G4PVPlacement(0,                       //no rotation
  //                  pos_thneuabs,                    //at position
  //                  logicThneuabs,             //its logical volume
  //                  "Thneuabs",                //its name
  //                  logicWorld,                //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking


  


  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
