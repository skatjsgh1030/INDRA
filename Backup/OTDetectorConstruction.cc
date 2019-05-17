#include "OTDetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iostream>
#include "TMath.h"

using std::stringstream;
using namespace std;
using namespace CLHEP;
		OTDetectorConstruction::OTDetectorConstruction()
: G4VUserDetectorConstruction()
{
}

OTDetectorConstruction::~OTDetectorConstruction()
{
}

G4VPhysicalVolume* OTDetectorConstruction::Construct()
{
		//		auto runManager = (G4RunManager *) G4RunManager::GetRunManager();

		G4NistManager* nist = G4NistManager::Instance();

		const G4double labTemp = CLHEP::STP_Temperature + 20.*kelvin;

		G4Element*  elC  = new G4Element("Carbon",   "C",  6,  12.011*g/mole);
		G4Element*  elH  = new G4Element("Hydrogen", "H",  1,  1.00794*g/mole);

		G4Material* Methylene = new G4Material("Methylene", 6.262e-04*g/CLHEP::cm3, 2, kStateGas, labTemp);
		Methylene -> AddElement(elC, 1);
		Methylene -> AddElement(elH, 2);

		// -----------------------------------------------------
		// World

		G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");
		G4double world_size = 2000*cm;

		G4Box* solidWorld =
				new G4Box("World",                       // its name
								5*world_size,                // half x
								5*world_size,                // half y
								5*world_size);               // half zNDDetectorConstruction::LANDDetectorConstruction()


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
								55,                     //copy number
								true);                 //overlaps checking


		// -----------------------------------------------------
		// Detector

		G4Material* scintillator_mat = nist -> FindOrBuildMaterial("G4_XYLENE");
		G4Material* veto_mat = nist -> FindOrBuildMaterial("G4_POLYETHYLENE");
		G4Material* tube_mat = nist -> FindOrBuildMaterial("G4_Pyrex_Glass");
		//G4Material* target_mat = nist -> FindOrBuildMaterial("Methylene");
		G4double detector_size = 1*mm;
		G4double scintillator_offset_z = (63.5/2)*mm;
		G4double veto_offset_z = (3005-420)*mm;
		G4double tube_offset_z = (63.5/2)*mm;


		G4double sizeX = 2000 * detector_size;
		G4double sizeY = 76.2 * detector_size;
		G4double sizeZ = 63.5 * detector_size;
		//target
	/*	G4Box* target =
				new G4Box("target",
								150*mm/2,
								150*mm/2,
								1*mm/2);
		G4LogicalVolume* logicalDetector3 =
				new G4LogicalVolume(target,
								target_mat,
								"target");

		new G4PVPlacement(0,
						G4ThreeVector(0,0,1*mm),
						logicalDetector3,
						"target",
						logicWorld,
						false,
						50,
						true);*/

		//hollow pyrex glass tube
		G4Box* outerTube =
				new G4Box("outerTube",
								sizeX/2,
								sizeY/2,
								sizeZ/2);

		G4Box* innerTube =
				new G4Box("innerTube",
								(2000*mm-6*mm)/2,
								(76.2*mm-6*mm)/2,
								(63.5*mm-6*mm)/2);

		G4SubtractionSolid*hollowTube =
				new G4SubtractionSolid("hollowTube",
								outerTube,
								innerTube);

		/*G4double AngleNN = -27.07*CLHEP::degree;
		  G4RotationMatrix* yRot_SS = new G4RotationMatrix;
		  yRot_SS -> rotateY(AngleNN);*/
		G4LogicalVolume* logicalDetector =
				new G4LogicalVolume(hollowTube,
								tube_mat,
								"hollowTube");
		new G4PVPlacement(0,
						G4ThreeVector(/*-1*tube_offset_z*TMath::Sin(AngleNN)*/0*mm,0*mm,tube_offset_z/**TMath::Cos(AngleNN)*/),
						logicalDetector,
						"hollowTube",
						logicWorld,
						false,
						4,
						true);
		//liquide Scintillator

		/*G4double AngleN = -27.07*CLHEP::degree;
		  G4RotationMatrix* yRot_S = new G4RotationMatrix;
		  yRot_S -> rotateY(AngleN);*/
		G4Box* lScintillator =
				new G4Box("lScintillator",
								(2000*mm-6*mm)/2,
								(76.2*mm-6*mm)/2,
								(63.5*mm-6*mm)/2);

		G4LogicalVolume* logicalDetector1 =
				new G4LogicalVolume(lScintillator,
								scintillator_mat,
								"lScintillator");
		/*	auto pvp = */new G4PVPlacement(/*yRot_S*/0,
						G4ThreeVector(/*-1*scintillator_offset_z*TMath::Sin(AngleN)*/0*mm,0*mm,scintillator_offset_z/**TMath::Cos(AngleN)*/),
						logicalDetector1,
						"lScintillator",
						logicWorld,
						false,
						2,
						true);
		//	runManager -> SetSensitiveDetector(pvp);
		G4Box* veto =
				new G4Box("veto",
								90*mm/2,
								2000*mm/2,
								10*mm/2);

		G4LogicalVolume* logicalDetector2 =
				new G4LogicalVolume(veto,
								veto_mat,
								"veto");

		/*	pvp = new*/ G4PVPlacement(0,//pi,theta,psi
						G4ThreeVector(/*-1*(veto_offset_z*TMath::Sin(AngelF))*/0*mm,0*mm,veto_offset_z/**TMath::Cos(AngelF)*/),
						logicalDetector2,
						"veto",
						logicWorld,
						false,
						1,
						true);

		//	runManager -> SetSensitiveDetector(pvp);


		logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());

		//G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
		//targetVisAtt -> SetForceWireframe(true);
		//logicalDetector3 -> SetVisAttributes (targetVisAtt);
		G4VisAttributes* vetoVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
		vetoVisAtt -> SetForceWireframe(true);
		logicalDetector2 -> SetVisAttributes (vetoVisAtt);
		G4VisAttributes* hollowTubeVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
		hollowTubeVisAtt -> SetForceWireframe(true);
		logicalDetector -> SetVisAttributes (hollowTubeVisAtt);

		return physWorld;
}

