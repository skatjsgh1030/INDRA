#include "INDRADetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4VPhysicalVolume.hh"
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
  INDRADetectorConstruction::INDRADetectorConstruction()
: G4VUserDetectorConstruction()
{
}

INDRADetectorConstruction::~INDRADetectorConstruction()
{
}

G4VPhysicalVolume* INDRADetectorConstruction::Construct()
{
  //----------------Material --------------------
  G4NistManager* nist = G4NistManager::Instance();

  const G4double labTemp = CLHEP::STP_Temperature + 20.*kelvin;

  G4Element*  elCs  = new G4Element("Caesium","Cs",  55,   132.90545*g/mole);
  G4Element*  elI   = new G4Element("Iodine", "I" , 53,   126.90447*g/mole);

  G4Material* CsI = new G4Material("CsI", 4.51*g/CLHEP::cm3, 2, kStateSolid, labTemp);
  CsI -> AddElement(elCs, 1);
  CsI -> AddElement(elI,  1);

  // ------------ Generate & Add Material Properties Table ------------
  /*/
	G4double photonEnergy[] =
	{	1.90738*eV, 1.98368*eV, 2.06633*eV,
	2.15617*eV, 2.25418*eV, 2.36152*eV,
	2.47960*eV, 2.61011*eV, 2.75511*eV,
	2.91718*eV, 3.0995*eV};
	const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  // CsI
  //fixed
  G4double refractiveIndex1[] =
  { 1.7789, 1.7820, 1.7856, 1.7897, 1.7945,
  1.8001, 1.8067, 1.8146, 1.8242, 1.8359,
  1.8506};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
  {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
  30.000*m, 28.500*m, 27.000*m,17.500*m, 14.500*m };

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
  { 1.00, 1.00, 1.00, 
  1.00, 1.00, 1.00, 
  1.00,1.00, 1.00,
  1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
  { 0.01, 1.00, 2.00,
  3.00, 4.00, 5.00,
  6.00,7.00, 6.00,
  5.00, 4.00 };

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
  ->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
  ->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
  ->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
  ->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",2./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  G4double energy_CsI_Cry[] = {
  1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
  1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
  1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
  1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
  1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
  2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
  2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
  2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
  2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
  2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
  3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
  3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
  3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
  4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
  5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
};

const G4int numentries_CsI_Cry = sizeof(energy_CsI_Cry)/sizeof(G4double);

//assume 100 times larger than the rayleigh scattering for now.
G4double mie_CsI_Cry[] = {
  167024.4*m, 158726.7*m, 150742  *m,
  143062.5*m, 135680.2*m, 128587.4*m,
  121776.3*m, 115239.5*m, 108969.5*m,
  102958.8*m, 97200.35*m, 91686.86*m,
  86411.33*m, 81366.79*m, 76546.42*m,
  71943.46*m, 67551.29*m, 63363.36*m,
  59373.25*m, 55574.61*m, 51961.24*m,
  48527.00*m, 45265.87*m, 42171.94*m,
  39239.39*m, 36462.50*m, 33835.68*m,
  31353.41*m, 29010.30*m, 26801.03*m,
  24720.42*m, 22763.36*m, 20924.88*m,
  19200.07*m, 17584.16*m, 16072.45*m,
  14660.38*m, 13343.46*m, 12117.33*m,
  10977.70*m, 9920.416*m, 8941.407*m,
  8036.711*m, 7202.470*m, 6434.927*m,
  5730.429*m, 5085.425*m, 4496.467*m,
  3960.210*m, 3473.413*m, 3032.937*m,
  2635.746*m, 2278.907*m, 1959.588*m,
  1675.064*m, 1422.710*m, 1200.004*m,
  1004.528*m, 833.9666*m, 686.1063*m
};

assert(sizeof(mie_CsI_Cry) == sizeof(energy_CsI_Cry));

// gforward, gbackward, forward backward ratio
G4double mie_CsI_Cry_const[3]={0.99,0.99,0.8};

myMPT1->AddProperty("MIEHG",energy_CsI_Cry,mie_CsI_Cry,numentries_CsI_Cry)
->SetSpline(true);
myMPT1->AddConstProperty("MIEHG_FORWARD",mie_CsI_Cry_const[0]);
myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_CsI_Cry_const[1]);
myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_CsI_Cry_const[2]);

G4cout << "CsI G4MaterialPropertiesTable" << G4endl;
myMPT1->DumpTable();

CsI->SetMaterialPropertiesTable(myMPT1);

// Set the Birks Constant for the Water scintillator

CsI->GetIonisation()->SetBirksConstant(0.126*mm/MeV);*/

// -----------------------------------------------------
// World

G4Material* world_mat = nist -> FindOrBuildMaterial("G4_Galactic");
//G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");		
G4double world_size = 5*m;

G4Box* solidWorld =
new G4Box("World",                       // its name
	3*world_size,                // half x
	3*world_size,                // half y
	3*world_size);               // half zNDDetectorConstruction::LANDDetectorConstruction()


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
	true);                 //overlaps checking


// -----------------------------------------------------
G4Material* scintillator_mat = nist -> FindOrBuildMaterial("CsI");
G4Material* Silicon_mat = nist -> FindOrBuildMaterial("G4_Si");

const int N_8 = 8;
const int N_16 = 16;
const int N_24 = 24;
const double gap = 0.3;
//--------------------Si---------
//---------------------------6------------------
double thetaSB_6 = 20.0;//degree
double thetaSF_6 = 14.0;
double distS_6 = 397.0;//mm
double thickS_6 = 0.3;
double thetaS_6 = (thetaSB_6 - thetaSF_6)/2;
double thetaSP_6 = (thetaSB_6 + thetaSF_6)/2;
double hXS1_6 = (distS_6/(TMath::Cos(thetaS_6*degree)))*TMath::Sin(thetaSF_6*degree);
double hXS2_6 = (distS_6/(TMath::Cos(thetaS_6*degree)))*TMath::Sin(thetaSB_6*degree);
double hXS3_6 = ((distS_6+thickS_6)/(TMath::Cos(thetaS_6*degree)))*TMath::Sin(thetaSF_6*degree);
double hXS4_6 = ((distS_6+thickS_6)/(TMath::Cos(thetaS_6*degree)))*TMath::Sin(thetaSB_6*degree);
//----------------------------7------------------
double thetaSB_7 = 27.0;//degree
double thetaSF_7 = 20.3;
double distS_7 = 397.0;//mm
double thickS_7 = 0.3;
double thetaS_7 = (thetaSB_7 - thetaSF_7)/2;
double thetaSP_7 = (thetaSB_7 + thetaSF_7)/2;
double hXS1_7 = (distS_7/(TMath::Cos(thetaS_7*degree)))*TMath::Sin(thetaSF_7*degree);
double hXS2_7 = (distS_7/(TMath::Cos(thetaS_7*degree)))*TMath::Sin(thetaSB_7*degree);
double hXS3_7 = ((distS_7+thickS_7)/(TMath::Cos(thetaS_7*degree)))*TMath::Sin(thetaSF_7*degree);
double hXS4_7 = ((distS_7+thickS_7)/(TMath::Cos(thetaS_7*degree)))*TMath::Sin(thetaSB_7*degree);
//-----------------------------8------------------------------
double thetaSB_8 = 35.0;//degree
double thetaSF_8 = 27.5;
double distS_8 = 267.0;//mm
double thickS_8 = 0.3;
double thetaS_8 = (thetaSB_8 - thetaSF_8)/2;
double thetaSP_8 = (thetaSB_8 + thetaSF_8)/2;
double hXS1_8 = (distS_8/(TMath::Cos(thetaS_8*degree)))*TMath::Sin(thetaSF_8*degree)*mm;
double hXS2_8 = (distS_8/(TMath::Cos(thetaS_8*degree)))*TMath::Sin(thetaSB_8*degree)*mm;
double hXS3_8 = ((distS_8+thickS_8)/(TMath::Cos(thetaS_8*degree)))*TMath::Sin(thetaSF_8*degree)*mm;
double hXS4_8 = ((distS_8+thickS_8)/(TMath::Cos(thetaS_8*degree)))*TMath::Sin(thetaSB_8*degree)*mm;
//-----------------------------9-----------------
double thetaSB_9 = 45.0;//degree
double thetaSF_9 = 35.5;
double distS_9 = 267.0;//mm
double thickS_9 = 0.3;
double thetaS_9 = (thetaSB_9 - thetaSF_9)/2;
double thetaSP_9 = (thetaSB_9 + thetaSF_9)/2;
double hXS1_9 = (distS_9/(TMath::Cos(thetaS_9*degree)))*TMath::Sin(thetaSF_9*degree)*mm;
double hXS2_9 = (distS_9/(TMath::Cos(thetaS_9*degree)))*TMath::Sin(thetaSB_9*degree)*mm;
double hXS3_9 = ((distS_9+thickS_9)/(TMath::Cos(thetaS_9*degree)))*TMath::Sin(thetaSF_9*degree)*mm;
double hXS4_9 = ((distS_9+thickS_9)/(TMath::Cos(thetaS_9*degree)))*TMath::Sin(thetaSB_9*degree)*mm;

//CsI crystal
//---------------------------6------------------
double thetaB_6 = 20.0;//degree
double thetaF_6 = 14.0;
double dist_6 = 400.0;//mm
double thick_6 = 97.0;
double theta_6 = (thetaB_6 - thetaF_6)/2;
double thetaP_6 = (thetaB_6 + thetaF_6)/2;
double hX1_6 = (dist_6/(TMath::Cos(theta_6*degree)))*TMath::Sin(thetaF_6*degree);
double hX2_6 = (dist_6/(TMath::Cos(theta_6*degree)))*TMath::Sin(thetaB_6*degree);
double hX3_6 = ((dist_6+thick_6)/(TMath::Cos(theta_6*degree)))*TMath::Sin(thetaF_6*degree);
double hX4_6 = ((dist_6+thick_6)/(TMath::Cos(theta_6*degree)))*TMath::Sin(thetaB_6*degree);
//----------------------------7------------------
double thetaB_7 = 27.0;//degree
double thetaF_7 = 20.3;
double dist_7 = 400.0;//mm
double thick_7 = 97.0;
double theta_7 = (thetaB_7 - thetaF_7)/2;
double thetaP_7 = (thetaB_7 + thetaF_7)/2;
double hX1_7 = (dist_7/(TMath::Cos(theta_7*degree)))*TMath::Sin(thetaF_7*degree);
double hX2_7 = (dist_7/(TMath::Cos(theta_7*degree)))*TMath::Sin(thetaB_7*degree);
double hX3_7 = ((dist_7+thick_7)/(TMath::Cos(theta_7*degree)))*TMath::Sin(thetaF_7*degree);
double hX4_7 = ((dist_7+thick_7)/(TMath::Cos(theta_7*degree)))*TMath::Sin(thetaB_7*degree);
//-----------------------------8------------------------------
double thetaB_8 = 35.0;//degree
double thetaF_8 = 27.5;
double dist_8 = 270.0;//mm
double thick_8 = 90.0;
double theta_8 = (thetaB_8 - thetaF_8)/2;
double thetaP_8 = (thetaB_8 + thetaF_8)/2;
double hX1_8 = (dist_8/(TMath::Cos(theta_8*degree)))*TMath::Sin(thetaF_8*degree)*mm;
double hX2_8 = (dist_8/(TMath::Cos(theta_8*degree)))*TMath::Sin(thetaB_8*degree)*mm;
double hX3_8 = ((dist_8+thick_8)/(TMath::Cos(theta_8*degree)))*TMath::Sin(thetaF_8*degree)*mm;
double hX4_8 = ((dist_8+thick_8)/(TMath::Cos(theta_8*degree)))*TMath::Sin(thetaB_8*degree)*mm;
//-----------------------------9-----------------
double thetaB_9 = 45.0;//degree
double thetaF_9 = 35.5;
double dist_9 = 270.0;//mm
double thick_9 = 90.0;
double theta_9 = (thetaB_9 - thetaF_9)/2;
double thetaP_9 = (thetaB_9 + thetaF_9)/2;
double hX1_9 = (dist_9/(TMath::Cos(theta_9*degree)))*TMath::Sin(thetaF_9*degree)*mm;
double hX2_9 = (dist_9/(TMath::Cos(theta_9*degree)))*TMath::Sin(thetaB_9*degree)*mm;
double hX3_9 = ((dist_9+thick_9)/(TMath::Cos(theta_9*degree)))*TMath::Sin(thetaF_9*degree)*mm;
double hX4_9 = ((dist_9+thick_9)/(TMath::Cos(theta_9*degree)))*TMath::Sin(thetaB_9*degree)*mm;
//---------------------------10------------------
double thetaB_10 = 57.0;//degree
double thetaF_10 = 45.5;
double dist_10 = 270.0;//mm
double thick_10 = 76.0;
double theta_10 = (thetaB_10 - thetaF_10)/2;
double thetaP_10 = (thetaB_10 + thetaF_10)/2;
double hX1_10 = (dist_10/(TMath::Cos(theta_10*degree)))*TMath::Sin(thetaF_10*degree);
double hX2_10 = (dist_10/(TMath::Cos(theta_10*degree)))*TMath::Sin(thetaB_10*degree);
double hX3_10 = ((dist_10+thick_10)/(TMath::Cos(theta_10*degree)))*TMath::Sin(thetaF_10*degree);
double hX4_10 = ((dist_10+thick_10)/(TMath::Cos(theta_10*degree)))*TMath::Sin(thetaB_10*degree);
//----------------------------
//---------------------------11------------------
double thetaB_11 = 70.0;//degree
double thetaF_11 = 57.5;
double dist_11 = 270.0;//mm
double thick_11 = 76.0;
double theta_11 = (thetaB_11 - thetaF_11)/2;
double thetaP_11 = (thetaB_11 + thetaF_11)/2;
double hX1_11 = (dist_11/(TMath::Cos(theta_11*degree)))*TMath::Sin(thetaF_11*degree);
double hX2_11 = (dist_11/(TMath::Cos(theta_11*degree)))*TMath::Sin(thetaB_11*degree);
double hX3_11 = ((dist_11+thick_11)/(TMath::Cos(theta_11*degree)))*TMath::Sin(thetaF_11*degree);
double hX4_11 = ((dist_11+thick_11)/(TMath::Cos(theta_11*degree)))*TMath::Sin(thetaB_11*degree);
//----------------------------
//---------------------------12------------------
double thetaB_12 = 88.0;//degree
double thetaF_12 = 70.0;
double dist_12 = 270.0;//mm
double thick_12 = 48.0;
double theta_12 = (thetaB_12 - thetaF_12)/2;
double thetaP_12 = (thetaB_12 + thetaF_12)/2;
double hX1_12 = (dist_12/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaF_12*degree);
double hX2_12 = (dist_12/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaB_12*degree);
double hX3_12 = ((dist_12+thick_12)/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaF_12*degree);
double hX4_12 = ((dist_12+thick_12)/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaB_12*degree);
//----------------------------
//---------------------------13------------------
double thetaB_13 = 92.0;//degree
double thetaF_13 = 109.5;
double dist_13 = 270.0;//mm
double thick_13 = 60.0;
double theta_13 = (thetaB_13 - thetaF_13)/2;
double thetaP_13 = (thetaB_13 + thetaF_13)/2;
double hX1_13 = (dist_13/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaF_13*degree);
double hX2_13 = (dist_13/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaB_13*degree);
double hX3_13 = ((dist_13+thick_13)/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaF_13*degree);
double hX4_13 = ((dist_13+thick_13)/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaB_13*degree);
//----------------------------
//---------------------------14------------------
double thetaB_14 = 110.0;//degree
double thetaF_14 = 125.0;
double dist_14 = 270.0;//mm
double thick_14 = 50.0;
double theta_14 = (thetaB_14 - thetaF_14)/2;
double thetaP_14 = (thetaB_14 + thetaF_14)/2;
double hX1_14 = (dist_14/(TMath::Cos(theta_14*degree)))*TMath::Sin(thetaF_14*degree);
double hX2_14 = (dist_14/(TMath::Cos(theta_14*degree)))*TMath::Sin(thetaB_14*degree);
double hX3_14 = ((dist_14+thick_14)/(TMath::Cos(theta_14*degree)))*TMath::Sin(thetaF_14*degree);
double hX4_14 = ((dist_14+thick_14)/(TMath::Cos(theta_14*degree)))*TMath::Sin(thetaB_14*degree);
//----------------------------
//---------------------------15------------------
double thetaB_15 = 126.0;//degree
double thetaF_15 = 140.0;
double dist_15 = 270.0;//mm
double thick_15 = 50.0;
double theta_15 = (thetaB_15 - thetaF_15)/2;
double thetaP_15 = (thetaB_15 + thetaF_15)/2;
double hX1_15 = (dist_15/(TMath::Cos(theta_15*degree)))*TMath::Sin(thetaF_15*degree);
double hX2_15 = (dist_15/(TMath::Cos(theta_15*degree)))*TMath::Sin(thetaB_15*degree);
double hX3_15 = ((dist_15+thick_15)/(TMath::Cos(theta_15*degree)))*TMath::Sin(thetaF_15*degree);
double hX4_15 = ((dist_15+thick_15)/(TMath::Cos(theta_15*degree)))*TMath::Sin(thetaB_15*degree);
//----------------------------
//---------------------------16------------------
double thetaB_16 = 142.0;//degree
double thetaF_16 = 155.0;
double dist_16 = 270.0;//mm
double thick_16 = 50.0;
double theta_16 = (thetaB_16 - thetaF_16)/2;
double thetaP_16 = (thetaB_16 + thetaF_16)/2;
double hX1_16 = (dist_16/(TMath::Cos(theta_16*degree)))*TMath::Sin(thetaF_16*degree);
double hX2_16 = (dist_16/(TMath::Cos(theta_16*degree)))*TMath::Sin(thetaB_16*degree);
double hX3_16 = ((dist_16+thick_16)/(TMath::Cos(theta_16*degree)))*TMath::Sin(thetaF_16*degree);
double hX4_16 = ((dist_16+thick_16)/(TMath::Cos(theta_16*degree)))*TMath::Sin(thetaB_16*degree);
//----------------------------
//---------------------------17------------------
double thetaB_17 = 157.0;//degree
double thetaF_17 = 176.0;
double dist_17 = 270.0;//mm
double thick_17 = 50.0;
double theta_17 = (thetaB_17 - thetaF_17)/2;
double thetaP_17 = (thetaB_17 + thetaF_17)/2;
double hX1_17 = (dist_17/(TMath::Cos(theta_17*degree)))*TMath::Sin(thetaF_17*degree);
double hX2_17 = (dist_17/(TMath::Cos(theta_17*degree)))*TMath::Sin(thetaB_17*degree);
double hX3_17 = ((dist_17+thick_17)/(TMath::Cos(theta_17*degree)))*TMath::Sin(thetaF_17*degree);
double hX4_17 = ((dist_17+thick_17)/(TMath::Cos(theta_17*degree)))*TMath::Sin(thetaB_17*degree);
//----------------------------
G4Trap* Si6[N_24];
G4LogicalVolume* logicalDetectorS6[N_24];
G4VPhysicalVolume* PhysicalDetectorS6[N_24];
G4VisAttributes* Si6VisAtt[N_24];

G4Trap* Si7[N_24];
G4LogicalVolume* logicalDetectorS7[N_24];
G4VPhysicalVolume* PhysicalDetectorS7[N_24];
G4VisAttributes* Si7VisAtt[N_24];

G4Trap* Si8[N_24];
G4LogicalVolume* logicalDetectorS8[N_24];
G4VPhysicalVolume* PhysicalDetectorS8[N_24];
G4VisAttributes* Si8VisAtt[N_24];

G4Trap* Si9[N_24];
G4LogicalVolume* logicalDetectorS9[N_24];
G4VPhysicalVolume* PhysicalDetectorS9[N_24];
G4VisAttributes* Si9VisAtt[N_24];


G4Trap* CsI6[N_24];
G4LogicalVolume* logicalDetector6[N_24];
G4VPhysicalVolume* PhysicalDetector6[N_24];
G4VisAttributes* CsI6VisAtt[N_24];

G4Trap* CsI7[N_24];
G4LogicalVolume* logicalDetector7[N_24];
G4VPhysicalVolume* PhysicalDetector7[N_24];
G4VisAttributes* CsI7VisAtt[N_24];

G4Trap* CsI8[N_24];
G4LogicalVolume* logicalDetector8[N_24];
G4VPhysicalVolume* PhysicalDetector8[N_24];
G4VisAttributes* CsI8VisAtt[N_24];

G4Trap* CsI9[N_24];
G4LogicalVolume* logicalDetector9[N_24];
G4VPhysicalVolume* PhysicalDetector9[N_24];
G4VisAttributes* CsI9VisAtt[N_24];

G4Trap* CsI10[N_24];
G4LogicalVolume* logicalDetector10[N_24];
G4VPhysicalVolume* PhysicalDetector10[N_24];
G4VisAttributes* CsI10VisAtt[N_24];

G4Trap* CsI11[N_24];
G4LogicalVolume* logicalDetector11[N_24];
G4VPhysicalVolume* PhysicalDetector11[N_24];
G4VisAttributes* CsI11VisAtt[N_24];

G4Trap* CsI12[N_24];
G4LogicalVolume* logicalDetector12[N_24];
G4VPhysicalVolume* PhysicalDetector12[N_24];
G4VisAttributes* CsI12VisAtt[N_24];

G4Trap* CsI13[N_24];
G4LogicalVolume* logicalDetector13[N_24];
G4VPhysicalVolume* PhysicalDetector13[N_24];
G4VisAttributes* CsI13VisAtt[N_24];

G4Trap* CsI14[N_16];
G4LogicalVolume* logicalDetector14[N_16];
G4VPhysicalVolume* PhysicalDetector14[N_16];
G4VisAttributes* CsI14VisAtt[N_16];

G4Trap* CsI15[N_16];
G4LogicalVolume* logicalDetector15[N_16];
G4VPhysicalVolume* PhysicalDetector15[N_16];
G4VisAttributes* CsI15VisAtt[N_16];

G4Trap* CsI16[N_8];
G4LogicalVolume* logicalDetector16[N_8];
G4VPhysicalVolume* PhysicalDetector16[N_8];
G4VisAttributes* CsI16VisAtt[N_8];

G4Trap* CsI17[N_8];
G4LogicalVolume* logicalDetector17[N_8];
G4VPhysicalVolume* PhysicalDetector17[N_8];
G4VisAttributes* CsI17VisAtt[N_8];
int Si6num[N_24] = {0,1,2,3,4,5,6,7,8,9,10,11,
0,1,2,3,4,5,6,7,8,9,10,11};
string Si6name[N_24] = {"Si600","Si601","Si602","Si603","Si604","Si605","Si606","Si607","Si608","Si609","Si610","Si611",
"Si600","Si601","Si602","Si603","Si604","Si605","Si606","Si607","Si608","Si609","Si610","Si611"};

int Si7num[N_24] = {12,13,14,15,16,17,18,19,20,21,22,23,
12,13,14,15,16,17,18,19,20,21,22,23};
string Si7name[N_24] = {"Si700","Si701","Si702","Si703","Si704","Si705","Si706","Si707","Si708","Si709","Si710","Si711",
"Si700","Si701","Si702","Si703","Si704","Si705","Si706","Si707","Si708","Si709","Si710","Si711"};

int Si8num[N_24] = {24,25,26,27,28,29,30,31,32,33,34,35,
24,25,26,27,28,29,30,31,32,33,34,35};
string Si8name[N_24] = {"Si800","Si801","Si802","Si803","Si804","Si805","Si806","Si807","Si808","Si809","Si810","Si811",
"Si800","Si801","Si802","Si803","Si804","Si805","Si806","Si807","Si808","Si809","Si810","Si811"};

int Si9num[N_24] = {36,37,38,39,40,41,42,43,44,45,46,47,
36,37,38,39,40,41,42,43,44,45,46,47};
string Si9name[N_24] = {"Si900","Si901","Si902","Si903","Si904","Si905","Si906","Si907","Si908","Si909","Si910","Si911",
"Si900","Si901","Si902","Si903","Si904","Si905","Si906","Si907","Si908","Si909","Si910","Si911"};


int CsI6num[N_24] = {0,1,2,3,4,5,6,7,8,9,10,11,
0,1,2,3,4,5,6,7,8,9,10,11};
string CsI6name[N_24] = {"CsI600","CsI601","CsI602","CsI603","CsI604","CsI605","CsI606","CsI607","CsI608","CsI609","CsI610","CsI611",
"CsI600","CsI601","CsI602","CsI603","CsI604","CsI605","CsI606","CsI607","CsI608","CsI609","CsI610","CsI611"};

int CsI7num[N_24] = {12,13,14,15,16,17,18,19,20,21,22,23,
12,13,14,15,16,17,18,19,20,21,22,23};
string CsI7name[N_24] = {"CsI700","CsI701","CsI702","CsI703","CsI704","CsI705","CsI706","CsI707","CsI708","CsI709","CsI710","CsI711",
"CsI700","CsI701","CsI702","CsI703","CsI704","CsI705","CsI706","CsI707","CsI708","CsI709","CsI710","CsI711"};

int CsI8num[N_24] = {24,25,26,27,28,29,30,31,32,33,34,35,
24,25,26,27,28,29,30,31,32,33,34,35};
string CsI8name[N_24] = {"CsI800","CsI801","CsI802","CsI803","CsI804","CsI805","CsI806","CsI807","CsI808","CsI809","CsI810","CsI811",
"CsI800","CsI801","CsI802","CsI803","CsI804","CsI805","CsI806","CsI807","CsI808","CsI809","CsI810","CsI811"};

int CsI9num[N_24] = {36,37,38,39,40,41,42,43,44,45,46,47,
36,37,38,39,40,41,42,43,44,45,46,47};
string CsI9name[N_24] = {"CsI900","CsI901","CsI902","CsI903","CsI904","CsI905","CsI906","CsI907","CsI908","CsI909","CsI910","CsI911",
"CsI900","CsI901","CsI902","CsI903","CsI904","CsI905","CsI906","CsI907","CsI908","CsI909","CsI910","CsI911"};

int CsI10num[N_24] = {36,37,38,39,40,41,42,43,44,45,46,47,
36,37,38,39,40,41,42,43,44,45,46,47};
string CsI10name[N_24] = {"CsI1000","CsI1001","CsI1002","CsI1003","CsI1004","CsI1005","CsI1006","CsI1007","CsI1008","CsI1009","CsI1010","CsI1011",
"CsI1000","CsI1001","CsI1002","CsI1003","CsI1004","CsI1005","CsI1006","CsI1007","CsI1008","CsI1009","CsI1010","CsI1011"};

int CsI11num[N_24] = {36,37,38,39,40,41,42,43,44,45,46,47,
36,37,38,39,40,41,42,43,44,45,46,47};
string CsI11name[N_24] = {"CsI1100","CsI1101","CsI1102","CsI1103","CsI1104","CsI1105","CsI1106","CsI1107","CsI1108","CsI1109","CsI1110","CsI1111",
"CsI1100","CsI1101","CsI1102","CsI1103","CsI1104","CsI1105","CsI1106","CsI1107","CsI1108","CsI1109","CsI1110","CsI1111"};

int CsI12num[N_24] = {36,37,38,39,40,41,42,43,44,45,46,47,
36,37,38,39,40,41,42,43,44,45,46,47};
string CsI12name[N_24] = {"CsI1200","CsI1201","CsI1202","CsI1203","CsI1204","CsI1205","CsI1206","CsI1207","CsI1208","CsI1209","CsI1210","CsI1211",
"CsI1200","CsI1201","CsI1202","CsI1203","CsI1204","CsI1205","CsI1206","CsI1207","CsI1208","CsI1209","CsI1210","CsI1211"};

int CsI13num[N_24] = {36,37,38,39,40,41,42,43,44,45,46,47,
36,37,38,39,40,41,42,43,44,45,46,47};
string CsI13name[N_24] = {"CsI1300","CsI1301","CsI1302","CsI1303","CsI1304","CsI1305","CsI1306","CsI1307","CsI1308","CsI1309","CsI1310","CsI1311",
"CsI1300","CsI1301","CsI1302","CsI1303","CsI1304","CsI1305","CsI1306","CsI1307","CsI1308","CsI1309","CsI1310","CsI1311"};

int CsI14num[N_16] = {36,37,38,39,40,41,42,43,44,45,46,47};
string CsI14name[N_16] = {"CsI1400","CsI1401","CsI1402","CsI1403","CsI1404","CsI1405","CsI1406","CsI1407","CsI1408","CsI1409","CsI1410","CsI1411",
"CsI1412","CsI1413","CsI1414","CsI1415"};

int CsI15num[N_16] = {36,37,38,39,40,41,42,43,44,45,46,47};
string CsI15name[N_16] = {"CsI1500","CsI1501","CsI1502","CsI1503","CsI1504","CsI1505","CsI1506","CsI1507","CsI1508","CsI1509","CsI1510","CsI1511",
"CsI1512","CsI1513","CsI1514","CsI1515"};

int CsI16num[N_8] = {36,37,38,39,40,41,42,43};
string CsI16name[N_8] = {"CsI1600","CsI1601","CsI1602","CsI1603","CsI1604","CsI1605","CsI1606","CsI1607"};

int CsI17num[N_8] = {36,37,38,39,40,41,42,43};
string CsI17name[N_8] = {"CsI1700","CsI1701","CsI1702","CsI1703","CsI1704","CsI1705","CsI1706","CsI1707"};

for(int NB = 0 ; NB < N_24 ; NB++)
{
  //-----------------------Ring 6 chio---------------------------------------34----
  CsI6[NB] = new G4Trap(CsI6name[NB],
	  thick_6/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_6*TMath::Tan(theta_6*degree)*mm,	//G4double pDy1,
	  hX1_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_6+thick_6)*TMath::Tan(theta_6*degree)*mm,	//G4double pDy2,
	  hX3_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_6;
  Rot_CsI_6 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_CsI_6 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_6*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_6*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_6 -> rotateX(1*thetaP_6*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_6 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_6*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_6*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_6 -> rotateX(thetaP_6*degree);}

  logicalDetector6[NB] = new G4LogicalVolume(CsI6[NB] ,scintillator_mat ,CsI6name[NB]);
  PhysicalDetector6[NB] = new G4PVPlacement(Rot_CsI_6,
	  G4ThreeVector((dist_6+(thick_6/2))*(TMath::Sin(thetaP_6*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_6+(thick_6/2))*(TMath::Sin(thetaP_6*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_6+(thick_6/2))*(TMath::Cos(thetaP_6*degree))*mm),
	  logicalDetector6[NB],
	  CsI6name[NB],
	  logicWorld,
	  false,
	  CsI6num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 7 chio---------------------------------------34----
  CsI7[NB] = new G4Trap(CsI7name[NB],
	  thick_7/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_7*TMath::Tan(theta_7*degree)*mm,	//G4double pDy1,
	  hX1_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_7+thick_7)*TMath::Tan(theta_7*degree)*mm,	//G4double pDy2,
	  hX3_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_7;
  Rot_CsI_7 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_CsI_7 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_7*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_7*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_7 -> rotateX(1*thetaP_7*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_7 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_7*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_7*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_7 -> rotateX(thetaP_7*degree);}

  logicalDetector7[NB] = new G4LogicalVolume(CsI7[NB] ,scintillator_mat ,CsI7name[NB]);
  PhysicalDetector7[NB] = new G4PVPlacement(Rot_CsI_7,
	  G4ThreeVector((dist_7+(thick_7/2))*(TMath::Sin(thetaP_7*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_7+(thick_7/2))*(TMath::Sin(thetaP_7*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_7+(thick_7/2))*(TMath::Cos(thetaP_7*degree))*mm),
	  logicalDetector7[NB],
	  CsI7name[NB],
	  logicWorld,
	  false,
	  CsI7num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 8 chio---------------------------------------34----
  CsI8[NB] = new G4Trap(CsI8name[NB],
	  thick_8/2*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_8*TMath::Tan(theta_8*degree)*mm,	//G4double pDy1,
	  hX1_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_8+thick_8)*TMath::Tan(theta_8*degree)*mm,	//G4double pDy2,
	  hX3_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)
  G4RotationMatrix* Rot_CsI_8;
  Rot_CsI_8 = new G4RotationMatrix;
  
  if(0 <= NB && NB <=11){
  Rot_CsI_8 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_8*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_8*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_8 -> rotateX(1*thetaP_8*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_8 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_8*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_8*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_8 -> rotateX(thetaP_8*degree);}

  logicalDetector8[NB] = new G4LogicalVolume(CsI8[NB] ,scintillator_mat ,CsI8name[NB]);
  PhysicalDetector8[NB] = new G4PVPlacement(Rot_CsI_8,
	  G4ThreeVector((dist_8+(thick_8/2))*(TMath::Sin(thetaP_8*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_8+(thick_8/2))*(TMath::Sin(thetaP_8*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_8+(thick_8/2))*(TMath::Cos(thetaP_8*degree))*mm),
	  logicalDetector8[NB],
	  CsI8name[NB],
	  logicWorld,
	  false,
	  CsI8num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 9 chio---------------------------------------34----
  CsI9[NB] = new G4Trap(CsI9name[NB],
	  thick_9/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_9*TMath::Tan(theta_9*degree)*mm,	//G4double pDy1,
	  hX1_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_9+thick_9)*TMath::Tan(theta_9*degree)*mm,	//G4double pDy2,
	  hX3_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_9;
  Rot_CsI_9 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_CsI_9 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_9*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_9*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_9 -> rotateX(1*thetaP_9*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_9 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_9*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_9*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_9 -> rotateX(thetaP_9*degree);}

  logicalDetector9[NB] = new G4LogicalVolume(CsI9[NB] ,scintillator_mat ,CsI9name[NB]);
  PhysicalDetector9[NB] = new G4PVPlacement(Rot_CsI_9,
	  G4ThreeVector((dist_9+(thick_9/2))*(TMath::Sin(thetaP_9*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_9+(thick_9/2))*(TMath::Sin(thetaP_9*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_9+(thick_9/2))*(TMath::Cos(thetaP_9*degree))*mm),
	  logicalDetector9[NB],
	  CsI9name[NB],
	  logicWorld,
	  false,
	  CsI9num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 6 Si---------------------------------------34----
  Si6[NB] = new G4Trap(Si6name[NB],
	  thickS_6/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  distS_6*TMath::Tan(thetaS_6*degree)*mm,	//G4double pDy1,
	  hXS1_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hXS2_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (distS_6+thickS_6)*TMath::Tan(thetaS_6*degree)*mm,	//G4double pDy2,
	  hXS3_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hXS4_6*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_Si_6;
  Rot_Si_6 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_Si_6 -> rotateZ(TMath::ATan(TMath::Sin(thetaSP_6*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_6*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_6 -> rotateX(1*thetaSP_6*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_Si_6 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaSP_6*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_6*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_6 -> rotateX(thetaSP_6*degree);}

  logicalDetectorS6[NB] = new G4LogicalVolume(Si6[NB] ,scintillator_mat ,Si6name[NB]);
  PhysicalDetectorS6[NB] = new G4PVPlacement(Rot_Si_6,
	  G4ThreeVector((distS_6+(thickS_6/2))*(TMath::Sin(thetaSP_6*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_6+(thickS_6/2))*(TMath::Sin(thetaSP_6*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_6+(thickS_6/2))*(TMath::Cos(thetaSP_6*degree))*mm),
	  logicalDetectorS6[NB],
	  Si6name[NB],
	  logicWorld,
	  false,
	  Si6num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 7 Si---------------------------------------34----
  Si7[NB] = new G4Trap(Si7name[NB],
	  thickS_7/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  distS_7*TMath::Tan(thetaS_7*degree)*mm,	//G4double pDy1,
	  hXS1_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hXS2_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (distS_7+thickS_7)*TMath::Tan(thetaS_7*degree)*mm,	//G4double pDy2,
	  hXS3_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hXS4_7*TMath::Tan(370.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_Si_7;
  Rot_Si_7 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_Si_7 -> rotateZ(TMath::ATan(TMath::Sin(thetaSP_7*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_7*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_7 -> rotateX(1*thetaSP_7*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_Si_7 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaSP_7*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_7*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_7 -> rotateX(thetaSP_7*degree);}

  logicalDetectorS7[NB] = new G4LogicalVolume(Si7[NB] ,scintillator_mat ,Si7name[NB]);
  PhysicalDetectorS7[NB] = new G4PVPlacement(Rot_Si_7,
	  G4ThreeVector((distS_7+(thickS_7/2))*(TMath::Sin(thetaSP_7*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_7+(thickS_7/2))*(TMath::Sin(thetaSP_7*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_7+(thickS_7/2))*(TMath::Cos(thetaSP_7*degree))*mm),
	  logicalDetectorS7[NB],
	  Si7name[NB],
	  logicWorld,
	  false,
	  Si7num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 8 Si---------------------------------------34----
  Si8[NB] = new G4Trap(Si8name[NB],
	  thickS_8/2*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  distS_8*TMath::Tan(thetaS_8*degree)*mm,	//G4double pDy1,
	  hXS1_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hXS2_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (distS_8+thickS_8)*TMath::Tan(thetaS_8*degree)*mm,	//G4double pDy2,
	  hXS3_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hXS4_8*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)
  G4RotationMatrix* Rot_Si_8;
  Rot_Si_8 = new G4RotationMatrix;
  
  if(0 <= NB && NB <=11){
  Rot_Si_8 -> rotateZ(TMath::ATan(TMath::Sin(thetaSP_8*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_8*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_8 -> rotateX(1*thetaSP_8*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_Si_8 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaSP_8*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_8*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_8 -> rotateX(thetaSP_8*degree);}

  logicalDetectorS8[NB] = new G4LogicalVolume(Si8[NB] ,scintillator_mat ,Si8name[NB]);
  PhysicalDetectorS8[NB] = new G4PVPlacement(Rot_Si_8,
	  G4ThreeVector((distS_8+(thickS_8/2))*(TMath::Sin(thetaSP_8*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_8+(thickS_8/2))*(TMath::Sin(thetaSP_8*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_8+(thickS_8/2))*(TMath::Cos(thetaSP_8*degree))*mm),
	  logicalDetectorS8[NB],
	  Si8name[NB],
	  logicWorld,
	  false,
	  Si8num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 9 Si---------------------------------------34----
  Si9[NB] = new G4Trap(Si9name[NB],
	  thickS_9/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  distS_9*TMath::Tan(thetaS_9*degree)*mm,	//G4double pDy1,
	  hXS1_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hXS2_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (distS_9+thickS_9)*TMath::Tan(thetaS_9*degree)*mm,	//G4double pDy2,
	  hXS3_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hXS4_9*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_Si_9;
  Rot_Si_9 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_Si_9 -> rotateZ(TMath::ATan(TMath::Sin(thetaSP_9*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_9*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_9 -> rotateX(1*thetaSP_9*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_Si_9 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaSP_9*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaSP_9*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_Si_9 -> rotateX(thetaSP_9*degree);}

  logicalDetectorS9[NB] = new G4LogicalVolume(Si9[NB] ,scintillator_mat ,Si9name[NB]);
  PhysicalDetectorS9[NB] = new G4PVPlacement(Rot_Si_9,
	  G4ThreeVector((distS_9+(thickS_9/2))*(TMath::Sin(thetaSP_9*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_9+(thickS_9/2))*(TMath::Sin(thetaSP_9*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (distS_9+(thickS_9/2))*(TMath::Cos(thetaSP_9*degree))*mm),
	  logicalDetectorS9[NB],
	  Si9name[NB],
	  logicWorld,
	  false,
	  Si9num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 10 chio---------------------------------------34----
  CsI10[NB] = new G4Trap(CsI10name[NB],
	  thick_10/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_10*TMath::Tan(theta_10*degree)*mm,	//G4double pDy1,
	  hX1_10*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_10*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_10+thick_10)*TMath::Tan(theta_10*degree)*mm,	//G4double pDy2,
	  hX3_10*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_10*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_10;
  Rot_CsI_10 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_CsI_10 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_10*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_10*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_10 -> rotateX(1*thetaP_10*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_10 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_10*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_10*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_10 -> rotateX(thetaP_10*degree);}

  logicalDetector10[NB] = new G4LogicalVolume(CsI10[NB] ,scintillator_mat ,CsI10name[NB]);
  PhysicalDetector10[NB] = new G4PVPlacement(Rot_CsI_10,
	  G4ThreeVector((dist_10+(thick_10/2))*(TMath::Sin(thetaP_10*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_10+(thick_10/2))*(TMath::Sin(thetaP_10*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_10+(thick_10/2))*(TMath::Cos(thetaP_10*degree))*mm),
	  logicalDetector10[NB],
	  CsI10name[NB],
	  logicWorld,
	  false,
	  CsI10num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 11 chio---------------------------------------34----
  CsI11[NB] = new G4Trap(CsI11name[NB],
	  thick_11/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_11*TMath::Tan(theta_11*degree)*mm,	//G4double pDy1,
	  hX1_11*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_11*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_11+thick_11)*TMath::Tan(theta_11*degree)*mm,	//G4double pDy2,
	  hX3_11*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_11*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_11;
  Rot_CsI_11 = new G4RotationMatrix;
  if(0 <= NB && NB <=11){
  Rot_CsI_11 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_11*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_11*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_11 -> rotateX(1*thetaP_11*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_11 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_11*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_11*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_11 -> rotateX(thetaP_11*degree);}

  logicalDetector11[NB] = new G4LogicalVolume(CsI11[NB] ,scintillator_mat ,CsI11name[NB]);
  PhysicalDetector11[NB] = new G4PVPlacement(Rot_CsI_11,
	  G4ThreeVector((dist_11+(thick_11/2))*(TMath::Sin(thetaP_11*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_11+(thick_11/2))*(TMath::Sin(thetaP_11*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_11+(thick_11/2))*(TMath::Cos(thetaP_11*degree))*mm),
	  logicalDetector11[NB],
	  CsI11name[NB],
	  logicWorld,
	  false,
	  CsI11num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 12 chio---------------------------------------34----
  CsI12[NB] = new G4Trap(CsI12name[NB],
	  thick_12/2*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_12*TMath::Tan(theta_12*degree)*mm,	//G4double pDy1,
	  hX1_12*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_12*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_12+thick_12)*TMath::Tan(theta_12*degree)*mm,	//G4double pDy2,
	  hX3_12*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_12*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)
  G4RotationMatrix* Rot_CsI_12;
  Rot_CsI_12 = new G4RotationMatrix;
  
  if(0 <= NB && NB <=11){
  Rot_CsI_12 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_12*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_12*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_12 -> rotateX(1*thetaP_12*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_12 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_12*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_12*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_12 -> rotateX(thetaP_12*degree);}

  logicalDetector12[NB] = new G4LogicalVolume(CsI12[NB] ,scintillator_mat ,CsI12name[NB]);
  PhysicalDetector12[NB] = new G4PVPlacement(Rot_CsI_12,
	  G4ThreeVector((dist_12+(thick_12/2))*(TMath::Sin(thetaP_12*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_12+(thick_12/2))*(TMath::Sin(thetaP_12*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_12+(thick_12/2))*(TMath::Cos(thetaP_12*degree))*mm),
	  logicalDetector12[NB],
	  CsI12name[NB],
	  logicWorld,
	  false,
	  CsI12num[NB],
	  true);
 //----------------------------------------------------------------------------------------------------------------------------
  //-----------------------Ring 13 chio---------------------------------------34----
  CsI13[NB] = new G4Trap(CsI13name[NB],
	  thick_13/2*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_13+thick_13)*TMath::Tan(theta_13*degree))*mm,	//G4double pDy2,
	  abs(hX3_13*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_13*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_13*TMath::Tan(theta_13*degree))*mm,	//G4double pDy1,
	  abs(hX1_13*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_13*TMath::Tan(360.0*degree/(2.0*N_24)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)
  G4RotationMatrix* Rot_CsI_13;
  Rot_CsI_13 = new G4RotationMatrix;
  
  if(0 <= NB && NB <=11){
  Rot_CsI_13 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_13*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_13*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_13 -> rotateX(180*degree+1*thetaP_13*degree);}
  
  if(12 <= NB && NB <=23){  
  Rot_CsI_13 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_13*degree)*TMath::Cos((360.0/(2*N_24)+(360.0/N_24)*NB)*degree)/(TMath::Sin(thetaP_13*degree)*TMath::Sin((360.0/(2*N_24)+(360.0/N_24)*NB)*degree))));
  Rot_CsI_13 -> rotateX(180*degree+thetaP_13*degree);}

  logicalDetector13[NB] = new G4LogicalVolume(CsI13[NB] ,scintillator_mat ,CsI13name[NB]);
  PhysicalDetector13[NB] = new G4PVPlacement(Rot_CsI_13,
	  G4ThreeVector((dist_13+(thick_13/2))*(TMath::Sin(thetaP_13*degree))*(TMath::Cos(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_13+(thick_13/2))*(TMath::Sin(thetaP_13*degree))*(TMath::Sin(360.0/(2*N_24)*degree+(360.0/N_24)*NB*degree))*mm, (dist_13+(thick_13/2))*(TMath::Cos(thetaP_13*degree))*mm),
	  logicalDetector13[NB],
	  CsI13name[NB],
	  logicWorld,
	  false,
	  CsI13num[NB],
	  true);
  //---------------------------------------------------------------------------------      
 logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());
  Si6VisAtt[NB] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  Si6VisAtt[NB] -> SetForceWireframe(false);
  logicalDetectorS6[NB] -> SetVisAttributes (Si6VisAtt[NB]);

  Si7VisAtt[NB] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  Si7VisAtt[NB] -> SetForceWireframe(false);
  logicalDetectorS7[NB] -> SetVisAttributes (Si7VisAtt[NB]);

  Si8VisAtt[NB] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  Si8VisAtt[NB] -> SetForceWireframe(false);
  logicalDetectorS8[NB] -> SetVisAttributes (Si8VisAtt[NB]);

  Si9VisAtt[NB] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  Si9VisAtt[NB] -> SetForceWireframe(false);
  logicalDetectorS9[NB] -> SetVisAttributes (Si9VisAtt[NB]);
 
 
  CsI6VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  CsI6VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector6[NB] -> SetVisAttributes (CsI6VisAtt[NB]);

  CsI7VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,1.0,0.5));
  CsI7VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector7[NB] -> SetVisAttributes (CsI7VisAtt[NB]);

  CsI8VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  CsI8VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector8[NB] -> SetVisAttributes (CsI8VisAtt[NB]);

  CsI9VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.5,1.0));
  CsI9VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector9[NB] -> SetVisAttributes (CsI9VisAtt[NB]);

  CsI10VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  CsI10VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector10[NB] -> SetVisAttributes (CsI10VisAtt[NB]);

  CsI11VisAtt[NB] = new G4VisAttributes(G4Colour(0.5,1.0,1.0));
  CsI11VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector11[NB] -> SetVisAttributes (CsI11VisAtt[NB]);

  CsI12VisAtt[NB] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  CsI12VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector12[NB] -> SetVisAttributes (CsI12VisAtt[NB]);

  CsI13VisAtt[NB] = new G4VisAttributes(G4Colour(0.2,1.0,1.0));
  CsI13VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector13[NB] -> SetVisAttributes (CsI13VisAtt[NB]);
  }
for(int NB = 0 ; NB < N_16 ; NB++)
{

  //-----------------------Ring 14 chio---------------------------------------34----
  CsI14[NB] = new G4Trap(CsI14name[NB],
	  thick_14/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_14+thick_14)*TMath::Tan(theta_14*degree))*mm,	//G4double pDy2,
	  abs(hX3_14*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_14*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_14*TMath::Tan(theta_14*degree))*mm,	//G4double pDy1,
	  abs(hX1_14*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_14*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_14;
  Rot_CsI_14 = new G4RotationMatrix;
  if(0 <= NB && NB <=7){
  Rot_CsI_14 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_14*degree)*TMath::Cos((360.0/(2*N_16)+(360.0/N_16)*NB)*degree)/(TMath::Sin(thetaP_14*degree)*TMath::Sin((360.0/(2*N_16)+(360.0/N_16)*NB)*degree))));
  Rot_CsI_14 -> rotateX(180*degree+1*thetaP_14*degree);}
  
  if(8 <= NB && NB <=15){  
  Rot_CsI_14 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_14*degree)*TMath::Cos((360.0/(2*N_16)+(360.0/N_16)*NB)*degree)/(TMath::Sin(thetaP_14*degree)*TMath::Sin((360.0/(2*N_16)+(360.0/N_16)*NB)*degree))));
  Rot_CsI_14 -> rotateX(180*degree+thetaP_14*degree);}

  logicalDetector14[NB] = new G4LogicalVolume(CsI14[NB] ,scintillator_mat ,CsI14name[NB]);
  PhysicalDetector14[NB] = new G4PVPlacement(Rot_CsI_14,
	  G4ThreeVector((dist_14+(thick_14/2))*(TMath::Sin(thetaP_14*degree))*(TMath::Cos(360.0/(2*N_16)*degree+(360.0/N_16)*NB*degree))*mm, (dist_14+(thick_14/2))*(TMath::Sin(thetaP_14*degree))*(TMath::Sin(360.0/(2*N_16)*degree+(360.0/N_16)*NB*degree))*mm, (dist_14+(thick_14/2))*(TMath::Cos(thetaP_14*degree))*mm),
	  logicalDetector14[NB],
	  CsI14name[NB],
	  logicWorld,
	  false,
	  CsI14num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 15 chio---------------------------------------34----
  CsI15[NB] = new G4Trap(CsI15name[NB],
	  thick_15/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_15+thick_15)*TMath::Tan(theta_15*degree))*mm,	//G4double pDy2,
	  abs(hX3_15*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_15*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_15*TMath::Tan(theta_15*degree))*mm,	//G4double pDy1,
	  abs(hX1_15*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_15*TMath::Tan(360.0*degree/(2.0*N_16)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_15;
  Rot_CsI_15 = new G4RotationMatrix;
  if(0 <= NB && NB <=7){
  Rot_CsI_15 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_15*degree)*TMath::Cos((360.0/(2*N_16)+(360.0/N_16)*NB)*degree)/(TMath::Sin(thetaP_15*degree)*TMath::Sin((360.0/(2*N_16)+(360.0/N_16)*NB)*degree))));
  Rot_CsI_15 -> rotateX(180*degree+1*thetaP_15*degree);}
  
  if(8 <= NB && NB <=15){  
  Rot_CsI_15 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_15*degree)*TMath::Cos((360.0/(2*N_16)+(360.0/N_16)*NB)*degree)/(TMath::Sin(thetaP_15*degree)*TMath::Sin((360.0/(2*N_16)+(360.0/N_16)*NB)*degree))));
  Rot_CsI_15 -> rotateX(180*degree+thetaP_15*degree);}

  logicalDetector15[NB] = new G4LogicalVolume(CsI15[NB] ,scintillator_mat ,CsI15name[NB]);
  PhysicalDetector15[NB] = new G4PVPlacement(Rot_CsI_15,
	  G4ThreeVector((dist_15+(thick_15/2))*(TMath::Sin(thetaP_15*degree))*(TMath::Cos(360.0/(2*N_16)*degree+(360.0/N_16)*NB*degree))*mm, (dist_15+(thick_15/2))*(TMath::Sin(thetaP_15*degree))*(TMath::Sin(360.0/(2*N_16)*degree+(360.0/N_16)*NB*degree))*mm, (dist_15+(thick_15/2))*(TMath::Cos(thetaP_15*degree))*mm),
	  logicalDetector15[NB],
	  CsI15name[NB],
	  logicWorld,
	  false,
	  CsI15num[NB],
	  true);
  //---------------------------------------------------------------------------------      
 logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());


  CsI14VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  CsI14VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector14[NB] -> SetVisAttributes (CsI14VisAtt[NB]);

  CsI15VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  CsI15VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector15[NB] -> SetVisAttributes (CsI15VisAtt[NB]);
}
for(int NB = 0 ; NB < N_8 ; NB++)
{

  //-----------------------Ring 16 chio---------------------------------------34----
  CsI16[NB] = new G4Trap(CsI16name[NB],
	  thick_16/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_16+thick_16)*TMath::Tan(theta_16*degree))*mm,	//G4double pDy2,
	  abs(hX3_16*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_16*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_16*TMath::Tan(theta_16*degree))*mm,	//G4double pDy1,
	  abs(hX1_16*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_16*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_16;
  Rot_CsI_16 = new G4RotationMatrix;
  if(0 <= NB && NB <=3){
  Rot_CsI_16 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_16*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_16*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_CsI_16 -> rotateX(180*degree+1*thetaP_16*degree);}
  
  if(4 <= NB && NB <=7){  
  Rot_CsI_16 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_16*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_16*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_CsI_16 -> rotateX(180*degree+thetaP_16*degree);}

  logicalDetector16[NB] = new G4LogicalVolume(CsI16[NB] ,scintillator_mat ,CsI16name[NB]);
  PhysicalDetector16[NB] = new G4PVPlacement(Rot_CsI_16,
	  G4ThreeVector((dist_16+(thick_16/2))*(TMath::Sin(thetaP_16*degree))*(TMath::Cos(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_16+(thick_16/2))*(TMath::Sin(thetaP_16*degree))*(TMath::Sin(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_16+(thick_16/2))*(TMath::Cos(thetaP_16*degree))*mm),
	  logicalDetector16[NB],
	  CsI16name[NB],
	  logicWorld,
	  false,
	  CsI16num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 17 chio---------------------------------------34----
  CsI17[NB] = new G4Trap(CsI17name[NB],
	  thick_17/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_17+thick_17)*TMath::Tan(theta_17*degree))*mm,	//G4double pDy2,
	  abs(hX3_17*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_17*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_17*TMath::Tan(theta_17*degree))*mm,	//G4double pDy1,
	  abs(hX1_17*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_17*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_CsI_17;
  Rot_CsI_17 = new G4RotationMatrix;
  if(0 <= NB && NB <=3){
  Rot_CsI_17 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_17*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_17*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_CsI_17 -> rotateX(180*degree+1*thetaP_17*degree);}
  
  if(4 <= NB && NB <=7){  
  Rot_CsI_17 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_17*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_17*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_CsI_17 -> rotateX(180*degree+thetaP_17*degree);}

  logicalDetector17[NB] = new G4LogicalVolume(CsI17[NB] ,scintillator_mat ,CsI17name[NB]);
  PhysicalDetector17[NB] = new G4PVPlacement(Rot_CsI_17,
	  G4ThreeVector((dist_17+(thick_17/2))*(TMath::Sin(thetaP_17*degree))*(TMath::Cos(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_17+(thick_17/2))*(TMath::Sin(thetaP_17*degree))*(TMath::Sin(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_17+(thick_17/2))*(TMath::Cos(thetaP_17*degree))*mm),
	  logicalDetector17[NB],
	  CsI17name[NB],
	  logicWorld,
	  false,
	  CsI17num[NB],
	  true);
  //---------------------------------------------------------------------------------      
 logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());


  CsI16VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.2,1.0));
  CsI16VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector16[NB] -> SetVisAttributes (CsI16VisAtt[NB]);

  CsI17VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,1.0,0.2));
  CsI17VisAtt[NB] -> SetForceWireframe(false);
  logicalDetector17[NB] -> SetVisAttributes (CsI17VisAtt[NB]);
}
// ------------- Surfaces --------------
//
// Csi crystal
/*/
  G4OpticalSurface* opTubeSurface = new G4OpticalSurface("TubeSurface");
  opTubeSurface->SetType(dielectric_metal);
  opTubeSurface->SetFinish(polished);
  opTubeSurface->SetModel(unified);

  new G4LogicalBorderSurface("TubeSurface",PhysicalDetector3,physWorld,opTubeSurface);

//G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>


//if (opticalSurface) opticalSurface->DumpInfo();

// Generate & Add Material Properties Table attached to the optical surfaces
//
const G4int num = 2;
G4double ephoton[num] = {2.034*eV, 4.136*eV};

//Optical alu rap Surface
G4double refractiveIndex[num] = {0.48, 1.55};
G4double specularLobe[num]    = {0.3, 0.3};
G4double specularSpike[num]   = {0.2, 0.2};
G4double backScatter[num]     = {0.2, 0.2};
G4double reflecitivity[num]   = {0.90, 0.92};

G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
myST1->AddProperty("REFLECITIVITY",         ephoton, reflecitivity,   num);
myST1->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
myST1->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);

G4cout << "Tube Surface G4MaterialPropertiesTable" << G4endl;
myST1->DumpTable();

opTubeSurface->SetMaterialPropertiesTable(myST1);

//runManager -> SetSensitiveDetector(pvp);*/

//Set detector's Color


return physWorld;
}

