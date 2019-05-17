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
const int N_12 = 12;
const double gap = 0.5;
//choi crystal
//---------------------------choi67------------------
double thetaB_67 = 27.0;//degree
double thetaF_67 = 14.0;
double dist_67 = 250.0;//mm
double thick_67 = 50.0;
double theta_67 = (thetaB_67 - thetaF_67)/2;
double thetaP_67 = (thetaB_67 + thetaF_67)/2;
double hX1_67 = (dist_67/(TMath::Cos(theta_67*degree)))*TMath::Sin(thetaF_67*degree);
double hX2_67 = (dist_67/(TMath::Cos(theta_67*degree)))*TMath::Sin(thetaB_67*degree);
double hX3_67 = ((dist_67+thick_67)/(TMath::Cos(theta_67*degree)))*TMath::Sin(thetaF_67*degree);
double hX4_67 = ((dist_67+thick_67)/(TMath::Cos(theta_67*degree)))*TMath::Sin(thetaB_67*degree);
//----------------------------choi89------------------
double thetaB_89 = 45.0;//degree
double thetaF_89 = 27.0;
double dist_89 = 120.0;//mm
double thick_89 = 50.0;
double theta_89 = (thetaB_89 - thetaF_89)/2;
double thetaP_89 = (thetaB_89 + thetaF_89)/2;
double hX1_89 = (dist_89/(TMath::Cos(theta_89*degree)))*TMath::Sin(thetaF_89*degree);
double hX2_89 = (dist_89/(TMath::Cos(theta_89*degree)))*TMath::Sin(thetaB_89*degree);
double hX3_89 = ((dist_89+thick_89)/(TMath::Cos(theta_89*degree)))*TMath::Sin(thetaF_89*degree);
double hX4_89 = ((dist_89+thick_89)/(TMath::Cos(theta_89*degree)))*TMath::Sin(thetaB_89*degree);
//-----------------------------choi1011------------------------------
double thetaB_1011 = 70.0;//degree
double thetaF_1011 = 45.5;
double dist_1011 = 120.0;//mm
double thick_1011 = 50.0;
double theta_1011 = (thetaB_1011 - thetaF_1011)/2;
double thetaP_1011 = (thetaB_1011 + thetaF_1011)/2;
double hX1_1011 = (dist_1011/(TMath::Cos(theta_1011*degree)))*TMath::Sin(thetaF_1011*degree)*mm;
double hX2_1011 = (dist_1011/(TMath::Cos(theta_1011*degree)))*TMath::Sin(thetaB_1011*degree)*mm;
double hX3_1011 = ((dist_1011+thick_1011)/(TMath::Cos(theta_1011*degree)))*TMath::Sin(thetaF_1011*degree)*mm;
double hX4_1011 = ((dist_1011+thick_1011)/(TMath::Cos(theta_1011*degree)))*TMath::Sin(thetaB_1011*degree)*mm;
//-----------------------------choi12-----------------
double thetaB_12 = 88.0;//degree
double thetaF_12 = 70.3;
double dist_12 = 120.0;//mm
double thick_12 = 50.0;
double theta_12 = (thetaB_12 - thetaF_12)/2;
double thetaP_12 = (thetaB_12 + thetaF_12)/2;
double hX1_12 = (dist_12/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaF_12*degree)*mm;
double hX2_12 = (dist_12/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaB_12*degree)*mm;
double hX3_12 = ((dist_12+thick_12)/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaF_12*degree)*mm;
double hX4_12 = ((dist_12+thick_12)/(TMath::Cos(theta_12*degree)))*TMath::Sin(thetaB_12*degree)*mm;
//---------------------------choi13------------------
double thetaB_13 = 92.0;//degree
double thetaF_13 = 110.0;
double dist_13 = 120.0;//mm
double thick_13 = 50.0;
double theta_13 = (thetaB_13 - thetaF_13)/2;
double thetaP_13 = (thetaB_13 + thetaF_13)/2;
double hX1_13 = (dist_13/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaF_13*degree);
double hX2_13 = (dist_13/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaB_13*degree);
double hX3_13 = ((dist_13+thick_13)/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaF_13*degree);
double hX4_13 = ((dist_13+thick_13)/(TMath::Cos(theta_13*degree)))*TMath::Sin(thetaB_13*degree);
//----------------------------
//---------------------------choi1415------------------
double thetaB_1415 = 111.0;//degree
double thetaF_1415 = 142.0;
double dist_1415 = 120.0;//mm
double thick_1415 = 50.0;
double theta_1415 = (thetaB_1415 - thetaF_1415)/2;
double thetaP_1415 = (thetaB_1415 + thetaF_1415)/2;
double hX1_1415 = (dist_1415/(TMath::Cos(theta_1415*degree)))*TMath::Sin(thetaF_1415*degree);
double hX2_1415 = (dist_1415/(TMath::Cos(theta_1415*degree)))*TMath::Sin(thetaB_1415*degree);
double hX3_1415 = ((dist_1415+thick_1415)/(TMath::Cos(theta_1415*degree)))*TMath::Sin(thetaF_1415*degree);
double hX4_1415 = ((dist_1415+thick_1415)/(TMath::Cos(theta_1415*degree)))*TMath::Sin(thetaB_1415*degree);
//----------------------------
//---------------------------choi1617------------------
double thetaB_1617 = 143.0;//degree
double thetaF_1617 = 176.0;
double dist_1617 = 120.0;//mm
double thick_1617 = 50.0;
double theta_1617 = (thetaB_1617 - thetaF_1617)/2;
double thetaP_1617 = (thetaB_1617 + thetaF_1617)/2;
double hX1_1617 = (dist_1617/(TMath::Cos(theta_1617*degree)))*TMath::Sin(thetaF_1617*degree);
double hX2_1617 = (dist_1617/(TMath::Cos(theta_1617*degree)))*TMath::Sin(thetaB_1617*degree);
double hX3_1617 = ((dist_1617+thick_1617)/(TMath::Cos(theta_1617*degree)))*TMath::Sin(thetaF_1617*degree);
double hX4_1617 = ((dist_1617+thick_1617)/(TMath::Cos(theta_1617*degree)))*TMath::Sin(thetaB_1617*degree);
//----------------------------

G4Trap* choi67[N_12];
G4LogicalVolume* logicalDetectorCHOI67[N_12];
G4VPhysicalVolume* PhysicalDetectorCHOI67[N_12];
G4VisAttributes* choi67VisAtt[N_12];

G4Trap* choi89[N_12];
G4LogicalVolume* logicalDetectorCHOI89[N_12];
G4VPhysicalVolume* PhysicalDetectorCHOI89[N_12];
G4VisAttributes* choi89VisAtt[N_12];

G4Trap* choi1011[N_12];
G4LogicalVolume* logicalDetectorCHOI1011[N_12];
G4VPhysicalVolume* PhysicalDetectorCHOI1011[N_12];
G4VisAttributes* choi1011VisAtt[N_12];

G4Trap* choi12[N_12];
G4LogicalVolume* logicalDetectorCHOI12[N_12];
G4VPhysicalVolume* PhysicalDetectorCHOI12[N_12];
G4VisAttributes* choi12VisAtt[N_12];

G4Trap* choi13[N_8];
G4LogicalVolume* logicalDetectorCHOI13[N_8];
G4VPhysicalVolume* PhysicalDetectorCHOI13[N_8];
G4VisAttributes* choi13VisAtt[N_8];

G4Trap* choi1415[N_8];
G4LogicalVolume* logicalDetectorCHOI1415[N_8];
G4VPhysicalVolume* PhysicalDetectorCHOI1415[N_8];
G4VisAttributes* choi1415VisAtt[N_8];

G4Trap* choi1617[N_8];
G4LogicalVolume* logicalDetectorCHOI1617[N_8];
G4VPhysicalVolume* PhysicalDetectorCHOI1617[N_8];
G4VisAttributes* choi1617VisAtt[N_8];

int choi67num[12] = {0,1,2,3,4,5,6,7,8,9,10,11};
string choi67name[12] = {"6700","6701","6702","6703","6704","6705","6706","6707","6708","6709","6710","6711"};
int choi89num[12] = {12,13,14,15,16,17,18,19,20,21,22,23};
string choi89name[12] = {"8900","8901","8902","8903","8904","8905","8906","8907","8908","8909","8910","8911"};
int choi1011num[12] = {24,25,26,27,28,29,30,31,32,33,34,35};
string choi1011name[12] = {"101100","101101","101102","101103","101104","101105","101106","101107","101108","101109","101110","101111"};
int choi12num[12] = {36,37,38,39,40,41,42,43,44,45,46,47};
string choi12name[12] = {"1200","1201","1202","1203","1204","1205","1206","1207","1208","1209","1210","1211"};
int choi13num[8] = {36,37,38,39,40,41,42,43};
string choi13name[8] = {"1200","1201","1202","1203","1204","1205","1206","1207"};
int choi1415num[8] = {36,37,38,39,40,41,42,43};
string choi1415name[8] = {"1200","1201","1202","1203","1204","1205","1206","1201"};
int choi1617num[8] = {36,37,38,39,40,41,42,43};
string choi1617name[8] = {"1200","1201","1202","1203","1204","1205","1206","1201"};

for(int NB = 0 ; NB < N_12 ; NB++)
{
  //-----------------------Ring 6,7 chio---------------------------------------34----
  choi67[NB] = new G4Trap(choi67name[NB],
	  thick_67/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_67*TMath::Tan(theta_67*degree)*mm,	//G4double pDy1,
	  hX1_67*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_67*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_67+thick_67)*TMath::Tan(theta_67*degree)*mm,	//G4double pDy2,
	  hX3_67*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_67*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_choi_67;
  Rot_choi_67 = new G4RotationMatrix;
  if(0 <= NB && NB <=5){
  Rot_choi_67 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_67*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_67*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_67 -> rotateX(1*thetaP_67*degree);}
  
  if(6 <= NB && NB <=11){  
  Rot_choi_67 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_67*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_67*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_67 -> rotateX(thetaP_67*degree);}

  logicalDetectorCHOI67[NB] = new G4LogicalVolume(choi67[NB] ,scintillator_mat ,choi67name[NB]);
  PhysicalDetectorCHOI67[NB] = new G4PVPlacement(Rot_choi_67,
	  G4ThreeVector((dist_67+(thick_67/2))*(TMath::Sin(thetaP_67*degree))*(TMath::Cos(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_67+(thick_67/2))*(TMath::Sin(thetaP_67*degree))*(TMath::Sin(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_67+(thick_67/2))*(TMath::Cos(thetaP_67*degree))*mm),
	  logicalDetectorCHOI67[NB],
	  choi67name[NB],
	  logicWorld,
	  false,
	  choi67num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 8,9 chio---------------------------------------34----
  choi89[NB] = new G4Trap(choi89name[NB],
	  thick_89/2*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_89*TMath::Tan(theta_89*degree)*mm,	//G4double pDy1,
	  hX1_89*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_89*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_89+thick_89)*TMath::Tan(theta_89*degree)*mm,	//G4double pDy2,
	  hX3_89*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_89*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)
  G4RotationMatrix* Rot_choi_89;
  Rot_choi_89 = new G4RotationMatrix;
  
  if(0 <= NB && NB <=5){
  Rot_choi_89 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_89*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_89*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_89 -> rotateX(1*thetaP_89*degree);}
  
  if(6 <= NB && NB <=11){  
  Rot_choi_89 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_89*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_89*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_89 -> rotateX(thetaP_89*degree);}

  logicalDetectorCHOI89[NB] = new G4LogicalVolume(choi89[NB] ,scintillator_mat ,choi89name[NB]);
  PhysicalDetectorCHOI89[NB] = new G4PVPlacement(Rot_choi_89,
	  G4ThreeVector((dist_89+(thick_89/2))*(TMath::Sin(thetaP_89*degree))*(TMath::Cos(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_89+(thick_89/2))*(TMath::Sin(thetaP_89*degree))*(TMath::Sin(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_89+(thick_89/2))*(TMath::Cos(thetaP_89*degree))*mm),
	  logicalDetectorCHOI89[NB],
	  choi89name[NB],
	  logicWorld,
	  false,
	  choi89num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 10,11 chio---------------------------------------34----
  choi1011[NB] = new G4Trap(choi1011name[NB],
	  thick_1011/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_1011*TMath::Tan(theta_1011*degree)*mm,	//G4double pDy1,
	  hX1_1011*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_1011*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_1011+thick_1011)*TMath::Tan(theta_1011*degree)*mm,	//G4double pDy2,
	  hX3_1011*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_1011*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_choi_1011;
  Rot_choi_1011 = new G4RotationMatrix;
  if(0 <= NB && NB <=5){
  Rot_choi_1011 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_1011*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_1011*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_1011 -> rotateX(1*thetaP_1011*degree);}
  
  if(6 <= NB && NB <=11){  
  Rot_choi_1011 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_1011*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_1011*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_1011 -> rotateX(thetaP_1011*degree);}

  logicalDetectorCHOI1011[NB] = new G4LogicalVolume(choi1011[NB] ,scintillator_mat ,choi1011name[NB]);
  PhysicalDetectorCHOI1011[NB] = new G4PVPlacement(Rot_choi_1011,
	  G4ThreeVector((dist_1011+(thick_1011/2))*(TMath::Sin(thetaP_1011*degree))*(TMath::Cos(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_1011+(thick_1011/2))*(TMath::Sin(thetaP_1011*degree))*(TMath::Sin(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_1011+(thick_1011/2))*(TMath::Cos(thetaP_1011*degree))*mm),
	  logicalDetectorCHOI1011[NB],
	  choi1011name[NB],
	  logicWorld,
	  false,
	  choi1011num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 12 chio---------------------------------------34----
  choi12[NB] = new G4Trap(choi12name[NB],
	  thick_12/2*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  dist_12*TMath::Tan(theta_12*degree)*mm,	//G4double pDy1,
	  hX1_12*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double  pDx1,  
	  hX2_12*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx2,
	  0*degree,		//G4double  pAlp1, 
	  (dist_12+thick_12)*TMath::Tan(theta_12*degree)*mm,	//G4double pDy2,
	  hX3_12*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx3,  
	  hX4_12*TMath::Tan(360.0*degree/(2.0*N_12)-gap*degree)*mm,	//G4double pDx4,
	  0*degree);	//G4double pAlp2)
  G4RotationMatrix* Rot_choi_12;
  Rot_choi_12 = new G4RotationMatrix;
  
  if(0 <= NB && NB <=5){
  Rot_choi_12 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_12*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_12*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_12 -> rotateX(1*thetaP_12*degree);}
  
  if(6 <= NB && NB <=11){  
  Rot_choi_12 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_12*degree)*TMath::Cos((360.0/(2*N_12)+(360.0/N_12)*NB)*degree)/(TMath::Sin(thetaP_12*degree)*TMath::Sin((360.0/(2*N_12)+(360.0/N_12)*NB)*degree))));
  Rot_choi_12 -> rotateX(thetaP_12*degree);}

  logicalDetectorCHOI12[NB] = new G4LogicalVolume(choi12[NB] ,scintillator_mat ,choi12name[NB]);
  PhysicalDetectorCHOI12[NB] = new G4PVPlacement(Rot_choi_12,
	  G4ThreeVector((dist_12+(thick_12/2))*(TMath::Sin(thetaP_12*degree))*(TMath::Cos(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_12+(thick_12/2))*(TMath::Sin(thetaP_12*degree))*(TMath::Sin(360.0/(2*N_12)*degree+(360.0/N_12)*NB*degree))*mm, (dist_12+(thick_12/2))*(TMath::Cos(thetaP_12*degree))*mm),
	  logicalDetectorCHOI12[NB],
	  choi12name[NB],
	  logicWorld,
	  false,
	  choi12num[NB],
	  true);
 //----------------------------------------------------------------------------------------------------------------------------
  choi67VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  choi67VisAtt[NB] -> SetForceWireframe(true);
  logicalDetectorCHOI67[NB] -> SetVisAttributes (choi67VisAtt[NB]);

  choi89VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  choi89VisAtt[NB] -> SetForceWireframe(true);
  logicalDetectorCHOI89[NB] -> SetVisAttributes (choi89VisAtt[NB]);

  choi1011VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  choi1011VisAtt[NB] -> SetForceWireframe(true);
  logicalDetectorCHOI1011[NB] -> SetVisAttributes (choi1011VisAtt[NB]);

  choi12VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  choi12VisAtt[NB] -> SetForceWireframe(true);
  logicalDetectorCHOI12[NB] -> SetVisAttributes (choi12VisAtt[NB]);
  }
for(int NB = 0 ; NB < N_8 ; NB++)
{

  //-----------------------Ring 13 chio---------------------------------------34----
  choi13[NB] = new G4Trap(choi13name[NB],
	  thick_13/2*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_13+thick_13)*TMath::Tan(theta_13*degree))*mm,	//G4double pDy2,
	  abs(hX3_13*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_13*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_13*TMath::Tan(theta_13*degree))*mm,	//G4double pDy1,
	  abs(hX1_13*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_13*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)
  G4RotationMatrix* Rot_choi_13;
  Rot_choi_13 = new G4RotationMatrix;
  
  if(0 <= NB && NB <=3){
  Rot_choi_13 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_13*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_13*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_choi_13 -> rotateX(180*degree+1*thetaP_13*degree);}
  
  if(4 <= NB && NB <=7){  
  Rot_choi_13 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_13*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_13*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_choi_13 -> rotateX(180*degree+thetaP_13*degree);}

  logicalDetectorCHOI13[NB] = new G4LogicalVolume(choi13[NB] ,scintillator_mat ,choi13name[NB]);
  PhysicalDetectorCHOI13[NB] = new G4PVPlacement(Rot_choi_13,
	  G4ThreeVector((dist_13+(thick_13/2))*(TMath::Sin(thetaP_13*degree))*(TMath::Cos(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_13+(thick_13/2))*(TMath::Sin(thetaP_13*degree))*(TMath::Sin(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_13+(thick_13/2))*(TMath::Cos(thetaP_13*degree))*mm),
	  logicalDetectorCHOI13[NB],
	  choi13name[NB],
	  logicWorld,
	  false,
	  choi13num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 14,15 chio---------------------------------------34----
  choi1415[NB] = new G4Trap(choi1415name[NB],
	  thick_1415/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_1415+thick_1415)*TMath::Tan(theta_1415*degree))*mm,	//G4double pDy2,
	  abs(hX3_1415*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_1415*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_1415*TMath::Tan(theta_1415*degree))*mm,	//G4double pDy1,
	  abs(hX1_1415*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_1415*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_choi_1415;
  Rot_choi_1415 = new G4RotationMatrix;
  if(0 <= NB && NB <=3){
  Rot_choi_1415 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_1415*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_1415*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_choi_1415 -> rotateX(180*degree+1*thetaP_1415*degree);}
  
  if(4 <= NB && NB <=7){  
  Rot_choi_1415 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_1415*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_1415*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_choi_1415 -> rotateX(180*degree+thetaP_1415*degree);}

  logicalDetectorCHOI1415[NB] = new G4LogicalVolume(choi1415[NB] ,scintillator_mat ,choi1415name[NB]);
  PhysicalDetectorCHOI1415[NB] = new G4PVPlacement(Rot_choi_1415,
	  G4ThreeVector((dist_1415+(thick_1415/2))*(TMath::Sin(thetaP_1415*degree))*(TMath::Cos(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_1415+(thick_1415/2))*(TMath::Sin(thetaP_1415*degree))*(TMath::Sin(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_1415+(thick_1415/2))*(TMath::Cos(thetaP_1415*degree))*mm),
	  logicalDetectorCHOI1415[NB],
	  choi1415name[NB],
	  logicWorld,
	  false,
	  choi1415num[NB],
	  true);
  //---------------------------------------------------------------------------------      
  //-----------------------Ring 16,17 chio---------------------------------------34----
  choi1617[NB] = new G4Trap(choi1617name[NB],
	  thick_1617/2.0*mm,	//G4double  pDz,   
	  0*degree,	//G4double pTheta,
	  0*degree,	//G4double  pPhi,  
	  abs((dist_1617+thick_1617)*TMath::Tan(theta_1617*degree))*mm,	//G4double pDy2,
	  abs(hX3_1617*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx3,  
	  abs(hX4_1617*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx4,
	  0*degree,		//G4double  pAlp1, 
	  abs(dist_1617*TMath::Tan(theta_1617*degree))*mm,	//G4double pDy1,
	  abs(hX1_1617*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double  pDx1,  
	  abs(hX2_1617*TMath::Tan(360.0*degree/(2.0*N_8)-gap*degree))*mm,	//G4double pDx2,
	  0*degree);	//G4double pAlp2)

  G4RotationMatrix* Rot_choi_1617;
  Rot_choi_1617 = new G4RotationMatrix;
  if(0 <= NB && NB <=3){
  Rot_choi_1617 -> rotateZ(TMath::ATan(TMath::Sin(thetaP_1617*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_1617*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_choi_1617 -> rotateX(180*degree+1*thetaP_1617*degree);}
  
  if(4 <= NB && NB <=7){  
  Rot_choi_1617 -> rotateZ(180*degree+TMath::ATan(TMath::Sin(thetaP_1617*degree)*TMath::Cos((360.0/(2*N_8)+(360.0/N_8)*NB)*degree)/(TMath::Sin(thetaP_1617*degree)*TMath::Sin((360.0/(2*N_8)+(360.0/N_8)*NB)*degree))));
  Rot_choi_1617 -> rotateX(180*degree+thetaP_1617*degree);}

  logicalDetectorCHOI1617[NB] = new G4LogicalVolume(choi1617[NB] ,scintillator_mat ,choi1617name[NB]);
  PhysicalDetectorCHOI1617[NB] = new G4PVPlacement(Rot_choi_1617,
	  G4ThreeVector((dist_1617+(thick_1617/2))*(TMath::Sin(thetaP_1617*degree))*(TMath::Cos(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_1617+(thick_1617/2))*(TMath::Sin(thetaP_1617*degree))*(TMath::Sin(360.0/(2*N_8)*degree+(360.0/N_8)*NB*degree))*mm, (dist_1617+(thick_1617/2))*(TMath::Cos(thetaP_1617*degree))*mm),
	  logicalDetectorCHOI1617[NB],
	  choi1617name[NB],
	  logicWorld,
	  false,
	  choi1617num[NB],
	  true);
  //---------------------------------------------------------------------------------      
 logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());
 logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());


  choi13VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  choi13VisAtt[NB] -> SetForceWireframe(true);
  logicalDetectorCHOI13[NB] -> SetVisAttributes (choi13VisAtt[NB]);

  choi1415VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  choi1415VisAtt[NB] -> SetForceWireframe(true);
  logicalDetectorCHOI1415[NB] -> SetVisAttributes (choi1415VisAtt[NB]);

  choi1617VisAtt[NB] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  choi1617VisAtt[NB] -> SetForceWireframe(true);
  logicalDetectorCHOI1617[NB] -> SetVisAttributes (choi1617VisAtt[NB]);
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

