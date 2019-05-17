#ifndef INDRADETECTORCONSTRUCTION_HH
#define INDRADETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

//using std::stringstream;
//using namespace std;
//using namespace CLHEP;

class G4VPhysicalVolume;
class G4LogicalVolume;

class INDRADetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    INDRADetectorConstruction();
    virtual ~INDRADetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
  

};

#endif
