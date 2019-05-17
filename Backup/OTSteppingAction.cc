#include "OTSteppingAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4ThreeVector.hh"

OTSteppingAction::OTSteppingAction()
: G4UserSteppingAction()
{
}

OTSteppingAction::~OTSteppingAction()
{
}

void OTSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
  G4int volumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();
  G4double edep = step -> GetTotalEnergyDeposit();
  G4double x = step -> GetDeltaPosition().x();  	
  G4double y = step -> GetDeltaPosition().y();  	
  G4double z = step -> GetDeltaPosition().z();
//G4double moment = step -> GetDeltaMomentum();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> FillNtupleIColumn(0, eventID);
  analysisManager -> FillNtupleIColumn(1, volumeID);
  analysisManager -> FillNtupleDColumn(2, edep);
  analysisManager -> FillNtupleDColumn(3, x);
  analysisManager -> FillNtupleDColumn(4, y);
  analysisManager -> FillNtupleDColumn(5, z);
  //analysisManager -> FillNtupleDColumn(6, moment);
  analysisManager -> AddNtupleRow();
}
