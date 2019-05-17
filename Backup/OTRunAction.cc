#include "OTRunAction.hh"
#include "g4root.hh"

OTRunAction::OTRunAction()
: G4UserRunAction()
{
}

OTRunAction::~OTRunAction()
{
  delete G4AnalysisManager::Instance();

}

void OTRunAction::BeginOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> OpenFile("Single_Neu100_02.05_T.root");

  analysisManager -> CreateNtuple("step", "step");
  analysisManager -> CreateNtupleIColumn("volumeID");
  analysisManager -> CreateNtupleDColumn("edep");
  analysisManager -> CreateNtupleDColumn("x");
  analysisManager -> CreateNtupleDColumn("y");
  analysisManager -> CreateNtupleDColumn("z");
  analysisManager -> FinishNtuple();
}

void OTRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> Write();
  analysisManager -> CloseFile();
}
