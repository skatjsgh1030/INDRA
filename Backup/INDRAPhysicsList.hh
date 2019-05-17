/// \file FAZIA/include/FAZIAPhysicsList.hh
/// \brief Definition of the FAZIAPhysicsList class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef FAZIAPhysicsList_h
#define FAZIAPhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class FAZIAPhysicsListMessenger;

//class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class FAZIAPhysicsList : public G4VUserPhysicsList
{
  public:

    FAZIAPhysicsList();
 //   virtual ~FAZIAPhysicsList();

  public:

    virtual void ConstructParticle();
    virtual void ConstructProcess();

    virtual void SetCuts();

    //these methods Construct physics processes and register them
    void ConstructDecay();
    void ConstructEM();
    void ConstructOp();

    //for the Messenger 
    void SetVerbose(G4int);
    //void SetNbOfPhotonsCerenkov(G4int);
 
  private:

    //FAZIAPhysicsListMessenger* fMessenger;

    static G4ThreadLocal G4int fVerboseLevel;
    static G4ThreadLocal G4int fMaxNumPhotonStep;

    //static G4ThreadLocal G4Cerenkov* fCerenkovProcess;
    static G4ThreadLocal G4Scintillation* fScintillationProcess;
    static G4ThreadLocal G4OpAbsorption* fAbsorptionProcess;
    static G4ThreadLocal G4OpRayleigh* fRayleighScatteringProcess;
    static G4ThreadLocal G4OpMieHG* fMieHGScatteringProcess;
    static G4ThreadLocal G4OpBoundaryProcess* fBoundaryProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* FAZIAPhysicsList_h */
