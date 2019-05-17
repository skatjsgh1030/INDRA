#include "FAZIAPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "TMath.h"
#include <iostream>
#define restM 105 //MeV

using namespace CLHEP;
using namespace std;

FAZIAPrimaryGeneratorAction::FAZIAPrimaryGeneratorAction(): G4VUserPrimaryGeneratorAction()
{  
}


FAZIAPrimaryGeneratorAction::~FAZIAPrimaryGeneratorAction()
{
  delete fParticleGun;
}


void FAZIAPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  costheta = RandGauss::shoot(0,0.398);  
  random = G4UniformRand();  
  phirand = G4UniformRand();
  Parti_ene = G4UniformRand();

  ///////////////////////////////////////////////////////////////////////////
  //ifstream file ("bess-dflux_p0_p1_mp_mm.dat");
  //double p0, p1, mup, mum ;

  //double Mu_momtMax = 404;

  //G4double muPtot = 253.948;
  //G4double muMtot = 218.106;
  G4double muPtotflux = 492.289;
  G4double muMtotflux = 405.911;
  G4double momt_avg[65];
  //G4double kineticE[65];
  G4double momfor[65];
  G4double momlat[65];
  G4double fluxP[65];
  G4double normfluxP[65];
  G4double fluxM[65];
  G4double normfluxM[65];

  for(G4int i = 0 ; i < 65 ; i++)
  {	
	//file >> p0 >> p1 >> mup  >> mum;

	if(i==0)
	{
	  momfor[i] = p0[i];
	  momlat[i] = p1[i];

	  momt_avg[i] = (momfor[i] + momlat[i])/2;

	  fluxP[i] = mup[i]*momt_avg[i];
	  fluxM[i] = mum[i]*momt_avg[i];
	}
	else
	{
	  momfor[i] = p0[i];
	  momlat[i] = p1[i];

	  momt_avg[i] = (momfor[i] + momlat[i])/2;

	  fluxP[i] = fluxP[i-1] + mup[i]*momt_avg[i];
	  fluxM[i] = fluxM[i-1] + mum[i]*momt_avg[i];
	}

	momt_avg[i] = (momfor[i] + momlat[i])/2;

	normfluxP[i] = fluxP[i]/muPtotflux;
	normfluxM[i] = fluxM[i]/muMtotflux;	
  }
  //////////////////////////////////////////////////////////////////////////////
  G4int n_particle = 1;

  G4double Phi_rand = 2*pi*phirand;
  //G4double Theta_rand = (60*CLHEP::degree)*((TMath::Cos(90*CLHEP::degree*costheta))*(TMath::Cos(90*CLHEP::degree*costheta)) );
  G4double Theta_rand = (pi/2)*costheta;
  G4int r_value = 20;

  G4double x = r_value*(TMath::Sin(Theta_rand))*(TMath::Cos(Phi_rand))*m;
  G4double y = r_value*(TMath::Sin(Theta_rand))*(TMath::Sin(Phi_rand))*m;
  G4double z = r_value*(TMath::Cos(Theta_rand))*m;

  fParticleGun = new G4ParticleGun(n_particle);
  G4ThreeVector gunposition = G4ThreeVector(x,y,z);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  G4ParticleDefinition* particle; 

  for(G4int j = 0 ; j < 65 ; j++)
  {
	//cout << random << " and " << fluxM[j]/fluxP[j] << endl;
	if(0. <= random && random < muMtotflux/muPtotflux)
	//if(0. <= random && random < fluxM[j]/fluxP[j])
	{
	  if(normfluxP[j] <= Parti_ene && Parti_ene < normfluxP[j+1])
	  {
		particle = particleTable -> FindParticle(particleName = "mu+");		

		fParticleGun -> SetParticlePosition(gunposition);  
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(-x,-y,-z));
		fParticleGun -> SetParticleMomentum(momt_avg[j]*GeV);
		fParticleGun -> GeneratePrimaryVertex(anEvent);
		std::cout <<"event" << eventNum++ << G4endl;
	  }
	  else {
		continue;
	  }
	}

	if(muMtotflux/muPtotflux <= random && random <= 1. )
	//if(fluxM[j]/fluxP[j] <= random && random <= 1. )
	{
	  if(normfluxM[j] <= Parti_ene && Parti_ene < normfluxM[j+1])
	  {
		particle = particleTable -> FindParticle(particleName = "mu-");		

		fParticleGun -> SetParticlePosition(gunposition);  
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(-x,-y,-z));
		fParticleGun -> SetParticleMomentum(momt_avg[j]*GeV);
		fParticleGun -> GeneratePrimaryVertex(anEvent);
		std::cout <<"event" << eventNum++ << G4endl;
	  }
	  else {
		continue;
	  }

	}
  }
  //fParticleGun -> GeneratePrimaryVertex(anEvent);
  //std::cout <<"event" << eventNum++ << G4endl;
}
