#ifndef FAZIAPRIMARYGENERATORACTION_HH
#define FAZIAPRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <iostream>

using namespace CLHEP;
using namespace std;


class FAZIAPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	FAZIAPrimaryGeneratorAction();    
	virtual ~FAZIAPrimaryGeneratorAction();

	// method from the base class
	virtual void GeneratePrimaries(G4Event*);         

	// method to access particle gun
	const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  private:
	G4ParticleGun*  fParticleGun; // pointer a to G4 gun class

	G4double random;
	G4double phirand;
	G4double costheta;
	G4double Parti_ene;
	G4double p0[65] = { 0.576,0.621,0.669,0.720,0.776,0.836,0.901,0.970,1.04 ,1.13 ,1.21 ,1.31 ,1.41 ,1.52 ,1.63 ,1.76 ,1.90 ,2.04 ,2.20 ,2.37 ,2.55 ,2.75 ,2.96 ,3.19 ,3.44 ,3.71 ,3.99 ,4.30 ,4.63 ,4.99 ,5.38 ,5.79 ,6.24 ,6.73 ,7.25 ,7.81 ,8.41 ,9.06 ,9.76 ,10.5 ,11.3 ,12.2 ,13.2 ,14.2 ,15.3 ,16.4 ,17.7 ,19.1 ,20.6 ,23.9 ,27.7 ,32.1 ,37.3 ,43.3 ,50.2 ,58.3 ,67.7 ,78.5 ,91.1 ,106. ,132. ,165. ,207. ,258. ,323.};
	G4double p1[65] = {0.621 ,0.669 ,0.720 ,0.776 ,0.836 ,0.901 ,0.970 ,1.04  ,1.13  ,1.21  ,1.31  ,1.41  ,1.52  ,1.63  ,1.76  ,1.90  ,2.04  ,2.20  ,2.37  ,2.55  ,2.75  ,2.96  ,3.19  ,3.44  ,3.71  ,3.99  ,4.30  ,4.63  ,4.99  ,5.38  ,5.79  ,6.24  ,6.73  ,7.25  ,7.81  ,8.41  ,9.06  ,9.76  ,10.5  ,11.3  ,12.2  ,13.2  ,14.2  ,15.3  ,16.4  ,17.7  ,19.1  ,20.6  ,23.9  ,27.7  ,32.1  ,37.3  ,43.3  ,50.2  ,58.3  ,67.7  ,78.5  ,91.1  ,106.  ,132.  ,165.  ,207.  ,258.  ,323.  ,404. };
	G4double mup[65] = { 13.4   ,13.3   ,13.2   ,12.7   ,12.6   ,12.3   ,12.1   ,11.6   ,11.0   ,10.7   ,10.2   ,9.78   ,9.38   ,8.72   ,8.59   ,7.85   ,7.41   ,7.03   ,6.38   ,6.01   ,5.45   ,5.02   ,4.62   ,4.17   ,3.79   ,3.37   ,3.02   ,2.74   ,2.41   ,2.11   ,1.87   ,1.63   ,1.44   ,1.21   ,1.06   ,.911   ,.807   ,.706   ,.567   ,.490   ,.438   ,.358   ,.307   ,.249   ,.213   ,.176   ,.149   ,.122   ,.0957  ,.0669  ,.0452  ,.0293  ,.0201  ,.0131  ,.00876 ,.00572 ,.00401 ,.00255 ,.00155 ,.000943,.000426,.000251,.000118,.0000558,.0000256};
	G4double mum[65] = {12.5 ,12.2 ,11.9 ,11.5 ,11.3 ,10.8 ,10.6 ,10.0 ,9.6 ,9.1 ,8.8 ,8.53 ,8.01 ,7.53 ,7.22 ,6.72 ,6.28 ,5.71 ,5.36 ,4.92 ,4.62 ,4.09 ,3.69 ,3.31 ,3.06 ,2.70 ,2.37 ,2.11 ,1.87 ,1.66 ,1.46 ,1.24 ,1.13 ,0.96 ,0.83 ,.693 ,.608 ,.515 ,.454 ,.385 ,.324 ,.274 ,.231 ,.197 ,.168,.140,.115,.096,.0736,.0504,.0348,.0226,.0155,.0104,.00691,.00454,.00288,.00193,.00117,.000695,.000365,.000147,.000075,.00000413,.00000245};
	/*///////////////////////////////////////////////////////////////////////////
	  double p0, p1, mup, mum; 


	  double momt_avg[65];
	  double momfor[65] = {   };
	  double momlat[65];
	  }*/

G4int eventNum = 0;
};

#endif

