#include "TH1D.h"
#include "TCanvas.h"
#define LigSpd 150000000 //[mm/s]
#include <iostream>
#include "TF1.h"
#include "TMath.h"
//#include "TGraph.h"

//Impulse Response
Double_t Impulse_Response(Double_t *t, Double_t *par)
{
		Double_t time = t[0];
		Double_t expo;
		if (time >= par[6])
		{   
				expo = par[7]*(par[0]*exp(-1*(time-par[6])/par[1])+par[2]*exp(-1*(time-par[6])/par[3])+par[4]*exp(-1*(time-par[6])/par[5]));
		}   

		else
		{   
				expo = 0*time;
		}   
		return expo;
}


void ConvolutionModi()
{
		TF1 *fit_R = new TF1("fit_R",Impulse_Response,0,512,8);
		TF1 *fit_L = new TF1("fit_L",Impulse_Response,0,512,8);
		//parameter of ImperseResponce Func
		fit_R ->SetParameter(0,1.51);
		fit_R ->SetParameter(1,24.81);
		fit_R ->SetParameter(2,-1.81);
		fit_R ->SetParameter(3,4.145);
		fit_R ->SetParameter(4,0.159);
		fit_R ->SetParameter(5,169);

		fit_L ->SetParameter(0,1.51);
		fit_L ->SetParameter(1,24.81);
		fit_L ->SetParameter(2,-1.81);
		fit_L ->SetParameter(3,4.145);
		fit_L ->SetParameter(4,0.159);
		fit_L ->SetParameter(5,169);
		//address of data File
		auto file = new TFile("/Users/namseonho/Documents/opentutorials_Geant4/opentutorials_Geant4-single/build/Single_Neu100_02.05_T.root");
		auto tree = (TTree *) file -> Get("data");

		TClonesArray *stepArray2 = nullptr;
		tree -> SetBranchAddress("step2",&stepArray2);
		// XR,XL_Dist = EnegyDeposition to position
		TH1D* XR_dist[50];
		TH1D* XL_dist[50];
		// After Convolution data
		TH1D* sum_Imp_ReR[50];
		TH1D* sum_Imp_ReL[50];

		TCanvas* graph[50];
		for (int i = 0; i < 50; i++)   // i = Event Number
		{
				XR_dist[i] = new TH1D (Form("XR_dist%d", i),"XR_Position To Edep",512,0,512);
				XL_dist[i] = new TH1D (Form("XL_dist%d", i),"XL_Position To Edep",512,0,512);
				sum_Imp_ReR[i] = new TH1D (Form("sum_Imp_ReR%d",i),"Sum_R_Edep",512,0,512); //sum_Imp_ReR[i] -> SetLineColor(kBlue);
				sum_Imp_ReL[i] = new TH1D (Form("sum_Imp_ReL%d",i),"Sum_L_Edep",512,0,512); //sum_Imp_ReL[i] -> SetLineColor(kRed);

				//To Remove NULL Canvases
				tree -> GetEntry(i);
				auto numSteps = step2 -> GetEntries();

				if(numSteps!=0)
				{
						graph[i] = new TCanvas(Form("graph%d", i),"",800,800);
						graph[i] -> Divide(2,2);
				}
		}
		double Dist_XR;
		double Dist_XL;


		double Dps_X;

		double Dps_XR;
		double Dps_XL;

		double ratio_X = 38401000;//512/(2000/(c/2)):x-axis ratio constant

		double EDep ;

		//auto numEvent = tree -> GetEntries();
		for (auto i = 0; i < 50; i++)
		{

				tree -> GetEntry(i);

				auto numSteps = stepArray2 -> GetEntries();


				/*double SumCon_XR;
				  double SumCon_XL;*/

				cout << "i : " << i <<endl;//check
				cout << "num steps : " << numSteps << endl;//check 

				//numSteps : final number of Step number which in a event
				for (auto iStep = 0; iStep < numSteps; iStep++)
				{
						double val_R[numSteps];
						double val_L[numSteps];
						auto step = (KBMCStep *) step -> At(iStep);

						EDep = step -> GetEdep();

						Dps_X = step -> GetX();


						Dist_XR = ((Dps_X + 1000)/LigSpd)*ratio_X;
						XR_dist[i] -> Fill(Dist_XR,EDep);

						Dist_XL = ((1000 - Dps_X)/LigSpd)*ratio_X;
						XL_dist[i] -> Fill(Dist_XL,EDep);
						//Set parameter of Impulse Response as Energy and deposit position
						fit_R ->SetParameter(6,Dist_XR);
						fit_R ->SetParameter(7,EDep);


						fit_L ->SetParameter(6,Dist_XL);
						fit_L ->SetParameter(7,EDep);
						//Convolution with Energy Peak and Impulse Response Func
						for(int j=0;j<512;j++)
						{

								val_R[iStep] = fit_R -> Eval(j);
								val_R[iStep] += val_R[iStep];

								val_L[iStep] = fit_L -> Eval(j);
								val_L[iStep] += val_L[iStep];

								if (val_R[iStep]<0)
								{
										val_R[iStep]=0;
										sum_Imp_ReR[i] -> Fill(j,val_R[iStep]/2);
								}
								else
								{
										sum_Imp_ReR[i] -> Fill(j,val_R[iStep]/2);
								}

								if (val_L[iStep]<0)
								{
										val_L[iStep]=0;
										sum_Imp_ReL[i] -> Fill(j,val_L[iStep]/2);
								}
								else
								{

										sum_Imp_ReL[i] -> Fill(j,val_L[iStep]/2);
								}
						}

				}

				if(numSteps != 0)
				{
						graph[i]->cd(1);
						XL_dist[i]->Draw();
						XL_dist[i] ->Sumw2(kFALSE);

						graph[i]->cd(2);
						XR_dist[i]->Draw();
						XR_dist[i] ->Sumw2(kFALSE);

						graph[i]->cd(4);
						// if(sum_Imp_ReR[i]->GetN() != 0)
						sum_Imp_ReR[i]->Draw();
						sum_Imp_ReR[i]->Sumw2(kFALSE);

						graph[i]->cd(3);
						//if(sum_Imp_ReL[i]->GetN() != 0)
						sum_Imp_ReL[i]->Draw();
						sum_Imp_ReL[i]->Sumw2(kFALSE);
				}
		}
}

