#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <random>
#include "math.h"
#include "string.h"
#include <vector>
#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TUnixSystem.h"
#include "TRandom3.h"
#endif
#include <stdio.h>
#include <iostream>
using namespace std;

const int maxTrack = 20000;
int main(int argc, char *argv[])
{
    if(argc < 3)    return -1;

    int particle1 = -321;
    int particle2 = 3334;
    int particle3 = -321;
    int particle4 = -3334;
    cout << "particle 1: " << particle1 << endl;
    cout << "particle 2: " << particle2 << endl;
    cout << "particle 3: " << particle3 << endl;
    cout << "particle 4: " << particle4 << endl;

    int fileStart = atoi(argv[1]);
    int fileEnd = atoi(argv[2]);
    TChain *hadronTree = new TChain("hadronTree");

    for (int i = fileStart; i < fileEnd; i++)
    {
        TString fileName = "../../../data06/afterART7.7b_";
        fileName += i;
        fileName += ".root";
        hadronTree->Add(fileName);
    }


    TH1D *samek_pair1 = new TH1D("samek_pair1", "same event k 1;k^{*}", 40, 0.0, 4); // particle1 particle2
    TH1D *samek_pair2 = new TH1D("samek_pair2", "same event k 2;k^{*}", 40, 0.0, 4); // particle1 particle2
    TH1D *samept_pair1 = new TH1D("samept_pair1", "same pt 1;#Deltap_{T}",40,0,4);
    TH1D *samept_pair2 = new TH1D("samept_pair2", "same pt 2;#Deltap_{T}",40,0,4);
    TH1D *samey_pair1 = new TH1D("samey_pair1", "same event y 1;#Deltay", 80, -4, 4); // particle1 particle2
    TH1D *samey_pair2 = new TH1D("samey_pair2", "same event y 2;#Deltay", 80, -4, 4); // particle1 particle2
    TH1D *samephi_pair1 = new TH1D("samephi_pair1", "same event phi 1;#Delta#phi", 63,0,6.3); // particle1 particle2
    TH1D *samephi_pair2 = new TH1D("samephi_pair2", "same event phi 2;#Delta#phi", 63,0,6.3); // particle1 particle2

    samek_pair1->Sumw2();
    samek_pair2->Sumw2();
    samept_pair1->Sumw2();
    samept_pair2->Sumw2();
    samey_pair1->Sumw2();
    samey_pair2->Sumw2();
    samephi_pair1->Sumw2();
    samephi_pair2->Sumw2();


    int mult, id[maxTrack];
    float x[maxTrack], y[maxTrack], z[maxTrack], px[maxTrack], py[maxTrack], pz[maxTrack], mass[maxTrack], t[maxTrack];
    float energy[maxTrack],b;

    hadronTree->SetBranchAddress("nMultiplicityTree", &mult);
    hadronTree->SetBranchAddress("b", &b);
    hadronTree->SetBranchAddress("id", &id);
    hadronTree->SetBranchAddress("px", &px);
    hadronTree->SetBranchAddress("py", &py);
    hadronTree->SetBranchAddress("pz", &pz);
    hadronTree->SetBranchAddress("energy", &energy);

    // total events
    const int nevents = hadronTree->GetEntries();
    cout << "nevent: " << nevents << endl;
    //-----------------------------------same_pt event-------------------------------------------------------------------
    int same_event_pair1 = 0;
    int same_event_pair2 = 0;


    for (Int_t nevent = 0; nevent < nevents; nevent++)
    {
        hadronTree->GetEntry(nevent);

        if(b > 3.4) continue;//0-5%
        int omega_num=0;
        int omegabar_num=0;
        for (Int_t j = i + 1; j < mult; j++)
        {
            if(id[j]==3334 ) omega_num++;
            if(id[j]==-3334) omegabar_num++;
        }
        if(omega_num ==0 && omegabar_num ==0)   continue;
        if(omega_num > 1)   continue; 

        bool p2_eve=false;
        bool p4_eve=false;

        for (Int_t i = 0; i < mult - 1; i++)
        {
            
            TLorentzVector p1;
            p1.SetPxPyPzE(px[i], py[i], pz[i], energy[i]);
            if( fabs(p1.PseudoRapidity()) >1 ) continue;

            for (Int_t j = i + 1; j < mult; j++)
            {
                if (id[i] == id[j])     continue;
                TLorentzVector p2;
                p2.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                if( fabs(p2.PseudoRapidity())>1 ) continue;
                // particle 1 and 2 
                if ((id[i] == particle1 && id[j] == particle2) || (id[i] == particle2 && id[j] == particle1))
                {

                    p1.SetPxPyPzE(px[i], py[i], pz[i], energy[i]);
                    p2.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                    TLorentzVector p3 = p1 + p2;
                    TLorentzVector p4 = p1 - p2;
                    p1.Boost(-p3.BoostVector());
                    p2.Boost(-p3.BoostVector());
                    p4.Boost(-p3.BoostVector());
                  
                    float delk = 0.5*(p1 - p2).Rho();
                    //float deltapt = (p1-p2).Perp();
                    float deltapt = p4.Perp();
                    float deltay = p1.Rapidity() - p2.Rapidity();
                    float phi = (p1-p2).Phi();
                    samept_pair1->Fill(deltapt, 1.0/0.1);
                    samek_pair1->Fill(delk, 1.0 / 0.1);
                    samey_pair1->Fill(deltay,1.0/0.1);
                    samephi_pair1->Fill(phi,1.0/0.1);
                    p2_eve=true;

                }
                //pariticle 3 & 4 
                if ((id[i] == particle3 && id[j] == particle4) || (id[i] == particle4 && id[j] == particle3))
                {

                    p1.SetPxPyPzE(px[i], py[i], pz[i], energy[i]);
                    p2.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                    TLorentzVector p3 = p1 + p2;
                    TLorentzVector p4 = p1 - p2;
                    p1.Boost(-p3.BoostVector());
                    p2.Boost(-p3.BoostVector());
                    p4.Boost(-p3.BoostVector());

                    float delk1 = 0.5*(p1 - p2).Rho();
                   // float deltapt = (p1-p2).Perp();
                    float deltapt = P4.Perp();
                    float deltay = p1.Rapidity() - p2.Rapidity();
                    float phi = p4.Phi();
                    samept_pair2->Fill(deltapt, 1.0/0.1);
                    samek_pair2->Fill(delk1, 1.0/0.1);
                    samey_pair2->Fill(deltay, 1.0/0.1);
                    samephi_pair2->Fill(phi, 1.0/0.1);
                    p4_eve=true;
                }
                
            }
           
        }
        if(p2_eve==true) same_event_pair1++;
        if(p4_eve==true) same_event_pair2++;
    }
    cout << "complete!" << endl;
    cout << "same pair 1 events: " << same_event_pair1 << endl;
    cout << "same pair 2 events: " << same_event_pair2 << endl;

//     samek_pair1->Scale(1.0/same_event_pair1);
//     samek_pair2->Scale(1.0/same_event_pair2);
//     samept_pair1->Scale(1.0/same_event_pair1);
//     samept_pair2->Scale(1.0/same_event_pair2);
//     samey_pair1->Scale(1.0/same_event_pair1);
//     samey_pair2->Scale(1.0/same_event_pair2);
//     samephi_pair1->Scale(1.0/same_event_pair1);
//     samephi_pair2->Scale(1.0/same_event_pair2);
    //k same canvas
    TCanvas *c1 = new TCanvas();
    samek_pair1->Draw("e");
    samek_pair1->SetLineColor(kRed);
    samek_pair2->Draw("same");
    samek_pair2->GetXaxis()->SetTitle("k^{*}");
    samek_pair2->GetYaxis()->SetTitle("count");
    samek_pair2->SetLineColor(kBlue);
    TLegend *leg = new TLegend();
    leg->AddEntry(samek_pair1,"k^{-}#Omega^{-}","l");
    leg->AddEntry(samek_pair2,"k^{-}#Omega^{+}","l");
    leg->Draw();
    c1->Draw();
    //Delta y in same canvas
    TCanvas *c2 = new TCanvas();
    samey_pair1->Draw("e");

    samey_pair1->SetLineColor(kRed);
    samey_pair2->Draw("same");

    samey_pair2->SetLineColor(kBlue);
    TLegend *leg2 = new TLegend();
    leg2->AddEntry(samey_pair1,"k^{-}#Omega^{-}","l");
    leg2->AddEntry(samey_pair2,"k^{-}#Omega^{+}","l");
    leg2->Draw();
    c2->Draw();

    //Delta pT in the same canvas.
    TCanvas *c3 = new TCanvas();
    samept_pair1->Draw("e");
    samept_pair1->SetLineColor(kRed);
    samept_pair2->Draw("same");

    samept_pair2->SetLineColor(kBlue);
    TLegend *leg3 = new TLegend();
    leg3->AddEntry(samept_pair1,"k^{-}#Omega^{-}","l");
    leg3->AddEntry(samept_pair2,"k^{-}#Omega^{+}","l");
    leg3->Draw();
    c3->Draw();


    //Delta phi in the same canvas.
    TCanvas *c4 = new TCanvas();
    samephi_pair1->Draw("e");
    samephi_pair1->SetLineColor(kRed);
    samephi_pair2->Draw("same");
    samephi_pair2->SetLineColor(kBlue);
    TLegend *leg4 = new TLegend();
    leg4->AddEntry(samephi_pair1,"k^{-}#Omega^{-}","l");
    leg4->AddEntry(samephi_pair2,"k^{-}#Omega^{+}","l");
    leg4->Draw();
    c4->Draw();



    //kstar deduction
    TH1D *deduction_kstar = new TH1D("deduction_kstar","k^{-}#bar{#Omega}^{+}-k^{-}#Omega^{-};k^{*};C(k^{*})",40,0,4);
    deduction_kstar->Sumw2();
    deduction_kstar->Add(samek_pair2,samek_pair1,1,-1);
    //Delta y deduction
    TH1D *deduction_y = new TH1D("deduction_y","k^{-}#bar{#Omega}^{+}-k^{-}#Omega^{-};#Deltay;C(#Deltay)",80,-4,4);
    deduction_y->Sumw2();
    deduction_y->Add(samey_pair2,samey_pair1,1,-1);

    //Delta pt deduction
    TH1D *deduction_pt = new TH1D("deduction_pt","k^{-}#bar{#Omega}^{+}-k^{-}#Omega^{-};#Deltap_{T};C(#Deltap_{T})",40,0,4);
    deduction_pt->Sumw2();
    deduction_pt->Add(samept_pair2,samept_pair1,1,-1);

    //phi
    TH1D *deduction_phi = new TH1D("deduction_phi","k^{-}#bar{#Omega}^{+}-k^{-}#Omega^{-};#Delta#phi",63,0,6.3);
    deduction_phi->Sumw2();
    deduction_phi->Add(samephi_pair2,samephi_pair1,1,-1);

    TString filename = "deduction_kmom_kmop_";
    filename += atoi(argv[3]);
    filename += ".root";

    TFile *result_corr = new TFile(filename, "RECREATE");
    c1->Write();
    c2->Write();
    c3->Write();    
    c4->Write();
    samept_pair1->Write();
    samept_pair2->Write();
    deduction_pt->Write();
    samek_pair1->Write();
    samek_pair2->Write();
    deduction_kstar->Write();
    samey_pair1->Write();
    samey_pair2->Write();
    deduction_y->Write();
    samephi_pair1->Write();
    samephi_pair2->Write();
    deduction_phi->Write();
    result_corr->Write();
    return 0;
}
