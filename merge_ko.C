#include "TFile.h"
#include <iostream>
#include <stdio.h>
#include "TH1D.h"
using namespace std;
void merge_ko(int event_pair1, int event_pair2, int event_pair3, int event_pair4){
    TFile *file = new TFile("kmom_kmop_dong_default.root","READ");
 
    TH1D *samek_pair1  = (TH1D*)file->Get("samek_pair1");
    TH1D *samek_pair2  = (TH1D*)file->Get("samek_pair2");
    TH1D *samept_pair1 = (TH1D*)file->Get("samept_pair1");
    TH1D *samept_pair2 = (TH1D*)file->Get("samept_pair2");
    TH1D *samey_pair1  = (TH1D*)file->Get("samey_pair1");
    TH1D *samey_pair2  = (TH1D*)file->Get("samey_pair2");

    samek_pair1->Scale(1.0/event_pair1);
    samek_pair2->Scale(1.0/event_pair2);
    samept_pair1->Scale(1.0/event_pair1);
    samept_pair2->Scale(1.0/event_pair2);
    samey_pair1->Scale(1.0/event_pair1);
    samey_pair2->Scale(1.0/event_pair2);

    
    
    //kstar deduction
    TH1D *kmop_kmom_kstar = new TH1D("kmop_kmom_kstar","k^{-}#bar{#Omega}^{+}-k^{-}#Omega^{-};k^{*};count",40,0,4);
    kmop_kmom_kstar->Sumw2();
    kmop_kmom_kstar->Add(samek_pair2,samek_pair1,1,-1);
    //Delta y deduction
    TH1D *kmop_kmom_y = new TH1D("kmop_kmom_y","k^{-}#bar{#Omega}^{+}-k^{-}#Omega^{-};#Deltay;count",80,-4,4);
    kmop_kmom_y->Sumw2();
    kmop_kmom_y->Add(samey_pair2,samey_pair1,1,-1);
    //Delta pt deduction
    TH1D *kmop_kmom_pt = new TH1D("kmop_kmom_pt","k^{-}#bar{#Omega}^{+}-k^{-}#Omega^{-};#Deltap_{T};count",40,0,4);
    kmop_kmom_pt->Sumw2();
    kmop_kmom_pt->Add(samept_pair2,samept_pair1,1,-1);
//------------------------------------------------------------------------

    TFile *file1 = new TFile("kpom_kpop_dong_default.root","READ");
    TH1D *samek_pair3  = (TH1D*)file1->Get("samek_pair1");
    TH1D *samek_pair4  = (TH1D*)file1->Get("samek_pair2");
    TH1D *samept_pair3 = (TH1D*)file1->Get("samept_pair1");
    TH1D *samept_pair4 = (TH1D*)file1->Get("samept_pair2");
    TH1D *samey_pair3  = (TH1D*)file1->Get("samey_pair1");
    TH1D *samey_pair4  = (TH1D*)file1->Get("samey_pair2");


    samek_pair3->Scale(1.0/event_pair3);
    samek_pair4->Scale(1.0/event_pair4);
    samept_pair3->Scale(1.0/event_pair3);
    samept_pair4->Scale(1.0/event_pair4);
    samey_pair3->Scale(1.0/event_pair3);
    samey_pair4->Scale(1.0/event_pair4);

    //kstar deduction
    TH1D *kpom_kpop_kstar = new TH1D("kpom_kpop_kstar","K^{+}#Omega^{-}-K^{+}#bar{#Omega}^{+};k^{*};count",40,0,4);
    kpom_kpop_kstar->Sumw2();
    kpom_kpop_kstar->Add(samek_pair3, samek_pair4, 1, -1);
    //Delta y deduction
    TH1D *kpom_kpop_y = new TH1D("kpom_kpop_y","K^{+}#Omega^{-}-K^{+}#bar{#Omega}^{+};#Deltay;count",80,-4,4);
    kpom_kpop_y->Sumw2();
    kpom_kpop_y->Add(samey_pair3, samey_pair4, 1, -1);
    //Delta pt deduction
    TH1D *kpom_kpop_pt = new TH1D("kpom_kpop_pt","K^{+}#Omega^{-}-K^{+}#bar{#Omega}^{+};#Deltap_{T};count",40,0,4);
    kpom_kpop_pt->Sumw2();
    kpom_kpop_pt->Add(samept_pair3, samept_pair4, 1, -1);

    TFile *koCorr = new TFile("kocorr_norm.root","RECREATE");
    kmop_kmom_kstar->Write();
    kmop_kmom_y->Write();
    kmop_kmom_pt->Write();

    kpom_kpop_kstar->Write();
    kpom_kpop_y->Write();
    kpom_kpop_pt->Write();
    koCorr->Write();
}
