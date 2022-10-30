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
#include "TUnixSystem.h"
#include "TRandom3.h"
#endif
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;
const int maxTrack = 20000;

int main()
{
    int particle1 = -321;
    int particle2 = 3334;
    int particle3 = -321;
    int particle4 = -3334;
    cout << "particle 1: " << particle1 << endl;
    cout << "particle 2: " << particle2 << endl;
    cout << "particle 3: " << particle3 << endl;
    cout << "particle 4: " << particle4 << endl;


    TChain *hadronTree = new TChain("hadronTree");
    for (int i = 1; i < 2385; i++)
    {
        TString fileName = "../data77/afterART7.7a_";
        fileName += i;
        fileName += ".root";
        hadronTree->Add(fileName);
    }
    for (int i = 1; i < 1328; i++)
    {
        TString fileName = "../data06/afterART7.7a_";
        fileName += i;
        fileName += ".root";
        hadronTree->Add(fileName);
    }

    TH1D *same_pt1 = new TH1D("same_pt1", "same_pt event", 40, 0.0, 4); // particle1 particle2
    TH1D *mix_pt1 = new TH1D("mix_pt1", "mix_pt event", 40, 0.0, 4);
    TH1D *same_k1 = new TH1D("same_k1", "same event kstar",40, 0.0, 4); // particle1 particle2
    TH1D *mix_k1 = new TH1D("mix_k1", "mix event kstar", 40, 0.0, 4);

    TH1D *same_pt2 = new TH1D("same_pt2", "same_pt event",40, 0.0, 4); // particle1 particle2
    TH1D *mix_pt2 = new TH1D("mix_pt2", "mix_pt event", 40, 0.0, 4);
    TH1D *same_k2 = new TH1D("same_k2", "same event kstar", 40, 0.0, 4); // particle1 particle2
    TH1D *mix_k2 = new TH1D("mix_k2", "mix event kstar", 40, 0.0, 4);

    same_pt1->Sumw2();
    mix_pt1->Sumw2();
    same_k1->Sumw2();
    mix_k1->Sumw2();
    same_pt2->Sumw2();
    mix_pt2->Sumw2();
    same_k2->Sumw2();
    mix_k2->Sumw2();


    int mult, id[maxTrack];
    float x[maxTrack], y[maxTrack], z[maxTrack], px[maxTrack], py[maxTrack], pz[maxTrack], mass[maxTrack], t[maxTrack];
    float energy[maxTrack], b;

    hadronTree->SetBranchAddress("nMultiplicityTree", &mult);
    hadronTree->SetBranchAddress("b", &b);
    hadronTree->SetBranchAddress("id", &id);
    hadronTree->SetBranchAddress("px", &px);
    hadronTree->SetBranchAddress("py", &py);
    hadronTree->SetBranchAddress("pz", &pz);
    hadronTree->SetBranchAddress("energy", &energy);
    
    // total events
    int nevents = hadronTree->GetEntries();
    cout << "nevent: " << nevents << endl;

    int same_pt_num1 = 0;
    int same_pt_num2 = 0;
    vector<int> pid1;
    vector<float> fpx1;
    vector<float> fpy1;
    vector<float> fpz1;
    vector<float> fenergy1;
    vector<int> pid2;
    vector<float> fpx2;
    vector<float> fpy2;
    vector<float> fpz2;
    vector<float> fenergy2;
    int event_pool1 = 0;
    int event_pool2 = 0;

    // int nch;
    // TFile *nchfile = new TFile("nch1.root","READ");
    // TTree *NchTree = (TTree*)nchfile->Get("NchTree");
    // NchTree->SetBranchAddress("nch",&nch);

    for (int k = 0; k < nevents; k++)
    {
        hadronTree->GetEntry(k);
        if (b > 3.4)  continue;
        // same event
        for (Int_t i = 0; i < mult - 1; i++)
        {
            TLorentzVector p1;
            p1.SetPxPyPzE(px[i], py[i], pz[i], energy[i]);
            if (fabs(p1.PseudoRapidity()) > 1.0)
                continue; // Eta -1,1
            for (Int_t j = i + 1; j < mult; j++)
            {
                if (id[i] == id[j])
                    continue;
                // particle 1 and 2
                if ((id[i] == particle1 && id[j] == particle2) || (id[i] == particle2 && id[j] == particle1))
                {
                    TLorentzVector p2;
                    p2.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                    if (fabs(p2.PseudoRapidity()) > 1.0)
                        continue; // eta -1,1
                    p1.SetPxPyPzE(px[i], py[i], pz[i], energy[i]);
                    TLorentzVector p3 = p1 + p2;
                    p1.Boost(-p3.BoostVector());
                    p2.Boost(-p3.BoostVector());
                    float delpt = (p1 - p2).Perp();
                    float delk = 0.5*(p1 - p2).Rho();
                    same_pt1->Fill(delpt, 1.0 / 0.1);
                    same_k1->Fill(delk, 1.0 / 0.1);
                    same_pt_num1++;
                }
                // particle 3 and 4
                if ((id[i] == particle3 && id[j] == particle4) || (id[i] == particle4 && id[j] == particle3))
                {
                    TLorentzVector p2;
                    p2.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                    if (fabs(p2.PseudoRapidity()) > 1.0)
                        continue; // Eta -1,1
                    p1.SetPxPyPzE(px[i], py[i], pz[i], energy[i]);
                    TLorentzVector p3 = p1 + p2;
                    p1.Boost(-p3.BoostVector());
                    p2.Boost(-p3.BoostVector());
                    float delpt = (p1 - p2).Perp();
                    float delk = 0.5*(p1 - p2).Rho();
                    same_pt2->Fill(delpt, 1.0 / 0.1);
                    same_k2->Fill(delk, 1.0 / 0.1);
                    same_pt_num2++;
                }
            }
        }

        /*---------------------event mixing--------------------------------------*/
        // ten events containing omega will be written into the container.
        bool hasOmega = false;
        bool hasAntiOmega = false;
        // check if this event has a omega/Anti Omega in Rap(-0.5,0.5)
        for (int j = 0; j < mult; j++)
        {
            if (id[j] != particle2)
                continue;
            TLorentzVector p1;
            p1.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
            if (fabs(p1.PseudoRapidity()) > 1.0)
                continue; // Eta -1,1
            hasOmega = true;
            break;
        }
        for (int j = 0; j < mult; j++)
        {
            if (id[j] != particle4)
                continue;
            TLorentzVector p1;
            p1.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
            if (fabs(p1.PseudoRapidity()) > 1.0)
                continue; // Eta -1,1
            hasAntiOmega = true;
            break;
        }

        // put omega and kaon into container.
        if (hasOmega == true)
        {
            for (int j = 0; j < mult; j++)
            {
                // has omega event
                if (id[j] != particle1 && id[j] != particle2)
                    continue;
                // cout << "This particle is omega" << endl;
                TLorentzVector particle;
                particle.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                if (fabs(particle.PseudoRapidity()) > 1.0)
                    continue;
                pid1.push_back(id[j]);
                fpx1.push_back(px[j]);
                fpy1.push_back(py[j]);
                fpz1.push_back(pz[j]);
                fenergy1.push_back(energy[j]);
            }
            event_pool1++;
        }
        // put omega and kaon into container.
        if (hasAntiOmega == true)
        {
            for (int j = 0; j < mult; j++)
            {
                // has omega event
                if (id[j] != particle3 && id[j] != particle4)
                    continue;
                // cout << "This particle is omega" << endl;
                TLorentzVector particle;
                particle.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                if (fabs(particle.PseudoRapidity()) > 1.0)
                    continue; // eta -1,1
                pid2.push_back(id[j]);
                fpx2.push_back(px[j]);
                fpy2.push_back(py[j]);
                fpz2.push_back(pz[j]);
                fenergy2.push_back(energy[j]);
            }
            event_pool2++;
        }

        // mix the first pair
        if (event_pool1 == 10)
        {
            // cout << "mixing particle number of pair 1: "<< pid1.size() <<endl;
            for (int m = 0; m < pid1.size() - 1; m++)
            {
                TLorentzVector p1;
                p1.SetPxPyPzE(fpx1[m], fpy1[m], fpz1[m], fenergy1[m]);
                for (int k = m + 1; k < pid1.size(); k++)
                {
                    if (pid1[m] == pid1[k])
                        continue;
                    TLorentzVector p2;
                    p2.SetPxPyPzE(fpx1[k], fpy1[k], fpz1[k], fenergy1[k]);
                    p1.SetPxPyPzE(fpx1[m], fpy1[m], fpz1[m], fenergy1[m]);
                    TLorentzVector p3 = p1 + p2;
                    p1.Boost(-p3.BoostVector());
                    p2.Boost(-p3.BoostVector());
                    float pt = (p1 - p2).Perp();
                    float dst = 0.5*(p1 - p2).Rho();
                    mix_pt1->Fill(pt, 1.0 / 0.1);
                    mix_k1->Fill(dst, 1.0 / 0.1);
                }
            }
            // reset the container.
            pid1.clear();
            fpx1.clear();
            fpy1.clear();
            fpz1.clear();
            fenergy1.clear();
            event_pool1 = 0;
        }
        // event mixing second pair
        if (event_pool2 == 10)
        {
            // cout << "mixing particle number of pair 2: "<< pid2.size() <<endl;
            for (int m = 0; m < pid2.size() - 1; m++)
            {
                TLorentzVector p1;
                p1.SetPxPyPzE(fpx2[m], fpy2[m], fpz2[m], fenergy2[m]);
                for (int k = m + 1; k < pid2.size(); k++)
                {
                    if (pid2[m] == pid2[k])
                        continue;
                    TLorentzVector p2;
                    p2.SetPxPyPzE(fpx2[k], fpy2[k], fpz2[k], fenergy2[k]);
                    p1.SetPxPyPzE(fpx2[m], fpy2[m], fpz2[m], fenergy2[m]);
                    TLorentzVector p3 = p1 + p2;
                    p1.Boost(-p3.BoostVector());
                    p2.Boost(-p3.BoostVector());
                    float pt = (p1 - p2).Perp();
                    float dst = 0.5*(p1 - p2).Rho();
                    mix_pt2->Fill(pt, 1.0 / 0.1);
                    mix_k2->Fill(dst, 1.0 / 0.1);
                }
            }
            // reset the container after mixing event.
            pid2.clear();
            fpx2.clear();
            fpy2.clear();
            fpz2.clear();
            fenergy2.clear();
            event_pool2 = 0;
        }
    }
    cout << "particle 1 and 2 in same event: " << same_pt_num1 << endl;
    cout << "particle 3 and 4 in same event: " << same_pt_num2 << endl;
    cout << "same and mixing event in same cycle complete!" << endl;
    // pt normalization kmom
    TF1 *constant_pt1 = new TF1("constant_pt1", "1", 0, 10);
    Int_t binNorm_pt1[2];
    binNorm_pt1[0] = same_pt1->FindBin(0.8);
    binNorm_pt1[1] = same_pt1->FindBin(1.3);
    Double_t factorN_pt1 = mix_pt1->Integral(binNorm_pt1[0], binNorm_pt1[1]) / same_pt1->Integral(binNorm_pt1[0], binNorm_pt1[1]);
    cout << "normalization factor: " << factorN_pt1 << endl;

    TCanvas *c1 = new TCanvas();
    same_pt1->Draw("e");
    same_pt1->GetXaxis()->SetTitle("#Deltap_{T}");
    same_pt1->GetYaxis()->SetTitle("count");
    c1->Draw();

    TCanvas *c2 = new TCanvas();
    mix_pt1->Draw("e");
    mix_pt1->GetXaxis()->SetTitle("#Deltap_{T}");
    mix_pt1->GetYaxis()->SetTitle("count");
    c2->Draw();

    TCanvas *c3 = new TCanvas();
    TH1D *norm_pt1 = (TH1D *)same_pt1->Clone();
    TH1D *mixClone_pt1 = (TH1D *)mix_pt1->Clone();

    norm_pt1->Sumw2();
    mixClone_pt1->Sumw2();

    mixClone_pt1->Divide(constant_pt1, factorN_pt1);
    norm_pt1->Divide(mixClone_pt1);
    norm_pt1->Draw("e");
    norm_pt1->GetXaxis()->SetTitle("#Deltap_{T}(GeV/c)");
    norm_pt1->GetYaxis()->SetTitle("c(#Deltap_{T})");
    norm_pt1->SetTitle("k^{-}#Omega^{-}");
    norm_pt1->SetName("norm_kmom_pt");
    c3->Draw();

    // k normalization kmom
    TF1 *constant_k1 = new TF1("constant_k1", "1", 0, 10);
    Int_t binNorm_k1[2];
    binNorm_k1[0] = same_k1->FindBin(0.8);
    binNorm_k1[1] = same_k1->FindBin(1.3);
    Double_t factorN1_k1 = mix_k1->Integral(binNorm_k1[0], binNorm_k1[1]) / same_k1->Integral(binNorm_k1[0], binNorm_k1[1]);

    TCanvas *c4 = new TCanvas();
    TH1D *norm_k1 = (TH1D *)same_k1->Clone();
    TH1D *mixClone_k1 = (TH1D *)mix_k1->Clone();

    norm_k1->Sumw2();
    mixClone_k1->Sumw2();

    mixClone_k1->Divide(constant_k1, factorN1_k1);
    norm_k1->Divide(mixClone_k1);
    norm_k1->Draw("e");
    norm_k1->GetXaxis()->SetTitle("k^{*}");
    norm_k1->GetYaxis()->SetTitle("C(k^{*})");
    // norm_k1->SetName("norm_kpom_k");
    norm_k1->SetName("norm_kmom_k");
    norm_k1->SetTitle("k^{-}#Omega^{-}");

    // pt normalization  kmop
    TF1 *constant_pt2 = new TF1("constant_pt2", "1", 0, 10);
    Int_t binNorm_pt2[2];
    binNorm_pt2[0] = same_pt2->FindBin(0.8);
    binNorm_pt2[1] = same_pt2->FindBin(1.3);
    Double_t factorN_pt2 = mix_pt2->Integral(binNorm_pt2[0], binNorm_pt2[1]) / same_pt2->Integral(binNorm_pt2[0], binNorm_pt2[1]);
    cout << "normalization factor: " << factorN_pt2 << endl;
    TCanvas *c11 = new TCanvas();
    same_pt2->Draw("e");
    same_pt2->GetXaxis()->SetTitle("#Deltap_{T}");
    same_pt2->GetYaxis()->SetTitle("count");
    c11->Draw();
    TCanvas *c22 = new TCanvas();
    mix_pt2->Draw("e");
    mix_pt2->GetXaxis()->SetTitle("#Deltap_{T}");
    mix_pt2->GetYaxis()->SetTitle("count");
    c22->Draw();
    TCanvas *c33 = new TCanvas();
    TH1D *norm_pt2 = (TH1D *)same_pt2->Clone();
    TH1D *mixClone_pt2 = (TH1D *)mix_pt2->Clone();

    norm_pt2->Sumw2();
    mixClone_pt2->Sumw2();

    mixClone_pt2->Divide(constant_pt2, factorN_pt2);
    norm_pt2->Divide(mixClone_pt2);
    norm_pt2->Draw("e");
    norm_pt2->GetXaxis()->SetTitle("#Deltap_{T}(GeV/c)");
    norm_pt2->GetYaxis()->SetTitle("C(#Deltap_{T})");
    norm_pt2->SetTitle("k^{-}#Omega^{+}");
    norm_pt2->SetName("norm_kmop_pt");
    c33->Draw();

    // k normalization kpop
    TF1 *constant_k2 = new TF1("constant_k2", "1", 0, 10);
    Int_t binNorm_k2[2];
    binNorm_k2[0] = same_k2->FindBin(0.8);
    binNorm_k2[1] = same_k2->FindBin(1.3);
    Double_t factorN1_k2 = mix_k2->Integral(binNorm_k2[0], binNorm_k2[1]) / same_k2->Integral(binNorm_k2[0], binNorm_k2[1]);

    TCanvas *c44 = new TCanvas();
    TH1D *norm_k2 = (TH1D *)same_k2->Clone();
    TH1D *mixClone_k2 = (TH1D *)mix_k2->Clone();

    norm_k2->Sumw2();
    mixClone_k2->Sumw2();


    mixClone_k2->Divide(constant_k2, factorN1_k2);
    norm_k2->Divide(mixClone_k2);
    norm_k2->Draw("e");
    norm_k2->GetXaxis()->SetTitle("k^{*}");
    norm_k2->GetYaxis()->SetTitle("C(k^{*})");
    norm_k2->SetName("norm_kmop_k");
    norm_k2->SetTitle("k^{-}#Omega^{+}");


    TH1D *kmom_deduction_kmop_pt = new TH1D("kmom_deduction_kmop_pt", "k^{-}#Omega^{-}-k^{-}#Omega^{+};#Deltap_{T};counts", 40, 0.0, 4);
    TH1D *kmom_deduction_kmop_k = new TH1D("kmom_deduction_kmop_k", "k^{-}#Omega^{-}-k^{-}#Omega^{+};k^{*};counts", 40, 0.0, 4);
    kmom_deduction_kmop_pt->Add(norm_pt1, norm_pt2, 1, -1);
    kmom_deduction_kmop_k->Add(norm_k1, norm_k2, 1, -1);

    TFile *result_corr = new TFile("kmom_kmop_mix.root", "RECREATE");
    same_pt1->Write();
    mix_pt1->Write();
    norm_pt1->Write();

    same_k1->Write();
    mix_k1->Write();
    norm_k1->Write();

    same_pt2->Write();
    mix_pt2->Write();
    norm_pt2->Write();

    same_k2->Write();
    mix_k2->Write();
    norm_k2->Write();

    kmom_deduction_kmop_pt->Write();
    kmom_deduction_kmop_k->Write();

    result_corr->Write();
    return 0;
}
