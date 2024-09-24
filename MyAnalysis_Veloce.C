#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <cmath>
#include <iostream>

double MyAnalysis::t(
    int &i) { // tolgo il tempo dallo SC al target, che lui conta in TWTOF
  return TWTOF->at(GLBtrackTWid->at(i)) * 1e-9 - dt;
}
double MyAnalysis::p(int &i) {
  return sqrt(pow(GLBtrackPx->at(i), 2) + pow(GLBtrackPy->at(i), 2) +
              pow(GLBtrackPz->at(i), 2));
  //* 1e3;
}
double MyAnalysis::beta(int &i) {
  return (GLBtrackLength->at(i) * 1e-2 / (c * t(i))); // traccia in cm
}
double MyAnalysis::E_k(int &i) { 
  return CAenergy->at(GLBtrackCAid->at(i)) +
         TWDe1Point->at(GLBtrackTWid->at(i)) * 1e-3 + //era GeV
         TWDe2Point->at(GLBtrackTWid->at(i)) * 1e-3;
}
int MyAnalysis::Z(int &i) { return TWChargePoint->at(GLBtrackTWid->at(i)); }

/*FUNZIONI PER L'ANALISI*/
void MyAnalysis::A_Z_solo(vector<TH1D *> *h_i_1, vector<TH1D *> *h_i_2,
                          vector<TH1D *> *h_i_3) {}
void MyAnalysis::compare_A_Z_solo(vector<TGraph *> *g_1_2_Z,
                                  vector<TGraph *> *g_2_3_Z,
                                  vector<TGraph *> *g_1_3_Z) {}

void MyAnalysis::Mass(vector<TH1D *> *h_i_1, vector<TH1D *> *h_i_2,
                         vector<TH1D *> *h_i_3, vector<TH2D *> *h_i_1_2,
                         vector<TH2D *> *h_i_2_3, vector<TH2D *> *h_i_1_3) {
  for (int j = 0; j < GLBtrackTWid->size(); j++) {
    if (GLBtrackTWid->at(j) < 0 ||
        GLBtrackCAid->at(j) < 0) { // tolgo i casi in cui indice è negativo
      continue;
    }
    /*RIEMPIMENTO ISTOGRAMMI*/
    A_1 = (p(j) / (u * beta(j) * (1 / sqrt(1 - pow(beta(j), 2)))));
    h_i_1->at(Z(j) - 1)->Fill(A_1);

    A_2 = (E_k(j) / (u * ((1 / sqrt(1 - pow((beta(j)), 2))) - 1)));
    h_i_2->at(Z(j) - 1)->Fill(A_2);

    A_3 = ((pow(p(j), 2) - pow(E_k(j), 2)) / (2 * u * E_k(j))); 
    h_i_3->at(Z(j) - 1)->Fill(A_3);
    
    h_i_1_2->at(Z(j) - 1)->Fill(A_1, A_2);
    h_i_2_3->at(Z(j) - 1)->Fill(A_2, A_3);
    h_i_1_3->at(Z(j) - 1)->Fill(A_1, A_3);
  }
}

void MyAnalysis::Frac_en_p(vector<TH1D *> *en, vector<TH1D *> *imp) {}

/*FUNZIONI FINALI*/
void MyAnalysis::Analysis() {
  
  BeforeLoop();
  Loop();
  AfterLoop();
  /*
  BeforeLoop_A();
  Loop_A();
  AfterLoop_A();
  */
}

void MyAnalysis::Loop() {}
void MyAnalysis::BeforeLoop() {}
void MyAnalysis::AfterLoop() {
  TCanvas *c_A_1 = new TCanvas("c_A_1", "A_1", 2000, 1000);
  TCanvas *c_A_2 = new TCanvas("c_A_2", "A_2", 2000, 1000);
  TCanvas *c_A_3 = new TCanvas("c_A_3", "A_3", 2000, 1000);
  TCanvas *c_A_1_2 = new TCanvas("c_A_1_2", "A_1_2", 2000, 1000);
  TCanvas *c_A_2_3 = new TCanvas("c_A_2_3", "A_2_3", 2000, 1000);
  TCanvas *c_A_1_3 = new TCanvas("c_A_1_3", "A_1_3", 2000, 1000);
  c_A_1->Divide(2, 3);
  c_A_2->Divide(2, 3);
  c_A_3->Divide(2, 3);
  c_A_1_2->Divide(2, 3);
  c_A_2_3->Divide(2, 3);
  c_A_1_3->Divide(2, 3);
  TFile *file_A = new TFile("A_histograms.root", "READ");
  TFile *file_comp_A = new TFile("comp_A_histograms.root", "READ");
  std::vector<TH1D*>* h_A_1 = new std::vector<TH1D*>(6, nullptr);
  std::vector<TH1D*>* h_A_2 = new std::vector<TH1D*>(6, nullptr); 
  std::vector<TH1D*>* h_A_3 = new std::vector<TH1D*>(6, nullptr); 
  std::vector<TH2D*>* h_comp_A_1_2 = new std::vector<TH2D*>(6, nullptr);
  std::vector<TH2D*>* h_comp_A_2_3 = new std::vector<TH2D*>(6, nullptr);
  std::vector<TH2D*>* h_comp_A_1_3 = new std::vector<TH2D*>(6, nullptr);
    
  for (int i = 0; i < 6; ++i) {
        std::string histName_1 = "h_A_1_Z_" + std::to_string(i+1);  
        std::string histName_2 = "h_A_2_Z_" + std::to_string(i+1);
        std::string histName_3 = "h_A_3_Z_" + std::to_string(i+1);
        std::string histName_1_2 = "h_A_1_2_Z_" + std::to_string(i+1);
        std::string histName_2_3 = "h_A_2_3_Z_" + std::to_string(i+1);
        std::string histName_1_3 = "h_A_1_3_Z_" + std::to_string(i+1);
        (*h_A_1)[i] = (TH1D*)file_A->Get(histName_1.c_str());  
        (*h_A_2)[i] = (TH1D*)file_A->Get(histName_2.c_str());
        (*h_A_3)[i] = (TH1D*)file_A->Get(histName_3.c_str());
        (*h_comp_A_1_2)[i] = (TH2D*)file_comp_A->Get(histName_1_2.c_str());
        (*h_comp_A_2_3)[i] = (TH2D*)file_comp_A->Get(histName_2_3.c_str());
        (*h_comp_A_1_3)[i] = (TH2D*)file_comp_A->Get(histName_1_3.c_str());
    }
    for (int i = 0; i < 6; i++) {
    c_A_1->cd(i + 1);
    (*h_A_1)[i]->Draw();
    c_A_2->cd(i + 1);
    (*h_A_2)[i]->Draw(); 
    c_A_3->cd(i + 1);
    (*h_A_3)[i]->Draw();
    c_A_1_2->cd(i + 1);
    (*h_comp_A_1_2)[i]->Draw("colz");
    c_A_2_3->cd(i + 1);
    (*h_comp_A_2_3)[i]->Draw("colz");
    c_A_1_3->cd(i + 1);
    (*h_comp_A_1_3)[i]->Draw("colz");
  }
  c_A_1->Update();
  c_A_2->Update();
  c_A_3->Update();
  c_A_1_2->Update();
  c_A_2_3->Update();
  c_A_1_3->Update();

  c_A_1->WaitPrimitive();
  c_A_2->WaitPrimitive();
  c_A_3->WaitPrimitive();
  c_A_1_2->WaitPrimitive();
  c_A_2_3->WaitPrimitive();
  c_A_1_3->WaitPrimitive();

  file_A->Close();
  file_comp_A->Close();
}




void MyAnalysis::Loop_A() {
  if (fChain == 0)
    return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 10000;
  //Long64_t nentries = 200000;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
  
    Mass(h_A_1_Z, h_A_2_Z, h_A_3_Z, h_compare_A_1_2_Z, h_compare_A_2_3_Z, h_compare_A_1_3_Z);
  }
    TFile *file_A = new TFile("A_histograms.root", "RECREATE");
    TFile *file_comp_A = new TFile("comp_A_histograms.root", "RECREATE");
    for (int i = 0; i < 6; i++)
    {
      file_A->cd();
      (*h_A_1_Z)[i]->Write();
      file_comp_A->cd();
      (*h_compare_A_1_2_Z)[i]->Write();
    }
    for (int i = 0; i < 6; i++)
    {
      file_A->cd();
      (*h_A_2_Z)[i]->Write();
      file_comp_A->cd();
      (*h_compare_A_2_3_Z)[i]->Write();
    }
    for (int i = 0; i < 6; i++)
    {
      file_A->cd();
      (*h_A_3_Z)[i]->Write();
      file_comp_A->cd();
      (*h_compare_A_1_3_Z)[i]->Write();
    }
    file_A->Close();
    file_comp_A->Close();
}
void MyAnalysis::BeforeLoop_A() {
  h_A_1_Z = new vector<TH1D *>(6);
  h_A_2_Z = new vector<TH1D *>(6);
  h_A_3_Z = new vector<TH1D *>(6);
  h_compare_A_1_2_Z = new vector<TH2D *>(6);
  h_compare_A_2_3_Z = new vector<TH2D *>(6);
  h_compare_A_1_3_Z = new vector<TH2D *>(6);
  for (int i = 0; i < 6; i++) {
    (*h_A_1_Z)[i] = new TH1D(([](int idx) {
                                std::ostringstream oss;
                                oss << "h_A_1_Z_" << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 20);
    (*h_A_1_Z)[i]->SetLineColor(kBlack);
    (*h_A_1_Z)[i]->SetXTitle("A1");   
    (*h_A_2_Z)[i] = new TH1D(([](int idx) {
                                std::ostringstream oss;
                                oss << "h_A_2_Z_" << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 20);
    (*h_A_2_Z)[i]->SetLineColor(kBlack);
    (*h_A_2_Z)[i]->SetXTitle("A2");
    (*h_A_3_Z)[i] = new TH1D(([](int idx) {
                                std::ostringstream oss;
                                oss << "h_A_3_Z_" << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 20);
    (*h_A_3_Z)[i]->SetLineColor(kBlack);
    (*h_A_3_Z)[i]->SetXTitle("A3");
    (*h_compare_A_1_2_Z)[i] = new TH2D(([](int idx) {
                                std::ostringstream oss;
                                oss << "h_A_1_2_Z_" << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 20, 1000, 0.0, 20);
    (*h_compare_A_1_2_Z)[i]->SetXTitle("A_1");
    (*h_compare_A_1_2_Z)[i]->SetYTitle("A_2");
    (*h_compare_A_2_3_Z)[i] = new TH2D(([](int idx) {
                                std::ostringstream oss;
                                oss << "h_A_2_3_Z_" << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 20, 1000, 0.0, 20);
    (*h_compare_A_2_3_Z)[i]->SetXTitle("A_2");
    (*h_compare_A_2_3_Z)[i]->SetYTitle("A_3");
    (*h_compare_A_1_3_Z)[i] = new TH2D(([](int idx) {
                                std::ostringstream oss;
                                oss << "h_A_1_3_Z_" << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 20, 1000, 0.0, 20);
    (*h_compare_A_1_3_Z)[i]->SetXTitle("A_1");
    (*h_compare_A_1_3_Z)[i]->SetYTitle("A_3");                                                                             
    if(i == 0) {
      (*h_A_1_Z)[i]->GetXaxis()->SetRangeUser(0,4);
      (*h_A_2_Z)[i]->GetXaxis()->SetRangeUser(0,4);
      (*h_A_3_Z)[i]->GetXaxis()->SetRangeUser(0,8);
      
      (*h_compare_A_1_2_Z)[i]->GetXaxis()->SetRangeUser(0,4);//1
      (*h_compare_A_1_2_Z)[i]->GetYaxis()->SetRangeUser(0,4);//2
      (*h_compare_A_2_3_Z)[i]->GetXaxis()->SetRangeUser(0,4);//2
      (*h_compare_A_2_3_Z)[i]->GetYaxis()->SetRangeUser(0,8);//3
      (*h_compare_A_1_3_Z)[i]->GetXaxis()->SetRangeUser(0,4);//1
      (*h_compare_A_1_3_Z)[i]->GetYaxis()->SetRangeUser(0,8);//3
    }
    else if(i == 1) {
      (*h_A_1_Z)[i]->GetXaxis()->SetRangeUser(1,6);
      (*h_A_2_Z)[i]->GetXaxis()->SetRangeUser(0,8);
      (*h_A_3_Z)[i]->GetXaxis()->SetRangeUser(0,12);

      (*h_compare_A_1_2_Z)[i]->GetXaxis()->SetRangeUser(1,6);//1
      (*h_compare_A_1_2_Z)[i]->GetYaxis()->SetRangeUser(0,8);//2
      (*h_compare_A_2_3_Z)[i]->GetXaxis()->SetRangeUser(0,8);//2
      (*h_compare_A_2_3_Z)[i]->GetYaxis()->SetRangeUser(0,12);//3
      (*h_compare_A_1_3_Z)[i]->GetXaxis()->SetRangeUser(1,6);//1
      (*h_compare_A_1_3_Z)[i]->GetYaxis()->SetRangeUser(0,12);//3
    }
    else if(i == 2) {
      (*h_A_1_Z)[i]->GetXaxis()->SetRangeUser(4, 10);
      (*h_A_2_Z)[i]->GetXaxis()->SetRangeUser(0,8);
      (*h_A_3_Z)[i]->GetXaxis()->SetRangeUser(4, 14);

      (*h_compare_A_1_2_Z)[i]->GetXaxis()->SetRangeUser(4, 10);//1
      (*h_compare_A_1_2_Z)[i]->GetYaxis()->SetRangeUser(0,8);//2
      (*h_compare_A_2_3_Z)[i]->GetXaxis()->SetRangeUser(0,8);//2
      (*h_compare_A_2_3_Z)[i]->GetYaxis()->SetRangeUser(4, 14);//3
      (*h_compare_A_1_3_Z)[i]->GetXaxis()->SetRangeUser(4, 10);//1
      (*h_compare_A_1_3_Z)[i]->GetYaxis()->SetRangeUser(4, 14);//3
    }
    else if(i == 3) {
      (*h_A_1_Z)[i]->GetXaxis()->SetRangeUser(5,12);
      (*h_A_2_Z)[i]->GetXaxis()->SetRangeUser(0,12);
      (*h_A_3_Z)[i]->GetXaxis()->SetRangeUser(5,12);

      (*h_compare_A_1_2_Z)[i]->GetXaxis()->SetRangeUser(5,12);//1
      (*h_compare_A_1_2_Z)[i]->GetYaxis()->SetRangeUser(0,12);//2
      (*h_compare_A_2_3_Z)[i]->GetXaxis()->SetRangeUser(0,12);//2
      (*h_compare_A_2_3_Z)[i]->GetYaxis()->SetRangeUser(5,12);//3
      (*h_compare_A_1_3_Z)[i]->GetXaxis()->SetRangeUser(5,12);//1
      (*h_compare_A_1_3_Z)[i]->GetYaxis()->SetRangeUser(5,12);//3
    }
    else if(i == 4) {
      (*h_A_1_Z)[i]->GetXaxis()->SetRangeUser(7,14);
      (*h_A_2_Z)[i]->GetXaxis()->SetRangeUser(0,12);
      (*h_A_3_Z)[i]->GetXaxis()->SetRangeUser(7,14);

      (*h_compare_A_1_2_Z)[i]->GetXaxis()->SetRangeUser(7, 14);//1
      (*h_compare_A_1_2_Z)[i]->GetYaxis()->SetRangeUser(0,12);//2
      (*h_compare_A_2_3_Z)[i]->GetXaxis()->SetRangeUser(0,12);//2
      (*h_compare_A_2_3_Z)[i]->GetYaxis()->SetRangeUser(7, 14);//3
      (*h_compare_A_1_3_Z)[i]->GetXaxis()->SetRangeUser(7, 14);//1
      (*h_compare_A_1_3_Z)[i]->GetYaxis()->SetRangeUser(7, 14);//3
    }
    else if(i == 5) {
      (*h_A_1_Z)[i]->GetXaxis()->SetRangeUser(10,14);
      (*h_A_2_Z)[i]->GetXaxis()->SetRangeUser(0,14);
      (*h_A_3_Z)[i]->GetXaxis()->SetRangeUser(10,16);

      (*h_compare_A_1_2_Z)[i]->GetXaxis()->SetRangeUser(10, 14);//1
      (*h_compare_A_1_2_Z)[i]->GetYaxis()->SetRangeUser(0,14);//2
      (*h_compare_A_2_3_Z)[i]->GetXaxis()->SetRangeUser(0,14);//2
      (*h_compare_A_2_3_Z)[i]->GetYaxis()->SetRangeUser(7, 16);//3
      (*h_compare_A_1_3_Z)[i]->GetXaxis()->SetRangeUser(10, 14);//1
      (*h_compare_A_1_3_Z)[i]->GetYaxis()->SetRangeUser(7, 16);//3
      
    }
  }
}
void MyAnalysis::AfterLoop_A() {
  TCanvas *c_A_1 = new TCanvas("c_A_1", "A_1", 2000, 1000);
  TCanvas *c_A_2 = new TCanvas("c_A_2", "A_2", 2000, 1000);
  TCanvas *c_A_3 = new TCanvas("c_A_3", "A_3", 2000, 1000);
  TCanvas *c_A_1_2 = new TCanvas("c_A_1_2", "A_1_2", 2000, 1000);
  TCanvas *c_A_2_3 = new TCanvas("c_A_2_3", "A_2_3", 2000, 1000);
  TCanvas *c_A_1_3 = new TCanvas("c_A_1_3", "A_1_3", 2000, 1000);
  c_A_1->Divide(2, 3);
  c_A_2->Divide(2, 3);
  c_A_3->Divide(2, 3);
  c_A_1_2->Divide(2, 3);
  c_A_2_3->Divide(2, 3);
  c_A_1_3->Divide(2, 3);

    /*  
  TFile *file_A = new TFile("A_histograms.root", "READ");
  file_A->ls();
  if (!file_A || file_A->IsZombie()) {
        std::cerr << "Errore: il file non può essere aperto!" << std::endl;
        return;
    }

  (*h_A_1_Z)[0] = (TH1D*)file_A->Get("h_A_1_Z_1");
  (*h_A_1_Z)[1] = (TH1D*)file_A->Get("h_A_1_Z_2");
  (*h_A_1_Z)[2] = (TH1D*)file_A->Get("h_A_1_Z_3");
  (*h_A_1_Z)[3] = (TH1D*)file_A->Get("h_A_1_Z_4");
  (*h_A_1_Z)[4] = (TH1D*)file_A->Get("h_A_1_Z_5");
  (*h_A_1_Z)[5] = (TH1D*)file_A->Get("h_A_1_Z_6");
  (*h_A_2_Z)[0] = (TH1D*)file_A->Get("h_A_2_Z_1");
  (*h_A_2_Z)[1] = (TH1D*)file_A->Get("h_A_2_Z_2");
  (*h_A_2_Z)[2] = (TH1D*)file_A->Get("h_A_2_Z_3");
  (*h_A_2_Z)[3] = (TH1D*)file_A->Get("h_A_2_Z_4");
  (*h_A_2_Z)[4] = (TH1D*)file_A->Get("h_A_2_Z_5");
  (*h_A_2_Z)[5] = (TH1D*)file_A->Get("h_A_2_Z_6");
  (*h_A_3_Z)[0] = (TH1D*)file_A->Get("h_A_3_Z_1");
  (*h_A_3_Z)[1] = (TH1D*)file_A->Get("h_A_3_Z_2");
  (*h_A_3_Z)[2] = (TH1D*)file_A->Get("h_A_3_Z_3");
  (*h_A_3_Z)[3] = (TH1D*)file_A->Get("h_A_3_Z_4");
  (*h_A_3_Z)[4] = (TH1D*)file_A->Get("h_A_3_Z_5");
  (*h_A_3_Z)[5] = (TH1D*)file_A->Get("h_A_3_Z_6");
  TFile *file_comp_A = new TFile("comp_A_histograms.root", "READ");
  (*h_compare_A_1_2_Z)[1] = (TH2D*)file_comp_A->Get("h_compare_A_1_2_Z_1");
  (*h_compare_A_1_2_Z)[2] = (TH2D*)file_comp_A->Get("h_compare_A_1_2_Z_2");
  (*h_compare_A_1_2_Z)[3] = (TH2D*)file_comp_A->Get("h_compare_A_1_2_Z_3");
  (*h_compare_A_1_2_Z)[0] = (TH2D*)file_comp_A->Get("h_compare_A_1_2_Z_4");
  (*h_compare_A_1_2_Z)[4] = (TH2D*)file_comp_A->Get("h_compare_A_1_2_Z_5");
  (*h_compare_A_1_2_Z)[5] = (TH2D*)file_comp_A->Get("h_compare_A_1_2_Z_6");
  (*h_compare_A_2_3_Z)[0] = (TH2D*)file_comp_A->Get("h_compare_A_2_3_Z_1");
  (*h_compare_A_2_3_Z)[1] = (TH2D*)file_comp_A->Get("h_compare_A_2_3_Z_2");
  (*h_compare_A_2_3_Z)[2] = (TH2D*)file_comp_A->Get("h_compare_A_2_3_Z_3");
  (*h_compare_A_2_3_Z)[3] = (TH2D*)file_comp_A->Get("h_compare_A_2_3_Z_4");
  (*h_compare_A_2_3_Z)[4] = (TH2D*)file_comp_A->Get("h_compare_A_2_3_Z_5");
  (*h_compare_A_2_3_Z)[5] = (TH2D*)file_comp_A->Get("h_compare_A_2_3_Z_6");
  (*h_compare_A_1_3_Z)[0] = (TH2D*)file_comp_A->Get("h_compare_A_1_3_Z_1");
  (*h_compare_A_1_3_Z)[1] = (TH2D*)file_comp_A->Get("h_compare_A_1_3_Z_2");
  (*h_compare_A_1_3_Z)[2] = (TH2D*)file_comp_A->Get("h_compare_A_1_3_Z_3");
  (*h_compare_A_1_3_Z)[3] = (TH2D*)file_comp_A->Get("h_compare_A_1_3_Z_4");
  (*h_compare_A_1_3_Z)[4] = (TH2D*)file_comp_A->Get("h_compare_A_1_3_Z_5");
  (*h_compare_A_1_3_Z)[5] = (TH2D*)file_comp_A->Get("h_compare_A_1_3_Z_6");
  */
  for (int i = 0; i < 6; i++) {
    c_A_1->cd(i + 1);
    (*h_A_1_Z)[i]->Draw();
    c_A_2->cd(i + 1);
    (*h_A_2_Z)[i]->Draw();
    c_A_3->cd(i + 1);
    (*h_A_3_Z)[i]->Draw();
    c_A_1_2->cd(i + 1);
    (*h_compare_A_1_2_Z)[i]->Draw("colz");
    c_A_2_3->cd(i + 1);
    (*h_compare_A_2_3_Z)[i]->Draw("colz");
    c_A_1_3->cd(i + 1);
    (*h_compare_A_1_3_Z)[i]->Draw("colz");
  }
  //file_A->Close();
  //file_comp_A->Close();

}

void MyAnalysis::BeforeLoop_Frac() {};
//void MyAnalysis::BeforeLoop_A() {};
void MyAnalysis::AfterLoop_Frac() {};
//void MyAnalysis::AfterLoop_A() {};

/*
  TF1 *Gauss = new TF1("Gauss", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0, 16);
  Gauss->SetParameters(5000, 1, 0.5);
  (*h_A_1_Z)[0]->Fit("Gauss");
 // (*h_A_2_Z)[0]->Fit("Gauss");
 // (*h_A_3_Z)[0]->Fit("Gauss");
  Gauss->SetParameters(3000, 4.1, 0.5);
  (*h_A_1_Z)[1]->Fit("Gauss");
 // (*h_A_2_Z)[1]->Fit("Gauss");
 // (*h_A_3_Z)[1]->Fit("Gauss");
  Gauss->SetParameters(200, 6.1, 0.5);
  (*h_A_1_Z)[2]->Fit("Gauss");
 // (*h_A_2_Z)[2]->Fit("Gauss");
 // (*h_A_3_Z)[2]->Fit("Gauss");
  Gauss->SetParameters(120, 7.2, 0.5);
  (*h_A_1_Z)[3]->Fit("Gauss");
 // (*h_A_2_Z)[3]->Fit("Gauss");
  //(*h_A_3_Z)[3]->Fit("Gauss");
  Gauss->SetParameters(200, 10.3, 0.5);
  (*h_A_1_Z)[4]->Fit("Gauss");
 // (*h_A_2_Z)[4]->Fit("Gauss");
 // (*h_A_3_Z)[4]->Fit("Gauss");
  Gauss->SetParameters(80000, 12.3, 0.5);
  (*h_A_1_Z)[5]->Fit("Gauss");
 // (*h_A_2_Z)[5]->Fit("Gauss");
 // (*h_A_3_Z)[5]->Fit("Gauss");

  c_A_1->Print("A1_fit.pdf");
  c_A_2->Print("A2_fit.pdf");
  c_A_3->Print("A3_fit.pdf");*/