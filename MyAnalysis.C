#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include <TCanvas.h>
#include <TH2.h>
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
// return ((CAenergy->at(GLBtrackCAid->at(i)) )* 1e3
// +TWDe1Point->at(GLBtrackTWid->at(i))+TWDe2Point->at(GLBtrackTWid->at(i)));


/*FUNZIONI PER L'ANALISI*/
void MyAnalysis::A_Z_solo(vector<TH1D *> *h_i_1, vector<TH1D *> *h_i_2,
                          vector<TH1D *> *h_i_3) {
  for (int j = 0; j < GLBtrackTWid->size();
       j++) { // ipotizzando stessa dimensione vettori
    if (GLBtrackTWid->at(j) < 0 ||
        GLBtrackCAid->at(j) < 0) { // tolgo i casi in cui indice è negativo
      continue;
    }
    A_1 = (p(j) / (u * beta(j) * (1 / sqrt(1 - pow(beta(j), 2)))));
    h_i_1->at(Z(j) - 1)->Fill(A_1);

    A_2 = (E_k(j) / (u * ((1 / sqrt(1 - pow((beta(j)), 2))) - 1)));
    h_i_2->at(Z(j) - 1)->Fill(A_2);

    A_3 = ((pow(p(j), 2) - pow(E_k(j), 2)) / (2 * E_k(j)));
    h_i_3->at(Z(j) - 1)->Fill(A_3);
  }
}
void MyAnalysis::compare_A_Z_solo(vector<TGraph *> *g_1_2_Z,
                                  vector<TGraph *> *g_2_3_Z,
                                  vector<TGraph *> *g_1_3_Z) {
  for (int j = 0; j < GLBtrackTWid->size();
       j++) { // ipotizzando stessa dimensione vettori
    if (GLBtrackTWid->at(j) < 0 ||
        GLBtrackCAid->at(j) < 0) { // tolgo i casi in cui indice è negativo
      continue;
    }
    A_1 = (p(j) / (u * beta(j) * (1 / sqrt(1 - pow(beta(j), 2)))));
    A_2 = (E_k(j) / (u * ((1 / sqrt(1 - pow((beta(j)), 2))) - 1)));
    A_3 = ((pow(p(j), 2) - pow(E_k(j), 2)) / (2 * E_k(j)));
    g_1_2_Z->at(Z(j) - 1)->SetPoint((*g_1_2_Z)[Z(j) - 1]->GetN(), A_1, A_2);
    g_2_3_Z->at(Z(j) - 1)->SetPoint((*g_2_3_Z)[Z(j) - 1]->GetN(), A_2, A_3);
    g_1_3_Z->at(Z(j) - 1)->SetPoint((*g_1_3_Z)[Z(j) - 1]->GetN(), A_1, A_3);
  }
}
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

    /*g_1_2_Z->at(Z(j) - 1)->SetPoint(j, A_1, A_2);NON SO PERCHE COSI NON VA
    g_2_3_Z->at(Z(j) - 1)->SetPoint(j, A_2, A_3);
    g_1_3_Z->at(Z(j) - 1)->SetPoint(j, A_1, A_3);*/
   // g_1_2_Z->at(Z(j) - 1)->SetPoint((*g_1_2_Z)[Z(j) - 1]->GetN(), A_1, A_2);
   // g_2_3_Z->at(Z(j) - 1)->SetPoint((*g_2_3_Z)[Z(j) - 1]->GetN(), A_2, A_3);
   // g_1_3_Z->at(Z(j) - 1)->SetPoint((*g_1_3_Z)[Z(j) - 1]->GetN(), A_1, A_3);
   
    h_i_1_2->at(Z(j) - 1)->Fill(A_1, A_2);
    h_i_2_3->at(Z(j) - 1)->Fill(A_2, A_3);
    h_i_1_3->at(Z(j) - 1)->Fill(A_1, A_3);
  }
}
void MyAnalysis::Frac_en_p(vector<TH1D *> *en, vector<TH1D *> *imp) {
  int TW_id;
  int MC_id_1;
  int MC_id_2;
  for (int j = 0; j < GLBtrackTWid->size();
       j++) { // ipotizzando stessa dimensione vettori
    TW_id = GLBtrackTWid->at(j);
    if ( TW_id < 0 ||
        GLBtrackCAid->at(j) < 0) { // tolgo i casi in cui indice è negativo
      continue;
    }
    //MC_id_1 = TATW_MCID_1->at(TW_ind);
    //MC_id_2 = TATW_MCID_2->at(TW_ind);
    MC_id_1 = (*TATW_MCID_1)[TW_id];
    MC_id_2 = (*TATW_MCID_2)[TW_id];
    if (MC_id_1 == MC_id_2) {
      en->at(Z(j) - 1)->Fill(E_k(j) / ( (sqrt(pow((*MC_InitMom_x)[ MC_id_1],2) + pow((*MC_InitMom_y)[ MC_id_1],2)
       + pow((*MC_InitMom_z)[ MC_id_1],2) + pow((*MC_Mass)[MC_id_1], 2))) - (*MC_Mass)[MC_id_1]) ); 
      imp->at(Z(j) - 1)->Fill(p(j) / sqrt(pow((*MC_InitMom_x)[ MC_id_1],2) + pow((*MC_InitMom_y)[ MC_id_1],2)
       + pow((*MC_InitMom_z)[ MC_id_1],2)) );
    } 
    else{ 
      continue;
    }
  }
}
void MyAnalysis::Loop() {
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
    // A_Z_solo(h_A_1_Z, h_A_2_Z, h_A_3_Z);
    // compare_A_Z_solo(g_compare_A_1_2_Z, g_compare_A_2_3_Z,
    // g_compare_A_1_3_Z);
    Mass(h_A_1_Z, h_A_2_Z, h_A_3_Z, h_compare_A_1_2_Z, h_compare_A_2_3_Z, h_compare_A_1_3_Z);
    //Frac_en_p(frac_en, frac_p);

  }
}
void MyAnalysis::Analysis() {
  BeforeLoop();
  //BeforeLoop_Frac();
  Loop();
  AfterLoop();
  //AfterLoop_Frac();
}
void MyAnalysis::BeforeLoop() {
  h_A_1_Z = new vector<TH1D *>(6);
  h_A_2_Z = new vector<TH1D *>(6);
  h_A_3_Z = new vector<TH1D *>(6);
  frac_en = new vector<TH1D *>(6);
  //g_compare_A_1_2_Z = new vector<TGraph *>(6);
  //g_compare_A_2_3_Z = new vector<TGraph *>(6);
  //g_compare_A_1_3_Z = new vector<TGraph *>(6);
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
    /*std::ostringstream title; 
    title << "Z = " << i + 1;
    (*g_compare_A_1_2_Z)[i] = new TGraph();
    (*g_compare_A_1_2_Z)[i]->SetTitle(title.str().c_str());
    (*g_compare_A_1_2_Z)[i]->GetXaxis()->SetTitle("A_1");
    (*g_compare_A_1_2_Z)[i]->GetYaxis()->SetTitle("A_2");
    (*g_compare_A_2_3_Z)[i] = new TGraph();
    (*g_compare_A_2_3_Z)[i]->SetTitle(title.str().c_str());
    (*g_compare_A_2_3_Z)[i]->GetXaxis()->SetTitle("A_2");
    (*g_compare_A_2_3_Z)[i]->GetYaxis()->SetTitle("A_3");
    (*g_compare_A_1_3_Z)[i] = new TGraph();
    (*g_compare_A_1_3_Z)[i]->SetTitle(title.str().c_str());
    (*g_compare_A_1_3_Z)[i]->GetXaxis()->SetTitle("A_1");
    (*g_compare_A_1_3_Z)[i]->GetYaxis()->SetTitle("A_3");*/

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
void MyAnalysis::AfterLoop() {
  TCanvas *c_A_1 = new TCanvas("c_A_1", "A_1", 2000, 1000);
  TCanvas *c_A_2 = new TCanvas("c_A_2", "A_2", 2000, 1000);
  TCanvas *c_A_3 = new TCanvas("c_A_3", "A_3", 2000, 1000);
  c_A_1->Divide(2, 3);
  c_A_2->Divide(2, 3);
  c_A_3->Divide(2, 3);
  TCanvas *c_compare_1_2 = new TCanvas("c_compare_1_2", "1_2", 2000, 1000);
  TCanvas *c_compare_2_3 = new TCanvas("c_compare_2_3", "2_3", 2000, 1000);
  TCanvas *c_compare_1_3 = new TCanvas("c_compare_1_3", "1_3", 2000, 1000);
  c_compare_1_2->Divide(2, 3);
  c_compare_2_3->Divide(2, 3);
  c_compare_1_3->Divide(2, 3);
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
}
////////////////////////////////////////////////////////////////////////
void MyAnalysis::BeforeLoop_Frac() {
  frac_en = new vector<TH1D *>(6);
  frac_p = new vector<TH1D *>(6);
  for (int i = 0; i < 6; i++) {
    (*frac_en)[i] = new TH1D(([](int idx) {
                                std::ostringstream oss;
                                oss << "E_Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 1.5);
    (*frac_en)[i]->SetLineColor(kBlack);
    (*frac_en)[i]->SetXTitle("Ekin/Ekingen");
    (*frac_p)[i] = new TH1D(([](int idx) {
                                std::ostringstream oss;
                                oss << "p_Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              ([](int idx) {
                                std::ostringstream oss;
                                oss << "Z = " << idx + 1;
                                return oss.str();
                              })(i)
                                  .c_str(),
                              1000, 0.0, 1.5);
    (*frac_p)[i]->SetLineColor(kBlack);
    (*frac_p)[i]->SetXTitle("P/Pgen");
  }
}

void MyAnalysis::AfterLoop_Frac() {
  TCanvas *c_en =
      new TCanvas("Confronto Energia", "Confronto Energia", 2000, 1000);
  c_en->Divide(2, 3);
  TCanvas *c_p =
      new TCanvas("Confronto Impulso", "Confronto Impulso", 2000, 1000);
  c_p->Divide(2, 3);
  for (int i = 0; i < 6; i++) {
    c_en->cd(i + 1);
    (*frac_en)[i]->Draw();
    c_p->cd(i + 1);
    (*frac_p)[i]->Draw();
  }
}
////////////////////////////////////////////////////////////////////////
void MyAnalysis::BeforeLoop_A() {
  h_A_1_Z = new vector<TH1D *>(6);
  h_A_2_Z = new vector<TH1D *>(6);
  h_A_3_Z = new vector<TH1D *>(6);
  frac_en = new vector<TH1D *>(6);
  g_compare_A_1_2_Z = new vector<TGraph *>(6);
  g_compare_A_2_3_Z = new vector<TGraph *>(6);
  g_compare_A_1_3_Z = new vector<TGraph *>(6);
  for (int i = 0; i < 6; i++) {
    h_A_1_Z->at(i) = new TH1D(([](int idx) {
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
    h_A_2_Z->at(i) = new TH1D(([](int idx) {
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
    h_A_3_Z->at(i) = new TH1D(([](int idx) {
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
    std::ostringstream title;
    title << "Z = " << i + 1;
    g_compare_A_1_2_Z->at(i) = new TGraph();
    g_compare_A_1_2_Z->at(i)->SetTitle(title.str().c_str());
    g_compare_A_1_2_Z->at(i)->GetXaxis()->SetTitle("A_1");
    g_compare_A_1_2_Z->at(i)->GetYaxis()->SetTitle("A_2");
    g_compare_A_2_3_Z->at(i) = new TGraph();
    g_compare_A_2_3_Z->at(i)->SetTitle(title.str().c_str());
    g_compare_A_2_3_Z->at(i)->GetXaxis()->SetTitle("A_2");
    g_compare_A_2_3_Z->at(i)->GetYaxis()->SetTitle("A_3");
    g_compare_A_1_3_Z->at(i) = new TGraph();
    g_compare_A_1_3_Z->at(i)->SetTitle(title.str().c_str());
    g_compare_A_1_3_Z->at(i)->GetXaxis()->SetTitle("A_1");
    g_compare_A_1_3_Z->at(i)->GetYaxis()->SetTitle("A_3");
  }
}
void MyAnalysis::AfterLoop_A() {
  TCanvas *c_A_1 = new TCanvas("c_A_1", "A_1", 2000, 1000);
  TCanvas *c_A_2 = new TCanvas("c_A_2", "A_2", 2000, 1000);
  TCanvas *c_A_3 = new TCanvas("c_A_3", "A_3", 2000, 1000);
  c_A_1->Divide(2, 3);
  c_A_2->Divide(2, 3);
  c_A_3->Divide(2, 3);
  TCanvas *c_compare_1_2 = new TCanvas("c_compare_1_2", "1_2", 2000, 1000);
  TCanvas *c_compare_2_3 = new TCanvas("c_compare_2_3", "2_3", 2000, 1000);
  TCanvas *c_compare_1_3 = new TCanvas("c_compare_1_3", "1_3", 2000, 1000);
  c_compare_1_2->Divide(2, 3);
  c_compare_2_3->Divide(2, 3);
  c_compare_1_3->Divide(2, 3);
  for (int i = 0; i < 6; i++) {
    c_A_1->cd(i + 1);
    h_A_1_Z->at(i)->Draw();
    c_A_2->cd(i + 1);
    h_A_2_Z->at(i)->Draw();
    c_A_3->cd(i + 1);
    h_A_3_Z->at(i)->Draw();
    c_compare_1_2->cd(i + 1);
    g_compare_A_1_2_Z->at(i)->Draw("AP");
    c_compare_2_3->cd(i + 1);
    g_compare_A_2_3_Z->at(i)->Draw("AP");
    c_compare_1_3->cd(i + 1);
    g_compare_A_1_3_Z->at(i)->Draw("AP");
  }
}
//en->at(Z(j) - 1)->Fill(E_k(j) / ( (*MC_Mass)[MC_id_1] * ( (1 / sqrt(1 - pow(((*MC_Track_Length)[MC_id_1] * 1e-2/ ((*MC_TOF)[MC_id_1] * c)), 2) )) -1)) );  
      /*std::cout << "1 : " << ( (sqrt(pow((*MC_InitMom_x)[ MC_id_1],2) + pow((*MC_InitMom_y)[ MC_id_1],2)
       + pow((*MC_InitMom_z)[ MC_id_1],2) + pow((*MC_Mass)[MC_id_1], 2))) - (*MC_Mass)[MC_id_1]) << "2 " << ( (*MC_Mass)[MC_id_1] * ( (1 / sqrt(1 - pow(((*MC_Track_Length)[MC_id_1] * 1e-2/ ((*MC_TOF)[MC_id_1] * c)), 2) )) -1)) << endl;*/
/*std::cout << "p: " << p(j) << "p_: " << sqrt(pow((*MC_InitMom_x)[ MC_id_1],2) + pow((*MC_InitMom_y)[ MC_id_1],2)
       + pow((*MC_InitMom_z)[ MC_id_1],2)) << endl;*/