#define MyAnalysis_cxx
#include "MyAnalysisPrima.h"
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
double MyAnalysis::E_k(int &i) { //GeV
  return CAenergy->at(GLBtrackCAid->at(i)) +
         TWDe1Point->at(GLBtrackTWid->at(i)) * 1e-3 +
         TWDe2Point->at(GLBtrackTWid->at(i)) * 1e-3;
}
int MyAnalysis::Z(int &i) { return TWChargePoint->at(GLBtrackTWid->at(i)); }
// return ((CAenergy->at(GLBtrackCAid->at(i)) )* 1e3
// +TWDe1Point->at(GLBtrackTWid->at(i))+TWDe2Point->at(GLBtrackTWid->at(i)));
void MyAnalysis::A_Z(vector<TH1D *> *h_i_1, vector<TH1D *> *h_i_2,
                     vector<TH1D *> *h_i_3, vector<vector<double> *> *v_i_1,
                     vector<vector<double> *> *v_i_2,
                     vector<vector<double> *> *v_i_3) {
  for (int j = 0; j < GLBtrackTWid->size();
       j++) { // ipotizzando stessa dimensione vettori
    if (GLBtrackTWid->at(j) < 0 ||
        GLBtrackCAid->at(j) < 0) { // tolgo i casi in cui indice è negativo
      continue;
    }
    A_1 = (p(j) / (u * beta(j) * (1 / sqrt(1 - pow(beta(j), 2)))));
    h_i_1->at(Z(j) - 1)->Fill(A_1);
    v_i_1->at(Z(j) - 1)->push_back(A_1);

    A_2 = (E_k(j) / (u * ((1 / sqrt(1 - pow((beta(j)), 2))) - 1)));
    h_i_2->at(Z(j) - 1)->Fill(A_2);
    v_i_2->at(Z(j) - 1)->push_back(A_2);

    A_3 = ((pow(p(j), 2) - pow(E_k(j), 2)) / (2 * E_k(j)));
    h_i_3->at(Z(j) - 1)->Fill(A_3);
    v_i_3->at(Z(j) - 1)->push_back(A_3);
  }

  return;
}
void MyAnalysis::compare_A_Z(vector<TGraph *> *g_1_2_Z,
                             vector<TGraph *> *g_2_3_Z,
                             vector<TGraph *> *g_1_3_Z) {
  A_Z(h_A_1_Z, h_A_2_Z, h_A_3_Z, v_A_1_Z, v_A_2_Z, v_A_3_Z);
  for (int i = 0; i < v_A_1_Z->size(); i++) {
    g_1_2_Z->at(i)->Set(v_A_1_Z->at(i)->size());
    g_2_3_Z->at(i)->Set(v_A_2_Z->at(i)->size());
    g_1_3_Z->at(i)->Set(v_A_3_Z->at(i)->size());
    for (int j = 0; j < v_A_1_Z->at(i)->size(); j++) {
      (g_1_2_Z->at(i))->SetPoint(j, (v_A_1_Z->at(i))->at(j),
                               (v_A_2_Z->at(i))->at(j));
      (g_2_3_Z->at(i))->SetPoint(j, (v_A_2_Z->at(i))->at(j),
                               (v_A_3_Z->at(i))->at(j));
      (g_1_3_Z->at(i))->SetPoint(j, (v_A_1_Z->at(i))->at(j),
                               (v_A_3_Z->at(i))->at(j));
    }
  }
}
void MyAnalysis::Loop() {
  if (fChain == 0)
    return;

  //Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries = 20000;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    //A_Z(h_A_1_Z, h_A_2_Z, h_A_3_Z, v_A_1_Z, v_A_2_Z, v_A_3_Z);
    compare_A_Z(g_compare_A_1_2_Z, g_compare_A_2_3_Z, g_compare_A_1_3_Z);
  }
}
void MyAnalysis::BeforeLoop() {
  h_A_1_Z = new vector<TH1D *>(6);
  h_A_2_Z = new vector<TH1D *>(6);
  h_A_3_Z = new vector<TH1D *>(6);
  v_A_1_Z = new vector<vector<double> *>(6);
  v_A_2_Z = new vector<vector<double> *>(6);
  v_A_3_Z = new vector<vector<double> *>(6);
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
                              100, 0.0, 20);
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
                              100, 0.0, 20);
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
                              100, 0.0, 20);
    v_A_1_Z->at(i) = new vector<double>(0);
    v_A_2_Z->at(i) = new vector<double>(0);
    v_A_3_Z->at(i) = new vector<double>(0);
    std::ostringstream title;
    title << "Z = " << i+1;
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
  for (int i = 0; i < v_A_1_Z->size(); i++) {
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
void MyAnalysis::Analysis() {
  BeforeLoop();
  Loop();
  AfterLoop();
    for(int j = 0; j<6; j++) {
    std::cout << g_compare_A_1_2_Z->at(j)->GetN() << endl;
  }
}

/*void MyAnalysis::A(TH1D *h_1, TH1D *h_2, TH1D *h_3, vector<double> *v_1,
                   vector<double> *v_2, vector<double> *v_3) {
  for (int j = 0; j < GLBtrackTWid->size();
       j++) { // ipotizzando stessa dimensione vettori
    if (GLBtrackTWid->at(j) < 0 ||
        GLBtrackCAid->at(j) < 0) { // tolgo i casi in cui indice è negativo
      continue;
    }
    A_1 = (p(j) / (u * beta(j) * (1 / sqrt(1 - pow(beta(j), 2)))));
    // h_1->Fill(A_1);
    v_1->push_back(A_1);

    A_2 = (E_k(j) / (u * ((1 / sqrt(1 - pow((beta(j)), 2))) - 1)));
    // h_2->Fill(A_2);
    v_2->push_back(A_2);

    A_3 = ((pow(p(j), 2) - pow(E_k(j), 2)) / (2 * E_k(j)));
    // h_3->Fill(A_3);
    v_3->push_back(A_3);
  }

  return;
}
void MyAnalysis::compare_A(TGraph *g_1, TGraph *g_2, TGraph *g_3) {
  A(h_A_1, h_A_2, h_A_3, v_A_1, v_A_2, v_A_3);
  if (v_A_1->size() == v_A_2->size() && v_A_2->size() == v_A_3->size() &&
      v_A_1->size() == v_A_3->size()) {
    g_1->Set(v_A_1->size());
    g_2->Set(v_A_1->size());
    g_3->Set(v_A_1->size());
    for (int j = 0; j < v_A_1->size(); j++) {
      g_1->SetPoint(j, v_A_1->at(j), v_A_2->at(j));
      g_2->SetPoint(j, v_A_2->at(j), v_A_3->at(j));
      g_3->SetPoint(j, v_A_1->at(j), v_A_3->at(j));
    }
  } else {
    std::cout << "A1: " << v_A_1->size() << " A2: " << v_A_2->size()
              << " A3: " << v_A_3->size() << endl;
  }
}*/

  /*if (v_i_1->size() != 6) {
    std::cout << "errore vec: " << v_i_1->size();
  }
  if (h_i_2->size() != 6) {
    std::cout << "errore vec: " << h_i_2->size();
  }
  for (int i = 0; i < v_i_1->size(); i++) {
    if (v_i_1->at(i)->size() != v_i_1->at(i)->size() ||
        v_i_1->at(i)->size() != v_i_3->at(i)->size()) {
      std::cout << "errore size: " << v_i_1->at(i)->size()
                << v_i_2->at(i)->size() << v_i_3->at(i)->size();
    }
  }*/
