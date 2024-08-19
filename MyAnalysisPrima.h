#ifndef MyAnalysis_h
#define MyAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <TGraph.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class MyAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ev;
   Int_t           vertex_n; //NUMERO DI VERTEX POINT IN UN EVENTO(i vertex point sono la proiezione delle tracce sul target)                         
   vector<float>   *vertex_y;// ogni oggetto è un vector, perche a ogni evento
                           // sono associate tot grandezze
   vector<float>   *vertex_z; //piccato in 0, ma il target ha una dimensione fisica che è circa 5/600 MICROMETRI 
   vector<float>   *vertex_x; 
   vector<int>     *vt_trk_n;// questo è un vettore anche se solitamente ha un solo  elemento, che è = 1 :
                         // 1 traccia ricostruita per ogni punto (perche ho solo 1 interazione) 
                         //se ho piu elementi=1: piu vertex point, ognuno una traccia associata
   vector<float>   *vt_trk_chi2;
   vector<TVector3> *vt_trk_slopez;
   vector<TVector3> *vt_trk_origin;
   vector<TVector3> *vt_trk_projTW;
   vector<int>     *vt_trk_clus_n;//solitamente ha un solo elemento, e in quel caso è SEMPRE=4, cioè quatro cluster per 1 traccia
                              //altrimenti ho piu elementi, e quelli anche molto piu di 4 
//la somma di tutti gli elementi di vt_trk_n è il numero di tracce ricostruite. quindi la size di vt_trk_clus_n è quella somma
   vector<int>     *trk_vtx_clus_MCId;
   vector<int>     *vt_trk_clus_tot_hits;
   vector<float>   *vt_trk_clus_x;// questo considera solo i cluster con cui ho
                                // ricostruito traccia
   vector<float>   *vt_trk_clus_y;
   vector<float>   *vt_trk_clus_z;//qui ho solo 4 picchi e nient'altro perche relativo ai 4 piani del vertex
   vector<int>     *vt_clus_n;// numero di cluster per ogni piano(4 piani), quindi PER OGNI
                  // EVENTO (tracce da 1 particella) ho SEMPRE 4 ELEMENTI NEL
                  // VETTORE e ho 1 cluster se ho sono riuscito a costruire la
                  // traccia (quindi elementi del vettore sono 1 1 1 1)
                  // altrimenti 0 (0 0 0 0) infine puo essere che ci sono piu
                  // interazioni quindi piu cluster per piano dato 1 evento
   vector<int>     *vtx_clus_MCId;
   vector<int>     *vt_clus_tot_hits;// 1 cluster = tot pixel accesi, quindi per ogni piano
                         // ti dice, dato un evento, quanti pixel hittati
   vector<float>   *vt_clus_x;// posizione x di ogni cluster dall'inizio(anche se
                            // non ricostruisco traccia)
   vector<float>   *vt_clus_y;
   vector<float>   *vt_clus_z;
   vector<float>   *vt_clus_MC;
   vector<int>     *msd_station_id;
   vector<int>     *msd_pt_n;
   vector<int>     *msd_pt_MCId;
   vector<double>  *msd_eloss1;
   vector<double>  *msd_eloss2;
   vector<TVector3> *msd_pos;
   vector<int>     *msd_station_clus_id;
   vector<int>     *msd_clus_n;
   vector<int>     *msd_clus_MCId;
   vector<double>  *msd_eloss;
   vector<TVector3> *msd_clus_pos;
   Int_t           TWPoints;  //punto : interazione delle due barre. ne ho piu di 1 se ho piu frammenti
   vector<int>     *TWChargePoint;
   vector<int>     *TATW_MCID_1;
   vector<int>     *TATW_MCID_2;
   vector<double>  *TWDe1Point;
   vector<double>  *TWDe2Point;
   vector<double>  *TWXPoint;//valori discreti a seconda della barra, quindi x e y molto piccati in 0
   vector<double>  *TWYPoint;
   vector<double>  *TWZPoint;
   Int_t           TWHit; 
   vector<double>  *TWDeHit; // energia rilasciata sul TW, se nell'evento non si ha un rilascio di energia sensato
                             //il valore è -99
                             //valore medio è in keV
   vector<int>     *TW_MCID_hit;
   vector<double>  *TWTOFHit;
   vector<double>  *TWbarHit;
   vector<double>  *TWlayerHit;
   vector<double>  *TWTOF; //tempo medio è 10 NANOSECONDI
   vector<double>  *TWTOF1;
   vector<double>  *TWTOF2;
   Int_t           CAclusN;
   vector<double>  *CAenergy;
   vector<double>  *CAposX;
   vector<double>  *CAposY;
   vector<double>  *CAposZ;
   Int_t           GLBtracks;
   vector<double>  *GLBtrackPx;
   vector<double>  *GLBtrackPy;
   vector<double>  *GLBtrackPz;
   vector<double>  *GLBtrackLength;
   vector<int>     *GLBtrackTWid; //indice da guardare nei vettori di TW
   vector<int>     *GLBtrackCAid;
   vector<int>     *MC_Dead_region;
   vector<int>     *MC_Generation_region;
   vector<int>     *MC_FlukaID;
   vector<int>     *MC_MotherID;
   vector<int>     *MC_BaryonN;
   vector<double>  *MC_Mass;
   vector<int>     *MC_Ptype;
   vector<int>     *MC_Charge;
   vector<double>  *MC_TOF;
   vector<double>  *MC_Track_Length;
   vector<float>   *MC_InitPos_x;
   vector<float>   *MC_InitPos_y;
   vector<float>   *MC_InitPos_z;
   vector<float>   *MC_FinalPos_x;
   vector<float>   *MC_FianlPos_y;
   vector<float>   *MC_FinalPos_z;
   vector<float>   *MC_InitMom_x;
   vector<float>   *MC_InitMom_y;
   vector<float>   *MC_InitMom_z;
   vector<float>   *MC_FinalMom_x;
   vector<float>   *MC_FinalMom_y;
   vector<float>   *MC_FinalMom_z;

   // List of branches  : SERVONO PER LINKARE GLI OGGETTI PRECEDENTI CON IL DATA
   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_vertex_n;   //!
   TBranch        *b_vertex_x;   //!
   TBranch        *b_vertex_y;   //!
   TBranch        *b_vertex_z;   //!
   TBranch        *b_vt_trk_n;   //!
   TBranch        *b_vt_trk_chi2;   //!
   TBranch        *b_vt_trk_slopez;   //!
   TBranch        *b_vt_trk_origin;   //!
   TBranch        *b_vt_trk_projTW;   //!
   TBranch        *b_vt_trk_clus_n;   //!
   TBranch        *b_trk_vtx_clus_MCId;   //!
   TBranch        *b_vt_trk_clus_tot_hits;   //!
   TBranch        *b_vt_trk_clus_x;   //!
   TBranch        *b_vt_trk_clus_y;   //!
   TBranch        *b_vt_trk_clus_z;   //!
   TBranch        *b_vt_clus_n;   //!
   TBranch        *b_vtx_clus_MCId;   //!
   TBranch        *b_vt_clus_tot_hits;   //!
   TBranch        *b_vt_clus_x;   //!
   TBranch        *b_vt_clus_y;   //!
   TBranch        *b_vt_clus_z;   //!
   TBranch        *b_vt_clus_MC;   //!
   TBranch        *b_msd_station_id;   //!
   TBranch        *b_msd_pt_n;   //!
   TBranch        *b_msd_pt_MCId;   //!
   TBranch        *b_msd_eloss1;   //!
   TBranch        *b_msd_eloss2;   //!
   TBranch        *b_msd_pos;   //!
   TBranch        *b_msd_station_clus_id;   //!
   TBranch        *b_msd_clus_n;   //!
   TBranch        *b_msd_clus_MCId;   //!
   TBranch        *b_msd_eloss;   //!
   TBranch        *b_msd_clus_pos;   //!
   TBranch        *b_TWPoints;   //!
   TBranch        *b_TWChargePoint;   //!
   TBranch        *b_TATW_MCID_1;   //!
   TBranch        *b_TATW_MCID_2;   //!
   TBranch        *b_TWDe1Point;   //!
   TBranch        *b_TWDe2Point;   //!
   TBranch        *b_TWXPoint;   //!
   TBranch        *b_TWYPoint;   //!
   TBranch        *b_TWZPoint;   //!
   TBranch        *b_TWHit;   //!
   TBranch        *b_TWDeHit;   //!
   TBranch        *b_TW_MCID_hit;   //!
   TBranch        *b_TWTOFHit;   //!
   TBranch        *b_TWbarHit;   //!
   TBranch        *b_TWlayerHit;   //!
   TBranch        *b_TWTOF;   //!
   TBranch        *b_TWTOF1;   //!
   TBranch        *b_TWTOF2;   //!
   TBranch        *b_CAclusN;   //!
   TBranch        *b_CAenergy;   //!
   TBranch        *b_CAposX;   //!
   TBranch        *b_CAposY;   //!
   TBranch        *b_CAposZ;   //!
   TBranch        *b_GLBtracks;   //!
   TBranch        *b_GLBtrackPx;   //!
   TBranch        *b_GLBtrackPy;   //!
   TBranch        *b_GLBtrackPz;   //!
   TBranch        *b_GLBtrackLength;   //!
   TBranch        *b_GLBtrackTWid;   //!
   TBranch        *b_GLBtrackCAid;   //!
   TBranch        *b_MC_Dead_region;   //!
   TBranch        *b_MC_Generation_region;   //!
   TBranch        *b_MC_FlukaID;   //!
   TBranch        *b_MC_MotherID;   //!
   TBranch        *b_MC_BaryonN;   //!
   TBranch        *b_MC_Mass;   //!
   TBranch        *b_MC_Ptype;   //!
   TBranch        *b_MC_Charge;   //!
   TBranch        *b_MC_TOF;   //!
   TBranch        *b_MC_Track_Length;   //!
   TBranch        *b_MC_InitPos_x;   //!
   TBranch        *b_MC_InitPos_y;   //!
   TBranch        *b_MC_InitPos_z;   //!
   TBranch        *b_MC_FinalPos_x;   //!
   TBranch        *b_MC_FianlPos_y;   //!
   TBranch        *b_MC_FinalPos_z;   //!
   TBranch        *b_MC_InitMom_x;   //!
   TBranch        *b_MC_InitMom_y;   //!
   TBranch        *b_MC_InitMom_z;   //!
   TBranch        *b_MC_FinalMom_x;   //!
   TBranch        *b_MC_FinalMom_y;   //!
   TBranch        *b_MC_FinalMom_z;   //!

   MyAnalysis(TTree *tree=0);
   virtual ~MyAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void BeforeLoop();
   virtual void Loop();
   virtual void AfterLoop();
   virtual void Analysis();
   virtual double t(int& i);
   virtual double p(int& i);
   virtual double beta(int& i);
   virtual double E_k(int& i);
   virtual int Z(int& i);
   //virtual void A(TH1D* h_1, TH1D* h_2, TH1D* h_3, vector<double> *v_1, vector<double> *v_2, vector<double> *v_3);
   //virtual void compare_A(TGraph *g_1,  TGraph *g_2,  TGraph *g_3);
   virtual void A_Z(vector<TH1D*> *h_i_1, vector<TH1D*> *h_i_2, vector<TH1D*> *h_i_3, vector<vector<double>*> *v_i_1, vector<vector<double>*> *v_i_2, vector<vector<double>*> *v_i_3);
   virtual void compare_A_Z(vector<TGraph*> *g_1_Z,  vector<TGraph*> *g_2_Z,  vector<TGraph*> *g_3_Z);
   
   /*virtual void A_Z_ex(TH1D *h_1, TH1D *h_2, TH1D *h_3, TH1D *h_4, TH1D *h_5,
TH1D *h_6, vector<double> *v_1, vector<double> *v_2, vector<double> *v_3,
 vector<double> *v_4, vector<double> *v_5, vector<double> *v_6);*/
   private:
   double u = 0.9314936148; //  [GeV/c^2]
   double E_k_in = 0.400; //  [Gev/u]  (è Carbonio12)
   double beta_in = sqrt(1-1/pow(1+E_k_in/u,2));
   double c = 299792458.0; 
   double dx = 0.45925; // [m] 
   double dt = dx/(beta_in*c); //  [s]
   
   double A_1;
   double A_2;
   double A_3;

   TH1D *h_A_1;
   TH1D *h_A_2;
   TH1D *h_A_3;
   vector<double> *v_A_1;
   vector<double> *v_A_2;
   vector<double> *v_A_3;
   TGraph *g_compare_A_1_2;
   TGraph *g_compare_A_2_3;
   TGraph *g_compare_A_1_3;

   vector<TH1D*> *h_A_1_Z;
   vector<TH1D*> *h_A_2_Z;
   vector<TH1D*> *h_A_3_Z;
   vector<vector<double>*> *v_A_1_Z;
   vector<vector<double>*> *v_A_2_Z;
   vector<vector<double>*> *v_A_3_Z;
   vector<TGraph*> *g_compare_A_1_2_Z;
   vector<TGraph*> *g_compare_A_2_3_Z;
   vector<TGraph*> *g_compare_A_1_3_Z; 

   /*TH1D *h_Z_1;
   TH1D *h_Z_2;
   TH1D *h_Z_3;
   TH1D *h_Z_4;
   TH1D *h_Z_5;
   TH1D *h_Z_6;
   vector<double> *v_Z_1;
   vector<double> *v_Z_2;
   vector<double> *v_Z_3;
   vector<double> *v_Z_4;
   vector<double> *v_Z_5;
   vector<double> *v_Z_6;*/

};

#endif

#ifdef MyAnalysis_cxx
MyAnalysis::MyAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("12C_400_C.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("12C_400_C.root");
      }
      f->GetObject("OuTree",tree);

   }
   Init(tree);
}

MyAnalysis::~MyAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vertex_x = 0;
   vertex_y = 0;
   vertex_z = 0;
   vt_trk_n = 0;
   vt_trk_chi2 = 0;
   vt_trk_slopez = 0;
   vt_trk_origin = 0;
   vt_trk_projTW = 0;
   vt_trk_clus_n = 0;
   trk_vtx_clus_MCId = 0;
   vt_trk_clus_tot_hits = 0;
   vt_trk_clus_x = 0;
   vt_trk_clus_y = 0;
   vt_trk_clus_z = 0;
   vt_clus_n = 0;
   vtx_clus_MCId = 0;
   vt_clus_tot_hits = 0;
   vt_clus_x = 0;
   vt_clus_y = 0;
   vt_clus_z = 0;
   vt_clus_MC = 0;
   msd_station_id = 0;
   msd_pt_n = 0;
   msd_pt_MCId = 0;
   msd_eloss1 = 0;
   msd_eloss2 = 0;
   msd_pos = 0;
   msd_station_clus_id = 0;
   msd_clus_n = 0;
   msd_clus_MCId = 0;
   msd_eloss = 0;
   msd_clus_pos = 0;
   TWChargePoint = 0;
   TATW_MCID_1 = 0;
   TATW_MCID_2 = 0;
   TWDe1Point = 0;
   TWDe2Point = 0;
   TWXPoint = 0;
   TWYPoint = 0;
   TWZPoint = 0;
   TWDeHit = 0;
   TW_MCID_hit = 0;
   TWTOFHit = 0;
   TWbarHit = 0;
   TWlayerHit = 0;
   TWTOF = 0;
   TWTOF1 = 0;
   TWTOF2 = 0;
   CAenergy = 0;
   CAposX = 0;
   CAposY = 0;
   CAposZ = 0;
   GLBtrackPx = 0;
   GLBtrackPy = 0;
   GLBtrackPz = 0;
   GLBtrackLength = 0;
   GLBtrackTWid = 0;
   GLBtrackCAid = 0;
   MC_Dead_region = 0;
   MC_Generation_region = 0;
   MC_FlukaID = 0;
   MC_MotherID = 0;
   MC_BaryonN = 0;
   MC_Mass = 0;
   MC_Ptype = 0;
   MC_Charge = 0;
   MC_TOF = 0;
   MC_Track_Length = 0;
   MC_InitPos_x = 0;
   MC_InitPos_y = 0;
   MC_InitPos_z = 0;
   MC_FinalPos_x = 0;
   MC_FianlPos_y = 0;
   MC_FinalPos_z = 0;
   MC_InitMom_x = 0;
   MC_InitMom_y = 0;
   MC_InitMom_z = 0;
   MC_FinalMom_x = 0;
   MC_FinalMom_y = 0;
   MC_FinalMom_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev, &b_ev);
   fChain->SetBranchAddress("vertex_n", &vertex_n, &b_vertex_n);
   fChain->SetBranchAddress("vertex_x", &vertex_x, &b_vertex_x);
   fChain->SetBranchAddress("vertex_y", &vertex_y, &b_vertex_y);
   fChain->SetBranchAddress("vertex_z", &vertex_z, &b_vertex_z);
   fChain->SetBranchAddress("vt_trk_n", &vt_trk_n, &b_vt_trk_n);
   fChain->SetBranchAddress("vt_trk_chi2", &vt_trk_chi2, &b_vt_trk_chi2);
   fChain->SetBranchAddress("vt_trk_slopez", &vt_trk_slopez, &b_vt_trk_slopez);
   fChain->SetBranchAddress("vt_trk_origin", &vt_trk_origin, &b_vt_trk_origin);
   fChain->SetBranchAddress("vt_trk_projTW", &vt_trk_projTW, &b_vt_trk_projTW);
   fChain->SetBranchAddress("vt_trk_clus_n", &vt_trk_clus_n, &b_vt_trk_clus_n);
   fChain->SetBranchAddress("trk_vtx_clus_MCId", &trk_vtx_clus_MCId, &b_trk_vtx_clus_MCId);
   fChain->SetBranchAddress("vt_trk_clus_tot_hits", &vt_trk_clus_tot_hits, &b_vt_trk_clus_tot_hits);
   fChain->SetBranchAddress("vt_trk_clus_x", &vt_trk_clus_x, &b_vt_trk_clus_x);
   fChain->SetBranchAddress("vt_trk_clus_y", &vt_trk_clus_y, &b_vt_trk_clus_y);
   fChain->SetBranchAddress("vt_trk_clus_z", &vt_trk_clus_z, &b_vt_trk_clus_z);
   fChain->SetBranchAddress("vt_clus_n", &vt_clus_n, &b_vt_clus_n);
   fChain->SetBranchAddress("vtx_clus_MCId", &vtx_clus_MCId, &b_vtx_clus_MCId);
   fChain->SetBranchAddress("vt_clus_tot_hits", &vt_clus_tot_hits, &b_vt_clus_tot_hits);
   fChain->SetBranchAddress("vt_clus_x", &vt_clus_x, &b_vt_clus_x);
   fChain->SetBranchAddress("vt_clus_y", &vt_clus_y, &b_vt_clus_y);
   fChain->SetBranchAddress("vt_clus_z", &vt_clus_z, &b_vt_clus_z);
   fChain->SetBranchAddress("vt_clus_MC", &vt_clus_MC, &b_vt_clus_MC);
   fChain->SetBranchAddress("msd_station_id", &msd_station_id, &b_msd_station_id);
   fChain->SetBranchAddress("msd_pt_n", &msd_pt_n, &b_msd_pt_n);
   fChain->SetBranchAddress("msd_pt_MCId", &msd_pt_MCId, &b_msd_pt_MCId);
   fChain->SetBranchAddress("msd_eloss1", &msd_eloss1, &b_msd_eloss1);
   fChain->SetBranchAddress("msd_eloss2", &msd_eloss2, &b_msd_eloss2);
   fChain->SetBranchAddress("msd_pos", &msd_pos, &b_msd_pos);
   fChain->SetBranchAddress("msd_station_clus_id", &msd_station_clus_id, &b_msd_station_clus_id);
   fChain->SetBranchAddress("msd_clus_n", &msd_clus_n, &b_msd_clus_n);
   fChain->SetBranchAddress("msd_clus_MCId", &msd_clus_MCId, &b_msd_clus_MCId);
   fChain->SetBranchAddress("msd_eloss", &msd_eloss, &b_msd_eloss);
   fChain->SetBranchAddress("msd_clus_pos", &msd_clus_pos, &b_msd_clus_pos);
   fChain->SetBranchAddress("TWPoints", &TWPoints, &b_TWPoints);
   fChain->SetBranchAddress("TWChargePoint", &TWChargePoint, &b_TWChargePoint);
   fChain->SetBranchAddress("TATW_MCID_1", &TATW_MCID_1, &b_TATW_MCID_1);
   fChain->SetBranchAddress("TATW_MCID_2", &TATW_MCID_2, &b_TATW_MCID_2);
   fChain->SetBranchAddress("TWDe1Point", &TWDe1Point, &b_TWDe1Point);
   fChain->SetBranchAddress("TWDe2Point", &TWDe2Point, &b_TWDe2Point);
   fChain->SetBranchAddress("TWXPoint", &TWXPoint, &b_TWXPoint);
   fChain->SetBranchAddress("TWYPoint", &TWYPoint, &b_TWYPoint);
   fChain->SetBranchAddress("TWZPoint", &TWZPoint, &b_TWZPoint);
   fChain->SetBranchAddress("TWHit", &TWHit, &b_TWHit);
   fChain->SetBranchAddress("TWDeHit", &TWDeHit, &b_TWDeHit);
   fChain->SetBranchAddress("TW_MCID_hit", &TW_MCID_hit, &b_TW_MCID_hit);
   fChain->SetBranchAddress("TWTOFHit", &TWTOFHit, &b_TWTOFHit);
   fChain->SetBranchAddress("TWbarHit", &TWbarHit, &b_TWbarHit);
   fChain->SetBranchAddress("TWlayerHit", &TWlayerHit, &b_TWlayerHit);
   fChain->SetBranchAddress("TWTOF", &TWTOF, &b_TWTOF);
   fChain->SetBranchAddress("TWTOF1", &TWTOF1, &b_TWTOF1);
   fChain->SetBranchAddress("TWTOF2", &TWTOF2, &b_TWTOF2);
   fChain->SetBranchAddress("CAclusN", &CAclusN, &b_CAclusN);
   fChain->SetBranchAddress("CAenergy", &CAenergy, &b_CAenergy);
   fChain->SetBranchAddress("CAposX", &CAposX, &b_CAposX);
   fChain->SetBranchAddress("CAposY", &CAposY, &b_CAposY);
   fChain->SetBranchAddress("CAposZ", &CAposZ, &b_CAposZ);
   fChain->SetBranchAddress("GLBtracks", &GLBtracks, &b_GLBtracks);
   fChain->SetBranchAddress("GLBtrackPx", &GLBtrackPx, &b_GLBtrackPx);
   fChain->SetBranchAddress("GLBtrackPy", &GLBtrackPy, &b_GLBtrackPy);
   fChain->SetBranchAddress("GLBtrackPz", &GLBtrackPz, &b_GLBtrackPz);
   fChain->SetBranchAddress("GLBtrackLength", &GLBtrackLength, &b_GLBtrackLength);
   fChain->SetBranchAddress("GLBtrackTWid", &GLBtrackTWid, &b_GLBtrackTWid);
   fChain->SetBranchAddress("GLBtrackCAid", &GLBtrackCAid, &b_GLBtrackCAid);
   fChain->SetBranchAddress("MC_Dead_region", &MC_Dead_region, &b_MC_Dead_region);
   fChain->SetBranchAddress("MC_Generation_region", &MC_Generation_region, &b_MC_Generation_region);
   fChain->SetBranchAddress("MC_FlukaID", &MC_FlukaID, &b_MC_FlukaID);
   fChain->SetBranchAddress("MC_MotherID", &MC_MotherID, &b_MC_MotherID);
   fChain->SetBranchAddress("MC_BaryonN", &MC_BaryonN, &b_MC_BaryonN);
   fChain->SetBranchAddress("MC_Mass", &MC_Mass, &b_MC_Mass);
   fChain->SetBranchAddress("MC_Ptype", &MC_Ptype, &b_MC_Ptype);
   fChain->SetBranchAddress("MC_Charge", &MC_Charge, &b_MC_Charge);
   fChain->SetBranchAddress("MC_TOF", &MC_TOF, &b_MC_TOF);
   fChain->SetBranchAddress("MC_Track_Length", &MC_Track_Length, &b_MC_Track_Length);
   fChain->SetBranchAddress("MC_InitPos_x", &MC_InitPos_x, &b_MC_InitPos_x);
   fChain->SetBranchAddress("MC_InitPos_y", &MC_InitPos_y, &b_MC_InitPos_y);
   fChain->SetBranchAddress("MC_InitPos_z", &MC_InitPos_z, &b_MC_InitPos_z);
   fChain->SetBranchAddress("MC_FinalPos_x", &MC_FinalPos_x, &b_MC_FinalPos_x);
   fChain->SetBranchAddress("MC_FianlPos_y", &MC_FianlPos_y, &b_MC_FianlPos_y);
   fChain->SetBranchAddress("MC_FinalPos_z", &MC_FinalPos_z, &b_MC_FinalPos_z);
   fChain->SetBranchAddress("MC_InitMom_x", &MC_InitMom_x, &b_MC_InitMom_x);
   fChain->SetBranchAddress("MC_InitMom_y", &MC_InitMom_y, &b_MC_InitMom_y);
   fChain->SetBranchAddress("MC_InitMom_z", &MC_InitMom_z, &b_MC_InitMom_z);
   fChain->SetBranchAddress("MC_FinalMom_x", &MC_FinalMom_x, &b_MC_FinalMom_x);
   fChain->SetBranchAddress("MC_FinalMom_y", &MC_FinalMom_y, &b_MC_FinalMom_y);
   fChain->SetBranchAddress("MC_FinalMom_z", &MC_FinalMom_z, &b_MC_FinalMom_z);
   Notify();
}

Bool_t MyAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyAnalysis_cxx
