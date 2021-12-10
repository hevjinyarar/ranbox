//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  2 13:41:50 2021 by ROOT version 6.24/04
// from TTree tree/tree
// found on file: bparked_Alessandro.root
//////////////////////////////////////////////////////////

#ifndef bpark_h
#define bpark_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class bpark {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         rho;
   vector<float>   *METraw_et;
   vector<float>   *METraw_phi;
   vector<float>   *METraw_sumEt;
   vector<float>   *MET_corrPx;
   vector<float>   *MET_corrPy;
   vector<float>   *MET_et;
   vector<float>   *MET_phi;
   vector<float>   *MET_puppi_et;
   vector<float>   *MET_puppi_phi;
   vector<float>   *MET_sumEt;
   vector<float>   *MET_JetEnUp;
   vector<float>   *MET_JetEnDown;
   vector<float>   *MET_JetResUp;
   vector<float>   *MET_JetResDown;
   vector<float>   *MET_UnclusteredEnUp;
   vector<float>   *MET_UnclusteredEnDown;
   Int_t           EVENT_event;
   Int_t           EVENT_run;
   Int_t           EVENT_lumiBlock;
   Int_t           PV_N;
   Bool_t          PV_filter;
   vector<float>   *PV_chi2;
   vector<float>   *PV_ndof;
   vector<float>   *PV_rho;
   vector<float>   *PV_z;
   vector<float>   *BeamSpot_x0;
   vector<float>   *BeamSpot_y0;
   vector<float>   *BeamSpot_z0;
   vector<int>     *IsBsTauTau;
   vector<int>     *BsTauTau_nCandidates;
   vector<float>   *BsTauTau_mu1_pt;
   vector<float>   *BsTauTau_mu1_eta;
   vector<float>   *BsTauTau_mu1_phi;
   vector<float>   *BsTauTau_mu1_mass;
   vector<int>     *BsTauTau_mu1_q;
   vector<int>     *BsTauTau_mu1_isLoose;
   vector<int>     *BsTauTau_mu1_isTight;
   vector<int>     *BsTauTau_mu1_isPF;
   vector<int>     *BsTauTau_mu1_isGlobal;
   vector<int>     *BsTauTau_mu1_isTracker;
   vector<int>     *BsTauTau_mu1_isSoft;
   vector<float>   *BsTauTau_mu1_vx;
   vector<float>   *BsTauTau_mu1_vy;
   vector<float>   *BsTauTau_mu1_vz;
   vector<float>   *BsTauTau_mu1_iso;
   vector<float>   *BsTauTau_mu1_dbiso;
   vector<float>   *BsTauTau_tau_pt;
   vector<float>   *BsTauTau_tau_eta;
   vector<float>   *BsTauTau_tau_phi;
   vector<float>   *BsTauTau_tau_mass;
   vector<float>   *BsTauTau_tau_rhomass1;
   vector<float>   *BsTauTau_tau_rhomass2;
   vector<int>     *BsTauTau_tau_q;
   vector<float>   *BsTauTau_tau_vx;
   vector<float>   *BsTauTau_tau_vy;
   vector<float>   *BsTauTau_tau_vz;
   vector<float>   *BsTauTau_tau_max_dr_3prong;
   vector<float>   *BsTauTau_tau_lip;
   vector<float>   *BsTauTau_tau_lips;
   vector<float>   *BsTauTau_tau_pvip;
   vector<float>   *BsTauTau_tau_pvips;
   vector<float>   *BsTauTau_tau_fl3d;
   vector<float>   *BsTauTau_tau_fls3d;
   vector<float>   *BsTauTau_tau_alpha;
   vector<float>   *BsTauTau_tau_vprob;
   vector<bool>    *BsTauTau_tau_isRight;
   vector<bool>    *BsTauTau_tau_isRight1;
   vector<bool>    *BsTauTau_tau_isRight2;
   vector<bool>    *BsTauTau_tau_isRight3;
   vector<float>   *BsTauTau_tau_dr1;
   vector<float>   *BsTauTau_tau_dr2;
   vector<float>   *BsTauTau_tau_dr3;
   vector<float>   *BsTauTau_tau_ptres1;
   vector<float>   *BsTauTau_tau_ptres2;
   vector<float>   *BsTauTau_tau_ptres3;
   vector<int>     *BsTauTau_tau_matched_ppdgId;
   vector<float>   *BsTauTau_tau_matched_gentaupt;
   vector<float>   *BsTauTau_tau_gentaupt;
   vector<float>   *BsTauTau_tau_sumofdnn;
   vector<int>     *BsTauTau_tau_pfidx1;
   vector<int>     *BsTauTau_tau_pfidx2;
   vector<int>     *BsTauTau_tau_pfidx3;
   vector<float>   *BsTauTau_tau_pi1_dnn;
   vector<float>   *BsTauTau_tau_pi2_dnn;
   vector<float>   *BsTauTau_tau_pi3_dnn;
   vector<float>   *BsTauTau_tau_pi1_pt;
   vector<float>   *BsTauTau_tau_pi1_eta;
   vector<float>   *BsTauTau_tau_pi1_phi;
   vector<float>   *BsTauTau_tau_pi1_mass;
   vector<int>     *BsTauTau_tau_pi1_charge;
   vector<float>   *BsTauTau_tau_pi1_dz;
   vector<float>   *BsTauTau_tau_muon_dr1;
   map<string,bool> *HLT_BPH_isFired;
   vector<float>   *BsTauTau_tau_pi2_pt;
   vector<float>   *BsTauTau_tau_pi2_eta;
   vector<float>   *BsTauTau_tau_pi2_phi;
   vector<float>   *BsTauTau_tau_pi2_mass;
   vector<int>     *BsTauTau_tau_pi2_charge;
   vector<float>   *BsTauTau_tau_pi2_dz;
   vector<float>   *BsTauTau_tau_muon_dr2;
   vector<float>   *BsTauTau_tau_pi3_pt;
   vector<float>   *BsTauTau_tau_pi3_eta;
   vector<float>   *BsTauTau_tau_pi3_phi;
   vector<float>   *BsTauTau_tau_pi3_mass;
   vector<int>     *BsTauTau_tau_pi3_charge;
   vector<float>   *BsTauTau_tau_pi3_dz;
   vector<float>   *BsTauTau_tau_muon_dr3;
   vector<float>   *BsTauTau_PV_vx;
   vector<float>   *BsTauTau_PV_vy;
   vector<float>   *BsTauTau_PV_vz;
   vector<float>   *BsTauTau_bbPV_vx;
   vector<float>   *BsTauTau_bbPV_vy;
   vector<float>   *BsTauTau_bbPV_vz;
   vector<float>   *BsTauTau_bbPV_refit_vx;
   vector<float>   *BsTauTau_bbPV_refit_vy;
   vector<float>   *BsTauTau_bbPV_refit_vz;
   vector<float>   *BsTauTau_genPV_vx;
   vector<float>   *BsTauTau_genPV_vy;
   vector<float>   *BsTauTau_genPV_vz;
   vector<float>   *BsTauTau_B_pt;
   vector<float>   *BsTauTau_B_eta;
   vector<float>   *BsTauTau_B_phi;
   vector<float>   *BsTauTau_B_mass;
   vector<float>   *BsTauTau_B_vprob;
   vector<float>   *BsTauTau_B_lip;
   vector<float>   *BsTauTau_B_lips;
   vector<float>   *BsTauTau_B_pvip;
   vector<float>   *BsTauTau_B_pvips;
   vector<float>   *BsTauTau_B_fl3d;
   vector<float>   *BsTauTau_B_fls3d;
   vector<float>   *BsTauTau_B_alpha;
   vector<float>   *BsTauTau_B_maxdoca;
   vector<float>   *BsTauTau_B_mindoca;
   vector<float>   *BsTauTau_muonpion_maxdoca;
   vector<float>   *BsTauTau_muonpion_mindoca;
   vector<float>   *BsTauTau_B_vx;
   vector<float>   *BsTauTau_B_vy;
   vector<float>   *BsTauTau_B_vz;
   vector<float>   *BsTauTau_B_iso;
   vector<int>     *BsTauTau_B_iso_ntracks;
   vector<float>   *BsTauTau_B_iso_mindoca;
   vector<int>     *BsTauTau_ngenmuons;
   vector<int>     *BsTauTau_isgen3;
   vector<int>     *BsTauTau_isgen3matched;
   vector<int>     *BsTauTau_nch;
   vector<int>     *BsTauTau_nch_after_dnn;
   vector<int>     *BsTauTau_nch_before_dnn;
   vector<int>     *BsTauTau_nch_qr;
   vector<int>     *BsTauTau_ngentau3;
   vector<int>     *BsTauTau_ngentau;
   vector<float>   *BsTauTau_gentaupt;
   vector<int>     *BsTauTau_gentaudm;

   // List of branches
   TBranch        *b_rho;   //!
   TBranch        *b_METraw_et;   //!
   TBranch        *b_METraw_phi;   //!
   TBranch        *b_METraw_sumEt;   //!
   TBranch        *b_MET_corrPx;   //!
   TBranch        *b_MET_corrPy;   //!
   TBranch        *b_MET_et;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_puppi_et;   //!
   TBranch        *b_MET_puppi_phi;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_MET_JetEnUp;   //!
   TBranch        *b_MET_JetEnDown;   //!
   TBranch        *b_MET_JetResUp;   //!
   TBranch        *b_MET_JetResDown;   //!
   TBranch        *b_MET_UnclusteredEnUp;   //!
   TBranch        *b_MET_UnclusteredEnDown;   //!
   TBranch        *b_EVENT_event;   //!
   TBranch        *b_EVENT_run;   //!
   TBranch        *b_EVENT_lumiBlock;   //!
   TBranch        *b_PV_N;   //!
   TBranch        *b_PV_filter;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_rho;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_BeamSpot_x0;   //!
   TBranch        *b_BeamSpot_y0;   //!
   TBranch        *b_BeamSpot_z0;   //!
   TBranch        *b_IsBsTauTau;   //!
   TBranch        *b_BsTauTau_nCandidates;   //!
   TBranch        *b_BsTauTau_mu1_pt;   //!
   TBranch        *b_BsTauTau_mu1_eta;   //!
   TBranch        *b_BsTauTau_mu1_phi;   //!
   TBranch        *b_BsTauTau_mu1_mass;   //!
   TBranch        *b_BsTauTau_mu1_q;   //!
   TBranch        *b_BsTauTau_mu1_isLoose;   //!
   TBranch        *b_BsTauTau_mu1_isTight;   //!
   TBranch        *b_BsTauTau_mu1_isPF;   //!
   TBranch        *b_BsTauTau_mu1_isGlobal;   //!
   TBranch        *b_BsTauTau_mu1_isTracker;   //!
   TBranch        *b_BsTauTau_mu1_isSoft;   //!
   TBranch        *b_BsTauTau_mu1_vx;   //!
   TBranch        *b_BsTauTau_mu1_vy;   //!
   TBranch        *b_BsTauTau_mu1_vz;   //!
   TBranch        *b_BsTauTau_mu1_iso;   //!
   TBranch        *b_BsTauTau_mu1_dbiso;   //!
   TBranch        *b_BsTauTau_tau_pt;   //!
   TBranch        *b_BsTauTau_tau_eta;   //!
   TBranch        *b_BsTauTau_tau_phi;   //!
   TBranch        *b_BsTauTau_tau_mass;   //!
   TBranch        *b_BsTauTau_tau_rhomass1;   //!
   TBranch        *b_BsTauTau_tau_rhomass2;   //!
   TBranch        *b_BsTauTau_tau_q;   //!
   TBranch        *b_BsTauTau_tau_vx;   //!
   TBranch        *b_BsTauTau_tau_vy;   //!
   TBranch        *b_BsTauTau_tau_vz;   //!
   TBranch        *b_BsTauTau_tau_max_dr_3prong;   //!
   TBranch        *b_BsTauTau_tau_lip;   //!
   TBranch        *b_BsTauTau_tau_lips;   //!
   TBranch        *b_BsTauTau_tau_pvip;   //!
   TBranch        *b_BsTauTau_tau_pvips;   //!
   TBranch        *b_BsTauTau_tau_fl3d;   //!
   TBranch        *b_BsTauTau_tau_fls3d;   //!
   TBranch        *b_BsTauTau_tau_alpha;   //!
   TBranch        *b_BsTauTau_tau_vprob;   //!
   TBranch        *b_BsTauTau_tau_isRight;   //!
   TBranch        *b_BsTauTau_tau_isRight1;   //!
   TBranch        *b_BsTauTau_tau_isRight2;   //!
   TBranch        *b_BsTauTau_tau_isRight3;   //!
   TBranch        *b_BsTauTau_tau_dr1;   //!
   TBranch        *b_BsTauTau_tau_dr2;   //!
   TBranch        *b_BsTauTau_tau_dr3;   //!
   TBranch        *b_BsTauTau_tau_ptres1;   //!
   TBranch        *b_BsTauTau_tau_ptres2;   //!
   TBranch        *b_BsTauTau_tau_ptres3;   //!
   TBranch        *b_BsTauTau_tau_matched_ppdgId;   //!
   TBranch        *b_BsTauTau_tau_matched_gentaupt;   //!
   TBranch        *b_BsTauTau_tau_gentaupt;   //!
   TBranch        *b_BsTauTau_tau_sumofdnn;   //!
   TBranch        *b_BsTauTau_tau_pfidx1;   //!
   TBranch        *b_BsTauTau_tau_pfidx2;   //!
   TBranch        *b_BsTauTau_tau_pfidx3;   //!
   TBranch        *b_BsTauTau_tau_pi1_dnn;   //!
   TBranch        *b_BsTauTau_tau_pi2_dnn;   //!
   TBranch        *b_BsTauTau_tau_pi3_dnn;   //!
   TBranch        *b_BsTauTau_tau_pi1_pt;   //!
   TBranch        *b_BsTauTau_tau_pi1_eta;   //!
   TBranch        *b_BsTauTau_tau_pi1_phi;   //!
   TBranch        *b_BsTauTau_tau_pi1_mass;   //!
   TBranch        *b_BsTauTau_tau_pi1_charge;   //!
   TBranch        *b_BsTauTau_tau_pi1_dz;   //!
   TBranch        *b_BsTauTau_tau_muon_dr1;   //!
   TBranch        *b_HLT_BPH_isFired;   //!
   TBranch        *b_BsTauTau_tau_pi2_pt;   //!
   TBranch        *b_BsTauTau_tau_pi2_eta;   //!
   TBranch        *b_BsTauTau_tau_pi2_phi;   //!
   TBranch        *b_BsTauTau_tau_pi2_mass;   //!
   TBranch        *b_BsTauTau_tau_pi2_charge;   //!
   TBranch        *b_BsTauTau_tau_pi2_dz;   //!
   TBranch        *b_BsTauTau_tau_muon_dr2;   //!
   TBranch        *b_BsTauTau_tau_pi3_pt;   //!
   TBranch        *b_BsTauTau_tau_pi3_eta;   //!
   TBranch        *b_BsTauTau_tau_pi3_phi;   //!
   TBranch        *b_BsTauTau_tau_pi3_mass;   //!
   TBranch        *b_BsTauTau_tau_pi3_charge;   //!
   TBranch        *b_BsTauTau_tau_pi3_dz;   //!
   TBranch        *b_BsTauTau_tau_muon_dr3;   //!
   TBranch        *b_BsTauTau_PV_vx;   //!
   TBranch        *b_BsTauTau_PV_vy;   //!
   TBranch        *b_BsTauTau_PV_vz;   //!
   TBranch        *b_BsTauTau_bbPV_vx;   //!
   TBranch        *b_BsTauTau_bbPV_vy;   //!
   TBranch        *b_BsTauTau_bbPV_vz;   //!
   TBranch        *b_BsTauTau_bbPV_refit_vx;   //!
   TBranch        *b_BsTauTau_bbPV_refit_vy;   //!
   TBranch        *b_BsTauTau_bbPV_refit_vz;   //!
   TBranch        *b_BsTauTau_genPV_vx;   //!
   TBranch        *b_BsTauTau_genPV_vy;   //!
   TBranch        *b_BsTauTau_genPV_vz;   //!
   TBranch        *b_BsTauTau_B_pt;   //!
   TBranch        *b_BsTauTau_B_eta;   //!
   TBranch        *b_BsTauTau_B_phi;   //!
   TBranch        *b_BsTauTau_B_mass;   //!
   TBranch        *b_BsTauTau_B_vprob;   //!
   TBranch        *b_BsTauTau_B_lip;   //!
   TBranch        *b_BsTauTau_B_lips;   //!
   TBranch        *b_BsTauTau_B_pvip;   //!
   TBranch        *b_BsTauTau_B_pvips;   //!
   TBranch        *b_BsTauTau_B_fl3d;   //!
   TBranch        *b_BsTauTau_B_fls3d;   //!
   TBranch        *b_BsTauTau_B_alpha;   //!
   TBranch        *b_BsTauTau_B_maxdoca;   //!
   TBranch        *b_BsTauTau_B_mindoca;   //!
   TBranch        *b_BsTauTau_muonpion_maxdoca;   //!
   TBranch        *b_BsTauTau_muonpion_mindoca;   //!
   TBranch        *b_BsTauTau_B_vx;   //!
   TBranch        *b_BsTauTau_B_vy;   //!
   TBranch        *b_BsTauTau_B_vz;   //!
   TBranch        *b_BsTauTau_B_iso;   //!
   TBranch        *b_BsTauTau_B_iso_ntracks;   //!
   TBranch        *b_BsTauTau_B_iso_mindoca;   //!
   TBranch        *b_BsTauTau_ngenmuons;   //!
   TBranch        *b_BsTauTau_isgen3;   //!
   TBranch        *b_BsTauTau_isgen3matched;   //!
   TBranch        *b_BsTauTau_nch;   //!
   TBranch        *b_BsTauTau_nch_after_dnn;   //!
   TBranch        *b_BsTauTau_nch_before_dnn;   //!
   TBranch        *b_BsTauTau_nch_qr;   //!
   TBranch        *b_BsTauTau_ngentau3;   //!
   TBranch        *b_BsTauTau_ngentau;   //!
   TBranch        *b_BsTauTau_gentaupt;   //!
   TBranch        *b_BsTauTau_gentaudm;   //!

   bpark(TTree *tree=0);
   virtual ~bpark();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef bpark_cxx
bpark::bpark(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("bparked_Alessandro.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("bparked_Alessandro.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("bparked_Alessandro.root:/ntuplizer");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

bpark::~bpark()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bpark::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bpark::LoadTree(Long64_t entry)
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

void bpark::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   METraw_et = 0;
   METraw_phi = 0;
   METraw_sumEt = 0;
   MET_corrPx = 0;
   MET_corrPy = 0;
   MET_et = 0;
   MET_phi = 0;
   MET_puppi_et = 0;
   MET_puppi_phi = 0;
   MET_sumEt = 0;
   MET_JetEnUp = 0;
   MET_JetEnDown = 0;
   MET_JetResUp = 0;
   MET_JetResDown = 0;
   MET_UnclusteredEnUp = 0;
   MET_UnclusteredEnDown = 0;
   PV_chi2 = 0;
   PV_ndof = 0;
   PV_rho = 0;
   PV_z = 0;
   BeamSpot_x0 = 0;
   BeamSpot_y0 = 0;
   BeamSpot_z0 = 0;
   IsBsTauTau = 0;
   BsTauTau_nCandidates = 0;
   BsTauTau_mu1_pt = 0;
   BsTauTau_mu1_eta = 0;
   BsTauTau_mu1_phi = 0;
   BsTauTau_mu1_mass = 0;
   BsTauTau_mu1_q = 0;
   BsTauTau_mu1_isLoose = 0;
   BsTauTau_mu1_isTight = 0;
   BsTauTau_mu1_isPF = 0;
   BsTauTau_mu1_isGlobal = 0;
   BsTauTau_mu1_isTracker = 0;
   BsTauTau_mu1_isSoft = 0;
   BsTauTau_mu1_vx = 0;
   BsTauTau_mu1_vy = 0;
   BsTauTau_mu1_vz = 0;
   BsTauTau_mu1_iso = 0;
   BsTauTau_mu1_dbiso = 0;
   BsTauTau_tau_pt = 0;
   BsTauTau_tau_eta = 0;
   BsTauTau_tau_phi = 0;
   BsTauTau_tau_mass = 0;
   BsTauTau_tau_rhomass1 = 0;
   BsTauTau_tau_rhomass2 = 0;
   BsTauTau_tau_q = 0;
   BsTauTau_tau_vx = 0;
   BsTauTau_tau_vy = 0;
   BsTauTau_tau_vz = 0;
   BsTauTau_tau_max_dr_3prong = 0;
   BsTauTau_tau_lip = 0;
   BsTauTau_tau_lips = 0;
   BsTauTau_tau_pvip = 0;
   BsTauTau_tau_pvips = 0;
   BsTauTau_tau_fl3d = 0;
   BsTauTau_tau_fls3d = 0;
   BsTauTau_tau_alpha = 0;
   BsTauTau_tau_vprob = 0;
   BsTauTau_tau_isRight = 0;
   BsTauTau_tau_isRight1 = 0;
   BsTauTau_tau_isRight2 = 0;
   BsTauTau_tau_isRight3 = 0;
   BsTauTau_tau_dr1 = 0;
   BsTauTau_tau_dr2 = 0;
   BsTauTau_tau_dr3 = 0;
   BsTauTau_tau_ptres1 = 0;
   BsTauTau_tau_ptres2 = 0;
   BsTauTau_tau_ptres3 = 0;
   BsTauTau_tau_matched_ppdgId = 0;
   BsTauTau_tau_matched_gentaupt = 0;
   BsTauTau_tau_gentaupt = 0;
   BsTauTau_tau_sumofdnn = 0;
   BsTauTau_tau_pfidx1 = 0;
   BsTauTau_tau_pfidx2 = 0;
   BsTauTau_tau_pfidx3 = 0;
   BsTauTau_tau_pi1_dnn = 0;
   BsTauTau_tau_pi2_dnn = 0;
   BsTauTau_tau_pi3_dnn = 0;
   BsTauTau_tau_pi1_pt = 0;
   BsTauTau_tau_pi1_eta = 0;
   BsTauTau_tau_pi1_phi = 0;
   BsTauTau_tau_pi1_mass = 0;
   BsTauTau_tau_pi1_charge = 0;
   BsTauTau_tau_pi1_dz = 0;
   BsTauTau_tau_muon_dr1 = 0;
   HLT_BPH_isFired = 0;
   BsTauTau_tau_pi2_pt = 0;
   BsTauTau_tau_pi2_eta = 0;
   BsTauTau_tau_pi2_phi = 0;
   BsTauTau_tau_pi2_mass = 0;
   BsTauTau_tau_pi2_charge = 0;
   BsTauTau_tau_pi2_dz = 0;
   BsTauTau_tau_muon_dr2 = 0;
   BsTauTau_tau_pi3_pt = 0;
   BsTauTau_tau_pi3_eta = 0;
   BsTauTau_tau_pi3_phi = 0;
   BsTauTau_tau_pi3_mass = 0;
   BsTauTau_tau_pi3_charge = 0;
   BsTauTau_tau_pi3_dz = 0;
   BsTauTau_tau_muon_dr3 = 0;
   BsTauTau_PV_vx = 0;
   BsTauTau_PV_vy = 0;
   BsTauTau_PV_vz = 0;
   BsTauTau_bbPV_vx = 0;
   BsTauTau_bbPV_vy = 0;
   BsTauTau_bbPV_vz = 0;
   BsTauTau_bbPV_refit_vx = 0;
   BsTauTau_bbPV_refit_vy = 0;
   BsTauTau_bbPV_refit_vz = 0;
   BsTauTau_genPV_vx = 0;
   BsTauTau_genPV_vy = 0;
   BsTauTau_genPV_vz = 0;
   BsTauTau_B_pt = 0;
   BsTauTau_B_eta = 0;
   BsTauTau_B_phi = 0;
   BsTauTau_B_mass = 0;
   BsTauTau_B_vprob = 0;
   BsTauTau_B_lip = 0;
   BsTauTau_B_lips = 0;
   BsTauTau_B_pvip = 0;
   BsTauTau_B_pvips = 0;
   BsTauTau_B_fl3d = 0;
   BsTauTau_B_fls3d = 0;
   BsTauTau_B_alpha = 0;
   BsTauTau_B_maxdoca = 0;
   BsTauTau_B_mindoca = 0;
   BsTauTau_muonpion_maxdoca = 0;
   BsTauTau_muonpion_mindoca = 0;
   BsTauTau_B_vx = 0;
   BsTauTau_B_vy = 0;
   BsTauTau_B_vz = 0;
   BsTauTau_B_iso = 0;
   BsTauTau_B_iso_ntracks = 0;
   BsTauTau_B_iso_mindoca = 0;
   BsTauTau_ngenmuons = 0;
   BsTauTau_isgen3 = 0;
   BsTauTau_isgen3matched = 0;
   BsTauTau_nch = 0;
   BsTauTau_nch_after_dnn = 0;
   BsTauTau_nch_before_dnn = 0;
   BsTauTau_nch_qr = 0;
   BsTauTau_ngentau3 = 0;
   BsTauTau_ngentau = 0;
   BsTauTau_gentaupt = 0;
   BsTauTau_gentaudm = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("METraw_et", &METraw_et, &b_METraw_et);
   fChain->SetBranchAddress("METraw_phi", &METraw_phi, &b_METraw_phi);
   fChain->SetBranchAddress("METraw_sumEt", &METraw_sumEt, &b_METraw_sumEt);
   fChain->SetBranchAddress("MET_corrPx", &MET_corrPx, &b_MET_corrPx);
   fChain->SetBranchAddress("MET_corrPy", &MET_corrPy, &b_MET_corrPy);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_puppi_et", &MET_puppi_et, &b_MET_puppi_et);
   fChain->SetBranchAddress("MET_puppi_phi", &MET_puppi_phi, &b_MET_puppi_phi);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("MET_JetEnUp", &MET_JetEnUp, &b_MET_JetEnUp);
   fChain->SetBranchAddress("MET_JetEnDown", &MET_JetEnDown, &b_MET_JetEnDown);
   fChain->SetBranchAddress("MET_JetResUp", &MET_JetResUp, &b_MET_JetResUp);
   fChain->SetBranchAddress("MET_JetResDown", &MET_JetResDown, &b_MET_JetResDown);
   fChain->SetBranchAddress("MET_UnclusteredEnUp", &MET_UnclusteredEnUp, &b_MET_UnclusteredEnUp);
   fChain->SetBranchAddress("MET_UnclusteredEnDown", &MET_UnclusteredEnDown, &b_MET_UnclusteredEnDown);
   fChain->SetBranchAddress("EVENT_event", &EVENT_event, &b_EVENT_event);
   fChain->SetBranchAddress("EVENT_run", &EVENT_run, &b_EVENT_run);
   fChain->SetBranchAddress("EVENT_lumiBlock", &EVENT_lumiBlock, &b_EVENT_lumiBlock);
   fChain->SetBranchAddress("PV_N", &PV_N, &b_PV_N);
   fChain->SetBranchAddress("PV_filter", &PV_filter, &b_PV_filter);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_rho", &PV_rho, &b_PV_rho);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("BeamSpot_x0", &BeamSpot_x0, &b_BeamSpot_x0);
   fChain->SetBranchAddress("BeamSpot_y0", &BeamSpot_y0, &b_BeamSpot_y0);
   fChain->SetBranchAddress("BeamSpot_z0", &BeamSpot_z0, &b_BeamSpot_z0);
   fChain->SetBranchAddress("IsBsTauTau", &IsBsTauTau, &b_IsBsTauTau);
   fChain->SetBranchAddress("BsTauTau_nCandidates", &BsTauTau_nCandidates, &b_BsTauTau_nCandidates);
   fChain->SetBranchAddress("BsTauTau_mu1_pt", &BsTauTau_mu1_pt, &b_BsTauTau_mu1_pt);
   fChain->SetBranchAddress("BsTauTau_mu1_eta", &BsTauTau_mu1_eta, &b_BsTauTau_mu1_eta);
   fChain->SetBranchAddress("BsTauTau_mu1_phi", &BsTauTau_mu1_phi, &b_BsTauTau_mu1_phi);
   fChain->SetBranchAddress("BsTauTau_mu1_mass", &BsTauTau_mu1_mass, &b_BsTauTau_mu1_mass);
   fChain->SetBranchAddress("BsTauTau_mu1_q", &BsTauTau_mu1_q, &b_BsTauTau_mu1_q);
   fChain->SetBranchAddress("BsTauTau_mu1_isLoose", &BsTauTau_mu1_isLoose, &b_BsTauTau_mu1_isLoose);
   fChain->SetBranchAddress("BsTauTau_mu1_isTight", &BsTauTau_mu1_isTight, &b_BsTauTau_mu1_isTight);
   fChain->SetBranchAddress("BsTauTau_mu1_isPF", &BsTauTau_mu1_isPF, &b_BsTauTau_mu1_isPF);
   fChain->SetBranchAddress("BsTauTau_mu1_isGlobal", &BsTauTau_mu1_isGlobal, &b_BsTauTau_mu1_isGlobal);
   fChain->SetBranchAddress("BsTauTau_mu1_isTracker", &BsTauTau_mu1_isTracker, &b_BsTauTau_mu1_isTracker);
   fChain->SetBranchAddress("BsTauTau_mu1_isSoft", &BsTauTau_mu1_isSoft, &b_BsTauTau_mu1_isSoft);
   fChain->SetBranchAddress("BsTauTau_mu1_vx", &BsTauTau_mu1_vx, &b_BsTauTau_mu1_vx);
   fChain->SetBranchAddress("BsTauTau_mu1_vy", &BsTauTau_mu1_vy, &b_BsTauTau_mu1_vy);
   fChain->SetBranchAddress("BsTauTau_mu1_vz", &BsTauTau_mu1_vz, &b_BsTauTau_mu1_vz);
   fChain->SetBranchAddress("BsTauTau_mu1_iso", &BsTauTau_mu1_iso, &b_BsTauTau_mu1_iso);
   fChain->SetBranchAddress("BsTauTau_mu1_dbiso", &BsTauTau_mu1_dbiso, &b_BsTauTau_mu1_dbiso);
   fChain->SetBranchAddress("BsTauTau_tau_pt", &BsTauTau_tau_pt, &b_BsTauTau_tau_pt);
   fChain->SetBranchAddress("BsTauTau_tau_eta", &BsTauTau_tau_eta, &b_BsTauTau_tau_eta);
   fChain->SetBranchAddress("BsTauTau_tau_phi", &BsTauTau_tau_phi, &b_BsTauTau_tau_phi);
   fChain->SetBranchAddress("BsTauTau_tau_mass", &BsTauTau_tau_mass, &b_BsTauTau_tau_mass);
   fChain->SetBranchAddress("BsTauTau_tau_rhomass1", &BsTauTau_tau_rhomass1, &b_BsTauTau_tau_rhomass1);
   fChain->SetBranchAddress("BsTauTau_tau_rhomass2", &BsTauTau_tau_rhomass2, &b_BsTauTau_tau_rhomass2);
   fChain->SetBranchAddress("BsTauTau_tau_q", &BsTauTau_tau_q, &b_BsTauTau_tau_q);
   fChain->SetBranchAddress("BsTauTau_tau_vx", &BsTauTau_tau_vx, &b_BsTauTau_tau_vx);
   fChain->SetBranchAddress("BsTauTau_tau_vy", &BsTauTau_tau_vy, &b_BsTauTau_tau_vy);
   fChain->SetBranchAddress("BsTauTau_tau_vz", &BsTauTau_tau_vz, &b_BsTauTau_tau_vz);
   fChain->SetBranchAddress("BsTauTau_tau_max_dr_3prong", &BsTauTau_tau_max_dr_3prong, &b_BsTauTau_tau_max_dr_3prong);
   fChain->SetBranchAddress("BsTauTau_tau_lip", &BsTauTau_tau_lip, &b_BsTauTau_tau_lip);
   fChain->SetBranchAddress("BsTauTau_tau_lips", &BsTauTau_tau_lips, &b_BsTauTau_tau_lips);
   fChain->SetBranchAddress("BsTauTau_tau_pvip", &BsTauTau_tau_pvip, &b_BsTauTau_tau_pvip);
   fChain->SetBranchAddress("BsTauTau_tau_pvips", &BsTauTau_tau_pvips, &b_BsTauTau_tau_pvips);
   fChain->SetBranchAddress("BsTauTau_tau_fl3d", &BsTauTau_tau_fl3d, &b_BsTauTau_tau_fl3d);
   fChain->SetBranchAddress("BsTauTau_tau_fls3d", &BsTauTau_tau_fls3d, &b_BsTauTau_tau_fls3d);
   fChain->SetBranchAddress("BsTauTau_tau_alpha", &BsTauTau_tau_alpha, &b_BsTauTau_tau_alpha);
   fChain->SetBranchAddress("BsTauTau_tau_vprob", &BsTauTau_tau_vprob, &b_BsTauTau_tau_vprob);
   fChain->SetBranchAddress("BsTauTau_tau_isRight", &BsTauTau_tau_isRight, &b_BsTauTau_tau_isRight);
   fChain->SetBranchAddress("BsTauTau_tau_isRight1", &BsTauTau_tau_isRight1, &b_BsTauTau_tau_isRight1);
   fChain->SetBranchAddress("BsTauTau_tau_isRight2", &BsTauTau_tau_isRight2, &b_BsTauTau_tau_isRight2);
   fChain->SetBranchAddress("BsTauTau_tau_isRight3", &BsTauTau_tau_isRight3, &b_BsTauTau_tau_isRight3);
   fChain->SetBranchAddress("BsTauTau_tau_dr1", &BsTauTau_tau_dr1, &b_BsTauTau_tau_dr1);
   fChain->SetBranchAddress("BsTauTau_tau_dr2", &BsTauTau_tau_dr2, &b_BsTauTau_tau_dr2);
   fChain->SetBranchAddress("BsTauTau_tau_dr3", &BsTauTau_tau_dr3, &b_BsTauTau_tau_dr3);
   fChain->SetBranchAddress("BsTauTau_tau_ptres1", &BsTauTau_tau_ptres1, &b_BsTauTau_tau_ptres1);
   fChain->SetBranchAddress("BsTauTau_tau_ptres2", &BsTauTau_tau_ptres2, &b_BsTauTau_tau_ptres2);
   fChain->SetBranchAddress("BsTauTau_tau_ptres3", &BsTauTau_tau_ptres3, &b_BsTauTau_tau_ptres3);
   fChain->SetBranchAddress("BsTauTau_tau_matched_ppdgId", &BsTauTau_tau_matched_ppdgId, &b_BsTauTau_tau_matched_ppdgId);
   fChain->SetBranchAddress("BsTauTau_tau_matched_gentaupt", &BsTauTau_tau_matched_gentaupt, &b_BsTauTau_tau_matched_gentaupt);
   fChain->SetBranchAddress("BsTauTau_tau_gentaupt", &BsTauTau_tau_gentaupt, &b_BsTauTau_tau_gentaupt);
   fChain->SetBranchAddress("BsTauTau_tau_sumofdnn", &BsTauTau_tau_sumofdnn, &b_BsTauTau_tau_sumofdnn);
   fChain->SetBranchAddress("BsTauTau_tau_pfidx1", &BsTauTau_tau_pfidx1, &b_BsTauTau_tau_pfidx1);
   fChain->SetBranchAddress("BsTauTau_tau_pfidx2", &BsTauTau_tau_pfidx2, &b_BsTauTau_tau_pfidx2);
   fChain->SetBranchAddress("BsTauTau_tau_pfidx3", &BsTauTau_tau_pfidx3, &b_BsTauTau_tau_pfidx3);
   fChain->SetBranchAddress("BsTauTau_tau_pi1_dnn", &BsTauTau_tau_pi1_dnn, &b_BsTauTau_tau_pi1_dnn);
   fChain->SetBranchAddress("BsTauTau_tau_pi2_dnn", &BsTauTau_tau_pi2_dnn, &b_BsTauTau_tau_pi2_dnn);
   fChain->SetBranchAddress("BsTauTau_tau_pi3_dnn", &BsTauTau_tau_pi3_dnn, &b_BsTauTau_tau_pi3_dnn);
   fChain->SetBranchAddress("BsTauTau_tau_pi1_pt", &BsTauTau_tau_pi1_pt, &b_BsTauTau_tau_pi1_pt);
   fChain->SetBranchAddress("BsTauTau_tau_pi1_eta", &BsTauTau_tau_pi1_eta, &b_BsTauTau_tau_pi1_eta);
   fChain->SetBranchAddress("BsTauTau_tau_pi1_phi", &BsTauTau_tau_pi1_phi, &b_BsTauTau_tau_pi1_phi);
   fChain->SetBranchAddress("BsTauTau_tau_pi1_mass", &BsTauTau_tau_pi1_mass, &b_BsTauTau_tau_pi1_mass);
   fChain->SetBranchAddress("BsTauTau_tau_pi1_charge", &BsTauTau_tau_pi1_charge, &b_BsTauTau_tau_pi1_charge);
   fChain->SetBranchAddress("BsTauTau_tau_pi1_dz", &BsTauTau_tau_pi1_dz, &b_BsTauTau_tau_pi1_dz);
   fChain->SetBranchAddress("BsTauTau_tau_muon_dr1", &BsTauTau_tau_muon_dr1, &b_BsTauTau_tau_muon_dr1);
   fChain->SetBranchAddress("HLT_BPH_isFired", &HLT_BPH_isFired, &b_HLT_BPH_isFired);
   fChain->SetBranchAddress("BsTauTau_tau_pi2_pt", &BsTauTau_tau_pi2_pt, &b_BsTauTau_tau_pi2_pt);
   fChain->SetBranchAddress("BsTauTau_tau_pi2_eta", &BsTauTau_tau_pi2_eta, &b_BsTauTau_tau_pi2_eta);
   fChain->SetBranchAddress("BsTauTau_tau_pi2_phi", &BsTauTau_tau_pi2_phi, &b_BsTauTau_tau_pi2_phi);
   fChain->SetBranchAddress("BsTauTau_tau_pi2_mass", &BsTauTau_tau_pi2_mass, &b_BsTauTau_tau_pi2_mass);
   fChain->SetBranchAddress("BsTauTau_tau_pi2_charge", &BsTauTau_tau_pi2_charge, &b_BsTauTau_tau_pi2_charge);
   fChain->SetBranchAddress("BsTauTau_tau_pi2_dz", &BsTauTau_tau_pi2_dz, &b_BsTauTau_tau_pi2_dz);
   fChain->SetBranchAddress("BsTauTau_tau_muon_dr2", &BsTauTau_tau_muon_dr2, &b_BsTauTau_tau_muon_dr2);
   fChain->SetBranchAddress("BsTauTau_tau_pi3_pt", &BsTauTau_tau_pi3_pt, &b_BsTauTau_tau_pi3_pt);
   fChain->SetBranchAddress("BsTauTau_tau_pi3_eta", &BsTauTau_tau_pi3_eta, &b_BsTauTau_tau_pi3_eta);
   fChain->SetBranchAddress("BsTauTau_tau_pi3_phi", &BsTauTau_tau_pi3_phi, &b_BsTauTau_tau_pi3_phi);
   fChain->SetBranchAddress("BsTauTau_tau_pi3_mass", &BsTauTau_tau_pi3_mass, &b_BsTauTau_tau_pi3_mass);
   fChain->SetBranchAddress("BsTauTau_tau_pi3_charge", &BsTauTau_tau_pi3_charge, &b_BsTauTau_tau_pi3_charge);
   fChain->SetBranchAddress("BsTauTau_tau_pi3_dz", &BsTauTau_tau_pi3_dz, &b_BsTauTau_tau_pi3_dz);
   fChain->SetBranchAddress("BsTauTau_tau_muon_dr3", &BsTauTau_tau_muon_dr3, &b_BsTauTau_tau_muon_dr3);
   fChain->SetBranchAddress("BsTauTau_PV_vx", &BsTauTau_PV_vx, &b_BsTauTau_PV_vx);
   fChain->SetBranchAddress("BsTauTau_PV_vy", &BsTauTau_PV_vy, &b_BsTauTau_PV_vy);
   fChain->SetBranchAddress("BsTauTau_PV_vz", &BsTauTau_PV_vz, &b_BsTauTau_PV_vz);
   fChain->SetBranchAddress("BsTauTau_bbPV_vx", &BsTauTau_bbPV_vx, &b_BsTauTau_bbPV_vx);
   fChain->SetBranchAddress("BsTauTau_bbPV_vy", &BsTauTau_bbPV_vy, &b_BsTauTau_bbPV_vy);
   fChain->SetBranchAddress("BsTauTau_bbPV_vz", &BsTauTau_bbPV_vz, &b_BsTauTau_bbPV_vz);
   fChain->SetBranchAddress("BsTauTau_bbPV_refit_vx", &BsTauTau_bbPV_refit_vx, &b_BsTauTau_bbPV_refit_vx);
   fChain->SetBranchAddress("BsTauTau_bbPV_refit_vy", &BsTauTau_bbPV_refit_vy, &b_BsTauTau_bbPV_refit_vy);
   fChain->SetBranchAddress("BsTauTau_bbPV_refit_vz", &BsTauTau_bbPV_refit_vz, &b_BsTauTau_bbPV_refit_vz);
   fChain->SetBranchAddress("BsTauTau_genPV_vx", &BsTauTau_genPV_vx, &b_BsTauTau_genPV_vx);
   fChain->SetBranchAddress("BsTauTau_genPV_vy", &BsTauTau_genPV_vy, &b_BsTauTau_genPV_vy);
   fChain->SetBranchAddress("BsTauTau_genPV_vz", &BsTauTau_genPV_vz, &b_BsTauTau_genPV_vz);
   fChain->SetBranchAddress("BsTauTau_B_pt", &BsTauTau_B_pt, &b_BsTauTau_B_pt);
   fChain->SetBranchAddress("BsTauTau_B_eta", &BsTauTau_B_eta, &b_BsTauTau_B_eta);
   fChain->SetBranchAddress("BsTauTau_B_phi", &BsTauTau_B_phi, &b_BsTauTau_B_phi);
   fChain->SetBranchAddress("BsTauTau_B_mass", &BsTauTau_B_mass, &b_BsTauTau_B_mass);
   fChain->SetBranchAddress("BsTauTau_B_vprob", &BsTauTau_B_vprob, &b_BsTauTau_B_vprob);
   fChain->SetBranchAddress("BsTauTau_B_lip", &BsTauTau_B_lip, &b_BsTauTau_B_lip);
   fChain->SetBranchAddress("BsTauTau_B_lips", &BsTauTau_B_lips, &b_BsTauTau_B_lips);
   fChain->SetBranchAddress("BsTauTau_B_pvip", &BsTauTau_B_pvip, &b_BsTauTau_B_pvip);
   fChain->SetBranchAddress("BsTauTau_B_pvips", &BsTauTau_B_pvips, &b_BsTauTau_B_pvips);
   fChain->SetBranchAddress("BsTauTau_B_fl3d", &BsTauTau_B_fl3d, &b_BsTauTau_B_fl3d);
   fChain->SetBranchAddress("BsTauTau_B_fls3d", &BsTauTau_B_fls3d, &b_BsTauTau_B_fls3d);
   fChain->SetBranchAddress("BsTauTau_B_alpha", &BsTauTau_B_alpha, &b_BsTauTau_B_alpha);
   fChain->SetBranchAddress("BsTauTau_B_maxdoca", &BsTauTau_B_maxdoca, &b_BsTauTau_B_maxdoca);
   fChain->SetBranchAddress("BsTauTau_B_mindoca", &BsTauTau_B_mindoca, &b_BsTauTau_B_mindoca);
   fChain->SetBranchAddress("BsTauTau_muonpion_maxdoca", &BsTauTau_muonpion_maxdoca, &b_BsTauTau_muonpion_maxdoca);
   fChain->SetBranchAddress("BsTauTau_muonpion_mindoca", &BsTauTau_muonpion_mindoca, &b_BsTauTau_muonpion_mindoca);
   fChain->SetBranchAddress("BsTauTau_B_vx", &BsTauTau_B_vx, &b_BsTauTau_B_vx);
   fChain->SetBranchAddress("BsTauTau_B_vy", &BsTauTau_B_vy, &b_BsTauTau_B_vy);
   fChain->SetBranchAddress("BsTauTau_B_vz", &BsTauTau_B_vz, &b_BsTauTau_B_vz);
   fChain->SetBranchAddress("BsTauTau_B_iso", &BsTauTau_B_iso, &b_BsTauTau_B_iso);
   fChain->SetBranchAddress("BsTauTau_B_iso_ntracks", &BsTauTau_B_iso_ntracks, &b_BsTauTau_B_iso_ntracks);
   fChain->SetBranchAddress("BsTauTau_B_iso_mindoca", &BsTauTau_B_iso_mindoca, &b_BsTauTau_B_iso_mindoca);
   fChain->SetBranchAddress("BsTauTau_ngenmuons", &BsTauTau_ngenmuons, &b_BsTauTau_ngenmuons);
   fChain->SetBranchAddress("BsTauTau_isgen3", &BsTauTau_isgen3, &b_BsTauTau_isgen3);
   fChain->SetBranchAddress("BsTauTau_isgen3matched", &BsTauTau_isgen3matched, &b_BsTauTau_isgen3matched);
   fChain->SetBranchAddress("BsTauTau_nch", &BsTauTau_nch, &b_BsTauTau_nch);
   fChain->SetBranchAddress("BsTauTau_nch_after_dnn", &BsTauTau_nch_after_dnn, &b_BsTauTau_nch_after_dnn);
   fChain->SetBranchAddress("BsTauTau_nch_before_dnn", &BsTauTau_nch_before_dnn, &b_BsTauTau_nch_before_dnn);
   fChain->SetBranchAddress("BsTauTau_nch_qr", &BsTauTau_nch_qr, &b_BsTauTau_nch_qr);
   fChain->SetBranchAddress("BsTauTau_ngentau3", &BsTauTau_ngentau3, &b_BsTauTau_ngentau3);
   fChain->SetBranchAddress("BsTauTau_ngentau", &BsTauTau_ngentau, &b_BsTauTau_ngentau);
   fChain->SetBranchAddress("BsTauTau_gentaupt", &BsTauTau_gentaupt, &b_BsTauTau_gentaupt);
   fChain->SetBranchAddress("BsTauTau_gentaudm", &BsTauTau_gentaudm, &b_BsTauTau_gentaudm);
   Notify();
}

Bool_t bpark::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bpark::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bpark::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef bpark_cxx
