//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  1 12:33:58 2022 by ROOT version 6.24/07
// from TTree tree/tree
// found on file: RelVal_DisplacedSUSY_NO_PU_Lxy_added.root
//////////////////////////////////////////////////////////

#ifndef RelVal_eff_h
#define RelVal_eff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class RelVal_eff {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *pt_glbMuon;
   vector<double>  *eta_glbMuon;
   vector<double>  *phi_glbMuon;
   vector<double>  *dxy_glbMuon;
   vector<double>  *dz_glbMuon;
   vector<double>  *d0_glbMuon;
   vector<double>  *pt_DisglbMuon;
   vector<double>  *eta_DisglbMuon;
   vector<double>  *phi_DisglbMuon;
   vector<double>  *dxy_DisglbMuon;
   vector<double>  *dz_DisglbMuon;
   vector<double>  *d0_DisglbMuon;
   vector<double>  *lxy_DisglbMuon;
   vector<double>  *pt_Sta;
   vector<double>  *eta_Sta;
   vector<double>  *phi_Sta;
   vector<double>  *dxy_Sta;
   vector<double>  *dz_Sta;
   vector<double>  *d0_Sta;
   vector<double>  *pt_DisSta;
   vector<double>  *eta_DisSta;
   vector<double>  *phi_DisSta;
   vector<double>  *dxy_DisSta;
   vector<double>  *dz_DisSta;
   vector<double>  *d0_DisSta;
   vector<double>  *lxy_DisSta;
   vector<double>  *pt_DisTrk;
   vector<double>  *eta_DisTrk;
   vector<double>  *phi_DisTrk;
   vector<double>  *dxy_DisTrk;
   vector<double>  *dz_DisTrk;
   vector<double>  *d0_DisTrk;
   vector<double>  *lxy_DisTrk;
   vector<double>  *pt_DisOutIn;
   vector<double>  *eta_DisOutIn;
   vector<double>  *phi_DisOutIn;
   vector<double>  *dxy_DisOutIn;
   vector<double>  *dz_DisOutIn;
   vector<double>  *d0_DisOutIn;
   vector<double>  *pt_DisEar;
   vector<double>  *eta_DisEar;
   vector<double>  *phi_DisEar;
   vector<double>  *dxy_DisEar;
   vector<double>  *dz_DisEar;
   vector<double>  *d0_DisEar;
   vector<double>  *pt_gen;
   vector<double>  *eta_gen;
   vector<double>  *phi_gen;
   vector<double>  *vx_gen;
   vector<double>  *vy_gen;
   vector<double>  *vz_gen;
   vector<double>  *vertex_genx;
   vector<double>  *vertex_geny;
   vector<double>  *dxy_gen;
   vector<double>  *lxy_gen;
   vector<double>  *pt_Muon;
   vector<double>  *eta_Muon;
   vector<double>  *phi_Muon;
   vector<double>  *mass_Muon;
   vector<double>  *dxy_Muon;
   vector<double>  *dz_Muon;
   vector<double>  *d0_Muon;
   vector<double>  *pt_displacedMuon;
   vector<double>  *eta_displacedMuon;
   vector<double>  *phi_displacedMuon;
   vector<double>  *mass_displacedMuon;
   vector<double>  *dxy_displacedMuon;
   vector<double>  *dz_displacedMuon;
   vector<double>  *d0_displacedMuon;
   Int_t           numberOfOISeeds;

   // List of branches
   TBranch        *b_pt_glbMuon;   //!
   TBranch        *b_eta_glbMuon;   //!
   TBranch        *b_phi_glbMuon;   //!
   TBranch        *b_dxy_glbMuon;   //!
   TBranch        *b_dz_glbMuon;   //!
   TBranch        *b_d0_glbMuon;   //!
   TBranch        *b_pt_DisglbMuon;   //!
   TBranch        *b_eta_DisglbMuon;   //!
   TBranch        *b_phi_DisglbMuon;   //!
   TBranch        *b_dxy_DisglbMuon;   //!
   TBranch        *b_dz_DisglbMuon;   //!
   TBranch        *b_d0_DisglbMuon;   //!
   TBranch        *b_lxy_DisglbMuon;   //!
   TBranch        *b_pt_Sta;   //!
   TBranch        *b_eta_Sta;   //!
   TBranch        *b_phi_Sta;   //!
   TBranch        *b_dxy_Sta;   //!
   TBranch        *b_dz_Sta;   //!
   TBranch        *b_d0_Sta;   //!
   TBranch        *b_pt_DisSta;   //!
   TBranch        *b_eta_DisSta;   //!
   TBranch        *b_phi_DisSta;   //!
   TBranch        *b_dxy_DisSta;   //!
   TBranch        *b_dz_DisSta;   //!
   TBranch        *b_d0_DisSta;   //!
   TBranch        *b_lxy_DisSta;   //!
   TBranch        *b_pt_DisTrk;   //!
   TBranch        *b_eta_DisTrk;   //!
   TBranch        *b_phi_DisTrk;   //!
   TBranch        *b_dxy_DisTrk;   //!
   TBranch        *b_dz_DisTrk;   //!
   TBranch        *b_d0_DisTrk;   //!
   TBranch        *b_lxy_DisTrk;   //!
   TBranch        *b_pt_DisOutIn;   //!
   TBranch        *b_eta_DisOutIn;   //!
   TBranch        *b_phi_DisOutIn;   //!
   TBranch        *b_dxy_DisOutIn;   //!
   TBranch        *b_dz_DisOutIn;   //!
   TBranch        *b_d0_DisOutIn;   //!
   TBranch        *b_pt_DisEar;   //!
   TBranch        *b_eta_DisEar;   //!
   TBranch        *b_phi_DisEar;   //!
   TBranch        *b_dxy_DisEar;   //!
   TBranch        *b_dz_DisEar;   //!
   TBranch        *b_d0_DisEar;   //!
   TBranch        *b_pt_gen;   //!
   TBranch        *b_eta_gen;   //!
   TBranch        *b_phi_gen;   //!
   TBranch        *b_vx_gen;   //!
   TBranch        *b_vy_gen;   //!
   TBranch        *b_vz_gen;   //!
   TBranch        *b_vertex_genx;   //!
   TBranch        *b_vertex_geny;   //!
   TBranch        *b_dxy_gen;   //!
   TBranch        *b_lxy_gen;   //!
   TBranch        *b_pt_Muon;   //!
   TBranch        *b_eta_Muon;   //!
   TBranch        *b_phi_Muon;   //!
   TBranch        *b_mass_Muon;   //!
   TBranch        *b_dxy_Muon;   //!
   TBranch        *b_dz_Muon;   //!
   TBranch        *b_d0_Muon;   //!
   TBranch        *b_pt_displacedMuon;   //!
   TBranch        *b_eta_displacedMuon;   //!
   TBranch        *b_phi_displacedMuon;   //!
   TBranch        *b_mass_displacedMuon;   //!
   TBranch        *b_dxy_displacedMuon;   //!
   TBranch        *b_dz_displacedMuon;   //!
   TBranch        *b_d0_displacedMuon;   //!
   TBranch        *b_numberOfOISeeds;   //!

   RelVal_eff(TTree *tree=0);
   virtual ~RelVal_eff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RelVal_eff_cxx
RelVal_eff::RelVal_eff(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RelVal_DisplacedSUSY_NO_PU_Lxy_added.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RelVal_DisplacedSUSY_NO_PU_Lxy_added.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("RelVal_DisplacedSUSY_NO_PU_Lxy_added.root:/ntuple");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

RelVal_eff::~RelVal_eff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RelVal_eff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RelVal_eff::LoadTree(Long64_t entry)
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

void RelVal_eff::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pt_glbMuon = 0;
   eta_glbMuon = 0;
   phi_glbMuon = 0;
   dxy_glbMuon = 0;
   dz_glbMuon = 0;
   d0_glbMuon = 0;
   pt_DisglbMuon = 0;
   eta_DisglbMuon = 0;
   phi_DisglbMuon = 0;
   dxy_DisglbMuon = 0;
   dz_DisglbMuon = 0;
   d0_DisglbMuon = 0;
   lxy_DisglbMuon = 0;
   pt_Sta = 0;
   eta_Sta = 0;
   phi_Sta = 0;
   dxy_Sta = 0;
   dz_Sta = 0;
   d0_Sta = 0;
   pt_DisSta = 0;
   eta_DisSta = 0;
   phi_DisSta = 0;
   dxy_DisSta = 0;
   dz_DisSta = 0;
   d0_DisSta = 0;
   lxy_DisSta = 0;
   pt_DisTrk = 0;
   eta_DisTrk = 0;
   phi_DisTrk = 0;
   dxy_DisTrk = 0;
   dz_DisTrk = 0;
   d0_DisTrk = 0;
   lxy_DisTrk = 0;
   pt_DisOutIn = 0;
   eta_DisOutIn = 0;
   phi_DisOutIn = 0;
   dxy_DisOutIn = 0;
   dz_DisOutIn = 0;
   d0_DisOutIn = 0;
   pt_DisEar = 0;
   eta_DisEar = 0;
   phi_DisEar = 0;
   dxy_DisEar = 0;
   dz_DisEar = 0;
   d0_DisEar = 0;
   pt_gen = 0;
   eta_gen = 0;
   phi_gen = 0;
   vx_gen = 0;
   vy_gen = 0;
   vz_gen = 0;
   vertex_genx = 0;
   vertex_geny = 0;
   dxy_gen = 0;
   lxy_gen = 0;
   pt_Muon = 0;
   eta_Muon = 0;
   phi_Muon = 0;
   mass_Muon = 0;
   dxy_Muon = 0;
   dz_Muon = 0;
   d0_Muon = 0;
   pt_displacedMuon = 0;
   eta_displacedMuon = 0;
   phi_displacedMuon = 0;
   mass_displacedMuon = 0;
   dxy_displacedMuon = 0;
   dz_displacedMuon = 0;
   d0_displacedMuon = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pt_glbMuon", &pt_glbMuon, &b_pt_glbMuon);
   fChain->SetBranchAddress("eta_glbMuon", &eta_glbMuon, &b_eta_glbMuon);
   fChain->SetBranchAddress("phi_glbMuon", &phi_glbMuon, &b_phi_glbMuon);
   fChain->SetBranchAddress("dxy_glbMuon", &dxy_glbMuon, &b_dxy_glbMuon);
   fChain->SetBranchAddress("dz_glbMuon", &dz_glbMuon, &b_dz_glbMuon);
   fChain->SetBranchAddress("d0_glbMuon", &d0_glbMuon, &b_d0_glbMuon);
   fChain->SetBranchAddress("pt_DisglbMuon", &pt_DisglbMuon, &b_pt_DisglbMuon);
   fChain->SetBranchAddress("eta_DisglbMuon", &eta_DisglbMuon, &b_eta_DisglbMuon);
   fChain->SetBranchAddress("phi_DisglbMuon", &phi_DisglbMuon, &b_phi_DisglbMuon);
   fChain->SetBranchAddress("dxy_DisglbMuon", &dxy_DisglbMuon, &b_dxy_DisglbMuon);
   fChain->SetBranchAddress("dz_DisglbMuon", &dz_DisglbMuon, &b_dz_DisglbMuon);
   fChain->SetBranchAddress("d0_DisglbMuon", &d0_DisglbMuon, &b_d0_DisglbMuon);
   fChain->SetBranchAddress("lxy_DisglbMuon", &lxy_DisglbMuon, &b_lxy_DisglbMuon);
   fChain->SetBranchAddress("pt_Sta", &pt_Sta, &b_pt_Sta);
   fChain->SetBranchAddress("eta_Sta", &eta_Sta, &b_eta_Sta);
   fChain->SetBranchAddress("phi_Sta", &phi_Sta, &b_phi_Sta);
   fChain->SetBranchAddress("dxy_Sta", &dxy_Sta, &b_dxy_Sta);
   fChain->SetBranchAddress("dz_Sta", &dz_Sta, &b_dz_Sta);
   fChain->SetBranchAddress("d0_Sta", &d0_Sta, &b_d0_Sta);
   fChain->SetBranchAddress("pt_DisSta", &pt_DisSta, &b_pt_DisSta);
   fChain->SetBranchAddress("eta_DisSta", &eta_DisSta, &b_eta_DisSta);
   fChain->SetBranchAddress("phi_DisSta", &phi_DisSta, &b_phi_DisSta);
   fChain->SetBranchAddress("dxy_DisSta", &dxy_DisSta, &b_dxy_DisSta);
   fChain->SetBranchAddress("dz_DisSta", &dz_DisSta, &b_dz_DisSta);
   fChain->SetBranchAddress("d0_DisSta", &d0_DisSta, &b_d0_DisSta);
   fChain->SetBranchAddress("lxy_DisSta", &lxy_DisSta, &b_lxy_DisSta);
   fChain->SetBranchAddress("pt_DisTrk", &pt_DisTrk, &b_pt_DisTrk);
   fChain->SetBranchAddress("eta_DisTrk", &eta_DisTrk, &b_eta_DisTrk);
   fChain->SetBranchAddress("phi_DisTrk", &phi_DisTrk, &b_phi_DisTrk);
   fChain->SetBranchAddress("dxy_DisTrk", &dxy_DisTrk, &b_dxy_DisTrk);
   fChain->SetBranchAddress("dz_DisTrk", &dz_DisTrk, &b_dz_DisTrk);
   fChain->SetBranchAddress("d0_DisTrk", &d0_DisTrk, &b_d0_DisTrk);
   fChain->SetBranchAddress("lxy_DisTrk", &lxy_DisTrk, &b_lxy_DisTrk);
   fChain->SetBranchAddress("pt_DisOutIn", &pt_DisOutIn, &b_pt_DisOutIn);
   fChain->SetBranchAddress("eta_DisOutIn", &eta_DisOutIn, &b_eta_DisOutIn);
   fChain->SetBranchAddress("phi_DisOutIn", &phi_DisOutIn, &b_phi_DisOutIn);
   fChain->SetBranchAddress("dxy_DisOutIn", &dxy_DisOutIn, &b_dxy_DisOutIn);
   fChain->SetBranchAddress("dz_DisOutIn", &dz_DisOutIn, &b_dz_DisOutIn);
   fChain->SetBranchAddress("d0_DisOutIn", &d0_DisOutIn, &b_d0_DisOutIn);
   fChain->SetBranchAddress("pt_DisEar", &pt_DisEar, &b_pt_DisEar);
   fChain->SetBranchAddress("eta_DisEar", &eta_DisEar, &b_eta_DisEar);
   fChain->SetBranchAddress("phi_DisEar", &phi_DisEar, &b_phi_DisEar);
   fChain->SetBranchAddress("dxy_DisEar", &dxy_DisEar, &b_dxy_DisEar);
   fChain->SetBranchAddress("dz_DisEar", &dz_DisEar, &b_dz_DisEar);
   fChain->SetBranchAddress("d0_DisEar", &d0_DisEar, &b_d0_DisEar);
   fChain->SetBranchAddress("pt_gen", &pt_gen, &b_pt_gen);
   fChain->SetBranchAddress("eta_gen", &eta_gen, &b_eta_gen);
   fChain->SetBranchAddress("phi_gen", &phi_gen, &b_phi_gen);
   fChain->SetBranchAddress("vx_gen", &vx_gen, &b_vx_gen);
   fChain->SetBranchAddress("vy_gen", &vy_gen, &b_vy_gen);
   fChain->SetBranchAddress("vz_gen", &vz_gen, &b_vz_gen);
   fChain->SetBranchAddress("vertex_genx", &vertex_genx, &b_vertex_genx);
   fChain->SetBranchAddress("vertex_geny", &vertex_geny, &b_vertex_geny);
   fChain->SetBranchAddress("dxy_gen", &dxy_gen, &b_dxy_gen);
   fChain->SetBranchAddress("lxy_gen", &lxy_gen, &b_lxy_gen);
   fChain->SetBranchAddress("pt_Muon", &pt_Muon, &b_pt_Muon);
   fChain->SetBranchAddress("eta_Muon", &eta_Muon, &b_eta_Muon);
   fChain->SetBranchAddress("phi_Muon", &phi_Muon, &b_phi_Muon);
   fChain->SetBranchAddress("mass_Muon", &mass_Muon, &b_mass_Muon);
   fChain->SetBranchAddress("dxy_Muon", &dxy_Muon, &b_dxy_Muon);
   fChain->SetBranchAddress("dz_Muon", &dz_Muon, &b_dz_Muon);
   fChain->SetBranchAddress("d0_Muon", &d0_Muon, &b_d0_Muon);
   fChain->SetBranchAddress("pt_displacedMuon", &pt_displacedMuon, &b_pt_displacedMuon);
   fChain->SetBranchAddress("eta_displacedMuon", &eta_displacedMuon, &b_eta_displacedMuon);
   fChain->SetBranchAddress("phi_displacedMuon", &phi_displacedMuon, &b_phi_displacedMuon);
   fChain->SetBranchAddress("mass_displacedMuon", &mass_displacedMuon, &b_mass_displacedMuon);
   fChain->SetBranchAddress("dxy_displacedMuon", &dxy_displacedMuon, &b_dxy_displacedMuon);
   fChain->SetBranchAddress("dz_displacedMuon", &dz_displacedMuon, &b_dz_displacedMuon);
   fChain->SetBranchAddress("d0_displacedMuon", &d0_displacedMuon, &b_d0_displacedMuon);
   fChain->SetBranchAddress("numberOfOISeeds", &numberOfOISeeds, &b_numberOfOISeeds);
   Notify();
}

Bool_t RelVal_eff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RelVal_eff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RelVal_eff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RelVal_eff_cxx
