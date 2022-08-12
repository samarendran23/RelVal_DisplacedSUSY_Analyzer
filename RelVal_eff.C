#define RelVal_eff_cxx
#include "RelVal_eff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void RelVal_eff::Loop()
{
//   In a ROOT session, you can do:
//      root> .L RelVal_eff.C
//      root> RelVal_eff t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch



   //TH1F *hNgen = new TH1F("hNgen","hNgen",5,0,5);
   //TH1F *hNdisglb = new TH1F("hNdisglb","hNdisglb",5,0,5);
   //TH1F *hNdisSta = new TH1F("hNdisSta","hNdisSta",5,0,5);


   TH1F *hPtgen = new TH1F("hPtgen","hPtgen",100,0,100);


   TH1F *hDisGlbpass_pt =new TH1F("hDisGlbpass_pt","hDisGlbpass_pt",10,0,1000);
   TH1F *hDisGlbtotal_pt =new TH1F("hDisGlbtotal_pt","hDisGlbtotal_pt",10,0,1000);

   TH1F *hDisStapass_pt =new TH1F("hDisStapass_pt","hDisStapass_pt",10,0,1000);
   TH1F *hDisStatotal_pt =new TH1F("hDisStatotal_pt","hDisStatotal_pt",10,0,1000);

   TH1F *hDisTrkpass_pt =new TH1F("hDisTrkpass_pt","hDisTrkpass_pt",10,0,1000);
   TH1F *hDisTrktotal_pt =new TH1F("hDisTrktotal_pt","hDisTrktotal_pt",10,0,1000);


   TH1F *hDisGlbpass_eta =new TH1F("hDisGlbpass_eta","hDisGlbpass_eta",10,-3,3);
   TH1F *hDisGlbtotal_eta =new TH1F("hDisGlbtotal_eta","hDisGlbtotal_eta",10,-3,3);

   TH1F *hDisStapass_eta =new TH1F("hDisStapass_eta","hDisStapass_eta",10,-3,3);
   TH1F *hDisStatotal_eta =new TH1F("hDisStatotal_eta","hDisStatotal_eta",10,-3,3);

   TH1F *hDisTrkpass_eta =new TH1F("hDisTrkpass_eta","hDisTrkpass_eta",10,-3,3);
   TH1F *hDisTrktotal_eta =new TH1F("hDisTrktotal_eta","hDisTrktotal_eta",10,-3,3);


   TH1F *hDisGlbpass_dxy =new TH1F("hDisGlbpass_dxy","hDisGlbpass_dxy",10,0.0,30.0);
   TH1F *hDisGlbtotal_dxy =new TH1F("hDisGlbtotal_dxy","hDisGlbtotal_dxy",10,0.0,30.0);

   TH1F *hDisStapass_dxy =new TH1F("hDisStapass_dxy","hDisStapass_dxy",10,0.0,30.0);
   TH1F *hDisStatotal_dxy =new TH1F("hDisStatotal_dxy","hDisStatotal_dxy",10,0.0,30.0);

   TH1F *hDisTrkpass_dxy =new TH1F("hDisTrkpass_dxy","hDisTrkpass_dxy",10,0.0,30.0);
   TH1F *hDisTrktotal_dxy =new TH1F("hDisTrktotal_dxy","hDisTrktotal_dxy",10,0.0,30.0);


   TH1F *hDisGlbpass_lxy =new TH1F("hDisGlbpass_lxy","hDisGlbpass_lxy",10,0.0,30.0);
   TH1F *hDisGlbtotal_lxy =new TH1F("hDisGlbtotal_lxy","hDisGlbtotal_lxy",10,0.0,30.0);

   TH1F *hDisStapass_lxy =new TH1F("hDisStapass_lxy","hDisStapass_lxy",10,0.0,30.0);
   TH1F *hDisStatotal_lxy =new TH1F("hDisStatotal_lxy","hDisStatotal_lxy",10,0.0,30.0);

   TH1F *hDisTrkpass_lxy =new TH1F("hDisTrkpass_lxy","hDisTrkpass_lxy",10,0.0,30.0);
   TH1F *hDisTrktotal_lxy =new TH1F("hDisTrktotal_lxy","hDisTrktotal_lxy",10,0.0,30.0);


   TH1F *hgen_pt =new TH1F("hgen_pt","hgen_pt",10,0,1000);
   TH1F *hgen_eta =new TH1F("hgen_eta","hgen_eta",10,-3,3);
   TH1F *hgen_dxy =new TH1F("hgen_dxy","hgen_dxy",10,0.0,30.0);
   TH1F *hgen_lxy =new TH1F("hgen_lxy","hgen_lxy",10,0.0,30.0);


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
   //hNgen->Fill(pt_gen->size());
   //hNdisglb->Fill(pt_DisglbMuon->size());
   //hNdisSta->Fill(pt_DisSta->size());

   for(int i=0;i<pt_gen->size();i++){//gen loop
   //hPtgen->Fill(pt_gen->at(i));

   hDisGlbtotal_pt->Fill(pt_gen->at(i));
   hDisStatotal_pt->Fill(pt_gen->at(i));
   hDisTrktotal_pt->Fill(pt_gen->at(i));

   hDisGlbtotal_eta->Fill(eta_gen->at(i));
   hDisStatotal_eta->Fill(eta_gen->at(i));
   hDisTrktotal_eta->Fill(eta_gen->at(i));

   hDisGlbtotal_dxy->Fill(dxy_gen->at(i));
   hDisStatotal_dxy->Fill(dxy_gen->at(i));
   hDisTrktotal_dxy->Fill(dxy_gen->at(i));

   hDisGlbtotal_lxy->Fill(lxy_gen->at(i));
   hDisStatotal_lxy->Fill(lxy_gen->at(i));
   hDisTrktotal_lxy->Fill(lxy_gen->at(i));

   double deltaR = 10;
   //double deltaR2 = 10;


  TLorentzVector genmu;

  genmu.SetPtEtaPhiM(pt_gen->at(i),eta_gen->at(i),phi_gen->at(i),0);

  for(int j=0;j<pt_DisglbMuon->size();j++){//glb

   TLorentzVector glbMu;
   glbMu.SetPtEtaPhiM(pt_DisglbMuon->at(j), eta_DisglbMuon->at(j),phi_DisglbMuon->at(j),0);
   deltaR=genmu.DeltaR(glbMu);


      if(deltaR<0.1){//dr

        hDisGlbpass_pt->Fill(pt_gen->at(i));
        hDisGlbpass_eta->Fill(eta_gen->at(i));
        hDisGlbpass_dxy->Fill(dxy_gen->at(i));
        hDisGlbpass_lxy->Fill(lxy_gen->at(i));

      } //dr


    }//glb 


  for(int k=0;k<pt_DisSta->size();k++){//sta loop

    TLorentzVector staMu;
    staMu.SetPtEtaPhiM(pt_DisSta->at(k), eta_DisSta->at(k),phi_DisSta->at(k),0);
    deltaR=genmu.DeltaR(staMu);

    ///double relptdiff=(pt_DisSta->at(k)-pt_gen->at(i))/pt_gen->at(i);


    //if((relptdiff<0.20))hDRgenSta->Fill(deltaR2);

    if(deltaR<0.1){//dr


   hDisStapass_pt->Fill(pt_gen->at(i));
   hDisStapass_eta->Fill(eta_gen->at(i));
   hDisStapass_dxy->Fill(dxy_gen->at(i));
   hDisStapass_lxy->Fill(lxy_gen->at(i));

    }//dr


   }//sta loop



  for(int l=0;l<pt_DisTrk->size();l++){//Trk loop

    TLorentzVector TrkMu;
    TrkMu.SetPtEtaPhiM(pt_DisTrk->at(l), eta_DisTrk->at(l),phi_DisTrk->at(l),0);
    deltaR=genmu.DeltaR(TrkMu);

    //double relptdiff=(pt_DisTrk->at(l)-pt_gen->at(i))/pt_gen->at(i);


    //if((relptdiff<0.20))hDRgenSta->Fill(deltaR2);

    if(deltaR<0.1){//dr


   hDisTrkpass_pt->Fill(pt_gen->at(i));
   hDisTrkpass_eta->Fill(eta_gen->at(i));
   hDisTrkpass_dxy->Fill(dxy_gen->at(i));
   hDisTrkpass_lxy->Fill(lxy_gen->at(i));

    }//dr


   }//trk loop


   }//gen loop



   }

  TCanvas *C=new TCanvas("C","C",1200,900);
  C->Divide(4,3);
  C->SetGrid();

  C->cd(1);
  TEfficiency *Eff_glb_pt = new TEfficiency(*hDisGlbpass_pt,*hDisGlbtotal_pt);
  Eff_glb_pt->SetTitle("p_{T} vs Efficiency;p_{T};#epsilon");
  Eff_glb_pt->SetMarkerColor(8);//97
  Eff_glb_pt->SetMarkerStyle(20);
  Eff_glb_pt->SetMarkerSize(1.0);
  Eff_glb_pt->Draw("AP");
  
  TLegend* l1 = new TLegend(.40,.35,.85,.45);
  l1->AddEntry(Eff_glb_pt,"displacedGlobalMuons","P");
  l1->Draw();


  C->cd(2);
  TEfficiency *Eff_glb_eta = new TEfficiency(*hDisGlbpass_eta,*hDisGlbtotal_eta);
  Eff_glb_eta->SetTitle("#eta vs Efficiency;#eta;#epsilon");
  Eff_glb_eta->SetMarkerColor(8);//97
  Eff_glb_eta->SetMarkerStyle(20);
  Eff_glb_eta->SetMarkerSize(1.0);
  Eff_glb_eta->Draw("AP");

  TLegend* l2 = new TLegend(.40,.35,.85,.45);
  l2->AddEntry(Eff_glb_eta,"displacedGlobalMuons","P");
  l2->Draw();

  C->cd(3);
  TEfficiency *Eff_glb_dxy = new TEfficiency(*hDisGlbpass_dxy,*hDisGlbtotal_dxy);
  Eff_glb_dxy->SetTitle("d_{xy} vs Efficiency;d_{xy};#epsilon");
  Eff_glb_dxy->SetMarkerColor(8);//97
  Eff_glb_dxy->SetMarkerStyle(20);
  Eff_glb_dxy->SetMarkerSize(1.0);
  Eff_glb_dxy->Draw("AP");

  TLegend* l3 = new TLegend(.40,.35,.85,.45);
  l3->AddEntry(Eff_glb_dxy,"displacedGlobalMuons","P");
  l3->Draw();


  C->cd(4);
  TEfficiency *Eff_glb_lxy = new TEfficiency(*hDisGlbpass_lxy,*hDisGlbtotal_lxy);
  Eff_glb_lxy->SetTitle("l_{xy} vs Efficiency;l_{xy};#epsilon");
  Eff_glb_lxy->SetMarkerColor(8);//97
  Eff_glb_lxy->SetMarkerStyle(20);
  Eff_glb_lxy->SetMarkerSize(1.0);
  Eff_glb_lxy->Draw("AP");

  TLegend* l4 = new TLegend(.40,.35,.85,.45);
  l4->AddEntry(Eff_glb_lxy,"displacedGlobalMuons","P");
  l4->Draw();


  C->cd(5);
  TEfficiency *Eff_sta_pt = new TEfficiency(*hDisStapass_pt,*hDisStatotal_pt);
  Eff_sta_pt->SetTitle("p_{T} vs Efficiency;p_{T};#epsilon");
  Eff_sta_pt->SetMarkerColor(7);//97
  Eff_sta_pt->SetMarkerStyle(20);
  Eff_sta_pt->SetMarkerSize(1.0);
  Eff_sta_pt->Draw("AP");

  TLegend* l5 = new TLegend(.40,.35,.85,.45);
  l5->AddEntry(Eff_sta_pt,"displacedSAMuons","P");
  l5->Draw();
  
  C->cd(6);
  TEfficiency *Eff_sta_eta = new TEfficiency(*hDisStapass_eta,*hDisStatotal_eta);
  Eff_sta_eta->SetTitle("#eta vs Efficiency;#eta;#epsilon");
  Eff_sta_eta->SetMarkerColor(7);//97
  Eff_sta_eta->SetMarkerStyle(20);
  Eff_sta_eta->SetMarkerSize(1.0);
  Eff_sta_eta->Draw("AP");

  TLegend* l6 = new TLegend(.40,.35,.85,.45);
  l6->AddEntry(Eff_sta_eta,"displacedSAMuons","P");
  l6->Draw();

  C->cd(7);
  TEfficiency *Eff_sta_dxy = new TEfficiency(*hDisStapass_dxy,*hDisStatotal_dxy);
  Eff_sta_dxy->SetTitle("d_{xy} vs Efficiency;d_{xy};#epsilon");
  Eff_sta_dxy->SetMarkerColor(7);//97
  Eff_sta_dxy->SetMarkerStyle(20);
  Eff_sta_dxy->SetMarkerSize(1.0);
  Eff_sta_dxy->Draw("AP");

  TLegend* l7 = new TLegend(.40,.35,.85,.45);
  l7->AddEntry(Eff_sta_dxy,"displacedSAMuons","P");
  l7->Draw();

  C->cd(8);
  TEfficiency *Eff_sta_lxy = new TEfficiency(*hDisStapass_lxy,*hDisStatotal_lxy);
  Eff_sta_lxy->SetTitle("l_{xy} vs Efficiency;l_{xy};#epsilon");
  Eff_sta_lxy->SetMarkerColor(7);//97
  Eff_sta_lxy->SetMarkerStyle(20);
  Eff_sta_lxy->SetMarkerSize(1.0);
  Eff_sta_lxy->Draw("AP");

  TLegend* l8 = new TLegend(.40,.35,.85,.45);
  l8->AddEntry(Eff_sta_lxy,"displacedSAMuons","P");
  l8->Draw();

  C->cd(9);
  TEfficiency *Eff_Trk_pt = new TEfficiency(*hDisTrkpass_pt,*hDisTrktotal_pt);
  Eff_Trk_pt->SetTitle("p_{T} vs Efficiency;p_{T};#epsilon");
  Eff_Trk_pt->SetMarkerColor(4);//97
  Eff_Trk_pt->SetMarkerStyle(20);
  Eff_Trk_pt->SetMarkerSize(1.0);
  Eff_Trk_pt->Draw("AP");

  TLegend* l9 = new TLegend(.40,.35,.85,.45);
  l9->AddEntry(Eff_Trk_pt,"displacedTRK","P");
  l9->Draw();

  C->cd(10);
  TEfficiency *Eff_Trk_eta = new TEfficiency(*hDisTrkpass_eta,*hDisTrktotal_eta);
  Eff_Trk_eta->SetTitle("#eta vs Efficiency;#eta;#epsilon");
  Eff_Trk_eta->SetMarkerColor(4);//97
  Eff_Trk_eta->SetMarkerStyle(20);
  Eff_Trk_eta->SetMarkerSize(1.0);
  Eff_Trk_eta->Draw("AP");

  TLegend* l10 = new TLegend(.40,.35,.85,.45);
  l10->AddEntry(Eff_Trk_eta,"displacedTRK","P");
  l10->Draw();


  C->cd(11);
  TEfficiency *Eff_Trk_dxy = new TEfficiency(*hDisTrkpass_dxy,*hDisTrktotal_dxy);
  Eff_Trk_dxy->SetTitle("d_{xy} vs Efficiency;d_{xy};#epsilon");
  Eff_Trk_dxy->SetMarkerColor(4);//97
  Eff_Trk_dxy->SetMarkerStyle(20);
  Eff_Trk_dxy->SetMarkerSize(1.0);
  Eff_Trk_dxy->Draw("AP");

  TLegend* l11 = new TLegend(.40,.35,.85,.45);
  l11->AddEntry(Eff_Trk_dxy,"displacedTRK","P");
  l11->Draw();

  C->cd(12);
  TEfficiency *Eff_Trk_lxy = new TEfficiency(*hDisTrkpass_lxy,*hDisTrktotal_lxy);
  Eff_Trk_lxy->SetTitle("l_{xy} vs Efficiency;l_{xy};#epsilon");
  Eff_Trk_lxy->SetMarkerColor(4);//97
  Eff_Trk_lxy->SetMarkerStyle(20);
  Eff_Trk_lxy->SetMarkerSize(1.0);
  Eff_Trk_lxy->Draw("AP");

  TLegend* l12 = new TLegend(.40,.35,.85,.45);
  l12->AddEntry(Eff_Trk_lxy,"displacedTRK","P");
  l12->Draw();


}
