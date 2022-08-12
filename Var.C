void Var()
{


    TFile *f=new TFile("RelVal_DisplacedSUSY_NO_PU_Lxy_added.root","UPDATE");

    TTree *myt=(TTree *)f->Get("ntuple/tree");


   TH1F *hDisGlbtotal_pt =new TH1F("hDisGlbtotal_pt","DisGlb p_{T}",10,0,1000);

   TH1F *hDisStatotal_pt =new TH1F("hDisStatotal_pt","DisSta p_{T}",10,0,1000);

   TH1F *hDisTrktotal_pt =new TH1F("hDisTrktotal_pt","DisTrk p_{T}",10,0,1000);

   TH1F *hgen_pt =new TH1F("hgen_pt","gen p_{T}",10,0,1000);


   TH1F *hDisGlbtotal_eta =new TH1F("hDisGlbtotal_eta","DisGlb #eta",30,-3,3);

   TH1F *hDisStatotal_eta =new TH1F("hDisStatotal_eta","DisSta #eta",30,-3,3);

   TH1F *hDisTrktotal_eta =new TH1F("hDisTrktotal_eta","DisTrk #eta",30,-3,3);

   TH1F *hgen_eta =new TH1F("hgen_eta","gen #eta",30,-3,3);


   TH1F *hDisGlbtotal_dxy =new TH1F("hDisGlbtotal_dxy","DisGlb d_{xy}",30,0.0,30.0);

   TH1F *hDisStatotal_dxy =new TH1F("hDisStatotal_dxy","DisSta d_{xy}",30,0.0,30.0);

   TH1F *hDisTrktotal_dxy =new TH1F("hDisTrktotal_dxy","DisTrk d_{xy}",30,0.0,30.0);

   TH1F *hgen_dxy =new TH1F("hgen_dxy","gen d_{xy}",30,0.0,30.0);


   TH1F *hDisGlbtotal_lxy =new TH1F("hDisGlbtotal_lxy","DisGlb l_{xy}",30,0.0,30.0);

   TH1F *hDisStatotal_lxy =new TH1F("hDisStatotal_lxy","DisSta l_{xy}",30,0.0,30.0);

   TH1F *hDisTrktotal_lxy =new TH1F("hDisTrktotal_lxy","DisTrk l_{xy}",30,0.0,30.0);

   TH1F *hgen_lxy =new TH1F("hgen_lxy","gen l_{xy}",30,0.0,30.0);


   myt->Project("hDisGlbtotal_pt","pt_DisglbMuon");
   myt->Project("hDisStatotal_pt","pt_DisSta");
   myt->Project("hDisTrktotal_pt","pt_DisTrk");
   myt->Project("hgen_pt","pt_gen");

   myt->Project("hDisGlbtotal_eta","eta_DisglbMuon");
   myt->Project("hDisStatotal_eta","eta_DisSta");
   myt->Project("hDisTrktotal_eta","eta_DisTrk");
   myt->Project("hgen_eta","eta_gen");

   myt->Project("hDisGlbtotal_dxy","dxy_DisglbMuon");
   myt->Project("hDisStatotal_dxy","dxy_DisSta");
   myt->Project("hDisTrktotal_dxy","dxy_DisTrk");
   myt->Project("hgen_dxy","dxy_gen");

   myt->Project("hDisGlbtotal_lxy","lxy_DisglbMuon");
   myt->Project("hDisStatotal_lxy","lxy_DisSta");
   myt->Project("hDisTrktotal_lxy","lxy_DisTrk");
   myt->Project("hgen_lxy","lxy_gen");


  TCanvas *c = new TCanvas("c","c",1200,1200);
  TStyle* mcStyle = new TStyle("mcStyle","Manuel's Root Styles");
  gStyle->SetOptStat(0);

  c->Divide(4,4);

  c->cd(1);
  hDisGlbtotal_pt->Draw("e1");

  c->cd(2);
  hDisGlbtotal_eta->Draw("e1");

  c->cd(3);
  hDisGlbtotal_dxy->Draw("e1");

  c->cd(4);
  hDisGlbtotal_lxy->Draw("e1");

  c->cd(5);
  hDisStatotal_pt->Draw("e1");

  c->cd(6);
  hDisStatotal_eta->Draw("e1");

  c->cd(7);
  hDisStatotal_dxy->Draw("e1");

  c->cd(8);
  hDisStatotal_lxy->Draw("e1");

  c->cd(9);
  hDisTrktotal_pt->Draw("e1");

  c->cd(10);
  hDisTrktotal_eta->Draw("e1");

  c->cd(11);
  hDisTrktotal_dxy->Draw("e1");

  c->cd(12);
  hDisTrktotal_lxy->Draw("e1");

  c->cd(13);
  hgen_pt->Draw("e1");

  c->cd(14);
  hgen_eta->Draw("e1");

  c->cd(15);
  hgen_dxy->Draw("e1");

  c->cd(16);
  hgen_lxy->Draw("e1");

}
