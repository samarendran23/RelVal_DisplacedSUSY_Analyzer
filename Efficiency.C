TCanvas *c[12];

const char*  collections[3] = {"displacedGlobalMuon", "displacedSAMuons", "displacedTracks"};

const char* Var_name[12] = {"pt_DisglbMuon","pt_DisSta","pt_DisTrk","eta_DisglbMuon","eta_DisSta","eta_DisTrk","dxy_DisglbMuon","dxy_DisSta","dxy_DisTrk","lxy_DisglbMuon","lxy_DisSta","lxy_DisTrk"};

double L_bound[12] = {0.0,0.0,0.0,-3.0,-3.0,-3.0,0.0,0.0,0.0,0.0,0.0,0.0};
double U_bound[12] = {1000.0,1000.0,1000.0,3.0,3.0,3.0,30.0,30.0,30.0,30.0,30.0,30.0};

void plot(TH1D *htotal, TH1D *hpass, const char* str, int i)
{


 cout<<"Hi I am Entering Plot Function !!"<<endl;

htotal->SetTitle("");
htotal->GetXaxis()->SetTitle(str);
htotal->GetYaxis()->SetTitle("Entries");

htotal->SetMarkerStyle(20);
htotal->SetMarkerColor(kRed);
htotal->SetMarkerSize(1.0);
htotal->SetFillStyle(3004);
htotal->SetFillColor(kRed);
htotal->SetLineWidth(1);
htotal->SetLineStyle(2);
htotal->SetLineColor(kRed);

hpass->SetMarkerStyle(21);
hpass->SetMarkerColor(kBlue);
hpass->SetMarkerSize(1.0);
hpass->SetFillStyle(3004);
hpass->SetFillColor(kBlue);
hpass->SetLineStyle(2);
hpass->SetLineWidth(1);
hpass->SetLineColor(kBlue);


TLegend *legend1=new TLegend(.50,.75,.85,.88);
legend1->SetFillColor(0);
legend1->SetTextSize(0.0);
legend1->AddEntry(htotal,"w/o Truth-Matching","FL");
legend1->AddEntry(hpass,"with Truth-Matching","FL");


c[i] = new TCanvas("c","c", 1400, 800);
gStyle->SetOptStat(0);

c[i]->Divide(2,1);

c[i]->cd(1);

htotal->Draw("E1 HIST");
hpass->Draw("E1 HIST SAME");
legend1->Draw();

c[i]->cd(2);

  TEfficiency *Eff = new TEfficiency(*hpass,*htotal);
  //Eff->SetTitle(;srt;"Fake Rate");
  Eff->SetMarkerColor(8);//97
  Eff->SetMarkerStyle(20);
  Eff->SetMarkerSize(1.4);
  Eff->Draw("AP");

  c[i]->SaveAs(Form("./Truth_Matchinh_Eff%d_%s.png",i,Var_name[i]));
  c[i]->SaveAs(Form("./Truth_Matchinh_Eff%d_%s.pdf",i,Var_name[i]));
  delete c[i];


}

void Efficiency()
{


 TFile *f=new TFile("truthmatch_added_wo_pu.root","UPDATE");

 TTree *myt=(TTree *)f->Get("ntuple/tree");

/* TCut TM_glb  = "istrue_DisglbMuon==1";

 TCut TM_sa  = "istrue_DisSta==1";

 TCut TM_trk  = "istrue_DisTrk==1";
*/

 TCut TM[3] = {"istrue_DisglbMuon==1","istrue_DisSta==1","istrue_DisTrk==1"};

const char* var_name[12] = {"pt_DisglbMuon","pt_DisSta","pt_DisTrk","eta_DisglbMuon","eta_DisSta","eta_DisTrk","dxy_DisglbMuon","dxy_DisSta","dxy_DisTrk","lxy_DisglbMuon","lxy_DisSta","lxy_DisTrk"};

const char* title[12] = {"pt_DisglbMuon","pt_DisSta","pt_DisTrk","eta_DisglbMuon","eta_DisSta","eta_DisTrk","dxy_DisglbMuon","dxy_DisSta","dxy_DisTrk","lxy_DisglbMuon","lxy_DisSta","lxy_DisTrk"};

double l_bound[12] = {0.0,0.0,0.0,-3.0,-3.0,-3.0,0.0,0.0,0.0,0.0,0.0,0.0};
double u_bound[12] = {1000.0,1000.0,1000.0,3.0,3.0,3.0,30.0,30.0,30.0,30.0,30.0,30.0};

int n = 10; // number of bins

TH1D* hpass[3] ;
TH1D* htotal[3] ;

const char* hist_name_pass[3] = {"hpass_glb" , "hpass_sa", "hpass_trk"};

const char* hist_name_total[3] = {"htotal_glb" , "htotal_sa", "htotal_trk"};


const char* hist_title_pass[3] = {"hpass_glb" , "hpass_sa", "hpass_trk"};

const char* hist_title_total[3] = {"htotal_glb" , "htotal_sa", "htotal_trk"};


for(int i=0;i<3;i++)
        {
                hpass[i] = new TH1D(hist_name_pass[i],hist_title_pass[i], n, l_bound[i], u_bound[i]);

                htotal[i] = new TH1D(hist_name_total[i],hist_title_total[i], n, l_bound[i], u_bound[i]);

                for(int j=0;i<12;i++){

                myt->Project(hist_name_pass[i], var_name[j], TM[i]);

                myt->Project(hist_name_total[i], var_name[j]);

                plot(hpass[i], htotal[i],title[j], j);

                delete hpass[i];
                delete htotal[i];

		}

        }


/*for(int i=0;i<3;i++)
        {
                h_total[i] = new TH1D(hist_name_total[i],hist_title_total[i], n, l_bound[i], u_bound[i]);

                for(int j=0;i<12;i++){

                myt->Project(h_total[i], var_name[j]);
 
        }


}*/


}
