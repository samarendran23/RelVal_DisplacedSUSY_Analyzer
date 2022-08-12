// -*- C++ -*-
//
// Package:    Analysis/Analyzer
// Class:      Analyzer
//
/**\class Analyzer Analyzer.cc Analysis/Analyzer/plugins/Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bibhuprasad Mahakud
//         Created:  Mon, 06 Sep 2021 18:04:25 GMT
// Edited : Samarendra Nayak
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoMuon/MuonIdentification/plugins/MuonIdProducer.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include<vector>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>

using namespace std;

using namespace edm;

using namespace reco;

const double PI = 3.141592653589793;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

using reco::GenParticleCollection;


class Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Analyzer(const edm::ParameterSet&);
      ~Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      double calPhi(double, double, double);

      double calEta(double, double, double);

      double calEtaPhiDistance (double, double, double, double, double, double);

      void calLxy (double, double, double, double, double, double, double*);


      // ----------member data ---------------------------
     // edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file

      edm::Service<TFileService> fs;

      //edm::EDGetTokenT<TrackCollection> Tracks_global_;

      edm::EDGetTokenT<TrackCollection> Tracks_Disglobal_;

      //edm::EDGetTokenT<TrackCollection> Tracks_standalone_;

      edm::EDGetTokenT<TrackCollection> Tracks_Disstandalone_;

      edm::EDGetTokenT<TrackCollection> displacedTracksToken_;

      //edm::EDGetTokenT<TrackCollection>  displacedOutInTracks_;

      //edm::EDGetTokenT<reco::MuonCollection>  earlyDisplacedMuonsToken_;

      edm::EDGetTokenT<GenParticleCollection> gen_Muon_;

      //edm::EDGetTokenT<std::vector<reco::Muon>> offlineMuonToken_;

      //edm::EDGetTokenT<std::vector<reco::Muon>> offlinedisplacedMuonToken_;


      //edm::EDGetTokenT<int> numOIseeds_;

      //edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
 
     //edm::EDGetTokenT<std::vector<MuonIdProducer>> Muonid_;


      reco::BeamSpot beamSpot_;

      double TruthMatchMuonMaxR_;



      TTree* tree;

      //unsigned int nprivtx;

      vector<double> *pt_glbMuon;

      vector<double> *eta_glbMuon;

      vector<double> *phi_glbMuon;

      vector<double> *dxy_glbMuon;

      vector<double> *dz_glbMuon;

      vector<double> *d0_glbMuon;

      //vector<double> *lxy_glbMuon;


      vector<bool> *istrue_DisglbMuon;

      vector<double> *px_DisglbMuon;

      vector<double> *py_DisglbMuon;

      vector<double> *pz_DisglbMuon;

      vector<double> *pt_DisglbMuon;

      vector<double> *eta_DisglbMuon;

      vector<double> *phi_DisglbMuon;

      vector<double> *dxy_DisglbMuon;

      vector<double> *dz_DisglbMuon;

      vector<double> *d0_DisglbMuon;

      vector<double> *lxy_DisglbMuon;


      vector<double> *pt_Sta;

      vector<double> *eta_Sta;

      vector<double> *phi_Sta;

      vector<double> *dxy_Sta;

      vector<double> *dz_Sta;

      vector<double> *d0_Sta;


      vector<bool> *istrue_DisSta;

      vector<double> *px_DisSta;

      vector<double> *py_DisSta;

      vector<double> *pz_DisSta;

      vector<double> *pt_DisSta;

      vector<double> *eta_DisSta;

      vector<double> *phi_DisSta;

      vector<double> *dxy_DisSta;

      vector<double> *dz_DisSta;

      vector<double> *d0_DisSta;

      vector<double> *lxy_DisSta;


      vector<bool> *istrue_DisTrk;

      vector<double> *px_DisTrk;

      vector<double> *py_DisTrk;

      vector<double> *pz_DisTrk;

      vector<double> *pt_DisTrk;

      vector<double> *eta_DisTrk;

      vector<double> *phi_DisTrk;

      vector<double> *dxy_DisTrk;

      vector<double> *dz_DisTrk;

      vector<double> *d0_DisTrk;

      vector<double> *lxy_DisTrk;


      vector<double> *pt_DisOutIn;

      vector<double> *eta_DisOutIn;

      vector<double> *phi_DisOutIn;

      vector<double> *dxy_DisOutIn;

      vector<double> *dz_DisOutIn;

      vector<double> *d0_DisOutIn;


      vector<double> *pt_DisEar;

      vector<double> *eta_DisEar;

      vector<double> *phi_DisEar;

      vector<double> *dxy_DisEar;

      vector<double> *dz_DisEar;

      vector<double> *d0_DisEar;


      double px_gen;

      double py_gen;

      double pz_gen;

      vector<double> *pt_gen;

      vector<double> *eta_gen;

      vector<double> *phi_gen;

      vector<double> *vx_gen;

      vector<double> *vy_gen;

      vector<double> *vz_gen;

      vector<double> *vertex_genx;

      vector<double> *vertex_geny;

      vector<double> *dxy_gen;

      vector<double> *lxy_gen;

      //vector<double> *Muon;

      //vector<double> *displacedMuon;


      vector<double> *pt_Muon;

      vector<double> *eta_Muon;

      vector<double> *phi_Muon;

      vector<double> *mass_Muon;

      vector<double> *dxy_Muon;

      vector<double> *dz_Muon;

      vector<double> *d0_Muon;


      vector<double> *pt_displacedMuon;

      vector<double> *eta_displacedMuon;

      vector<double> *phi_displacedMuon;

      vector<double> *mass_displacedMuon;

      vector<double> *dxy_displacedMuon;

      vector<double> *dz_displacedMuon;

      vector<double> *d0_displacedMuon;

 
      int numberOfOISeeds;

      //int Tracks;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig)
 :
 // tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),

 //Tracks_global_(consumes<TrackCollection>(iConfig.getParameter<edm::InputTag>("globaltracks"))),

 Tracks_Disglobal_(consumes<TrackCollection>(iConfig.getParameter<edm::InputTag>("displacedglobaltracks"))),

 //Tracks_standalone_(consumes<TrackCollection>(iConfig.getParameter<edm::InputTag>("standalonetracks"))),

 Tracks_Disstandalone_(consumes<TrackCollection>(iConfig.getParameter<edm::InputTag>("displacedstandalonetracks"))),

 displacedTracksToken_(consumes<TrackCollection>(iConfig.getParameter<edm::InputTag>("displacedtracks"))),

 //displacedOutInTracks_(consumes<TrackCollection>(iConfig.getParameter<edm::InputTag>("displacedoutintracks"))),

 //earlyDisplacedMuonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("earlydisplacedmuons"))),

 gen_Muon_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticle"))),

 //offlineMuonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),

 //offlinedisplacedMuonToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("displacedmuons"))),

 //numOIseeds_(consumes<int>(iConfig.getParameter<edm::InputTag>("numoiseeds"))),

 //tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks")))

 //Muonid_(consumes<std::vector<MuonIdProducer>>(iConfig.getParameter<edm::InputTag>("muonid")))

 TruthMatchMuonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchMuonMaxR")),



  pt_glbMuon(0),

  eta_glbMuon(0),

  phi_glbMuon(0),

  dxy_glbMuon(0),

  dz_glbMuon(0),

  d0_glbMuon(0),

  //lxy_glbMuon(0),


  istrue_DisglbMuon(0),

  px_DisglbMuon(0),

  py_DisglbMuon(0),

  pz_DisglbMuon(0),

  pt_DisglbMuon(0),

  eta_DisglbMuon(0),

  phi_DisglbMuon(0),

  dxy_DisglbMuon(0),

  dz_DisglbMuon(0),

  d0_DisglbMuon(0),

  lxy_DisglbMuon(0),


  pt_Sta(0),

  eta_Sta(0),

  phi_Sta(0),

  dxy_Sta(0),

  dz_Sta(0),

  d0_Sta(0),


  istrue_DisSta(0),

  px_DisSta(0),

  py_DisSta(0),

  pz_DisSta(0),

  pt_DisSta(0),

  eta_DisSta(0),

  phi_DisSta(0),

  dxy_DisSta(0),

  dz_DisSta(0),

  d0_DisSta(0),

  lxy_DisSta(0),


  istrue_DisTrk(0),

  px_DisTrk(0),

  py_DisTrk(0),

  pz_DisTrk(0),

  pt_DisTrk(0),

  eta_DisTrk(0),

  phi_DisTrk(0),

  dxy_DisTrk(0),

  dz_DisTrk(0),

  d0_DisTrk(0),

  lxy_DisTrk(0),


  pt_DisOutIn(0),

  eta_DisOutIn(0),

  phi_DisOutIn(0),

  dxy_DisOutIn(0),

  dz_DisOutIn(0),

  d0_DisOutIn(0),


  pt_DisEar(0),

  eta_DisEar(0),

  phi_DisEar(0),

  dxy_DisEar(0),

  dz_DisEar(0),

  d0_DisEar(0),


  px_gen(0),

  py_gen(0),

  pz_gen(0),

  pt_gen(0),

  eta_gen(0),

  phi_gen(0),

  vx_gen(0),

  vy_gen(0),

  vz_gen(0),

  vertex_genx(0),

  vertex_geny(0),

  dxy_gen(0),

  lxy_gen(0),

  //Muon(0),

  //displacedMuon(0)


  pt_Muon(0),

  eta_Muon(0),

  phi_Muon(0),

  mass_Muon(0),

  dxy_Muon(0),

  dz_Muon(0),

  d0_Muon(0),


  pt_displacedMuon(0),

  eta_displacedMuon(0),

  phi_displacedMuon(0),

  mass_displacedMuon(0),

  dxy_displacedMuon(0),

  dz_displacedMuon(0),

  d0_displacedMuon(0)

  //numberOfOISeeds(0)

  //Tracks(0)

{
   //now do what ever initialization is needed

   tree = fs->make<TTree>("tree","tree");

   tree->Branch("pt_glbMuon",&pt_glbMuon);

   tree->Branch("eta_glbMuon",&eta_glbMuon);

   tree->Branch("phi_glbMuon",&phi_glbMuon);

   tree->Branch("dxy_glbMuon",&dxy_glbMuon);

   tree->Branch("dz_glbMuon",&dz_glbMuon);

   tree->Branch("d0_glbMuon",&d0_glbMuon);

   //tree->Branch("lxy_glbMuon",&lxy_glbMuon);


   tree->Branch("istrue_DisglbMuon",&istrue_DisglbMuon);

   tree->Branch("px_DisglbMuon",&px_DisglbMuon);

   tree->Branch("py_DisglbMuon",&py_DisglbMuon);

   tree->Branch("pz_DisglbMuon",&pz_DisglbMuon);

   tree->Branch("pt_DisglbMuon",&pt_DisglbMuon);

   tree->Branch("eta_DisglbMuon",&eta_DisglbMuon);

   tree->Branch("phi_DisglbMuon",&phi_DisglbMuon);

   tree->Branch("dxy_DisglbMuon",&dxy_DisglbMuon);

   tree->Branch("dz_DisglbMuon",&dz_DisglbMuon);

   tree->Branch("d0_DisglbMuon",&d0_DisglbMuon);

   tree->Branch("lxy_DisglbMuon",&lxy_DisglbMuon);


   tree->Branch("pt_Sta",&pt_Sta);

   tree->Branch("eta_Sta",&eta_Sta);

   tree->Branch("phi_Sta",&phi_Sta);

   tree->Branch("dxy_Sta",&dxy_Sta);

   tree->Branch("dz_Sta",&dz_Sta);

   tree->Branch("d0_Sta",&d0_Sta);


   tree->Branch("istrue_DisSta",&istrue_DisSta);

   tree->Branch("px_DisSta",&px_DisSta);

   tree->Branch("py_DisSta",&py_DisSta);

   tree->Branch("pz_DisSta",&pz_DisSta);

   tree->Branch("pt_DisSta",&pt_DisSta);

   tree->Branch("eta_DisSta",&eta_DisSta);

   tree->Branch("phi_DisSta",&phi_DisSta);

   tree->Branch("dxy_DisSta",&dxy_DisSta);

   tree->Branch("dz_DisSta",&dz_DisSta);

   tree->Branch("d0_DisSta",&d0_DisSta);

   tree->Branch("lxy_DisSta",&lxy_DisSta);


   tree->Branch("istrue_DisTrk",&istrue_DisTrk);

   tree->Branch("px_DisTrk",&px_DisTrk);

   tree->Branch("py_DisTrk",&py_DisTrk);

   tree->Branch("pz_DisTrk",&pz_DisTrk);

   tree->Branch("pt_DisTrk",&pt_DisTrk);

   tree->Branch("eta_DisTrk",&eta_DisTrk);

   tree->Branch("phi_DisTrk",&phi_DisTrk);

   tree->Branch("dxy_DisTrk",&dxy_DisTrk);

   tree->Branch("dz_DisTrk",&dz_DisTrk);

   tree->Branch("d0_DisTrk",&d0_DisTrk);

   tree->Branch("lxy_DisTrk",&lxy_DisTrk);


   tree->Branch("pt_DisOutIn",&pt_DisOutIn);

   tree->Branch("eta_DisOutIn",&eta_DisOutIn);

   tree->Branch("phi_DisOutIn",&phi_DisOutIn);

   tree->Branch("dxy_DisOutIn",&dxy_DisOutIn);

   tree->Branch("dz_DisOutIn",&dz_DisOutIn);

   tree->Branch("d0_DisOutIn",&d0_DisOutIn);


   tree->Branch("pt_DisEar",&pt_DisEar);

   tree->Branch("eta_DisEar",&eta_DisEar);

   tree->Branch("phi_DisEar",&phi_DisEar);

   tree->Branch("dxy_DisEar",&dxy_DisEar);

   tree->Branch("dz_DisEar",&dz_DisEar);

   tree->Branch("d0_DisEar",&d0_DisEar);


   tree->Branch("px_gen",&px_gen);

   tree->Branch("py_gen",&py_gen);

   tree->Branch("pz_gen",&pz_gen);

   tree->Branch("pt_gen",&pt_gen);

   tree->Branch("eta_gen",&eta_gen);

   tree->Branch("phi_gen",&phi_gen);

   tree->Branch("vx_gen",&vx_gen);

   tree->Branch("vy_gen",&vy_gen);

   tree->Branch("vz_gen",&vz_gen);

   tree->Branch("vertex_genx",&vertex_genx);

   tree->Branch("vertex_geny",&vertex_geny);

   tree->Branch("dxy_gen",&dxy_gen);

   tree->Branch("lxy_gen",&lxy_gen);


   //tree->Branch("Muon",&Muon);

   //tree->Branch("displacedMuon",&displacedMuon);
   

   tree->Branch("pt_Muon",&pt_Muon);

   tree->Branch("eta_Muon",&eta_Muon);

   tree->Branch("phi_Muon",&phi_Muon);

   tree->Branch("mass_Muon",&mass_Muon);

   tree->Branch("dxy_Muon",&dxy_Muon);

   tree->Branch("dz_Muon",&dz_Muon);

   tree->Branch("d0_Muon",&d0_Muon);


   tree->Branch("pt_displacedMuon",&pt_displacedMuon);

   tree->Branch("eta_displacedMuon",&eta_displacedMuon);

   tree->Branch("phi_displacedMuon",&phi_displacedMuon);

   tree->Branch("mass_displacedMuon",&mass_displacedMuon);

   tree->Branch("dxy_displacedMuon",&dxy_displacedMuon);

   tree->Branch("dz_displacedMuon",&dz_displacedMuon);

   tree->Branch("d0_displacedMuon",&d0_displacedMuon);


   tree->Branch("numberOfOISeeds",&numberOfOISeeds);

   //tree->Branch("Tracks",&Tracks);


}


Analyzer::~Analyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void
Analyzer::calLxy (double Vx, double Vy, double Vz,
                  double Wx, double Wy, double Wz,
                  double* varLxy)
{
  *varLxy = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
}

double
Analyzer::calPhi (double Px, double Py, double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}


double
Analyzer::calEta (double Px, double Py, double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
 	 return 0.5*log((P + Pz) / (P - Pz));
}

double
Analyzer::calEtaPhiDistance (double Px1, double Py1, double Pz1,
				double Px2, double Py2, double Pz2)
{
  double phi1 = calPhi (Px1,Py1,Pz1);
  double eta1 = calEta (Px1,Py1,Pz1);
  double phi2 = calPhi (Px2,Py2,Pz2);
  double eta2 = calEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}


// ------------ method called for each event  ------------
void
Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



//   Handle<TrackCollection> globaltracks_;
//   iEvent.getByToken(Tracks_global_, globaltracks_);

   Handle<TrackCollection> displacedglobaltracks_;
   iEvent.getByToken(Tracks_Disglobal_, displacedglobaltracks_);

//   Handle<TrackCollection> standalonetracks_;
//   iEvent.getByToken(Tracks_standalone_, standalonetracks_);

   Handle<TrackCollection> displacedstandalonetracks_;
   iEvent.getByToken(Tracks_Disstandalone_, displacedstandalonetracks_);

   Handle<TrackCollection> displacedtracks_;
   iEvent.getByToken(displacedTracksToken_, displacedtracks_);

//   Handle<TrackCollection> displacedoutintracks_;
//   iEvent.getByToken(displacedOutInTracks_, displacedoutintracks_);

//   edm::Handle<reco::MuonCollection> earlydisplacedmuons_;
//   iEvent.getByToken(earlyDisplacedMuonsToken_, earlydisplacedmuons_);

   Handle<GenParticleCollection> genparticle_;
   iEvent.getByToken(gen_Muon_, genparticle_);

//   edm::Handle<std::vector<reco::Muon> > muons_;  
//   iEvent.getByToken(offlineMuonToken_, muons_);

//   edm::Handle<std::vector<reco::Muon> > displacedmuons_;
//   iEvent.getByToken(offlinedisplacedMuonToken_, displacedmuons_);


   /*Handle<int> numoiseeds_;
   iEvent.getByToken(numOIseeds_,numoiseeds_);*/

   //numberOfOISeeds =-10;

  // edm::Handle<std::vector<reco::Track> > tracks;
  // iEvent.getByToken(tracksToken_, tracks);

  /*edm::Handle<std::vector<MuonIdProducer> > muonid;
  iEvent.getByToken(Muonid_, muonid);*/



  pt_glbMuon->clear();

  eta_glbMuon->clear();

  phi_glbMuon->clear();

  dxy_glbMuon->clear();

  dz_glbMuon->clear();

  d0_glbMuon->clear();

  //lxy_glbMuon->clear();


  istrue_DisglbMuon->clear();

  pt_DisglbMuon->clear();

  px_DisglbMuon->clear();

  py_DisglbMuon->clear();

  pz_DisglbMuon->clear();

  eta_DisglbMuon->clear();

  phi_DisglbMuon->clear();

  dxy_DisglbMuon->clear();

  dz_DisglbMuon->clear();

  d0_DisglbMuon->clear();

  lxy_DisglbMuon->clear();


  pt_Sta->clear();

  eta_Sta->clear();

  phi_Sta->clear();

  dxy_Sta->clear();

  dz_Sta->clear();

  d0_Sta->clear();


  istrue_DisSta->clear();

  px_DisSta->clear();

  py_DisSta->clear();

  pz_DisSta->clear();

  pt_DisSta->clear();

  eta_DisSta->clear();

  phi_DisSta->clear();

  dxy_DisSta->clear();

  dz_DisSta->clear();

  d0_DisSta->clear();

  lxy_DisSta->clear();


  istrue_DisTrk->clear();

  px_DisTrk->clear();

  py_DisTrk->clear();

  pz_DisTrk->clear();

  pt_DisTrk->clear();

  eta_DisTrk->clear();

  phi_DisTrk->clear();

  dxy_DisTrk->clear();

  dz_DisTrk->clear();

  d0_DisTrk->clear();

  lxy_DisTrk->clear();


  pt_DisOutIn->clear();

  eta_DisOutIn->clear();

  phi_DisOutIn->clear();

  dxy_DisOutIn->clear();

  dz_DisOutIn->clear();

  d0_DisOutIn->clear();


  pt_DisEar->clear();

  eta_DisEar->clear();

  phi_DisEar->clear();

  dxy_DisEar->clear();

  dz_DisEar->clear();

  d0_DisEar->clear();


  pt_gen->clear();

  eta_gen->clear();

  phi_gen->clear();

  vx_gen->clear();

  vy_gen->clear();

  vz_gen->clear();

  vertex_genx->clear();

  vertex_geny->clear();

  dxy_gen->clear();

  lxy_gen->clear();


  //Muon->clear();

  //displacedMuon->clear();


  pt_Muon->clear();

  eta_Muon->clear();

  phi_Muon->clear();

  mass_Muon->clear();

  dxy_Muon->clear();

  dz_Muon->clear();

  d0_Muon->clear();


  pt_displacedMuon->clear();

  eta_displacedMuon->clear();

  phi_displacedMuon->clear();

  mass_displacedMuon->clear();

  dxy_displacedMuon->clear();

  dz_displacedMuon->clear();

  d0_displacedMuon->clear();


    //Begin loop over GlobalMuonTrac  
/*    for(TrackCollection::const_iterator glbTrack = globaltracks_->begin();
        glbTrack != globaltracks_->end();
        ++glbTrack) {

      pt_glbMuon->push_back(glbTrack->pt());

      eta_glbMuon->push_back(glbTrack->eta());

      phi_glbMuon->push_back(glbTrack->phi());

      dxy_glbMuon->push_back(glbTrack->dxy());

      dz_glbMuon->push_back(glbTrack->dz());

      d0_glbMuon->push_back(glbTrack->d0());

      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }
    //End loop over GlobalMuonTracks
*/

   for(TrackCollection::const_iterator disglbTrack = displacedglobaltracks_->begin();
        disglbTrack != displacedglobaltracks_->end();
        ++disglbTrack) {

      px_DisglbMuon->push_back(disglbTrack->px());
 
      py_DisglbMuon->push_back(disglbTrack->py());
 
      pz_DisglbMuon->push_back(disglbTrack->pz());

      pt_DisglbMuon->push_back(disglbTrack->pt());

      eta_DisglbMuon->push_back(disglbTrack->eta());

      phi_DisglbMuon->push_back(disglbTrack->phi());

      dxy_DisglbMuon->push_back(disglbTrack->dxy());

      dz_DisglbMuon->push_back(disglbTrack->dz());

      d0_DisglbMuon->push_back(disglbTrack->d0());


     double LXY_DisglbMuon;

     calLxy (disglbTrack->vertex().x(), disglbTrack->vertex().y(), 0.0,
             beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
             &LXY_DisglbMuon);

     lxy_DisglbMuon->push_back(LXY_DisglbMuon);

    double size_glb = pt_DisglbMuon->size();

    cout<< "SIZE GLOBAL " <<size_glb<<endl;

      //lxy_DisglbMuon->push_back(disglbTrack->lxy());

      }
    
    //Begin loop over StandaloneMuonTracks    
/*    for(TrackCollection::const_iterator saTrack = standalonetracks_->begin();
        saTrack != standalonetracks_->end();
        ++saTrack) {

      pt_Sta->push_back(saTrack->pt());

      eta_Sta->push_back(saTrack->eta());

      phi_Sta->push_back(saTrack->phi());

      dxy_Sta->push_back(saTrack->dxy());

      dz_Sta->push_back(saTrack->dz());

      d0_Sta->push_back(saTrack->d0());

     }
    // End loop over StandaloneMuonTracks
*/

    for(TrackCollection::const_iterator DissaTrack = displacedstandalonetracks_->begin();
        DissaTrack != displacedstandalonetracks_->end();
        ++DissaTrack) {

      px_DisSta->push_back(DissaTrack->px());

      py_DisSta->push_back(DissaTrack->py());

      pz_DisSta->push_back(DissaTrack->pz());
      
      pt_DisSta->push_back(DissaTrack->pt());

      eta_DisSta->push_back(DissaTrack->eta());

      phi_DisSta->push_back(DissaTrack->phi());

      dxy_DisSta->push_back(DissaTrack->dxy());

      dz_DisSta->push_back(DissaTrack->dz());

      d0_DisSta->push_back(DissaTrack->d0());

     double LXY_DisSta;

     calLxy (DissaTrack->vertex().x(), DissaTrack->vertex().y(), 0.0,
             beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
             &LXY_DisSta);

    lxy_DisSta->push_back(LXY_DisSta);

    double size_sa = pt_DisSta->size();

    cout<< "SIZE SA" <<size_sa<<endl;


     }


    for(TrackCollection::const_iterator DisTrack = displacedtracks_->begin();
        DisTrack != displacedtracks_->end();
        ++DisTrack) {

      px_DisTrk->push_back(DisTrack->px());

      py_DisTrk->push_back(DisTrack->py());

      pz_DisTrk->push_back(DisTrack->pz());

      pt_DisTrk->push_back(DisTrack->pt());

      eta_DisTrk->push_back(DisTrack->eta());

      phi_DisTrk->push_back(DisTrack->phi());

      dxy_DisTrk->push_back(DisTrack->dxy());

      dz_DisTrk->push_back(DisTrack->dz());

      d0_DisTrk->push_back(DisTrack->d0());


     double LXY_DisTrk;

     calLxy (DisTrack->vertex().x(), DisTrack->vertex().y(), 0.0,
             beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
             &LXY_DisTrk);

     lxy_DisTrk->push_back(LXY_DisTrk);

    double size_trk = pt_DisTrk->size();

    cout<< "SIZE TRACK" <<size_trk<<endl;

     }



/*    for(TrackCollection::const_iterator outInTrack = displacedoutintracks_->begin();
        outInTrack != displacedoutintracks_->end();
        ++outInTrack) {

        pt_DisOutIn->push_back(outInTrack->pt());

        eta_DisOutIn->push_back(outInTrack->eta());

        phi_DisOutIn->push_back(outInTrack->phi());

        dxy_DisOutIn->push_back(outInTrack->dxy());

        dz_DisOutIn->push_back(outInTrack->dz());

        d0_DisOutIn->push_back(outInTrack->d0());


          }
*/
     
/*    for (reco::MuonCollection::const_iterator muo = earlydisplacedmuons_->begin(); muo != earlydisplacedmuons_->end(); ++muo){

        pt_DisEar->push_back(muo->pt());

        eta_DisEar->push_back(muo->eta());

        phi_DisEar->push_back(muo->phi());*/

        /*dxy_DisEar->push_back(muo->dxy());

        dz_DisEar->push_back(muo->dz());

        d0_DisEar->push_back(muo->d0());*/



//    }


    //Begin loop over genMuon
    for(GenParticleCollection::const_iterator genTrack = genparticle_->begin();
        genTrack != genparticle_->end();
        ++genTrack) {

     const reco::GenParticle& genPart = (*genTrack);
     if (genPart.status()==1 && std::abs(genPart.pdgId()) == 13 && genPart.isPromptFinalState()){

      px_gen = genTrack->px();

      py_gen = genTrack->py();

      pz_gen = genTrack->pz();

      pt_gen->push_back(genTrack->pt());

      eta_gen->push_back(genTrack->eta());

      phi_gen->push_back(genTrack->phi());

      vx_gen->push_back(genTrack->vx());

      vy_gen->push_back(genTrack->vy());

      vz_gen->push_back(genTrack->vz());

      double vertex_Genx = genTrack->vertex().x();

      vertex_genx->push_back(vertex_Genx);

      double vertex_Geny = genTrack->vertex().y();

      vertex_geny->push_back(vertex_Geny);

      double dxy_Gen = TMath::Abs((genTrack->p4().y()*genTrack->vertex().x() - genTrack->p4().x()*genTrack->vertex().y())/genTrack->pt());

      dxy_gen->push_back(dxy_Gen);

     double LXY_gen;

     calLxy (genTrack->vertex().x(), genTrack->vertex().y(), 0.0,
             beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
             &LXY_gen);

     lxy_gen->push_back(LXY_gen);

    double size_gen = genparticle_->size();

    cout<< "SIZE GEN " <<size_gen<<endl;

}//muon condition
}
//End loop over genMuon



  //std::cout<<"event ---"<<std::endl;

/*  for(std::vector<reco::Muon>::const_iterator itMu = muons_->begin();

    itMu != muons_->end();
    ++itMu) {
*/
    
   //itMu->sharedHits();
   //std::cout<<"Match:"<<itMu->MuonIdProducer::sharedHits()<<std::endl;

/*      pt_Muon->push_back(itMu->pt());

      eta_Muon->push_back(itMu->eta());

      phi_Muon->push_back(itMu->phi());

      mass_Muon->push_back(itMu->mass());*/

      /*dxy_Muon->push_back(itMu->dxy());

      dz_Muon->push_back(itMu->dz());

      d0_Muon->push_back(itMu->d0());*/


   //std::cout<<"muon pt:  "<<itMu->pt()<<std::endl;

   /*if(itMu->isGlobalMuon())std::cout<<"is a global muon"<<std::endl;   
   if(itMu->isDisplacedMuon())std::cout<<"is a Displaced muon"<<std::endl;
   if(itMu->isTrackerMuon())std::cout<<"is a Tracker muon"<<std::endl;
   if(itMu->isPFMuon())std::cout<<"is a PF muon"<<std::endl;*/


/*   }


  for(std::vector<reco::Muon>::const_iterator itdisplacedMu = displacedmuons_->begin();

    itdisplacedMu != displacedmuons_->end();
    ++itdisplacedMu) { 


     const reco::Muon& disMuon = (*itdisplacedMu);
     if (disMuon.isGlobalMuon()==1){ 			 ////suggested by Celia


      pt_displacedMuon->push_back(itdisplacedMu->pt());

      eta_displacedMuon->push_back(itdisplacedMu->eta());

      phi_displacedMuon->push_back(itdisplacedMu->phi());

      mass_displacedMuon->push_back(itdisplacedMu->mass());*/

      /*dxy_displacedMuon->push_back(itdisplacedMu->dxy());

      dz_displacedMuon->push_back(itdisplacedMu->dz());

      d0_displacedMuon->push_back(itdisplacedMu->d0());*/


    //std::cout<<"displacedMuon pt:  "<<itdisplacedMu->pt()<<std::endl;


//   }
//}

/*  for(const int nseeds = numoiseeds_->begin();
       nseeds != numoiseeds_->end();
        ++nseeds) {*/

      //int numberOfOISeeds = *numoiseeds_;

      //numberOfOISeeds.push_back(numberOfOISeeds_);
//}


 
  double deltaEtaPhi;


//Truth Match with Displaced Global 
  for (vector<int>::size_type i = 0; i < pt_DisglbMuon->size(); i++) {

  deltaEtaPhi = calEtaPhiDistance(px_gen, py_gen, pz_gen,
				    px_DisglbMuon->at(i), py_DisglbMuon->at(i), pz_DisglbMuon->at(i));

   cout<< "deltaEtaPhi for Displaced Global  "<<deltaEtaPhi<<endl;

    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istrue_DisglbMuon->push_back(true);
    } else {
      istrue_DisglbMuon->push_back(false);
    }
}

//Truth Match with Displaced SA
  for (vector<int>::size_type i = 0; i < pt_DisSta->size(); i++) {

  deltaEtaPhi = calEtaPhiDistance(px_gen, py_gen, pz_gen,
                                    px_DisSta->at(i), py_DisSta->at(i), pz_DisSta->at(i));

   cout<< "deltaEtaPhi for Displaced SA  "<<deltaEtaPhi<<endl;

    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istrue_DisSta->push_back(true);
    } else {
      istrue_DisSta->push_back(false);
    }
}

//Truth Match with Displaced Track
  for (vector<int>::size_type i = 0; i < pt_DisTrk->size(); i++) {

  deltaEtaPhi = calEtaPhiDistance(px_gen, py_gen, pz_gen,
                                    px_DisTrk->at(i), py_DisTrk->at(i), pz_DisTrk->at(i));

   cout<< "deltaEtaPhi for Displaced Track  "<<deltaEtaPhi<<endl;

    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istrue_DisTrk->push_back(true);
    } else {
      istrue_DisTrk->push_back(false);
    }
}


   tree->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void
Analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
Analyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
