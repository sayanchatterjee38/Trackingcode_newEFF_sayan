#include <memory>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TTree.h>
#include "TParallelCoord.h"
#include "TParallelCoordVar.h"
#include "THnSparse.h"
#include "TAxis.h"


#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
// RecoJets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

//#include "HITrackingStudies/HITrackingStudies/interface/HITrackCorrectionTreeHelper.h"
#include "TrackingCode/HIRun2015Ana/interface/HITrackCorrectionTreeHelper.h"

//added by sayan ~ GenParticleCollection
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


class HITrackCorrectionAnalyzer : public edm::EDAnalyzer 
{
public:
  explicit HITrackCorrectionAnalyzer(const edm::ParameterSet&);
  ~HITrackCorrectionAnalyzer();
  
  static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void initHistos(const edm::Service<TFileService> & fs);
  bool multCuts(const reco::Track & track, const reco::Vertex & vertex);
  bool passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex, float trkMVA);
  bool caloMatched(const reco::Track & track, const edm::Event& iEvent, unsigned it );

  // ----------member data ---------------------------


  std::map<std::string,TTree*> trkTree_;
  std::map<std::string,TH2F*> trkCorr2D_;
  std::map<std::string,TH2F*> ptcent2D_;
  std::map<std::string,TH2F*> etacent2D_;
  std::map<std::string,TH2F*> phicent2D_;
  std::map<std::string,TH3F*> etaphicent3D_;
  std::map<std::string,TH3F*> trkCorr3D_;
  std::map<std::string,THnSparseD*> trkCorr4D_;
  TH3F * momRes_;
  TH1F * vtxZ_;
  /*TH1F * pthat_;
  TH1F * ptreco_;
  TH1F * pttrue_;
  TH1F * ptfake_;
  TH1F * etareco_;
  TH1F * etatrue_;
  TH1F * etafake_;
  TH1F * phireco_;      
  TH1F * phitrue_;
  TH1F * phifake_; 
  
  //sayan ~ 23.10.2021          
  TH1F * ptgen_sayan_TPC;
  TH1F * ptgen_sayan_GPC;
  TH1F * etagen_sayan_TPC;
  TH1F * etagen_sayan_GPC;
  */

  TH1F* pt_reco;
  TH1F* pt_recotosim;
  TH1F* pt_sec;
  TH1F* pt_fake;
  TH1F* pt_sim;
  TH1F* pt_beforedocalomatch;
  TH1F* pt_afterdocalomatch;
  TH1F* pt_eff;
  TH1F* pt_mul;

  TH1F* eta_reco;
  TH1F* eta_recotosim;
  TH1F* eta_sec;
  TH1F* eta_fake;
  TH1F* eta_sim;
  TH1F* eta_beforedocalomatch;
  TH1F* eta_afterdocalomatch;
  TH1F* eta_eff;
  TH1F* eta_mul;

  TH1F* phi_reco;
  TH1F* phi_recotosim;
  TH1F* phi_sec;
  TH1F* phi_fake;
  TH1F* phi_sim;
  TH1F* phi_beforedocalomatch;
  TH1F* phi_afterdocalomatch;
  TH1F* phi_eff;
  TH1F* phi_mul;

  TF1 * vtxWeightFunc_;

  HITrackCorrectionTreeHelper treeHelper_;


  edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;
  edm::EDGetTokenT<TrackingParticleCollection> tpFakSrc_;
  edm::EDGetTokenT<TrackingParticleCollection> tpEffSrc_;
  edm::EDGetTokenT<reco::RecoToSimCollection> associatorMapRTS_;
  edm::EDGetTokenT<reco::SimToRecoCollection> associatorMapSTR_;      
  
  std::vector<double> ptBins_;
  std::vector<double> etaBins_;
  std::vector<double> occBins_;
  std::vector<double> phiBins_;
  
  bool doCaloMatched_;
  double reso_;
  double crossSection_;
  
  std::vector<double> vtxWeightParameters_;
  std::vector<int> algoParameters_;
  bool doVtxReweighting_;

  bool applyVertexZCut_;
  double vertexZMax_;
  
  bool applyTrackCuts_;
  std::string qualityString_;
  double dxyErrMax_;
  double dzErrMax_;
  double ptErrMax_;
  int    nhitsMin_;
  double chi2nMax_;      

  bool doMomRes_;
  
  bool fillNTuples_;

  bool useCentrality_;
  edm::EDGetTokenT<int> centralitySrc_;
  edm::EDGetTokenT<reco::CaloJetCollection> jetSrc_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandSrc_;
  edm::EDGetTokenT<std::vector<float>> mvaSrc_;
  
  edm::EDGetTokenT<reco::GenParticleCollection> trackTagsgen_;

};

HITrackCorrectionAnalyzer::HITrackCorrectionAnalyzer(const edm::ParameterSet& iConfig):
  treeHelper_(),
  vertexSrc_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
  trackSrc_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("trackSrc"))),
  tpFakSrc_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpFakSrc"))),
  tpEffSrc_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpEffSrc"))),
  associatorMapRTS_(consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("associatorMap"))),
  associatorMapSTR_(consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("associatorMap"))),
  ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
  etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
  occBins_(iConfig.getParameter<std::vector<double> >("occBins")),
  phiBins_(iConfig.getParameter<std::vector<double> >("phiBins")),
  doCaloMatched_(iConfig.getParameter<bool>("doCaloMatched")),
  reso_(iConfig.getParameter<double>("reso")),
  crossSection_(iConfig.getParameter<double>("crossSection")),
  vtxWeightParameters_(iConfig.getParameter<std::vector<double> >("vtxWeightParameters")),
  algoParameters_(iConfig.getParameter<std::vector<int> >("algoParameters")),
  doVtxReweighting_(iConfig.getParameter<bool>("doVtxReweighting")),
  applyVertexZCut_(iConfig.getParameter<bool>("applyVertexZCut")),
  vertexZMax_(iConfig.getParameter<double>("vertexZMax")),
  applyTrackCuts_(iConfig.getParameter<bool>("applyTrackCuts")),
  qualityString_(iConfig.getParameter<std::string>("qualityString")),
  dxyErrMax_(iConfig.getParameter<double>("dxyErrMax")),
  dzErrMax_(iConfig.getParameter<double>("dzErrMax")),
  ptErrMax_(iConfig.getParameter<double>("ptErrMax")),
  nhitsMin_(iConfig.getParameter<int>("nhitsMin")),
  chi2nMax_(iConfig.getParameter<double>("chi2nMax")),
  doMomRes_(iConfig.getParameter<bool>("doMomRes")),
  fillNTuples_(iConfig.getParameter<bool>("fillNTuples")),
  useCentrality_(iConfig.getParameter<bool>("useCentrality")),
  centralitySrc_(consumes<int>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
  jetSrc_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  pfCandSrc_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSrc"))),
  mvaSrc_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("mvaSrc"))),
  //added by sayan ~ 23.10.2021
  trackTagsgen_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("tracksgen")))
  
{
  edm::Service<TFileService> fs;
  initHistos(fs);
  
  vtxWeightFunc_ = new TF1("vtxWeight","gaus(0)/gaus(3)",-50.,50.);
  // vtxWeightParameters should have size 6,
  // one really should throw an error if not
  if( (int)vtxWeightParameters_.size() == 6 )
    {
      for( unsigned int i=0;i<vtxWeightParameters_.size(); i++)
	vtxWeightFunc_->SetParameter(i,vtxWeightParameters_[i]);
   }
  
  if( fillNTuples_ )
    {
      trkTree_["rec"] = fs->make<TTree>("recTree","recTree");
      trkTree_["rec"]->Branch("recValues",&treeHelper_.b,treeHelper_.hiTrackLeafString.Data());
      trkTree_["sim"] = fs->make<TTree>("simTree","simTree");
      trkTree_["sim"]->Branch("simValues",&treeHelper_.b,treeHelper_.hiTrackLeafString.Data());
    }
}

HITrackCorrectionAnalyzer::~HITrackCorrectionAnalyzer()
{
  delete vtxWeightFunc_;
}

void
HITrackCorrectionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // obtain collections of simulated particles 
   edm::Handle<TrackingParticleCollection>  TPCollectionHeff, TPCollectionHfake;
   iEvent.getByToken(tpEffSrc_,TPCollectionHeff);
   iEvent.getByToken(tpFakSrc_,TPCollectionHfake);

   // obtain collections of Gen Particles  ~ added by sayan on 23.10.2021
   edm::Handle< reco::GenParticleCollection > tracksgen;
   iEvent.getByToken(trackTagsgen_, tracksgen);

   // obtain association map between tracks and simulated particles
   reco::RecoToSimCollection recSimColl;
   reco::SimToRecoCollection simRecColl;
   edm::Handle<reco::SimToRecoCollection > simtorecoCollectionH;
   edm::Handle<reco::RecoToSimCollection > recotosimCollectionH;
   iEvent.getByToken(associatorMapSTR_,simtorecoCollectionH);
   simRecColl= *(simtorecoCollectionH.product());
   iEvent.getByToken(associatorMapRTS_,recotosimCollectionH);
   recSimColl= *(recotosimCollectionH.product());

   //calo jets
   Handle<reco::CaloJetCollection> JetCollection;
   iEvent.getByToken(jetSrc_, JetCollection);
   if( !JetCollection.isValid() ) return; 
   double leadingJet = 0.;
   for(unsigned irecojet = 0; irecojet < JetCollection->size(); irecojet++ ){
      const reco::CaloJet & JetCand = (*JetCollection)[irecojet];
      if( JetCand.pt() > leadingJet ) leadingJet = JetCand.pt();//finding leading pT jets
   }

   //pthat_->Fill( leadingJet, crossSection_ );

   // obtain reconstructed tracks
   Handle<edm::View<reco::Track> > tcol;
   iEvent.getByToken(trackSrc_, tcol);

   // obtain primary vertices
   Handle<reco::VertexCollection> vertex;
   iEvent.getByToken(vertexSrc_, vertex);
  
   //mva selections //sayan ~ 26.09.2021                                                             
   edm::Handle<std::vector<float>> mvaoutput;
   iEvent.getByToken(mvaSrc_, mvaoutput);

   // sort the vertcies by number of tracks in descending order
   reco::VertexCollection vsorted = *vertex;
   std::sort( vsorted.begin(), vsorted.end(), HITrackCorrectionAnalyzer::vtxSort );

   // skip events with no PV, this should not happen
   if( vsorted.size() == 0) return;

   // skip events failing vertex cut
   if( applyVertexZCut_)
   {
     if( fabs(vsorted[0].z()) > vertexZMax_ ) return;
   }

   // determine vertex reweighting factor
   double w = 1.0;
   w = w * crossSection_;
   if ( doVtxReweighting_ )
     w *= vtxWeightFunc_->Eval(vsorted[0].z());

   vtxZ_->Fill(vsorted[0].z(),w);

   // determine event multipliticy
   int multiplicity =0;
   for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){
     edm::RefToBase<reco::Track> track(tcol, i);
     reco::Track* tr=const_cast<reco::Track*>(track.get());
     if( multCuts(*tr, vsorted[0]) )
       multiplicity++;
   }

   // determine centrality if set
   // note if there is no centrality information multiplicity 
   // will be used in place of the centrality
   double cbin = multiplicity;
   double occ = multiplicity;
   if( useCentrality_ )
   {
     edm::Handle<int> centralityBin;
     iEvent.getByToken(centralitySrc_, centralityBin);
     cbin = *centralityBin;
     occ = cbin;
   } 

   // ---------------------
   // loop through reco tracks to fill fake, reco, and secondary histograms
   // ---------------------

   int it1 = 0;

   for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){ 
     edm::RefToBase<reco::Track> track(tcol, i);
     reco::Track* tr=const_cast<reco::Track*>(track.get());

     float trkMVA1 = (*mvaoutput)[it1];
     it1++;

     // skip tracks that fail cuts, using vertex with most tracks as PV       
     if( ! passesTrackCuts(*tr, vsorted[0], trkMVA1) ) continue;
     if( ! caloMatched(*tr, iEvent, i) ) continue;
     
     if(tr->charge()==0)continue; //for all charges 
     //if(tr->charge()<=0)continue; //for positive
     //if(tr->charge()>=0)continue; //for negative

     Double_t x4d_recsecfak[4] = {tr->eta(), tr->pt(), tr->phi(), occ};

     ptcent2D_["hrecpT"]->Fill( tr->pt(), occ, w);
     etacent2D_["hreceta"]->Fill( tr->eta(), occ, w);
     phicent2D_["hrecphi"]->Fill( tr->phi(), occ, w);

     
     trkCorr2D_["hrec"]->Fill(tr->eta(), tr->pt(), w);
     trkCorr3D_["hrec3D"]->Fill(tr->eta(), tr->pt(), occ, w);

     etaphicent3D_["hrecetaphi3D"]->Fill(tr->eta(),tr->phi(), occ, w);
     
     trkCorr4D_["hrec4D"]->Fill(x4d_recsecfak, w);

     pt_reco->Fill(tr->pt());
     eta_reco->Fill(tr->eta());
     phi_reco->Fill(tr->phi());
     
     // look for match to simulated particle, use first match if it exists
     std::vector<std::pair<TrackingParticleRef, double> > tp;
     const TrackingParticle *mtp=0;
     if(recSimColl.find(track) != recSimColl.end())
     {
       tp = recSimColl[track];
       mtp = tp.begin()->first.get();  
       
       //added to look at the distribution of reco which is matching to simulated
       pt_recotosim->Fill(tr->pt());
       eta_recotosim->Fill(tr->eta());
       phi_recotosim->Fill(tr->phi());
       
       if( fillNTuples_) treeHelper_.Set(*mtp, *tr, vsorted[0], tp.size(), cbin); 
       if( mtp->status() < 0 ) 
       {

	 ptcent2D_["hsecpT"]->Fill( tr->pt(), occ, w);
	 etacent2D_["hseceta"]->Fill( tr->eta(), occ, w);
	 phicent2D_["hsecphi"]->Fill( tr->phi(), occ, w);

	 
         trkCorr2D_["hsec"]->Fill(tr->eta(), tr->pt(), w);     
         trkCorr3D_["hsec3D"]->Fill(tr->eta(), tr->pt(), occ, w);    
         pt_sec->Fill(tr->pt());
     	 eta_sec->Fill(tr->eta());
     	 phi_sec->Fill(tr->phi()); 

	 etaphicent3D_["hsecetaphi3D"]->Fill(tr->eta(),tr->phi(), occ, w);
	 
	 trkCorr4D_["hsec4D"]->Fill(x4d_recsecfak, w);

       }
     }
     else
     {

       ptcent2D_["hfakpT"]->Fill( tr->pt(), occ, w);
       etacent2D_["hfaketa"]->Fill( tr->eta(), occ, w);
       phicent2D_["hfakphi"]->Fill( tr->phi(), occ, w);
       
       if( fillNTuples_) treeHelper_.Set(*tr, vsorted[0], cbin); 
       trkCorr2D_["hfak"]->Fill(tr->eta(), tr->pt(), w);
       trkCorr3D_["hfak3D"]->Fill(tr->eta(), tr->pt(), occ, w);
       
       pt_fake->Fill(tr->pt());
       eta_fake->Fill(tr->eta());
       phi_fake->Fill(tr->phi()); 

       etaphicent3D_["hfaketaphi3D"]->Fill(tr->eta(),tr->phi(), occ, w);
       
       trkCorr4D_["hfak4D"]->Fill(x4d_recsecfak, w);


     }
     if( fillNTuples_) trkTree_["rec"]->Fill(); 
   }

   // ---------------------
   // loop through sim particles to fill matched, multiple,  and sim histograms 
   // ---------------------

   for(TrackingParticleCollection::size_type i=0; i<TPCollectionHeff->size(); i++) 
   {      
     TrackingParticleRef tpr(TPCollectionHeff, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
         
     if(tp->status() < 0 || tp->charge()==0 ) continue; //only charged primaries
     //if(tp->status() < 0 || tp->charge()<=0 ) continue; //only positive charged primaries
     //if(tp->status() < 0 || tp->charge()>=0 ) continue; //only negative charged primaries
     
     if (tp->pt()< 0.5) continue;
     if (fabs(tp->eta()) > 2.4) continue;

     pt_sim->Fill(tp->pt());
     eta_sim->Fill(tp->eta());
     phi_sim->Fill(tp->phi());

     ptcent2D_["hsimpT"]->Fill( tp->pt(), occ, w);
     etacent2D_["hsimeta"]->Fill( tp->eta(), occ, w);
     phicent2D_["hsimphi"]->Fill( tp->phi(), occ, w);

     
     trkCorr2D_["hsim"]->Fill(tp->eta(),tp->pt(), w);
     trkCorr3D_["hsim3D"]->Fill(tp->eta(),tp->pt(), occ, w);

     etaphicent3D_["hsimetaphi3D"]->Fill(tp->eta(),tp->phi(), occ, w);
     
     Double_t x4d_sim[4] = {tp->eta(), tp->pt(), tp->phi(), occ};
     trkCorr4D_["hsim4D"]->Fill(x4d_sim, w);

     // find number of matched reco tracks that pass cuts
     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
     size_t nrec=0;
     if(simRecColl.find(tpr) != simRecColl.end())
     {
       rt = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) simRecColl[tpr];
       std::vector<std::pair<edm::RefToBase<reco::Track>, double> >::const_iterator rtit;
       
       int it2 = 0; //bysayan ~ 26.09.2021

       for (rtit = rt.begin(); rtit != rt.end(); ++rtit)
       {
         const reco::Track* tmtr = rtit->first.get();

	 float trkMVA2 = (*mvaoutput)[it2];
         it2++;

         if( ! passesTrackCuts(*tmtr, vsorted[0], trkMVA2) ) continue;
         //cross check
	 pt_beforedocalomatch -> Fill(tmtr->pt());
	 eta_beforedocalomatch -> Fill(tmtr->eta());
	 phi_beforedocalomatch -> Fill(tmtr->phi());
	 
	 unsigned index = -1;
         if( doCaloMatched_ ){ 
          for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){ 
             edm::RefToBase<reco::Track> track(tcol, i);
             reco::Track* tr=const_cast<reco::Track*>(track.get());
             index++;
             if( tmtr->pt() == tr->pt() && tmtr->eta() == tr->eta() && tmtr->phi() == tr->phi() && tmtr->numberOfValidHits() == tr->numberOfValidHits() ) break;//simple match to find the corresponding index number (i-th track) in the track collection
          }
         if( ! caloMatched(*tmtr, iEvent, index) ) continue;
         }  
         nrec++;
         if( doMomRes_ ) momRes_->Fill( tp->eta(), tp->pt(), tmtr->pt(), w);
       }
     }
     //cross check
     pt_afterdocalomatch -> Fill(tp->pt());
     eta_afterdocalomatch -> Fill(tp->eta());
     phi_afterdocalomatch -> Fill(tp->phi());


     if( nrec>0 && fillNTuples_ ) treeHelper_.Set(*tp, *(rt.begin()->first.get()), vsorted[0], rt.size(), cbin);
     if( nrec==0 && fillNTuples_ ) treeHelper_.Set(*tp, cbin);

     if(nrec>0) ptcent2D_["heffpT"]->Fill( tp->pt(), occ, w);
     if(nrec>0) etacent2D_["heffeta"]->Fill( tp->eta(), occ, w);
     if(nrec>0) phicent2D_["heffphi"]->Fill( tp->phi(), occ, w);
     if(nrec>1) ptcent2D_["hmulpT"]->Fill( tp->pt(), occ, w);
     if(nrec>1) etacent2D_["hmuleta"]->Fill( tp->eta(), occ, w);
     if(nrec>1) phicent2D_["hmulphi"]->Fill( tp->phi(), occ, w);
     
     if(nrec>0) trkCorr2D_["heff"]->Fill(tp->eta(),tp->pt(), w);
     if(nrec>0) trkCorr3D_["heff3D"]->Fill(tp->eta(),tp->pt(), occ, w);
     if(nrec>1) trkCorr2D_["hmul"]->Fill(tp->eta(),tp->pt(), w);
     if(nrec>1) trkCorr3D_["hmul3D"]->Fill(tp->eta(),tp->pt(), occ, w);
     if( fillNTuples_) trkTree_["sim"]->Fill(); 

     if(nrec>0) etaphicent3D_["heffetaphi3D"]->Fill(tp->eta(),tp->phi(), occ, w);
     if(nrec>1) etaphicent3D_["hmuletaphi3D"]->Fill(tp->eta(),tp->phi(), occ, w);

     
     Double_t x4d_effmul[4] = {tp->eta(), tp->pt(), tp->phi(), occ};
     if(nrec>0) trkCorr4D_["heff4D"]->Fill(x4d_effmul, w);
     if(nrec>1) trkCorr4D_["hmul4D"]->Fill(x4d_effmul, w);

     //cross check
     if(nrec>0) pt_eff->Fill(tp->pt());
     if(nrec>1) pt_mul->Fill(tp->pt());

     if(nrec>0) eta_eff->Fill(tp->eta());
     if(nrec>1) eta_mul->Fill(tp->eta());

     if(nrec>0) phi_eff->Fill(tp->phi());
     if(nrec>1) phi_mul->Fill(tp->phi());
}


   //added for gen particles
   for( reco::GenParticleCollection::const_iterator itTrk = tracksgen->begin(); itTrk != tracksgen->end(); ++itTrk )
     {
       // Get eta, pt, phi and charge of the track                                                                                                                                                        
       //  if( itTrk->status() != 1 || itTrk->charge() == 0  ) continue; //only charged primaries
       //if( itTrk->status() != 1 || itTrk->charge() <= 0  ) continue; //only positive charged primaries
       if( itTrk->status() != 1 || itTrk->charge() >= 0  ) continue; //only negative charged primaries
              
       if ( itTrk->pt() < 0.5 ) continue;
       if ( fabs(itTrk->eta()) > 2.4) continue;

       //   ptgen_sayan_GPC->Fill(itTrk->pt());
       //etagen_sayan_GPC->Fill(itTrk->eta());
     }


}

bool
HITrackCorrectionAnalyzer::multCuts(const reco::Track & track, const reco::Vertex & vertex)
{

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   vtxPoint=vertex.position();
   vzErr=vertex.zError();
   vxErr=vertex.xError();
   vyErr=vertex.yError();

   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);

   if(track.quality(reco::TrackBase::qualityByName(qualityString_)) != 1) return false;
   if(fabs(dxy/dxysigma) >= 3.0) return false;
   if(fabs(dz/dzsigma) >= 3.0) return false;
   if(fabs(track.ptError()) / track.pt() >= 0.1) return false;
   if( track.pt() < 0.5 || fabs(track.eta()) > 2.4 ) return false;

   return true;

}

bool
HITrackCorrectionAnalyzer::passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex, float trkMVA)
{
   if ( ! applyTrackCuts_ ) return true;

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   vtxPoint=vertex.position();
   vzErr=vertex.zError();
   vxErr=vertex.xError();
   vyErr=vertex.yError();

   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
   
   double chi2n = track.normalizedChi2();
   double nlayers = track.hitPattern().trackerLayersWithMeasurement();
   chi2n = chi2n/nlayers;
   int nhits = track.numberOfValidHits();
   int algo  = track.algo(); 

   if( !track.quality(reco::TrackBase::highPurity) ) return false;
   if(fabs(dxy/dxysigma) >= dxyErrMax_) return false;
   if(fabs(dz/dzsigma) >= dzErrMax_) return false;
   if(fabs(track.ptError()) / track.pt() >= ptErrMax_) return false;
   if(nhits < nhitsMin_ ) return false;
   if(chi2n >= chi2nMax_ ) return false;  
   if(algo == 6 && trkMVA < 0.98) return false;
   
   if ( track.pt() < 0.5 ) return false;
   if ( fabs(track.eta()) > 2.4) return false;
   
   return true;
}

bool 
HITrackCorrectionAnalyzer::caloMatched( const reco::Track & track, const edm::Event& iEvent, unsigned it )
{
  if( ! doCaloMatched_ ) return true;
  
  // obtain pf candidates
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandSrc_, pfCandidates);
  if( !pfCandidates.isValid() ) return false;

  double ecalEnergy = 0.;
  double hcalEnergy = 0.;

  for( unsigned ic = 0; ic < pfCandidates->size(); ic++ ) {//calo matching loops

      const reco::PFCandidate& cand = (*pfCandidates)[ic];

      int type = cand.particleId();

      // only charged hadrons and leptons can be asscociated with a track
      if(!(type == reco::PFCandidate::h ||     //type1
      type == reco::PFCandidate::e ||     //type2
      type == reco::PFCandidate::mu      //type3
      )) continue;

      reco::TrackRef trackRef = cand.trackRef();
      if( it == trackRef.key() ) {
        // cand_index = ic;
        ecalEnergy = cand.ecalEnergy();
        hcalEnergy = cand.hcalEnergy();              
        break;
      } 
  }

  //if((track.pt()-reso_*track.ptError())*TMath::CosH( track.eta() )>15 && (track.pt()-reso_*track.ptError())*TMath::CosH( track.eta() ) > hcalEnergy+ecalEnergy ) return false;
  if( track.pt() < 20 || ( (hcalEnergy+ecalEnergy)/( track.pt()*TMath::CosH(track.eta() ) ) > reso_ && (hcalEnergy+ecalEnergy)/(TMath::CosH(track.eta())) > (track.pt() - 80.0) )  ) return true;
  else {
    return false;
  }
}


void
HITrackCorrectionAnalyzer::initHistos(const edm::Service<TFileService> & fs)
{

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  std::vector<std::string> hNames2D = { "hsim", "hrec", "hmul", "hfak", "heff", "hsec" };

  for( auto name : hNames2D )
  {
     trkCorr2D_[name] = fs->make<TH2F>(name.c_str(),";#eta;p_{T}",
				       etaBins_.size()-1, &etaBins_[0],
				       ptBins_.size()-1, &ptBins_[0]);
  }
  

  //pt cent
  std::vector<std::string> hNamespT2D = { "hsimpT", "hrecpT", "hmulpT", "hfakpT", "heffpT", "hsecpT" };
  for( auto name : hNamespT2D )
    {
      ptcent2D_[name] = fs->make<TH2F>(name.c_str(),";p_{T};occ",
				       ptBins_.size()-1, &ptBins_[0],
				       occBins_.size()-1, &occBins_[0]);
    }
  
  //eta cent
  std::vector<std::string> hNameseta2D = { "hsimeta", "hreceta", "hmuleta", "hfaketa", "heffeta", "hseceta" };
  for( auto name : hNameseta2D )
    {
      etacent2D_[name] = fs->make<TH2F>(name.c_str(),";#eta;occ",
					etaBins_.size()-1, &etaBins_[0],
					occBins_.size()-1, &occBins_[0]);
    }
  
  //phi cent
  std::vector<std::string> hNamesphi2D = { "hsimphi", "hrecphi", "hmulphi", "hfakphi", "heffphi", "hsecphi" };
  for( auto name : hNamesphi2D )
    {
      phicent2D_[name] = fs->make<TH2F>(name.c_str(),";#phi;occ",
					phiBins_.size()-1, &phiBins_[0],
					occBins_.size()-1, &occBins_[0]);
    }

  //eta pT cent
  std::vector<std::string> hNames3D = { "hsim3D", "hrec3D", "hmul3D", "hfak3D", "heff3D", "hsec3D" };

  for( auto name : hNames3D )
  {
     trkCorr3D_[name] = fs->make<TH3F>(name.c_str(),";#eta;p_{T};occ",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]);
  }

  //eta phi cent
  std::vector<std::string> hNamesetaphi3D = { "hsimetaphi3D", "hrecetaphi3D", "hmuletaphi3D", "hfaketaphi3D", "heffetaphi3D", "hsecetaphi3D" };

  for( auto name : hNamesetaphi3D )
  {
     etaphicent3D_[name] = fs->make<TH3F>(name.c_str(),";#eta;#phi;occ",
                           etaBins_.size()-1, &etaBins_[0],
                           phiBins_.size()-1, &phiBins_[0],
                           occBins_.size()-1, &occBins_[0]);
  }


  //4d
  std::vector<std::string> hNames4D = { "hsim4D", "hrec4D", "hmul4D", "hfak4D", "heff4D", "hsec4D" };
  
  const Int_t ndims =4;    //eta, pT, phi, centBin
  Int_t bins[ndims] = {18, 50, 63, 8};
  Double_t xmin[ndims]={0., 0., -3.15, 0.};
  Double_t xmax[ndims] = {10., 10., 3.15, 10.};
  

  for( auto name : hNames4D )
    {
      trkCorr4D_[name] = fs->make<THnSparseD>(name.c_str(),";#eta;p_{T};#phi;occ",  ndims, bins, xmin, xmax);
                           
      const int nVarBins0 = 18;
      Double_t varBins0[nVarBins0+1] = { -2.4, -2.0, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4 };
      trkCorr4D_[name]->GetAxis(0)->Set(nVarBins0, varBins0);
    
      const int nVarBins1 = 50;
      Double_t varBins1[nVarBins1+1] = {  0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0, 180.0, 300.0, 500.0 };
      trkCorr4D_[name]->GetAxis(1)->Set(nVarBins1, varBins1);
      
      //const int nVarBins2 = 22;
      //Double_t varBins2[nVarBins2+1] = { -3.3, -3.0, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3 };
      //trkCorr4D_[name]->GetAxis(2)->Set(nVarBins2, varBins2);
      
      //const int nVarBins2 = 27;
      //Double_t varBins2[nVarBins2+1] = { -3.3, -3.0, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3 };
      //trkCorr4D_[name]->GetAxis(2)->Set(nVarBins2, varBins2);
      
      const int nVarBins3 = 8;
      Double_t varBins3[nVarBins3+1] = { 0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 140.0, 200.0 };
      trkCorr4D_[name]->GetAxis(3)->Set(nVarBins3, varBins3);
      
      trkCorr4D_[name]->Sumw2();
      
    }
  
  

  vtxZ_ = fs->make<TH1F>("vtxZ","Vertex z position",100,-30,30);
  
  /*pthat_ = fs->make<TH1F>("pthat", "p_{T}(GeV)", 8000,0,800);
  ptreco_ = fs->make<TH1F>("ptreco_", "p_{T}(GeV)", 1000,0,100);
  pttrue_ = fs->make<TH1F>("pttrue_", "p_{T}(GeV)", 1000,0,100);
  ptfake_ = fs->make<TH1F>("ptfake_", "p_{T}(GeV)", 1000,0,100);
  etareco_ = fs->make<TH1F>("etareco_", "#eta", 200,-2.5,2.5);
  etatrue_ = fs->make<TH1F>("etatrue_", "#eta", 200,-2.5,2.5);
  etafake_ = fs->make<TH1F>("etafake_", "#eta", 200,-2.5,2.5);
  phireco_ = fs->make<TH1F>("phireco_", "#phi", 200,-3.2,3.2);
  phitrue_ = fs->make<TH1F>("phitrue_", "#phi", 200,-3.2,3.2);
  phifake_ = fs->make<TH1F>("phifake_", "#phi", 200,-3.2,3.2); 

  */
  
  pt_reco = fs->make<TH1F>("pt_reco", "p_{T}(GeV)", 1000,0,100);
  pt_recotosim = fs->make<TH1F>("pt_recotosim", "p_{T}(GeV)", 1000,0,100);
  pt_sec = fs->make<TH1F>("pt_sec", "p_{T}(GeV)", 1000,0,100);
  pt_fake = fs->make<TH1F>("pt_fake", "p_{T}(GeV)", 1000,0,100);
  pt_sim = fs->make<TH1F>("pt_sim", "p_{T}(GeV)", 1000,0,100);
  pt_beforedocalomatch = fs->make<TH1F>("pt_beforedocalomatch", "p_{T}(GeV)", 1000,0,100);
  pt_afterdocalomatch = fs->make<TH1F>("pt_afterdocalomatch", "p_{T}(GeV)", 1000,0,100);
  pt_eff = fs->make<TH1F>("pt_eff", "p_{T}(GeV)", 1000,0,100);
  pt_mul = fs->make<TH1F>("pt_mul", "p_{T}(GeV)", 1000,0,100);

  eta_reco = fs->make<TH1F>("eta_reco", "#eta", 48, -2.4, 2.4);
  eta_recotosim = fs->make<TH1F>("eta_recotosim", "#eta", 48, -2.4, 2.4);
  eta_sec = fs->make<TH1F>("eta_sec", "#eta", 48, -2.4, 2.4);
  eta_fake = fs->make<TH1F>("eta_fake", "#eta", 48, -2.4, 2.4);
  eta_sim = fs->make<TH1F>("eta_sim", "#eta", 48, -2.4, 2.4);
  eta_beforedocalomatch = fs->make<TH1F>("eta_beforedocalomatch", "#eta", 48, -2.4, 2.4);
  eta_afterdocalomatch = fs->make<TH1F>("eta_afterdocalomatch", "#eta", 48, -2.4, 2.4);
  eta_eff = fs->make<TH1F>("eta_eff", "#eta", 48, -2.4, 2.4);
  eta_mul = fs->make<TH1F>("eta_mul", "#eta", 48,-2.4, 2.4);

  phi_reco = fs->make<TH1F>("phi_reco", "#phi", 63, -3.15, 3.15);
  phi_recotosim = fs->make<TH1F>("phi_recotosim", "#phi", 63, -3.15, 3.15);
  phi_sec = fs->make<TH1F>("phi_sec", "#phi", 63, -3.15, 3.15);
  phi_fake = fs->make<TH1F>("phi_fake", "#phi", 63, -3.15, 3.15);
  phi_sim = fs->make<TH1F>("phi_sim", "#phi", 63, -3.15, 3.15);
  phi_beforedocalomatch = fs->make<TH1F>("phi_beforedocalomatch", "#phi", 63, -3.15, 3.15);
  phi_afterdocalomatch = fs->make<TH1F>("phi_afterdocalomatch", "#phi", 63, -3.15, 3.15);
  phi_eff = fs->make<TH1F>("phi_eff", "#phi", 63, -3.15, 3.15);
  phi_mul = fs->make<TH1F>("phi_mul", "#phi", 63, -3.15, 3.15);
  
  

  /*
  //23.10.2021 ~ sayan
  ptgen_sayan_TPC = fs->make<TH1F>("ptgen_sayan_TPC", "p_{T}(GeV)", 75,0.5,8.0);
  ptgen_sayan_GPC = fs->make<TH1F>("ptgen_sayan_GPC", "p_{T}(GeV)", 75,0.5,8.0);
  etagen_sayan_TPC = fs->make<TH1F>("etagen_sayan_TPC", "#eta", etaBins_.size()-1, &etaBins_[0]);
  etagen_sayan_GPC = fs->make<TH1F>("etagen_sayan_GPC", "#eta", etaBins_.size()-1, &etaBins_[0]);
*/
  std::vector<double> ptBinsFine;
  for( unsigned int bin = 0; bin<ptBins_.size()-1; bin++)
  {
    double bStart = ptBins_[bin];
    double bWid = ptBins_[bin+1] - ptBins_[bin];
    for( int i=0;i<5;i++)
      ptBinsFine.push_back( bStart + (double)i * bWid / 5. );
  }
  ptBinsFine.push_back(ptBins_[ptBins_.size()-1]);

  momRes_ = fs->make<TH3F>("momRes","momentum resolution sim vs reco",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBinsFine.size()-1, &ptBinsFine[0],
                           ptBinsFine.size()-1, &ptBinsFine[0]);
}

bool
HITrackCorrectionAnalyzer::vtxSort( const reco::Vertex &  a, const reco::Vertex & b )
{
  if( a.tracksSize() != b.tracksSize() )
    return  a.tracksSize() > b.tracksSize() ? true : false ;
  else
    return  a.chi2() < b.chi2() ? true : false ;  
}

void
HITrackCorrectionAnalyzer::beginJob()
{
}

void
HITrackCorrectionAnalyzer::endJob()
{
}

DEFINE_FWK_MODULE(HITrackCorrectionAnalyzer);
