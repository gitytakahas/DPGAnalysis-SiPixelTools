// File: PixClusterAna.cc
// Description: T0 test the pixel clusters. 
// Author: Danek Kotlinski 
// Creation Date:  Initial version. 3/06
// Modify to work with CMSSW354, 11/03/10 d.k.
// Modify to work with CMSSW620, SLC6, CMSSW700, 10/10/13 d.k.
// Change to ByToken (clusters only) so other data access does not work for the moment.
// Make for 4 layers.
//--------------------------------------------
#include <memory>
#include <string>
#include <iostream>

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

//#include "DataFormats/Common/interface/EDProduct.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
//#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

// For L1
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"

// For HLT
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"

// for resyncs 
#include "DataFormats/Scalers/interface/Level1TriggerScalers.h"

// For luminisoty
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrector.h"

// For PVs
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


// To use root histos
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// For ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

using namespace std;

#define PV


//=======================================================================

// decode a simplified ROC address
int rocId(int col, int row) {
  int rocRow = row/80;
  int rocCol = col/52;
  int rocId = rocCol + rocRow*8;
  return rocId;
}


//==========================================================================================

class PixClusterAna : public edm::EDAnalyzer {
 public:
  
  explicit PixClusterAna(const edm::ParameterSet& conf);  
  virtual ~PixClusterAna();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  virtual void beginRun(const edm::EventSetup& iSetup);
  virtual void beginJob();
  virtual void endJob();
  
 private:
  edm::ParameterSet conf_;
  edm::InputTag src_;
  //const static bool PRINT = false;
  bool PRINT;
  int select1, select2;
  int countEvents, countAllEvents;
  double sumClusters, sumPixels, countLumi;
  //float rocHits[3][10]; // layer, ids, store hits in 30 rocs
  //float moduleHits[3][10]; // layer, ids, store hits in 30 modules
  bool phase1_;
  bool eventFlag[40];

  // Needed for the ByToken method
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > myClus;
  edm::EDGetTokenT<LumiSummary> LumiToken;
  edm::EDGetTokenT<edm::ConditionsInLumiBlock> ConditionsInLumiBlockToken;
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> L1TrigReadoutToken;
  edm::EDGetTokenT<edm::TriggerResults> TrigResultsToken;
  edm::EDGetTokenT<reco::VertexCollection> VertexCollectionToken;


  //TFile* hFile;

  TFile* myfile;
  
  TH1F* h_total;

  TTree* event_tree;
  Int_t b_bxid;
  Int_t b_run;
  Int_t b_event;
  Int_t b_lumi;
  Int_t b_orbit;
  Int_t b_npv;
  Int_t b_nmodule;

  Int_t b_countEvents;
  Int_t b_numberOfPixels;
  Int_t b_numberOfClusters;
  Int_t b_maxClusPerDet;
  Int_t b_maxPixPerClu;
  Int_t b_numOfClustersPerLay1;
  Int_t b_numOfClustersPerLay2;
  Int_t b_numOfClustersPerLay3;
  Int_t b_numOfClustersPerLay4;
  Int_t b_numOfClustersPerDisk1;
  Int_t b_numOfClustersPerDisk2;
  Int_t b_numOfClustersPerDisk3;
  Int_t b_numOfClustersPerDisk4;
  Int_t b_numOfPixPerLay1;
  Int_t b_numOfPixPerLay2;
  Int_t b_numOfPixPerLay3;
  Int_t b_numOfPixPerLay4;
  Int_t b_numOfPixPerDisk1;
  Int_t b_numOfPixPerDisk2;
  Int_t b_numOfPixPerDisk3;
  Int_t b_numOfPixPerDisk4;
  Int_t b_numberOfDetUnits1;
  Int_t b_numberOfDetUnits2;
  Int_t b_numberOfDetUnits3;
  Int_t b_numberOfDetUnits4;
  Int_t b_numOfPixBarrel;
  Int_t b_numOfPixForward;
  Int_t b_numOfClusterBarrel;
  Int_t b_numOfClusterForward;
  Int_t b_pixb;
  Int_t b_numberOfNoneEdgePixels;
  Int_t b_maxPixPerDet;


  TTree* pixel_tree;
  Int_t p_bxid;
  Int_t p_run;
  Int_t p_event;
  Int_t p_lumi;
  Int_t p_orbit;
  Int_t p_npv;
  Int_t p_nmodule;
  Int_t p_layer;
  Int_t p_module;
  Int_t p_ladder;
  Int_t p_sector;
  Int_t p_subid;
  Int_t p_adc;
  Int_t p_pixx;
  Int_t p_pixy;
  Float_t p_tmp1;
  Float_t p_tmp2;
  Int_t p_roc;
  Int_t p_link;
  Int_t p_rocInCol;
  Float_t p_ch;
  Int_t p_size;
  Int_t p_sizex;
  Int_t p_sizey;
  Float_t p_x;
  Float_t p_y;
  Int_t p_minPixelRow;
  Int_t p_maxPixelRow;
  Int_t p_minPixelCol;
  Int_t p_maxPixelCol;
  Float_t p_lx;
  Float_t p_ly;
  Float_t p_zPos;
  Float_t p_rPos;

  TTree* cluster_tree;
  Int_t c_module;
  Int_t c_layer;
  Int_t c_ladder;
  Int_t c_size;
  Int_t c_sizex;
  Int_t c_sizey;
  Int_t c_subid;
  Int_t c_ch;
  Float_t c_x;
  Float_t c_y;
  //  Float_t c_z;
  Float_t c_zPos;
  Float_t c_rPos;
  Float_t c_phiPos;
  Int_t c_sector;


  TH2F *hmap[4][20];

  TH1D *hclusPerDet1,*hclusPerDet2,*hclusPerDet3,*hclusPerDet4;
  TH1D *hpixPerDet1,*hpixPerDet2,*hpixPerDet3,*hpixPerDet4;
  //TH1D *hpixPerLink1,*hpixPerLink2,*hpixPerLink3,*hpixPerLink4;

  TH1D *hlt1, *hlt2, *hlt3;

  TH1D *hpixPerDet11,*hpixPerDet12,*hpixPerDet13,*hpixPerDet14;
  TH1D *hpixPerDet21,*hpixPerDet22,*hpixPerDet23,*hpixPerDet24;
  TH1D *hpixPerDet31,*hpixPerDet32,*hpixPerDet33,*hpixPerDet34;
  TH1D *hpixPerDet41,*hpixPerDet42,*hpixPerDet43,*hpixPerDet44;

};

/////////////////////////////////////////////////////////////////
// Contructor, empty.
PixClusterAna::PixClusterAna(edm::ParameterSet const& conf) 
  : conf_(conf), src_(conf.getParameter<edm::InputTag>( "src" )) { 
  PRINT = conf.getUntrackedParameter<bool>("Verbosity",false);
  select1 = conf.getUntrackedParameter<int>("Select1",0);
  select2 = conf.getUntrackedParameter<int>("Select2",0);
  //src_ =  conf.getParameter<edm::InputTag>( "src" );
  if(PRINT) cout<<" Construct "<<endl;

  // For the ByToken method
  myClus = consumes<edmNew::DetSetVector<SiPixelCluster> >(conf.getParameter<edm::InputTag>( "src" ));
  LumiToken                  = consumes <LumiSummary>(edm::InputTag("lumiProducer"));
  ConditionsInLumiBlockToken = consumes <edm::ConditionsInLumiBlock> (edm::InputTag("conditionsInEdm"));
  L1TrigReadoutToken         = consumes <L1GlobalTriggerReadoutRecord>(edm::InputTag("gtDigis"));
  TrigResultsToken           = consumes <edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  VertexCollectionToken      = consumes <reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  phase1_ = conf.getUntrackedParameter<bool>("phase1",false);


}
// Virtual destructor needed.
PixClusterAna::~PixClusterAna() { }  

// ------------ method called at the begining   ------------
void PixClusterAna::beginRun(const edm::EventSetup& iSetup) {
  cout << "beginRun -  PixelClusterTest " <<endl;
}

// ------------ method called at the begining   ------------
void PixClusterAna::beginJob() {
  cout << "Initialize PixClusterAna " <<endl;

  //  edm::Service<TFileService> fs;

  myfile = new TFile("Myroot.root","recreate");
  event_tree = new TTree("event_tree", "event_tree");
  event_tree->Branch("bxid", &b_bxid, "bxid/I");
  event_tree->Branch("run", &b_run, "run/I");
  event_tree->Branch("event", &b_event, "event/I");
  event_tree->Branch("lumi", &b_lumi, "lumi/I");
  event_tree->Branch("orbit", &b_orbit, "orbit/I");
  event_tree->Branch("npv", &b_npv, "npv/I");
  event_tree->Branch("nmodule", &b_nmodule, "nmodule/I");

  event_tree->Branch("countEvents", &b_countEvents, "countEvents/I");
  event_tree->Branch("numberOfPixels", &b_numberOfPixels, "numberOfPixels/I");
  event_tree->Branch("numberOfClusters", &b_numberOfClusters, "numberOfClusters/I");
  event_tree->Branch("maxClusPerDet", &b_maxClusPerDet, "maxClusPerDet/I");
  event_tree->Branch("maxPixPerClu", &b_maxPixPerClu, "maxPixPerClu/I");

  event_tree->Branch("numOfClustersPerLay1", &b_numOfClustersPerLay1, "numOfClustersPerLay1/I");
  event_tree->Branch("numOfClustersPerLay2", &b_numOfClustersPerLay2, "numOfClustersPerLay2/I");
  event_tree->Branch("numOfClustersPerLay3", &b_numOfClustersPerLay3, "numOfClustersPerLay3/I");
  event_tree->Branch("numOfClustersPerLay4", &b_numOfClustersPerLay4, "numOfClustersPerLay4/I");
  event_tree->Branch("numOfClustersPerDisk1", &b_numOfClustersPerDisk1, "numOfClustersPerDisk1/I");
  event_tree->Branch("numOfClustersPerDisk2", &b_numOfClustersPerDisk2, "numOfClustersPerDisk2/I");
  event_tree->Branch("numOfClustersPerDisk3", &b_numOfClustersPerDisk3, "numOfClustersPerDisk3/I");
  event_tree->Branch("numOfClustersPerDisk4", &b_numOfClustersPerDisk4, "numOfClustersPerDisk4/I");

  event_tree->Branch("numOfPixPerLay1", &b_numOfPixPerLay1, "numOfPixPerLay1/I");
  event_tree->Branch("numOfPixPerLay2", &b_numOfPixPerLay2, "numOfPixPerLay2/I");
  event_tree->Branch("numOfPixPerLay3", &b_numOfPixPerLay3, "numOfPixPerLay3/I");
  event_tree->Branch("numOfPixPerLay4", &b_numOfPixPerLay4, "numOfPixPerLay4/I");
  event_tree->Branch("numOfPixPerDisk1", &b_numOfPixPerDisk1, "numOfPixPerDisk1/I");
  event_tree->Branch("numOfPixPerDisk2", &b_numOfPixPerDisk2, "numOfPixPerDisk2/I");
  event_tree->Branch("numOfPixPerDisk3", &b_numOfPixPerDisk3, "numOfPixPerDisk3/I");
  event_tree->Branch("numOfPixPerDisk4", &b_numOfPixPerDisk4, "numOfPixPerDisk4/I");

  event_tree->Branch("numberOfDetUnits1", &b_numberOfDetUnits1, "numberOfDetUnits1/I");
  event_tree->Branch("numberOfDetUnits2", &b_numberOfDetUnits2, "numberOfDetUnits2/I");
  event_tree->Branch("numberOfDetUnits3", &b_numberOfDetUnits3, "numberOfDetUnits3/I");
  event_tree->Branch("numberOfDetUnits4", &b_numberOfDetUnits4, "numberOfDetUnits4/I");

  event_tree->Branch("numOfPixBarrel", &b_numOfPixBarrel, "numOfPixBarrel/I");
  event_tree->Branch("numOfPixForward", &b_numOfPixForward, "numOfPixForward/I");
  event_tree->Branch("numOfClusterBarrel", &b_numOfClusterBarrel, "numOfClusterBarrel/I");
  event_tree->Branch("numOfClusterForward", &b_numOfClusterForward, "numOfClusterForward/I");
  event_tree->Branch("numberOfNoneEdgePixels", &b_numberOfNoneEdgePixels, "numberOfNoneEdgePixels/I");
  event_tree->Branch("maxPixPerDet", &b_maxPixPerDet, "maxPixPerDet/I");


  
  pixel_tree = new TTree("pixel_tree", "pixel_tree");
  pixel_tree->Branch("bxid", &p_bxid, "bxid/I");
  pixel_tree->Branch("run", &p_run, "run/I");
  pixel_tree->Branch("event", &p_event, "event/I");
  pixel_tree->Branch("lumi", &p_lumi, "lumi/I");
  pixel_tree->Branch("orbit", &p_orbit, "orbit/I");
  pixel_tree->Branch("npv", &p_npv, "npv/I");
  pixel_tree->Branch("nmodule", &p_nmodule, "nmodule/I");
  pixel_tree->Branch("layer", &p_layer, "layer/I");
  pixel_tree->Branch("module", &p_module, "module/I");
  pixel_tree->Branch("ladder", &p_ladder, "ladder/I");
  pixel_tree->Branch("sector", &p_sector, "sector/I");
  pixel_tree->Branch("subid", &p_subid, "subid/I");
  pixel_tree->Branch("adc", &p_adc, "adc/I");
  pixel_tree->Branch("pixx", &p_pixx, "pixx/I");
  pixel_tree->Branch("pixy", &p_pixy, "pixy/I");
  pixel_tree->Branch("tmp1", &p_tmp1, "tmp1/F");
  pixel_tree->Branch("tmp2", &p_tmp2, "tmp2/F");
  pixel_tree->Branch("roc", &p_roc, "roc/I");
  pixel_tree->Branch("link", &p_link, "link/I");
  pixel_tree->Branch("rocInCol", &p_rocInCol, "rocInCol/I");

  pixel_tree->Branch("ch", &p_ch, "ch/F");
  pixel_tree->Branch("size", &p_size, "size/I");
  pixel_tree->Branch("sizex", &p_sizex, "sizex/I");
  pixel_tree->Branch("sizey", &p_sizey, "sizey/I");
  pixel_tree->Branch("x", &p_x, "x/F");
  pixel_tree->Branch("y", &p_y, "y/F");
  pixel_tree->Branch("minPixelRow", &p_minPixelRow, "minPixelRow/I");
  pixel_tree->Branch("maxPixelRow", &p_maxPixelRow, "maxPixelRow/I");
  pixel_tree->Branch("minPixelCol", &p_minPixelCol, "minPixelCol/I");
  pixel_tree->Branch("maxPixelCol", &p_maxPixelCol, "maxPixelCol/I");
  pixel_tree->Branch("lx", &p_lx, "lx/F");
  pixel_tree->Branch("ly", &p_ly, "ly/F");
  pixel_tree->Branch("zPos", &p_zPos, "zPos/F");
  pixel_tree->Branch("rPos", &p_rPos, "rPos/F");

  cluster_tree = new TTree("cluster_tree", "cluster_tree");
  cluster_tree->Branch("module", &c_module, "module/I");
  cluster_tree->Branch("layer", &c_layer, "layer/I");
  cluster_tree->Branch("ladder", &c_ladder, "ladder/I");
  cluster_tree->Branch("size", &c_size, "size/I");
  cluster_tree->Branch("sizex", &c_sizex, "sizex/I");
  cluster_tree->Branch("sizey", &c_sizey, "sizey/I");
  cluster_tree->Branch("ch", &c_ch, "ch/I");
  cluster_tree->Branch("subid", &c_subid, "subid/I");
  cluster_tree->Branch("x", &c_x, "x/F");
  cluster_tree->Branch("y", &c_y, "y/F");
  //  cluster_tree->Branch("z", &c_z, "z/F");
  cluster_tree->Branch("zPos", &c_zPos, "zPos/F");
  cluster_tree->Branch("rPos", &c_rPos, "rPos/F");
  cluster_tree->Branch("phiPos", &c_phiPos, "phiPos/F");
  cluster_tree->Branch("sector", &c_sector, "sector/I");


  h_total = new TH1F("total", "total", 1,0,1);
	
  for(int ilayer=0; ilayer<4; ilayer++){
    for(int ii=0; ii<20; ii++){
      TString hname = "hmap_L";
      hname += ilayer+1;
      hname += "_id";
      hname += ii;

//      if(ilayer==0){
//	hmap[ilayer][ii] = new TH2F(hname, hname, 72, -4.5, 4.5, 26, -6.5, 6.5);
//      }else if(ilayer==1){
//	hmap[ilayer][ii] = new TH2F(hname, hname, 72, -4.5, 4.5, 58, -14.5, 14.5);
//      }else if(ilayer==2){
//	hmap[ilayer][ii] = new TH2F(hname, hname, 72, -4.5, 4.5, 90, -22.5, 22.5);
//      }else if(ilayer==3){
//	hmap[ilayer][ii] = new TH2F(hname, hname, 72, -4.5, 4.5, 130, -32.5, 32.5);
//      }

      hmap[ilayer][ii] = new TH2F(hname, hname, 56, -28, 28, 30, -3.14, 3.14);


//      if(ilayer==0){
//      }else if(ilayer==1){
//	hmap[ilayer][ii] = new TH2F(hname, hname, 20, -28, 28, 12, -3.14, 3.14);
//      }
    }
  }

  //=====================================================================

//  int sizeH=200;
//  float lowH = -0.5;
//  float highH = 199.5;
//
//  hclusPerDet1 = fs->make<TH1D>( "hclusPerDet1", "Clus per det l1",
//			    sizeH, lowH, highH);
//  hclusPerDet2 = fs->make<TH1D>( "hclusPerDet2", "Clus per det l2",
//			    sizeH, lowH, highH);
//  hclusPerDet3 = fs->make<TH1D>( "hclusPerDet3", "Clus per det l3",
//			    sizeH, lowH, highH);
//  hclusPerDet4 = fs->make<TH1D>( "hclusPerDet4", "Clus per det l4",
//			    sizeH, lowH, highH);
//
//  sizeH=1000;
//  highH = 1999.5;
//  hpixPerDet1 = fs->make<TH1D>( "hpixPerDet1", "Pix per det l1",
//			    sizeH, lowH, highH);
//  hpixPerDet2 = fs->make<TH1D>( "hpixPerDet2", "Pix per det l2",
//			    sizeH, lowH, highH);
//  hpixPerDet3 = fs->make<TH1D>( "hpixPerDet3", "Pix per det l3",
//			    sizeH, lowH, highH);
//  hpixPerDet4 = fs->make<TH1D>( "hpixPerDet4", "Pix per det l4",
//			    sizeH, lowH, highH);
//
//  hpixPerDet11 = fs->make<TH1D>( "hpixPerDet11", "Pix per det l1 - ring 1",
//			    sizeH, lowH, highH);
//  hpixPerDet12 = fs->make<TH1D>( "hpixPerDet12", "Pix per det l1 - ring 2",
//			    sizeH, lowH, highH);
//  hpixPerDet13 = fs->make<TH1D>( "hpixPerDet13", "Pix per det l1 - ring 3",
//			    sizeH, lowH, highH);
//  hpixPerDet14 = fs->make<TH1D>( "hpixPerDet14", "Pix per det l1 - ring 4",
//			    sizeH, lowH, highH);
//  hpixPerDet21 = fs->make<TH1D>( "hpixPerDet21", "Pix per det l2 - ring 1",
//			    sizeH, lowH, highH);
//  hpixPerDet22 = fs->make<TH1D>( "hpixPerDet22", "Pix per det l2 - ring 2",
//			    sizeH, lowH, highH);
//  hpixPerDet23 = fs->make<TH1D>( "hpixPerDet23", "Pix per det l2 - ring 3",
//			    sizeH, lowH, highH);
//  hpixPerDet24 = fs->make<TH1D>( "hpixPerDet24", "Pix per det l2 - ring 4",
//			    sizeH, lowH, highH);
//  hpixPerDet31 = fs->make<TH1D>( "hpixPerDet31", "Pix per det l3 - ring 1",
//			    sizeH, lowH, highH);
//  hpixPerDet32 = fs->make<TH1D>( "hpixPerDet32", "Pix per det l3 - ring 2",
//			    sizeH, lowH, highH);
//  hpixPerDet33 = fs->make<TH1D>( "hpixPerDet33", "Pix per det l3 - ring 3",
//			    sizeH, lowH, highH);
//  hpixPerDet34 = fs->make<TH1D>( "hpixPerDet34", "Pix per det l3 - ring 4",
//			    sizeH, lowH, highH);
//  hpixPerDet41 = fs->make<TH1D>( "hpixPerDet41", "Pix per det l4 - ring 1",
//			    sizeH, lowH, highH);
//  hpixPerDet42 = fs->make<TH1D>( "hpixPerDet42", "Pix per det l4 - ring 2",
//			    sizeH, lowH, highH);
//  hpixPerDet43 = fs->make<TH1D>( "hpixPerDet43", "Pix per det l4 - ring 3",
//			    sizeH, lowH, highH);
//  hpixPerDet44 = fs->make<TH1D>( "hpixPerDet44", "Pix per det l4 - ring 4",
//			    sizeH, lowH, highH);
//
//  sizeH=300;
//  highH = 299.5; // charge limit in kelec
//
//  hlt1 = fs->make<TH1D>("hlt1","hlt1",512,-0.5,511.5);
//  hlt2 = fs->make<TH1D>("hlt2","hlt2",512,-0.5,511.5);
//  hlt3 = fs->make<TH1D>("hlt3","hlt3",512,-0.5,511.5);

  countEvents=0;
  countAllEvents=0;
  sumClusters=0., sumPixels=0.;
  countLumi = 0.;

  for(int i=0; i<40;++i) eventFlag[i]=true;

  //for(int ilayer = 0; ilayer<3; ++ilayer) 
  //for(int id=0;id<10;++id) {rocHits[ilayer][id]=0.; moduleHits[ilayer][id]=0.;}

}
// ------------ method called to at the end of the job  ------------
void PixClusterAna::endJob(){
  double norm = 1;
  cout << " End PixClusterAna, events all/with hits=  " 
       << countAllEvents<<"/"<<countEvents<<endl;

  if(countEvents>0) {
    sumClusters = sumClusters/float(countEvents);
    sumPixels = sumPixels/float(countEvents);
    norm = 1./float(countEvents);
  }


  cout <<" clus/pix per full event "<<sumClusters<<"/"<<sumPixels<<endl;
  cout<<" 2D plots are rescaled by the number of full events "<<countEvents<<endl;
  cout << "norm = " << norm << endl;
  
  //countLumi /= 1000.;
  //double c1=0, c2=0;
  //if(c3>0.) {c1=sumClusters/c3; c2=sumPixels/c3;}
  //cout<<" Lumi = "<<countLumi<<" still the /10 bug? "<<"clu and pix per lumi unit"<<c1<<" "<<c2<<endl;
  
  if(countEvents==0) return;

  //Divide the size histos

  // Rescale all 2D plots
//  hDetsMap1->Scale(norm);
//  hDetsMap2->Scale(norm);
//  hDetsMap3->Scale(norm);
//  hDetsMap4->Scale(norm);
//  hDetMap1->Scale(norm);
//  hDetMap2->Scale(norm);
//  hDetMap3->Scale(norm);
//  hDetMap4->Scale(norm);

  h_total->SetBinContent(1, countEvents);
  
  h_total->Write();

  for(int ilayer=0; ilayer<4; ilayer++){
    for(int ii=0; ii<20; ii++){
      hmap[ilayer][ii]->Write();
    }
  }


  

  myfile->Write();
  myfile->Close();

}
//////////////////////////////////////////////////////////////////
// Functions that gets called by framework every event
void PixClusterAna::analyze(const edm::Event& e, 
			      const edm::EventSetup& es) {
  using namespace edm;
  //const int MAX_CUT = 1000000; unused
  //  bool select = false;
  //static double pixsum = 0., clussum=0.;
  //static int nsum = 0, lsold=-999, lumiold=0.;

  // Get event setup 
  edm::ESHandle<TrackerGeometry> geom;
  es.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);

  countAllEvents++;
  int run       = e.id().run();
  int event     = e.id().event();

  int lumiBlock = e.luminosityBlock();
  int bx        = e.bunchCrossing();
  int orbit     = e.orbitNumber();

  //  std::cout << b_bxid << " " << bx << std::endl;

  //  if(abs(b_bxid-bx) == 1){
  //    std::cout << "Continuous events !!!!!!!!!!!!!!!!" << std::endl;
  //  }


  b_bxid = bx;
  b_run = run; 
  b_event = event; 
  b_lumi = lumiBlock;
  b_orbit = orbit;

  int numPVsGood = 0;
  if(run>165000) { // skip for earlier runs, crashes
    edm::Handle<reco::VertexCollection> vertices;
    e.getByToken(VertexCollectionToken, vertices);

    if( !vertices.failedToGet() && vertices.isValid() ) {
      for( reco::VertexCollection::const_iterator iVertex = vertices->begin();
	   iVertex != vertices->end(); ++iVertex ) {
	
	if(! iVertex->isValid() ) continue;
	if( iVertex->isFake() ) continue; 
	numPVsGood++;
	
	if( PRINT ){
	  cout << "vertex";
	  cout << ": x " << iVertex->x();
	  cout << ", y " << iVertex->y();
	  cout << ", z " << iVertex->z();
	  cout << ", ndof " << iVertex->ndof();
	  cout << endl;
	} // print 
      } // for loop
    } // if vertex    
  } // if run
    
  b_npv = numPVsGood;



  // Get Cluster Collection from InputTag
  edm::Handle< edmNew::DetSetVector<SiPixelCluster> > clusters;
  //e.getByLabel( src_ , clusters);
  // New By Token method
  e.getByToken( myClus , clusters);


  const edmNew::DetSetVector<SiPixelCluster>& input = *clusters;     

 // number of modules with pix
  int numOf = input.size();
  b_nmodule = numOf;


  if(PRINT) cout<<"run "<<run<<" event "<<event<<" bx "<<bx<<" lumi "<<lumiBlock<<" orbit "<<orbit<<" "
		<<numOf;


  // Analyse HLT
  //bool passHLT1=false,passHLT2=false,passHLT3=false,passHLT4=false,passHLT5=false;
  const int hltSize = 1024;
  //  bool hlt[hltSize];
  //  for(int i=0;i<hltSize;++i) hlt[i]=false;

  edm::TriggerNames TrigNames;
  edm::Handle<edm::TriggerResults> HLTResults;
  // Extract the HLT results
  //e.getByLabel(edm::InputTag("TriggerResults","","HLT"),HLTResults);
  e.getByToken(TrigResultsToken, HLTResults);

  if ((HLTResults.isValid() == true) && (HLTResults->size() > 0)) {
    //TrigNames.init(*HLTResults);
    const edm::TriggerNames & TrigNames = e.triggerNames(*HLTResults);

    if(TrigNames.triggerNames().size() > hltSize) cout<<" extend the hlt array above "<<hltSize<<endl;

    for (unsigned int i = 0; i < TrigNames.triggerNames().size(); i++) {  // loop over trigger
      if(countAllEvents==1) cout<<i<<" "<<TrigNames.triggerName(i)<<endl;

      if ( 
	   (HLTResults->wasrun(TrigNames.triggerIndex(TrigNames.triggerName(i))) == true) &&
	   (HLTResults->accept(TrigNames.triggerIndex(TrigNames.triggerName(i))) == true) &&
	   (HLTResults->error( TrigNames.triggerIndex(TrigNames.triggerName(i))) == false) ) {

	//	std::cout << "pass" << std::endl;

	//	hlt[i]=true;
	//	hlt1->Fill(float(i));

      } // if hlt
    } // loop 
  } // if valid


  // Select events if select1>0
  if(select1>0) {
    // skip events with few pixel dets
    if(select1==1) { if(numOf<select2) return; } 
    // select events only for a defined bx
    else if(select1==2) { if(bx!=select2) return; } 
    else if(select1==3) { if(  !( (bx==39)||(bx==201)||(bx==443)||(bx==499)||(bx==1083)||(bx==1337)||(bx==1492)||(bx==1977)||(bx==2231)||(bx==2287)||(bx==2871)||(bx==3224)||(bx==3280) )   ) return; } 
    else if(select1==4) { if( ( (bx==1)||(bx==39)||(bx==201)||(bx==443)||(bx==499)||(bx==1083)||(bx==1337)||(bx==1492)||(bx==1977)||(bx==2231)||(bx==2287)||(bx==2871)||(bx==3224)||(bx==3280) )   ) return; } 
    // select specific event
    else if(select1==10) { if(event!=select2) return; } 
    //....
  }
  

  //  for (unsigned int i=0;i<256;i++) if(hlt[i]==true) hlt2->Fill(float(i));


  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoH;
  es.get<TrackerTopologyRcd>().get(tTopoH);
  const TrackerTopology *tTopo=tTopoH.product();

  //---------------------------------------

  countEvents++;
  int numberOfClusters = 0;
  int numberOfPixels = 0;
  int numberOfNoneEdgePixels = 0;
  int numberOfDetUnits1 = 0;
  int numOfClustersPerDet1=0;        
  int numOfClustersPerLay1=0;        
  int numberOfDetUnits2 = 0;
  int numOfClustersPerDet2=0;        
  int numOfClustersPerLay2=0;        
  int numberOfDetUnits3 = 0;
  int numOfClustersPerDet3=0;        
  int numOfClustersPerLay3=0;        
  int numberOfDetUnits4 = 0;
  int numOfClustersPerDet4=0;        
  int numOfClustersPerLay4=0;        

  int numOfPixPerLay1=0;     
  int numOfPixPerLay2=0;     
  int numOfPixPerLay3=0;     
  int numOfPixPerLay4=0;     

  int numOfPixPerDet1=0;  
  int numOfPixPerDet2=0;  
  int numOfPixPerDet3=0;  
  int numOfPixPerDet4=0;  
      
  //int numOfPixPerLink11=0;  
  //int numOfPixPerLink12=0;  
  //int numOfPixPerLink21=0;  
  //int numOfPixPerLink22=0;  
  //int numOfPixPerLink3=0;  

  int maxClusPerDet=0;
  int maxPixPerDet=0;
  unsigned int maxPixPerClu=0;

  int numOfClustersPerDisk1=0;  
  int numOfClustersPerDisk2=0;  
  int numOfClustersPerDisk3=0;  
  int numOfClustersPerDisk4=0;  
  int numOfPixPerDisk1=0;  
  int numOfPixPerDisk2=0;  
  int numOfPixPerDisk3=0;  
  int numOfPixPerDisk4=0;  

  //static int module1[416][160] = {{0}};
  //static int module2[416][160] = {{0}};
  //static int module3[416][160] = {{0}};

  
  // get vector of detunit ids
  //--- Loop over detunits.
  edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=input.begin();
  for ( ; DSViter != input.end() ; DSViter++) {
    //bool valid = false;
    unsigned int detid = DSViter->detId();
    // Det id
    DetId detId = DetId(detid);       // Get the Detid object
    unsigned int detType=detId.det(); // det type, pixel=1
    unsigned int subid=detId.subdetId(); //subdetector type, barrel=1
 
    if(PRINT)
      cout<<"Det: "<<detId.rawId()<<" "<<detId.null()<<" "<<detType<<" "<<subid<<endl;    


    if(detType!=1) continue; // look only at pixels

    // Get the geom-detector
    const PixelGeomDetUnit * theGeomDet =
      dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );
    double detZ = theGeomDet->surface().position().z();
    double detPhi = theGeomDet->surface().position().phi();
    double detR = theGeomDet->surface().position().perp();

    
    const PixelTopology * topol = &(theGeomDet->specificTopology());

    // barrel ids
    unsigned int layerC=0;
    unsigned int ladderC=0;
    unsigned int zindex=0;
    int shell  = 0; // shell id // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
    int sector = 0; // 1-8
    int ladder = 0; // 1-22
    int layer  = 0; // 1-3
    int module = 0; // 1-4
    bool half  = false; // 

    // Endcap ids
    unsigned int disk=0; //1,2,3
    unsigned int blade=0; //1-24
    unsigned int zindexF=0; //
    unsigned int side=0; //size=1 for -z, 2 for +z
    unsigned int panel=0; //panel=1

    edmNew::DetSet<SiPixelCluster>::const_iterator clustIt;

    // Subdet id, pix barrel=1, forward=2
    if(subid==2) {  // forward

      disk=tTopo->pxfDisk(detid); //1,2,3
      blade=tTopo->pxfBlade(detid); //1-24
      zindex=tTopo->pxfModule(detid); //
      side=tTopo->pxfSide(detid); //size=1 for -z, 2 for +z
      panel=tTopo->pxfPanel(detid); //panel=1
      PixelEndcapName pen(detid,tTopo,phase1_);

      if(PRINT) cout<<" forward det, disk "<<disk<<", blade "
 		    <<blade<<", module "<<zindexF<<", side "<<side<<", panel "
 		    <<panel<<" pos = "<<detZ<<" "<<detR<<endl;
 

    } else if (subid==1) {  // barrel

      layerC=tTopo->pxbLayer(detid);
      ladderC=tTopo->pxbLadder(detid);
      zindex=tTopo->pxbModule(detid);
      //PixelBarrelName pbn(detid);
      PixelBarrelName pbn(detid,tTopo,phase1_);

      // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
      PixelBarrelName::Shell sh = pbn.shell(); //enum
      sector = pbn.sectorName();
      ladder = pbn.ladderName();
      layer  = pbn.layerName();
      module = pbn.moduleName();
      half  = pbn.isHalfModule();
      shell = int(sh);
      // change the module sign for z<0
      if(shell==1 || shell==2) module = -module;
      // change ladeer sign for Outer )x<0)
      if(shell==1 || shell==3) ladder = -ladder;
      
      if(PRINT) { 
	cout<<" Barrel layer, ladder, module "
	    <<layerC<<" "<<ladderC<<" "<<zindex<<" "
	    <<sh<<"("<<shell<<") "<<sector<<" "<<layer<<" "<<ladder<<" "
	    <<module<<" "<<half<< endl;
	//cout<<" Barrel det, thick "<<detThick<<" "
	//  <<" layer, ladder, module "
	//  <<layer<<" "<<ladder<<" "<<zindex<<endl;
	//cout<<" col/row, pitch "<<cols<<" "<<rows<<" "
	//  <<pitchX<<" "<<pitchY<<endl;
      }      
      
    } // if subid
    
    if(PRINT) {
      cout<<"List clusters : "<<endl;
      cout<<"Num Charge Size SizeX SizeY X Y Xmin Xmax Ymin Ymax Edge"
	  <<endl;
    }

    // Loop over clusters
    for (clustIt = DSViter->begin(); clustIt != DSViter->end(); clustIt++) {
      sumClusters++;
      numberOfClusters++;

      std::cout << "numberofclusteres = " << numberOfClusters << std::endl;

      float ch = float(clustIt->charge())/1000.; // convert ke to electrons
      int size = clustIt->size();
      int sizeX = clustIt->sizeX(); //x=row=rfi, 
      int sizeY = clustIt->sizeY(); //y=col=z_global
      float x = clustIt->x(); // row, cluster position in pitch units, as float (int+0.5)
      float y = clustIt->y(); // col, analog average
      //      float z = clustIt->z(); // col, analog average
      // Returns int index of the cluster min/max  
      int minPixelRow = clustIt->minPixelRow(); //x
      int maxPixelRow = clustIt->maxPixelRow();
      int minPixelCol = clustIt->minPixelCol(); //y
      int maxPixelCol = clustIt->maxPixelCol();

      //unsigned int geoId = clustIt->geographicalId(); // always 0?!
      // edge method moved to topologu class
      bool edgeHitX = (topol->isItEdgePixelInX(minPixelRow)) || 
	(topol->isItEdgePixelInX(maxPixelRow)); 
      bool edgeHitY = (topol->isItEdgePixelInY(minPixelCol)) || 
	(topol->isItEdgePixelInY(maxPixelCol)); 

      bool edgeHitX2 = false; // edge method moved 
      bool edgeHitY2 = false; // to topologu class
            
      if(PRINT) cout<<numberOfClusters<<" "<<ch<<" "<<size<<" "<<sizeX<<" "<<sizeY<<" "
		    <<x<<" "<<y<<" "<<minPixelRow<<" "<<maxPixelRow<<" "<<minPixelCol<<" "
		    <<maxPixelCol<<" "<<edgeHitX<<" "<<edgeHitY<<endl;

      // get global z position of teh cluster
      LocalPoint lp = topol->localPosition(MeasurementPoint(x,y));
      //      GlobalPoint clustgp = theGeomDet->surface().toGlobal(clustlp);
      float lx = lp.x(); // local cluster position in cm
      float ly = lp.y();


      float zPos = detZ - ly;
      float rPos = detR + lx;

      // Get the pixels in the Cluster
      const vector<SiPixelCluster::Pixel>& pixelsVec = clustIt->pixels();
      if(PRINT) cout<<" Pixels in this cluster (i/x/y/char)"<<endl;
      //bool bigInX=false, bigInY=false;
      // Look at pixels in this cluster. ADC is calibrated, in electrons
      bool edgeInX = false; // edge method moved 
      bool edgeInY = false; // to topologu class
      //bool cluBigInX = false; // does this clu include a big pixel
      //bool cluBigInY = false; // does this clu include a big pixel
      //int noisy = 0;

      if(subid==1 && layer>0){
	hmap[layer-1][countEvents]->Fill(zPos, detPhi, ch);
      }
      


      if(pixelsVec.size()>maxPixPerClu) maxPixPerClu = pixelsVec.size();
 
      for (unsigned int i = 0;  i < pixelsVec.size(); ++i) { // loop over pixels
	//bool isBig=false;
	sumPixels++;
	numberOfPixels++;
	float pixx = pixelsVec[i].x; // index as float=iteger, row index
	float pixy = pixelsVec[i].y; // same, col index
	float adc = (float(pixelsVec[i].adc)/1000.);

	int roc = rocId(int(pixy),int(pixx));  // 0-15, column, row
	int link = int(roc/8); // link 0 & 1
	int rocInCol = roc%8; // 0-7

	//int chan = PixelChannelIdentifier::pixelToChannel(int(pixx),int(pixy));

	bool bigInX = topol->isItBigPixelInX(int(pixx));
	bool bigInY = topol->isItBigPixelInY(int(pixy));
	if( !(bigInX || bigInY) ) {numberOfNoneEdgePixels++;}
	//else {isBig=true;}

	// Pixel histos

	float tmp1= float(module) - 0.5 + (0.125/2.) + (float(rocInCol) * 0.125); 
	float tmp2= float(ladder) - 0.5 + (0.5/2.)   + (float(link) * 0.5); 

	
	p_bxid = bx;
	p_run = run; 
	p_event = event; 
	p_lumi = lumiBlock;
	p_orbit = orbit;
	p_npv = numPVsGood;
	p_nmodule = numOf;
	p_layer = layer;
	p_module = module;
	p_ladder = ladder;
	p_sector = sector;
	p_subid = subid;
	p_adc = adc;
	p_pixx = pixx;
	p_pixy = pixy;
	p_tmp1 = tmp1;
	p_tmp2 = tmp2;
	p_roc = roc;
	p_link = link;
	p_rocInCol = rocInCol;

	p_ch = ch;
	p_size = size;
	p_sizex = sizeX;
	p_sizey = sizeY;
	p_x = x;
	p_y = y;
	p_minPixelRow = minPixelRow;
	p_maxPixelRow = maxPixelRow;
	p_minPixelCol = minPixelCol;
	p_maxPixelCol = maxPixelCol;
	p_lx = lx;
	p_ly = ly;
	p_zPos = zPos;
	p_rPos = rPos;

	pixel_tree->Fill();


	if (subid==1){


	  if(layer==1) {

	    //	    std::cout << "clusterid, layer, tmp1, tmp2 " << numberOfClusters << " " << layer << " " << tmp1 << " " << tmp2 << std::endl;


	    //cout<<" module "<<layer<<" "<<ladder<<" "<<module<<" "
	    //	<<pixx<<" "<<pixy<<endl;

	    numOfPixPerDet1++;
	    numOfPixPerLay1++;     




	  } else if(layer==2) {

	    numOfPixPerDet2++;
	    numOfPixPerLay2++;   

	  } else if(layer==3) {

	    numOfPixPerDet3++;
	    numOfPixPerLay3++;

	  } else if(layer==4) {

	    numOfPixPerDet4++;
	    numOfPixPerLay4++;
 

	  }  // if layer

	} else if (subid==2) {  // endcap
	  // pixels

	  if(disk==1) { // disk1 -+z
	    if(side==1) numOfPixPerDisk2++;      // d1,-z
	    else if(side==2) numOfPixPerDisk3++; // d1, +z
	    else cout<<" unknown side "<<side<<endl;

	  } else if(disk==2) { // disk2 -+z
	    
	    if(side==1) numOfPixPerDisk1++;      // d2, -z
	    else if(side==2) numOfPixPerDisk4++; // d2, +z
	    else cout<<" unknown side "<<side<<endl;

	  } else if (disk==3){

	  }else cout<<" unknown disk "<<disk<<endl;

	} // end if subdet (pixel loop)

	// Find if this cluster includes an edge pixel
	edgeInX = topol->isItEdgePixelInX(int(pixx));
	edgeInY = topol->isItEdgePixelInY(int(pixy));
	
	if(PRINT) cout<<i<<" "<<pixx<<" "<<pixy<<" "<<adc<<endl;
	//" "<<bigInX<<" "<<bigInY<<" "<<edgeInX<<" "<<edgeInY<<" "<<isBig<<endl;
	
	if(edgeInX) edgeHitX2=true;
	if(edgeInY) edgeHitY2=true; 
	//if(bigInX) cluBigInX=true;
	//if(bigInY) cluBigInY=true;

      } // pixel loop




      c_module = module;
      c_layer = layer;
      c_ladder = ladder;
      c_size = size;
      c_sizex = sizeX;
      c_sizey = sizeY;
      c_ch = ch;
      c_subid = subid;
      c_x = x;
      c_y = y;
      //      c_z = z;
      c_zPos = zPos;
      c_rPos = rPos;
      c_phiPos = detPhi;
      c_sector = sector;
      cluster_tree->Fill();

      
      // Cluster histos
      if (subid==1 ) {  // barrel
	//if (subid==1) {  // barrel

	if(layer==1) {  // layer 1
	  
	  numOfClustersPerDet1++;
	  numOfClustersPerLay1++;

	} else if(layer==2) {

	  
	  numOfClustersPerDet2++;
	  numOfClustersPerLay2++;

	} else if(layer==3) {



	  numOfClustersPerDet3++;
	  numOfClustersPerLay3++;


	} else if(layer==4) {


	  numOfClustersPerDet4++;
	  numOfClustersPerLay4++;


	} // end if layer

      } else if (subid==2) {  // endcap

	//cout<<disk<<" "<<side<<endl;
	if(disk==1) { // disk1 -+z
	  if(side==1) {
	    numOfClustersPerDisk2++;      // d1,-z
	  } else if(side==2) {
	    numOfClustersPerDisk3++; // d1, +z
	  } else cout<<" unknown side "<<side<<endl;

	} else if(disk==2) { // disk2 -+z

	  if(side==1) {
	    numOfClustersPerDisk1++;      // d2, -z
	  } else if(side==2) { 
	    numOfClustersPerDisk4++; // d2, +z
	  } else cout<<" unknown side "<<side<<endl;

	} else if(disk==3) { // disk3 -+z

	} else cout<<" unknown disk "<<disk<<endl;

      } // end barrel/forward cluster loop
      
      if( (edgeHitX != edgeHitX2) && sizeX<64) 
	cout<<" wrong egdeX "<<edgeHitX<<" "<<edgeHitX2<<endl;
      if( (edgeHitY != edgeHitY2) && sizeY<64)  
	cout<<" wrong egdeY "<<edgeHitY<<" "<<edgeHitY2<<endl;

    } // clusters 

    
    if(numOfClustersPerDet1>maxClusPerDet) maxClusPerDet = numOfClustersPerDet1;
    if(numOfClustersPerDet2>maxClusPerDet) maxClusPerDet = numOfClustersPerDet2;
    if(numOfClustersPerDet3>maxClusPerDet) maxClusPerDet = numOfClustersPerDet3;

    if(PRINT) {
      if(layer==1) 
	cout<<"Lay1: number of clusters per det = "<<numOfClustersPerDet1<<endl;
      else if(layer==2) 
	cout<<"Lay2: number of clusters per det = "<<numOfClustersPerDet2<<endl;
      else if(layer==3) 
	cout<<"Lay3: number of clusters per det = "<<numOfClustersPerDet3<<endl;
    } // end if PRINT
  
    if (subid==1) {  // barrel

      // Det histos
      if(layer==1) {
	
	++numberOfDetUnits1;
	//	hclusPerDet1->Fill(float(numOfClustersPerDet1));
	//	hpixPerDet1->Fill(float(numOfPixPerDet1));
	if(numOfPixPerDet1>maxPixPerDet) maxPixPerDet = numOfPixPerDet1;  
	//if(numOfPixPerLink11>798 || numOfPixPerLink12>798) select=true; 
	//hpixPerLink1->Fill(float(numOfPixPerLink11));
	//hpixPerLink1->Fill(float(numOfPixPerLink12));

//	if(abs(module)==1)      hpixPerDet11->Fill(float(numOfPixPerDet1));
//	else if(abs(module)==2) hpixPerDet12->Fill(float(numOfPixPerDet1));
//	else if(abs(module)==3) hpixPerDet13->Fill(float(numOfPixPerDet1));
//	else if(abs(module)==4) hpixPerDet14->Fill(float(numOfPixPerDet1));

	//	numOfClustersPerDet1=0;        
	numOfPixPerDet1=0;
	//numOfPixPerLink11=0;        
	//numOfPixPerLink12=0;        

      } else if(layer==2) {

	++numberOfDetUnits2;
	//	hclusPerDet2->Fill(float(numOfClustersPerDet2));
	//	hpixPerDet2->Fill(float(numOfPixPerDet2));
	if(numOfPixPerDet2>maxPixPerDet) maxPixPerDet = numOfPixPerDet2;  
	//hpixPerLink2->Fill(float(numOfPixPerLink21));
	//hpixPerLink2->Fill(float(numOfPixPerLink22));


//	if(abs(module)==1)      hpixPerDet21->Fill(float(numOfPixPerDet2));
//	else if(abs(module)==2) hpixPerDet22->Fill(float(numOfPixPerDet2));
//	else if(abs(module)==3) hpixPerDet23->Fill(float(numOfPixPerDet2));
//	else if(abs(module)==4) hpixPerDet24->Fill(float(numOfPixPerDet2));

	//	numOfClustersPerDet2=0;
	numOfPixPerDet2=0;        
	//numOfPixPerLink21=0;        
	//numOfPixPerLink22=0;        

      } else if(layer==3) {

	++numberOfDetUnits3;
	//	hclusPerDet3->Fill(float(numOfClustersPerDet3));
	//	hpixPerDet3->Fill(float(numOfPixPerDet3));
	if(numOfPixPerDet3>maxPixPerDet) maxPixPerDet = numOfPixPerDet3;  

//	if(abs(module)==1)      hpixPerDet31->Fill(float(numOfPixPerDet3));
//	else if(abs(module)==2) hpixPerDet32->Fill(float(numOfPixPerDet3));
//	else if(abs(module)==3) hpixPerDet33->Fill(float(numOfPixPerDet3));
//	else if(abs(module)==4) hpixPerDet34->Fill(float(numOfPixPerDet3));

// 	if(      module==-2&& ladder==13)  hpixPerDet100->Fill(float(numOfPixPerDet3));
// 	else if( module==+1&& ladder==11)  hpixPerDet101->Fill(float(numOfPixPerDet3));
// 	else if( module==-1&& ladder==9 )  hpixPerDet102->Fill(float(numOfPixPerDet3));
// 	else if( module==-4&& ladder==11)  hpixPerDet103->Fill(float(numOfPixPerDet3));
// 	else if( module==-2&& ladder==11)  hpixPerDet104->Fill(float(numOfPixPerDet3));
// 	else if( module==-2&& ladder==19)  hpixPerDet105->Fill(float(numOfPixPerDet3));

// 	if( numOfPixPerDet3>191 ) {
// 	  cout<<" Layer 3 module "<<ladder<<" "<<module<<" "<<numOfPixPerDet3<<endl;
// 	  select = true;
// 	}

//	numOfClustersPerDet3=0;
	numOfPixPerDet3=0;        
	//numOfPixPerLink3=0;        

      } else if(layer==4) {

	++numberOfDetUnits4;
	//	hclusPerDet4->Fill(float(numOfClustersPerDet4));
	//	hpixPerDet4->Fill(float(numOfPixPerDet4));
	if(numOfPixPerDet4>maxPixPerDet) maxPixPerDet = numOfPixPerDet4;  

	//	numOfClustersPerDet4=0;
	numOfPixPerDet4=0;        
	//numOfPixPerLink4=0;        

      } // layer
      
    } // end barrel/forward
    
  } // detunits loop

  
  if( PRINT) {  // 
    cout<<"run "<<run<<" event "<<event<<" bx "<<bx<<" lumi "<<lumiBlock<<" orbit "<<orbit<<" "<<countEvents<<endl;   
    cout<<"Num of pix "<<numberOfPixels<<" num of clus "<<numberOfClusters<<" max clus per det "
	<<maxClusPerDet<<" max pix per clu "<<maxPixPerClu<<" count "
	<<countEvents<<endl;
    cout<<"Number of clusters per Lay1,2,3: "<<numOfClustersPerLay1<<" "
	<<numOfClustersPerLay2<<" "<<numOfClustersPerLay3<<endl;
    cout<<"Number of pixels per Lay1,2,3: "<<numOfPixPerLay1<<" "
	<<numOfPixPerLay2<<" "<<numOfPixPerLay3<<endl;
    cout<<"Number of dets with clus in Lay1,2,3: "<<numberOfDetUnits1<<" "
	<<numberOfDetUnits2<<" "<<numberOfDetUnits3<<endl;
  } // if PRINT
  
  b_countEvents = countEvents;
  b_numberOfPixels = numberOfPixels;
  b_numberOfClusters = numberOfClusters;
  b_maxClusPerDet = maxClusPerDet;
  b_maxPixPerClu = maxPixPerClu;

  b_numOfClustersPerLay1 = numOfClustersPerLay1;
  b_numOfClustersPerLay2 = numOfClustersPerLay2;
  b_numOfClustersPerLay3 = numOfClustersPerLay3;
  b_numOfClustersPerLay4 = numOfClustersPerLay4;
  b_numOfClustersPerDisk1 = numOfClustersPerDisk1;
  b_numOfClustersPerDisk2 = numOfClustersPerDisk2;
  b_numOfClustersPerDisk3 = numOfClustersPerDisk3;
  b_numOfClustersPerDisk4 = numOfClustersPerDisk4;

  b_numOfPixPerLay1 = numOfPixPerLay1;
  b_numOfPixPerLay2 = numOfPixPerLay2;
  b_numOfPixPerLay3 = numOfPixPerLay3;
  b_numOfPixPerLay4 = numOfPixPerLay4;
  b_numOfPixPerDisk1 = numOfPixPerDisk1;
  b_numOfPixPerDisk2 = numOfPixPerDisk2;
  b_numOfPixPerDisk3 = numOfPixPerDisk3;
  b_numOfPixPerDisk4 = numOfPixPerDisk4;

  b_numberOfDetUnits1 = numberOfDetUnits1;
  b_numberOfDetUnits2 = numberOfDetUnits2;
  b_numberOfDetUnits3 = numberOfDetUnits3;
  b_numberOfDetUnits4 = numberOfDetUnits4;

  int pixb = numOfPixPerLay1+numOfPixPerLay2+numOfPixPerLay3;
  int pixf = numOfPixPerDisk1 + numOfPixPerDisk2 + numOfPixPerDisk3 + numOfPixPerDisk4; 
  int clusb = numOfClustersPerLay1+numOfClustersPerLay2+numOfClustersPerLay3+numOfClustersPerLay4;
  int clusf = numOfClustersPerDisk1+numOfClustersPerDisk2+numOfClustersPerDisk3+numOfClustersPerDisk4;

  b_numOfPixBarrel = pixb;
  b_numOfPixForward = pixf;
  b_numOfClusterBarrel = clusb;
  b_numOfClusterForward = clusf;
  b_pixb = pixb;
  b_numberOfNoneEdgePixels = numberOfNoneEdgePixels;
  b_maxPixPerDet = maxPixPerDet;


  event_tree->Fill();  

  //  for (unsigned int i=0;i<256;i++) if(hlt[i]) hlt3->Fill(float(i));
    
      
} // end 

//define this as a plug-in
DEFINE_FWK_MODULE(PixClusterAna);
