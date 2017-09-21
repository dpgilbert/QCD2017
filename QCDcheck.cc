#include "QCDcheck.h"

using namespace std;
using namespace mt2;
using namespace duplicate_removal;

void QCDcheck::loop (const TChain * chain, const char* sample, const char* outdir, const int maxEvt) {
  TString outname = Form("%s/%s.root",outdir,sample);
  TFile * outfile = new TFile(outname.Data(), "RECREATE");

  const char* lumi17 = "11.379";

  const int nmt2_bins = 12;
  const float mt2_min = 0;
  const float mt2_max = 600;

  const int nht_bins = 120;
  const float ht_min = 300;
  const float ht_max = 4500;

  const int njpt_bins = 40;
  const float jpt_min = 30;
  const float jpt_max = 2030;

  const int nmet_bins = 50;
  const float met_min = 0;
  const float met_max = 1000;

  const int nphi_bins = 80;
  const float phi_min = -3.14159;
  const float phi_max = 3.14159;

  // Book Histograms
  // MHT anomaly cross-check histograms 
  // plot MT2-binned histogram down to 0
  // HT-binned histogram with no MT2 cut
  // HT-binned histogram with MT2 > 200 (already have this one, make for new data)
  // jet 1 pt wihtout MT2 cut
  // jet 2 pt without MT2 cut
  // met_pt without MT2 cut, with HT > 1000 so we can get down to the MET > 30 region
  // met_phi without MT2 cut, with HT > 1000

  TString h_mt2_title = Form("M_{T2} Distribution for #Delta#phi_{min} < 0.3;M_{T2} (GeV);Count @ %s fb^{-1}", lumi17);
  TH1F h_mt2("h_mt2", h_mt2_title.Data(), nmt2_bins, mt2_min, mt2_max);
  TString h_ht_title = Form("H_{T} Distribution for #Delta#phi_{min} < 0.3;H_{T} (GeV);Count @ %s fb^{-1}", lumi17);
  TH1F h_ht ("h_ht", h_ht_title.Data(), nht_bins, ht_min, ht_max);
  TString h_ht_mt2cut_title = Form("H_{T} Distribution for #Delta#phi_{min} < 0.3; H_{T} (GeV);Count @ %s fb^{-1}", lumi17);
  TH1F h_ht_mt2cut ("h_ht_mt2cut", h_ht_mt2cut_title.Data(), nht_bins, ht_min, ht_max);
  TString h_j1pt_title = Form("Jet 1 p_{T} Distribution for #Delta#phi_{min} < 0.3;Jet 1 p_{T} (GeV);Count @ %s fb^{-1}", lumi17);
  TH1F h_j1pt("h_j1pt", h_j1pt_title.Data(), njpt_bins, jpt_min, jpt_max);
  TString h_j2pt_title = Form("Jet 2 p_{T} Distribution for #Delta#phi_{min} < 0.3;Jet 2 p_{T} (GeV);Count @ %s fb^{-1}", lumi17);
  TH1F h_j2pt("h_j2pt", h_j2pt_title.Data(), njpt_bins, jpt_min, jpt_max);
  TString h_met_pt_title = Form("MET Distribution for #Delta#phi_{min} < 0.3;MET (GeV);Count @ %s fb^{-1}",lumi17);
  TH1F h_met_pt("h_met_pt", h_met_pt_title.Data(), nmet_bins, met_min, met_max);
  TString h_met_phi_title = Form("#phi_{MET} Distribution for #Delta#phi_{min} < 0.3;#phi_{MET};Count @ %s fb^{-1}", lumi17);
  TH1F h_met_phi("h_met_phi", h_met_phi_title.Data(), nphi_bins, phi_min, phi_max);
  TString h_jet1_phi_title = Form("Jet 1 #phi# Distribution;#phi_{Jet 1};Count @ %s fb^{-1}", lumi17);
  TH1F h_jet1_phi("h_jet1_phi", h_jet1_phi_title.Data(), nphi_bins, phi_min, phi_max);

  TH1F h_mt2_hi("h_mt2_hi", h_mt2_title.Data(), nmt2_bins, mt2_min, mt2_max);
  TH1F h_ht_hi ("h_ht_hi", h_ht_title.Data(), nht_bins, ht_min, ht_max);
  TH1F h_ht_mt2cut_hi ("h_ht_mt2cut_hi", h_ht_mt2cut_title.Data(), nht_bins, ht_min, ht_max);
  TH1F h_j1pt_hi("h_j1pt_hi", h_j1pt_title.Data(), njpt_bins, jpt_min, jpt_max);
  TH1F h_j2pt_hi("h_j2pt_hi", h_j2pt_title.Data(), njpt_bins, jpt_min, jpt_max);
  TH1F h_met_pt_hi("h_met_pt_hi", h_met_pt_title.Data(), nmet_bins, met_min, met_max);
  TH1F h_met_phi_hi("h_met_phi_hi", h_met_phi_title.Data(), nphi_bins, phi_min, phi_max);
  TH1F h_jet1_phi_hi("h_jet1_phi_hi", h_jet1_phi_title.Data(), nphi_bins, phi_min, phi_max);

  TH1::SetDefaultSumw2();

  TTree * tree = (TTree*)chain->Clone("mt2");
  mt2tree mt2_tree;
  mt2_tree.Init(tree);

  int available = tree->GetEntries();
  int max_event = maxEvt < 0 ? available : min(maxEvt, available);
  
  cout << "Running on " << max_event << " events." << endl;
  cout << "Tree contains " << available << " events." << endl;

  //  const char* json_file = "../babymaker/jsons/Cert_294927-299420_13TeV_PromptReco_Collisions17_JSON_snt.txt";
  const char* json_file = "../babymaker/jsons/Cert_294927-301567_13TeV_PromptReco_Collisions17_JSON_snt.txt";
  set_goodrun_file(json_file);
  
  eventFilter metFilterTxt;
  cout<<"Loading bad event files ..."<<endl;
  // updated lists for full dataset
  metFilterTxt.loadBadEventList("/nfs-6/userdata/mt2utils/csc2015_Dec01.txt");
  metFilterTxt.loadBadEventList("/nfs-6/userdata/mt2utils/ecalscn1043093_Dec01.txt");
  metFilterTxt.loadBadEventList("/nfs-6/userdata/mt2utils/badResolutionTrack_Jan13.txt");
  metFilterTxt.loadBadEventList("/nfs-6/userdata/mt2utils/muonBadTrack_Jan13.txt");
  cout<<" ... finished!"<<endl;
  
  int nDup = 0;
  outfile->cd();

  for (int event = 0; event < max_event; event++) {
    if ( (event+1) % 100000 == 0) cout << "Processing event " << event+1 << ", " << (100. * event) / max_event << "% complete" << endl;

    mt2_tree.GetEntry(event);

    const int run_ = mt2_tree.run;
    const int evt_ = mt2_tree.evt;
    const int lumi_ = mt2_tree.lumi;

    DorkyEventIdentifier id (run_, evt_, lumi_);
    if (is_duplicate(id)) {
      nDup++;
      continue;
    }

    if (!goodrun(run_, lumi_)) continue;
    if (mt2_tree.nVert == 0) continue;
    
    // Filters
    if (!mt2_tree.Flag_globalSuperTightHalo2016Filter) continue; 
    //    if (!mt2_tree.Flag_badMuonFilterV2) continue;
    if (!mt2_tree.Flag_eeBadScFilter) continue; 

    if (!mt2_tree.Flag_goodVertices) continue;
    if (!mt2_tree.Flag_HBHENoiseFilter) continue;
    if (!mt2_tree.Flag_HBHENoiseIsoFilter) continue;
    if (!mt2_tree.Flag_EcalDeadCellTriggerPrimitiveFilter) continue;
    //    if (!mt2_tree.Flag_badChargedHadronFilterV2) continue; 
    
    const float met_pt_ = mt2_tree.met_pt;
    const float ht_ = mt2_tree.ht;

    if (isinf(met_pt_) || isnan(met_pt_) || isinf(ht_) || isnan(ht_)) continue;

    if (ht_ < 300) continue;    
    if (met_pt_ < 30 || (ht_ < 1200 && met_pt_ < 300)) continue;

    const int njet_ = mt2_tree.nJet30;
    if (njet_ < 2) continue;

    if (mt2_tree.jet_pt[0] > 13000.) continue;

    bool clean = mt2_tree.met_miniaodPt / mt2_tree.met_caloPt <= 5 && mt2_tree.nJet200MuFrac50DphiMet == 0;

    if (!clean) continue;
    
    const float w_ = 1; 
    const float met_phi_ = mt2_tree.met_phi;
    const int nlepveto_ = mt2_tree.nMuons10 + mt2_tree.nElectrons10 + mt2_tree.nPFLep5LowMT + mt2_tree.nPFHad10LowMT;
    const float mt2_ = mt2_tree.mt2;
    const float jet1_pt_ = mt2_tree.jet1_pt;
    const float jet2_pt_ = mt2_tree.jet2_pt;
    const float deltaPhiMin_ = mt2_tree.deltaPhiMin;
    const float diffMetMhtOverMet_ = mt2_tree.diffMetMht / met_pt_;
    
    if (nlepveto_ > 0 || diffMetMhtOverMet_ > 0.5) continue;
    
    // QCD
    bool passQCDtrig = (mt2_tree.HLT_PFHT1050 || mt2_tree.HLT_PFJet500 || mt2_tree.HLT_PFHT800_PFMET75_PFMHT75 || mt2_tree.HLT_PFHT500_PFMET100_PFMHT100 || mt2_tree.HLT_PFMET120_PFMHT120 || mt2_tree.HLT_PFMETNoMu120_PFMHTNoMu120) && (mt2_tree.met_miniaodPt / mt2_tree.met_caloPt) <= 5.0 && mt2_tree.nJet200MuFrac50DphiMet == 0;
    
    if (! passQCDtrig) continue;

    if (deltaPhiMin_ > 0.3) continue;

    // fill histos
    h_mt2.Fill(mt2_, w_);
    h_ht.Fill(ht_, w_);
    if (mt2_ > 200) h_ht_mt2cut.Fill(ht_, w_); 
    h_j1pt.Fill(jet1_pt_, w_);
    h_j2pt.Fill(jet2_pt_, w_);
    if (ht_ > 1200) {
      h_met_pt.Fill(met_pt_, w_);
      h_met_phi.Fill(met_phi_, w_);
      h_jet1_phi.Fill(mt2_tree.jet_phi[0], w_);
    }    

    if (met_pt_ > 300) {
      h_mt2_hi.Fill(mt2_, w_);
      h_ht_hi.Fill(ht_, w_);
      if (mt2_ > 200) h_ht_mt2cut_hi.Fill(ht_, w_); 
      h_j1pt_hi.Fill(jet1_pt_, w_);
      h_j2pt_hi.Fill(jet2_pt_, w_);
      if (ht_ > 1200) {
	h_met_pt_hi.Fill(met_pt_, w_);
	h_met_phi_hi.Fill(met_phi_, w_);
	h_jet1_phi_hi.Fill(mt2_tree.jet_phi[0], w_);
      }
    }    
  } // event loop

  h_mt2.Write();
  h_ht.Write();
  h_ht_mt2cut.Write();
  h_j1pt.Write();
  h_j2pt.Write();
  h_met_pt.Write();
  h_met_phi.Write();
  h_jet1_phi.Write();

  h_mt2_hi.Write();
  h_ht_hi.Write();
  h_ht_mt2cut_hi.Write();
  h_j1pt_hi.Write();
  h_j2pt_hi.Write();
  h_met_pt_hi.Write();
  h_met_phi_hi.Write();
  h_jet1_phi_hi.Write();

  outfile->Close();

  cout << "Number of duplicates was: " << nDup << endl;
  
} // loop function


int main (int argc, char ** argv) {
  int maxEvt = -1;
  if (argc > 4) maxEvt = atoi(argv[4]);

  TChain* ch = new TChain("mt2");
  
  TString infile_name = Form("%s/%s*.root", argv[1], argv[2]);
  ch->Add(infile_name);

  QCDcheck * l = new QCDcheck();
  l->loop(ch, argv[2], argv[3], maxEvt);
  return 0;
}

QCDcheck::QCDcheck() {}

QCDcheck::~QCDcheck() {};
