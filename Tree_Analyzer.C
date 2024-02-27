
//#include "call_libraries.h"  // call libraries from ROOT and C++
#include "JetCorrector.h" // reader for JEC
#include <TROOT.h>       // ROOT system header
#include <TFile.h>       // ROOT file I/O
#include <TTree.h>       // ROOT tree
#include <TH1.h>  // ROOT 1D histograms
#include <TH1D.h>
#include <TH2.h>         // ROOT 2D histograms
#include <TCanvas.h>     // ROOT canvas for plotting
#include <TGraph.h>      // ROOT graph
#include <TGraph2D.h>    // ROOT 2D graph
#include <TF1.h>         // ROOT 1D functions
#include <TMath.h>       // ROOT math functions
#include <TStyle.h>      // ROOT style settings
#include <TLegend.h>     // ROOT legends
#include <TPad.h>        // ROOT pads for arranging plots
#include <TMathText.h>   // ROOT text rendering
#include <TChain.h>      // ROOT chain for combining TTree(s)
#include <TString.h>     // ROOT string class
#include <iostream>      // C++ standard I/O
using namespace std;


void Tree_Analyzer(Int_t kFile = 270, bool isMC = 1)
{
  int is_MC=1;

  vector<string> Files;
Files.push_back("Pbgoing_files/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt");
 
JetCorrector JEC(Files);

  TH1::SetDefaultSumw2();

  	
	
  
  TH1D* hEvents = new TH1D("hEvents", "", 10, 0, 10);
  TH1D* hpthat = new TH1D("hpthat", "", 200, 0, 600.);
  TH1D* hpthatW = new TH1D("hpthatW", "", 200, 0, 600.);
  TH1D* hZvtx = new TH1D("hZvtx", "", 200, -20, 20);

  
  TH1D* hgenpt= new TH1D("genpt","",200,0.,1000.);
  TH1D* hgeneta=new TH1D("geneta","",50,-3.,3.);
  TH1D* hgenphi = new TH1D("genphi", "", 64, -TMath::Pi(), TMath::Pi());

  
  TH1D* Eta_leading=new TH1D("Eta_leading","",20,-5,5);
  TH1D* Eta_Subleading=new TH1D("Eta_Subleading","",20,-5,5);
  
  TH1D* Phi_leading=new TH1D("Phi_leading","",20,-TMath::Pi(),TMath::Pi());
  TH1D* Phi_Subleading=new TH1D("Phi_Subleading","",20,-TMath::Pi(),TMath::Pi());
  
  TH1D* pT_leading=new TH1D("pT_leading","",30,0,1000);
  TH1D* pT_Subleading=new TH1D("pT_Subleading","",20,0,1000);

  TH1D* Eta_leading_reco=new TH1D("Eta_leading_reco","",20,-5,5);
  TH1D* Eta_Subleading_reco=new TH1D("Eta_Subleading_reco","",20,-5,5);

  TH1D* Phi_leading_reco=new TH1D("Phi_leading_reco","",20,-TMath::Pi(),TMath::Pi());
  TH1D* Phi_Subleading_reco=new TH1D("Phi_Subleading_reco","",20,-TMath::Pi(),TMath::Pi());

  TH1D* pT_leading_reco=new TH1D("pT_leading_reco","",50,0,1000);
  TH1D* pT_Subleading_reco=new TH1D("pT_Subleading_reco","",50,0,1000);
  
  
  TH1D* pT_dijet=new TH1D("pT_dijet","",20,0,1000);
  TH1D* Phi_dijet=new TH1D("Phi_dijet","",20,(2*TMath::Pi()/3),2*TMath::Pi());
  TH1D* X_dijet=new TH1D("X_dijet","",40,0.,1.);

  TH1D* pT_dijet_reco=new TH1D("pT_dijet_reco","",20,0,1000);
  TH1D* Phi_dijet_reco=new TH1D("Phi_dijet_reco","",20,(2*TMath::Pi()/3),2*TMath::Pi());
  TH1D* X_dijet_reco=new TH1D("X_dijet_reco","",40,0.,1.);


  TH1D* pT_jet=new TH1D("pT_jet","",20,0,1000);
  TH1D* Phi_jet=new TH1D("Phi_jet","",20,(2*TMath::Pi()/3),2*TMath::Pi());
  TH1D* X_jet=new TH1D("X_jet","",40,0.,1.);

  TH1D* pT_jet_reco=new TH1D("pT_jet_reco","",50,0,1000);
  TH1D* pT_jet_corr=new TH1D("pT_jet_corr","",20,0,1000);
  TH1D* Phi_jet_reco=new TH1D("Phi_jet_reco","",20,(2*TMath::Pi()/3),2*TMath::Pi());
  TH1D* X_jet_reco=new TH1D("X_jet_reco","",40,0.,1.);

  

  Double_t binEdges1[] = {-2.915,-2.6333,-2.07,-1.7883,-1.5067,-1.225,-0.9433,-0.6617,-0.38,-0.0983,0.183,0.465,0.7467,1.0283,1.31,1.5917,1.8733,2.4367,3.0};
  Double_t binEdges2[] = {-2.6333,-2.07,-1.7883,-1.5067,-1.225,-0.9433,-0.6617,-0.38,-0.0983,0.183,0.465,0.7467,1.0283,1.31,1.5917,1.8733,2.4367,3.0};
  Double_t binEdges3[] = {-2.07,-1.7883,-1.5067,-1.225,-0.9433,-0.6617,-0.38,-0.0983,0.183,0.465,0.7467,1.0283,1.31,1.5917,1.8733,2.4367,3.0};

  TH1D* Eta_dijet=new TH1D("Eta_dijet","",18,-5,5);
  TH1D* Eta_dijet_reco= new TH1D("Eta_dijet_reco","",18,-5,5);
  
  TH1D* Eta_dijet_25 = new TH1D("Eta_dijet_25", "", 18, -5,5);
  TH1D* Eta_dijet_50 = new TH1D("Eta_dijet_50", "", 18,-5,5);
  TH1D* Eta_dijet_70 = new TH1D("Eta_dijet_70", "", 18, -5,5);
  TH1D* Eta_dijet_90 = new TH1D("Eta_dijet_90", "", 18, -5,5);
  TH1D* Eta_dijet_110 = new TH1D("Eta_dijet_110", "", 18,-5,5);
  TH1D* Eta_dijet_130 = new TH1D("Eta_dijet_130", "", 18,-5,5);

  TH1D* Eta_dijet_reco_25 = new TH1D("Eta_dijet_reco_25", "", 18, -5,5);
  TH1D* Eta_dijet_reco_50 = new TH1D("Eta_dijet_reco_50", "", 18,-5,5);
  TH1D* Eta_dijet_reco_70 = new TH1D("Eta_dijet_reco_70", "", 18, -5,5);
  TH1D* Eta_dijet_reco_90 = new TH1D("Eta_dijet_reco_90", "", 18, -5,5);
  TH1D* Eta_dijet_reco_110 = new TH1D("Eta_dijet_reco_110", "", 18,-5,5);
  TH1D* Eta_dijet_reco_130 = new TH1D("Eta_dijet_reco_130", "", 18,-5,5);

  TH1D* JES_25 = new TH1D("JES_25","",40,-0.2,0.2);
  TH1D* JES_50 = new TH1D("JES_50","",40,-0.2,0.2);
  TH1D* JES_100 = new TH1D("JES_100","",40,-0.2,0.2);
  TH1D* JES_130 = new TH1D("JES_130","",40,-0.2,0.2);
  TH1D* JES_160= new TH1D("JES_160","",40,-0.2,0.2);
  TH1D* JES_200= new TH1D("JES_200","",40,-0.2,0.2);

  TH1D* Eta_jet_25 = new TH1D("Eta_jet_25", "", 18, -5,5);
  TH1D* Eta_jet_50 = new TH1D("Eta_jet_50", "", 18,-5,5);
  TH1D* Eta_jet_70 = new TH1D("Eta_jet_70", "", 18, -5,5);
  TH1D* Eta_jet_90 = new TH1D("Eta_jet_90", "", 18, -5,5);
  TH1D* Eta_jet_110 = new TH1D("Eta_jet_110", "", 18,-5,5);
  TH1D* Eta_jet_130 = new TH1D("Eta_jet_130", "", 18,-5,5);

  TH1D* Eta_jet_reco_25 = new TH1D("Eta_jet_reco_25", "", 18, -5,5);
  TH1D* Eta_jet_reco_50 = new TH1D("Eta_jet_reco_50", "", 18,-5,5);
  TH1D* Eta_jet_reco_70 = new TH1D("Eta_jet_reco_70", "", 18, -5,5);
  TH1D* Eta_jet_reco_90 = new TH1D("Eta_jet_reco_90", "", 18, -5,5);
  TH1D* Eta_jet_reco_110 = new TH1D("Eta_jet_reco_110", "", 18,-5,5);
  TH1D* Eta_jet_reco_130 = new TH1D("Eta_jet_reco_130", "", 18,-5,5);

  TH1D* Eta_jet_corr_25 = new TH1D("Eta_jet_corr_25", "", 18, -5,5);
  TH1D* Eta_jet_corr_50 = new TH1D("Eta_jet_corr_50", "", 18,-5,5);
  TH1D* Eta_jet_corr_70 = new TH1D("Eta_jet_corr_70", "", 18, -5,5);
  TH1D* Eta_jet_corr_90 = new TH1D("Eta_jet_corr_90", "", 18, -5,5);
  TH1D* Eta_jet_corr_110 = new TH1D("Eta_jet_corr_110", "", 18,-5,5);
  TH1D* Eta_jet_corr_130 = new TH1D("Eta_jet_corr_130", "", 18,-5,5);

 
  
  


  
    TH1D* hjtCorrpt = new TH1D("hjtCorrpt", "", 200, 0., 1000.);
  TH1D* hjtEta = new TH1D("hjtEta", "", 50, -2.5, 2.5);
  TH1D* hjtPhi = new TH1D("hjtPhi", "", 64, -TMath::Pi(), TMath::Pi());
  TH1D* hjtWTAEta = new TH1D("hjtWTAEta", "", 50, -2.5, 2.5);
  TH1D* hjtWTAPhi = new TH1D("hjtWTAPhi", "", 64, -TMath::Pi(), TMath::Pi());

  const int DR_nbins = 18;
  double DR_bins_edge[DR_nbins+1] = {0.0, 0.0035, 0.007, 0.011, 0.016, 0.022, 0.029, 0.037, 0.046, 0.056, 0.067, 0.080, 0.095, 0.115, 0.14, 0.17, 0.22, 0.30, 0.4};

  TH1D* hDeltaR_axis = new TH1D("hDeltaR_axis", "", DR_nbins, DR_bins_edge);

  const int nEtaBins = 201;
  const int nPhiBins = 199;
  TH1D* hDeltaEta_axis = new TH1D("hDeltaEta_axis", "", nEtaBins, -0.5, 0.5);
  TH1D* hDeltaPhi_axis = new TH1D("hDeltaPhi_axis", "", nPhiBins, -0.5, 0.5);

  std::vector<std::string> githubTxtUrls = {
          "Pbgoing_files/MB_PD1_Pbgoing.txt",
	  "Pbgoing_files/MB_PD2_Pbgoing.txt",
	  //  "Pbgoing_files/MB_PD3_Pbgoing.txt"
	  //  "Pbgoing_files/MB_PD4_Pbgoing.txt"
    //  "MB_PD5_pgoing.txt",
    //  "MB_PD6_pgoing.txt",
       //  "MB_PD7_pgoing.txt",
     //  "MB_PD8_pgoing.txt"
  };

  //  float pthat_values[]={120.,220.,370.,540.,15.,280.,80.,460.,170.,30.,50.};

  // int nfiles=11;
  int nfile=0;
  for (const auto& txtUrl : githubTxtUrls)
    { nfile++;
      
      
      std::ifstream file(txtUrl);
      std::stringstream ss;
      if (!file.is_open()) {
	std::cerr << "Error opening file from GitHub: " << txtUrl << std::endl;
	continue;
      }
      ss<<file.rdbuf();
      file.close();
      std::string line;
     
      while (std::getline(ss,line)){
       
	std::cout<<"File no."<<" "<<nfile<<" "<<"Root file address"<<line<<endl;
       
       
      
      TFile *f =TFile::Open(line.c_str(),"READ");
      cout<<"file opened";
      if (!f) continue;    



      
      
      TTree *hea_tree = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
      if (!hea_tree) continue;
      hea_tree->SetBranchStatus("*", 0);

      Float_t vz;
      Float_t weight;
      Float_t pthat;

      hea_tree->SetBranchStatus("vz", 1);
      hea_tree->SetBranchAddress("vz", &vz);

      /*      if(is_MC)
	{
	  hea_tree->SetBranchStatus("weight", 1);
	  hea_tree->SetBranchAddress("weight", &weight);
	  
	  hea_tree->SetBranchStatus("pthat", 1);
	  hea_tree->SetBranchAddress("pthat", &pthat);
	  } */

      TTree *hlt_tree = (TTree*)f->Get("hltanalysis/HltTree");
      if (!hlt_tree) continue;
      hlt_tree->SetBranchStatus("*", 0);

      // Int_t HLT_HIAK4CaloJet80_v1;

      //  hlt_tree->SetBranchStatus("HLT_HIAK4CaloJet80_v1", 1);
      //  hlt_tree->SetBranchAddress("HLT_HIAK4CaloJet80_v1", &HLT_HIAK4CaloJet80_v1);

      TTree *ski_tree = (TTree*)f->Get("skimanalysis/HltTree");
      if (!ski_tree) continue;
      ski_tree->SetBranchStatus("*", 0);

      Int_t HBHENoiseFilterResultRun2Loose;
      Int_t pPAprimaryVertexFilter;
      Int_t pBeamScrapingFilter;


      //     TTree *_tree = (TTree*)f->Get("skimanalysis/HltTree");
            Int_t phfCoincFilter;
       Int_t pVertexFilterCutdz1p0;
      

      ski_tree->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
      ski_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose);

      ski_tree->SetBranchStatus("pPAprimaryVertexFilter", 1);
      ski_tree->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter);

      ski_tree->SetBranchStatus("pBeamScrapingFilter", 1);
      ski_tree->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter);


            ski_tree->SetBranchStatus("phfCoincFilter",1);
        ski_tree->SetBranchAddress("phfCoincFilter",&phfCoincFilter);


       ski_tree->SetBranchStatus("pVertexFilterCutdz1p0",1);
       ski_tree->SetBranchAddress("pVertexFilterCutdz1p0",&pVertexFilterCutdz1p0);
      TTree *jet_tree = (TTree*)f->Get("ak4PFJetAnalyzer/t");
      if (!jet_tree) continue;
      jet_tree->SetBranchStatus("*", 0);

      const int nmax = 99999;
      Int_t nref;
      Float_t rawpt[nmax];
      Float_t trackMax[nmax];
      Float_t jtpt[nmax];
      Float_t jteta[nmax];
      Float_t jtphi[nmax];
      Float_t WTAeta[nmax];
      Float_t WTAphi[nmax];

      Float_t genpt[nmax];
      Float_t geneta[nmax];
      Float_t genphi[nmax];
      
      jet_tree->SetBranchStatus("nref", 1);
      jet_tree->SetBranchAddress("nref", &nref);

      jet_tree->SetBranchStatus("rawpt", 1);
      jet_tree->SetBranchAddress("rawpt", &rawpt);

      jet_tree->SetBranchStatus("trackMax", 1);
      jet_tree->SetBranchAddress("trackMax", &trackMax);

      //jet_tree->SetBranchStatus("genpt",1);
      // jet_tree->SetBranchAddress("genpt",&genpt);

      // jet_tree->SetBranchStatus("geneta",1);
      //jet_tree->SetBranchAddress("geneta",&geneta);

      //jet_tree->SetBranchStatus("genphi",1);
      // jet_tree->SetBranchAddress("genphi",&genphi);

      // jet_tree->SetBranchStatus("jtpt", 1);
      // jet_tree->SetBranchAddress("jtpt", &jtpt);

      jet_tree->SetBranchStatus("jteta", 1);
      jet_tree->SetBranchAddress("jteta", &jteta);

      jet_tree->SetBranchStatus("jtphi", 1);
      jet_tree->SetBranchAddress("jtphi", &jtphi);

      jet_tree->SetBranchStatus("WTAeta", 1);
      jet_tree->SetBranchAddress("WTAeta", &WTAeta);

      jet_tree->SetBranchStatus("WTAphi", 1);
      jet_tree->SetBranchAddress("WTAphi", &WTAphi);
      
      
      int nevents = hea_tree->GetEntries(); // number of events
      
      cout << "Total number of events in those files: "<< nevents << endl;

     
	  
      

     
      
        for(int i = 0; i < nevents; i++) //event loop start
      // for(int i = 0; i < 1000; i++) //event loop start
	{
	  hea_tree->GetEntry(i);
	  hlt_tree->GetEntry(i);
	  ski_tree->GetEntry(i);
	  jet_tree->GetEntry(i);
	 
	  //  if(i%10000 == 0)
	    {
	      //  std::cout<<"File location"<<" "<<nfile<<" "<<i<<" "<< " events running"<<std::endl;
	    }
	  	  if(vz <= -15. || vz >= 15.) continue; // vertex cut
	  
	  hEvents->AddBinContent(1,1);
	  hEvents->AddBinContent(2,1);
	  
	  
	  // cout<<"pTHat is  "<<pT_hat<<"  evtweight is"<<ptHatw<<endl;

  hEvents->AddBinContent(3,1);

  
	  if(HBHENoiseFilterResultRun2Loose != 1 || pPAprimaryVertexFilter != 1 || pBeamScrapingFilter != 1) continue; //apply the skimmed event filters

	  	  if (phfCoincFilter !=1 || pVertexFilterCutdz1p0!=1) continue;
	  hEvents->AddBinContent(4,1);
	  
	  //  if(HLT_HIAK4CaloJet80_v1 != 1) continue; // apply jet trigger
	  
	  hEvents->AddBinContent(5,1);
	  // cout<<"Check 1"<<endl;

	  if(nref <= 0) continue; // if there is no jets in an event
	  // cout<<i<<" "<<"nref value"<<" "<<nref<<" "<<"Check 2"<<endl;

	  hEvents->AddBinContent(6,1);
	  hEvents->AddBinContent(7,1);


	  int ptHatw=1;
	  hZvtx->Fill(vz,ptHatw);
	 

  

	 

	  //  cout<<"Check 3"<<endl;
	  //   std::vector<std::tuple<Float_t, Float_t, Float_t>> GenjetTriplets;
	  std::vector<std::tuple<Float_t, Float_t, Float_t>> RecojetTriplets;

	  //  std::vector<std::tuple<Float_t, Float_t, Float_t>>* GenjetTriplets= new  std::vector<std::tuple<Float_t, Float_t, Float_t>>();
	  //std::vector<std::tuple<Float_t, Float_t, Float_t>>* RecojetTriplets= new std::vector<std::tuple<Float_t, Float_t, Float_t>>();

	  
	  for (int j = 0; j < nref; j++) //Jet loop start
	     {
	       //   cout<<i<<" "<<trackMax[j]/rawpt[j]<<endl;
	       if(trackMax[j]/rawpt[j] < 0.01 || rawpt[j]<20) {RecojetTriplets.push_back(std::make_tuple(0.,0.,0.));
		 //	 GenjetTriplets.push_back(std::make_tuple(0,0,0));
		 continue;} // Cut for jets for very low maxium pT track
	       if(trackMax[j]/rawpt[j] > 0.98){RecojetTriplets.push_back(std::make_tuple(0.,0.,0.));
		 // GenjetTriplets.push_back(std::make_tuple(0,0,0))
		 ;continue;} // Cut for jets where all the pT is taken by one track
	      // cout<<"Size of jt_eta"<<sizeof(jteta)<<endl;
	      // cout<<"nref is"<<nref<<endl;
	      Float_t jt_rawpt = rawpt[j];
	      Float_t jt_eta = jteta[j];
              Float_t jt_phi = jtphi[j];
	      
	      
	      //    if (jt_rawpt < 20) continue;
	      // Float_t gen_pT= genpt[j];
	      //  Float_t gen_eta=geneta[j];
	      // Float_t gen_phi=genphi[j];

	      
	      // std::cout<<genpt[j]<<std::endl;
	       JEC.SetJetPT(jt_rawpt); 
	       JEC.SetJetEta(jt_eta); 
	       JEC.SetJetPhi(jt_phi);
	      
	       Float_t jt_corr_pt = JEC.GetCorrectedPT();
	       //  cout<<jt_corr_pt<<endl;
	      //   cout<<jt_rawpt<<"  "<<jt_corr_pt<<endl;      
	      // Float_t jt_corr_pt = jtpt[j];
	      Float_t jt_WTAeta = WTAeta[j];
	      Float_t jt_WTAphi = WTAphi[j];

	      /* float jet_jes= (jt_corr_pt-gen_pT)/gen_pT;

	          if (gen_pT>25 && gen_pT<50) JES_25->Fill(jet_jes,ptHatw);
		  if (gen_pT>50 && gen_pT<100) JES_50->Fill(jet_jes,ptHatw);
		  if (gen_pT>100&& gen_pT<130) JES_100->Fill(jet_jes,ptHatw);
		  if (gen_pT>130 && gen_pT<160) JES_130->Fill(jet_jes,ptHatw);
		  if (gen_pT>160 && gen_pT<200) JES_160->Fill(jet_jes,ptHatw);
		  if (gen_pT>200)  JES_200->Fill(jet_jes,ptHatw); */

		  // cout<<"Check 5"<<endl;
		  RecojetTriplets.push_back(std::make_tuple(jt_corr_pt,jteta[j],jtphi[j]));
		  //  GenjetTriplets.push_back(std::make_tuple(genpt[j], geneta[j], genphi[j]));
		  //                  cout<<j<<" "<<gen_pT<<" "<<jt_corr_pt<<endl;
		  //  cout<<j<<" "<<std::get<0>(RecojetTriplets[j])<<"  "<<std::get<1>(RecojetTriplets[j])<<" "<<std::get<2>(RecojetTriplets[j])<<endl;
	      

		


	      if (abs(jt_eta)<5)
		{
		 
		 
                  pT_jet_reco->Fill(jt_corr_pt,ptHatw);

		  if (jt_corr_pt>25 && jt_corr_pt<50) {Eta_jet_reco_25->Fill(jt_eta,ptHatw);}
		  if (jt_corr_pt>50 && jt_corr_pt<70) {Eta_jet_reco_50->Fill(jt_eta,ptHatw);}
		  if (jt_corr_pt>70 && jt_corr_pt<90) {Eta_jet_reco_70->Fill(jt_eta,ptHatw);}
		  if (jt_corr_pt>90 && jt_corr_pt<110) {Eta_jet_reco_90->Fill(jt_eta,ptHatw);}
		  if (jt_corr_pt>110 && jt_corr_pt<130) {Eta_jet_reco_110->Fill(jt_eta,ptHatw);}
		  if (jt_corr_pt>130 ){ Eta_jet_reco_130->Fill(jt_eta,ptHatw);}





		  //  hgenpt->Fill(gen_pT,ptHatw);
		  //  hgeneta->Fill(gen_eta,ptHatw);
		  //  hgenphi->Fill(gen_phi,ptHatw);
		}

	     
	      //cout<<jt_corr_pt<<endl;    
	      // hjtCorrpt->Fill(jt_corr_pt, ptHatw);
	      hjtEta->Fill(jt_eta, ptHatw);
	      hjtPhi->Fill(jt_phi, ptHatw);
	      hjtWTAEta->Fill(jt_WTAeta, ptHatw);
	      hjtWTAPhi->Fill(jt_WTAphi, ptHatw);

	      double DEta = jt_eta - jt_WTAeta;
	      //  double DPhi = TVector2::Phi_mpi_pi(jt_phi - jt_WTAphi);
	      // double DR = TMath::Sqrt(pow(DEta,2) + pow(DPhi,2));

	      // hDeltaR_axis->Fill(DR, ptHatw);
	      // hDeltaEta_axis->Fill(DEta, ptHatw);
	      // hDeltaPhi_axis->Fill(DPhi, ptHatw);

	    } // jet loop end
	  // cout<<"Check 6"<<endl;
	  //   cout<<std::get<0>(RecojetTriplets[0])<<"  "<<std::get<0>(RecojetTriplets[1])<<" "<<std::get<0>(RecojetTriplets[2])<<endl;
	  

	  // Gen-jets ****************************************************************************************************************************************************
	  // Leading and Subleading jet criteria
	  //pT1>30, pT2>20, -3<eta<3, Delta phi>2 pi/3
	  // 	  std::sort(GenjetTriplets.begin(), GenjetTriplets.end(), std::greater<std::tuple<Float_t, Float_t, Float_t>>());
	  std::sort(RecojetTriplets.begin(), RecojetTriplets.end(), std::greater<std::tuple<Float_t, Float_t, Float_t>>());
	  
		if (std::get<0>(RecojetTriplets[0]) > 30  &&  std::get<0>(RecojetTriplets[1]) > 20         &&         std::abs(std::get<2>(RecojetTriplets[0]) - std::get<2>(RecojetTriplets[1])) >( 2 * TMath::Pi() / 3))
		  {



		    Float_t Leading_pT=std::get<0>(RecojetTriplets[0]);
		    Float_t Subleading_pT=std::get<0>(RecojetTriplets[1]);
		    Float_t Leading_Eta=std::get<1>(RecojetTriplets[0]);
		    Float_t Subleading_Eta=std::get<1>(RecojetTriplets[1]);
		    // cout<<Leading_pT<<" "<<Subleading_pT<<" "<<Leading_Eta<<" "<<Subleading_Eta<<" "<<std::abs(std::get<2>(jetTriplets[0]) - std::get<2>(jetTriplets[1]))<<endl;
		    Float_t pTdijet=(Leading_pT+Subleading_pT)/2;
		    Float_t Etadijet=(Leading_Eta+Subleading_Eta)/2;
		    Float_t Phidijet=std::abs(std::get<2>(RecojetTriplets[0]) - std::get<2>(RecojetTriplets[1]));
		    Float_t x=Subleading_pT/Leading_pT;
		    //              std::cout<<Leading_Eta<<" "<<Subleading_Eta<<" "<<Etadijet<<std::endl;

		    if (abs(Etadijet)<5)
		      {

			Eta_leading_reco->Fill(Leading_Eta,ptHatw);
			Eta_Subleading_reco->Fill(Subleading_Eta,ptHatw);
			pT_leading_reco->Fill(Leading_pT,ptHatw);
			pT_Subleading_reco->Fill(Subleading_pT,ptHatw);
			Eta_dijet_reco->Fill(Etadijet,ptHatw);
			pT_dijet_reco->Fill(pTdijet,ptHatw);
			Phi_dijet_reco->Fill(Phidijet,ptHatw);
			X_dijet_reco->Fill(x,ptHatw);

			//Fill histograms for different dijet pT range
			if (pTdijet>25 && pTdijet<50)    Eta_dijet_reco_25->Fill(Etadijet,ptHatw);

			if (pTdijet>50 && pTdijet<70)  Eta_dijet_reco_50->Fill(Etadijet,ptHatw);

			if (pTdijet>70 && pTdijet<90) Eta_dijet_reco_70->Fill(Etadijet,ptHatw);

			if (pTdijet>90 && pTdijet<110) Eta_dijet_reco_90->Fill(Etadijet,ptHatw);

			if (pTdijet>110 && pTdijet<130) Eta_dijet_reco_110->Fill(Etadijet,ptHatw);

	 		if (pTdijet>130 ) Eta_dijet_reco_130->Fill(Etadijet,ptHatw);

		      }
		    
	  	  
		      }
	  



		    





		// cout<<"Came here"<<" "<<i<<endl;
	} // events loop end
	  delete hea_tree;      
		  delete hlt_tree;      
		  delete ski_tree;
	  delete jet_tree;
		

	  f->Close();
      }// while loop ends
     
        } // file loop end

 
  TFile *fout = new TFile("/afs/cern.ch/user/v/vpant/pPb_8.16TeV_Analysis/pPb_8_16_data_results.root", "recreate");
  // cout<<"Came here 2"<<endl;
  // hjtCorrpt->Draw();
 hgenpt->Write();
 hgeneta->Write();
 // hgenphi->Write();
 Eta_leading->Write();
 Eta_Subleading->Write();
 // Phi_leading->Write();
 // Phi_Subleading->Write();
 pT_leading->Write();
 pT_Subleading->Write();
 pT_jet_corr->Write();
 Eta_dijet->Write();
 pT_dijet->Write();
 Phi_dijet->Write();
 X_dijet->Write();
 pT_leading_reco->Write();
 pT_Subleading_reco->Write();
 Eta_dijet_reco->Write();
 pT_dijet_reco->Write();
 Phi_dijet_reco->Write();
 X_dijet_reco->Write();

 pT_jet->Write();
 pT_jet_reco->Write();
 
 JES_25->Write();
 JES_50->Write();
 JES_100->Write();
 JES_130->Write();
 JES_160->Write();
 JES_200->Write();
 
 Eta_dijet_25->Write();
 Eta_dijet_50->Write();
 Eta_dijet_70->Write();
 Eta_dijet_90->Write();
 Eta_dijet_110->Write();
 Eta_dijet_130->Write();

 Eta_jet_25->Write();
 Eta_jet_50->Write();
 Eta_jet_70->Write();
 Eta_jet_90->Write();
 Eta_jet_110->Write();
 Eta_jet_130->Write();

 Eta_dijet_reco_25->Write();
 Eta_dijet_reco_50->Write();
 Eta_dijet_reco_70->Write();
 Eta_dijet_reco_90->Write();
 Eta_dijet_reco_110->Write();
 Eta_dijet_reco_130->Write();

 Eta_jet_reco_25->Write();
 Eta_jet_reco_50->Write();
 Eta_jet_reco_70->Write();
 Eta_jet_reco_90->Write();
 Eta_jet_reco_110->Write();
 Eta_jet_reco_130->Write();

 
  hEvents->Write();
  hZvtx->Write();
  if(is_MC)
    {
      hpthat->Write();
      hpthatW->Write();
    }

  hjtCorrpt->Write();
  hjtEta->Write();
  hjtPhi->Write();
  hjtWTAEta->Write();
  hjtWTAPhi->Write();
  hDeltaR_axis->Write();
  hDeltaEta_axis->Write();
  hDeltaPhi_axis->Write();

  fout->Write();
  fout->Close();
  delete fout;
  cout<<"finished";


  

      
} // program end


