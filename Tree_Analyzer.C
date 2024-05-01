//#include "call_libraries.h"  // call libraries from ROOT and C++
#include "JetCorrector.h" // reader for JEC
#include <TROOT.h>       // ROOT system header
#include <TFile.h>       // ROOT file I/O
#include <TTree.h>       // ROOT tree
#include <TH1.h>  // ROOT 1D histograms
#include <TH1D.h>
#include "TH3.h"

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


void Tree_Analyzer(TString input_file, TString outputFileName,int is_MC,Float_t pthat_value)
{
  
  int is_Pbgoing=1;
vector<string> Files;
 int Eta_correction=1;
 // cout<<"here";
 
 if (is_MC && is_Pbgoing) Files.push_back("/afs/cern.ch/user/v/vpant/pPb_8.16TeV_Analysis/MC_files/Autumn16_HI_pPb_pgoing_Unembedded_MC_L2Relative_AK4PF.txt");

 if (is_MC && !is_Pbgoing)
   {Files.push_back("/afs/cern.ch/user/v/vpant/pPb_8.16TeV_Analysis/MC_files/Autumn16_HI_pPb_Pbgoing_Unembedded_MC_L2Relative_AK4PF.txt");
     Eta_correction=-1;
   }
 
 if (!is_MC && is_Pbgoing) Files.push_back("Pbgoing_files/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt");

  if (!is_MC && !is_Pbgoing)
    {Files.push_back("Pbgoing_files/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt");
      Eta_correction=-1;
    }
  
JetCorrector JEC(Files);

TH1::SetDefaultSumw2();
TH2::SetDefaultSumw2();
TH3::SetDefaultSumw2();

  TFile* file1= new TFile("hvztx_weighting.root","READ");

   TH1D* Ratio_vztx = (TH1D*) file1->Get("Ratio");
  
  // Read the input file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
      
	if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(file_chain.c_str());}
	inputfile.close();

	// Read the trees to be added in the Chain
	TChain *hlt_tree = new TChain("hltanalysis/HltTree");
	TChain *jet_tree = new TChain("ak4PFJetAnalyzer/t");
	TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
       TChain *ski_tree = new TChain("skimanalysis/HltTree");

	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
	  cout<<"list operator is "<<*listIterator<<endl;          
	  TFile *testfile = TFile::Open(*listIterator);
		if (!testfile) cout<<"file not found";
                if(!testfile || testfile->IsZombie() || testfile->TestBit(TFile::kRecovered)) cout << "File: " << *listIterator << " failed to open" << endl;
                if(!testfile || testfile->IsZombie() || testfile->TestBit(TFile::kRecovered)) continue;
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		hlt_tree->Add(*listIterator);
		hea_tree->Add(*listIterator);
		jet_tree->Add(*listIterator);
		ski_tree->Add(*listIterator);
		//	if(is_MC){gen_tree->Add(*listIterator);}
	}
  	
	
  
  TH1D* hEvents = new TH1D("hEvents", "", 10, 0, 10);
  TH1D* hpthat = new TH1D("hpthat", "", 200, 0, 600.);
  TH1D* hpthatW = new TH1D("hpthatW", "", 200, 0, 600.);
  TH1D* hZvtx = new TH1D("hZvtx", "", 80, -20, 20);

  
  
  TH1D* pT_jet_reco=new TH1D("pT_jet_reco","",50,0,1000);
  TH1D* pT_jet_corr=new TH1D("pT_jet_corr","",20,0,1000);
  TH1D* Phi_jet_reco=new TH1D("Phi_jet_reco","",20,(2*TMath::Pi()/3),2*TMath::Pi());
  TH1D* X_jet_reco=new TH1D("X_jet_reco","",40,0.,1.);

  TH2D* hpteta_dijet = new TH2D("hpteta_dijet","",1000,0,1000,100,-5,5);
  TH2D* hpteta_jet = new TH2D("hpteta_jet","",1000,0,1000,100,-5,5);

  
  TH2D* hgenpteta_dijet = new TH2D("hgenpteta_dijet","",1000,0,1000,100,-5,5);
  TH2D* hgenpteta_jet = new TH2D("hgenpteta_jet","",1000,0,1000,100,-5,5);
    
  
  TH3D* heta_diff_pt_jetL_eta = new TH3D("heta_diff_pt_jetL_eta","",500,0,1000,1000,-0.4,0.4,100,-5,5);
  TH3D* heta_diff_pt_jetSL_eta = new TH3D("heta_diff_pt_jetSL_eta","",500,0,1000,1000,-0.4,0.4,100,-5,5);

  
  
  TH3D* hpt_diff_pt_jet_eta = new TH3D("hpt_diff_pt_jet_eta","",500,0,1000,1000,0.,2.,100,-5.,5.);
  //  TH2D* hphi_diff_pt_jet = new TH2D("hphi_diff_pt_jet","",1000,0,1000,4000,-TMath::Pi(),TMath::Pi());
  TH3D* hpt_subL_L_jet_pT=new TH3D("hpt_subL_L_jet_pT","",500,0,1000,500,0.,1000.,200,0.,1000);    
  

 TH2D* hpteta_hlt_Calo_40_jet= new TH2D("hpteta_hlt_Calo_40_jet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_Calo_60_jet= new TH2D("hpteta_hlt_Calo_60_jet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_Calo_80_jet= new TH2D("hpteta_hlt_Calo_80_jet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_Calo_100_jet= new TH2D("hpteta_hlt_Calo_100_jet","",1000,0,1000,200,-5,5);

  TH2D* hpteta_hlt_Calo_40_dijet= new TH2D("hpteta_hlt_Calo_40_dijet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_Calo_60_dijet= new TH2D("hpteta_hlt_Calo_60_dijet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_Calo_80_dijet= new TH2D("hpteta_hlt_Calo_80_dijet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_Calo_100_dijet= new TH2D("hpteta_hlt_Calo_100_dijet","",1000,0,1000,200,-5,5);
 

  TH2D* hpteta_hlt_PF_60_jet= new TH2D("hpteta_hlt_PF_60_jet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_PF_80_jet= new TH2D("hpteta_hlt_PF_80_jet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_PF_100_jet= new TH2D("hpteta_hlt_PF_100_jet","",1000,0,1000,200,-5,5);


  
  TH2D* hpteta_hlt_PF_60_dijet= new TH2D("hpteta_hlt_PF_60_dijet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_PF_80_dijet= new TH2D("hpteta_hlt_PF_80_dijet","",1000,0,1000,200,-5,5);
  TH2D* hpteta_hlt_PF_100_dijet= new TH2D("hpteta_hlt_PF_100_dijet","",1000,0,1000,200,-5,5);
  


  
   


       
       
      
        
  int nevents = hea_tree->GetEntries(); // number of events
	cout << "Total number of events in those files: "<< nevents << endl;
	cout << endl;
	cout << "-------------------------------------------------" << endl;
        
	// Start loop over events
	
	for(int i = 0; i < nevents; i++)
	  {
	    if (i%10000==0) cout<<i<<" events passed"<<endl;    
      hea_tree->SetBranchStatus("*", 0);

      Float_t vz;
      Float_t weight;
      Float_t pthat;

      hea_tree->SetBranchStatus("vz", 1);
      hea_tree->SetBranchAddress("vz", &vz);

            if(is_MC)
	{
	  hea_tree->SetBranchStatus("weight", 1);
	  hea_tree->SetBranchAddress("weight", &weight);
	  
	  hea_tree->SetBranchStatus("pthat", 1);
	  hea_tree->SetBranchAddress("pthat", &pthat);
	  } 
	    // cout<<"Reached 1"<<endl;
     
      hlt_tree->SetBranchStatus("*", 0);
      //*************************************Calo jets**********************************************
        Int_t HLT_PAAK4CaloJet40_Eta5p1_v3=0;
TBranch* branch = hlt_tree->GetBranch("HLT_PAAK4CaloJet40_Eta5p1_v3");

 if(branch){
        hlt_tree->SetBranchStatus("HLT_PAAK4CaloJet40_Eta5p1_v3", 1);
        hlt_tree->SetBranchAddress("HLT_PAAK4CaloJet40_Eta5p1_v3", &HLT_PAAK4CaloJet40_Eta5p1_v3);
 }

	Int_t HLT_PAAK4CaloJet60_Eta5p1_v3;

        hlt_tree->SetBranchStatus("HLT_PAAK4CaloJet60_Eta5p1_v3", 1);
        hlt_tree->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v3", &HLT_PAAK4CaloJet60_Eta5p1_v3);


	Int_t HLT_PAAK4CaloJet80_Eta5p1_v3;

        hlt_tree->SetBranchStatus("HLT_PAAK4CaloJet80_Eta5p1_v3", 1);
        hlt_tree->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v3", &HLT_PAAK4CaloJet80_Eta5p1_v3);


	Int_t HLT_PAAK4CaloJet100_Eta5p1_v3;

        hlt_tree->SetBranchStatus("HLT_PAAK4CaloJet100_Eta5p1_v3", 1);
        hlt_tree->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v3", &HLT_PAAK4CaloJet100_Eta5p1_v3);

	//*********************************** PF jets ************************************************8

	Int_t HLT_PAAK4PFJet60_Eta5p1_v4;

        hlt_tree->SetBranchStatus("HLT_PAAK4PFJet60_Eta5p1_v4", 1);
        hlt_tree->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4", &HLT_PAAK4PFJet60_Eta5p1_v4);


        Int_t HLT_PAAK4PFJet80_Eta5p1_v3;

        hlt_tree->SetBranchStatus("HLT_PAAK4PFJet80_Eta5p1_v3", 1);
        hlt_tree->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3", &HLT_PAAK4PFJet80_Eta5p1_v3);


        Int_t HLT_PAAK4PFJet100_Eta5p1_v3;

        hlt_tree->SetBranchStatus("HLT_PAAK4PFJet100_Eta5p1_v3", 1);
        hlt_tree->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3", &HLT_PAAK4PFJet100_Eta5p1_v3);


	
    
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
     
      jet_tree->SetBranchStatus("*", 0);

      const int nmax = 1000;
      Int_t nref=0;
      Int_t ngen=0;
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
      Float_t refeta[nmax];
      Float_t refpt[nmax];
      Float_t refphi[nmax];
      
      Float_t NHF[nmax];  
      Float_t NEMF[nmax]; 
      Float_t CHF[nmax];  
      Float_t MUF[nmax];  
      Float_t CEMF[nmax];
      Int_t CHM[nmax]; 
      Int_t CM[nmax]; 
      Int_t NHM[nmax];
      Int_t NM[nmax];
      Int_t MM[nmax];
      
           jet_tree->SetBranchStatus("jtPfNHF", 1);
      jet_tree->SetBranchAddress("jtPfNHF", &NHF);    

      jet_tree->SetBranchStatus("jtPfNEF", 1);
      jet_tree->SetBranchAddress("jtPfNEF", &NEMF);

      jet_tree->SetBranchStatus("jtPfCHF", 1);
      jet_tree->SetBranchAddress("jtPfCHF", &CHF);

      jet_tree->SetBranchStatus("jtPfMUF", 1);
      jet_tree->SetBranchAddress("jtPfMUF", &MUF);

      jet_tree->SetBranchStatus("jtPfCEF", 1);
      jet_tree->SetBranchAddress("jtPfCEF", &CEMF);

      jet_tree->SetBranchStatus("jtPfCHM", 1);
      jet_tree->SetBranchAddress("jtPfCHM", &CHM);

       jet_tree->SetBranchStatus("jtPfCEM", 1);
      jet_tree->SetBranchAddress("jtPfCEM", &CM);

         jet_tree->SetBranchStatus("jtPfNHM", 1);
      jet_tree->SetBranchAddress("jtPfNHM", &NHM);
      
       jet_tree->SetBranchStatus("jtPfNEM", 1);
      jet_tree->SetBranchAddress("jtPfNEM", &NM);
      
      jet_tree->SetBranchStatus("jtPfMUM", 1);
      jet_tree->SetBranchAddress("jtPfMUM", &MM);
       
      
      jet_tree->SetBranchStatus("nref", 1);
      jet_tree->SetBranchAddress("nref", &nref);

      jet_tree->SetBranchStatus("rawpt", 1);
      jet_tree->SetBranchAddress("rawpt", &rawpt);

      jet_tree->SetBranchStatus("trackMax", 1);
      jet_tree->SetBranchAddress("trackMax", &trackMax);

      if (is_MC)
	{
      jet_tree->SetBranchStatus("genpt",1);
      jet_tree->SetBranchAddress("genpt",&genpt);

      jet_tree->SetBranchStatus("geneta",1);
      jet_tree->SetBranchAddress("geneta",&geneta);

      jet_tree->SetBranchStatus("genphi",1);
      jet_tree->SetBranchAddress("genphi",&genphi);

      jet_tree->SetBranchStatus("refeta", 1);
      jet_tree->SetBranchAddress("refeta", &refeta);

       jet_tree->SetBranchStatus("refpt", 1);
      jet_tree->SetBranchAddress("refpt", &refpt);

      jet_tree->SetBranchStatus("ngen", 1);
      jet_tree->SetBranchAddress("ngen", &ngen);

      jet_tree->SetBranchStatus("refphi",1);
      jet_tree->SetBranchAddress("refphi",&refphi);

	}
      // cout<<"Reached 2"<<endl;
      jet_tree->SetBranchStatus("jteta", 1);
      jet_tree->SetBranchAddress("jteta", &jteta);

      jet_tree->SetBranchStatus("jtphi", 1);
      jet_tree->SetBranchAddress("jtphi", &jtphi);

      jet_tree->SetBranchStatus("WTAeta", 1);
      jet_tree->SetBranchAddress("WTAeta", &WTAeta);

      jet_tree->SetBranchStatus("WTAphi", 1);
      jet_tree->SetBranchAddress("WTAphi", &WTAphi);
      
      
      
	  hea_tree->GetEntry(i);
	  hlt_tree->GetEntry(i);
	  ski_tree->GetEntry(i);
	  jet_tree->GetEntry(i);
	 
	  //  if(i%10000 == 0)
	    {
	      //  std::cout<<"File location"<<" "<<nfile<<" "<<i<<" "<< " events running"<<std::endl;
	    }
	  	  if(vz <= -15. || vz >= 15.) continue; // vertex cut
	          hEvents->Fill(1);

  
	  if(HBHENoiseFilterResultRun2Loose != 1 || pPAprimaryVertexFilter != 1 || pBeamScrapingFilter != 1) continue; //apply the skimmed event filters

	  	  if (phfCoincFilter !=1 || pVertexFilterCutdz1p0!=1) continue;
		  // hEvents->AddBinContent(4,1);
	  
	  //  if(HLT_HIAK4CaloJet80_v1 != 1) continue; // apply jet trigger
	  
		  // hEvents->AddBinContent(5,1);
	  // cout<<"Check 1"eck<<endl;

	  if(nref <= 0) continue; // if there is no jets in an event
	  // cout<<i<<" "<<"nref value"<<" "<<nref<<" "<<"Check 2"<<endl;

	  // hEvents->AddBinContent(6,1);
	  // hEvents->AddBinContent(7,1);


	  Float_t ptHatw=1.0;
	  float pT_hat=0.;

	   if(is_MC)
            {  pT_hat=pthat;
	      
              float evtweight=1.0;
              if(pT_hat > 15.0 && pT_hat <= 30.){evtweight = 1.0404701e-06 * 961104;}
              else if(pT_hat > 30. && pT_hat <= 50.){evtweight = 7.7966624e-08 * 952110;}
              else if(pT_hat > 50. && pT_hat <= 80.){evtweight = 1.0016052e-08 * 952554;}
              else if(pT_hat > 80. && pT_hat <= 120.){evtweight = 1.3018269e-09 * 996844;}
              else if(pT_hat > 120. && pT_hat <= 170.){evtweight = 2.2648493e-10 * 964681;}
              else if(pT_hat > 170. && pT_hat <= 220.){evtweight = 4.0879112e-11 * 999260;}
              else if(pT_hat > 220. && pT_hat <= 280.){evtweight = 1.1898939e-11 * 964336;}
              else if(pT_hat > 280. && pT_hat <= 370.){evtweight = 3.3364433e-12 * 995036;}
              else if(pT_hat > 370. && pT_hat <= 460.){evtweight = 7.6612402e-13 * 958160;}
              else if(pT_hat > 460. && pT_hat <= 540.){evtweight = 2.1341026e-13 * 981427;}
              else if(pT_hat > 540.){evtweight = 7.9191586e-14 * 1000000;}

              ptHatw=evtweight/nevents;
	        Int_t bin = Ratio_vztx->FindBin(vz);
	        Double_t ratio_vztx = Ratio_vztx->GetBinContent(bin);
	       ptHatw*=ratio_vztx;
          if (pthat_value==15.0 && !(pT_hat>15.0 && pT_hat<=30.)) continue;
          if (pthat_value==30.0 && !(pT_hat>30.0 && pT_hat<=50.)) continue;
          if (pthat_value==50.0 && !(pT_hat>50.0 && pT_hat<=80.)) continue;
          if (pthat_value==80.0 && !(pT_hat>80.0 && pT_hat<=120.)) continue;
          if (pthat_value==120.0 && !(pT_hat>120.0 && pT_hat<=170.)) continue;
          if (pthat_value==170.0 && !(pT_hat>170.0 && pT_hat<=220.)) continue;
          if (pthat_value==220.0 && !(pT_hat>220.0 && pT_hat<=280.)) continue;
          if (pthat_value==280.0 && !(pT_hat>280.0 && pT_hat<=370.)) continue;

          if (pthat_value==370.0 && !(pT_hat>370.0 && pT_hat<=460.)) continue;

          if (pthat_value==460.0 && !(pT_hat>460.0 && pT_hat<=540.)) continue;

          if (pthat_value==540.0 && !(pT_hat>540.)) continue;
	    }
	   // cout<<"Reached 3"<<endl;

if(is_MC)
  {
              hpthat->Fill(pT_hat);
              hpthatW->Fill(pT_hat,ptHatw);
            }
	   hZvtx->Fill(vz,ptHatw);
	 
	   //           cout<<ptHatw<<endl;
  
	   //	   cout<<ptHatw<<endl;
	 

	   // cout<<"Check 3"<<endl;
	  //   std::vector<std::tuple<Float_t, Float_t, Float_t>> GenjetTriplets;
	   std::vector<std::tuple<Float_t, Float_t, Float_t,Float_t, Float_t, Float_t>> RecojetTriplets;
          std::vector<std::tuple<Float_t, Float_t, Float_t>> GenjetTriplets;
	  

	  RecojetTriplets.push_back(std::make_tuple(0.,0.,0.,0.,0.,0.));
                 GenjetTriplets.push_back(std::make_tuple(0.,0.,0.));

	int jetid=0;
	int Eta_corr=Eta_correction;
	  for (int j = 0; j < nref; j++)  // Reco Jet loop start
	    {

	       Int_t NumConst = CM[j]+NM[j];
               Int_t NumNeutralParticle =NM[j];

               
	      
	       //if(abs(jteta[j])<=2.4 && trackMax[j]/rawpt[j] < 0.01) {
	       // continue;} // Cut for jets for very low maxium pT track
	      
	       //if(trackMax[j]/rawpt[j] > 0.98){
	       //continue;} // Cut for jets where all the pT is taken by one track

	      Float_t jt_rawpt = rawpt[j];
	      Float_t jt_eta = jteta[j];
              Float_t jt_phi = jtphi[j];
              
	      
	       if (abs(jt_eta)<=2.7)  jetid=((NHF[j]<0.90 && NEMF[j]<0.90 && NumConst>1 && MUF[j]<0.8) && ((abs(jt_eta)<=2.4 && CHF[j]>0 && CHM[j]>0 && CEMF[j]<0.90) || abs(jt_eta)>2.4) && abs(jt_eta)<=2.7);
	       
   
	      if (jt_eta>2.7 &&	jt_eta<=3) jetid=((NHF[j]<0.98 && NEMF[j]>0.01 && NumNeutralParticle>2 && abs(jt_eta)>2.7 && abs(jt_eta)<=3.0 ));

	      //	      if (abs(jt_eta)>3) jetid=((NEMF[j]<0.90 && NumNeutralParticle>10 && abs(jt_eta)>3.0 ));
	         
	      
			    if (!jetid) continue;
	    
	       JEC.SetJetPT(jt_rawpt); 
	       JEC.SetJetEta(jt_eta); 
	       JEC.SetJetPhi(jt_phi);
              Float_t jt_corr_pt = JEC.GetCorrectedPT();	       

	      jt_eta*=Eta_corr;
	      
	       Float_t ref_eta;
	       Float_t ref_pt;
	       Float_t ref_phi;

	       
	      if (is_MC)
                {
		  ref_phi=refphi[j];
                ref_eta=refeta[j];
                ref_pt=refpt[j];
		ref_eta*=Eta_corr;
		RecojetTriplets.push_back(std::make_tuple(jt_corr_pt,jteta[j],jtphi[j],ref_pt,ref_eta,ref_phi));
                }
	      else	  RecojetTriplets.push_back(std::make_tuple(jt_corr_pt,jt_eta,jt_phi,0.,0.,0.));
		  
		  //	   cout<<"Reached 4 a "<<endl;
	       if (abs(jt_eta)<5)
		{
		 
		 
                  
                  hpteta_jet->Fill(jt_corr_pt,jt_eta,ptHatw);

		  if (HLT_PAAK4CaloJet40_Eta5p1_v3) hpteta_hlt_Calo_40_jet->Fill(jt_corr_pt,jt_eta,ptHatw);
                  if (HLT_PAAK4CaloJet60_Eta5p1_v3) hpteta_hlt_Calo_60_jet->Fill(jt_corr_pt,jt_eta,ptHatw);
                  if (HLT_PAAK4CaloJet80_Eta5p1_v3) hpteta_hlt_Calo_80_jet->Fill(jt_corr_pt,jt_eta,ptHatw);
                  if (HLT_PAAK4CaloJet100_Eta5p1_v3) hpteta_hlt_Calo_100_jet->Fill(jt_corr_pt,jt_eta,ptHatw);

                  if (HLT_PAAK4PFJet60_Eta5p1_v4) hpteta_hlt_PF_60_jet->Fill(jt_corr_pt,jt_eta,ptHatw);
      	      	  if (HLT_PAAK4PFJet80_Eta5p1_v3) hpteta_hlt_PF_80_jet->Fill(jt_corr_pt,jt_eta,ptHatw);
                  if (HLT_PAAK4PFJet100_Eta5p1_v3) hpteta_hlt_PF_100_jet->Fill(jt_corr_pt,jt_eta,ptHatw);

                   if (is_MC)
                    { Float_t eta_diff=jt_eta-ref_eta;
                      Float_t pt_diff=jt_corr_pt/ref_pt;
                      Float_t phi_diff=TMath::ACos(TMath::Cos(jt_phi))- TMath::ACos(TMath::Cos(ref_phi));
		      //   cout<<phi_diff<<endl;
		   
		      // if (ref_pt>0) heta_diff_pt_jet->Fill(jt_corr_pt,eta_diff,ptHatw);
                      if (ref_pt>0) hpt_diff_pt_jet_eta->Fill(ref_pt,pt_diff,ref_eta,ptHatw);
		      //  if (ref_pt>0) hphi_diff_pt_jet->Fill(jt_corr_pt,phi_diff,ptHatw);
		        
		      } 

		   //cout<<"Reached 4"<<endl;
		}

	   
	    } // reco jet loop end





	  	  for (int j = 0; j < ngen; j++) //gen Jet loop start
	    {
	      /*	      if (abs(jt_eta)<=2.7)  jetid=((NHF[j]<0.90 && NEMF[j]<0.90 && NumConst>1 && MUF[j]<0.8) && ((abs(jt_eta)<=2.4 && CHF[j]> \
0 && CHM[j]>0 && CEMF[j]<0.90) || abs(jt_eta)>2.4) && abs(jt_eta)<=2.7);


              if (jt_eta>2.7 && jt_eta<=3) jetid=((NHF[j]<0.98 && NEMF[j]>0.01 && NumNeutralParticle>2 && abs(jt_eta)>2.7 && abs(jt_eta)<=3.0 ));

              if (abs(jt_eta)>3) jetid=((NEMF[j]<0.90 && NumNeutralParticle>10 && abs(jt_eta)>3.0 ));


	      if (!jetid) continue;
	      */
	      Float_t gen_pT;
                Float_t gen_eta;
               Float_t gen_phi;

	       
	      if (is_MC)
		{
	       gen_pT= genpt[j];
	    
	        gen_eta=geneta[j];
	        gen_phi=genphi[j];
		gen_eta*=Eta_corr;
		}

	      
	      // if (gen_pT<30) continue;
	      if (abs(gen_eta)<5)
		{
		 
		 
                  
              

		  if (is_MC)
		    { hgenpteta_jet->Fill(gen_pT,gen_eta,ptHatw);
                      		      
          	      GenjetTriplets.push_back(std::make_tuple(genpt[j], geneta[j], genphi[j]));
		      //   cout<<eta_diff<<"  "<<pt_diff<<endl;   
		    
		      // cout<<"Reached 5"<<endl;
		}
		}
	   
	    } //gen jet loop end
	  

		  // cout<<"Check 6"<<endl;
	  //   cout<<std::get<0>(RecojetTriplets[0])<<"  "<<std::get<0>(RecojetTriplets[1])<<" "<<std::get<0>(RecojetTriplets[2])<<endl;
	  

	  // Gen-jets ****************************************************************************************************************************************************
	  // Leading and Subleading jet criteria
	  //pT1>30, pT2>20, -3<eta<3, Delta phi>2 pi/3
	  if (is_MC)
	    {std::sort(GenjetTriplets.begin(), GenjetTriplets.end(), std::greater<std::tuple<Float_t, Float_t, Float_t>>());
	      // cout<<std::get<0>(GenjetTriplets[0])<<"  "<<std::get<0>(GenjetTriplets[1])<<"   "<<TMath::ACos(TMath::Cos((std::get<2>(GenjetTriplets[0]) - std::get<2>(GenjetTriplets[1]))))<<endl;
	      if (std::get<0>(GenjetTriplets[0]) > 50  &&  std::get<0>(GenjetTriplets[1]) > 30   && TMath::ACos(TMath::Cos((std::get<2>(GenjetTriplets[0]) - std::get<2>(GenjetTriplets[1]))) >( 2 * TMath::Pi() / 3)))
                    {



                      Float_t Leading_pT=std::get<0>(GenjetTriplets[0]);
            Float_t Subleading_pT=std::get<0>(GenjetTriplets[1]);
            Float_t Leading_Eta=std::get<1>(GenjetTriplets[0]);
            Float_t Subleading_Eta=std::get<1>(GenjetTriplets[1]);
            // cout<<Leading_pT<<" "<<Subleading_pT<<" "<<Leading_Eta<<" "<<Subleading_Eta<<" "<<std::abs(std::get<2>(jetTriplets[0]) - std::get<2>(j\
etTriplets[1]))<<endl;                                                                                                                                
                Float_t pTdijet=(Leading_pT+Subleading_pT)/2;
                Float_t Etadijet=(Leading_Eta+Subleading_Eta)/2;
                Float_t Phidijet=TMath::ACos(TMath::Cos((std::get<2>(GenjetTriplets[0]) - std::get<2>(GenjetTriplets[1]))));
                Float_t x=Subleading_pT/Leading_pT;


                if (abs(Etadijet)<5)
                {
hgenpteta_dijet->Fill(pTdijet,Etadijet,ptHatw);



	    }
		    }
	    }
	  //  cout<<"Reached 5"<<endl;
	      // if (RecojetTriplets.size()<2) continue;

	  //	  cout<<std::get<0>(RecojetTriplets[1])<<" "<<std::get<3>(RecojetTriplets[1])<<" ";
	  std::sort(RecojetTriplets.begin(), RecojetTriplets.end(), std::greater<std::tuple<Float_t, Float_t, Float_t,Float_t,Float_t,Float_t>>());
	  // cout<<"Out of jet loop"<<endl;
	  // cout<<std::get<0>(RecojetTriplets[0])<<"  "<<std::get<0>(RecojetTriplets[1])<<endl;

	  
	  if (std::get<0>(RecojetTriplets[0]) > 50  &&  std::get<0>(RecojetTriplets[1]) > 30         &&         TMath::ACos(TMath::Cos((std::get<2>(RecojetTriplets[0]) - std::get<2>(RecojetTriplets[1]))) >( 2 * TMath::Pi() / 3)))
		  {

		    //               cout<<"Inside comparison statement"<<endl;

		    Float_t Leading_pT=std::get<0>(RecojetTriplets[0]);
		    Float_t Subleading_pT=std::get<0>(RecojetTriplets[1]);

		    Float_t Leading_pT_ref=std::get<3>(RecojetTriplets[0]);
                    Float_t Subleading_pT_ref=std::get<3>(RecojetTriplets[1]);

		    
		    
		    Float_t Leading_Eta=std::get<1>(RecojetTriplets[0]);
		    Float_t Subleading_Eta=std::get<1>(RecojetTriplets[1]);

		    Float_t Leading_Eta_ref=std::get<4>(RecojetTriplets[0]);
                    Float_t Subleading_Eta_ref=std::get<4>(RecojetTriplets[1]);

		    Float_t Eta_diff_L=Leading_Eta - Leading_Eta_ref;
		    Float_t Eta_diff_SL=Subleading_Eta - Subleading_Eta_ref;

		    
		    
		    Float_t pTdijet=(Leading_pT+Subleading_pT)/2;
		    Float_t pTdijet_ref=(Leading_pT_ref+Subleading_pT_ref)/2;
		    
		    Float_t Etadijet=(Leading_Eta+Subleading_Eta)/2;
		    Float_t Phidijet=TMath::ACos(TMath::Cos((std::get<2>(RecojetTriplets[0]) - std::get<2>(RecojetTriplets[1]))));
		    Float_t x=Subleading_pT/Leading_pT;
		    //	    cout<<std::get<0>(RecojetTriplets[0])<<" "<<std::get<3>(RecojetTriplets[0])<<endl;
		    if (is_MC)
		      {   if (std::get<3>(RecojetTriplets[0]) > 0)  heta_diff_pt_jetL_eta->Fill(pTdijet_ref,Eta_diff_L,Leading_Eta_ref,ptHatw);
			// 	cout<<pTdijet<<" "<<Eta_diff_L<<" "<<Leading_Eta_ref<<endl;
			if (std::get<3>(RecojetTriplets[1]) >0)	heta_diff_pt_jetSL_eta->Fill(pTdijet_ref,Eta_diff_SL,Subleading_Eta_ref,ptHatw);
		      }
		    if (abs(Etadijet)<5)
		      {
                        hpt_subL_L_jet_pT->Fill(Leading_pT,Subleading_pT,pTdijet,ptHatw);
			hpteta_dijet->Fill(pTdijet,Etadijet,ptHatw);
		      

		  if (HLT_PAAK4CaloJet40_Eta5p1_v3) hpteta_hlt_Calo_40_dijet->Fill(pTdijet,Etadijet,ptHatw);
	  	  if (HLT_PAAK4CaloJet60_Eta5p1_v3) hpteta_hlt_Calo_60_dijet->Fill(pTdijet,Etadijet,ptHatw);
		  if (HLT_PAAK4CaloJet80_Eta5p1_v3) hpteta_hlt_Calo_80_dijet->Fill(pTdijet,Etadijet,ptHatw);
		  if (HLT_PAAK4CaloJet100_Eta5p1_v3) hpteta_hlt_Calo_100_dijet->Fill(pTdijet,Etadijet,ptHatw);

		  if (HLT_PAAK4PFJet60_Eta5p1_v4) hpteta_hlt_PF_60_dijet->Fill(pTdijet,Etadijet,ptHatw);
                  if (HLT_PAAK4PFJet80_Eta5p1_v3) hpteta_hlt_PF_80_dijet->Fill(pTdijet,Etadijet,ptHatw);
		  if (HLT_PAAK4PFJet100_Eta5p1_v3) hpteta_hlt_PF_100_dijet->Fill(pTdijet,Etadijet,ptHatw);
		      }
		  }
	  



		    





	  // cout<<"Came here"<<" "<<i<<endl;
	} // events loop end
	  delete hea_tree;      
		  delete hlt_tree;      
		  delete ski_tree;
	  delete jet_tree;
		

	  
     
 

	  TFile *fout = new TFile(outputFileName, "RECREATE");
        
 
  // cout<<"Came here 2"<<endl;
  // hjtCorrpt->Draw();
 
 

  
 hpteta_dijet->Write();
 hpteta_jet->Write();
 hpt_subL_L_jet_pT->Write();
 if(is_MC)
   { hgenpteta_dijet->Write();
     hgenpteta_jet->Write();
         heta_diff_pt_jetL_eta->Write();
       heta_diff_pt_jetSL_eta->Write();
       hpt_diff_pt_jet_eta->Write();
   }
 
 hpteta_hlt_Calo_40_jet->Write();
 hpteta_hlt_Calo_60_jet->Write();
 hpteta_hlt_Calo_80_jet->Write();
 hpteta_hlt_Calo_100_jet->Write();

 hpteta_hlt_Calo_40_dijet->Write();
 hpteta_hlt_Calo_60_dijet->Write();
 hpteta_hlt_Calo_80_dijet->Write();
 hpteta_hlt_Calo_100_dijet->Write();
 
 // hpteta_hlt_PF_40_jet->Write();
 hpteta_hlt_PF_60_jet->Write();
 hpteta_hlt_PF_80_jet->Write();
 hpteta_hlt_PF_100_jet->Write();

 // hpteta_hlt_PF_40_dijet->Write();
 hpteta_hlt_PF_60_dijet->Write();
 hpteta_hlt_PF_80_dijet->Write();
 hpteta_hlt_PF_100_dijet->Write();


 hEvents->Write();
  hZvtx->Write();
   if(is_MC)
    {
      hpthat->Write();
      hpthatW->Write();
    }

  
  fout->Write();
  fout->Close();
  delete fout;
  cout<<"finished";


 

      
} // program end

int main(int argc, char** argv){
                                TString firstArgument(argv[1]);
                                TString outfile(argv[2]);
				int mc = atoi(argv[3]);
                                Float_t pthat_value= atoi(argv[4]);
				//		int is_Pbgoing=atoi(argv[5]);
				//cout<<is_Pbgoing;
				Tree_Analyzer(firstArgument,outfile,mc,pthat_value);
}


