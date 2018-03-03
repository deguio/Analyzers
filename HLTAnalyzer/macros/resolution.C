//calculate rate vs cut for scouting triggers
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

void setGlobalStyle()
{
  // For the statistics box:                                                                                                                                                                      
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");                                                                                                               
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  tdrStyle->SetPadGridX(true);
  tdrStyle->SetPadGridY(true);

  tdrStyle->SetMarkerSize(0.5);
}

inline
float f_deltaRR (float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = std::abs(phi1-phi2); if (dphi>float(M_PI)) dphi-=float(2*M_PI);  
    return deta*deta + dphi*dphi;
  }
inline
float f_deltaPhi (float phi1, float phi2)
{
  float dphi = std::abs(phi1-phi2); if (dphi>float(M_PI)) dphi-=float(2*M_PI);
  return dphi;
}

/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction)
{
  float integralMax = fraction * histo->Integral();
  
  int N = histo -> GetNbinsX();
  int M1 = 0;
  int M2 = 0;
  for(int bin1 = 0; bin1 < N; ++bin1)
    {
      if( histo->GetBinContent(bin1+1) > 0. && M1 == 0 ) M1 = bin1-1;
      if( histo->GetBinContent(bin1+1) > 0. ) M2 = bin1+2;
    }
  
  std::map<int,float> binCenters;
  std::map<int,float> binContents;
  std::map<int,float> binIntegrals;
  for(int bin1 = M1; bin1 < M2; ++bin1)
    {
      binCenters[bin1] = histo->GetBinCenter(bin1+1);
      binContents[bin1] = histo->GetBinContent(bin1+1);
    
      for(int bin2 = M1; bin2 <= bin1; ++bin2)
	binIntegrals[bin1] += binContents[bin2];
    }
  
  float min = 0.;
  float max = 0.;
  float delta = 999999.;
  for(int bin1 = M1; bin1 < M2; ++bin1)
    {
      for(int bin2 = bin1+1; bin2 < M2; ++bin2)
	{
	  if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;
      
	  float tmpMin = histo -> GetBinCenter(bin1+1);
	  float tmpMax = histo -> GetBinCenter(bin2+1);
      
	  if( (tmpMax-tmpMin) < delta )
	    {
	      delta = (tmpMax - tmpMin);
	      min = tmpMin;
	      max = tmpMax;
	    }
      
	  break;
	}
    }
  
  TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
  for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
    {
      if( smallHisto->GetBinCenter(bin) < min )
	smallHisto -> SetBinContent(bin,0);
    
      if( smallHisto->GetBinCenter(bin) > max )
	smallHisto -> SetBinContent(bin,0);
    }
  smallHisto -> SetFillColor(kYellow);
  //smallHisto->Write();
  
  float mean = smallHisto -> GetMean();
  float meanErr = smallHisto -> GetMeanError();  
  
  ret[0] = mean;
  ret[1] = meanErr;
  ret[2] = min;
  ret[3] = max;
}

int resolution(std::string fileName, int thr, std::string refmet)
{
  setGlobalStyle();
  std::cout << "USING THRESHOLD: " << thr << std::endl;
  TFile* outFile = new TFile(("OUT_"+fileName+".root").c_str(),"RECREATE");
  
  
  //###############################
  unsigned int PFTAU50MET90 = 0;
  unsigned int PFMET120PFMHT120CALOBTAG = 1;
  unsigned int PFJET80PFMET120PFMHT120 = 2;
  unsigned int PFMET120PFMHT120  = 3;
  unsigned int PFMET100PFMHT100  = 4;

  bool matched = true;
  //###############################

  TChain* tt = new TChain("MyAnalysis/HLTree");
  tt->Add((fileName+".root").c_str());

  //set branches
  TBranch* b_lumi;

  //offline
  TBranch* b_caloMET;
  TBranch* b_caloMETPhi;
  TBranch* b_PFMET;
  TBranch* b_PFMETNoMu;
  //online
  TBranch* b_hltMetClean;
  TBranch* b_hltMetCleanNoMu;
  TBranch* b_hltMetCleanPhi;
  TBranch* b_hltMht;
  TBranch* b_hltPFMET;
  TBranch* b_hltPFMHTTightID;
  TBranch* b_hltPFMHTNoMuTightID;
  TBranch* b_hltPFMETNoMu;

  TBranch* b_hltAccept;

  int lumi_ = 0;
  //offline
  float caloMET_ = 0;
  float caloMETPhi_ = 0;
  float PFMET_ = 0;
  float PFMETNoMu_ = 0;
  //online
  float hltMetClean_ = 0;
  float hltMetCleanNoMu_ = 0;
  float hltMetCleanPhi_ = 0;
  float hltMht_ = 0;
  float hltPFMET_ = 0;
  float hltPFMHTTightID_ = 0;
  float hltPFMHTNoMuTightID_ = 0;
  float hltPFMETNoMu_ = 0;

  std::vector<int>* hltAccept_ = 0;

  tt->SetBranchAddress("lumi", &lumi_, &b_lumi);
  //offline
  tt->SetBranchAddress("caloMet", &caloMET_, &b_caloMET);
  tt->SetBranchAddress("caloMetPhi", &caloMETPhi_, &b_caloMETPhi);
  tt->SetBranchAddress("PFMet", &PFMET_, &b_PFMET);
  tt->SetBranchAddress("PFMetNoMu", &PFMETNoMu_, &b_PFMETNoMu);
  //online
  tt->SetBranchAddress("hltMetClean", &hltMetClean_, &b_hltMetClean);
  tt->SetBranchAddress("hltMetCleanNoMu", &hltMetCleanNoMu_, &b_hltMetCleanNoMu);
  tt->SetBranchAddress("hltMetCleanPhi", &hltMetCleanPhi_, &b_hltMetCleanPhi);
  tt->SetBranchAddress("hltMht", &hltMht_, &b_hltMht);
  tt->SetBranchAddress("hltPFMET", &hltPFMET_, &b_hltPFMET);
  tt->SetBranchAddress("hltPFMHTTightID", &hltPFMHTTightID_, &b_hltPFMHTTightID);
  tt->SetBranchAddress("hltPFMHTNoMuTightID", &hltPFMHTNoMuTightID_, &b_hltPFMHTNoMuTightID);
  tt->SetBranchAddress("hltPFMETNoMu", &hltPFMETNoMu_, &b_hltPFMETNoMu);

  tt->SetBranchAddress("hltAccept", &hltAccept_, &b_hltAccept);

  int nentries = tt->GetEntries();
  std::cout << "Number of entries: " << nentries << std::endl;

  //book graphs and plots
  Double_t binsx[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,150,250,500};
  //Double_t binsx[] = {0,5};
  Int_t  binnumx = sizeof(binsx)/sizeof(Double_t) - 1;

  float min = -1.5;
  float max = 8.;
  int nBins = 1000;
  Double_t binsy[nBins];
  int pos=0;
  for(float ii=min;ii<max;ii+=(max-min)/nBins)
    {
      binsy[pos] = ii;
      ++pos;
    }
  Int_t  binnumy = sizeof(binsy)/sizeof(Double_t) - 1;

  TF1* f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])-[0]*TMath::Erf((-x-[1])/[2])",0,500);
  f1->SetParameters(0.5,120,40);  
  f1->FixParameter(0,0.5);
  f1->SetLineWidth(2.);
  f1->SetLineColor(kRed);

  TEfficiency* CaloMET_eff = new TEfficiency("CaloMET_eff","CaloMET_eff",50,0.,250.);
  CaloMET_eff->SetLineWidth(2);
  CaloMET_eff->SetTitle("turnOn;MET [GeV]");
  TH1F* CaloMET_num = new TH1F("CaloMET_num","CaloMET_num",25,0.,250.);
  TH1F* CaloMET_den = new TH1F("CaloMET_den","CaloMET_den",25,0.,250.);

  TEfficiency* PFMET_eff = new TEfficiency("PFMET_eff","PFMET_eff",50,0.,250.);
  PFMET_eff->SetLineWidth(2);
  PFMET_eff->SetTitle("turnOn;MET [GeV]");

  TH1F* PFMET_res = new TH1F("PFMET_res","PFMET_res",nBins,min,max);
  PFMET_res->GetXaxis()->SetTitle("(hltPFMET - PFMET)/PFMET");
  TProfile* PFMET_res_vs_PFMET = new TProfile("PFMET_res_vs_PFMET","PFMET_res_vs_PFMET",binnumx,binsx);
  PFMET_res_vs_PFMET->GetYaxis()->SetTitle("(hltPFMET - PFMET)/PFMET");
  PFMET_res_vs_PFMET->GetXaxis()->SetTitle("PFMET");
  TH2F* PFMET_resVsPFMET = new TH2F("PFMET_resVsPFMET","PFMET_resVsPFMET",binnumx,binsx,binnumy,binsy);


  TH1F* CaloMET_res = new TH1F("CaloMET_res","CaloMET_res",nBins,min,max);
  CaloMET_res->GetXaxis()->SetTitle("(hltCaloMET - PFMET)/PFMET");
  TProfile* CaloMET_res_vs_PFMET = new TProfile("CaloMET_res_vs_PFMET","CaloMET_res_vs_PFMET",binnumx,binsx);
  CaloMET_res_vs_PFMET->GetYaxis()->SetTitle("(hltCaloMET - PFMET)/PFMET");
  CaloMET_res_vs_PFMET->GetXaxis()->SetTitle("PFMET");
  TH2F* CaloMET_resVsPFMET = new TH2F("CaloMET_resVsPFMET","CaloMET_resVsPFMET",binnumx,binsx,binnumy,binsy);
  
  
  TGraph* caloMet_rate = new TGraph();
  caloMet_rate->SetName("caloMet_rate");
  caloMet_rate->SetMinimum(0);
  TGraph* PFMet_rate = new TGraph();
  PFMet_rate->SetName("PFMet_rate");
  PFMet_rate->SetMinimum(0);

  
  //study rate for calibration
  int bin=0;
  for(int range=80; range<140; ++range)
    {
      int count_caloMet = 0;
      int count_PFMet = 0;
      for (Long64_t jentry=0; jentry<nentries;++jentry)
   	{
   	  tt->GetEntry(jentry);
   	  if(jentry%1000 == 0)
   	    printProgress((float)jentry/(float)nentries);
  	  
   	  if(hltMetClean_ > range)
   	    ++count_caloMet;
   	  if(hltPFMET_ > range)
   	    ++count_PFMet;
   	}
      caloMet_rate->SetPoint(bin,range,count_caloMet);
      PFMet_rate->SetPoint(bin,range,count_PFMet);
      ++bin;
    }


  //study efficiency
  for (Long64_t jentry=0; jentry<nentries;++jentry)
    {
      tt->GetEntry(jentry);
      if(jentry%5000 == 0)
	printProgress((float)jentry/(float)nentries);

      //### select MET type
      float REFMET_ = -1;
      if(refmet == "calo")
	REFMET_ = caloMET_;
      if(refmet == "pf")
	REFMET_ = PFMET_;  
      if(refmet == "pfnomu")
	REFMET_ = PFMETNoMu_;  

    

      //### fill eff plots using HLT bit and thresholds
      // if(hltAccept_->at(PFMET120PFMHT120CALOBTAG)==1)
      // 	PFMET_eff->Fill(hltPFMHTTightID_ > thr
      // 			&& hltPFMET_ > thr,
      // 			REFMET_);
      
       // if(hltAccept_->at(PFTAU50MET90)==1)
       //  	CaloMET_eff->Fill(hltMetClean_ > thr,
       //  			  caloMET_);

      //### MET only turnon
      //if(f_deltaPhi(hltMetCleanPhi_,caloMETPhi_) < 0.5)

      //CaloMET_eff->Fill(hltMetClean_ > thr, REFMET_);
      CaloMET_eff->Fill(hltMetClean_ > thr, REFMET_);
      PFMET_eff->Fill(hltPFMET_ > thr, REFMET_);


      //fill res plots
      float metres = (hltPFMET_ - REFMET_)/REFMET_;
      PFMET_res->Fill(metres);
      PFMET_res_vs_PFMET->Fill(REFMET_, (metres));
      PFMET_resVsPFMET->Fill(REFMET_, (metres));

      metres = (hltMetClean_ - REFMET_)/REFMET_;
      CaloMET_res->Fill(metres);
      CaloMET_res_vs_PFMET->Fill(REFMET_, (metres));
      CaloMET_resVsPFMET->Fill(REFMET_, (metres));
    }

  std::cout << std::endl;


  // PFMET_eff->Fit(f1,"r");
  // CaloMET_eff->Fit(f1,"r");
  
  TLegend* leg0 = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg0->AddEntry(PFMET_eff,"PFMET","P");

  TCanvas* c1 = new TCanvas();
  c1->Divide(2,1);
  c1->cd(1);
  PFMET_eff->Draw();
  c1->cd(2);
  CaloMET_eff->Draw();
  //leg0->Draw("sames");

  TCanvas* c2 = new TCanvas();
  c2->Divide();
  c2->cd(1);
  caloMet_rate->Draw("AP");
  c2->cd(2);
  PFMet_rate->Draw("AP");

  TCanvas* c3 = new TCanvas();
  c3->Divide(2,1);
  c3->cd(1);
  PFMET_res->Draw();
  c3->cd(2);
  PFMET_res_vs_PFMET->Draw();

  TCanvas* c4 = new TCanvas();
  c4->Divide(2,1);
  c4->cd(1);
  CaloMET_res->Draw();
  c4->cd(2);
  CaloMET_res_vs_PFMET->Draw();

  TCanvas* c5 = new TCanvas();
  CaloMET_resVsPFMET->Draw("COLZ");
  TCanvas* c55 = new TCanvas();
  PFMET_resVsPFMET->Draw("COLZ");


  TGraph* caloResGraphVsPFMET = new TGraph();
  caloResGraphVsPFMET->SetName("caloResGraphVsPFMET");
  caloResGraphVsPFMET->GetXaxis()->SetTitle("PFMET");
  caloResGraphVsPFMET->GetYaxis()->SetTitle("#sigma[(hltCaloMET-PFMET)/PFMET]");
  caloResGraphVsPFMET->SetMinimum(0);
  caloResGraphVsPFMET->SetMaximum(4);
  TGraph* PFResGraphVsPFMET = new TGraph();
  PFResGraphVsPFMET->SetName("PFResGraphVsPFMET");
  PFResGraphVsPFMET->GetXaxis()->SetTitle("PFMET");
  PFResGraphVsPFMET->GetYaxis()->SetTitle("#sigma[(hltPFMET-PFMET)/PFMET]");
  PFResGraphVsPFMET->SetMinimum(0);
  PFResGraphVsPFMET->SetMaximum(4);
/////  //compute resolution
/////  for(int ii = 1; ii < CaloMET_resVsPFMET->GetXaxis()->GetNbins()+1; ++ii)
/////    {
/////      std::cout << "Bin " << ii << " of " << CaloMET_resVsPFMET->GetXaxis()->GetNbins() << std::endl;
/////
/////      TH1F* tmpcalo = (TH1F*)CaloMET_resVsPFMET->ProjectionY("tmpcalo",ii,ii);
/////      //tmpcalo->Write();
/////      float parscalo[4];
/////      FindSmallestInterval(parscalo, tmpcalo, 0.68);
/////      caloResGraphVsPFMET->SetPoint(ii-1,CaloMET_resVsPFMET->GetXaxis()->GetBinCenter(ii), (parscalo[3]-parscalo[2]));
/////
/////      TH1F* tmppf = (TH1F*)PFMET_resVsPFMET->ProjectionY("tmppf",ii,ii);
/////      //tmppf->Write();
/////      float parspf[4];
/////      FindSmallestInterval(parspf, tmppf, 0.68);
/////      PFResGraphVsPFMET->SetPoint(ii-1,PFMET_resVsPFMET->GetXaxis()->GetBinCenter(ii), (parspf[3]-parspf[2]));
/////    }
  TCanvas* c6 = new TCanvas();
  c6->Divide(2,1);
  c6->cd(1);
  caloResGraphVsPFMET->Draw("AP");
  c6->cd(2);
  PFResGraphVsPFMET->Draw("AP");


  //### OUTPUTFILES ###
  // std::string folder;
  // if(matched)
  //   folder = "matched";
  // else
  //   folder = "notMatched";
  // if(wideJets)
  //   folder += "WideJets/";
  // else
  //   folder += "CaloJets/";
  // c1->Print((folder+"turnOn.pdf").c_str(),"pdf");
  // c2->Print((folder+"L1Rates.pdf").c_str(),"pdf");
  // c3->Print((folder+"mjj_res.pdf").c_str(),"pdf");
  // c31->Print((folder+"mjj_diff.pdf").c_str(),"pdf");
  // c4->Print((folder+"deltaR1.pdf").c_str(),"pdf");
  // c41->Print((folder+"deltaR2.pdf").c_str(),"pdf");
  // c5->Print((folder+"pt1_res.pdf").c_str(),"pdf");
  // c6->Print((folder+"pt2_res.pdf").c_str(),"pdf");
  // c7->Print((folder+"deltaR1_vsMjj.pdf").c_str(),"pdf");
  // c8->Print((folder+"deltaPhi_vsMjj.pdf").c_str(),"pdf");
  // c9->Print((folder+"deltaMjj_vsMjj.pdf").c_str(),"pdf");


  PFMET_eff->Write();
  CaloMET_eff->Write();

  PFMET_res_vs_PFMET->Write();
  CaloMET_res_vs_PFMET->Write();

  PFMET_res->Write();
  CaloMET_res->Write();

  caloMet_rate->Write("rate");
  PFMet_rate->Write("rate");
  caloResGraphVsPFMET->Write(caloResGraphVsPFMET->GetName());
  PFResGraphVsPFMET->Write(PFResGraphVsPFMET->GetName());
  
  return 0;


}
