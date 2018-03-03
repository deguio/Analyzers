//simple macro to compare MAHIOFF and MAHION data

void loopdir(TFile* file, std::vector<TH1*>& histovec, std::vector<TGraph*>& graphvec, std::vector<TEfficiency*>& effvec)
{
  TIter next(file->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TH1"))
      {
   	TH1 *h = (TH1*)key->ReadObj();
   	histovec.push_back(h);
      }
    else if (cl->InheritsFrom("TGraph"))
      {
	TGraph *g = (TGraph*)key->ReadObj();
        graphvec.push_back(g);
      }
    else if (cl->InheritsFrom("TEfficiency"))
      {
	TEfficiency *e = (TEfficiency*)key->ReadObj();
        effvec.push_back(e);
      }
  }
}

void setMAHIOFFstyle(TAttMarker* mark, TAttLine* line)
{
  mark->SetMarkerSize(1);  
}

void setMAHIONstyle(TAttMarker* mark, TAttLine* line)
{
  line->SetLineColor(2);
  mark->SetMarkerColor(2);
}

void setGlobalStyle()
{
  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(000000); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

}



void comparePlots(std::string plotsMAHIOFF, std::string plotsMAHION)
{
  setGlobalStyle();
  int setLog[] = {0,0,0,0,0,0,0,0,0,0};
  
  //open files
  TFile* f_plotsMAHIOFF = TFile::Open(plotsMAHIOFF.c_str());
  TFile* f_plotsMAHION = TFile::Open(plotsMAHION.c_str());

  std::vector<TH1*> h_MAHIOFF;
  std::vector<TH1*> h_MAHION;
  std::vector<TGraph*> g_MAHIOFF;
  std::vector<TGraph*> g_MAHION;
  std::vector<TEfficiency*> e_MAHIOFF;
  std::vector<TEfficiency*> e_MAHION;

  loopdir(f_plotsMAHIOFF, h_MAHIOFF, g_MAHIOFF, e_MAHIOFF);
  loopdir(f_plotsMAHION, h_MAHION, g_MAHION, e_MAHION);

  for(int itr=0; itr<h_MAHIOFF.size(); ++itr)
    {
      std::cout << "histoName: " << h_MAHIOFF.at(itr)->GetName() << std::endl;

      setMAHIOFFstyle(h_MAHIOFF.at(itr),h_MAHIOFF.at(itr));
      setMAHIONstyle(h_MAHION.at(itr),h_MAHION.at(itr));

      //legend
      TLegend* leg = new TLegend(0.82, 0.78, 1.03, 0.89);
      //leg->SetFillColor(kWhite);
      leg->AddEntry(h_MAHIOFF.at(itr),"MAHIOFF Data","L");
      leg->AddEntry(h_MAHION.at(itr),"MAHION Data","P");

      TCanvas* c1 = new TCanvas();
      c1 -> cd();
      TPad* p1 = new TPad("p1","p1",0., 0.25, 1., 1.);
      if(setLog[itr] == 1)
	p1->SetLogy();
      TPad* p2 = new TPad("p2","p2",0., 0., 1., 0.25);
      p1 -> Draw();
      p2 -> Draw();

      //draw first pad
      p1 -> cd();
      p1 -> SetGridx();
      p1 -> SetGridy();

      //scale and compute ratio
      h_MAHIOFF.at(itr)->Sumw2();
      h_MAHION.at(itr)->Sumw2();
      h_MAHIOFF.at(itr)->Scale(h_MAHION.at(itr)->Integral()/h_MAHIOFF.at(itr)->Integral());

      h_MAHIOFF.at(itr)->Draw("HISTO,E");
      h_MAHION.at(itr)->Draw("P,sames");
      leg->Draw("same");


      TH1F* ratioHisto;
      if(h_MAHIOFF.at(itr)->GetXaxis()->GetXbins()->GetArray() != nullptr)
	ratioHisto = new TH1F("tmp","tmp",h_MAHIOFF.at(itr)->GetNbinsX(),h_MAHIOFF.at(itr)->GetXaxis()->GetXbins()->GetArray());
      else
	ratioHisto = new TH1F("tmp","tmp",h_MAHIOFF.at(itr)->GetNbinsX(),
			      h_MAHIOFF.at(itr)->GetBinLowEdge(1),
			      h_MAHIOFF.at(itr)->GetBinLowEdge(h_MAHIOFF.at(itr)->GetNbinsX())+h_MAHIOFF.at(itr)->GetBinWidth(1));
      ratioHisto->Sumw2();
      
      for(int bin = 1; bin <= h_MAHIOFF.at(itr)->GetNbinsX(); ++bin)
	{
	  if(h_MAHIOFF.at(itr)->GetBinContent(bin) == 0. || h_MAHION.at(itr)->GetBinContent(bin) == 0.) continue;

	  float valMAHIOFF = h_MAHIOFF.at(itr)->GetBinContent(bin);
	  float valMAHION = h_MAHION.at(itr)->GetBinContent(bin);
	  float sigmaMAHIOFF = h_MAHIOFF.at(itr)->GetBinError(bin);
	  float sigmaMAHION = h_MAHION.at(itr)->GetBinError(bin);

	  ratioHisto -> SetBinContent(bin, valMAHIOFF/valMAHION);
	  ratioHisto -> SetBinError(bin, sqrt((sigmaMAHIOFF/valMAHION)*(sigmaMAHIOFF/valMAHION) +
					      (valMAHIOFF*sigmaMAHION/valMAHION/valMAHION)*(valMAHIOFF*sigmaMAHION/valMAHION/valMAHION)));
	}
      
      p2 -> cd();
      p2 -> SetGridx();
      p2 -> SetGridy();

      // For the axis labels:
      //tdrStyle->SetLabelSize(0.1, "XYZ");
  
      ratioHisto -> GetYaxis() -> SetRangeUser(0., 2.);
      ratioHisto -> DrawCopy("P");
      
      TF1* line = new TF1("line", "1.", -1000000., 1000000.);
      line -> SetLineWidth(2);
      line -> SetLineColor(kRed);
      line -> Draw("same");

      delete ratioHisto;
    }

  for(int itr=0; itr<g_MAHIOFF.size(); ++itr)
    {
      std::cout << "graphName: " << g_MAHIOFF.at(itr)->GetName() << std::endl;

      setMAHIOFFstyle(g_MAHIOFF.at(itr),g_MAHIOFF.at(itr));
      setMAHIONstyle(g_MAHION.at(itr),g_MAHION.at(itr));

      //legend
      TLegend* leg = new TLegend(0.82, 0.78, 1.03, 0.89);
      //leg->SetFillColor(kWhite);                                                                                                                                                       
      leg->AddEntry(g_MAHIOFF.at(itr),"MAHIOFF Data","L");
      leg->AddEntry(g_MAHION.at(itr),"MAHION Data","P");

      TCanvas* g1 = new TCanvas();
      g1 -> cd();
      g1->SetGridx();
      g1->SetGridy();

      g_MAHIOFF.at(itr)->Draw("AP");
      g_MAHION.at(itr)->Draw("P,sames");
      leg->Draw("same");

    }

  for(int itr=0; itr<e_MAHIOFF.size(); ++itr)
    {
      std::cout << "teffName: " << e_MAHIOFF.at(itr)->GetName() << std::endl;

      setMAHIOFFstyle(e_MAHIOFF.at(itr),e_MAHIOFF.at(itr));
      setMAHIONstyle(e_MAHION.at(itr),e_MAHION.at(itr));

      //legend
      TLegend* leg = new TLegend(0.82, 0.78, 1.03, 0.89);
      //leg->SetFillColor(kWhite);                                                                                                                                                       
      leg->AddEntry(e_MAHIOFF.at(itr),"MAHIOFF Data","L");
      leg->AddEntry(e_MAHION.at(itr),"MAHION Data","P");

      TCanvas* e1 = new TCanvas();
      e1 -> cd();
      e1->SetGridx();
      e1->SetGridy();

      e_MAHIOFF.at(itr)->Draw("AP");
      e_MAHION.at(itr)->Draw("P,sames");
      leg->Draw("same");

    }
}
