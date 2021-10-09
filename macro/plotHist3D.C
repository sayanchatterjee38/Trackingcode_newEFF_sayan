#include "RiceStyle.h"

void plotHist2D() {

  RiceStyle();

  gStyle->SetOptStat(0);

  TFile *f = new TFile("postrk.root");

  char ndir[256] = "HITrackCorrections";
  
  // sim-to-reco hists
  TH3F *hSim = (TH3F*) f->Get(Form("%s/hsim3D",ndir));  
  TH3F *hEff = (TH3F*) f->Get(Form("%s/heff3D",ndir)); 
  TH3F *hMul = (TH3F*) f->Get(Form("%s/hmul3D",ndir)); 

  // reco-to-sim hists
  TH3F *hRec = (TH3F*) f->Get(Form("%s/hrec3D",ndir)); 
  TH3F *hFak = (TH3F*) f->Get(Form("%s/hfak3D",ndir)); 
  TH3F *hSec = (TH3F*) f->Get(Form("%s/hsec3D",ndir)); 
  
  // ratio histograms
  TH3F *Eff3D = (TH3F*) hEff->Clone("Eff3D");
  TH3F *rMul = (TH3F*) hMul->Clone("rMul");
  TH3F *Fak3D = (TH3F*) hFak->Clone("Fak3D");
  TH3F *rSec = (TH3F*) hSec->Clone("rSec");

  // reco efficiency fraction
  TCanvas *c2 = new TCanvas("c2","Reco Efficiency Fraction",600,500);
  Eff3D->Divide(hEff,hSim,1,1,"B");
  Eff3D->SetMaximum(1.0); Eff3D->SetMinimum(0.0);
  Eff3D->SetTitle("Absolute Efficiency");
  Eff3D->Draw("surf1");

  // fake reco fraction
  TCanvas *c4 = new TCanvas("c4","Fake Reco Fraction",600,500);
  Fak3D->Divide(hFak,hRec,1,1,"B");
  Fak3D->SetMaximum(0.25); Fak3D->SetMinimum(0.0);
  Fak3D->SetTitle("Fake Reconstruction Fraction");
  Fak3D->Draw("surf1");

  TFile *f1 = new TFile("newEff_postrk.root","RECREATE");
  Eff3D->Write();
  Fak3D->Write();
  f1->Close();
  
}



