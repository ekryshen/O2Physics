#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "map"
#include "vector"
#include "TEllipse.h"
#include "TGraphErrors.h"
#include "ActsFatras/EventData/Barcode.hpp"

using namespace std;

void draw(){
  gStyle->SetOptStat(0);
  gStyle->SetPalette(0);
  gStyle->SetPadLeftMargin(0.135);
  gStyle->SetPadRightMargin(0.03);
  TFile* fHits = new TFile("hits.root");
  TTree* tHits = (TTree*) fHits->Get("hits");
  float tz = 0;
  float tx = 0;
  float ty = 0;
  UInt_t event_id = 0;
  ULong64_t particle_id = 0;
  ULong64_t geometry_id = 0;
  tHits->SetBranchAddress("event_id",&event_id);
  tHits->SetBranchAddress("particle_id",&particle_id);
  tHits->SetBranchAddress("tx",&tx);
  tHits->SetBranchAddress("ty",&ty);
  tHits->SetBranchAddress("tz",&tz);
  
  double xMax = 1700;
  double yMax = 1700;
  double zMax = 1700;

  TH3F* h = new TH3F("h","h;x;y;z",200,-xMax,xMax,200,-yMax,yMax,200,-zMax,zMax);
  TH2F* hXY = new TH2F("hXY","hXY;x;y",1700,-xMax,xMax,1700,-yMax,yMax);
  TH2F* hXZ = new TH2F("hXZ","hXZ;x;z",1700,-xMax,xMax,1700,-zMax,zMax);
  TH2F* hYZ = new TH2F("hYZ","hYZ;y;z",1700,-yMax,yMax,1700,-zMax,zMax);
  for (int entry=0;entry< tHits->GetEntries();entry++){
    tHits->GetEntry(entry);
//    printf("%f, %f, %f\n", tx, ty, tz);
    h->Fill(tx,ty,tz);
    hXY->Fill(tx,ty);
    hXZ->Fill(tx,tz);
    hYZ->Fill(ty,tz);
  }

  TCanvas* cc = new TCanvas("cc","cc",1040,1100);
  cc->Divide(2,2,0.001,0.001);
  cc->cd(1);
  hXY->Draw();
  cc->cd(2);
  hXZ->Draw();
  cc->cd(3);
  hYZ->Draw();
  cc->cd(4);
  h->Draw();

  TCanvas* cXY = new TCanvas("cXY","cXY",1040,1100);
  hXY->GetXaxis()->SetRangeUser(-1300,1300);
  hXY->GetYaxis()->SetRangeUser(-1300,1300);
  hXY->Draw();

}
