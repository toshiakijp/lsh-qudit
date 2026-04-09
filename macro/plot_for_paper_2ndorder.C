TLegend *tl;
TLatex *ttx;
TText *txP;
TF1 *f1;
TF1 *f2;
double dt=1.25;
TCanvas *c0;
const int nSite=4;
TGraph *grTH[3];
TLine *tli;

bool fEmulator=true;
TGraphErrors *grQnx_dep1[nSite+1];
TGraphErrors *grQnx_dep3[nSite+1];
TGraphErrors *grQnx_dep2_fit[nSite+1];
TGraphErrors *grQnx_dep3_fit[nSite+1];
TH1D *hAll;
TH1D *hPass;
TH2D *hFrame;
TH2D *hRatioFrame;
TPad *pad1[2];
const int nTmax = 5;
const double TLists[nTmax] = {0.2, 0.4, 0.6, 0.8, 1.0};
const int MyCol[4] = {EColor::kRed, EColor::kOrange, EColor::kGreen+1, EColor::kBlue};

const int nObs=5;
const int nStep=1;
//const int nStep=4;
TH1D *hRaw[nObs];
TH1D *hPost[nObs];

TGraphAsymmErrors *grDiscard[nTmax];
TGraphErrors *grTrotter_Raw[nTmax][3];
TGraphErrors *grTrotter_Sim[nTmax][3];
TGraphErrors *grTrotter_Post[nTmax][3];
TGraphErrors *grTrotter_Dep2[nTmax][3];
TGraphErrors *grTrotterEM_Raw[3];
TGraphErrors *grTrotterEM_Sim[3];
TGraphErrors *grTrotterEM_Post[3];
TGraphErrors *grTrotterEM_Dep2[3];
double vRaw[nObs][nStep+1];
double eRaw[nObs][nStep+1];

double vQnx_dep1[nObs][nStep+1];
double eQnx_dep1[nObs][nStep+1];
double vQnx_dep2[nTmax][nObs][nStep+1];
double eQnx_dep2[nTmax][nObs][nStep+1];
double vQnx_dep2_fit[nObs][nStep+1];
double eQnx_dep2_fit[nObs][nStep+1];
double vQnx_dep3[nObs][nStep+1];
double eQnx_dep3[nObs][nStep+1];
double vQnx_dep3_fit[nObs][nStep+1];
double eQnx_dep3_fit[nObs][nStep+1];
double vDiscard[nTmax][nStep+1];

double vPost[5][nStep+1];
double ePost[nObs][nStep+1];

double vPdep1[nObs][nStep+1];
double ePdep1[nObs][nStep+1];

double vPdep2[nTmax][nObs][nStep+1];
double ePdep2[nTmax][nObs][nStep+1];
double vPdep2_fit[nObs][nStep+1];
double ePdep2_fit[nObs][nStep+1];

double vPdep3[nObs][nStep+1];
double ePdep3[nObs][nStep+1];
double vPdep3_fit[nObs][nStep+1];
double ePdep3_fit[nObs][nStep+1];

TGraphErrors *grPdep1[nObs];
TGraphErrors *grPdep3[nObs];

TGraphErrors *grSim[nTmax][nObs];
TGraphErrors *grRaw[nStep][nObs];
TGraphErrors *grPost[nStep][nObs];
TGraphErrors *grPdep2[nTmax][nObs];
TGraphErrors *grQnx_dep2[nStep][nObs];

TGraphErrors *grRaw_Ratio[nStep][nObs];
TGraphErrors *grPost_Ratio[nStep][nObs];
TGraphErrors *grQnx_dep2_Ratio[nStep][nObs];

TH1D *h1;
ifstream ifs;

double N_L0,N_L1,N_L2,N_Lfix;
double N_R0,N_R1,N_R2,N_Rfix;
int anc,counts;
double vIn[nObs];

const double ccc[nObs] = {3./16., 10./16., 2.0, 10./16., 0.5};

const string Tag="sf1";
//const string Tag="initfull";
//const string Tag="full";
//const string Tag="old";
//const string Tag="precomp";
//const string Tag="test";
const string PDFName=Form("../PDF/PaperPlot_H2Emu_2ndOrder_%s.pdf", Tag.c_str());
//const string PDFName=Form("PaperPlot_H2_%s.pdf", Tag.c_str());
//const string PDFName=Form("PaperPlot_H1_%s.pdf", Tag.c_str());
int InitSetup(void);
int LoadMainResults(void);
int LoadZeroInit(void);
int LoadNumCalc(void);

int plot_for_paper_2ndorder(void){
  // Initial Setup
  InitSetup();
  LoadNumCalc();
  LoadZeroInit();
  
  // Load Main Results
  LoadMainResults();  

  hFrame = new TH2D("hFrame", ";Number of Trotter Steps;Discard Rate", 50, -0.5, nStep+0.5, 50, 0.0, 1.0);
  hFrame->GetYaxis()->SetNdivisions(408);
  hFrame->GetXaxis()->SetNdivisions(7);
  c0->Print(Form("%s[", PDFName.c_str()), "pdf");
  
  for(int j=0;j<3;j++){
    cout << "j = " << j << endl;
    delete hFrame; hFrame=NULL;
    if(j==0){
      hFrame = new TH2D("hFrame", ";Time;Electric Energy Density", 50, -0.05, 1.05, 50, 0.0, 1.0);
    }else{
      hFrame = new TH2D("hFrame", ";Time;Electric Energy Density", 50, -0.05, 1.05, 50, 0.0, 2.2);
    }
    hFrame->GetYaxis()->SetNdivisions(408);
    hFrame->GetXaxis()->SetNdivisions(208);
    hFrame->Draw();
    tl = new TLegend(0.58-0.35-0.01, 0.68+0.1-0.55, 0.92-0.45-0.1, 0.92-0.4);
    tl->AddEntry(grTH[j],  "Theoretical calc.", "l");

    for(int iStep=0;iStep<nStep;iStep++){
      tl->AddEntry(grRaw[iStep][j],  Form("Raw (%d-step)", iStep+1), "P");
      tl->AddEntry(grQnx_dep2[iStep][j], Form("EM Depolarizing (%d-step)", iStep+1), "P");
      tl->AddEntry(grPost[iStep][j], Form("EM Post-selection (%d-step)", iStep+1), "P");
    }
    tl->SetBorderSize(0);
    tl->SetTextSize(0.03);

    ttx = new TLatex(0.85, 0.965, Form("link %d", j));
    ttx->SetNDC();
    ttx->SetTextSize(0.040);
    ttx->Draw();

    tl->Draw();

    grTH[j]->Draw("same");
    for(int iStep=0;iStep<nStep;iStep++){
      grRaw[iStep][j]->Draw("sameP");
      grQnx_dep2[iStep][j]->Draw("sameP");
      grPost[iStep][j]->Draw("sameP");
    }
    
    if(fEmulator) txP->Draw();
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
    delete ttx; ttx=NULL;
    delete tl; tl=NULL;
  }// for j

  for(int j=0;j<3;j++){
    cout << "j = " << j << endl;
    delete hFrame; hFrame=NULL;
    if(j==0){
      //hFrame = new TH2D("hFrame", ";Time;Electric Energy Density", 50, -0.05, 1.05, 50, 0.0, 1.0);
      hFrame = new TH2D("hFrame", ";;Electric Energy Density", 50, -0.05, 1.05, 50, 0.0, 1.0);
    }else{
      //hFrame = new TH2D("hFrame", ";Time;Electric Energy Density", 50, -0.05, 1.05, 50, 0.0, 2.2);
      hFrame = new TH2D("hFrame", ";;Electric Energy Density", 50, -0.05, 1.05, 50, 0.0, 2.2);
    }
    delete hRatioFrame; hRatioFrame=NULL;
    hRatioFrame = new TH2D("hRatioFrame", Form(";Time;Difference [%%]"), 50, -0.05, 1.05, 50, -30, 30);
    hFrame->GetYaxis()->SetTitleOffset(0.7);
    hFrame->GetYaxis()->SetTitleSize(0.07);
    hFrame->GetYaxis()->SetLabelSize(0.07);
    hRatioFrame->GetXaxis()->SetLabelSize(0.15);
    hRatioFrame->GetXaxis()->SetTitleSize(0.15);
    hRatioFrame->GetXaxis()->SetTitleOffset(1.0);
    hRatioFrame->GetYaxis()->SetLabelSize(0.13);
    hRatioFrame->GetYaxis()->SetTitleSize(0.13);
    hRatioFrame->GetYaxis()->SetTitleOffset(0.38);
    hRatioFrame->GetYaxis()->SetNdivisions(505);
    hRatioFrame->GetXaxis()->SetNoExponent(true);

    hFrame->GetXaxis()->SetLabelSize(0);
    c0->cd();
    pad1[0] = new TPad("pad1_0", "pad1_0", 0.0, 0.35, 1.0, 1.0);
    pad1[1] = new TPad("pad1_1", "pad1_1", 0.0, 0.0, 1.0, 0.35);
    pad1[0]->SetBottomMargin(0.03);
    pad1[0]->SetLeftMargin(0.12);
    pad1[1]->SetTopMargin(0);
    pad1[1]->SetLeftMargin(0.12);
    pad1[1]->SetBottomMargin(0.35);
    pad1[1]->SetGridy();

    hFrame->GetYaxis()->SetNdivisions(408);
    hFrame->GetXaxis()->SetNdivisions(208);
    hRatioFrame->GetXaxis()->SetNdivisions(208);
    c0->Clear();
    pad1[0]->Draw();
    pad1[0]->cd();
    hFrame->Draw();
    tl = new TLegend(0.58-0.35-0.05, 0.68+0.1-0.7, 0.92-0.35-0.15, 0.92-0.5);
    tl->AddEntry(grTH[j],  "Theoretical calc.", "l");

    for(int iStep=0;iStep<nStep;iStep++){
      tl->AddEntry(grRaw[iStep][j],  Form("Raw"), "P");
      tl->AddEntry(grQnx_dep2[iStep][j], Form("EM Depolarizing"), "P");
      tl->AddEntry(grPost[iStep][j], Form("EM Post-selection"), "P");
      //tl->AddEntry(grRaw[iStep][j],  Form("Raw (%d-step)", iStep+1), "P");
      //tl->AddEntry(grQnx_dep2[iStep][j], Form("EM Depolarizing (%d-step)", iStep+1), "P");
      //tl->AddEntry(grPost[iStep][j], Form("EM Post-selection (%d-step)", iStep+1), "P");

      grRaw_Ratio[iStep][j]      = (TGraphErrors *)grRaw[iStep][j]->Clone();
      grPost_Ratio[iStep][j]     = (TGraphErrors *)grPost[iStep][j]->Clone();
      grQnx_dep2_Ratio[iStep][j] = (TGraphErrors *)grQnx_dep2[iStep][j]->Clone();

      grRaw_Ratio[iStep][j]     ->Clear();
      grPost_Ratio[iStep][j]    ->Clear();
      grQnx_dep2_Ratio[iStep][j]->Clear();
      for(int iT=0;iT<nTmax;iT++){
	grRaw_Ratio[iStep][j]->SetPoint(iT, TLists[iT], 100.0*(grRaw[iStep][j]->GetPointY(iT)/grTH[j]->Eval(TLists[iT]) - 1.0));
	grRaw_Ratio[iStep][j]->SetPointError(iT, 0.0, 100.0*(grRaw[iStep][j]->GetErrorY(iT)/grTH[j]->Eval(TLists[iT])));
	grPost_Ratio[iStep][j]->SetPoint(iT, TLists[iT], 100.0*(grPost[iStep][j]->GetPointY(iT)/grTH[j]->Eval(TLists[iT]) - 1.0));
	grPost_Ratio[iStep][j]->SetPointError(iT, 0.0, 100.0*(grPost[iStep][j]->GetErrorY(iT)/grTH[j]->Eval(TLists[iT])));
	grQnx_dep2_Ratio[iStep][j]->SetPoint(iT, TLists[iT], 100.0*(grQnx_dep2[iStep][j]->GetPointY(iT)/grTH[j]->Eval(TLists[iT]) - 1.0));
	grQnx_dep2_Ratio[iStep][j]->SetPointError(iT, 0.0, 100.0*(grQnx_dep2[iStep][j]->GetErrorY(iT)/grTH[j]->Eval(TLists[iT])));

      }// for iT

      for(int iT=0;iT<nTmax;iT++) cout << Form("Raw   (%.1f): %.1f", TLists[iT], 100.0*(grRaw[iStep][j]->GetPointY(iT)/grTH[j]->Eval(TLists[iT]) - 1.0)) << endl;
      for(int iT=0;iT<nTmax;iT++) cout << Form("Dep   (%.1f): %.1f", TLists[iT], 100.0*(grQnx_dep2[iStep][j]->GetPointY(iT)/grTH[j]->Eval(TLists[iT]) - 1.0)) << endl;
      for(int iT=0;iT<nTmax;iT++) cout << Form("Post  (%.1f): %.1f", TLists[iT], 100.0*(grPost[iStep][j]->GetPointY(iT)/grTH[j]->Eval(TLists[iT]) - 1.0)) << endl;
    }

    tl->SetBorderSize(0);
    tl->SetTextSize(0.05);

    ttx = new TLatex(0.85, 0.965, Form("link %d", j));
    ttx->SetNDC();
    ttx->SetTextSize(0.05);
    ttx->Draw();

    tl->Draw();

    grTH[j]->Draw("same");
    for(int iStep=0;iStep<nStep;iStep++){
      grRaw[iStep][j]->Draw("sameP");
      grQnx_dep2[iStep][j]->Draw("sameP");
      grPost[iStep][j]->Draw("sameP");
    }
    
    if(fEmulator) txP->Draw();

    c0->cd();
    pad1[1]->SetFillColorAlpha(0, 0.0);
    pad1[1]->Draw();
    pad1[1]->cd();
    hRatioFrame->Draw();
    for(int iStep=0;iStep<nStep;iStep++){
      grRaw_Ratio[iStep][j]->Draw("sameP");
      grQnx_dep2_Ratio[iStep][j]->Draw("sameP");
      grPost_Ratio[iStep][j]->Draw("sameP");
    }

    c0->SaveAs(Form("../PDF/Nsite4_2ndorder_link%d.pdf", j));
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
    delete ttx; ttx=NULL;
    delete tl; tl=NULL;
  }// for j

  c0->Print(Form("%s]", PDFName.c_str()), "pdf");
    
  return 0;
}

int InitSetup(void){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetPaperSize(20,26);
  // set margin sizes
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  Double_t tsize=0.05;
  Int_t font=42; // Helvetica
  gStyle->SetTextFont(font);

  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");
  
  gStyle->SetLabelSize(tsize,"x");
  gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(tsize,"y");
  gStyle->SetTitleSize(tsize,"y");
  gStyle->SetLabelSize(tsize,"z");
  gStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  gStyle->SetErrorX(0.0001);
  gStyle->SetEndErrorSize(0.);

  txP = new TText(0.15, 0.965, "Emulator Results");
  txP->SetNDC();
  txP->SetTextSize(0.040);
  txP->SetTextColor(kRed);
  
  c0 = new TCanvas("c0", "c0", 800, 600);
  hAll  = new TH1D("hAll", "", nStep, 0.5, nStep+0.5);
  hPass = new TH1D("hPass", "", nStep, 0.5, nStep+0.5);

  f1 = new TF1("f1", "[0]+x*[1]+x*x*[2]", 0.0, 1.0);
  f1->SetLineWidth(1);
  f1->SetNpx(1000);

  for(int j=0;j<3;j++){
    grTH[j] = new TGraph();
    grTH[j]->SetLineStyle(2);
    grTH[j]->SetLineWidth(2);
  }// for i

  for(int iStep=0;iStep<nStep;iStep++){
    for(int j=0;j<nObs;j++){
      grRaw[iStep][j] = new TGraphErrors();
      grPost[iStep][j] = new TGraphErrors();

      grQnx_dep2[iStep][j] = new TGraphErrors();
      grQnx_dep2[iStep][j]->SetMarkerStyle(21+iStep*4);
      grQnx_dep2[iStep][j]->SetMarkerColor(kBlue);
      grQnx_dep2[iStep][j]->SetLineColor(kBlue);

      grRaw[iStep][j]->SetLineColor(kBlack);
      grPost[iStep][j]->SetLineColor(kRed);
      grRaw[iStep][j]->SetMarkerColor(kBlack);
      grPost[iStep][j]->SetMarkerColor(kRed);
      grRaw[iStep][j]->SetMarkerStyle(20+iStep*4);
      grPost[iStep][j]->SetMarkerStyle(22+iStep*4);
    }// for obs
  }// for iT

  for(int i=0;i<nObs;i++){
    hRaw[i] = new TH1D(Form("hRaw_%d", i), "", 1000, 0, 10);
    hPost[i] = new TH1D(Form("hPost_%d", i), "", 1000, 0, 10);
  }

  return 0;
}

int LoadMainResults(void){
  cout << "Load main results" << endl;
  for(int iT=0;iT<nTmax;iT++){
    hPass->Reset();
    hAll->Reset();
    for(int iStep=0;iStep<nStep;iStep++){
      for(int j=0;j<nObs;j++){
	hRaw[j]->Reset();
	hPost[j]->Reset();
      }

      string filename = Form("../savefiles/H2-Emulator_2ndorderT%d_%s_%dstep.txt", (int)(10*TLists[iT]), Tag.c_str(), iStep+1);
      ifs.open(filename.c_str());
      cout << "   - Load File: " << filename << endl;
    
      while(ifs>>N_L0>>N_L1>>N_L2>>N_Lfix>>N_R0>>N_R1>>N_R2>>N_Rfix>>anc>>counts){
	vIn[0] = (N_L0/2.0)*(N_L0/2.0 + 1);
	vIn[1] = (N_L1/2.0)*(N_L1/2.0 + 1);
	vIn[2] = (N_L2/2.0)*(N_L2/2.0 + 1);
	vIn[3] = (N_Lfix/2.0)*(N_Lfix/2.0 + 1);
	vIn[4] = (N_Rfix/2.0)*(N_Rfix/2.0 + 1);
      
	for(int j=0;j<nObs;j++){
	  for(int i=0;i<counts;i++){	  
	    hRaw[j]->Fill(vIn[j]);
	    if(j==0)
	      hAll->Fill(iStep+1);
	  
	    if(anc==0 && N_Rfix==0 && N_L0==N_R0 && N_L1==N_R1 && N_L2==N_R2){
	      hPost[j]->Fill(vIn[j]);
	    }else{
	      if(j==0)
		hPass->Fill(iStep+1);
	    }
	  } // for counts
	}// for j
      }// while
      ifs.close();
    
      for(int j=0;j<nObs;j++){
	vRaw[j][iStep] = hRaw[j]->GetMean();
	eRaw[j][iStep] = hRaw[j]->GetMeanError();
	grRaw[iStep][j]->SetPoint(iT, TLists[iT], vRaw[j][iStep]);
	grRaw[iStep][j]->SetPointError(iT,   0.0, eRaw[j][iStep]);

	vPost[j][iStep] = hPost[j]->GetMean();
	ePost[j][iStep] = hPost[j]->GetMeanError();
	grPost[iStep][j]->SetPoint(iT, TLists[iT], vPost[j][iStep]);
	grPost[iStep][j]->SetPointError(iT,   0.0, ePost[j][iStep]);
	
	double vObs = vRaw[j][iStep];
	double eObs = eRaw[j][iStep];
	double vP = vPdep2[iT][j][iStep];
	double eP = ePdep2[iT][j][iStep];
	vQnx_dep2[iT][j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
	eQnx_dep2[iT][j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				       (TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));
	grQnx_dep2[iStep][j]->SetPoint(iT, TLists[iT], vQnx_dep2[iT][j][iStep]);
	grQnx_dep2[iStep][j]->SetPointError(iT, 0.0, eQnx_dep2[iT][j][iStep]);	
      }// for obs
      vDiscard[iT][iStep] = 1.0 - hPost[0]->GetEntries()/hRaw[0]->GetEntries();
    }// for iStep
    grDiscard[iT] = new TGraphAsymmErrors(hPass, hAll);
  }// for iT
    
  cout << " Discard rate for   1-step   2-step" << endl;
  for(int iT=0;iT<nTmax;iT++){
    cout << "        T = " << Form("%.1f", TLists[iT]) << " :    ";
    for(int iStep=0;iStep<nStep;iStep++){
      cout << Form("%.3f    ", vDiscard[iT][iStep]);
    }// for i step
    cout << endl;
  }// for iT
  cout << endl;

  return 0;
}

int LoadZeroInit(void){
  for(int iT=0;iT<nTmax;iT++){
    cout << "Estimate depolarizing noise probability for T = " << TLists[iT] << endl;
    for(int iStep=0;iStep<nStep;iStep++){
      for(int j=0;j<nObs;j++){
	hRaw[j]->Reset();
      }

      string filename = Form("../savefiles/H2-Emulator_2ndorderT%d_%s_zeroinit_%dstep.txt", (int)(10*TLists[iT]), Tag.c_str(), iStep+1);
      ifs.open(filename.c_str());
      cout << "   - Load File: " << filename << endl;
    
      while(ifs>>N_L0>>N_L1>>N_L2>>N_Lfix>>N_R0>>N_R1>>N_R2>>N_Rfix>>anc>>counts){
	vIn[0] = (N_L0/2.0)*(N_L0/2.0 + 1);
	vIn[1] = (N_L1/2.0)*(N_L1/2.0 + 1);
	vIn[2] = (N_L2/2.0)*(N_L2/2.0 + 1);
	vIn[3] = (N_Lfix/2.0)*(N_Lfix/2.0 + 1);
	vIn[4] = (N_Rfix/2.0)*(N_Rfix/2.0 + 1);
      
	for(int j=0;j<nObs;j++){
	  for(int i=0;i<counts;i++){	  
	    hRaw[j]->Fill(vIn[j]);
	  } // for counts
	}// for j
      }// while
      ifs.close();

      for(int j=0;j<nObs;j++){
	vRaw[j][iStep] = hRaw[j]->GetMean();
	eRaw[j][iStep] = hRaw[j]->GetMeanError();

	vPdep2[iT][j][iStep]= 1.0 - (vRaw[j][iStep] - ccc[j])/(0.0 - ccc[j]);
	ePdep2[iT][j][iStep]= TMath::Abs(eRaw[j][iStep]/(0.0 - ccc[j]));

	//grPdep2[iT][j]->SetPoint(iStep, iStep+1, vPdep2[iT][j][iStep]);
	//grPdep2[iT][j]->SetPointError(iStep, 0.0, ePdep2[iT][j][iStep]);
      }// for obj
    }// for iStep
  }// for iT
  
  for(int iT=0;iT<nTmax;iT++){
    cout << " T = " << TLists[iT] << endl;
    cout << "  p for    H_E(0)   H_E(1)   H_E(2)" << endl;
    for(int iStep=0;iStep<nStep;iStep++){
      cout << Form("  %d-step   ", iStep+1);
      for(int j=0;j<nObs;j++){
	cout << Form("%.3f    ", vPdep2[iT][j][iStep]);
      }// for obj
      cout << endl;
    }// for iStep
  }// for iT
  
  return 0;
}

int LoadNumCalc(){
  string filename = "../savefiles/theoretical_calc.txt";
  ifs.open(filename.c_str());
  cout << "Load File: " << filename << endl;

  double tp;
  double expval[3];
  int iP=0;
  while(ifs>>tp>>expval[0]>>expval[1]>>expval[2]){
    for(int i=0;i<3;i++){
      grTH[i]->SetPoint(iP, tp, expval[i]);
    }// for i
    iP++;
  }// while
  ifs.close();
  
  return 0;
}

