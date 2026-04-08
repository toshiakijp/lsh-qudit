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
TGraphErrors *grQnx_dep2[nSite+1];
TGraphErrors *grQnx_dep3[nSite+1];
TGraphErrors *grQnx_dep2_fit[nSite+1];
TGraphErrors *grQnx_dep3_fit[nSite+1];
TH1D *hAll;
TH1D *hPass;
TH2D *hFrame;
const int nTmax = 5;
const double TLists[nTmax] = {0.2, 0.4, 0.6, 0.8, 1.0};
const int MyCol[4] = {EColor::kRed, EColor::kOrange, EColor::kGreen+1, EColor::kBlue};

const int nObs=5;
const int nStep=6;
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
TGraphErrors *grPdep2[nTmax][nObs];
TGraphErrors *grPdep3[nObs];

TGraphErrors *grSim[nTmax][nObs];
TGraphErrors *grRaw[nTmax][nObs];
TGraphErrors *grPost[nTmax][nObs];

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
const string PDFName=Form("../PDF/plot_H2Emu_TrotterEM_%s.pdf", Tag.c_str());
//const string PDFName=Form("PaperPlot_H2_%s.pdf", Tag.c_str());
//const string PDFName=Form("PaperPlot_H1_%s.pdf", Tag.c_str());
int InitSetup(void);
int LoadMainResults(void);
int LoadZeroInit(void);
int LoadNumCalc(void);
int LoadSimResults(void);
int ApplyTrotterEM(void);

int plot_TrotterEM(void){
  // Initial Setup
  InitSetup();
  LoadNumCalc();
  LoadSimResults();
  LoadZeroInit();
  
  // Load Main Results
  LoadMainResults();  

  // Apply Trotter Error Mitigation
  ApplyTrotterEM();
  
  hFrame = new TH2D("hFrame", ";Number of Trotter Steps;Discard Rate", 50, -0.5, nStep+0.5, 50, 0.0, 1.0);
  hFrame->GetYaxis()->SetNdivisions(408);
  hFrame->GetXaxis()->SetNdivisions(7);
  c0->Print(Form("%s[", PDFName.c_str()), "pdf");

  for(int iT=0;iT<nTmax;iT++){
    hFrame->Draw();
    grDiscard[iT]->SetMarkerStyle(8);
    grDiscard[iT]->SetMarkerSize(2.0);
    //grDiscard[iT]->SetLineColor(kBlack);
    grDiscard[iT]->SetLineColor(80+iT*5-1);
    grDiscard[iT]->SetMarkerColor(80+iT*5-1);
    grDiscard[iT]->SetLineWidth(2);
    for(int i=0;i<nStep;i++){
      grDiscard[iT]->SetPointEXhigh(i, 0.0);
      grDiscard[iT]->SetPointEXlow(i, 0.0);
    }
    grDiscard[iT]->Draw("sameP");

    ttx = new TLatex(0.80, 0.965, Form("T_{max} = %.1f", TLists[iT]));
    ttx->SetNDC();
    ttx->SetTextSize(0.040);
    ttx->Draw();

    if(fEmulator) txP->Draw();
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
  }

  hFrame->Draw();
  tl = new TLegend(0.58-0.38, 0.68+0.02-0.1, 0.92-0.48, 0.92-0.03);
  tl->SetBorderSize(0);
  tl->SetTextSize(0.045);
  for(int iT=0;iT<nTmax;iT++){
    grDiscard[iT]->SetMarkerSize(1.0);
    grDiscard[iT]->Draw("sameP");
    tl->AddEntry(grDiscard[iT], Form("T_{max} = %.1f", TLists[iT]), "P");
  }
  tl->Draw();
  if(fEmulator) txP->Draw();
  c0->Print(Form("%s", PDFName.c_str()), "pdf");

  
  delete hFrame; hFrame=NULL;
  hFrame = new TH2D("hFrame", ";Number of Trotter Steps;Depolarizing Noise Probability", 50, -0.5, nStep+0.5, 50, 0.0, 1.0);
  hFrame->GetYaxis()->SetNdivisions(408);
  hFrame->GetXaxis()->SetNdivisions(7);

  for(int iT=0;iT<nTmax;iT++){
    hFrame->Draw();
    tl = new TLegend(0.58-0.38, 0.68+0.02, 0.92-0.48, 0.92-0.03);
    for(int j=0;j<3;j++){
      grPdep2[iT][j]->Draw("sameEP");
      tl->AddEntry(grPdep2[iT][j], Form("#it{p} for #it{h_{E}}(%d)", j), "P");
      tl->SetBorderSize(0);
      tl->SetTextSize(0.045);
    }
    ttx = new TLatex(0.80, 0.965, Form("T_{max} = %.1f", TLists[iT]));
    ttx->SetNDC();
    ttx->SetTextSize(0.040);
    ttx->Draw();
    tl->Draw();
    if(fEmulator) txP->Draw();
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
  }

  
  for(int j=0;j<3;j++){
    cout << "j = " << j << endl;
    delete hFrame; hFrame=NULL;
    hFrame = new TH2D("hFrame", ";Time;Electric Energy Density", 50, -0.05, 1.05, 50, 0.0, 4.0);
    hFrame->GetYaxis()->SetNdivisions(408);
    hFrame->GetXaxis()->SetNdivisions(416);
    hFrame->Draw();
    tl = new TLegend(0.58-0.35-0.05, 0.68+0.1, 0.92-0.45-0.1, 0.92);
    tl->AddEntry(grTH[j],  "Theoretical calc.", "l");
    tl->AddEntry(grTrotterEM_Raw[j],  "Raw + Trotter-EM", "P");
    tl->AddEntry(grTrotterEM_Dep2[j], "Depolarizing-EM + Trotter-EM", "P");
    tl->AddEntry(grTrotterEM_Post[j], "Post-EM + Trotter-EM", "P");
    tl->SetBorderSize(0);
    tl->SetTextSize(0.03);

    ttx = new TLatex(0.85, 0.965, Form("link %d", j));
    ttx->SetNDC();
    ttx->SetTextSize(0.040);
    ttx->Draw();

    tl->Draw();

    grTH[j]->Draw("same");
    grTrotterEM_Raw[j]->Draw("sameEP");
    grTrotterEM_Sim[j]->Draw("sameEP");
    grTrotterEM_Dep2[j]->Draw("sameEP");
    grTrotterEM_Post[j]->Draw("sameEP");
    
    if(fEmulator) txP->Draw();
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
    delete ttx; ttx=NULL;
    delete tl; tl=NULL;
  }

  for(int iT=0;iT<nTmax;iT++){
    for(int j=0;j<3;j++){
      hFrame = new TH2D("hFrame", ";time step;Electric Energy Density", 50, 0.0, TLists[iT], 50, 0.0, 3.0);
      //hFrame->GetYaxis()->SetNdivisions(408);
      hFrame->GetXaxis()->SetNdivisions(504);
      hFrame->Draw();
      tl = new TLegend(0.58-0.38+0.38, 0.68+0.02, 0.92-0.48+0.38, 0.92-0.03);
      tl->SetBorderSize(0);
      tl->SetTextSize(0.045);

      ttx = new TLatex(0.85-0.1, 0.965, Form("T_{max} = %.1f, link %d", TLists[iT],  j));
      ttx->SetNDC();
      ttx->SetTextSize(0.040);
      ttx->Draw();

      tli = new TLine(0.0, grTH[j]->Eval(TLists[iT]), TLists[iT], grTH[j]->Eval(TLists[iT]));
      tli->SetLineColor(kGray);
      tli->SetLineStyle(2);
      tli->Draw("same");
      tl->AddEntry(tli, Form("Theoretical calc."), "l");

      grTrotter_Raw[iT][j]->Draw("sameP");
      f2 = (TF1 *)grTrotter_Raw[iT][j]->GetFunction("f1")->Clone("f2");
      f2->SetRange(0.0, TLists[iT]/6.0);
      f2->SetLineStyle(2);
      f2->Draw("same");
      tl->AddEntry(grTrotter_Raw[iT][j], Form("Raw"), "p");

      
      grTrotter_Sim[iT][j]->Draw("sameP");
      f2 = (TF1 *)grTrotter_Sim[iT][j]->GetFunction("f1")->Clone("f2");
      f2->SetRange(0.0, TLists[iT]/6.0);
      f2->SetLineStyle(2);
      f2->Draw("same");
      tl->AddEntry(grTrotter_Sim[iT][j], Form("Noiseless sim."), "p");
      
      grTrotter_Dep2[iT][j]->Draw("sameP");
      f2 = (TF1 *)grTrotter_Dep2[iT][j]->GetFunction("f1")->Clone("f2");
      f2->SetRange(0.0, TLists[iT]/6.0);
      f2->SetLineStyle(2);
      f2->Draw("same");
      tl->AddEntry(grTrotter_Dep2[iT][j], Form("EM Depolarizing"), "p");

      grTrotter_Post[iT][j]->Draw("sameP");
      f2 = (TF1 *)grTrotter_Post[iT][j]->GetFunction("f1")->Clone("f2");
      f2->SetRange(0.0, TLists[iT]/6.0);
      f2->SetLineStyle(2);
      f2->Draw("same");
      tl->AddEntry(grTrotter_Post[iT][j], Form("EM Post-selection"), "p");

      tl->Draw();
      if(fEmulator) txP->Draw();
      c0->Print(Form("%s", PDFName.c_str()), "pdf");

      delete hFrame; hFrame=NULL;
    }// for j
  }// for iT
  
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
    grTrotterEM_Raw[j] = new TGraphErrors();
    grTrotterEM_Sim[j] = new TGraphErrors();
    grTrotterEM_Post[j] = new TGraphErrors();
    grTrotterEM_Dep2[j] = new TGraphErrors();

    grTrotterEM_Sim[j]->SetLineColor(kGreen+1);
    grTrotterEM_Sim[j]->SetMarkerColor(kGreen+1);
    grTrotterEM_Sim[j]->SetMarkerStyle(33);
    grTrotterEM_Raw[j]->SetLineColor(kBlack);
    grTrotterEM_Raw[j]->SetMarkerColor(kBlack);
    grTrotterEM_Raw[j]->SetMarkerStyle(20);
    grTrotterEM_Post[j]->SetLineColor(kRed);
    grTrotterEM_Post[j]->SetMarkerColor(kRed);
    grTrotterEM_Post[j]->SetMarkerStyle(22);
    grTrotterEM_Dep2[j]->SetLineColor(kBlue);
    grTrotterEM_Dep2[j]->SetMarkerColor(kBlue);
    grTrotterEM_Dep2[j]->SetMarkerStyle(21);
  }// for i

  for(int iT=0;iT<nTmax;iT++){
    for(int i=0;i<3;i++){
      grTrotter_Raw[iT][i] = new TGraphErrors();
      grTrotter_Sim[iT][i] = new TGraphErrors();
      grTrotter_Post[iT][i] = new TGraphErrors();
      grTrotter_Dep2[iT][i] = new TGraphErrors();

      grTrotter_Raw[iT][i]->SetLineColor(kBlack);
      grTrotter_Raw[iT][i]->SetMarkerColor(kBlack);
      grTrotter_Raw[iT][i]->SetMarkerStyle(20);
      grTrotter_Sim[iT][i]->SetLineColor(kGreen+1);
      grTrotter_Sim[iT][i]->SetMarkerColor(kGreen+1);
      grTrotter_Sim[iT][i]->SetMarkerStyle(33);
      grTrotter_Post[iT][i]->SetLineColor(kRed);
      grTrotter_Post[iT][i]->SetMarkerColor(kRed);
      grTrotter_Post[iT][i]->SetMarkerStyle(22);
      grTrotter_Dep2[iT][i]->SetLineColor(kBlue);
      grTrotter_Dep2[iT][i]->SetMarkerColor(kBlue);
      grTrotter_Dep2[iT][i]->SetMarkerStyle(21);
    }// for i
    
    for(int j=0;j<nObs;j++){
      grRaw[iT][j] = new TGraphErrors();
      grSim[iT][j] = new TGraphErrors();
      grPost[iT][j] = new TGraphErrors();

      grPdep2[iT][j] = new TGraphErrors(); 
      grPdep2[iT][j]->SetMarkerStyle(21+j);
      if(j!=4){
	grPdep2[iT][j]->SetLineColor(MyCol[j]);
	grPdep2[iT][j]->SetMarkerColor(MyCol[j]);
      }else{
	grPdep2[iT][j]->SetLineColor(kBlack);
	grPdep2[iT][j]->SetMarkerColor(kBlack);
      }
      
      grRaw[iT][j]->SetLineColor(kBlack);
      grPost[iT][j]->SetLineColor(kRed);
      grRaw[iT][j]->SetMarkerColor(kBlack);
      grPost[iT][j]->SetMarkerColor(kRed);
      grRaw[iT][j]->SetMarkerStyle(20);
      grPost[iT][j]->SetMarkerStyle(22);            
      grSim[iT][j]->SetLineColor(kGreen+1);
      grSim[iT][j]->SetMarkerColor(kGreen+1);
      grSim[iT][j]->SetMarkerStyle(33);
    }// for obs
  }// for iT


  for(int i=0;i<nObs;i++){
    hRaw[i] = new TH1D(Form("hRaw_%d", i), "", 1000, 0, 10);
    hPost[i] = new TH1D(Form("hPost_%d", i), "", 1000, 0, 10);
    
    grQnx_dep2[i] = new TGraphErrors();
    grQnx_dep2[i]->SetMarkerStyle(22);
    grQnx_dep2[i]->SetMarkerColor(MyCol[1]);
    grQnx_dep2[i]->SetLineColor(MyCol[1]);

    grQnx_dep2_fit[i] = new TGraphErrors();
    grQnx_dep2_fit[i]->SetMarkerStyle(21);
    grQnx_dep2_fit[i]->SetMarkerColor(MyCol[1]);
    grQnx_dep2_fit[i]->SetLineColor(MyCol[1]);

    grQnx_dep3[i] = new TGraphErrors();
    grQnx_dep3[i]->SetMarkerStyle(22);
    grQnx_dep3[i]->SetMarkerColor(MyCol[2]);
    grQnx_dep3[i]->SetLineColor(MyCol[2]);

    grQnx_dep3_fit[i] = new TGraphErrors();
    grQnx_dep3_fit[i]->SetMarkerStyle(21);
    grQnx_dep3_fit[i]->SetMarkerColor(MyCol[2]);
    grQnx_dep3_fit[i]->SetLineColor(MyCol[2]);

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

      string filename = Form("../savefiles/H2-Emulator_trotterT%d_%s_%dstep.txt", (int)(10*TLists[iT]), Tag.c_str(), iStep+1);
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
	grRaw[iT][j]->SetPoint(iStep, dt*(iStep+1), vRaw[j][iStep]);
	grRaw[iT][j]->SetPointError(iStep, 0.0, eRaw[j][iStep]);

	vPost[j][iStep] = hPost[j]->GetMean();
	ePost[j][iStep] = hPost[j]->GetMeanError();
	grPost[iT][j]->SetPoint(iStep, dt*(iStep+1), vPost[j][iStep]);
	grPost[iT][j]->SetPointError(iStep, 0.0, ePost[j][iStep]);
	
	double vObs = vRaw[j][iStep];
	double eObs = eRaw[j][iStep];
	double vP = vPdep2[iT][j][iStep];
	double eP = ePdep2[iT][j][iStep];
	vQnx_dep2[iT][j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
	eQnx_dep2[iT][j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				       (TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));
	
	//if(j<3 && iStep>0){
	int offset=1;
	if(j<3){
	  grTrotter_Raw[iT][j]->SetPoint(iStep-1+offset, TLists[iT]/(iStep+1), vObs);
	  grTrotter_Raw[iT][j]->SetPointError(iStep-1+offset, 0.0, eObs);

	  grTrotter_Post[iT][j]->SetPoint(iStep-1+offset, TLists[iT]/(iStep+1), vPost[j][iStep]);
	  grTrotter_Post[iT][j]->SetPointError(iStep-1+offset, 0.0, ePost[j][iStep]);
	  
	  grTrotter_Dep2[iT][j]->SetPoint(iStep-1+offset, TLists[iT]/(iStep+1), vQnx_dep2[iT][j][iStep]);
	  grTrotter_Dep2[iT][j]->SetPointError(iStep-1+offset, 0.0, eQnx_dep2[iT][j][iStep]);
	}
      }// for obs
      vDiscard[iT][iStep] = 1.0 - hPost[0]->GetEntries()/hRaw[0]->GetEntries();
    }// for iStep
    grDiscard[iT] = new TGraphAsymmErrors(hPass, hAll);
  }// for iT
    
  cout << " Discard rate for   1-step   2-step   3-step   4-step   5-step   6-step" << endl;
  for(int iT=0;iT<nTmax;iT++){
    cout << "        T = " << Form("%.1f", TLists[iT]) << " :    ";
    for(int iStep=0;iStep<nStep;iStep++){
      cout << Form("%.3f    ", vDiscard[iT][iStep]);
    }// for i step
    cout << endl;
  }// for iT
  cout << endl;

  /*
  cout << " Apply depolarizing noise error mitigation for " << endl;
  cout << "    (nTmax=" << nTmax << ") * (nStep=" << nStep << ") * (nObs=3)" << endl;
  for(int iT=0;iT<nTmax;iT++){
    for(int iStep=0;iStep<nStep;iStep++){
      for(int j=0;j<nObs;j++){
	double vObs = vRaw[j][iStep+1];
	double eObs = eRaw[j][iStep+1];
	double vP = vPdep2[iT][j][iStep];
	double eP = ePdep2[iT][j][iStep];
      
	vQnx_dep2[iT][j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
	eQnx_dep2[iT][j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				       (TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));
      }// for j
    }// for iStep
  }// for iT
  */
  return 0;
}

int LoadZeroInit(void){
  for(int iT=0;iT<nTmax;iT++){
    cout << "Estimate depolarizing noise probability for T = " << TLists[iT] << endl;
    for(int iStep=0;iStep<nStep;iStep++){
      for(int j=0;j<nObs;j++){
	hRaw[j]->Reset();
      }

      string filename = Form("../savefiles/H2-Emulator_trotterT%d_%s_zeroinit_%dstep.txt", (int)(10*TLists[iT]), Tag.c_str(), iStep+1);
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

	grPdep2[iT][j]->SetPoint(iStep, iStep+1, vPdep2[iT][j][iStep]);
	grPdep2[iT][j]->SetPointError(iStep, 0.0, ePdep2[iT][j][iStep]);
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

int ApplyTrotterEM(void){
  for(int iT=0;iT<nTmax;iT++){
    for(int j=0;j<3;j++){
      double yval, yerr;

      f1->SetParameters(1.0, 0.0, -1.0);
      
      grTrotter_Sim[iT][j]->Fit(f1, "Q", "", TLists[iT]/6.0, TLists[iT]/2.0);
      grTrotter_Sim[iT][j]->GetFunction("f1")->SetLineColor(kGreen+1);
      yval = f1->GetParameter(0);
      yerr = f1->GetParError(0);

      grTrotterEM_Sim[j]->SetPoint(iT, TLists[iT], yval);
      grTrotterEM_Sim[j]->SetPointError(iT, 0.0, yerr);

      grTrotter_Raw[iT][j]->Fit(f1, "Q", "", TLists[iT]/6.0, TLists[iT]/2.0);
      grTrotter_Raw[iT][j]->GetFunction("f1")->SetLineColor(kBlack);
      yval = f1->GetParameter(0);
      yerr = f1->GetParError(0);

      grTrotterEM_Raw[j]->SetPoint(iT, TLists[iT], yval);
      grTrotterEM_Raw[j]->SetPointError(iT, 0.0, yerr);

      
      f1->SetParameters(1.0, 0.0, -1.0);      
      grTrotter_Post[iT][j]->Fit(f1, "Q", "", TLists[iT]/6.0, TLists[iT]/2.0);
      grTrotter_Post[iT][j]->GetFunction("f1")->SetLineColor(kRed);
      yval = f1->GetParameter(0);
      yerr = f1->GetParError(0);

      grTrotterEM_Post[j]->SetPoint(iT, TLists[iT], yval);
      grTrotterEM_Post[j]->SetPointError(iT, 0.0, yerr);
      
      f1->SetParameters(1.0, 0.0, -1.0);

      grTrotter_Dep2[iT][j]->Fit(f1, "Q", "", TLists[iT]/6.0, TLists[iT]/2.0);
      grTrotter_Dep2[iT][j]->GetFunction("f1")->SetLineColor(kBlue);
      yval = f1->GetParameter(0);
      yerr = f1->GetParError(0);
      
      grTrotterEM_Dep2[j]->SetPoint(iT, TLists[iT], yval);
      grTrotterEM_Dep2[j]->SetPointError(iT, 0.0, yerr);
    }// for j obj
  }// for iT
  
  return 0;
}

int LoadSimResults(void){
  cout << "Load noiseless results" << endl;
  for(int iT=0;iT<nTmax;iT++){
    hPass->Reset();
    hAll->Reset();
    for(int iStep=0;iStep<nStep;iStep++){
      for(int j=0;j<nObs;j++){
	hRaw[j]->Reset();
      }

      string filename = Form("../savefiles/Sim_trotterT%d_%dstep.txt", (int)(10*TLists[iT]), iStep+1);
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
	grSim[iT][j]->SetPoint(iStep, dt*(iStep+1), vRaw[j][iStep]);
	grSim[iT][j]->SetPointError(iStep, 0.0, eRaw[j][iStep]);
	int offset=1;
	if(j<3){
	  grTrotter_Sim[iT][j]->SetPoint(iStep-1+offset, TLists[iT]/(iStep+1), vRaw[j][iStep]);
	  grTrotter_Sim[iT][j]->SetPointError(iStep-1+offset, 0.0, eRaw[j][iStep]);
	}
      }// for obs
    }// for iStep
  }// for iT

  return 0;
}
