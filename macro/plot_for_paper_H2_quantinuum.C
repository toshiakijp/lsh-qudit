TLegend *tl;
TText *ttx;
TF1 *f1;
double dt=1.25;
TCanvas *c0;
const int nSite=4;
//TGraph *grSim[nSite];
TGraphErrors *grQnx_dep1[nSite+1];
TGraphErrors *grQnx_dep2[nSite+1];
TGraphErrors *grQnx_dep3[nSite+1];
TGraphErrors *grQnx_dep2_fit[nSite+1];
TGraphErrors *grQnx_dep3_fit[nSite+1];
TGraph *grQnx_deplin[nSite];
TGraph *grQnx_post[nSite];
TGraphAsymmErrors *grDiscard;
TH1D *hAll;
TH1D *hPass;
TH2D *hFrame;
const int MyCol[4] = {EColor::kRed, EColor::kOrange, EColor::kGreen+1, EColor::kBlue};

const int nObs=5;
const int nStep=6;
//const int nStep=4;
TH1D *hSim[nObs];
TH1D *hRaw[nObs];
TH1D *hPost[nObs];

double vRaw[nObs][nStep+1];
double eRaw[nObs][nStep+1];

double vQnx_dep1[nObs][nStep+1];
double eQnx_dep1[nObs][nStep+1];
double vQnx_dep2[nObs][nStep+1];
double eQnx_dep2[nObs][nStep+1];
double vQnx_dep2_fit[nObs][nStep+1];
double eQnx_dep2_fit[nObs][nStep+1];
double vQnx_dep3[nObs][nStep+1];
double eQnx_dep3[nObs][nStep+1];
double vQnx_dep3_fit[nObs][nStep+1];
double eQnx_dep3_fit[nObs][nStep+1];
double vDiscard[nStep+1];

//double vSim[5][nStep+1];
const double vSim[nObs][6+1]={{0.75, 0.058365000000000014, 0.4886174999999999, 0.4945875000000001, 0.3296774999999999, 0.3161175000000001, 0.37761749999999983},
				  {2.0, 0.35956000000000016, 0.5294749999999998, 1.0853350000000004, 0.42379999999999995, 0.4748275000000001, 0.6441724999999998},
				  {2.0, 0.23769250000000003, 0.5223649999999999, 1.4649300000000005, 0.3910199999999999, 0.6453099999999999, 0.9654649999999999},
				  {0.75, 0.7500000000000002, 0.75, 0.7500000000000002, 0.75, 0.7499999999999999, 0.7499999999999999},
				  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

double vPost[5][nStep+1];
double ePost[nObs][nStep+1];

double vPdep1[nObs][nStep+1];
double ePdep1[nObs][nStep+1];

double vPdep2[nObs][nStep+1];
double ePdep2[nObs][nStep+1];
double vPdep2_fit[nObs][nStep+1];
double ePdep2_fit[nObs][nStep+1];

double vPdep3[nObs][nStep+1];
double ePdep3[nObs][nStep+1];
double vPdep3_fit[nObs][nStep+1];
double ePdep3_fit[nObs][nStep+1];

TGraphErrors *grPdep1[nObs];
TGraphErrors *grPdep2[nObs];
TGraphErrors *grPdep3[nObs];

TGraphErrors *grSim[5];
TGraphErrors *grRaw[nObs];
TGraphErrors *grPost[5];

TH1D *h1;
ifstream ifs;

double N_L0,N_L1,N_L2,N_Lfix;
double N_R0,N_R1,N_R2,N_Rfix;
int anc,counts;
double vIn[nObs];

const double ccc[nObs] = {3./16., 10./16., 2.0, 10./16., 0.5};

const string Tag="initfull_sf1";
//const string Tag="initfull";
//const string Tag="full";
//const string Tag="old";
//const string Tag="precomp";
//const string Tag="test";
const string PDFName=Form("../PDF/PaperPlot_H2_%s_test.pdf", Tag.c_str());
//const string PDFName=Form("PaperPlot_H2_%s.pdf", Tag.c_str());
//const string PDFName=Form("PaperPlot_H1_%s.pdf", Tag.c_str());
int InitSetup(void);
int LoadMainResults(void);
int LoadZeroInit(void);
int LoadPairResults(void);

int plot_for_paper_H2_quantinuum(void){
  // Initial Setup
  InitSetup();
  
  LoadZeroInit();
  LoadPairResults();
  
  // Load Main Results
  LoadMainResults();  

  hFrame = new TH2D("hFrame", ";Number of Trotter Steps;Discard Rate", 50, -0.5, nStep+0.5, 50, 0.0, 1.0);
  hFrame->GetYaxis()->SetNdivisions(408);
  hFrame->GetXaxis()->SetNdivisions(7);
  c0->Print(Form("%s[", PDFName.c_str()), "pdf");

  hFrame->Draw();
  grDiscard->SetMarkerStyle(8);
  grDiscard->SetMarkerSize(2.0);
  grDiscard->SetLineColor(kBlack);
  grDiscard->SetLineWidth(2);
  for(int i=0;i<nStep+1;i++){
    grDiscard->SetPointEXhigh(i, 0.0);
    grDiscard->SetPointEXlow(i, 0.0);
  }
  grDiscard->Draw("sameP");
  c0->Print(Form("%s", PDFName.c_str()), "pdf");
  
  cout << "Discard rate: " << endl;
  for(int i=0;i<nStep+1;i++){
    cout << i << " step: " << grDiscard->GetPointY(i) << " + " << grDiscard->GetErrorYhigh(i) << " - " << grDiscard->GetErrorYlow(i) << endl;
  }
  /*
  delete hFrame; hFrame=NULL;
  hFrame = new TH2D("hFrame", ";Number of Trotter Steps;Depolarizing Noise Probability", 50, -0.5, nStep+0.5, 50, 0.0, 1.0);
  hFrame->GetYaxis()->SetNdivisions(408);
  hFrame->GetXaxis()->SetNdivisions(5);
  hFrame->Draw();
  grPdep1[3]->Draw("sameEP");
  c0->Print(Form("%s", PDFName.c_str()), "pdf");
  c0->Print(Form("%s]", PDFName.c_str()), "pdf");
  */
  delete hFrame; hFrame=NULL;
  hFrame = new TH2D("hFrame", ";Number of Trotter Steps;Depolarizing Noise Probability", 50, -0.5, nStep+0.5, 50, 0.0, 1.0);
  hFrame->GetYaxis()->SetNdivisions(408);
  hFrame->GetXaxis()->SetNdivisions(7);
  hFrame->Draw();
  for(int j=0;j<5;j++){
    grPdep1[j]->Draw("sameEP");
  }
  c0->Print(Form("%s", PDFName.c_str()), "pdf");

  hFrame->Draw();
  tl = new TLegend(0.58-0.38, 0.68+0.02, 0.92-0.48, 0.92-0.03);
  tl->AddEntry(grPdep1[3], "global #it{p}", "P");
  tl->SetBorderSize(0);
  tl->SetTextSize(0.045);
  grPdep1[3]->Draw("sameEP");
  tl->Draw();
  c0->Print(Form("%s", PDFName.c_str()), "pdf");

  hFrame->Draw();
  tl = new TLegend(0.58-0.38, 0.68+0.02, 0.92-0.48, 0.92-0.03);
  for(int j=0;j<3;j++){
    grPdep2[j]->Draw("sameEP");
    tl->AddEntry(grPdep2[j], Form("#it{p} for #it{h_{E}}(%d)", j), "P");
    tl->SetBorderSize(0);
    tl->SetTextSize(0.045);
  }
  tl->Draw();
  c0->Print(Form("%s", PDFName.c_str()), "pdf");

  hFrame->Draw();
  tl = new TLegend(0.58-0.38, 0.68+0.02, 0.92-0.48, 0.92-0.03);
  for(int j=0;j<2;j++){
    grPdep3[j]->Draw("sameEP");
    tl->AddEntry(grPdep3[j], Form("#it{p} for #it{h_{E}}(%d)", j), "P");
    tl->SetBorderSize(0);
    tl->SetTextSize(0.045);
  }
  tl->Draw();
  c0->Print(Form("%s", PDFName.c_str()), "pdf");

  cout << "hoge" << endl;
  for(int iStep=0;iStep<nStep;iStep++){
    for(int j=0;j<nObs;j++){
      double vObs = vRaw[j][iStep+1];
      double eObs = eRaw[j][iStep+1];
      double vP = vPdep2[j][iStep];
      double eP = ePdep2[j][iStep];
      
      vQnx_dep2[j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
      eQnx_dep2[j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				 (TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));

      grQnx_dep2[j]->SetPoint(iStep, dt*(iStep+1), vQnx_dep2[j][iStep]);
      grQnx_dep2[j]->SetPointError(iStep, 0.0, eQnx_dep2[j][iStep]);
    }// for j

    for(int j=0;j<nObs;j++){
      double vObs = vRaw[j][iStep+1];
      double eObs = eRaw[j][iStep+1];
      double vP = vPdep2_fit[j][iStep];
      double eP = ePdep2_fit[j][iStep];
      
      vQnx_dep2_fit[j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
      eQnx_dep2_fit[j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				(TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));

      grQnx_dep2_fit[j]->SetPoint(iStep, dt*(iStep+1), vQnx_dep2_fit[j][iStep]);
      grQnx_dep2_fit[j]->SetPointError(iStep, 0.0, eQnx_dep2_fit[j][iStep]);
    }// for j

    for(int j=0;j<nObs;j++){
      double vObs = vRaw[j][iStep+1];
      double eObs = eRaw[j][iStep+1];
      double vP = vPdep3_fit[j][iStep];
      double eP = ePdep3_fit[j][iStep];
      
      if(j==2){
	vP = vPdep3_fit[0][iStep];
	eP = ePdep3_fit[0][iStep];
	vQnx_dep3_fit[j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
	eQnx_dep3_fit[j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				       (TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));

	grQnx_dep3_fit[j]->SetPoint(iStep, dt*(iStep+1), vQnx_dep3_fit[j][iStep]);
	grQnx_dep3_fit[j]->SetPointError(iStep, 0.0, eQnx_dep3_fit[j][iStep]);
      }else{
	vQnx_dep3_fit[j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
	eQnx_dep3_fit[j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				       (TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));

	grQnx_dep3_fit[j]->SetPoint(iStep, dt*(iStep+1), vQnx_dep3_fit[j][iStep]);
	grQnx_dep3_fit[j]->SetPointError(iStep, 0.0, eQnx_dep3_fit[j][iStep]);
	
      }
    }// for j
    
  }// for iStep
  
  for(int j=0;j<4;j++){
    cout << "j = " << j << endl;
    delete hFrame; hFrame=NULL;
    hFrame = new TH2D("hFrame", ";Time;Electric Energy Density", 50, -0.05, 7.55, 50, 0.0, 2.0);
    hFrame->GetYaxis()->SetNdivisions(408);
    hFrame->GetXaxis()->SetNdivisions(416);
    hFrame->Draw();
    tl = new TLegend(0.58-0.35, 0.68, 0.92-0.45, 0.92);
    tl->AddEntry(grSim[j]  , "Noiseless", "l");
    tl->AddEntry(grRaw[j], "Raw", "P");
    tl->AddEntry(grQnx_dep1[j], "EM Depolarizing (global)", "P");
    tl->AddEntry(grQnx_dep2[j], "EM Depolarizing (zero)", "P");
    tl->AddEntry(grQnx_dep2_fit[j], "EM Depolarizing (zero-fit)", "P");
    tl->AddEntry(grQnx_dep3_fit[j], "EM Depolarizing (pair-fit)", "P");
    tl->AddEntry(grPost[j], "EM Post-selection", "P");
    tl->SetBorderSize(0);
    tl->SetTextSize(0.03);

    ttx = new TText(0.85, 0.965, Form("link %d", j));
    ttx->SetNDC();
    ttx->SetTextSize(0.040);
    ttx->Draw();

    tl->Draw();

    grSim[j]->SetLineWidth(2);
    grSim[j]->Draw("sameEL");
    grRaw[j]->Draw("sameEP");
    grQnx_dep1[j]->Draw("sameEP");
    grQnx_dep2[j]->Draw("sameEP");
    grQnx_dep2_fit[j]->Draw("sameEP");
    grQnx_dep3_fit[j]->Draw("sameEP");
    grPost[j]->Draw("sameEP");

    
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
    delete ttx; ttx=NULL;
    delete tl; tl=NULL;
  }

  for(int j=0;j<4;j++){
    cout << "j = " << j << endl;
    delete hFrame; hFrame=NULL;
    hFrame = new TH2D("hFrame", ";Time;Electric Energy Density", 50, -0.05, 7.55, 50, 0.0, 2.0);
    hFrame->GetYaxis()->SetNdivisions(408);
    hFrame->GetXaxis()->SetNdivisions(416);
    hFrame->Draw();
    tl = new TLegend(0.58-0.3, 0.68, 0.92-0.4, 0.92);
    tl->AddEntry(grSim[j]  , "Noiseless", "l");
    tl->AddEntry(grRaw[j], "Raw", "P");
    tl->SetBorderSize(0);
    tl->SetTextSize(0.045);
    //tl->SetTextSize(0.04);

    ttx = new TText(0.85, 0.965, Form("link %d", j));
    ttx->SetNDC();
    ttx->SetTextSize(0.040);
    ttx->Draw();

    tl->Draw();

    grSim[j]->SetLineWidth(2);
    grSim[j]->Draw("sameEL");
    grRaw[j]->Draw("sameEP");
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
    delete ttx; ttx=NULL;
    delete tl; tl=NULL;
  }  
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

  c0 = new TCanvas("c0", "c0", 800, 600);
  hAll  = new TH1D("hAll", "", nStep+1, -0.5, nStep+0.5);
  hPass = new TH1D("hPass", "", nStep+1, -0.5, nStep+0.5);
  f1 = new TF1("f1", "[0]*x", -0.5, 6.5);
  f1->SetNpx(1000);
  
  for(int i=0;i<nObs;i++){
    hSim[i] = new TH1D(Form("hSim_%d", i), "", 1000, 0, 10);
    hRaw[i] = new TH1D(Form("hRaw_%d", i), "", 1000, 0, 10);
    hPost[i] = new TH1D(Form("hPost_%d", i), "", 1000, 0, 10);

    grSim[i] = new TGraphErrors();
    grRaw[i] = new TGraphErrors();
    grPost[i] = new TGraphErrors();

    grPdep1[i] = new TGraphErrors();
    grPdep1[i]->SetMarkerStyle(8);
    if(i!=4){
      grPdep1[i]->SetLineColor(MyCol[i]);
      grPdep1[i]->SetMarkerColor(MyCol[i]);
    }else{
      grPdep1[i]->SetLineColor(kBlack);
      grPdep1[i]->SetMarkerColor(kBlack);
    }

    grPdep2[i] = new TGraphErrors(); 
    //grPdep2[i]->SetMarkerStyle(8);
    grPdep2[i]->SetMarkerStyle(21+i);
   if(i!=4){
      grPdep2[i]->SetLineColor(MyCol[i]);
      grPdep2[i]->SetMarkerColor(MyCol[i]);
    }else{
      grPdep2[i]->SetLineColor(kBlack);
      grPdep2[i]->SetMarkerColor(kBlack);
    }

    grPdep3[i] = new TGraphErrors();
    //grPdep3[i]->SetMarkerStyle(8);
    grPdep3[i]->SetMarkerStyle(21+i);
    if(i!=4){
      grPdep3[i]->SetLineColor(MyCol[i]);
      grPdep3[i]->SetMarkerColor(MyCol[i]);
    }else{
      grPdep3[i]->SetLineColor(kBlack);
      grPdep3[i]->SetMarkerColor(kBlack);
    }
    
    
    grSim[i]->SetLineColor(kBlack);
    grSim[i]->SetLineStyle(2);
    grRaw[i]->SetLineColor(kBlack);
    grPost[i]->SetLineColor(kRed);
    grRaw[i]->SetMarkerColor(kBlack);
    grPost[i]->SetMarkerColor(kRed);
    grRaw[i]->SetMarkerStyle(20);
    grPost[i]->SetMarkerStyle(22);

    grQnx_dep1[i] = new TGraphErrors();
    grQnx_dep1[i]->SetMarkerStyle(21);
    grQnx_dep1[i]->SetMarkerColor(kBlue);
    grQnx_dep1[i]->SetLineColor(kBlue);

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

  for(int iObs=0;iObs<nObs;iObs++){
    for(int iStep=0;iStep<nStep+1;iStep++){
      grSim[iObs]->SetPoint(iStep, dt*iStep, vSim[iObs][iStep]);
    }// for iStep
  }// for iObs

  return 0;
}

int LoadMainResults(void){
  hPass->Reset();
  hAll->Reset();
  for(int iStep=0;iStep<nStep+1;iStep++){
    for(int j=0;j<nObs;j++){
      hRaw[j]->Reset();
      hPost[j]->Reset();
    }

    //string filename = Form("H1-Emulator_%s_%dstep.txt", Tag.c_str(), iStep);
    //string filename = Form("H2-Emulator_%s_%dstep.txt", Tag.c_str(), iStep);
    //string filename = Form("H2-Emulator_%s_%dstep.txt", Tag.c_str(), iStep);
    string filename = Form("../savefiles/H2-Emulator_%s_%dstep.txt", Tag.c_str(), iStep);
    //string filename = Form("%s_%dstep.txt", Tag.c_str(), iStep);
    ifs.open(filename.c_str());
    cout << "Load File: " << filename << endl;
    
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
	    hAll->Fill(iStep);
	  
	  if(anc==0 && N_Rfix==0 && N_L0==N_R0 && N_L1==N_R1 && N_L2==N_R2){
	    hPost[j]->Fill(vIn[j]);
	  }else{
	    if(j==0)
	      hPass->Fill(iStep);
	  }
	} // for counts
      }// for j
    }// while
    ifs.close();
    
    for(int j=0;j<nObs;j++){
      vRaw[j][iStep] = hRaw[j]->GetMean();
      eRaw[j][iStep] = hRaw[j]->GetMeanError();
      grRaw[j]->SetPoint(iStep, dt*iStep, vRaw[j][iStep]);
      grRaw[j]->SetPointError(iStep, 0.0, eRaw[j][iStep]);

      vPost[j][iStep] = hPost[j]->GetMean();
      ePost[j][iStep] = hPost[j]->GetMeanError();
      grPost[j]->SetPoint(iStep, dt*iStep, vPost[j][iStep]);
      grPost[j]->SetPointError(iStep, 0.0, ePost[j][iStep]);

      vPdep1[j][iStep]= 1.0 - (vRaw[j][iStep] - ccc[j])/(vSim[j][iStep] - ccc[j]);
      ePdep1[j][iStep]= TMath::Abs(eRaw[j][iStep]/(vSim[j][iStep] - ccc[j]));
    }
    for(int j=0;j<nObs;j++){
      double vObs = vRaw[j][iStep];
      double eObs = eRaw[j][iStep];
      double vP = vPdep1[3][iStep];
      double eP = ePdep1[3][iStep];
      cout << "iStep = " << iStep << ", j = " << j << ", vObs =  " << vObs << endl;
      vQnx_dep1[j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
      eQnx_dep1[j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				(TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));

      grQnx_dep1[j]->SetPoint(iStep, dt*iStep, vQnx_dep1[j][iStep]);
      grQnx_dep1[j]->SetPointError(iStep, 0.0, eQnx_dep1[j][iStep]);
    }
    //cout << iStep << " trotter step: ";
    //cout << 100.0*(1.0 - hCond1[0]->GetEntries()/hCond0[0]->GetEntries()) << ", ";
    //cout << 100.0*(1.0 - hCond2[0]->GetEntries()/hCond0[0]->GetEntries()) << ", ";
    //cout << 100.0*(1.0 - hCond3[0]->GetEntries()/hCond0[0]->GetEntries()) << endl;
    
    vDiscard[iStep] = 1.0 - hPost[0]->GetEntries()/hRaw[0]->GetEntries();
  }// for iStep
  cout << "Discard rate: ";
  for(int i=0;i<nStep+1;i++){
    cout << vDiscard[i] << ", ";
  }
  cout << endl;

  grDiscard = new TGraphAsymmErrors(hPass, hAll);
  //grDiscard = new TGraphAsymmErrors(hCond3[0], hAll);

  hPass->Reset();
  hAll->Reset();

  cout << "Depolarizing noise probability" << endl;
  for(int i=0;i<nStep+1;i++){
    for(int j=0;j<nObs;j++){
      cout << Form("%.3lf +- %.3lf,   ", vPdep1[j][i], ePdep1[j][i]);
      grPdep1[j]->SetPoint(i, i, vPdep1[j][i]);
      grPdep1[j]->SetPointError(i, 0.0, ePdep1[j][i]);
    }
    cout << endl;
  }

  
  return 0;
}


int LoadZeroInit(void){
  for(int iStep=0;iStep<nStep;iStep++){
    for(int j=0;j<nObs;j++){
      hRaw[j]->Reset();
    }

    //string filename = Form("H1-Emulator_%s_zeroinit_%dstep.txt", Tag.c_str(), iStep+1);
    //string filename = Form("H2-Emulator_%s_zeroinit_%dstep.txt", Tag.c_str(), iStep+1);
    string filename = Form("../savefiles/H2-Emulator_%s_zeroinit_%dstep.txt", Tag.c_str(), iStep+1);
    //string filename = Form("%s_%dstep.txt", Tag.c_str(), iStep);
    ifs.open(filename.c_str());
    cout << "Load File: " << filename << endl;
    
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

      vPdep2[j][iStep]= 1.0 - (vRaw[j][iStep] - ccc[j])/(0.0 - ccc[j]);
      ePdep2[j][iStep]= TMath::Abs(eRaw[j][iStep]/(0.0 - ccc[j]));
    }
    for(int j=0;j<nObs;j++){
      double vObs = vRaw[j][iStep];
      double eObs = eRaw[j][iStep];
      double vP = vPdep2[j][iStep];
      double eP = ePdep2[j][iStep];
      cout << "iStep = " << iStep+1 << ", j = " << j << ", vObs =  " << vObs << endl;
      /*
      vQnx_dep2[j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
      eQnx_dep2[j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				(TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));

      grQnx_dep2[j]->SetPoint(iStep, dt*(iStep+1), vQnx_dep2[j][iStep]);
      grQnx_dep2[j]->SetPointError(iStep, 0.0, eQnx_dep2[j][iStep]);
      */
      grPdep2[j]->SetPoint(iStep, iStep+1, vPdep2[j][iStep]);
      grPdep2[j]->SetPointError(iStep, 0.0, ePdep2[j][iStep]);
    }
  }// for iStep

  for(int j=0;j<nObs;j++){
    grPdep2[j]->Fit(f1, "0");
    for(int iStep=0;iStep<nStep;iStep++){
      vPdep2_fit[j][iStep] = (iStep+1)*grPdep2[j]->GetFunction("f1")->GetParameter(0);
      ePdep2_fit[j][iStep] = (iStep+1)*grPdep2[j]->GetFunction("f1")->GetParError(0);
    }// for j
  }// for iStep

  return 0;
}
int LoadPairResults(void){
  cout << "Load Pair Results" << endl;
  for(int iStep=0;iStep<nStep/2;iStep++){
    for(int j=0;j<nObs;j++){
      hRaw[j]->Reset();
    }

    //string filename = Form("H1-Emulator_%s_pair_%dstep.txt", Tag.c_str(), 2*(iStep+1));
    //string filename = Form("H2-Emulator_%s_pair_%dstep.txt", Tag.c_str(), 2*(iStep+1));
    string filename = Form("../savefiles/H2-Emulator_%s_pair_%dstep.txt", Tag.c_str(), 2*(iStep+1));
    //string filename = Form("%s_%dstep.txt", Tag.c_str(), iStep);
    ifs.open(filename.c_str());
    cout << "Load File: " << filename << endl;
    
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

      vPdep3[j][iStep]= 1.0 - (vRaw[j][iStep] - ccc[j])/(vSim[j][0] - ccc[j]);
      ePdep3[j][iStep]= TMath::Abs(eRaw[j][iStep]/(vSim[j][0] - ccc[j]));
      cout <<  "iStep = " << 2*(iStep+1) << ", j = " << j << ", vRaw = " << vRaw[j][iStep] << ", vSim = " << vSim[j][0] << ", ccc = " << ccc[j] << endl;
    }
    for(int j=0;j<nObs;j++){
      double vObs = vRaw[j][iStep];
      double eObs = eRaw[j][iStep];
      double vP = vPdep3[j][iStep];
      double eP = ePdep3[j][iStep];
      cout << "iStep = " << 2*(iStep+1) << ", j = " << j << ", vObs =  " << vObs << ", vP = " << vP << ", eP = " << eP << endl;
      /*
      vQnx_dep3[j][iStep] = (vObs - ccc[j])/(1.0 - vP) + ccc[j];
      eQnx_dep3[j][iStep] = sqrt(TMath::Power((eObs)/(1.0 - vP), 2) +
				(TMath::Power(vObs - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));

      grQnx_dep3[j]->SetPoint(iStep, dt*2*(iStep+1), vQnx_dep3[j][iStep]);
      grQnx_dep3[j]->SetPointError(iStep, 0.0, eQnx_dep3[j][iStep]);
      */
      grPdep3[j]->SetPoint(iStep, 2*(iStep+1), vPdep3[j][iStep]);
      grPdep3[j]->SetPointError(iStep, 0.0, ePdep3[j][iStep]);
    }
  }// for iStep

  for(int j=0;j<nObs;j++){
    grPdep3[j]->Fit(f1, "0");
    for(int iStep=0;iStep<nStep;iStep++){
      vPdep3_fit[j][iStep] = (iStep+1)*grPdep3[j]->GetFunction("f1")->GetParameter(0);
      ePdep3_fit[j][iStep] = (iStep+1)*grPdep3[j]->GetFunction("f1")->GetParError(0);
    }// for j
  }// for iStep

  
  return 0;
}
