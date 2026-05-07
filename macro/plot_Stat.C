bool fEmulator=true;
bool fUpdateFigure=true;
TLegend *tl;
TText *ttx;
TText *txP;
TLine *myline;
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
TGraph *grRemain;
TH2D *hFrame;
TH2D *hRatioFrame;
TPad *pad1[2];
const int MyCol[4] = {EColor::kRed, EColor::kOrange, EColor::kGreen+1, EColor::kBlue};

const int nObs=5;
const int nStep=6;
//const int nStep=4;
TH1D *hSim[nObs];
TH1D *hRaw[nStep][nObs];
TH1D *hZero[nStep][nObs];
TH1D *hPost[nObs];

TGraphErrors *grObs_Raw[nStep][nObs];
TGraphErrors *grObs_Dep[nStep][nObs];

TH1D *hDiff[3];
TH1D *hDiff_Sigma[3];
TH1D *hDiff_Ratio[3];

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
TGraphErrors *grRaw_Ratio[nObs];
TGraphErrors *grPost_Ratio[nObs];
TGraphErrors *grQnx_dep2_Ratio[nObs];

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
//const string PDFName=Form("../PDF/PaperPlot_H2_%s_zero100.pdf", Tag.c_str());
//const string PDFName=Form("../PDF/PaperPlot_H2_%s_main100_zero100.pdf", Tag.c_str());
//const string PDFName=Form("../PDF/PaperPlot_H2_%s_main100.pdf", Tag.c_str());
//const string PDFName=Form("../PDF/PaperPlot_H2_%s_var.pdf", Tag.c_str());
const string PDFName=Form("../PDF/plot_Stat.pdf");
//const string PDFName=Form("../PDF/PaperPlot_H1_%s.pdf", Tag.c_str());
int InitSetup(void);
int LoadMainResults(void);
int LoadZeroInit(void);

ifstream ifsZero[nStep];
ifstream ifsMain[nStep];

const int nShotsStep=100;
const int nShotsMax=6000;

//string filename = Form("../savefiles/H2-Emulator_%s_zeroinit_%dstep.txt", Tag.c_str(), iStep+1);

int plot_Stat(void){
  // Initial Setup
  InitSetup();

  for(int iStep=0;iStep<nStep;iStep++){
    //ifsZero[iStep].open(Form("../savefiles/H2-Emulator_shots_zeroinit_0_%dstep.txt", iStep+1));
    //ifsMain[iStep].open(Form("../savefiles/H2-Emulator_shots_0_%dstep.txt", iStep+1));
    ifsZero[iStep].open(Form("../savefiles/H2-Emulator_shots_zeroinit_0_%dstep.txt", iStep+1));
    ifsMain[iStep].open(Form("../savefiles/H2-Emulator_shots_0_%dstep.txt", iStep+1));
  }
  int iFile=0;
  for(int iP=0;iP<(nShotsMax/nShotsStep);iP++){
    if(iP>((1000/nShotsStep)*(iFile+1)-1)){
      iFile++;
      for(int iStep=0;iStep<nStep;iStep++){
	ifsZero[iStep].close();
	ifsMain[iStep].close();
	ifsZero[iStep].open(Form("../savefiles/H2-Emulator_shots_zeroinit_%d_%dstep.txt", iFile, iStep+1));
	ifsMain[iStep].open(Form("../savefiles/H2-Emulator_shots_%d_%dstep.txt", iFile, iStep+1));
      }// for iStep   
    }// if
      
    for(int iShots=0;iShots<nShotsStep;iShots++){
      for(int iStep=0;iStep<nStep;iStep++){
	//cout << "iP = " << iP << ", iShots = " << iShots << endl;
	// zero init
	ifsZero[iStep]>>N_L0>>N_L1>>N_L2>>N_Lfix>>N_R0>>N_R1>>N_R2>>N_Rfix>>anc>>counts;
	//cout << counts << endl;
	vIn[0] = (N_L0/2.0)*(N_L0/2.0 + 1);
	vIn[1] = (N_L1/2.0)*(N_L1/2.0 + 1);
	vIn[2] = (N_L2/2.0)*(N_L2/2.0 + 1);
	vIn[3] = (N_Lfix/2.0)*(N_Lfix/2.0 + 1);
	vIn[4] = (N_Rfix/2.0)*(N_Rfix/2.0 + 1);
      
	for(int j=0;j<nObs;j++){
	  hZero[iStep][j]->Fill(vIn[j]);
	}// for j
	
	// main
	ifsMain[iStep]>>N_L0>>N_L1>>N_L2>>N_Lfix>>N_R0>>N_R1>>N_R2>>N_Rfix>>anc>>counts;
	vIn[0] = (N_L0/2.0)*(N_L0/2.0 + 1);
	vIn[1] = (N_L1/2.0)*(N_L1/2.0 + 1);
	vIn[2] = (N_L2/2.0)*(N_L2/2.0 + 1);
	vIn[3] = (N_Lfix/2.0)*(N_Lfix/2.0 + 1);
	vIn[4] = (N_Rfix/2.0)*(N_Rfix/2.0 + 1);
      
	for(int j=0;j<nObs;j++){
	  hRaw[iStep][j]->Fill(vIn[j]);
	}// for j
      }// for iStep
    }// for iShots

    for(int iStep=0;iStep<nStep;iStep++){
      //double vP = 1.0 - (hRaw[iStep][3]->GetMean() - ccc[3])/(vSim[3][iStep+1] - ccc[3]);
      //double eP = TMath::Abs(hRaw[iStep][3]->GetMeanError()/(vSim[3][iStep+1] - ccc[3]));
      for(int j=0;j<nObs;j++){
	double vRaw = hRaw[iStep][j]->GetMean();
	double eRaw = hRaw[iStep][j]->GetMeanError();
	//cout << "vRaw = " << vRaw << ", eRaw = " << eRaw << endl;
	
	grObs_Raw[iStep][j]->SetPoint(iP, nShotsStep*(iP+1), vRaw);
	grObs_Raw[iStep][j]->SetPointError(iP, 0.0, eRaw);

	double vP = 1.0 - (hZero[iStep][j]->GetMean() - ccc[j])/(0.0 - ccc[j]);
	double eP = TMath::Abs(hZero[iStep][j]->GetMeanError()/(0.0 - ccc[j]));
	//cout << "vP = " << vP << ", eP = " << eP << endl;

	double vDep =  (vRaw - ccc[j])/(1.0 - vP) + ccc[j];
	double eDep = sqrt(TMath::Power((eRaw)/(1.0 - vP), 2) +
			   (TMath::Power(vRaw - ccc[j], 2)*eP*eP)/(TMath::Power(1.0 - vP, 4)));
	
	grObs_Dep[iStep][j]->SetPoint(iP, nShotsStep*(iP+1), vDep);
	grObs_Dep[iStep][j]->SetPointError(iP, 0.0, eDep);	
      }// for j
    }// for iStep
  }// for iP
  
  hFrame = new TH2D("hFrame", ";Number of Shots;Electric Energy Density", 50, 0.0, nShotsMax+10, 50, 0.0, 2.0);
  hFrame->GetYaxis()->SetNdivisions(408);
  hFrame->GetXaxis()->SetNdivisions(7);
  c0->Print(Form("%s[", PDFName.c_str()), "pdf");

  for(int iStep=0;iStep<nStep;iStep++){
    if(iStep==0) hFrame->GetYaxis()->SetRangeUser(0.0, 0.7);
    if(iStep==1) hFrame->GetYaxis()->SetRangeUser(0.0, 1.0);
    if(iStep==2) hFrame->GetYaxis()->SetRangeUser(0.0, 2.0);
    if(iStep==3) hFrame->GetYaxis()->SetRangeUser(0.0, 1.0);
    if(iStep==4) hFrame->GetYaxis()->SetRangeUser(0.0, 1.5);
    if(iStep==5) hFrame->GetYaxis()->SetRangeUser(0.0, 1.5);
    hFrame->Draw();
    ttx = new TText(0.75, 0.965, Form("%d-Trotter-step", iStep+1));
    ttx->SetNDC();
    ttx->SetTextSize(0.050);
    ttx->Draw();

    if(fEmulator) txP->Draw();
    tl = new TLegend(0.5, 0.68+0.02+0.15, 0.92, 0.92-0.03);
    tl->SetNColumns(3);
    tl->SetBorderSize(0);
    tl->SetTextSize(0.045);
    for(int j=0;j<3;j++){
      myline = new TLine(0.0, vSim[j][iStep+1], nShotsMax, vSim[j][iStep+1]);
      myline->SetLineStyle(2);
      //myline->SetLineColor(60+4*j);
      myline->SetLineColor(MyCol[j]);
      myline->Draw("same");
      grObs_Raw[iStep][j]->SetMarkerStyle(20+j);
      grObs_Raw[iStep][j]->SetMarkerColor(91+4*j);
      grObs_Raw[iStep][j]->SetLineColor(91+4*j);
      //grObs_Raw[iStep][j]->Draw("sameP");
      grObs_Dep[iStep][j]->SetMarkerStyle(24+j);
      //grObs_Dep[iStep][j]->SetMarkerColor(60+4*j);
      //grObs_Dep[iStep][j]->SetLineColor(60+4*j);
      grObs_Dep[iStep][j]->SetMarkerColor(MyCol[j]);
      grObs_Dep[iStep][j]->SetLineColor(MyCol[j]);
      grObs_Dep[iStep][j]->Draw("sameP");
      tl->AddEntry(grObs_Dep[iStep][j], Form("link-%d", j), "l");
    }
    tl->Draw();
    c0->Print(Form("%s", PDFName.c_str()), "pdf");
  }// for iStep  
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
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
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

  txP = new TText(0.18, 0.965, "Emulator results");
  txP->SetNDC();
  txP->SetTextSize(0.040);
  txP->SetTextColor(kRed);

  c0 = new TCanvas("c0", "c0", 800, 600);

  for(int iStep=0;iStep<nStep;iStep++){
    for(int i=0;i<nObs;i++){
      hZero[iStep][i] = new TH1D(Form("hZero_%dstep_%d", iStep, i), "", 1000, 0, 10);
      hRaw[iStep][i] = new TH1D(Form("hRaw_%dstep_%d", iStep, i), "", 1000, 0, 10);

      grObs_Raw[iStep][i] = new TGraphErrors();
      grObs_Dep[iStep][i] = new TGraphErrors();
    }
  }// for iStep

  return 0;
}

