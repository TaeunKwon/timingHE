#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()

#include "TSystem.h"
#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TGaxis.h"
#include <algorithm>


#include <string>

#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;


///Define color vector
TColor *one = new TColor(2001, 0.906, 0.153, 0.094);
TColor *two = new TColor(2002, 0.906, 0.533, 0.094);
TColor *three = new TColor(2003, 0.086, 0.404, 0.576);
TColor *four = new TColor(2004, 0.071, 0.694, 0.18);
TColor *five = new TColor(2005, 0.388, 0.098, 0.608);
TColor *six = new TColor(2006, 0.906, 0.878, 0.094);
const int colors[] = {1, 2001, 2003, 2004, 2005, 2006, 2002, 2, 3, 4, 6, 7, 5, 1, 8, 9, 29, 38, 46};


void plot_TH1F(vector<TH1F*> h1_vec, TString name, int wedge);
void plot_vecTH1F(vector< vector<TH1F*> > h1_vec, TString name, int wedge);
void plot_TH2F(vector<TH2F*> h2_vec, TString name, int wedge);
void plot_pulse(vector<TH1F*> h1_vec, vector<TH1F*> h1_ped, TString name, int wedge);
void plot_pulse(vector<TH1F*> h1_vec, TString name, int wedge);

void printSummaryScatter(TGraphErrors g, TString name, TString xlabel, TString ylabel);
void printSimpleTH1F(TH1F * h, TString name);
void printSummaryTH2(TH2F * h2, TString name, float min = 0, float max = 0);
void printSummaryTH1vec(vector<TH1F*> h1, TString name, float minn = 0, float maxx = 0, bool logY = false);
void writeCorrectionTable(TString name, int ifC = 2);
void h1cosmetic(TH1F *hist);
void setSummaryStyle();
string getPhi(TString title);
int getWedge(int iPhi) {
    if (iPhi >= 71 || iPhi <= 2) return 1;
    else return floor((iPhi + 5) / 4.);
}

const int minieta = -29;
const int nieta = 59;

//map IPhi 1 to 72 to histogram indices 0 to 71 (index = IPhi - miniphi)
const int miniphi = 1;
const int niphi = 72;

//map Depth 1 to 7 to histogram indices 0 to 6 (index = depth - mindepth)
const int mindepth = 1;
const int ndepth = 7;

const int nWedges = 18;
const int nfC = 4;
const int fCbins[] = {0, 5000, 7000, 11000};
TString plotDir;
TString summaryDir;
TString webplotDir;

//Prepare vectors to store various properties (one vector for each charge bin)
//Each vector will have one entry per channel, with matching indices
vector<int> entries[nfC];
vector<int> etas[nfC];
vector<int> phis[nfC];
vector<int> depths[nfC];

vector<float> meanTDCs[nfC];
vector<float> meanErrTDCs[nfC];
vector<float> rmsTDCs[nfC];
vector<float> medianTDCs[nfC];
vector<float> medianErrTDCs[nfC];

vector<float> meanTS2s[nfC];
vector<float> meanErrTS2s[nfC];
vector<float> meanTS4s[nfC];
vector<float> meanErrTS4s[nfC];
vector<float> meanEs[nfC];
vector<float> meanErrEs[nfC];

//Eta-phi maps, one per depth and charge bin
TH2F * h2_entries[nfC][ndepth];
TH2F * h2_meanTDC[nfC][ndepth];
TH2F * h2_rmsTDC[nfC][ndepth];
TH2F * h2_medianTDC[nfC][ndepth];
TH2F * h2_mean_medianTDC[nfC][ndepth];
TH2F * h2_sigma[nfC][ndepth];
TH2F * h2_meanErrTDC[nfC][ndepth];
TH2F * h2_meanTS2[nfC][ndepth];
TH2F * h2_meanTS4[nfC][ndepth];
TH2F * h2_meanE[nfC][ndepth];
TH2F * h2_correction[nfC][ndepth];

//Simple histogram filled with value for each channel
//ndepth+1 reserved for combination of ALL depths
TH1F * h1_entries[nfC][ndepth + 1];
TH1F * h1_meanTDC[nfC][ndepth + 1];
TH1F * h1_rmsTDC[nfC][ndepth + 1];
TH1F * h1_medianTDC[nfC][ndepth + 1];
TH1F * h1_mean_medianTDC[nfC][ndepth + 1];
TH1F * h1_sigma[nfC][ndepth + 1];
TH1F * h1_meanErrTDC[nfC][ndepth + 1];
TH1F * h1_meanTS2[nfC][ndepth + 1];
TH1F * h1_meanTS4[nfC][ndepth + 1];
TH1F * h1_meanE[nfC][ndepth + 1];
TH1F * h1_correction[nfC][ndepth + 1];

//Mean shape vs eta, averaged over all phi, in slices of depth
TH1F * h1_entries_vs_eta[nfC][ndepth];
TH1F * h1_meanTDC_vs_eta[nfC][ndepth];
TH1F * h1_rmsTDC_vs_eta[nfC][ndepth];
TH1F * h1_medianTDC_vs_eta[nfC][ndepth];
TH1F * h1_mean_medianTDC_vs_eta[nfC][ndepth];
TH1F * h1_sigma_vs_eta[nfC][ndepth];
TH1F * h1_meanErrTDC_vs_eta[nfC][ndepth];
TH1F * h1_meanTS2_vs_eta[nfC][ndepth];
TH1F * h1_meanTS4_vs_eta[nfC][ndepth];
TH1F * h1_meanE_vs_eta[nfC][ndepth];

//Mean shape vs phi, averaged over all eta, in slices of depth. Separated for HE- and HE+
//HEM
TH1F * h1_entries_vs_phi_HEM[nfC][ndepth];
TH1F * h1_meanTDC_vs_phi_HEM[nfC][ndepth];
TH1F * h1_rmsTDC_vs_phi_HEM[nfC][ndepth];
TH1F * h1_medianTDC_vs_phi_HEM[nfC][ndepth];
TH1F * h1_mean_medianTDC_vs_phi_HEM[nfC][ndepth];
TH1F * h1_sigma_vs_phi_HEM[nfC][ndepth];
TH1F * h1_meanErrTDC_vs_phi_HEM[nfC][ndepth];
TH1F * h1_meanTS2_vs_phi_HEM[nfC][ndepth];
TH1F * h1_meanTS4_vs_phi_HEM[nfC][ndepth];
TH1F * h1_meanE_vs_phi_HEM[nfC][ndepth];

//HEP
TH1F * h1_entries_vs_phi_HEP[nfC][ndepth];
TH1F * h1_meanTDC_vs_phi_HEP[nfC][ndepth];
TH1F * h1_rmsTDC_vs_phi_HEP[nfC][ndepth];
TH1F * h1_medianTDC_vs_phi_HEP[nfC][ndepth];
TH1F * h1_mean_medianTDC_vs_phi_HEP[nfC][ndepth];
TH1F * h1_sigma_vs_phi_HEP[nfC][ndepth];
TH1F * h1_meanErrTDC_vs_phi_HEP[nfC][ndepth];
TH1F * h1_meanTS2_vs_phi_HEP[nfC][ndepth];
TH1F * h1_meanTS4_vs_phi_HEP[nfC][ndepth];
TH1F * h1_meanE_vs_phi_HEP[nfC][ndepth];

float getMedian(TH1F * h) {
    Double_t median;
    Double_t fiftypercent = 0.5;
    h->GetQuantiles(1, &median, &fiftypercent);
    return median;
}

void plotDistributions(TString histFile, TString tag, int wedgeSel = 0, int mode = 0) {

    /// Wedge sel: do one wedge at a time, for testing or quick checks

    ///mode:
    //0: all plots
    //1: Only per-channel distributions and summary hists
    //2: Only 1D summaries, collapsed in eta / phi


//vector< <vector<TH2F*> > ADCvsTDC_fibers;
    plotDir = Form("plots/%s/", tag.Data());
    summaryDir = Form("summary/%s/", tag.Data());
    webplotDir = Form("webplots/%s/", tag.Data());
    gSystem->mkdir(plotDir, kTRUE);
    gSystem->mkdir(summaryDir, kTRUE);
    gSystem->mkdir(webplotDir, kTRUE);
    gSystem->mkdir(webplotDir + "summary/", kTRUE);
    TFile * inFile = TFile::Open(histFile);

//1D shapes: vector for each wedge, eta [depth], but integrated over phi
//1D charge, time: vector of vectors for each wedge,eta ([depth][phi])
//2D time slew, integrated over phi
    for (int ifC = 0; ifC < nfC; ifC++) {
        for (int depth = 0; depth < ndepth + 1; depth++) {
            h1_entries[ifC][depth] = new TH1F(Form("h1_entries_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Log10(Entries);Channels", depth + 1, ifC) , 40, 0, 7);
            h1_meanTDC[ifC][depth] = new TH1F(Form("h1_meanTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Mean TDC time [ns];Channels", depth + 1, ifC) , 32, 72.75, 88.75);
            h1_rmsTDC[ifC][depth] = new TH1F(Form("h1_rmsTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;RMS in TDC time [ns];Channels", depth + 1, ifC) , 40, 0, 12);
            h1_medianTDC[ifC][depth] = new TH1F(Form("h1_medianTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Median TDC time [ns];Channels", depth + 1, ifC) , 32, 72.75, 88.75);
            h1_mean_medianTDC[ifC][depth] = new TH1F(Form("h1_mean_medianTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Mean - median TDC time [ns];Channels", depth + 1, ifC) , 30, -4, 4);
            h1_sigma[ifC][depth] = new TH1F(Form("h1_sigma_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i; Pull from 81 ns [sigma];Channels", depth + 1, ifC) , 30, -5, 5);
            h1_meanErrTDC[ifC][depth] = new TH1F(Form("h1_meanErrTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Error on mean TDC time [ns];Channels", depth + 1, ifC) , 40, 0, 1.5);
            h1_meanTS2[ifC][depth] = new TH1F(Form("h1_meanTS2_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Mean TS2 fraction [%];Channels", depth + 1, ifC) , 30, 0, 2);
            h1_meanTS4[ifC][depth] = new TH1F(Form("h1_meanTS4_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Mean TS4 fraction [%];Channels", depth + 1, ifC) , 40, 20, 80);
            h1_meanE[ifC][depth] = new TH1F(Form("h1_meanE_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Mean total charge [fC];Channels", depth + 1, ifC) , 40, 5000, 15000);
            h1_correction[ifC][depth] = new TH1F(Form("h1_correction_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;Correction to 81 [ns];Channels", depth + 1, ifC) , 40, -15.5, 4.5);
            if (depth < ndepth) {
                h2_entries[ifC][depth] = new TH2F(Form("h2_entries_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Entries", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_meanTDC[ifC][depth] = new TH2F(Form("h2_meanTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Mean TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_rmsTDC[ifC][depth] = new TH2F(Form("h2_rmsTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;RMS in TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_medianTDC[ifC][depth] = new TH2F(Form("h2_medianTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Median TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_mean_medianTDC[ifC][depth] = new TH2F(Form("h2_mean_medianTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Mean - median TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_sigma[ifC][depth] = new TH2F(Form("h2_sigma_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi; Pull from 81 ns [sigma]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_meanErrTDC[ifC][depth] = new TH2F(Form("h2_meanErrTDC_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Error on mean TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_meanTS2[ifC][depth] = new TH2F(Form("h2_meanTS2_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Mean TS2 fraction [%]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_meanTS4[ifC][depth] = new TH2F(Form("h2_meanTS4_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Mean TS4 fraction [%]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_meanE[ifC][depth] = new TH2F(Form("h2_meanE_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi;Mean total charge [fC]", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);
                h2_correction[ifC][depth] = new TH2F(Form("h2_correction_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;iPhi; Correction to 81 [ns];", depth + 1, ifC) , 59, -29.5, 29.5, 73, 0.5, 73.5);


                h1_entries_vs_eta[ifC][depth] = new TH1F(Form("h1_entries_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Entries", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_meanTDC_vs_eta[ifC][depth] = new TH1F(Form("h1_meanTDC_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Mean TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_rmsTDC_vs_eta[ifC][depth] = new TH1F(Form("h1_rmsTDC_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;RMS in TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_medianTDC_vs_eta[ifC][depth] = new TH1F(Form("h1_medianTDC_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Median TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_mean_medianTDC_vs_eta[ifC][depth] = new TH1F(Form("h1_mean_medianTDC_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Mean - median TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_sigma_vs_eta[ifC][depth] = new TH1F(Form("h1_sigma_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta; Pull from 81 ns [sigma]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_meanErrTDC_vs_eta[ifC][depth] = new TH1F(Form("h1_meanErrTDC_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Error on mean TDC time [ns]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_meanTS2_vs_eta[ifC][depth] = new TH1F(Form("h1_meanTS2_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Mean TS2 fraction [%]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_meanTS4_vs_eta[ifC][depth] = new TH1F(Form("h1_meanTS4_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Mean TS4 fraction [%]", depth + 1, ifC) , 59, -29.5, 29.5);
                h1_meanE_vs_eta[ifC][depth] = new TH1F(Form("h1_meanE_vs_eta_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iEta;Mean total charge [fC]", depth + 1, ifC) , 59, -29.5, 29.5);


                // Use excessive x-axis range to make space for legend
                h1_entries_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_entries_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Entries", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanTDC_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_meanTDC_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_rmsTDC_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_rmsTDC_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;RMS TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_medianTDC_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_medianTDC_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Median TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_mean_medianTDC_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_mean_medianTDC_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean - median TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_sigma_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_sigma_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi; Pull from 81 ns [sigma]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanErrTDC_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_meanErrTDC_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Error on mean TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanTS2_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_meanTS2_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean TS2 fraction [%]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanTS4_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_meanTS4_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean TS4 fraction [%]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanE_vs_phi_HEM[ifC][depth] = new TH1F(Form("h1_meanE_vs_phi_HEM_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean total charge [fC]", depth + 1, ifC) , 96, 0.5, 96.5);

                h1_entries_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_entries_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Entries", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanTDC_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_meanTDC_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_rmsTDC_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_rmsTDC_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;RMS in TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_medianTDC_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_medianTDC_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Median TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_mean_medianTDC_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_mean_medianTDC_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean - median TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_sigma_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_sigma_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi; Pull from 81 ns [sigma]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanErrTDC_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_meanErrTDC_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Error on mean TDC time [ns]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanTS2_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_meanTS2_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean TS2 fraction [%]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanTS4_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_meanTS4_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean TS4 fraction [%]", depth + 1, ifC) , 96, 0.5, 96.5);
                h1_meanE_vs_phi_HEP[ifC][depth] = new TH1F(Form("h1_meanE_vs_phi_HEP_fC%i_depth%i", ifC, depth + 1), Form("Depth %i, charge %i;iPhi;Mean total charge [fC]", depth + 1, ifC) , 96, 0.5, 96.5);
            }
        }
    }
    //Start looping over all channels to fill all histograms and lists
    if (mode == 0 || mode == 1) {
        for (int ifC = 0; ifC < nfC; ifC++) {
            if (ifC != 0 && ifC != 2) continue; // only do inclusive charge (for charge plots) and measurement bin (for time plots) to save time.
            for (int iWedge = 1; iWedge <= nWedges; iWedge++) {
                if (wedgeSel != 0 && iWedge != wedgeSel) continue;

                //Get phi values corresponding to this wedge
                gSystem->mkdir(Form("%s/HE%i", plotDir.Data(), iWedge), kTRUE);
                gSystem->mkdir(Form("%s/HE%i", webplotDir.Data(), iWedge), kTRUE);
                cout << "Wedge " << iWedge << endl;
                vector<int> iphis_this_wedge;
                if (iWedge == 1) {
                    iphis_this_wedge.push_back(71); iphis_this_wedge.push_back(72); iphis_this_wedge.push_back(1); iphis_this_wedge.push_back(2);
                    //   iphis_this_wedge.push_back(73); // iphi 73: fake iphi placeholder for integral over iphi
                }
                else {
                    for (int i = 0; i < 4; i++) {iphis_this_wedge.push_back(iWedge * 4 - 5 + i);}
                }
                //Loop over all 
                for (int ieta = 0; ieta < nieta; ieta++) {
                    int labelEta = ieta + minieta;
                    if (abs(labelEta) <= 15)continue; //These are not part of HE

                    //Prepare pointers to hold histograms for this iEta
                    vector<TH1F*> pulses, peds;
                    vector<TH2F*> chgTS2vsTDC, chgTS4vsTDC, chgvsTDC;
                    vector< vector<TH1F*> > energy, energyADC, timeTDC, chgTS2, chgTS4;
                    cout << "Eta " << labelEta << endl;

                    TString detector = "HEP";
                    if (labelEta < 0) detector = "HEM";

                    //Form name for this channel
                    TString group_name = detector + Form("%i_ieta%i_fC%i", iWedge, labelEta, ifC);
                    TString cos_name;
                    if (ifC == 0) cos_name = detector + Form("%i, iEta%i, Q > %i fC", iWedge, labelEta, fCbins[1]);
                    else if (ifC > 0 && ifC + 1 < nfC) cos_name = detector + Form("%i, iEta%i, %i < Q < %i fC", iWedge, labelEta, fCbins[ifC], fCbins[ifC + 1]);
                    else cos_name = Form("HE%i, iEta%i, Q > %i fC", iWedge, labelEta, fCbins[ifC]);
                    // Loop over depths
                    for (int depth = 0; depth < ndepth; depth++) {
                        int labelDepth = depth + mindepth;
                        if (labelDepth == 7 && abs(labelEta) < 26) continue; // skip channels that don't exist to save time.
                        if (labelDepth == 6 && abs(labelEta) <= 18) continue;
                        vector<TH1F*> energy_thisdepth, energyADC_thisdepth, timeTDC_thisdepth, chgTS2_thisdepth, chgTS4_thisdepth;
                        TH1F* pulses_combined, *peds_combined;
                        TH2F* chgTS2vsTDC_combined, *chgTS4vsTDC_combined, *chgvsTDC_combined;

                        for (int i = 0; i < iphis_this_wedge.size(); i++) {
                            int labelPhi = iphis_this_wedge[i];

                            if (abs(labelEta) >= 21 && labelPhi % 2 == 0) continue;
                            /// 1D plots: overlay 4 different phi histograms for each each eta,depth canvas
                            TH1F * energy_thisphi = (TH1F*)inFile->Get(Form("h1_energy_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                            TH1F * energyADC_thisphi = (TH1F*)inFile->Get(Form("h1_energyADC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                            TH1F * timeTDC_thisphi = (TH1F*)inFile->Get(Form("h1_TDC_time_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                            TH1F * chgTS2_thisphi = (TH1F*)inFile->Get(Form("h1_chgfracTS2_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                            TH1F * chgTS4_thisphi = (TH1F*)inFile->Get(Form("h1_chgfracTS4_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));


                            if (energy_thisphi) { //if hist is found, 
                                //store histograms in vector for plotting with appropriate group
                                energy_thisdepth.push_back( energy_thisphi );
                                energyADC_thisdepth.push_back( energyADC_thisphi );
                                float thisE = energy_thisphi->GetMean();
                                float thisErrE = energy_thisphi->GetMeanError();

                                //Add values for this channel to lists
                                meanEs[ifC].push_back(thisE);
                                meanErrEs[ifC].push_back(thisErrE);

                                //Set values for 2D histograms
                                h2_meanE[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisE);
                                //Fill 1D histogram
                                h1_meanE[ifC][depth]->Fill(thisE);

                            }

                            if (timeTDC_thisphi) {//if hist is found, 
                                //store histogram in vector for plotting with appropriate group
                                timeTDC_thisdepth.push_back( timeTDC_thisphi );
                                float thisEntries = timeTDC_thisphi->GetEntries();
                                float thisT = timeTDC_thisphi->GetMean();
                                float thisTerr = timeTDC_thisphi->GetMeanError();
                                float thisRMS = timeTDC_thisphi->GetRMS();
                                float thisMedian = getMedian(timeTDC_thisphi);

                                //Add values for this channel to lists
                                entries[ifC].push_back(thisEntries);
                                etas[ifC].push_back(labelEta);
                                phis[ifC].push_back(labelPhi);
                                depths[ifC].push_back(labelDepth);

                                meanTDCs[ifC].push_back(thisT);
                                rmsTDCs[ifC].push_back(thisRMS);
                                meanErrTDCs[ifC].push_back(thisTerr);
                                medianTDCs[ifC].push_back(thisMedian);
                                medianErrTDCs[ifC].push_back(timeTDC_thisphi->GetBinWidth(1) / sqrt(12)); // use any bin width resolution as uncertainty on median

                                //Set values for 2D histograms
                                h2_entries[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisEntries);
                                h2_meanTDC[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisT);
                                h2_rmsTDC[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisRMS);
                                h2_medianTDC[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisMedian);
                                h2_mean_medianTDC[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisT - thisMedian);
                                h2_meanErrTDC[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisTerr);
                                h2_sigma[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, (81. - thisT) / thisTerr);
                                h2_correction[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, 81. - thisT);

                                //Fill 1D histograms
                                h1_entries[ifC][depth]->Fill( log10(thisEntries));
                                h1_meanTDC[ifC][depth]->Fill( thisT);
                                h1_rmsTDC[ifC][depth]->Fill( thisRMS);
                                h1_medianTDC[ifC][depth]->Fill( thisMedian);
                                h1_mean_medianTDC[ifC][depth]->Fill( thisT - thisMedian);
                                h1_meanErrTDC[ifC][depth]->Fill( thisTerr);
                                h1_sigma[ifC][depth]->Fill( (81. - thisT) / thisTerr);
                                h1_correction[ifC][depth]->Fill(81. - thisT);


                            }
                            //Same pattern again
                            if (chgTS2_thisphi) {
                                chgTS2_thisdepth.push_back( chgTS2_thisphi );
                                float thisMeanTS2 = chgTS2_thisphi->GetMean();
                                float thisMeanErrTS2 = chgTS2_thisphi->GetMeanError();
                                meanTS2s[ifC].push_back(thisMeanTS2);
                                meanErrTS2s[ifC].push_back(thisMeanErrTS2);

                                h2_meanTS2[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisMeanTS2);
                                h1_meanTS2[ifC][depth]->Fill(thisMeanTS2);

                            }
                            //Same pattern again
                            if (chgTS4_thisphi) {
                                chgTS4_thisdepth.push_back( chgTS4_thisphi );
                                float thisMeanTS4 = chgTS4_thisphi->GetMean();
                                float thisMeanErrTS4 = chgTS4_thisphi->GetMeanError();
                                meanTS4s[ifC].push_back(thisMeanTS4);
                                meanErrTS4s[ifC].push_back(thisMeanErrTS4);

                                h2_meanTS4[ifC][depth]->SetBinContent(labelEta + 30, labelPhi, thisMeanTS4);
                                h1_meanTS4[ifC][depth]->Fill(thisMeanTS4);

                            }
                            if (mode == 0) {
                                /// 2D plots, and pulse shapes: 
                                // instead of overlaying 4 plots, average 4 histograms together in each eta,depth canvas (overlay not practical)
                             
                                if (i == 0) {//first phi channel: clone histogram
                                    pulses_combined = (TH1F*)inFile->Get(Form("h1_fC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    peds_combined = (TH1F*)inFile->Get(Form("h1_ped_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    chgTS2vsTDC_combined = (TH2F*)inFile->Get(Form("h2_chgfracTS2vsTDC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    chgTS4vsTDC_combined = (TH2F*)inFile->Get(Form("h2_chgfracTS4vsTDC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    chgvsTDC_combined = (TH2F*)inFile->Get(Form("h2_TotFCvsTDC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));

                                }
                                else {//Later phi channels: add to clone

                                    //Get
                                    TH1F * pulse_thisphi = (TH1F*)inFile->Get(Form("h1_fC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    TH1F * peds_thisphi = (TH1F*)inFile->Get(Form("h1_ped_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    TH2F * chgTS2vsTDC_thisphi = (TH2F*)inFile->Get(Form("h2_chgfracTS2vsTDC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    TH2F * chgTS4vsTDC_thisphi = (TH2F*)inFile->Get(Form("h2_chgfracTS4vsTDC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));
                                    TH2F * chgvsTDC_thisphi = (TH2F*)inFile->Get(Form("h2_TotFCvsTDC_ieta%i_iphi%i_idepth%i_fC%i", labelEta, labelPhi, labelDepth, ifC));

                                    //Add
                                    if (pulses_combined && pulse_thisphi) pulses_combined->Add(pulse_thisphi);
                                    if (peds_combined && peds_thisphi) peds_combined->Add(peds_thisphi);
                                    if (chgTS2vsTDC_combined && chgTS2vsTDC_thisphi) chgTS2vsTDC_combined->Add(chgTS2vsTDC_thisphi);
                                    if (chgTS4vsTDC_combined && chgTS4vsTDC_thisphi) chgTS4vsTDC_combined->Add(chgTS4vsTDC_thisphi);
                                    if (chgvsTDC_combined && chgvsTDC_thisphi) chgvsTDC_combined->Add(chgvsTDC_thisphi);
                                }
                            }
                        }//loop over phi in this wedge

                        // Add summed histograms to list of depths for this wedge,eta
                        if (mode == 0) {
                            if (pulses_combined) {
                                pulses_combined->SetTitle(cos_name + Form(", Depth %i, all Phi", labelDepth));
                                pulses.push_back(pulses_combined);
                            }
                            if (peds_combined) {
                                peds_combined->SetTitle(cos_name + Form(", Depth %i, all Phi", labelDepth));
                                peds.push_back(peds_combined);
                            }

                            if (chgTS2vsTDC_combined) {
                                chgTS2vsTDC_combined->SetTitle(cos_name + Form(", Depth %i, all Phi", labelDepth));
                                chgTS2vsTDC.push_back(chgTS2vsTDC_combined);
                            }
                            if (chgTS4vsTDC_combined) {
                                chgTS4vsTDC_combined->SetTitle(cos_name + Form(", Depth %i, all Phi", labelDepth));
                                chgTS4vsTDC.push_back(chgTS4vsTDC_combined);
                            }
                            if (chgvsTDC_combined) {
                                chgvsTDC_combined->SetTitle(cos_name + Form(", Depth %i, all Phi", labelDepth));
                                chgvsTDC.push_back(chgvsTDC_combined);
                            }
                            if (energy_thisdepth.size() > 0) energy.push_back(energy_thisdepth);
                            if (energyADC_thisdepth.size() > 0) energyADC.push_back(energyADC_thisdepth);
                            if (timeTDC_thisdepth.size() > 0) timeTDC.push_back(timeTDC_thisdepth);
                            if (chgTS2_thisdepth.size() > 0) chgTS2.push_back(chgTS2_thisdepth);
                            if (chgTS4_thisdepth.size() > 0) chgTS4.push_back(chgTS4_thisdepth);
                        }
                    }//loop over depth

                    /// Plot the per-wedge single channel distributions
                    if (mode == 0) {
                        if (ifC == 2) {//Only plot measurement charge bin (save time)
                            if (pulses.size() > 0) {
                                plot_pulse(pulses, peds, "pulseshape_" + group_name, iWedge);
                                plot_TH2F(chgTS2vsTDC, "chgTS2vsTDC_" + group_name, iWedge);
                                plot_TH2F(chgTS4vsTDC, "chgTS4vsTDC_" + group_name, iWedge);
                                plot_vecTH1F(timeTDC, "timeTDC_" + group_name, iWedge);
                                plot_vecTH1F(chgTS2, "chgTS2_" + group_name, iWedge);
                                plot_vecTH1F(chgTS4, "chgTS4_" + group_name, iWedge);
                            }
                        }
                        //Plots with total charge as an axis: only use inclusive version
                        if (ifC == 0) {
                            if (energy.size() > 0) {
                                plot_vecTH1F(energy, "energy_" + group_name, iWedge);
                                plot_vecTH1F(energyADC, "energyADC_" + group_name, iWedge);
                                plot_TH2F(chgvsTDC, "chgvsTDC_" + group_name, iWedge);
                            }
                        }
                    }
                }//loop over eta
            }//loop over wedge
        }//loop over fC
    } //mode !=2


    //Loop again to get 1D averaged summary plots as function of eta, phi
    if (mode == 0 || mode == 2) {
        for (int ifC = 0; ifC < nfC; ifC++) {
            if (ifC != 2) continue;
            for (int depth = 0; depth < ndepth; depth++) {
                int labelDepth = depth + mindepth;

                //Loop over eta to fill histograms as function of eta
                for (int ieta = 0; ieta < nieta; ieta++) {
                    int labelEta = ieta + minieta;
                    if (abs(labelEta) <= 15) continue;
                    /// iPhi 73 corresponds to histograms that combine all channels in phi
                    TH1F * energy_allphi = (TH1F*)inFile->Get(Form("h1_energy_ieta%i_iphi%i_idepth%i_fC%i", labelEta, 73, labelDepth, ifC));
                    TH1F * timeTDC_allphi = (TH1F*)inFile->Get(Form("h1_TDC_time_ieta%i_iphi%i_idepth%i_fC%i", labelEta, 73, labelDepth, ifC));
                    TH1F * chgTS2_allphi = (TH1F*)inFile->Get(Form("h1_chgfracTS2_ieta%i_iphi%i_idepth%i_fC%i", labelEta, 73, labelDepth, ifC));
                    TH1F * chgTS4_allphi = (TH1F*)inFile->Get(Form("h1_chgfracTS4_ieta%i_iphi%i_idepth%i_fC%i", labelEta, 73, labelDepth, ifC));

                    if (energy_allphi) {
                        float thisE = energy_allphi->GetMean();
                        float thisErrE = energy_allphi->GetMeanError();
                        float thisEntries = timeTDC_allphi->GetEntries();
                        float thisT = timeTDC_allphi->GetMean();
                        float thisRMS = timeTDC_allphi->GetRMS();
                        float thisTerr = timeTDC_allphi->GetMeanError();
                        float thisMedian = getMedian(timeTDC_allphi);
                        float thisMeanTS2 = chgTS2_allphi->GetMean();
                        float thisMeanErrTS2 = chgTS2_allphi->GetMeanError();
                        float thisMeanTS4 = chgTS4_allphi->GetMean();
                        float thisMeanErrTS4 = chgTS4_allphi->GetMeanError();

                        h1_meanE_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisE);
                        h1_meanE_vs_eta[ifC][depth]->SetBinError(labelEta + 30, thisErrE);

                        h1_entries_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisEntries);
                        h1_entries_vs_eta[ifC][depth]->SetBinError(labelEta + 30, pow(thisEntries, 0.5));
                        h1_meanTDC_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisT);
                        h1_meanTDC_vs_eta[ifC][depth]->SetBinError(labelEta + 30, thisTerr);
                        h1_rmsTDC_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisRMS);
                        h1_medianTDC_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisMedian);
                        h1_medianTDC_vs_eta[ifC][depth]->SetBinError(labelEta + 30, timeTDC_allphi->GetBinWidth(1) / sqrt(12));
                        h1_mean_medianTDC_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisT - thisMedian);
                        h1_mean_medianTDC_vs_eta[ifC][depth]->SetBinError(labelEta + 30, thisTerr);
                        h1_meanErrTDC_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisTerr);
                        h1_sigma_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, (81. - thisT) / thisTerr);

                        h1_meanTS2_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisMeanTS2);
                        h1_meanTS2_vs_eta[ifC][depth]->SetBinError(labelEta + 30, thisMeanErrTS2);

                        h1_meanTS4_vs_eta[ifC][depth]->SetBinContent(labelEta + 30, thisMeanTS4);
                        h1_meanTS4_vs_eta[ifC][depth]->SetBinError(labelEta + 30, thisMeanErrTS4);

                    }
                }
                //loop over phi to fill hists as function of phi
                for (int iphi = 0; iphi < niphi; iphi++) {
                    int labelPhi = iphi + miniphi;
                    //Eta = 30 corresponds to histograms combined in eta for all eta in HEM 
                    TH1F * energy_alletaHEM = (TH1F*)inFile->Get(Form("h1_energy_ieta%i_iphi%i_idepth%i_fC%i", 30, labelPhi, labelDepth, ifC));
                    TH1F * timeTDC_alletaHEM = (TH1F*)inFile->Get(Form("h1_TDC_time_ieta%i_iphi%i_idepth%i_fC%i", 30, labelPhi, labelDepth, ifC));
                    TH1F * chgTS2_alletaHEM = (TH1F*)inFile->Get(Form("h1_chgfracTS2_ieta%i_iphi%i_idepth%i_fC%i", 30, labelPhi, labelDepth, ifC));
                    TH1F * chgTS4_alletaHEM = (TH1F*)inFile->Get(Form("h1_chgfracTS4_ieta%i_iphi%i_idepth%i_fC%i", 30, labelPhi, labelDepth, ifC));

                    //Eta = 30 corresponds to histograms combined in eta for all eta in HEP
                    TH1F * energy_alletaHEP = (TH1F*)inFile->Get(Form("h1_energy_ieta%i_iphi%i_idepth%i_fC%i", 31, labelPhi, labelDepth, ifC));
                    TH1F * timeTDC_alletaHEP = (TH1F*)inFile->Get(Form("h1_TDC_time_ieta%i_iphi%i_idepth%i_fC%i", 31, labelPhi, labelDepth, ifC));
                    TH1F * chgTS2_alletaHEP = (TH1F*)inFile->Get(Form("h1_chgfracTS2_ieta%i_iphi%i_idepth%i_fC%i", 31, labelPhi, labelDepth, ifC));
                    TH1F * chgTS4_alletaHEP = (TH1F*)inFile->Get(Form("h1_chgfracTS4_ieta%i_iphi%i_idepth%i_fC%i", 31, labelPhi, labelDepth, ifC));

                    if (energy_alletaHEM) {
                        //Fill HEM values
                        float thisE = energy_alletaHEM->GetMean();
                        float thisErrE = energy_alletaHEM->GetMeanError();
                        float thisEntries = timeTDC_alletaHEM->GetEntries();
                        float thisT = timeTDC_alletaHEM->GetMean();
                        float thisRMS = timeTDC_alletaHEM->GetRMS();
                        float thisTerr = timeTDC_alletaHEM->GetMeanError();
                        float thisMedian = getMedian(timeTDC_alletaHEM);
                        float thisMeanTS2 = chgTS2_alletaHEM->GetMean();
                        float thisMeanErrTS2 = chgTS2_alletaHEM->GetMeanError();
                        float thisMeanTS4 = chgTS4_alletaHEM->GetMean();
                        float thisMeanErrTS4 = chgTS4_alletaHEM->GetMeanError();

                        h1_meanE_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisE);
                        h1_meanE_vs_phi_HEM[ifC][depth]->SetBinError(labelPhi, thisErrE);

                        h1_entries_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisEntries);
                        h1_entries_vs_phi_HEM[ifC][depth]->SetBinError(labelPhi, pow(thisEntries, 0.5));
                        h1_meanTDC_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisT);
                        h1_meanTDC_vs_phi_HEM[ifC][depth]->SetBinError(labelPhi, thisTerr);
                        h1_rmsTDC_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisRMS);
                        h1_medianTDC_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisMedian);
                        h1_medianTDC_vs_phi_HEM[ifC][depth]->SetBinError(labelPhi, timeTDC_alletaHEM->GetBinWidth(1) / sqrt(12));
                        h1_mean_medianTDC_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisT - thisMedian);
                        h1_mean_medianTDC_vs_phi_HEM[ifC][depth]->SetBinError(labelPhi, thisTerr);
                        h1_meanErrTDC_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisTerr);
                        h1_sigma_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, (81. - thisT) / thisTerr);

                        h1_meanTS2_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisMeanTS2);
                        h1_meanTS2_vs_phi_HEM[ifC][depth]->SetBinError(labelPhi, thisMeanErrTS2);
                        h1_meanTS4_vs_phi_HEM[ifC][depth]->SetBinContent(labelPhi, thisMeanTS4);
                        h1_meanTS4_vs_phi_HEM[ifC][depth]->SetBinError(labelPhi, thisMeanErrTS4);

                    }
                    if (energy_alletaHEP) {
                        //Fill HEP values
                        float thisE = energy_alletaHEP->GetMean();
                        float thisErrE = energy_alletaHEP->GetMeanError();
                        float thisEntries = timeTDC_alletaHEP->GetEntries();
                        float thisT = timeTDC_alletaHEP->GetMean();
                        float thisRMS = timeTDC_alletaHEP->GetRMS();
                        float thisTerr = timeTDC_alletaHEP->GetMeanError();
                        float thisMedian = getMedian(timeTDC_alletaHEP);
                        float thisMeanTS2 = chgTS2_alletaHEP->GetMean();
                        float thisMeanErrTS2 = chgTS2_alletaHEP->GetMeanError();
                        float thisMeanTS4 = chgTS4_alletaHEP->GetMean();
                        float thisMeanErrTS4 = chgTS4_alletaHEP->GetMeanError();
                        h1_meanE_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisE);
                        h1_meanE_vs_phi_HEP[ifC][depth]->SetBinError(labelPhi, thisErrE);

                        h1_entries_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisEntries);
                        h1_entries_vs_phi_HEP[ifC][depth]->SetBinError(labelPhi, pow(thisEntries, 0.5));
                        h1_meanTDC_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisT);
                        h1_meanTDC_vs_phi_HEP[ifC][depth]->SetBinError(labelPhi, thisTerr);
                        h1_rmsTDC_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisRMS);
                        h1_medianTDC_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisMedian);
                        h1_medianTDC_vs_phi_HEP[ifC][depth]->SetBinError(labelPhi, timeTDC_alletaHEP->GetBinWidth(1) / sqrt(12));
                        h1_mean_medianTDC_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisT - thisMedian);
                        h1_mean_medianTDC_vs_phi_HEP[ifC][depth]->SetBinError(labelPhi, thisTerr);
                        h1_meanErrTDC_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisTerr);
                        h1_sigma_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, (81. - thisT) / thisTerr);

                        h1_meanTS2_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisMeanTS2);
                        h1_meanTS2_vs_phi_HEP[ifC][depth]->SetBinError(labelPhi, thisMeanErrTS2);
                        h1_meanTS4_vs_phi_HEP[ifC][depth]->SetBinContent(labelPhi, thisMeanTS4);
                        h1_meanTS4_vs_phi_HEP[ifC][depth]->SetBinError(labelPhi, thisMeanErrTS4);

                    }
                }
            }
        }
    }
    if (mode == 0 || mode == 1) writeCorrectionTable(Form("test_%s_fC%i.csv", tag.Data(), 2), 2);
    setSummaryStyle();

    for (int ifC = 0; ifC < nfC; ifC++)
    {
        if (ifC != 2) continue;
        if (mode == 0 || mode == 1) {
            //These make scatter plots based on the vector<floats> filled once per channel
            printSummaryScatter(TGraphErrors(meanTDCs[ifC].size(), meanTDCs[ifC].data(), meanTS2s[ifC].data(), meanErrTDCs[ifC].data(), meanErrTS2s[ifC].data()), Form("TS2frac_vs_meanT_fC%i", ifC), "Mean TDC time [ns]", "Mean TS2 charge fraction [%]");
            printSummaryScatter(TGraphErrors(meanTDCs[ifC].size(), meanTDCs[ifC].data(), meanTS4s[ifC].data(), meanErrTDCs[ifC].data(), meanErrTS4s[ifC].data()), Form("TS4frac_vs_meanT_fC%i", ifC), "Mean TDC time [ns]", "Mean TS4 charge fraction [%]");
            printSummaryScatter(TGraphErrors(meanTDCs[ifC].size(), meanTDCs[ifC].data(), medianTDCs[ifC].data(), meanErrTDCs[ifC].data(), medianErrTDCs[ifC].data()), Form("medianT_vs_meanT_fC%i", ifC), "Mean TDC time [ns]", "Median TDC time [ns]");
        }
        vector<TH1F*> vh1_meanTDC_vs_eta, vh1_meanTS2_vs_eta, vh1_meanTS4_vs_eta, vh1_rmsTDC_vs_eta, vh1_meanE_vs_eta ;
        vector<TH1F*> vh1_meanTDC_vs_phi_HEM, vh1_meanTS2_vs_phi_HEM, vh1_meanTS4_vs_phi_HEM, vh1_rmsTDC_vs_phi_HEM, vh1_meanE_vs_phi_HEM;
        vector<TH1F*> vh1_meanTDC_vs_phi_HEP, vh1_meanTS2_vs_phi_HEP, vh1_meanTS4_vs_phi_HEP, vh1_rmsTDC_vs_phi_HEP, vh1_meanE_vs_phi_HEP;
        for (int idepth = 0; idepth < ndepth; idepth++) {
            //if (idepth != 0) continue;
            if (mode == 0 || mode == 1) {
                //Print all the eta-phi maps
                printSummaryTH2(h2_meanTDC[ifC][idepth], Form("meanTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_rmsTDC[ifC][idepth], Form("rmsTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_meanErrTDC[ifC][idepth], Form("meanErrTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_medianTDC[ifC][idepth], Form("medianTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_mean_medianTDC[ifC][idepth], Form("mean_medianTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_entries[ifC][idepth], Form("entries_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_meanTS2[ifC][idepth], Form("meanTS2_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_meanTS4[ifC][idepth], Form("meanTS4_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_correction[ifC][idepth], Form("correction_idepth%i_fC%i", idepth + 1, ifC));
                printSummaryTH2(h2_meanE[0][idepth], Form("meanE_idepth%i_fC%i", idepth + 1, 0));

                //print 1D summary of eta-phi maps
                printSimpleTH1F(h1_meanTDC[ifC][idepth], Form("meanTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_rmsTDC[ifC][idepth], Form("rmsTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_meanErrTDC[ifC][idepth], Form("meanErrTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_medianTDC[ifC][idepth], Form("medianTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_mean_medianTDC[ifC][idepth], Form("mean_medianTDC_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_entries[ifC][idepth], Form("entries_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_meanTS2[ifC][idepth], Form("meanTS2_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_meanTS4[ifC][idepth], Form("meanTS4_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_correction[ifC][idepth], Form("correction_idepth%i_fC%i", idepth + 1, ifC));
                printSimpleTH1F(h1_meanE[0][idepth], Form("meanE_idepth%i_fC%i", idepth + 1, 0));

                h1_meanTDC[ifC][ndepth]->Add( h1_meanTDC[ifC][idepth]);
                h1_rmsTDC[ifC][ndepth]->Add( h1_rmsTDC[ifC][idepth]);
                h1_meanErrTDC[ifC][ndepth]->Add(h1_meanErrTDC[ifC][idepth]);
                h1_medianTDC[ifC][ndepth]->Add(h1_medianTDC[ifC][idepth]);
                h1_mean_medianTDC[ifC][ndepth]->Add(h1_mean_medianTDC[ifC][idepth]);
                h1_entries[ifC][ndepth]->Add(h1_entries[ifC][idepth]);
                h1_meanTS2[ifC][ndepth]->Add(h1_meanTS2[ifC][idepth]);
                h1_meanTS4[ifC][ndepth]->Add(h1_meanTS4[ifC][idepth]);
                h1_correction[ifC][ndepth]->Add(h1_correction[ifC][idepth]);

            }

            //Prepare vector of 1D summary plots vs eta or phi to overlap different depths on one canvas
            vh1_meanTDC_vs_eta.push_back(h1_meanTDC_vs_eta[ifC][idepth]);
            vh1_rmsTDC_vs_eta.push_back(h1_rmsTDC_vs_eta[ifC][idepth]);
            vh1_meanTS2_vs_eta.push_back(h1_meanTS2_vs_eta[ifC][idepth]);
            vh1_meanTS4_vs_eta.push_back(h1_meanTS4_vs_eta[ifC][idepth]);
            vh1_meanE_vs_eta.push_back(h1_meanE_vs_eta[0][idepth]);

            vh1_meanTDC_vs_phi_HEM.push_back(h1_meanTDC_vs_phi_HEM[ifC][idepth]);
            vh1_rmsTDC_vs_phi_HEM.push_back(h1_rmsTDC_vs_phi_HEM[ifC][idepth]);
            vh1_meanTS2_vs_phi_HEM.push_back(h1_meanTS2_vs_phi_HEM[ifC][idepth]);
            vh1_meanTS4_vs_phi_HEM.push_back(h1_meanTS4_vs_phi_HEM[ifC][idepth]);
            vh1_meanE_vs_phi_HEM.push_back(h1_meanE_vs_phi_HEM[0][idepth]);


            vh1_meanTDC_vs_phi_HEP.push_back(h1_meanTDC_vs_phi_HEP[ifC][idepth]);
            vh1_rmsTDC_vs_phi_HEP.push_back(h1_rmsTDC_vs_phi_HEP[ifC][idepth]);
            vh1_meanTS2_vs_phi_HEP.push_back(h1_meanTS2_vs_phi_HEP[ifC][idepth]);
            vh1_meanTS4_vs_phi_HEP.push_back(h1_meanTS4_vs_phi_HEP[ifC][idepth]);
            vh1_meanE_vs_phi_HEP.push_back(h1_meanE_vs_phi_HEP[0][idepth]);

        }
        if (mode == 0 || mode == 2) {

            //summary plots vs eta, for depth 1-3
            printSummaryTH1vec(vh1_meanTDC_vs_eta, Form("meanTDC_vs_eta_fC%i", ifC), 78, 84);
            vector<TH1F*> vh1_meanTDC_vs_eta_depth123;
            vh1_meanTDC_vs_eta_depth123.push_back(vh1_meanTDC_vs_eta[0]);
            vh1_meanTDC_vs_eta_depth123.push_back(vh1_meanTDC_vs_eta[1]);
            vh1_meanTDC_vs_eta_depth123.push_back(vh1_meanTDC_vs_eta[2]);
            printSummaryTH1vec(vh1_meanTDC_vs_eta_depth123, Form("meanTDC_vs_eta_depth123_fC%i", ifC), 78, 84);

            //summary plots vs eta, for depth 4-7
            vector<TH1F*> vh1_meanTDC_vs_eta_depth4567;
            vh1_meanTDC_vs_eta_depth4567.push_back(vh1_meanTDC_vs_eta[3]);
            vh1_meanTDC_vs_eta_depth4567.push_back(vh1_meanTDC_vs_eta[4]);
            vh1_meanTDC_vs_eta_depth4567.push_back(vh1_meanTDC_vs_eta[5]);
            vh1_meanTDC_vs_eta_depth4567.push_back(vh1_meanTDC_vs_eta[6]);
            printSummaryTH1vec(vh1_meanTDC_vs_eta_depth4567, Form("meanTDC_vs_eta_depth4567_fC%i", ifC), 78, 84);

            printSummaryTH1vec(vh1_rmsTDC_vs_eta, Form("rmsTDC_vs_eta_fC%i", ifC), 0, 3);
            printSummaryTH1vec(vh1_meanTS2_vs_eta, Form("meanTS2_vs_eta_fC%i", ifC), 0, 2);
            printSummaryTH1vec(vh1_meanTS2_vs_eta, Form("meanTS2_vs_eta_fC%i", ifC), 0.001, 3, true);
            printSummaryTH1vec(vh1_meanTS4_vs_eta, Form("meanTS4_vs_eta_fC%i", ifC), 0, 55);
            printSummaryTH1vec(vh1_meanE_vs_eta, Form("meanE_vs_eta_fC%i", 0), 5000, 15000);

            //Summary vs phi, for HEM and HEP separately
            printSummaryTH1vec(vh1_meanTDC_vs_phi_HEM, Form("meanTDC_vs_phi_HEM_fC%i", ifC), 75, 87);
            printSummaryTH1vec(vh1_rmsTDC_vs_phi_HEM, Form("rmsTDC_vs_phi_HEM_fC%i", ifC), 0, 3);
            printSummaryTH1vec(vh1_meanTS2_vs_phi_HEM, Form("meanTS2_vs_phi_HEM_fC%i", ifC), 0, 2);
            printSummaryTH1vec(vh1_meanTS2_vs_phi_HEM, Form("meanTS2_vs_phi_HEM_fC%i", ifC), 0.001, 3, true);
            printSummaryTH1vec(vh1_meanTS4_vs_phi_HEM, Form("meanTS4_vs_phi_HEM_fC%i", ifC), 0, 55);
            printSummaryTH1vec(vh1_meanE_vs_phi_HEM, Form("meanE_vs_phi_HEM_fC%i", 0), 5000, 15000);

            printSummaryTH1vec(vh1_meanTDC_vs_phi_HEP, Form("meanTDC_vs_phi_HEP_fC%i", ifC), 75, 87);
            printSummaryTH1vec(vh1_rmsTDC_vs_phi_HEP, Form("rmsTDC_vs_phi_HEP_fC%i", ifC), 0, 3);
            printSummaryTH1vec(vh1_meanTS2_vs_phi_HEP, Form("meanTS2_vs_phi_HEP_fC%i", ifC), 0, 2);
            printSummaryTH1vec(vh1_meanTS2_vs_phi_HEP, Form("meanTS2_vs_phi_HEP_fC%i", ifC), 0.001, 3, true);
            printSummaryTH1vec(vh1_meanTS4_vs_phi_HEP, Form("meanTS4_vs_phi_HEP_fC%i", ifC), 0, 55);
            printSummaryTH1vec(vh1_meanE_vs_phi_HEP, Form("meanE_vs_phi_HEP_fC%i", 0), 5000, 15000);
        }

        if (mode == 0 || mode == 1) {

            //These print simple 1D summary plots filled once per channel
            printSimpleTH1F(h1_meanTDC[ifC][ndepth], Form("meanTDC_fC%i", ifC));
            printSimpleTH1F(h1_rmsTDC[ifC][ndepth], Form("rmsTDC_fC%i", ifC));
            printSimpleTH1F(h1_meanErrTDC[ifC][ndepth], Form("meanErrTDC_fC%i", ifC));
            printSimpleTH1F(h1_medianTDC[ifC][ndepth], Form("medianTDC_fC%i", ifC));
            printSimpleTH1F(h1_mean_medianTDC[ifC][ndepth], Form("mean_medianTDC_fC%i", ifC));
            printSimpleTH1F(h1_entries[ifC][ndepth], Form("entries_fC%i", ifC));
            printSimpleTH1F(h1_meanTS2[ifC][ndepth], Form("meanTS2_fC%i", ifC));
            printSimpleTH1F(h1_meanTS4[ifC][ndepth], Form("meanTS4_fC%i", ifC));
            printSimpleTH1F(h1_correction[ifC][ndepth], Form("correction_fC%i", ifC));

        }
    }
}

//This writes the correction table, can be useful for looking up channel information as well
void writeCorrectionTable(TString name, int ifC) {
    // Open file
    gSystem->Exec(Form("rm -f %s", name.Data()));
    ofstream fout(name.Data(), ios_base::app | ios_base::out);
    if (fout.is_open())
    {
        fout << "nDigis,iEta,iPhi,Depth,Mean TDC time[ns], TDC time RMS[ns], Uncertainty,Adjustment[ns],Adjustment[phase units]\n";
        for (int ichan = 0; ichan < etas[ifC].size(); ichan++) {

            int adjustment = 0;
            //Only adjust channels with 5 hits, and with mean more than 1 sigma away from 81 ns
            if (entries[ifC][ichan] > 5 && abs(81. - meanTDCs[ifC][ichan]) / meanErrTDCs[ifC][ichan] > 1 ) adjustment = round(2.*(81. - meanTDCs[ifC][ichan])); //0.5 ns granularity, rounded to integer
            fout << entries[ifC][ichan] << "," << etas[ifC][ichan] << "," << phis[ifC][ichan] << "," << depths[ifC][ichan] << "," << std::fixed << std::setprecision(2) << meanTDCs[ifC][ichan] << "," << rmsTDCs[ifC][ichan] << "," << std::setprecision(3) << meanErrTDCs[ifC][ichan] << "," << std::setprecision(1) << 81. - meanTDCs[ifC][ichan] << "," << adjustment << "\n";

        }
        fout.close();
    }

}

void printSummaryScatter(TGraphErrors g, TString name, TString xlabel, TString ylabel) {
    TCanvas c;
    c.SetGrid();
    //g.SetTitle("Run "+run)
    g.GetXaxis()->SetTitle(xlabel);
    g.GetYaxis()->SetTitle(ylabel);
    //g.GetYaxis()->SetTitleOffset(1.1);
    g.GetXaxis()->SetTitleOffset(1.1);
    g.SetLineWidth(0.5);
    //g.SetTitleSize(0.07);
    g.SetMarkerStyle(20);
    g.SetMarkerSize(0.4);
    //g.SetMarkerSize(0.4);
    g.SetMarkerColor(602);
    g.Draw("AEPZ");
    cout << g.GetN() << " channels" << endl;
    c.Print(summaryDir + "g_" + name + ".pdf");
    c.Print(webplotDir + "summary/g_" + name + ".png");
    c.SetLogy();
    c.Print(summaryDir + "g_" + name + "_log.pdf");
    c.Print(webplotDir + "summary/g_" + name + "_log.png");

    c.SetLogy(false);
    //Also print with no error bars-- too dense.
    g.Draw("APZ");

    c.Print(summaryDir + "g_no_err" + name + ".pdf");
    c.Print(webplotDir + "summary/g_no_err" + name + ".png");
    c.SetLogy();
    c.Print(summaryDir + "g_no_err_" + name + "_log.pdf");
    c.Print(webplotDir + "summary/g_no_err_" + name + "_log.png");
}

void printSimpleTH1F(TH1F * h, TString name) {
    TCanvas c;
    c.SetGrid();

    h1cosmetic(h);

    h->Draw("hist");

    c.Print(summaryDir + "h1_" + name + ".pdf");
    c.Print(webplotDir + "summary/h1_" + name + ".png");
    c.SetLogy();
    c.Print(summaryDir + "h1_" + name + "_log.pdf");
    c.Print(webplotDir + "summary/h1_" + name + "_log.png");
}




void printSummaryTH2(TH2F * h2, TString name, float min, float max) {
    //if h2->GetNbinsX() != 32: c = ROOT.TCanvas()
    //else:
    TCanvas c("c", "c", 1000, 600);
    c.SetRightMargin(0.15);
    c.SetLeftMargin(0.1);

    h2->SetLabelSize(0.05, "z");
    h2->SetLabelSize(0.05, "x");
    h2->SetLabelSize(0.05, "y");
    h2->SetStats(0);

    h2->SetTitleOffset(0.8, "y");
    h2->SetTitleOffset(0.82, "z");
    //   h2->SetTitleOffset(1.03,"z");
    //   h2->SetTitleOffset(0.9,"y");
    //   h2->SetTitleOffset(0.9,"x");
    //   h2->GetXaxis()->SetTitleSize(0.05);
    //   h2->GetYaxis()->SetTitleSize(0.05);
    h2->SetTitleSize(0.06, "X");
    h2->SetTitleSize(0.06, "Y");
    h2->GetZaxis()->SetTitleSize(0.06);
    //   h2->SetTitleSize(0.07);
    if (h2->GetMinimum() >= 0) h2->SetMinimum(0.95 * h2->GetMinimum(0.000001));
//    if (min != 0) h2->SetMinimum(0.95 * h2->GetMinimum(0.000001));
    if (max != 0) {
        h2->SetMaximum(max);
        h2->SetMinimum(min);
    }
    h2->Draw("colz");
    c.Print(summaryDir + "h2_" + name + ".pdf");
    c.Print(webplotDir + "summary/h2_" + name + ".png");

    c.SetLogz();
    c.Print(summaryDir + "h2_" + name + "_log.pdf");
    c.Print(webplotDir + "summary/h2_" + name + "_log.png");
}

void printSummaryTH1vec(vector<TH1F*> h1, TString name, float minn, float maxx, bool logY) {
    TCanvas c;
    c.SetGrid();
    c.SetLogy(logY);

    TLegend *leg;
    if (name.Contains("vs_eta")) leg = new TLegend(0.336, 0.14, 0.674, 0.43);
    //else if(name.Contains("meanTS2_vs_phi") && !logY) leg = new TLegend(0.12, 0.6, 0.9, 0.43);
    else leg = new TLegend(0.74, 0.14, 0.9, 0.65);

    TString drawstyle = "hist e";
    if (name.Contains("vs_phi")) drawstyle = "e";

    if (!name.Contains("vs_phi")) leg->SetNColumns(2);
    float min = 10000;
    float max = 0;

    float mult = 1.15;

    int nDepths = h1.size();

    //if "meanT_" in name: mult=1.05
    for (int idep = 0; idep < nDepths; idep++) {
        if (mult * h1[idep]->GetMaximum() > max) max =  mult * h1[idep]->GetMaximum();
        if (0.99 * h1[idep]->GetMinimum(0.000001) < min) min =  0.99 * h1[idep]->GetMinimum(0.000001);
    }
    // #  h1[maxdep].SetMinimum(min)
    //  # h1[maxdep].Draw("hist e")
    for (int idep = 0; idep < nDepths; idep++) {
        h1[idep]->SetStats(0);
        h1cosmetic(h1[idep]);
        h1[idep]->SetLineColor(colors[idep]);
        h1[idep]->SetLineWidth(2);

        if (maxx == 0) {
            h1[idep]->SetMinimum(min);
            if (max > 0) h1[idep]->SetMaximum(max);
        }
        else {
            h1[idep]->SetMinimum(minn);
            h1[idep]->SetMaximum(maxx);
        }

        if(!name.Contains("4567")) leg->AddEntry(h1[idep], Form("Depth %i", (1 + idep)), "el");
        else leg->AddEntry(h1[idep], Form("Depth %i", (4 + idep)), "el");
        
        if (idep == 0) h1[idep]->Draw(drawstyle);
        else h1[idep]->Draw(drawstyle + " same");
    }
    leg->Draw("same");
    if (! logY) {
        c.Print(summaryDir + "h1_" + name + ".pdf");
        c.Print(webplotDir + "summary/h1_" + name + ".png");
    }
    else {
        c.Print(summaryDir + "h1_" + name + "_log.pdf");
        c.Print(webplotDir + "summary/h1_" + name + "_log.png");
    }

}

void plot_TH2F(vector<TH2F*> h2_vec, TString name, int wedge)
{

    TCanvas *c = new TCanvas("c", "c", 1400, 600);
    if (h2_vec.size() < 7) c->Divide(3, 2);
    else c->Divide(4, 2);
    for (unsigned int ichan = 0; ichan < h2_vec.size(); ichan++)
    {
        //  cout<<"pad number "<<iphi<<endl;
        TPad *pad(NULL);
        pad = static_cast<TPad *>(c->cd(ichan + 1));
        pad->SetLeftMargin(0.2);
        pad->SetBottomMargin(0.14);
        //pad->SetTopMargin(0.16);
        pad->SetLogz();
        // pad->SetGrid();
        // TLatex *peak = new TLatex(0.7,0.84,Form("Peak TS: %i",));
        //peak->SetNDC();
        //peak->SetTextSize(textSize);


        h2_vec[ichan]->SetStats(0);
        h2_vec[ichan]->SetNdivisions(505, "X");
        //gStyle->SetTitleFontSize(0.1);
        h2_vec[ichan]->GetXaxis()->SetTitleSize(0.07);
        h2_vec[ichan]->GetXaxis()->SetLabelSize(0.06);
        h2_vec[ichan]->GetYaxis()->SetTitleSize(0.07);
        h2_vec[ichan]->GetYaxis()->SetLabelSize(0.06);
        h2_vec[ichan]->GetYaxis()->SetTitleOffset(1.15);
        //h2_vec[ichan]->SetMinimum(0);
        h2_vec[ichan]->Draw("colz");

    }

    c->Print(Form("%s/HE%i/%s.pdf", plotDir.Data(), wedge, name.Data()));
    c->Print(Form("%s/HE%i/%s.png", webplotDir.Data(), wedge, name.Data()));
    delete c;
}

void plot_TH1F(vector<TH1F*> h1_vec, TString name, int wedge)
{
    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    if (h1_vec.size() < 7) c->Divide(3, 2);
    else c->Divide(4, 2);
    for (unsigned int ichan = 0; ichan < h1_vec.size(); ichan++)
    {
        //  cout<<"pad number "<<iphi<<endl;
        TPad *pad(NULL);
        pad = static_cast<TPad *>(c->cd(ichan + 1));
        pad->SetLeftMargin(0.13);
        pad->SetBottomMargin(0.14);
        // pad->SetTopMargin(0.16);
        pad->SetGrid();
        pad->SetLogy();

        // TLatex *peak = new TLatex(0.7,0.84,Form("Peak TS: %i",));
        //peak->SetNDC();
        //peak->SetTextSize(textSize);
        h1_vec[ichan]->SetStats(0);
        //gStyle->SetTitleFontSize(0.1);

        h1_vec[ichan]->GetXaxis()->SetTitleSize(0.07);
        h1_vec[ichan]->GetXaxis()->SetLabelSize(0.06);
        h1_vec[ichan]->GetYaxis()->SetTitleSize(0.07);
        h1_vec[ichan]->GetYaxis()->SetLabelSize(0.06);
        h1_vec[ichan]->GetYaxis()->SetTitleOffset(1.15);
        //h1_vec[ichan]->SetMinimum(0);

        h1_vec[ichan]->SetNdivisions(505, "X");

        TString prec = "0";
        TString n = Form("Entries: %." + prec + "f", h1_vec[ichan]->GetEntries());
        TString m = Form("Mean: %." + prec + "f", h1_vec[ichan]->GetMean());
        TString r = Form("RMS: %." + prec + "f", h1_vec[ichan]->GetRMS());
        TString o = Form("Overflow: %.0f", h1_vec[ichan]->GetBinContent(h1_vec[ichan]->GetNbinsX() + 1));
        TString u = Form("Underflow: %.0f", h1_vec[ichan]->GetBinContent(0));


        if (h1_vec[ichan]->GetBinContent(0) > 0) { h1_vec[ichan]->SetBinContent(1, h1_vec[ichan]->GetBinContent(0) + h1_vec[ichan]->GetBinContent(1)); h1_vec[ichan]->SetBinContent(0, 0);}
        if (h1_vec[ichan]->GetBinContent(h1_vec[ichan]->GetNbinsX() + 1) > 0) {
            h1_vec[ichan]->SetBinContent(h1_vec[ichan]->GetNbinsX(), h1_vec[ichan]->GetBinContent(h1_vec[ichan]->GetNbinsX() + 1) + h1_vec[ichan]->GetBinContent(h1_vec[ichan]->GetNbinsX()));
            h1_vec[ichan]->SetBinContent(h1_vec[ichan]->GetNbinsX() + 1, 0);
        }

        if (h1_vec[ichan]->Integral() > 0) h1_vec[ichan]->Scale(1. / h1_vec[ichan]->Integral());
        //  cout <<"Mean before: "<<m<<" mean after: "<<h1_vec[ichan]->GetMean();
        h1_vec[ichan]->Draw("hist");
        float textSize = 0.055;
        float xpos = 0.57;



        TLatex *ent = new TLatex(xpos, 0.78, n);
        TLatex *mean = new TLatex(xpos, 0.72, m);
        TLatex *rms = new TLatex(xpos, 0.66, r);
        TLatex *of = new TLatex(xpos, 0.6, o);
        TLatex *uf = new TLatex(xpos, 0.54, u);
        ent->SetNDC();
        mean->SetNDC();
        of->SetNDC();
        rms->SetNDC();
        uf->SetNDC();
        ent->SetTextSize(textSize);
        mean->SetTextSize(textSize);
        of->SetTextSize(textSize);
        rms->SetTextSize(textSize);
        of->SetTextColor(kRed);
        uf->SetTextSize(textSize);
        uf->SetTextColor(kRed);
        mean->Draw();
        rms->Draw();
        ent->Draw();
        //  if(h1_vec[ichan]->GetBinContent(h1_vec[ichan]->GetNbinsX()+1) > 0) of->Draw();
        //if(h1_vec[ichan]->GetBinContent(0) > 0) uf->Draw();




    }

    c->Print(Form("%s/HE%i/%s.pdf", plotDir.Data(), wedge, name.Data()));
    c->Print(Form("%s/HE%i/%s.png", webplotDir.Data(), wedge, name.Data()));
    delete c;
}


void plot_vecTH1F(vector< vector<TH1F*> > h1_vec, TString name, int wedge)
{
    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    if (h1_vec.size() < 7) c->Divide(3, 2);
    else c->Divide(4, 2);
    for (unsigned int ichan = 0; ichan < h1_vec.size(); ichan++)
    {
        //  cout<<"pad number "<<iphi<<endl;
        TPad *pad(NULL);
        pad = static_cast<TPad *>(c->cd(ichan + 1));
        pad->SetLeftMargin(0.16);
        pad->SetBottomMargin(0.14);
        // pad->SetTopMargin(0.16);
        pad->SetGrid();
        pad->SetLogy();

        // TLatex *peak = new TLatex(0.7,0.84,Form("Peak TS: %i",));
        //peak->SetNDC();
        //peak->SetTextSize(textSize);

        ///Brutal hack to remove phi label from top title and apply to legend entries instead
        string title_intersection;
        TString title;

        if (h1_vec[ichan].size() > 1) {
            string t1 = h1_vec[ichan][0]->GetTitle();
            string t2 = h1_vec[ichan][1]->GetTitle();
            //  cout<<"t1 "+t1+" , t2 "+t2<<endl;
            std::set_intersection(t1.begin(), t1.end(), t2.begin(), t2.end(), std::back_inserter(title_intersection));
            title = title_intersection;
        }
        else if (h1_vec[ichan].size() > 0) title = h1_vec[ichan][0]->GetTitle();
        title.ReplaceAll("iPhi 7D", "D");
        title.ReplaceAll("iPhi 6D", "D");
        title.ReplaceAll("iPhi 5D", "D");
        title.ReplaceAll("iPhi 4D", "D");
        title.ReplaceAll("iPhi 3D", "D");
        title.ReplaceAll("iPhi 2D", "D");
        title.ReplaceAll("iPhi 1D", "D");
        title.ReplaceAll("iPhi D", "D");

        title.ReplaceAll("  ", " ");
        //title.ReplaceAll("iPhi 67", "All iPhi");
        // title.ReplaceAll("Depth 8", "All depths");
        float ystart = 0.9 - 0.075 * h1_vec[ichan].size();
        TLegend *leg = new TLegend(0.16, ystart, 0.9, 0.9);

        for (unsigned int ihist = 0; ihist < h1_vec[ichan].size(); ihist++) {
            h1_vec[ichan][ihist]->SetStats(0);
            //gStyle->SetTitleFontSize(0.1);

            h1_vec[ichan][ihist]->GetXaxis()->SetTitleSize(0.07);
            h1_vec[ichan][ihist]->GetXaxis()->SetLabelSize(0.06);
            h1_vec[ichan][ihist]->GetYaxis()->SetTitleSize(0.07);
            h1_vec[ichan][ihist]->GetYaxis()->SetLabelSize(0.06);
            h1_vec[ichan][ihist]->GetYaxis()->SetTitleOffset(1.19);
            //  h1_vec[ichan][ihist]->SetMaximum(20);
            //  h1_vec[ichan][ihist]->SetMinimum(0.0001);

            TString thistitle = h1_vec[ichan][ihist]->GetTitle();
            string phi = getPhi(thistitle);
            //    cout<<"Got phi: "<<phi<<endl;
            //std::set_difference(st.begin(), st.end(), title_intersection.begin(), title_intersection.end(), std::back_inserter(phi));
            TString tphi = "iPhi " + phi;
            // tphi.ReplaceAll(" ", "");
            // tphi.ReplaceAll(",", "");
            // tphi = "iPhi 6" + tphi;
            //tphi.ReplaceAll("iPhi 67", "All iPhi");
            h1_vec[ichan][ihist]->SetTitle(title);
            h1_vec[ichan][ihist]->SetNdivisions(505, "X");

            TString prec = "0";
            if (name.Contains("TS2")) prec = "1";
            TString legentry = tphi + Form(" (N: %.0f, #mu: %." + prec + "f #pm %.1f)", h1_vec[ichan][ihist]->GetEntries(), h1_vec[ichan][ihist]->GetMean(), h1_vec[ichan][ihist]->GetMeanError());
            if (h1_vec[ichan].size() == 1) legentry = Form("N: %.0f, #mu: %." + prec + "f #pm %.1f", h1_vec[ichan][ihist]->GetEntries(), h1_vec[ichan][ihist]->GetMean(), h1_vec[ichan][ihist]->GetMeanError());

            if (h1_vec[ichan][ihist]->GetBinContent(0) > 0) { h1_vec[ichan][ihist]->SetBinContent(1, h1_vec[ichan][ihist]->GetBinContent(0) + h1_vec[ichan][ihist]->GetBinContent(1)); h1_vec[ichan][ihist]->SetBinContent(0, 0);}
            if (h1_vec[ichan][ihist]->GetBinContent(h1_vec[ichan][ihist]->GetNbinsX() + 1) > 0) {
                h1_vec[ichan][ihist]->SetBinContent(h1_vec[ichan][ihist]->GetNbinsX(), h1_vec[ichan][ihist]->GetBinContent(h1_vec[ichan][ihist]->GetNbinsX() + 1) + h1_vec[ichan][ihist]->GetBinContent(h1_vec[ichan][ihist]->GetNbinsX()));
                h1_vec[ichan][ihist]->SetBinContent(h1_vec[ichan][ihist]->GetNbinsX() + 1, 0);
            }

            if (h1_vec[ichan][ihist]->Integral() > 0) h1_vec[ichan][ihist]->Scale(1. / h1_vec[ichan][ihist]->Integral());
            h1_vec[ichan][ihist]->SetAxisRange(0.0005, 250., "Y");
            //  cout <<"Mean before: "<<m<<" mean after: "<<h1_vec[ichan][ihist]->GetMean();
            // cout<<"colors[ihist] "<<colors[ihist]<<endl;
            h1_vec[ichan][ihist]->SetLineColor(colors[ihist]);
            //h1_vec[ichan][ihist]->SetFillColor(colors[ihist]);
            leg->AddEntry(h1_vec[ichan][ihist], legentry, "l");
            if (ihist == 0)
                h1_vec[ichan][ihist]->Draw("hist");
            else
                h1_vec[ichan][ihist]->Draw("hist same");

        }
        leg->Draw();

    }

    c->Print(Form("%s/HE%i/%s.pdf", plotDir.Data(), wedge, name.Data()));
    c->Print(Form("%s/HE%i/%s.png", webplotDir.Data(), wedge, name.Data()));
    delete c;
}



//Plot average pulse shape for 6 channels
void plot_pulse(vector<TH1F*> h1_vec, TString name, int wedge)
{
    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    if (h1_vec.size() < 7) c->Divide(3, 2);
    else c->Divide(4, 2);
    int max = h1_vec.size();
    for (int ichan = 0; ichan < max; ichan++)
    {
        //  cout<<"pad number "<<ichan<<endl;
        TPad *pad(NULL);
        pad = static_cast<TPad *>(c->cd(ichan + 1));
        pad->SetLeftMargin(0.2);
        pad->SetBottomMargin(0.2);
        pad->SetGrid();
        // TLatex *peak = new TLatex(0.7,0.84,Form("Peak TS: %i",));
        //peak->SetNDC();
        //peak->SetTextSize(textSize);

        TString title = h1_vec[ichan]->GetTitle();
        h1_vec[ichan]->SetTitle(title);

        h1_vec[ichan]->SetStats(0);
        //h1_vec[ichan]->SetTitleSize(0.1);
        h1_vec[ichan]->GetXaxis()->SetTitleSize(0.07);
        h1_vec[ichan]->GetXaxis()->SetLabelSize(0.06);
        h1_vec[ichan]->GetYaxis()->SetTitleSize(0.07);
        h1_vec[ichan]->GetYaxis()->SetLabelSize(0.06);
        h1_vec[ichan]->GetYaxis()->SetTitleOffset(1.17);
        h1_vec[ichan]->SetMinimum(0);
        h1_vec[ichan]->Draw("EP");

    }

    c->Print(Form("%s/HE%i/%s.pdf", plotDir.Data(), wedge, name.Data()));
    c->Print(Form("%s/HE%i/%s.png", webplotDir.Data(), wedge, name.Data()));
    delete c;
}

//Plot average pulse shape for 6 channels
void plot_pulse(vector<TH1F*> h1_vec, vector<TH1F*> h1_ped, TString name, int wedge)
{
    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    if (h1_vec.size() < 7) c->Divide(3, 2);
    else c->Divide(4, 2);
    int max = h1_vec.size();
    for (int ichan = 0; ichan < max; ichan++)
    {
        //  cout<<"pad number "<<ichan<<endl;
        TPad *pad(NULL);
        pad = static_cast<TPad *>(c->cd(ichan + 1));
        pad->SetLeftMargin(0.2);
        pad->SetBottomMargin(0.2);
        pad->SetGrid();
        // TLatex *peak = new TLatex(0.7,0.84,Form("Peak TS: %i",));
        //peak->SetNDC();
        //peak->SetTextSize(textSize);

        //// Normalize shapes ///
        /// Shapes in histogram are just the sum of all pulses //
        int ndigis = h1_vec[ichan]->GetEntries() / 8;
        if (ndigis > 0) h1_vec[ichan]->Scale(1. / ndigis);
        ndigis = h1_ped[ichan]->GetEntries() / 8;
        if (ndigis > 0) h1_ped[ichan]->Scale(1. / ndigis);

        TString title = h1_vec[ichan]->GetTitle();
        h1_vec[ichan]->SetTitle(title);
        h1_ped[ichan]->SetTitle(title);
        h1_vec[ichan]->SetStats(0);
        //gStyle->SetTitleFontSize(0.1);
        h1_vec[ichan]->GetXaxis()->SetTitleSize(0.07);
        h1_vec[ichan]->GetXaxis()->SetLabelSize(0.06);
        h1_vec[ichan]->GetYaxis()->SetTitleSize(0.07);
        h1_vec[ichan]->GetYaxis()->SetLabelSize(0.06);
        h1_vec[ichan]->GetYaxis()->SetTitleOffset(1.17);
        h1_vec[ichan]->SetMinimum(0);
        h1_vec[ichan]->Draw("EP");
        h1_ped[ichan]->SetLineColor(kRed);
        h1_ped[ichan]->SetLineStyle(3);
        h1_ped[ichan]->Draw("hist same");

    }

    c->Print(Form("%s/HE%i/%s.pdf", plotDir.Data(), wedge, name.Data()));
    c->Print(Form("%s/HE%i/%s.png", webplotDir.Data(), wedge, name.Data()));
    delete c;
}

void h1cosmetic(TH1F * hist) {


    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.05);
    // hist->GetYaxis()->SetTitleOffset(1.15);


}


void setSummaryStyle() {
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLabelSize(0.05, "xyz");
    //#gStyle->SetTitleFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetTitleFont(42, "t");
    //#gStyle->SetTitleSize(0.05);
    gStyle->SetTitleSize(0.06, "xyz");
    gStyle->SetTitleSize(0.06, "t");
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetTitleOffset(1, "Y");
    //#gStyle->SetLegendTextSize(0.05);
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(13);

}

//This is like the useful python split function
vector<string> split(const string & str, const string & delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos - prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

string getPhi(TString title) {
    //cout << "splitting " << title << endl;
    vector<string> phidepth = split(title.Data(), "iPhi ");
    //cout<<phidepth[0]<<" "<<phidepth[1]<<endl;
    vector<string> phi = split(phidepth[1], ",");
    //cout<<phi[0]<<endl;
    return phi[0];
}


# ifndef __CINT__  // the following code will be invisible for the interpreter
int main(int argc, char **argv)
{
    if (argc == 3) plotDistributions(argv[1], argv[2]);
    else if (argc == 4) plotDistributions(argv[1], argv[2], stoi(argv[3]));
    else if (argc == 5) plotDistributions(argv[1], argv[2], stoi(argv[3]), stoi(argv[4]));
    else cout << "Please give input histogram root file, and dataset tag for output name." << endl;

}
# endif