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

//map IEta -29 to 29 to histogram indices 0 to 58 (index = IEta - minieta)
//will skip IEta -15 to 15 (indices 14 to 44)
const int minieta = -29;
const int nieta = 59;

//map IPhi 1 to 72 to histogram indices 0 to 71 (index = IPhi - miniphi)
const int miniphi = 1;
const int niphi = 72;

//map Depth 1 to 7 to histogram indices 0 to 6 (index = depth - mindepth)
const int mindepth = 1;
const int ndepth = 7;

const int nfiber = 24;
//Number of charge bins
//bin 0: inclusive
//bin 1: 5000-7000 fC
//bin 2: 7000-11000 fC
//bin 3: 11000+
const int nfC = 4;
const int fCbins[] = {0, 5000, 7000, 11000};

//Store average pedestal for each channel and capacitor
float peds[nieta][niphi][ndepth][4];

//Compile like this
// g++ -o makeTimingHistsHE makeTimingHistsHE.C  `root-config --cflags --glibs`
void makeTimingHistsHE(TString inputDir, TString outputTag, int pedestalRun);
# ifndef __CINT__  // the following code will be invisible for the interpreter
int main(int argc, char **argv)
{
    if (argc == 4) makeTimingHistsHE(argv[1], argv[2], stoi(argv[3]));
    else cout << "Please give input directory, tag for output name, and pedestal run number." << endl;

}
# endif

//Find sum of ADC from max-1 to max+2
float findTotalADC(vector<int> ADC) {
    int max_idx = distance(ADC.begin(), max_element(ADC.begin(), ADC.end()));
    //cout<<"Max index is "<<max_idx<<endl;
    int min = max_idx - 1;
    int max = max_idx + 2;
    if (min < 0) min = 0;
    if (max > 7) max = 7;
    float total = 0.;
    for (int i = min; i <= max; i++) {
        //cout<<"i is "<<i<<", fC is "<<ADC.at(i)<<endl;
        total += ADC.at(i);
    }
    // cout<<"total is "<<total<<endl;
    return total;

}

//Find sum of fC from max-1 to max+2
float findTotalfC(vector<float> fC_ped_subt) {
    int max_idx = distance(fC_ped_subt.begin(), max_element(fC_ped_subt.begin(), fC_ped_subt.end()));
    //cout<<"Max index is "<<max_idx<<endl;
    int min = max_idx - 1;
    int max = max_idx + 2;
    if (min < 0) min = 0;
    if (max > 7) max = 7;
    float total = 0.;
    for (int i = min; i <= max; i++) {
        //cout<<"i is "<<i<<", fC is "<<fC_ped_subt.at(i)<<endl;
        total += fC_ped_subt.at(i);
    }
    // cout<<"total is "<<total<<endl;
    return total;

}


void loadPedestals(int pedestalRun) {

    TString pedname = Form("pedestals/ped_%i.txt", pedestalRun);
    cout << "Using pedestal table " << pedname << endl;
    std::ifstream pedFile(pedname);
    std::string pedline;
   // set<TString> pedchannels;
  //  set<int> intpedchans;
    while (getline(pedFile, pedline)) {
        if (!pedline.length()) continue;
        std::istringstream iss(pedline);
        int eta, phi, dep;
        string det, ID; 
        float cap1, cap2, cap3, cap4;
        iss >> eta >> phi >> dep >> det >> cap1 >> cap2 >> cap3 >> cap4 >> ID;
        //cout << Form("eta %i, phi%i, dep %i, cap1 %.1f, cap2 %.1f, cap3 %.1f, cap4 %.1f", eta, phi, dep, cap1, cap2, cap3, cap4)<< endl;
        peds[eta - minieta][phi - miniphi][dep - mindepth][0] = cap1;
        peds[eta - minieta][phi - miniphi][dep - mindepth][1] = cap2;
        peds[eta - minieta][phi - miniphi][dep - mindepth][2] = cap3;
        peds[eta - minieta][phi - miniphi][dep - mindepth][3] = cap4;
        //pedchannels.insert(Form("iEta %i, iPhi %i, Depth %i", eta, phi, dep));
        //intpedchans.insert(dep + 10 * (phi) + 1000 * (eta));

    }
}


void makeTimingHistsHE(TString inputDir, TString outputTag, int pedestalRun) {

    gStyle->SetGridColor(16);


    loadPedestals(pedestalRun);

    TChain *tree = new TChain("hcalTupleTree/tree");
    tree->Add(inputDir + "/*.root");

    int nevents = tree->GetEntries();
    int nevents_with_digis = tree->GetEntries("Length$(QIE11DigiIEta)>0");
    cout << Form("Number of events: %i", nevents) << endl;
    cout << Form("Number of events with digis: %i", nevents_with_digis) << endl;


    // event,ls and run
    UInt_t   event = 0;
    tree->SetBranchAddress("event", &event);
    UInt_t   ls = 0;
    tree->SetBranchAddress("ls", &ls);
    UInt_t   run = 0;
    tree->SetBranchAddress("run", &run);
    UInt_t   bx = 0;
    tree->SetBranchAddress("bx", &bx);

    // QIE11

    //vector<int>   *QIE11DigiSubdet = 0;
    //tree->SetBranchAddress("QIE11DigiSubdet", &QIE11DigiSubdet);
    vector<int>   *QIE11DigiIEta = 0;
    tree->SetBranchAddress("QIE11DigiIEta", &QIE11DigiIEta);
    vector<int>   *QIE11DigiIPhi = 0;
    tree->SetBranchAddress("QIE11DigiIPhi", &QIE11DigiIPhi);
    vector<int>   *QIE11DigiDepth = 0;
    tree->SetBranchAddress("QIE11DigiDepth", &QIE11DigiDepth);
    vector<vector<int> >   *QIE11DigiCapID = 0;
    tree->SetBranchAddress("QIE11DigiCapID", &QIE11DigiCapID);
    vector<vector<float> >   *QIE11DigiFC = 0;
    tree->SetBranchAddress("QIE11DigiFC", &QIE11DigiFC);
    //vector<vector<int> >   *QIE11DigiTDC = 0;
    //tree->SetBranchAddress("QIE11DigiTDC", &QIE11DigiTDC);
    vector<vector<int> >   *QIE11DigiADC = 0;
    tree->SetBranchAddress("QIE11DigiADC", &QIE11DigiADC);
    //vector<vector<int> >   *QIE11DigiSOI = 0;
    //tree->SetBranchAddress("QIE11DigiSOI", &QIE11DigiSOI);

    vector<float>   *QIE11DigiTimeTDC = 0;
    tree->SetBranchAddress("QIE11DigiTimeTDC", &QIE11DigiTimeTDC);
    vector<float>   *QIE11DigiTimeFC = 0;
    tree->SetBranchAddress("QIE11DigiTimeFC", &QIE11DigiTimeFC);
    vector<float>   *QIE11DigiTotFC = 0;
    tree->SetBranchAddress("QIE11DigiTotFC", &QIE11DigiTotFC);

    vector<int>   *QIE11DigiNTDC = 0;
    tree->SetBranchAddress("QIE11DigiNTDC", &QIE11DigiNTDC);

    cout << "Finished loading branches." << endl;

    TString outFilename = "hists/hists_" + outputTag + ".root";
    TFile* outFile = new TFile(outFilename, "RECREATE");

    //Pulse shapes
    TH1F *h1_fC[nieta + 2][niphi + 1][ndepth][nfC]; //Hists to store average fC in each time slice
    TH1F *h1_ped[nieta + 2][niphi + 1][ndepth][nfC]; //Hists to store average pedestal in each time slice
    TH1F *h1_energy[nieta + 2][niphi + 1][ndepth][nfC]; //energy distributions (charge distributions)
    TH1F *h1_energyADC[nieta + 2][niphi + 1][ndepth][nfC]; //energy distributions (charge distributions)
    TH1F *h1_chg_time[nieta + 2][niphi + 1][ndepth][nfC]; //Filled with charge-weighted average (uses peak finding)
    TH1F *h1_TDC_time[nieta + 2][niphi + 1][ndepth][nfC]; //in ns
    TH1F *h1_nTDC[nieta + 2][niphi + 1][ndepth][nfC]; //Number of times TDC fires in one digi
    TH1F *h1_chgfracTS2[nieta + 2][niphi + 1][ndepth][nfC];
    TH1F *h1_chgfracTS4[nieta + 2][niphi + 1][ndepth][nfC];
    TH2F *h2_TotFCvsTDC[nieta + 2][niphi + 1][ndepth][nfC];
    TH2F *h2_chgfracTS2vsTDC[nieta + 2][niphi + 1][ndepth][nfC]; //ns
    TH2F *h2_chgfracTS4vsTDC[nieta + 2][niphi + 1][ndepth][nfC]; //ns

    /// index nieta: integration over eta, HEM
    /// index nieta+1: integration over eta, HEP
    /// index niphi: integration over phi




    //Range for charge distributions
    int enbins = 40;
    int eminx = 5000;
    int emaxx = 25000;

    cout << "Booking histograms." << endl;
    for (int ieta = 0; ieta < nieta + 2; ieta++)
    {
        //////// HE only ////////
        if (ieta >= 14 && ieta <= 44) continue;

        for (int iphi = 0; iphi < niphi + 1; iphi++)
        {

            for (int idepth = 0; idepth < ndepth; idepth++) {

                for (int ifC = 0; ifC < nfC; ifC++) {

                    TString title;
                    if (ifC == 0) title = Form("iEta %i, iPhi %i, Depth %i, Q > %i fC;", ieta + minieta, iphi + miniphi, idepth + mindepth, fCbins[1]);
                    else if (ifC + 1 < nfC) title = Form("iEta %i, iPhi %i, Depth %i, %i < Q < %i fC;", ieta + minieta, iphi + miniphi, idepth + mindepth, fCbins[ifC], fCbins[ifC + 1]);
                    else  title = Form("iEta %i, iPhi %i, Depth %i, Q > %i fC;", ieta + minieta, iphi + miniphi, idepth + mindepth, fCbins[ifC]);
                    TString name = Form("_ieta%i_iphi%i_idepth%i_fC%i", ieta + minieta, iphi + miniphi, idepth + mindepth, ifC);


                    //Pulse shapes- sum of ALL digis in each channel. Will have to normalize when plotting
                    h1_fC[ieta][iphi][idepth][ifC] = new TH1F("h1_fC" + name, "h1_fC" + name, 8, -0.5, 7.5);
                    h1_fC[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_fC[ieta][iphi][idepth][ifC]->SetTitle(title + "TS; Average charge [fC]");

                    h1_ped[ieta][iphi][idepth][ifC] = new TH1F("h1_ped" + name, "h1_ped" + name, 8, -0.5, 7.5);
                    h1_ped[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_ped[ieta][iphi][idepth][ifC]->SetTitle(title + "TS; Average charge [fC]");

                    h1_TDC_time[ieta][iphi][idepth][ifC] = new TH1F("h1_TDC_time" + name, "h1_TDC_time" + name, 50, 40, 140);
                    h1_TDC_time[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_TDC_time[ieta][iphi][idepth][ifC]->SetTitle(title + "TDC time [ns]; Fraction of hits");

                    h1_nTDC[ieta][iphi][idepth][ifC] = new TH1F("h1_nTDC" + name, "h1_nTDC" + name, 50, 40, 140);
                    h1_nTDC[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_nTDC[ieta][iphi][idepth][ifC]->SetTitle(title + "Number of TDC fires; Fraction of hits");

                    h1_chg_time[ieta][iphi][idepth][ifC] = new TH1F("h1_chg_time" + name, "h1_chg_time" + name, 50, 40, 140);
                    h1_chg_time[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_chg_time[ieta][iphi][idepth][ifC]->SetTitle(title + "Charge-averaged time [ns]; Fraction of hits");

                    //energy = total charge
                    h1_energy[ieta][iphi][idepth][ifC] = new TH1F("h1_energy" + name, "h1_energy" + name, enbins, eminx, emaxx);
                    h1_energy[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_energy[ieta][iphi][idepth][ifC]->SetTitle(title + "Total charge, 4 TS around peak [fC]; Fraction of hits");
                    h1_energy[ieta][iphi][idepth][ifC]->StatOverflows(kTRUE);

                    h1_energyADC[ieta][iphi][idepth][ifC] = new TH1F("h1_energyADC" + name, "h1_energyADC" + name, enbins, 0, 500);
                    h1_energyADC[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_energyADC[ieta][iphi][idepth][ifC]->SetTitle(title + "Total charge, 4 TS around peak [fC]; Fraction of hits");
                    h1_energyADC[ieta][iphi][idepth][ifC]->StatOverflows(kTRUE);

                    h1_chgfracTS2[ieta][iphi][idepth][ifC] = new TH1F("h1_chgfracTS2" + name, "h1_chgfracTS2" + name, 30, 0, 15);
                    h1_chgfracTS2[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_chgfracTS2[ieta][iphi][idepth][ifC]->StatOverflows(kTRUE);
                    h1_chgfracTS2[ieta][iphi][idepth][ifC]->SetTitle(title + "TS2 charge fraction [%]; Fraction of hits");

                    h1_chgfracTS4[ieta][iphi][idepth][ifC] = new TH1F("h1_chgfracTS4" + name, "h1_chgfracTS4" + name, 35, 0, 70);
                    h1_chgfracTS4[ieta][iphi][idepth][ifC]->Sumw2();
                    h1_chgfracTS4[ieta][iphi][idepth][ifC]->StatOverflows(kTRUE);
                    h1_chgfracTS4[ieta][iphi][idepth][ifC]->SetTitle(title + "TS4 charge fraction [%]; Fraction of hits");

                    h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC] = new TH2F("h2_chgfracTS2vsTDC" + name, "h2_chgfracTS2vsTDC" + name, 50, 40, 140, 30, 0, 15);
                    h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]->Sumw2();
                    h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]->SetTitle(title + "TDC time [ns];TS2 charge fraction [%]");   //ns

                    h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC] = new TH2F("h2_chgfracTS4vsTDC" + name, "h2_chgfracTS4vsTDC" + name, 50, 40, 140, 35, 0, 70);
                    h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]->Sumw2();
                    h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]->SetTitle(title + "TDC time [ns];TS4 charge fraction [%]");   //ns

                    h2_TotFCvsTDC[ieta][iphi][idepth][ifC] = new TH2F("h2_TotFCvsTDC" + name, "h2_TotFCvsTDC" + name, 50, 40, 140, enbins, eminx, emaxx);
                    h2_TotFCvsTDC[ieta][iphi][idepth][ifC]->Sumw2();
                    h2_TotFCvsTDC[ieta][iphi][idepth][ifC]->SetTitle(title + "TDC time [ns]; Total charge, 4 TS around peak [fC]");


                }
            }
        }
    }

    //Start event loop
    cout << "Starting event loop" << endl;
    for (int ievent = 0; ievent < nevents; ievent++)
    {
        //if (ievent % 1 == 0) cout << "Processed " << ievent << " events" << endl;
        if (ievent % 100000 == 0) cout << "Processed " << ievent << " events" << endl;
        tree->GetEntry(ievent);

        //For each event, loop over digis
        for (unsigned int idigi = 0; idigi < QIE11DigiIEta->size(); idigi++) {
            //Get indices for histogram vectors
            int ieta =  QIE11DigiIEta->at(idigi) - minieta;
            int iphi =  QIE11DigiIPhi->at(idigi) - miniphi;
            int idepth = QIE11DigiDepth->at(idigi) - mindepth;


            vector<float> QIE11DigiFC_ped_subt;
            //For each digi, loop over all time slices
            for(unsigned int its = 0; its<QIE11DigiFC->at(idigi).size(); its++){
              float fC = QIE11DigiFC->at(idigi).at(its);
              int capid = QIE11DigiCapID->at(idigi).at(its); 
              float fC_ped_subt = fC - peds[ieta][iphi][idepth][capid];
              if(peds[ieta][iphi][idepth][capid]== 0) cout<<"Missing pedestal value, eta phi depth: "<<QIE11DigiIEta->at(idigi)<<" "<<QIE11DigiIPhi->at(idigi)<<" "<<QIE11DigiDepth->at(idigi)<<endl;
              QIE11DigiFC_ped_subt.push_back(fC_ped_subt);
            }



            float totFC = findTotalfC(QIE11DigiFC_ped_subt); //= QIE11DigiTotFC->at(idigi);
            float timeTDC = QIE11DigiTimeTDC->at(idigi);
            int nTDC = QIE11DigiNTDC->at(idigi);
            float timeFC = QIE11DigiTimeFC->at(idigi);
            float chgfracTS2 = 100.*QIE11DigiFC_ped_subt.at(2) / totFC;
            float chgfracTS4 = 100.*QIE11DigiFC_ped_subt.at(4) / totFC;

            //Find charge bin
            int ifC = 0;
            if (totFC >= fCbins[1] && totFC <= fCbins[2]) ifC = 1;
            else if (totFC > fCbins[2] && totFC <= fCbins[3]) ifC = 2;
            else if (totFC > fCbins[3]) ifC = 3;

            if (ifC == 0) continue;

            h1_energy[ieta][iphi][idepth][ifC]->Fill(totFC);
            h1_energyADC[ieta][iphi][idepth][ifC]->Fill(findTotalADC(QIE11DigiADC->at(idigi)));
            h1_chg_time[ieta][iphi][idepth][ifC]->Fill(timeFC);
            h1_chgfracTS2[ieta][iphi][idepth][ifC]->Fill(chgfracTS2);
            h1_chgfracTS4[ieta][iphi][idepth][ifC]->Fill(chgfracTS4);

            //Fill average pulse shape (now that we know the charge bin)
            for (unsigned int its = 0; its < QIE11DigiFC->at(idigi).size(); its++) {
                // cout<<ieta<<" "<<iphi<<" "<<idepth<<" "<<ifC<<" "<<its<<" "<<QIE11DigiFC->at(idigi).at(its)<<endl;
                //Don't subtract pedestal for the pulse shape histogram- instead overlay with measured pedestals as validation
                int capid = QIE11DigiCapID->at(idigi).at(its); 
                h1_fC[ieta][iphi][idepth][ifC]->Fill(its, QIE11DigiFC->at(idigi).at(its)); 
                h1_ped[ieta][iphi][idepth][ifC]->Fill(its, peds[ieta][iphi][idepth][capid]); 
                
            }

            //TDC info
            h1_nTDC[ieta][iphi][idepth][ifC]->Fill(nTDC);
            if (nTDC == 1) {
                h1_TDC_time[ieta][iphi][idepth][ifC]->Fill(timeTDC);
                h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]->Fill(timeTDC, chgfracTS2);
                h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]->Fill(timeTDC, chgfracTS4);
                h2_TotFCvsTDC[ieta][iphi][idepth][ifC]->Fill(timeTDC, totFC);
            }
        }
    }
    //End event loop
    cout << "Finished event loop" << endl;

    //Perform integrations for inclusive fC histograms
    for (int ieta = 0; ieta < nieta; ieta++) {
        if (ieta >= 14 && ieta <= 44) continue;
        for (int iphi = 0; iphi < niphi; iphi++) {
            for (int idepth = 0; idepth < ndepth; idepth++) {
                for (int ifC = 1; ifC < nfC; ifC++) {

                    //Bin 0 is not filled; reserved for inclusive charge bin
                    h1_energy[ieta][iphi][idepth][0]->Add(h1_energy[ieta][iphi][idepth][ifC]);
                    h1_energyADC[ieta][iphi][idepth][0]->Add(h1_energyADC[ieta][iphi][idepth][ifC]);
                    h1_fC[ieta][iphi][idepth][0]->Add(h1_fC[ieta][iphi][idepth][ifC]);
                    h1_ped[ieta][iphi][idepth][0]->Add(h1_ped[ieta][iphi][idepth][ifC]);
                    h1_TDC_time[ieta][iphi][idepth][0]->Add(h1_TDC_time[ieta][iphi][idepth][ifC]);
                    h1_nTDC[ieta][iphi][idepth][0]->Add(h1_nTDC[ieta][iphi][idepth][ifC]);
                    h1_chg_time[ieta][iphi][idepth][0]->Add(h1_chg_time[ieta][iphi][idepth][ifC]);
                    h1_chgfracTS2[ieta][iphi][idepth][0]->Add(h1_chgfracTS2[ieta][iphi][idepth][ifC]);
                    h1_chgfracTS4[ieta][iphi][idepth][0]->Add(h1_chgfracTS4[ieta][iphi][idepth][ifC]);

                    h2_chgfracTS2vsTDC[ieta][iphi][idepth][0]->Add(h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]);
                    h2_chgfracTS4vsTDC[ieta][iphi][idepth][0]->Add(h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]);
                    h2_TotFCvsTDC[ieta][iphi][idepth][0]->Add(h2_TotFCvsTDC[ieta][iphi][idepth][ifC]);

                }

                //Perform integrations in phi and eta
                // "nth" index is reserved for combined
                //integrate eta separately for HEM and HEP (nieta and nieta+1)
                for (int ifC = 1; ifC < nfC; ifC++) {
                    h1_energy[ieta][niphi][idepth][ifC]->Add(h1_energy[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_energy[nieta][iphi][idepth][ifC]->Add(h1_energy[ieta][iphi][idepth][ifC]);
                    else h1_energy[nieta + 1][iphi][idepth][ifC]->Add(h1_energy[ieta][iphi][idepth][ifC]);

                    h1_energyADC[ieta][niphi][idepth][ifC]->Add(h1_energyADC[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_energyADC[nieta][iphi][idepth][ifC]->Add(h1_energyADC[ieta][iphi][idepth][ifC]);
                    else h1_energyADC[nieta + 1][iphi][idepth][ifC]->Add(h1_energyADC[ieta][iphi][idepth][ifC]);

                    h1_fC[ieta][niphi][idepth][ifC]->Add(h1_fC[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_fC[nieta][iphi][idepth][ifC]->Add(h1_fC[ieta][iphi][idepth][ifC]);
                    else h1_fC[nieta + 1][iphi][idepth][ifC]->Add(h1_fC[ieta][iphi][idepth][ifC]);

                    h1_ped[ieta][niphi][idepth][ifC]->Add(h1_ped[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_ped[nieta][iphi][idepth][ifC]->Add(h1_ped[ieta][iphi][idepth][ifC]);
                    else h1_ped[nieta + 1][iphi][idepth][ifC]->Add(h1_ped[ieta][iphi][idepth][ifC]);

                    h1_TDC_time[ieta][niphi][idepth][ifC]->Add(h1_TDC_time[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_TDC_time[nieta][iphi][idepth][ifC]->Add(h1_TDC_time[ieta][iphi][idepth][ifC]);
                    else h1_TDC_time[nieta + 1][iphi][idepth][ifC]->Add(h1_TDC_time[ieta][iphi][idepth][ifC]);

                    h1_nTDC[ieta][niphi][idepth][ifC]->Add(h1_nTDC[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_nTDC[nieta][iphi][idepth][ifC]->Add(h1_nTDC[ieta][iphi][idepth][ifC]);
                    else h1_nTDC[nieta + 1][iphi][idepth][ifC]->Add(h1_nTDC[ieta][iphi][idepth][ifC]);

                    h1_chg_time[ieta][niphi][idepth][ifC]->Add(h1_chg_time[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_chg_time[nieta][iphi][idepth][ifC]->Add(h1_chg_time[ieta][iphi][idepth][ifC]);
                    else h1_chg_time[nieta + 1][iphi][idepth][ifC]->Add(h1_chg_time[ieta][iphi][idepth][ifC]);

                    h1_chgfracTS2[ieta][niphi][idepth][ifC]->Add(h1_chgfracTS2[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_chgfracTS2[nieta][iphi][idepth][ifC]->Add(h1_chgfracTS2[ieta][iphi][idepth][ifC]);
                    else h1_chgfracTS2[nieta + 1][iphi][idepth][ifC]->Add(h1_chgfracTS2[ieta][iphi][idepth][ifC]);

                    h1_chgfracTS4[ieta][niphi][idepth][ifC]->Add(h1_chgfracTS4[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h1_chgfracTS4[nieta][iphi][idepth][ifC]->Add(h1_chgfracTS4[ieta][iphi][idepth][ifC]);
                    else h1_chgfracTS4[nieta + 1][iphi][idepth][ifC]->Add(h1_chgfracTS4[ieta][iphi][idepth][ifC]);

                    h2_chgfracTS2vsTDC[ieta][niphi][idepth][0]->Add(h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h2_chgfracTS2vsTDC[nieta][iphi][idepth][0]->Add(h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]);
                    else h2_chgfracTS2vsTDC[nieta + 1][iphi][idepth][0]->Add(h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]);

                    h2_chgfracTS4vsTDC[ieta][niphi][idepth][0]->Add(h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h2_chgfracTS4vsTDC[nieta][iphi][idepth][0]->Add(h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]);
                    else h2_chgfracTS4vsTDC[nieta + 1][iphi][idepth][0]->Add(h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]);

                    h2_TotFCvsTDC[ieta][niphi][idepth][0]->Add(h2_TotFCvsTDC[ieta][iphi][idepth][ifC]);
                    if (ieta + minieta < 0) h2_TotFCvsTDC[nieta][iphi][idepth][0]->Add(h2_TotFCvsTDC[ieta][iphi][idepth][ifC]);
                    else h2_TotFCvsTDC[nieta + 1][iphi][idepth][0]->Add(h2_TotFCvsTDC[ieta][iphi][idepth][ifC]);
                }
            }
        }
    }

    for (int ieta = 0; ieta < nieta + 2; ieta++) {
        if (ieta >= 14 && ieta <= 44) continue;
        for (int iphi = 0; iphi < niphi + 1; iphi++) {
            for (int idepth = 0; idepth < ndepth; idepth++) {

                //Save to disk (if non-empty)
                for (int ifC = 0; ifC < nfC; ifC++) {
                    int nevents_this_chan = h1_energy[ieta][iphi][idepth][ifC]->GetEntries();
                    if (nevents_this_chan > 0) {
                        h1_fC[ieta][iphi][idepth][ifC]->Write();
                        h1_ped[ieta][iphi][idepth][ifC]->Write();
                        h1_energy[ieta][iphi][idepth][ifC]->Write();
                        h1_energyADC[ieta][iphi][idepth][ifC]->Write();
                        h1_TDC_time[ieta][iphi][idepth][ifC]->Write();
                        h1_nTDC[ieta][iphi][idepth][ifC]->Write();
                        h1_chg_time[ieta][iphi][idepth][ifC]->Write();
                        h1_chgfracTS2[ieta][iphi][idepth][ifC]->Write();
                        h1_chgfracTS4[ieta][iphi][idepth][ifC]->Write();
                        h2_chgfracTS2vsTDC[ieta][iphi][idepth][ifC]->Write();
                        h2_chgfracTS4vsTDC[ieta][iphi][idepth][ifC]->Write();
                        h2_TotFCvsTDC[ieta][iphi][idepth][ifC]->Write();
                    }
                }
            }
        }
    }
    outFile->Close();
}