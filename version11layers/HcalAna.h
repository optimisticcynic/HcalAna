#ifndef HcalAna_H
#define HcalAna_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include "TLorentzVector.h"
#include <TDirectory.h>

#include "AnaInput.h"

// Define the size of array variables. - need to be the same as ntuple maker
#define MAXMU 10
#define MAXJET 20
#define MAXGEN 20
#define mcuts 3
#define MAXBI 50

// definitions of input data strut sizes
const int MAX_N_ISO = 5;
const int MAX_N_LEAVES = 100;
const int MAX_N_VALS = 5;

class HcalAna : public TObject {

public:

   HcalAna( string datacardfile = "DataCard.txt");     
   ~HcalAna();     
   
   void Analysis();
   void makingisoplotj1();
   void HistoWrite( string theFolder , TFile* file ) ;
   void OpenHistogram() ;
   void DrawHistogram( ) ;
   void Draw_H1( vector <TH1D*> h1, string plotname , vector <TH1D*> h2 , string splotname= "" ) ;
   void Draw_H2( vector <TH1D*> h1, string plotname , vector <TH1D*> h2 , string splotname= "" ) ;
   void efficiency(  TH1D* h1 ,  TH1D* h2 ,  double &maxb, double &mins, double &effb, double &blocks);
   void massTableofDoom(int par, int k, float& mesht1, float& mesht2, float& meshht1, float& meshht2);

private:

   AnaInput*     Input ;
   vector <double> mucuts, jcuts;
   vector <int> parameter1, parameter2, parameter3, sparameter1, sparameter2, sparameter3, parameterh1, parameterh2, parameterh3, sparameterh1, sparameterh2, sparameterh3;
   vector <int> nparameter1, nparameter2, nparameter3, snparameter1, snparameter2, snparameter3, nparameterh1, nparameterh2, nparameterh3, snparameterh1, snparameterh2, snparameterh3;
   ////////////////////// JOE WAS HERE
   template<class T>
   struct InputDataStruct
   {
   	  //This is a struct to hold the iso values from the files 
   private:
   
   	  T data[MAX_N_ISO][MAX_N_LEAVES][MAX_N_VALS];
   	  
   public:
   	  //accesssor function 
   	  T operator()(int i, int j, int k)
   	  {
   	  	 if(i < 1 || j < 0 || k < 0 || i > MAX_N_ISO || j >= MAX_N_LEAVES || k >= MAX_N_VALS)
   	  	 {
   	  	     std::cout << "INVALID InputDataStruct Values!!!!!" << std::endl;
   	  	     exit(1);
   	  	 }
   	  	 return data[i - 1][j][k];
   	  }
   	  
   	  T* getArrayPtr(int i)
   	  {
   	  	 if(i < 1 || i > MAX_N_ISO)
   	  	 {
   	  	     std::cout << "INVALID InputDataStruct Valuess!!!!!" << std::endl;
   	  	     exit(1);
   	  	 }
   	     return (T*)data[i - 1];
   	  }
   };
   
   InputDataStruct<float> giso;
   /*float giso1[100][5];//deposted iso reality
   float giso2[100][5];
   float giso3[100][5];
   float giso4[100][5];
   float giso5[100][5];
   */


   InputDataStruct<int> gisoh;
   
   int gisoh1[100][5];//deposted iso hit reality
   int gisoh2[100][5];
   int gisoh3[100][5];
    int gisoh4[100][5];
   int gisoh5[100][5];
   
   
   InputDataStruct<float> riso;
   
   float riso1[100][5];//deposted iso reconstructed
   float riso2[100][5];
   float riso3[100][5];
   float riso4[100][5];
   float riso5[100][5];
   
   
   InputDataStruct<int> risoh;
   int risoh1[100][5];//deposted iso reconstructed
   int risoh2[100][5];
   int risoh3[100][5];
   int risoh4[100][5];
   int risoh5[100][5];
   
   string hfolder ;
   string plotType ;
   string comp1;
   string comp2;
   vector <string> comp;

   //string parx;
   int ProcessEvents ;
   string datafileName ;
   char hfName[16] ;
   
   map<string, float[100][5]>adding;
   
   // Giving variables to hold value from ntuple
   int nMu, nJets, nGen ;
   int bid[MAXBI], mom[MAXBI];
   int hisn;
   int offby;
   int tnum;
   int tgen;
   
   float  genE[10];
   
   float muPx[MAXMU], muPy[MAXMU], muPz[MAXMU], muE[MAXMU], jetE[MAXJET],jetPx[MAXJET], jetPy[MAXJET], jetPz[MAXJET];
   float bPx[MAXBI], bPy[MAXBI], bPz[MAXBI], bE[MAXBI];
   TTree* tr ;



   TFile *theFile ;//file
   TDirectory *jets_directory;
   TDirectory *w_directory;
   TDirectory *t_directory;
   TDirectory *rec_directory;
   TDirectory *compare_directory;
   
   // Define the histograms
   TH1D* h_nMu  ;
   TH1D* h_muPt ;
   TH1D* h_nJets;
   TH1D* h_Jets;
   TH1D* h_JetsPt;
   TH1D* h_gwz;
   TH1D* h_gjets;

   TH1D* g_isolj1;
   TH1D* g_isolj12;//holds real muons isolation energy j jet w wboson decay t total
   TH1D* g_isolj123;
   
   TH1D* g_isoljh1;
   TH1D* g_isoljh12;//holds real muons isolation energy j jet w wboson decay t total
   TH1D* g_isoljh123;
   
   TH1D* g_isolw1;//energy deposit w decay
   TH1D* g_isolw12;
   TH1D* g_isolw123;
   TH1D* g_isolwh1;//hits w decay iso
   TH1D* g_isolwh12;
   TH1D* g_isolwh123;
   
   TH1D* g_isolt1;//total muons
   TH1D* g_isolt12;
   TH1D* g_isolt123;
   TH1D* g_isolth1;//hits
   TH1D* g_isolth12;
   TH1D* g_isolth123;
   
   TH1D* r_isol1;//holds reconstructed muons isolation energy
   TH1D* r_isol12;
   TH1D* r_isol123;
   TH1D* r_isolh1;//holds reconstructed muons isolation energy
   TH1D* r_isolh12;
   TH1D* r_isolh123;
   
   TH1D* g_crashtp;
   TH1D* g_crashjp;
   TH1D* g_crashwp;
   TH1D* g_crashte;
   TH1D* g_crashje;
   TH1D* g_crashwe;
   
   TH1D* g_realte;
   TH1D* g_realje;
   TH1D* g_realwe;
   TH1D* g_realtp;
   TH1D* g_realjp;
   TH1D* g_realwp;
   
   vector <TH1D*> mesh1;//energy deposted
   vector <TH1D*> mesh2;
   
   vector <TH1D*> meshh1;//this is for hits
   vector <TH1D*> meshh2;
   
   //map<string,TH1D*> histos;
   //ClassDef(HcalAna, 1);
};

//#if !defined(__CINT__)
//    ClassImp(HcalAna);
#endif

