#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <math.h>
#include <cmath> 
#include <map>
#include "Math/GenVector/VectorUtil.h"
//#include "hDraw.h"

#include "HcalAna.h"

HcalAna::HcalAna(string datacardfile, int u, int y, int t)
{
    Input = new AnaInput(datacardfile);

    Input->GetParameters("PlotType", &plotType);
    Input->GetParameters("Path", &hfolder);
    Input->GetParameters("ProcessEvents", &ProcessEvents);
    Input->GetParameters("TheData", &datafileName);
    //Input->GetParameters("HistoName",     &hfName ) ;
    Input->GetParameters("MuonCuts", &mucuts);
    Input->GetParameters("JetCuts", &jcuts);
    cout << "getting parameter" << endl;
    Input->GetParameters("parameter1", &parameter1);
    Input->GetParameters("parameter2", &parameter2);
    Input->GetParameters("parameter3", &parameter3);
    Input->GetParameters("sparameter1", &sparameter1);
    Input->GetParameters("sparameter2", &sparameter2);
    Input->GetParameters("sparameter3", &sparameter3);


    Input->GetParameters("parameterh1", &parameterh1);
    Input->GetParameters("parameterh2", &parameterh2);
    Input->GetParameters("parameterh3", &parameterh3);
    Input->GetParameters("sparameterh1", &sparameterh1);
    Input->GetParameters("sparameterh2", &sparameterh2);
    Input->GetParameters("sparameterh3", &sparameterh3);


    Input->GetParameters("nparameter1", &nparameter1);
    Input->GetParameters("nparameter2", &nparameter2);
    Input->GetParameters("nparameter3", &nparameter3);
    Input->GetParameters("snparameter1", &snparameter1);
    Input->GetParameters("snparameter2", &snparameter2);
    Input->GetParameters("snparameter3", &snparameter3);


    Input->GetParameters("nparameterh1", &nparameterh1);
    Input->GetParameters("nparameterh2", &nparameterh2);
    Input->GetParameters("nparameterh3", &nparameterh3);
    Input->GetParameters("snparameterh1", &snparameterh1);
    Input->GetParameters("snparameterh2", &snparameterh2);
    Input->GetParameters("snparameterh3", &snparameterh3);
    Input->GetParameters("histonum", &hisn);
    Input->GetParameters("classification", &rec);
    Input->GetParameters("type", &type);
    Input->GetParameters("innercone", &asub);
    Input->GetParameters("maxdr", &drval);

    maxdr = drval;
    string clever(8, '0');
    sprintf(&clever[0], "_dr_%f", drval);
    if(sparameter1[0] == 4) rec = 1;



    if(asub == 0) sub = "";
    cout << "asub is " << asub << endl;
    if(asub == 1)
    {
        sub = "_sub" + clever;
        cout << "sub is " << sub << endl;
    }
    if(asub == 2) sub = "_subd" + clever;
    if(asub == 3) sub = "_subu" + clever;


    if(rec == 1) ofile = "reco" + sub;
    else ofile = "gen" + sub;
    string typ1, typ2;
    if(parameter1[0] == 1) typ1 = "wdecay";
    if(parameter1[0] == 2) typ1 = "jet";
    if(parameter1[0] == 3) typ1 = "total";
    if(parameter1[0] == 4) typ1 = "rec";
    if(sparameter1[0] == 1) typ2 = "wdecay";
    if(sparameter1[0] == 2) typ2 = "jet";
    if(sparameter1[0] == 3) typ2 = "total";
    if(sparameter1[0] == 4) typ2 = "rec";


    ofile = ofile + typ1 + "vrs" + typ2;

    if(type == 13) ifile = "muons";
    if(type == 11) ifile = "electron";

    isreco = (rec == 1);
    if(u >= 0)parameter1[1] = u;
    if(y >= 0)parameter2[1] = y;
    if(t >= 0)parameter3[1] = t;

    ///Hey file name hopefully

    //////more lazy work  
    snparameterh1[1] = snparameter1[1] = nparameterh1[1] = nparameter1[1];
    snparameterh2[1] = snparameter2[1] = nparameterh2[1] = nparameter2[1];
    snparameterh3[1] = snparameter3[1] = nparameterh3[1] = nparameter3[1];

    sparameterh1[1] = sparameter1[1] = parameterh1[1] = parameter1[1];
    sparameterh2[1] = sparameter2[1] = parameterh2[1] = parameter2[1];
    sparameterh3[1] = sparameter3[1] = parameterh3[1] = parameter3[1];

    nparameterh1[0] = nparameter1[0] = parameterh1[0] = parameter1[0];
    snparameterh1[0] = snparameter1[0] = sparameterh1[0] = sparameter1[0];

    parameter3[0] = parameter2[0] = parameter1[0];
    sparameter3[0] = sparameter2[0] = sparameter1[0];
    parameterh3[0] = parameterh2[0] = parameterh1[0];
    sparameterh3[0] = sparameterh2[0] = sparameterh1[0];

    nparameter3[0] = nparameter2[0] = nparameter1[0];
    snparameter3[0] = snparameter2[0] = snparameter1[0];
    nparameterh3[0] = nparameterh2[0] = nparameterh1[0];
    snparameterh3[0] = snparameterh2[0] = snparameterh1[0];


    ///////////////// LAZY AS WAY OF GETTING OUT OF WORK


    if(parameter1[1] == 0) parameter1[0] = 0;
    if(parameter2[1] == 0) parameter2[0] = 0;
    if(parameter3[1] == 0) parameter3[0] = 0;

    if(sparameter1[1] == 0) sparameter1[0] = 0;
    if(sparameter2[1] == 0) sparameter2[0] = 0;
    if(sparameter3[1] == 0) sparameter3[0] = 0;

    if(parameterh1[1] == 0) parameterh1[0] = 0;
    if(parameterh2[1] == 0) parameterh2[0] = 0;
    if(parameterh3[1] == 0) parameterh3[0] = 0;

    if(sparameterh1[1] == 0) sparameterh1[0] = 0;
    if(sparameterh2[1] == 0) sparameterh2[0] = 0;
    if(sparameterh3[1] == 0) sparameterh3[0] = 0;

    if(nparameter1[1] == 0) nparameter1[0] = 0;
    if(nparameter2[1] == 0) nparameter2[0] = 0;
    if(nparameter3[1] == 0) nparameter3[0] = 0;

    if(snparameter1[1] == 0) snparameter1[0] = 0;
    if(snparameter2[1] == 0) snparameter2[0] = 0;
    if(snparameter3[1] == 0) snparameter3[0] = 0;

    if(nparameterh1[1] == 0) nparameterh1[0] = 0;
    if(nparameterh2[1] == 0) nparameterh2[0] = 0;
    if(nparameterh3[1] == 0) nparameterh3[0] = 0;

    if(snparameterh1[1] == 0) snparameterh1[0] = 0;
    if(snparameterh2[1] == 0) snparameterh2[0] = 0;
    if(snparameterh3[1] == 0) snparameterh3[0] = 0;


    ////////////

    cout << "parameter1[0] " << parameter1[0] << endl;
    cout << "parameterh1[1] " << parameterh1[1] << endl;
    cout << "parameter2[0] " << parameter2[0] << endl;
    cout << "parameterh2[1] " << parameterh2[1] << endl;
    cout << "parameter3[0] " << parameter3[0] << endl;
    cout << "parameterh3[1] " << parameterh3[1] << endl;


    Input->GetParameters("comp1", &comp1);
    Input->GetParameters("comp2", &comp2);

    comp.push_back(comp1);

    // Set the address for the branches

    tr = Input->GetTree(datafileName, "HcalUpgrade");
    cout << datafileName << std::endl;


    tr->SetBranchAddress("lepPx", lepPx);
    tr->SetBranchAddress("lepPy", lepPy);
    tr->SetBranchAddress("lepPz", lepPz);
    tr->SetBranchAddress("lepE", lepE);


    tr->SetBranchAddress("nJets", &nJets); //so number of jets
    tr->SetBranchAddress("jetE", jetE);
    tr->SetBranchAddress("jetPx", jetPx);
    tr->SetBranchAddress("jetPy", jetPy);
    tr->SetBranchAddress("jetPz", jetPz);


    tr->SetBranchAddress("nGen", &nGen); //find number of events byprodects that are detected
    tr->SetBranchAddress("pdgId", &bid); //baby ID
    tr->SetBranchAddress("momId", &mom); //moms ID

    tr->SetBranchAddress("genPx", bPx); //babies shit
    tr->SetBranchAddress("genPy", bPy);
    tr->SetBranchAddress("genPz", bPz);
    tr->SetBranchAddress("genE", bE);
    tr->SetBranchAddress("gen_hcalE", gen_hcalE);
    tr->SetBranchAddress("gen_hcaldR", gen_hcaldR);

    //So this should be the isolation shit reality muons
    if(isreco)
    {
        tr->SetBranchAddress("lepIso1", iso.getArrayPtr(1));
        tr->SetBranchAddress("lepIso2", iso.getArrayPtr(2));
        tr->SetBranchAddress("lepIso3", iso.getArrayPtr(3));
        tr->SetBranchAddress("lepIso4", iso.getArrayPtr(4));
        tr->SetBranchAddress("lepIso5", iso.getArrayPtr(5));

        tr->SetBranchAddress("lepIhit1", isoh.getArrayPtr(1));
        tr->SetBranchAddress("lepIhit2", isoh.getArrayPtr(2));
        tr->SetBranchAddress("lepIhit3", isoh.getArrayPtr(3));
        tr->SetBranchAddress("lepIhit4", isoh.getArrayPtr(4));
        tr->SetBranchAddress("lepIhit5", isoh.getArrayPtr(5));
    } else
    {
        tr->SetBranchAddress("genIso1", iso.getArrayPtr(1)); //energy
        tr->SetBranchAddress("genIso2", iso.getArrayPtr(2));
        tr->SetBranchAddress("genIso3", iso.getArrayPtr(3));
        tr->SetBranchAddress("genIso4", iso.getArrayPtr(4));
        tr->SetBranchAddress("genIso5", iso.getArrayPtr(5));

        //number of hits for reals
        tr->SetBranchAddress("genIhit1", isoh.getArrayPtr(1));
        tr->SetBranchAddress("genIhit2", isoh.getArrayPtr(2));
        tr->SetBranchAddress("genIhit3", isoh.getArrayPtr(3));
        tr->SetBranchAddress("genIhit4", isoh.getArrayPtr(4));
        tr->SetBranchAddress("genIhit5", isoh.getArrayPtr(5));
    }


    tr->SetBranchAddress("genE", genE);


    //iso due to reconstructed
    //tr->SetBranchAddress("muIso1",      riso1);


    tr->SetBranchAddress("nLeptons", &nlep);
    tr->SetBranchAddress("lepE", lepE);

    /////////////////////
    ///////////////
    // Input->GetParameters("aroundwego", &aro);
    // grabdata(datacardfile);





    if(u > 0)aro = 1;
    else aro = 0;
    if(aro == 1) ofile = "all_" + ofile;




    sprintf(hfName, "%d%d%d_%d%d%d", parameter1[1], parameter2[1], parameter3[1], nparameter1[1], nparameter2[1], nparameter3[1]);

    if(aro == 1) comp2 = string(hfName);
    comp.push_back(comp2);

    string parx("null");


    // Read the ntuple tree 
}

HcalAna::~HcalAna()
{

    delete Input;

}

float eta(float x, float y, float z)
{
    float theta, feta, tht; //tht- tan of half of theta feta is final eta
    theta = atan(fabs(sqrt(y * y + x * x) / z));
    tht = tan(theta / 2);
    feta = -log(tht);
    return feta;
}









// analysis template

void HcalAna::Analysis()
{


    if(parameter1[0] == 0)parameter1[1] = 0;
    if(parameter2[0] == 0)parameter2[1] = 0;
    if(parameter3[0] == 0)parameter3[1] = 0;
    if(nparameter1[0] == 0)nparameter1[1] = 0;
    if(nparameter2[0] == 0)nparameter2[1] = 0;
    if(nparameter3[0] == 0)nparameter3[1] = 0;


    sprintf(hfName, "%d%d%d_%d%d%d", parameter1[1], parameter2[1], parameter3[1], nparameter1[1], nparameter2[1], nparameter3[1]);
    // Prepare file and folder for histograms
    cout << "so here should have other shit " << ofile << endl;
    gSystem->mkdir(ofile.c_str());
    gSystem->mkdir((ofile + "/" + ifile).c_str());

    hfolder = "SLHC_" + string(hfName);

    gSystem->mkdir((ofile + "/" + ifile + "/" + hfolder).c_str());
    // create histogram files 
    TString Path_fName = ofile + "/" + ifile + "/" + hfolder + "/" + hfName + ".root"; //probably dont need to worry about 
    cout << Path_fName << endl << endl;
    theFile = new TFile(Path_fName, "RECREATE");
    theFile->cd();
    jets_directory = theFile->mkdir("Jets", "Jets");
    w_directory = theFile->mkdir("wbos", "wbos");
    t_directory = theFile->mkdir("gtot", "gtot");
    rec_directory = theFile->mkdir("reconstructed", "reconstructed");
    compare_directory = theFile->mkdir("compare", "compare");

    // booking histograms 
    //histos["h_nMu"] = new TH1D("h_nMu",   " number of reco Muons ", 9, 0., 9. ) ;
    //t_directory->cd();

    g_realte = new TH1D("g_realte", "the eta of the generated particels of the total", 50, -.8, .8); // eta of total 
    h_lepPt = new TH1D("h_lepPt", " muon Pt ", 50, 30, 400);
    h_nJets = new TH1D("h_nJets", "should be total number of generated oer event", 15, 0, 15); // again added by me should make new histograms
    //h_Jets    = new TH1D("h_Jets ",   " highest energy  of a jet per event", 50, 0, 5000 );
    g_isolt1 = new TH1D("g_isolt1", "the look of the iso due to muons of the first layer", 50, 0, 300); //jet
    g_isolt12 = new TH1D("g_isolt12", "the look of the iso due to muons of the first and second layer", 50, 0, 300);
    g_isolt123 = new TH1D("g_isolt123", "the look of the iso due to muons of all three layers", 50, 0, 300);
    g_isolth1 = new TH1D("g_isolth1", "the number of hits in the isolated regeon around the muon on the first layer", 50, 0, 50); //w decay hits
    g_isolth12 = new TH1D("g_isolth12", "he number of hits in the isolated regeon around the muon on the first and second layers", 50, 0, 50);
    g_isolth123 = new TH1D("g_isolth123", "the number of hits in the isolated regeon around the muon on the three layers", 50, 0, 50);

    g_realtp = new TH1D("g_realtp", "the pt of the  real of the total", 50, 0, 500); //momentum of total



    jets_directory->cd();
    h_gjets = new TH1D("gjets", "number of muons created by jets", 5, 0, 5);

    g_isolj1 = new TH1D("g_isolj1", "the energy deposted  on the first layer for the isolated muons", 50, 0, 300); //deposted stuff isolated g mean denerated
    g_isolj12 = new TH1D("g_isolj12", "the energy deposted  on the first and second layer for the isolated muons", 50, 0, 300);
    g_isolj123 = new TH1D("g_isolj123", "the energy deposted  on the all three layers for the isolated muons", 50, 0, 100);

    h_JetsPt = new TH1D("h_JetsPt ", " highest momentum  of a jet per event", 50, 30, 400);
    g_realje = new TH1D("g_realje", "the eta of the muon due to a jet", 50, -.8, .8); //eta of jets

    g_realjp = new TH1D("g_realjp", "the pt of the  of the jet", 50, 0, 500); //momentum of jets
    g_isoljh1 = new TH1D("g_isoljh1", "the number of hits in the isolated regeon around the muon on the first layer due to jets", 60, 0, 60); //number of hits
    g_isoljh12 = new TH1D("g_isoljh12", "the number of hits in the isolated regeon around the muon on the first and second layers due to jets", 60, 0, 60);
    g_isoljh123 = new TH1D("g_isoljh123", "the number of hits in the isolated regeon around the muon on the three layers due to jets", 60, 0, 60);

    rec_directory->cd();

    h_nMu = new TH1D("h_nMu", " number of reco Muons ", 9, 0., 9.);

    r_isol1 = new TH1D("r_isol1", "the engery deposted iso of the reconstructed muon on the first layer", 50, 0, 300); //reconstructed iso
    r_isol12 = new TH1D("r_isol12", "the engery deposted iso of the reconstructed muon on the first and second layer", 50, 0, 300);
    r_isol123 = new TH1D("r_isol123", "the engery deposted iso of the reconstructed muon on the  all three layers", 50, 0, 300);
    r_isolh1 = new TH1D("r_isolh1", "the number of hits around the isolated muon on the first layer", 60, 0, 60); //reconstructed iso hits
    r_isolh12 = new TH1D("r_isolh12", "the number of hits around the isolated muon on the first and second layer", 60, 0, 60);
    r_isolh123 = new TH1D("r_isolh123", "the number of hits around the isolated muon on all three layers", 60, 0, 60);




    w_directory->cd();

    //w 
    g_realwp = new TH1D("g_realwp", " the pt of the  due to w", 50, 0, 500); //momentum of w decayed muons
    g_realwe = new TH1D("g_realwe", " the eta of the zero crash due to w", 50, -.8, .8); //eta of w boson
    g_isolw1 = new TH1D("g_isolw1", "the look of the iso muons energy depost due to w decay on the first layer", 50, 0, 100); //w decay iso deposit
    g_isolw12 = new TH1D("g_isolw12", "the look of the iso  muons energy deposit caused by w decay on the first and second layer", 50, 0, 300);
    g_isolw123 = new TH1D("g_isolw123", "the look of the iso muons energy depost from w decays on all three layers", 50, 0, 100);
    h_gwz = new TH1D("gwz", "number of muons created by w decay", 5, 0, 5); //number of generated by w decay 

    g_isolwh1 = new TH1D("g_isolwh1", "the number of hits in the isolated regeon around the muon on the first layer due to w", 50, 0, 50); //w decay hits
    g_isolwh12 = new TH1D("g_isolwh12", "he number of hits in the isolated regeon around the muon on the first and second layers due to w", 50, 0, 50);
    g_isolwh123 = new TH1D("g_isolwh123", "the number of hits in the isolated regeon around the muon on the three layers due to w", 50, 0, 50);

    //hey don't forget you need this if you are using multiple efficiancies at the same time
    compare_directory -> cd();
    for(int i = 0; i < 6; i++)
    {

        char mesh1_name[5];
        char mesh2_name[5];
        char meshh1_name[5];
        char meshh2_name[5];



        sprintf(mesh1_name, "mesh1%d", i);
        sprintf(mesh2_name, "mesh2%d", i);
        sprintf(meshh1_name, "meshh1%d", i);
        sprintf(meshh2_name, "meshh2%d", i);


        if(parameter1[3] == 1)
        {
            mesh1.push_back(new TH1D(mesh1_name, "so this is whatever you chose to add together deposit energy", 400, 0, 200));
        }

        if(sparameter1[3] == 1)
        {
            mesh2.push_back(new TH1D(mesh2_name, "so this is whatever you chose to add together deposit energy", 400, 0, 200));
        }


        meshh1.push_back(new TH1D(meshh1_name, "so this is whate.00001ver you chose to add together deposit hits", 50, 0, 50));
        meshh2.push_back(new TH1D(meshh2_name, "so this is whatever you chose to add together deposit hits", 50, 0, 50));

    }

    float pjets = 0;
    int count;
    int totalN = tr->GetEntries();
    cout << " Total Number of entries : " << totalN << endl;




    // for(int x=1; x<=5; x++)
    //{

    for(int i = 0; i < totalN; i++) //particular event 
    {
        for(int h = 0; h < 20; h++)
        {
            test[h] = 0;
        }



        pjets = 0;
        count = 0;
        //cout<<endl<<"leaf we are on"<<i<<endl;                  
        if(ProcessEvents > 0 && i >= (ProcessEvents)) break;
        // Get the entry from the ntuple  
        tr->GetEntry(i);



        //So here is supposed to look at the generated shit
        int evngwz = 0; //event Group from w
        int evngnz = 0; // event group not from w


        float mesht1 = 0; //so this gives us the total deposted enenergy by the particle
        float mesht2 = 0;
        float meshht1 = 0; //so this gives us the total deposted enenergy by the particle
        float meshht2 = 0;


        //cout<<" =============================== "<<endl ;
        for(int k = 0; k < nGen; k++)//DUMBASS THIS DOES REALality nGen
        {



            for(int p = 0; p < (isreco ? nlep : 1); p++)
            {
                if(!isreco) p = k;
                if(test[p] == 1) continue;
                mesht1 = 0; //so this gives us the total deposted enenergy by the particle
                mesht2 = 0;
                meshht1 = 0; //so this gives us the total hits
                meshht2 = 0;


                TLorentzVector gP4;
                gP4.SetPxPyPzE(bPx[k], bPy[k], bPz[k], bE[k]);
                TLorentzVector rP4;
                rP4.SetPxPyPzE(lepPx[p], lepPy[p], lepPz[p], lepE[p]);
                //regedr=sqrt((rP4.Eta()-gP4.Eta())*(rP4.Eta()-gP4.Eta())+(rP4.Theta()-gP4.Theta())*(rP4.Theta()-gP4.Theta()));

                // regedr=DeltaR(gP4,rP4);
                if(rec == 1) regedr = gP4.DeltaR(rP4);
                else regedr = 0;
                //printf("Gen (%d) Reco (%d) -> dR = %.5f ", k, p, regedr ) ;
                if(abs(bid[k]) == type && abs(mom[k]) == 24 && gP4.Pt() > mucuts[0]&& (isreco ? ((regedr) < .1) : true) && fabs(gP4.Eta()) < mucuts[2])//sees if it comes from a w boson
                {

                    /* if(gisoh1[k][0]+gisoh1[k][1]+gisoh1[k][2]==0)
                       {
                       g_crashwe->Fill(gP4.Eta());
                       g_crashwp->Fill(gP4.Pt());
                       }
                     */

                    g_realwe->Fill(gP4.Eta());
                    g_realwp->Fill(gP4.Pt());

                    massTableofDoom(1, p, mesht1, mesht2, meshht1, meshht2);
                    ////new
                }


                if(abs(bid[k]) == type && abs(mom[k]) != 24 && gP4.Pt() > mucuts[0]&&(isreco ? (regedr < .1) : true) && eta(bPx[k], bPy[k], bPz[k]) < mucuts[2])// sees if it is a muon form a jet
                {

                    /*if(gisoh1[k][0]+gisoh1[k][1]+gisoh1[k][2]==0)  delete me ones this is working

                      {
                      g_crashje->Fill(gP4.Eta());
                      g_crashjp->Fill(gP4.Pt());
                      }*/
                    g_realje->Fill(gP4.Eta());
                    g_realjp->Fill(gP4.Pt());

                    g_isolj1->Fill(iso(1, k, 0)); // energy
                    g_isolj12->Fill(iso(1, k, 0) + iso(1, k, 1));
                    g_isolj123->Fill(iso(1, k, 0) + iso(1, k, 1) + iso(1, k, 2));
                    g_isoljh1->Fill(isoh(1, k, 0)); //hits
                    g_isoljh12->Fill(isoh(1, k, 0) + isoh(1, k, 1));
                    g_isoljh123->Fill(isoh(1, k, 0) + isoh(1, k, 1) + isoh(1, k, 2));
                    evngnz++;


                    massTableofDoom(2, p, mesht1, mesht2, meshht1, meshht2);


                }



                if(abs(bid[k]) == type && gP4.Pt() > mucuts[0] && eta(bPx[k], bPy[k], bPz[k]) < mucuts[2]&&(isreco ? (regedr < .1) : true))
                {
                    //THIS IS TOTAL
                    //
                    massTableofDoom(3, p, mesht1, mesht2, meshht1, meshht2);

                    test[p] = 1;

                    if(abs(bid[k]) == type && gP4.Pt() > mucuts[0] && gP4.Eta() < mucuts[2]&&(isreco ? regedr < .001 : true))
                    {

                        if((parameter1[0] == 1 && abs(mom[k]) == 24) || (parameter1[0] == 2 && abs(mom[k]) != 24) || (parameter1[0] == 3))
                        {
                            if(parameter1[3] == 2)
                            {
                                mesh1.at(hisn)->Fill(mesht1 / gP4.Pt());
                            }
                            if(parameter1[3] == 1)
                            {
                                mesh1.at(hisn)->Fill(mesht1 + .00001);
                            }
                        }
                        if((sparameter1[0] == 1 && abs(mom[k]) == 24) || (sparameter1[0] == 2 && abs(mom[k]) != 24) || (sparameter1[0] == 3))
                        {
                            if(sparameter1[3] == 2)
                            {
                                mesh2.at(hisn)->Fill(mesht2 / gP4.Pt());
                            }
                            if(sparameter1[3] == 1)
                            {
                                mesh2.at(hisn)->Fill(mesht2 + .0000001);
                            }
                        }

                        if((parameterh1[0] == 1 && abs(mom[k]) == 24) || (parameterh1[0] == 2 && abs(mom[k]) != 24) || (parameterh1[0] == 3))
                        {
                            //this should be a graph of hits
                            if(parameterh1[3] == 2)
                            {
                                meshh1.at(hisn)->Fill(meshht1 / lepE[p]);
                            }
                            if(parameterh1[3] == 1)
                            {

                                meshh1.at(hisn)->Fill(meshht1);
                            }
                        }

                        if((sparameterh1[0] == 1 && abs(mom[k]) == 24) || (sparameterh1[0] == 2 && abs(mom[k]) != 24) || (sparameterh1[0] == 3))
                        {
                            if(sparameterh1[3] == 2)
                            {
                                meshh2.at(hisn)->Fill(meshht2 / lepE[p]);
                            }
                            if(sparameterh1[3] == 1)
                            { // if (meshht2<1) cout<<"huh this aint right"<<endl;
                                meshh2.at(hisn)->Fill(meshht2);

                            }

                        }


                    }
                }


            }

        }


        //so this should be used to count number of 

        /////////////////////////////// 
        for(int p = 0; p < (isreco ? nlep : 0); p++)
        {

            TLorentzVector rP4;
            rP4.SetPxPyPzE(lepPx[p], lepPy[p], lepPz[p], lepE[p]);
            TLorentzVector gP4;
            gP4.SetPxPyPzE(bPx[p], bPy[p], bPz[p], bE[p]);

            if(rP4.Pt() > mucuts[0] && eta(lepPx[p], lepPy[p], lepPz[p]) < mucuts[2])
            {
                massTableofDoom(4, p, mesht1, mesht2, meshht1, meshht2);
                if(sparameter1[0] == 4 && test[p] != 1)
                {
                    if(sparameter1[3] == 1)mesh2.at(hisn)->Fill(mesht2 + .0000001);
                    if(sparameter1[3] == 2)mesh2.at(hisn)->Fill(mesht2 / rP4.Pt() + .0000001);
                    meshh2.at(hisn)->Fill(meshht2);
                }

            }

        }
        /////////////////




        /////////////////////////
        //cout<<"number of baby "<<test<<endl;

        //int a=evngwz+evngnz;


        h_gwz->Fill(evngwz);
        h_gjets->Fill(evngnz);
        //h_nJets->Fill(a);//DELETE ME this is total number of real particles

        //cout<<"number of seen muons "<<count<<endl;
        h_nMu->Fill(count);





        mesht1 = 0;
        mesht2 = 0;


    } // end of event looping



    // record file in histogram format
    theFile->cd();
    HistoWrite("", theFile);
    // cout<<"so we are off by="<<offby<<endl<<"total number of reco muons="<<tnum<<endl<<"total number of gen="<<tgen<<endl;
    cout << " historams written ! " << endl;
    theFile->Close();

    mesh1.clear();
    mesh2.clear();
    meshh1.clear();
    meshh2.clear();

}




// **************************************** 
// *           Draw Histograms            *
// **************************************** 

void HcalAna::DrawHistogram()
{

    OpenHistogram();


    cout << "crash yet?" << endl;

    Draw_H1(mesh1, "first comparisons you wanted to make", mesh2, "second comparisons you wanted to make");
    //Draw_H1( mesh2, "second comparisons you wanted to make");

    Draw_H2(meshh1, "first hit comparisons you wanted to make ", meshh2, "second comparisons you wanted to make hits");
    // Draw_H2( meshh2, "second comparisons you wanted to make");


    theFile->Close();
}

void HcalAna::Draw_H1(vector <TH1D*> h1, string plotname, vector <TH1D*> h2 /* = NULL */, string splotname /*""*/)
{
    double beg = 0, end = 0;
    cout << "Draw_H1 start\n";

    TGraph * effi[2];
    for(int k = 0; k < 2; k++)
    {
        cout << "Make " << k << " plot\n";
        effi[k] = new TGraph(10);
        hfolder = ofile + "/" + ifile + "/SLHC_" + comp.at(k) + "/";
        TCanvas* c1 = new TCanvas("c1", "", 800, 700);
        c1->SetFillColor(10);
        c1->SetFillColor(10);
        ////c1->SetLogy();

        gStyle->SetOptStat("eoumi");
        gStyle->SetOptStat("");

        h1.at(k)->SetLineColor(2);

        c1->SetLogy();
        string typ1, typ2, dr1(1, 0), dr2(1, 0), lay1(1, 0), lay2(1, 0), lay3(1, 0), slay1(1, 0), slay2(1, 0), slay3(1, 0);
        if(parameter1[0] == 1) typ1 = "wdecay";
        if(parameter1[0] == 2) typ1 = "jet";
        if(parameter1[0] == 3) typ1 = "total";
        if(parameter1[0] == 4) typ1 = "rec";
        if(sparameter1[0] == 1) typ2 = "wdecay";
        if(sparameter1[0] == 2) typ2 = "jet";
        if(sparameter1[0] == 3) typ2 = "total";
        if(sparameter1[0] == 4) typ2 = "rec";


        if(parameter1[0] != 0)sprintf(&dr1.at(0), "%d", parameter1[1]);
        if(sparameter1[0] != 0)sprintf(&dr2.at(0), "%d", sparameter1[1]);
        //int x=1+parameter1[2];
        if(parameter1[0] != 0)sprintf(&lay1.at(0), "%d", (1 + parameter1[2]));
        if(parameter2[0] != 0)sprintf(&lay2.at(0), "%d", (1 + parameter2[2]));
        if(parameter3[0] != 0)sprintf(&lay3.at(0), "%d", (1 + parameter3[2]));

        if(sparameter1[0] != 0)sprintf(&slay1.at(0), "%d", (1 + sparameter1[2]));
        if(sparameter2[0] != 0)sprintf(&slay2.at(0), "%d", (1 + sparameter2[2]));
        if(sparameter3[0] != 0)sprintf(&slay3.at(0), "%d", (1 + sparameter3[2]));




        //  dr1=parameter1[1];

        if(h2.at(k) != NULL)
        {


            float max1 = h1.at(k)->GetAt(h1.at(k)->GetMaximumBin());

            float max2 = h2.at(k)->GetAt(h2.at(k)->GetMaximumBin());

            h1.at(k)->SetMaximum(1.1 * max(max1, max2));

            h1.at(k)->Draw();

            h2.at(k)->SetLineColor(8);

            h2.at(k)->Draw("same");

            TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9); // lets make a legend to tell everything apart

            leg->SetHeader("energy deposted");

            leg-> SetFillColor(kWhite);

            leg-> SetShadowColor(kWhite);

            leg->AddEntry(h1.at(k), (typ1 + " dr of " + comp.at(k)).c_str(), "l");

            leg->AddEntry(h2.at(k), (typ2 + " dr of " + comp.at(k)).c_str(), "l");

            leg->Draw();



            double maxb[20], mins[20], effb[20], blocks[20];
            maxb[0] = 1;
            mins[0] = .8;
            double maxb2, mins2;
            maxb2 = maxb[0];
            mins2 = mins[0];

            effb[0] = 1;
            blocks[0] = 1; //slight misnomer this is how much of the noise gets through

            //TGraph* effi = new TGraph(10,mins,blocks);


            for(int i = 0; i < 10; i++)
            {


                efficiency(h1.at(k), h2.at(k), maxb[i], mins[i], effb[i], blocks[i]);
                //( TH1D* h1 , string plotname, TH1D* h2 , string splotname, double maxb, double mins, double *effb, double *blocks )
                mins2 = .02 + mins2;
                maxb2 = maxb2 - .02;

                effi[k]->SetPoint(i, mins[i], blocks[i]);

                mins[i + 1] = mins2;
                maxb[i + 1] = maxb2;

            }
            beg = mins[0];
            end = mins[9];


            // first attempt at getting this bloody thing to graph
            TCanvas *c2 = new TCanvas("c2", "sooooooooo working? no?", 200, 10, 700, 500);

            //Int_t n = 20;


            effi[k]->SetName("efficiency");
            effi[k]->SetTitle("efficiency");
            effi[k]->GetXaxis()->SetTitle("efficiency");
            effi[k]->GetYaxis()->SetTitle("Noise");



            effi[k]->Draw("AC*");
            c2->Update();

            TString plotname2 = hfolder + plotname + "effic" + "." + plotType;
            cout << "plotname2:" << plotname2 << endl;
            c2->Print(plotname2);
            delete c2;




        } else h1.at(k)->Draw();

        c1->SetLogy();

        c1->Update();




        TString plotname1 = hfolder + plotname + splotname + "." + plotType;
        c1->Print(plotname1);

        delete c1;

    }

    effi[0]->SetMarkerColor(2);
    effi[1]->SetMarkerColor(3);
    gStyle->SetOptFit();
    TCanvas *c3 = new TCanvas("c3", "multigraph", 700, 500);

    TMultiGraph *mg = new TMultiGraph("compo", "combo");

    c3->SetFillColor(10);
    c3->SetFillColor(10);

    mg->Add(effi[0]);
    mg->Add(effi[1]);

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("efficiency");
    mg->GetYaxis()->SetTitle("backround servival");

    // Change the axis limits
    //c3->BuildLegend();
    ////////
    TLegend* leg = new TLegend(0.15, .7, 0.35, .9); // lets make a legend to tell everything apart
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    //leg->AddEntry(h1.at(0),(typ1 + " dr of " + dr1 + " layers " + lay1 + lay2 + lay3 ).c_str(), "l");
    leg->AddEntry(effi[0], (comp1.c_str()), "p");
    leg->AddEntry(effi[1], (comp2.c_str()), "p");
    leg->Draw();

    c3->Update();

    gSystem->mkdir((ofile + "/" + ifile + "/compare"+comp1).c_str());

    TString plotname3 = ofile + "/" + ifile + "/compare"+comp1+"/" + comp1 + sub + "_" + comp2 + sub + "." + plotType;
    c3->Print(plotname3);
    /////////
    string helpplease = ofile + "/" + ifile + "/compare"+comp1+"/" + comp1 + sub + "_t" + comp2 + sub + ".root";
    TFile * gt = new TFile(helpplease.c_str(), "RECREATE");
    gt->cd();
    mg->Write();

    delete c3;

    TCanvas *c4 = new TCanvas("c4", "effratio", 700, 500);
    TF1 * fit0 = new TF1("expfit0", "[0]+exp([1]+x*[2])");
    effi[0]->Fit(fit0, "NQ");
    TF1 * fit1 = new TF1("expfit1", "[0]+exp([1]+x*[2])");
    effi[1]->Fit(fit1, "NQ");

    TF1 * ratio = new TF1("ratio", "([0]+exp([1]+x*[2]))/([3]+exp([4]+x*[5]))", beg, end);
    for(int i = 0; i < 3; i++)
        ratio->SetParameter(i, fit0->GetParameter(i));
    for(int i = 0; i < 3; i++)
        ratio->SetParameter(i + 3, fit1->GetParameter(i));

    gt->cd();
    string tit = comp1 + "/" + comp2 + "_ratio";
    ratio->SetTitle(tit.c_str());
    ratio->GetXaxis()->SetTitle("efficiency");
    ratio->Write();
    c4->cd();
    ratio->Draw();
    if(aro == 1) gSystem->mkdir((ofile + "/" + ifile + "/better" + comp1).c_str());
    if(aro != 1)
    {
        TString plotname4 = ofile + "/" + ifile + "/compare"+comp1+"/ratio" + comp1 + sub + "_" + comp2 + sub + "." + plotType;

        c4->Print(plotname4);
    }
    if((ratio->Integral(beg, end)>(end - beg)) && aro == 1)
    {

        TString plotname4 = ofile + "/" + ifile + "/better" + comp1 + "/" + comp1 + sub + "_" + comp2 + sub + "." + plotType;
        c4->Print(plotname4);

    }


    delete c4;

}

void HcalAna::Draw_H2(vector <TH1D*> h1, string plotname, vector <TH1D*> h2 /* = NULL */, string splotname /*""*/)
{
    double beg = .8, end = 1;


    TGraph * effi[2];
    for(int k = 0; k < 2; k++)
    {
        effi[k] = new TGraph(10);

        hfolder = ofile + "/" + ifile + "/SLHC_" + comp.at(k) + "/";

        TCanvas* c1 = new TCanvas("c1", "", 800, 700);
        c1->SetFillColor(10);
        c1->SetFillColor(10);


        gStyle->SetOptStat("eoumi");
        gStyle->SetOptStat("");

        //h1->GetXaxis()->SetTitle( "x title" );
        //h1->GetYaxis()->SetTitle( "y title" );
        h1.at(k)->SetLineColor(2);

        c1->cd();
        c1->SetLogy();
        string typ1, typ2, dr1(1, 0), dr2(1, 0), lay1(1, 0), lay2(1, 0), lay3(1, 0), slay1(1, 0), slay2(1, 0), slay3(1, 0);
        if(parameterh1[0] == 1) typ1 = "wdecay";
        if(parameterh1[0] == 2) typ1 = "jet";
        if(parameterh1[0] == 3) typ1 = "total";
        if(parameterh1[0] == 4) typ1 = "rec";
        if(sparameterh1[0] == 1) typ2 = "wdecay";
        if(sparameterh1[0] == 2) typ2 = "jet";
        if(sparameterh1[0] == 3) typ2 = "total";
        if(sparameterh1[0] == 4) typ2 = "rec";


        if(parameterh1[0] != 0)sprintf(&dr1.at(0), "%d", parameterh1[1]);
        if(sparameterh1[0] != 0)sprintf(&dr2.at(0), "%d", sparameterh1[1]);
        //int x=1+parameter1[2];
        if(parameterh1[0] != 0)sprintf(&lay1.at(0), "%d", (1 + parameterh1[2]));
        if(parameterh2[0] != 0)sprintf(&lay2.at(0), "%d", (1 + parameterh2[2]));
        if(parameterh3[0] != 0)sprintf(&lay3.at(0), "%d", (1 + parameterh3[2]));

        if(sparameterh1[0] != 0)sprintf(&slay1.at(0), "%d", (1 + sparameterh1[2]));
        if(sparameterh2[0] != 0)sprintf(&slay2.at(0), "%d", (1 + sparameterh2[2]));
        if(sparameterh3[0] != 0)sprintf(&slay3.at(0), "%d", (1 + sparameterh3[2]));




        //  dr1=parameter1[1];

        if(h2.at(0) != NULL)
        {


            float max1 = h1.at(k)->GetAt(h1.at(k)->GetMaximumBin());
            float max2 = h2.at(k)->GetAt(h2.at(k)->GetMaximumBin());

            h1.at(k)->SetMaximum(1.1 * max(max1, max2));

            h1.at(k)->Draw();
            h2.at(k)->SetLineColor(8);
            h2.at(k)->Draw("same");
            TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9); // lets make a legend to tell everything apart
            leg->SetHeader("combo of hits");
            leg-> SetFillColor(kWhite);
            leg-> SetShadowColor(kWhite);
            //leg->AddEntry(h1.at(0),(typ1 calElectron/+ " dr of " + dr1 + " layers " + lay1 + lay2 + lay3 ).c_str(), "l");
            leg->AddEntry(h1.at(k), (typ1 + " dr of " + comp.at(k)).c_str(), "l");
            leg->AddEntry(h2.at(k), (typ2 + " dr of " + comp.at(k)).c_str(), "l");
            leg->Draw();

            //HEY HEYHEY this should be making eff graphs



            double maxb[20], mins[20], effb[20], blocks[20];
            maxb[0] = 1;
            mins[0] = .8;
            double maxb2, mins2;
            maxb2 = maxb[0];
            mins2 = mins[0];

            effb[0] = 1;
            blocks[0] = 1; //slight misnomer this is how much of the noise gets through

            //TGraph* effi = new TGraph(10,mins,blocks);


            for(int i = 0; i < 10; i++)
            {


                efficiency(h1.at(k), h2.at(k), maxb[i], mins[i], effb[i], blocks[i]);
                // why is this here? ( TH1D* h1 , string plotname, TH1D* h2 , string splotname, double maxb, double mins, double *effb, double *blocks )
                mins2 = .02 + mins2;
                maxb2 = maxb2 - .02;
                effi[k]->SetPoint(i, mins[i], blocks[i]);

                mins[i + 1] = mins2;
                maxb[i + 1] = maxb2;

            }
            beg = mins[0];
            end = mins[10];
            // first attempt at getting this bloody thing to graph
            TCanvas *c2 = new TCanvas("c2", "sooooooooo working? hits", 200, 10, 700, 500);


            effi[k]->SetTitle("efficiencyh");
            effi[k]->GetXaxis()->SetTitle("efficiency");
            effi[k]->GetYaxis()->SetTitle("Noise");
            effi[k]->Draw("AC*");

            c2->Update();

            TString plotname2 = hfolder + plotname + "effich" + "." + plotType;
            c2->Print(plotname2);
            delete c2;
            // yep all this is trying to get this to graph
        } else
            h1.at(k)->Draw();
        c1->SetLogy();

        c1->Update();
        /*
           if(h2.at(0)!=NULL)
           {//what the hell is this and what does it do? 
        //may just be example

        TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);// lets make a legend to tell everything apart
        leg->SetHeader("The Legend Title");
        leg->AddEntry(h1,"first entry", "f");
        leg->AddEntry(h2,"second entry", "f");
        leg->Draw();
        } */

        //keep this guy



        string hhfName = comp.at(k);
        TString Path_fName = hfolder + hhfName + ".root";


        string plotname1 = hfolder + plotname + "h" + splotname + "." + plotType;
        c1->Print(plotname1.c_str());
        delete c1;

    }
    /////////
    effi[0]->SetMarkerColor(2);
    effi[1]->SetMarkerColor(3);
    gStyle->SetOptFit();
    TCanvas *c3 = new TCanvas("c3", "multigraph", 700, 500);



    TMultiGraph *mg = new TMultiGraph("comboH", "comboH");

    c3->SetFillColor(10);
    c3->SetFillColor(10);
    mg->Add(effi[0]);
    mg->Add(effi[1]);

    mg->Draw("AP"); //apl


    mg->GetXaxis()->SetTitle("efficiency");
    mg->GetYaxis()->SetTitle("backround servival");

    // Change the axis limits
    //c3->BuildLegend();
    ////////
    TLegend* leg = new TLegend(0.15, .7, 0.35, .9); // lets make a legend to tell everything apart
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    //leg->AddEntry(h1.at(0),(typ1 + " dr of " + dr1 + " layers " + lay1 + lay2 + lay3 ).c_str(), "l");
    leg->AddEntry(effi[0], (comp1.c_str()), "p");
    leg->AddEntry(effi[1], (comp2.c_str()), "p");
    leg->Draw();

    c3->Update();

    gSystem->mkdir((ofile + "/" + ifile + "/compare"+comp1+"/").c_str());

    TString plotname3 = ofile + "/" + ifile + "/compare"+comp1+"/" + comp1 + sub + "h_" + comp2 + sub + "." + plotType;
    c3->Print(plotname3);
    string helpplease = ofile + "/" + ifile + "/compare"+comp1+"/" + comp1 + sub + "h_t" + comp2 + sub + ".root";
    TFile * gt = new TFile(helpplease.c_str(), "RECREATE");
    gt->cd();
    mg->Write();

    delete c3;


    ////////
    TCanvas * c4 = new TCanvas("c4h", "effratioh", 700, 500);
    TF1 * fit0 = new TF1("expfit0", "[0]+exp([1]+x*[2])");
    effi[0]->Fit(fit0, "NQ");
    TF1 * fit1 = new TF1("expfit1", "[0]+exp([1]+x*[2])");
    effi[1]->Fit(fit1, "NQ");

    TF1 * ratio = new TF1("ratioh", "([0]+exp([1]+x*[2]))/([3]+exp([4]+x*[5]))", beg, end);

    for(int i = 0; i < 3; i++)
        ratio->SetParameter(i, fit0->GetParameter(i));
    for(int i = 0; i < 3; i++)
        ratio->SetParameter(i + 3, fit1->GetParameter(i));

    gt->cd();
    string tit = comp1 + sub + "/" + comp2 + sub + "_ratioh";
    ratio->SetTitle(tit.c_str());
    //ratio->SetTitle("WHAT THE FUCK");
    ratio->Write();
    ratio->GetXaxis()->SetTitle("efficiency");
    c4->cd();
    ratio->Draw();
    //TString plotname4 = ofile+"/"+ifile+"/compare/ratioh"+comp1+sub+"_"+comp2+sub+"."+plotType ;
    //c4->Print( plotname4 );

    if(aro == 1) gSystem->mkdir((ofile + "/" + ifile + "/better" + comp1).c_str());
    if(aro != 1)
    {
        TString plotname4 = ofile + "/" + ifile + "/compare"+comp1+"/ratioh" + comp1 + sub + "_" + comp2 + sub + "." + plotType;

        c4->Print(plotname4);
    }
    cout << "aro:" << aro << endl;
    cout<<"ratio:"<<ratio->Integral(beg, end)<<"  (end - beg)"<<(end - beg)<<endl;
    if((ratio->Integral(beg, end)>(end - beg)) && aro == 1)
    {

        TString plotname4 = ofile + "/" + ifile + "/better" + comp1 + "/h" + comp1 + sub + "_" + comp2 + sub + "." + plotType;
        c4->Print(plotname4);

    }

    delete c4;
    /////////////
}

void HcalAna::efficiency(TH1D* h1, TH1D* h2, double &maxb, double &mins, double &effb, double &blocks)
{
    double H1t = h1->Integral(); // total area under the curve for the h1
    double H2t = h2->Integral(); //total area under the curve for the h2
    double H2 = h2->Integral(1, 1); //this should be the area under the curve for h2 to whatever bin we are up to
    double H1 = h1->Integral(1, 1); // H1 to h1 is same as h2 to H2
    double H1p;
    double H2p;
    //int nbins= h1->GcoutetNbinsX(); // this will give us the number of bins if I ever need it.
    int binb = 1, bins = 1; //binb should be the bin that gives us the highest bin we go to 
    int i = 0;
    H2p = H2 / H2t;
    H1p = H1 / H1t;
    do
    {
        i++;
        H2 = h2->Integral(1, i);
        H1 = h1->Integral(1, i);

        H2p = H2 / H2t;
        H1p = H1 / H1t;

        if(H2p < maxb) binb = i; //so this is probably the problem
        if(H1p < mins) bins = i + 1;
        //delete these down to while
        H2 = h2->Integral(1, bins);
        H1 = h1->Integral(1, bins);
        blocks = H2 / H2t;
        H1p = H1 / H1t;
        //cout<<"so the of signal is "<<H1<<" backround is "<<H2<<endl;
    } while(H1p < mins);
    //cout<<H2p<<" so that was H2p and H1p is"<< H1p<<endl;
    H2 = h2->Integral(1, bins);
    H1 = h1->Integral(1, binb);
    blocks = H2 / H2t; //slight misnomer this is how much of the noise gets through
    //cout<<H1<<"so this is H1 "<< binb << "and this is the bin causing the issues"<<endl; 
    effb = H1 / H1t;
    maxb = (h2->Integral(1, binb)) / H2t;
    mins = (h1->Integral(1, bins)) / H1t;
    //cout<< "sooooooooooooo lets hope we get nice numbers blocks "<< blocks<<" and of course eff "<< effb<<endl;
}

void HcalAna::OpenHistogram()
{//not sure look into this

    hfolder = ofile + "/" + ifile + "/SLHC_" + comp.at(0) + "/";
    string h1fName = string(comp.at(0));
    TString Path_fName = hfolder + comp.at(0) + ".root";
    cout << " Opening : " << Path_fName << endl;

    theFile = (TFile*) TFile::Open(Path_fName, "READ");
    //hFile->cd() ;
    cout << " file opened ! " << endl;



    TString c_dir = "compare/";

    char mesh1_name[5];
    char mesh2_name[5];
    char meshh1_name[5];
    char meshh2_name[5];



    sprintf(mesh1_name, "mesh1%d", hisn);
    sprintf(mesh2_name, "mesh2%d", hisn);
    sprintf(meshh1_name, "meshh1%d", hisn);
    sprintf(meshh2_name, "meshh2%d", hisn);


    //sprintf(hfName,"%d%d%d_%d%d%d" ,parameter1[1],parameter2[1],parameter3[1],nparameter1[1],nparameter2[1],nparameter3[1]);


    mesh1.push_back((TH1D*) theFile->Get(c_dir + mesh1_name));
    mesh2.push_back((TH1D*) theFile->Get(c_dir + mesh2_name));

    meshh1.push_back((TH1D*) theFile->Get(c_dir + meshh1_name));
    meshh2.push_back((TH1D*) theFile->Get(c_dir + meshh2_name));
    //////////////////

    hfolder = ofile + "/" + ifile + "/SLHC_" + comp.at(1) + "/";
    cout << "so this is the folder " << hfolder << endl;
    //hfName=comp.at(1)
    TString Path_fName1 = hfolder + comp.at(1) + ".root";
    cout << " Opening : " << Path_fName1 << endl;

    theFile = (TFile*) TFile::Open(Path_fName1, "READ");
    //hFile->cd() ;
    cout << " file opened ! " << endl;





    sprintf(mesh1_name, "mesh1%d", hisn);
    sprintf(mesh2_name, "mesh2%d", hisn);
    sprintf(meshh1_name, "meshh1%d", hisn);
    sprintf(meshh2_name, "meshh2%d", hisn);


    //sprintf(hfName,"%d%d%d_%d%d%d" ,parameter1[1],parameter2[1],parameter3[1],nparameter1[1],nparameter2[1],nparameter3[1]);


    // create histogram files 

    mesh1.push_back((TH1D*) theFile->Get(c_dir + mesh1_name));
    mesh2.push_back((TH1D*) theFile->Get(c_dir + mesh2_name));

    meshh1.push_back((TH1D*) theFile->Get(c_dir + meshh1_name));
    meshh2.push_back((TH1D*) theFile->Get(c_dir + meshh2_name));

    cout << " Done opening. " << endl;
}

void HcalAna::HistoWrite(string theFolder, TFile * file)
{

    if(theFolder.size() > 0) file->cd(theFolder.c_str()); //actually saving files


    //g_realte->t_direcory->Write();






    cout << theFolder.c_str() << endl << endl;

    /*
       jets_directory->cd();
       g_isolj1->Write();
       g_isolj12->Write();
       g_isolj123->Write();
       g_isoljh1->Write();
       g_isoljh12->Write();
       g_isoljh123->Write();
       g_realje->Write();
       g_realjp->Write();
       h_gjets->Write();//number of muons created by jets


       w_directory->cd();
       g_isolw1->Write();
       g_isolw12->Write();
       g_isolw123->Write();
       g_isolwh1->Write();
       g_isolwh12->Write();
       g_isolwh123->Write();
       h_gwz->Write();
       g_realwe->Write();
       g_realwp->Write();

       t_directory->cd();
       g_isolt1->Write();
       g_isolt12->Write();
       g_isolt123->Write();
       g_isolth1->Write();
       g_isolth12->Write();
       g_isolth123->Write();
       g_realtp->Write();
       g_realte->Write();
       h_lepPt->Write()  ;
       h_nJets->Write() ;//note this is thdr1.at(i)e jets pt not the muons

       rec_directory->cd();
       h_nMu->Write()   ;//number of muons recon
       r_isol1->Write();
       r_isol12->Write();
       r_isol123->Write();
       r_isolh1->Write();
       r_isolh12->Write();
       r_isolh123->Write();
     */
    compare_directory->cd();



    mesh1.at(hisn)->Write();
    mesh2.at(hisn)->Write();
    meshh1.at(hisn)->Write();
    meshh2.at(hisn)->Write();


}

void HcalAna::massTableofDoom(int par, int k, float& mesht1, float& mesht2, float& meshht1, float& meshht2)
{


    if(parameter1[0] == par)
    {
        mesht1 = (iso(parameter1[1], k, parameter1[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) mesht1 -= gen_hcalE[k][0];
    }
    if(parameter2[0] == par)
    {
        mesht1 += (iso(parameter2[1], k, parameter2[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) mesht1 -= gen_hcalE[k][1];
    }
    if(parameter3[0] == par)
    {
        mesht1 += (iso(parameter3[1], k, parameter3[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) mesht1 -= gen_hcalE[k][2];
    }
    //subtracting
    if(nparameter1[0] == par)
    {
        mesht1 -= (iso(nparameter1[1], k, nparameter1[2]));


    }
    if(nparameter2[0] == par)
    {
        mesht1 -= (iso(nparameter2[1], k, nparameter2[2]));
    }
    if(nparameter3[0] == par)
    {
        mesht1 -= (iso(nparameter3[1], k, nparameter3[2]));

    }



    ////////////
    if(sparameter1[0] == par)
    {
        mesht2 = (iso(sparameter1[1], k, sparameter1[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) mesht2 -= gen_hcalE[k][0];
    }
    if(sparameter2[0] == par)
    {
        mesht2 += (iso(sparameter2[1], k, sparameter2[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) mesht2 -= gen_hcalE[k][1];
    }
    if(sparameter3[0] == par)
    {
        mesht2 += (iso(sparameter3[1], k, sparameter3[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) mesht2 -= gen_hcalE[k][2];
    }
    //////////subtracting
    if(snparameter1[0] == par)
    {
        mesht2 -= (iso(snparameter1[1], k, snparameter1[2]));
    }
    if(snparameter2[0] == par)
    {
        mesht2 -= (iso(snparameter2[1], k, snparameter2[2]));
    }
    if(snparameter3[0] == par)
    {
        mesht2 -= (iso(snparameter3[1], k, snparameter3[2]));
    }
    //////
    //okay hits w boson
    if(parameterh1[0] == par)
    {
        meshht1 = (isoh(parameterh1[1], k, parameterh1[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) meshht1 -= 1;
    }
    if(parameterh2[0] == par)
    {
        meshht1 += (isoh(parameterh2[1], k, parameterh2[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) meshht1 -= 1;
    }
    if(parameterh3[0] == par)
    {
        meshht1 += (isoh(parameterh3[1], k, parameterh3[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) meshht1 -= 1;
    }
    //subtracting


    if(nparameterh1[0] == par)
    {
        meshht1 -= (isoh(nparameterh1[1], k, nparameterh1[2]));

    }
    if(nparameterh2[0] == par)
    {
        meshht1 -= (isoh(nparameterh2[1], k, nparameterh2[2]));
    }
    if(nparameterh3[0] == par)
    {
        meshht1 -= (isoh(nparameterh3[1], k, nparameterh3[2]));

    }
    ///////////
    if(sparameterh1[0] == par)
    {
        meshht2 = (isoh(sparameterh1[1], k, sparameterh1[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) meshht2 -= 1;
    }
    if(sparameterh2[0] == par)
    {
        meshht2 += (isoh(sparameterh2[1], k, sparameterh2[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) meshht2 -= 1;
    }
    if(sparameterh3[0] == par)
    {
        meshht2 += (isoh(sparameterh3[1], k, sparameterh3[2]));
        if(asub == 1 && gen_hcaldR[k][1] <= drval) meshht2 -= 1;
    }
    ////////////////subtracting
    if(snparameterh1[0] == par)
    {
        meshht2 -= (isoh(snparameterh1[1], k, snparameterh1[2]));
    }
    if(snparameterh2[0] == par)
    {
        meshht2 -= (isoh(snparameterh2[1], k, snparameterh2[2]));
    }
    if(snparameterh3[0] == par)
    {
        meshht2 -= (isoh(snparameterh3[1], k, snparameterh3[2]));
    }


    if(asub == 2 && parameter1[0] == par)
    {
        if(gen_hcaldR[k][0] < gen_hcaldR[k][1] && gen_hcaldR[k][0] < gen_hcaldR[k][2] && gen_hcaldR[k][0] < maxdr)mesht1 -= gen_hcalE[k][0];
        if(gen_hcaldR[k][1] < gen_hcaldR[k][0] && gen_hcaldR[k][1] < gen_hcaldR[k][2] && gen_hcaldR[k][1] < maxdr)mesht1 -= gen_hcalE[k][1];
        if(gen_hcaldR[k][2] < gen_hcaldR[k][1] && gen_hcaldR[k][2] < gen_hcaldR[k][0] && gen_hcaldR[k][2] < maxdr)mesht1 -= gen_hcalE[k][2];

        if(gen_hcaldR[k][0] < gen_hcaldR[k][1] && gen_hcaldR[k][0] == gen_hcaldR[k][2] && gen_hcalE[k][0] < gen_hcalE[k][2] && gen_hcaldR[k][0] < maxdr)mesht1 -= gen_hcalE[k][0];

        if(gen_hcaldR[k][1] < gen_hcaldR[k][0] && gen_hcaldR[k][1] == gen_hcaldR[k][2] && gen_hcalE[k][1] < gen_hcalE[k][2] && gen_hcaldR[k][1] < maxdr)mesht1 -= gen_hcalE[k][1];

        if(gen_hcaldR[k][1] < gen_hcaldR[k][2] && gen_hcaldR[k][1] == gen_hcaldR[k][0] && gen_hcalE[k][1] < gen_hcalE[k][0] && gen_hcaldR[k][1] < maxdr)mesht1 -= gen_hcalE[k][1];


        if(gen_hcaldR[k][0] < gen_hcaldR[k][2] && gen_hcaldR[k][0] == gen_hcaldR[k][1] && gen_hcalE[k][0] < gen_hcalE[k][1] && gen_hcaldR[k][0] < maxdr)mesht1 -= gen_hcalE[k][0];

        if(gen_hcaldR[k][2] < gen_hcaldR[k][0] && gen_hcaldR[k][1] == gen_hcaldR[k][2] && gen_hcalE[k][2] < gen_hcalE[k][1] && gen_hcaldR[k][2] < maxdr)mesht1 -= gen_hcalE[k][2];

        if(gen_hcaldR[k][2] < gen_hcaldR[k][1] && gen_hcaldR[k][2] == gen_hcaldR[k][0] && gen_hcalE[k][2] < gen_hcalE[k][0] && gen_hcaldR[k][2] < maxdr)mesht1 -= gen_hcalE[k][2];
    }


    if(asub == 3 && parameter1[0] == par)
    {
        if(gen_hcaldR[k][0] < gen_hcaldR[k][1] && gen_hcaldR[k][0] < gen_hcaldR[k][2] && gen_hcaldR[k][0] < maxdr)mesht1 -= gen_hcalE[k][0];
        if(gen_hcaldR[k][1] < gen_hcaldR[k][0] && gen_hcaldR[k][1] < gen_hcaldR[k][2] && gen_hcaldR[k][2] < maxdr)mesht1 -= gen_hcalE[k][1];
        if(gen_hcaldR[k][2] < gen_hcaldR[k][1] && gen_hcaldR[k][2] < gen_hcaldR[k][0] && gen_hcaldR[k][2] < maxdr)mesht1 -= gen_hcalE[k][2];

        if(gen_hcaldR[k][0] < gen_hcaldR[k][1] && gen_hcaldR[k][0] == gen_hcaldR[k][2] && gen_hcalE[k][0] > gen_hcalE[k][2] && gen_hcaldR[k][0] < maxdr)mesht1 -= gen_hcalE[k][0];

        if(gen_hcaldR[k][1] < gen_hcaldR[k][0] && gen_hcaldR[k][1] == gen_hcaldR[k][2] && gen_hcalE[k][1] > gen_hcalE[k][2] && gen_hcaldR[k][1] < maxdr)mesht1 -= gen_hcalE[k][1];

        if(gen_hcaldR[k][1] < gen_hcaldR[k][2] && gen_hcaldR[k][1] == gen_hcaldR[k][0] && gen_hcalE[k][1] > gen_hcalE[k][0] && gen_hcaldR[k][1] < maxdr)mesht1 -= gen_hcalE[k][1];


        if(gen_hcaldR[k][0] < gen_hcaldR[k][2] && gen_hcaldR[k][0] == gen_hcaldR[k][1] && gen_hcalE[k][0] > gen_hcalE[k][1] && gen_hcaldR[k][0] < maxdr)mesht1 -= gen_hcalE[k][0];

        if(gen_hcaldR[k][2] < gen_hcaldR[k][0] && gen_hcaldR[k][1] == gen_hcaldR[k][2] && gen_hcalE[k][2] > gen_hcalE[k][1] && gen_hcaldR[k][2] < maxdr)mesht1 -= gen_hcalE[k][2];

        if(gen_hcaldR[k][2] < gen_hcaldR[k][1] && gen_hcaldR[k][2] == gen_hcaldR[k][0] && gen_hcalE[k][2] > gen_hcalE[k][0] && gen_hcaldR[k][2] < maxdr)mesht1 -= gen_hcalE[k][2];
    }


    if((asub == 2 || asub == 3) && parameterh1[0] == par)// should be subtracting one for hits... I hope
    {
        if(gen_hcaldR[k][0] < maxdr || gen_hcaldR[k][1] < maxdr || gen_hcaldR[k][2] < maxdr)meshht1 -= 1;
    }

    if((asub == 2 || asub == 3) && sparameterh1[0] == par)// should be subtracting one for hits... I hope
    {

        if(gen_hcaldR[k][0] < maxdr || gen_hcaldR[k][1] < maxdr || gen_hcaldR[k][2] < maxdr)meshht2 -= 1;
    }

}

//void HcalAna::grabdata(string datacardfile)
