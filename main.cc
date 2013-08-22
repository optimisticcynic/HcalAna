#include <iostream> 
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMinuit.h>

#include "AnaInput.h"
#include "HcalAna.h"
#include "getopt.h"

using namespace std;

int main(int argc, char* argv[])
{
    bool runAll = false;
    string datacardfile = "DataCard.txt";
    int opt = 0;
    while((opt = getopt(argc, argv, "bd:")) != -1)
    {
        switch(opt)
        {
            case 'b':
                runAll = true;
                break;
            case 'd':
                datacardfile = optarg;
                break;
        }
    }

    //cout << "runALL: " << runAll << endl << "datacard: " << datacardfile << endl;


    //string datacardfile = (argc > 1) ? argv[1] : "DataCard.txt";

    AnaInput *Input = new AnaInput(datacardfile);

    int module = -1;
    Input->GetParameters("Module", &module);
    cout << "  Run module " << module << endl;

    if(runAll)
    {
        for(int u = 1; u <=5; u++)
        {
            for(int y = 0; y <=5; y++)
            {
                for(int t = 0; t <=5; t++)
                {


                    HcalAna *hAna = new HcalAna(datacardfile, u, y, t);

                    if((module == 1) || (module == 3)) hAna->Analysis();
                    if((module == 2) || (module == 3)) hAna->DrawHistogram(); //probably going to want to change this back to one

                    delete hAna;
                }
            }
        }
    } else
    {
        HcalAna *hAna = new HcalAna(datacardfile);
        if((module == 1) || (module == 3)) hAna->Analysis();
        if((module == 2) || (module == 3)) hAna->DrawHistogram(); //probably going to want to change this back to one

        delete hAna;
    }

    return 0;
}

