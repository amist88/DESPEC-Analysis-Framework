#include "CalibParameter.h"
#include <iostream>

using namespace std;
 
CalibParameter::CalibParameter()
: TGo4Parameter("AnaOnlineCalPar")
{ Int_t i;
    //plastic VME
    for (i = 0; i < 32; i++){
       AplasQDC[i] = 0;
       BplasQDC[i] = 0;
        } 
        for (i = 0; i < 32; i++){
          AplasTDC_Chref_dT[i] = 0;
          BplasTDC_SC41dT[i] = 0;
        }
     //FATIMA
    for (i = 0; i < 50; i++){
       DetIDFAT = 0; 
       Extra1fat[i] = 0;   
       Afat[i] = 0;
       Bfat[i] = 0;
       Cfat[i] = 0;
       Dfat[i] = 0;
       
       TFatTDC_Chref_dT[i] = 0;
       TFatTDC_SC41dT[i] = 0;
       TDCfatID = 0;
        
     //GALILEO
     for (i = 0; i < 36; i++){
       DetIDGal = 0;  
       AGal[i] = 0;
       BGal[i] = 0;
       CGal[i] = 0;
        } 
    } 
}
  //----------------------------------------------
CalibParameter::CalibParameter(const Text_t *name)
: TGo4Parameter(name)
{
  Int_t i;
  ifstream    file;
  
      file.open("Configuration_Files/Plastic_VME_Energy_Calibration.dat");
  if (file.fail()) {
    cout << "ERROR:  Could not open file: Plastic_VME_Calibration.dat ! (params set to 1 and 0)\n";
    for (i = 0; i < 32; i++){
       AplasQDC[i] = 1.;
       BplasQDC[i] = 0.;
      
    }
  }
  
   else {
    cout << "CalibParameter - reading calibration from: Plastic_VME_Energy_Calibration.dat\n";
    for (i = 0; i < 32; i++){
       if(IsData(file)) file >> AplasQDC[i]>> BplasQDC[i] ;
       if (file.fail()) cout << "ERROR reading Plastic_VME_Energy_Calibration.dat\n";
    }
  }
  file.close();
  //-----------------------------------------------------------------------------------------------------------
   file.open("Configuration_Files/Plastic_VME_Ref_Time_Calibration.dat");
  if (file.fail()) {
    cout << "ERROR:  Could not open file: Plastic_VME_Ref_Time_Calibration.dat ! (param set to 0)\n";
    for (i = 0; i < 32; i++){
       AplasTDC_Chref_dT[i] = 0;
       BplasTDC_SC41dT[i] = 0.;
    }
  }
  
   else {
    cout << "CalibParameter - reading calibration from: Plastic_VME_Ref_Time_Calibration.dat\n";
    for (i = 0; i < 32; i++){
      if(IsData(file)) file >> AplasTDC_Chref_dT[i] >> BplasTDC_SC41dT[i];
       if (file.fail()) cout << "ERROR reading Plastic_VME_Ref_Time_Calibration.dat\n";
    }
  }
  file.close();
  //-----------------------------------------------------------------------------------------------------//
   file.open("Configuration_Files/FATIMA_Energy_Calibration.txt");
  if (file.fail()) {
    cout << "ERROR:  Could not open file: FATIMA_Energy_Calibration.txt ! (params set to 0 1 and 0)\n";
    for (i = 0; i < 50; i++){
       DetIDFAT = i; 
       Extra1fat[i] = 0;    
       Afat[i] = 0;
       Bfat[i] = 0;
       Cfat[i] = 1;
       Dfat[i] = 0;
    }
  }
  
   else {
    cout << "CalibParameter - reading calibration from: FATIMA_Energy_Calibration.txt\n";
    for (i = 0; i < 50; i++){
      if(IsData(file)) file >> DetIDFAT>> Extra1fat[i]   >> Afat[i] >> Bfat[i] >>Cfat[i] >>Dfat[i] ;
       if (file.fail()) cout << "ERROR reading FATIMA_Energy_Calibration.txt\n ";
    }
  }
  file.close();
  //------------------------------------------------------------------------------//
  file.open("Configuration_Files/FATIMA_Ref_Time_Calibration.dat");
  if (file.fail()) {
    cout << "ERROR:  Could not open file: FATIMA_Ref_Time_Calibration.dat ! (param set 0)\n";
    for (i = 0; i < 50; i++){
       TDCfatID = i;
       TFatTDC_Chref_dT[i] = 0;
       TFatTDC_SC41dT[i] = 0;
    }
  }
  
   else {
    cout << "CalibParameter - reading calibration from: FATIMA_Ref_Time_Calibration.dat\n";
    for (i = 0; i < 50; i++){
      if(IsData(file)) file >> TDCfatID >> TFatTDC_Chref_dT[i]>>TFatTDC_SC41dT[i];
       if (file.fail()) cout << "ERROR reading FATIMA_Ref_Time_Calibration.dat\n";
    }
  }
  file.close();
  
  //------------------------------------------------------------------------------//
  file.open("Configuration_Files/GALILEO_Energy_Calibration.txt");
  if (file.fail()) {
    cout << "ERROR:  Could not open file: GALILEO_Energy_Calibration.txt ! (params set to 1 and 0)\n";
    for (i = 0; i < 36; i++){
       DetIDGal=i;
       AGal[i] = 0.;
       BGal[i] = 1.;
       CGal[i] = 0.;
    }
  }
  
   else {
    cout << "CalibParameter - reading calibration from: GALILEO_Energy_Calibration.txt\n";
    for (i = 0; i < 36; i++){
      if(IsData(file)) file >> DetIDGal >>AGal[i] >> BGal[i] >> CGal[i] ;
       if (file.fail()) cout << "ERROR reading GALILEO_Energy_Calibration.txt\n";
    }
  }
  file.close();
}
//------------------------------------------------------------------------------//

CalibParameter::~CalibParameter()
{}
//------------------------------------------------------------------------------

Int_t CalibParameter::PrintParameter(Text_t *buf, Int_t)
{
  //Int_t  i;
  //cout << "\n AnaOnline Calibration Parameter: " << GetName() << endl;

  //Plastic main
//   cout << "Plastic Energy calibration:\n";
//   for (i = 0; i < 32; i++){
//       cout << "strip " << i << "  \t A = " << AplasQDC[i] << "\t B = " << BplasQDC[i] << endl;
//     }
//     
//     cout << "Plastic Timing calibration:\n";
//   for (i = 0; i < 32; i++){
//       cout << "strip " << i << "  \t A = " << AplasTDC_Chref_dT[i] << endl;
//     }
//     
//     //FATIMA
//     cout << "FATIMA Energy calibration:\n";
//   for (i = 0; i < 50; i++){
//       cout << "detector " << i << "  \t A = " << Afat[i] << "\t B = " << Bfat[i] << "\t C = " << Cfat[i]<< "\t D = " << Dfat[i]<<  endl;
//     }
// //     
// //       cout << "FATIMA RefPlas-Fat calibration:\n";
// //   for (i = 0; i < 50; i++){
// //       cout << "strip " << i << "  \t A = " << Tfat[i] << endl;
//   
//         cout << "GALILEO Energy calibration:\n";
//         for (i = 0; i < 36; i++){
//         cout << "strip " << i << "  \t A = " << AGal[i] << "  \t B = " << BGal[i] << "  \t C = " << CGal[i]<< endl;
//     }
    return 1;
}

//------------------------------------------------------------------------------

Bool_t CalibParameter::UpdateFrom(TGo4Parameter *pp)
{
  if(pp->InheritsFrom("CalibParameter"))
  {
    Int_t i;
    CalibParameter *from = (CalibParameter *) pp;
    
    for (i = 0; i < 32; i++){
       AplasQDC[i]  = from->AplasQDC[i];
       BplasQDC[i]  = from->BplasQDC[i];

       cout << "CalibParameter - Parameter : " << GetName() << " UPDATED\n";
     }
     
      for (i = 0; i < 32; i++){
       AplasTDC_Chref_dT[i]  = from->AplasTDC_Chref_dT[i];
       BplasTDC_SC41dT[i] = from ->BplasTDC_SC41dT[i];
       cout << "CalibParameter - Parameter : " << GetName() << " UPDATED\n";
     }
     
      for (i = 0; i < 50; i++){
       DetIDFAT = from->DetIDFAT; 
       Extra1fat[i] = from->Extra1fat[i];
       
       Afat[i] = from->Afat[i];
       Bfat[i] = from->Bfat[i];
       Cfat[i] = from->Cfat[i];
       Dfat[i] = from->Dfat[i];
       TDCfatID = from-> TDCfatID;
       TFatTDC_Chref_dT[i]  = from->TFatTDC_Chref_dT[i];
       TFatTDC_SC41dT[i]  = from->TFatTDC_SC41dT[i];
       
      
      cout << "CalibParameter - Parameter : " << GetName() << " UPDATED\n";
    }
    
//     for (i = 0; i < 50; i++){
//        Tfat[i]  = from->Tfat[i];
//        cout << "CalibParameter - Parameter : " << GetName() << " UPDATED\n";
//      }
  
  
    for (i = 0; i < 36; i++){
       DetIDGal = from->DetIDGal;
       AGal[i]  = from->AGal[i];
       BGal[i]  = from->BGal[i];
       CGal[i]  = from->CGal[i];
       
       cout << "CalibParameter - Parameter : " << GetName() << " UPDATED\n";
     }
  }
     
      else {
      cout << "ERROR: CalibParameter - Wrong parameter object: " << pp->ClassName() << endl;}
 
        return kTRUE;
         
}
    int CalibParameter::IsData(ifstream &f) {
        char dum;
        char dumstr[300];
        int retval = 0;

        /* 'operator >>' does not read End-of-line, therefore check if read 
            character is not EOL (10) */
        do {
            dum = f.get();
            if (dum == '#')    // comment line => read whole line and throw it away
            f.getline(dumstr,300);
        }
        while ((dum == '#') || ((int)dum == 10)); 

        f.unget();   // go one character back
        retval = 1;
        return retval;
    }

ClassImp(CalibParameter)
