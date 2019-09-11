// $Id: EventAnlProc.cxx 754 2011-05-18 11:04:52Z adamczew $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum f�r Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
    //-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#include "EventAnlProc.h"

#include <cstdlib>
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"

#include "TGo4WinCond.h"
#include "TGo4Analysis.h"

#include "TGo4Picture.h"

#include "EventAnlStore.h"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "TSCNParameter.h"

#include "TAidaConfiguration.h"
//-----------------------------------------------------------
EventAnlProc::EventAnlProc() :
TGo4EventProcessor(),
fParam1(0)
{
}
//-----------------------------------------------------------
EventAnlProc::EventAnlProc(const char* name) :
TGo4EventProcessor(name)
{
  //Clear up for AIDA
  implantEvents = 0;
  decayEvents = 0;
  pulserEvents = 0;
  nonsenseEvents = 0;




  cout << "**** EventAnlProc: Create" << endl;

  TFile* f;
  f = new TFile("Parameters.root");
  if(f->IsOpen()) //The File Parameter.root exist
  {
    cout<<"Reading initial parameters from file Parameters.root"<<endl;
    fParam1=(TSCNParameter*)f->Get("SCNParameter");
    f->Close();
  }
  else //Data file with parameters
  {
    ifstream myfile;
    myfile.open ("Parameters.dat", ios::in);
    if(myfile.is_open()) //I have the file
    {
      fParam1 = (TSCNParameter*)  GetParameter("SCNParameter");
      string intro;
      myfile>>intro;
      cout<<intro<<endl;
      cout<<"Readind data file Parameters.dat"<<endl;
      int a;
      double b;
      for(int i=0;i<SCN_NUM_CHAN;i++)
      {
        myfile>>a;
        fParam1->SetPedestal(i,a);
      }
      myfile>>intro;
      for(int i=0;i<SCN_NUM_CHAN;i++)
      {
        myfile>>b;
        fParam1->SetFactor(i,b);
      }
      myfile.close();
    }
    else  //No Parameters at all. Go to default constructor
    {
      cout<<"No parameters file. Using default constructor"<<endl;
      fParam1 = (TSCNParameter*)  GetParameter("SCNParameter");
    }
  }
  fCal = new CalibParameter("CalibPar");
  AddParameter(fCal);
  // fCal = (MoDSSCalibParameter*) GetParameter("CalibPar");
  if (fCal) fCal->PrintParameter(0,0);
  else cout << "**** ERRR - CalibPar doesn't exist - program will crash.\n";

    fCorrel = new CorrelParameter("CorrelPar");
  AddParameter(fCorrel);
  if (fCorrel) fCorrel->PrintParameter(0,0);
  else cout << "**** ERRR - CorrelPar doesn't exist - program will crash.\n";

  //FINGER Polygon ToT vs strip gate window test
    Double_t ToTvalues[8]={0,8,23,25,35,45,25,0};
    Double_t Stripvalues[8]={17,22,27,32,37,42,47,50};
    TCutG* ToTvsStripcut = new TCutG("initialcut",8,ToTvalues,Stripvalues);
    fCond_FingToTvsStrip = new TGo4PolyCond("FING_TOTvsStrip");
    fCond_FingToTvsStrip -> SetValues(ToTvsStripcut);
    AddAnalysisCondition(fCond_FingToTvsStrip);
    fCond_FingToTvsStrip -> Enable();
    delete ToTvsStripcut;
    
    
    read_setup_parameters();
    get_used_Systems();
}
//-----------------------------------------------------------
EventAnlProc::~EventAnlProc()
{
  cout << "**** EventAnlProc: Delete" << endl;

}
//-----------------------------------------------------------
// void EventAnlProc::load_FingerID_File(){
//
//      const char* format = "%d %d %d";
//     ifstream data("Configuration_Files/Finger_allocation.txt");
//     if(data.fail()){
//         cerr << "Could not find Finger_allocation config file!" << endl;
//         exit(0);
//     }
// //     int id[5] = {0,0,0,0,0};
//     //int i = 0;
//     int tamid = 0;
//     int tamch = 0;
//     int fingid = 0;
//     string line;
//     //char s_tmp[100];
//     while(data.good()){
//
//         getline(data,line,'\n');
//         if(line[0] == '#') continue;
//         sscanf(line.c_str(),format,&tamid,&tamch,&fingid);
//
//
// //        fingID[tamid][tamch] = fingid;
//
//     }
// }

Bool_t EventAnlProc::BuildEvent(TGo4EventElement* dest)
{

  Bool_t isValid=kFALSE; // validity of output event

  EventUnpackStore* pInput  = (EventUnpackStore*) GetInputEvent();
  EventAnlStore* pOutput = (EventAnlStore*) dest;



  if((pInput==0) || !pInput->IsValid()){ // input invalid
    pOutput->SetValid(isValid); // invalid
    return isValid; // must be same is for SetValid
  }
  isValid=kTRUE;

   ///general inputs from the unpacker
    event_number = pInput->fevent_number;
  
  for (int i = 0; i<7; i++){
            if(pInput->fProcID[i]>-1){
      PrcID_Conv[i] = pInput->fProcID[i];
            //Used_Systems[i] = pInput->fUsed_Systems[i];

    }
  }
   static bool create =false;
  //Create histograms
  if (!create)
  {
    if (Used_Systems[1]) Make_Aida_Histos();
    if (Used_Systems[2]) Make_Plastic_VME_Histos();
    if (Used_Systems[3]) Make_Fatima_Histos();
    if (Used_Systems[4]) Make_Galileo_Histos();
    if (Used_Systems[5]) Make_Finger_Histos();
    if (Used_Systems[2] && Used_Systems[3]) Make_Fat_Plas_Histos();
    if (Used_Systems[5] && Used_Systems[2]) Make_Fing_Plas_Histos();
        create = true;
        }
                 /** Now extract the data from the stored Unpacker array (root tree)**/
    ///--------------------------------------/**FRS Input**/------------------------------------------///
    
      if (Used_Systems[0] && PrcID_Conv==0){
        ///MUSIC
       for(int i =0; i<2; ++i){
            FRS_dE[i] = pInput->fFRS_dE[i];
            FRS_dE_cor[i] = pInput->fFRS_dE_cor[i];
           }
        ///SCI
       for(int l=0;l<12;++l){
            FRS_sci_l[l] = pInput->fFRS_sci_l[l];
            FRS_sci_r[l] = pInput->fFRS_sci_r[l];
            FRS_sci_e[l] = pInput->fFRS_sci_e[l];
            FRS_sci_tx[l] = pInput->fFRS_sci_tx[l];
            FRS_sci_x[l] = pInput->fFRS_sci_x[l];
           }
        ///SCI TOF
        FRS_sci_tofll2 = pInput->fFRS_sci_tofll2;
        FRS_sci_tofll3 = pInput->fFRS_sci_tofll3;
        FRS_sci_tof2 = pInput->fFRS_sci_tof2;
        FRS_sci_tofrr2 = pInput->fFRS_sci_tofrr2;
        FRS_sci_tofrr3 = pInput->fFRS_sci_tofrr3;
        FRS_sci_tof3 = pInput->fFRS_sci_tof3;
        ///ID 2 4
        FRS_ID_x2 = pInput->fFRS_ID_x2;
        FRS_ID_y2 = pInput->fFRS_ID_y2;
        FRS_ID_a2 = pInput->fFRS_ID_a2;
        FRS_ID_b2 = pInput->fFRS_ID_b2;

        FRS_ID_x4 = pInput->fFRS_ID_x4;
        FRS_ID_y4 = pInput->fFRS_ID_y4;
        FRS_ID_a4 = pInput->fFRS_ID_a4;
        FRS_ID_b4 = pInput->fFRS_ID_b4;
            ///SCI dT
        FRS_sci_dt_21l_21r = pInput->fFRS_sci_dt_21l_21r;
        FRS_sci_dt_41l_41r = pInput->fFRS_sci_dt_41l_41r;
        FRS_sci_dt_42l_42r = pInput->fFRS_sci_dt_42l_42r;
        FRS_sci_dt_43l_43r = pInput->fFRS_sci_dt_43l_43r;

        FRS_sci_dt_21l_41l = pInput->fFRS_sci_dt_21l_41l;
        FRS_sci_dt_21r_41r = pInput->fFRS_sci_dt_21r_41r;

        FRS_sci_dt_21l_42l = pInput->fFRS_sci_dt_21l_42l;
        FRS_sci_dt_21r_42r = pInput->fFRS_sci_dt_21r_42r;
            ///ID Beta Rho
        for(int i =0; i<2; ++i){
        FRS_ID_brho[i] = pInput->fFRS_ID_brho[i];
        FRS_ID_rho[i] = pInput->fFRS_ID_rho[i];
        }
        FRS_beta = pInput->fFRS_beta;
        FRS_beta3 = pInput->fFRS_beta3;
        FRS_gamma  = pInput->fFRS_gamma;
            ///ID Z AoQ
        FRS_AoQ = pInput->fFRS_AoQ;
        FRS_AoQ_corr = pInput->fFRS_AoQ_corr;

        FRS_z = pInput->fFRS_z;
        FRS_z2 = pInput->fFRS_z2;
        FRS_z3 = pInput->fFRS_z3;
            ///ID Timestamp
        FRS_timestamp = pInput->fFRS_timestamp;
        FRS_ts = pInput->fFRS_ts;
        FRS_ts2 = pInput->fFRS_ts2;
            }
    

   ///-------------------------------- /**AIDA Input**/ --------------------------------///
      //  if (Used_Systems[1]&&  PrcID_Conv[1]==1) {ProcessAida(pInput);}
        ProcessAida(pInput);
        Aida_Fired = 0;
        for(int i=0; i<10000; i++) WR_Aida_Det_diff[i] = 0;

        Aida_Fired = pInput->fAIDAHits;

        for(int i=0; i< Aida_Fired; i++){
            WR_Aida_Det_diff[i]=pInput->fWR_Aida_Det_diff[i]; //in ns

        }
   ///-------------------------------- /**bPlastic VME Input**/ --------------------------------///

  bPlasQDCFired = 0;

  bPlasTDCFired = 0;
  bPlasTDC_ref = 0;
  bPlas_TDC_diff_sum = 0;

  for (int i=0; i<32; i++){
    bPlasQDC[i] = 0;
    bPlasQDCID[i] = 0;
    bPlas_TDC_Multiplicity[i] = 0;
  }
  for (int j=0; j<50; j++){
    bPlasTDCID[j] = 0;
    
  for (int k=0; k<32; k++){
    bPlasTDC_TS[j][k] = 0;
    bPlas_TDC_diff[k] = 0;
    }
  }
  
  if (Used_Systems[2]&& PrcID_Conv[2]==2){
    bPlasQDCFired =  pInput->fbPlas_VME_firedQDC;
    //QDC
    for (int i = 0; i<bPlasQDCFired; i++){
      bPlasQDCID[i] = pInput->fbPlas_VME_QDC_ID[i];
      bPlasQDC[i] = pInput->fbPlas_VME_QDC_E[i];
    }
    //TDC
    bPlasTDCFired = pInput->fbPlas_VME_firedTDC;
    
    for (int i = 0; i<bPlasTDCFired; i++){
      bPlasTDCID[i] = pInput->fbPlas_VME_TDC_ID[i];
      bPlasTDC_TS[i][bPlasTDCID[i]] = pInput->fbPlas_VME_TDC_TS[i][bPlasTDCID[i]];
          //  bPlas_TDC_Multiplicity[bPlasTDCID[i]]++ ;

      if(bPlasTDCID[i] == 0){
        bPlasTDC_ref= bPlasTDC_TS[i][bPlasTDCID[i]];
            }

        }
          Do_Plastic_VME_Histos(pOutput);
 ///--------------------------------------/**Scalar Input**/------------------------------------------///

    ScalarFired = pInput->fScalar_fired;
    for (int i = 0; i<ScalarFired; i++){
      ScalarID = pInput->fScalar_ID;
    }
  }

  ///--------------------------------------/**Fatima Input**/------------------------------------------///
        FatQDCFired = 0;
        FatTDCFired = 0;
        SC41 = 0;
        SC41_ns =0;
        Fat_CHA_0_TDC = 0;
        Fat_WR = 0;
  for (int i=0; i<50; i++){
         FatQDC[i] = 0;
         FatQDC_T[i] =0;
  }
  for (int j=0; j<50; j++){
          FatTDCID[j] = -1;
          FatQDCID[j] = -1;
          FatTDC_Multipl[j] = 0;

          for (int k=0; k<50; k++){
      FatTDC_TS[j][k] = 0;

    }
  }
  if (Used_Systems[3]&& PrcID_Conv[3]==3){
    FatQDCFired =  pInput->fFat_firedQDC;
    Fat_WR = pInput->fFat_WR;
    //QDC
    for (int i = 0; i<FatQDCFired; i++){
      FatQDCID[i] = pInput->fFat_QDC_ID[i];
      FatQDC[i] = pInput->fFat_QDC_E[i];
      FatQDC_T[i] = pInput->fFat_QDC_T[i];
                                       }

    //TDC
    FatTDCFired =  pInput->fFat_firedTDC;
    for (int i = 0; i<FatTDCFired; i++){
      FatTDCID[i] = pInput->fFat_TDC_ID[i];
            FatTDC_TS[i][FatTDCID[i]] = (pInput->fFat_TDC_TS[i][FatTDCID[i]])*0.025;
            FatTDC_Multipl[FatTDCID[i]] = pInput-> fFat_TDC_Multiplicity[FatTDCID[i]];
            ///Reference channel
             if( FatTDC_TS[i][0]>0  && FatTDCID[i] ==0){
                           Fat_CHA_0_TDC =   FatTDC_TS[i][0];
                           pOutput->pFat_Ch0_TDC = Fat_CHA_0_TDC;
                         }
            
            //SC41 Trigger
      if (FatTDCID[i]==40){
            SC41 = pInput ->fSC41[i]; //in 25ps
            SC41_ns = SC41*0.025; //SC41 signal in ns
      }
    }

        Do_Fatima_Histos(pOutput);
  }

 ///--------------------------------------/**Galileo Input**/------------------------------------------///
   GalFired = -1;
   //Gal_WR = 0;
   for(int g = 0; g<20; g++){
      GalID[g] = -1;
      GalE[g] = -1;
      GalT[g] =-1;
   }
   if (Used_Systems[4] && PrcID_Conv[4]==4){
        
      
    GalFired =  pInput->fGal_fired;
    GalPileup = pInput->fGal_Pileup;
    Gal_WR = pInput->fGal_WR;
    //for(int f=0;f<1000;f++){
   
    //}

    for (int i = 0; i<GalFired; i++){
      GalID[i] = pInput->fGal_ID[i];
      GalE[GalID[i]] = pInput->fGal_E[GalID[i]];
      GalT[GalID[i]] = pInput->fGal_T[GalID[i]];
    }
    

        Do_Galileo_Histos(pOutput);
  }

 ///--------------------------------------/**Finger Input**/------------------------------------------///
        Fing_firedTamex = -1;
        for (int i = 0; i<4; i++)
        {
            Fing_leadHits[i] = -1;
            Fing_trailHits[i] = -1;
            Fing_iterator[i] = -1;
            Fing_trig[i] = 0;

        for (int j = 0; j<32; j++){
            Fing_tamex_ch[i][j] = -1;
            Fing_leadChan[i][j] = -1;
            Fing_leadT[i][j] = 0;
            Fing_trailChan[i][j] = -1;
            Fing_trailT[i][j] = 0;
            Fing_TOT[i][j] = 0;
            Fing_TOT_added[i][j] = 0;
            Fing_chID[i][j] = -1;
            Fing_lead_coarse[i][j] =  -1;
            Fing_lead_fine[i][j]  =-1;
            Fing_trail_coarse[i][j] = -1;
            Fing_trail_fine[i][j] = -1;
            }
        }
  if (Used_Systems[5] && PrcID_Conv[5]==5){
          Fing_firedTamex = pInput->ffing_tamexhits;
          maxToT = Fing_TOT[0][0];
          maxToT_added = Fing_TOT_added[0][0];
    
          for (int i=0; i< Fing_firedTamex; i++){
            Fing_leadHits[i] = pInput-> ffing_leadHits[i];
            Fing_trailHits[i] = pInput-> ffing_trailHits[i];
            Fing_iterator[i] = pInput-> ffing_iterator[i];
            Fing_trig[i] = pInput->ffing_Trig[i];

                for (int j =0; j<Fing_iterator[i]; j++){
                    Fing_tamex_ch[i][j] =  pInput->ffing_tamexCh[i][j]; 
                    Fing_lead_coarse[i][j] =  pInput->ffing_lead_coarse[i][j];
                    Fing_lead_fine[i][j]  = pInput->ffing_lead_fine[i][j];
                    Fing_trail_coarse[i][j] =  pInput->ffing_trail_coarse[i][j];
                    Fing_trail_fine[i][j] =  pInput->ffing_trail_fine[i][j];
                    
                    Fing_chID[i][j] = pInput->ffing_chID[i][j];
                  //  cout << "  Fing_chID[i][j] " << Fing_chID[i][j]<< " i " << i << " j " << j << " ffing_Lead_T[i][j] " << pInput->ffing_Lead_T[i][j] <<endl;
                    if(Fing_chID[i][j] % 2 == 1){
                        Fing_leadChan[i][j] = pInput->ffing_Lead_Phys_Chan[i][j];
                        Fing_leadT[i][j] = pInput->ffing_Lead_T[i][j];
                       
                        }                                                                    

                else{
                     Fing_trailChan[i][j] = pInput->ffing_Trail_Phys_Chan[i][j];
                     Fing_trailT[i][j] = pInput->ffing_Trail_T[i][j];
                    
                }
                ///Note: ToT value here is only for the 'up' PMTs
                    Fing_TOT[i][j] = pInput->ffing_TOT[i][j];
                    /// Up PMT ToT + Down PMT ToT
                    Fing_TOT_added[i][j] =   pInput->ffing_TOT_added[i][j];
            //    cout << "1) ev "<< event_number <<" Fing_TOT[i][j] "<< Fing_TOT[i][j] << " i " << i << " j " << j << endl;
                //for(int l=0;l<Fing_leadChan[i][j]; l++){
                        if(Fing_chID[i][j] % 2 == 1){
                             
                          if(maxToT<Fing_TOT[i][j] && Fing_leadChan[i][j]>0 ){
                            
                                maxToT = Fing_TOT[i][j];
                                maxToTChan = Fing_leadChan[i][j];
                        }
                        //Get max ToT for PMT pairings
                        if(maxToT_added<Fing_TOT_added[i][j]&& Fing_leadChan[i][j]>0 ){
                            maxToT_added = Fing_TOT_added[i][j];
                            maxToT_added_Chan = Fing_leadChan[i][j];
                            
                        }
             //  cout <<"2) ev " << event_number <<" maxtot " <<maxToT <<endl;
                    
                }
            }
          }
    //    if(maxToT>-1){
    //   cout <<"1a) ev " << event_number <<  " maxToT  " << maxToT <<" maxToTChan " << maxToTChan<< endl;
 
            Do_Finger_Histos(pOutput);
  }

///------------------------/**Setup for some temporary correlations**/----------------------------------///
  if (Used_Systems[2]&& PrcID_Conv[2]==2 && Used_Systems[3]&& PrcID_Conv[3]==3){
            Do_Fat_Plas_Histos(pOutput);
            }
  if (Used_Systems[5]&& PrcID_Conv[5]==5 && Used_Systems[2]&& PrcID_Conv[2]==2){
    Do_Fing_Plas_Histos();
            }
                        /** End of Unpack Tree input**/

  BuildEvent2(dest);
  pOutput->SetValid(isValid);
  return isValid;

 }  //End of BuildEvent

    Bool_t EventAnlProc::BuildEvent2(TGo4EventElement* dest){
//        EventUnpackStore* pInput  = (EventUnpackStore*) GetInputEvent();
//         int event_number1 = pInput -> fevent_number;
//         int event_counter = pInput -> fArray_count;

       //cout <<"event " << event_number << " WR " <<  pInput ->fWR_main_array[event_counter] << " event_counter " << event_counter << endl;
      return(1);
        }
///End of Input from Unpacker ///
 ///-----------------------------------------------------------------------------------------------------------------///
   void EventAnlProc::get_used_Systems(){
    for(int i = 0;i < 7;i++) Used_Systems[i] = false;

    ifstream data("Configuration_Files/Used_Systems.txt");
    if(data.fail()){
        cerr << "Could not find Used_Systems config file!" << endl;
        exit(0);
    }
    int i = 0;
    int id = 0;
    string line;
    char s_tmp[100];
    while(data.good()){
        getline(data,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),"%s %d",s_tmp,&id);
        Used_Systems[i] = (id == 1);
        i++;
    }
    string DET_NAME[6] = {"FRS","AIDA","PLASTIC","FATIMA","GALILEO","FINGER"};

    cout << "\n=====================================================" << endl;
    cout << "USED SYSTEMS" << endl;
    cout << "-----------------------------------------------------" << endl;
    for(int j = 0;j < 6;++j){
        if(Used_Systems[j]) cout << DET_NAME[j] << endl;
    }
    cout << "=====================================================" << endl;


}
///-----------------------------------------------------------------------------------------------------------------///
  void EventAnlProc::read_setup_parameters(){

    // unused // const char* format = "%s %d";

    ifstream file("Configuration_Files/Detector_System_Setup_File.txt");

    if(file.fail()){
        cerr << "Could not find File for setup parameters!" << endl;
        exit(0);
    }

    string line;
    string var_name;
    // unused //int dummy_var;
    //file.ignore(256,'GENERAL_CONFIGURATION');

    file.ignore(256,':');
    file >> FAT_exclusion_dist;//dummy_var;

    file.ignore(256,':');
    file >> FAT_nearest_neighbour_exclusion;//dummy_var;

    file.ignore(256,':');
    file >> same_ring_exclusion;//dummy_var;

    file.ignore(256,':');
    file >> output_position_matrix;//dummy_var;

    cout<<endl;
    cout<<endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"Setup Parameters List Analysis Proc: "<<endl;
    if(FAT_exclusion_dist > 0) cout<<"FATIMA Detectors Excluded if Linear Difference Exceeds "<<FAT_exclusion_dist<<" mm"<<endl;
    else if(FAT_exclusion_dist == 0) cout<<"'Nearest Neighbour Exclusion': Disabled (Distance set to 0)"<<endl;
    cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<endl;
    cout<<endl;
    /*while(file.good()){
        getline(file,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format,&var_name,&dummy_var);

        cout<<"Hello Again?"<<endl;

        if (var_name == "White_Rabbit_Enabled:" && dummy_var == 1)  WHITE_RABBIT_USED = true;
        else if (var_name == "White_Rabbit_Enabled:" && dummy_var == 0)  WHITE_RABBIT_USED = false;

        if (var_name == "FATIMA_Gain_Match_Enabled:" && dummy_var == 1)  FAT_gain_match_used = true;
        else if (var_name == "FATIMA_Gain_Match_Enabled:" && dummy_var == 0) FAT_gain_match_used  = false;

    }*/

}
/**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------    FRS   ---------------------------------------------**/
 /**----------------------------------------------------------------------------------------------**/
// void EventAnlProc::Make_FRS_Histos(){
//      hFRS_z1_z2 = MakeTH2('D',"FRS/z1vsz2","FRS Z1 vs Z2", 500, 0, 5000,  500, 0, 5000);
//      hFRS_AoQ_s2X = MakeTH2('D',"FRS/z1vsz2","FRS Z1 vs Z2", 500, 0, 5000,  500, 0, 5000);
//     
//     
// }
// 
// void EventAnlProc::Do_FRS_Histos(EventAnlStore* pOutput){
//     hFRS_z1_z2->Fill(FRS_z,FRS_z2);
//     
//     
// }

/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------    AIDA   ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/

void EventAnlProc::Make_Aida_Histos(){
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
  implants_strip_xy.resize(conf->DSSDs());
  implants_pos_xy.resize(conf->DSSDs());
  implants_e.resize(conf->DSSDs());
  implants_e_xy.resize(conf->DSSDs());
  implants_time_delta.resize(conf->DSSDs());
  implants_strip_1d.resize(conf->DSSDs());
  implants_per_event.resize(conf->DSSDs());
  decays_strip_xy.resize(conf->DSSDs());
  decays_pos_xy.resize(conf->DSSDs());
  decays_e.resize(conf->DSSDs());
  decays_e_xy.resize(conf->DSSDs());
  decays_time_delta.resize(conf->DSSDs());
  decays_strip_1d.resize(conf->DSSDs());
  decays_per_event.resize(conf->DSSDs());
  implants_channels.resize(conf->DSSDs());
  decays_channels.resize(conf->DSSDs());
  for (int i = 0; i < conf->DSSDs(); ++i)
  {
    implants_strip_xy[i] = MakeTH2('I', Form("AIDA/Implants/DSSD%d_implants_strip_XY", i+1), Form("DSSD %d implant hit pattern", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
    implants_pos_xy[i] = MakeTH2('D', Form("AIDA/Implants/DSSD%d_implants_pos_XY", i+1), Form("DSSD %d implant position", i+1), 128, -37.8, 37.8, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    implants_e[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_energy", i+1), Form("DSSD %d implant energy", i+1), 1000, 0, 10000, "Implant Energy/MeV");
    implants_e_xy[i] = MakeTH2('F', Form("AIDA/Implants/DSSD%d_implants_energy_XY", i+1), Form("DSSD %d implant front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 10000, "X Energy", "Y Energy");
    implants_time_delta[i] = MakeTH1('F', Form("AIDA/Implants/DSSD%d_implants_time_delta", i+1), Form("DSSD %d implant front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    implants_strip_1d[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_strip_1d", i+1), Form("DSSD %d implant 1D hit pattern", i+1), 256, 0, 256, "Strip number");
    implants_per_event[i] = MakeTH1('I', Form("AIDA/Implants/DSSD%d_implants_per_event", i+1), Form("DSSD %d implants per event", i+1), 100, 0, 100, "Number of implants");

    decays_strip_xy[i] = MakeTH2('I', Form("AIDA/Decays/DSSD%d_decays_strip_XY", i+1), Form("DSSD %d decay hit pattern", i+1), 128, 0, 128, 128, 0, 128, "X strip", "Y strip");
    decays_pos_xy[i] = MakeTH2('D', Form("AIDA/Decays/DSSD%d_decays_pos_XY", i+1), Form("DSSD %d decay position", i+1), 128, -37.8, 37.8, 128, -37.8, 37.8, "X position/mm", "Y position/mm");
    decays_e[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_energy", i+1), Form("DSSD %d decay energy", i+1), 1000, 0, 20000, "Decay Energy/keV");
    decays_e_xy[i] = MakeTH2('F', Form("AIDA/Decays/DSSD%d_decays_energy_XY", i+1), Form("DSSD %d decay front energy vs back energy", i+1), 1000, 0, 10000, 1000, 0, 20000, "X Energy", "Y Energy");
    decays_time_delta[i] = MakeTH1('F', Form("AIDA/Decays/DSSD%d_decays_time_delta", i+1), Form("DSSD %d decay front vs back time", i+1), 1000, -10000, 10000, "Time Difference/ns");
    decays_strip_1d[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_strip_1d", i+1), Form("DSSD %d decay 1D hit pattern", i+1), 256, 0, 256, "Strip number");
    decays_per_event[i] = MakeTH1('I', Form("AIDA/Decays/DSSD%d_decays_per_event", i+1), Form("DSSD %d decays per event", i+1), 100, 0, 100, "Number of decays");

    implants_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_implants_channels", i+1), Form("DSSD %d number of implant channels", i+1), 256, 0, 256);
    decays_channels[i] = MakeTH1('I', Form("AIDA/DSSD%d_decays_channels", i+1), Form("DSSD %d number of decay channels", i+1), 769, 0, 769);
  }
}

//////////Process AIDA//////////////////////////
void EventAnlProc::ProcessAida(EventUnpackStore* pInputMain) {
 // int Aida_hits =0;
//       double bPlasQDCGainMatch_AIDA[32] ={0};
  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();

  //      cout <<"event_number " << event_number << " pInput->ImplantEvents.size() " << pInput->ImplantEvents.size()<< ": pInput->DecayEvents.size() " << pInput->DecayEvents.size()<<endl;
  //         Aida_hits = pInput->AIDAHits;
  //         cout << " evt.Channel " << evt.Channel << endl;

  for(AidaUnpackData& pInputD : pInputMain->Aida)
  {

    AidaUnpackData* pInput = &pInputD;

    if (pInput->ImplantEvents.size() > 1)
    {
      //         cout << " pInput->ImplantEvents.size() " << pInput->ImplantEvents.size() <<  endl;

      implantEvents++;


      // Cluster events on adjecent strips into one
      std::vector<AidaCluster> clusters = EventsToClusters(pInput->ImplantEvents);
      //
      //     // Match front-back clusters which define physical hits on the detector
      std::vector<std::pair<AidaCluster, AidaCluster>> hits;
      //
      //

      for (auto& i : clusters)
      {
        if(i.DSSD == -1 || i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        for (auto& j : clusters)
        {
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          /// Gates (set in TAidaConfiguration)
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyH() && i.IsGoodTime(j, conf->FrontBackWindow()))
          {
            hits.push_back({i, j});
          }
        }
      }

      int channels[768] = {0};
      for (auto& i : pInput->ImplantEvents)
      {
        channels[i.Module * 64 + i.Channel]++;
      }
      int channelM = 0;
      for (int i = 0; i < 768; ++i){
        if (channels[i]) ++channelM;
      }
      implants_channels[0]->Fill(channelM);

      //     // Generate stored data for hits and plot the histograms
      for (auto& i : hits)
      {
        AidaHit hit = ClusterPairToHit(i);
        //
        pInput->Implants.push_back(hit);

     /*   for (int i=0; i<bPlasQDCID; i++){
    //cout <<" bPlasQDCGainMatch_AIDA[i] " << bPlasQDCGainMatch_AIDA[bPlasQDCID]<<endl;

            for (int j =0; j<Aida_Fired; j++){

//             if(bPlasQDCGainMatch_AIDA[bPlasQDCID]>0){
//                                             }
    if(WR_Aida_Det_diff[j]/1000>-10 || WR_Aida_Det_diff[j]/1000<10){
     bPlasQDCGainMatch_AIDA[i] = fCal->AplasQDC[i]*bPlasQDC[i] + fCal->BplasQDC[i];
        //   cout <<" bPlasQDCGainMatch_AIDA[i][j] " << bPlasQDCGainMatch_AIDA[i][j]<<endl;
        //  cout << "WR_Aida_Det_diff[j] " <<WR_Aida_Det_diff[j] << " bPlasQDCGainMatch_AIDA[i] "<< bPlasQDCGainMatch_AIDA[bPlasQDCID] <<endl;

                }
            }
        } */

        implants_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
        implants_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        implants_e[hit.DSSD - 1]->Fill(hit.Energy);
        implants_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        //implants_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        implants_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        int channel = i.first.Strip;
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + 128;
        implants_strip_1d[hit.DSSD - 1]->Fill(channel);
      }
      //
      if (hits.size() > 0)
      {
        implants_per_event[0]->Fill(hits.size());
      }
    }
    else if (pInput->DecayEvents.size() > 1)
    {
      decayEvents++;

      int channels[768] = {0};
      for (auto& i : pInput->DecayEvents)
      {
        channels[i.Module * 64 + i.Channel]++;
      }
      int channelM = 0;
      for (int i = 0; i < 768; ++i)
      if (channels[i]) ++channelM;
      decays_channels[0]->Fill(channelM);
      //     cout <<"i " <<i << " channelM " << channelM<< endl;
      if (channelM > 400)
      {
        decayEvents--;
        pulserEvents++;
        return;
      }

      // Clean up huge event buffers - for now we just destroy them
      if (pInput->DecayEvents.size() > 400)
      {
        decayEvents--;
        nonsenseEvents++;
        //pInput->SetValid(kFALSE);
        return;
      }

      std::vector<AidaCluster> clusters = EventsToClusters(pInput->DecayEvents);

      std::vector<std::pair<AidaCluster, AidaCluster>> hits;

      for (auto& i : clusters)
      {

        if(i.DSSD == -1 || i.Side != conf->DSSD(i.DSSD -1).XSide) continue;

        //if(i.Energy < 100) continue;
        for (auto& j : clusters)
        {
          if(j.DSSD != i.DSSD || j.Side != conf->DSSD(j.DSSD -1).YSide) continue;
          //if(j.Energy < 100) continue;
          // Gates
          if (abs(i.Energy - j.Energy) < conf->FrontBackEnergyL() && i.IsGoodTime(j, conf->FrontBackWindow()))
          hits.push_back({i, j});
        }
      }

      for (auto& i : hits)
      {
        AidaHit hit = ClusterPairToHit(i);

        pInput->Decays.push_back(hit);
        decays_strip_xy[hit.DSSD - 1]->Fill(hit.StripX, hit.StripY);
        decays_pos_xy[hit.DSSD - 1]->Fill(hit.PosX, hit.PosY);
        decays_e[hit.DSSD - 1]->Fill(hit.Energy);
        decays_e_xy[hit.DSSD - 1]->Fill(hit.EnergyFront, hit.EnergyBack);
        //decays_time_delta[hit.DSSD - 1]->Fill(hit.TimeFront - hit.TimeBack);
        decays_time_delta[hit.DSSD - 1]->Fill(hit.FastTimeFront - hit.FastTimeBack);

        int channel = i.first.Strip;
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
        channel = i.second.Strip + 128;
        decays_strip_1d[hit.DSSD - 1]->Fill(channel);
      }

      if (clusters.size() > 0)
      {
        decays_per_event[0]->Fill(clusters.size());
      }
      //   cout <<" decayEvents " << decayEvents << endl;
    }
    else
    {
      nonsenseEvents++;
      //pInput->SetValid(kFALSE);
    }

  }
}


std::vector<AidaCluster> EventAnlProc::EventsToClusters(std::vector<AidaEvent> const& events)
{
  std::vector<AidaCluster> clusters;
  for (auto& i : events)
  {
    // Don't cluster invalid events
    if (i.DSSD == -1) continue;

    bool added = false;

    // Try to add the event to an existing cluster
    for (auto& j : clusters)
    {
      if(j.IsAdjacent(i) && j.IsGoodTime(i))
      {
        j.AddEvent(i);
        added = true;
        break;
      }
    }

    // Otherwise make a new cluster for the event
    if (!added)
    {
      AidaCluster c_test;
      c_test.AddEvent(i);
      clusters.push_back(c_test);
    }
  }
  return clusters;
}

AidaHit EventAnlProc::ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const& i)
{
  AidaHit hit;
  hit.DSSD = i.first.DSSD;

  hit.StripX = i.first.Strip;
  hit.StripY = i.second.Strip;
  hit.PosX = 75.6 * i.first.Strip / 128. - 37.75;
  hit.PosY = 75.6 * i.second.Strip / 128. - 37.75;

  hit.StripXMin = i.first.StripMin;
  hit.StripXMax = i.first.StripMax;
  hit.StripYMin = i.second.StripMin;
  hit.StripYMax = i.second.StripMax;

  hit.Energy = (i.first.Energy + i.second.Energy) / 2;
  hit.EnergyFront = i.first.Energy;
  hit.EnergyBack = i.second.Energy;

  hit.Time = std::min(i.first.Time, i.second.Time);
  hit.TimeFront = i.first.Time;
  hit.TimeBack = i.second.Time;
  hit.FastTime = std::min(i.first.FastTime, i.second.FastTime);
  hit.FastTimeFront = i.first.FastTime;
  hit.FastTimeBack = i.second.FastTime;

  return hit;
}
///End of Aida ///
/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  bPLASTIC VME (+Scalar)  ----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/

void EventAnlProc::Make_Plastic_VME_Histos(){
  for (int i=0; i<32; i++){
    hPLAS_QDCCalib1[i] =  MakeTH1('D', Form("bPlastic/Energy/QDC1Calib/QDC1Calib_Ch.%2d",i), Form("QDC1 Calib Ch. %2d",i), 20000, 0., 20000.);
   // hPLAS_TDCCalib1[i] =  MakeTH1('D', Form("bPlastic/Timing/TDC1Calib/TDC1Calib_Ch.%2d",i), Form("TDC1 Ch. %2d",i), 2E4, 0, 2E5);
    hPLAS_TimeDiffSiPM_Ch_Raw[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SiPM1-SiPMRaw/TDCSiPM_dT_Ch.%02d",i),Form("SiPM dT(ns) Ch1 - Ch%0d",i),10000,-10000,10000);
    hPLAS_TimeDiffSiPM_Ch_Calib[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SiPM1-SiPMCalib/TDCSiPM1_dT_Ch.%02d",i),Form("TimeDiff Gainmatched Ch1 - Ch%0d",i),10000,-10000,10000);
    hPLAS_TimeDiffSiPM_Ch_Calib_Egated[i] = MakeTH1('D',Form("bPlastic/Timing/EGated/TDCdt_SiPM1-SiPM/TDCSiPM1_dT_EGated_Ch.%02d",i),Form("TimeDiff Gainmatched Ch1 - Ch%0d",i),10000,-10000,10000);

    //hPLAS_TimeDiff_Ch_Raw[i] = MakeTH1('D',Form("bPlastic/Timing/Raw/TDCdt_ref-plas/TDCSiPM1_dt_Raw_Ch.%02d",i),Form("TimeDiff Ch0 - Ch%0d",i),1E5,-2E5,2E5);
   // hPLAS_TimeDiffSiPM_Ch[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SiPM1-SiPm/TDCSiPMdt_Ch.%02d",i),Form("TimeDiff Ch1 - Ch%0d",i),1E5,-2E5,2E5);
    ///new
    hPLAS_TDC_FiredRatio[i] = MakeTH1('D', Form("bPlastic/Timing/TDC1_FiredRatio/TDC_FiredRatio_Ch.%2d",i), Form("TDC1 Ratio CalibTDC/FiredTDC Ch. %2d",i), 4000, 0, 4000);
    hPLAS_TimeDiff_SC41_Raw[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SC41-SiPMRaw/TDCdT_SC41_Plas_RawCh.%02d",i),Form("SC41-SiPM dT(ns) Raw Ch%0d",i),10000,-10000,10000);
    hPLAS_TimeDiff_SC41_Calib[i] = MakeTH1('D',Form("bPlastic/Timing/TDCdt_SC41-SiPMCalib/TDCdT_SC41_Plas_Ch.%02d",i),Form("SC41-SiPM dT Calib  Ch%0d",i),10000,-10000,10000);
    hPLAS_TimeDiff_SC41_Calib_Egated[i] = MakeTH1('D',Form("bPlastic/Timing/EGated/TDCdt_SC41-SiPM/TDCdT_SC41_Plas_EGated_Ch.%02d",i),Form("SC41-SiPM dT Calib  Ch%0d",i),10000,-10000,10000);

    hPLAS_TDC_multich[i] = MakeTH1('D', Form("bPlastic/Stats/TDC_MultiCh/TDCMch%2d",i), Form("TDC channel Multi %2d",i), 50, 0, 50);
    hPLAS_CoincE1E2[i] = MakeTH2('D',Form("bPlastic/Energy/CoincEnergy/Coinc_Energy_Energy_Ch.%02d",i), Form("Coinc_Energy_Energy_Sum_Ch.%0d",i), 500, 0, 5000,  500, 0, 5000);


  }
    //hPLAS_E1E2TimeGated[i] = MakeTH2('D',Form("bPlastic/CoincEEGated/CoincEEGated%02d",i), Form("CoincGated%0d",i), 400, 0, 4000.0,  400, 0, 4000.0);

    hPLAS_CoincE_dTSiPM_Ch[7] =   MakeTH2('D',Form("bPlastic/Timing/Coinc_Energy_SiPM_dT/Coinc_E_dTCh1-Ch.%02d",7), Form("Energy vs. SiPM1-SiPMCh.%0d",7), 2500, 0, 5000, 2000,-10000,10000);

    //Sum spectra
    hPLAS_QDCCalib1Sum = MakeTH1('D',"bPlastic/Energy/QDCSum","QDC Calibrated Sum",2000,0,20000);
    hPLAS_TDCCalib1Sum = MakeTH1('D',"bPlastic/Timing/TDCSum","TDC Calibrated Sum (ns)",2E4, 0, 2E5);
    hPLAS_TimeDiffSiPM_Ch_Sum = MakeTH1('D',"bPlastic/Timing/TDCdT_SiPM1-SiPM","SiPM1 - SiPM Ch.x dT(ns) (calibrated)",2000,-1000,1000);
    hPLAS_TimeDiffSiPM_Ch_Sum_M1 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SiPM1-SiPM_Multi1","SiPM1 - SiPM Ch.x dT(ns) Multiplicity1 (calibrated)",2000,-1000,1000);
    hPLAS_TimeDiffSiPM_Ch_Sum_M2 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SiPM1-SiPM_Multi2","SiPM1 - SiPM Ch.x dT(ns) Multiplicity2 (calibrated)",2000,-1000,1000);
    hPLAS_TimeDiffSiPM_Ch_Sum_M3P = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SiPM1-SiPM_Multi3","SiPM1 - SiPM Ch.x dT(ns) Multiplicity>3 (calibrated)",2000,-1000,1000);

    hPLAS_TimeDiffSiPM_Ch_Sum_Egated = MakeTH1('D',"bPlastic/Energy/EGated/TDCdT_SiPM1-SiPM_EGated","Energy gated SiPM1 - SiPM  dT(ns) Ch.x (calibrated)",2000,-1000,1000);

    hPLAS_TimeDiff_SC41_Sum = MakeTH1('D',"bPlastic/Timing/TDCdT_SC41-SiPM","SC41 - SiPM dT(ns) (calibrated)",2000,-1000,1000);
    hPLAS_TimeDiff_SC41_Sum_M1 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SC41-SiPM_Multi1","SC41 - SiPM  Ch.x dT(ns) Multiplicity1 (calibrated)",2000,-1000,1000);
    hPLAS_TimeDiff_SC41_Sum_M2 = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SC41-SiPM_Multi2","SC41 - SiPM  Ch.x dT(ns) Multiplicity2 (calibrated)",2000,-1000,1000);
    hPLAS_TimeDiff_SC41_Sum_M3P = MakeTH1('D',"bPlastic/Timing/Multiplicity/TDCdT_SC41-SiPM_Multi3+","SC41 - SiPM  Ch.x dT(ns) Multiplicity >3 (calibrated)",2000,-1000,1000);


    hPLAS_TimeDiff_SC41_Sum_Egated  = MakeTH1('D',"bPlastic/Timing/EGated/TDCdT_SC41-SiPM_EGated","Energy gated SC41 - SiPM dT(ns) Ch.x (calibrated)",2000,-1000,1000);
    hPLAS_TDC_FiredRatio_Sum = MakeTH1('D',"bPlastic/Timing/TDC_FiredRatio","TDC1 Ratio CalibTDC/FiredTDC",4000, 0, 4000);

    hPLAS_QDC1_hits  = MakeTH1('D',"bPlastic/Stats/QDC1_hits","bPlastic hit pattern QDC1",32,0,32);
    hPLAS_TDC_hits  = MakeTH1('D',"bPlastic/Stats/TDC_hits","bPlastic hit pattern TDC",32,0,32);
    hPLAS_TDC_multi = MakeTH1('D',"bPlastic/Stats/TDC_Multi","bPlastic TDC Multiplicity",50,0,50);

    hPLAS_CoincE1E2_Sum = MakeTH2('D',"bPlastic/Energy/Coinc_Energy_Energy_Sum","bPlastic Energy-Energy Sum", 500, 0, 5000,  500, 0, 5000);
    hPLAS_CoincE_dTSiPM_Sum = MakeTH2('D',"bPlastic/Energy/CoincEnergy_SiPM1-SiPMCh.x","Energy vs. SiPMCh.1-SiPM_all", 5000, 0, 5000, 100,-100,100);

    hScalar_hit_pattern = MakeTH1('D',"Scalar/HitPat","Scalar Hit pattern",32,0,32);
}
//------------------------------------------------------------------------------------------------------------------------//
void EventAnlProc::Do_Plastic_VME_Histos(EventAnlStore* pOutput){

  int bPlasTDCIDMain;
    int bPlasQDCID_i, bPlasQDCID_j;
    double bPlasTDC_TS_Raw[32],    bPlasTDC_T_Calib[32] ;
    double bPlasTDC_TS_Firedratio[32];
    double bPlas_SiPM_dT_Raw[32], bPlas_SiPM_dT_Calib[32];
    double bPlas_SC41_dT_Raw[32], bPlas_SC41_dT_Calib[32];
    double bPlasQDCGainMatch_i[32],bPlasQDCGainMatch_j[32];
    double bPlas_TDC_Cha1;
  bPlasTDCIDMain = -1;
    bPlas_TDC_Cha1 = 0;
    bPlasQDCID_i=0;
    bPlasQDCID_j=0;
  for (int i=0; i< 32; i++){

        bPlasTDC_TS_Raw[i]=0;
        bPlasTDC_T_Calib[i]=0;
        bPlasQDCGainMatch_i[i] = 0;
        bPlasQDCGainMatch_j[i] = 0;
        bPlas_SiPM_dT_Raw[i] = 0;
        bPlas_SiPM_dT_Calib[i] = 0;
        bPlas_SC41_dT_Raw[i] = 0;
        bPlas_SC41_dT_Calib[i] = 0;
  }
  /**------------------bPlastic Energy -----------------------------------------**/
         pOutput->pbPlas_QDCFired = bPlasQDCFired;
    for (int i=0; i<bPlasQDCFired; i++){
         bPlasQDCID_i = bPlasQDCID[i];
         pOutput->pbPlas_QDCID[i] = bPlasQDCID[i];
         bPlasQDCGainMatch_i[bPlasQDCID_i] = fCal->AplasQDC[bPlasQDCID_i]*bPlasQDC[bPlasQDCID_i] + fCal -> BplasQDC[bPlasQDCID_i]; //gain matching
          hPLAS_QDC1_hits->Fill(bPlasQDCID_i);
        
          if(bPlasQDCGainMatch_i[bPlasQDCID_i]>30){
             pOutput->pbPlas_QDCGainMatch_i[bPlasQDCID_i] = bPlasQDCGainMatch_i[bPlasQDCID_i];
             hPLAS_QDCCalib1[bPlasQDCID_i]->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i]);
             hPLAS_QDCCalib1Sum->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i]);
                             }
               ///Energy-Energy Matrix
                for(int j=0; j<bPlasQDCFired; j++){
                  bPlasQDCID_j = bPlasQDCID[j];
                   ///Dont loop on the first hit again: (check)
                  if(bPlasQDCID_i<bPlasQDCID_j){
                  bPlasQDCGainMatch_j[bPlasQDCID_j] = fCal->AplasQDC[bPlasQDCID_j]*bPlasQDC[bPlasQDCID_j] + fCal -> BplasQDC[bPlasQDCID_j]; //gain matching
                  hPLAS_CoincE1E2[bPlasQDCID_i] ->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i],bPlasQDCGainMatch_j[bPlasQDCID_j]);
                  hPLAS_CoincE1E2_Sum->Fill(bPlasQDCGainMatch_i[bPlasQDCID_i],bPlasQDCGainMatch_j[bPlasQDCID_j]);

               }
            }
          }
     /**----------------------------bPlastic Timing -----------------------------------------**/
     ///Channel 0 is SC41 trigger
        pOutput->pbPlas_TDCFired = bPlasTDCFired;
  for (int i=0; i<bPlasTDCFired;i++){

          bPlasTDCIDMain = bPlasTDCID[i];
          pOutput->pbPlas_TDCID[i] = bPlasTDCID[i];
          bPlasTDC_TS_Raw[bPlasTDCIDMain] = bPlasTDC_TS[i][bPlasTDCIDMain];
          bPlas_TDC_Multiplicity[bPlasTDCIDMain]++;
          pOutput-> pbPlas_TDC_Multiplicity[bPlasTDCIDMain] = bPlas_TDC_Multiplicity[bPlasTDCIDMain];
          hPLAS_TDC_hits->Fill(bPlasTDCIDMain);
          hPLAS_TDC_multi->Fill(bPlas_TDC_Multiplicity[bPlasTDCIDMain]);
       //  cout <<"FAT Event " << event_number <<" bPlasTDC_TS_Raw[bPlasTDCIDMain]" <<bPlasTDC_TS_Raw[bPlasTDCIDMain]  << " bPlasTDCID[i] " << bPlasTDCID[i] << endl;

          //cout <<"ev " << event_number << " fired " <<  bPlasTDCFired << " i " << i <<" bPlasTDCID[i] " << bPlasTDCID[i]  <<endl;

          ///Calibrate raw TDC
        //  bPlasTDC_T_Calib[bPlasTDCIDMain] =  bPlasTDC_TS_Raw[bPlasTDCIDMain] + fCal->AplasTDC_Raw[bPlasTDCIDMain];
          pOutput->pbPlasTDC_T[bPlasTDCIDMain] =   bPlasTDC_TS_Raw[bPlasTDCIDMain]; //Output bPlas Raw TDC
            ///Get the first hit of the reference channel (Cha.1)
           if(bPlasTDC_T_Calib[1]>0 && bPlas_TDC_Multiplicity[bPlasTDCIDMain]==1 &&bPlasTDCIDMain ==1){
                           bPlas_TDC_Cha1 =   bPlasTDC_T_Calib[1];
                                   }

          //hPLAS_TDCCalib1[bPlasTDCIDMain] -> Fill(bPlasTDC_T_Calib[bPlasTDCIDMain]);
          hPLAS_TDCCalib1Sum -> Fill(bPlasTDC_T_Calib[bPlasTDCIDMain]);

          bPlasTDC_TS_Firedratio[bPlasTDCIDMain] = bPlasTDC_T_Calib[bPlasTDCIDMain]/bPlasTDCFired;
          hPLAS_TDC_FiredRatio[bPlasTDCIDMain] -> Fill(bPlasTDC_TS_Firedratio[bPlasTDCIDMain]);
          hPLAS_TDC_FiredRatio_Sum -> Fill( bPlasTDC_TS_Firedratio[bPlasTDCIDMain]);

          ///Ref SiPMCh.1 - SiPM Ch.x
       //   if(bPlas_TDC_Multiplicity[bPlasTDCIDMain]>1 && bPlasTDCIDMain==1){

       //   }
          if(bPlas_TDC_Cha1>0 &&bPlasTDC_TS_Raw[bPlasTDCIDMain]>0)
        //  cout<<"event " << event_number<<" bPlas_TDC_Cha1 " <<bPlas_TDC_Cha1 << " bPlasTDC_TS_Raw[bPlasTDCIDMain] " << bPlasTDC_TS_Raw[bPlasTDCIDMain]<< endl;
           bPlas_SiPM_dT_Raw[bPlasTDCIDMain] = (bPlas_TDC_Cha1 - bPlasTDC_TS_Raw[bPlasTDCIDMain]);
           bPlas_SiPM_dT_Calib[bPlasTDCIDMain] = ((bPlas_TDC_Cha1 - bPlasTDC_TS_Raw[bPlasTDCIDMain] )+ fCal->AplasTDC_Chref_dT[bPlasTDCIDMain]);
           pOutput ->  pbPlas_SiPM_dT_Calib[bPlasTDCIDMain] =   bPlas_SiPM_dT_Calib[bPlasTDCIDMain];

           if(bPlasTDCIDMain!=1){
            hPLAS_TimeDiffSiPM_Ch_Raw[bPlasTDCIDMain] -> Fill(bPlas_SiPM_dT_Raw[bPlasTDCIDMain]);
            hPLAS_TimeDiffSiPM_Ch_Calib[bPlasTDCIDMain] -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
            hPLAS_TimeDiffSiPM_Ch_Sum -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);

                      ///Fill for Multiplicites
if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==1) hPLAS_TimeDiffSiPM_Ch_Sum_M1 ->Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==2) hPLAS_TimeDiffSiPM_Ch_Sum_M2 ->Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] >2)  hPLAS_TimeDiffSiPM_Ch_Sum_M3P ->Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
        }
          ///SC41 - SiPM Ch.x
          if(bPlasTDC_TS_Raw[0]>0 && bPlasTDC_TS_Raw[bPlasTDCIDMain]>0){

           bPlas_SC41_dT_Raw[bPlasTDCIDMain] = (bPlasTDC_TS_Raw[0] - bPlasTDC_TS_Raw[bPlasTDCIDMain]);
           bPlas_SC41_dT_Calib[bPlasTDCIDMain] =  bPlas_SC41_dT_Raw[bPlasTDCIDMain] + fCal->BplasTDC_SC41dT[bPlasTDCIDMain];

            pOutput-> pbPlas_SC41_dT[bPlasTDCIDMain] = bPlas_SC41_dT_Calib[bPlasTDCIDMain]; //Output SC41-bPlas dT
            if(bPlasTDCIDMain!=0){
                hPLAS_TimeDiff_SC41_Raw[bPlasTDCIDMain] ->Fill(bPlas_SC41_dT_Raw[bPlasTDCIDMain]);
                hPLAS_TimeDiff_SC41_Calib[bPlasTDCIDMain] ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
                hPLAS_TimeDiff_SC41_Sum ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
                ///Fill for Multiplicites
if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==1) hPLAS_TimeDiff_SC41_Sum_M1 ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] ==2) hPLAS_TimeDiff_SC41_Sum_M2 ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
if (bPlas_TDC_Multiplicity[bPlasTDCIDMain] >2)  hPLAS_TimeDiff_SC41_Sum_M3P ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);

            }
          }
  /**-----------------------bPlastic Energy gated Timing (gates defined in Correlations.dat)-----------------------------------------**/
           /// NOTE: TDC ID = QDC ID +1 (TDC Ch.0 =SC41)
            for (int j=0; j<bPlasQDCFired; j++){

                ///Energy-Time matrices
                if(bPlasQDCID[j] == bPlasTDCIDMain+1 && bPlasTDC_T_Calib[8]>0){
               hPLAS_CoincE_dTSiPM_Ch[7] -> Fill(bPlasQDCGainMatch_i[7], bPlasTDC_T_Calib[8]);
               hPLAS_CoincE_dTSiPM_Sum -> Fill(bPlasQDCGainMatch_i[bPlasQDCID[j]], bPlasTDC_T_Calib[bPlasTDCIDMain+1]);

                ///Energy Gate
     if( bPlasQDCGainMatch_i[bPlasQDCID[j]] > fCorrel->GbPlas_Egate_low &&  bPlasQDCGainMatch_i[bPlasQDCID[j]] < fCorrel->GbPlas_Egate_high ){
                ///Energy gated SiPM - SiPM x.
          if(bPlasTDCIDMain!=1){
           hPLAS_TimeDiffSiPM_Ch_Calib_Egated[bPlasTDCIDMain] -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
           hPLAS_TimeDiffSiPM_Ch_Sum_Egated -> Fill(bPlas_SiPM_dT_Calib[bPlasTDCIDMain]);
                 }

          if(bPlasTDC_TS_Raw[0]>0 && bPlasTDC_TS_Raw[bPlasTDCIDMain]>0){
                ///Energy gated SC41 - SiPM Ch.x
            if(bPlasTDCIDMain!=0){
           hPLAS_TimeDiff_SC41_Calib_Egated[bPlasTDCIDMain] ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
           hPLAS_TimeDiff_SC41_Sum_Egated ->Fill(bPlas_SC41_dT_Calib[bPlasTDCIDMain]);
              }
            }
          }
        }
      }
  }
}

 /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------------  FATIMA  ----------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::FAT_det_pos_setup(){

    FAT_positions   = new double*[36];
    FAT_neighbour_check = new bool*[36];
    FAT_angle_diffs = new double*[36];

    for(int i = 0; i < 36; ++i){
    FAT_positions[i] = new double[3];
    FAT_angle_diffs[i] = new double[36];
    FAT_neighbour_check[i] = new bool[36];
    for (int j = 0; j < 3; ++j) FAT_positions[i][j] = -1;
    for (int k = 0; k < 36; ++k){

        FAT_neighbour_check[i][k] = true;

        FAT_angle_diffs[i][k] = -1;

            }
    }

    const char* format = "%d %lf %lf %lf";

    ifstream file("Configuration_Files/FATIMA_Detector_Positions.txt");

    if(file.fail()){
        cerr << "Could not find FATIMA Detector Positions File!" << endl;
        exit(0);
    }

    string line;
    int pos_num;
    double r, theta, phi;

    while(file.good()){
        getline(file,line,'\n');
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format, &pos_num, &r, &theta, &phi);

        FAT_positions[pos_num][0] = r;
        FAT_positions[pos_num][1] = theta;
        FAT_positions[pos_num][2] = phi;

    }

    if(FAT_nearest_neighbour_exclusion){

    for(int i = 0; i < 36; ++i){

        if(i%12 == 11) FAT_neighbour_check[i][(i-11)] = false; // Same Ring Rignt
        else FAT_neighbour_check[i][(i+1)] = false; // Same Ring Right
        if(i%12 == 0) FAT_neighbour_check[i][i+11] = false; // Same Ring Left
        else FAT_neighbour_check[i][(i-1)] = false; // Same Ring Left

        if(!same_ring_exclusion){

        if(i < 12){

            FAT_neighbour_check[i][(i+12)] = false; // Middle Ring Left

            if(i == 12) FAT_neighbour_check[i][(i+1)] = false; // Middle Ring Left for 11

            else FAT_neighbour_check[i][(i+13)] = false; // Middle Ring Right
        }
        if(i > 11 && i < 24){

            FAT_neighbour_check[i][(i+12)] = false; // Upper Outer Ring
            FAT_neighbour_check[i][(i-12)] = false; // Lower Outer Ring

            if(i == 12){
             FAT_neighbour_check[i][(i+23)] = false; // Upper Outer Ring
             FAT_neighbour_check[i][(i-1)] = false; // Lower Outer Ring
            }
            else{
            FAT_neighbour_check[i][(i+11)] = false; // Upper Outer Ring
            FAT_neighbour_check[i][(i-13)] = false; // Lower Outer Ring
            }
        }
        if(i > 23){

            FAT_neighbour_check[i][(i-12)] = false; // Middle Ring Left

            if(i == 35) FAT_neighbour_check[i][(i-23)] = false; // Middle Ring Left for 35

            else FAT_neighbour_check[i][(i-11)] = false; // Middle Ring Right
                }
            }
        }
    }

    ofstream output_position_matrix_file;
    output_position_matrix_file.open ("Configuration_Files/FATIMA_Exclusion_Matrix.txt");
    cout<<endl;
    cout << "============================================================" << endl;
    cout << "A Matrix of excluded detector pairings can be found in" << endl;
    cout << "'Configuration_Files/FATIMA_Exclusion_Matrix.txt'"<<endl;
    cout << "============================================================" << endl;
    cout<<endl;

    if (output_position_matrix) output_position_matrix_file <<"        "<<"0 "<<"1 "<<"2 "<<"3 "<<"4 "<<"5 "<<"6 "<<"7 "<<"8 "
            <<"9 "<<"10 "<<"11 "<<"12 "<<"13 "<<"14 "<<"15 "<<"16 "<<"17 "
            <<"18 "<<"19 "<<"20 "<<"21 "<<"22 "<<"23 "<<"24 "<<"25 "<<"26 "
            <<"27 "<<"28 "<<"29 "<<"30 "<<"31 "<<"32 "<<"33 "<<"34 "<<"35 "<<endl;

    for(int i = 0; i < 36; ++i){

    if (i >= 10 && output_position_matrix) output_position_matrix_file <<"Det "<<i<<": ";
    if (i < 10  && output_position_matrix) output_position_matrix_file <<"Det "<<i<<" : ";

    for (int k = 0; k < 36; ++k){

        if(k > 9 && output_position_matrix) output_position_matrix_file<<" ";

        double dist = distance_between_detectors( FAT_positions[i][0],  FAT_positions[i][1],  FAT_positions[i][2],
                              FAT_positions[k][0],  FAT_positions[k][1],  FAT_positions[k][2]);

        double angle = angle_between_detectors(FAT_positions[i][0], FAT_positions[k][0], dist);

        FAT_angle_diffs[i][k] = angle;

        if((dist < FAT_exclusion_dist && (((i < 12 && k < 12) ||
                        (i < 24 && i > 11 && k < 24 && k > 11) ||
                        (i > 23 && k > 23)) || !same_ring_exclusion )) || i == k ){


         FAT_neighbour_check[i][k] = false;

        }


        if (output_position_matrix && !FAT_neighbour_check[i][k]) output_position_matrix_file<<"X ";

        else if(output_position_matrix && FAT_neighbour_check[i][k]) output_position_matrix_file<<"0 ";

    }

    if (output_position_matrix) output_position_matrix_file<<endl;

    }

    output_position_matrix_file.close();

}
//-----------------------------------------------------------------------------------------------------------------------------//

double EventAnlProc::distance_between_detectors(double _r, double _theta, double _phi, double r_, double theta_, double phi_){

    _theta = _theta * M_PI/180.0;
    theta_ = theta_ * M_PI/180.0;

    _phi = _phi * M_PI/180.0;
    phi_ = phi_ * M_PI/180.0;

    double dist = sqrt(_r*_r + r_*r_ - 2.0*_r*r_*(sin(_theta)*sin(theta_)*cos(_phi - phi_) + cos(_theta)*cos(theta_)));

    return dist;


}
//-----------------------------------------------------------------------------------------------------------------------------//
double EventAnlProc::angle_between_detectors(double _r, double r_, double dist_){


    double angle_diff = acos((_r*_r + r_*r_ - dist_*dist_)/(2.0*_r*r_));

    angle_diff = angle_diff * 180.0/M_PI;

    return angle_diff;

}
//-----------------------------------------------------------------------------------------------------------------------------//

void EventAnlProc::Make_Fatima_Histos(){


  for (int i=0; i<50; i++){
    hFAT_QDCCalib1[i] =  MakeTH1('D', Form("FATIMA/Energy/EnergyCalib/LaBr_ECalib_Ch.%2d",i), Form("QDC Calib Ch. %2d",i), 4000,0,4000);
    hFAT_QDCdt[i]   = MakeTH1('D', Form("FATIMA/Timing/QDCdt/QDCdt%2d",i), Form("QDCdT Ch.%2d",i), 3201,-40,40);
    //hFAT_TDCCalib1[i] =  MakeTH1('D', Form("FATIMA/Timing/TDCCalib/LaBr_Tcalib%2d",i), Form("TDC channel Calib %2d",i), 1E5,0,2E5);
    hFAT_TDCdt_refSC41[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC/TDCdT_SC41_LaBr%02d", i), Form("TDC dtSC41 All Multip SC41- LaBr%02d", i),4000,-1000,1000);
//     hFAT_TDCdt_refSC41_M1[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC_M1/TDCdT_SC41_M1_LaBr%02d", i), Form("TDC dtSC41 Multip 1 SC41- LaBr%02d", i),25000,-50000,50000);
//     hFAT_TDCdt_refSC41_M2[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC_M2/TDCdT_SC41_M2_LaBr%02d", i), Form("TDC dtSC41 Multip 2 SC41- LaBr%02d", i),25000,-50000,50000);
//     hFAT_TDCdt_refSC41_M3[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC_M3/TDCdT_SC41_M2+_LaBr%02d", i), Form("TDC dtSC41 Multip 2+ SC41- LaBr%02d", i),25000,-50000,50000);
//
    hFAT_TDCdt_refSC41_gated[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_SC41-FatTDC_EGated/TDCdT_Egated_SC41_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4000,-1000,1000);
/*    hFAT_TDCdt_refSC41_M1_gated[i] = MakeTH1('D', Form("FATIMA/Timing/EGated/TDCdt_SC41-FatTDC/TDCdT_Egated_SC41_M1_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4000,-1000,1000);
    hFAT_TDCdt_refSC41_M2_gated[i] = MakeTH1('D', Form("FATIMA/Timing/EGated/TDCdt_SC41-FatTDC/TDCdT_Egated_SC41_M2_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4000,-1000,1000);
    hFAT_TDCdt_refSC41_M3_gated[i] = MakeTH1('D', Form("FATIMA/Timing/EGated/TDCdt_SC41-FatTDC/TDCdT_Egated_SC41_M2+_LaBr%02d", i), Form("TDC Gamma gated dtSC41 SC41- LaBr%02d", i),4000,-1000,1000);
   */

    hFAT_TDCdt_refCha[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_AllM/TDCdT_Cha_AllM_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel All Multip LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
    hFAT_TDCdt_refCha_gated[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_EGated/TDCdT_Cha_EGated_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel All Multip LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
    //     hFAT_TDCdt_refCha_M1[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_M1/TDCdT_Cha_M1_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel Multip 1 LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
//     hFAT_TDCdt_refCha_M2[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_M2/TDCdT_Cha_M2_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel Multip 2 LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
//     hFAT_TDCdt_refCha_M3[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdT_TDC0-TDC_M2+/TDCdT_Cha_M2+_LaBr%02d_LaBr%02d", 0, i), Form("TDC dt Channel Multip 2+ LaBr%02d - LaBr%02d",0 , i),4000,-1000,1000);
//
    hFAT_TDC_Multipl_ch[i] = MakeTH1('D', Form("FATIMA/Stats/TDC_MultiplCh/TDCM_Ch_LaBr%2d",i), Form("TDC channel Multi Fatima %2d",i), 50, 0, 50);

  }
    hFAT_QDC_vs_TDC_SiPMdT_Ch[7] = MakeTH2('D',Form("FATIMA/Timing/Energy_vs._Time_SiPMdT_Ch/Energy_vs._Time_SiPMdT_Ch.%02d", 7),Form("Fatima Energy vs SiPMCh.0-SiPMCh.%02d", 7),4000,0,4000, 4000,-1000,1000);
    hFAT_QDC_vs_TDC_SC41dT_Ch[7] = MakeTH2('D',Form("FATIMA/Timing/Energy_vs._Time_SC41_Ch/Energy_vs._Time_SC41_Ch.%02d", 7),Form("Fatima Energy vs SC41-SiPMCh.%02d", 7),4000,0,4000, 4000,-1000,1000);

    hFAT_QDCCalib1Sum = MakeTH1('D', "FATIMA/Energy/Fat_EnergySum", "LaBr Energy (all detectors)",4000,0,4000);
    hFAT_hits_QDC       = MakeTH1('D', "FATIMA/Stats/QDC_FAThits", "bPlastic hit pattern QDC1",50,0,50);
    hFAT_E_Mat_Sum = MakeTH2('D', "FATIMA/Energy/Gam-GamSum", "FATIMA Gamma-Gamma (all detectors)",4000,0,4000, 4000,0,4000);
    hFAT_hits_TDC       = MakeTH1('D', "FATIMA/Stats/TDC_FAThits", "FATIMA TDC statistics",50,0,50);
    hFAT_TDC_Multipl_PerChan       = MakeTH1('D', "FATIMA/Stats/TDC_FAT_Multiplicity_perCh", "FATIMA TDC Multiplicity (hits per channel)",50,0,50);
    hFAT_TDC_Multipl       = MakeTH1('D', "FATIMA/Stats/TDC_FAT_Multiplicity", "FATIMA TDC Multiplicity",50,0,50);

    hFAT_TDCdt_refSC41_Sum       = MakeTH1('D', "FATIMA/Timing/TDCdt_refSC41_Sum", "TDC dT Ref SC41(all detectors)",10000,-1000,1000);
    hFAT_TDCdt_refSC41_Sum_gated       = MakeTH1('D', "FATIMA/Timing/TDCdt_refSC41_Sum_EGated", "TDC dT (all detectors) Energy gated", 10000,-1000,1000);
    hFAT_TDCdt_refCha_Sum       = MakeTH1('D', "FATIMA/Timing/TDCdt_ref0_AllM_Sum", "TDC dT LaBr0 - LaBr Multip All (all detectors)", 4000,-1000,1000);
    hFAT_TDCdt_refCha_Sum_M1       = MakeTH1('D', "FATIMA/Timing/TDCdt_ref0_M1_Sum", "TDC dT LaBr0 - LaBr Multip 1 (all detectors)", 4000,-1000,1000);
    hFAT_TDCdt_refCha_Sum_M2       = MakeTH1('D', "FATIMA/Timing/TDCdt_ref0_M2_Sum", "TDC dT LaBr0 - LaBr Multip 2 (all detectors)", 4000,-1000,1000);
    hFAT_TDCdt_refCha_Sum_M3       = MakeTH1('D', "FATIMA/Timing/TDCdt_ref0_M2+_Sum", "TDC dT LaBr0 - LaBr Multip 2+ (all detectors)", 4000,-1000,1000);

    hFAT_TDCdt_refCha_Sum_gated     = MakeTH1('D', "FATIMA/Timing/TDCdt_ref0_Sum_EGated","TDC dT LaBr0 Gamma gated (all detectors)",4000,-1000,1000);
    hFAT_QDC_vs_TDC_SiPMdT = MakeTH2('D',"FATIMA/Energy_vs._Time_SiPMdT","Energy_vs._Time_SiPMdT",4000,0,4000, 1000,-1000,1000);
    hFAT_QDC_vs_TDC_SC41dT = MakeTH2('D',"FATIMA/Energy_vs._Time_S41dT","Energy_vs._Time_S41dT",4000,0,4000,  1000,-1000,1000);

}

void EventAnlProc::Do_Fatima_Histos(EventAnlStore* pOutput){
    double Fat_QDC_i[50], Fat_QDC_j[50];
    double Fat_QDC_GainMatch[50], Fat_QDCGainMatch_j[50];
    double FATgate1_low, FATgate1_high;
    double Fat_TDC_T_Main[50], Fat_SC41_dT_Raw[50], Fat_SC41_dT_Calib[50],  Fat_Ch_dT[50], Fat_Ch_dT_Calib[50];
    int Fat_QDC_IDMain_i, Fat_QDC_IDMain_j, Fat_TDC_IDMain;
    double Fat_QDC_dt, Fat_QDCtime1, Fat_QDCtime2;
    int Fat_TDC_Incr;
    int  Fat_TDC_Multipl_perCh[50] ={0};
    

    Fat_QDC_IDMain_i = -1;
    Fat_QDC_IDMain_j = -1;
    Fat_TDC_IDMain = -1;
    Fat_QDC_dt = 0;
    Fat_QDCtime1 = 0;
    Fat_QDCtime2 = 0;
  

    for(int i=0; i<50; i++){
        Fat_QDC_i[i] = -1;
        Fat_QDC_j[i] = -1;
        Fat_QDC_GainMatch[i] = 0;
        Fat_QDCGainMatch_j[i] = 0;
        Fat_Ch_dT_Calib[i] = 0;


        Fat_TDC_T_Main[i] = 0;
        Fat_SC41_dT_Raw[i] = 0;
        Fat_SC41_dT_Calib[i] = 0;
        Fat_Ch_dT[i] = 0;
  }
      //Fatima Energy gates
    FATgate1_low  = fCorrel->GFat_Egate_low;
    FATgate1_high = fCorrel->GFat_Egate_high;
   // Fat_E_gate1 = FATgate1_low + (FATgate1_high - FATgate1_low)/2.;

    /**------------------------------FATIMA Energy -----------------------------------------**/
        pOutput->pFat_QDCFired = FatQDCFired;
        pOutput->pFat_WR = Fat_WR;
      
    for (int i=0; i<FatQDCFired; i++){
     
        Fat_QDC_IDMain_i = FatQDCID[i]; //Channel ID
        pOutput->pFat_QDCID[i] = FatQDCID[i];
        hFAT_hits_QDC->Fill(Fat_QDC_IDMain_i);

        Fat_QDC_i[Fat_QDC_IDMain_i] = FatQDC[i];  //Raw energy
        Fat_QDCtime1 = FatQDC_T[i];

          ///FATIMA Calibrated Energy Singles
        //Fat_QDC_GainMatch[Fat_QDC_IDMain_i] = fCal->Afat[Fat_QDC_IDMain_i]* pow(Fat_QDC_i[Fat_QDC_IDMain_i],3) + fCal->Bfat[Fat_QDC_IDMain_i]* pow(FatQDC[i],2) + fCal->Cfat[Fat_QDC_IDMain_i]*FatQDC[i] + fCal->Dfat[Fat_QDC_IDMain_i];
        Fat_QDC_GainMatch[Fat_QDC_IDMain_i] = Fat_QDC_i[Fat_QDC_IDMain_i];
        hFAT_QDCCalib1[Fat_QDC_IDMain_i]->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i]);
        pOutput->pFat_QDCGainMatch[Fat_QDC_IDMain_i] = Fat_QDC_GainMatch[Fat_QDC_IDMain_i];
        hFAT_QDCCalib1Sum->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i]);

               ///Gamma-Gamma Fatima
          for (int j=0; j<FatQDCFired; j++){
            Fat_QDC_IDMain_j= FatQDCID[j];
            Fat_QDC_j[Fat_QDC_IDMain_j] = FatQDC[j];
            Fat_QDCtime2 = FatQDC_T[j];

            ///Dont loop on the first hit again:
            if(Fat_QDC_IDMain_i < Fat_QDC_IDMain_j){
                Fat_QDC_j[Fat_QDC_IDMain_j] = FatQDC[j];
               // Fat_QDCGainMatch_j[Fat_QDC_IDMain_j] = fCal->Afat[Fat_QDC_IDMain_j]* pow(Fat_QDC_j[Fat_QDC_IDMain_j],3) + fCal->Bfat[Fat_QDC_IDMain_j]* pow(Fat_QDC_j[j],2) + fCal->Cfat[Fat_QDC_IDMain_j]*Fat_QDC_j[j] + fCal->Dfat[Fat_QDC_IDMain_j];
                Fat_QDC_dt = Fat_QDCtime1 - Fat_QDCtime2;
                hFAT_QDCdt[Fat_QDC_IDMain_i] ->Fill(Fat_QDC_dt);
                //Fill Energy-Energy matrix (NOTE: turned off making 2D matrices for each channel for now to speed things up)
               // hFAT_Chan_E_Mat[Fat_QDC_IDMain_i]->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i], Fat_QDCGainMatch_j[Fat_QDC_IDMain_j]);
                hFAT_E_Mat_Sum->Fill(Fat_QDC_GainMatch[Fat_QDC_IDMain_i], Fat_QDCGainMatch_j[Fat_QDC_IDMain_j]);
          }
       }
  }

 /**---------------------------------FATIMA TIMING -----------------------------------------**/
          hFAT_TDC_Multipl->Fill(FatTDCFired);
          pOutput -> pFat_TDCFired = FatTDCFired;
         for (int i=0; i<FatTDCFired; i++){

           ///FAT TDC ID and Raw TDC data
          Fat_TDC_IDMain = FatTDCID[i];
          pOutput -> pFat_TDCID[i] = FatTDCID[i];
          Fat_TDC_T_Main[Fat_TDC_IDMain] = FatTDC_TS[i][Fat_TDC_IDMain];
          pOutput ->  pFat_TDC_T[Fat_TDC_IDMain] =  Fat_TDC_T_Main[Fat_TDC_IDMain];
     
          ///FAT Multiplicity
           Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]++;
           pOutput -> pFat_TDC_Multipl_perCh[Fat_TDC_IDMain] =  Fat_TDC_Multipl_perCh[Fat_TDC_IDMain];
           hFAT_TDC_Multipl_ch[i] -> Fill(Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]);
           hFAT_TDC_Multipl_PerChan -> Fill(Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]);

         

           ///Hit pattern
           hFAT_hits_TDC->Fill(Fat_TDC_IDMain);
            ///SC41 (TDC Ch.40) - FAT SiPM Ch.x
           //if( Fat_TDC_T_Main[Fat_TDC_IDMain]>fCorrel->GFat_TRawgate_low && Fat_TDC_T_Main[Fat_TDC_IDMain]<fCorrel->GFat_TRawgate_high){
            Fat_SC41_dT_Raw[Fat_TDC_IDMain] = (SC41_ns -  Fat_TDC_T_Main[Fat_TDC_IDMain]);
            Fat_SC41_dT_Calib[Fat_TDC_IDMain]  = Fat_SC41_dT_Raw[Fat_TDC_IDMain] + fCal-> TFatTDC_SC41dT[Fat_TDC_IDMain];
            pOutput ->  pFat_SC41_dT_Calib[Fat_TDC_IDMain] =  Fat_SC41_dT_Calib[Fat_TDC_IDMain];

            ///SC41 - Fatima TDC Only take the first hit per channel (Sultan)
           Fat_TDC_Incr = 0;
            for(int j=0; j<=i; j++){
             if (Fat_TDC_IDMain != FatTDCID[j] ) Fat_TDC_Incr++;
            }
            //if(Fat_TDC_Incr == i  && Fat_TDC_IDMain < 40 ){
             for (int j = 0; j< FatQDCFired; j++){
            if( FatQDCID[j] == Fat_TDC_IDMain){

              hFAT_TDCdt_refSC41[Fat_TDC_IDMain] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
              hFAT_TDCdt_refSC41_Sum -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
              //
              //pOutput->pFat_SC41_dT_Calib[Fat_TDC_IDMain] = Fat_SC41_dT_Calib[Fat_TDC_IDMain];
            }
      
              /// Fatima Time SiPM 0 - SiPM Ch.x (Ch. 0 used as the reference)
              if(Fat_TDC_IDMain < 40 && Fat_CHA_0_TDC>0  && Fat_TDC_T_Main[Fat_TDC_IDMain] > 0&& FatQDCID[j] == Fat_TDC_IDMain){
                //cout << "evn " << event_number <<" Fat_CHA_0_TDC " << Fat_CHA_0_TDC << " Fat_TDC_IDMain " << Fat_TDC_IDMain <<"Fat_TDC_T_Main[0] " << Fat_TDC_T_Main[0] << " Fat_TDC_T_Main[Fat_TDC_IDMain] " << Fat_TDC_T_Main <<" multi " << Fat_TDC_Multipl_perCh[Fat_TDC_IDMain] << endl;

                    Fat_Ch_dT[Fat_TDC_IDMain] =  (Fat_CHA_0_TDC - Fat_TDC_T_Main[Fat_TDC_IDMain]);
                 //   if(Fat_QDC_GainMatch[FatQDCID[j]] > FATgate1_low && Fat_QDC_GainMatch[FatQDCID[j]] < FATgate1_high){

                 //   if(Fat_Ch_dT[Fat_TDC_IDMain]>114&&Fat_Ch_dT[Fat_TDC_IDMain]<139) hFAT_test[Fat_TDC_IDMain]->Fill(Fat_TDC_T_Main[Fat_TDC_IDMain]);
                    Fat_Ch_dT_Calib[Fat_TDC_IDMain] =  Fat_Ch_dT[Fat_TDC_IDMain] + fCal-> TFatTDC_Chref_dT[Fat_TDC_IDMain];
                    pOutput->pFat_Ch_dT[Fat_TDC_IDMain] =  Fat_Ch_dT_Calib[Fat_TDC_IDMain] ;
                    if(Fat_TDC_IDMain!=0){
                        hFAT_TDCdt_refCha[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                        hFAT_TDCdt_refCha_Sum  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);

                    /// Histogram TDC cha 0 - TDC Ch.x
                    if(Fat_TDC_IDMain!=0 && Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==1 ){
                    //hFAT_TDCdt_refCha_M1[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                    hFAT_TDCdt_refCha_Sum_M1  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                        }
                    if(Fat_TDC_IDMain!=0 && Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==2 ){
                   // hFAT_TDCdt_refCha_M2[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                    hFAT_TDCdt_refCha_Sum_M2  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                        }
                    if(Fat_TDC_IDMain!=0 && Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]>2 ){
                   // hFAT_TDCdt_refCha_M3[Fat_TDC_IDMain]->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                    hFAT_TDCdt_refCha_Sum_M3  ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                            }
                        }
                   // }
              }
                        ///Energy Time matrix (just one channel (ch. 7) for now to speed things up)

                       if(FatQDCID[j] == Fat_TDC_IDMain && Fat_QDC_GainMatch[FatQDCID[i]]>0  ){
                    hFAT_QDC_vs_TDC_SiPMdT_Ch[7] ->Fill( Fat_QDC_GainMatch[7],Fat_Ch_dT[Fat_TDC_IDMain]);
                    hFAT_QDC_vs_TDC_SC41dT_Ch[7] ->Fill( Fat_QDC_GainMatch[7],Fat_SC41_dT_Calib[Fat_TDC_IDMain]);

                    hFAT_QDC_vs_TDC_SiPMdT ->Fill( Fat_QDC_GainMatch[FatQDCID[j]],Fat_Ch_dT[Fat_TDC_IDMain]);
                    hFAT_QDC_vs_TDC_SC41dT ->Fill( Fat_QDC_GainMatch[FatQDCID[j]],Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
                            }


            ///Gamma energy gates
             if(Fat_QDC_GainMatch[FatQDCID[j]] > FATgate1_low && Fat_QDC_GainMatch[FatQDCID[j]] < FATgate1_high){
                 /// Fatima Time SiPM 0 - SiPM Ch.x Energy gated
                if(Fat_TDC_IDMain!=0) {
                   hFAT_TDCdt_refCha_gated[Fat_TDC_IDMain] ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
                   hFAT_TDCdt_refCha_Sum_gated ->Fill(Fat_Ch_dT[Fat_TDC_IDMain]);
            ///SC41 - Fatima TDC Energy Gated
            hFAT_TDCdt_refSC41_gated[Fat_TDC_IDMain] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
         //  if (Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==1) hFAT_TDCdt_refSC41_M1_gated[i] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
         //    if (Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==2) hFAT_TDCdt_refSC41_M2_gated[i] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
          //   if (Fat_TDC_Multipl_perCh[Fat_TDC_IDMain]==3) hFAT_TDCdt_refSC41_M3_gated[i] -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);
            hFAT_TDCdt_refSC41_Sum_gated -> Fill(Fat_SC41_dT_Calib[Fat_TDC_IDMain]);

                   }
                }
             }
           //}
        }
    }

/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  GALILEO  ---------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
 void EventAnlProc::Make_Galileo_Histos(){
    hGAL_ESum = MakeTH1('I',"GALILEO/Sum/GALILEO_ESum","GALILEO Energy Sum",5000,0,5000);
    hGAL_Hit_Pat = MakeTH1('I',"GALILEO/Stats/GALILEO_Hit_Pat","GALILEO Hit Pattern",32,0,32);
    hGAL_Multi = MakeTH1('I',"GALILEO/Stats/GALILEO_Multiplicity","GALILEO Multiplicity",50,0,50);
    hGAL_Chan_E_Mat = MakeTH2('D',"GALILEO/GALILEO_E_Mat","GALILEO Energy-Energy Matrix",5001,0,10000,5001,0,10000);
    hGAL_Chan_E_M1= MakeTH1('I',"GALILEO/Stats/GALILEO_multiplicity_1","GALILEO Channel Energy",5000,0,5000);
    hGAL_Chan_E_M2= MakeTH1('I',"GALILEO/Stats/GALILEO_multiplicity_2","GALILEO Channel Energy",5000,0,5000);
    hGAL_AddbackSum = MakeTH1('I',"GALILEO/Sum/GALILEO_Addback","GALILEO Addback Energy Sum",5000,0,5000);
    
    for (int j=0; j<32; j++)
    {
    hGAL_Chan_E[j] = MakeTH1('D',Form("GALILEO/GALILEO_Energy/GALILEO_E%2d",j), Form("GALILEO Channel Energy Channel %2d",j),5000,0,5000);
    hGAL_FatdT[j] = MakeTH1('I',Form("Correlations/Fatima_Galilieo/Fat_GAldT%2d",j),Form("GALILEO Fatima dT Ch. %2d",j),2000,-1000,1000);
   // hGAL_Chan_E2[j] = MakeTH1('D',Form("GALILEO/GALILEO_Energy2/GALILEO_E2%2d",j), Form("GALILEO Channel Energy Channel %2d",j),5000,0,5000);
    //hGAL_Chan_Egate[j] = MakeTH1('D',Form("GALILEO/gated energy/GALILEO_Egate%2d",j), Form("GALILEO Channel Energy Channel %2d",j),5000,0,5000);
    }
    for (int k=0; k<32; k++){
    
    hGAL_Chan_Time_Diff[k] = MakeTH1('D',Form("GALILEO/Time_diff/GALILEO_Chan_Time_Diff%2d",k), Form("GALILEO Channel Time Difference for %2d",k),100,-1000,1000);
     
    hGAL_Time_Diff_vs_Energy[k] = MakeTH2('D',Form("GALILEO/GALILEO_dT_vs_Energy_Spectra/GALILEO_dT_vs_E%2d",k), Form("GALILEO Time Difference Vs Channel Energy Channel %2d",k),5000,0,5000,100,-1000,1000);
        } 
    }

void EventAnlProc::Do_Galileo_Histos(EventAnlStore* pOutput){
        double  GalE1[32],GalE2_i[32],Gal_time_diff,sumM2;
        double GalE_Cal_i[32],GalE2_Cal_i[32];
        double GalE_Cal_Sum_i,sum1,sum2;
        long GalT_i[32];
        double AddbackSum;
        long Fat_Gal_dT;
        GalE_Cal_Sum_i = 0;
        Gal_time_diff = 0;
        sum1 = 0;
        sum2 = 0;
        sumM2 = 0;
        AddbackSum = 0;
        for (int i =0; i<32; i++){
            GalE_Cal_i[i] = 0;
            GalE2_Cal_i[i] = 0;
            GalE1[i] = 0;
            GalE2_i[i] = 0;
            GalT_i[i] = 0;
        
        }
//         for(int j=0;j<1000;j++){
//             Tdiff_WR_gal[j]=0;
//         }

        //Galileo multiplicity
        hGAL_Multi -> Fill(GalFired);
        pOutput-> pGalFired = GalFired;
        pOutput->pGal_WR = Gal_WR;
        
//         int time_gate=1000;
//         if(Gal_WR>0){
//         for (int p = 0; p<time_gate; p++){
//         Tdif" f_WR_gal[p] = Gal_WR;
       // if (Used_Systems[3] && PrcID_Conv[3]==3 ){
       //     if (Fat_WR>0){
  //  cout << "Gal_WR "<<Gal_WR<< " Fat_WR "<<Fat_WR<<" Fat_WR-Gal_WR " << Fat_WR-Gal_WR <<endl;}}
    for (int i = 0; i < GalFired; i++){
         
        pOutput-> pGalID[i] = GalID[i];
        GalT_i[i]=GalT[GalID[i]]; //Galileo First hit Time
        pOutput-> pGalT[GalID[i]] = GalT[GalID[i]]; 
        GalE1[GalID[i]] = GalE[GalID[i]]; //Galileo raw energy
        Fat_Gal_dT =  (Fat_WR-Gal_WR);
        hGAL_FatdT[GalID[i]]->Fill(Fat_Gal_dT);
        //Gal_Multipl++;
        ///Galileo calibrated energy
        GalE_Cal_i[GalID[i]] = fCal->AGal[GalID[i]]* pow( GalE1[GalID[i]],2) + fCal->BGal[GalID[i]]*  GalE1[GalID[i]] + fCal->CGal[GalID[i]];
        pOutput-> pGalE_Cal_i[GalID[i]] = GalE_Cal_i[GalID[i]];
        ///Galileo Hit pattern
        hGAL_Hit_Pat->Fill(GalID[i]);
        
        ///Galileo energy sum all channels
        GalE_Cal_Sum_i =  GalE_Cal_i[GalID[i]];

        if(GalE_Cal_i[GalID[i]] > 0){
            ///Galileo fill energy spectra
            hGAL_Chan_E[GalID[i]]->Fill(GalE_Cal_i[GalID[i]]);
            hGAL_ESum ->Fill(GalE_Cal_Sum_i);
            AddbackSum += GalE_Cal_Sum_i; 
    }
         // cout<<"1)event " << event_number << " addback="<<AddbackSum<<" GalE_Cal_i[GalID[i]] " << GalE_Cal_i[GalID[i]] << endl;
           ///Galileo Gamma-Gamma, (additional fired crystals)
        for(int k =i+1; k < GalFired; k++){
            GalE2_i[k] = GalE[GalID[k]];
            GalE_Cal_i[GalID[k]] = fCal->AGal[GalID[k]]* pow(GalE2_i[k],2) + fCal->BGal[GalID[k]]* GalE2_i[k] + fCal->CGal[GalID[k]];
            GalT_i[k]=GalT[GalID[k]];
            Gal_time_diff = GalT_i[i] - GalT_i[k];
            hGAL_Chan_Time_Diff[GalID[k]]->Fill(Gal_time_diff);
            pOutput->pGal_dT = Gal_time_diff;
      // cout<<"2)event " << event_number << " addback="<<AddbackSum<<" GalE_Cal_i[GalID[i]] " << GalE_Cal_i[GalID[i]] <<" GalE_Cal_i[GalID[k]] " << GalE_Cal_i[GalID[k]] << endl;

            //Time gate for Multiplicity 2 events
            if (Gal_time_diff>-50 && Gal_time_diff<50)
            {
                sumM2 += GalE_Cal_Sum_i;
                hGAL_Chan_E_M2->Fill(sumM2);
            }

       /// Gamma-Gamma Energy matrix
        hGAL_Chan_E_Mat->Fill(GalE_Cal_i[GalID[i]],GalE_Cal_i[GalID[k]]);
            }
        }
        if(AddbackSum>0){
           hGAL_AddbackSum->Fill(AddbackSum);
           pOutput->pGalE_Addback = AddbackSum;
        }
    }
/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  FINGER  ----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_Finger_Histos(){

  for (int i =0; i<52; i++){
    hFING_lead_lead[i] = MakeTH1('D',Form("FINGER/lead-lead/lead-leadCh.%02d",i),Form("lead-leadCh.%02d",i),1000, -50000, 50000.);
    hFING_ToT[i] = MakeTH1('D',Form("FINGER/TOT/TOTCh%02d",i),Form("TOT Ch(%2d)",i), 2001, -100000., 100000.);
    hFING_trig_lead[i] = MakeTH1('D',Form("FINGER/trig-lead/trig-leadCh.%02d",i),Form("trig-leadCh.%02d",i), 500, -25000., 25000.);
    //hFING_ToT_lead_lead[i] = MakeTH2('D',Form("FINGER/ToT-Trig-Lead/StatCh%02d",i),Form("ToT_LeadCh%02d",i), 500, 0., 50., 200, 300., 400.);
    hFING_MaxToT[i] = MakeTH1('D',Form("FINGER/MAXTOT/MAXTOTCh%02d",i),Form("MAXTOT Ch(%2d)",i), 1000, 0., 10000.);
    hFING_fcoarse[i] = MakeTH1('D',Form("FINGER/PADI_Coarse/CoarseCh%02d",i),Form("Coarse%2d",i), 1E6, 0., 1.2E7);
    hFING_ffine[i] = MakeTH1('D',Form("FINGER/PADI_Fine/FineCh%02d",i),Form("Fine%2d",i), 12625, -500., 50000.);
    hFING_trail_trail[i]= MakeTH1('D',Form("FINGER/trail-trail/trail-trailCh.%02d",i),Form("trail-trailCh.%02d",i),1000, -50000., 50000.);
      
}

  hFING_Hit_Pat = MakeTH1('I',"FINGER/Stats/FINGER_Hit_Pat","FINGER Hit Pattern",52,0,52);
  hFING_ToT_StripID = MakeTH2('I',"FINGER/TOT_vs_PMT","ToT vs Strip number", 2001, -100000., 100000., 52, 0, 52);
  hFING_MaxToT_StripID = MakeTH2('I',"FINGER/MaxTOT_vs_PMT","MaxToT vs Strip number", 2001, -100000., 100000., 52, 0, 52);
  hFING_Pos = MakeTH2('D',"FINGER/position","Time ratio vs Strips",51,1,51, 1000, -10., 10.);
  hFING_Pos_ToT = MakeTH2('D',"FINGER/positionToT","ToT ratio vs Strips",51,1,51, 5000, -1., 1.);
  hFING_Pos_ToT_Max = MakeTH2('D',"FINGER/positionToTMax","ToT ratio vs Strips Max",51,1,51, 5000, -1., 1.);
 
  hFING_ToT_StripID_Exp= MakeTH2('I',"FINGER/TOT_vs_PMT_Exp","ToT exponential vs Strip number", 1000, 0., 10000., 52, 0, 52);
  hFING_MaxToTExp_StripID = MakeTH2('I',"FINGER/MaxTOTExp_vs_PMT","MaxToT exponential vs Strip number", 1000, 0., 100000., 52, 0, 52);

  hFING_ToT_StripID_UpDown = MakeTH2('I',"FINGER/TOT_vs_PMT_sumpmt","ToT vs Strip number sum PMT", 2001, -100000., 100000., 52, 0, 52);
     
  hFING_ToT_StripID_UpDown_max = MakeTH2('I',"FINGER/TOT_vs_PMT_sumpmt_MAX","ToT vs Strip number sum PMT MAX", 2001, -100000., 100000., 52, 0, 52);
}

void EventAnlProc::Do_Finger_Histos(EventAnlStore* pOutput){
        double Fing_LeadDiff[100],Fing_TrailDiff[100],Fing_SC41_diff[100],  Fing_LeadPlus[100];
        double Fing_TOT_Chan[52];
        double FingToT_E[4][32];
        double FingToT_E_Max;
        double Fing_pos_ToT_max[52];
        double Fing_ratio;
        double Fing_ToT_UpDown[4][32];
        double Fing_pos_ToT[4][32];
        //  double upData = 0;                                                          //////E.Sahin    pos vs strip calc
        double downData = 0;
        double total_time = 0;
        //  int pmtUp;
        // int pmtDown;
        int firedPmtNumber = -1;
        int totalFiredPMT = -1;

        FingToT_E_Max = 0;
        for (int i=0; i<100; i++){
            Fing_LeadDiff[i] = 0;
            Fing_TrailDiff[i] = 0;
            Fing_SC41_diff[i] = 0;
            Fing_LeadPlus[i] = 0;   
            }
    
        for(int k=0; k<4;k++){
            for(int l=0;l<32;l++){
            Fing_ToT_UpDown[k][l] = -9999;
            FingToT_E[k][l] = -9999;
            Fing_pos_ToT[k][l] = -9999;
            }      
        }
   for (int k=0; k<52; k++){
        Fing_TOT_Chan[k] = -1;
       
   }
        
  for(int a=0;a<50;a++){
    pmtSetPerEvent[a] = -9999;
    dataSetPerEvent[a] = -9999;
    totaltimeEvent[a] = -9999;
  }
  pOutput -> pFing_firedTamex = Fing_firedTamex;
    for(int i =0; i<Fing_firedTamex; i++){

    //Lead hits
        pOutput -> pFing_iterator[i] = Fing_iterator[i];
    
        for (int j =0; j<Fing_iterator[i]; j++){
            pOutput->pFing_LeadChan[i][j] = Fing_leadChan[i][j];
     // Trail - Trail
            if(Fing_chID[i][j] % 2 == 0){
       //      cout <<"2ev " << event_number << " Fing_trailT[i][j]  " << Fing_trailT[i][j] <<" Fing_trailChan[i][j] "<<Fing_trailChan[i][j] << endl;
                Fing_TrailDiff[Fing_trailChan[i][j]] = (Fing_trailT[i][j] - Fing_trailT[i][j+2]);
                hFING_trail_trail[Fing_trailChan[i][j]] -> Fill( Fing_TrailDiff[Fing_trailChan[i][j]] );
       
        }
        
        
            if(Fing_leadChan[i][j]>-1){
                hFING_fcoarse[Fing_leadChan[i][j]]->Fill(Fing_lead_coarse[i][j]);
                hFING_ffine[ Fing_leadChan[i][j]] ->Fill(Fing_lead_fine[i][j]);
            // cout << "1) ev " << event_number << " Fing_lead_coarse[i][j] " << Fing_lead_coarse[i][j]<< " Fing_lead_fine[i][j] " << Fing_lead_fine[i][j] <<" i " << i << " j "<<j <<endl;
            }
            if(Fing_trailChan[i][j]>-1){
                hFING_fcoarse[Fing_trailChan[i][j]]->Fill(Fing_trail_coarse[i][j]);
                hFING_ffine[ Fing_trailChan[i][j]]   ->Fill(Fing_trail_fine[i][j]);
              // cout << "2) ev " << event_number << " Fing_trail_coarse[i][j] " << Fing_trail_coarse[i][j]<<" Fing_trail_fine[i][j] " << Fing_trail_fine[i][j]<< " i " << i << " j "<<j << endl;
            }
            //cout << "1)event " << event_number << "Fing_TOT_added[i][j] " << Fing_TOT_added[i][j] <<" Fing_leadChan[i][j] "<<Fing_leadChan[i][j]<<endl;
            
    //cout << "2)event " << event_number << "maxToT_added " << maxToT_added <<" maxToT_added_Chan "<<maxToT_added_Chan<<endl;
            if (j<32&& Fing_leadChan[i][j]>-1){

                hFING_Hit_Pat ->Fill(Fing_leadChan[i][j]);
                pOutput -> pFing_leadT[i][j] = Fing_leadT[i][j];
    
            ///Lead - Lead Time
                 if(Fing_leadT[i][j]>0&&Fing_leadT[i][j+2]>0){ 
                if(Fing_leadChan[i][j] % 2 == 0){//Even Strips (for top bottom PMT matching)
                Fing_LeadDiff[Fing_leadChan[i][j]] = (Fing_leadT[i][j] - Fing_leadT[i][j+2]);
                Fing_LeadPlus[Fing_leadChan[i][j]] = (Fing_leadT[i][j] + Fing_leadT[i][j+2]);
                }
                 if(Fing_leadChan[i][j] % 2 == 1){ //Odd Strips
                Fing_LeadDiff[Fing_leadChan[i][j]] = (Fing_leadT[i][j+2] - Fing_leadT[i][j]);
                Fing_LeadPlus[Fing_leadChan[i][j]] = (Fing_leadT[i][j+2] + Fing_leadT[i][j]);
                }
            
                Fing_LeadPlus[Fing_leadChan[i][j]] = (Fing_leadT[i][j] + Fing_leadT[i][j+2]);
                hFING_lead_lead[Fing_leadChan[i][j]] ->Fill(Fing_LeadDiff[Fing_leadChan[i][j]]);
                pOutput -> pFing_LeadDiff[Fing_leadChan[i][j]]  = Fing_LeadDiff[Fing_leadChan[i][j]];
                pOutput -> pFing_LeadPlus[Fing_leadChan[i][j]] =  Fing_LeadPlus[Fing_leadChan[i][j]];
              
            //    if(Fing_leadChan[i][j]==maxToTChan){
                    hFING_LeadLead_StripID->Fill(Fing_LeadDiff[Fing_leadChan[i][j]],Fing_leadChan[i][j]);
               // }
            }

        ///Trigger - Lead
        if(Fing_leadT[i][j]>0 && Fing_leadT[0][0]>0){

          Fing_SC41_diff[Fing_leadChan[i][j]] =  (Fing_leadT[0][0] - Fing_leadT[i][j]);
          pOutput-> pFing_SC41_diff[Fing_leadChan[i][j]] = Fing_SC41_diff[Fing_leadChan[i][j]];
          // cout << "ev " << event_number  << " i "<<i << " j " <<j << " Fing_leadT[0][0] " << Fing_leadT[0][0] <<" Fing_leadT[i][j] " << Fing_leadT[i][j] <<" Fing_SC41_diff " << Fing_SC41_diff[Fing_leadChan[i][j]]<<  " Fing_leadChan[i][j] " << Fing_leadChan[i][j] <<endl;
          hFING_trig_lead[Fing_leadChan[i][j]] ->Fill(Fing_SC41_diff[Fing_leadChan[i][j]]);
        }
        ///Time/Threshold for PMT up and PMT down sum
           if(Fing_leadChan[i][j]>17){
                Fing_ToT_UpDown[i][j] = (Fing_trailT[i][j+1]-Fing_leadT[i][j]) + (Fing_trailT[i][j-1]- Fing_leadT[i][j-2]);
          //  }
//      cout <<"1) event " << event_number << " Fing_trailT[i][j+1] " << Fing_trailT[i][j+1] <<" Fing_leadT[i][j] " << Fing_leadT[i][j] <<  " Fing_trailT[i][j-1] " << Fing_trailT[i][j-1] <<" Fing_leadT[i][j-2] " << Fing_leadT[i][j-2]  << " Fing_ToT_UpDown[i][j] " << Fing_ToT_UpDown[i][j] <<" Fing_leadChan[i][j] " << Fing_leadChan[i][j] << " Fing_tamex_ch[i][j] " << Fing_tamex_ch[i][j]<<  " i " << i << " j " << j <<endl;
            hFING_ToT_StripID_UpDown -> Fill(Fing_ToT_UpDown[i][j],Fing_leadChan[i][j]);
            ///ToT for the Maximum value in a given event
            if((Fing_tamex_ch[i][j] && Fing_tamex_ch[i][j]+1) >0 && maxToT_added_Chan>17){
               hFING_ToT_StripID_UpDown_max -> Fill(maxToT_added,maxToT_added_Chan);
        }
        ///Get position from ToT 
        if(Fing_TOT[i][j]>0  ){
            Fing_pos_ToT_max[maxToT_added_Chan] = (Fing_TOT[i][j]/maxToT_added);
            Fing_pos_ToT[i][j] = (Fing_TOT[i][j]/Fing_ToT_UpDown[i][j]);
     
          if(Fing_pos_ToT[i][j]>0){
            hFING_Pos_ToT -> Fill(Fing_leadChan[i][j], Fing_pos_ToT[i][j]);
            if( maxToT_added_Chan>17){
            hFING_Pos_ToT_Max ->Fill(Fing_leadChan[i][j], Fing_pos_ToT_max[Fing_leadChan[i][j]]);
            }
        }
          //  cout << " event " << event_number <<" Fing_TOT[i][j] " << Fing_TOT[i][j] <<" Fing_ToT_UpDown[i][j] " << Fing_ToT_UpDown[i][j] << " Fing_pos_ToT[i][j] " << Fing_pos_ToT[i][j] << " i " << i << " j " << j <<endl;  
        }
            if(Fing_leadChan[i][j]>0){         
                ///ToT (Only up) to Energy conversion, Tau needs to be checked 
                FingToT_E[i][j] = 2*exp(Fing_pos_ToT[i][j]/1000000);
               //cout <<"FingToT_E " << FingToT_E[i][j] <<" Fing_TOT[i][j] "<< Fing_TOT[i][j]<< " i " << i << " j " << j <<endl;  
                hFING_ToT_StripID_Exp ->Fill( FingToT_E[i][j], Fing_leadChan[i][j]);
                hFING_ToT[Fing_leadChan[i][j]]->Fill(Fing_TOT[i][j]);
                hFING_ToT_StripID ->Fill(Fing_TOT[i][j], Fing_leadChan[i][j]);
                pOutput-> pFing_TOT[i][j] = Fing_TOT[i][j]; 
                               
                    }
           }
        firedPmtNumber = Fing_leadChan[i][j];
        dataSetPerEvent[Fing_leadChan[i][j]] = Fing_LeadDiff[Fing_leadChan[i][j]];
        totaltimeEvent[Fing_leadChan[i][j]] = Fing_LeadPlus[Fing_leadChan[i][j]]/1E3;
        
      if(Fing_leadChan[i][j]>-1){
      downData = dataSetPerEvent[Fing_leadChan[i][j]];
      total_time = totaltimeEvent[Fing_leadChan[i][j]];
//      cout <<"2event " << event_number <<" downData " <<downData << " dataSetPerEvent[fingchan] "<<dataSetPerEvent[Fing_leadChan[i][j]] << " total_time " <<  total_time << " Fing_leadChan[i][j] " << Fing_leadChan[i][j]<<endl;
      
      
      
      ///Needs testing 
      pOutput -> pFing_downData = downData;
      pOutput -> pFing_total_time = total_time;
      ///Lead-lead/lead+lead
      if(Fing_leadChan[i][j]>0){
       Fing_ratio = Fing_LeadDiff[Fing_leadChan[i][j]]/total_time;
      // cout <<"Fing_ratio "<<Fing_ratio <<" Fing_LeadDiff[Fing_leadChan[i][j]] " << Fing_LeadDiff[Fing_leadChan[i][j]] <<" total_time " << total_time << endl;
       hFING_Pos->Fill(Fing_leadChan[i][j], Fing_ratio);
        //pmtSetPerEvent[iter] =  Fing_leadChan[i][j];
       //if (event_number==2157836){
 //cout <<"1) Fing_ratio " <<Fing_ratio << " Fing_LeadDiff[Fing_leadChan[i][j]] "<< Fing_LeadDiff[Fing_leadChan[i][j]]<<" Fing_LeadPlus[Fing_leadChan[i][j]] " << Fing_LeadPlus[Fing_leadChan[i][j]] <<endl;  
        totalFiredPMT++;
            }
        }
    }
    if(maxToT>0){          
                 hFING_MaxToT[maxToTChan]->Fill(maxToT);
                 hFING_MaxToT_StripID->Fill(maxToT, maxToTChan);
                  ///MaxToT Energy conversion 
                 FingToT_E_Max = 2*exp(maxToT/1000000);
                 hFING_MaxToTExp_StripID->Fill(FingToT_E_Max, maxToTChan);
    }
    
  //////E.Sahin June 2019   pos vs strip calc
  //if (event_number==2157836){
  
//     for(int m=0; m<Fing_firedTamex; m++){
//         for(int n=0; n<Fing_iterator[m]; n++ ){
// //
//    // if (st%2==0){
// //      upData = dataSetPerEvent[st];
// 
// //       pmtUp = st;
// //       pmtDown = st-1;
//       //  cout << st << endl;
//   //  }
// 
//     //else{
//     // upData = dataSetPerEvent[st-1];
// //       downData = dataSetPerEvent[st];
// //       total_time = totaltimeEvent[st];
// //       pOutput -> pFing_downData = downData;
// //       pOutput -> pFing_total_time = total_time;
// //      cout  <<"3event " << event_number <<" downData " <<downData <<  " total_time " <<  total_time << " st " << st<<endl;
// 
// //       pmtUp = st-1;
// //       pmtDown = st;
// 
//       // cout << "upData  " << upData << "  downData  " << downData << " st " << st << endl;
// 
//     //}
//   //  if(downData>-100){
// 
// //      double sum = upData + downData;
//      // double ratio = upData / sum;
//     
//       //////////Histogram filling///do not forget ratio
//           }
//         }
      }
     }
   }
   

    /**----------------------------------------------------------------------------------------------**/
    /**--------------------------- FATIMA - Plastic Online Correlations -----------------------------**/
    /**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_Fat_Plas_Histos(){
  for (int i =0; i<50; i++){
 // hFAT_TDCdt_refSC41calib[i] = MakeTH1('D', Form("FATIMA/Timing/TDCdt_refSC41/Calib_ref/TDCdt_LaBr%02d_LaBr%02d_calib", 40, i), Form("TDC dt LaBr%02d LaBr%02d", 40, i),27500,-5000,50000);
        hFat_bplas_Corr_SC41_Dets[i] =  MakeTH1('D', Form("Correlations/FATIMA_PLASTIC_General/Fat_SC41-Plas_SC41_dT/Fat-bPlas_Tdiff_SC41_FatCh.%2d",i), Form("Time difference SC41_Fat_bplas: Fatima Ch. %2d",i), 1000,-1000,1000);
        hFat_bplas_Corr_SC41_Dets_gated[i] =  MakeTH1('D', Form("Correlations/FATIMA_PLASTIC_Gated/Fat_SC41-Plas_SC41_dT_Gated/Fat-bPlas_Tdiff_SC41_Fat_Egate_Ch.%2d",i), Form("Time difference SC41_Fat_bplas gated on fatima energy. Ch. %2d",i), 5000,-10000,10000);
        hFat_bplas_Corr_SiPM_Dets[i] =  MakeTH1('D', Form("Correlations/FATIMA_PLASTIC_General/Fat_Ref-ChdT--Plas_Ref-ChdT/Fat-ChdT_bPlasChdT_FatCh.%2d",i), Form("Time difference: Fatima (Ch.0-Ch.x) - bplas (Ch.1-Ch.x average) Per Fatima Ch. %2d",i), 1000,-1000,1000);
        hFat_bplas_Corr_RawTDCFAT_RawTDCbP_Dets[i] =  MakeTH1('D', Form("Correlations/FATIMA_PLASTIC_General/FatTDC_Ch--PlasTDC_Avg/Fat-TDC_bPlasTDCAvg_FatCh.%2d",i), Form("Time difference: Fatima Raw TDC - bplas TDC Average Per Fatima Ch. %2d",i), 1000,-1000,1000);
        hFat_bplas_Corr_SiPM_Gated_Dets[i] =  MakeTH1('D', Form("Correlations/FATIMA_PLASTIC_Gated/FatTDC_Ch--PlasTDC_Avg/Fat-TDC_bPlasTDCAvg_EGated_FatCh.%2d",i), Form("Time difference: Fatima Raw TDC - bplas TDC Average Fatima EGated Per Fatima Ch. %2d",i), 1000,-1000,1000);

    }
 // hFAT_TDCdt_refSC41calib_sum = MakeTH1('I',"FATIMA/Timing/TDCdt_refSC41/TDCdt_calib_Sum","Average SC41-Plastic",27500,-5000,50000);
  hPLAS_TimeDiff_SC41_avg = MakeTH1('I',"bPlastic/Timing/TDCdT_SC41-SiPM_Avg","SC41 - SiPM dT(ns) Average (calibrated)",2000,-1000,1000);
  hPLAS_TimeDiff_SiPM_avg = MakeTH1('I',"bPlastic/Timing/TDCdT_SiPM1-SiPM_Avg","SiPM1 - SiPM Ch.x dT(ns) Average (calibrated)",2000,-1000,1000);
  hPLAS_bPlasTDC_avg= MakeTH1('I',"bPlastic/Timing/TDC_Average","Average TDC",5000,-5000,5000);
  hFat_bplas_Corr_EFatEbPlas = MakeTH2('D',"Correlations/FATIMA_PLASTIC_General/EnergyvsEnergy","Coinc Fatima bPlastic Energy",2000, 0.0, 4000.0, 2000, 0.0, 4000);

}
  /**---------------------------------------------------------------------------------------------------------**/
 void EventAnlProc::Do_Fat_Plas_Histos(EventAnlStore* pOutput){
     double bPlas_SC41_dT_bPlasFatCorr[32], Fat_SC41_dT_bPlasFatCorr[50];
     double bPlas_SiPM_dT_bPlasFatCorr[32];
     double bPlas_GainMatch_bPlasFatCorr[32];
     int Fat_QDCID_bPlasFatCorr;
     double Fat_GainMatch_bPlasFatCorr[50];
     double bPlas_SC41_dT_bPlasFatCorr_sum, bPlas_SiPM_dT_bPlasFatCorr_sum;
     double bPlas_SC41_dT_bPlasFatCorr_avg, bPlas_SiPM_dT_bPlasFatCorr_avg;
     int bPlasTDChits;
     double  Fat_SC41_bPlas_bPlasFatCorr[50];
     double Fat_PM_bPlasFatCorr[50],  Fat_PMdT_bPlasdT_bPlasFatCorr[50];
     double Fat_bPlas_rawTDC_bPlasFatCorr[50];
     int bPlas_TDC_Incr;
     double bPlas_TDC_sum, bPlas_TDC_avg;
     int    Fat_bPlasTDCIDMain[50];
     int Fat_TDC_Multipl[50];
     double Fat_minus_plasticTDC[50];

       bPlas_SC41_dT_bPlasFatCorr_sum = 0;
       bPlas_SC41_dT_bPlasFatCorr_avg = 0;
       bPlas_SiPM_dT_bPlasFatCorr_avg = 0;
       bPlas_SiPM_dT_bPlasFatCorr_sum = 0;
       bPlas_TDC_sum = 0;
       Fat_QDCID_bPlasFatCorr = 0;
       bPlas_TDC_avg = 0;

       for (int i =0; i<32; i++){
           bPlas_SC41_dT_bPlasFatCorr[i] = 0;
           bPlas_SiPM_dT_bPlasFatCorr[i] = 0;
           bPlas_GainMatch_bPlasFatCorr[i] = 0;
       }
       for (int i =0; i<50; i++){
         Fat_SC41_bPlas_bPlasFatCorr[i] = 0;   
         Fat_GainMatch_bPlasFatCorr[i] = 0;
         Fat_minus_plasticTDC[i] = 0;
         Fat_PM_bPlasFatCorr[i] = 0;
         Fat_PMdT_bPlasdT_bPlasFatCorr[i] = 0;
         Fat_bPlas_rawTDC_bPlasFatCorr[i] = 0;
         Fat_bPlasTDCIDMain[i] = 0;
         Fat_TDC_Multipl[i] = 0;
         Fat_SC41_dT_bPlasFatCorr[i] = 0;
  }

     /**----------------------------------------------------------------------------------------------**/
            ///Get FATIMA Energy ///
  for (int i=0; i<FatQDCFired; i++){
         Fat_QDCID_bPlasFatCorr = FatQDCID[i];
       ///Calibrated Fatima Energy
         Fat_GainMatch_bPlasFatCorr[Fat_QDCID_bPlasFatCorr] = pOutput->pFat_QDCGainMatch[Fat_QDCID_bPlasFatCorr];
  }
            ///Get bPLASTIC Time///
            bPlasTDChits = 0;
    for (int i = 0; i<bPlasTDCFired; i++){

            /// Plastic TDC
            int bPlasTDCID_bPlasFatCorr = bPlasTDCID[i];
          //  bPlas_TDCRaw[bPlasTDCID_bPlasFatCorr] =
             ///Take only first hit of a given channel
            bPlas_TDC_Incr = 0;
            for(int j=0; j<=i; j++){
             if (bPlasTDCID_bPlasFatCorr != bPlasTDCID[j] ) bPlas_TDC_Incr++;

            }
            if(bPlas_TDC_Incr == i ){
            bPlas_SC41_dT_bPlasFatCorr[bPlasTDCID_bPlasFatCorr] = pOutput-> pbPlas_SC41_dT[bPlasTDCID_bPlasFatCorr];
            bPlas_SiPM_dT_bPlasFatCorr[bPlasTDCID_bPlasFatCorr] = pOutput-> pbPlas_SiPM_dT_Calib[bPlasTDCID_bPlasFatCorr];

            }
           /// Sum of S41-bPlas TDC dT
           if(bPlasTDCID_bPlasFatCorr!=0 &&bPlas_SC41_dT_bPlasFatCorr[bPlasTDCID_bPlasFatCorr]>0 ){
              bPlas_SC41_dT_bPlasFatCorr_sum += bPlas_SC41_dT_bPlasFatCorr[bPlasTDCID_bPlasFatCorr];
              bPlas_SiPM_dT_bPlasFatCorr_sum += bPlas_SiPM_dT_bPlasFatCorr[bPlasTDCID_bPlasFatCorr];
              bPlas_TDC_sum += bPlasTDC_TS[i][bPlasTDCID_bPlasFatCorr];
                  bPlasTDChits++; //needed for averaging later
                     }
      }
                ///Get FATIMA TDC
  for (int i=0; i<FatTDCFired; i++){
    Fat_bPlasTDCIDMain[i] = FatTDCID[i]; //FAT TDC ID

    // Fat_bPlasTDC_TS_Raw[i] = FatTDC_TS[i][Fat_bPlasTDCIDMain[i]]; //Fat TDC Data
    Fat_TDC_Multipl[Fat_bPlasTDCIDMain[i]]++; //Fat TDC Multiplicity
                             ///Call in SC41-Fat and Fat ref-Fat Ch.
                Fat_SC41_dT_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] = pOutput->pFat_SC41_dT_Calib[Fat_bPlasTDCIDMain[i]];
               // if(Fat_SC41_dT_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] && Fat_bPlasTDC_TS_Raw[i]>0 && Fat_bPlasTDCIDMain[i]<41){


//                     for (int k=0; k<bPlasTDCFired;k++){
//
//                         bPlas_Fat_TDC_IDMain[k] = bPlasTDCID[k]; //Plastic TDC ID
//                         bPlas_TDC_Multipl[k][bPlas_Fat_TDC_IDMain[k]]++; // Plastic TDC Multiplicity
//                         }
                ///Average over Plastic TDC events (only if the channel fires once)
                    bPlas_SC41_dT_bPlasFatCorr_avg = bPlas_SC41_dT_bPlasFatCorr_sum/bPlasTDChits;
                    bPlas_SiPM_dT_bPlasFatCorr_avg = bPlas_SiPM_dT_bPlasFatCorr_sum/bPlasTDChits;

                    hPLAS_TimeDiff_SC41_avg ->Fill(bPlas_SC41_dT_bPlasFatCorr_avg);
                    hPLAS_TimeDiff_SiPM_avg ->Fill(bPlas_SiPM_dT_bPlasFatCorr_avg);
                ///For each Fatima: (SC41-Fatima) - (SC41-Plastic) Taking only the first Fat hit/Channel
                    int Fat_TDC_Incr_bPlasFatCorr = 0;
                    for(int j=0; j<=i; j++){
                        if (Fat_bPlasTDCIDMain[i]!= FatTDCID[j] ) Fat_TDC_Incr_bPlasFatCorr++;

                    }

            if(Fat_bPlasTDCIDMain[i]<40 ){

                ///Raw FAT TDC - Raw Plastic TDC1
                    bPlas_TDC_avg =  bPlas_TDC_sum/bPlasTDChits;

                   hPLAS_bPlasTDC_avg -> Fill(bPlas_TDC_avg);
              
                   if(Fat_GainMatch_bPlasFatCorr[Fat_QDCID_bPlasFatCorr]>fCorrel->GFat_Egate_low && Fat_GainMatch_bPlasFatCorr[Fat_QDCID_bPlasFatCorr]<fCorrel->GFat_Egate_high ){

                   if(FatTDC_TS[i][Fat_bPlasTDCIDMain[i]]>fCorrel->GFat_TRawgate_low && FatTDC_TS[i][Fat_bPlasTDCIDMain[i]]<fCorrel->GFat_TRawgate_high){
                    Fat_bPlas_rawTDC_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] = (FatTDC_TS[i][Fat_bPlasTDCIDMain[i]] - bPlas_TDC_avg);
                    //cout<<"ev " << event_number << " bPlas_TDC_sum " << bPlas_TDC_sum << " bPlasTDChits " << bPlasTDChits << " bPlas_TDC_avg " << bPlas_TDC_avg << "  FatTDC_TS " << FatTDC_TS[i][Fat_bPlasTDCIDMain[i]] << " Fat_bPlasTDCIDMain[i] " << Fat_bPlasTDCIDMain[i]<< "  Fat_bPlas_rawTDC_bPlasFatCorr " <<  Fat_bPlas_rawTDC_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] << " i " << i << endl;
                    hFat_bplas_Corr_RawTDCFAT_RawTDCbP_Dets[Fat_bPlasTDCIDMain[i]]->Fill(Fat_bPlas_rawTDC_bPlasFatCorr[Fat_bPlasTDCIDMain[i]]);
                    // cout <<"event " << event_number << " Fat_bPlasTDC_TS_Raw "<< Fat_bPlasTDC_TS_Raw[Fat_bPlasTDCIDMain[i]] <<" bPlasTDC_TS " << bPlasTDC_TS[k][bPlasTDCID[k]]<<" Fat_bPlas_rawTDC_bPlasFatCorr " << Fat_bPlas_rawTDC_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] << " Fat_bPlasTDCIDMain[i] " <<Fat_bPlasTDCIDMain[i] <<" bPlasTDCID[j] " << bPlasTDCID[k]<< endl;
                   //}
           

                ///PM-Fatima dT - SiPM-bPlas dT average
                Fat_PM_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] =  pOutput->pFat_Ch_dT[Fat_bPlasTDCIDMain[i]];
                Fat_PMdT_bPlasdT_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] = (Fat_PM_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] - bPlas_SiPM_dT_bPlasFatCorr_avg);
                
                hFat_bplas_Corr_SiPM_Dets[Fat_bPlasTDCIDMain[i]] -> Fill(Fat_PMdT_bPlasdT_bPlasFatCorr[Fat_bPlasTDCIDMain[i]]);

                if(bPlas_SC41_dT_bPlasFatCorr_avg>0){
                ///SC41-Fatima - SC41-bPlas dT average
                    Fat_SC41_bPlas_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] =  Fat_SC41_dT_bPlasFatCorr[Fat_bPlasTDCIDMain[i]] - bPlas_SC41_dT_bPlasFatCorr_avg; ///Fatima - Plastic SC41
                    hFat_bplas_Corr_SC41_Dets[Fat_bPlasTDCIDMain[i]]->Fill(Fat_SC41_bPlas_bPlasFatCorr[Fat_bPlasTDCIDMain[i]]);


                    ///Fatima Energy gate
                     hFat_bplas_Corr_SC41_Dets_gated[Fat_bPlasTDCIDMain[i]]->Fill(Fat_minus_plasticTDC[i]);
                     hFat_bplas_Corr_SiPM_Gated_Dets[Fat_bPlasTDCIDMain[i]] -> Fill(Fat_PMdT_bPlasdT_bPlasFatCorr[Fat_bPlasTDCIDMain[i]]);
                        }
                    }
                    for (int k = 0; k<bPlasTDCFired; k++){
                    ///Fatima Energy - bPlastic Energy
                  bPlas_GainMatch_bPlasFatCorr[ bPlasTDCID[k]-1] = pOutput->pbPlas_QDCGainMatch_i[ bPlasTDCID[k]-1];
                  hFat_bplas_Corr_EFatEbPlas->Fill(Fat_GainMatch_bPlasFatCorr[Fat_QDCID_bPlasFatCorr], bPlas_GainMatch_bPlasFatCorr[ bPlasTDCID[k]-1] );
                    }                    
                }
            }
        }
    }

/**----------------------------------------------------------------------------------------------**/
/**--------------------------------------  FINGER - Plastic (temporary for testing purposes) ----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventAnlProc::Make_Fing_Plas_Histos(){
  for (int i =0; i<32; i++){
    h_FingToT_PlasT[i] = MakeTH2('D',Form("Correlations/FINGER_bPlas/FINGER_bPlas_ToT_bPlasT%02d",i), Form ("FING_ToT_bPlastT_%2d",i),52,0,52,2500,-5000,5000);
  }
  h_FingStrip_PlasID = MakeTH2('D',"Correlations/FINGER_bPlas/FINGER_Strip_TDCID","FINGER Strip vs PlasTDCID",52,0,52,35,0,35);
  h_FingToT_PlasE = MakeTH2('D',"Correlations/FINGER_bPlas/FINGER_bPlas_ToT_QDCE","FINGER ToT vs PlasQDC Energy",2000,0,200, 20000,0,20000);
}
void EventAnlProc::Do_Fing_Plas_Histos(){

  double bPlasFingQDCGainMatch[32];
  int bPlasFingTDCIDMain = 0;
  double bPlasFingTDC_TS_Main[32]; 
  double bPlasFing_SiPM_diff[32]; 
  int bPlasFing_FingStripID = 0;

  for(int i = 0; i<32; i++){
      bPlasFingQDCGainMatch[i] = 0;
      bPlasFingTDC_TS_Main[i] = 0;
      bPlasFing_SiPM_diff[i] = 0;
  }      

  //Test: bPlas QDC vs Finger ToT
     for (int i=0; i<bPlasQDCFired; i++){
    // cout << "1) bPlasQDC[i] " << bPlasQDC[i] << endl;
    bPlasFingQDCGainMatch[i] = fCal->AplasQDC[i]*bPlasQDC[i] + fCal->BplasQDC[i]; //gain matching
    for(int j =0; j<Fing_firedTamex; j++){
      for (int k =0; k<Fing_iterator[k]; k++){

        if(bPlasFingQDCGainMatch[i]>0){

          h_FingToT_PlasE ->Fill(Fing_TOT[j][k],bPlasFingQDCGainMatch[i]);
        }
       
      }
    }
  }
  ///bPlas TDC SiPM1-SiPM Ch.x vs Finger ToT (all)
  for (int i=0; i<bPlasTDCFired;i++){
    bPlasFingTDCIDMain = bPlasTDCID[i];
    bPlasFingTDC_TS_Main[bPlasFingTDCIDMain] = bPlasTDC_TS[i][bPlasFingTDCIDMain];
    bPlasFing_SiPM_diff[bPlasFingTDCIDMain] = (bPlasFingTDC_TS_Main[1] - bPlasFingTDC_TS_Main[bPlasFingTDCIDMain]); //in ns

    //Finger part
    for(int j =0; j<Fing_firedTamex; j++){
      for (int k =0; k<Fing_iterator[k]; k++){
        if(k<32){
        // if(fCond_FingToTvsStrip->Test(Fing_TOT[j][k], Fing_leadChan[i][j]))
                        
          h_FingToT_PlasT[bPlasFingTDCIDMain]->Fill(Fing_leadChan[j][k], bPlasFing_SiPM_diff[bPlasFingTDCIDMain]);
          h_FingStrip_PlasID ->Fill(Fing_leadChan[j][k],bPlasFingTDCIDMain);
                        }
                    }
    }
  }
}
