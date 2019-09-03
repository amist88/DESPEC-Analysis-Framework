#include "EventUnpackProc.h"
#include "EventUnpackStore.h"
#include "Riostream.h"

// Root Includes //
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TArc.h"
#include "TTree.h"

#include <time.h>
#include <math.h>
#include <iomanip>

// Go4 Includes //
#include "TGo4UserException.h"
#include "TGo4Picture.h"
#include "TGo4MbsEvent.h"

#include "TGo4MbsSubEvent.h"

// General Includes //
#include <fstream>
#include <vector>

#include "Detector_System.cxx"
#include "AIDA_Detector_System.h"
#include "FATIMA_Detector_System.h"
#include "PLASTIC_Detector_System.h"
#include "PLASTIC_VME_Detector_System.h"
#include "FINGER_Detector_System.h"
#include "GALILEO_Detector_System_TEST.h"
#include "FRS_Detector_System.h"

#include "TAidaConfiguration.h"

#include "CalibParameter.h"
#include "CorrelParameter.h"

#include "Data_Stream.cxx"
#include "White_Rabbit.h"

#include <string>


using namespace std;


//***********************************************************
EventUnpackProc::EventUnpackProc() :TGo4EventProcessor("Proc")
{
  cout << "**** EventUnpackProc: Create instance " << endl;


}
//***********************************************************
// standard factory
EventUnpackProc::EventUnpackProc(const char* name) : TGo4EventProcessor(name)
{


  cout << "**** EventUnpackProc: Create" << endl;

  //  input_data_path_old = "old";
   // WR_out.open ("WR_diff_store_270289.txt");
  WR_used = false;

  //used_systems
  get_used_Systems();
  get_WR_Config();
    //read_setup_parameters();

  //  FAT_det_pos_setup();

  //create White Rabbit obj
  WR = new White_Rabbit();

  fCal = new CalibParameter("CalibPar");
  AddParameter(fCal);
  if (fCal) fCal->PrintParameter(0,0);
  else cout << "**** ERRR - CalibPar doesn't exist - program will crash.\n";

  fCorrel = new CorrelParameter("CorrelPar");
  AddParameter(fCorrel);
  if (fCorrel) fCorrel->PrintParameter(0,0);
  else cout << "**** ERRR - CorrelPar doesn't exist - program will crash.\n";

  //create Detector Systems
  Detector_Systems = new Detector_System*[7];

  // all non used systems intialized as NULL
  //-> calling uninitialized system will cause an error !
  Detector_Systems[0] = !Used_Systems[0] ? nullptr : new FRS_Detector_System();
  Detector_Systems[1] = !Used_Systems[1] ? nullptr : new AIDA_Detector_System();
  Detector_Systems[2] = !Used_Systems[2] ? nullptr : new PLASTIC_VME_Detector_System();
  Detector_Systems[3] = !Used_Systems[3] ? nullptr : new FATIMA_Detector_System();
  Detector_Systems[4] = !Used_Systems[4] ? nullptr : new GALILEO_Detector_System();
  Detector_Systems[5] = !Used_Systems[5] ? nullptr : new FINGER_Detector_System();

  //     if(Used_Systems[2] && VME_TAMEX==true){
  //          Detector_Systems[5] = new PLASTIC_VME_Detector_System();
  //     }
  //          else Detector_Systems[5] = nullptr;
  // Detector_Systems[5] = (!Used_Systems[2] && VME_TAMEX==false) ? nullptr : new PLASTIC_VME_Detector_System();

  for(int i = 0;i < 6;++i) if(!Used_Systems[i]) Detector_Systems[i] = nullptr;

    PLASTIC_CALIBRATION = Used_Systems[2] ? Check_Cal_Plastic() : false; //Not needed anymore

  //Only create histograms if system is used
  if(Used_Systems[0]) Make_FRS_Histos();

  if(Used_Systems[1]) Make_AIDA_Histos();

  if(Used_Systems[3]) Make_FATIMA_Histos();

  if(Used_Systems[2] && !PLASTIC_CALIBRATION && VME_TAMEX==false) Make_Plastic_Histos();

  if(Used_Systems[2] && !PLASTIC_CALIBRATION && VME_TAMEX==true) Make_Plastic_VME_Histos();

  if(Used_Systems[4]) Make_GALILEO_Histos();

  //if(Used_Systems[5]) Make_Finger_Histos();

  get_interest_arrays();
  //Skip event building if plastic calibration is enabled
  // SKIP_EVT_BUILDING = SKIP_EVT_BUILDING || PLASTIC_CALIBRATION;

  if(!SKIP_EVT_BUILDING){
    EvtBuilder = new EventBuilder*[1];
    EvtBuilder[0] = new Time_EventBuilder(amount_interest,length_interest,interest_array);

  }

  checkTAMEXorVME();
  checkPADI_or_PADIWA();

  //Raw_Event object to handle data
  //RAW = new Raw_Event(PADI_OR_PADIWA);

  RAW = new Raw_Event();

  load_PrcID_File();

  load_FingerID_File();

  read_setup_parameters();
//
//     FAT_det_pos_setup();
    WR_count = 0;
    count = 0;
    array_count = 0;
    iterator = 0;
    val_it = 0;

  //  FAT_gain_match_done = false;

  Cout_counter = 0;

  //Clear for AIDA
  lastTime = 0;
  ID = 0;
  totalEvents = 0;
  startTime = 0;
  stopTime = 0;
  fAida.ImplantEvents.clear();
  fAida.DecayEvents.clear();
  fAida.Implants.clear();
  fAida.Decays.clear();
  // Setup AIDA arrays
  if(Used_Systems[1])
  {
    TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
    adcLastTimestamp.resize(conf->FEEs());
    adcCounts.resize(conf->FEEs());
  }

}

//----------------------------------------------------------


EventUnpackProc::~EventUnpackProc()
{

  if(!SKIP_EVT_BUILDING){
    cout << "------------------" << endl;
    cout << "Deleting Event Builder" << endl;
    delete EvtBuilder[0];
    EvtBuilder[0] = nullptr;
    delete[] EvtBuilder;
  }
  string DET_NAME[7] = {"FRS","AIDA","PLASTIC_TAMEX","FATIMA","GALILEO","FINGER","PLASTIC_VME"};
  //Detector_Systems[3]->write();
  cout << "------------------" << endl;
  for(int i = 0;i < 6;++i){
    if(Detector_Systems[i]){
      delete Detector_Systems[i];
      Detector_Systems[i] = nullptr;
      //cout<<"Detector_System " << DET_NAME[i] << " deleted" << endl;
    }
  }
  cout << "------------------" << endl;
  for(int i = 0;i < 10;++i) if(interest_array[i]) delete[] interest_array[i];
  delete[] interest_array;
  delete[] length_interest;
  delete[] Detector_Systems;


  delete RAW;//Clear for AIDA
  delete WR;
  cout << "**** EventUnpackProc: Delete instance" << endl;
}

//----------------------------------------------------------
Bool_t EventUnpackProc::BuildEvent(TGo4EventElement* dest)
{
  bool skip;
  static TString oldfname = "";
  TGo4MbsEvent *fMbsEvent = dynamic_cast<TGo4MbsEvent*>    (GetInputEvent("Unpack"));
  s_filhe* fileheader=fMbsEvent->GetMbsSourceHeader();
  //s_bufhe* head = GetMbsBufferHeader();

  EventUnpackStore* fOutput = (EventUnpackStore*) dest;
  TGo4MbsEvent* fInput = (TGo4MbsEvent*) GetInputEvent();

  if(fOutput==0){
    cout << "UnpackProc: no unpack output event !";
    return false;
  }
  input_data_path = fileheader->filhe_file;

  count++;

  if (count % 100000 == 0){

    cout << "\r";
    cout << "Event " << count << " Reached!!!"<<"    Data File Number : "<<data_file_number;
    cout <<"\t\t\t\t";
    cout.flush();
  }

  Bool_t isValid=kFALSE; // validity of output event //

  if (fInput==0) // Ensures that there is data in the event //
  {
    cout << "EventUnpackProc: no input event !"<< endl;
    fOutput->SetValid(isValid);
    return isValid;
  }
  isValid=kTRUE;
  event_number=fInput->GetCount();
  fOutput-> fevent_number = event_number;
  // fOutput->flestore=1;
  //  cout << "fOutput->fevent_number " << fOutput->fevent_number<< endl;


  fInput->ResetIterator();
  TGo4MbsSubEvent* psubevt(0);


  // ------------------------------------------------------ //
  // |                                                    | //
  // |               START OF EVENT ANALYSIS              | //
  // |                                                    | //
  // ------------------------------------------------------ //
  //Temporary to fix crashes caused by TAMEX 06.04.19 
  //Fixed 12.07.19 but what happened needs to be investigated
  
   // if (event_number==270283 || event_number==270282 ||  event_number==270281||  event_number==270280 ||  event_number==270279 ||  event_number==270278 ||  event_number==270277 ||  event_number==270276 ||  event_number==270275 ||  event_number==270274 ||  event_number==270273 ){
      // if (event_number==2144767 ){
          
//     if ((event_number==974467 && input_data_path == "/media/sdc1/April_online/f0007_higheri_0001.lmd")
//     ||((event_number==3927530 ||event_number==3998732) && input_data_path == "/media/sdc1/April_online/f000-3_0004.lmd")
//     ||((event_number==4188593 || event_number==4587527 || event_number==4615226) && input_data_path == "/media/sdc1/April_online/f000-3_0005.lmd")
//     ||((event_number==3921332 && input_data_path=="/media/sdc1/April_online/f005_diffuse_0001.lmd"))
//     ||((event_number==5701008 ||5784742 ||7009167) && input_data_path=="/media/sdc1/April_online/f0006_maskTPC0001.lmd")
//     ||((event_number==7009167||event_number==7037211||event_number==8118970 ||event_number== 8127984) && input_data_path=="/media/sdc1/April_online/f0006_maskphase20001.lmd")
//     ||((event_number==7009167) && input_data_path=="/media/sdc1/April_online/f0006_maskphase20001.lmd"))
// 
//     {
//       skip =true;
//     }
     if (event_number>0){
        //  cout<<"done" <<endl;
  //  if(skip ==false){ //needed if some event causes a crash
   // cout << "Event_num " << event_number << endl;
      int subevent_iter = 0;

      Int_t PrcID_Conv = 0;

      Int_t* pdata = nullptr;
      Int_t lwords = 0;
      Int_t PrcID = 0;
      Int_t sub_evt_length = 0;
      WR_tmp = 0;
      WR_d=0;
      AIDA_Loop = 0;
      WR_main=0;
    for (int i=0; i<8196; i++){
        WR_AIDA[i]=0;
        WR_diff[i] = 0;
      }

      while ((psubevt = fInput->NextSubEvent()) != 0) // subevent loop //
      {
        subevent_iter++;
        pdata = psubevt->GetDataField();
        lwords = psubevt->GetIntLen();
        PrcID = psubevt->GetProcid();
        PrcID_Conv = get_Conversion(PrcID);
        fOutput -> fProcID[PrcID_Conv] = PrcID_Conv;

        if(array_count==10000) array_count = 0; //reset the array after x events

        for (int i =0; i<7;i++){
          fOutput->fUsed_Systems[i] = Used_Systems[i];
     
        }

        sub_evt_length  = (psubevt->GetDlen() - 2) / 2;

    ///------------------------------WHITE RABBIT: UNDER DEVELOPMENT --------------------------------------////
        if(WHITE_RABBIT_USED){
          sub_evt_length = sub_evt_length - 5;

            //Pulls it straight from White_Rabbit class
          WR_d = WR->get_Detector_id();
            WR_tmp = WR->get_White_Rabbit(pdata);
            pdata = WR->get_pdata();
           
             // For Galileo and Fatima
           if(WR_d==3) fOutput->fFat_WR = WR_tmp;
           if(WR_d==4) fOutput->fGal_WR = WR_tmp;


            WR_main = WR_tmp;
     
            //For AIDA White Rabbit time differences
           AIDA_Loop = RAW->get_AIDA_HITS();
            for(int i=0; i<AIDA_Loop; i++) {

                WR_AIDA[i] = RAW-> get_AIDA_WR(i);

                if( WR_AIDA[i]>0 && WR_main>0){

              WR_diff[i] = (WR_AIDA[i] - WR_main); // in ns
           //   fOutput->fWR_Aida_Det_diff[i] = WR_diff[i]; ///--->This guy fills alot of data in the unpacktree.
              WR_count++;
                                    }
                            }
                     }
          //    cout <<"1) event_number " << event_number << " WR_diff[i] "<< WR_diff[i]/1000 << " i " << i <<endl;
            //cout <<"event " << event_number <<" WR_main " << WR_main <<" WR_AIDA[i] " << WR_AIDA[i] << "WR_diff[i]" <<WR_diff[i]<< endl;
            // cout <<"event " << event_number  << " WR_diff[i] " <<WR_diff[i]<< endl;
          // WR_out<< "event " << event_number <<" WR_main " << WR_main <<" WR_AIDA[i] " << WR_AIDA[i] << " WR_diff[i]" <<WR_diff[i]<< endl;
               
                //cout <<"1) event_number " << event_number << " WR_diff[i] "<< WR_diff[i]/1000 << " i " << i <<endl;

//                  if((WR_diff[i]/1000)<-10 || (WR_diff[i]/1000)>10) { //in mus
//                  WR_count=0;
//                  for(int i=0; i<100000; i++)
//                  WR_array_AIDA[i] = 0;
//                  WR_array[i] = 0 ;
//                  }
               //  cout <<"1) event_number " << event_number << " WR_diff[i] "<< WR_diff[i]/1000 << " i " << i <<endl;


//                 }
//                  if((WR_diff[i]/1000)>-10 || (WR_diff[i]/1000)<10) { //in mus
//                  //  WR_array_AIDA[WR_count] = WR_AIDA[i];
//                    //WR_array[WR_count] = WR_main;
//
//            // cout <<"ev_num " << event_number << " WR_count " << WR_count<< "  WR_diff[i]/1000 " <<  WR_diff[i]/1000<<endl;;
//             }


//         for (int i=0; i<WR_count; i++){
//             cout <<"1) WR_array_AIDA[WR_count] "<< WR_array_AIDA[WR_count] << endl;
//
//             for (int j=0; j<WR_count; j++){
//                 if(i!=j){
//             cout << "2) WR_array_AIDA[WR_count] " << WR_array_AIDA[WR_count] << endl;
//                 }
//             }
//         }
///-----------------------------------------------------------------------------------------------------------///
        //if necessary, directly print MBS for wanted Detector_System
        if(PrcID_Conv == AIDA && false) print_MBS(pdata,lwords);
        if(PrcID_Conv == FATIMA && false) print_MBS(pdata,lwords);
        if(PrcID_Conv == PLASTIC && false) print_MBS(pdata,lwords);
        if(PrcID_Conv == GALILEO && false) print_MBS(pdata,lwords);
        // if(PrcID_Conv == FINGER && false) print_MBS(pdata,lwords);

        //=================================================================
        //UNPACKING
        ///send subevent to respective unpacker

        if(Detector_Systems[PrcID_Conv] !=0){
        Detector_Systems[PrcID_Conv]->Process_MBS(psubevt);

        Detector_Systems[PrcID_Conv]->Process_MBS(pdata);

        ///get mbs stream data from unpacker (pointer copy solution)
        pdata = Detector_Systems[PrcID_Conv]->get_pdata();

        ///get data from subevent and send to RAW
        Detector_Systems[PrcID_Conv]->get_Event_data(RAW);
        }


        //=================================================================

        //cals_done = Detector_Systems[PrcID_Conv]->calibration_done();

        //=================================================================
        //Event Building
    //   if(!SKIP_EVT_BUILDING ) EvtBuilder[0]->set_Event(RAW);
        //=================================================================

      //  if(cals_done) break;

        //=================================================================
        //HISTOGRAM FILLING (only singles)
        FILL_HISTOGRAMS(PrcID_Conv);
        //=================================================================

        pdata = nullptr;

        ///--------------------------------------------------------------------------------------------///
                                /** Unpack Tree for each detector subsystem**/
        ///--------------------------------------------------------------------------------------------///
                                                /** Output FRS **/
        ///--------------------------------------------------------------------------------------------///
    if (Used_Systems[0] && PrcID_Conv==0){
        ///MUSIC
           for(int i =0; i<2; ++i){
            fOutput->fFRS_dE[i] = RAW->get_FRS_dE(i);
            fOutput->fFRS_dE_cor[i] = RAW->get_FRS_dE_corr(i);
           }
        ///SCI
           for(int l=0;l<12;++l){
            fOutput->fFRS_sci_l[l] = RAW->get_FRS_sci_l(l);
            fOutput->fFRS_sci_r[l] = RAW->get_FRS_sci_r(l);
            fOutput->fFRS_sci_e[l] = RAW->get_FRS_sci_e(l);
            fOutput->fFRS_sci_tx[l] = RAW->get_FRS_sci_tx(l);
            fOutput->fFRS_sci_x[l] = RAW->get_FRS_sci_x(l);
           }
            ///SCI TOF
        fOutput->fFRS_sci_tofll2 = RAW->get_FRS_tofll2();
        fOutput->fFRS_sci_tofll3 = RAW->get_FRS_tofll3();
        fOutput->fFRS_sci_tof2 = RAW->get_FRS_tof2();
        fOutput->fFRS_sci_tofrr2 = RAW->get_FRS_tofrr2();
        fOutput->fFRS_sci_tofrr3 = RAW->get_FRS_tofrr3();
        fOutput->fFRS_sci_tof3 = RAW->get_FRS_tof3();
        ///ID 2 4
        fOutput->fFRS_ID_x2 = RAW->get_FRS_x2();
        fOutput->fFRS_ID_y2 = RAW->get_FRS_y2();
        fOutput->fFRS_ID_a2 = RAW->get_FRS_a2();
        fOutput->fFRS_ID_b2 = RAW->get_FRS_b2();

        fOutput->fFRS_ID_x4 = RAW->get_FRS_x4();
        fOutput->fFRS_ID_y4 = RAW->get_FRS_y4();
        fOutput->fFRS_ID_a4 = RAW->get_FRS_a4();
        fOutput->fFRS_ID_b4 = RAW->get_FRS_b4();
            ///SCI dT
        fOutput->fFRS_sci_dt_21l_21r = RAW->get_FRS_dt_21l_21r();
        fOutput->fFRS_sci_dt_41l_41r = RAW->get_FRS_dt_41l_41r();
        fOutput->fFRS_sci_dt_42l_42r = RAW->get_FRS_dt_42l_42r();
        fOutput->fFRS_sci_dt_43l_43r = RAW->get_FRS_dt_43l_43r();

        fOutput->fFRS_sci_dt_21l_41l = RAW->get_FRS_dt_21l_41l();
        fOutput->fFRS_sci_dt_21r_41r = RAW->get_FRS_dt_21r_41r();

        fOutput->fFRS_sci_dt_21l_42l = RAW->get_FRS_dt_21l_42l();
        fOutput->fFRS_sci_dt_21r_42r = RAW->get_FRS_dt_21r_42r();
            ///ID Beta Rho
        for(int i =0; i<2; ++i){
            fOutput->fFRS_ID_brho[i] = RAW->get_FRS_brho(i);
            fOutput->fFRS_ID_rho[i] = RAW->get_FRS_rho(i);
        }
        fOutput->fFRS_beta = RAW->get_FRS_beta();
        fOutput->fFRS_beta3 = RAW->get_FRS_beta3();
        fOutput->fFRS_gamma = RAW->get_FRS_gamma();
            ///ID Z AoQ
        fOutput->fFRS_AoQ = RAW->get_FRS_AoQ();
        fOutput->fFRS_AoQ_corr = RAW->get_FRS_AoQ_corr();

        fOutput->fFRS_z = RAW->get_FRS_z();
        fOutput->fFRS_z2 = RAW->get_FRS_z2();
        fOutput->fFRS_z3 = RAW->get_FRS_z3();
            ///ID Timestamp
        fOutput->fFRS_timestamp = RAW->get_FRS_timestamp();
        fOutput->fFRS_ts = RAW->get_FRS_ts();
        fOutput->fFRS_ts2 = RAW->get_FRS_ts2();
         }
         ///--------------------------------------------------------------------------------------------///
                                            /** Output AIDA **/
        ///--------------------------------------------------------------------------------------------///

        if (Used_Systems[1] && PrcID_Conv==1){
          TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
          AIDA_Hits = RAW->get_AIDA_HITS();

          AidaEvent evt;
          fOutput->fAIDAHits = AIDA_Hits;
          for(int i = 0; i<AIDA_Hits; i++){

            AIDA_Energy[i] = RAW->get_AIDA_Energy(i);

            AIDA_FEE[i] = RAW-> get_AIDA_FEE_ID(i);
            AIDA_ChID[i] = RAW-> get_AIDA_CHA_ID(i);
            AIDA_Time[i] = RAW-> get_AIDA_WR(i);
            AIDA_HighE_veto[i] = RAW-> get_AIDA_HighE_VETO(i);
            AIDA_Side[i] = RAW-> get_AIDA_SIDE(i);
            AIDA_Strip[i] = RAW-> get_AIDA_STRIP(i);
            AIDA_evtID[i] = RAW-> get_AIDA_EVTID(i);

            evt.Channel = AIDA_ChID[i];
            evt.Module = AIDA_FEE[i];
            evt.Time = AIDA_Time[i];
            evt.HighEnergy =  AIDA_HighE_veto[i];
            evt.DSSD = conf->FEE(AIDA_FEE[i]).DSSD;
            evt.Side = AIDA_Side[i];
            evt.Strip = AIDA_Strip[i];
            evt.ID = AIDA_evtID[i];
            evt.Energy = AIDA_Energy[i];
            evt.FastTime = RAW->get_AIDA_FastTime(i);

            if (!startTime) startTime = evt.Time;
            stopTime = evt.Time;
            /// Build events from everything until there's a gap of 2000 �s (event window)

            /// If lastTime is 0 it's the first event
            /// New event for timewarps
            if (lastTime > 0 && (evt.Time - lastTime > conf->EventWindow()))
            {
              // if event happened too late, redo the event again with a new out_event
              lastTime = 0;
              ResetMultiplexer();

              totalEvents++;
        
              fOutput->Aida.push_back(fAida);
              fAida.ImplantEvents.clear();
              fAida.DecayEvents.clear();
              fAida.AIDATime = 0;
            }

            lastTime = evt.Time;
            CorrectTimeForMultiplexer(evt);
            

            if (evt.HighEnergy==1)
            {

              fAida.ImplantEvents.push_back(evt);

            }
            else
            {
              fAida.DecayEvents.push_back(evt);

            }
           
            if (fOutput->AIDATime == 0)
            {
              fAida.AIDATime = evt.Time;

            }
          }
        }
        ///--------------------------------------------------------------------------------------------///
                                        /** Output bPlastic + SCALAR **/
        ///--------------------------------------------------------------------------------------------///

        //  int Plas_VME_QDC_ID;
        int Plas_VME_TDC_ID_sing, Plas_VME_TDC_ID[50], bPlasTDCmulti[50];
        Plas_VME_TDC_ID_sing=-1;

        for (int i =0; i<50; i++){
          Plas_VME_TDC_ID[i] = -1;
          bPlasTDCmulti[i] = 0;
        }

        if (Used_Systems[2] && PrcID_Conv==2){
          fOutput->fbPlas_VME_firedQDC = RAW->get_plastic_VME_QDC_fired();

          for (int i = 0; i < RAW->get_plastic_VME_QDC_fired();  i++){
            // Plas_VME_QDC_ID = RAW->get_plastic_VME_QDC_cha(i);
            // bPlasQDCmulti[Plas_VME_QDC_ID]++;

            fOutput-> fbPlas_VME_QDC_ID[i] = RAW->get_plastic_VME_QDC_cha(i);
            fOutput-> fbPlas_VME_QDC_E[i] = RAW->get_plastic_VME_QDC_dat1(i);
            // fOutput-> fbPlas_VME_QDC_Multiplicity = bPlasQDCmulti[Plas_VME_QDC_ID];
                       
            }
          fOutput->fbPlas_VME_firedTDC = RAW->get_plastic_VME_TDC_fired();

          for (int j=0; j < RAW->get_plastic_VME_TDC_fired(); j++){
            Plas_VME_TDC_ID[j] = RAW->get_plastic_VME_TDC_cha(j); //Plastic TDC ID
            Plas_VME_TDC_ID_sing = RAW->get_plastic_VME_TDC_cha(j);

            bPlasTDCmulti[Plas_VME_TDC_ID_sing]++;
            fOutput->fbPlas_VME_TDC_ID[j] = Plas_VME_TDC_ID[j];
            fOutput->fbPlas_VME_TDC_TS[j][Plas_VME_TDC_ID[j]] = RAW->get_plastic_VME_TDC_dat(Plas_VME_TDC_ID[j])*0.025;//25ps //Plastic TDC Data
            fOutput->fbPlas_VME_TDC_Multiplicity[Plas_VME_TDC_ID_sing] =  bPlasTDCmulti[Plas_VME_TDC_ID_sing];

          }

          fOutput->fScalar_fired = RAW->get_scalar_iterator();
          for (int g=0; g < RAW->get_scalar_iterator(); g++){
            fOutput-> fScalar_ID = RAW->get_scalar_chan(g);
           
          }
        }
        ///--------------------------------------------------------------------------------------------///
                                                /**Output FATIMA **/
        ///--------------------------------------------------------------------------------------------///

        int Fat_QDC_ID;
        int Fat_TDC_ID_sing;
        int Fat_TDC_ID[48];
        int Fat_TDC_multi[48];

        for (int i = 0; i<48; i++){
          Fat_TDC_multi[i] = 0;
         
        }
        if (Used_Systems[3]&& PrcID_Conv==3){
          fOutput->fFat_firedQDC = RAW->get_FAT_QDCs_fired();

          for (int i=0; i<RAW->get_FAT_QDCs_fired(); i++){

            fOutput->fFat_QDC_ID[i] =  RAW->get_FAT_QDC_id(i);
            fOutput->fFat_QDC_E[i] = RAW->get_FAT_QLong_Raw(i);
            fOutput->fFat_QDC_T[i] = RAW->get_FAT_t_qdc(i);
          }
          fOutput->fFat_firedTDC = RAW->get_FAT_TDCs_fired();
          for (int j=0; j<RAW->get_FAT_TDCs_fired(); j++){
            Fat_TDC_ID[j] =  RAW->get_FAT_TDC_id(j);
            Fat_TDC_ID_sing = RAW->get_FAT_TDC_id(j);

            Fat_TDC_multi[Fat_TDC_ID_sing]++;

            fOutput->fFat_TDC_ID[j] =  RAW->get_FAT_TDC_id(j);
            fOutput->fFat_TDC_TS[j][Fat_TDC_ID[j]] = RAW->get_FAT_TDC_timestamp(j);
            fOutput->fFat_TDC_Multiplicity[Fat_TDC_ID_sing] =  Fat_TDC_multi[Fat_TDC_ID_sing];

            /** SC41 Trigger **/
            if(Fat_TDC_ID[j]==40){
            fOutput->fSC41[j] =  RAW->get_FAT_TDC_timestamp(j);
            }
          }
        }
        ///--------------------------------------------------------------------------------------------///
                                            /**Output GALILEO **/
        ///--------------------------------------------------------------------------------------------///
        if (Used_Systems[4]&& PrcID_Conv==4){
          fOutput->fGal_fired = RAW->get_GALILEO_am_Fired();
          for (int i=0; i<RAW->get_GALILEO_am_Fired(); i++){
            int Gal_ID = RAW->get_GALILEO_Det_ids(i);
            fOutput->fGal_ID[i] =  RAW->get_GALILEO_Det_ids(i);
            fOutput->fGal_E[Gal_ID] = RAW->get_GALILEO_Chan_E(i)/1000;
            fOutput->fGal_T[Gal_ID] = RAW->get_GALILEO_Chan_T(i);
            fOutput->fGal_Pileup = RAW->get_GALILEO_Pileup(i);
          }
        }
        ///--------------------------------------------------------------------------------------------///
                                        /** Output FINGER **/
  ///--------------------------------------------------------------------------------------------///
        if (Used_Systems[5]&& PrcID_Conv==5){

          int Phys_Channel_Lead[4][16] = {0,0};
          int Phys_Channel_Trail[4][16] = {0,0};
          int trailHits[4] = {0};
          int leadHits[4] = {0};
          //int MaxHits = 0;
          int fingfired[4] = {0};

          fOutput->ffing_tamexhits = RAW->get_FINGER_tamex_hits();

          for (int i=0; i<RAW->get_FINGER_tamex_hits(); i++){
            leadHits[i] =  RAW->get_FINGER_lead_hits(i);
            trailHits[i] = RAW->get_FINGER_trail_hits(i);
            fingfired[i] = RAW->get_FINGER_am_Fired(i);
            fOutput->ffing_leadHits[i] = leadHits[i];
            fOutput->ffing_trailHits[i] = trailHits[i];
            fOutput->ffing_iterator[i] = fingfired[i];
            fOutput->ffing_Trig[i] = RAW->get_FINGER_trigger_T(i);

            for(int j = 0;j < fingfired[i];j++){

              fOutput->ffing_chID[i][j] = RAW->get_FINGER_CH_ID(i,j);
              fOutput->ffing_lead_coarse[i][j] = RAW->get_FINGER_coarse_lead(i,j);
              fOutput->ffing_lead_fine[i][j] = RAW->get_FINGER_fine_lead(i,j);
              fOutput->ffing_trail_coarse[i][j] = RAW->get_FINGER_coarse_trail(i,j);
              fOutput->ffing_trail_fine[i][j] = RAW->get_FINGER_fine_trail(i,j);
            
              if(fOutput->ffing_chID[i][j] % 2 == 0){
                Phys_Channel_Lead[i][j] = fingID[i][j/2]; //From allocation file
                fOutput->ffing_Lead_Phys_Chan[i][j] = Phys_Channel_Lead[i][j];  //Lead Chan.
                fOutput->ffing_Lead_T[i][j] = RAW->get_FINGER_lead_T(i,j); //Lead Time
          }

              else{
                Phys_Channel_Trail[i][j] = RAW->get_FINGER_physical_channel(i,j);
                fOutput->ffing_Trail_Phys_Chan[i][j] = Phys_Channel_Trail[i][j]; //Trail Chan.
                fOutput->ffing_Trail_T[i][j] = RAW->get_FINGER_trail_T(i,j); //Trail Time
              }

              fOutput->ffing_TOT[i][j] = RAW->get_FINGER_TOT(i,j);

            }
          }
        }
        ///--------------------------------------------------------------------------------------------///

      } //End of subevent loop


      fOutput->SetValid(isValid);

      pdata = nullptr;
    //} //End of Skip
  }
  //
  return isValid;

}

void EventUnpackProc::FILL_HISTOGRAMS(int PrcID_Conv){
  switch(PrcID_Conv){
    case 0:
    Fill_FRS_Histos();
    break;
    case 1:
    //   Process_AIDA_Event(fOutput);
    Fill_AIDA_Histos();
    break;
    case 2:
    if(!PLASTIC_CALIBRATION && VME_TAMEX==false) Fill_Plastic_Histos();
    if(!PLASTIC_CALIBRATION && VME_TAMEX==true) Fill_Plastic_VME_Histos();
    break;
    case 3:
    Fill_FATIMA_Histos();
    break;
    case 4:
    Fill_GALILEO_Histos();
    break;
    case 5:
    //Fill_Finger_Histos();
    break;
    default:
    cerr << "PrcID_Conv " << PrcID_Conv << " not known" << endl;
    exit(0);
  }
}


//-----------------------------------------------------------------------------------------------------------------------------//
void EventUnpackProc::ResetMultiplexer()
{
  for (int i = 0; i < 12; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      adcLastTimestamp[i][j] = 0;
      adcCounts[i][j] = 0;
    }
  }
}



void EventUnpackProc::CorrectTimeForMultiplexer(AidaEvent &evt)
{
  int fee = evt.Module;
  int adc = evt.Channel / 16;
  int64_t time = evt.Time;

  if ((time - adcLastTimestamp[fee][adc] > 2500) && adcLastTimestamp[fee][adc] != 0)
  adcCounts[fee][adc] = 0;

  adcLastTimestamp[fee][adc] = time;

  evt.Time = time - (2000 * adcCounts[fee][adc]++);
}

//-----------------------------------------------------------------------------------------------------------------------------//


void EventUnpackProc::load_PrcID_File(){
  ifstream data("Configuration_Files/PrcID_to_Det_Sys.txt");
  if(data.fail()){
    cerr << "Could not find PrcID config file!" << endl;
    exit(0);
  }
  int id[5] = {0,0,0,0,0};
  int i = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d %d %d %d %d",s_tmp,&id[0],&id[1],&id[2],&id[3],&id[4]);
    for(int j = 0; j < 5; ++j) PrcID_Array[i][j] = id[j];
    i++;
  }
}
//---------------------------------------------------------------------------------------------------
void EventUnpackProc::load_FingerID_File(){

  const char* format = "%d %d %d";
  ifstream data("Configuration_Files/Finger_allocation.txt");
  if(data.fail()){
    cerr << "Could not find Finger_allocation config file!" << endl;
    exit(0);
  }
  //     int id[5] = {0,0,0,0,0};
  //int i = 0;
  int tamid = 0;
  int tamch = 0;
  int fingid = 0;
  string line;
  //char s_tmp[100];
  while(data.good()){

    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),format,&tamid,&tamch,&fingid);
    fingID[tamid][tamch] = fingid;
  }
}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::read_setup_parameters(){

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
  file >> WHITE_RABBIT_USED;//dummy_var;


  cout<<endl;
  cout<<endl;
  cout<<"////////////////////////////////////////////////////////////////////////"<<endl;
    cout<<"Setup Parameters List Unpack Proc: "<<endl;
  if(WHITE_RABBIT_USED) cout<<"White Rabbit: Enabled"<<endl;
  else if(!WHITE_RABBIT_USED) cout<<"White Rabbit: Disabled"<<endl;
   // if(FAT_exclusion_dist > 0) cout<<"FATIMA Detectors Excluded if Linear Difference Exceeds "<<FAT_exclusion_dist<<" mm"<<endl;
   // else if(FAT_exclusion_dist == 0) cout<<"'Nearest Neighbour Exclusion': Disabled (Distance set to 0)"<<endl;
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

//-----------------------------------------------------------------------------------------------------------------------------//
Int_t EventUnpackProc::get_Conversion(Int_t PrcID){

  for(int i = 0;i < 6;++i){
    for(int j = 0;j < 5;++j){
      if(PrcID == PrcID_Array[i][j]) return i;
    }
  }
  cerr << "ProcID " << PrcID << " not known!" << endl;
  exit(0);

}

void EventUnpackProc::get_used_Systems(){
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

//     cout << "\n=====================================================" << endl;
//     cout << "USED SYSTEMS" << endl;
//     cout << "-----------------------------------------------------" << endl;
//     for(int j = 0;j < 6;++j){
//         if(Used_Systems[j]) cout << DET_NAME[j] << endl;
//     }
//     cout << "=====================================================" << endl;


}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::get_WR_Config(){
  ifstream data("Configuration_Files/White_Rabbit.txt");
  if(data.fail()){
    cerr << "Could not find White_Rabbit config file!" << endl;
    exit(0);
  }

  int id = 0;
  string line;
  char s_tmp[100];
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),"%s %d",s_tmp,&id);
    WR_used = (id == 1);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------//
//Implement in future
void EventUnpackProc::get_interest_arrays(){
  //SKIP_EVT_BUILDING = false; //this causes the code to crash
  amount_interest = 0;

  length_interest = new int[10];
  interest_array = new int*[10];
  for(int i = 0;i < 10;++i){
    length_interest[i] = 0;
    interest_array[i] = new int[6];
    for(int j = 0;j < 6;++j) interest_array[i][j] = -1;
  }

  string DET_NAME[6] = {"FRS","AIDA","PLASTIC","FATIMA","GALILEO","FINGER"};

  const char* format = "%d %d %d %d %d %d";

  int tmp_values[6];

  ifstream data("Configuration_Files/Coincidences_of_Interest.txt");

  //Print statements if loading of data fails
  if(data.fail()){
    string input_string;
    cerr << endl;
    cerr << "No Detector_System coincidence file found!" << endl;
    cerr << "Do you want to omit the Time_EventBuilding?  (y/n)\t\t ";
    getline(cin,input_string);
    cerr << endl;
    if(input_string == "y"){
      cerr << "------------------------------------------------------" << endl;
      cerr << "!Time_EventBuilding will be skipped! ONLY SINGLES" << endl;
      cerr << "------------------------------------------------------" << endl;
      //SKIP_EVT_BUILDING = true;
      return;
    }
    else{
      cerr << "Not skipping Time_EventBuilding.\n";
      cerr << "PLEASE CREATE Configuration_Files/Coincidences_of_Interest.txt file!" << endl;
      cerr << "\nEXITING PROGRAM NOW" << endl;
      exit(0);
    }
  }

  string line;

  cout << "\n=====================================================" << endl;
  cout << "Coincidences of interest are: " << endl;
  cout << "-----------------------------------------------------" << endl;
  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    for(int i = 0;i < 6;++i) tmp_values[i] = -1;

    sscanf(line.c_str(),format,&tmp_values[0],&tmp_values[1]
    ,&tmp_values[2],&tmp_values[3]
    ,&tmp_values[4],&tmp_values[5]);

    if(tmp_values[0] != -1){
      if(tmp_values[1] == -1){
        cerr << endl;
        cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cerr << "SINGLE COINCIDENCE DETECTED! => EDIT Coincidences_of_Interest.txt FILE!" << endl;
        cerr << "-> EXITING PROGRAM" << endl;
        cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cerr << endl;
        exit(0);
      }

      for(int i = 0;i < 6;++i){
        if(tmp_values[i] != -1){
          interest_array[amount_interest][length_interest[amount_interest]] = tmp_values[i];
          length_interest[amount_interest]++;
          cout << "interest_array " << interest_array[amount_interest][length_interest[amount_interest]] << " amount_interest "<< amount_interest <<" length_interest "<<length_interest << endl;
        }
      }
      cout << "-> ";
      for(int i = 0;i < length_interest[amount_interest]-1;++i) cout << DET_NAME[interest_array[amount_interest][i]] << " + ";
      cout << DET_NAME[interest_array[amount_interest][length_interest[amount_interest]-1]];
      cout << endl;
      amount_interest++;

    }
  }
  cout << "=====================================================" << endl;
}


  //-----------------------------------------------------------------------------------------------------------------------------//
  // ################################################################## //
  // ################################################################## //
  // ################# Raw Histogram Filling Section ################## //
  // ################################################################## //
  // ################################################################## //
  /**----------------------------------------------------------------------------------------------**/
  /**---------------------------------------------  FRS  ------------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/

  void EventUnpackProc::Make_FRS_Histos(){


    hsci_tofll2 = MakeTH1('D',"FRS/hsci_tofll2","hsci_tofll2",1500,0.,62000.); //From TAC
    hsci_tofll3 = MakeTH1('D',"FRS/hsci_tofll3","hsci_tofll3",1500,0.,62000.);
    hsci_tof2 = MakeTH1('D',"FRS/hsci_tof2","hsci_tof2",1000,0.,62000.);
    hsci_tofrr2 = MakeTH1('D',"FRS/hsci_tofrr2","hsci_tofrr2",1500,0.,62000.);
    hsci_tofrr3 = MakeTH1('D',"FRS/hsci_tofrr3","hsci_tofrr3",1500,0.,62000.);
    hsci_tof3 = MakeTH1('D',"FRS/hsci_tof3","hsci_tof3",1000,0.,62000.);

    hsci_dt_21l_21r = MakeTH1('D',"FRS/hsci_dt_21l_21r","hsci_dt_21l_21r",5001,0,5000);
    hsci_dt_41l_41r = MakeTH1('D',"FRS/hsci_dt_41l_41r","hsci_dt_41l_41r",5001,0,5000);
    hsci_dt_42l_42r = MakeTH1('D',"FRS/hsci_dt_42l_42r","hsci_dt_42l_42r",5001,0,5000);

    hsci_dt_21l_41l = MakeTH1('D',"FRS/hsci_dt_21l_41l","hsci_dt_21l_41l",5001,0,5000); //from Multihit TDCS
    hsci_dt_21r_41r = MakeTH1('D',"FRS/hsci_dt_21r_41r","hsci_dt_21r_41r",5001,0,5000);

    hsci_dt_21l_42l = MakeTH1('D',"FRS/hsci_dt_21l_42l","hsci_dt_21l_42l",5001,0,5000);
    hsci_dt_21r_42r = MakeTH1('D',"FRS/hsci_dt_21r_42r","hsci_dt_21r_42r",5001,0,5000);


    hID_x2 = MakeTH1('D',"ID_x2","ID_x2",3000,0,5000);
    hID_y2 = MakeTH1('D',"ID_y2","ID_y2",3000,0,5000);
    hID_a2 = MakeTH1('D',"ID_a2","ID_a2",3000,0,5000);
    hID_b2 = MakeTH1('D',"ID_b2","ID_b2",3000,0,5000);

    hID_x4 = MakeTH1('D',"ID_x4","ID_x4",3000,0,5000);
    hID_y4 = MakeTH1('D',"ID_y4","ID_y4",3000,0,5000);
    hID_a4 = MakeTH1('D',"ID_a4","ID_a4",3000,0,5000);
    hID_b4 = MakeTH1('D',"ID_b4","ID_b4",3000,0,5000);

    hbeta = MakeTH1('D',"FRS/beta","beta",100,0.,1.); //velocity
    hbeta3 = MakeTH1('D',"FRS/beta3","beta3",100,0.,1.);
    hgamma = MakeTH1('D',"FRS/gamma","gamma",100,0.,1.);

    hAoQ = MakeTH1('D',"FRS/AoQ","AoQ",200,1.4,3.0); // 200,1.4,3.0
    hAoQ_corr = MakeTH1('D',"FRS/AoQ_corr","AoQ_corr",200,1.4,3.0); // 200,1.4,3.0

    //We neeed MUSIC 41 and 42 check.
    hz = MakeTH1('D',"FRS/z","z",100,0.,93.); 
    hz2 = MakeTH1('D',"FRS/z2","z2",100,0.,93.);
    hz3 = MakeTH1('D',"FRS/z3","z3",100,10.,93.);

    htimestamp = MakeTH1('D',"FRS/timestamp","timestamp",30,0.,300.);
    hts = MakeTH1('D',"FRS/ts","ts",30,0.,300.);
    hts2 = MakeTH1('D',"FRS/ts2","ts2",30,0.,300.);

  }
  //-----------------------------------------------------------------------------------------------------------------------------//
  void EventUnpackProc::Fill_FRS_Histos(){

    Float_t dE, dE_cor;

    Float_t sci_l, sci_r, sci_e, sci_tx, sci_x;

    Float_t sci_tofll2, sci_tofll3, sci_tof2, sci_tofrr2, sci_tofrr3, sci_tof3;

    Float_t ID_x2, ID_y2, ID_a2, ID_b2;

    Float_t ID_x4, ID_y4, ID_a4, ID_b4;


    Int_t sci_dt_21l_21r, sci_dt_41l_41r, sci_dt_42l_42r, sci_dt_43l_43r;

    Int_t sci_dt_21l_41l, sci_dt_21r_41r, sci_dt_21l_42l, sci_dt_21r_42r;

    Float_t ID_brho, ID_rho;

    Float_t beta, beta3, gamma;

    Float_t AoQ, AoQ_corr;

    Float_t z, z2, z3;

    Float_t timestamp, ts, ts2;
    
    for(int i =0; i<2; ++i){
        dE = RAW->get_FRS_dE(i);
        dE_cor = RAW->get_FRS_dE_corr(i);
    }
    
    for(int l=0;l<12;++l){
        sci_l = RAW->get_FRS_sci_l(l);
        sci_r = RAW->get_FRS_sci_r(l);
        sci_e = RAW->get_FRS_sci_e(l);
        sci_tx = RAW->get_FRS_sci_tx(l);
        sci_x = RAW->get_FRS_sci_x(l);
    }
    sci_tofll2 = RAW->get_FRS_tofll2();
    sci_tofll3 = RAW->get_FRS_tofll3();
    sci_tof2 = RAW->get_FRS_tof2();
    sci_tofrr2 = RAW->get_FRS_tofrr2();
    sci_tofrr3 = RAW->get_FRS_tofrr3();
    sci_tof3 = RAW->get_FRS_tof3();

    ID_x2 = RAW->get_FRS_x2();
    ID_y2 = RAW->get_FRS_y2();
    ID_a2 = RAW->get_FRS_a2();
    ID_b2 = RAW->get_FRS_b2();

    ID_x4 = RAW->get_FRS_x4();
    ID_y4 = RAW->get_FRS_y4();
    ID_a4 = RAW->get_FRS_a4();
    ID_b4 = RAW->get_FRS_b4();

    sci_dt_21l_21r = RAW->get_FRS_dt_21l_21r();
    sci_dt_41l_41r = RAW->get_FRS_dt_41l_41r();
    sci_dt_42l_42r = RAW->get_FRS_dt_42l_42r();
    sci_dt_43l_43r = RAW->get_FRS_dt_43l_43r();

    sci_dt_21l_41l = RAW->get_FRS_dt_21l_41l();
    sci_dt_21r_41r = RAW->get_FRS_dt_21r_41r();

    sci_dt_21l_42l = RAW->get_FRS_dt_21l_42l();
    sci_dt_21r_42r = RAW->get_FRS_dt_21r_42r();
    
    for(int k =0; k<2; ++k){
        ID_brho = RAW->get_FRS_brho(k);
        ID_rho = RAW->get_FRS_rho(k);
    }
    beta = RAW->get_FRS_beta();
    beta3 = RAW->get_FRS_beta3();
    gamma = RAW->get_FRS_gamma();

    AoQ = RAW->get_FRS_AoQ();
    AoQ_corr = RAW->get_FRS_AoQ_corr();

    z = RAW->get_FRS_z();
    z2 = RAW->get_FRS_z2();
    z3 = RAW->get_FRS_z3();

    timestamp = RAW->get_FRS_timestamp();
    ts = RAW->get_FRS_ts();
    ts2 = RAW->get_FRS_ts2();

    // ---------------------------------------- //

    if(sci_tofll2) hsci_tofll2->Fill(sci_tofll2);
    if(sci_tofll3) hsci_tofll3->Fill(sci_tofll3);
    if(sci_tof2) hsci_tof2->Fill(sci_tof2);
    if(sci_tofrr2) hsci_tofrr2->Fill(sci_tofrr2);
    if(sci_tofrr3) hsci_tofrr3->Fill(sci_tofrr3);
    if(sci_tof3) hsci_tof3->Fill(sci_tof3);

    if(ID_x2) hID_x2->Fill(ID_x2);
    if(ID_y2) hID_y2->Fill(ID_y2);
    if(ID_a2) hID_a2->Fill(ID_a2);
    if(ID_b2) hID_b2->Fill(ID_b2);

    if(ID_x4) hID_x4->Fill(ID_x4);
    if(ID_y4) hID_y4->Fill(ID_y4);
    if(ID_a4) hID_a4->Fill(ID_a4);
    if(ID_b4) hID_b4->Fill(ID_b4);

    if(sci_dt_21l_21r) hsci_dt_21l_21r->Fill(sci_dt_21l_21r);
    if(sci_dt_41l_41r) hsci_dt_41l_41r->Fill(sci_dt_41l_41r);
    if(sci_dt_42l_42r) hsci_dt_42l_42r->Fill(sci_dt_42l_42r);

    if(sci_dt_21l_41l) hsci_dt_21l_41l->Fill(sci_dt_21l_41l);
    if(sci_dt_21r_41r) hsci_dt_21r_41r->Fill(sci_dt_21r_41r);

    if(sci_dt_21l_42l) hsci_dt_21l_42l->Fill(sci_dt_21l_42l);
    if(sci_dt_21r_42r) hsci_dt_21r_42r->Fill(sci_dt_21r_42r);



    if(beta) hbeta->Fill(beta);
    if(beta3) hbeta3->Fill(beta3);
    if(gamma) hgamma->Fill(gamma);

    if(AoQ) hAoQ->Fill(AoQ);
    if(AoQ_corr) hAoQ_corr->Fill(AoQ_corr);

    if(z) hz->Fill(z);
    if(z2) hz2->Fill(z2);
    if(z3) hz3->Fill(z3);

    if(timestamp) htimestamp->Fill(timestamp);
    if(ts) hts->Fill(ts);
    if(ts2) hts2->Fill(ts2);

  }

  /**----------------------------------------------------------------------------------------------**/
  /**-------------------------------------------  AIDA   ------------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/

  void EventUnpackProc::Make_AIDA_Histos(){

    TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
    hAIDA_ADC.resize(conf->FEEs());

    for (int i = 0; i < conf->FEEs(); i++)
    {
      for (int j = 0; j < 64; j++)
      {
        hAIDA_ADC[i][j][0] = MakeTH1('I',
          Form("AIDA/Unpacker/FEE%d/Fee%d_L_Channel%02d", i+1, i+1, j+1),
          Form("FEE %d Channel %2d (Low Energy)", i+1, j+1),
          2000, -32768, 32767
        );
      }
    }

    for (int i = 0; i < conf->FEEs(); i++)
    {
      for (int j = 0; j < 64; j++)
      {
        hAIDA_ADC[i][j][1] = MakeTH1('I',
          Form("AIDA/Unpacker/FEE%d/Fee%d_H_Channel%02d", i+1, i+1, j+1),
          Form("FEE %d Channel %2d (High Energy)", i+1, j+1),
          2000, -32768, 32767
        );
      }
    }
}

void EventUnpackProc::Fill_AIDA_Histos() {
  AIDA_Hits = RAW->get_AIDA_HITS();

  for(int i = 0; i<AIDA_Hits; i++) {
    int fee = RAW-> get_AIDA_FEE_ID(i);
    int chan = RAW-> get_AIDA_CHA_ID(i);
    int adc = RAW->get_AIDA_ADC(i);
    int veto = RAW->get_AIDA_HighE_VETO(i) ? 1 : 0;

    hAIDA_ADC[fee][chan][veto]->Fill(adc - 32767);
  }
}


  /**----------------------------------------------------------------------------------------------**/
  /**---------------------------------------  bPLASTIC TAMEX  -------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/
  void EventUnpackProc::Make_Plastic_Histos(){

    TOT_TOT = new TH1***[100];
    TOT_Single = new TH1**[100];
    TRAIL_TRAIL = new TH1***[100];
    LEAD_LEAD = new TH1***[100];

    for(int i = 0;i < 100;++i){
      TOT_Single[i] = new TH1*[100];
      TOT_TOT[i] = new TH1**[100];
      TRAIL_TRAIL[i] = new TH1**[100];
      LEAD_LEAD[i] = new TH1**[100];

      for(int j = 0;j < 100;++j){
        TOT_TOT[i][j] = new TH1*[100];
        TRAIL_TRAIL[i][j] = new TH1*[100];
        LEAD_LEAD[i][j] = new TH1*[100];

        for(int k = 0;k < 100;++k){
          TOT_TOT[i][j][k] = nullptr;
          TRAIL_TRAIL[i][j][k] = nullptr;
          LEAD_LEAD[i][j][k]   = nullptr;
        }

        TOT_Single[i][j] = nullptr;

      }
    }
  }
  //-----------------------------------------------------------------------------------------------------------------------------//

  void EventUnpackProc::Fill_Plastic_Histos(){

    //get amount of fired Tamex modules
    int TamexHits = RAW->get_PLASTIC_tamex_hits();

    //int Physical_hits = 0;
    int leadHits = 0,leadHitsCh = 0;
    int trailHits = 0,trailHitsCh = 0;
    int Phys_Channel[2] = {0,0};
    double Lead[2] = {0,0};
    double Trail[2] = {0,0};
    double TOT[2] = {0,0};

    double Diff = 0;

    int MaxHits = 0;

    for(int i = 0;i < TamexHits;++i){

      leadHits = RAW->get_PLASTIC_lead_hits(i);
      trailHits = RAW->get_PLASTIC_trail_hits(i);
      //cout << "trail " <<trailHits << endl;
      MaxHits = (leadHits >= trailHits) ? leadHits : trailHits;

      for(int j = 0;j < leadHits;++j){
        Phys_Channel[0] = RAW->get_PLASTIC_physical_channel(i,j);
        Lead[0] = RAW->get_PLASTIC_lead_T(i,Phys_Channel[0]);

        //Leading - Leading
        for(int k = 0;k < leadHits;++k){
          if(k != j){
            Phys_Channel[1] = RAW->get_PLASTIC_physical_channel(i,k);
            Lead[1] = RAW->get_PLASTIC_lead_T(i,Phys_Channel[1]);
            Diff = Lead[0] - Lead[1];

            if(!LEAD_LEAD[i][Phys_Channel[0]][Phys_Channel[1]]){
              LEAD_LEAD[i][Phys_Channel[0]][Phys_Channel[1]] = MakeTH1('D',
              Form("PLASTIC/lead_minus_lead_all_chans/lead_minus_lead_board_%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),
              Form("lead_minus_lead_board_%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),10000, -500., 500.);
            }
            LEAD_LEAD[i][Phys_Channel[0]][Phys_Channel[1]]->Fill(Diff);
          }
        }
      }
      for(int j = 0;j < trailHits;++j){
        Phys_Channel[0] = RAW->get_PLASTIC_physical_channel(i,j);
        Trail[0] = RAW->get_PLASTIC_trail_T(i,Phys_Channel[0]);

        //Trailing - Trailing
        for(int k = 0;k < trailHits;++k){
          if(k != j){
            Phys_Channel[1] = RAW->get_PLASTIC_physical_channel(i,k);
            Trail[1] = RAW->get_PLASTIC_trail_T(i,Phys_Channel[1]);
            Diff = Trail[0] - Trail[1];

            if(!TRAIL_TRAIL[i][Phys_Channel[0]][Phys_Channel[1]]){
              TRAIL_TRAIL[i][Phys_Channel[0]][Phys_Channel[1]] = MakeTH1('D',
              Form("PLASTIC/trail_minus_trail/trail_minus_trail_board_%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),
              Form("trail_minus_trail_board%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),10000, -500., 500.);
            }

            TRAIL_TRAIL[i][Phys_Channel[0]][Phys_Channel[1]]->Fill(Diff);
          }
        }
      }
      for(int j = 0;j < MaxHits;++j){

        Phys_Channel[0] = RAW->get_PLASTIC_physical_channel(i,j);

        leadHitsCh = RAW->get_PLASTIC_physical_lead_hits(i,Phys_Channel[0]);
        trailHitsCh = RAW->get_PLASTIC_physical_trail_hits(i,Phys_Channel[0]);

        if(leadHitsCh == trailHitsCh){
          TOT[0] = RAW->get_PLASTIC_TOT(i,Phys_Channel[0]);
          //Trailing - Trailing
          for(int k = 0;k < MaxHits;++k){
            if(k != j){
              Phys_Channel[1] = RAW->get_PLASTIC_physical_channel(i,k);

              leadHitsCh = RAW->get_PLASTIC_physical_lead_hits(i,Phys_Channel[1]);
              trailHitsCh = RAW->get_PLASTIC_physical_trail_hits(i,Phys_Channel[1]);

              if(leadHitsCh == trailHitsCh){
                TOT[1] = RAW->get_PLASTIC_TOT(i,Phys_Channel[1]);
                Diff = TOT[0] - TOT[1];

                if(!TOT_TOT[i][Phys_Channel[0]][Phys_Channel[1]]){
                  TOT_TOT[i][Phys_Channel[0]][Phys_Channel[1]] = MakeTH1('D',
                  Form("PLASTIC/TOT/TOT_Diffs/TOT_board_%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),
                  Form("TOT_board%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),10000, -500., 500.);
                }

                TOT_TOT[i][Phys_Channel[0]][Phys_Channel[1]]->Fill(Diff);
              }
            }
          }
          if(!TOT_Single[i][Phys_Channel[0]]){
            TOT_Single[i][Phys_Channel[0]] = MakeTH1('D',Form("PLASTIC/TOT/TOTs/TOT_board_%d_ch%d",i,Phys_Channel[0]),
            Form("TOT_board%d_ch%d",i,Phys_Channel[0]),10000, -500., 500.);
          }
          TOT_Single[i][Phys_Channel[0]]->Fill(TOT[0]);
        }
      }
    }
  }
  /**----------------------------------------------------------------------------------------------**/
  /**--------------------------------------  bPLASTIC VME (+Scalar)  ----------------------------------------**/
  /**----------------------------------------------------------------------------------------------**/
  void EventUnpackProc::Make_Plastic_VME_Histos(){
    for (int i = 0; i<32; i++){

  hPLAS_QDCRaw1[i] =  MakeTH1('D', Form("bPlastic/Energy/Raw/QDC1Raw/QDC1_Ch.%2d",i), Form("QDC1 Ch. %2d",i), 2500, 0., 5000.);
  hPLAS_QDCRaw2[i] =  MakeTH1('D', Form("bPlastic/Energy/Raw/QDC2Raw/QDC2_Ch.%2d",i), Form("QDC2 Ch. %2d",i), 2500, 0., 5000.);
  hPLAS_TDCRaw[i] =  MakeTH1('D', Form("bPlastic/Timing/Raw/TDCRaw/TDC_Ch.%2d",i), Form("TDC Ch. %2d",i), 1000, 0, 10000);
    }
    hScalar_hit_pattern = MakeTH1('D',"Scalar/HitPat","Scalar Hit pattern",32,0,32);
}


  //-----------------------------------------------------------------------------------------------------------------------------//
  void EventUnpackProc::Fill_Plastic_VME_Histos(){

    double QDC1[32], QDC2[32];
    double bplasTDC[32];
    int bPlas_QCDID;
    int bplasTDCID, bPlasIT;

    int bPlasmulti[32];
    int Scalar_iterator;
    int Scalar_Chan;
   // double Scalar_Data[32];

    //Reset arrays and variables://
    Scalar_Chan = 0;
    bplasTDCID=0;
    bPlasIT = 0;

    for (int i =0; i<32; i++){
      QDC1[i] = 0;
      QDC2[i] = 0;
      bplasTDC[i] = 0;
    //  Scalar_Data[i] = 0;
    
    }

    /**------------------PLASTIC Energy -----------------------------------------**/

    for (int jj = 0; jj < RAW->get_plastic_VME_QDC_fired();  jj++){
      bPlas_QCDID = RAW->get_plastic_VME_QDC_cha(jj);
      QDC1[bPlas_QCDID] = RAW->get_plastic_VME_QDC_dat1(jj);
      hPLAS_QDCRaw1[bPlas_QCDID]->Fill(QDC1[bPlas_QCDID]);
      hPLAS_QDCRaw2[bPlas_QCDID]->Fill(QDC2[bPlas_QCDID]);

    
    }
    /**------------------PLASTIC Time -----------------------------------------**/
   
    bPlasIT = RAW->get_plastic_VME_TDC_fired(); //Fired TDC channels

    for (int i=0; i<bPlasIT; i++){

      bplasTDCID = RAW->get_plastic_VME_TDC_cha(i); //Plastic TDC ID
      bplasTDC[bplasTDCID] = RAW->get_plastic_VME_TDC_dat(bplasTDCID)*0.025;//1ns //Plastic TDC Data
     
      bPlasmulti[bplasTDCID]++;


    if(bplasTDC[bplasTDCID]>0){
     
            hPLAS_TDCRaw[bplasTDCID]->Fill(bplasTDC[bplasTDCID]);  //1ns

             }  
            }
              //Scalar
        Scalar_iterator = RAW->get_scalar_iterator();
        for (int g=0; g<Scalar_iterator; g++){
            Scalar_Chan = RAW->get_scalar_chan(g);
            hScalar_hit_pattern->Fill(Scalar_Chan);
            }
         }

/**----------------------------------------------------------------------------------------------**/
/**-----------------------------------------  FATIMA   ------------------------------------------**/
/**----------------------------------------------------------------------------------------------**/
void EventUnpackProc::Make_FATIMA_Histos(){


  for (int i=0; i<FAT_MAX_DET; i++){
    hFAT_Eraw[i] = MakeTH1('D', Form("FATIMA/Energy/Raw/E_Raw_LaBr%02d", i),
    Form("LaBr%02d energy (raw)", i),2000,0,40000);

    hFAT_Traw[i] = MakeTH1('D', Form("FATIMA/Timing/Raw/Traw_LaBr%02d", i),
                    Form("LaBr%02d energy", i),5000,0,5000);

            }
        }
//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::Fill_FATIMA_Histos(){
  double FAT_E[50],FAT_T[50]; 
  double En_i;  
  int detQDC, detTDC;
  detTDC = 0;
  detQDC = 0;
  for (int k=0; k<50; k++){
    FAT_T[k] = 0;
    FAT_E[k] = 0;
  }



  /**------------------FATIMA Energy -----------------------------------------**/
  
  for (int i=0; i<RAW->get_FAT_QDCs_fired(); i++){ /** Loops over only channels in the QDC **/

    detQDC = RAW->get_FAT_QDC_id(i); /**FAT ID QDC*/
    En_i = RAW->get_FAT_QLong_Raw(i); /**Raw FAT Energy*/

    hFAT_Eraw[detQDC]->Fill(En_i);


  }
  /**------------------FATIMA TIMING -----------------------------------------**/
  for (int i=0; i<RAW->get_FAT_TDCs_fired(); i++){ /** Loops over only channels in the TDC 1-4 **/

    detTDC = (RAW->get_FAT_TDC_id(i));
    FAT_T[detTDC] = (RAW->get_FAT_TDC_timestamp(i));
         hFAT_Traw[detTDC]->Fill(FAT_T[detTDC]*0.025); //in ns

         }
}



/**----------------------------------------------------------------------------------------------**/
/**----------------------------------------   GALILEO   -----------------------------------------**/
/**----------------------------------------------------------------------------------------------**/


void EventUnpackProc::Make_GALILEO_Histos(){
  for (int j; j<32; j++){
        hGAL_Raw_E[j] = MakeTH1('D',Form("GALILEO/Raw/GALILEO_Energy_Spectra/GALILEO_Raw_E%2d",j),
                            Form("GALILEO Channel Energy Channel Raw %2d",j),20000,0,20000);

                    }
                }
//-----------------------------------------------------------------------------------------------------------------------------//
void EventUnpackProc::Fill_GALILEO_Histos(){

    double tmpGAL[32];
    int  GALILEO_hits, GalID;

     /**------------------GALILEO Raw Energy -----------------------------------------**/
      GALILEO_hits = RAW->get_GALILEO_am_Fired();
         for(int i=0; i<GALILEO_hits; i++){
        GalID = RAW->get_GALILEO_Det_ids(i);
        tmpGAL[GalID] = RAW->get_GALILEO_Chan_E(i)/1000;
        hGAL_Raw_E[GalID]->Fill(tmpGAL[GalID]);
         }
   }



//         /**----------------------------------------------------------------------------------------------**/
//         /**---------------------------  FINGER TAMEX ------------------------ ---------------------------**/
//         /**----------------------------------------------------------------------------------------------**/
// void EventUnpackProc::Make_Finger_Histos(){
// 
//   for( int i=0;i< 48;i++){
// 
//     // hlead_lead[i] = MakeTH1('D',Form("FINGER/lead-lead/lead-leadCh.%02d",i),Form("lead-leadCh.%02d",i), 200001, -2500., 2500.);
// 
//     //hFin_PadiCoarse[i] = MakeTH1('D',Form("FINGER/PADI_Coarse/CoarseCh%02d",i),Form("Coarse%2d",i), 1001, -500., 500.);
//     // hFin_PadiFine[i] = MakeTH1('D',Form("FINGER/PADI_Fine/FineCh%02d",i),Form("Fine%2d",i), 1001, -500., 500.);
//     //hFin_TimeDiffTrailing[i] = MakeTH1('D',Form("FINGER/TimeDiffTrailing/TrailingDiffCh%02d",i),Form("trailing Ch(%2d)-Ch(%2d)",i,i+1), 20001, -100., 100.);
//     // hFin_ToT[i] = MakeTH1('D',Form("FINGER/TOT/TOTCh%02d",i),Form("TOT Ch(%2d)-Ch(%2d)",i,i+1), 800, 0., 200.);
//   }
//   TOT_TOT = new TH1***[100];
//   TOT_Single = new TH1**[100];
//   TRAIL_TRAIL = new TH1***[100];
//   LEAD_LEAD = new TH1***[100];
// 
//   for(int i = 0;i < 100;++i){
//     TOT_Single[i] = new TH1*[100];
//     TOT_TOT[i] = new TH1**[100];
//     TRAIL_TRAIL[i] = new TH1**[100];
//     LEAD_LEAD[i] = new TH1**[100];
// 
//     for(int j = 0;j < 100;++j){
//       TOT_TOT[i][j] = new TH1*[100];
//       TRAIL_TRAIL[i][j] = new TH1*[100];
//       LEAD_LEAD[i][j] = new TH1*[100];
// 
//       for(int k = 0;k < 100;++k){
//         TOT_TOT[i][j][k] = nullptr;
//         TRAIL_TRAIL[i][j][k] = nullptr;
//         LEAD_LEAD[i][j][k]   = nullptr;
//       }
// 
//       TOT_Single[i][j] = nullptr;
// 
//     }
//   }
// }
// //-----------------------------------------------------------------------------------------------------------------------------//
// 
// void EventUnpackProc::Fill_Finger_Histos(){
//   if(event_number>0){
//     //get amount of fired Tamex modules
//     int TamexHits = RAW->get_FINGER_tamex_hits();
// 
//     //int Physical_hits = 0;
//     int leadHits[3]={0,0};
//     int trailHits = 0;
//     // int Phys_Channel[2] = {0,0};
//     int Phys_chan[100][100]={0,0};
//     int chID[100][100]={0,0};
//     bool TamexLeadisFired;
//     double Coarse_T_lead[100][100]={0,0};
//     //     double Lead[2] = {0,0};
//     //     double Trail[2] = {0,0};
//     //double TOT[2] = {0,0};
//     double Fin_ToT[100][100]={0,0};
//     double LeadDiff[48]={0,0};
//     int MaxHits = 0;
// 
//     //  cout << "TamexHits " << TamexHits << endl;
//     for(int i = 0;i < TamexHits;i++){
//       //  cout << "Fin Trig " << RAW->get_FINGER_trigger_T(i) << endl;
// 
//       leadHits[i] = RAW->get_FINGER_lead_hits(i);
//       if(leadHits[i]>-1)TamexLeadisFired = true;
// 
//       //  cout << "iterator " <<  RAW->get_FINGER_am_Fired(i) << endl;
//       trailHits = RAW->get_FINGER_trail_hits(i);
//       // cout << "0) TamexHits "<< TamexHits << " leadHits " << leadHits <<endl;
//       //cout << "trail " <<trailHits << endl;
//       MaxHits = (leadHits[i] >= trailHits) ? leadHits[i] : trailHits;
// 
//       for(int j = 0;j < RAW->get_FINGER_am_Fired(i);j++){
//         chID[i][j] = RAW->get_FINGER_CH_ID(i,j);
// 
//         Phys_chan[i][j] = (RAW->get_FINGER_physical_channel(i,j));
// 
// 
// 
//         // cout << "4)ev " << event_number << " Phys_chan[i][j] " << Phys_chan[i][j] << "  chID[i][j] " <<  chID[i][j] <<" Coarse_T_lead[i][j] " <<Coarse_T_lead[i][j]<< " i " << i << " j " << j << endl;
// 
//         //    cout << "a) i " << i << " j " << j << " Phys_chan[i][j] " << Phys_chan[i][j] << "  chID[i][j] " <<  chID[i][j] << " Coarse_T_lead " << Coarse_T_lead[i][j] << endl;
// 
//         Fin_ToT[i][j] = RAW->get_FINGER_TOT(i,j);
//         // cout << "1) Phys_chan[i][j] "<< Phys_chan[i][j] << " Coarse_T_lead[i][j] " << Coarse_T_lead[i][j] <<" i " << i << " j " <<  j << " chID[i][j]  " <<chID[i][j] << endl;
// 
// 
// 
//         // Phys_Channel[0] = RAW->get_FINGER_physical_channel(i,j);
// 
//         //Lead[0] = RAW->get_FINGER_lead_T(i,Phys_Channel[0]);
//         //  cout <<  "Phys_chan[i][j] "<<  Phys_chan[i][j] << " Coarse_T_lead[i][j] " << Coarse_T_lead[i][j] <<" i " << i << " j " << j <<  endl;
// 
//       }
//     }
// 
//     ///SET PADI or PADIWA!////
//     if (PADI_OR_PADIWA == true){
//       //  if(TamexLeadisFired==true){
// 
//       for(int i = 0;i < TamexHits;i++){
//         for(int j = 0;j <  RAW->get_FINGER_am_Fired(i);j++){
//           Coarse_T_lead[i][j] = RAW->get_FINGER_lead_T(i,j);
//           //Leading - Leading
//           if(i==0) Phys_chan[i][j] = Phys_chan[i][j];
//           if(i==1) Phys_chan[i][j] = Phys_chan[i][j]+16;
// 
//           if(i==2) Phys_chan[i][j] = Phys_chan[i][j]+32;
//           if(Coarse_T_lead[i][j]!=0 && Coarse_T_lead[i][j+2]!=0){
//             LeadDiff[Phys_chan[i][j]] =  Coarse_T_lead[i][j] -  Coarse_T_lead[i][j+2];
//             // cout<<"i " << i << " j " << j<<" ChID "<<chID[i][j]<<" Phys_chan[i][j] " << Phys_chan[i][j] <<" Coarse_T_lead[i][j] " <<Coarse_T_lead[i][j] <<" Coarse_T_lead[i][j+2] " <<Coarse_T_lead[i][j+2] <<" LeadDiff[r] " << LeadDiff[Phys_chan[i][j]]<<endl;
//             //for (int l=0; l<16*3;l++){
//             // cout<<"event " << event_number <<" i " << i << " j " << j<<" ChID "<<chID[i][j]<<" Phys_chan[i][j] " << Phys_chan[i][j] <<endl;
//             //   cout <<"i " << i << " j " << j << "chID[i][j] " << chID[i][j] <<"  Phys_chan[i][j]  " <<  Phys_chan[i][j]    << " get_FINGER_physical_lead_hits " << RAW->get_FINGER_physical_lead_hits(i,j)<< " leadHits " << leadHits[i] << endl;
//             //  if (Phys_chan[i][j]==19&& LeadDiff[Phys_chan[i][j]]>120){
//             // hlead_lead[Phys_chan[i][j]]->Fill(LeadDiff[Phys_chan[i][j]]);
//             //cout<<"event " << event_number <<" i " << i << " j " << j<<" ChID "<<chID[i][j]<<" Phys_chan[i][j] " << Phys_chan[i][j] <<" Coarse_T_lead[i][j] " <<Coarse_T_lead[i][j] <<" Coarse_T_lead[i][j+2] " <<Coarse_T_lead[i][j+2] <<" LeadDiff[r] " << LeadDiff[Phys_chan[i][j]]<<endl;
// 
//             //   }
// 
//             //  hFin_ToT[Phys_chan[i][j]]->Fill(Fin_ToT[i][j]);
// 
//           }
//         }
//       }
//     }
//   }
// }
//         for(int j = 0;j < trailHits;++j){
//             Phys_Channel[0] = RAW->get_FINGER_physical_channel(i,j);
//             Trail[0] = RAW->get_FINGER_trail_T(i,Phys_Channel[0]);
//
//             //Trailing - Trailing
//             for(int k = 0;k < trailHits;++k){
//                 if(k != j){
//                     Phys_Channel[1] = RAW->get_FINGER_physical_channel(i,k);
//                     Trail[1] = RAW->get_FINGER_trail_T(i,Phys_Channel[1]);
//                     Diff = Trail[0] - Trail[1];
//
//                     if(!TRAIL_TRAIL[i][Phys_Channel[0]][Phys_Channel[1]]){
//                         TRAIL_TRAIL[i][Phys_Channel[0]][Phys_Channel[1]] = MakeTH1('D',
//                                       Form("FINGER/trail_minus_trail/trail_minus_trail_board_%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),
//                                       Form("trail_minus_trail_board%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),10000, -500., 500.);
//                     }
//
//                     TRAIL_TRAIL[i][Phys_Channel[0]][Phys_Channel[1]]->Fill(Diff);
//                 }
//             }
//         }
//}
//         for(int j = 0;j < MaxHits;++j){
//
//             Phys_Channel[0] = RAW->get_FINGER_physical_channel(i,j);
//
//             leadHitsCh = RAW->get_FINGER_physical_lead_hits(i,Phys_Channel[0]);
//             trailHitsCh = RAW->get_FINGER_physical_trail_hits(i,Phys_Channel[0]);
//
//
//             if(leadHitsCh == trailHitsCh){
//                 TOT[0] = RAW->get_FINGER_TOT(i,Phys_Channel[0]);
//                 //Trailing - Trailing
//                 for(int k = 0;k < MaxHits;++k){
//                     if(k != j){
//                         Phys_Channel[1] = RAW->get_FINGER_physical_channel(i,k);
//
//                         leadHitsCh = RAW->get_FINGER_physical_lead_hits(i,Phys_Channel[1]);
//                         trailHitsCh = RAW->get_FINGER_physical_trail_hits(i,Phys_Channel[1]);
//
//                         if(leadHitsCh == trailHitsCh){
//                             TOT[1] = RAW->get_FINGER_TOT(i,Phys_Channel[1]);
//                             Diff = TOT[0] - TOT[1];
//
//                             if(!TOT_TOT[i][Phys_Channel[0]][Phys_Channel[1]]){
//                                 TOT_TOT[i][Phys_Channel[0]][Phys_Channel[1]] = MakeTH1('D',
//                                       Form("FINGER/TOT/TOT_Diffs/TOT_board_%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),
//                                       Form("TOT_board%d_from_ch%d_to_%d",i,Phys_Channel[0],Phys_Channel[1]),10000, -500., 500.);
//                             }
//
//                             TOT_TOT[i][Phys_Channel[0]][Phys_Channel[1]]->Fill(Diff);
//                         }
//                     }
//                 }
//                 if(!TOT_Single[i][Phys_Channel[0]]){
//                     TOT_Single[i][Phys_Channel[0]] = MakeTH1('D',Form("FINGER/TOT/TOTs/TOT_board_%d_ch%d",i,Phys_Channel[0]),
//                                                      Form("TOT_board%d_ch%d",i,Phys_Channel[0]),10000, -500., 500.);
//                 }
//                 TOT_Single[i][Phys_Channel[0]]->Fill(TOT[0]);
//             }
//         }
//     }
// }

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::checkPADI_or_PADIWA(){

  std::ifstream PADIFILE("Configuration_Files/PADI_or_PADIWA.txt");

  std::string line;

  if(PADIFILE.fail()){
    std::cerr << "Could not find Configuration_Files/PADI_or_PADIWA.txt file" << std::endl;
    exit(1);
  }
  bool P_or_PW = false;
  while(std::getline(PADIFILE,line)){
    if(line[0] == '#') continue;

    if(line == "PADI") P_or_PW = true;
    if(line == "PADIWA") P_or_PW = false;

    if(line != "PADIWA" && line != "PADI"){
      std::cerr << line << " module of PLASTIC not known!" <<std::endl;
      exit(1);
    }
  }

  PADI_OR_PADIWA = P_or_PW;

}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::checkTAMEXorVME(){

  std::ifstream PL_FILE("Configuration_Files/TAMEX_or_VME.txt");

  std::string line;

  if(PL_FILE.fail()){
    std::cerr << "Could not find Configuration_Files/TAMEX_or_VME.txt file" << std::endl;
    exit(1);
  }
  bool T_or_V = false;
  while(std::getline(PL_FILE,line)){
    if(line[0] == '#') continue;

    if(line == "VME") T_or_V = true;
    if(line == "TAMEX") T_or_V = false;

    if(line != "VME" && line != "TAMEX"){
      std::cerr << line << " module of PLASTIC not known!" <<std::endl;
      exit(1);
    }
  }

  VME_TAMEX = T_or_V;

}

//-----------------------------------------------------------------------------------------------------------------------------//

bool EventUnpackProc::Check_Cal_Plastic(){
  ifstream data("Configuration_Files/PLASTIC_CALIB_FILE.txt");
  if(data.fail()){
    cerr << "Could not find Calibration type file for PLASTIC" << endl;
    exit(0);
  }
  string line;
  const char* format = "%s %d";
  char s[100];
  int val;
  bool CALIBRATE = false;

  while(data.good()){
    getline(data,line,'\n');
    if(line[0] == '#') continue;
    sscanf(line.c_str(),format,&s,&val);
    if(string(s) == string("ONLINE")) CALIBRATE = (val == 1);
  }

  return CALIBRATE;

}

//-----------------------------------------------------------------------------------------------------------------------------//

void EventUnpackProc::print_MBS(int* pdata,int lwords){
  cout << "---------------------\n";
  for(int i = 0;i < lwords;++i){
    cout << hex << *(pdata + i) << " ";
    if(i % 5 == 0 && i > 0) cout << endl;
  }
  cout << "\n---------------------\n";
}
//-----------------------------------------------------------------------------------------------------------------------------//
//                                                            END                                                              //
//-----------------------------------------------------------------------------------------------------------------------------//
