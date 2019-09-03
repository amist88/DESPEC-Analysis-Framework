#include "Raw_Event.h"
#include "AidaMover.h" 
#include <iostream>

using namespace std;

//---------------------------------------------------------------

//Raw_Event::Raw_Event(bool PADI_OR_PADIWA) : PLASTIC_Data(),PLASTIC_VME_Data(){
    Raw_Event::Raw_Event(){


    Event_Type = -1;

//     PLASTIC_Data.PADI_OR_PADIWA = PADI_OR_PADIWA;

}

//---------------------------------------------------------------

Raw_Event::~Raw_Event(){}

//---------------------------------------------------------------

// #################################################################


void Raw_Event::set_DATA_MUSIC(Float_t* FRS_dE,Float_t* FRS_dE_cor){

    for(int i; i<3; ++i){

    dE[i] = FRS_dE[i];
    dE_cor[i] = FRS_dE_cor[i];

    }

}
void Raw_Event::set_DATA_SCI(Float_t* FRS_sci_l,Float_t* FRS_sci_r,Float_t* FRS_sci_e,Float_t* FRS_sci_tx,Float_t* FRS_sci_x){

    for(int i; i<12; ++i){

    sci_l[i] = FRS_sci_l[i];
    sci_r[i] = FRS_sci_r[i];
    sci_e[i] = FRS_sci_e[i];
    sci_tx[i] = FRS_sci_tx[i];
    sci_x[i] = FRS_sci_x[i];

    }

}

void Raw_Event::set_DATA_SCI_dT(Int_t FRS_dt_21l_21r, Int_t FRS_dt_41l_41r,
                Int_t FRS_dt_21l_41l, Int_t FRS_dt_21r_41r,
                Int_t FRS_dt_42l_42r, Int_t FRS_dt_43l_43r,
                Int_t FRS_dt_21l_42l, Int_t FRS_dt_21r_42r,
                Int_t FRS_dt_81l_81r, Int_t FRS_dt_21l_81l, Int_t FRS_dt_21r_81r){


    dt_21l_21r = FRS_dt_21l_21r;
    dt_41l_41r = FRS_dt_41l_41r;
    dt_21l_41l = FRS_dt_21l_41l;
    dt_21r_41r = FRS_dt_21r_41r;
    dt_42l_42r = FRS_dt_42l_42r;
    dt_43l_43r = FRS_dt_43l_43r;
    dt_21l_42l = FRS_dt_21l_42l;
    dt_21r_42r = FRS_dt_21r_42r;
    dt_81l_81r = FRS_dt_81l_81r;
    dt_21l_81l = FRS_dt_21l_81l;
    dt_21r_81r = FRS_dt_21r_81r;




}
void Raw_Event::set_DATA_SCI_ToF(Float_t FRS_sci_tofll2,Float_t FRS_sci_tofll3,Float_t FRS_sci_tof2,Float_t FRS_sci_tofrr2,Float_t FRS_sci_tofrr3,Float_t FRS_sci_tof3){

    sci_tofll2 = FRS_sci_tofll2;
    sci_tofll3 = FRS_sci_tofll3;
    sci_tof2   = FRS_sci_tof2;
    sci_tofrr2 = FRS_sci_tofrr2;
    sci_tofrr3 = FRS_sci_tofrr3;
    sci_tof3   = FRS_sci_tof3;

}
void Raw_Event::set_DATA_ID_2_4(Float_t FRS_ID_x2,Float_t FRS_ID_y2,Float_t FRS_ID_a2,Float_t FRS_ID_b2,Float_t FRS_ID_x4,Float_t FRS_ID_y4,Float_t FRS_ID_a4,Float_t FRS_ID_b4){

    ID_x2 = FRS_ID_x2;
    ID_y2 = FRS_ID_y2;
    ID_a2 = FRS_ID_a2;
    ID_b2 = FRS_ID_b2;
    ID_x4 = FRS_ID_x4;
    ID_y4 = FRS_ID_y4;
    ID_a4 = FRS_ID_a4;
    ID_b4 = FRS_ID_b4;

}
void Raw_Event::set_DATA_ID_Beta_Rho(Float_t* FRS_ID_brho,Float_t* FRS_ID_rho,Float_t FRS_beta,Float_t FRS_beta3,Float_t FRS_gamma){

    for(int i; i<2; ++i){

    ID_brho[i] = FRS_ID_brho[i];
    ID_rho[i] = FRS_ID_rho[i];

    }

    beta = FRS_beta;
    beta3 = FRS_beta3;
    gamma = FRS_gamma;

}
void Raw_Event::set_DATA_ID_Z_AoQ(Float_t FRS_AoQ,Float_t FRS_AoQ_corr,Float_t FRS_z,Float_t FRS_z2,Float_t FRS_z3){

    AoQ = FRS_AoQ;
    AoQ_corr = FRS_AoQ_corr;
    z = FRS_z;
    z2 = FRS_z2;
    z3 = FRS_z3;

}
void Raw_Event::set_DATA_ID_Timestamp(Float_t FRS_timestamp,Float_t FRS_ts,Float_t FRS_ts2){

    timestamp = FRS_timestamp;
    ts = FRS_ts;
    ts2 = FRS_ts2;

}

// #################################################################

// int Raw_Event::get_Event_type(){
//     return Event_Type;
// }
//
//
// void Raw_Event::set_DATA_FATIMA(int QDC_FIRED,int TDC_FIRED,std::vector<double> &Ql_Raw,std::vector<double> &Qs_Raw,
//                                 std::vector<double> &Ql,std::vector<ULong64_t> &TDC,std::vector<double> &TDC_ns,
//                                 std::vector<ULong64_t> &QDC_c,std::vector<double> &QDC_f,std::vector<int> &det_ids_QDC,
//                                 std::vector<int> &det_ids_TDC){
//
//     FATIMA_Data.SetDATA(QDC_FIRED,TDC_FIRED,Ql_Raw,Qs_Raw,Ql,TDC,TDC_ns,QDC_c,QDC_f,det_ids_QDC,det_ids_TDC);
//  cout << FATIMA_Data.SetDATA(&Ql) << endl;
//     Event_Type = 3;
//
// }



int Raw_Event::get_Event_type(){
    return Event_Type;

}
//
//
void Raw_Event::set_DATA_FATIMA(int QDC_FIRED,int TDC_FIRED,
                                double* Ql_Raw,double* Qs_Raw,
                                double* Ql,
                                ULong64_t* TDC, double* TDC_ns,
                                ULong64_t* QDC_c, double* QDC_f,
                                int* det_ids_QDC,int* det_ids_TDC){

    this->FAT_QDCs_FIRED = QDC_FIRED;
    this->FAT_TDCs_FIRED = TDC_FIRED;

    int dets_fired = 0;
    for (int i=0; i<QDC_FIRED; i++) {
        this->FAT_QDC_id[i] = det_ids_QDC[i];
        //if (det_ids_QDC[i] == 35) cout<<"I am in QDC"<<endl;
        this->FAT_QLong[i]  = Ql[i];
        this->FAT_QLong_Raw[i]  = Ql_Raw[i];
        this->FAT_QShort_Raw[i] = Qs_Raw[i];
        this->FAT_QDC_t_coarse[i] = QDC_c[i];
        this->FAT_QDC_t_fine[i] = QDC_f[i];
        for (int j=0; j<TDC_FIRED; j++) {
            if (det_ids_QDC[i] == det_ids_TDC[j]) {
                this->FAT_id[dets_fired] = det_ids_QDC[i];
            //    if (det_ids_TDC[j] == 35) cout<<"I am in TDC"<<endl;

                this->FAT_E[dets_fired] = Ql[i];
                this->FAT_ratio[dets_fired] = (double) Qs_Raw[i]/Ql_Raw[i];
                this->FAT_t[dets_fired] = TDC_ns[j];
                this->FAT_t_qdc[dets_fired] = QDC_c[i];
                dets_fired++;
            }
        }
    }
    this->FAT_DET_FIRED = dets_fired;

    for (int i=0; i<TDC_FIRED; i++) {
        this->FAT_TDC_id[i]        = det_ids_TDC[i];
        this->FAT_TDC_timestamp[i] = TDC[i];


    }
       Event_Type = 3;
}
// void Raw_Event::set_AIDA_Event(){}
 ////////////////////////////////// AIDA ////////////////////////
 void Raw_Event::set_DATA_AIDA(double* Aida_E, int* Aida_feeID, int* Aida_ChID, ULong64_t* Aida_WR, int Aida_hits, bool* Aida_veto, int* Aida_side, int* Aida_strip, int* Aida_EvtID, ULong64_t* Aida_Fast, int* Aida_ADC){

    this->  AIDA_Hits = Aida_hits;
    for (int i =0; i< Aida_hits; i++){
        this->  AIDA_Energy[i] = Aida_E[i];
        this->  AIDA_FEE[i] = Aida_feeID[i];
        this->  AIDA_CHA_ID[i] = Aida_ChID[i];
        this->  AIDA_WR[i] = Aida_WR[i];

        this->  AIDA_HighE_VETO[i] = Aida_veto[i];
        this->  AIDA_SIDE[i] = Aida_side[i];
        this->  AIDA_STRIP[i] =  Aida_strip[i];
        this->  AIDA_EVT_ID[i] = Aida_EvtID[i];

        this-> AIDA_FastTime[i] = Aida_Fast[i];

        this->AIDA_ADC[i] = Aida_ADC[i];

       }
   Event_Type = 1;
 }



//--------------------------------------FINGER-----------------------------------------------------------------//
void Raw_Event::set_DATA_FINGER(int* it,double** Edge_Coarse,double** Edge_fine,UInt** ch_ed,double* Coarse_Trigger,double* Fine_Trigger,int amount_hit_tamex, int** Lead_Arr){

    this->amount_hit_tamex = amount_hit_tamex;
    //reset lead and trail hits
    for(int i = 0;i < amount_hit_tamex;i++){
        for(int j = 0;j < 32;j++){
            leading_hits_ch[i][j] = 0;
            trailing_hits_ch[i][j] = 0;
            leading_array[i][j] = 0;
        }
    }

    //loop over all 4 tamex modules
    for(int i = 0;i < amount_hit_tamex;i++){
        iterator[i] = it[i];
        trigger_coarse[i] = Coarse_Trigger[i];
        trigger_fine[i] = Fine_Trigger[i];
        //cout << "trigger_coarse[i] " << trigger_coarse[i]<< " trigger_fine[i] "<<trigger_fine[i]<< endl;
        fired_tamex[i] = (iterator[i] > 0);
        leading_hits[i] = 0;
        trailing_hits[i] = 0;
          
        for(int j = 0;j < iterator[i];j++){
            ch_ID[i][j] = ch_ed[i][j];
            leading_array[i][j] = Lead_Arr[i][j];
            //  cout<<"1 ch_ID[i][j] " <<ch_ID[i][j]<< " coarse_T_edge_lead[i][j] " <<  coarse_T_edge_lead[i][j]<< " leading_array[i][j] " << leading_array[i][j] << endl;
            if(ch_ID[i][j] % 2 == 0){
               //cout<<"2 ch_ID[i][j] " <<ch_ID[i][j]<< " coarse_T_edge_lead[i][j] " <<  coarse_T_edge_lead[i][j]<< " leading_array[i][j] " << leading_array[i][j]<< " i " << i << " j " << j << endl;
       
                coarse_T_edge_lead[i][j] = (double) Edge_Coarse[i][j];
                fine_T_edge_lead[i][j] = (double) Edge_fine[i][j];
                phys_channel[i][j] = (ch_ID[i][j]+1)/2;
                leading_hits[i]++;
                leading_hits_ch[i][phys_channel[i][j]]++;
           //     cout << "2)phys_channel[i][j] " << phys_channel[i][j] << " i " << i << " j " << j<< " ch_ID[i][j] " << ch_ID[i][j]<< endl;
                // cout << "2) eve" << " i " << i << " j " << j <<" ChID "<< ch_ID[i][j]<<" Phys chan " << phys_channel[i][j]<<" coarse_T_edge_lead[i][j] " << coarse_T_edge_lead[i][j] << " fine_T_edge_lead[i][j] " << fine_T_edge_lead[i][j] << endl;
         
            }
            else{
                coarse_T_edge_trail[i][j] = (double)  Edge_Coarse[i][j];
                fine_T_edge_trail[i][j] =(double)  Edge_fine[i][j];
               // cout <<"RAW trail fine " << Edge_fine[i][j] << " i " << i << " j " << j <<endl;
                trailing_hits[i]++;
                phys_channel[i][j] = (ch_ID[i][j])/2;
                trailing_hits_ch[i][phys_channel[i][j]]++;

            }
        }
    }

    Event_Type = 2;
}
//---------------PLASTIC TAMEX IMPLEMENT LATER ------------------------------------------------
// void Raw_Event::set_DATA_PLASTIC(int* it,double** Edge_Coarse,double** Edge_fine,UInt** ch_ed,double* Coarse_Trigger,double* Fine_Trigger,int amount_hit_tamex){
//
//     this->amount_hit_tamex = amount_hit_tamex;
//     //reset lead and trail hits
//     for(int i = 0;i < amount_hit_tamex;++i){
//         for(int j = 0;j < 17;j++){
//             leading_hits_ch[i][j] = 0;
//             trailing_hits_ch[i][j] = 0;
//         }
//     }
//
//     //loop over all 4 tamex modules
//     for(int i = 0;i < amount_hit_tamex;i++){
//         iterator[i] = it[i];
//         trigger_coarse[i] = Coarse_Trigger[i];
//         trigger_fine[i] = Fine_Trigger[i];
//         fired_tamex[i] = (iterator[i] > 0);
//         leading_hits[i] = 0;
//         trailing_hits[i] = 0;
//
//         for(int j = 0;j < iterator[i];++j){
//             ch_ID[i][j] = ch_ed[i][j];
//
//             if(ch_ID[i][j] % 2 == 1){
//
//                 coarse_T_edge_lead[i][j] = (double) Edge_Coarse[i][j]*5.0;
//                 fine_T_edge_lead[i][j] = (double) Edge_fine[i][j]*5.0;
//
//                 phys_channel[i][j] = (ch_ID[i][j]+1)/2;
//                 leading_hits[i]++;
//                 leading_hits_ch[i][phys_channel[i][j]]++;
//                  }
//             else{
//
//                 coarse_T_edge_trail[i][j] = (double)  Edge_Coarse[i][j]*5.0;
//                 fine_T_edge_trail[i][j] =(double)  Edge_fine[i][j]*5.0;
//
//                 trailing_hits[i]++;
//                 phys_channel[i][j] = (ch_ID[i][j])/2;
//                 trailing_hits_ch[i][phys_channel[i][j]]++;
//           }
//         }
//     }
//
//     Event_Type = 2;
// }


void Raw_Event::set_DATA_PLASTIC_VME(int QDC_ite, int TDC_ite, double* VME_QDC_Dat1, double* VME_QDC_Dat2, int* VME_QDC_Chan, double* VME_TDC_Dat, int* VME_TDC_Chan){
     QDC_IT=0;
    this-> QDC_IT = QDC_ite;
    this->TDC_IT = TDC_ite;
    for (int i = 0; i<QDC_IT; i++){
    this->VME_QDC_DAT1[i] = VME_QDC_Dat1[i];
    this->VME_QDC_DAT2[i] = VME_QDC_Dat2[i];
    this->VME_QDC_CHA[i] = VME_QDC_Chan[i];
      }

    for(int j=0; j<TDC_IT; j++) {

    this->VME_TDC_CHA[j] = VME_TDC_Chan[j];
    this->VME_TDC_DAT[VME_TDC_Chan[j]] = VME_TDC_Dat[VME_TDC_Chan[j]];
    // cout << "TDC_ite " <<j<< " VME_TDC_Dat[VME_TDC_Chan[j]] " << VME_TDC_Dat[VME_TDC_Chan[j]]<< " VME_TDC_Chan[j] " <<VME_TDC_Chan[j] <<endl;

         }
       Event_Type=2;
        }

    void Raw_Event::set_DATA_SCALAR(int Scalar_ite,  double* Scalar_Dat, int* Scalar_Chan){

    this->SCALAR_ITERATOR = Scalar_ite;

    for(int i=0; i<SCALAR_ITERATOR; i++) {
    this->SCALAR_DATA[i] = Scalar_Dat[i];
    this->SCALAR_CHAN[i] = Scalar_Chan[i];
  //   cout << "TDC_ite " <<TDC_ite<< " VME_TDC_Dat[j] " << VME_TDC_Dat[j]<< " VME_TDC_Chan[j] " <<VME_TDC_Chan[j] <<endl;

         }

    Event_Type=2;
        }
//---------------------------------------------------------------

bool Raw_Event::PLASTIC_CheckVME(){
    return VME_Event;
}

//---------------------------------------------------------------

void Raw_Event::set_DATA_GALILEO(int GAL_FIRED,ULong64_t* sum_time,int* pileup,int* hit_pattern,ULong64_t* chan_time,double* chan_en, int* FEBEX_det_ids){
    this->GAL_FIRED = GAL_FIRED;

    for(int i = 0;i < GAL_FIRED;++i){

        GALILEO_Det_Nums[i] = FEBEX_det_ids[i];
        GALILEO_sum_time[i] = sum_time[i];
        GALILEO_pileup[i] = pileup[i];
        GALILEO_hit_pattern[i] = hit_pattern[i];
        GALILEO_chan_time[i] = chan_time[i];
        GALILEO_chan_energy[i] = chan_en[i];
      // cout << "2) i " <<i<<" GALILEO_chan_energy[i] "<< GALILEO_chan_energy[i]<< endl;

    }

    Event_Type = 4;
}
//TEMPORARY GETTERS FOR FRS, FATIMA, PLASTIC, and GALILEO

// #############################################################

// FRS

Float_t Raw_Event::get_FRS_dE(int i){return dE[i];}
Float_t Raw_Event::get_FRS_dE_corr(int i){return dE_cor[i];}

Float_t Raw_Event::get_FRS_sci_l(int i){return sci_l[i];}
Float_t Raw_Event::get_FRS_sci_r(int i){return sci_r[i];}
Float_t Raw_Event::get_FRS_sci_e(int i){return sci_e[i];}
Float_t Raw_Event::get_FRS_sci_tx(int i){return sci_tx[i];}
Float_t Raw_Event::get_FRS_sci_x(int i){return sci_x[i];}


Int_t Raw_Event::get_FRS_dt_21l_21r(){return dt_21l_21r;}
Int_t Raw_Event::get_FRS_dt_41l_41r(){return dt_41l_41r;}
Int_t Raw_Event::get_FRS_dt_21l_41l(){return dt_21l_41l;}
Int_t Raw_Event::get_FRS_dt_21r_41r(){return dt_21r_41r;}
Int_t Raw_Event::get_FRS_dt_42l_42r(){return dt_42l_42r;}
Int_t Raw_Event::get_FRS_dt_43l_43r(){return dt_43l_43r;}
Int_t Raw_Event::get_FRS_dt_21l_42l(){return dt_21l_42l;}
Int_t Raw_Event::get_FRS_dt_21r_42r(){return dt_21r_42r;}
Int_t Raw_Event::get_FRS_dt_81l_81r(){return dt_81l_81r;}
Int_t Raw_Event::get_FRS_dt_21l_81l(){return dt_21l_81l;}
Int_t Raw_Event::get_FRS_dt_21r_81r(){return dt_21r_81r;}

Float_t Raw_Event::get_FRS_tofll2(){return sci_tofll2;}
Float_t Raw_Event::get_FRS_tofll3(){return sci_tofll3;}
Float_t Raw_Event::get_FRS_tof2(){return sci_tof2;}
Float_t Raw_Event::get_FRS_tofrr2(){return sci_tofrr2;}
Float_t Raw_Event::get_FRS_tofrr3(){return sci_tofrr3;}
Float_t Raw_Event::get_FRS_tof3(){return sci_tof3;}

Float_t Raw_Event::get_FRS_x2(){return ID_x2;}
Float_t Raw_Event::get_FRS_y2(){return ID_y2;}
Float_t Raw_Event::get_FRS_a2(){return ID_a2;}
Float_t Raw_Event::get_FRS_b2(){return ID_b2;}

Float_t Raw_Event::get_FRS_x4(){return ID_x4;}
Float_t Raw_Event::get_FRS_y4(){return ID_y4;}
Float_t Raw_Event::get_FRS_a4(){return ID_a4;}
Float_t Raw_Event::get_FRS_b4(){return ID_b4;}

Float_t Raw_Event::get_FRS_brho(int i){return ID_brho[i];}
Float_t Raw_Event::get_FRS_rho(int i){return ID_rho[i];}

Float_t Raw_Event::get_FRS_beta(){return beta;}
Float_t Raw_Event::get_FRS_beta3(){return beta3;}
Float_t Raw_Event::get_FRS_gamma(){return gamma;}

Float_t Raw_Event::get_FRS_AoQ(){return AoQ;}
Float_t Raw_Event::get_FRS_AoQ_corr(){return AoQ_corr;}
Float_t Raw_Event::get_FRS_z(){return z;}
Float_t Raw_Event::get_FRS_z2(){return z2;}
Float_t Raw_Event::get_FRS_z3(){return z3;}

Float_t Raw_Event::get_FRS_timestamp(){return timestamp;}
Float_t Raw_Event::get_FRS_ts(){return ts;}
Float_t Raw_Event::get_FRS_ts2(){return ts2;}

// #######################################################

//White Rabbit

//---------------------------------------------------------------

void Raw_Event::set_WR(ULong64_t WR){this->WR = WR;}

ULong64_t Raw_Event::get_WR(){return WR;}
//---------------------------------------------------------------



     double Raw_Event::get_AIDA_Energy(int i){return AIDA_Energy[i];}
     int Raw_Event::get_AIDA_FEE_ID(int i){return AIDA_FEE[i];}
     int Raw_Event::get_AIDA_CHA_ID(int i){return AIDA_CHA_ID[i];}
     int Raw_Event::get_AIDA_HITS(){return AIDA_Hits;}
     ULong64_t Raw_Event::get_AIDA_WR(int i){return AIDA_WR[i];}
     bool Raw_Event::get_AIDA_HighE_VETO(int i){return AIDA_HighE_VETO[i];}
     int Raw_Event::get_AIDA_SIDE(int i){return AIDA_SIDE[i];}
     int Raw_Event::get_AIDA_STRIP(int i){return AIDA_STRIP[i];}
     int Raw_Event::get_AIDA_EVTID(int i){return AIDA_EVT_ID[i];}
     ULong64_t Raw_Event::get_AIDA_FastTime(int i){return AIDA_FastTime[i];}
     int Raw_Event::get_AIDA_ADC(int i){return AIDA_ADC[i];}



//FATIMA

//---------------------------------------------------------------

      int Raw_Event::get_FAT_det_fired(){return FAT_DET_FIRED;}
      int Raw_Event::get_FAT_id(int i){return FAT_id[i];}
   double Raw_Event::get_FAT_E(int i){return FAT_E[i];}
   double Raw_Event::get_FAT_ratio(int i){return FAT_ratio[i];}
   double Raw_Event::get_FAT_t(int i){return FAT_t[i];}
   double Raw_Event::get_FAT_t_qdc(int i){return FAT_t_qdc[i];}

      int Raw_Event::get_FAT_QDCs_fired(){return FAT_QDCs_FIRED;}
      int Raw_Event::get_FAT_QDC_id(int i){return FAT_QDC_id[i];}
   double Raw_Event::get_FAT_QLong(int i){return FAT_QLong[i];}
   double Raw_Event::get_FAT_QShort_Raw(int i){return FAT_QShort_Raw[i];}
   double Raw_Event::get_FAT_QLong_Raw(int i){return FAT_QLong_Raw[i];}
   ULong64_t Raw_Event::get_FAT_QDC_t_Coarse(int i){return FAT_QDC_t_coarse[i];}
   double Raw_Event::get_FAT_QDC_t_Fine(int i){return FAT_QDC_t_fine[i];}

      int Raw_Event::get_FAT_TDCs_fired(){return FAT_TDCs_FIRED;}
      int Raw_Event::get_FAT_TDC_id(int i){return FAT_TDC_id[i];}
   double Raw_Event::get_FAT_TDC_timestamp(int i){return FAT_TDC_timestamp[i];}


//FATIMA_DataStruct* Raw_Event::PassFATIMA(){return &FATIMA_Data;}


//---------------------------------------------------------------
//FINGER

//---------------------------------------------------------------

    int Raw_Event::get_FINGER_tamex_hits(){return amount_hit_tamex;}

    int Raw_Event::get_FINGER_am_Fired(int i){return iterator[i];}

    double Raw_Event::get_FINGER_trigger_T(int i){return (trigger_coarse[i] - trigger_fine[i])*5000;}

    int Raw_Event::get_FINGER_CH_ID(int i,int j){return ch_ID[i][j];}

    double Raw_Event::get_FINGER_lead_T(int i,int j){
        //cout << "SEND l" << coarse_T_edge_lead[i][j] << " " << fine_T_edge_lead[i][j]  << " " <<  coarse_T_edge_lead[i][j]*5 - fine_T_edge_lead[i][j] << endl;
        return (coarse_T_edge_lead[i][j] - fine_T_edge_lead[i][j])*5000;
                    }

    double Raw_Event::get_FINGER_coarse_lead(int i,int j){

        return coarse_T_edge_lead[i][j]*5000;
   }
                    
   double Raw_Event::get_FINGER_fine_lead(int i,int j){
        return fine_T_edge_lead[i][j]*5000;
   }
   double Raw_Event::get_FINGER_coarse_trail(int i,int j){
        return coarse_T_edge_trail[i][j]*5000;
   }
   double Raw_Event::get_FINGER_fine_trail(int i,int j){
        return fine_T_edge_trail[i][j]*5000;
   }
   
    double Raw_Event::get_FINGER_trail_T(int i,int j){
        //cout << "SEND t" << coarse_T_edge_trail[i][j] << " " << fine_T_edge_trail[i][j] << endl;
        return (coarse_T_edge_trail[i][j] - fine_T_edge_trail[i][j])*5000; 
                    }

    double Raw_Event::get_FINGER_TOT(int i,int j){
        // i is board ID, j is physical channel
        double T_lead = (coarse_T_edge_lead[i][j] - fine_T_edge_lead[i][j])*5000;
        double T_trail = (coarse_T_edge_trail[i][j+1] - fine_T_edge_trail[i][j+1])*5000;
               return T_trail - T_lead;
                    }

    int Raw_Event::get_FINGER_trail_hits(int i){return trailing_hits[i];}

    int Raw_Event::get_FINGER_lead_hits(int i){return leading_hits[i];}

    int Raw_Event::get_FINGER_physical_channel(int i,int j){return phys_channel[i][j];}

    int Raw_Event::get_FINGER_physical_lead_hits(int i,int j){return leading_hits_ch[i][j];}

    int Raw_Event::get_FINGER_physical_trail_hits(int i,int j){return trailing_hits_ch[i][j];}


//---------------------------------------------------------------



// PLASTIC VME

 int Raw_Event::get_plastic_VME_QDC_fired(){return QDC_IT;}
 int Raw_Event::get_plastic_VME_TDC_fired(){return TDC_IT;}
 double Raw_Event::get_plastic_VME_QDC_dat1(int i){return VME_QDC_DAT1[i];}
 double Raw_Event::get_plastic_VME_QDC_dat2(int i){return VME_QDC_DAT2[i];}
 int Raw_Event::get_plastic_VME_QDC_cha(int i){return VME_QDC_CHA[i];}
 double Raw_Event::get_plastic_VME_TDC_dat(int i){return VME_TDC_DAT[i];}
 int Raw_Event::get_plastic_VME_TDC_cha(int i){return VME_TDC_CHA[i];}

 int Raw_Event::get_scalar_iterator(){return SCALAR_ITERATOR;}
 double Raw_Event::get_scalar_data(int i){return SCALAR_DATA[i];}
 int Raw_Event::get_scalar_chan(int i){return SCALAR_CHAN[i];}
/*
PLASTIC_VME_DataStruct* Raw_Event::PassPLASTIC_VME(){ return &PLASTIC_VME_Data;}

//---------------------------------------------------------------

PLASTIC_DataStruct* Raw_Event::PassPLASTIC(){ return &PLASTIC_Data;}*/

//---------------------------------------------------------------

//GALILEO


int Raw_Event::get_GALILEO_am_Fired(){return GAL_FIRED;}

ULong64_t Raw_Event::get_GALILEO_Sum_T(int i){return GALILEO_sum_time[i];}

int Raw_Event::get_GALILEO_Pileup(int i){return GALILEO_pileup[i];}

int Raw_Event::get_GALILEO_Hit_Pattern(int i){return GALILEO_hit_pattern[i];}

ULong64_t Raw_Event::get_GALILEO_Chan_T(int i){return GALILEO_chan_time[i];}

double Raw_Event::get_GALILEO_Chan_E(int i){return GALILEO_chan_energy[i];}

int Raw_Event::get_GALILEO_Det_ids(int i){return GALILEO_Det_Nums[i];}
