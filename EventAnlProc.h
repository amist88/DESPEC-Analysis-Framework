// $Id: EventAnlProc.h 755 2011-05-20 08:04:11Z linev $
//Mistry 10.04.19
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

#ifndef EVENTANLPROCESSOR_H
#define EVENTANLPROCESSOR_H

#include "TGo4EventProcessor.h"

//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "CalibParameter.h"
#include "CorrelParameter.h"
#include "AIDA_Event.h"
#include "TCutG.h"
#include "TGraph.h"
#include "Go4ConditionsBase/TGo4WinCond.h"
#include "Go4ConditionsBase/TGo4PolyCond.h"

class EventAnlStore;
class TSCNParameter;

class EventAnlProc : public TGo4EventProcessor {
   public:
      EventAnlProc();
      EventAnlProc(const char * name);
      virtual ~EventAnlProc();
      CalibParameter *fCal;
      CorrelParameter *fCorrel;
      virtual Bool_t BuildEvent(TGo4EventElement* dest);
      virtual Bool_t BuildEvent2(TGo4EventElement* dest);
      TGo4PolyCond *fCond_FingToTvsStrip;

      
      void ProcessAida(EventUnpackStore* pInput);

      int PrcID_Conv[7];
      int Used_Systems[7];
      int event_number;

      void get_used_Systems();
      
      Float_t  FRS_dE[2],  FRS_dE_cor[2];
      Float_t  FRS_sci_l[12],  FRS_sci_r[12],  FRS_sci_e[12],  FRS_sci_tx[12],  FRS_sci_x[12];
      Float_t  FRS_sci_tofll2,  FRS_sci_tofll3, FRS_sci_tof2, FRS_sci_tofrr2, FRS_sci_tofrr3, FRS_sci_tof3;
      Float_t FRS_ID_x2, FRS_ID_y2, FRS_ID_a2, FRS_ID_b2;
      Float_t FRS_ID_x4, FRS_ID_y4, FRS_ID_a4, FRS_ID_b4;
      Int_t FRS_sci_dt_21l_21r, FRS_sci_dt_41l_41r, FRS_sci_dt_42l_42r, FRS_sci_dt_43l_43r;
      Int_t FRS_sci_dt_21l_41l, FRS_sci_dt_21r_41r, FRS_sci_dt_21l_42l, FRS_sci_dt_21r_42r;
      Float_t FRS_ID_brho[2], FRS_ID_rho[2];
      Float_t FRS_beta, FRS_beta3, FRS_gamma;
      Float_t FRS_AoQ, FRS_AoQ_corr;
      Float_t FRS_z, FRS_z2, FRS_z3;
      Float_t FRS_timestamp, FRS_ts, FRS_ts2;  
      
      int Aida_Fired;
      Long64_t WR_Aida_Det_diff[10000];
      
      int bPlasQDCFired;
      int bPlasQDCID[32];
      double bPlasQDC[32];
      double bPlasQDCGain[32];

      int bPlasTDCFired;
      int bPlasTDCID[50];
      double bPlasTDC_TS[50][32];
      double bPlasTDC_ref;
      int bPlas_TDC_Multiplicity[32];
      double bPlas_TDC_diff[50];
      double bPlas_TDC_diff_sum;
      int ScalarFired;
      int ScalarID;

      int FatQDCFired;
      int FatQDCID[50];
      double FatQDC[50];
      double FatQDC_T[50];

      int FatTDCFired;
      int FatTDCID[50];
      double FatTDC_TS[50][50];
      int FatTDC_Multipl[50];
      double Fat_CHA_0_TDC;
      Long64_t Fat_WR;
      double SC41, SC41_ns;

      int GalFired;
      int GalID[32];
      double GalE[32];
      double GalT[32];
      int GalPileup;
      Long64_t Gal_WR;

      int Fing_firedTamex;
      int Fing_leadHits[4];
      int Fing_trailHits[4];
      int Fing_iterator[4];
      double Fing_trig[4];
      int Fing_leadChan[4][32];
      int Fing_trailChan[4][32];
      double Fing_lead_coarse[4][32];
      double Fing_lead_fine[4][32];
      double Fing_trail_coarse[4][32];
      double Fing_trail_fine[4][32];
      
      double Fing_leadT[4][32];
      double Fing_trailT[4][32];
      double Fing_TOT[4][32];
      int Fing_chID[4][32];
      double dataSetPerEvent[50];
      double pmtSetPerEvent[50];
      double totaltimeEvent[50];
      double maxToT;
      int maxToTChan;

     void Make_FRS_Histos();
     void Make_Aida_Histos();
     void Make_Plastic_VME_Histos();
     void Make_Fatima_Histos();
     void Make_Galileo_Histos();
     void Make_Finger_Histos();

     void Make_Fat_Plas_Histos();
     void Make_Fing_Plas_Histos();

     void Do_FRS_Histos(EventAnlStore* pOutput);
     void Do_Plastic_VME_Histos(EventAnlStore* pOutput);
     void Do_Fatima_Histos(EventAnlStore* pOutput);
     void Do_Galileo_Histos(EventAnlStore* pOutput);
     void Do_Finger_Histos(EventAnlStore* pOutput);
     void Do_Fing_Plas_Histos();
     void Do_Fat_Plas_Histos(EventAnlStore* pOutput);
      // TH1 *GermaniumCal;

     void read_setup_parameters();
     void FAT_det_pos_setup();
     double distance_between_detectors(double, double, double, double, double, double);
     double angle_between_detectors(double, double, double);
     
            long lastTime;
            int ID;
            AidaEvent old;
            AidaEvent evt;

            long startTime;
            long stopTime;

      /* Multiplexer correction */


            int totalEvents;
            int implantEvents;
            int decayEvents;
            int pulserEvents;
            int nonsenseEvents;

            std::vector<TH2*> implants_strip_xy;
            std::vector<TH2*> implants_pos_xy;
            std::vector<TH1*> implants_e;
            std::vector<TH2*> implants_e_xy;
            std::vector<TH1*> implants_time_delta;
            std::vector<TH1*> implants_strip_1d;
            std::vector<TH1*> implants_per_event;
            std::vector<TH1*> implants_channels;

            std::vector<TH2*> decays_strip_xy;
            std::vector<TH2*> decays_pos_xy;
            std::vector<TH1*> decays_e;
            std::vector<TH2*> decays_e_xy;
            std::vector<TH1*> decays_time_delta;
            std::vector<TH1*> decays_strip_1d;
            std::vector<TH1*> decays_per_event;
            std::vector<TH1*> decays_channels;
//       std::vector<AidaEvent> ImplantEvents;
//       std::vector<AidaEvent> DecayEvents;
//
//       std::vector<AidaHit> Implants;
//       std::vector<AidaHit> Decays;
             std::vector<AidaCluster> EventsToClusters(std::vector<AidaEvent> const&);
            AidaHit ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const&);

            TH1 *hPLAS_QDCCalib1[32];
            TH1 *hPLAS_QDCCalib1Sum;
            TH1 *hPLAS_QDC1_hits;
            TH2 *hPLAS_CoincE1E2[32];
            TH2 *hPLAS_CoincE1E2_Sum;
            TH2 *hPLAS_E1E2TimeGated[32];
            
            TH1 *hPLAS_TDCCalib1[32];
            TH1 *hPLAS_TDCCalib1Sum;
            TH1 *hPLAS_TDC_hits;
            TH1 *hPLAS_TDC_multi;
            TH1 *hPLAS_TDC_FiredRatio[32];
            TH1 *hPLAS_TDC_FiredRatio_Sum;
            TH1 *hPLAS_TimeDiffSiPM_Ch_Raw[32];
            TH1 *hPLAS_TimeDiffSiPM_Ch_Calib[32];
            TH1 *hPLAS_TimeDiffSiPM_Ch_Calib_Egated[32];
            TH1 *hPLAS_TimeDiffSiPM_Ch_Sum;
            TH1 *hPLAS_TimeDiffSiPM_Ch_Sum_M1;
            TH1 *hPLAS_TimeDiffSiPM_Ch_Sum_M2;
            TH1 *hPLAS_TimeDiffSiPM_Ch_Sum_M3P;
            
            TH1 *hPLAS_TimeDiffSiPM_Ch_Sum_Egated;
            TH1 *hPLAS_TimeDiff_SC41_Raw[32];
            TH1 *hPLAS_TimeDiff_SC41_Calib[32];
            TH1 *hPLAS_TimeDiff_SC41_Calib_Egated[32];
            TH1 *hPLAS_TimeDiff_SC41_Sum;
            TH1 *hPLAS_TimeDiff_SC41_Sum_M1;
            TH1 *hPLAS_TimeDiff_SC41_Sum_M2;
            TH1 *hPLAS_TimeDiff_SC41_Sum_M3P;
            TH1 *hPLAS_TimeDiff_SC41_Sum_Egated;
            TH1 *hPLAS_TimeDiff_SC41_avg;
            TH1 *hPLAS_TimeDiff_SiPM_avg;
            TH2 *hPLAS_CoincE_dTSiPM_Ch[32];
            TH2 *hPLAS_CoincE_dTSiPM_Sum;
            TH1 *hPLAS_bPlasTDC_avg;
            TH1 *hPLAS_TDC_multich[32];
            
            TH1 *hScalar_hit_pattern;
            
            TH1 *hFAT_QDCCalib1[50];
            TH1 *hFAT_QDCdt[50];      
            TH1 *hFAT_QDCCalib1Sum;
            TH1 *hFAT_hits_QDC;
            TH2 *hFAT_Chan_E_Mat[50];
            TH2 *hFAT_E_Mat_Sum;
            
            TH1 *hFAT_TDCCalib1[50];
            TH1 *hFAT_TDC_Multipl_ch[50];
            TH1 *hFAT_TDC_Multipl_PerChan;
            TH1 *hFAT_TDC_Multipl;
            TH1 *hFAT_hits_TDC;
            TH1 *hFAT_TDCdt_refSC41[50];
            TH1 *hFAT_TDCdt_refSC41_M1[50];
            TH1 *hFAT_TDCdt_refSC41_M2[50];
            TH1 *hFAT_TDCdt_refSC41_M3[50];
            
            TH1 *hFAT_TDCdt_refSC41_gated[50];
            TH1 *hFAT_TDCdt_refSC41_M1_gated[50];
            TH1 *hFAT_TDCdt_refSC41_M2_gated[50];
            TH1 *hFAT_TDCdt_refSC41_M3_gated[50];
            TH1 *hFAT_TDCdt_refSC41calib[50];
            
            TH1 *hFAT_TDCdt_refSC41_Sum;
            TH1 *hFAT_TDCdt_refSC41_Sum_calib;
            TH1 *hFAT_TDCdt_refCha[50];
            TH1 *hFAT_TDCdt_refCha_M1[50];
            TH1 *hFAT_TDCdt_refCha_M2[50];
            TH1 *hFAT_TDCdt_refCha_M3[50];
            TH1 *hFAT_TDCdt_refCha_gated[50];
            TH1 *hFAT_TDCdt_refCha_Sum;
            TH1 *hFAT_TDCdt_refCha_Sum_M1;
            TH1 *hFAT_TDCdt_refCha_Sum_M2;
            TH1 *hFAT_TDCdt_refCha_Sum_M3;
            TH1 *hFAT_TDCdt_refCha_Sum_gated;
            TH1 *hFAT_TDCdt_refSC41_Sum_gated;
            TH2 *hFAT_QDC_vs_TDC_SiPMdT_Ch[50];
            TH2 *hFAT_QDC_vs_TDC_SiPMdT;
            TH2 *hFAT_QDC_vs_TDC_SC41dT_Ch[50];
            TH2 *hFAT_QDC_vs_TDC_SC41dT;
            TH1 *hFAT_test[50];
            
            TH1 *hGAL_Chan_E[32];
            //TH1 *hGAL_Chan_E2;
            TH1 *hGAL_Chan_Egate;
            TH1 *hGAL_Chan_E_M1;
            TH1 *hGAL_Chan_E_M2;
            TH1 *hGAL_AddbackSum;
            TH1 *hGAL_Chan_Time_Diff[32];
            TH1 *hGAL_Time_Diff_vs_Energy[32];
            TH1 *hGAL_ESum;
            TH1 *hGAL_Hit_Pat;
            TH1 *hGAL_Multi;
            TH1 *hGAL_Pileup;
            TH2 *hGAL_Chan_E_Mat;
            
            TH1 *hFING_Hit_Pat;
            TH1 *hFING_lead_lead[52];
            TH1 *hFING_trail_trail[52];
            TH1 *hFING_ToT[52];
            TH1 *hFING_MaxToT[52];
            TH1 *hFING_trig_lead[52];
            TH2 *hFING_ToT_StripID;
            TH2 *hFING_MaxToT_StripID;
            TH2 *hFING_MaxToTExp_StripID;
            TH2 *hFING_ToT_lead_lead[52];
            TH2 *hFING_Pos;                 //E.Sahin 24.05
            TH2 *hFING_ToT_StripID_Exp;
            TH1 *hFING_fcoarse[52]; 
            TH1 *hFING_ffine[52];
          
            
            TH1 *hFat_minus_plasticTDC[50];
            TH1 *hFat_minus_plasticTDCFatGate[50];
            TH1 *hFat_bplas_Corr_RawTDCFAT_RawTDCbP_Dets[50];
            
            TH1 *hFat_bplas_Corr_SC41_Dets[50];
            TH1 *hFat_bplas_Corr_SC41_Dets_gated[50];
            TH1 *hFat_bplas_Corr_SiPM_Dets[50];
            TH2 *hFat_bplas_Corr_EFatEbPlas;
            TH1 *hFat_bplas_Corr_SiPM_Gated_Dets[50];
            
            TH2 *h_FingStrip_PlasID;
            TH2 *h_FingToT_PlasE; 
            TH2 *h_FingToT_PlasT[32];
            
            TSCNParameter *fParam1;
            TGo4Picture *picPMT;
private :
            bool FAT_dist_corr_used; // Read from General Setup File
            int FAT_exclusion_dist; // Read from General Setup File
            int FAT_num_TDC_modules; // Read from General Setup File
            bool FAT_nearest_neighbour_exclusion; // Read from General Setup File
            bool same_ring_exclusion; // Read from General Setup File
            bool output_position_matrix; // Read from General Setup File
            
            double** FAT_positions;
            double** FAT_angle_diffs;
            bool** FAT_neighbour_check;

      ClassDef(EventAnlProc, 1)
	};
#endif //EVENTANLPROCESSOR_H
