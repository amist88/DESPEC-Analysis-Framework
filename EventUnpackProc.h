// $Id: TSCNUnpackProc.h 755 2011-05-20 08:04:11Z linev $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum für Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef EVENTUNPACKPROC_H
#define EVENTUNPACKPROC_H

#include "Riostream.h"

// Root Includes //
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TArc.h"
#include "TTree.h"


// Go4 Includes //
#include "TGo4UserException.h"
#include "TGo4Picture.h"
#include "TGo4MbsEvent.h"


#include "EventUnpackStore.h"
#include "CalibParameter.h"
#include "CorrelParameter.h"

#include "Detector_System.cxx"
#include "AIDA_Detector_System.h"
#include "FATIMA_Detector_System.h"
#include "PLASTIC_Detector_System.h"
#include "PLASTIC_VME_Detector_System.h"
#include "GALILEO_Detector_System_TEST.h"
#include "DESPECAnalysis.h"
//#include "FRS_Detector_System.h"

#include "Data_Stream.cxx"
#include "White_Rabbit.h"
#include "AIDA_Event.h"
#include "AIDA_Decay_Event_Store.h"
#include "PLASTIC_Data_Stream.h"

#include "Raw_Event.h"

#include "EventBuilder.cxx"
#include "Time_EventBuilder.h"
#include "AidaMover.h"
#include <string>

#include <array>
#include <vector>


/////////////////////////////////////

#include "AIDA_Headers.h"
#include "AIDA_Event.h"
#include "AIDA_Data_Types.h"
#include "TGo4MbsEvent.h"
#include "Data_Stream.cxx"
#include "EventBuilder.cxx"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"
#include "AIDA_Decay_Event_Store.h"
//#include "AIDA_Processor.h"

#include "Detector_System.cxx"

///////////////////////////////

using namespace std;

#include "TGo4EventProcessor.h"

	class EventUnpackStore;

	//#include "EventUnpackStore.h"

	class EventUnpackProc : public TGo4EventProcessor
	{
		public:

            EventUnpackProc();
			EventUnpackProc(const char* name);
            virtual ~EventUnpackProc();
			//FATIMA variables
			double FATgate1_low, FATgate1_high;
			double FATgate2_low, FATgate2_high;

            Bool_t BuildEvent(TGo4EventElement* dest);
             void Process_AIDA_Event(EventUnpackStore* event);

             ofstream WR_out;

//            std::vector<AidaCluster> EventsToClusters(std::vector<AidaEvent> const&);
//             AidaHit ClusterPairToHit(std::pair<AidaCluster, AidaCluster> const&);

            //parameters
            CalibParameter *fCal;
            CorrelParameter *fCorrel;
            EventUnpackStore *fOutput;
            AidaUnpackData    fAida;

		protected:
//                 TGo4MbsEvent *fInput;


			TH1* hsci_tofll2;

			TH1* hsci_tofll3;
			TH1* hsci_tof2;
			TH1* hsci_tofrr2;
			TH1* hsci_tofrr3;
			TH1* hsci_tof3;

			TH1* hID_x2;
			TH1* hID_y2;
			TH1* hID_a2;
			TH1* hID_b2;

			TH1* hID_x4;
			TH1* hID_y4;
			TH1* hID_a4;
			TH1* hID_b4;

			TH1* hsci_dt_21l_21r;
			TH1* hsci_dt_41l_41r;
			TH1* hsci_dt_42l_42r;
			TH1* hsci_dt_43l_43r;
			TH1* hsci_dt_81l_81r;

			TH1* hsci_dt_21l_41l;
			TH1* hsci_dt_21r_41r;

			TH1* hsci_dt_21l_42l;
			TH1* hsci_dt_21r_42r;

			TH1* hsci_dt_21l_81l;
			TH1* hsci_dt_21r_81r;


			TH1* hbeta;
			TH1* hbeta3;
			TH1* hgamma;

			TH1* hAoQ;
			TH1* hAoQ_corr;

			TH1* hz;
			TH1* hz2;
			TH1* hz3;

			TH1* htimestamp;
			TH1* hts;
			TH1* hts2;

            //AIDA Histograms
//             TH1* hAIDA_EnergyImpSum;
//             TH1* hAIDA_EnergyDecSum;
//             TH1* hAIDA_EnergyDecHitsX;
//             TH1* hAIDA_EnergyImpHitsX;
//             TH1* hAIDA_EnergyDecHitsY;
//             TH1* hAIDA_EnergyImpHitsY;

			//Fatima histograms
			//-general
			TH1* hFAT_Esum;
			TH2* hFAT_gg;
            TH1* hFAT_TDCdt_ref_sum;
			TH1* hFAT_TDCdt_refcalib_sum;
            TH1* hFAT_TDCdt_ref0_sum;
			TH1* hFAT_QDCdtsum;
			TH1* hFAT_TDCdtsum_ref_gated;
            TH1* hFAT_TDCdtsum_ref0_gated;
			TH1* hFAT_QDCdtsum_ref_gated;   //for now...
			TH1* hFAT_Angular_Diff_ref_gated; // Histogram of Gated Angular Differences
			//-statistics
			TH1* hFAT_hits;		     //number of hits per detector id
			TH1* hFAT_hits_QDC;
			TH1* hFAT_hits_TDC;
            TH1* hFAT_TDC_multi;

			TH2* hFAT_QDC_TDC_hitmap; //hits of qdc and tdc in same event
			TH2* hFAT_correlations;   //det-det coincidence map
			//-energy
			TH1* hFAT_E[48];

            TH1*  hFAT_Traw[48];
			TH1*  hFAT_Eraw[48];

            TH1* hFAT_TDC_multich[48];

			TH2** hFAT_E_ratio;
			TH2** hFAT_gg_ref;
			//-timing
			TH1* hFAT_TDCdt_ref[48];
            TH1* hFAT_TDCdt_refcalib[48];
            TH1* hFAT_TDCdt_ref0[48];
			TH1** hFAT_QDCdt_ref;
			TH2** hFAT_TDC_QDC_dt;
			TH1*  hFAT_TDCdt_ref_gated[48];
            TH1*  hFAT_TDCdt_ref0_gated[48];
			TH2** hFAT_E_TDCdt_ref_gated;

			//Other histograms
			TH1* WR_HIST;
			TH1* WR_HIST2;
			TH1* C_t;
			TH1** tamex_Mult_lead;
			TH1** tamex_Mult_trail;

			TH1*** mat;
			TH1* all;
			TH1* all2;

			TH1* WR_F;

            TH1 *hPLAS_QDCRaw1[32];
            TH1 *hPLAS_QDCRaw2[32];
            TH1 *hPLAS_TDCRaw[32];
            TH1 *hPLAS_TimeDiff[32]; // Time difference Histogram from TDC channels
            TH1 *hPLAS_TimeDiffSiPM[32];
            TH1 *hPLAS_TimeDiffCalib[32]; // Time difference Histogram from TDC channels Calibrated
            TH2 *hPLAS_CoincQ1A1; //Coincidences with germanium
            TH2 *hPLAS_CoincQ1Q2Add; //Coincidences Addback
            TH2 *hPLAS_CoincQ1Q2[32];
            TH2 *hPLAS_EETimeGated[32];
            TH1 *hPLAS_TimeD1;
            TH1 *hPLAS_TimeD2;
            TH1 *hPLAS_QDC1_hits;
            TH1 *hPLAS_TDC_hits;
            TH1 *hPLAS_TDC_multich[32];
            TH1 *hPLAS_QDC1_multi;
            TH1 *hPLAS_TDC_multi;
            TH2 *hPLAS_CoincQ1T1[32]; //Coincidences with germaniun
            TH2 *hPLAS_TempCoincQ1T1[32]; //Coincidences with germaniun
            TH1 *hPLAS_QDCGain[32];
            TH1 *hPLAS_TDCGain[32];
            TH1 *hPLAS_TDC_FATgated[32];
            TH1 *hPLAS_QDCRaw1_FATgated[32];
            TH1 *hPLAS_TimeDiff_FATgated[32];
            TH1 *hPLAS_QDCCalib1_FATgated[32];

            TH1 *hScalar_hit_pattern;

            TH2 *hCoinc_FAT_PLAS_TDetc;
            TH2 *hCoinc_FAT_PLAS_ESum;
            TH2 *hCoinc_FAT_PLAS_T;
            TH2 *hCoinc_FAT_PLAS_E[32];
            TH1 *hCoinc_FAT_PLAS_TDiff[32];
            TH1 *hCoinc_FAT_PLAS_TDiff_gg[32];
          //  TH1 *hCoinc_FATID_PLAS_TDiff[34];

            TH1 *hCoinc_FAT_PLAS_TDiffSum;
            TH1 *hCoinc_FAT_PLAS_TDiffSum_gg;

            TH1 *hPLAS_TDCSumCalib;
            //TH1 *hPLAS_QDCSum;
            TH1 *hPLAS_GatedHist;

			//Plastic histograms
			//TH1*** tamex_Mult_Ch_lead;
			//TH1*** tamex_Mult_Ch_trail;
			//TH2** tamex_mult_mat_lead;
			//TH2** tamex_mult_mat_trail;
			//TH1*** Trail_LEAD;
			TH1**** TRAIL_TRAIL;
			TH1**** LEAD_LEAD;

			//TH1*** Coarse;
			//TH1** DIFF_ARR;
			TH1**** TOT_TOT;
			TH1*** TOT_Single;
			//TH1*** LEAD_LEAD_Total;

			TH1 *hFin_ToT[33];
			TH1 *hlead_lead[48];

			TH1* FAT_TDC_Diff; // ****NEWLY ADDED****

			//for the SIS modules

			// GALILEO Histograms //
			TH1 *hGAL_Raw_E[32];

		private:


//             long adcLastTimestamp[12][4];
//             int adcCounts[12][4];
//             void ResetMultiplexer();
//             void CorrectTimeForMultiplexer(AidaEvent& evt);




             int AIDA_Hits=0;
                double AIDA_Energy[10000] = {0};
                int AIDA_FEE[10000] = {0};
                int AIDA_ChID[10000] = {0};
                ULong64_t AIDA_Time[10000] = {0};
                bool AIDA_HighE_veto[10000] = {false};
                int AIDA_Side[10000] = {0};
                int AIDA_Strip[10000] = {0};
                int AIDA_evtID[10000] = {0};

			const int FATIMA_reference_det = 0;
			const int FAT_MAX_DET = 48;

			const int FRS = 0;
			const int AIDA = 1;
			const int PLASTIC = 2;
			const int FATIMA = 3;
			const int GALILEO = 4;
			const int FINGER = 5;

			int FAT_REF_DET;
			

			float E_gate1,E_gate2;

			Bool_t ffill;
			Int_t fshift;
			ULong64_t White_Rabbbit_old;

			Int_t PrcID_Array[10][5];
			bool Used_Systems[10];

            int tamID[4][100];
            int tamCH[4][100];
            int fingID[4][16];

			bool SKIP_EVT_BUILDING;

			bool PADI_OR_PADIWA,VME_TAMEX;

        ///For AIDA

            long lastTime;
            int ID;
            AidaEvent old;
            AidaEvent evt;

            long startTime;
            long stopTime;

      /* Multiplexer correction */
						std::vector<std::array<uint64_t, 4>> adcLastTimestamp;
						std::vector<std::array<int, 4>> adcCounts;
            void ResetMultiplexer();
            void CorrectTimeForMultiplexer(AidaEvent& evt);

            int totalEvents;
            int implantEvents;
            int decayEvents;
            int pulserEvents;
            int nonsenseEvents;

						// AIDA histograms
						std::vector<std::array<std::array<TH1*, 2>, 64>>	hAIDA_ADC;
						TH2* hAIDA_ADC_unaligned;
						TH2* hAIDA_ADC_aligned;

      ///End AIDA
			double vals[100000];
			int val_it;
            int event_number;
			string input_data_path;
			string input_data_path_old;
            string test;

			bool cals_done,WR_used;
			bool WHITE_RABBIT_USED; // Read from General Setup File
			bool FAT_make_raw_histograms;
			
			
			int file_pwd, file_end;
			std::string gain_match_filename;
			int data_file_number = 0;

			Detector_System** Detector_Systems;
			Data_Stream** data_stream;
			White_Rabbit* WR;
			Raw_Event* RAW;
           // AIDA_Event* Aida_inp;
			EventBuilder** EvtBuilder;
            AidaEvent CallTheAida();

			int amount_interest;
			int* length_interest;
			int** interest_array;
            int FAT_GamGatelow;
            int FAT_GamGatehigh;

			//Event_Builder** EvtBuilder;

			double fatima_E_save[4];
			int am_FATIMA_hits;
			int num_full_FAT_evts;
			

			Int_t get_Conversion(Int_t);
			void get_used_Systems();
			void get_WR_Config();

			void read_setup_parameters();

            void load_FingerID_File();

			void load_PrcID_File();
			void get_interest_arrays();

			void Make_FRS_Histos();
			void Fill_FRS_Histos();

            void Make_AIDA_Histos();
            void Fill_AIDA_Histos();

			void Make_Plastic_Histos();
			void Fill_Plastic_Histos();

            void Make_Plastic_VME_Histos();
            void Fill_Plastic_VME_Histos();

			void Make_FATIMA_Histos();
			void Fill_FATIMA_Histos();

			void Make_GALILEO_Histos();
			void Fill_GALILEO_Histos();

            void Make_Finger_Histos();
            void Fill_Finger_Histos();

			void FILL_HISTOGRAMS(int);


			bool Check_Cal_Plastic();

			void checkPADI_or_PADIWA();
			void checkTAMEXorVME();

			bool PLASTIC_CALIBRATION;


			
// 			void FAT_det_pos_setup();
// 			double distance_between_detectors(double, double, double, double, double, double);
// 			double angle_between_detectors(double, double, double);

			void print_MBS(int*,int);

			int count;
            int array_count;
			int called[2];
			int iterator;
			int Cout_counter;

			ULong64_t WR_tmp;
            ULong64_t WR_main;
            int WR_count;
            int WR_d;
            ULong64_t WR_AIDA[10000];
            Long64_t  WR_diff[10000];
            int AIDA_Loop;
			ClassDef(EventUnpackProc,1)
	};

#endif 
