#ifndef RAW_EVENT_H
#define RAW_EVENT_H

#include <math.h>
#include <float.h>
#include <vector>


//#include "PLASTIC_DataStruct.h"
//#include "PLASTIC_VME_DataStruct.h"
//#include "FATIMA_DataStruct.h"
#include "AIDA_Event.h"
#include "AIDA_Decay_Event_Store.h"

//#include "Rtypes.h"

typedef unsigned long ULong;
typedef unsigned int UInt;
typedef unsigned long long ULong64_t;
typedef float Float_t;
typedef int Int_t;


class Raw_Event{

private:


	// ##########################################################

	//FRS
	// MUSIC PARAMETERS //

	Float_t dE[3];      // set_DATA_MUSIC
	Float_t dE_cor[3];  // set_DATA_MUSIC

	// SCINTILLATOR PARAMETERS //

	Float_t sci_l[12];  // set_DATA_SCI
	Float_t sci_r[12];  // set_DATA_SCI
	Float_t sci_e[12];  // set_DATA_SCI
	Float_t sci_tx[12]; // set_DATA_SCI
	Float_t sci_x[12];  // set_DATA_SCI



	Int_t dt_21l_21r;
	Int_t dt_41l_41r;
	Int_t dt_21l_41l;
	Int_t dt_21r_41r;
	Int_t dt_42l_42r;
	Int_t dt_43l_43r;
	Int_t dt_21l_42l;
	Int_t dt_21r_42r;
	Int_t dt_81l_81r;
	Int_t dt_21l_81l;
	Int_t dt_21r_81r;


	Float_t sci_tofll2; // set_DATA_SCI_ToF
	Float_t sci_tofll3; // set_DATA_SCI_ToF
	Float_t sci_tof2;   // set_DATA_SCI_ToF
	Float_t sci_tofrr2; // set_DATA_SCI_ToF
	Float_t sci_tofrr3; // set_DATA_SCI_ToF
	Float_t sci_tof3;   // set_DATA_SCI_ToF

	// ID PARAMETERS //

	Float_t ID_x2;      // set_DATA_ID_2_4
	Float_t ID_y2;      // set_DATA_ID_2_4
	Float_t ID_a2;      // set_DATA_ID_2_4
	Float_t ID_b2;      // set_DATA_ID_2_4

	Float_t ID_x4;      // set_DATA_ID_2_4
	Float_t ID_y4;      // set_DATA_ID_2_4
	Float_t ID_a4;      // set_DATA_ID_2_4
	Float_t ID_b4;      // set_DATA_ID_2_4

	Float_t ID_brho[2]; // set_DATA_ID_Beta_Rho
	Float_t ID_rho[2];  // set_DATA_ID_Beta_Rho

	Float_t beta;       // set_DATA_ID_Beta_Rho
	Float_t beta3;      // set_DATA_ID_Beta_Rho
	Float_t gamma;      // set_DATA_ID_Beta_Rho

	Float_t AoQ;        // set_DATA_ID_Z_AoQ
	Float_t AoQ_corr;   // set_DATA_ID_Z_AoQ

	Float_t z;          // set_DATA_ID_Z_AoQ
	Float_t z2;         // set_DATA_ID_Z_AoQ
	Float_t z3;         // set_DATA_ID_Z_AoQ

	Float_t timestamp;  // set_DATA_ID_Timestamp
	Float_t ts;         // set_DATA_ID_Timestamp
	Float_t ts2;        // set_DATA_ID_Timestamp


	// ##########################################################

	//White Rabbit
	ULong64_t WR;

    //AIDA
    double AIDA_Energy[10000];
    int AIDA_FEE[10000];
    int AIDA_CHA_ID[10000];
    ULong64_t AIDA_WR[10000];
    int AIDA_Hits;
    bool AIDA_HighE_VETO[10000];
    int AIDA_SIDE[10000];
    int AIDA_STRIP[10000];
    int AIDA_EVT_ID[10000];
	ULong64_t AIDA_FastTime[10000];
	int AIDA_ADC[10000];



//     double AIDA_DecayEnergy[12][64][2];
//     ULong64_t AIDA_DecayTimestamp[24][64];
//     int AIDA_DecayID[24];
//     int AIDA_DecayFEE;
//     int AIDA_Decayhits;
//
//     double AIDA_ImpEnergy[24][64];
//     ULong64_t AIDA_ImpTimestamp[24][64];
//     int AIDA_ImpID[24];
//     int AIDA_ImpFEE;
//     int AIDA_Imphits;

    int    FAT_DET_FIRED;         //number of completed detectors in evt
    int    FAT_id[100];           //id according to allocation file
    double FAT_E[100];            //energy (calibrated nad gain matched)
    double FAT_ratio[100];        //Qshort/Qlong
    double FAT_t[100];            //tdc time  (ns, shifted)
    double FAT_t_qdc[100];        //qdc time  (ns, shifted)
    //For trouble debugging:
    int    FAT_QDCs_FIRED;
    int    FAT_QDC_id[100];
    double FAT_QLong[100];        //calibrated
    double FAT_QLong_Raw[100];
    double FAT_QShort_Raw[100];
    ULong  FAT_QDC_t_coarse[100];
    double FAT_QDC_t_fine[100];  //qdc time (ns)
    //
    int    FAT_TDCs_FIRED;
    int    FAT_TDC_id[100];
    double FAT_TDC_timestamp[100];//tdc time raw
    // for vector array
	//FATIMA_DataStruct FATIMA_Data;


	//PLASTIC
//  	PLASTIC_DataStruct PLASTIC_Data;
//  	PLASTIC_VME_DataStruct PLASTIC_VME_Data;

//FINGER
    int amount_hit_tamex;
    int iterator[4];
    double trigger_coarse[32];
    double trigger_fine[32];
    int leading_array[4][32];
    int leading_hits[4];
    int trailing_hits[4];
    int phys_channel[4][32];
    int leading_hits_ch[4][32];
    int trailing_hits_ch[4][32];
    UInt ch_ID[4][32];

    double coarse_T_edge_lead[4][32];
    double coarse_T_edge_trail[4][32];
    double fine_T_edge_lead[4][32];
    double fine_T_edge_trail[4][32];

    bool fired_tamex[4];


	bool VME_Event;

     //plastic VME
    int QDC_IT;
    int TDC_IT;
    double VME_QDC_DAT1[32];
    double VME_QDC_DAT2[32];
    int VME_QDC_CHA[32];
    double VME_TDC_DAT[50];
    int VME_TDC_CHA[50];

    int SCALAR_ITERATOR;
    double SCALAR_DATA[50];
    int SCALAR_CHAN[50];

	//GALILEO
	int GAL_FIRED;
	int GALILEO_Det_Nums[36];
	ULong GALILEO_sum_time[36];
	int GALILEO_pileup[36];
	int GALILEO_hit_pattern[36];
	ULong GALILEO_chan_time[36];
	double GALILEO_chan_energy[36];
    int DET[36];


	bool ch51;

	int Event_Type;


public:
	//Raw_Event(bool);
    Raw_Event();
	~Raw_Event();

	// ##########################################################

	// FRS STUFF //
	void set_DATA_MUSIC(Float_t*,Float_t*);
	void set_DATA_SCI(Float_t*,Float_t*,Float_t*,Float_t*,Float_t*);
	void set_DATA_SCI_dT(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);
	void set_DATA_SCI_ToF(Float_t,Float_t,Float_t,Float_t,Float_t,Float_t);
	void set_DATA_ID_2_4(Float_t,Float_t,Float_t,Float_t,Float_t,Float_t,Float_t,Float_t);
	void set_DATA_ID_Beta_Rho(Float_t*,Float_t*,Float_t,Float_t,Float_t);
	void set_DATA_ID_Z_AoQ(Float_t,Float_t,Float_t,Float_t,Float_t);
	void set_DATA_ID_Timestamp(Float_t,Float_t,Float_t);
	// FRS STUFF //


	// ##########################################################
//      void set_DATA_AIDA_DECAY(double***, int, int*, ULong64_t**, int);
//      void set_DATA_AIDA_IMP(double**, int, int*, ULong64_t**, int);

     void set_DATA_AIDA(double*, int*, int*, ULong64_t*, int, bool*, int*, int*, int*, ULong64_t*, int*);
     //void Nset_DATA_AIDA(AidaEvent*);

     void set_AIDA_Event(int);

     void set_DATA_FATIMA(int,int,
                         double*,double*,double*,
                         ULong64_t*,double*,
                         ULong64_t*,double*,
                         int*,int*);

    void set_DATA_FINGER(int*,double**,double**,UInt**,double*,double*,int,int**);

    void set_DATA_PLASTIC(int*,double**,double**,UInt**,double*,double*,int);

    void set_DATA_PLASTIC_VME(int, int, double*, double*, int*, double*, int*);

    void set_DATA_SCALAR(int, double*, int*);

    void set_DATA_GALILEO(int,ULong64_t*,int*,int*,ULong64_t*,double*,int*);



// 	void set_DATA_PLASTIC(std::vector<int> &it,std::vector<std::vector<double> > &Edge_Coarse,
//                           std::vector<std::vector<double> > &Edge_fine, std::vector<std::vector<UInt> > &ch_ed,
//                           std::vector<double> &Coarse_Trigger,std::vector<double> &Fine_Trigger,int amount_hit_tamex);
//
// 	void set_DATA_PLASTIC_VME(int TDC_iterator,std::vector<double> &VME_QDC_Data,std::vector<double> &VME_QDC_Channels,
// 						  std::vector<double> &VME_TDC_Data,std::vector<double> &VME_TDC_Channels);

	//void set_DATA_GALILEO(int,ULong64_t*,int*,int*,ULong64_t*,double*,int*);

	int get_Event_type();

	bool PLASTIC_CheckVME();


	// ####################################################


	//temporary FRS getters
	Float_t get_FRS_dE(int);
	Float_t get_FRS_dE_corr(int);

	Float_t get_FRS_sci_l(int);
	Float_t get_FRS_sci_r(int);
	Float_t get_FRS_sci_e(int);
	Float_t get_FRS_sci_tx(int);
	Float_t get_FRS_sci_x(int);

	Int_t get_FRS_dt_21l_21r();
	Int_t get_FRS_dt_41l_41r();
	Int_t get_FRS_dt_21l_41l();
	Int_t get_FRS_dt_21r_41r();
	Int_t get_FRS_dt_42l_42r();
	Int_t get_FRS_dt_43l_43r();
	Int_t get_FRS_dt_21l_42l();
	Int_t get_FRS_dt_21r_42r();
	Int_t get_FRS_dt_81l_81r();
	Int_t get_FRS_dt_21l_81l();
	Int_t get_FRS_dt_21r_81r();

	Float_t get_FRS_tofll2();
	Float_t get_FRS_tofll3();
	Float_t get_FRS_tof2();
	Float_t get_FRS_tofrr2();
	Float_t get_FRS_tofrr3();
	Float_t get_FRS_tof3();

	Float_t get_FRS_x2();
	Float_t get_FRS_y2();
	Float_t get_FRS_a2();
	Float_t get_FRS_b2();

	Float_t get_FRS_x4();
	Float_t get_FRS_y4();
	Float_t get_FRS_a4();
	Float_t get_FRS_b4();

	Float_t get_FRS_brho(int);
	Float_t get_FRS_rho(int);

	Float_t get_FRS_beta();
	Float_t get_FRS_beta3();
	Float_t get_FRS_gamma();

	Float_t get_FRS_AoQ();
	Float_t get_FRS_AoQ_corr();
	Float_t get_FRS_z();
	Float_t get_FRS_z2();
	Float_t get_FRS_z3();

	Float_t get_FRS_timestamp();
	Float_t get_FRS_ts();
	Float_t get_FRS_ts2();

	// ####################################################
    double get_AIDA_Energy(int i);
    int get_AIDA_FEE_ID(int i);
    int get_AIDA_CHA_ID(int i);
    ULong64_t get_AIDA_WR(int i);
    int get_AIDA_HITS();
    bool get_AIDA_HighE_VETO(int i);
    int get_AIDA_SIDE(int i);
    int get_AIDA_STRIP(int i);
    int get_AIDA_EVTID(int i);
	ULong64_t get_AIDA_FastTime(int i);
	int get_AIDA_ADC(int i);


//     double get_AIDA_DecayEnergy(int i,int j,int k);
//     ULong64_t get_AIDA_DecayTimestamp(int i,int j);
//     int get_AIDA_DecayID(int i);
//     int get_AIDA_Decayhits();
//     int get_AIDA_DecayFEEID();
//
//     double get_AIDA_ImpEnergy(int i,int j);
//     ULong64_t get_AIDA_ImpTimestamp(int i,int j);
//     int get_AIDA_ImpID(int i);
//     int get_AIDA_Imphits();
//     int get_AIDA_ImpFEEID();

	//temporary FATIMA getters
	  int get_FAT_det_fired();
	  int get_FAT_id(int i);
      double get_FAT_E(int i);
      double get_FAT_ratio(int i);
      double get_FAT_t(int i);
      double get_FAT_t_qdc(int i);

      int get_FAT_QDCs_fired();
      int get_FAT_QDC_id(int i);
      double get_FAT_QLong(int i);
      double get_FAT_QShort_Raw(int i);
      double get_FAT_QLong_Raw(int i);
      ULong64_t get_FAT_QDC_t_Coarse(int i);
      double get_FAT_QDC_t_Fine(int i);

      int get_FAT_TDCs_fired();
	  int get_FAT_TDC_id(int i);
      double get_FAT_TDC_timestamp(int i);

   // FATIMA_DataStruct* PassFATIMA();

    //temporary FINGER getters
    int get_FINGER_am_Fired(int);
    double get_FINGER_trigger_T(int);
    int get_FINGER_CH_ID(int,int);
    double get_FINGER_lead_T(int,int);
    double get_FINGER_trail_T(int,int);
    int get_FINGER_trail_hits(int);
    int get_FINGER_lead_hits(int);
    int get_FINGER_physical_channel(int,int);
    int get_FINGER_physical_lead_hits(int,int);
    int get_FINGER_physical_trail_hits(int,int);
    double get_FINGER_coarse_lead(int,int);
    double get_FINGER_fine_lead(int,int);
    double get_FINGER_coarse_trail(int,int);
    double get_FINGER_fine_trail(int,int);
    double get_FINGER_TOT(int,int);
    double get_FINGER_TOT_added(int,int);
    int get_FINGER_tamex_hits();

	//temporary PLASTIC getters
	int get_PLASTIC_am_Fired(int);
	double get_PLASTIC_trigger_T(int);
	int get_PLASTIC_CH_ID(int,int);
	double get_PLASTIC_lead_T(int,int);
	double get_PLASTIC_trail_T(int,int);
	int get_PLASTIC_trail_hits(int);
	int get_PLASTIC_lead_hits(int);
	int get_PLASTIC_physical_channel(int,int);
	int get_PLASTIC_physical_lead_hits(int,int);
	int get_PLASTIC_physical_trail_hits(int,int);
	double get_PLASTIC_coarse_lead(int,int);
	double get_PLASTIC_TOT(int,int);
	int get_PLASTIC_tamex_hits();

// 	PLASTIC_VME_DataStruct* PassPLASTIC_VME();
// 	PLASTIC_DataStruct* PassPLASTIC();
        int get_plastic_VME_QDC_fired();
        int get_plastic_VME_TDC_fired();
        double get_plastic_VME_QDC_dat1(int);
        double get_plastic_VME_QDC_dat2(int);
        int get_plastic_VME_QDC_cha(int);
        double get_plastic_VME_TDC_dat(int);
        int get_plastic_VME_TDC_cha(int);

      int get_scalar_iterator();
      double get_scalar_data(int);
      int get_scalar_chan(int);



	//temporary GALILEO getters
	int get_GALILEO_am_Fired();
	ULong64_t get_GALILEO_Sum_T(int);
	int get_GALILEO_Pileup(int);
	int get_GALILEO_Hit_Pattern(int);
	ULong64_t get_GALILEO_Chan_T(int);
	double get_GALILEO_Chan_E(int);
    int get_GALILEO_Det_ids(int);


	//White Rabbit setter and getter
	void set_WR(ULong64_t);
	ULong64_t get_WR();


};


#endif
