#ifndef AIDA_DETECTOR_SYSTEM_H
#define AIDA_DETECTOR_SYSTEM_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <stdlib.h>
#include <cstdlib>

#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TArc.h"
#include "TTree.h"
#include "TH1.h"
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



typedef unsigned long long ULong64_t;
// struct UnpackProcStats
// {
//   uint64_t Start;
//   uint64_t Stop;
//   int MBSEvents;
//   int AIDAEvents;
//   int TimeWarps;
//   struct
//   {
//     int ADC;
//     int Info;
//     int WAVE;
//     int Unknown;
//   } Type;
//   struct
//   {
//     int Pause;
//     int Resume;
//     int Sync48;
//     int Sync63;
//     int Discriminator;
//     int Unknown;
//   } InfoType;
// };

class AIDA_Detector_System : public Detector_System
{
private:

  int Hits;
  double* Energy;
  int* FEE_ID;
  int* Ch_ID;
  ULong64_t* Time;
  ULong64_t* FastTime;
  bool* HighE_veto;
  int* Side;
  int* Strip;
  int* EvtID;
  int* AdcData;

  Int_t* pdata;
  Int_t* pdata_start;
  Int_t lwords;
  Int_t sub_evt_length;

  std::vector<std::array<int, 64>> adcOffsets;
  std::vector<uint64_t> oldtime_d;
  std::vector<uint64_t> upperTime48_d;
  std::vector<uint64_t> upperTime63_d;
  std::vector<std::array<int, 64>> discriminatorTimes;

  Int_t* pMbsData;
  Int_t* pMbsEnd;


  uint64_t oldtime;
  uint64_t upperTime48;
  uint64_t upperTime63;

  UnpackProcStats stats;
  uint64_t GetTime(uint64_t lowTS);
  uint64_t GetTime_d(int fee, uint64_t lowTS);

  bool BuildAIDAEvent(TGo4MbsSubEvent* psubevt, Int_t word1, Int_t word2);
  bool BuildAIDAADCEvent(TGo4MbsSubEvent* psubevt, Int_t word1, Int_t word2);
  bool BuildAIDAInfoEvent(TGo4MbsSubEvent* psubevt, Int_t word1, Int_t word2);

public:
  AIDA_Detector_System();
  ~AIDA_Detector_System();

  ///Nic additions

  void CorrectMultiplexer(ADCDataItem&);
  virtual void PrintStatistics();

  void Process_MBS(TGo4MbsSubEvent* psubevt);


  void Process_MBS(int*){};
  void get_Event_data(Raw_Event*);
  int* get_pdata();

  void write(){return;};
  bool calibration_done(){return false;}

  void set_Gain_Match_Filename(std::string){return;};
};



#endif
