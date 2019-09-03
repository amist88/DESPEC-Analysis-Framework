#include "AIDA_Detector_System.h"
#include "AIDA_Data_Types.h"
#include "AIDA_Event.h"
#include "TGo4MbsEvent.h"
#include <fstream>
#include "TGo4MbsSubEvent.h"
#include "TROOT.h"
#include <cstdlib>
#include "TH1.h"
#include "TTree.h"
#include <time.h>
#include "Timestamp.h"
#include "TAidaConfiguration.h"
#include "TGo4Log.h"
#include <vector>

using namespace std;

static const int FeeToStrip[64] = {
  62, 63, 59, 60, 61, 56, 57, 58, 52, 53, 54, 55, 49, 50, 51, 45,
  46, 47, 48, 42, 43, 44, 38, 39, 40, 41, 35, 36, 37, 31, 32, 33,
  34, 28, 29, 30, 24, 25, 26, 27, 21, 22, 23, 17, 18, 19, 20, 14,
  15, 16, 10, 11, 12,  7,  3,  0,  8,  4,  1,  9,  5,  2, 13,  6
};
//---------------------------------------------------------------

AIDA_Detector_System::AIDA_Detector_System()
{
  TGo4Log::Info("DESPEC Analysis: Creating AIDA Configuration");
  TAidaConfiguration::Create("Configuration_Files/AIDA.txt");

  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();

  adcOffsets.resize(conf->FEEs());
  discriminatorTimes.resize(conf->FEEs());
  oldtime_d.resize(conf->FEEs());
  upperTime48_d.resize(conf->FEEs());
  upperTime63_d.resize(conf->FEEs());


  Energy = new double[10000];
  FEE_ID =  new int[10000];
  Ch_ID = new int[10000];
  Time = new ULong64_t [10000];
  FastTime = new ULong64_t [10000];
  HighE_veto = new bool[10000];
  Side = new int[10000];
  Strip = new int[10000];
  EvtID = new int[10000];
  AdcData = new int[10000];

  Hits = 0;
  for (int i = 0; i < 10000; i++)
  {
    Energy[i] = 0;
    FEE_ID[i] = 0;
    Ch_ID[i] = 0;
    Time[i] = 0;
    FastTime[i] = 0;
    HighE_veto[i] = 0;
    Side[i] = 0;
    Strip[i] = 0;
    EvtID[i] = 0;
    AdcData[i] = 0;
  }

  for (int i = 0; i < conf->FEEs(); ++i)
  {
    oldtime_d[i] = 0;
    upperTime48_d[i] = 0;
    upperTime63_d[i] = 0;
    for (int j = 0; j < 64; ++j)
    {
      adcOffsets[i][j] = 0;
      discriminatorTimes[i][j] = 0;
    }
  }

  upperTime48 = 0;
  upperTime63 = 0;
  //timewarpAmount = 0;
  memset(&stats, 0, sizeof(UnpackProcStats));

  std::ifstream fs("Configuration_Files/AIDA_offsets.txt");
  fs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  while (fs)
  {
    if (fs.peek() == '#')
    {
      fs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    int fee, channel;
    double offset;
    fs >> fee >> channel >> offset;
    if (!fs) break;
    adcOffsets[fee-1][channel-1] = offset;
    fs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  TGo4Log::Info("AIDA: Loaded ADC Offsets");
  fs.close();
}

//---------------------------------------------------------------

AIDA_Detector_System::~AIDA_Detector_System(){


  delete[] Energy;
  delete[] FEE_ID;
  delete[] Ch_ID;
  delete[] Time;
  delete[] FastTime;
  delete[] HighE_veto;
  delete[] Side;
  delete[] Strip;
  delete[] EvtID;
  delete[] AdcData;

}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
void AIDA_Detector_System::get_Event_data(Raw_Event* RAW){
  //
  RAW->set_DATA_AIDA(Energy, FEE_ID, Ch_ID, Time, Hits, HighE_veto, Side, Strip, EvtID, FastTime, AdcData);
}
//---------------------------------------------------------------

void AIDA_Detector_System::Process_MBS(TGo4MbsSubEvent* psubevt)
{
  pdata = psubevt->GetDataField();
  Int_t lwords = psubevt->GetIntLen();
  pdata_start = pdata;

  sub_evt_length  = ((psubevt->GetDlen() - 2) / 2);

  if (*pdata++ == 0x200)
  {

    uint64_t wr = 0;
    uint64_t bits = *pdata++ & 0xFFFF;
    wr |= bits;
    bits = *pdata++ & 0xFFFF;
    wr |= (bits << 16);
    bits = *pdata++ & 0xFFFF;
    wr |= (bits << 32);
    bits = *pdata++ & 0xFFFF;
    wr |= (bits << 48);

    upperTime63 = std::max(upperTime63, wr >> 48);
    upperTime48 = std::max(upperTime48, (wr >> 28) & 0x000FFFFF);

    if(!stats.Start) stats.Start = wr;

    if ((lwords - 5 ) & 0x1)
    {
      std::cout << "AIDA data size is incorrect, ignoring data" << std::endl;
    }
    Hits = 0;
    while ((pdata - pdata_start) < (sub_evt_length)) // subevent loop
    {
      stats.MBSEvents++;

      Int_t word1 = *pdata++;
      Int_t word2 = *pdata++;

      //           printf( "word1 %d (0x%x)\n" , word1);
      //           printf( "word2 %d (0x%x)\n" , word2);

      BuildAIDAEvent(psubevt, word1, word2);
    }//end  while

  }
}

bool AIDA_Detector_System::BuildAIDAEvent(TGo4MbsSubEvent* psubevt, Int_t word1, Int_t word2)
{
  int type = (word1 >> 30) & 0x3;
  int lowTS = (word2 & 0x0FFFFFFF);
  stats.Stop = GetTime(lowTS);

  stats.AIDAEvents++;

  if (type == 3)
  {
    return BuildAIDAADCEvent(psubevt, word1, word2);
  }
  else if (type == 2)
  {
    stats.Type.Info++;
    return BuildAIDAInfoEvent(psubevt, word1, word2);
  }
  else
  {
    std::cout << "Unknown event type " << type << std::endl;
    return false;
  }
}

bool AIDA_Detector_System::BuildAIDAADCEvent(TGo4MbsSubEvent* psubevt, Int_t word1, Int_t word2)
{
  AidaEvent evt;

  // evt.ID = 3;
  EvtID[Hits] = 3;

  int channelID = (word1 >> 16) & 0xFFF;
  int feeID = (channelID >> 6) & 0x3F;
  channelID &= 0x3F;
  int lowTS = (word2 & 0x0FFFFFFF);
  int data = (word1 & 0xFFFF);
  int veto = (word1 >> 28) & 0x1;

  evt.Data = data;

  ///For RAW event (to unpack procedure stage)
  FEE_ID[Hits] = feeID;
  Ch_ID[Hits] = channelID;
  Time[Hits] = GetTime_d(feeID, lowTS);
  HighE_veto[Hits] = (veto == 1);
  FastTime[Hits] = discriminatorTimes[feeID][channelID];
  AdcData[Hits] = data;

  //cout <<" Ch_ID[hi] " <<Ch_ID[Hits]<< " Time[hi] " << Time[Hits] << " hits " <<Hits<< endl;
  // }
  /* Real DSSD Mapping */

  TAidaConfiguration const* conf = TAidaConfiguration::GetInstance();
  evt.DSSD = conf->FEE(feeID).DSSD;

  //Side, Front or back determines the polarity
  if (evt.DSSD > 0)
  {
    evt.Side = conf->FEE(feeID).Side;

    Side[Hits] = conf->FEE(feeID).Side;

    if (conf->FEE(feeID).High)
    {
      evt.Strip = 127 - FeeToStrip[channelID];
      Strip[Hits] = 127 - FeeToStrip[channelID];
    }
    else
    {
      evt.Strip = FeeToStrip[channelID];
      Strip[Hits] = FeeToStrip[channelID];
    }

    evt.Intensity = (evt.Data - 32768) * evt.Side;
    // if (evt.HighEnergy)
    // test->Fill(evt.Strip + (evt.Side == -1 ? 0 : 1) * 128);

    if (veto == 1)
    {
      evt.Energy = (evt.Intensity - adcOffsets[feeID][channelID]) * 0.7; // Energy in MeV
      Energy[Hits] = (evt.Intensity - adcOffsets[feeID][channelID]) * 0.7; // Energy in MeV
    }
    else
    {
      evt.Energy = (evt.Intensity- adcOffsets[feeID][channelID]) * 0.7; // Energy in keV
      Energy[Hits] = (evt.Intensity- adcOffsets[feeID][channelID]) * 0.7; // Energy in keV
    }
    // delete data without an offset (bad channels)		    }
    if (!adcOffsets[feeID][channelID])
    {
      //Side[Hits] = -10;
      Energy[Hits] = 0;
    }
  }

  Hits++;
  return true;
}


//Don't need to output this
bool AIDA_Detector_System::BuildAIDAInfoEvent(TGo4MbsSubEvent* psubevt, Int_t word1, Int_t word2)
{
  int moduleNumber = (word1 >> 24) & 0x3F;
  int infoCode = (word1 >> 20) & 0x000F;
  int infoField = (word1 & 0x000FFFFF);
  int lowTS = (word2 & 0x0FFFFFFF);

  // We just decode Timestamp events here

  // SYNC (PAUSE/RESUME too(?)) for bits 48-28 of the timestamp
  if (infoCode == 2) {
    stats.InfoType.Pause++;
  }
  else if (infoCode == 3) {
    stats.InfoType.Resume++;
  }
  else if (infoCode == 4)
  {
    stats.InfoType.Sync48++;
    upperTime48 = infoField;
    upperTime48_d[moduleNumber] = infoField;
    if (infoCode == 4 && stats.Start == 0) stats.Start = GetTime(lowTS);
  }
  // WR-SYNC for bits 63-48 of the timestamp
  else if (infoCode == 5)
  {
    stats.InfoType.Sync63++;
    upperTime63 = infoField;
    upperTime63_d[moduleNumber] = infoField;
  }
  else if (infoCode == 6)
  {
    stats.InfoType.Discriminator++;
    int adc = ((infoField >> 16) & 0xF) - 1;
    int hits = infoField & 0xFFFF;
    int64_t time = GetTime(lowTS);
    for (int i = 0; i < 16; ++i)
    {
      if (hits & (1 << i))
      {
        discriminatorTimes[moduleNumber][adc*16 + i] = time;
      }
    }
  }
  else
  {
    stats.InfoType.Unknown++;
  }

  return false;
}

uint64_t AIDA_Detector_System::GetTime(uint64_t lowTS)
{
  return (upperTime63 << 48) | (upperTime48 << 28) | lowTS;
}

uint64_t AIDA_Detector_System::GetTime_d(int fee, uint64_t lowTS)
{
  return GetTime(lowTS);
  return (upperTime63_d[fee] << 48) | (upperTime48_d[fee] << 28) | lowTS;
}

void AIDA_Detector_System::PrintStatistics()
{
  std::cout << "======================================" << std::endl;
  std::cout << "AIDA Unpacker Analysis" << std::endl;
  std::cout << "--------------------------------------" << std::endl;
  std::cout << "First event: " << wr_to_string(stats.Start);
  std::cout << "Last event : " << wr_to_string(stats.Stop);
  double duration = wr_to_duration(stats.Start, stats.Stop);
  std::cout << "--------------------------------------" << std::endl;
  std::cout << "MBS Events   : " << stats.MBSEvents << " (" << frequency(stats.MBSEvents, duration) << ")" <<  std::endl;
  std::cout << "AIDA Words   : " << stats.AIDAEvents << " (" << frequency(stats.AIDAEvents, duration) << ")" <<  std::endl;
  std::cout << "ADC Words    : " << stats.Type.ADC << " (" << frequency(stats.Type.ADC, duration) << ")" <<  std::endl;
  std::cout << "Info Words   : " << stats.Type.Info << " (" << frequency(stats.Type.Info, duration) << ")" <<  std::endl;
  std::cout << "SYNC48 Words : " << stats.InfoType.Sync48 << " (" << frequency(stats.InfoType.Sync48, duration) << ")" <<  std::endl;
  std::cout << "SYNC63 Words : " << stats.InfoType.Sync63 << " (" << frequency(stats.InfoType.Sync63, duration) << ")" <<  std::endl;
  std::cout << "PAUSE Words  : " << stats.InfoType.Pause << " (" << frequency(stats.InfoType.Pause, duration) << ")" <<  std::endl;
  std::cout << "RESUME Words : " << stats.InfoType.Resume << " (" << frequency(stats.InfoType.Resume, duration) << ")" <<  std::endl;
  std::cout << "Discrim Words: " << stats.InfoType.Discriminator << " (" << frequency(stats.InfoType.Discriminator, duration) << ")" <<  std::endl;
  std::cout << "Unknown Words: " << stats.InfoType.Unknown << " (" << frequency(stats.InfoType.Unknown, duration) << ")" <<  std::endl;
  std::cout << "Timewarps    : " << stats.TimeWarps << " (" << frequency(stats.TimeWarps, duration) << ")" <<  std::endl;
  std::cout << "======================================" << std::endl;
}


int* AIDA_Detector_System::get_pdata(){return pdata;}

///End of Nic Unpack here
//---------------------------------------------------------------
