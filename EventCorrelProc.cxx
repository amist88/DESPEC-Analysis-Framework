// $Id: EventCorrelProc.cxx 754 2011-05-18 11:04:52Z adamczew $
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

#include "EventCorrelProc.h"

#include <cstdlib>
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "TGo4WinCond.h"
#include "TGo4CondArray.h"
#include "TGo4Analysis.h"

#include "TGo4Picture.h"

#include "EventCorrelStore.h"
#include "EventUnpackStore.h"
#include "EventAnlStore.h"
#include "TSCNParameter.h"

//#include "BrentMethod.c"
//#include "Linearizator.h" //pos reconstruction
//#include "GoldenMethod.c"

#define range 6


//-----------------------------------------------------------
EventCorrelProc::EventCorrelProc() :
  TGo4EventProcessor()/*,
  fParam1(0),
  fTimeDiff(0),
  fGatedHist(0),
  fCoincQ1A1(0),
  fCoincQ1T1(0),
  fconHis1(0)*/
{
}
//-----------------------------------------------------------
EventCorrelProc::EventCorrelProc(const char* name) :
   TGo4EventProcessor(name)
{
  /*cout << "**** EventCorrelProc: Create" << endl;
  fParam1 = (TSCNParameter*)  GetParameter("SCNParameter");
 

  
  fTimeDiff = MakeTH1('D',"TimeDiff","TimeDiff",20000,-10000,10000);
  fGatedHist = MakeTH1('D',"GatedHist","GatedHist",4000,1,4000);

  fCoincQ1A1 = MakeTH2('D',"CoincE","CoincQDC1ADC1", 200, 0.0, 2000.0, 200, 0.0, 2000.0);     
  fCoincQ1T1 = MakeTH2('D',"CoincET","CoincQDC1TDC1", 2000, -2000.0, 2000.0, 3000, 0.0, 3000.0);     
  fconHis1   = MakeWinCond("cHis1", -200, 200, "TimeDiff"); 
  fWinCon1   = MakeWinCond("wincon1", -100, 100, 0, 3000);

  fconHis1->Enable();
  fWinCon1->Enable();*/


 }
//-----------------------------------------------------------
EventCorrelProc::~EventCorrelProc()
{
  cout << "**** EventCorrelProc: Delete" << endl;
}
//-----------------------------------------------------------
Bool_t EventCorrelProc::BuildEvent(TGo4EventElement* dest)
{
  Bool_t isValid=kFALSE; // validity of output event
  
  EventAnlStore* inp_evt  = (EventAnlStore*) GetInputEvent();
  EventCorrelStore* out_evt = (EventCorrelStore*) dest;

 
  
  // see comments in UnpackProc
  if((inp_evt==0) || !inp_evt->IsValid()){ // input invalid
    out_evt->SetValid(isValid); // invalid
    return isValid; // must be same is for SetValid
  }
  isValid=kTRUE;
  

  Double_t tdc_event[SCN_NUM_CHAN];
  Double_t tdc_event_diff[SCN_NUM_CHAN];
  // unused // Double_t timediffdata;



//   for(Int_t j=0;j<SCN_NUM_CHAN;j++)
//     {
// 
//       out_evt->Data1[j] = inp_evt->frDataC1[j];
//       out_evt->Data2[j] = inp_evt->frDataC2[j];
//       out_evt->Data3[j] = inp_evt->frDataC3[j];
//    }



  

 
 out_evt->SetValid(isValid);
  return isValid;
}

 /**----------------------------------------------------------------------------------------------**/
 /**----------------------------------     FRS-AIDA (Implanted ion)   -------------------------**/
 /**----------------------------------------------------------------------------------------------**/

 /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------- FRS-AIDA-bPlastic  (Beta Decay)  -------------------**/
 /**----------------------------------------------------------------------------------------------**/
 
  /**----------------------------------------------------------------------------------------------**/
  /**----------------------------     FRS-AIDA-bPlastic-FATIMA/GALILEO  (Isomers/Beta decay spectroscopy) --**/
  /**----------------------------------------------------------------------------------------------**/

 /**----------------------------------------------------------------------------------------------**/
 /**--------------------------------  FRS-AIDA-bPlastic-FATIMA  (Beta Decay Timing)  -----**/
 /**----------------------------------------------------------------------------------------------**/