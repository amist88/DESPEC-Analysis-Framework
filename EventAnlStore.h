// $Id: EventAnlStore.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef TSCNCALEVENT_H
#define TSCNCALEVENT_H

#include "TGo4EventElement.h"
//#include "TSCNUnpackEvent.h"
#include "EventUnpackStore.h"

class EventAnlStore : public TGo4EventElement {
   public:
      EventAnlStore() : TGo4EventElement() {}
      EventAnlStore(const char* name) : TGo4EventElement(name) {}
      virtual ~EventAnlStore() {}

      virtual void  Clear(Option_t *t="");

//       Float_t frDataC1[SCN_NUM_CHAN]; //For QDC 1
//       Float_t frDataC2[SCN_NUM_CHAN]; //For QDC 2
//       Double_t frDataC3[SCN_NUM_CHAN]; //For TDC

      ///Plastic AnlProc Outputs
      Int_t    pbPlas_QDCFired; 
      Int_t    pbPlas_QDCID[32];
      Double_t pbPlas_QDCGainMatch_i[32];
      Int_t    pbPlas_TDCFired;
      Int_t    pbPlas_TDC_Multiplicity[32]; 
      Int_t    pbPlas_TDCID[50];
      Double_t pbPlasTDC_T[32];
      Double_t pbPlas_SC41_dT[32];
      Double_t pbPlas_SiPM_dT_Calib[32];
     
      Long64_t pFat_WR;
      Int_t    pFat_QDCFired;
      Int_t    pFat_QDCID[50];
      Double_t pFat_QDCGainMatch[50];
      Int_t    pFat_TDCFired;
      Int_t    pFat_TDCID[50];
      Double_t pFat_TDC_T[50];
      Double_t pFat_SC41_dT_Calib[50];
      Int_t    pFat_TDC_Multipl_perCh[50];
      Double_t pFat_Ch_dT[50];
      Double_t pFat_Ch0_TDC;
      
      Int_t pGalFired;
      Long64_t pGal_WR;
      Int_t pGalID[32];
      Long64_t pGalT[32];
      Double_t pGalE_Cal_i[32];
      Long64_t pGal_dT;
      Double_t pGalE_Addback;
      
      Int_t pFing_firedTamex; 
      Int_t pFing_iterator[4];
      Int_t pFing_LeadChan[4][32];
      Double_t pFing_leadT[4][32];
      Double_t pFing_LeadDiff[100];
      Double_t pFing_LeadPlus[100];
      Double_t pFing_SC41_diff[100];
      Double_t pFing_TOT[4][32];     
      Double_t pFing_downData;
      Double_t pFing_total_time;

   ClassDef(EventAnlStore,1)
};
#endif //TSCNCALEVENT_H



