// $Id: TSCNAnlEvent.h 755 2011-05-20 08:04:11Z linev $
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

#ifndef TSCNANLEVENT_H
#define TSCNANLEVENT_H

#include "TGo4EventElement.h"
//#include "TSCNUnpackEvent.h"
#include "EventAnlStore.h"

class EventCorrelStore : public TGo4EventElement {
   public:
      EventCorrelStore() : TGo4EventElement() {}
      EventCorrelStore(const char* name) : TGo4EventElement(name) {}
      virtual ~EventCorrelStore() {}

      virtual void  Clear(Option_t *t="");



      Double_t Data1[SCN_NUM_CHAN]; 
      Double_t Data2[SCN_NUM_CHAN]; 
      Double_t Data3[SCN_NUM_CHAN]; 


      

   ClassDef(EventCorrelStore,1)
};
#endif //TSCNANLEVENT_H



