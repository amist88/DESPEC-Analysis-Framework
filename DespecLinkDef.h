// $Id: SCNLinkDef.h 500 2009-11-04 18:42:41Z goofy $
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

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TSCNParameter+;

#pragma link C++ class DESPECAnalysis+;
#pragma link C++ class EventUnpackProc+;
#pragma link C++ class EventUnpackStore+;
#pragma link C++ class EventUnpackFact+;


#pragma link C++ class EventAnlStore+;
#pragma link C++ class EventAnlProc+;
#pragma link C++ class CalibParameter+;
#pragma link C++ class CorrelParameter+;
#pragma link C++ class EventCorrelStore+;
#pragma link C++ class EventCorrelProc+;

#pragma link C++ class AIDA_Event+;
#pragma link C++ class AidaEvent+;
#pragma link C++ class AidaHit+;
/*#pragma link C++ class TXRSAnalysis+;
#pragma link C++ class TXRSBasicProc+;
#pragma link C++ class TXRSUnpackEvent+;
#pragma link C++ class TXRSUnpackProc+;
#pragma link C++ class TXRSSortEvent+;
#pragma link C++ class TXRSSortProc+;
#pragma link C++ class TXRSCalibrEvent+;
#pragma link C++ class TXRSCalibrProc+;
#pragma link C++ class TXRSAnlEvent+;
#pragma link C++ class TXRSAnlProc+;
#pragma link C++ class TXRSUserEvent+;
#pragma link C++ class TXRSUserProc+;*/


#pragma link C++ class TXRSParameter+;
#pragma link C++ class TModParameter+;
#pragma link C++ class Map1+;
#pragma link C++ class TMWParameter+;
#pragma link C++ class TTPCParameter+;
#pragma link C++ class TMUSICParameter+;
#pragma link C++ class TSCIParameter+;
#pragma link C++ class TIDParameter+;
#pragma link C++ class TSIParameter+;
#pragma link C++ class TMRTOFMSParameter+;
#pragma link C++ class AidaMover+;
#endif



