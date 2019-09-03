DESPEC Analysis code v4.2 22.08.19 AM
General
Code is on the despec account in /u/despec/DESPEC-Online-April_v4.2/
Do ‘make all’ to compile. If you make changes in the code then just do ‘make'
•	Systems included: AIDA, bPlastic (VME +Scalar), FATIMA, GALILEO, Finger (TAMEX)
•	Systems selected in /Configuration_Files/Used_Systems.txt – Can be turned on/off as desired
•	Proc IDs set in PrcID_to_Det_sys.txt (bplastic was changed from 20 (offline tests) to 80 for the online run)
•	Energy and Time gates for FATIMA and Plastic can be set in Configuration_Files/Correlations.dat
•	Gain matching for energy (time) can be found in /Configuration_Files/
Plastic_VME_Energy_Calibration.dat
Plastic_VME_Time_Calibration.dat (Raw times)
FATIMA_Energy_Calibration.dat
FATIMA_Ref_Time_Calibration.dat (Reference-Fatima time) 
GALILEO_Energy_Calibration.txt
/Calibration_FINGER/ (each channel own calib file)

Unpacker (EventUnpackProc)
•	Each data unpack step is done in the respective Detector_System class. The data is combined in the event builder named: EventUnpackProc.cxx  The raw histograms are here. 
•	All getters for the various detector systems (i.e. RAW defined in Raw_Event.cxx) should be there (Energy, Time etc).
•	fOutput sends variables to the root tree (UnpackxTree) defined in EventUnpackStore

•	In the auto saved histogram .root file after running via Go4, the histograms are placed into folders according to their respective detector system plus one folder for bplastic-fatima correlations. The ‘Info’ tab in the Go4 histogram browser gives specific details.

UnpackxTree
•	FATIMA, Plastic and galileo outputs (given below) can be filled into root trees. (use the 
–store rootfilename.root when running go4analysis in batch mode)
•	Once the data has been stored to the root tree, using the stored unpack tree as the source for the AnlProc and disabling the UnpackProc will speed up the code(needs testing)
For example to see the Fatima raw energy of channel 6 then open the saved root file (i.e. rootfilename.root) and type 
UnpackxTree->Draw(“fFat_QDC_E >>(20000,0,40000)”,”fFat_TDC_ Multiplicity >0”)

Analysis Procedure (EventAnlProc)
•	The outputs from the UnpackxTree tree are sent to the EventAnlProc on lines 186-360 and variables are defined here to use throughout the Analysis procedure. 
•	Here, the main histograms are filled and (for now) correlations are done. You can find out information on the various histograms in the info tab of the histogram browser in Go4. You can change the name of the histograms in the .cxx file as you wish. 
•	Calibrations (except for finger) are applied in EventAnlProc.cxx. 
•	The Finger fine calibration can be done for the first few events by activating FORCE 1 and ‘ONLINE 1 in the file /Configuration_Files/FINGER_CALIB_FILE’ Then the data can be run and once finished set back to 0, 0. 

AnalysisxTree
The analysis procedure (calibrated) data from EventAnlProc.cxx can be now stored into root trees. This can be activated using:
go4analysis -file /data.local/d001f012_0001.lmd  -step 1 -store analysistree.root -asf histos.root
			      
                                             Filename                    Step	  Name of root tree	Histo autosave file
						(0 unpack, 1 Analysis)			
  
---------------------------------------------------------------------------------------------------------------------
UnpackxTree Leaves: Activate output with –step 1 –store ‘filename.root’ in batch mode

fevent_number = Event number
fProcID[ID] =  process ID
fUsed_Systems[ID] = Used system defined in ‘/Configuration_Files/Used_Systems.txt’
  
        AIDA
All variables (expect hits) stored in one leaf of the UnpackStore branch and separated in the EventAnlProc

        bPlastic            
        fbPlas_VME_firedQDC = Plastic fired QDC channels;    
        fbPlas_VME_QDC_E[‘fired QDC’] = Plastic raw QDC data;
        fbPlas_VME_QDC_ID[‘fired QDC’] = Plastic QDC ID;
         
        fbPlas_VME_firedTDC = Plastic fired TDC channels ;    
        fbPlas_VME_TDC_ID[‘fired TDC’] = Plastic TDC ID ;
        fbPlas_VME_TDC_TS[‘fired TDC’][‘TDC ID’] = Plastic raw TDC data;
        fbPlas_VME_TDC_Multiplicity[‘TDC ID’] = Plastic TDC multiplicity;
         
        fScalar_fired = fired Scalar;
        fScalar_ID = Scalar ID;

        Fatima
fFat_firedQDC = Fatima fired QDC channels
fFat_QDC_ID[‘fired QDC’] = Fatima QDC ID;
fFat_QDC_E[‘fired QDC’] = Fatima raw QDC Data;

fFat_firedTDC = fired Fatima TDC;
fFat_TDC_ID[‘fired TDC’] = Fatima TDC ID;
fFat_TDC_Multiplicity[‘TDC ID’] = Fatima TDC multiplicity;
fFat_TDC_TS[‘fired TDC’][‘TDC ID’] = Fatima TDC data;       
                          
            fS41[fired TDC’] = S41 trigger;



UnpackxTree cont.
Galileo   
   fGal_fired = Galileo fired channels
   fGal_ID[‘fired gal’] = Gal ID ;
   fGal_E[‘fired gal’] = Gal energy ;
   fGal_T[‘fired gal’] = Gal time;
   fGal_Pileup = Gal pileup events (needs testing)
        
Finger
ffing_tamexhits   = Fired TAMEX hits    
ffing_leadHits[‘tamex hits’] = lead hits
ffing_trailHits[‘tamex hits’] = trail hits
ffing_iterator[‘tamex hits’]  = loops over leading/trailing edges
ffing_chID[‘tamex hits][ffing iterator] = Finger Tamex ID 

(Channels defined in Configuration_Files/Finger_allocation.txt)
ffing_Lead_Phys_Chan[‘tamex hits’][fing iterator] = Finger Leading ID 
ffing_Lead_T[‘tamex hits’][fing iterator] = Finger Leading Time
ffing_Trail_Phys_Chan[‘tamex hits’][fing iterator]  = Finger Trailing ID 
ffing_Trail_T[‘tamex hits][fing iterator] = Finger Trailing Time
ffing_TOT[‘tamex hits’][fing iterator]   = Time/Threshold (computed in Raw_Event.cxx)

---------------------------------------------------------------------------------------------------------------------
AnalysisxTree Leaves: Activate output with –step 1 –store ‘filename.root’ in batch mode
bPlastic 
      pbPlas_QDCFired = Plastic fired QDC channels;     
      pbPlas_QDCID[‘fired QDC’] = Plastic QDC ID;
      pbPlas_QDCGainMatch_i[‘fired QDC’] = Plastic gainmatched QDC data;

      pbPlas_TDCFired = Plastic fired TDC channels ;    
      pbPlas_TDC_Multiplicity[‘TDC ID’] = Plastic TDC multiplicity;
      pbPlas_TDCID[‘fired TDC’] = Plastic TDC ID ;
      pbPlasTDC_T[‘fired TDC’] = Raw  Plastic TDC     
      pbPlas_SC41_dT[‘fired TDC’] = Gainmatched Time difference Sci41-bPlas Ch.x data (second  col of Plastic_VME_Ref_Time_calibration.dat)
      pbPlas_SiPM_dT_Calib[‘fired TDC’] = Gainmatched time difference bPlas Ch.1 - bPlas Ch.x data (first col of Plastic_VME_Ref_Time_calibration.dat)


   Fatima  
      pFat_WR = Fatima/bPlastic white rabbit
      pFat_QDCFired = Fatima fired QDC channels
      pFat_QDCID[‘fired QDC’] = Fatima QDC ID;
      pFat_QDCGainMatch[‘fired QDC’] = Fatima gainmatched QDC Data;

      pFat_TDCFired = fired Fatima TDC;
      pFat_TDCID[‘fired TDC’] = Fatima TDC ID;
      pFat_TDC_T[‘pFatTDCID’] = Fatima Raw TDC Data;

      pFat_SC41_dT_Calib[‘pFatTDCID’] = Gainmatched Time difference Sci41-Fatima Ch.x data ; (third col of Fatima_Ref_Time_calibration.dat)

      pFat_TDC_Multipl_perCh[‘pFatTDCID’] Fatima Multiplicity (per channel)
      pFat_Ch_dT[‘pFatTDCID’] = Gainmatched time difference Fatima Ch.0 - Fatima Ch.x; (second  col of Fatima_Ref_Time_calibration.dat)
      pFat_Ch0_TDC =  Raw Fatima TDC data Ch.0;
     
Galileo    
      pGalFired = Galileo fired channels (multiplicity);
      pGal_WR = Galileo White rabbit;
      pGalID[‘GalID’] = Gal ID ;
      pGalT[‘GalID’] = Gal time;
      pGalE_Cal_i[‘GalID’] = Gal calibrated energy ;
      pGal_dT = Time difference between two fired channels;
      pGalE_Addback = Added back energy sum (>1 hit/event);
      
Finger
      pFing_firedTamex = Finger Fired TAMEX module    ; 
      pFing_iterator[‘tamexID’] = Finger lead/trail hits;
      pFing_LeadChan[‘tamexID’][‘l/t hits’] = Finger leading channel number (mapped);
      pFing_leadT[‘tamexID’][ ‘l/t hits’] = Finger Leading Time;
      pFing_LeadDiff[‘l/t hits’] = Finger Lead-Lead Time;
      pFing_LeadPlus[‘l/t hits’] Finger Lead+Lead Time;;
      pFing_SC41_diff[‘l/t hits’] = Finger Sci41-Lead Time;;
      pFing_TOT[‘tamexID’][ ‘l/t hits’] = Finger Time/Threshold;     
      pFing_downData ;
      pFing_total_time;
