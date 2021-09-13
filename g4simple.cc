// standard library includes
#include <iostream>
#include <regex>
#include <string>
#include <utility>
#include <vector>

// external includes
#include "G4EventManager.hh"
#include "G4GDMLParser.hh"
#include "G4GeneralParticleSource.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysListFactory.hh"
#include "G4PVPlacement.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgrMessenger.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"
#include "G4UnitsTable.hh"
#include "G4UserSteppingAction.hh"
#include "G4VisExecutive.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

// experiment includes
#include "g4hdf5.hh"

class G4SimpleSteppingAction : public G4UserSteppingAction, public G4UImessenger {
  protected:
    G4UIcommand       * fVolIDCmd;
    G4UIcmdWithAString* fOutputFormatCmd;
    G4UIcmdWithAString* fOutputOptionCmd;
    G4UIcmdWithABool  * fRecordAllStepsCmd;
    G4UIcmdWithAString* fSilenceOutputCmd;
    G4UIcmdWithAString* fAddOutputCmd;

    enum EOption {kStepWise, kEventWise};
    EOption fOption;
    bool fRecordAllSteps;

    std::vector<std::pair<std::regex, std::string>> fPatternPairs;
    G4int                                           fNEvents;
    G4int                                           fEventNumber;
    std::vector<G4int>                              fPID;
    std::vector<G4int>                              fTrackID;
    std::vector<G4int>                              fParentID;
    std::vector<G4int>                              fStepNumber;
    std::vector<G4double>                           fKE;
    std::vector<G4double>                           fEDep;
    std::vector<G4double>                           fX;
    std::vector<G4double>                           fY;
    std::vector<G4double>                           fZ;
    std::vector<G4double>                           fLX;
    std::vector<G4double>                           fLY;
    std::vector<G4double>                           fLZ;
    std::vector<G4double>                           fPdX;
    std::vector<G4double>                           fPdY;
    std::vector<G4double>                           fPdZ;
    std::vector<G4double>                           fT;
    std::vector<G4int>                              fVolID;
    std::vector<G4int>                              fIRep;
    std::map<G4VPhysicalVolume*, int>               fVolIDMap;

    // bools for setting which fields to write to output
    G4bool fWEv, fWPid, fWTS, fWKE, fWEDep, fWR, fWLR, fWP, fWT, fWV;

  public:
    G4SimpleSteppingAction() :
        fNEvents    (0   ),
        fEventNumber(0   ),
        fWEv        (true),
        fWPid       (true),
        fWTS        (true),
        fWKE        (true),
        fWEDep      (true),
        fWR         (true),
        fWLR        (true),
        fWP         (true),
        fWT         (true),
        fWV         (true) {
      ResetVars();

      fVolIDCmd = new G4UIcommand("/g4simple/setVolID", this);
      fVolIDCmd->SetParameter(new G4UIparameter("pattern"    , 's', false));
      fVolIDCmd->SetParameter(new G4UIparameter("replacement", 's', false));
      fVolIDCmd->SetGuidance(
        "Volumes with name matching [pattern] will be given volume ID "
        "based on the [replacement] rule. Replacement rule must produce an integer. "
        "Patterns which replace to 0 or -1 are forbidden and will be omitted.");

      fOutputOptionCmd = new G4UIcmdWithAString("/g4simple/setOutputOption", this);
      std::string candidates = "stepwise eventwise";
      fOutputOptionCmd->SetCandidates(candidates.c_str());
      fOutputOptionCmd->SetGuidance("Set output option:");
      fOutputOptionCmd->SetGuidance("  stepwise: one row per step");
      fOutputOptionCmd->SetGuidance("  eventwise: one row per event");
      fOption = kStepWise;

      fRecordAllStepsCmd = new G4UIcmdWithABool("/g4simple/recordAllSteps", this);
      fRecordAllStepsCmd->SetParameterName("recordAllSteps", true);
      fRecordAllStepsCmd->SetDefaultValue(true);
      fRecordAllStepsCmd->SetGuidance("Write out every single step, not just those in sensitive volumes.");
      fRecordAllSteps = false;

      fSilenceOutputCmd = new G4UIcmdWithAString("/g4simple/silenceOutput", this);
      fSilenceOutputCmd->SetGuidance("Silence output fields");
      fAddOutputCmd = new G4UIcmdWithAString("/g4simple/addOutput", this);
      fAddOutputCmd->SetGuidance("Add output fields");
      candidates = "event pid track_step kinetic_energy energy_deposition position local_position momentum time volume all";
      fSilenceOutputCmd->SetCandidates(candidates.c_str());
      fAddOutputCmd    ->SetCandidates(candidates.c_str());
    }

    G4VAnalysisManager* GetAnalysisManager() {
      return G4Hdf5::G4AnalysisManager::Instance();
    }

    ~G4SimpleSteppingAction() { 
      G4VAnalysisManager* man = GetAnalysisManager();
      if (man->IsOpenFile()) {
        if (fOption == kEventWise && fPID.size() > 0) {
          WriteRow();
        }
        man->Write();
        man->CloseFile();
      }

      delete man;
      delete fVolIDCmd;
      delete fOutputFormatCmd;
      delete fOutputOptionCmd;
      delete fRecordAllStepsCmd;
    } 

    void SetNewValue(G4UIcommand* command, G4String newValues) {
      if (command == fVolIDCmd) {
        std::istringstream iss(newValues);
        std::string pattern;
        std::string replacement;
        iss >> pattern >> replacement;
        fPatternPairs.push_back(std::pair<std::regex, std::string>(std::regex(pattern), replacement));
      }
      if (command == fOutputFormatCmd) {
        // also set recommended options
        // override option by subsequent call to /g4simple/setOutputOption
        fOption = kStepWise;
        GetAnalysisManager(); // call once to make all of the /analysis commands available
      }
      if (command == fOutputOptionCmd) {
        if (newValues == "stepwise" ) fOption = kStepWise ;
        if (newValues == "eventwise") fOption = kEventWise;
      }
      if (command == fRecordAllStepsCmd) {
        fRecordAllSteps = fRecordAllStepsCmd->GetNewBoolValue(newValues);
      }
      if (command == fSilenceOutputCmd) {
        G4bool all = (newValues == "all");
        if (all || newValues == "event"            ) fWEv   = false;
        if (all || newValues == "pid"              ) fWPid  = false;
        if (all || newValues == "track_step"       ) fWTS   = false;
        if (all || newValues == "kinetic_energy"   ) fWKE   = false;
        if (all || newValues == "energy_deposition") fWEDep = false;
        if (all || newValues == "position"         ) fWR    = false;
        if (all || newValues == "local_position"   ) fWLR   = false;
        if (all || newValues == "momentum"         ) fWP    = false;
        if (all || newValues == "time"             ) fWT    = false;
        if (all || newValues == "volume"           ) fWV    = false;
      }
      if (command == fAddOutputCmd) {
        G4bool all = (newValues == "all");
        if (all || newValues == "event"            ) fWEv   = true;
        if (all || newValues == "pid"              ) fWPid  = true;
        if (all || newValues == "track_step"       ) fWTS   = true;
        if (all || newValues == "kinetic_energy"   ) fWKE   = true;
        if (all || newValues == "energy_deposition") fWEDep = true;
        if (all || newValues == "position"         ) fWR    = true;
        if (all || newValues == "local_position"   ) fWLR   = true;
        if (all || newValues == "momentum"         ) fWP    = true;
        if (all || newValues == "time"             ) fWT    = true;
        if (all || newValues == "volume"           ) fWV    = true;
      }
    }

    void ResetVars() {
      fPID       .clear();
      fTrackID   .clear();
      fParentID  .clear();
      fStepNumber.clear();
      fKE        .clear();
      fEDep      .clear();
      fX         .clear();
      fY         .clear();
      fZ         .clear();
      fLX        .clear();
      fLY        .clear();
      fLZ        .clear();
      fPdX       .clear();
      fPdY       .clear();
      fPdZ       .clear();
      fT         .clear();
      fVolID     .clear();
      fIRep      .clear();
    }

    G4int GetVolID(G4StepPoint* stepPoint) {
      G4VPhysicalVolume* vpv = stepPoint->GetPhysicalVolume();
      G4int id = fVolIDMap[vpv];
      if (id == 0 && fPatternPairs.size() > 0) {
        std::string name = (vpv == NULL) ? "NULL" : vpv->GetName();

        for (auto& pp : fPatternPairs) {
          if (regex_match(name, pp.first)) {
            std::string replaced = regex_replace(name,pp.first,pp.second);
	          std::cout << "Setting ID for " << name << " to " << replaced << std::endl;
            int id_new = stoi(replaced);
            if (id_new == 0) {
              std::cout << "Volume " << name << ": Can't use ID = 0" << std::endl;
            } else {
              id = id_new;
            }

            break;
          }
        }

        fVolIDMap[vpv] = id;
      }

      return id;
    }

    void PushData(const G4Step* step, G4bool usePreStep=false, G4bool zeroEdep=false) {
      // g4simple output convention:
      // Each two rows form a pre-post step point pair.
      // In G4 one is always "in" the vol of the prestep point in G4
      // However in g4simple the volID along with the Edep of the step get
      // recorded along with the poststep point info.
      // This means that boundary crossings are recorded in g4simple output at
      // the first step AFTER hitting the boundary.
      G4StepPoint* stepPoint = step->GetPreStepPoint();
      fVolID.push_back(GetVolID(stepPoint));
      if (!usePreStep) stepPoint = stepPoint = step->GetPostStepPoint();
      fPID       .push_back(step->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
      fTrackID   .push_back(step->GetTrack()->GetTrackID());
      fParentID  .push_back(step->GetTrack()->GetParentID());
      fStepNumber.push_back(step->GetTrack()->GetCurrentStepNumber() - int(usePreStep));
      fKE.push_back(stepPoint->GetKineticEnergy());
      if (usePreStep || zeroEdep) fEDep.push_back(0);
      else fEDep.push_back(step->GetTotalEnergyDeposit());
      G4ThreeVector pos = stepPoint->GetPosition();
      fX.push_back(pos.x());
      fY.push_back(pos.y());
      fZ.push_back(pos.z());
      G4TouchableHandle vol = stepPoint->GetTouchableHandle();
      G4ThreeVector lPos = vol->GetHistory()->GetTopTransform().TransformPoint(pos);
      fLX.push_back(lPos.x());
      fLY.push_back(lPos.y());
      fLZ.push_back(lPos.z());
      G4ThreeVector momDir = stepPoint->GetMomentumDirection();
      fPdX .push_back(momDir.x());
      fPdY .push_back(momDir.y());
      fPdZ .push_back(momDir.z());
      fT   .push_back(stepPoint->GetGlobalTime());
      fIRep.push_back(vol->GetReplicaNumber());

      if (fOption == kStepWise) WriteRow();
    }

    void WriteRow() {
      G4VAnalysisManager* man = GetAnalysisManager();
      int iCol = 0;
      if (fWEv) man->FillNtupleIColumn(iCol++, fNEvents    );
      if (fWEv) man->FillNtupleIColumn(iCol++, fEventNumber);
      if (fOption == kStepWise) {
        size_t i = fPID.size() - 1;
        if (fWPid ) man->FillNtupleIColumn(iCol++, fPID       [i]);
        if (fWTS  ) man->FillNtupleIColumn(iCol++, fTrackID   [i]);
        if (fWTS  ) man->FillNtupleIColumn(iCol++, fParentID  [i]);
        if (fWTS  ) man->FillNtupleIColumn(iCol++, fStepNumber[i]);
        if (fWKE  ) man->FillNtupleDColumn(iCol++, fKE        [i]);
        if (fWEDep) man->FillNtupleDColumn(iCol++, fEDep      [i]);
        if (fWR   ) man->FillNtupleDColumn(iCol++, fX         [i]);
        if (fWR   ) man->FillNtupleDColumn(iCol++, fY         [i]);
        if (fWR   ) man->FillNtupleDColumn(iCol++, fZ         [i]);
        if (fWLR  ) man->FillNtupleDColumn(iCol++, fLX        [i]);
        if (fWLR  ) man->FillNtupleDColumn(iCol++, fLY        [i]);
        if (fWLR  ) man->FillNtupleDColumn(iCol++, fLZ        [i]);
        if (fWP   ) man->FillNtupleDColumn(iCol++, fPdX       [i]);
        if (fWP   ) man->FillNtupleDColumn(iCol++, fPdY       [i]);
        if (fWP   ) man->FillNtupleDColumn(iCol++, fPdZ       [i]);
        if (fWT   ) man->FillNtupleDColumn(iCol++, fT         [i]);
        if (fWV   ) man->FillNtupleIColumn(iCol++, fVolID     [i]);
        if (fWV   ) man->FillNtupleIColumn(iCol++, fIRep      [i]);
      }

      // for event-wise, manager copies data from vectors over
      // automatically in the next line
      man->AddNtupleRow();
    }

    void UserSteppingAction(const G4Step* step) {
      // This is the main function where we decide what to pull out and write
      // to an output file
      G4VAnalysisManager* man = GetAnalysisManager();

      // Open up a file if one is not open already
      if (!man->IsOpenFile()) {
        // need to create the ntuple before opening the file in order to avoid
        // writing error in csv, xml, and hdf5
        man->CreateNtuple("g4sntuple", "steps data");
        if (fWEv) man->CreateNtupleIColumn("nEvents");
        if (fWEv) man->CreateNtupleIColumn("event"  );
        if (fOption == kEventWise) {
          if (fWPid ) man->CreateNtupleIColumn("pid"     , fPID       );
          if (fWTS  ) man->CreateNtupleIColumn("trackID" , fTrackID   );
          if (fWTS  ) man->CreateNtupleIColumn("parentID", fParentID  );
          if (fWTS  ) man->CreateNtupleIColumn("step"    , fStepNumber);
          if (fWKE  ) man->CreateNtupleDColumn("KE"      , fKE        );
          if (fWEDep) man->CreateNtupleDColumn("Edep"    , fEDep      );
          if (fWR   ) man->CreateNtupleDColumn("x"       , fX         );
          if (fWR   ) man->CreateNtupleDColumn("y"       , fY         );
          if (fWR   ) man->CreateNtupleDColumn("z"       , fZ         );
          if (fWLR  ) man->CreateNtupleDColumn("lx"      , fLX        );
          if (fWLR  ) man->CreateNtupleDColumn("ly"      , fLY        );
          if (fWLR  ) man->CreateNtupleDColumn("lz"      , fLZ        );
          if (fWP   ) man->CreateNtupleDColumn("pdx"     , fPdX       );
          if (fWP   ) man->CreateNtupleDColumn("pdy"     , fPdY       );
          if (fWP   ) man->CreateNtupleDColumn("pdz"     , fPdZ       );
          if (fWT   ) man->CreateNtupleDColumn("t"       , fT         );
          if (fWV   ) man->CreateNtupleIColumn("volID"   , fVolID     );
          if (fWV   ) man->CreateNtupleIColumn("iRep"    , fIRep      );
        } else if (fOption == kStepWise) {
          if (fWPid ) man->CreateNtupleIColumn("pid"     );
          if (fWTS  ) man->CreateNtupleIColumn("trackID" );
          if (fWTS  ) man->CreateNtupleIColumn("parentID");
          if (fWTS  ) man->CreateNtupleIColumn("step"    );
          if (fWKE  ) man->CreateNtupleDColumn("KE"      );
          if (fWEDep) man->CreateNtupleDColumn("Edep"    );
          if (fWR   ) man->CreateNtupleDColumn("x"       );
          if (fWR   ) man->CreateNtupleDColumn("y"       );
          if (fWR   ) man->CreateNtupleDColumn("z"       );
          if (fWLR  ) man->CreateNtupleDColumn("lx"      );
          if (fWLR  ) man->CreateNtupleDColumn("ly"      );
          if (fWLR  ) man->CreateNtupleDColumn("lz"      );
          if (fWP   ) man->CreateNtupleDColumn("pdx"     );
          if (fWP   ) man->CreateNtupleDColumn("pdy"     );
          if (fWP   ) man->CreateNtupleDColumn("pdz"     );
          if (fWT   ) man->CreateNtupleDColumn("t"       );
          if (fWV   ) man->CreateNtupleIColumn("volID"   );
          if (fWV   ) man->CreateNtupleIColumn("iRep"    );
        } else {
          std::cout << "ERROR: Unknown output option " << fOption << std::endl;
          return;
        }

        man->FinishNtuple();
	      man->SetFileName("g4simple");
        std::cout << "Opening file " << man->GetFileName() << std::endl;
        man->OpenFile();

        ResetVars();
        fNEvents = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
        fVolIDMap.clear();
      }

      // Get the event number for recording
      fEventNumber = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
      static G4int lastEventID = fEventNumber;
      if (fEventNumber != lastEventID) {
        if (fOption == kEventWise && fPID.size() > 0) {
          WriteRow();
        }

        ResetVars();
        lastEventID = fEventNumber;
      }

      // If writing out all steps, just write and return.
      G4bool usePreStep = true;
      if (fRecordAllSteps) {
        if (step->GetTrack()->GetCurrentStepNumber() == 1) {
          PushData(step, usePreStep);
        }

        PushData(step);
        return;
      }

      // Not writing out all steps:
      // First record primary event info from pre-step of first step of first track
      if (step->GetTrack()->GetTrackID()           == 1 &&
          step->GetTrack()->GetCurrentStepNumber() == 1) {
        PushData(step, usePreStep);
      }

      // Below here: writing out only steps in sensitive volumes (volID != 0)
      G4int preID  = GetVolID(step->GetPreStepPoint ());
      G4int postID = GetVolID(step->GetPostStepPoint());

      // Record step data if in a sensitive volume and Edep > 0
      if (preID != 0 && step->GetTotalEnergyDeposit() > 0) {
        // Record pre-step data for the first step of a E-depositing track in a
        // sens vol. Note: if trackID == 1, we already recorded it
        if (step->GetTrack()->GetCurrentStepNumber() == 1 &&
            step->GetTrack()->GetTrackID()           >  1) {
          PushData(step, usePreStep);
        }
        // Record post-step data for all sens vol steps
        PushData(step);
        return; // don't need to re-write poststep below if we already wrote it out
      }

      // Record the step point when an energy-depositing particle first enters a
      // sensitive volume: it's the post-step-point of the step where the phys
      // vol pointer changes.
      // Have to do this last because we might have already written it out
      // during the last step of the previous volume if it was also sensitive.
      G4bool zeroEdep = true;
      if (step->GetPreStepPoint()->GetPhysicalVolume() != step->GetPostStepPoint()->GetPhysicalVolume()) {
        if (postID != 0 && step->GetTotalEnergyDeposit() > 0) {
          PushData(step, usePreStep = false, zeroEdep);
        }
      }
    }
};

class G4SimplePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    void GeneratePrimaries(G4Event* event) {
      fParticleGun.GeneratePrimaryVertex(event);
    }

  private:
    G4GeneralParticleSource fParticleGun;
};

class G4SimpleDetectorConstruction : public G4VUserDetectorConstruction {
  public:
    G4SimpleDetectorConstruction() {}

    virtual G4VPhysicalVolume* Construct() {
      G4NistManager* nistManager = G4NistManager::Instance();

      G4Material *mVacuum = new G4Material(
        "Vacuum",
        1.0,
        1.01 * g / mole,
        universe_mean_density,
        kStateGas,
        3.0e-18 * pascal,
        2.73 * kelvin);

      const G4int nEntries_vacuum = 2;
      G4double photonEnergy_vacuum   [nEntries_vacuum] = {1.0 * eV, 5.0 * eV};
      G4double refractiveIndex_vacuum[nEntries_vacuum] = {1.0, 1.0};
      G4MaterialPropertiesTable* table_vacuum = new G4MaterialPropertiesTable();
      table_vacuum->AddProperty("RINDEX", photonEnergy_vacuum, refractiveIndex_vacuum, nEntries_vacuum);

      mVacuum->SetMaterialPropertiesTable(table_vacuum);

      G4Material* mCe = nistManager->FindOrBuildMaterial("G4_Ce");
      G4Material* mO  = nistManager->FindOrBuildMaterial("G4_O" );
      G4Material* mSi = nistManager->FindOrBuildMaterial("G4_Si");
      G4Material* mLu = nistManager->FindOrBuildMaterial("G4_Lu");
      G4Material* mY  = nistManager->FindOrBuildMaterial("G4_Y" );

      G4Material* mLYSO = new G4Material("LYSO", 7.4 * g / cm3, 5);

      double pCe = 0.01;
      double pO  = (1 - pCe) * 5.0 / 8.0;
      double pSi = (1 - pCe) * 1.0 / 8.0;
      double pLu = (1 - pCe) * 2.0 / 8.0 * (1 - 0.05);
      double pY  = (1 - pCe) * 2.0 / 8.0 * (    0.05);

      mLYSO->AddMaterial(mCe, pCe);
      mLYSO->AddMaterial(mO , pO );
      mLYSO->AddMaterial(mSi, pSi);
      mLYSO->AddMaterial(mLu, pLu);
      mLYSO->AddMaterial(mY , pY );

      G4MaterialPropertiesTable* table_LYSO = new G4MaterialPropertiesTable();

      const G4int nEntries_LYSO = 117;

      G4double wavelength_LYSO[nEntries_LYSO] = {
        361.14468864468864 * nm,
        364.97252747252750 * nm,
        368.45238095238096 * nm,
        371.75824175824175 * nm,
        373.98851148851150 * nm,
        376.05006105006106 * nm,
        377.37345987345986 * nm,
        378.54395604395610 * nm,
        379.70390720390720 * nm,
        380.24521774521776 * nm,
        380.67053317053320 * nm,
        381.29304029304030 * nm,
        381.91941391941396 * nm,
        382.31379731379730 * nm,
        383.06776556776560 * nm,
        382.83577533577534 * nm,
        383.76373626373630 * nm,
        385.38461538461536 * nm,
        386.92307692307690 * nm,
        389.12271062271066 * nm,
        388.63553113553115 * nm,
        389.92804814233386 * nm,
        389.67948717948720 * nm,
        390.72344322344327 * nm,
        390.37545787545787 * nm,
        391.26475376475380 * nm,
        392.11538461538464 * nm,
        392.02838827838830 * nm,
        392.60256410256410 * nm,
        393.13034188034190 * nm,
        393.99450549450546 * nm,
        395.49581371009940 * nm,
        395.07326007326010 * nm,
        394.96886446886447 * nm,
        396.37820512820514 * nm,
        396.94851444851446 * nm,
        397.68315018315020 * nm,
        398.51296139757676 * nm,
        399.51007326007330 * nm,
        400.62520812520813 * nm,
        403.94688644688640 * nm,
        407.77472527472526 * nm,
        411.60256410256410 * nm,
        415.43040293040290 * nm,
        419.20470555085940 * nm,
        421.74053724053726 * nm,
        424.40345368916803 * nm,
        426.40775890775890 * nm,
        428.30586080586080 * nm,
        430.04578754578756 * nm,
        431.23888016745160 * nm,
        432.82967032967030 * nm,
        434.11721611721610 * nm,
        435.04412254412256 * nm,
        436.43606393606400 * nm,
        437.75946275946274 * nm,
        439.15140415140417 * nm,
        440.54334554334554 * nm,
        441.78238428238430 * nm,
        443.21123321123320 * nm,
        444.62953712953714 * nm,
        445.67032967032970 * nm,
        446.68581418581425 * nm,
        448.10939060939063 * nm,
        449.01098901098900 * nm,
        450.11294261294260 * nm,
        451.15689865689865 * nm,
        452.20085470085473 * nm,
        453.39560439560444 * nm,
        454.40476190476190 * nm,
        455.20015698587133 * nm,
        456.46103896103900 * nm,
        457.50183150183150 * nm,
        458.87057387057393 * nm,
        460.40081713158640 * nm,
        461.94444444444446 * nm,
        463.99369149369150 * nm,
        465.22130647130650 * nm,
        467.28021978021980 * nm,
        469.02014652014657 * nm,
        470.91824841824850 * nm,
        472.68981018981015 * nm,
        474.90912933220625 * nm,
        477.02380952380950 * nm,
        479.13848971541280 * nm,
        481.59401709401710 * nm,
        484.30830280830276 * nm,
        488.15934065934070 * nm,
        491.98717948717956 * nm,
        495.81501831501840 * nm,
        499.66962524654830 * nm,
        503.47069597069594 * nm,
        507.29853479853480 * nm,
        511.12637362637360 * nm,
        514.95421245421240 * nm,
        518.78205128205130 * nm,
        522.60989010989010 * nm,
        526.43772893772890 * nm,
        530.26556776556780 * nm,
        534.09340659340660 * nm,
        537.92124542124540 * nm,
        541.74908424908430 * nm,
        545.57692307692310 * nm,
        549.40476190476190 * nm,
        553.23260073260080 * nm,
        557.06043956043960 * nm,
        560.88827838827840 * nm,
        564.71611721611730 * nm,
        568.54395604395610 * nm,
        572.37179487179490 * nm,
        576.19963369963380 * nm,
        580.02747252747260 * nm,
        583.85531135531140 * nm,
        587.68315018315020 * nm,
        591.51098901098910 * nm,
        595.33882783882790 * nm,
        598.81868131868130 * nm};

      G4double emission_LYSO  [nEntries_LYSO] = {
        0.17844558542233813,
        0.80967482130273540,
        1.87553832902671050,
        2.83657649106154960,
        4.15232054766939000,
        5.37652270210411100,
        6.66289137219371200,
        7.87806078503754750,
        9.08499446290144900,
        10.4731143103236230,
        11.7998646487018560,
        13.3931339977851700,
        15.2973421926910330,
        17.1419342930970940,
        19.7467290102949140,
        18.2290513104466640,
        21.7534145441122180,
        24.1860465116279070,
        26.9767441860465030,
        31.7724252491694300,
        30.6256921373200440,
        34.7136529030216800,
        33.3084163898117400,
        37.1934820439803800,
        35.8333333333333300,
        38.8959640703826750,
        41.4837117017349600,
        40.0283776301218100,
        43.4186046511627950,
        44.9949243263196760,
        46.6168327796234700,
        50.7235801297263100,
        49.0891472868217000,
        47.9529346622369840,
        52.5170727205610750,
        53.8467454165128460,
        55.0975144579795640,
        56.4332566658147950,
        57.9877260981912000,
        59.3705325682069760,
        60.1117487163998660,
        60.2360817477096400,
        60.0782744387395400,
        59.6670190274841400,
        58.9136638555243050,
        57.7755629383536200,
        56.5436639772187850,
        55.2627604953186240,
        54.0442967884828300,
        52.8414544112218400,
        51.7718715393133950,
        50.6204626553463650,
        49.3311184939091800,
        48.2379442263163100,
        47.1189469445283300,
        45.9066998892580140,
        44.6573920265780800,
        43.3817829457364200,
        42.0882412161481900,
        40.6902916205241640,
        39.3404812242021500,
        37.9795127353266800,
        36.6319339575153600,
        35.1877579784556500,
        33.8475913621262400,
        32.6479635781961260,
        31.3621262458471800,
        30.0762889134982170,
        28.7372646733111740,
        27.4344776670358070,
        26.3573801613668760,
        25.1598207993556850,
        23.8715393133997800,
        22.5512181616832800,
        21.1572535991140700,
        19.7054263565891570,
        18.4847575981296900,
        17.2734403839055100,
        15.9706533776301280,
        14.7555370985603530,
        13.5681566495520090,
        12.4013389711064120,
        11.2275321577647280,
        10.0493724621631630,
        8.99799812590511300,
        7.91196013289037100,
        6.91602067183463000,
        6.04122621564484100,
        5.16133091714488050,
        4.43446088794927100,
        3.69324473975638060,
        3.21504077318032000,
        2.83247759991947130,
        2.55033725963959060,
        2.33992751434612960,
        2.11995368972114300,
        1.88563374609887550,
        1.68957011980269560,
        1.51263465216953820,
        1.35482734319944600,
        1.19702003422933960,
        1.03921272525924730,
        0.89575153528642200,
        0.73316218665057420,
        0.55622671901743100,
        0.40798348937885010,
        0.36494513238699255,
        0.31234269606362375,
        0.25017618040874370,
        0.22148394241418146,
        0.22148394241418146,
        0.20235578375113050,
        0.14018926809625043,
        0.11627906976745805,
        0.07802275244138457,
        0.06367663344408925,
        0.06367663344408925};

      G4double photonEnergy_LYSO    [nEntries_LYSO];
      G4double refractiveIndex_LYSO [nEntries_LYSO];
      G4double absorptionLength_LYSO[nEntries_LYSO];

      double maxEmission_LYSO = *std::max_element(std::begin(emission_LYSO), std::end(emission_LYSO));

      for (int i = 0; i < nEntries_LYSO; ++i){
        photonEnergy_LYSO    [i]  = 0.001240 * MeV * nm / wavelength_LYSO[i];
        emission_LYSO        [i] /= maxEmission_LYSO;
        refractiveIndex_LYSO [i]  = 1.83;
        absorptionLength_LYSO[i]  = 20.0;
      }

      table_LYSO->AddProperty("FASTCOMPONENT", photonEnergy_LYSO, emission_LYSO        , nEntries_LYSO);
      table_LYSO->AddProperty("SLOWCOMPONENT", photonEnergy_LYSO, emission_LYSO        , nEntries_LYSO);
      table_LYSO->AddProperty("RINDEX"       , photonEnergy_LYSO, refractiveIndex_LYSO , nEntries_LYSO);
      table_LYSO->AddProperty("ABSLENGTH"    , photonEnergy_LYSO, absorptionLength_LYSO, nEntries_LYSO);

      table_LYSO->AddConstProperty("FASTTIMECONSTANT"  ,    42.0 * ns );
      table_LYSO->AddConstProperty("SLOWTIMECONSTANT"  ,    42.0 * ns );
      table_LYSO->AddConstProperty("SCINTILLATIONYIELD", 32300.0 / MeV);
      table_LYSO->AddConstProperty("RESOLUTIONSCALE"   ,     1.0      );
      table_LYSO->AddConstProperty("YIELDRATIO"        ,     1.0      );

      mLYSO->SetMaterialPropertiesTable(table_LYSO);

      G4OpticalSurface* mTedlar = new G4OpticalSurface("Tedlar");

      mTedlar->SetType      (dielectric_dielectric);
      mTedlar->SetModel     (unified              );
      mTedlar->SetFinish    (groundbackpainted    );
      mTedlar->SetSigmaAlpha(0.07379              );

      const G4int nEntries_tedlar = 24;

      G4double photonEnergy_tedlar[nEntries_tedlar] = {
        1.38 * eV,
        1.55 * eV,
        1.77 * eV,
        1.82 * eV,
        1.88 * eV,
        1.94 * eV,
        2.00 * eV,
        2.07 * eV,
        2.14 * eV,
        2.21 * eV,
        2.30 * eV,
        2.38 * eV,
        2.48 * eV,
        2.58 * eV,
        2.70 * eV,
        2.82 * eV,
        2.95 * eV,
        3.10 * eV,
        3.31 * eV,
        3.54 * eV,
        3.81 * eV,
        4.13 * eV,
        4.51 * eV,
        4.96 * eV};

      G4double reflectivity_tedlar[nEntries_tedlar] = {
        0.0669,
        0.0741,
        0.0708,
        0.0708,
        0.0708,
        0.0713,
        0.0718,
        0.0728,
        0.0735,
        0.0741,
        0.0749,
        0.0760,
        0.0776,
        0.0788,
        0.0805,
        0.0821,
        0.0835,
        0.0831,
        0.0679,
        0.0601,
        0.0605,
        0.0631,
        0.0635,
        0.0637};

      G4double refractiveIndex_tedlar[nEntries_tedlar];
      G4double specularSpike_tedlar  [nEntries_tedlar];
      G4double specularLobe_tedlar   [nEntries_tedlar];
      G4double backscatter_tedlar    [nEntries_tedlar];

      std::fill_n(refractiveIndex_tedlar, nEntries_tedlar, 1.0);
      std::fill_n(specularSpike_tedlar  , nEntries_tedlar, 0.0);
      std::fill_n(specularLobe_tedlar   , nEntries_tedlar, 1.0);
      std::fill_n(backscatter_tedlar    , nEntries_tedlar, 0.0);

      G4MaterialPropertiesTable* table_tedlar = new G4MaterialPropertiesTable();

      table_tedlar->AddProperty("RINDEX"               , photonEnergy_tedlar, refractiveIndex_tedlar, nEntries_tedlar);
      table_tedlar->AddProperty("SPECULARSPIKECONSTANT", photonEnergy_tedlar, specularSpike_tedlar  , nEntries_tedlar);
      table_tedlar->AddProperty("SPECULARLOBECONSTANT" , photonEnergy_tedlar, specularLobe_tedlar   , nEntries_tedlar);
      table_tedlar->AddProperty("BACKSCATTERCONSTANT"  , photonEnergy_tedlar, backscatter_tedlar    , nEntries_tedlar);
      table_tedlar->AddProperty("REFLECTIVITY"         , photonEnergy_tedlar, reflectivity_tedlar   , nEntries_tedlar);

      mTedlar->SetMaterialPropertiesTable(table_tedlar);

      G4OpticalSurface* mTedlarReverse = new G4OpticalSurface("TedlarReverse");

      mTedlarReverse->SetType  (dielectric_dielectric);
      mTedlarReverse->SetModel (unified              );
      mTedlarReverse->SetFinish(groundfrontpainted   );

      G4MaterialPropertiesTable* table_tedlarReverse = new G4MaterialPropertiesTable();

      table_tedlarReverse->AddProperty("REFLECTIVITY", photonEnergy_tedlar, reflectivity_tedlar, nEntries_tedlar);

      mTedlarReverse->SetMaterialPropertiesTable(table_tedlarReverse);

      int nXtalRows = 1;
      int nXtalCols = 1;

      double wrappingGap =  0.00 * cm; // cm
      double xtalWidth   =  2.50 * cm; // cm
      double xtalHeight  =  2.50 * cm; // cm
      double xtalDepth   = 14.00 * cm; // cm

      G4ThreeVector xhat(1.0, 0.0, 0.0);
      G4ThreeVector yhat(0.0, 1.0, 0.0);
      G4ThreeVector zhat(0.0, 0.0, 1.0);

      G4Box* world_S = new G4Box("world_S", 300 * cm, 300 * cm, 300 * cm);
      G4LogicalVolume* world_L = new G4LogicalVolume(world_S, mVacuum, "world_L");
      world_P = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), world_L, "world_P", 0, false, 0, true);

      G4ThreeVector xtalOrigin(0.0, 0.0, 0.0);

      for (int irow = 0; irow < nXtalRows; ++irow) {
        G4ThreeVector offset =
          - double(nXtalCols - 1) / 2.0 * (xtalWidth + wrappingGap) * xhat
          - (double(nXtalRows - 1) / 2.0 - double(irow)) * (xtalHeight + wrappingGap) * yhat;
        G4ThreeVector xtalPos = xtalOrigin + offset;

        for (int icol = 0; icol < nXtalCols; ++icol) {
          G4Box* xtal_S = new G4Box(
            "xtal_S",
            xtalWidth  / 2.0,
            xtalHeight / 2.0,
            xtalDepth  / 2.0);

          G4LogicalVolume  * xtal_L = new G4LogicalVolume(xtal_S, mLYSO, "xtal_L");
          G4VPhysicalVolume* xtal_P = new G4PVPlacement  (
            0, xtalPos, xtal_L, "xtal_P", world_L, false, irow * nXtalCols + icol, true);

          new G4LogicalBorderSurface("LYSOSurface"       , xtal_P , world_P, mTedlar       );
          new G4LogicalBorderSurface("LYSOSurfaceReverse", world_P, xtal_P , mTedlarReverse);

          xtalPos += (xtalWidth + wrappingGap) * xhat;
        }
      }

      return world_P;
    }

  private:
    G4VPhysicalVolume* world_P;
};

class G4SimpleRunManager : public G4RunManager, public G4UImessenger {
  private:
    G4UIdirectory     * fDirectory;
    G4UIcmdWithAString* fPhysListCmd;
    G4UIcmdWithABool  * fRandomSeedCmd;
    G4UIcmdWithAString* fListVolsCmd;

  public:
    G4SimpleRunManager() {
      fDirectory = new G4UIdirectory("/g4simple/");
      fDirectory->SetGuidance("Parameters for g4simple MC");

      fPhysListCmd = new G4UIcmdWithAString("/g4simple/setReferencePhysList", this);
      fPhysListCmd->SetGuidance("Set reference physics list to be used");

      fRandomSeedCmd = new G4UIcmdWithABool("/g4simple/setRandomSeed", this);
      fRandomSeedCmd->SetParameterName("useURandom", true);
      fRandomSeedCmd->SetDefaultValue(false);
      fRandomSeedCmd->SetGuidance("Seed random number generator with a read from /dev/random");
      fRandomSeedCmd->SetGuidance("Set useURandom to true to read instead from /dev/urandom (faster but less random)");

      fListVolsCmd = new G4UIcmdWithAString("/g4simple/listPhysVols", this);
      fListVolsCmd->SetParameterName("pattern", true);
      fListVolsCmd->SetGuidance("List name of all instantiated physical volumes");
      fListVolsCmd->SetGuidance("Optionally supply a regex pattern to only list matching volume names");
      fListVolsCmd->AvailableForStates(G4State_Idle, G4State_GeomClosed, G4State_EventProc);
    }

    ~G4SimpleRunManager() {
      delete fDirectory;
      delete fPhysListCmd;
      delete fRandomSeedCmd;
      delete fListVolsCmd;
    }

    void SetNewValue(G4UIcommand* command, G4String newValues) {
      if (command == fPhysListCmd) {
        SetUserInitialization((new G4PhysListFactory)->GetReferencePhysList(newValues));
        SetUserAction(new G4SimplePrimaryGeneratorAction); // must come after phys list
        SetUserAction(new G4SimpleSteppingAction); // must come after phys list
      } else if (command == fRandomSeedCmd) {
        bool useURandom = fRandomSeedCmd->GetNewBoolValue(newValues);
        std::string path = useURandom ? "/dev/urandom" : "/dev/random";

        std::ifstream devrandom(path.c_str());
        if (!devrandom.good()) {
          std::cout << "setRandomSeed: couldn't open " << path << ". Your seed is not set." << std::endl;
          return;
        }

        long seed[2];
        devrandom.read((char*)(seed    ), sizeof(long));
        devrandom.read((char*)(seed + 1), sizeof(long));

        CLHEP::HepRandom::setTheSeeds(seed, 2);
        std::cout << "CLHEP::HepRandom seeds set to: " << seed[0] << ' ' << seed[1] << std::endl;
        devrandom.close();
      } else if (command == fListVolsCmd) {
        std::regex pattern(newValues);
        bool doMatching = (newValues != "");
        G4PhysicalVolumeStore* volumeStore = G4PhysicalVolumeStore::GetInstance();
        std::cout << "Physical volumes";
        if (doMatching) std::cout << " matching pattern " << newValues;
        std::cout << ":" << std::endl;
        for (size_t i = 0; i < volumeStore->size(); ++i) {
          std::string name = volumeStore->at(i)->GetName();
	        int iRep = volumeStore->at(i)->GetCopyNo();
          if (!doMatching || regex_match(name, pattern)) {
            std::cout << name << ' ' << iRep << std::endl;
          }
        }
      }
    }
};

int main(int argc, char** argv) {
  if (argc > 2) {
    std::cout << "usage: " << argv[0] << " [macro]" << std::endl;
    return 1;
  }

  G4VisManager      * visManager = new G4VisExecutive;
  G4SimpleRunManager* runManager = new G4SimpleRunManager;

  visManager->Initialize();
  runManager->SetUserInitialization(new G4SimpleDetectorConstruction());

  if (argc == 1) {
    (new G4UIterminal(new G4UItcsh))->SessionStart();
  } else {
    G4UImanager::GetUIpointer()->ApplyCommand(G4String("/control/execute ") + argv[1]);
  }

  delete visManager;
  delete runManager;

  return 0;
}
