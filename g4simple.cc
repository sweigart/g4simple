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
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysListFactory.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
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

    void SetNewValue(G4UIcommand *command, G4String newValues) {
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

    void UserSteppingAction(const G4Step *step) {
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

        // look for filename set by macro command: /analysis/setFileName [name]
	      if (man->GetFileName() == "") man->SetFileName("g4simpleout");
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
    G4SimpleDetectorConstruction(G4VPhysicalVolume* world = 0) {
      fWorld = world;
    }

    virtual G4VPhysicalVolume* Construct() {
      return fWorld;
    }

  private:
    G4VPhysicalVolume* fWorld;
};

class G4SimpleRunManager : public G4RunManager, public G4UImessenger {
  private:
    G4UIdirectory     * fDirectory;
    G4UIcmdWithAString* fPhysListCmd;
    G4UIcommand       * fDetectorCmd;
    G4UIcmdWithABool  * fRandomSeedCmd;
    G4UIcmdWithAString* fListVolsCmd;

  public:
    G4SimpleRunManager() {
      fDirectory = new G4UIdirectory("/g4simple/");
      fDirectory->SetGuidance("Parameters for g4simple MC");

      fPhysListCmd = new G4UIcmdWithAString("/g4simple/setReferencePhysList", this);
      fPhysListCmd->SetGuidance("Set reference physics list to be used");

      fDetectorCmd = new G4UIcmdWithAString("/g4simple/setDetectorGDML", this);
      fDetectorCmd->SetGuidance("Provide GDML filename specifying the detector construction");

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
      delete fDetectorCmd;
      delete fRandomSeedCmd;
      delete fListVolsCmd;
    }

    void SetNewValue(G4UIcommand* command, G4String newValues) {
      if (command == fPhysListCmd) {
        SetUserInitialization((new G4PhysListFactory)->GetReferencePhysList(newValues));
        SetUserAction(new G4SimplePrimaryGeneratorAction); // must come after phys list
        SetUserAction(new G4SimpleSteppingAction); // must come after phys list
      } else if (command == fDetectorCmd) {
        G4GDMLParser parser;
        parser.Read(newValues, true);
        SetUserInitialization(new G4SimpleDetectorConstruction(parser.GetWorldVolume())); // <====== !!!!!!!!!
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

  G4SimpleRunManager* runManager = new G4SimpleRunManager;
  G4VisManager      * visManager = new G4VisExecutive;
  visManager->Initialize();

  if (argc == 1) {
    (new G4UIterminal(new G4UItcsh))->SessionStart();
  } else {
    G4UImanager::GetUIpointer()->ApplyCommand(G4String("/control/execute ") + argv[1]);
  }

  delete visManager;
  delete runManager;

  return 0;
}
