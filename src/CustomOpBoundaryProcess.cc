
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

#include "CustomOpBoundaryProcess.hh"
#include "G4GeometryTolerance.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ParallelWorldProcess.hh"

#include "G4SystemOfUnits.hh"

CustomOpBoundaryProcess::CustomOpBoundaryProcess(const G4String& processName,
                                               G4ProcessType type)
             : G4VDiscreteProcess(processName, type)
{
	SetVerboseLevel(0);////my
        if ( verboseLevel > 0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        SetProcessSubType(fOpBoundary);
	theStatus = Undefined;
	theModel = glisur;
	theFinish = polished;
        theReflectivity =  1.;
        theEfficiency   =  0.;
        theTransmittance = 0.;

        prob_sl = 0.;
        prob_ss = 0.;
        prob_bs = 0.;

        PropertyPointer  = NULL;
        PropertyPointer1 = NULL;
        PropertyPointer2 = NULL;

        Material1 = NULL;
        Material2 = NULL;

        OpticalSurface = NULL;

        kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

        iTE = iTM = 0;
        thePhotonMomentum = 0.;
        Rindex1 = Rindex2 = cost1 = cost2 = sint1 = sint2 = 0.;

        idx = idy = 0;
        DichroicVector = NULL;
		to_decider_ = NULL;
}

CustomOpBoundaryProcess::~CustomOpBoundaryProcess(){}


// PostStepDoIt
// ------------
//

G4VParticleChange*
CustomOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
        theStatus = Undefined;

        aParticleChange.Initialize(aTrack);
        aParticleChange.ProposeVelocity(aTrack.GetVelocity());

        // Get hyperStep from  G4ParallelWorldProcess
        //  NOTE: PostSetpDoIt of this process should be
        //        invoked after G4ParallelWorldProcess!

        const G4Step* pStep = &aStep;
		to_decider_ = pStep;
        const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();
        
        if (hStep) pStep = hStep;

        G4bool isOnBoundary =
                (pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);

        if (isOnBoundary) {
           Material1 = pStep->GetPreStepPoint()->GetMaterial();
           Material2 = pStep->GetPostStepPoint()->GetMaterial();
        } else {
           theStatus = NotAtBoundary;
           if ( verboseLevel > 0) BoundaryProcessVerbose();
           return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
        }

        G4VPhysicalVolume* thePrePV  =
                               pStep->GetPreStepPoint() ->GetPhysicalVolume();
        G4VPhysicalVolume* thePostPV =
                               pStep->GetPostStepPoint()->GetPhysicalVolume();

        if ( verboseLevel > 0 ) {
           G4cout << " Photon at Boundary! " << G4endl;
           if (thePrePV)  G4cout << " thePrePV:  " << thePrePV->GetName()  << G4endl;
           if (thePostPV) G4cout << " thePostPV: " << thePostPV->GetName() << G4endl;
        }

	if (aTrack.GetStepLength()<=kCarTolerance/2){
	        theStatus = StepTooSmall;
                if ( verboseLevel > 0) BoundaryProcessVerbose();
	        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

	thePhotonMomentum = aParticle->GetTotalMomentum();
        OldMomentum       = aParticle->GetMomentumDirection();
	OldPolarization   = aParticle->GetPolarization();

        if ( verboseLevel > 0 ) {
           G4cout << " Old Momentum Direction: " << OldMomentum     << G4endl;
           G4cout << " Old Polarization:       " << OldPolarization << G4endl;
        }

        G4ThreeVector theGlobalPoint = pStep->GetPostStepPoint()->GetPosition();

        G4bool valid;
        //  Use the new method for Exit Normal in global coordinates,
        //    which provides the normal more reliably.

        // ID of Navigator which limits step

        G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
        std::vector<G4Navigator*>::iterator iNav =
                G4TransportationManager::GetTransportationManager()->
                                         GetActiveNavigatorsIterator();
        theGlobalNormal =
                   (iNav[hNavId])->GetGlobalExitNormal(theGlobalPoint,&valid);
        if (valid) {
          theGlobalNormal = -theGlobalNormal;
        }
        else 
        {
		  G4ExceptionDescription ed;
          ed << " G4OpBoundaryProcess/PostStepDoIt(): "
                 << " The Navigator reports that it returned an invalid normal"
                 << G4endl;
          G4Exception("G4OpBoundaryProcess::PostStepDoIt", "OpBoun01",
                      EventMustBeAborted,ed,
                      "Invalid Surface Normal - Geometry must return valid surface normal");
        }

        if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4OPTICAL_DEBUG
           G4ExceptionDescription ed;
           ed << " G4OpBoundaryProcess/PostStepDoIt(): "
              << " theGlobalNormal points in a wrong direction. "
              << G4endl;
           ed << "    The momentum of the photon arriving at interface (oldMomentum)"
              << " must exit the volume cross in the step. " << G4endl;
           ed << "  So it MUST have dot < 0 with the normal that Exits the new volume (globalNormal)." << G4endl;
           ed << "  >> The dot product of oldMomentum and global Normal is " << OldMomentum*theGlobalNormal << G4endl;
           ed << "     Old Momentum  (during step)     = " << OldMomentum << G4endl;
           ed << "     Global Normal (Exiting New Vol) = " << theGlobalNormal << G4endl;
           ed << G4endl;
           G4Exception("G4OpBoundaryProcess::PostStepDoIt", "OpBoun02",
                       EventMustBeAborted,  // Or JustWarning to see if it happens repeatedbly on one ray
                       ed,
                      "Invalid Surface Normal - Geometry must return valid surface normal pointing in the right direction");
#else
           theGlobalNormal = -theGlobalNormal;
#endif
        }

	G4MaterialPropertiesTable* aMaterialPropertiesTable;
        G4MaterialPropertyVector* Rindex;

	aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();
        if (aMaterialPropertiesTable) {
		Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	}
	else {
	        theStatus = NoRINDEX;
                if ( verboseLevel > 0) BoundaryProcessVerbose();
                aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        if (Rindex) {
		Rindex1 = Rindex->Value(thePhotonMomentum);
	}
	else {
	        theStatus = NoRINDEX;
                if ( verboseLevel > 0) BoundaryProcessVerbose();
                aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
		aParticleChange.ProposeTrackStatus(fStopAndKill);
		return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        theReflectivity =  1.;
        theEfficiency   =  0.;
        theTransmittance = 0.;

        theModel = glisur;
        theFinish = polished;

        G4SurfaceType type = dielectric_dielectric;

        Rindex = NULL;
        OpticalSurface = NULL;

        G4LogicalSurface* Surface = NULL;

        Surface = G4LogicalBorderSurface::GetSurface(thePrePV, thePostPV);

        if (Surface == NULL){
          G4bool enteredDaughter= (thePostPV->GetMotherLogical() ==
                                   thePrePV ->GetLogicalVolume());
	  if(enteredDaughter){
	    Surface = 
              G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume());
	    if(Surface == NULL)
	      Surface =
                G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume());
	  }
	  else {
	    Surface =
              G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume());
	    if(Surface == NULL)
	      Surface =
                G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume());
	  }
	}

	if (Surface) OpticalSurface = 
           dynamic_cast <G4OpticalSurface*> (Surface->GetSurfaceProperty());

	if (OpticalSurface) {

           type      = OpticalSurface->GetType();
	   theModel  = OpticalSurface->GetModel();
	   theFinish = OpticalSurface->GetFinish();

	   aMaterialPropertiesTable = OpticalSurface->
					GetMaterialPropertiesTable();

           if (aMaterialPropertiesTable) {

              if (theFinish == polishedbackpainted ||
                  theFinish == groundbackpainted ) {
                  Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	          if (Rindex) {
                     Rindex2 = Rindex->Value(thePhotonMomentum);
                  }
                  else {
		     theStatus = NoRINDEX;
                     if ( verboseLevel > 0) BoundaryProcessVerbose();
                     aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                     aParticleChange.ProposeTrackStatus(fStopAndKill);
                     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
                  }
              }

              PropertyPointer =
                      aMaterialPropertiesTable->GetProperty("REFLECTIVITY");
              PropertyPointer1 =
                      aMaterialPropertiesTable->GetProperty("REALRINDEX");
              PropertyPointer2 =
                      aMaterialPropertiesTable->GetProperty("IMAGINARYRINDEX");

              iTE = 1;
              iTM = 1;

              if (PropertyPointer) {

                 theReflectivity =
                          PropertyPointer->Value(thePhotonMomentum);

              } else if (PropertyPointer1 && PropertyPointer2) {

                 CalculateReflectivity();

              }

              PropertyPointer =
              aMaterialPropertiesTable->GetProperty("EFFICIENCY");
              if (PropertyPointer) {
                      theEfficiency =
                      PropertyPointer->Value(thePhotonMomentum);
              }

              PropertyPointer =
              aMaterialPropertiesTable->GetProperty("TRANSMITTANCE");
              if (PropertyPointer) {
                      theTransmittance =
                      PropertyPointer->Value(thePhotonMomentum);
              }

	      if ( theModel == unified ) {
	        PropertyPointer =
		aMaterialPropertiesTable->GetProperty("SPECULARLOBECONSTANT");
	        if (PropertyPointer) {
                         prob_sl =
			 PropertyPointer->Value(thePhotonMomentum);
                } else {
                         prob_sl = 0.0;
                }

	        PropertyPointer =
		aMaterialPropertiesTable->GetProperty("SPECULARSPIKECONSTANT");
	        if (PropertyPointer) {
                         prob_ss =
			 PropertyPointer->Value(thePhotonMomentum);
                } else {
                         prob_ss = 0.0;
                }

	        PropertyPointer =
		aMaterialPropertiesTable->GetProperty("BACKSCATTERCONSTANT");
	        if (PropertyPointer) {
                         prob_bs =
			 PropertyPointer->Value(thePhotonMomentum);
                } else {
                         prob_bs = 0.0;
                }
	      }
	   }
           else if (theFinish == polishedbackpainted ||
                    theFinish == groundbackpainted ) {
                      aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                      aParticleChange.ProposeTrackStatus(fStopAndKill);
                      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
           }
        }

        if (type == dielectric_dielectric ) {
           if (theFinish == polished || theFinish == ground ) {

	      if (Material1 == Material2){
		 theStatus = SameMaterial;
                 if ( verboseLevel > 0) BoundaryProcessVerbose();
		 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	      }
              aMaterialPropertiesTable =
                     Material2->GetMaterialPropertiesTable();
              if (aMaterialPropertiesTable)
                 Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
              if (Rindex) {
                 Rindex2 = Rindex->Value(thePhotonMomentum);
              }
              else {
		 theStatus = NoRINDEX;
                 if ( verboseLevel > 0) BoundaryProcessVerbose();
                 aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                 aParticleChange.ProposeTrackStatus(fStopAndKill);
                 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	      }
           }
        }

	if (type == dielectric_metal) {

	  DielectricMetal(&aStep);

          // Uncomment the following lines if you wish to have 
          //         Transmission instead of Absorption
          // if (theStatus == Absorption) {
          //    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
          // }

	}
        else if (type == dielectric_LUT) {

          DielectricLUT();

        }
        else if (type == dielectric_dichroic) {

          DielectricDichroic();

        }
	else if (type == dielectric_dielectric) {

          if ( theFinish == polishedbackpainted ||
               theFinish == groundbackpainted ) {
             DielectricDielectric(&aStep);
          }
	  else {
             if ( !G4BooleanRand(theReflectivity) ) {
                DoAbsorption();
             }
             else {
                if ( theFinish == polishedfrontpainted ) {
                   DoReflection();
                }
                else if ( theFinish == groundfrontpainted ) {
                   theStatus = LambertianReflection;
                   DoReflection();
                }
                else {
					DielectricDielectric(&aStep);
                }
             }
          }
        }
	else {

	  G4cerr << " Error: G4BoundaryProcess: illegal boundary type " << G4endl;
	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

	}

        NewMomentum = NewMomentum.unit();
        NewPolarization = NewPolarization.unit();

        if ( verboseLevel > 0) {
	   G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;
	   G4cout << " New Polarization:       " << NewPolarization << G4endl;
           BoundaryProcessVerbose();
        }

	aParticleChange.ProposeMomentumDirection(NewMomentum);
	aParticleChange.ProposePolarization(NewPolarization);

        if ( theStatus == FresnelRefraction ) {
           G4MaterialPropertyVector* groupvel =
           Material2->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
           G4double finalVelocity = groupvel->Value(thePhotonMomentum);
           aParticleChange.ProposeVelocity(finalVelocity);
        }

        if ( theStatus == Detection ) InvokeSD(pStep);

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void CustomOpBoundaryProcess::BoundaryProcessVerbose() const
{
	if (verboseLevel > 0)
	{
		if (theStatus == Undefined)
			G4cout << " *** Undefined *** " << G4endl;
		if (theStatus == FresnelRefraction)
			G4cout << " *** FresnelRefraction *** " << G4endl;
		if (theStatus == FresnelReflection)
			G4cout << " *** FresnelReflection *** " << G4endl;
		if (theStatus == TotalInternalReflection)
			G4cout << " *** TotalInternalReflection *** " << G4endl;
		if (theStatus == LambertianReflection)
			G4cout << " *** LambertianReflection *** " << G4endl;
		if (theStatus == LobeReflection)
			G4cout << " *** LobeReflection *** " << G4endl;
		if (theStatus == SpikeReflection)
			G4cout << " *** SpikeReflection *** " << G4endl;
		if (theStatus == BackScattering)
			G4cout << " *** BackScattering *** " << G4endl;
		if (theStatus == PolishedLumirrorAirReflection)
			G4cout << " *** PolishedLumirrorAirReflection *** " << G4endl;
		if (theStatus == PolishedLumirrorGlueReflection)
			G4cout << " *** PolishedLumirrorGlueReflection *** " << G4endl;
		if (theStatus == PolishedAirReflection)
			G4cout << " *** PolishedAirReflection *** " << G4endl;
		if (theStatus == PolishedTeflonAirReflection)
			G4cout << " *** PolishedTeflonAirReflection *** " << G4endl;
		if (theStatus == PolishedTiOAirReflection)
			G4cout << " *** PolishedTiOAirReflection *** " << G4endl;
		if (theStatus == PolishedTyvekAirReflection)
			G4cout << " *** PolishedTyvekAirReflection *** " << G4endl;
		if (theStatus == PolishedVM2000AirReflection)
			G4cout << " *** PolishedVM2000AirReflection *** " << G4endl;
		if (theStatus == PolishedVM2000GlueReflection)
			G4cout << " *** PolishedVM2000GlueReflection *** " << G4endl;
		if (theStatus == EtchedLumirrorAirReflection)
			G4cout << " *** EtchedLumirrorAirReflection *** " << G4endl;
		if (theStatus == EtchedLumirrorGlueReflection)
			G4cout << " *** EtchedLumirrorGlueReflection *** " << G4endl;
		if (theStatus == EtchedAirReflection)
			G4cout << " *** EtchedAirReflection *** " << G4endl;
		if (theStatus == EtchedTeflonAirReflection)
			G4cout << " *** EtchedTeflonAirReflection *** " << G4endl;
		if (theStatus == EtchedTiOAirReflection)
			G4cout << " *** EtchedTiOAirReflection *** " << G4endl;
		if (theStatus == EtchedTyvekAirReflection)
			G4cout << " *** EtchedTyvekAirReflection *** " << G4endl;
		if (theStatus == EtchedVM2000AirReflection)
			G4cout << " *** EtchedVM2000AirReflection *** " << G4endl;
		if (theStatus == EtchedVM2000GlueReflection)
			G4cout << " *** EtchedVM2000GlueReflection *** " << G4endl;
		if (theStatus == GroundLumirrorAirReflection)
			G4cout << " *** GroundLumirrorAirReflection *** " << G4endl;
		if (theStatus == GroundLumirrorGlueReflection)
			G4cout << " *** GroundLumirrorGlueReflection *** " << G4endl;
		if (theStatus == GroundAirReflection)
			G4cout << " *** GroundAirReflection *** " << G4endl;
		if (theStatus == GroundTeflonAirReflection)
			G4cout << " *** GroundTeflonAirReflection *** " << G4endl;
		if (theStatus == GroundTiOAirReflection)
			G4cout << " *** GroundTiOAirReflection *** " << G4endl;
		if (theStatus == GroundTyvekAirReflection)
			G4cout << " *** GroundTyvekAirReflection *** " << G4endl;
		if (theStatus == GroundVM2000AirReflection)
			G4cout << " *** GroundVM2000AirReflection *** " << G4endl;
		if (theStatus == GroundVM2000GlueReflection)
			G4cout << " *** GroundVM2000GlueReflection *** " << G4endl;
		if (theStatus == Absorption)
			G4cout << " *** Absorption *** " << G4endl;
		if (theStatus == Detection)
			G4cout << " *** Detection *** " << G4endl;
		if (theStatus == NotAtBoundary)
			G4cout << " *** NotAtBoundary *** " << G4endl;
		if (theStatus == SameMaterial)
			G4cout << " *** SameMaterial *** " << G4endl;
		if (theStatus == StepTooSmall)
			G4cout << " *** StepTooSmall *** " << G4endl;
		if (theStatus == NoRINDEX)
			G4cout << " *** NoRINDEX *** " << G4endl;
		if (theStatus == Dichroic)
			G4cout << " *** Dichroic Transmission *** " << G4endl;
	}
}

G4ThreeVector
CustomOpBoundaryProcess::GetFacetNormal(const G4ThreeVector& Momentum,
			            const G4ThreeVector&  Normal ) const
{
        G4ThreeVector FacetNormal;

	if (theModel == unified || theModel == LUT) {

	/* This function code alpha to a random value taken from the
           distribution p(alpha) = g(alpha; 0, sigma_alpha)*std::sin(alpha),
           for alpha > 0 and alpha < 90, where g(alpha; 0, sigma_alpha)
           is a gaussian distribution with mean 0 and standard deviation
           sigma_alpha.  */

	   G4double alpha;

	   G4double sigma_alpha = 0.0;
	   if (OpticalSurface) sigma_alpha = OpticalSurface->GetSigmaAlpha();

           if (sigma_alpha == 0.0) return FacetNormal = Normal;

	   G4double f_max = std::min(1.0,4.*sigma_alpha);

	   do {
	      do {
	         alpha = G4RandGauss::shoot(0.0,sigma_alpha);
	      } while (G4UniformRand()*f_max > std::sin(alpha) || alpha >= halfpi );

	      G4double phi = G4UniformRand()*twopi;

	      G4double SinAlpha = std::sin(alpha);
	      G4double CosAlpha = std::cos(alpha);
              G4double SinPhi = std::sin(phi);
              G4double CosPhi = std::cos(phi);

              G4double unit_x = SinAlpha * CosPhi;
              G4double unit_y = SinAlpha * SinPhi;
              G4double unit_z = CosAlpha;

	      FacetNormal.setX(unit_x);
	      FacetNormal.setY(unit_y);
	      FacetNormal.setZ(unit_z);

	      G4ThreeVector tmpNormal = Normal;

	      FacetNormal.rotateUz(tmpNormal);
	   } while (Momentum * FacetNormal >= 0.0);
	}
	else {

	   G4double  polish = 1.0;
	   if (OpticalSurface) polish = OpticalSurface->GetPolish();

           if (polish < 1.0) {
              do {
                 G4ThreeVector smear;
                 do {
                    smear.setX(2.*G4UniformRand()-1.0);
                    smear.setY(2.*G4UniformRand()-1.0);
                    smear.setZ(2.*G4UniformRand()-1.0);
                 } while (smear.mag()>1.0);
                 smear = (1.-polish) * smear;
                 FacetNormal = Normal + smear;
              } while (Momentum * FacetNormal >= 0.0);
              FacetNormal = FacetNormal.unit();
	   }
           else {
              FacetNormal = Normal;
           }
	}
	return FacetNormal;
}

void CustomOpBoundaryProcess::DielectricMetal(const G4Step* aStep)
{
	G4double net_refl_prob = 1;
    G4int n = 0;
#ifdef DEBUG_MC_NODES
	G4cout << "COBP: dielectric_metal()" << G4endl;
#endif
	do {

           n++;

           //MY if( !G4BooleanRand(theReflectivity) && n == 1 ) {

           //  // Comment out DoAbsorption and uncomment theStatus = Absorption;
           //  // if you wish to have Transmission instead of Absorption

           //  DoAbsorption();
           //  // theStatus = Absorption;
           //  break;

           //}
		   if (n == 1)
		   {
			   net_refl_prob *= theReflectivity;
		   }
           /*MY else */{

             if (PropertyPointer1 && PropertyPointer2) {
                if ( n > 1 ) {
                   CalculateReflectivity();
				   net_refl_prob *= theReflectivity;
                   /*MY if ( !G4BooleanRand(theReflectivity) ) {
                      DoAbsorption();
                      break;
                   }*/
                }
             }

             if ( theModel == glisur || theFinish == polished ) {
				
                DoReflection();

             } else { //not realized in my models

                if ( n == 1 ) ChooseReflection();
                                                                                
                if ( theStatus == LambertianReflection ) {
                   DoReflection();
                }
                else if ( theStatus == BackScattering ) {
                   NewMomentum = -OldMomentum;
                   NewPolarization = -OldPolarization;
                }
                else {

                   if(theStatus==LobeReflection){
                     if ( PropertyPointer1 && PropertyPointer2 ){
                     } else {
                        theFacetNormal =
                            GetFacetNormal(OldMomentum,theGlobalNormal);
                     }
                   }

                   G4double PdotN = OldMomentum * theFacetNormal;
                   NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
                   G4double EdotN = OldPolarization * theFacetNormal;

                   G4ThreeVector A_trans, A_paral;

                   if (sint1 > 0.0 ) {
                      A_trans = OldMomentum.cross(theFacetNormal);
                      A_trans = A_trans.unit();
                   } else {
                      A_trans  = OldPolarization;
                   }
                   A_paral   = NewMomentum.cross(A_trans);
                   A_paral   = A_paral.unit();

                   if(iTE>0&&iTM>0) {
                     NewPolarization = 
                           -OldPolarization + (2.*EdotN)*theFacetNormal;
                   } else if (iTE>0) {
                     NewPolarization = -A_trans;
                   } else if (iTM>0) {
                     NewPolarization = -A_paral;
                   }

                }

             }

             OldMomentum = NewMomentum;
             OldPolarization = NewPolarization;

	   }

	} while (NewMomentum * theGlobalNormal < 0.0);
	CustomRunManager* manman = (CustomRunManager*)G4RunManager::GetRunManager();
	G4int chosen_proc = manman->select_photon_BP(to_decider_,G4ThreeVector(0,0,0), NewMomentum);
	if ((RM_CHOOSE_KILL == chosen_proc) || (RM_CHOOSE_DEFL==chosen_proc))
	{
		manman->SetPhEvProb(0);
		manman->SetPhEvType(RM_PHOTON_NO_TRACK);
		manman->process_end(aStep);
	}
	if ((RM_CHOOSE_REFL == chosen_proc) || (RM_CHOOSE_BOTH == chosen_proc))
	{
		manman->SetPhEvProb(net_refl_prob);
		manman->SetPhEvType(RM_PHOTON_REFL);
		manman->process_end(aStep);
	}
}

void CustomOpBoundaryProcess::DielectricLUT()
{
        G4int thetaIndex, phiIndex;
        G4double AngularDistributionValue, thetaRad, phiRad, EdotN;
        G4ThreeVector PerpendicularVectorTheta, PerpendicularVectorPhi;

		theStatus = CustomOpBoundaryProcessStatus(G4int(theFinish) +
                           (G4int(NoRINDEX)-G4int(groundbackpainted)));

        G4int thetaIndexMax = OpticalSurface->GetThetaIndexMax();
        G4int phiIndexMax   = OpticalSurface->GetPhiIndexMax();

        do {
           if ( !G4BooleanRand(theReflectivity) ) // Not reflected, so Absorbed
              DoAbsorption();
           else {
              // Calculate Angle between Normal and Photon Momentum
              G4double anglePhotonToNormal = 
                                          OldMomentum.angle(-theGlobalNormal);
              // Round it to closest integer
              G4int angleIncident = G4int(std::floor(180/pi*anglePhotonToNormal+0.5));

              // Take random angles THETA and PHI, 
              // and see if below Probability - if not - Redo
              do {
                 thetaIndex = G4RandFlat::shootInt(thetaIndexMax-1);
                 phiIndex = G4RandFlat::shootInt(phiIndexMax-1);
                 // Find probability with the new indeces from LUT
                 AngularDistributionValue = OpticalSurface -> 
                   GetAngularDistributionValue(angleIncident,
                                               thetaIndex,
                                               phiIndex);
              } while ( !G4BooleanRand(AngularDistributionValue) );

              thetaRad = (-90 + 4*thetaIndex)*pi/180;
              phiRad = (-90 + 5*phiIndex)*pi/180;
              // Rotate Photon Momentum in Theta, then in Phi
              NewMomentum = -OldMomentum;
              PerpendicularVectorTheta = NewMomentum.cross(theGlobalNormal);
              if (PerpendicularVectorTheta.mag() > kCarTolerance ) {
                 PerpendicularVectorPhi = 
                                  PerpendicularVectorTheta.cross(NewMomentum);
              }
              else {
                 PerpendicularVectorTheta = NewMomentum.orthogonal();
                 PerpendicularVectorPhi =
                                  PerpendicularVectorTheta.cross(NewMomentum);
              }
              NewMomentum =
                 NewMomentum.rotate(anglePhotonToNormal-thetaRad,
                                    PerpendicularVectorTheta);
              NewMomentum = NewMomentum.rotate(-phiRad,PerpendicularVectorPhi);
              // Rotate Polarization too:
              theFacetNormal = (NewMomentum - OldMomentum).unit();
              EdotN = OldPolarization * theFacetNormal;
              NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
           }
        } while (NewMomentum * theGlobalNormal <= 0.0);
}

void CustomOpBoundaryProcess::DielectricDichroic()
{
        // Calculate Angle between Normal and Photon Momentum
        G4double anglePhotonToNormal = OldMomentum.angle(-theGlobalNormal);

        // Round it to closest integer
        G4double angleIncident = std::floor(180/pi*anglePhotonToNormal+0.5);

        if (!DichroicVector) {
           if (OpticalSurface) DichroicVector = OpticalSurface->GetDichroicVector();
        }


        if (DichroicVector) {
           G4double wavelength = h_Planck*c_light/thePhotonMomentum;
           theTransmittance =
             DichroicVector->Value(wavelength/nm,angleIncident,idx,idy)*perCent;
//            G4cout << "wavelength: " << std::floor(wavelength/nm) 
//                                     << "nm" << G4endl;
//            G4cout << "Incident angle: " << angleIncident << "deg" << G4endl;
//            G4cout << "Transmittance: " 
//                   << std::floor(theTransmittance/perCent) << "%" << G4endl;
        } else {
           G4ExceptionDescription ed;
           ed << " G4OpBoundaryProcess/DielectricDichroic(): "
              << " The dichroic surface has no G4Physics2DVector"
              << G4endl;
           G4Exception("G4OpBoundaryProcess::DielectricDichroic", "OpBoun03",
                       FatalException,ed,
                       "A dichroic surface must have an associated G4Physics2DVector");
        }

        if ( !G4BooleanRand(theTransmittance) ) // Not transmitted, so reflect
           DoReflection();
        else {
           theStatus = Dichroic;
           NewMomentum = OldMomentum;
           NewPolarization = OldPolarization;
        }
}

void CustomOpBoundaryProcess::DielectricDielectric(const G4Step* aStep)
{
#ifdef DEBUG_MC_NODES
	G4cout << "COBP: dielectric_dielectric()" << G4endl;
#endif
	G4bool Inside = false;
	G4bool Swap = false;

	leap:

    G4bool Through = false;
	G4bool Done = false;

	do {

	   if (Through) {
	      Swap = !Swap;
	      Through = false;
	      theGlobalNormal = -theGlobalNormal;
	      G4SwapPtr(Material1,Material2);
	      G4SwapObj(&Rindex1,&Rindex2);
	   }
	
	   if ( theFinish == polished ) {
	        theFacetNormal = theGlobalNormal;
	   }
	   else {
	        theFacetNormal =
	                     GetFacetNormal(OldMomentum,theGlobalNormal);
	   }

	   G4double PdotN = OldMomentum * theFacetNormal;
	   G4double EdotN = OldPolarization * theFacetNormal;

	   cost1 = - PdotN;
	   if (std::abs(cost1) < 1.0-kCarTolerance){
	      sint1 = std::sqrt(1.-cost1*cost1);
	      sint2 = sint1*Rindex1/Rindex2;     // *** Snell's Law ***
	   }
	   else {
	      sint1 = 0.0;
	      sint2 = 0.0;
	   }
	   //MY addition
	   G4ThreeVector Refl_momentum;
	   G4ThreeVector Refl_polar;
	   G4ThreeVector Defl_momentum;
	   G4ThreeVector Defl_polar;
	   G4double transmittance_ = 0; //if -1, then total internal refl is proceded

	   if (sint2 >= 1.0) 
	   {
		   transmittance_ = -1;
	      // Simulate total internal reflection
	      if (Swap) Swap = !Swap;
          theStatus = TotalInternalReflection;
		  if (theModel == unified && theFinish != polished)
			  ChooseReflection();
	      if ( theStatus == LambertianReflection ) 
				DoReflection();
		  else
		  {
			  /*this is the only case at the moment. If more complicated methods will be employed
			  the algorithm must be modifyed: consider all possible probabilities and using decider (from RunManager) create new appropriate MC nodes
			  for every type of process with corresponding distributions. At the moment (start of this project) decider checks whether
			  simple total internal reflection (and deflectioin below) should be processed further
			  (considering out momentum and geometry). If no, current process is killed*/
			  if (theStatus == BackScattering) 
			  {
				  NewMomentum = -OldMomentum;
				  NewPolarization = -OldPolarization;
				 
			  }
			  else 
			  {
				  PdotN = OldMomentum * theFacetNormal;
				  NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
				  EdotN = OldPolarization * theFacetNormal;
				  NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
			  }
			  CustomRunManager* manman = (CustomRunManager*)G4RunManager::GetRunManager();
			  manman->SetPhEvType(RM_PHOTON_TOT_REFL);
			  G4int chosen_proc = manman->select_photon_BP(to_decider_, G4ThreeVector(0, 0, 0), NewMomentum);
			  if ((chosen_proc == RM_CHOOSE_BOTH) || (chosen_proc == RM_CHOOSE_REFL))
				  manman->SetPhEvProb(1);
			  else manman->SetPhEvProb(0); //effectively kills the event
			  manman->process_end(aStep);
		  }
	   }
	   else if (sint2 < 1.0) {

	      // Calculate amplitude for transmission (Q = P x N)

	      if (cost1 > 0.0) {
	         cost2 =  std::sqrt(1.-sint2*sint2);
	      }
	      else {
	         cost2 = -std::sqrt(1.-sint2*sint2);
	      }

	      G4ThreeVector A_trans, A_paral, E1pp, E1pl;
	      G4double E1_perp, E1_parl;

	      if (sint1 > 0.0) {
	         A_trans = OldMomentum.cross(theFacetNormal);
                 A_trans = A_trans.unit();
	         E1_perp = OldPolarization * A_trans;
                 E1pp    = E1_perp * A_trans;
                 E1pl    = OldPolarization - E1pp;
                 E1_parl = E1pl.mag();
              }
	      else {
	         A_trans  = OldPolarization;
	         // Here we Follow Jackson's conventions and we set the
	         // parallel component = 1 in case of a ray perpendicular
	         // to the surface
	         E1_perp  = 0.0;
	         E1_parl  = 1.0;
	      }

              G4double s1 = Rindex1*cost1;
              G4double E2_perp = 2.*s1*E1_perp/(Rindex1*cost1+Rindex2*cost2);
              G4double E2_parl = 2.*s1*E1_parl/(Rindex2*cost1+Rindex1*cost2);
              G4double E2_total = E2_perp*E2_perp + E2_parl*E2_parl;
              G4double s2 = Rindex2*cost2*E2_total;

              G4double TransCoeff;

              if (theTransmittance > 0) TransCoeff = theTransmittance;
              else if (cost1 != 0.0) TransCoeff = s2/s1;
              else TransCoeff = 0.0;

	      G4double E2_abs, C_parl, C_perp;
		  transmittance_ = TransCoeff;
	      //if ( !G4BooleanRand(TransCoeff) ) {
			  {
	         //Now - calculate both refl and defl. depricated: Simulate reflection

                 //if (Swap) Swap = !Swap;

		 theStatus = FresnelReflection;

		 if ( theModel == unified && theFinish != polished )
						     ChooseReflection();

		 if ( theStatus == LambertianReflection ) {
		    DoReflection();
		 }
		 else if ( theStatus == BackScattering ) {
		    NewMomentum = -OldMomentum;
		    NewPolarization = -OldPolarization;
		 }
		 else {

                    PdotN = OldMomentum * theFacetNormal;
	            NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

	            if (sint1 > 0.0) {   // incident ray oblique

		       E2_parl   = Rindex2*E2_parl/Rindex1 - E1_parl;
		       E2_perp   = E2_perp - E1_perp;
		       E2_total  = E2_perp*E2_perp + E2_parl*E2_parl;
                       A_paral   = NewMomentum.cross(A_trans);
                       A_paral   = A_paral.unit();
		       E2_abs    = std::sqrt(E2_total);
		       C_parl    = E2_parl/E2_abs;
		       C_perp    = E2_perp/E2_abs;

                       NewPolarization = C_parl*A_paral + C_perp*A_trans;

	            }

	            else {               // incident ray perpendicular

	               if (Rindex2 > Rindex1) {
		          NewPolarization = - OldPolarization;
	               }
	               else {
	                  NewPolarization =   OldPolarization;
	               }

	            }
	         }
		  }
		 Refl_momentum = NewMomentum.unit();
		 Refl_polar = NewPolarization.unit();
	      { // photon gets transmitted

	     // Simulate transmission/refraction
		 theStatus = FresnelRefraction;

	         if (sint1 > 0.0) {      // incident ray oblique

		    G4double alpha = cost1 - cost2*(Rindex2/Rindex1);
		    NewMomentum = OldMomentum + alpha*theFacetNormal;
		    NewMomentum = NewMomentum.unit();
		    PdotN = -cost2;
                    A_paral = NewMomentum.cross(A_trans);
                    A_paral = A_paral.unit();
		    E2_abs     = std::sqrt(E2_total);
		    C_parl     = E2_parl/E2_abs;
		    C_perp     = E2_perp/E2_abs;

                    NewPolarization = C_parl*A_paral + C_perp*A_trans;

	         }
	         else {                  // incident ray perpendicular

		    NewMomentum = OldMomentum;
		    NewPolarization = OldPolarization;

	         }
			 Defl_momentum = NewMomentum.unit();
			 Defl_polar = NewPolarization.unit();
	      }
	   }
	   if (transmittance_ != -1) //useless, this is in 'else' place anyway//not total internal reflections case () //p.s. ugly again
	   {
#ifdef DEBUG_MC_NODES
		   if (!(Defl_momentum * theGlobalNormal > 0.0))
			   G4cout<<"W! Optics Boundary: deflection wrong output momentum, 'done' place"<<G4endl;
		   if (!(Refl_momentum * theGlobalNormal < 0.0))
			   G4cout << "W! Optics Boundary: reflection wrong output momentum, 'done' place" << G4endl;
#endif
		   CustomRunManager* manman = (CustomRunManager*)G4RunManager::GetRunManager();
		   G4int chosen_proc = manman->select_photon_BP(to_decider_,Defl_momentum, Refl_momentum);
		   switch (chosen_proc)
		   {
		   case RM_CHOOSE_KILL:
		   {
			   if (Swap) Swap = !Swap;
			   manman->SetPhEvProb(0);
			   manman->SetPhEvType(RM_PHOTON_NO_TRACK);
			   manman->process_end(aStep);
			   theStatus = FresnelRefraction;
			   NewMomentum = Defl_momentum;
			   NewPolarization = Defl_polar;
			   break;
		   }
		   case RM_CHOOSE_DEFL:
		   {
			   Inside = !Inside;
			   Through = true;
			   manman->SetPhEvProb(transmittance_);
			   manman->SetPhEvType(RM_PHOTON_DEFL);
			   manman->process_end(aStep);
			   theStatus = FresnelRefraction;
			   NewMomentum = Defl_momentum;
			   NewPolarization = Defl_polar;
			   break;
		   }
		   case RM_CHOOSE_REFL:
		   {
			   if (Swap) Swap = !Swap;
			   manman->SetPhEvProb(1.0 - transmittance_);
			   manman->SetPhEvType(RM_PHOTON_REFL);
			   manman->process_end(aStep);
			   theStatus = FresnelReflection;
			   NewMomentum = Refl_momentum;
			   NewPolarization = Refl_polar;
			   break;
		   }
		   case RM_CHOOSE_BOTH: //fork simulation here cause both reflected and deflected photons are of interest
		   {
			   Inside = !Inside;
			   Through = true;
			   manman->SetPhEvProb(transmittance_);
			   manman->SetPhEvType(RM_PHOTON_DEFL);
			   manman->process_end(aStep);
			   theStatus = FresnelRefraction;
			   NewMomentum = Defl_momentum;
			   NewPolarization = Defl_polar;
			   manman->spawn_new_MC_node(to_decider_,1.0-transmittance_, Refl_momentum,Refl_polar);
			   break;
		   }
		   }

	   }
	   OldMomentum = NewMomentum.unit();
	   OldPolarization = NewPolarization.unit();
	  /* if (theStatus == FresnelRefraction) {
	      Done = (NewMomentum * theGlobalNormal <= 0.0);
	   } 
	   else {
	      Done = (NewMomentum * theGlobalNormal >= 0.0);
	   }
*/ //My: I don't see the point in this, how 'Done' may be false? Is here several processes with photon involved?
	   Done = true;
	} while (!Done);
	//My: below are not reworked nor analysed
	if (Inside && !Swap) {
          if( theFinish == polishedbackpainted ||
              theFinish == groundbackpainted ) {

	      if( !G4BooleanRand(theReflectivity) ) {
		DoAbsorption();
              }
	      else {
		if (theStatus != FresnelRefraction ) {
		   theGlobalNormal = -theGlobalNormal;
	        }
	        else {
		   Swap = !Swap;
		   G4SwapPtr(Material1,Material2);
		   G4SwapObj(&Rindex1,&Rindex2);
	        }
		if ( theFinish == groundbackpainted )
					theStatus = LambertianReflection;

	        DoReflection();

	        theGlobalNormal = -theGlobalNormal;
		OldMomentum = NewMomentum;

	        goto leap;
	      }
	  }
	}
}

// GetMeanFreePath
// ---------------
//
G4double CustomOpBoundaryProcess::GetMeanFreePath(const G4Track&,
                                              G4double ,
                                              G4ForceCondition* condition)
{
	*condition = Forced;

	return DBL_MAX;
}

G4double CustomOpBoundaryProcess::GetIncidentAngle()
{
    G4double PdotN = OldMomentum * theFacetNormal;
    G4double magP= OldMomentum.mag();
    G4double magN= theFacetNormal.mag();
    G4double incidentangle = pi - std::acos(PdotN/(magP*magN));

    return incidentangle;
}

G4double CustomOpBoundaryProcess::GetReflectivity(G4double E1_perp,
                                              G4double E1_parl,
                                              G4double incidentangle,
                                              G4double RealRindex,
                                              G4double ImaginaryRindex)
{

  G4complex Reflectivity, Reflectivity_TE, Reflectivity_TM;
  G4complex N(RealRindex, ImaginaryRindex);
  G4complex CosPhi;

  G4complex u(1,0);           //unit number 1

  G4complex numeratorTE;      // E1_perp=1 E1_parl=0 -> TE polarization
  G4complex numeratorTM;      // E1_parl=1 E1_perp=0 -> TM polarization
  G4complex denominatorTE, denominatorTM;
  G4complex rTM, rTE;

  // Following two equations, rTM and rTE, are from: "Introduction To Modern
  // Optics" written by Fowles

  CosPhi=std::sqrt(u-((std::sin(incidentangle)*std::sin(incidentangle))/(N*N)));

  numeratorTE   = std::cos(incidentangle) - N*CosPhi;
  denominatorTE = std::cos(incidentangle) + N*CosPhi;
  rTE = numeratorTE/denominatorTE;

  numeratorTM   = N*std::cos(incidentangle) - CosPhi;
  denominatorTM = N*std::cos(incidentangle) + CosPhi;
  rTM = numeratorTM/denominatorTM;

  // This is my calculaton for reflectivity on a metalic surface
  // depending on the fraction of TE and TM polarization
  // when TE polarization, E1_parl=0 and E1_perp=1, R=abs(rTE)^2 and
  // when TM polarization, E1_parl=1 and E1_perp=0, R=abs(rTM)^2

  Reflectivity_TE =  (rTE*conj(rTE))*(E1_perp*E1_perp)
                    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity_TM =  (rTM*conj(rTM))*(E1_parl*E1_parl)
                    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity    = Reflectivity_TE + Reflectivity_TM;

  do {
     if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TE))
       {iTE = -1;}else{iTE = 1;}
     if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TM))
       {iTM = -1;}else{iTM = 1;}
  } while(iTE<0&&iTM<0);

  return real(Reflectivity);

}

void CustomOpBoundaryProcess::CalculateReflectivity()
{
  G4double RealRindex =
           PropertyPointer1->Value(thePhotonMomentum);
  G4double ImaginaryRindex =
           PropertyPointer2->Value(thePhotonMomentum);

  // calculate FacetNormal
  if ( theFinish == ground ) {
     theFacetNormal =
               GetFacetNormal(OldMomentum, theGlobalNormal);
  } else {
     theFacetNormal = theGlobalNormal;
  }

  G4double PdotN = OldMomentum * theFacetNormal;
  cost1 = -PdotN;

  if (std::abs(cost1) < 1.0 - kCarTolerance) {
     sint1 = std::sqrt(1. - cost1*cost1);
  } else {
     sint1 = 0.0;
  }

  G4ThreeVector A_trans, A_paral, E1pp, E1pl;
  G4double E1_perp, E1_parl;

  if (sint1 > 0.0 ) {
     A_trans = OldMomentum.cross(theFacetNormal);
     A_trans = A_trans.unit();
     E1_perp = OldPolarization * A_trans;
     E1pp    = E1_perp * A_trans;
     E1pl    = OldPolarization - E1pp;
     E1_parl = E1pl.mag();
  }
  else {
     A_trans  = OldPolarization;
     // Here we Follow Jackson's conventions and we set the
     // parallel component = 1 in case of a ray perpendicular
     // to the surface
     E1_perp  = 0.0;
     E1_parl  = 1.0;
  }

  //calculate incident angle
  G4double incidentangle = GetIncidentAngle();

  //calculate the reflectivity depending on incident angle,
  //polarization and complex refractive

  theReflectivity =
             GetReflectivity(E1_perp, E1_parl, incidentangle,
                                                 RealRindex, ImaginaryRindex);
}

G4bool CustomOpBoundaryProcess::InvokeSD(const G4Step* pStep)
{
  G4Step aStep = *pStep;

  aStep.AddTotalEnergyDeposit(thePhotonMomentum);

  G4VSensitiveDetector* sd = aStep.GetPostStepPoint()->GetSensitiveDetector();
  if (sd) return sd->Hit(&aStep);
  else return false;
}
