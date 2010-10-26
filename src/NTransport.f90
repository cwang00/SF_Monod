

MODULE NTransport

   IMPLICIT NONE
   REAL, PARAMETER :: MYPAR = 3.D7
   INTEGER, PARAMETER :: TotalNumberOfSpecies = 11
   REAL, PARAMETER :: NO3MoleculeWeight = 62
   REAL, PARAMETER :: NH4MoleculeWeight = 18
   REAL, PARAMETER :: CH2OMoleculeWeight = 30
   REAL, PARAMETER :: O2MoleculeWeight = 32 
   REAL, PARAMETER :: CO2MoleculeWeight = 44
   REAL, PARAMETER :: HCO3MoleculeWeight = 61
   REAL, PARAMETER :: CO3MoleculeWeight = 60
   REAL, PARAMETER :: HMoleculeWeight = 1
   REAL, PARAMETER :: N2MoleculeWeight = 28
   REAL, PARAMETER :: OHMoleculeWeight = 17
   REAL, PARAMETER :: CaMoleculeWeight = 40

   
   TYPE NTransParameters
     CHARACTER(10), DIMENSION(TotalNumberOfSpecies)::SpeciesNames
     REAL :: Biomass1Conc_mg_l
     REAL :: Biomass2Conc_mg_l
     REAL :: Biomass1YieldCoeff_M_biomass_M_substrate
     REAL :: Biomass2YieldCoeff_M_biomass_M_substrate
     REAL :: EmpericalBiomass1InhibitionConst_mg_l
     REAL :: EmpericalBiomass2InhibitionConst_mg_l
     REAL :: CH2OHalfSaturationConst_mg_l
     REAL :: O2HalfSaturationConst_mg_l
     REAL :: NH4HalfSaturationConst_mg_l
     REAL :: NO3HalfSaturationConst_mg_l
     REAL :: O2InhibitionCoef_mg_l
     REAL :: CH2ODistCoeff_cm3_g
     REAL :: NH4DistCoeff_cm3_g
     REAL :: Biomass1SpecDecayConst_1_day
     REAL :: Biomass2SpecDecayConst_1_day
     REAL :: O2HenryConst
     REAL :: N2HenryConst
     REAL :: CO2HenryConst
     REAL :: SubstrateFreeAirDiffusionCoeff_m2_day
     REAL :: EAcceptorFreeAirDiffusionCoeff_m2_day
     REAL :: OrgCarbonOxMaxRate_1_day
     REAL :: DenitOxMaxRate_1_day
     REAL :: NitOxMaxRate_1_day
     REAL :: K_f_CO2_1_day
     REAL :: K_b_CO2_1_mday
     REAL :: K_f_H_1_mday
     REAL :: K_b_H_1_day
     REAL :: K_f_HCO3_1_day
     REAL :: K_b_HCO3_1_mday
     REAL :: K_f_H2O_1_day
     REAL :: K_b_H2O_1_mday
     REAL :: K_f_1_cm_day
     REAL :: K_b_1_cm_mday
     REAL :: K_f_2_cm_day
     REAL :: K_b_2_cm_m2day
     REAL :: K_f_3_mol_cm2day
     REAL :: K_b_3_cm_mday
     REAL :: Acc_1_cm
     REAL, DIMENSION(TotalNumberOfSpecies)::backgroundConc
     REAL :: MYNPARS
   END TYPE NTransParameters

   REAL, ALLOCATABLE :: NitrifierConc(:,:,:)
   REAL, ALLOCATABLE :: DenitrifierConc(:,:,:)

   TYPE( NTransParameters ) :: npars

   CONTAINS

     SUBROUTINE ReadNTransPars( nparfile, nx, ny, nz )
      IMPLICIT NONE

      CHARACTER(LEN=100), INTENT(IN) :: nparfile
      INTEGER :: I, J, K, idx, nx, ny, nz
      OPEN(98, FILE = nparfile, STATUS = 'old' ) 

      DO I = 1, TotalNumberOfSpecies
         READ(98,*)  idx, npars%SpeciesNames( idx )
      END DO
     READ(98,*) npars%Biomass1Conc_mg_l
     READ(98,*) npars%Biomass2Conc_mg_l
     READ(98,*) npars%Biomass1YieldCoeff_M_biomass_M_substrate
     READ(98,*) npars%Biomass2YieldCoeff_M_biomass_M_substrate
     READ(98,*) npars%EmpericalBiomass1InhibitionConst_mg_l
     READ(98,*) npars%EmpericalBiomass2InhibitionConst_mg_l
     READ(98,*) npars%CH2OHalfSaturationConst_mg_l
     READ(98,*) npars%O2HalfSaturationConst_mg_l
     READ(98,*) npars%NH4HalfSaturationConst_mg_l
     READ(98,*) npars%NO3HalfSaturationConst_mg_l
     READ(98,*) npars%O2InhibitionCoef_mg_l
     READ(98,*) npars%CH2ODistCoeff_cm3_g
     READ(98,*) npars%NH4DistCoeff_cm3_g
     READ(98,*) npars%Biomass1SpecDecayConst_1_day
     READ(98,*) npars%Biomass2SpecDecayConst_1_day
     READ(98,*) npars%O2HenryConst
     READ(98,*) npars%N2HenryConst
     READ(98,*) npars%CO2HenryConst
     READ(98,*) npars%SubstrateFreeAirDiffusionCoeff_m2_day
     READ(98,*) npars%EAcceptorFreeAirDiffusionCoeff_m2_day
     READ(98,*) npars%OrgCarbonOxMaxRate_1_day
     READ(98,*) npars%DenitOxMaxRate_1_day
     READ(98,*) npars%NitOxMaxRate_1_day
     READ(98,*) npars%K_f_CO2_1_day
     READ(98,*) npars%K_b_CO2_1_mday
     READ(98,*) npars%K_f_H_1_mday
     READ(98,*) npars%K_b_H_1_day
     READ(98,*) npars%K_f_HCO3_1_day
     READ(98,*) npars%K_b_HCO3_1_mday
     READ(98,*) npars%K_f_H2O_1_day
     READ(98,*) npars%K_b_H2O_1_mday
     READ(98,*) npars%K_f_1_cm_day
     READ(98,*) npars%K_b_1_cm_mday
     READ(98,*) npars%K_f_2_cm_day
     READ(98,*) npars%K_b_2_cm_m2day
     READ(98,*) npars%K_f_3_mol_cm2day
     READ(98,*) npars%K_b_3_cm_mday
     READ(98,*) npars%Acc_1_cm
     DO I = 1, TotalNumberOfSpecies
        READ(98,*)  npars%backgroundConc( I )
     END DO
     CLOSE(98)
  
     ALLOCATE( NitrifierConc(nx, ny, nz ), DenitrifierConc(nx, ny, nz) )

     DO I = 1, nx
       DO J = 1, ny
         DO K = 1, nz
          NitrifierConc( I, J, K ) = npars%Biomass1Conc_mg_l
          DenitrifierConc( I, J, K ) = npars%Biomass2Conc_mg_l
         END DO  
       END DO
     END DO

  END SUBROUTINE ReadNTransPars
      
  INTEGER FUNCTION SpeciesNameToInd( n )
   IMPLICIT NONE
   INTEGER :: I
   CHARACTER(*) :: n 
   CHARACTER(10) :: leftn

   SpeciesNameToInd = -1 
   DO I = 1, TotalNumberOfSpecies
     leftn = TRIM( npars%SpeciesNames( I ))
     IF ( leftn .EQ. n ) THEN
       SpeciesNameToInd = I
       EXIT 
     END IF
   END DO
   IF ( SpeciesNameToInd .LT. 1 ) THEN
       WRITE(*,*) 'ERROR: Can'' find Species ', n
       STOP
   ENDIF
  END FUNCTION SpeciesNameToInd

  REAL FUNCTION DeltaX1( D_NH4, X1, Y1, kd1 )
   IMPLICIT NONE
   REAL :: D_NH4, X1, Y1, kd1

   DeltaX1 = -1 * D_NH4 * Y1 - X1 * kd1
  END FUNCTION 

  REAL FUNCTION DeltaX2( D_CH2O, X2, Y2, kd2 )
   IMPLICIT NONE
   REAL :: D_CH2O, X2, Y2, kd2

   DeltaX2 = -1 * D_CH2O * Y2 - X2 * kd2
  END FUNCTION 

  REAL FUNCTION BiomassGrowthInhibition( Kx, X )
   IMPLICIT NONE
   REAL :: Kx, X
   BiomassGrowthInhibition = Kx / ( Kx + X )
  END FUNCTION

  REAL FUNCTION F_X1( X1 )
   IMPLICIT NONE
   REAL :: X1
   F_X1 = BiomassGrowthInhibition(                                     &
              npars%EmpericalBiomass1InhibitionConst_mg_l, X1 )

  END FUNCTION
  
  REAL FUNCTION F_X2( X2 )
   IMPLICIT NONE
   REAL :: X2
   F_X2 = BiomassGrowthInhibition(                                     &
              npars%EmpericalBiomass2InhibitionConst_mg_l, X2 )

  END FUNCTION

 ! see K.T.B. MacQuarrie, E.A. Sudicky/ Journal of Contaminant Hydrology
 ! 47 (2001) 53-84
  REAL FUNCTION R1( Kmax_ox, X2, CH2O, O2, K_CH2O, K_O2 )
   IMPLICIT NONE
   REAL :: Kmax_ox, X2, CH2O, O2, K_CH2O, K_O2
   R1 = Kmax_ox * X2 * F_X2 ( X2 ) * CH2O / ( CH2O + K_CH2O ) *  &
        O2 / (O2 + K_O2 )
  END FUNCTION

  REAL FUNCTION R2( Kmax_nit, X1, NH4, O2, K_NH4, K_O2 )
   IMPLICIT NONE
   REAL :: Kmax_nit, X1, NH4, O2, K_NH4, K_O2
   R2 = Kmax_nit * X1 * F_X1 ( X1 ) * NH4 / ( NH4 + K_NH4 ) *  &
        O2 / (O2 + K_O2 )
  END FUNCTION

  REAL FUNCTION R3( Kmax_denit, X2, CH2O, NO3, O2, K_NO3, K_CH2O, K_O2I )
   IMPLICIT NONE
   REAL :: Kmax_denit, X2, CH2O, NO3, O2, K_NO3, K_CH2O, K_O2I
   R3 = Kmax_denit * X2 * F_X2( X2 ) * CH2O / ( CH2O + K_CH2O ) *  &
        NO3 / (NO3 + K_NO3 ) * K_O2I / (O2 + K_O2I )
  END FUNCTION

  REAL FUNCTION D_CH2O( CH2O, O2, NO3, X2 )
   IMPLICIT NONE
   REAL :: CH2O, O2, NO3, X2
!   IF( CH2O .LT. 0.0 ) CH2O = 0.0
!   IF( O2 .LT. 0.0 ) O2 = 0.0
!   IF( NO3 .LT. 0.0 ) NO3 = 0.0
!   IF( X2 .LT. 0.0 ) X2 = 0.0
   D_CH2O = -1.0 * R1( npars%OrgCarbonOxMaxRate_1_day,  X2, CH2O, O2,     &
                 npars%CH2OHalfSaturationConst_mg_l,                      &
                 npars%O2HalfSaturationConst_mg_l )  -                    &
                 R3( npars%DenitOxMaxRate_1_day,  X2,                     &
                     CH2O, NO3,O2,                                        &
                     npars%NO3HalfSaturationConst_mg_l,                   &
                     npars%CH2OHalfSaturationConst_mg_l,                  &
                     npars%O2InhibitionCoef_mg_l )
!   WRITE(*,*) 'D_CH2O = ', D_CH2O
  END FUNCTION

  REAL FUNCTION D_NH4( NH4, O2, X1 )
   IMPLICIT NONE
   REAL :: NH4, O2, X1

!   IF( NH4 .LT. 0.0 ) NH4 = 0.0
!   IF( O2 .LT. 0.0 ) O2 = 0.0
!   IF( X1 .LT. 0.0 ) X1 = 0.0
   D_NH4 = -R2( npars%NitOxMaxRate_1_day,  X1, NH4, O2,                   &
                      npars%NH4HalfSaturationConst_mg_l,                  &
                      npars%O2HalfSaturationConst_mg_l )
!   WRITE(*,*) 'D_NH4 = ', D_NH4
  END FUNCTION

  REAL FUNCTION D_NO3( NO3, NH4, O2, X1, X2, CH2O ) 
   IMPLICIT NONE
   REAL :: NO3, NH4, O2, X1, X2, CH2O
!
!   IF( NO3 .LT. 0.0 ) NO3 = 0.0
!   IF( NH4 .LT. 0.0 ) NH4 = 0.0
!   IF( O2 .LT. 0.0 ) O2 = 0.0
!   IF( X1 .LT. 0.0 ) X1 = 0.0
!   IF( X2 .LT. 0.0 ) X2 = 0.0
!   IF( CH2O .LT. 0.0 ) CH2O = 0.0
   D_NO3 = NO3MoleculeWeight / NH4MoleculeWeight *                           &
             R2( npars%NitOxMaxRate_1_day,  X1, NH4, O2,                     &
                      npars%NH4HalfSaturationConst_mg_l,                     &
                      npars%O2HalfSaturationConst_mg_l )                     &
           -4 * NO3MoleculeWeight / ( 5 * CH2OMoleculeWeight ) *              &
                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                       &
                     CH2O, NO3, O2,                                          &
                     npars%NO3HalfSaturationConst_mg_l,                      &
                     npars%CH2OHalfSaturationConst_mg_l,                     &
                     npars%O2InhibitionCoef_mg_l )
!   WRITE(*,*) 'D_NO3 = ', D_NO3
  END FUNCTION

  REAL FUNCTION D_O2( NH4, O2, X1, X2, CH2O ) 
   IMPLICIT NONE
   REAL :: NH4, O2, X1, X2, CH2O

!   IF( NH4 .LT. 0.0 ) NH4 = 0.0
!   IF( O2 .LT. 0.0 ) O2 = 0.0
!   IF( X1 .LT. 0.0 ) X1 = 0.0
!   IF( X2 .LT. 0.0 ) X2 = 0.0
!   IF( CH2O .LT. 0.0 ) CH2O = 0.0
   D_O2 = -2 * O2MoleculeWeight / NH4MoleculeWeight *                        &
             R2( npars%NitOxMaxRate_1_day,  X1, NH4, O2,                     &
                      npars%NH4HalfSaturationConst_mg_l,                     &
                      npars%O2HalfSaturationConst_mg_l )                     &
                - O2MoleculeWeight / CH2OMoleculeWeight *                    &
                R1( npars%OrgCarbonOxMaxRate_1_day,  X2, CH2O, O2,           &
                    npars%CH2OHalfSaturationConst_mg_l,                      &
                 npars%O2HalfSaturationConst_mg_l )
!   WRITE(*,*) 'D_O2 = ', D_O2
  END FUNCTION

  REAL FUNCTION D_CO2( O2, X2, NO3, CH2O, CO2, HCO3, H, OH, Ca ) 
   IMPLICIT NONE
   REAL :: O2, X2, NO3, CH2O, CO2, HCO3, H, OH, Ca

!   WRITE(*,*) '----------D_CO2------------'
!   WRITE(*,*) 'O2 = ', O2
!   WRITE(*,*) 'X2 = ', X2
!   WRITE(*,*) 'NO3 = ', NO3 
!   WRITE(*,*) 'CH2O = ', CH2O 
!   WRITE(*,*) 'CO2 = ', CO2 
!   WRITE(*,*) 'HCO3 = ', HCO3 
!   WRITE(*,*) 'H = ', H 
!   WRITE(*,*) 'OH = ', OH 
!   WRITE(*,*) 'Ca = ', Ca 

!   IF( O2 .LT. 0.0 ) O2 = 0.0
!   IF( X2 .LT. 0.0 ) X2 = 0.0
!   IF( NO3 .LT. 0.0 ) NO3 = 0.0
!   IF( CH2O .LT. 0.0 ) CH2O = 0.0
!   IF( CO2 .LT. 0.0 ) CO2 = 0.0
!   IF( HCO3 .LT. 0.0 ) HCO3 = 0.0
!   IF( H .LT. 0.0 ) H = 0.0
!   IF( OH .LT. 0.0 ) OH = 0.0
!   IF( Ca .LT. 0.0 ) Ca = 0.0
!   WRITE(*,*) '----------D_CO2 reset ------------'
!   WRITE(*,*) 'O2 = ', O2
!   WRITE(*,*) 'X2 = ', X2
!   WRITE(*,*) 'NO3 = ', NO3 
!   WRITE(*,*) 'CH2O = ', CH2O 
!   WRITE(*,*) 'CO2 = ', CO2 
!   WRITE(*,*) 'HCO3 = ', HCO3 
!   WRITE(*,*) 'H = ', H 
!   WRITE(*,*) 'OH = ', OH 
!   WRITE(*,*) 'Ca = ', Ca 
!   WRITE( *, *) '1 = ',                                                      &
!      CO2MoleculeWeight / ( 5 * CH2OMoleculeWeight ) *                       &
!                R1( npars%OrgCarbonOxMaxRate_1_day,  X2, CH2O, O2,           &
!                    npars%CH2OHalfSaturationConst_mg_l,                      &
!                 npars%O2HalfSaturationConst_mg_l )
!   WRITE( *, *) '2 = ',                                                      &
!             CO2MoleculeWeight / CH2OMoleculeWeight *                        &
!                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                       &
!                     CH2O, NO3, O2,                                          &
!                     npars%NO3HalfSaturationConst_mg_l,                      &
!                     npars%CH2OHalfSaturationConst_mg_l,                     &
!                     npars%O2InhibitionCoef_mg_l ) 
!   WRITE( *, *) '3 = ',                                                      &
!             - npars%K_f_CO2_1_day * CO2
!   WRITE( *, *) '4 = ',                                                      &
!        npars%K_b_CO2_1_mday * ( HCO3 / 1000 / HCO3MoleculeWeight)           &
!         * ( H /1000 / HMoleculeWeight ) * CO2MoleculeWeight * 1000
!   WRITE( *, *) '5 = ',                                                      &
!               -npars%K_f_H_1_mday * CO2 *                                   & 
!               ( OH / 1000 / OHMoleculeWeight )
!   WRITE( *, *) '6 = ',                                                      &
!          npars%K_b_H_1_day * HCO3 / HCO3MoleculeWeight * CO2Moleculeweight
!   WRITE( *, *) '7 = ',                                                      &
!               -npars%Acc_1_cm * npars%K_f_2_cm_day * CO2
!   WRITE( *, *) '8 = ',                                                      &
!                npars%Acc_1_cm * npars%K_b_2_cm_m2day *                      &
!                ( HCO3 /1000 / HCO3MoleculeWeight)**2 *                      &
!                 ( Ca /1000 /CaMoleculeWeight )                              &
!                * CO2MoleculeWeight * 1000

     D_CO2 = CO2MoleculeWeight / CH2OMoleculeWeight *                        &
                R1( npars%OrgCarbonOxMaxRate_1_day,  X2, CH2O, O2,           &
                    npars%CH2OHalfSaturationConst_mg_l,                      &
                 npars%O2HalfSaturationConst_mg_l ) +                        &
             CO2MoleculeWeight / ( 5 * CH2OMoleculeWeight ) *                &
                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                       &
                     CH2O, NO3, O2,                                          &
                     npars%NO3HalfSaturationConst_mg_l,                      &
                     npars%CH2OHalfSaturationConst_mg_l,                     &
                     npars%O2InhibitionCoef_mg_l )                           &
             - npars%K_f_CO2_1_day * CO2 +                                   & 
        npars%K_b_CO2_1_mday * ( HCO3 / 1000 / HCO3MoleculeWeight)           &
         * ( H /1000 / HMoleculeWeight ) * CO2MoleculeWeight * 1000  -       &
               npars%K_f_H_1_mday * CO2 *                                    & 
               ( OH / 1000 / OHMoleculeWeight ) +                            &
          npars%K_b_H_1_day * HCO3 / HCO3MoleculeWeight * CO2Moleculeweight  &
               -npars%Acc_1_cm * npars%K_f_2_cm_day * CO2 +                  &
                npars%Acc_1_cm * npars%K_b_2_cm_m2day *                      &
                ( HCO3 /1000 / HCO3MoleculeWeight)**2 *                      &
                 ( Ca /1000 /CaMoleculeWeight )                              &
                * CO2MoleculeWeight * 1000
!   WRITE(*,*) 'D_CO2 = ', D_CO2
  END FUNCTION

  REAL FUNCTION D_HCO3( O2, X2, NO3, CH2O, CO2, HCO3, H, OH, Ca, CO3 ) 
   IMPLICIT NONE
   REAL :: O2, X2, NO3, CH2O, CO2, HCO3, H,OH, Ca, CO3

!   WRITE(*,*) '----------D_HCO3------------'
!   WRITE(*,*) 'O2 = ', O2
!   WRITE(*,*) 'X2 = ', X2
!   WRITE(*,*) 'NO3 = ', NO3
!   WRITE(*,*) 'CH2O = ', CH2O
!   WRITE(*,*) 'CO2 = ', CO2
!   WRITE(*,*) 'HCO3 = ', HCO3
!   WRITE(*,*) 'H = ', H
!   WRITE(*,*) 'OH = ', OH
!   WRITE(*,*) 'Ca = ', Ca
!   WRITE(*,*) 'CO3 = ', CO3
!
!   IF( O2 .LT. 0.0 ) O2 = 0.0
!   IF( X2 .LT. 0.0 ) X2 = 0.0
!   IF( NO3 .LT. 0.0 ) NO3 = 0.0
!   IF( CH2O .LT. 0.0 ) CH2O = 0.0
!   IF( CO2 .LT. 0.0 ) CO2 = 0.0
!   IF( HCO3 .LT. 0.0 ) HCO3 = 0.0
!   IF( H .LT. 0.0 ) H = 0.0
!   IF( OH .LT. 0.0 ) OH = 0.0
!   IF( Ca .LT. 0.0 ) Ca = 0.0
!   IF( CO3 .LT. 0.0 ) CO3 = 0.0
!   WRITE(*,*) '----------D_HCO3 reset------------'
!   WRITE(*,*) 'O2 = ', O2
!   WRITE(*,*) 'X2 = ', X2
!   WRITE(*,*) 'NO3 = ', NO3
!   WRITE(*,*) 'CH2O = ', CH2O
!   WRITE(*,*) 'CO2 = ', CO2
!   WRITE(*,*) 'HCO3 = ', HCO3
!   WRITE(*,*) 'H = ', H
!   WRITE(*,*) 'OH = ', OH
!   WRITE(*,*) 'Ca = ', Ca
!   WRITE(*,*) 'CO3 = ', CO3
!
!   WRITE(*,*) '1 = ',                                                        &
!            HCO3MoleculeWeight / CH2OMoleculeWeight *                        &
!                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                       &
!                     CH2O, NO3, O2,                                          &
!                     npars%NO3HalfSaturationConst_mg_l,                      &
!                     npars%CH2OHalfSaturationConst_mg_l,                     &
!                     npars%O2InhibitionCoef_mg_l )
!   WRITE(*,*) '2 = ',                                                        &
!       HCO3MoleculeWeight / CO2MoleculeWeight * npars%K_f_CO2_1_day * CO2
!   WRITE(*,*) '3 = ',                                                        &
!        -npars%K_b_CO2_1_mday * HCO3 * ( H / 1000 / HMoleculeWeight )
!   WRITE(*,*) '4 = ',                                                        &
!      npars%K_f_H_1_mday * ( CO2 /1000 / CO2MoleculeWeight )               &
!            * ( OH / 1000 /OHMoleculeWeight ) * HCO3MoleculeWeight * 1000 
!   WRITE(*,*) '5 = ', -npars%K_b_H_1_day * HCO3 
!   WRITE(*,*) '6 = ',                                                        &
!        HCO3MoleculeWeight / CO3MoleculeWeight * npars%K_f_HCO3_1_day * CO3
!   WRITE(*,*) '7 = ',                                                        &
!      - npars%K_b_HCO3_1_mday * HCO3 * ( OH / 1000 / OHMoleculeWeight)
!   WRITE(*,*) '8 = ',                                                        &
!          npars%Acc_1_cm * HCO3MoleculeWeight / HMoleculeWeight *           &
!                                   npars%K_f_1_cm_day * H
!   WRITE(*,*) '9 = ',                                                        &
!         - npars%Acc_1_cm * npars%K_b_1_cm_mday *                            &
!             HCO3 * ( Ca / 1000 / CaMoleculeWeight ) 
!   WRITE(*,*) '10 = ',                                                        &
!          npars%Acc_1_cm * HCO3MoleculeWeight / CO2MoleculeWeight *         &
!                     npars%K_f_2_cm_day *  CO2
!   WRITE(*,*) '11 = ',                                                        &
!         - npars%Acc_1_cm * npars%K_b_2_cm_m2day * HCO3MoleculeWeight * 1000 &
!            * ( HCO3 / 1000 / HCO3MoleculeWeight) ** 2  *                    &
!              ( Ca / 1000 / CaMoleculeWeight )
!

   D_HCO3 = ( 4 * HCO3MoleculeWeight ) / ( 5 *  CH2OMoleculeWeight ) *       &
                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                       &
                     CH2O, NO3, O2,                                          &
                     npars%NO3HalfSaturationConst_mg_l,                      &
                     npars%CH2OHalfSaturationConst_mg_l,                     &
                     npars%O2InhibitionCoef_mg_l )                           &
      + HCO3MoleculeWeight / CO2MoleculeWeight * npars%K_f_CO2_1_day * CO2 - &
        npars%K_b_CO2_1_mday * HCO3 * ( H / 1000 / HMoleculeWeight )         &
      + npars%K_f_H_1_mday * ( CO2 /1000 / CO2MoleculeWeight )               &
            * ( OH / 1000 /OHMoleculeWeight ) * HCO3MoleculeWeight * 1000 -  &
         npars%K_b_H_1_day * HCO3                                            &
       + HCO3MoleculeWeight / CO3MoleculeWeight * npars%K_f_HCO3_1_day * CO3 &
      - npars%K_b_HCO3_1_mday * HCO3 * ( OH / 1000 / OHMoleculeWeight)       &
         + npars%Acc_1_cm * HCO3MoleculeWeight / HMoleculeWeight *           &
                                   npars%K_f_1_cm_day * H                    &
         - npars%Acc_1_cm * npars%K_b_1_cm_mday *                            &
             HCO3 * ( Ca / 1000 / CaMoleculeWeight )                         &
         + npars%Acc_1_cm * HCO3MoleculeWeight / CO2MoleculeWeight *         &
                     npars%K_f_2_cm_day *  CO2                               &
         - npars%Acc_1_cm * npars%K_b_2_cm_m2day * HCO3MoleculeWeight * 1000 &
            * ( HCO3 / 1000 / HCO3MoleculeWeight) ** 2  *                    &
              ( Ca / 1000 / CaMoleculeWeight )
!   WRITE(*,*) 'D_HCO3 = ', D_HCO3
  END FUNCTION

  REAL FUNCTION D_H( O2, NH4, X1, CO2, HCO3, H, OH, Ca ) 
   IMPLICIT NONE
   REAL :: O2, NH4, X1, CO2, HCO3, H, OH, Ca
!   IF( O2 .LT. 0.0 ) O2 = 0.0
!   IF( NH4 .LT. 0.0 ) NH4 = 0.0
!   IF( X1 .LT. 0.0 ) X1 = 0.0
!   IF( CO2 .LT. 0.0 ) CO2 = 0.0
!   IF( HCO3 .LT. 0.0 ) HCO3 = 0.0
!   IF( H .LT. 0.0 ) H = 0.0
!   IF( OH .LT. 0.0 ) OH = 0.0
!   IF( Ca .LT. 0.0 ) Ca = 0.0
   D_H = HMoleculeWeight / NH4MoleculeWeight *                               &
             R2( npars%NitOxMaxRate_1_day,  X1, NH4, O2,                     &
                      npars%NH4HalfSaturationConst_mg_l,                     &
                      npars%O2HalfSaturationConst_mg_l )  * 2                &
   + HMoleculeWeight / CO2MoleculeWeight * npars%K_f_CO2_1_day * CO2         &
  - npars%K_b_CO2_1_mday * ( HCO3 / 1000 / HCO3MoleculeWeight ) *  H         &
         - npars%Acc_1_cm * npars%K_f_1_cm_day * H +                         &
         npars%Acc_1_cm * npars%K_b_1_cm_mday *                              &
            ( HCO3 / 1000 / HCO3MoleculeWeight ) *                           &
            ( Ca / 1000 / CaMoleculeWeight ) *                               &
               HMoleculeWeight * 1000                                        &
         + npars%K_f_H2O_1_day * HMoleculeWeight * 1000 -                    &
           npars%K_b_H2O_1_mday * ( OH / 1000 / OHMoleculeWeight ) * H 
!   WRITE(*,*) 'D_H = ', D_H
  END FUNCTION

  REAL FUNCTION D_OH( H, OH, CO2, HCO3, CO3  ) 
   IMPLICIT NONE
   REAL :: OH, H, CO2, HCO3, CO3
!   IF( H .LT. 0.0 ) H = 0.0
!   IF( OH .LT. 0.0 ) OH = 0.0
!   IF( CO2 .LT. 0.0 ) CO2 = 0.0
!   IF( CO3 .LT. 0.0 ) CO3 = 0.0
!   IF( HCO3 .LT. 0.0 ) HCO3 = 0.0
   D_OH = npars%K_f_H2O_1_day * OHMoleculeWeight * 1000 -                    &
          npars%K_b_H2O_1_mday * OH * ( H / 1000 / HMoleculeWeight )         &
          -npars%K_f_H_1_mday * ( CO2 / 1000 / CO2MoleculeWeight ) * OH      &
          +npars%K_b_H_1_day * HCO3 / HCO3MoleculeWeight * OHMoleculeweight  &
       + OHMoleculeWeight / CO3MoleculeWeight * npars%K_f_HCO3_1_day * CO3   &
       - npars%K_b_HCO3_1_mday * ( HCO3 / 1000 /HCO3MoleculeWeight ) * OH

!   WRITE(*,*) 'D_OH = ', D_OH
  END FUNCTION

  REAL FUNCTION D_CO3( CO3, HCO3, OH, Ca ) 
   IMPLICIT NONE
   REAL :: CO3, HCO3, OH, Ca

!   WRITE(*,*) '----------D_CO3------------'
!   WRITE(*,*) 'CO3 = ', CO3
!   WRITE(*,*) 'HCO3 = ', HCO3
!   WRITE(*,*) 'OH = ', OH
!   WRITE(*,*) 'Ca = ', Ca
!
!   IF( CO3 .LT. 0.0 ) CO3 = 0.0
!   IF( HCO3 .LT. 0.0 ) HCO3 = 0.0
!   IF( OH .LT. 0.0 ) OH = 0.0
!   IF( Ca .LT. 0.0 ) Ca = 0.0
!   WRITE(*,*) '----------D_CO3-reset-----------'
!   WRITE(*,*) 'CO3 = ', CO3
!   WRITE(*,*) 'HCO3 = ', HCO3
!   WRITE(*,*) 'OH = ', OH
!   WRITE(*,*) 'Ca = ', Ca
!   WRITE(*,*) '1 = ', -npars%K_f_HCO3_1_day * CO3
!   WRITE(*,*) '2 = ', npars%K_b_HCO3_1_mday *                                &
!                                 ( HCO3 /1000 / HCO3MoleculeWeight ) *       &
!         ( OH / 1000 / OHMoleculeWeight ) * CO3MoleculeWeight * 1000 
!   WRITE(*,*), '3 = ',                                                       &
!        npars%Acc_1_cm * npars%K_f_3_mol_cm2day * CO3MoleculeWeight * 1000
!   WRITE(*,*), '4 = ',                                                       &
!        -  npars%Acc_1_cm * npars%K_b_3_cm_mday * CO3 *                      &
!           ( Ca / 1000 /CaMoleculeWeight )

   D_CO3 = -npars%K_f_HCO3_1_day * CO3                                       &
       + npars%K_b_HCO3_1_mday * ( HCO3 /1000 / HCO3MoleculeWeight ) *       &
         ( OH / 1000 / OHMoleculeWeight ) * CO3MoleculeWeight * 1000  +      &
        npars%Acc_1_cm * npars%K_f_3_mol_cm2day * CO3MoleculeWeight * 1000   &
        -  npars%Acc_1_cm * npars%K_b_3_cm_mday * CO3 *                      &
           ( Ca / 1000 /CaMoleculeWeight )
!   WRITE(*,*) 'D_CO3 = ', D_CO3
  END FUNCTION

  REAL FUNCTION D_Ca( CO3, HCO3, H, Ca, CO2 ) 
   IMPLICIT NONE
   REAL :: CO3, HCO3, H, Ca, CO2
!   IF( CO3 .LT. 0.0) CO3=0.0
!   IF( HCO3 .LT. 0.0) HCO3=0.0
!   IF( H .LT. 0.0) H=0.0
!   IF( Ca .LT. 0.0) Ca=0.0
!   IF( CO2 .LT. 0.0) CO2=0.0
   D_Ca = npars%Acc_1_cm * npars%K_f_1_cm_day *                              &
         ( H / 1000 / HMoleculeWeight ) *  CaMoleculeWeight * 1000           &
        - npars%Acc_1_cm * npars%K_b_1_cm_mday *                             &
          ( HCO3 / 1000 / HCO3MoleculeWeight ) * Ca                          &
         + npars%Acc_1_cm * npars%K_f_2_cm_day *                             &
          ( CO2 / 1000 / CO2MoleculeWeight ) *  CaMoleculeWeight * 1000 -    &
           npars%Acc_1_cm * npars%K_b_2_cm_m2day *                           &
           ( HCO3 / 1000 / HCO3MoleculeWeight ) ** 2 * Ca  +                 &
         npars%Acc_1_cm * npars%K_f_3_mol_cm2day * CaMoleculeWeight * 1000   &
      - npars%Acc_1_cm * npars%K_b_3_cm_mday *                               &
        ( CO3 / 1000 / CO3MoleculeWeight ) * Ca
!   WRITE(*,*) 'D_Ca = ', D_Ca
  END FUNCTION

  REAL FUNCTION D_N2( CH2O, O2, NO3, X2 )
   IMPLICIT NONE
   REAL :: CH2O, O2, NO3, X2

!   IF( CH2O .LT. 0.0) CH2O=0.0
!   IF( O2 .LT. 0.0) O2=0.0
!   IF( NO3 .LT. 0.0) NO3=0.0
!   IF( X2 .LT. 0.0) X2=0.0

   D_N2 =  2 * N2MoleculeWeight / ( 5 * CH2OMoleculeWeight ) *            &
                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                    &
                     CH2O, NO3, O2,                                       &
                     npars%NO3HalfSaturationConst_mg_l,                   &
                     npars%CH2OHalfSaturationConst_mg_l,                  &
                     npars%O2InhibitionCoef_mg_l )
!   WRITE(*,*) 'D_N2 = ', D_N2
  END FUNCTION

  REAL FUNCTION ReactRateAtCell( conc, xloc, yloc, zloc, specName  )
   IMPLICIT NONE
   INTEGER :: xloc, yloc, zloc
   REAL*4, DIMENSION(:,:,:,:) :: conc
   CHARACTER(10) :: specName

   IF( TRIM( specName ) .eq. 'CH2O' ) THEN
       ReactRateAtCell = D_CH2O(                                              &
                        conc( SpeciesNameToInd('CH2O'),  & 
                              xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(               &
                                   'CH2O') ),              &
                        conc( SpeciesNameToInd( 'O2'),    &
                              xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(               &
                                  'O2') ),                &
                        conc( SpeciesNameToInd('NO3'),   &
                              xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(               &
                                 'NO3') ),               &
                              DenitrifierConc( xloc, yloc, zloc ) )

     ELSE IF( TRIM(specName) .eq. 'NH4' ) THEN

       ReactRateAtCell = D_NH4(                                            &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc )    +                           &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4') ),            &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2') ),             &
                    NitrifierConc( xloc, yloc, zloc ) )

     ELSE IF( TRIM(specName) .eq. 'NO3' ) THEN
       ReactRateAtCell = D_NO3(                                            &
                    conc( SpeciesNameToInd( 'NO3'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NO3') ),            &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4') ),            &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2') ),             &
                    NitrifierConc( xloc, yloc, zloc ),                     &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O') )            &
                         )
                         
     ELSE IF( TRIM(specName) .eq. 'O2' ) THEN
       ReactRateAtCell = D_O2(                                             &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4') ),            &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2') ),             &
                    NitrifierConc( xloc, yloc, zloc ),                     &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O') )            &
                         )

     ELSE IF( TRIM(specName) .eq. 'CO2' ) THEN
       ReactRateAtCell =  D_CO2(                                           &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2') ),             &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'NO3'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NO3') ),            &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O') ),           &
                    conc( SpeciesNameToInd( 'CO2'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO2') ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3') ),           &
                    conc( SpeciesNameToInd( 'H'),      &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'H') ),              &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH') ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2') )              &
                         )
     ELSE IF( TRIM(specName) .eq. 'HCO3' ) THEN
       ReactRateAtCell =  D_HCO3(                                          &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2') ),             &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'NO3'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NO3') ),            &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O') ),           &
                    conc( SpeciesNameToInd( 'CO2'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO2') ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3') ),           &
                    conc( SpeciesNameToInd( 'H'),      &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'H') ),              &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH') ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2') ),             &
                    conc( SpeciesNameToInd( 'CO3'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO3') )             &
                         )

     ELSE IF( TRIM(specName) .eq. 'H' ) THEN
       ReactRateAtCell =  D_H(                                             &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2') ),             &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4') ),            &
                    NitrifierConc( xloc, yloc, zloc ),                     &
                    conc( SpeciesNameToInd( 'CO2'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO2') ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3') ),           &
                    conc( SpeciesNameToInd( 'H'),      &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'H') ),              &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH') ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2') )              &
                         )
       
     ELSE IF( TRIM(specName) .eq. 'CO3' ) THEN
       ReactRateAtCell =  D_CO3(                                           &
                    conc( SpeciesNameToInd( 'CO3'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO3') ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3') ),           &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH') ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2') )              &
                         )

     ELSE IF( TRIM(specName) .eq. 'Ca2' ) THEN
       ReactRateAtCell =  D_Ca(                                          &
                    conc( SpeciesNameToInd( 'CO3'),  &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO3') ),          &
                    conc( SpeciesNameToInd( 'HCO3'), &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'HCO3') ),         &
                    conc( SpeciesNameToInd( 'H'),    &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'H') ),            &
                    conc( SpeciesNameToInd( 'Ca2'),   &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'Ca2') ),           &
                    conc( SpeciesNameToInd( 'CO2'),  &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO2') )           &
                         )

     ELSE IF( TRIM(specName) .eq. 'OH' ) THEN
       ReactRateAtCell =  D_OH(                                          &
                    conc( SpeciesNameToInd( 'H'),    &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'H') ),            &
                    conc( SpeciesNameToInd( 'OH'),   &
                         xloc, yloc, zloc )    +                         &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'OH') ),           &
                    conc( SpeciesNameToInd( 'CO2'),  &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO2') ),           &
                    conc( SpeciesNameToInd( 'HCO3'), &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'HCO3') ),         &
                    conc( SpeciesNameToInd( 'CO3'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO3') )                               &
                         )
     ELSE IF( TRIM(specName) .eq. 'N2' ) THEN
       ReactRateAtCell =  D_N2(                                          &
                    conc( SpeciesNameToInd( 'CH2O'), &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CH2O') ),         &
                    conc( SpeciesNameToInd( 'O2'),   &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'O2') ),           &
                    conc( SpeciesNameToInd( 'NO3'),  &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'NO3') ),          &
                    DenitrifierConc( xloc, yloc, zloc ) )
     ELSE 
       WRITE(*,*) 'ERROR: WRONG Species Name ', specName
       STOP
     END IF
  END FUNCTION

      SUBROUTINE VODEReactionRatesAtCell( conc, xloc, yloc, zloc, dt, rate  )

      INTEGER NEQ, ITASK, IOPT, LRW, LIW, MF, IOUT, ISTATE, ITOL, IPAR, &
              IWORK, I
      REAL ATOL, RPAR, RTOL, RWORK, T, TOUT, Y, dt
      DIMENSION Y(13), ATOL(13), RWORK(477), IWORK(43)
      REAL, DIMENSION(:) :: rate 
      INTEGER :: xloc, yloc, zloc
      REAL*4, DIMENSION(:,:,:,:) :: conc

      NEQ = TotalNumberOfSpecies + 2
      Y(1) = conc( SpeciesNameToInd( 'NH4'), xloc, yloc, zloc ) +        &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'NH4') )     ! NH4 
      Y(2) = conc( SpeciesNameToInd( 'NO3'), xloc, yloc, zloc )  +       &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'NO3') )
      Y(3) = conc( SpeciesNameToInd( 'CH2O'), xloc, yloc, zloc ) +       &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CH2O') )
      Y(4) = conc( SpeciesNameToInd( 'O2'), xloc, yloc, zloc )    +      &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'O2') )
      Y(5) = conc( SpeciesNameToInd( 'CO2'),  xloc, yloc, zloc )  +      &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO2') )
      Y(6) = conc( SpeciesNameToInd( 'HCO3'),  xloc, yloc, zloc ) +      &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'HCO3') )
      Y(7) = conc( SpeciesNameToInd( 'H'),  xloc, yloc, zloc )    +      &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'H') )
      Y(8) = conc( SpeciesNameToInd( 'CO3'),  xloc, yloc, zloc ) +       &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO3') )
      Y(9) = conc( SpeciesNameToInd( 'Ca2'),  xloc, yloc, zloc ) +       &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'Ca2') )
      Y(10) = conc( SpeciesNameToInd( 'N2'),  xloc, yloc, zloc ) +       &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'N2') )
      Y(11) = conc( SpeciesNameToInd( 'OH'),  xloc, yloc, zloc ) +       &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'OH') )
      Y(12) = NitrifierConc( xloc, yloc, zloc )      ! X1
      Y(13) = DenitrifierConc( xloc, yloc, zloc )    ! X2

      DO I = 1, TotalNumberOfSpecies + 2 
       rate(I) = Y(I)
      END DO

      T = 0.0D0
      TOUT = dt
      ITOL = 2
      RTOL = 1.D-4
      ATOL(1) = 1.D-6
      ATOL(2) = 1.D-6
      ATOL(3) = 1.D-6
      ATOL(4) = 1.D-6
      ATOL(5) = 1.D-6
      ATOL(6) = 1.D-6
      ATOL(7) = 1.D-6
      ATOL(8) = 1.D-6
      ATOL(9) = 1.D-6
      ATOL(10) = 1.D-6
      ATOL(11) = 1.D-6
      ATOL(12) = 1.D-6
      ATOL(13) = 1.D-6
      ITASK = 1
      ISTATE = 1
      IOPT = 1
      LRW = 477  ! 20 + NYH * (MAXORD + 1) + 3 * NEQ + LWM
                ! NYH = NEQ = 13
                ! MAXORD = 5 ( METH = 2 )
                ! NEQ = 13
                ! LWM = 2 * NEQ**2  + 2 ( MITER =1 and MF .gt. 0 )
                ! LRW = 20 + 13 * ( 5 + 1 ) + 3 * 13 + 2 * 13 ** 2 + 2
                !     = 477
      LIW = 43  ! 30 + NEQ
      MF = 21
      IWORK(5) = 0 
      IWORK(6) = 1000
      IWORK(7) = 0
      IWORK(8) = 0
      IWORK(9) = 0
      IWORK(10) = 0
      RWORK(5) = 0.0
      RWORK(6) = 0.0
      RWORK(7) = 0.0
      RWORK(8) = 0.0
      RWORK(9) = 0.0
      RWORK(10) = 0.0
      CALL SVODE(FREACT,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,   &
                 IOPT,RWORK,LRW,IWORK,LIW,JACRACT,MF,RPAR,IPAR)
!      WRITE(6,20)T,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9), &
!                     Y(10),Y(11), Y(12), Y(13)
!  20  FORMAT(' At t =',D12.4,'   y =',13D14.6)
      IF (ISTATE .LT. 0) GO TO 80
!      WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),         &
!                 IWORK(20),IWORK(21),IWORK(22)
!  60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,                 &
!            '   No. J-s =',I4,'   No. LU-s =',I4/                 & 
!            '  No. nonlinear iterations =',I4/                    &
!            '  No. nonlinear convergence failures =',I4/          &
!            '  No. error test failures =',I4/)
      DO I = 1, TotalNumberOfSpecies + 2 
       rate(I) = Y(I) - rate(I)
      END DO
      RETURN
  80  WRITE(6,90)ISTATE
  90  FORMAT(///' Error halt: ISTATE =',I3)
      STOP
      END SUBROUTINE

      SUBROUTINE STIFFEX()
      INTEGER NEQ, ITASK, IOPT, LRW, LIW, MF, IOUT, ISTATE, ITOL, IPAR, &
              IWORK
      REAL ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
      DIMENSION Y(3), ATOL(3), RWORK(67), IWORK(33)
      npars%MYNPARS = MYPAR

      NEQ = 3
      Y(1) = 1.0D0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      T = 0.0D0
      TOUT = 0.4D0
      ITOL = 2
      RTOL = 1.D-4
      ATOL(1) = 1.D-8
      ATOL(2) = 1.D-14
      ATOL(3) = 1.D-6
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LRW = 67
      LIW = 33
      MF = 21
      DO 40 IOUT = 1,12
        CALL SVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,   &
                 IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
        WRITE(6,20)T,Y(1),Y(2),Y(3)
  20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
        IF (ISTATE .LT. 0) GO TO 80
  40    TOUT = TOUT*10.
      WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),         &
                 IWORK(20),IWORK(21),IWORK(22)
  60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,                 &
            '   No. J-s =',I4,'   No. LU-s =',I4/                 & 
            '  No. nonlinear iterations =',I4/                    &
            '  No. nonlinear convergence failures =',I4/          &
            '  No. error test failures =',I4/)
      STOP
  80  WRITE(6,90)ISTATE
  90  FORMAT(///' Error halt: ISTATE =',I3)
      STOP
      END SUBROUTINE

      REAL  FUNCTION FUNC( A, B )
      REAL A, B
       FUNC = 1.D4 * A * B 
      END FUNCTION
 
      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
      INTEGER  NEQ, IPAR
      REAL  RPAR, T, Y, YDOT
      DIMENSION Y(NEQ), YDOT(NEQ)
!      YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
      YDOT(1) = -.04D0*Y(1) + FUNC( Y(2), Y(3) ) 
   
!      YDOT(3) = 3.D7*Y(2)*Y(2)
      YDOT(3) = ( npars%MYNPARS)*Y(2)*Y(2)
      YDOT(2) = -YDOT(1) - YDOT(3)
      RETURN
      END SUBROUTINE
 
      SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      INTEGER  NEQ, NRPD, IPAR, ML, MU
      REAL PD, RPAR, T, Y
      DIMENSION Y(NEQ), PD(NRPD,NEQ)
      PD(1,1) = -.04D0
      PD(1,2) = 1.D4*Y(3)
      PD(1,3) = 1.D4*Y(2)
      PD(2,1) = .04D0
      PD(2,3) = -PD(1,3)
      PD(3,2) = 6.D7*Y(2)
      PD(2,2) = -PD(1,2) - PD(3,2)
      RETURN
      END SUBROUTINE

      SUBROUTINE RACT()
      INTEGER NEQ, ITASK, IOPT, LRW, LIW, MF, IOUT, ISTATE, ITOL, IPAR, &
              IWORK
      REAL ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
      DIMENSION Y(13), ATOL(13), RWORK(477), IWORK(43)
      npars%MYNPARS = MYPAR

      NEQ = 13
      Y(1) = 40                       ! NH4 
      Y(2) = 5.76                     ! NO3
      Y(3) = 82.0                     ! CH2O
      Y(4) = 6.0 + 6.0                ! O2
      Y(5) = 57.2  + 23.32            ! CO2
      Y(6) = 445.3 + 207.4            ! HCO3
      Y(7) = 6.8e-5  + 5.9E-5         ! H
      Y(8) = 0.24   + 0.132           ! CO3
      Y(9) = 37.2  + 68.0             ! Ca
      Y(10) = 0.0                     ! N2
      Y(11) =  1.122e-3 + 1.292e-3    ! OH 
      Y(12) = 0.256                   ! X1
      Y(13) = 0.256                   ! X2

      T = 0.0D0
      TOUT = 0.4D0
      ITOL = 2
      RTOL = 1.D-4
      ATOL(1) = 1.D-8
      ATOL(2) = 1.D-14
      ATOL(3) = 1.D-6
      ATOL(4) = 1.D-6
      ATOL(5) = 1.D-6
      ATOL(6) = 1.D-6
      ATOL(7) = 1.D-6
      ATOL(8) = 1.D-6
      ATOL(9) = 1.D-6
      ATOL(10) = 1.D-6
      ATOL(11) = 1.D-6
      ATOL(12) = 1.D-6
      ATOL(13) = 1.D-6
      ITASK = 1
      ISTATE = 1
      IOPT = 1
      LRW = 477  ! 20 + NYH * (MAXORD + 1) + 3 * NEQ + LWM
                ! NYH = NEQ = 13
                ! MAXORD = 5 ( METH = 2 )
                ! NEQ = 13
                ! LWM = 2 * NEQ**2  + 2 ( MITER =1 and MF .gt. 0 )
                ! LRW = 20 + 13 * ( 5 + 1 ) + 3 * 13 + 2 * 13 ** 2 + 2
                !     = 477
      LIW = 43  ! 30 + NEQ
      MF = 21
      IWORK(5) = 0 
      IWORK(6) = 1000
      IWORK(7) = 0
      IWORK(8) = 0
      IWORK(9) = 0
      IWORK(10) = 0
      RWORK(5) = 0.0
      RWORK(6) = 0.0
      RWORK(7) = 0.0
      RWORK(8) = 0.0
      RWORK(9) = 0.0
      RWORK(10) = 0.0
      DO 40 IOUT = 1,12
        CALL SVODE(FREACT,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,   &
                 IOPT,RWORK,LRW,IWORK,LIW,JACRACT,MF,RPAR,IPAR)
        WRITE(6,20)T,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9), &
                     Y(10),Y(11), Y(12), Y(13)
  20    FORMAT(' At t =',D12.4,'   y =',13D14.6)
        IF (ISTATE .LT. 0) GO TO 80
  40    TOUT = TOUT*10.
      WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),         &
                 IWORK(20),IWORK(21),IWORK(22)
  60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,                 &
            '   No. J-s =',I4,'   No. LU-s =',I4/                 & 
            '  No. nonlinear iterations =',I4/                    &
            '  No. nonlinear convergence failures =',I4/          &
            '  No. error test failures =',I4/)
      STOP
  80  WRITE(6,90)ISTATE
  90  FORMAT(///' Error halt: ISTATE =',I3)
      STOP
      END SUBROUTINE
      SUBROUTINE FREACT (NEQ, T, Y, YDOT, RPAR, IPAR)
      INTEGER  NEQ, IPAR
      REAL  RPAR, T, Y, YDOT
      DIMENSION Y(NEQ), YDOT(NEQ)
!     
!     Y(1) = NH4
!     Y(2) = NO3
!     Y(3) = CH2O
!     Y(4) = O2
!     Y(5) = CO2
!     Y(6) = HCO3
!     Y(7) = H
!     Y(8) = CO3
!     Y(9) = Ca2
!     Y(10) = N2
!     Y(11) = OH
!     Y(12) = X1
!     Y(13) = X2
!     
               ! D_NH4( NH4, O2, X1 )
      YDOT(1) = D_NH4( Y(1), Y(4), Y(12) )
               ! D_NO3( NO3, NH4, O2, X1, X2, CH2O ) 
      YDOT(2) = D_NO3( Y(2), Y(1), Y(4), Y(12), Y(13), Y(3) ) 
               ! D_CH2O( CH2O, O2, NO3, X2 )
      YDOT(3) = D_CH2O( Y(3), Y(4), Y(2), Y(13) )
               ! D_O2( NH4, O2, X1, X2, CH2O ) 
      YDOT(4) = D_O2( Y(1), Y(4), Y(12), Y(13), Y(3) ) 
              ! D_CO2( O2, X2, NO3, CH2O, CO2, HCO3, H, OH, Ca ) 
      YDOT(5) = D_CO2( Y(4), Y(13), Y(2), Y(3), Y(5), Y(6), Y(7), Y(11), Y(9) ) 
   !                O2, X2,   NO3, CH2O, CO2, HCO3, H, OH,   Ca, CO3  
  YDOT(6) = D_HCO3(Y(4),Y(13),Y(2),Y(3),Y(5),Y(6),Y(7),Y(11),Y(9),Y(8) ) 
   !             O2,    NH4, X1,    CO2,  HCO3, H,     OH,    Ca  
  YDOT(7) = D_H( Y(4), Y(1), Y(12), Y(5), Y(6), Y(7), Y(11), Y(9) ) 
   !               CO3,  HCO3, OH,    Ca  
  YDOT(8) = D_CO3( Y(8), Y(6), Y(11), Y(9) ) 
   !              CO3, HCO3,  H,    Ca,   CO2  
  YDOT(9) = D_Ca( Y(8), Y(6), Y(7), Y(9), Y(5) ) 
   !                CH2O, O2, NO3, X2 
  YDOT(10) = D_N2( Y(3), Y(4), Y(2), Y(13) )
   !               H,    OH,    CO2,  HCO3, CO3  
  YDOT(11) = D_OH( Y(7), Y(11), Y(5), Y(6), Y(8) ) 
   !                  D_NH4, X1, Y1, kd1 
  YDOT(12) = DeltaX1( YDOT(1), Y(12),                                      &
                         npars%Biomass1YieldCoeff_M_biomass_M_substrate,   &
                         npars%Biomass1SpecDecayConst_1_day )
   ! DeltaX2( D_CH2O, X2, Y2, kd2 )
  YDOT(13) = DeltaX2( YDOT(3), Y(13),                                         &
                         npars%Biomass2YieldCoeff_M_biomass_M_substrate,   &
                         npars%Biomass2SpecDecayConst_1_day )
      END SUBROUTINE

      SUBROUTINE JACRACT (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      INTEGER  NEQ, NRPD, IPAR, ML, MU
      REAL PD, RPAR, T, Y
      DIMENSION Y(NEQ), PD(NRPD,NEQ)
      REAL NH4, NO3, CH2O, O2, CO2, HCO3, H, CO3, Ca, N2, OH, X1, X2

     NH4 = Y(1)
     NO3 = Y(2)
     CH2O = Y(3)
     O2 = Y(4)
     CO2 = Y(5)
     HCO3 = Y(6)
     H = Y(7)
     CO3 = Y(8)
     Ca = Y(9)
     N2 = Y(10)
     OH = Y(11)
     X1 = Y(12)
     X2 = Y(13)

PD(1,1) = 1.0*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))-1.0&
*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
PD(1,2) = 0
PD(1,3) = 0
PD(1,4) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))-1.0&
*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
PD(1,5) = 0
PD(1,6) = 0
PD(1,7) = 0
PD(1,8) = 0
PD(1,9) = 0
PD(1,10) = 0
PD(1,11) = 0
PD(1,12) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)-1.&
0*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))
PD(1,13) = 0
PD(2,1) = 3.444444444444445*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.&
0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
PD(2,2) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
**2*(O2+1.0)*(X2+0.5))-8.266666666666666*CH2O*X2/((CH2O+10)*(NO3+&
0.5)*(O2+1.0)*(X2+0.5))
PD(2,3) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
.5)*(O2+1.0)*(X2+0.5))-8.266666666666666*NO3*X2/((CH2O+10)*(NO3+0&
.5)*(O2+1.0)*(X2+0.5))
PD(2,4) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
*(O2+1.0)**2*(X2+0.5))+3.444444444444445*NH4*X1/((NH4+0.1)*(O2+0.&
1)*(X1+1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(&
X1+1.0))
PD(2,5) = 0
PD(2,6) = 0
PD(2,7) = 0
PD(2,8) = 0
PD(2,9) = 0
PD(2,10) = 0
PD(2,11) = 0
PD(2,12) = 3.444444444444445*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+&
1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2&
)
PD(2,13) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
)*(O2+1.0)*(X2+0.5)**2)-8.266666666666666*CH2O*NO3/((CH2O+10)*(NO&
3+0.5)*(O2+1.0)*(X2+0.5))
PD(3,1) = 0
PD(3,2) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2+1.0)*(&
X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5))
PD(3,3) = -5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)&
)+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5))-5.0*&
O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+10)**2*&
(O2+0.1)*(X2+0.5))
PD(3,4) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)**2*(&
X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2&
/((CH2O+10)*(O2+0.1)**2*(X2+0.5))
PD(3,5) = 0
PD(3,6) = 0
PD(3,7) = 0
PD(3,8) = 0
PD(3,9) = 0
PD(3,10) = 0
PD(3,11) = 0
PD(3,12) = 0
PD(3,13) = -5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0&
.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO3*X2/((&
CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((CH2O+10&
)*(O2+0.1)*(X2+0.5)**2)
PD(4,1) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)&
*(X1+1.0))-3.555555555555555*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
PD(4,2) = 0
PD(4,3) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)**2*(O2+0.1&
)*(X2+0.5))-5.333333333333333*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))
PD(4,4) = -5.333333333333333*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2&
+0.5))+5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.&
5))-3.555555555555555*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))+3.5555&
55555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
PD(4,5) = 0
PD(4,6) = 0
PD(4,7) = 0
PD(4,8) = 0
PD(4,9) = 0
PD(4,10) = 0
PD(4,11) = 0
PD(4,12) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(&
X1+1.0)**2)-3.555555555555555*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0)&
)
PD(4,13) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)*&
(X2+0.5)**2)-5.333333333333333*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.&
5))
PD(5,1) = 0
PD(5,2) = 1.466666666666667*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2&
+1.0)*(X2+0.5))-1.466666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
)**2*(O2+1.0)*(X2+0.5))
PD(5,3) = 1.466666666666667*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+&
1.0)*(X2+0.5))-1.466666666666667*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
.5)*(O2+1.0)*(X2+0.5))+7.333333333333333*O2*X2/((CH2O+10)*(O2+0.1&
)*(X2+0.5))-7.333333333333333*CH2O*O2*X2/((CH2O+10)**2*(O2+0.1)*(&
X2+0.5))
PD(5,4) = -1.466666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
)*(O2+1.0)**2*(X2+0.5))+7.333333333333333*CH2O*X2/((CH2O+10)*(O2+&
0.1)*(X2+0.5))-7.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)**&
2*(X2+0.5))
PD(5,5) = -43200.0*OH-2.608416
PD(5,6) = 2.561246976619189e-5*Ca*HCO3+4917.901639344263*H+.&
1414754098360656*CO2Moleculeweight
PD(5,7) = 4917.901639344263*HCO3
PD(5,8) = 0
PD(5,9) = 1.2806234883095945e-5*HCO3**2
PD(5,10) = 0
PD(5,11) = -43200.0*CO2
PD(5,12) = 0
PD(5,13) = 1.466666666666667*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(&
O2+1.0)*(X2+0.5))+7.333333333333333*CH2O*O2/((CH2O+10)*(O2+0.1)*(&
X2+0.5))-1.466666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1&
.0)*(X2+0.5)**2)-7.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)&
*(X2+0.5)**2)
PD(6,1) = 0
PD(6,2) = 8.133333333333333*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2&
+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
)**2*(O2+1.0)*(X2+0.5))
PD(6,3) = 8.133333333333333*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+&
1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
.5)*(O2+1.0)*(X2+0.5))
PD(6,4) = -8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
)*(O2+1.0)**2*(X2+0.5))
PD(6,5) = 59890.90909090908*OH+3.616213090909091
PD(6,6) = -5421.176470588235*OH-3.550819672131148e-5*Ca*HCO3&
-6818.0*H-7.3055e-6*Ca-8.63
PD(6,7) = 1782.542-6818.0*HCO3
PD(6,8) = 10980.0
PD(6,9) = -1.775409836065574e-5*HCO3**2-7.3055e-6*HCO3
PD(6,10) = 0
PD(6,11) = 59890.90909090908*CO2-5421.176470588235*HCO3
PD(6,12) = 0
PD(6,13) = 8.133333333333333*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(&
O2+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0&
.5)*(O2+1.0)*(X2+0.5)**2)
PD(7,1) = .1111111111111111*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.&
0))-.1111111111111111*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
PD(7,2) = 0
PD(7,3) = 0
PD(7,4) = .1111111111111111*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1&
.0))-.1111111111111111*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
PD(7,5) = .05890909090909091
PD(7,6) = 1.1976229508196721e-7*Ca-111.7704918032787*H
PD(7,7) = -1447058.823529412*OH-111.7704918032787*HCO3-29.22&
2
PD(7,8) = 0
PD(7,9) = 1.1976229508196721e-7*HCO3
PD(7,10) = 0
PD(7,11) = -1447058.823529412*H
PD(7,12) = .1111111111111111*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+&
1.0))-.1111111111111111*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2&
)
PD(7,13) = 0
PD(8,1) = 0
PD(8,2) = 0
PD(8,3) = 0
PD(8,4) = 0
PD(8,5) = 0
PD(8,6) = 5332.304725168757*OH
PD(8,7) = 0
PD(8,8) = -0.14364*Ca-10800.0
PD(8,9) = -0.14364*CO3
PD(8,10) = 0
PD(8,11) = 5332.304725168757*HCO3
PD(8,12) = 0
PD(8,13) = 0
PD(9,1) = 0
PD(9,2) = 0
PD(9,3) = 0
PD(9,4) = 0
PD(9,5) = .01492363636363636
PD(9,6) = -2.3284063423810806e-5*Ca*HCO3-4.790491803278688e-&
6*Ca
PD(9,7) = 1168.88
PD(9,8) = -0.09576*Ca
PD(9,9) = -1.1642031711905403e-5*HCO3**2-4.790491803278688e-&
6*HCO3-0.09576*CO3
PD(9,10) = 0
PD(9,11) = 0
PD(9,12) = 0
PD(9,13) = 0
PD(10,1) = 0
PD(10,2) = 1.866666666666667*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O&
2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
5)**2*(O2+1.0)*(X2+0.5))
PD(10,3) = 1.866666666666667*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2&
+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)**2*(NO3+&
0.5)*(O2+1.0)*(X2+0.5))
PD(10,4) = -1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
5)*(O2+1.0)**2*(X2+0.5))
PD(10,5) = 0
PD(10,6) = 0
PD(10,7) = 0
PD(10,8) = 0
PD(10,9) = 0
PD(10,10) = 0
PD(10,11) = 0
PD(10,12) = 0
PD(10,13) = 1.866666666666667*CH2O*NO3/((CH2O+10)*(NO3+0.5)*&
(O2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+&
0.5)*(O2+1.0)*(X2+0.5)**2)
PD(11,1) = 0
PD(11,2) = 0
PD(11,3) = 0
PD(11,4) = 0
PD(11,5) = -16690.90909090909*OH
PD(11,6) = .1414754098360656*OHMoleculeweight-1510.819672131&
148*OH
PD(11,7) = -2.46e+7*OH
PD(11,8) = 3060.0
PD(11,9) = 0
PD(11,10) = 0
PD(11,11) = -1510.819672131148*HCO3-2.46e+7*H-16690.90909090&
909*CO2
PD(11,12) = 0
PD(11,13) = 0
PD(12,1) = 0.17*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH4&
*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
PD(12,2) = 0
PD(12,3) = 0
PD(12,4) = 0.17*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH&
4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
PD(12,5) = 0
PD(12,6) = 0
PD(12,7) = 0
PD(12,8) = 0
PD(12,9) = 0
PD(12,10) = 0
PD(12,11) = 0
PD(12,12) = 0.17*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*N&
H4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)
PD(12,13) = 0
PD(13,1) = 0
PD(13,2) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2&
+1.0)*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5&
)))
PD(13,3) = -0.5*(-5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(&
X2+0.5))+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5&
))-5.0*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+&
10)**2*(O2+0.1)*(X2+0.5)))
PD(13,4) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.&
0)**2*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2&
O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.5)))
PD(13,5) = 0
PD(13,6) = 0
PD(13,7) = 0
PD(13,8) = 0
PD(13,9) = 0
PD(13,10) = 0
PD(13,11) = 0
PD(13,12) = 0
PD(13,13) = -0.5*(-5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0&
)*(X2+0.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO&
3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((&
CH2O+10)*(O2+0.1)*(X2+0.5)**2))

!     Y(1) = NH4
!     Y(2) = NO3
!     Y(3) = CH2O
!     Y(4) = O2
!     Y(5) = CO2
!     Y(6) = HCO3
!     Y(7) = H
!     Y(8) = CO3
!     Y(9) = Ca2
!     Y(10) = N2
!     Y(11) = OH
!     Y(12) = X1
!     Y(13) = X2
      END SUBROUTINE

END MODULE NTransport
