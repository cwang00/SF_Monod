

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
     REAL :: Soil_bulk_den_g_cm3
!     REAL, DIMENSION(TotalNumberOfSpecies)::backgroundConc
     REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: backgroundConc
     REAL :: MYNPARS
   END TYPE NTransParameters

   REAL, DIMENSION(:,:,:), ALLOCATABLE :: NitrifierConc
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: DenitrifierConc

   TYPE( NTransParameters ) :: npars

   CONTAINS

     SUBROUTINE ReadNTransPars( nparfile, nx, ny, nz )
      IMPLICIT NONE

      CHARACTER(LEN=100), INTENT(IN) :: nparfile
      CHARACTER(LEN=128) :: bgConcType
      CHARACTER(LEN=500) :: bgConcPFBFile
      INTEGER :: I, J, K, L, M, N, idx, nx, ny, nz
      REAL*8 :: bgConcValue, dx, dy, dz
      REAL*8 :: temp(nx,ny,nz)
      INTEGER*4 :: nx1, ny1, nz1
interface

    SUBROUTINE pf_read(x,filename,nx,ny,nz,dx2,dy2,dz2)
    real*8  :: x(:,:,:)
    character*500 :: filename
    integer*4 :: nx
    integer*4 :: ny
    integer*4 :: nz
    real*8  :: dx2
    real*8  :: dy2
    real*8  :: dz2
    END SUBROUTINE pf_read
end interface

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
     READ(98,*) npars%Soil_bulk_den_g_cm3

     ALLOCATE( npars%backgroundConc(TotalNumberOfSpecies,nx, ny, nz ) )

     DO I = 1, TotalNumberOfSpecies

        READ(98,*) bgConcType
        IF ( TRIM( bgConcType  ) .eq. 'const' ) THEN
           READ(98, *) bgConcValue
           npars%backgroundConc(I, :, :, : ) = bgConcValue
        ELSE IF ( TRIM( bgConcType ) .eq. 'PFBFile' ) THEN
           READ(98, *) bgConcPFBFile
           CALL PF_READ( temp, bgConcPFBFile, nx1, ny1, nz1, &
                           dx, dy, dz )
            DO L = 1, nx
             DO M = 1, ny
              DO N = 1, nz
               npars%backgroundConc(I,L,M,N) = temp(L,M,N)
              END DO
             END DO
           END DO
        ELSE
                WRITE(*,*) 'ERROR: UNKNOWN Background concentration type ', &
                            bgConcType 
                STOP
        ENDIF
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
           -4 * NO3MoleculeWeight / ( 5 * CH2OMoleculeWeight ) *             &
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

     D_CO2 = 1.0 / ( CH2OMoleculeWeight * 1000.0 ) *                         &
                R1( npars%OrgCarbonOxMaxRate_1_day,  X2, CH2O, O2,           &
                    npars%CH2OHalfSaturationConst_mg_l,                      &
                 npars%O2HalfSaturationConst_mg_l ) +                        &
                1.0 / ( 5 * CH2OMoleculeWeight * 1000.0 ) *                  &
                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                       &
                     CH2O, NO3, O2,                                          &
                     npars%NO3HalfSaturationConst_mg_l,                      &
                     npars%CH2OHalfSaturationConst_mg_l,                     &
                     npars%O2InhibitionCoef_mg_l )                           &
             - npars%K_f_CO2_1_day * CO2 +                                   & 
        npars%K_b_CO2_1_mday * ( HCO3 ) * ( H )  -                           &
               npars%K_f_H_1_mday * CO2 * ( OH ) +                           &
          npars%K_b_H_1_day * HCO3                                           &
               -npars%Acc_1_cm * npars%K_f_2_cm_day * CO2 +                  &
                npars%Acc_1_cm * npars%K_b_2_cm_m2day *                      &
                ( HCO3 )**2 * ( Ca ) 
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

!   D_HCO3 = ( 4 * HCO3MoleculeWeight ) / ( 5 *  CH2OMoleculeWeight ) *       &
   D_HCO3 =  4  / ( 5 *  CH2OMoleculeWeight * 1000.0 ) *       &
                 R3( npars%DenitOxMaxRate_1_day,  X2 ,                       &
                     CH2O, NO3, O2,                                          &
                     npars%NO3HalfSaturationConst_mg_l,                      &
                     npars%CH2OHalfSaturationConst_mg_l,                     &
                     npars%O2InhibitionCoef_mg_l )                           &
      + npars%K_f_CO2_1_day * CO2 -                                          &
        npars%K_b_CO2_1_mday * HCO3 * ( H )                                  &
      + npars%K_f_H_1_mday * ( CO2 )                                         &
            * ( OH  )  -                                                     &
         npars%K_b_H_1_day * HCO3                                            &
       + npars%K_f_HCO3_1_day * CO3                                          &
      - npars%K_b_HCO3_1_mday * HCO3 * ( OH )                                &
         + npars%Acc_1_cm *                                                  &
                                   npars%K_f_1_cm_day * H                    &
         - npars%Acc_1_cm * npars%K_b_1_cm_mday *                            &
             HCO3 * ( Ca )                                                   &
         + npars%Acc_1_cm *                                                  &
                     npars%K_f_2_cm_day *  CO2                               &
         - npars%Acc_1_cm * npars%K_b_2_cm_m2day                             &
            * ( HCO3 ) ** 2  *                                               &
              ( Ca )
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
   D_H = 2.0 / ( NH4MoleculeWeight * 1000.0) *                               &
             R2( npars%NitOxMaxRate_1_day,  X1, NH4, O2,                     &
                      npars%NH4HalfSaturationConst_mg_l,                     &
                      npars%O2HalfSaturationConst_mg_l )  * 2                &
   + npars%K_f_CO2_1_day * CO2         &
  - npars%K_b_CO2_1_mday * ( HCO3) *  H         &
         - npars%Acc_1_cm * npars%K_f_1_cm_day * H +                         &
         npars%Acc_1_cm * npars%K_b_1_cm_mday *                              &
            ( HCO3 ) *                           &
            ( Ca )                               &
         + npars%K_f_H2O_1_day  -                    &
           npars%K_b_H2O_1_mday * ( OH ) * H 
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
   D_OH = npars%K_f_H2O_1_day  -                    &
          npars%K_b_H2O_1_mday * OH * ( H )         &
          -npars%K_f_H_1_mday * ( CO2 ) * OH      &
          +npars%K_b_H_1_day * HCO3   &
       + npars%K_f_HCO3_1_day * CO3   &
       - npars%K_b_HCO3_1_mday * ( HCO3 ) * OH

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
       + npars%K_b_HCO3_1_mday * ( HCO3 ) *       &
         ( OH ) +      &
        npars%Acc_1_cm * npars%K_f_3_mol_cm2day   &
        -  npars%Acc_1_cm * npars%K_b_3_cm_mday * CO3 *                      &
           ( Ca )
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
         ( H )            &
        - npars%Acc_1_cm * npars%K_b_1_cm_mday *                             &
          ( HCO3 ) * Ca                          &
         + npars%Acc_1_cm * npars%K_f_2_cm_day *                             &
          ( CO2 ) -    &
           npars%Acc_1_cm * npars%K_b_2_cm_m2day *                           &
           ( HCO3 ) ** 2 * Ca  +                 &
         npars%Acc_1_cm * npars%K_f_3_mol_cm2day    &
      - npars%Acc_1_cm * npars%K_b_3_cm_mday *                               &
        ( CO3 ) * Ca
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
                                   'CH2O'), xloc, yloc, zloc ),              &
                        conc( SpeciesNameToInd( 'O2'),    &
                              xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(               &
                                  'O2'), xloc, yloc, zloc ),                &
                        conc( SpeciesNameToInd('NO3'),   &
                              xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(               &
                                 'NO3'), xloc, yloc, zloc ),               &
                              DenitrifierConc( xloc, yloc, zloc ) )

     ELSE IF( TRIM(specName) .eq. 'NH4' ) THEN

       ReactRateAtCell = D_NH4(                                            &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc )    +                           &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2'), xloc, yloc, zloc ),             &
                    NitrifierConc( xloc, yloc, zloc ) )

     ELSE IF( TRIM(specName) .eq. 'NO3' ) THEN
       ReactRateAtCell = D_NO3(                                            &
                    conc( SpeciesNameToInd( 'NO3'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NO3'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2'), xloc, yloc, zloc ),             &
                    NitrifierConc( xloc, yloc, zloc ),                     &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O'), xloc, yloc, zloc )            &
                         )
                         
     ELSE IF( TRIM(specName) .eq. 'O2' ) THEN
       ReactRateAtCell = D_O2(                                             &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2'), xloc, yloc, zloc ),             &
                    NitrifierConc( xloc, yloc, zloc ),                     &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O'), xloc, yloc, zloc )            &
                         )

     ELSE IF( TRIM(specName) .eq. 'CO2' ) THEN
       ReactRateAtCell =  D_CO2(                                           &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2'), xloc, yloc, zloc ),             &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'NO3'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NO3'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'CO2'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO2'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'H'),      &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'H'), xloc, yloc, zloc ),              &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH'), xloc, yloc, zloc ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2'), xloc, yloc, zloc )              &
                         )
     ELSE IF( TRIM(specName) .eq. 'HCO3' ) THEN
       ReactRateAtCell =  D_HCO3(                                          &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2'), xloc, yloc, zloc ),             &
                    DenitrifierConc( xloc, yloc, zloc ),                   &
                    conc( SpeciesNameToInd( 'NO3'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NO3'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'CH2O'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CH2O'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'CO2'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO2'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'H'),      &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'H'), xloc, yloc, zloc ),              &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH'), xloc, yloc, zloc ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2'), xloc, yloc, zloc ),             &
                    conc( SpeciesNameToInd( 'CO3'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO3'), xloc, yloc, zloc )             &
                         )

     ELSE IF( TRIM(specName) .eq. 'H' ) THEN
       ReactRateAtCell =  D_H(                                             &
                    conc( SpeciesNameToInd( 'O2'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'O2'), xloc, yloc, zloc ),             &
                    conc( SpeciesNameToInd( 'NH4'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'NH4'), xloc, yloc, zloc ),            &
                    NitrifierConc( xloc, yloc, zloc ),                     &
                    conc( SpeciesNameToInd( 'CO2'),    &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO2'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'H'),      &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'H'), xloc, yloc, zloc ),              &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH'), xloc, yloc, zloc ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2'), xloc, yloc, zloc )              &
                         )
       
     ELSE IF( TRIM(specName) .eq. 'CO3' ) THEN
       ReactRateAtCell =  D_CO3(                                           &
                    conc( SpeciesNameToInd( 'CO3'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO3'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'HCO3'),   &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'HCO3'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'OH'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'OH'), xloc, yloc, zloc ),             &
                    conc( SpeciesNameToInd( 'Ca2'),     &
                         xloc, yloc, zloc ) +                              &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'Ca2'), xloc, yloc, zloc )              &
                         )

     ELSE IF( TRIM(specName) .eq. 'Ca2' ) THEN
       ReactRateAtCell =  D_Ca(                                          &
                    conc( SpeciesNameToInd( 'CO3'),  &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO3'), xloc, yloc, zloc ),          &
                    conc( SpeciesNameToInd( 'HCO3'), &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'HCO3'), xloc, yloc, zloc ),         &
                    conc( SpeciesNameToInd( 'H'),    &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'H'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'Ca2'),   &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'Ca2'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'CO2'),  &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO2'), xloc, yloc, zloc )           &
                         )

     ELSE IF( TRIM(specName) .eq. 'OH' ) THEN
       ReactRateAtCell =  D_OH(                                          &
                    conc( SpeciesNameToInd( 'H'),    &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'H'), xloc, yloc, zloc ),            &
                    conc( SpeciesNameToInd( 'OH'),   &
                         xloc, yloc, zloc )    +                         &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'OH'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'CO2'),  &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CO2'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'HCO3'), &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'HCO3'), xloc, yloc, zloc ),         &
                    conc( SpeciesNameToInd( 'CO3'),    &
                         xloc, yloc, zloc )  +                             &
                        npars%backgroundConc( SpeciesNameToInd(            &
                                  'CO3'), xloc, yloc, zloc )               &
                         )
     ELSE IF( TRIM(specName) .eq. 'N2' ) THEN
       ReactRateAtCell =  D_N2(                                          &
                    conc( SpeciesNameToInd( 'CH2O'), &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'CH2O'), xloc, yloc, zloc ),         &
                    conc( SpeciesNameToInd( 'O2'),   &
                         xloc, yloc, zloc ) +                            &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'O2'), xloc, yloc, zloc ),           &
                    conc( SpeciesNameToInd( 'NO3'),  &
                         xloc, yloc, zloc )  +                           &
                        npars%backgroundConc( SpeciesNameToInd(          &
                                  'NO3'), xloc, yloc, zloc ),          &
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

!      WRITE(*,*) 'LOCX = ', xloc, ' LOCY = ', yloc, ' LOCZ = ', zloc
      NEQ = TotalNumberOfSpecies + 2
      Y(1) = conc( SpeciesNameToInd( 'NH4'), xloc, yloc, zloc ) ! +        &
                  !       npars%backgroundConc( SpeciesNameToInd(          &
!                                  'NH4') )     ! NH4 
      Y(2) = conc( SpeciesNameToInd( 'NO3'), xloc, yloc, zloc ) !  +       &
                   !      npars%backgroundConc( SpeciesNameToInd(          &
                    !               'NO3') )
      Y(3) = conc( SpeciesNameToInd( 'CH2O'), xloc, yloc, zloc ) ! +       &
                   !    npars%backgroundConc( SpeciesNameToInd(          &
                   !              'CH2O') )
      Y(4) = conc( SpeciesNameToInd( 'O2'), xloc, yloc, zloc )  ! +      &
                   !    npars%backgroundConc( SpeciesNameToInd(          &
                   !              'O2') )
      Y(5) =  conc( SpeciesNameToInd( 'CO2'),  xloc, yloc, zloc ) 
!      Y(5) = ( conc( SpeciesNameToInd( 'CO2'),  xloc, yloc, zloc )  +    &
!               npars%backgroundConc( SpeciesNameToInd(  'CO2') ) ) /     &
!               ( CO2MoleculeWeight * 1000.0 ) ! mg/l to M
      Y(6) =  conc( SpeciesNameToInd( 'HCO3'),  xloc, yloc, zloc )
!      Y(6) = ( conc( SpeciesNameToInd( 'HCO3'),  xloc, yloc, zloc ) +    &
!                  npars%backgroundConc( SpeciesNameToInd( 'HCO3') ) ) /  &
!                  ( HCO3MoleculeWeight * 1000.0 ) ! mg/l to M
      Y(7) =  conc( SpeciesNameToInd( 'H'),  xloc, yloc, zloc )
!      Y(7) = ( conc( SpeciesNameToInd( 'H'),  xloc, yloc, zloc )    +    &
!                    npars%backgroundConc( SpeciesNameToInd( 'H') ) ) /   &
!                  ( HMoleculeWeight * 1000.0 ) ! mg/l to M
      Y(8) =  conc( SpeciesNameToInd( 'CO3'),  xloc, yloc, zloc ) 
!      Y(8) = ( conc( SpeciesNameToInd( 'CO3'),  xloc, yloc, zloc ) +     &
!               npars%backgroundConc( SpeciesNameToInd( 'CO3') ) )   /    &
!                  ( CO3MoleculeWeight * 1000.0 ) ! mg/l to M
      Y(9) =  conc( SpeciesNameToInd( 'Ca2'),  xloc, yloc, zloc )
!      Y(9) = ( conc( SpeciesNameToInd( 'Ca2'),  xloc, yloc, zloc ) +     &
!               npars%backgroundConc( SpeciesNameToInd(  'Ca2') ) )  /    &
!                  ( CaMoleculeWeight * 1000.0 ) ! mg/l to M
      Y(10) = conc( SpeciesNameToInd( 'N2'),  xloc, yloc, zloc ) ! +       &
!                   npars%backgroundConc( SpeciesNameToInd( 'N2') )
      Y(11) =  conc( SpeciesNameToInd( 'OH'),  xloc, yloc, zloc ) 
!      Y(11) = ( conc( SpeciesNameToInd( 'OH'),  xloc, yloc, zloc ) +     &
!                  npars%backgroundConc( SpeciesNameToInd( 'OH') ) ) /    &
!                  ( OHMoleculeWeight * 1000.0 ) ! mg/l to M
      Y(12) = NitrifierConc( xloc, yloc, zloc )      ! X1
      Y(13) = DenitrifierConc( xloc, yloc, zloc )    ! X2

      DO I = 1, TotalNumberOfSpecies + 2 
       rate(I) = Y(I)
      END DO

      T = 0.0D0
      TOUT = dt
      ITOL = 2
      !RTOL = 1.D-4
      RTOL = 0.0
      !ATOL(1) = 1.D-4
      !ATOL(2) = 1.D-4
      !ATOL(3) = 1.D-4
      !ATOL(4) = 1.D-4
      !ATOL(5) = 1.D-6
      !ATOL(6) = 1.D-6
      !ATOL(7) = 1.D-10
      !ATOL(8) = 1.D-8
      !ATOL(9) = 1.D-7
      !ATOL(10) = 1.D-8
      !ATOL(11) = 1.D-10
      !ATOL(12) = 1.D-6
      !ATOL(13) = 1.D-6
      ATOL(1) = 1.D-3
      ATOL(2) = 1.D-3
      ATOL(3) = 1.D-3
      ATOL(4) = 1.D-3
      ATOL(5) = 1.D-5
      ATOL(6) = 1.D-5
      ATOL(7) = 1.D-8
      ATOL(8) = 1.D-7
      ATOL(9) = 1.D-4
      ATOL(10) = 1.D-4
      ATOL(11) = 1.D-7
      ATOL(12) = 1.D-4
      ATOL(13) = 1.D-4

      ITASK = 1
      ISTATE = 1
!      IOPT = 1
      IOPT = 0
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
!      IWORK(6) = 50000
      IWORK(6) = 20000
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
!      WRITE(*,*) ( Y(I), I=1, 13 )

      DO 
       CALL SVODE(FREACT,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,   &
                 IOPT,RWORK,LRW,IWORK,LIW,JACRACT,MF,RPAR,IPAR)
!      WRITE(6,20)T,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9), &
!                     Y(10),Y(11), Y(12), Y(13)
!  20  FORMAT(' At t =',D12.4,'   y =',13D14.6)
        IF (ISTATE .LT. -1 ) GO TO 80
        IF ( ISTATE .GT. 1 ) EXIT
        IF ( ISTATE .EQ. -1 ) ISTATE = 2  
       ENDDO  
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
!
!     Covert M to mg/l for some of the species
!
!      rate( 5 ) = rate( 5 ) * ( CO2MoleculeWeight * 1000.0 ) ! M to mg/l
!      rate( 6 ) = rate( 6 ) * ( HCO3MoleculeWeight * 1000.0 ) ! M to mg/l
!      rate( 7 ) = rate( 7 ) * ( HMoleculeWeight * 1000.0 ) ! M to mg/l
!      rate( 8 ) = rate( 8 ) * ( CO3MoleculeWeight * 1000.0 ) ! M to mg/l
!      rate( 9 ) = rate( 9 ) * ( CaMoleculeWeight * 1000.0 ) ! M to mg/l
!      rate( 11 ) = rate( 11 ) * ( OHMoleculeWeight * 1000.0 ) ! M to mg/l

      NitrifierConc( xloc, yloc, zloc ) = Y(12)      ! X1
      DenitrifierConc( xloc, yloc, zloc ) = Y(13)    ! X2
      RETURN
  80  WRITE(6,90)ISTATE, xloc, yloc, zloc
      rate = 0.0 

  90  FORMAT(///' Error halt: ISTATE =',I3, ' Cell (', I4, ', ', I4, ', ', I4, ')' )
!      STOP
      RETURN
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
      IWORK(6) = 5000
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

PD(1,1) = NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*np&
ars%NitOxMaxRate_1_day*O2*X1/((npars%NH4HalfSaturationConst_mg_l+&
NH4)**2*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+npars%Emperical&
Biomass1InhibitionConst_mg_l))-npars%EmpericalBiomass1InhibitionC&
onst_mg_l*npars%NitOxMaxRate_1_day*O2*X1/((npars%NH4HalfSaturatio&
nConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+npars%&
EmpericalBiomass1InhibitionConst_mg_l))
PD(1,2) = 0
PD(1,3) = 0
PD(1,4) = NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*np&
ars%NitOxMaxRate_1_day*O2*X1/((npars%NH4HalfSaturationConst_mg_l+&
NH4)*(O2+npars%O2HalfSaturationConst_mg_l)**2*(X1+npars%Emperical&
Biomass1InhibitionConst_mg_l))-NH4*npars%EmpericalBiomass1Inhibit&
ionConst_mg_l*npars%NitOxMaxRate_1_day*X1/((npars%NH4HalfSaturati&
onConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+npars&
%EmpericalBiomass1InhibitionConst_mg_l))
PD(1,5) = 0
PD(1,6) = 0
PD(1,7) = 0
PD(1,8) = 0
PD(1,9) = 0
PD(1,10) = 0
PD(1,11) = 0
PD(1,12) = NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*n&
pars%NitOxMaxRate_1_day*O2*X1/((npars%NH4HalfSaturationConst_mg_l&
+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+npars%EmpericalBi&
omass1InhibitionConst_mg_l)**2)-NH4*npars%EmpericalBiomass1Inhibi&
tionConst_mg_l*npars%NitOxMaxRate_1_day*O2/((npars%NH4HalfSaturat&
ionConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+npar&
s%EmpericalBiomass1InhibitionConst_mg_l))
PD(1,13) = 0
PD(2,1) = NO3MoleculeWeight*npars%EmpericalBiomass1Inhibitio&
nConst_mg_l*npars%NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeight*(np&
ars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationCo&
nst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_mg_l))-NO3Mo&
leculeWeight*NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*npar&
s%NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeight*(npars%NH4HalfSatur&
ationConst_mg_l+NH4)**2*(O2+npars%O2HalfSaturationConst_mg_l)*(X1&
+npars%EmpericalBiomass1InhibitionConst_mg_l))
PD(2,2) = (-4.0)*NO3MoleculeWeight*CH2O*npars%DenitOxMaxRate&
_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibi&
tionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturatio&
nConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npa&
rs%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCo&
nst_mg_l))+4.0*NO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_1_&
day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibitio&
nCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationCo&
nst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)**2*(O2+npa&
rs%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCo&
nst_mg_l))
PD(2,3) = (-4.0)*NO3MoleculeWeight*NO3*npars%DenitOxMaxRate_&
1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibit&
ionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturation&
Const_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npar&
s%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCon&
st_mg_l))+4.0*NO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_1_d&
ay*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibition&
Coef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationCon&
st_mg_l+CH2O)**2*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npar&
s%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCon&
st_mg_l))
PD(2,4) = 4.0*NO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRat&
e_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhib&
itionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturati&
onConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+np&
ars%O2InhibitionCoef_mg_l)**2*(X2+npars%EmpericalBiomass2Inhibiti&
onConst_mg_l))+NO3MoleculeWeight*NH4*npars%EmpericalBiomass1Inhib&
itionConst_mg_l*npars%NitOxMaxRate_1_day*X1/(NH4MoleculeWeight*(n&
pars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationC&
onst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_mg_l))-NO3M&
oleculeWeight*NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*npa&
rs%NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeight*(npars%NH4HalfSatu&
rationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)**2*(X&
1+npars%EmpericalBiomass1InhibitionConst_mg_l))
PD(2,5) = 0
PD(2,6) = 0
PD(2,7) = 0
PD(2,8) = 0
PD(2,9) = 0
PD(2,10) = 0
PD(2,11) = 0
PD(2,12) = NO3MoleculeWeight*NH4*npars%EmpericalBiomass1Inhi&
bitionConst_mg_l*npars%NitOxMaxRate_1_day*O2/(NH4MoleculeWeight*(&
npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2HalfSaturation&
Const_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_mg_l))-NO3&
MoleculeWeight*NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*np&
ars%NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeight*(npars%NH4HalfSat&
urationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+&
npars%EmpericalBiomass1InhibitionConst_mg_l)**2)
PD(2,13) = (-4.0)*NO3MoleculeWeight*CH2O*NO3*npars%DenitOxMa&
xRate_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2I&
nhibitionCoef_mg_l/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturat&
ionConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+n&
pars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2Inhibition&
Const_mg_l))+4.0*NO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_&
1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibit&
ionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturation&
Const_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npar&
s%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCon&
st_mg_l)**2)
PD(3,1) = 0
PD(3,2) = CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%Emperica&
lBiomass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/((np&
ars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3HalfSaturationCo&
nst_mg_l+NO3)**2*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%Emper&
icalBiomass2InhibitionConst_mg_l))-CH2O*npars%DenitOxMaxRate_1_da&
y*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2InhibitionC&
oef_mg_l*X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3&
HalfSaturationConst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X&
2+npars%EmpericalBiomass2InhibitionConst_mg_l))
PD(3,3) = -NO3*npars%DenitOxMaxRate_1_day*npars%EmpericalBio&
mass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/((npars%&
CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_&
mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBio&
mass2InhibitionConst_mg_l))+CH2O*NO3*npars%DenitOxMaxRate_1_day*n&
pars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2InhibitionCoef&
_mg_l*X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)**2*(npars%NO3&
HalfSaturationConst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X&
2+npars%EmpericalBiomass2InhibitionConst_mg_l))-1.0*npars%Emperic&
alBiomass2InhibitionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*&
X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSat&
urationConst_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_&
l))+1.0*CH2O*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%Or&
gCarbonOxMaxRate_1_day*O2*X2/((npars%CH2OHalfSaturationConst_mg_l&
+CH2O)**2*(O2+npars%O2HalfSaturationConst_mg_l)*(X2+npars%Emperic&
alBiomass2InhibitionConst_mg_l))
PD(3,4) = CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%Emperica&
lBiomass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/((np&
ars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3HalfSaturationCo&
nst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)**2*(X2+npars%Emper&
icalBiomass2InhibitionConst_mg_l))-1.0*CH2O*npars%EmpericalBiomas&
s2InhibitionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*X2/((npars%&
CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSaturationCons&
t_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l))+1.0*CH2&
O*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%OrgCarbonOxMa&
xRate_1_day*O2*X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+&
npars%O2HalfSaturationConst_mg_l)**2*(X2+npars%EmpericalBiomass2I&
nhibitionConst_mg_l))
PD(3,5) = 0
PD(3,6) = 0
PD(3,7) = 0
PD(3,8) = 0
PD(3,9) = 0
PD(3,10) = 0
PD(3,11) = 0
PD(3,12) = 0
PD(3,13) = -CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%Emperi&
calBiomass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l/((npa&
rs%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3HalfSaturationCon&
st_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%Emperical&
Biomass2InhibitionConst_mg_l))-1.0*CH2O*npars%EmpericalBiomass2In&
hibitionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2/((npars%CH2O&
HalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSaturationConst_mg&
_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l))+CH2O*NO3*np&
ars%DenitOxMaxRate_1_day*npars%EmpericalBiomass2InhibitionConst_m&
g_l*npars%O2InhibitionCoef_mg_l*X2/((npars%CH2OHalfSaturationCons&
t_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars%O2&
InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_m&
g_l)**2)+1.0*CH2O*npars%EmpericalBiomass2InhibitionConst_mg_l*npa&
rs%OrgCarbonOxMaxRate_1_day*O2*X2/((npars%CH2OHalfSaturationConst&
_mg_l+CH2O)*(O2+npars%O2HalfSaturationConst_mg_l)*(X2+npars%Emper&
icalBiomass2InhibitionConst_mg_l)**2)
PD(4,1) = 2*O2MoleculeWeight*NH4*npars%EmpericalBiomass1Inhi&
bitionConst_mg_l*npars%NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeigh&
t*(npars%NH4HalfSaturationConst_mg_l+NH4)**2*(O2+npars%O2HalfSatu&
rationConst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_mg_l&
))-2*O2MoleculeWeight*npars%EmpericalBiomass1InhibitionConst_mg_l&
*npars%NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeight*(npars%NH4Half&
SaturationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(&
X1+npars%EmpericalBiomass1InhibitionConst_mg_l))
PD(4,2) = 0
PD(4,3) = O2MoleculeWeight*CH2O*npars%EmpericalBiomass2Inhib&
itionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*X2/(CH2OMolecul&
eWeight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)**2*(O2+npars%O2&
HalfSaturationConst_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCo&
nst_mg_l))-O2MoleculeWeight*npars%EmpericalBiomass2InhibitionCons&
t_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*X2/(CH2OMoleculeWeight*(&
npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSaturati&
onConst_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l))
PD(4,4) = -O2MoleculeWeight*CH2O*npars%EmpericalBiomass2Inhi&
bitionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*X2/(CH2OMoleculeW&
eight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfS&
aturationConst_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_m&
g_l))+O2MoleculeWeight*CH2O*npars%EmpericalBiomass2InhibitionCons&
t_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*X2/(CH2OMoleculeWeight*(&
npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSaturati&
onConst_mg_l)**2*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l)&
)-2*O2MoleculeWeight*NH4*npars%EmpericalBiomass1InhibitionConst_m&
g_l*npars%NitOxMaxRate_1_day*X1/(NH4MoleculeWeight*(npars%NH4Half&
SaturationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(&
X1+npars%EmpericalBiomass1InhibitionConst_mg_l))+2*O2MoleculeWeig&
ht*NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*npars%NitOxMax&
Rate_1_day*O2*X1/(NH4MoleculeWeight*(npars%NH4HalfSaturationConst&
_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)**2*(X1+npars%Emp&
ericalBiomass1InhibitionConst_mg_l))
PD(4,5) = 0
PD(4,6) = 0
PD(4,7) = 0
PD(4,8) = 0
PD(4,9) = 0
PD(4,10) = 0
PD(4,11) = 0
PD(4,12) = 2*O2MoleculeWeight*NH4*npars%EmpericalBiomass1Inh&
ibitionConst_mg_l*npars%NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeig&
ht*(npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2HalfSatura&
tionConst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_mg_l)*&
*2)-2*O2MoleculeWeight*NH4*npars%EmpericalBiomass1InhibitionConst&
_mg_l*npars%NitOxMaxRate_1_day*O2/(NH4MoleculeWeight*(npars%NH4Ha&
lfSaturationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)&
*(X1+npars%EmpericalBiomass1InhibitionConst_mg_l))
PD(4,13) = O2MoleculeWeight*CH2O*npars%EmpericalBiomass2Inhi&
bitionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*X2/(CH2OMolecu&
leWeight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2Ha&
lfSaturationConst_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCons&
t_mg_l)**2)-O2MoleculeWeight*CH2O*npars%EmpericalBiomass2Inhibiti&
onConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2/(CH2OMoleculeWeigh&
t*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSatur&
ationConst_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l)&
)
PD(5,1) = 0
PD(5,2) = CH2O*npars%DenitOxMaxRate_1_day*npars%EmpericalBio&
mass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/(CH2OMol&
eculeWeight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3H&
alfSaturationConst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X2&
+npars%EmpericalBiomass2InhibitionConst_mg_l))/5000.0-CH2O*NO3*np&
ars%DenitOxMaxRate_1_day*npars%EmpericalBiomass2InhibitionConst_m&
g_l*npars%O2InhibitionCoef_mg_l*X2/(CH2OMoleculeWeight*(npars%CH2&
OHalfSaturationConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_&
l+NO3)**2*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBio&
mass2InhibitionConst_mg_l))/5000.0
PD(5,3) = NO3*npars%DenitOxMaxRate_1_day*npars%EmpericalBiom&
ass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/(CH2OMole&
culeWeight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3Ha&
lfSaturationConst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X2+&
npars%EmpericalBiomass2InhibitionConst_mg_l))/5000.0-CH2O*NO3*npa&
rs%DenitOxMaxRate_1_day*npars%EmpericalBiomass2InhibitionConst_mg&
_l*npars%O2InhibitionCoef_mg_l*X2/(CH2OMoleculeWeight*(npars%CH2O&
HalfSaturationConst_mg_l+CH2O)**2*(npars%NO3HalfSaturationConst_m&
g_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiom&
ass2InhibitionConst_mg_l))/5000.0+npars%EmpericalBiomass2Inhibiti&
onConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*X2/(CH2OMoleculeWe&
ight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSa&
turationConst_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg&
_l))/1000.0-CH2O*npars%EmpericalBiomass2InhibitionConst_mg_l*npar&
s%OrgCarbonOxMaxRate_1_day*O2*X2/(CH2OMoleculeWeight*(npars%CH2OH&
alfSaturationConst_mg_l+CH2O)**2*(O2+npars%O2HalfSaturationConst_&
mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l))/1000.0
PD(5,4) = -CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%Emperic&
alBiomass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/(CH&
2OMoleculeWeight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars&
%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l&
)**2*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l))/5000.0+CH2&
O*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%OrgCarbonOxMa&
xRate_1_day*X2/(CH2OMoleculeWeight*(npars%CH2OHalfSaturationConst&
_mg_l+CH2O)*(O2+npars%O2HalfSaturationConst_mg_l)*(X2+npars%Emper&
icalBiomass2InhibitionConst_mg_l))/1000.0-CH2O*npars%EmpericalBio&
mass2InhibitionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*X2/(C&
H2OMoleculeWeight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+n&
pars%O2HalfSaturationConst_mg_l)**2*(X2+npars%EmpericalBiomass2In&
hibitionConst_mg_l))/1000.0
PD(5,5) = -npars%K_f_H_1_mday*OH-npars%K_f_CO2_1_day-npars%A&
cc_1_cm*npars%K_f_2_cm_day
PD(5,6) = npars%K_b_H_1_day+H*npars%K_b_CO2_1_mday+2*Ca*HCO3&
*npars%Acc_1_cm*npars%K_b_2_cm_m2day
PD(5,7) = HCO3*npars%K_b_CO2_1_mday
PD(5,8) = 0
PD(5,9) = HCO3**2*npars%Acc_1_cm*npars%K_b_2_cm_m2day
PD(5,10) = 0
PD(5,11) = -CO2*npars%K_f_H_1_mday
PD(5,12) = 0
PD(5,13) = CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%Emperic&
alBiomass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l/(CH2OM&
oleculeWeight*(npars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO&
3HalfSaturationConst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(&
X2+npars%EmpericalBiomass2InhibitionConst_mg_l))/5000.0+CH2O*npar&
s%EmpericalBiomass2InhibitionConst_mg_l*npars%OrgCarbonOxMaxRate_&
1_day*O2/(CH2OMoleculeWeight*(npars%CH2OHalfSaturationConst_mg_l+&
CH2O)*(O2+npars%O2HalfSaturationConst_mg_l)*(X2+npars%EmpericalBi&
omass2InhibitionConst_mg_l))/1000.0-CH2O*NO3*npars%DenitOxMaxRate&
_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibi&
tionCoef_mg_l*X2/(CH2OMoleculeWeight*(npars%CH2OHalfSaturationCon&
st_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars%O&
2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_&
mg_l)**2)/5000.0-CH2O*npars%EmpericalBiomass2InhibitionConst_mg_l&
*npars%OrgCarbonOxMaxRate_1_day*O2*X2/(CH2OMoleculeWeight*(npars%&
CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSaturationCons&
t_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l)**2)/1000&
.0
PD(6,1) = 0
PD(6,2) = 4.0*HCO3MoleculeWeight*CH2O*npars%DenitOxMaxRate_1&
_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibiti&
onCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationC&
onst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars&
%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCons&
t_mg_l))+(-4.0)*HCO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_&
1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibit&
ionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturation&
Const_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)**2*(O2+n&
pars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2Inhibition&
Const_mg_l))
PD(6,3) = 4.0*HCO3MoleculeWeight*NO3*npars%DenitOxMaxRate_1_&
day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibitio&
nCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationCo&
nst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars%&
O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst&
_mg_l))+(-4.0)*HCO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_1&
_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibiti&
onCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationC&
onst_mg_l+CH2O)**2*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+np&
ars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionC&
onst_mg_l))
PD(6,4) = (-4.0)*HCO3MoleculeWeight*CH2O*NO3*npars%DenitOxMa&
xRate_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2I&
nhibitionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSatu&
rationConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O&
2+npars%O2InhibitionCoef_mg_l)**2*(X2+npars%EmpericalBiomass2Inhi&
bitionConst_mg_l))
PD(6,5) = npars%K_f_H_1_mday*OH+npars%K_f_CO2_1_day+npars%Ac&
c_1_cm*npars%K_f_2_cm_day
PD(6,6) = -npars%K_b_HCO3_1_mday*OH-npars%K_b_H_1_day-H*npar&
s%K_b_CO2_1_mday-2*Ca*HCO3*npars%Acc_1_cm*npars%K_b_2_cm_m2day-Ca&
*npars%Acc_1_cm*npars%K_b_1_cm_mday
PD(6,7) = npars%Acc_1_cm*npars%K_f_1_cm_day-HCO3*npars%K_b_C&
O2_1_mday
PD(6,8) = npars%K_f_HCO3_1_day
PD(6,9) = -HCO3**2*npars%Acc_1_cm*npars%K_b_2_cm_m2day-HCO3*&
npars%Acc_1_cm*npars%K_b_1_cm_mday
PD(6,10) = 0
PD(6,11) = CO2*npars%K_f_H_1_mday-HCO3*npars%K_b_HCO3_1_mday

PD(6,12) = 0
PD(6,13) = 4.0*HCO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxR&
ate_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inh&
ibitionCoef_mg_l/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturatio&
nConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npa&
rs%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCo&
nst_mg_l))+(-4.0)*HCO3MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRat&
e_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhib&
itionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturati&
onConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+np&
ars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionC&
onst_mg_l)**2)
PD(7,1) = npars%EmpericalBiomass1InhibitionConst_mg_l*npars%&
NitOxMaxRate_1_day*O2*X1/(NH4MoleculeWeight*(npars%NH4HalfSaturat&
ionConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+npar&
s%EmpericalBiomass1InhibitionConst_mg_l))/250.0-NH4*npars%Emperic&
alBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_1_day*O2*X1/(NH&
4MoleculeWeight*(npars%NH4HalfSaturationConst_mg_l+NH4)**2*(O2+np&
ars%O2HalfSaturationConst_mg_l)*(X1+npars%EmpericalBiomass1Inhibi&
tionConst_mg_l))/250.0
PD(7,2) = 0
PD(7,3) = 0
PD(7,4) = NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*np&
ars%NitOxMaxRate_1_day*X1/(NH4MoleculeWeight*(npars%NH4HalfSatura&
tionConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+npa&
rs%EmpericalBiomass1InhibitionConst_mg_l))/250.0-NH4*npars%Emperi&
calBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_1_day*O2*X1/(N&
H4MoleculeWeight*(npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npar&
s%O2HalfSaturationConst_mg_l)**2*(X1+npars%EmpericalBiomass1Inhib&
itionConst_mg_l))/250.0
PD(7,5) = npars%K_f_CO2_1_day
PD(7,6) = Ca*npars%Acc_1_cm*npars%K_b_1_cm_mday-H*npars%K_b_&
CO2_1_mday
PD(7,7) = -npars%K_b_H2O_1_mday*OH-npars%Acc_1_cm*npars%K_f_&
1_cm_day-HCO3*npars%K_b_CO2_1_mday
PD(7,8) = 0
PD(7,9) = HCO3*npars%Acc_1_cm*npars%K_b_1_cm_mday
PD(7,10) = 0
PD(7,11) = -H*npars%K_b_H2O_1_mday
PD(7,12) = NH4*npars%EmpericalBiomass1InhibitionConst_mg_l*n&
pars%NitOxMaxRate_1_day*O2/(NH4MoleculeWeight*(npars%NH4HalfSatur&
ationConst_mg_l+NH4)*(O2+npars%O2HalfSaturationConst_mg_l)*(X1+np&
ars%EmpericalBiomass1InhibitionConst_mg_l))/250.0-NH4*npars%Emper&
icalBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_1_day*O2*X1/(&
NH4MoleculeWeight*(npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npa&
rs%O2HalfSaturationConst_mg_l)*(X1+npars%EmpericalBiomass1Inhibit&
ionConst_mg_l)**2)/250.0
PD(7,13) = 0
PD(8,1) = 0
PD(8,2) = 0
PD(8,3) = 0
PD(8,4) = 0
PD(8,5) = 0
PD(8,6) = npars%K_b_HCO3_1_mday*OH
PD(8,7) = 0
PD(8,8) = -npars%K_f_HCO3_1_day-Ca*npars%Acc_1_cm*npars%K_b_3_cm_mday
PD(8,9) = -CO3*npars%Acc_1_cm*npars%K_b_3_cm_mday
PD(8,10) = 0
PD(8,11) = HCO3*npars%K_b_HCO3_1_mday
PD(8,12) = 0
PD(8,13) = 0
PD(9,1) = 0
PD(9,2) = 0
PD(9,3) = 0
PD(9,4) = 0
PD(9,5) = npars%Acc_1_cm*npars%K_f_2_cm_day
PD(9,6) = -2*Ca*HCO3*npars%Acc_1_cm*npars%K_b_2_cm_m2day-Ca*&
npars%Acc_1_cm*npars%K_b_1_cm_mday
PD(9,7) = npars%Acc_1_cm*npars%K_f_1_cm_day
PD(9,8) = -Ca*npars%Acc_1_cm*npars%K_b_3_cm_mday
PD(9,9) = -CO3*npars%Acc_1_cm*npars%K_b_3_cm_mday-HCO3**2*np&
ars%Acc_1_cm*npars%K_b_2_cm_m2day-HCO3*npars%Acc_1_cm*npars%K_b_1&
_cm_mday
PD(9,10) = 0
PD(9,11) = 0
PD(9,12) = 0
PD(9,13) = 0
PD(10,1) = 0
PD(10,2) = 2.0*N2MoleculeWeight*CH2O*npars%DenitOxMaxRate_1_&
day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibitio&
nCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationCo&
nst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars%&
O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst&
_mg_l))+(-2.0)*N2MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_1_d&
ay*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibition&
Coef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationCon&
st_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)**2*(O2+npar&
s%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCon&
st_mg_l))
PD(10,3) = 2.0*N2MoleculeWeight*NO3*npars%DenitOxMaxRate_1_d&
ay*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibition&
Coef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationCon&
st_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars%O&
2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_&
mg_l))+(-2.0)*N2MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_1_da&
y*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2InhibitionC&
oef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationCons&
t_mg_l+CH2O)**2*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars&
%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCons&
t_mg_l))
PD(10,4) = (-2.0)*N2MoleculeWeight*CH2O*NO3*npars%DenitOxMax&
Rate_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2In&
hibitionCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSatur&
ationConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2&
+npars%O2InhibitionCoef_mg_l)**2*(X2+npars%EmpericalBiomass2Inhib&
itionConst_mg_l))
PD(10,5) = 0
PD(10,6) = 0
PD(10,7) = 0
PD(10,8) = 0
PD(10,9) = 0
PD(10,10) = 0
PD(10,11) = 0
PD(10,12) = 0
PD(10,13) = 2.0*N2MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRa&
te_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhi&
bitionCoef_mg_l/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturation&
Const_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npar&
s%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCon&
st_mg_l))+(-2.0)*N2MoleculeWeight*CH2O*NO3*npars%DenitOxMaxRate_1&
_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inhibiti&
onCoef_mg_l*X2/(5.0*CH2OMoleculeWeight*(npars%CH2OHalfSaturationC&
onst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars&
%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2InhibitionCons&
t_mg_l)**2)
PD(11,1) = 0
PD(11,2) = 0
PD(11,3) = 0
PD(11,4) = 0
PD(11,5) = -npars%K_f_H_1_mday*OH
PD(11,6) = npars%K_b_H_1_day-npars%K_b_HCO3_1_mday*OH
PD(11,7) = -npars%K_b_H2O_1_mday*OH
PD(11,8) = npars%K_f_HCO3_1_day
PD(11,9) = 0
PD(11,10) = 0
PD(11,11) = -CO2*npars%K_f_H_1_mday-HCO3*npars%K_b_HCO3_1_md&
ay-H*npars%K_b_H2O_1_mday
PD(11,12) = 0
PD(11,13) = 0
PD(12,1) = npars%Biomass1YieldCoeff_M_biomass_M_substrate*np&
ars%EmpericalBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_1_da&
y*O2*X1/((npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2Half&
SaturationConst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_&
mg_l))-NH4*npars%Biomass1YieldCoeff_M_biomass_M_substrate*npars%E&
mpericalBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_1_day*O2*&
X1/((npars%NH4HalfSaturationConst_mg_l+NH4)**2*(O2+npars%O2HalfSa&
turationConst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_mg&
_l))
PD(12,2) = 0
PD(12,3) = 0
PD(12,4) = NH4*npars%Biomass1YieldCoeff_M_biomass_M_substrat&
e*npars%EmpericalBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_&
1_day*X1/((npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2Hal&
fSaturationConst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst&
_mg_l))-NH4*npars%Biomass1YieldCoeff_M_biomass_M_substrate*npars%&
EmpericalBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_1_day*O2&
*X1/((npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2HalfSatu&
rationConst_mg_l)**2*(X1+npars%EmpericalBiomass1InhibitionConst_m&
g_l))
PD(12,5) = 0
PD(12,6) = 0
PD(12,7) = 0
PD(12,8) = 0
PD(12,9) = 0
PD(12,10) = 0
PD(12,11) = 0
PD(12,12) = NH4*npars%Biomass1YieldCoeff_M_biomass_M_substra&
te*npars%EmpericalBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate&
_1_day*O2/((npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2Ha&
lfSaturationConst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionCons&
t_mg_l))-NH4*npars%Biomass1YieldCoeff_M_biomass_M_substrate*npars&
%EmpericalBiomass1InhibitionConst_mg_l*npars%NitOxMaxRate_1_day*O&
2*X1/((npars%NH4HalfSaturationConst_mg_l+NH4)*(O2+npars%O2HalfSat&
urationConst_mg_l)*(X1+npars%EmpericalBiomass1InhibitionConst_mg_&
l)**2)-npars%Biomass1SpecDecayConst_1_day
PD(12,13) = 0
PD(13,1) = 0
PD(13,2) = -npars%Biomass2YieldCoeff_M_biomass_M_substrate*(&
CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%EmpericalBiomass2Inhibi&
tionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/((npars%CH2OHalfSat&
urationConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)**&
2*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2Inh&
ibitionConst_mg_l))-CH2O*npars%DenitOxMaxRate_1_day*npars%Emperic&
alBiomass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/((n&
pars%CH2OHalfSaturationConst_mg_l+CH2O)*(npars%NO3HalfSaturationC&
onst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%Emperic&
alBiomass2InhibitionConst_mg_l)))
PD(13,3) = -npars%Biomass2YieldCoeff_M_biomass_M_substrate*(&
-NO3*npars%DenitOxMaxRate_1_day*npars%EmpericalBiomass2Inhibition&
Const_mg_l*npars%O2InhibitionCoef_mg_l*X2/((npars%CH2OHalfSaturat&
ionConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O2+n&
pars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2Inhibition&
Const_mg_l))+CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%EmpericalB&
iomass2InhibitionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/((npar&
s%CH2OHalfSaturationConst_mg_l+CH2O)**2*(npars%NO3HalfSaturationC&
onst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_mg_l)*(X2+npars%Emperic&
alBiomass2InhibitionConst_mg_l))-1.0*npars%EmpericalBiomass2Inhib&
itionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*X2/((npars%CH2O&
HalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSaturationConst_mg&
_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l))+1.0*CH2O*np&
ars%EmpericalBiomass2InhibitionConst_mg_l*npars%OrgCarbonOxMaxRat&
e_1_day*O2*X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)**2*(O2+n&
pars%O2HalfSaturationConst_mg_l)*(X2+npars%EmpericalBiomass2Inhib&
itionConst_mg_l)))
PD(13,4) = -npars%Biomass2YieldCoeff_M_biomass_M_substrate*(&
CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%EmpericalBiomass2Inhibi&
tionConst_mg_l*npars%O2InhibitionCoef_mg_l*X2/((npars%CH2OHalfSat&
urationConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(&
O2+npars%O2InhibitionCoef_mg_l)**2*(X2+npars%EmpericalBiomass2Inh&
ibitionConst_mg_l))-1.0*CH2O*npars%EmpericalBiomass2InhibitionCon&
st_mg_l*npars%OrgCarbonOxMaxRate_1_day*X2/((npars%CH2OHalfSaturat&
ionConst_mg_l+CH2O)*(O2+npars%O2HalfSaturationConst_mg_l)*(X2+npa&
rs%EmpericalBiomass2InhibitionConst_mg_l))+1.0*CH2O*npars%Emperic&
alBiomass2InhibitionConst_mg_l*npars%OrgCarbonOxMaxRate_1_day*O2*&
X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2+npars%O2HalfSat&
urationConst_mg_l)**2*(X2+npars%EmpericalBiomass2InhibitionConst_&
mg_l)))
PD(13,5) = 0
PD(13,6) = 0
PD(13,7) = 0
PD(13,8) = 0
PD(13,9) = 0
PD(13,10) = 0
PD(13,11) = 0
PD(13,12) = 0
PD(13,13) = -npars%Biomass2YieldCoeff_M_biomass_M_substrate*&
(-CH2O*NO3*npars%DenitOxMaxRate_1_day*npars%EmpericalBiomass2Inhi&
bitionConst_mg_l*npars%O2InhibitionCoef_mg_l/((npars%CH2OHalfSatu&
rationConst_mg_l+CH2O)*(npars%NO3HalfSaturationConst_mg_l+NO3)*(O&
2+npars%O2InhibitionCoef_mg_l)*(X2+npars%EmpericalBiomass2Inhibit&
ionConst_mg_l))-1.0*CH2O*npars%EmpericalBiomass2InhibitionConst_m&
g_l*npars%OrgCarbonOxMaxRate_1_day*O2/((npars%CH2OHalfSaturationC&
onst_mg_l+CH2O)*(O2+npars%O2HalfSaturationConst_mg_l)*(X2+npars%E&
mpericalBiomass2InhibitionConst_mg_l))+CH2O*NO3*npars%DenitOxMaxR&
ate_1_day*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%O2Inh&
ibitionCoef_mg_l*X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)*(n&
pars%NO3HalfSaturationConst_mg_l+NO3)*(O2+npars%O2InhibitionCoef_&
mg_l)*(X2+npars%EmpericalBiomass2InhibitionConst_mg_l)**2)+1.0*CH&
2O*npars%EmpericalBiomass2InhibitionConst_mg_l*npars%OrgCarbonOxM&
axRate_1_day*O2*X2/((npars%CH2OHalfSaturationConst_mg_l+CH2O)*(O2&
+npars%O2HalfSaturationConst_mg_l)*(X2+npars%EmpericalBiomass2Inh&
ibitionConst_mg_l)**2))-npars%Biomass2SpecDecayConst_1_day

!PD(1,1) = 1.0*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))-1.0&
!*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,2) = 0
!PD(1,3) = 0
!PD(1,4) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))-1.0&
!*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,5) = 0
!PD(1,6) = 0
!PD(1,7) = 0
!PD(1,8) = 0
!PD(1,9) = 0
!PD(1,10) = 0
!PD(1,11) = 0
!PD(1,12) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)-1.&
!0*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,13) = 0
!PD(2,1) = 3.444444444444445*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.&
!0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
!PD(2,2) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
!**2*(O2+1.0)*(X2+0.5))-8.266666666666666*CH2O*X2/((CH2O+10)*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5))
!PD(2,3) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))-8.266666666666666*NO3*X2/((CH2O+10)*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))
!PD(2,4) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
!*(O2+1.0)**2*(X2+0.5))+3.444444444444445*NH4*X1/((NH4+0.1)*(O2+0.&
!1)*(X1+1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(&
!X1+1.0))
!PD(2,5) = 0
!PD(2,6) = 0
!PD(2,7) = 0
!PD(2,8) = 0
!PD(2,9) = 0
!PD(2,10) = 0
!PD(2,11) = 0
!PD(2,12) = 3.444444444444445*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+&
!1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2&
!)
!PD(2,13) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)*(O2+1.0)*(X2+0.5)**2)-8.266666666666666*CH2O*NO3/((CH2O+10)*(NO&
!3+0.5)*(O2+1.0)*(X2+0.5))
!PD(3,1) = 0
!PD(3,2) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2+1.0)*(&
!X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5))
!PD(3,3) = -5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)&
!)+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5))-5.0*&
!O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+10)**2*&
!(O2+0.1)*(X2+0.5))
!PD(3,4) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)**2*(&
!X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2&
!/((CH2O+10)*(O2+0.1)**2*(X2+0.5))
!PD(3,5) = 0
!PD(3,6) = 0
!PD(3,7) = 0
!PD(3,8) = 0
!PD(3,9) = 0
!PD(3,10) = 0
!PD(3,11) = 0
!PD(3,12) = 0
!PD(3,13) = -5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0&
!.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO3*X2/((&
!CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((CH2O+10&
!)*(O2+0.1)*(X2+0.5)**2)
!PD(4,1) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)&
!*(X1+1.0))-3.555555555555555*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(4,2) = 0
!PD(4,3) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)**2*(O2+0.1&
!)*(X2+0.5))-5.333333333333333*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))
!PD(4,4) = -5.333333333333333*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2&
!+0.5))+5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.&
!5))-3.555555555555555*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))+3.5555&
!55555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
!PD(4,5) = 0
!PD(4,6) = 0
!PD(4,7) = 0
!PD(4,8) = 0
!PD(4,9) = 0
!PD(4,10) = 0
!PD(4,11) = 0
!PD(4,12) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(&
!X1+1.0)**2)-3.555555555555555*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0)&
!)
!PD(4,13) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)*&
!(X2+0.5)**2)-5.333333333333333*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.&
!5))
!PD(5,1) = 0
!PD(5,2) = 3.3333333333333335e-5*CH2O*X2/((CH2O+10)*(NO3+0.5)&
!*(O2+1.0)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*&
!(NO3+0.5)**2*(O2+1.0)*(X2+0.5))
!PD(5,3) = 3.3333333333333335e-5*NO3*X2/((CH2O+10)*(NO3+0.5)*&
!(O2+1.0)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)**&
!2*(NO3+0.5)*(O2+1.0)*(X2+0.5))+1.666666666666667e-4*O2*X2/((CH2O+&
!10)*(O2+0.1)*(X2+0.5))-1.666666666666667e-4*CH2O*O2*X2/((CH2O+10)&
!**2*(O2+0.1)*(X2+0.5))
!PD(5,4) = -3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*(NO3&
!+0.5)*(O2+1.0)**2*(X2+0.5))+1.666666666666667e-4*CH2O*X2/((CH2O+1&
!0)*(O2+0.1)*(X2+0.5))-1.666666666666667e-4*CH2O*O2*X2/((CH2O+10)*&
!(O2+0.1)**2*(X2+0.5))
!PD(5,5) = -7.344e+8*OH-2592.016416
!PD(5,6) = 866.4*Ca*HCO3+6.818e+9*H+8.63
!PD(5,7) = 6.818e+9*HCO3
!PD(5,8) = 0
!PD(5,9) = 433.2*HCO3**2
!PD(5,10) = 0
!PD(5,11) = -7.344e+8*CO2
!PD(5,12) = 0
!PD(5,13) = 3.3333333333333335e-5*CH2O*NO3/((CH2O+10)*(NO3+0.&
!5)*(O2+1.0)*(X2+0.5))+1.666666666666667e-4*CH2O*O2/((CH2O+10)*(O2&
!+0.1)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*(NO3&
!+0.5)*(O2+1.0)*(X2+0.5)**2)-1.666666666666667e-4*CH2O*O2*X2/((CH2&
!O+10)*(O2+0.1)*(X2+0.5)**2)
!PD(6,1) = 0
!PD(6,2) = 8.133333333333333*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2&
!+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)**2*(O2+1.0)*(X2+0.5))
!PD(6,3) = 8.133333333333333*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+&
!1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))
!PD(6,4) = -8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)*(O2+1.0)**2*(X2+0.5))
!PD(6,5) = 7.344e+8*OH+2592.016416
!PD(6,6) = -9.22e+14*OH-866.4*Ca*HCO3-6.818e+9*H-0.29222*Ca-8&
!.63
!PD(6,7) = 29.222-6.818e+9*HCO3
!PD(6,8) = 1.11e+11
!PD(6,9) = -433.2*HCO3**2-0.29222*HCO3
!PD(6,10) = 0
!PD(6,11) = 7.344e+8*CO2-9.22e+14*HCO3
!PD(6,12) = 0
!PD(6,13) = 8.133333333333333*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(&
!O2+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5)**2)
!PD(7,1) = 2.2222222222222223e-4*O2*X1/((NH4+0.1)*(O2+0.1)*(X&
!1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X&
!1+1.0))
!PD(7,2) = 0
!PD(7,3) = 0
!PD(7,4) = 2.2222222222222223e-4*NH4*X1/((NH4+0.1)*(O2+0.1)*(&
!X1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(&
!X1+1.0))
!PD(7,5) = 2592.0
!PD(7,6) = 0.29222*Ca-6.818e+9*H
!PD(7,7) = -2.46e+16*OH-6.818e+9*HCO3-29.222
!PD(7,8) = 0
!PD(7,9) = 0.29222*HCO3
!PD(7,10) = 0
!PD(7,11) = -2.46e+16*H
!PD(7,12) = 2.2222222222222223e-4*NH4*O2/((NH4+0.1)*(O2+0.1)*&
!(X1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1&
!+1.0)**2)
!PD(7,13) = 0
!PD(8,1) = 0
!PD(8,2) = 0
!PD(8,3) = 0
!PD(8,4) = 0
!PD(8,5) = 0
!PD(8,6) = 9.22e+14*OH
!PD(8,7) = 0
!PD(8,8) = -0.0057456*Ca-1.11e+11
!PD(8,9) = -0.0057456*CO3
!PD(8,10) = 0
!PD(8,11) = 9.22e+14*HCO3
!PD(8,12) = 0
!PD(8,13) = 0
!PD(9,1) = 0
!PD(9,2) = 0
!PD(9,3) = 0
!PD(9,4) = 0
!PD(9,5) = 0.016416
!PD(9,6) = -866.4*Ca*HCO3-0.29222*Ca
!PD(9,7) = 29.222
!PD(9,8) = -0.0057456*Ca
!PD(9,9) = -433.2*HCO3**2-0.29222*HCO3-0.0057456*CO3
!PD(9,10) = 0
!PD(9,11) = 0
!PD(9,12) = 0
!PD(9,13) = 0
!PD(10,1) = 0
!PD(10,2) = 1.866666666666667*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O&
!2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
!5)**2*(O2+1.0)*(X2+0.5))
!PD(10,3) = 1.866666666666667*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2&
!+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)**2*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5))
!PD(10,4) = -1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
!5)*(O2+1.0)**2*(X2+0.5))
!PD(10,5) = 0
!PD(10,6) = 0
!PD(10,7) = 0
!PD(10,8) = 0
!PD(10,9) = 0
!PD(10,10) = 0
!PD(10,11) = 0
!PD(10,12) = 0
!PD(10,13) = 1.866666666666667*CH2O*NO3/((CH2O+10)*(NO3+0.5)*&
!(O2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5)**2)
!PD(11,1) = 0
!PD(11,2) = 0
!PD(11,3) = 0
!PD(11,4) = 0
!PD(11,5) = -7.344e+8*OH
!PD(11,6) = 8.63-9.22e+14*OH
!PD(11,7) = -2.46e+16*OH
!PD(11,8) = 1.11e+11
!PD(11,9) = 0
!PD(11,10) = 0
!PD(11,11) = -9.22e+14*HCO3-2.46e+16*H-7.344e+8*CO2
!PD(11,12) = 0
!PD(11,13) = 0
!PD(12,1) = 0.17*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH4&
!*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
!PD(12,2) = 0
!PD(12,3) = 0
!PD(12,4) = 0.17*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH&
!4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
!PD(12,5) = 0
!PD(12,6) = 0
!PD(12,7) = 0
!PD(12,8) = 0
!PD(12,9) = 0
!PD(12,10) = 0
!PD(12,11) = 0
!PD(12,12) = 0.17*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*N&
!H4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)
!PD(12,13) = 0
!PD(13,1) = 0
!PD(13,2) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2&
!+1.0)*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5&
!)))
!PD(13,3) = -0.5*(-5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(&
!X2+0.5))+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5&
!))-5.0*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+&
!10)**2*(O2+0.1)*(X2+0.5)))
!PD(13,4) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.&
!0)**2*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2&
!O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.5)))
!PD(13,5) = 0
!PD(13,6) = 0
!PD(13,7) = 0
!PD(13,8) = 0
!PD(13,9) = 0
!PD(13,10) = 0
!PD(13,11) = 0
!PD(13,12) = 0
!PD(13,13) = -0.5*(-5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0&
!)*(X2+0.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO&
!3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((&
!CH2O+10)*(O2+0.1)*(X2+0.5)**2))
!
! 
! Chen's parameters
!PD(1,1) = 1.0*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))-1.0&
!*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,2) = 0
!PD(1,3) = 0
!PD(1,4) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))-1.0&
!*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,5) = 0
!PD(1,6) = 0
!PD(1,7) = 0
!PD(1,8) = 0
!PD(1,9) = 0
!PD(1,10) = 0
!PD(1,11) = 0
!PD(1,12) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)-1.&
!0*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,13) = 0
!PD(2,1) = 3.444444444444445*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.&
!0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
!PD(2,2) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
!**2*(O2+1.0)*(X2+0.5))-8.266666666666666*CH2O*X2/((CH2O+10)*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5))
!PD(2,3) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))-8.266666666666666*NO3*X2/((CH2O+10)*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))
!PD(2,4) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
!*(O2+1.0)**2*(X2+0.5))+3.444444444444445*NH4*X1/((NH4+0.1)*(O2+0.&
!1)*(X1+1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(&
!X1+1.0))
!PD(2,5) = 0
!PD(2,6) = 0
!PD(2,7) = 0
!PD(2,8) = 0
!PD(2,9) = 0
!PD(2,10) = 0
!PD(2,11) = 0
!PD(2,12) = 3.444444444444445*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+&
!1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2&
!)
!PD(2,13) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)*(O2+1.0)*(X2+0.5)**2)-8.266666666666666*CH2O*NO3/((CH2O+10)*(NO&
!3+0.5)*(O2+1.0)*(X2+0.5))
!PD(3,1) = 0
!PD(3,2) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2+1.0)*(&
!X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5))
!PD(3,3) = -5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)&
!)+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5))-5.0*&
!O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+10)**2*&
!(O2+0.1)*(X2+0.5))
!PD(3,4) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)**2*(&
!X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2&
!/((CH2O+10)*(O2+0.1)**2*(X2+0.5))
!PD(3,5) = 0
!PD(3,6) = 0
!PD(3,7) = 0
!PD(3,8) = 0
!PD(3,9) = 0
!PD(3,10) = 0
!PD(3,11) = 0
!PD(3,12) = 0
!PD(3,13) = -5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0&
!.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO3*X2/((&
!CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((CH2O+10&
!)*(O2+0.1)*(X2+0.5)**2)
!PD(4,1) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)&
!*(X1+1.0))-3.555555555555555*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(4,2) = 0
!PD(4,3) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)**2*(O2+0.1&
!)*(X2+0.5))-5.333333333333333*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))
!PD(4,4) = -5.333333333333333*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2&
!+0.5))+5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.&
!5))-3.555555555555555*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))+3.5555&
!55555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
!PD(4,5) = 0
!PD(4,6) = 0
!PD(4,7) = 0
!PD(4,8) = 0
!PD(4,9) = 0
!PD(4,10) = 0
!PD(4,11) = 0
!PD(4,12) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(&
!X1+1.0)**2)-3.555555555555555*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0)&
!)
!PD(4,13) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)*&
!(X2+0.5)**2)-5.333333333333333*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.&
!5))
!PD(5,1) = 0
!PD(5,2) = 3.3333333333333335e-5*CH2O*X2/((CH2O+10)*(NO3+0.5)&
!*(O2+1.0)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*&
!(NO3+0.5)**2*(O2+1.0)*(X2+0.5))
!PD(5,3) = 3.3333333333333335e-5*NO3*X2/((CH2O+10)*(NO3+0.5)*&
!(O2+1.0)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)**&
!2*(NO3+0.5)*(O2+1.0)*(X2+0.5))+1.666666666666667e-4*O2*X2/((CH2O+&
!10)*(O2+0.1)*(X2+0.5))-1.666666666666667e-4*CH2O*O2*X2/((CH2O+10)&
!**2*(O2+0.1)*(X2+0.5))
!PD(5,4) = -3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*(NO3&
!+0.5)*(O2+1.0)**2*(X2+0.5))+1.666666666666667e-4*CH2O*X2/((CH2O+1&
!0)*(O2+0.1)*(X2+0.5))-1.666666666666667e-4*CH2O*O2*X2/((CH2O+10)*&
!(O2+0.1)**2*(X2+0.5))
!PD(5,5) = -7.344e+8*OH-2.608416
!PD(5,6) = 86640.0*Ca*HCO3+6818000.0*H+8.63
!PD(5,7) = 6818000.0*HCO3
!PD(5,8) = 0
!PD(5,9) = 43320.0*HCO3**2
!PD(5,10) = 0
!PD(5,11) = -7.344e+8*CO2
!PD(5,12) = 0
!PD(5,13) = 3.3333333333333335e-5*CH2O*NO3/((CH2O+10)*(NO3+0.&
!5)*(O2+1.0)*(X2+0.5))+1.666666666666667e-4*CH2O*O2/((CH2O+10)*(O2&
!+0.1)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*(NO3&
!+0.5)*(O2+1.0)*(X2+0.5)**2)-1.666666666666667e-4*CH2O*O2*X2/((CH2&
!O+10)*(O2+0.1)*(X2+0.5)**2)
!PD(6,1) = 0
!PD(6,2) = 8.133333333333333*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2&
!+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)**2*(O2+1.0)*(X2+0.5))
!PD(6,3) = 8.133333333333333*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+&
!1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))
!PD(6,4) = -8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)*(O2+1.0)**2*(X2+0.5))
!PD(6,5) = 7.344e+8*OH+2.608416
!PD(6,6) = -9.216e+7*OH-86640.0*Ca*HCO3-6818000.0*H-0.29222*C&
!a-8.63
!PD(6,7) = 29.222-6818000.0*HCO3
!PD(6,8) = 10800.0
!PD(6,9) = -43320.0*HCO3**2-0.29222*HCO3
!PD(6,10) = 0
!PD(6,11) = 7.344e+8*CO2-9.216e+7*HCO3
!PD(6,12) = 0
!PD(6,13) = 8.133333333333333*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(&
!O2+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5)**2)
!PD(7,1) = 2.2222222222222223e-4*O2*X1/((NH4+0.1)*(O2+0.1)*(X&
!1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X&
!1+1.0))
!PD(7,2) = 0
!PD(7,3) = 0
!PD(7,4) = 2.2222222222222223e-4*NH4*X1/((NH4+0.1)*(O2+0.1)*(&
!X1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(&
!X1+1.0))
!PD(7,5) = 2.592
!PD(7,6) = 0.29222*Ca-6818000.0*H
!PD(7,7) = -2.46e+10*OH-6818000.0*HCO3-29.222
!PD(7,8) = 0
!PD(7,9) = 0.29222*HCO3
!PD(7,10) = 0
!PD(7,11) = -2.46e+10*H
!PD(7,12) = 2.2222222222222223e-4*NH4*O2/((NH4+0.1)*(O2+0.1)*&
!(X1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1&
!+1.0)**2)
!PD(7,13) = 0
!PD(8,1) = 0
!PD(8,2) = 0
!PD(8,3) = 0
!PD(8,4) = 0
!PD(8,5) = 0
!PD(8,6) = 9.216e+7*OH
!PD(8,7) = 0
!PD(8,8) = -5745.599999999999*Ca-10800.0
!PD(8,9) = -5745.599999999999*CO3
!PD(8,10) = 0
!PD(8,11) = 9.216e+7*HCO3
!PD(8,12) = 0
!PD(8,13) = 0
!PD(9,1) = 0
!PD(9,2) = 0
!PD(9,3) = 0
!PD(9,4) = 0
!PD(9,5) = 0.016416
!PD(9,6) = -86640.0*Ca*HCO3-0.29222*Ca
!PD(9,7) = 29.222
!PD(9,8) = -5745.599999999999*Ca
!PD(9,9) = -43320.0*HCO3**2-0.29222*HCO3-5745.599999999999*CO&
!3
!PD(9,10) = 0
!PD(9,11) = 0
!PD(9,12) = 0
!PD(9,13) = 0
!PD(10,1) = 0
!PD(10,2) = 1.866666666666667*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O&
!2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
!5)**2*(O2+1.0)*(X2+0.5))
!PD(10,3) = 1.866666666666667*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2&
!+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)**2*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5))
!PD(10,4) = -1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
!5)*(O2+1.0)**2*(X2+0.5))
!PD(10,5) = 0
!PD(10,6) = 0
!PD(10,7) = 0
!PD(10,8) = 0
!PD(10,9) = 0
!PD(10,10) = 0
!PD(10,11) = 0
!PD(10,12) = 0
!PD(10,13) = 1.866666666666667*CH2O*NO3/((CH2O+10)*(NO3+0.5)*&
!(O2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5)**2)
!PD(11,1) = 0
!PD(11,2) = 0
!PD(11,3) = 0
!PD(11,4) = 0
!PD(11,5) = -7.344e+8*OH
!PD(11,6) = 8.63-9.216e+7*OH
!PD(11,7) = -2.46e+10*OH
!PD(11,8) = 10800.0
!PD(11,9) = 0
!PD(11,10) = 0
!PD(11,11) = -9.216e+7*HCO3-2.46e+10*H-7.344e+8*CO2
!PD(11,12) = 0
!PD(11,13) = 0
!PD(12,1) = 0.17*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH4&
!*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
!PD(12,2) = 0
!PD(12,3) = 0
!PD(12,4) = 0.17*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH&
!4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
!PD(12,5) = 0
!PD(12,6) = 0
!PD(12,7) = 0
!PD(12,8) = 0
!PD(12,9) = 0
!PD(12,10) = 0
!PD(12,11) = 0
!PD(12,12) = 0.17*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*N&
!H4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)
!PD(12,13) = 0
!PD(13,1) = 0
!PD(13,2) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2&
!+1.0)*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5&
!)))
!PD(13,3) = -0.5*(-5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(&
!X2+0.5))+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5&
!))-5.0*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+&
!10)**2*(O2+0.1)*(X2+0.5)))
!PD(13,4) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.&
!0)**2*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2&
!O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.5)))
!PD(13,5) = 0
!PD(13,6) = 0
!PD(13,7) = 0
!PD(13,8) = 0
!PD(13,9) = 0
!PD(13,10) = 0
!PD(13,11) = 0
!PD(13,12) = 0
!PD(13,13) = -0.5*(-5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0&
!)*(X2+0.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO&
!3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((&
!CH2O+10)*(O2+0.1)*(X2+0.5)**2))
!
! MacQuarrie's parameter
!
!PD(1,1) = 1.0*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))-1.0&
!*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,2) = 0
!PD(1,3) = 0
!PD(1,4) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))-1.0&
!*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,5) = 0
!PD(1,6) = 0
!PD(1,7) = 0
!PD(1,8) = 0
!PD(1,9) = 0
!PD(1,10) = 0
!PD(1,11) = 0
!PD(1,12) = 1.0*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)-1.&
!0*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(1,13) = 0
!PD(2,1) = 3.444444444444445*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.&
!0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
!PD(2,2) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
!**2*(O2+1.0)*(X2+0.5))-8.266666666666666*CH2O*X2/((CH2O+10)*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5))
!PD(2,3) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))-8.266666666666666*NO3*X2/((CH2O+10)*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))
!PD(2,4) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)&
!*(O2+1.0)**2*(X2+0.5))+3.444444444444445*NH4*X1/((NH4+0.1)*(O2+0.&
!1)*(X1+1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(&
!X1+1.0))
!PD(2,5) = 0
!PD(2,6) = 0
!PD(2,7) = 0
!PD(2,8) = 0
!PD(2,9) = 0
!PD(2,10) = 0
!PD(2,11) = 0
!PD(2,12) = 3.444444444444445*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+&
!1.0))-3.444444444444445*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2&
!)
!PD(2,13) = 8.266666666666666*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)*(O2+1.0)*(X2+0.5)**2)-8.266666666666666*CH2O*NO3/((CH2O+10)*(NO&
!3+0.5)*(O2+1.0)*(X2+0.5))
!PD(3,1) = 0
!PD(3,2) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2+1.0)*(&
!X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5))
!PD(3,3) = -5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)&
!)+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5))-5.0*&
!O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+10)**2*&
!(O2+0.1)*(X2+0.5))
!PD(3,4) = 5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)**2*(&
!X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2&
!/((CH2O+10)*(O2+0.1)**2*(X2+0.5))
!PD(3,5) = 0
!PD(3,6) = 0
!PD(3,7) = 0
!PD(3,8) = 0
!PD(3,9) = 0
!PD(3,10) = 0
!PD(3,11) = 0
!PD(3,12) = 0
!PD(3,13) = -5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0&
!.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO3*X2/((&
!CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((CH2O+10&
!)*(O2+0.1)*(X2+0.5)**2)
!PD(4,1) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)&
!*(X1+1.0))-3.555555555555555*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))
!PD(4,2) = 0
!PD(4,3) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)**2*(O2+0.1&
!)*(X2+0.5))-5.333333333333333*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))
!PD(4,4) = -5.333333333333333*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2&
!+0.5))+5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.&
!5))-3.555555555555555*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))+3.5555&
!55555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
!PD(4,5) = 0
!PD(4,6) = 0
!PD(4,7) = 0
!PD(4,8) = 0
!PD(4,9) = 0
!PD(4,10) = 0
!PD(4,11) = 0
!PD(4,12) = 3.555555555555555*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(&
!X1+1.0)**2)-3.555555555555555*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0)&
!)
!PD(4,13) = 5.333333333333333*CH2O*O2*X2/((CH2O+10)*(O2+0.1)*&
!(X2+0.5)**2)-5.333333333333333*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.&
!5))
!PD(5,1) = 0
!PD(5,2) = 3.3333333333333335e-5*CH2O*X2/((CH2O+10)*(NO3+0.5)&
!*(O2+1.0)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*&
!(NO3+0.5)**2*(O2+1.0)*(X2+0.5))
!PD(5,3) = 3.3333333333333335e-5*NO3*X2/((CH2O+10)*(NO3+0.5)*&
!(O2+1.0)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)**&
!2*(NO3+0.5)*(O2+1.0)*(X2+0.5))+1.666666666666667e-4*O2*X2/((CH2O+&
!10)*(O2+0.1)*(X2+0.5))-1.666666666666667e-4*CH2O*O2*X2/((CH2O+10)&
!**2*(O2+0.1)*(X2+0.5))
!PD(5,4) = -3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*(NO3&
!+0.5)*(O2+1.0)**2*(X2+0.5))+1.666666666666667e-4*CH2O*X2/((CH2O+1&
!0)*(O2+0.1)*(X2+0.5))-1.666666666666667e-4*CH2O*O2*X2/((CH2O+10)*&
!(O2+0.1)**2*(X2+0.5))
!PD(5,5) = -7.344e+8*OH-2592.016416
!PD(5,6) = 866.4*Ca*HCO3+6.818e+9*H+8.63
!PD(5,7) = 6.818e+9*HCO3
!PD(5,8) = 0
!PD(5,9) = 433.2*HCO3**2
!PD(5,10) = 0
!PD(5,11) = -7.344e+8*CO2
!PD(5,12) = 0
!PD(5,13) = 3.3333333333333335e-5*CH2O*NO3/((CH2O+10)*(NO3+0.&
!5)*(O2+1.0)*(X2+0.5))+1.666666666666667e-4*CH2O*O2/((CH2O+10)*(O2&
!+0.1)*(X2+0.5))-3.3333333333333335e-5*CH2O*NO3*X2/((CH2O+10)*(NO3&
!+0.5)*(O2+1.0)*(X2+0.5)**2)-1.666666666666667e-4*CH2O*O2*X2/((CH2&
!O+10)*(O2+0.1)*(X2+0.5)**2)
!PD(6,1) = 0
!PD(6,2) = 8.133333333333333*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2&
!+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)**2*(O2+1.0)*(X2+0.5))
!PD(6,3) = 8.133333333333333*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+&
!1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5))
!PD(6,4) = -8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5&
!)*(O2+1.0)**2*(X2+0.5))
!PD(6,5) = 7.344e+8*OH+2592.016416
!PD(6,6) = -9.22e+14*OH-866.4*Ca*HCO3-6.818e+9*H-0.29222*Ca-8&
!.63
!PD(6,7) = 29.222-6.818e+9*HCO3
!PD(6,8) = 1.11e+11
!PD(6,9) = -433.2*HCO3**2-0.29222*HCO3
!PD(6,10) = 0
!PD(6,11) = 7.344e+8*CO2-9.22e+14*HCO3
!PD(6,12) = 0
!PD(6,13) = 8.133333333333333*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(&
!O2+1.0)*(X2+0.5))-8.133333333333333*CH2O*NO3*X2/((CH2O+10)*(NO3+0&
!.5)*(O2+1.0)*(X2+0.5)**2)
!PD(7,1) = 2.2222222222222223e-4*O2*X1/((NH4+0.1)*(O2+0.1)*(X&
!1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X&
!1+1.0))
!PD(7,2) = 0
!PD(7,3) = 0
!PD(7,4) = 2.2222222222222223e-4*NH4*X1/((NH4+0.1)*(O2+0.1)*(&
!X1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(&
!X1+1.0))
!PD(7,5) = 2592.0
!PD(7,6) = 0.29222*Ca-6.818e+9*H
!PD(7,7) = -2.46e+16*OH-6.818e+9*HCO3-29.222
!PD(7,8) = 0
!PD(7,9) = 0.29222*HCO3
!PD(7,10) = 0
!PD(7,11) = -2.46e+16*H
!PD(7,12) = 2.2222222222222223e-4*NH4*O2/((NH4+0.1)*(O2+0.1)*&
!(X1+1.0))-2.2222222222222223e-4*NH4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1&
!+1.0)**2)
!PD(7,13) = 0
!PD(8,1) = 0
!PD(8,2) = 0
!PD(8,3) = 0
!PD(8,4) = 0
!PD(8,5) = 0
!PD(8,6) = 9.22e+14*OH
!PD(8,7) = 0
!PD(8,8) = -0.0057456*Ca-1.11e+11
!PD(8,9) = -0.0057456*CO3
!PD(8,10) = 0
!PD(8,11) = 9.22e+14*HCO3
!PD(8,12) = 0
!PD(8,13) = 0
!PD(9,1) = 0
!PD(9,2) = 0
!PD(9,3) = 0
!PD(9,4) = 0
!PD(9,5) = 0.016416
!PD(9,6) = -866.4*Ca*HCO3-0.29222*Ca
!PD(9,7) = 29.222
!PD(9,8) = -0.0057456*Ca
!PD(9,9) = -433.2*HCO3**2-0.29222*HCO3-0.0057456*CO3
!PD(9,10) = 0
!PD(9,11) = 0
!PD(9,12) = 0
!PD(9,13) = 0
!PD(10,1) = 0
!PD(10,2) = 1.866666666666667*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O&
!2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
!5)**2*(O2+1.0)*(X2+0.5))
!PD(10,3) = 1.866666666666667*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2&
!+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)**2*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5))
!PD(10,4) = -1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+0.&
!5)*(O2+1.0)**2*(X2+0.5))
!PD(10,5) = 0
!PD(10,6) = 0
!PD(10,7) = 0
!PD(10,8) = 0
!PD(10,9) = 0
!PD(10,10) = 0
!PD(10,11) = 0
!PD(10,12) = 0
!PD(10,13) = 1.866666666666667*CH2O*NO3/((CH2O+10)*(NO3+0.5)*&
!(O2+1.0)*(X2+0.5))-1.866666666666667*CH2O*NO3*X2/((CH2O+10)*(NO3+&
!0.5)*(O2+1.0)*(X2+0.5)**2)
!PD(11,1) = 0
!PD(11,2) = 0
!PD(11,3) = 0
!PD(11,4) = 0
!PD(11,5) = -7.344e+8*OH
!PD(11,6) = 8.63-9.22e+14*OH
!PD(11,7) = -2.46e+16*OH
!PD(11,8) = 1.11e+11
!PD(11,9) = 0
!PD(11,10) = 0
!PD(11,11) = -9.22e+14*HCO3-2.46e+16*H-7.344e+8*CO2
!PD(11,12) = 0
!PD(11,13) = 0
!PD(12,1) = 0.17*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH4&
!*O2*X1/((NH4+0.1)**2*(O2+0.1)*(X1+1.0))
!PD(12,2) = 0
!PD(12,3) = 0
!PD(12,4) = 0.17*NH4*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*NH&
!4*O2*X1/((NH4+0.1)*(O2+0.1)**2*(X1+1.0))
!PD(12,5) = 0
!PD(12,6) = 0
!PD(12,7) = 0
!PD(12,8) = 0
!PD(12,9) = 0
!PD(12,10) = 0
!PD(12,11) = 0
!PD(12,12) = 0.17*NH4*O2/((NH4+0.1)*(O2+0.1)*(X1+1.0))-0.17*N&
!H4*O2*X1/((NH4+0.1)*(O2+0.1)*(X1+1.0)**2)
!PD(12,13) = 0
!PD(13,1) = 0
!PD(13,2) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)**2*(O2&
!+1.0)*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5&
!)))
!PD(13,3) = -0.5*(-5.0*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(&
!X2+0.5))+5.0*CH2O*NO3*X2/((CH2O+10)**2*(NO3+0.5)*(O2+1.0)*(X2+0.5&
!))-5.0*O2*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*O2*X2/((CH2O+&
!10)**2*(O2+0.1)*(X2+0.5)))
!PD(13,4) = -0.5*(5.0*CH2O*NO3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.&
!0)**2*(X2+0.5))-5.0*CH2O*X2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2&
!O*O2*X2/((CH2O+10)*(O2+0.1)**2*(X2+0.5)))
!PD(13,5) = 0
!PD(13,6) = 0
!PD(13,7) = 0
!PD(13,8) = 0
!PD(13,9) = 0
!PD(13,10) = 0
!PD(13,11) = 0
!PD(13,12) = 0
!PD(13,13) = -0.5*(-5.0*CH2O*NO3/((CH2O+10)*(NO3+0.5)*(O2+1.0&
!)*(X2+0.5))-5.0*CH2O*O2/((CH2O+10)*(O2+0.1)*(X2+0.5))+5.0*CH2O*NO&
!3*X2/((CH2O+10)*(NO3+0.5)*(O2+1.0)*(X2+0.5)**2)+5.0*CH2O*O2*X2/((&
!CH2O+10)*(O2+0.1)*(X2+0.5)**2))

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
