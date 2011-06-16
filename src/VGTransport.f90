

MODULE VGTransport

   IMPLICIT NONE
   REAL, PARAMETER :: MYPAR = 3.D7
   INTEGER, PARAMETER :: VGTotalNumberOfSpecies = 1

   
   TYPE VGTransParameters
     CHARACTER(10), DIMENSION(VGTotalNumberOfSpecies)::SpeciesNames
     REAL :: firstorderdecay_1_day 
     REAL :: zeroorderdecay_1_day
     REAL, DIMENSION(VGTotalNumberOfSpecies)::backgroundConc
   END TYPE VGTransParameters

   TYPE( VGTransParameters ) :: vgpars

   CONTAINS

     SUBROUTINE ReadVGTransPars( vgparfile, nx, ny, nz )
      IMPLICIT NONE

      CHARACTER(LEN=100), INTENT(IN) :: vgparfile
      INTEGER :: I, J, K, idx, nx, ny, nz
      OPEN(98, FILE = vgparfile, STATUS = 'old' ) 

      DO I = 1, VGTotalNumberOfSpecies
         READ(98,*)  idx, vgpars%SpeciesNames( idx )
      END DO
     READ(98,*) vgpars%firstorderdecay_1_day
     READ(98,*) vgpars%zeroorderdecay_1_day
     DO I = 1, VGTotalNumberOfSpecies
        READ(98,*)  vgpars%backgroundConc( I )
     END DO
     CLOSE(98)
  
  END SUBROUTINE ReadVGTransPars
      
  SUBROUTINE VGReactRateAtCell( conc, xloc, yloc, zloc, dt, rate  )
   IMPLICIT NONE
   INTEGER :: xloc, yloc, zloc
   REAL*4, DIMENSION(:,:,:,:) :: conc
   REAL, DIMENSION(:) :: rate
   REAL :: dt

!   rate(1) = ( - conc( 1, xloc, yloc, zloc ) *                        &
!                         vgpars%firstorderdecay_1_day    +            &
!                        vgpars%zeroorderdecay_1_day ) * dt
   rate(1) = conc( 1, xloc, yloc, zloc ) *      &
               ( EXP( -vgpars%firstorderdecay_1_day * dt ) - 1.0 )

  END SUBROUTINE

END MODULE VGTransport
