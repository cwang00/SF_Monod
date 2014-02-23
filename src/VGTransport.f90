

MODULE VGTransport

   IMPLICIT NONE
   REAL, PARAMETER :: MYPAR = 3.D7
   INTEGER, PARAMETER :: VGTotalNumberOfSpecies = 1

   
   TYPE VGTransParameters
     CHARACTER(10), DIMENSION(VGTotalNumberOfSpecies)::SpeciesNames
     REAL :: firstorderdecay_1_day 
     REAL :: zeroorderdecay_1_day
!     REAL, DIMENSION(VGTotalNumberOfSpecies)::backgroundConc
     REAL, ALLOCATABLE::backgroundConc(:,:,:,:)
   END TYPE VGTransParameters

   TYPE( VGTransParameters ) :: vgpars

   CONTAINS

     SUBROUTINE ReadVGTransPars( vgparfile, nx, ny, nz )
      IMPLICIT NONE

      CHARACTER(LEN=100), INTENT(IN) :: vgparfile
      CHARACTER(LEN=128) :: bgConcType
      CHARACTER(LEN=500) :: bgConcPFBFile
      INTEGER :: I, J, K,L,M,N, idx, nx, ny, nz
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
      OPEN(98, FILE = vgparfile, STATUS = 'old' ) 

      DO I = 1, VGTotalNumberOfSpecies
         READ(98,*)  idx, vgpars%SpeciesNames( idx )
      END DO
     READ(98,*) vgpars%firstorderdecay_1_day
     READ(98,*) vgpars%zeroorderdecay_1_day

     ALLOCATE( vgpars%backgroundConc(VGTotalNumberOfSpecies,nx, ny, nz ) )

     DO I = 1, VGTotalNumberOfSpecies
        READ(98,*) bgConcType
        IF ( TRIM( bgConcType  ) .eq. 'const' ) THEN
           READ(98, *) bgConcValue
           vgpars%backgroundConc(I, :, :, : ) = bgConcValue
        ELSE IF ( TRIM( bgConcType ) .eq. 'PFBFile' ) THEN
           READ(98, *) bgConcPFBFile
           CALL PF_READ( temp, bgConcPFBFile, &
                      nx1, ny1, nz1, &
                           dx, dy, dz )
            DO L = 1, nx
             DO M = 1, ny
              DO N = 1, nz
               vgpars%backgroundConc(I,L,M,N) = temp(L,M,N)
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
