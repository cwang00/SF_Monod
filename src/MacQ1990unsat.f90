

MODULE MacQ1990unsat

   IMPLICIT NONE
   INTEGER, PARAMETER :: MacQ1990unsatTotalNumberOfSpecies = 3

   
   TYPE MacQ1990unsatParameters
     CHARACTER(10), DIMENSION(MacQ1990unsatTotalNumberOfSpecies)::SpeciesNames
     REAL :: kmax ! maximum rate of substrate utilization
     REAL :: Ks   ! substrate half-saturation constant
     REAL :: Ka   ! electron acceptor half-saturation constant
     REAL :: Rs   ! substrate retardation
     REAL :: Rm   ! electron acceptor retardation
     REAL :: Y    ! biomass yield coefficient
     REAL :: r1   ! substrate stoichiometric coefficient
     REAL :: r2   ! electron acceptor stoichiometric coefficient
     REAL :: decay ! biomass decay coefficient
     REAL :: Kb    ! biomass inhibition coefficient
     REAL :: SH    ! substrate Henry's constant
     REAL :: EH    ! electorn acceptor Henry's constant
     REAL :: SD0    ! substrate free air diffusion coefficient
     REAL :: ED0    ! electorn free air diffusion coefficient
!     REAL, DIMENSION(MacQ1990unsatTotalNumberOfSpecies)::backgroundConc
     REAL, ALLOCATABLE::backgroundConc(:,:,:,:)
   END TYPE MacQ1990unsatParameters

   REAL, ALLOCATABLE :: UBiomassConc(:,:,:)

   TYPE( MacQ1990unsatParameters ) :: umpars

   CONTAINS

     SUBROUTINE ReadMacQ1990unsatTransPars( mparfile, nx, ny, nz )
      IMPLICIT NONE

      CHARACTER(LEN=100), INTENT(IN) :: mparfile
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
      OPEN(98, FILE = mparfile, STATUS = 'old' ) 

      DO I = 1, MacQ1990unsatTotalNumberOfSpecies
         READ(98,*)  idx, umpars%SpeciesNames( idx )
      END DO
     READ(98,*) umpars%kmax
     READ(98,*) umpars%Ks
     READ(98,*) umpars%Ka
     READ(98,*) umpars%Rs
     READ(98,*) umpars%Rm
     READ(98,*) umpars%Y
     READ(98,*) umpars%r1
     READ(98,*) umpars%r2
     READ(98,*) umpars%decay
     READ(98,*) umpars%Kb
     READ(98,*) umpars%SH
     READ(98,*) umpars%EH
     READ(98,*) umpars%SD0
     READ(98,*) umpars%ED0

     ALLOCATE( umpars%backgroundConc(MacQ1990unsatTotalNumberOfSpecies,nx, ny, nz ) )

     DO I = 1, MacQ1990unsatTotalNumberOfSpecies
        READ(98,*) bgConcType
        IF ( TRIM( bgConcType  ) .eq. 'const' ) THEN
           READ(98, *) bgConcValue
           umpars%backgroundConc(I, :, :, : ) = bgConcValue
        ELSE IF ( TRIM( bgConcType ) .eq. 'PFBFile' ) THEN
           READ(98, *) bgConcPFBFile
           CALL PF_READ( temp, bgConcPFBFile, &
                        nx1, ny1, nz1, &
                           dx, dy, dz )
            DO L = 1, nx
             DO M = 1, ny
              DO N = 1, nz
               umpars%backgroundConc(I,L,M,N) = temp(L,M,N)
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
     ALLOCATE( UBiomassConc(nx, ny, nz) )

     DO I = 1, nx
       DO J = 1, ny
         DO K = 1, nz
          UBiomassConc( I, J, K ) = umpars%backgroundConc(3, I,J,K)
         END DO  
       END DO
     END DO
  
  END SUBROUTINE ReadMacQ1990unsatTransPars

  REAL FUNCTION F_X( X )
   IMPLICIT NONE
   REAL :: X
   F_X = umpars%Kb / ( umpars%Kb + X )
  END FUNCTION

  REAL FUNCTION D_S( S, A, X )
   IMPLICIT NONE
   REAL :: S, A, X 
   D_S = -1.0 * umpars%r1 * X * ( umpars%kmax / umpars%Rs ) *  &
            umpars%Kb / ( umpars%Kb + X ) * ( S / ( umpars%Ks + S ) ) * ( A / ( umpars%Ka + A ) )
  END FUNCTION

  REAL FUNCTION D_A( S, A, X )
   IMPLICIT NONE
   REAL :: S, A, X 
   D_A = -1.0 * umpars%r2 * X * ( umpars%kmax / umpars%Rm )  &
            * umpars%Kb / ( umpars%Kb + X ) * ( S / ( umpars%Ks + S ) ) * ( A / ( umpars%Ka + A ) )
  END FUNCTION

  REAL FUNCTION D_X( S, A, X )
   IMPLICIT NONE
   REAL :: S, A, X 
   D_X = umpars%Y * X * umpars%kmax * umpars%Kb / ( umpars%Kb + X ) *             &
         ( S / ( umpars%Ks + S ) ) * ( A / ( umpars%Ka + A ) ) -  &
          umpars%decay * X
  END FUNCTION
      
  SUBROUTINE MacQ1990unsatReactRateAtCell( conc, xloc, yloc, zloc, dt, rate  )
   IMPLICIT NONE
      INTEGER NEQ, ITASK, IOPT, LRW, LIW, MF, IOUT, ISTATE, ITOL, IPAR, &
              IWORK, I
      REAL ATOL, RPAR, RTOL, RWORK, T, TOUT, Y, dt
      DIMENSION Y(3), ATOL(3), RWORK(67), IWORK(33)
      REAL, DIMENSION(:) :: rate 
      INTEGER :: xloc, yloc, zloc
      REAL*4, DIMENSION(:,:,:,:) :: conc

      NEQ = MacQ1990unsatTotalNumberOfSpecies 

      Y(1) = conc( 1, xloc, yloc, zloc ) ! +        &
  !                       umpars%backgroundConc( 1 )     ! Substrate
      Y(2) = conc( 2, xloc, yloc, zloc ) ! +       &
  !                       umpars%backgroundConc( 2 ) ! electron acceptor
   !   Y(3) = conc( 3, xloc, yloc, zloc ) +       &
   !                     umpars%backgroundConc( 3 ) ! biomass
      Y(3) = UBiomassConc( xloc, yloc, zloc ) 

      DO I = 1, MacQ1990unsatTotalNumberOfSpecies
       rate(I) = Y(I)
      END DO

      T = 0.0D0
      TOUT = dt
      ITOL = 2
      RTOL = 1.D-4
      ATOL(1) = 1.D-6
      ATOL(2) = 1.D-6
      ATOL(3) = 1.D-6
      ITASK = 1
      ISTATE = 1
      IOPT = 1
      LRW = 67  ! 20 + NYH * (MAXORD + 1) + 3 * NEQ + LWM
                ! NYH = NEQ = 3
                ! MAXORD = 5 ( METH = 2 )
                ! NEQ = 3
                ! LWM = 2 * NEQ**2  + 2 ( MITER =1 and MF .gt. 0 )
                ! LRW = 20 + 3 * ( 5 + 1 ) + 3 * 3 + 2 * 3 ** 2 + 2
                !     = 67
      LIW = 33  ! 30 + NEQ
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

      DO 
       CALL SVODE(MACQ1990REACT,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,   &
                 IOPT,RWORK,LRW,IWORK,LIW,MACQ1990JAC,MF,RPAR,IPAR)
!      WRITE(6,20)T,Y(1),Y(2),Y(3)
!  20  FORMAT(' At t =',D12.4,'   y =',3D14.6)
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
      DO I = 1, MacQ1990unsatTotalNumberOfSpecies
         rate(I) = Y(I) - rate(I)
      END DO
      UBiomassConc( xloc, yloc, zloc ) = Y(3)
      RETURN
  80  WRITE(6,90)ISTATE, xloc, yloc, zloc
      rate = 0.0 

  90  FORMAT(///' Error halt: ISTATE =',I3, ' Cell (', I4, ', ', I4, ', ', I4, ')' )
!      STOP
      RETURN

  END SUBROUTINE

  SUBROUTINE MACQ1990REACT (NEQ, T, Y, YDOT, RPAR, IPAR)
     INTEGER  NEQ, IPAR
     REAL  RPAR, T, Y, YDOT
     DIMENSION Y(NEQ), YDOT(NEQ)
     YDOT(1) = D_S( Y(1), Y(2), Y(3) )
     YDOT(2) = D_A( Y(1), Y(2), Y(3) )
     YDOT(3) = D_X( Y(1), Y(2), Y(3) )
  END SUBROUTINE

  SUBROUTINE MACQ1990JAC (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      INTEGER  NEQ, NRPD, IPAR, ML, MU
      REAL PD, RPAR, T, Y
      DIMENSION Y(NEQ), PD(NRPD,NEQ)
      REAL S, A, X
      S = Y(1)
      A = Y(2)
      X = Y(3)
 
      PD(1,1) = umpars%kmax*umpars%r1*A*umpars%Kb*S*X/((umpars%Ka            &
        +A)*umpars%Rs*(S+umpars%Ks)**2*(X+umpars%Kb))-umpars%kmax*umpars%r1*  &
        A*umpars%Kb*X/((umpars%Ka+A)*umpars%Rs*(S+umpars%Ks)*(X+umpars%Kb))
      PD(1,2) = umpars%kmax*umpars%r1*A*umpars%Kb*S*X/((umpars%Ka            &
        +A)**2*umpars%Rs*(S+umpars%Ks)*(X+umpars%Kb))-umpars%kmax*umpars%r1* &
        umpars%Kb*S*X/((umpars%Ka+A)*umpars%Rs*(S+umpars%Ks)*(X+umpars%Kb))
      PD(1,3) = umpars%kmax*umpars%r1*A*umpars%Kb*S*X/((umpars%Ka            &
        +A)*umpars%Rs*(S+umpars%Ks)*(X+umpars%Kb)**2)-umpars%kmax*umpars%r1* &
        A*umpars%Kb*S/((umpars%Ka+A)*umpars%Rs*(S+umpars%Ks)*(X+umpars%Kb))
      PD(2,1) = umpars%kmax*umpars%r2*A*umpars%Kb*S*X/((umpars%Ka            &
        +A)*umpars%Rm*(S+umpars%Ks)**2*(X+umpars%Kb))-umpars%kmax*umpars%r2*  &
        A*umpars%Kb*X/((umpars%Ka+A)*umpars%Rm*(S+umpars%Ks)*(X+umpars%Kb))
      PD(2,2) = umpars%kmax*umpars%r2*A*umpars%Kb*S*X/((umpars%Ka            &
        +A)**2*umpars%Rm*(S+umpars%Ks)*(X+umpars%Kb))-umpars%kmax*umpars%r2*  &
        umpars%Kb*S*X/((umpars%Ka+A)*umpars%Rm*(S+umpars%Ks)*(X+umpars%Kb))
      PD(2,3) = umpars%kmax*umpars%r2*A*umpars%Kb*S*X/((umpars%Ka            &
        +A)*umpars%Rm*(S+umpars%Ks)*(X+umpars%Kb)**2)-umpars%kmax*umpars%r2*  &
        A*umpars%Kb*S/((umpars%Ka+A)*umpars%Rm*(S+umpars%Ks)*(X+umpars%Kb))
      PD(3,1) = umpars%kmax*A*umpars%Kb*umpars%Y*X/((umpars%Ka+A)            &
        *(S+umpars%Ks)*(X+umpars%Kb))-umpars%kmax*A*umpars%Kb*umpars%Y*S*X/(  &
        (umpars%Ka+A)*(S+umpars%Ks)**2*(X+umpars%Kb))
      PD(3,2) = umpars%kmax*umpars%Kb*umpars%Y*S*X/((umpars%Ka+A)            &
        *(S+umpars%Ks)*(X+umpars%Kb))-umpars%kmax*A*umpars%Kb*umpars%Y*S*X/( &
        (umpars%Ka+A)**2*(S+umpars%Ks)*(X+umpars%Kb))
      PD(3,3) = umpars%kmax*A*umpars%Kb*umpars%Y*S/((umpars%Ka+A)           &
        *(S+umpars%Ks)*(X+umpars%Kb))-umpars%kmax*A*umpars%Kb*umpars%Y*S*X/( &
        (umpars%Ka+A)*(S+umpars%Ks)*(X+umpars%Kb)**2)-umpars%decay
  END SUBROUTINE

END MODULE MacQ1990unsat
