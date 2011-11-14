

MODULE Chen1992

   IMPLICIT NONE
   REAL, PARAMETER :: MYPAR = 3.D7
   INTEGER, PARAMETER :: Chen1992TotalNumberOfSpecies = 5

   
   TYPE Chen1992Parameters
     CHARACTER(10), DIMENSION(Chen1992TotalNumberOfSpecies)::SpeciesNames
     REAL :: kmax_tol ! maximum rate of substrate utilization
     REAL :: kmax_ben ! maximum rate of substrate utilization
     REAL :: Ktol   ! toluene half-saturation constant
     REAL :: Kben   ! benzene half-saturation constant
     REAL :: Kdo   ! electron acceptor half-saturation constant
     REAL :: Rtol   ! substrate retardation
     REAL :: Rben   ! substrate retardation
     REAL :: Rm   ! electron acceptor retardation
     REAL :: Y    ! biomass yield coefficient
     REAL :: r_tol_do    ! Utilizaton ratio, X 
     REAL :: r_ben_do    ! Utilizaton ratio, X 
     REAL :: decay ! biomass decay coefficient
     REAL, DIMENSION(Chen1992TotalNumberOfSpecies) ::backgroundConc
   END TYPE Chen1992Parameters

   REAL, ALLOCATABLE :: Biomass1Conc(:,:,:), Biomass2Conc(:,:,:)

   TYPE( Chen1992Parameters ) :: cpars

   CONTAINS

     SUBROUTINE ReadChen1992TransPars( mparfile, nx, ny, nz )
      IMPLICIT NONE

      CHARACTER(LEN=100), INTENT(IN) :: mparfile
      INTEGER :: I, J, K, idx, nx, ny, nz
      OPEN(98, FILE = mparfile, STATUS = 'old' ) 

      DO I = 1, Chen1992TotalNumberOfSpecies
         READ(98,*)  idx, cpars%SpeciesNames( idx )
      END DO
     READ(98,*) cpars%kmax_tol
     READ(98,*) cpars%kmax_ben
     READ(98,*) cpars%Ktol
     READ(98,*) cpars%Kben
     READ(98,*) cpars%Kdo
     READ(98,*) cpars%Rtol
     READ(98,*) cpars%Rben
     READ(98,*) cpars%Rm
     READ(98,*) cpars%Y
     READ(98,*) cpars%r_tol_do ! Utilization ratio, X1
     READ(98,*) cpars%r_ben_do ! Utilization ratio, X2
     READ(98,*) cpars%decay
     DO I = 1, Chen1992TotalNumberOfSpecies
        READ(98,*)  cpars%backgroundConc( I )
     END DO
     CLOSE(98)
     ALLOCATE( Biomass1Conc(nx, ny, nz), Biomass2Conc(nx, ny, nz) )

     DO I = 1, nx
       DO J = 1, ny
         DO K = 1, nz
          Biomass1Conc( I, J, K ) = cpars%backgroundConc(4)
          Biomass2Conc( I, J, K ) = cpars%backgroundConc(5)
         END DO  
       END DO
     END DO
  
  END SUBROUTINE ReadChen1992TransPars

  REAL FUNCTION D_tol( tol, A, X )
   IMPLICIT NONE
   REAL :: tol, A, X 
   D_tol = -1.0 * X * ( cpars%kmax_tol / cpars%Rtol ) *  &
             ( tol / ( cpars%Ktol + tol ) ) * ( A / ( cpars%Kdo + A ) )
  END FUNCTION

  REAL FUNCTION D_ben( ben, A, X )
   IMPLICIT NONE
   REAL :: ben, A, X 
   D_ben = -1.0 * X * ( cpars%kmax_ben / cpars%Rben ) *  &
             ( ben / ( cpars%Kben + ben ) ) * ( A / ( cpars%Kdo + A ) )
  END FUNCTION

  REAL FUNCTION D_A( tol, ben , A, X1, X2 )
   IMPLICIT NONE
   REAL :: tol, ben, A, X1, X2 
   D_A = -1.0 * cpars%r_tol_do * X1 * ( cpars%kmax_tol / cpars%Rtol )  &
            * ( tol / ( cpars%Ktol + tol ) ) * ( A / ( cpars%Kdo + A ) ) &
         -1.0 * cpars%r_ben_do * X2 * ( cpars%kmax_ben / cpars%Rben ) *  &
             ( ben / ( cpars%Kben + ben ) ) * ( A / ( cpars%Kdo + A ) )
  END FUNCTION

  REAL FUNCTION D_X1( tol, A, X1 )
   IMPLICIT NONE
   REAL :: tol, A, X1 
   D_X1 = cpars%Y * X1 * cpars%kmax_tol *              &
         ( tol / ( cpars%Ktol + tol ) ) * ( A / ( cpars%Kdo + A ) ) -  &
          cpars%decay * X1
  END FUNCTION

  REAL FUNCTION D_X2( ben, A, X2 )
   IMPLICIT NONE
   REAL :: ben, A, X2 
   D_X2 = cpars%Y * X2 * cpars%kmax_ben *              &
         ( ben / ( cpars%Kben + ben ) ) * ( A / ( cpars%Kdo + A ) ) -  &
          cpars%decay * X2
  END FUNCTION
      
  SUBROUTINE Chen1992ReactRateAtCell( conc, xloc, yloc, zloc, dt, rate  )
   IMPLICIT NONE
      INTEGER NEQ, ITASK, IOPT, LRW, LIW, MF, IOUT, ISTATE, ITOL, IPAR, &
              IWORK, I
      REAL ATOL, RPAR, RTOL, RWORK, T, TOUT, Y, dt
      ! Num. of RWORK elements = 22 + 9 * NEQ + 2*NEQ**2 for MF = 21 or 22
      !                        = 22 + 9 * 5 + 2 * 5 ** 2 
      !                        = 117
      ! Num. of IWORK elements = 30 + NEQ  for MF = 21, 22 24 or 25
      !                        = 30 + 5
      !                        = 35
      DIMENSION Y(5), ATOL(5), RWORK(117), IWORK(35)
      REAL, DIMENSION(:) :: rate 
      INTEGER :: xloc, yloc, zloc
      REAL*4, DIMENSION(:,:,:,:) :: conc

      NEQ = Chen1992TotalNumberOfSpecies 

      Y(1) = conc( 1, xloc, yloc, zloc ) +        &
                        cpars%backgroundConc( 1 )     ! toluene
      Y(2) = conc( 2, xloc, yloc, zloc ) +        &
                        cpars%backgroundConc( 2 )     ! benzene
      Y(3) = conc( 3, xloc, yloc, zloc )  +       &
                        cpars%backgroundConc( 3 ) ! electron acceptor
   !   Y(3) = conc( 3, xloc, yloc, zloc ) +       &
   !                     cpars%backgroundConc( 3 ) ! biomass
      Y(4) = Biomass1Conc( xloc, yloc, zloc ) 
      Y(5) = Biomass2Conc( xloc, yloc, zloc ) 

      DO I = 1, Chen1992TotalNumberOfSpecies
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
      ITASK = 1
      ISTATE = 1
      IOPT = 1
      LRW = 117  ! 20 + NYH * (MAXORD + 1) + 3 * NEQ + LWM
                ! NYH = NEQ = 5
                ! MAXORD = 5 ( METH = 2 )
                ! NEQ = 5
                ! LWM = 2 * NEQ**2  + 2 ( MITER =1 and MF .gt. 0 )
                ! LRW = 20 + 5 * ( 5 + 1 ) + 3 * 5 + 2 * 5 ** 2 + 2
                !     = 117
      LIW = 35  ! 30 + NEQ
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
       CALL SVODE(Chen1992REACT,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,   &
                 IOPT,RWORK,LRW,IWORK,LIW,Chen1992JAC,MF,RPAR,IPAR)
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
      DO I = 1, Chen1992TotalNumberOfSpecies
         rate(I) = Y(I) - rate(I)
      END DO
      Biomass1Conc( xloc, yloc, zloc ) = Y(4)
      Biomass2Conc( xloc, yloc, zloc ) = Y(5)
      RETURN
  80  WRITE(6,90)ISTATE, xloc, yloc, zloc
      rate = 0.0 

  90  FORMAT(///' Error halt: ISTATE =',I3, ' Cell (', I4, ', ', I4, ', ', I4, ')' )
!      STOP
      RETURN

  END SUBROUTINE

  SUBROUTINE Chen1992REACT (NEQ, T, Y, YDOT, RPAR, IPAR)
     INTEGER  NEQ, IPAR
     REAL  RPAR, T, Y, YDOT
     DIMENSION Y(NEQ), YDOT(NEQ)
     YDOT(1) = D_tol( Y(1), Y(3), Y(4) )
     YDOT(2) = D_ben( Y(2), Y(3), Y(5) )
     YDOT(3) = D_A( Y(1), Y(2), Y(3), Y(4), Y(5) )
     YDOT(4) = D_X1( Y(1), Y(3), Y(4) )
     YDOT(5) = D_X2( Y(2), Y(3), Y(5) )

  END SUBROUTINE

  SUBROUTINE Chen1992JAC (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      INTEGER  NEQ, NRPD, IPAR, ML, MU
      REAL PD, RPAR, T, Y
      DIMENSION Y(NEQ), PD(NRPD,NEQ)
      REAL tol, ben, A, X1, X2
      tol = Y(1)
      ben = Y(2)
      A   = Y(3)
      X1  = Y(4)
      X2  = Y(5)

      PD(1,1) = cpars%kmax_tol*tol*A*X1/((cpars%Kdo+A)*(cpars%&
         &Ktol+tol)**2*cpars%Rtol)-cpars%kmax_tol*A*X1/((cpars%Kdo+A)*(cp&
         &ars%Ktol+tol)*cpars%Rtol)
      PD(1,2) = 0
      PD(1,3) = cpars%kmax_tol*tol*A*X1/((cpars%Kdo+A)**2*(cpa&
         &rs%Ktol+tol)*cpars%Rtol)-cpars%kmax_tol*tol*X1/((cpars%Kdo+A)*(&
         &cpars%Ktol+tol)*cpars%Rtol)
      PD(1,4) = -cpars%kmax_tol*tol*A/((cpars%Kdo+A)*(cpars%Kt&
         &ol+tol)*cpars%Rtol)
      PD(1,5) = 0
      PD(2,1) = 0
      PD(2,2) = ben*cpars%kmax_ben*A*X2/((cpars%Kben+ben)**2*(&
         &cpars%Kdo+A)*cpars%Rben)-cpars%kmax_ben*A*X2/((cpars%Kben+ben)*&
         &(cpars%Kdo+A)*cpars%Rben)
      PD(2,3) = ben*cpars%kmax_ben*A*X2/((cpars%Kben+ben)*(cpa&
         &rs%Kdo+A)**2*cpars%Rben)-ben*cpars%kmax_ben*X2/((cpars%Kben+ben&
         &)*(cpars%Kdo+A)*cpars%Rben)
      PD(2,4) = 0
      PD(2,5) = -ben*cpars%kmax_ben*A/((cpars%Kben+ben)*(cpars&
         &%Kdo+A)*cpars%Rben)
      PD(3,1) = cpars%kmax_tol*cpars%r_tol_do*tol*A*X1/((cpars&
         &%Kdo+A)*(cpars%Ktol+tol)**2*cpars%Rtol)-cpars%kmax_tol*cpars%r_&
         &tol_do*A*X1/((cpars%Kdo+A)*(cpars%Ktol+tol)*cpars%Rtol)
      PD(3,2) = ben*cpars%kmax_ben*cpars%r_ben_do*A*X2/((cpars%&
         &Kben+ben)**2*(cpars%Kdo+A)*cpars%Rben)-cpars%kmax_ben*cpars%r_&
         &ben_do*A*X2/((cpars%Kben+ben)*(cpars%Kdo+A)*cpars%Rben)
      PD(3,3) = -ben*cpars%kmax_ben*cpars%r_ben_do*X2/((cpars%&
         &Kben+ben)*(cpars%Kdo+A)*cpars%Rben)+ben*cpars%kmax_ben*cpars%r_&
         &ben_do*A*X2/((cpars%Kben+ben)*(cpars%Kdo+A)**2*cpars%Rben)-cpar&
         &s%kmax_tol*cpars%r_tol_do*tol*X1/((cpars%Kdo+A)*(cpars%Ktol+tol&
         &)*cpars%Rtol)+cpars%kmax_tol*cpars%r_tol_do*tol*A*X1/((cpars%Kd&
         &o+A)**2*(cpars%Ktol+tol)*cpars%Rtol)
      PD(3,4) = -cpars%kmax_tol*cpars%r_tol_do*tol*A/((cpars%K&
         &do+A)*(cpars%Ktol+tol)*cpars%Rtol)
      PD(3,5) = -ben*cpars%kmax_ben*cpars%r_ben_do*A/((cpars%K&
         &ben+ben)*(cpars%Kdo+A)*cpars%Rben)
      PD(4,1) = cpars%kmax_tol*A*cpars%Y*X1/((cpars%Kdo+A)*(cp&
         &ars%Ktol+tol))-cpars%kmax_tol*tol*A*cpars%Y*X1/((cpars%Kdo+A)*(&
         &cpars%Ktol+tol)**2)
      PD(4,2) = 0
      PD(4,3) = cpars%kmax_tol*tol*cpars%Y*X1/((cpars%Kdo+A)*(&
         &cpars%Ktol+tol))-cpars%kmax_tol*tol*A*cpars%Y*X1/((cpars%Kdo+A)&
         &**2*(cpars%Ktol+tol))
      PD(4,4) = cpars%kmax_tol*tol*A*cpars%Y/((cpars%Kdo+A)*(c&
         &pars%Ktol+tol))-cpars%decay
      PD(4,5) = 0
      PD(5,1) = 0
      PD(5,2) = cpars%kmax_ben*A*cpars%Y*X2/((cpars%Kben+ben)*&
         &(cpars%Kdo+A))-ben*cpars%kmax_ben*A*cpars%Y*X2/((cpars%Kben+ben&
         &)**2*(cpars%Kdo+A))
      PD(5,3) = ben*cpars%kmax_ben*cpars%Y*X2/((cpars%Kben+ben&
         &)*(cpars%Kdo+A))-ben*cpars%kmax_ben*A*cpars%Y*X2/((cpars%Kben+b&
         &en)*(cpars%Kdo+A)**2)
      PD(5,4) = 0
      PD(5,5) = ben*cpars%kmax_ben*A*cpars%Y/((cpars%Kben+ben)&
         &*(cpars%Kdo+A))-cpars%decay

  END SUBROUTINE

END MODULE Chen1992
