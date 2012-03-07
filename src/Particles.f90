

MODULE Particles

   USE NTransport
   USE VGTransport
   USE MacQ1990Transport
   USE MacQ1990unsat
   USE Chen1992

   IMPLICIT NONE

   INTERFACE
         REAL*8 FUNCTION ran1( idum )
           INTEGER*4 idum
         END FUNCTION ran1
   END INTERFACE

   TYPE :: Ptr_to_array
   
      INTEGER*4, DIMENSION(:), POINTER :: arr 

   END TYPE Ptr_to_array

   TYPE ( Ptr_to_array ), DIMENSION(:,:,:,:), ALLOCATABLE :: part_numbers 
   INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: number_of_parts
   INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: max_number_of_parts
!   REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmpp, tmplastprint
!   INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmpipwell
!   INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: tmpip, tmpiprp
!   INTEGER*4, DIMENSION(:), ALLOCATABLE :: tmpirp
!   INTEGER, DIMENSION(:), ALLOCATABLE :: partNumber
   INTEGER*4, DIMENSION(:), ALLOCATABLE :: removed
   INTEGER totalremoved
   REAL*8, DIMENSION(TotalNumberOfSpecies) :: maxParticleMass

   CONTAINS

     SUBROUTINE countParticles( pp,ip, totalparticles, nx, ny, nz, ns, &
                                 delv, tnext )

!       TYPE( Ptr_to_array ), DIMENSION(:,:,:,:) :: part_num
!       INTEGER, DIMENSION(:,:,:,:) :: number_of_parts
       REAL*8, DIMENSION(:,:) :: pp
       INTEGER*4, DIMENSION(:,:) :: ip
       INTEGER*4 :: nx, ny, nz, ns, totalparticles, i,j,k,l, p, x, y, z
       REAL*8, DIMENSION(3) :: delv 
       INTEGER*4, DIMENSION(ns, nx, ny, nz) :: pcount
       REAL*8 :: tnext

       DO l = 1, ns
         DO i = 1, nx
           DO j = 1, ny
             DO k = 1, nz
              number_of_parts(l, i, j, k ) = 0
              pcount(l, i, j, k ) = 0
!              NULLIFY( part_numbers( l, i, j, k)%arr )
             END DO
           END DO
         END DO
       END DO

       !count number of particles for each cell for each species
       DO p = 1, totalparticles
         x = IDINT( pp( p, 1 ) / delv( 1 ) ) + 1
         y = IDINT( pp( p, 2 ) / delv( 2 ) ) + 1
         z = IDINT( pp( p, 3 ) / delv( 3 ) ) + 1
         IF ( ( ip(p,2) .EQ. 1 ) .AND. ( x .GT. 0 ) .AND. ( y .GT. 0 ) .AND. &
              ( z .GT. 0 )  .AND. ( x .LE. nx ) .AND. ( y .LE. ny ) .AND.    &
              ( z .LE. nz ) .AND. ( pp( p, 7 ) .LE. tnext ) ) THEN
           number_of_parts( ip( p, 1 ), x, y, z ) =                           &
              number_of_parts( ip( p, 1 ), x, y, z ) + 1
         ENDIF
       END DO

       ! allocate space
       DO l = 1, ns
         DO i = 1, nx
          DO j = 1, ny
             DO k = 1, nz
              IF( ALLOCATED( part_numbers( l, i, j, k )%arr )) THEN
                   DEALLOCATE( part_numbers( l, i, j, k)%arr )              
                    NULLIFY( part_numbers( l, i, j, k)%arr )
               ENDIF 
               IF ( number_of_parts( l, i,j,k) > 0 ) THEN
                 ALLOCATE(part_numbers(l, i, j, k)%arr(                      &
                                          number_of_parts(l, i, j, k)))
              ENDIF
             END DO
           END DO
         END DO
       END DO

       DO p = 1, totalparticles
         x = IDINT( pp( p, 1 ) / delv( 1 ) ) + 1
         y = IDINT( pp( p, 2 ) / delv( 2 ) ) + 1
         z = IDINT( pp( p, 3 ) / delv( 3 ) ) + 1
         IF ( ( ip(p,2) .EQ. 1 ) .AND. ( x .GT. 0 ) .AND. ( y .GT. 0 ) .AND. &
             ( z .GT. 0 ) .AND. ( x .LE. nx ) .AND. ( y .LE. ny ) .AND.      &
             ( z .LE. nz ) .AND. ( pp( p, 7 ) .LE. tnext ) ) THEN
           pcount( ip(p, 1), x, y, z ) = pcount( ip(p, 1), x, y, z ) + 1
           part_numbers( ip(p, 1), x, y, z )%arr( pcount( ip(p, 1), x, y,z ))&
                                                                         = p
         ENDIF
         !
         ! recycle particles
         IF ( ( ip(p,2) .EQ. 1 ) .AND. (                              &
             ( x .LE. 0 .OR. y .LE. 0 .OR. z .LE.  0 ) ) .AND.        &
              ( pp(p, 7) .LE. tnext ) ) THEN
           IF( ALL(removed(1:totalremoved) .NE. p ) ) THEN
             totalremoved = totalremoved + 1
             removed( totalremoved ) = p 
           ENDIF
         ENDIF
       END DO
     END SUBROUTINE countParticles


     SUBROUTINE countCellPart( n, ploc,pp, ip, nx, ny, nz, tnext, &
                                part_dens )

       REAL*8, DIMENSION(:,:) :: pp
       INTEGER*4, DIMENSION(:,:) :: ip
       INTEGER*4 :: nx, ny, nz, ns, n
       INTEGER*4, DIMENSION(:) :: ploc, part_dens
       REAL*8 :: tnext
       INTEGER, DIMENSION(:), ALLOCATABLE :: temp

       IF ( ( ip(n,2) .EQ. 1 ) .AND. ( ploc(1) .GT. 0 ) .AND.     &
                 ( ploc(2) .GT. 0 ) .AND. ( ploc(3)  .GT. 0 )     &
                  .AND. ( ploc(1) .LE. nx ) .AND. ( ploc(2) .LE. ny ) &
                  .AND. ( ploc(3) .LE. nz ) .AND.                     &
                  ( pp( n, 7 ) .LE. tnext ) ) THEN

           number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) =     &
             number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) + 1

         IF ( number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) .GT.  &
            max_number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) )  ) THEN

            max_number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) = &
             max_number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) + &
                  part_dens( ip(n,1) )

              ALLOCATE( temp(                                              &
                number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) - 1 ) )
                temp = part_numbers( ip(n, 1), ploc(1), ploc(2), ploc(3) )%arr
              DEALLOCATE(                                                   &
                   part_numbers( ip(n, 1), ploc(1), ploc(2), ploc(3) )%arr )
              ALLOCATE(                                                     &
                   part_numbers( ip(n, 1), ploc(1), ploc(2), ploc(3) )%arr( &
                 max_number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) &
                     ) )
              part_numbers( ip(n, 1), ploc(1), ploc(2), ploc(3) )%arr(   &
              1:number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) - 1 ) &
                = temp 
              DEALLOCATE( temp )
         END IF

           part_numbers( ip(n, 1), ploc(1), ploc(2), ploc(3) )%arr(    &
             number_of_parts( ip(n, 1), ploc(1), ploc(2), ploc(3) ) ) = n
             
        END IF

         !
         ! recycle particles
         IF ( ( ip(n,2) .EQ. 1 ) .AND. (                              &
             ( ploc(1) .LE. 0 .OR. ploc(2) .LE. 0 .OR. ploc(3) .LE. 0 ) ) &
              .AND. ( pp(n, 7) .LE. tnext ) ) THEN
           IF( ALL(removed(1:totalremoved) .NE. n ) ) THEN
             totalremoved = totalremoved + 1
             removed( totalremoved ) = n 
           ENDIF
         ENDIF

     END SUBROUTINE countCellPart

     SUBROUTINE allocateParticles_Memory(  nx, ny, nz, ns, npmax, part_dens )
        INTEGER*4 :: nx, ny, nz, ns, npmax, i, j, k, l
        INTEGER*4, DIMENSION(:) ::  part_dens(:)
        ALLOCATE( number_of_parts( ns, nx, ny, nz ) )
        ALLOCATE( max_number_of_parts( ns, nx, ny, nz ) )
!       ALLOCATE( tmpp( npmax, 10 ) )
!       ALLOCATE( tmplastprint( npmax, 3 ) )
!       ALLOCATE( tmpipwell( npmax, 20 ) )
!       ALLOCATE( tmpip( npmax, 10 ) )
!       ALLOCATE( tmpirp( npmax ) )
!       ALLOCATE( tmpiprp( npmax, 2 ) )
       ALLOCATE( removed( npmax ) )
!       ALLOCATE( partNumber(npmax) )
       ALLOCATE( part_numbers( ns, nx, ny, nz ) )
       DO I = 1, nx
         DO J = 1, ny
           DO K = 1, nz
             DO l = 1, ns
              max_number_of_parts(l, I, J, K) = part_dens(l)
              ALLOCATE( part_numbers( l, I, J, K )%arr( part_dens(l) ) )
             END DO
           END DO
         END DO
       END DO
       
       totalremoved = 0
!       rmcnter1 = 0
!       rmcnter2 = 0
     END SUBROUTINE allocateParticles_Memory

     SUBROUTINE deallocateParticles_Memory(  nx, ny, nz, ns )

        INTEGER*4 :: nx, ny, nz, ns, l, i, j, k

!         IF( ALLOCATED( number_of_parts ) ) THEN
             DEALLOCATE( number_of_parts )
             DEALLOCATE( max_number_of_parts )
!       DEALLOCATE( tmpp )
!       DEALLOCATE( tmplastprint )
!       DEALLOCATE( tmpipwell )
!       DEALLOCATE( tmpip )
!       DEALLOCATE( tmpirp )
!       DEALLOCATE( tmpiprp )
       DEALLOCATE( removed )
!       DEALLOCATE( partNumber )
!         ENDIF 

       
       IF( ALLOCATED( part_numbers ) ) THEN
         DO l = 1, ns
           DO i = 1, nx
             DO j = 1, ny
               DO k = 1, nz
                 IF( ALLOCATED( part_numbers( l, i, j, k )%arr )) THEN
                    DEALLOCATE( part_numbers( l, i, j, k)%arr )              
                 ENDIF 
              END DO
             END DO
           END DO
         END DO
         DEALLOCATE( part_numbers )
       ENDIF

     END SUBROUTINE deallocateParticles_Memory


     SUBROUTINE addRemoveParticles( pp, ip, ipwell, irp, iprp, lastprint, &
                              np, npmax, delv,  ns, nx, ny, nz, conc, dt_day, &
                              porosity, sat, modelname  )
       REAL*8, DIMENSION(:,:) :: pp, lastprint
       REAL*8, DIMENSION(:,:,:) :: porosity, sat
       INTEGER, DIMENSION(:,:) :: ipwell
       INTEGER*4,  DIMENSION(:,:) :: ip, iprp
       INTEGER*4,  DIMENSION(:) ::  irp
       REAL*8, DIMENSION(3) :: delv
       INTEGER :: ns, nx, ny, nz, i, j, k, l, m, n, idx, p, ns_moving
       REAL*4, DIMENSION(:,:,:,:) :: conc
       REAL, DIMENSION(13) :: rates
       REAL*4 :: rate, minmass,maxmass, mass, backgroundmass, particlemass
       REAL*8 :: dt_day
       INTEGER*4 ::  npmax, np
       CHARACTER*20 substratename, modelname
       LOGICAL :: found, recycle
       REAL*4 :: bconc

       IF ( modelname == 'Chen1992' ) THEN
               ns_moving = Chen1992TotalNumberOfSpecies - 2
       ELSE IF( modelname == 'MacQ1990' ) THEN
               ns_moving = MacQ1990TotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ1990unsat' ) THEN
               ns_moving = MacQ1990unsatTotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ' ) THEN
               ns_moving = TotalNumberOfSpecies 
       ELSE IF( modelname == 'VG' ) THEN
               ns_moving = VGTotalNumberOfSpecies 
       ELSE IF( modelname == 'noreact' ) THEN
               ns_moving = 0
       ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
       ENDIF


    DO i = 3, nx - 2
      DO j = 3, ny - 2
        DO k = 3, nz - 2
        
          IF ( SUM( number_of_parts(1:ns_moving,i,j,k) ) .GT. 0 ) THEN 
            IF ( modelname == 'Chen1992' ) THEN
              CALL Chen1992ReactRateAtCell( conc, i, j, k, REAL(dt_day), rates )
            ELSE IF( modelname == 'MacQ1990' ) THEN
              CALL MacQ1990ReactRateAtCell( conc, i, j, k, REAL(dt_day), rates )
            ELSE IF( modelname == 'MacQ1990unsat' ) THEN
              CALL MacQ1990unsatReactRateAtCell( conc, i, j, k, REAL(dt_day), rates )
            ELSE IF( modelname == 'MacQ' ) THEN
              IF ( number_of_parts( 1, i,j,k ) .GT. 0 .AND.    &
                   number_of_parts( 2, i,j,k ) .GT. 0 .AND.    &
                   number_of_parts( 3, i,j,k ) .GT. 0 ) THEN
                CALL VODEReactionRatesAtCell(                          &
                      conc, i, j, k, REAL(dt_day), rates )
              ELSE
                 rates = 0.0
              ENDIF
            ELSE IF( modelname == 'VG' ) THEN
              CALL VGReactRateAtCell( conc, i, j, k, REAL(dt_day), rates )
            ELSE IF( modelname == 'noreact' ) THEN
                 rates = 0.0
            ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
            ENDIF

!            DO l = 1, ns
            DO l = 1, ns_moving 
!              CALL partNumbersAtCellSpec( partNumber,                        &
!                                         l, i, j, k, pp, ip, np,             &
!                                         number_of_parts(l,i,j,k), delv )

!              IF( conc( l,i,j,k) .LT. 0 .AND. conc(l,i,j,k) +           &
!                 npars%backgroundConc(l) .LT. 0 ) THEN
!                 WRITE(*,*) l, i, j, k
!              ENDIF 

!              rate = ReactRateAtCell( conc, i, j, k, npars%SpeciesNames( l ) )
!              IF ( TRIM(npars%SpeciesNames( l) ) .eq. 'NH4' ) THEN
!                 rate = rates( 1 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'NO3' ) THEN
!                 rate = rates( 2 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'CH2O' ) THEN
!                 rate = rates( 3 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'O2' ) THEN
!                 rate = rates( 4 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'CO2' ) THEN
!                 rate = rates( 5 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'HCO3' ) THEN
!                 rate = rates( 6 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'H' ) THEN
!                 rate = rates( 7 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'CO3' ) THEN
!                 rate = rates( 8 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'Ca2' ) THEN
!                 rate = rates( 9 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'N2' ) THEN
!                 rate = rates( 10 )
!              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'OH' ) THEN
!                 rate = rates( 11 )
!              ELSE IF ( TRIM(vgpars%SpeciesNames( l ) ) .eq. 'VG' ) THEN
!                 rate = rates( 1 )
!              ELSE
!
!                 PRINT*, 'Unknow species naem ', npars%SpeciesNames(l)
!
!              ENDIF
               rate = rates( l )
              mass = rate * delv(1) * delv(2) * delv(3) * 1000.0           &
                        * sat( i,j , k)  * porosity( i,j,k)
               ! update concentrations
               conc( l, i, j, k)  = conc(l, i, j, k ) + rate
!               IF ( conc( l, i, j, k ) .LT. 0.0 ) conc( l, i, j, k ) = 0.0
              ! 
              !gases in unsaturated soil
              !
             IF( modelname == 'MacQ' ) THEN
              IF ( TRIM(npars%SpeciesNames( l) ) .eq. 'N2' ) THEN
                mass = mass * sat(i,j,k) /                                 & 
                ( npars%N2HenryConst * ( 1 - sat( i,j,k) ) + sat(i,j,k) )
              ELSE IF ( TRIM(npars%SpeciesNames( l) ) .eq. 'CO2' ) THEN
                mass = mass * sat(i,j,k) /                                 & 
                ( npars%CO2HenryConst * ( 1 - sat( i,j,k) ) + sat(i,j,k) )
             ELSE IF ( TRIM(npars%SpeciesNames( l) ) .eq. 'O2' ) THEN
                mass = mass * sat(i,j,k) /                                 & 
                ( npars%O2HenryConst * ( 1 - sat( i,j,k) ) + sat(i,j,k) )
              ENDIF
             ENDIF

             IF ( modelname =='MacQ1990unsat' ) THEN
              IF ( l .eq. 1 ) THEN
                mass = mass * sat(i,j,k) /                                 &
                ( umpars%SH * ( 1 - sat( i,j,k) ) + sat(i,j,k) )
              ELSE IF ( l .eq. 2 ) THEN
                mass = mass * sat(i,j,k) /                                 &
                ( umpars%EH * ( 1 - sat( i,j,k) ) + sat(i,j,k) )
              ENDIF
             ENDIF

              IF ( modelname == 'Chen1992' ) THEN
               bconc = cpars%backgroundConc(l)
              ELSE IF( modelname == 'MacQ1990' ) THEN
               bconc = mpars%backgroundConc(l )
              ELSE IF( modelname == 'MacQ1990unsat' ) THEN
               bconc = umpars%backgroundConc(l )
              ELSE IF( modelname == 'MacQ' ) THEN
               bconc = npars%backgroundConc(l )
              ELSE IF( modelname == 'VG' ) THEN
               bconc = vgpars%backgroundConc(l)
              ELSE IF( modelname == 'noreact' ) THEN
                      bconc = 0.0
              ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
              ENDIF
              backgroundmass = bconc                                        &
                           * delv(1) * delv(2) * delv(3) * 1000.0           &
                        * sat( i,j , k)  * porosity( i,j,k)

              IF( number_of_parts( l, i, j, k ) == 0 .AND. mass > 0.0) THEN
!                minmass = mass / 10.0
                maxmass = maxParticleMass( l )
                IF (maxmass == 0.0 ) THEN
!                   maxmass = mass / 10.0
                    maxmass = MAXVAL( maxParticleMass )
                ENDIF

                DO WHILE( mass > 0.0 )
!                 mass = mass - minmass
                 mass = mass - maxmass
!                 IF ( rmcnter2 .GT. rmcnter1 ) THEN
!                      idx = removed( rmcnter2 )
!                      rmcnter2 = rmcnter2 - 1
!                      totalremoved = totalremoved - 1
                 IF ( totalremoved .GT. 0 ) THEN
                      idx = removed( totalremoved )
                      totalremoved = totalremoved - 1
                 ELSE 
                   idx = np + 1
                   np = np + 1
                 ENDIF
                 IF ( np <= npmax ) THEN
                  pp(idx,1) = ( i - 1 ) * delv( 1 ) + delv( 1 ) * ran1( 0 )
                  pp(idx,2) = ( j - 1 ) * delv( 2 ) + delv( 2 ) * ran1( 0 )
                  pp(idx,3) = ( k - 1 ) * delv( 3 ) + delv( 3 ) * ran1( 0 )
!                  pp(idx,4) = minmass
                  IF ( mass > 0 ) THEN
                    pp(idx,4) = maxmass
                  ELSE
                    pp(idx,4) = mass + maxmass 
                  ENDIF
                  pp(idx,5) = DBLE(idx)
                  pp(idx,7) = timeCellAllSpec( i, j, k, pp, modelname ) 
!                  pp(idx,7) = timeCellAllSpec( i, j, k, pp, np, ip, delv ) 
                  ip(idx,1) = l
                  ip(idx,2) = 1

                 ELSE
                   PRINT*, '2 max num particles exceeded. &
                           &increase npmax parameter ',npmax
                   STOP
                 ENDIF
                ENDDO

              ELSE IF( number_of_parts( l, i, j, k ) == 0 .AND. &
                                                          mass < 0.0 ) THEN
                IF( ABS( mass ) .GT. backgroundmass ) THEN
                  mass = - backgroundmass
                ENDIF
!                minmass = ABS(mass) / 10.0
                maxmass = maxParticleMass( l )
                IF (maxmass == 0.0 ) THEN
!                   maxmass = mass / 10.0
                    maxmass = MAXVAL( maxParticleMass )
                ENDIF
                DO WHILE( mass < 0.0 )
!                 mass = mass + minmass
                 mass = mass + maxmass
!                 backgroundmass = backgroundmass - minmass
                 backgroundmass = backgroundmass - maxmass
!                 IF ( rmcnter2 .GT. rmcnter1 ) THEN
!                      idx = removed( rmcnter2 )
!                      rmcnter2 = rmcnter2 - 1
!                      totalremoved = totalremoved - 1
                 IF ( totalremoved .GT. 0 ) THEN
                      idx = removed( totalremoved )
                      totalremoved = totalremoved - 1
                 ELSE 
                   idx = np + 1
                   np = np + 1
                 ENDIF
                 IF ( np <= npmax ) THEN
                  pp(idx,1) = ( i - 1 ) * delv( 1 ) + delv( 1 ) * ran1( 0 )
                  pp(idx,2) = ( j - 1 ) * delv( 2 ) + delv( 2 ) * ran1( 0 )
                  pp(idx,3) = ( k - 1 ) * delv( 3 ) + delv( 3 ) * ran1( 0 )
                  IF (mass < 0.0 ) THEN
                    pp(idx,4) = -maxmass
                  !pp(idx,4) = -minmass
                  ELSE
                    pp(idx,4) = mass - maxmass
                  ENDIF
                  pp(idx,5) = DBLE(idx)
                  pp(idx,7) = timeCellAllSpec( i, j, k, pp, modelname ) 
                  !pp(idx,7) = timeCellAllSpec( i, j, k, pp, np, ip, delv ) 
                  ip(idx,1) = l
                  ip(idx,2) = 1

                 ELSE
                   PRINT*, '2 max num particles exceeded. &
                           &increase npmax parameter ',npmax
                   STOP
                 ENDIF

                ENDDO

              ELSE IF( number_of_parts(l, i, j, k ) > 0 .AND. mass > 0.0 ) THEN

                m = 0
                DO WHILE( mass .GT. 0.0 .AND.              &
                     m .LT.  number_of_parts( l, i,j,k ) )

                  m = m + 1 
                  IF( pp( part_numbers( l, i,j,k)%arr( m ), 4) .LE. 0.0 ) THEN
                      mass = mass + pp( part_numbers(l,i,j,k)%arr( m ), 4 )
                      totalremoved = totalremoved + 1
!                      IF ( rmcnter1 .GT. 0 ) THEN
!                         rmcnter1 = rmcnter1 - 1
!                         removed( rmcnter1 ) = partNumber(m) 
!                      ELSE
!                         rmcnter2 = rmcnter2 + 1
!                         removed( rmcnter2 ) = partNumber(m) 
!                      ENDIF
                         removed( totalremoved ) = part_numbers(l,i,j,k)%arr(m) 
                  ENDIF
                ENDDO

!                minmass = ABS( pp( partNumber( 1 ), 4 ) )
!                minmass = mass / 10.0
                maxmass = maxParticleMass( l )
!                IF( minmass .GT. mass ) THEN
!                   minmass = mass
!                ENDIF
                DO WHILE( mass .GT. 0.0 )
!                 mass = mass - minmass
                 mass = mass - maxmass
!                 IF ( rmcnter2 .GT. rmcnter1 ) THEN
!                   idx = removed( rmcnter2 )
!                   rmcnter2 = rmcnter2 - 1 
!                   totalremoved = totalremoved - 1
                 IF ( totalremoved .GT. 0 ) THEN
                   idx = removed( totalremoved )
                   totalremoved = totalremoved - 1
                   recycle = .TRUE.
                 ELSE 
                   idx = np + 1
                   np = np + 1
                   recycle = .FALSE.
                 ENDIF

                 IF ( np <= npmax / 2 .OR. recycle ) THEN

                    pp(idx,1) = ( i - 1 ) * delv( 1 ) + delv( 1 ) * ran1( 0 )
                    pp(idx,2) = ( j - 1 ) * delv( 2 ) + delv( 2 ) * ran1( 0 )
                    pp(idx,3) = ( k - 1 ) * delv( 3 ) + delv( 3 ) * ran1( 0 )
                    IF ( mass > 0.0 ) THEN
!                      pp(idx,4) = minmass
                      pp(idx,4) = maxmass
                    ELSE
                      pp(idx,4) = mass + maxmass
                    ENDIF
                    pp(idx,5) = DBLE(idx)
                    pp(idx,7) = pp( part_numbers(l,i,j,k)%arr(1), 7 ) 

                    ip(idx,1) = l
                    ip(idx,2) = 1
                    ip(idx,3) = ip( part_numbers(l,i,j,k)%arr(1), 3 ) 
                    ip(idx,4) = ip( part_numbers(l,i,j,k)%arr(1), 4 )
                    ip(idx,5) = ip( part_numbers(l,i,j,k)%arr(1), 5 )
                    ip(idx,6) = ip( part_numbers(l,i,j,k)%arr(1), 6 )
                    ip(idx,7) = ip( part_numbers(l,i,j,k)%arr(1), 7 )
                    ip(idx,8) = ip( part_numbers(l,i,j,k)%arr(1), 8 )
                    ip(idx,9) = ip( part_numbers(l,i,j,k)%arr(1), 9 )
                    ip(idx,10) =ip( part_numbers(l,i,j,k)%arr(1), 10 )

                 ELSE
!                   PRINT*, '1 max num particles exceeded. &
!                           &increase npmax parameter ',npmax
!                   STOP
                  IF( number_of_parts(l,i,j,k) > totalremoved ) THEN
                     DO  m = 1, number_of_parts(l,i,j,k)
                       found = .false.
                       DO p = 1, totalremoved
                         IF( part_numbers(l,i,j,k)%arr( m ) .EQ.  & 
                                removed( p ) )                      THEN
                                found = .true.
                                EXIT
                         ENDIF
                       ENDDO
                       IF( .NOT. found ) THEN
                         pp( part_numbers(l,i,j,k)%arr(m), 4 ) =             &
                             pp( part_numbers(l,i,j,k)%arr(m), 4 ) +         &
                                      ( mass + maxmass ) /                  &
                                  ( number_of_parts(l,i,j,k) - totalremoved ) 
                       ENDIF
                     ENDDO    
                  ELSE
                   PRINT*, 'WARNING: total removed greater than or equal to &
                   &number of parts' 
                   STOP
                  ENDIF
                  mass = -1.0 !exit do while
                  np = np - 1
                 ENDIF
                 ENDDO ! DO WHILE

              ELSE IF( number_of_parts( l, i, j, k ) > 0 .AND. &
                                                          mass < 0.0 ) THEN
                particlemass = 0.0
                DO n = 1, number_of_parts(l, i,j,k )
                   particlemass = particlemass + pp( part_numbers(l,i,j,k)%arr( n ), 4 )
                ENDDO

                m = 0
                DO WHILE( mass .LT. 0.0 .AND.              &
                     m .LT.  number_of_parts( l, i,j,k ) .AND.     &
                    particlemass + backgroundmass .GT. 0.0 )

                  m = m + 1 
                  IF( pp( part_numbers(l,i,j,k)%arr( m ), 4) .GT. 0.0 ) THEN
                      mass = mass + pp( part_numbers(l,i,j,k)%arr( m ), 4 )
                      particlemass = particlemass - pp( part_numbers(l,i,j,k)%arr(m), 4 )
                      totalremoved = totalremoved + 1
!                      IF ( rmcnter1 .GT. 0 ) THEN
!                         rmcnter1 = rmcnter1 - 1
!                         removed( rmcnter1 ) = partNumber(m) 
!                      ELSE
!                         rmcnter2 = rmcnter2 + 1
!                         removed( rmcnter2 ) = partNumber(m) 
!                      ENDIF
                         removed( totalremoved ) = part_numbers(l,i,j,k)%arr(m) 
                  ENDIF
                ENDDO

!                minmass = ABS( pp( partNumber( 1 ), 4 ) )
!                minmass = ABS( mass ) / 10.0
                maxmass = maxParticleMass( l )
!                IF( minmass .GT. ABS( mass ) ) THEN
!                   minmass = ABS( mass )
!                ENDIF
                DO WHILE( mass < 0.0 .AND.    &
                      particlemass + backgroundmass .GT. 0.0  )
!                 mass = mass + minmass
                 mass = mass + maxmass
!                 particlemass = particlemass - minmass
                 particlemass = particlemass - maxmass
!                 IF ( rmcnter2 .GT. rmcnter1 ) THEN
!                   idx = removed( rmcnter2 )
!                   rmcnter2 = rmcnter2 - 1 
!                   totalremoved = totalremoved - 1
                 IF ( totalremoved .GT. 0 ) THEN
                   idx = removed( totalremoved )
                   totalremoved = totalremoved - 1
                   recycle = .TRUE.
                 ELSE 
                   idx = np + 1
                   np = np + 1
                   recycle = .FALSE.
                 ENDIF
                 IF ( np <= npmax / 2 .OR. recycle) THEN
                  pp(idx,1) = ( i - 1 ) * delv( 1 ) + delv( 1 ) * ran1( 0 )
                  pp(idx,2) = ( j - 1 ) * delv( 2 ) + delv( 2 ) * ran1( 0 )
                  pp(idx,3) = ( k - 1 ) * delv( 3 ) + delv( 3 ) * ran1( 0 )
                  IF ( mass < 0.0 ) THEN
!                    pp(idx,4) = -minmass
                    pp(idx,4) = -maxmass
                  ELSE
                    pp(idx,4) = mass - maxmass
                  ENDIF
                  pp(idx,5) = DBLE(idx)
                  pp(idx,7) = timeCellAllSpec( i, j, k, pp, modelname ) 
                  !pp(idx,7) = timeCellAllSpec( i, j, k, pp, np, ip, delv ) 
                  ip(idx,1) = l
                  ip(idx,2) = 1

                 ELSE
!                   PRINT*, '2 max num particles exceeded. &
!                           &increase npmax parameter ',npmax
!                   STOP

                  IF( number_of_parts(l,i,j,k) > totalremoved ) THEN
                     DO  m = 1, number_of_parts(l,i,j,k)
                       found = .false.
                       DO p = 1, totalremoved
                         IF( part_numbers(l,i,j,k)%arr( m ) .EQ.  & 
                                removed( p ) )                      THEN
                                found = .true.
                                EXIT
                         ENDIF
                       ENDDO
                       IF( .NOT. found ) THEN
                        IF ( abs( mass - maxmass ) < backgroundmass ) THEN 
                          pp( part_numbers(l,i,j,k)%arr(m), 4 ) =             &
                               pp( part_numbers(l,i,j,k)%arr(m), 4 ) +        &
                                      ( mass - maxmass ) /                   &
                                  ( number_of_parts(l,i,j,k) - totalremoved ) 
                        ELSE
                          pp( part_numbers(l,i,j,k)%arr(m), 4 ) =             &
                               pp( part_numbers(l,i,j,k)%arr(m), 4 ) -        &
                                      backgroundmass /                       &
                                  ( number_of_parts(l,i,j,k) - totalremoved ) 
                        ENDIF

                       ENDIF
                     ENDDO    
                  ELSE
                   PRINT*, 'WARNING: total removed greater than or equal to &
                   &number of parts' 
                   STOP
                  ENDIF
                  mass = 1.0 !exit do while
                  np = np - 1
                 ENDIF
                ENDDO
            ENDIF ! IF( number_of_parts( l, i, j, k ) == 0 .AND. mass > 0.0)
           END DO
        ENDIF !IF ( SUM( number_of_parts(1:ns,i,j,k) ) .GT. 0 ) THEN 
!        IF( modelname == "MacQ" ) THEN
!         NitrifierConc( i, j, k ) = NitrifierConc( i, j, k ) + rates( 12 )
!         DenitrifierConc( i, j, k ) = DenitrifierConc( i, j, k ) +  rates( 13 )
!        ENDIF
         IF ( modelname == 'Chen1992' ) THEN
             conc(4, i, j, k) = Biomass1Conc( i, j, k)
                conc(5, i, j, k) = Biomass2Conc( i, j, k)
         ELSE IF( modelname == 'MacQ1990' ) THEN
                conc(3, i, j, k) = BiomassConc( i, j, k)
         ELSE IF( modelname == 'MacQ1990unsat' ) THEN
                conc(3, i, j, k) = UBiomassConc( i, j, k)
         ELSE IF( modelname == 'MacQ' ) THEN
         ELSE IF( modelname == 'VG' ) THEN
         ELSE IF( modelname == 'noreact' ) THEN
         ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
         ENDIF

      END DO
    END DO
  END DO

!       CALL removeParticles( pp, ip, ipwell, irp, iprp, lastprint, np, &
!                             npmax, removed, totalremoved )
    CALL markRemovedParticles( pp, removed, totalremoved )


! Calculate the two biomasses
!       DO i = 3, nx - 2
!         DO j = 3, ny - 2
!           DO k = 3, nz - 2
!            IF ( SUM( number_of_parts(1:ns,i,j,k) ) .GT. 0 ) THEN 
!              substratename = 'NH4'
!              rate = ReactRateAtCell( conc, i, j, k, substratename )
!
!              NitrifierConc( i, j, k ) = NitrifierConc( i, j, k ) +        &
!
!                DeltaX1( rate, NitrifierConc( i, j, k ),                   &
!                         npars%Biomass1YieldCoeff_M_biomass_M_substrate,   &
!                         npars%Biomass1SpecDecayConst_1_day ) * dt_day
!
!              substratename = 'CH2O'
!              rate = ReactRateAtCell( conc, i, j, k, substratename )
!              DenitrifierConc( i, j, k ) = DenitrifierConc( i, j, k ) +    &
!                 DeltaX2( rate, DenitrifierConc( i, j, k ),                &
!                         npars%Biomass2YieldCoeff_M_biomass_M_substrate,   &
!                         npars%Biomass2SpecDecayConst_1_day ) * dt_day
!        ENDIF! IF ( SUM( number_of_parts(1:ns,i,j,k) ) .GT. 0 ) THEN 
!      END DO
!    END DO
!  END DO

     END SUBROUTINE addRemoveParticles

     REAL*8 FUNCTION minMassCellSpec( s, x, y, z, pp, np, delv, ip )

       INTEGER :: s, x, y, z, i, ix, iy, iz
       REAL*8, DIMENSION(:), ALLOCATABLE :: mass
       REAL*8, DIMENSION(:,:) :: pp
       INTEGER*4,  DIMENSION(:,:) :: ip
       REAL*8, DIMENSION(3) :: delv
       INTEGER*4 :: np
       LOGICAL :: found

       DO i = 1, np 
         ix = IDINT( pp( i, 1 ) / delv( 1 ) ) + 1
         iy = IDINT( pp( i, 2 ) / delv( 2 ) ) + 1
         iz = IDINT( pp( i, 3 ) / delv( 3 ) ) + 1
         IF ( ip(i,1) .EQ. s .AND.              &
              ix .EQ. x .AND. iy .EQ. y .AND. iz .EQ. z ) THEN
          found = .TRUE.
          minMassCellSpec = pp( i, 4 )
          EXIT 
         ENDIF  
       END DO
     
       IF ( .NOT. found ) THEN
!         PRINT*, 'NO particles are found in cell( ', x, ', ', y, ', ', z, ')' 
!         PRINT*, 'for species ', s 
!         PRINT*, 'when creating new particles, set time of particle to zero!' 

       ENDIF
     END FUNCTION minMassCellSpec

!     REAL*8 FUNCTION minMassCellAllSpec( ns, x, y, z, pp )
!
!       INTEGER :: ns, s, x, y, z, i, cnt
!       REAL*8, DIMENSION(:), ALLOCATABLE :: mass
!       REAL*8, DIMENSION(:,:) :: pp
!     
!       cnt = 0
!
!!       DO s = 1, TotalNumberOfSpecies 
!       DO s = 1, ns 
!         cnt = cnt + number_of_parts( s, x, y, z )
!       ENDDO
!
!       ALLOCATE( mass( cnt) )
!
!       cnt = 0
!!       DO s = 1, TotalNumberOfSpecies 
!       DO s = 1, ns 
!         DO i = 1, number_of_parts( s, x, y, z ) 
!           cnt = cnt + 1
!           mass( cnt ) = pp( part_numbers( s, x, y, z )%arr( i ), 4 ) 
!         ENDDO
!       END DO
!
!       minMassCellAllSpec = MINVAL( mass )
!
!       DEALLOCATE( mass )
!
!     END FUNCTION minMassCellAllSpec

!     REAL*8 FUNCTION minMassAllCellSpec( s, nx, ny, nz, pp )
!
!       INTEGER :: s, nx, ny, nz, i, j, k, l, cnt
!       REAL*8, DIMENSION(:), ALLOCATABLE :: mass
!       REAL*8, DIMENSION(:,:) :: pp
!     
!       cnt = 0
!
!!       DO s = 1, TotalNumberOfSpecies 
!       DO i = 1, nx 
!         DO j = 1, ny 
!           DO k = 1, ny 
!            cnt = cnt + number_of_parts( s, i, j, k )
!           ENDDO
!         ENDDO
!       ENDDO
!
!       ALLOCATE( mass( cnt) )
!
!       cnt = 0
!!       DO s = 1, TotalNumberOfSpecies 
!       DO i = 1, nx 
!         DO j = 1, ny 
!           DO k = 1, ny 
!            DO l = 1, number_of_parts( s, i, j, k ) 
!              cnt = cnt + 1
!              mass( cnt ) = pp( part_numbers( s, i, j, k )%arr( l ), 4 ) 
!            ENDDO
!           ENDDO
!         ENDDO
!       ENDDO
!
!       minMassAllCellSpec = MINVAL( mass )
!
!       DEALLOCATE( mass )
!
!     END FUNCTION minMassAllCellSpec

!     REAL*8 FUNCTION minMassAllCellAllSpec( pp, np )
!
!       INTEGER*4 :: i, np, cnt
!       REAL*8, DIMENSION(:), ALLOCATABLE :: mass
!       REAL*8, DIMENSION(:,:) :: pp
!     
!       cnt = 0
!
!!       DO s = 1, TotalNumberOfSpecies 
!       DO i = 1, np 
!         IF (pp(i, 4) > 0.0 ) THEN
!           cnt = cnt + 1
!         ENDIF
!       ENDDO
!
!       ALLOCATE( mass( cnt) )
!
!       cnt = 0
!!       DO s = 1, TotalNumberOfSpecies 
!       DO i = 1, np 
!         IF (pp(i, 4) > 0.0 ) THEN
!           cnt = cnt + 1
!           mass( cnt ) = pp( i, 4 ) 
!         ENDIF 
!       END DO
!
!       minMassAllCellAllSpec = MINVAL( mass )
!
!       DEALLOCATE( mass )
!
!     END FUNCTION minMassAllCellAllSpec

     SUBROUTINE maxMassAllCellSpec( pp, np, ip )
       INTEGER :: I
       INTEGER*4 ::  np
       REAL*8, DIMENSION(:,:) :: pp
       INTEGER*4,  DIMENSION(:,:) :: ip

       maxParticleMass = 0.0

       DO I = 1, np
         IF ( ABS( pp( I, 4 ) ) > maxParticleMass( ip(I, 1) ) ) THEN
            maxParticleMass( ip(I, 1 ) ) = ABS( pp(I,4) )
         ENDIF
       END DO

       RETURN
     END SUBROUTINE

     REAL*8 FUNCTION timeCellAllSpec( x, y, z, pp, modelname )
!     REAL*8 FUNCTION timeCellAllSpec( x, y, z, pp, np, ip, delv )

!       INTEGER ::  x, y, z, i, ix, iy, iz
       INTEGER ::  x, y, z, i, ns_moving
       REAL*8, DIMENSION(:,:) :: pp
       LOGICAL found
!       INTEGER*4 ::  np
!       INTEGER*4,  DIMENSION(:,:) :: ip
!       REAL*8, DIMENSION(3) :: delv
       CHARACTER (LEN=20) :: modelname
     
      found = .FALSE. 
      timeCellAllSpec = 0.0

!      DO i = 1, np
!         ix = IDINT( pp( i, 1 ) / delv( 1 ) ) + 1
!         iy = IDINT( pp( i, 2 ) / delv( 2 ) ) + 1
!         iz = IDINT( pp( i, 3 ) / delv( 3 ) ) + 1
!         IF ( ip(i,2) .EQ. 1 .AND.              &
!              ix .EQ. x .AND. iy .EQ. y .AND. iz .EQ. z ) THEN
!          found = .TRUE.
!          timeCellAllSpec = pp( i, 7 )
!          EXIT 
!         ENDIF  
!      END DO

       IF ( modelname == 'Chen1992' ) THEN
               ns_moving = Chen1992TotalNumberOfSpecies - 2
       ELSE IF( modelname == 'MacQ1990' ) THEN
               ns_moving = MacQ1990TotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ1990unsat' ) THEN
               ns_moving = MacQ1990unsatTotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ' ) THEN
               ns_moving = TotalNumberOfSpecies 
       ELSE IF( modelname == 'VG' ) THEN
               ns_moving = VGTotalNumberOfSpecies 
       ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
       ENDIF
      DO i = 1, ns_moving
        IF ( number_of_parts(i, x, y, z ) > 0 ) THEN
           found = .TRUE.
           timeCellAllSpec = pp( part_numbers(i, x, y, z)%arr(1), 7 )
           EXIT
        ENDIF
      ENDDO

      IF ( .NOT. found ) THEN
 !        PRINT*, 'NO particles are found in cell( ', x, ', ', y, ', ', z, ')' 
 !        PRINT*, 'when creating new particles, set time of particle to zero!' 
      ENDIF

     END FUNCTION timeCellAllSpec

     SUBROUTINE partNumbersAtCellSpec(  partNumber,                        &
                                            s, x, y, z, pp, ip, np,         &
                                            num_parts, delv )
       INTEGER :: s, x, y, z, i, cnt, num_parts, ix, iy, iz
       REAL*8, DIMENSION(:,:) :: pp
       INTEGER*4,  DIMENSION(:,:) :: ip
       INTEGER*4 :: np
       REAL*8, DIMENSION(3) :: delv
       INTEGER, DIMENSION(:) :: partNumber

!       cnt = 0
!       DO i = 1, np
!         ix = IDINT( pp( i, 1 ) / delv( 1 ) ) + 1
!         iy = IDINT( pp( i, 2 ) / delv( 2 ) ) + 1
!         iz = IDINT( pp( i, 3 ) / delv( 3 ) ) + 1
!         IF ( ip(i,2) .EQ. 1 .AND. ip(i,1) .EQ. s .AND.             &
!              ix .EQ. x .AND. iy .EQ. y .AND. iz .EQ. z ) THEN
!           cnt = cnt + 1
!           partNumber( cnt ) = i
!         ENDIF  
!       ENDDO
      DO I = 1, num_parts
        partNumber( I ) = part_numbers(s,x,y,z)%arr(I)
      ENDDO
     END SUBROUTINE partNumbersAtCellSpec

!     SUBROUTINE removeParticles( pp, ip, ipwell, irp, iprp, lastprint, np, &
!                                                         npmax, removedparts, &
!                                                         rmpartnum  )
!       REAL*8, DIMENSION(:,:) :: pp, lastprint
!       INTEGER, DIMENSION(:,:) :: ipwell
!       INTEGER*4,  DIMENSION(:,:) :: ip, iprp
!       INTEGER*4,  DIMENSION(:) ::  irp
!       INTEGER,  DIMENSION(:) :: removedparts 
!
!       INTEGER :: i, j, k, rmpartnum, cnt
!       INTEGER*4 ::  np, npmax
!
!       !copy pp
!
!       tmpp = pp
!       tmplastprint = lastprint
!       tmpipwell = ipwell
!       tmpip = ip
!       tmpirp = irp
!       tmpiprp = iprp
!
!       cnt = 0
!       IF ( rmpartnum > 0 ) THEN
!         DO i = 1, np
!           IF ( .NOT. ANY( removedparts(1:rmpartnum) .EQ. i ) ) THEN
!             cnt = cnt + 1
!             DO j =1, 10
!               pp( cnt, j ) = tmpp( i, j )
!               ip( cnt, j ) = tmpip( i, j )
!             ENDDO
!             pp( cnt, 5 ) = cnt
!             DO j =1, 20
!               ipwell( cnt, j ) = tmpipwell( i, j )
!             ENDDO
!             DO j =1, 3 
!               lastprint( cnt, j ) = tmplastprint( i, j )
!             ENDDO
!             DO j =1, 2 
!               iprp( cnt, j ) = tmpiprp( i, j )
!             ENDDO
!             irp( cnt ) = tmpirp( i )
!           ENDIF
!         ENDDO
!
!!         DO i = cnt + 1, np
!!             DO j =1, 10
!!!               pp( i, j ) = 0.0
!!               ip( i, j ) = 0
!!!             ENDDO
!!             DO j =1, 20
!!               ipwell( i, j ) = 0
!!             ENDDO
!!             DO j =1, 3 
!!               lastprint( i, j ) = 0.0
!!             ENDDO
!!             DO j =1, 2 
!!               iprp( i, j ) = 0 
!!             ENDDO
!!             irp( i ) = 0
!!         ENDDO
!         np = cnt
!       ENDIF
!
!     END SUBROUTINE removeParticles

     SUBROUTINE markRemovedParticles( pp, removedparts, totalremoved  )
       REAL*8, DIMENSION(:,:) :: pp
       INTEGER,  DIMENSION(:) :: removedparts 

       INTEGER :: i, totalremoved

       DO i = 1, totalremoved
         pp( removedparts( i ), 1 ) = -999
         pp( removedparts( i ), 2 ) = -999
         pp( removedparts( i ), 3 ) = -999
       ENDDO

     END SUBROUTINE markRemovedParticles

     SUBROUTINE  addGasDispersion( n, pp, ip, delv, tloc, fbl, porst, satutn, &
                                   dxx, dyy, dzz, modelname )
      INTEGER*4 :: n, tp1, tp2, tp3
      INTEGER*4 :: i, j, k
      INTEGER*4, DIMENSION(:,:) :: ip
      INTEGER*4, DIMENSION(:) :: tloc
      REAL*8, DIMENSION(:) :: fbl, delv
      REAL*8, DIMENSION(:,:) :: pp
      REAL*8, DIMENSION(:,:,:) :: porst, satutn
      REAL*8 :: dxx, dyy, dzz, D, H
      REAL*8, DIMENSION(4,4,4) :: pbl, sbl
      CHARACTER (LEN=20) :: modelname

      IF ( modelname == 'MacQ' ) THEN
       IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' .OR.        &
           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' .OR.       &
           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2'  ) THEN

        DO  i = 1, 3
          DO  j = 1, 3
            DO  k = 1, 3

              tp1 = tloc(1) -2 + i
              tp2 = tloc(2) -2 + j
              tp3 = tloc(3) -2 + k
             
              pbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * porst( tp1, tp2, tp3 ) +                   &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * porst( tp1, tp2 + 1, tp3 )      +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           porst( tp1, tp2, tp3 + 1 )  +                &
                           fbl(2) * fbl(3) * porst( tp1, tp2 + 1, tp3 + 1 )
              sbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * satutn( tp1, tp2, tp3 ) +                  &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * satutn( tp1, tp2 + 1, tp3 )     +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           satutn( tp1, tp2, tp3 + 1 )  +               &
                           fbl(2) * fbl(3) * satutn( tp1, tp2 + 1, tp3 + 1 )
            ENDDO
          ENDDO
        ENDDO

        IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' ) THEN
          H = npars%N2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2' ) THEN
          H = npars%O2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' ) THEN
          H = npars%CO2HenryConst
          D = npars%SubstrateFreeAirDiffusionCoeff_m2_day

        ENDIF

        dxx = dxx + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 3, 2, 2 ) - pbl(1,2,2) ) / ( 2 * delv( 1 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 3, 2, 2 ) - sbl( 1,2,2) ) / ( 2 * delv( 1 ) ) ) * &
                 H * D

        dyy = dyy + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 2, 3, 2 ) - pbl(2,1,2) ) / ( 2 * delv( 2 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 2, 3, 2 ) - sbl( 2,1,2) ) / ( 2 * delv( 2 ) ) ) * &
                 H * D
        dzz = dzz + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 2, 2, 3 ) - pbl(2,2,1) ) / ( 2 * delv( 3 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 2, 2, 3 ) - sbl( 2,2,1) ) / ( 2 * delv( 3 ) ) ) * &
                 H * D

       ENDIF
      ENDIF

     END SUBROUTINE addGasDispersion

     SUBROUTINE  addGasDispersionX( n, pp, ip, delv, tloc, fbl, porst, satutn, &
                                   dxx, modelname )
      INTEGER*4 :: n, tp1, tp2, tp3
      INTEGER*4 :: i, j, k
      INTEGER*4, DIMENSION(:,:) :: ip
      INTEGER*4, DIMENSION(:) :: tloc
      REAL*8, DIMENSION(:) :: fbl, delv
      REAL*8, DIMENSION(:,:) :: pp
      REAL*8, DIMENSION(:,:,:) :: porst, satutn
      REAL*8 :: dxx, D, H
      REAL*8, DIMENSION(4,4,4) :: pbl, sbl

      CHARACTER (LEN=20) :: modelname

      IF ( modelname == 'MacQ' ) THEN
!      IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' .OR.        &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' .OR.       &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2'  ) THEN
       IF( ip(n, 1 ) .EQ. 10 .OR. ip(n, 1) .EQ. 5 .OR. ip(n,1) .EQ. 4 ) THEN  

        DO  i = 1, 3
          DO  j = 1, 3
            DO  k = 1, 3

              tp1 = tloc(1) -2 + i
              tp2 = tloc(2) -2 + j
              tp3 = tloc(3) -2 + k
             
              pbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * porst( tp1, tp2, tp3 ) +                   &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * porst( tp1, tp2 + 1, tp3 )      +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           porst( tp1, tp2, tp3 + 1 )  +                &
                           fbl(2) * fbl(3) * porst( tp1, tp2 + 1, tp3 + 1 )
              sbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * satutn( tp1, tp2, tp3 ) +                  &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * satutn( tp1, tp2 + 1, tp3 )     +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           satutn( tp1, tp2, tp3 + 1 )  +               &
                           fbl(2) * fbl(3) * satutn( tp1, tp2 + 1, tp3 + 1 )
            ENDDO
          ENDDO
        ENDDO

!        IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' ) THEN
        IF( ip(n, 1 ) .EQ. 10 ) THEN  
          H = npars%N2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

!        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2' ) THEN
        ELSE IF( ip(n, 1 ) .EQ. 4 ) THEN  
          H = npars%O2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

!        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' ) THEN
        ELSE IF( ip(n, 1 ) .EQ. 5 ) THEN  
          H = npars%CO2HenryConst
          D = npars%SubstrateFreeAirDiffusionCoeff_m2_day

        ENDIF

        dxx = dxx + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 3, 2, 2 ) - pbl(1,2,2) ) / ( 2 * delv( 1 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 3, 2, 2 ) - sbl( 1,2,2) ) / ( 2 * delv( 1 ) ) ) * &
                 H * D
        ENDIF

      ELSE IF( modelname == 'MacQ1990unsat' ) THEN

       IF ( ip(n, 1) .EQ. 1. .OR. ip(n, 1) .EQ. 2 ) THEN
        DO  i = 1, 3
          DO  j = 1, 3
            DO  k = 1, 3

              tp1 = tloc(1) -2 + i
              tp2 = tloc(2) -2 + j
              tp3 = tloc(3) -2 + k
             
              pbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * porst( tp1, tp2, tp3 ) +                   &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * porst( tp1, tp2 + 1, tp3 )      +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           porst( tp1, tp2, tp3 + 1 )  +                &
                           fbl(2) * fbl(3) * porst( tp1, tp2 + 1, tp3 + 1 )
              sbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * satutn( tp1, tp2, tp3 ) +                  &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * satutn( tp1, tp2 + 1, tp3 )     +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           satutn( tp1, tp2, tp3 + 1 )  +               &
                           fbl(2) * fbl(3) * satutn( tp1, tp2 + 1, tp3 + 1 )
            ENDDO
          ENDDO
        ENDDO

        IF ( ip(n, 1) .EQ. 1 ) THEN
          H = umpars%SH
          D = umpars%SD0
        ELSE IF ( ip(n, 1) .EQ. 2 ) THEN
          H = umpars%EH
          D = umpars%ED0
        ENDIF

        dxx = dxx + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 3, 2, 2 ) - pbl(1,2,2) ) / ( 2 * delv( 1 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 3, 2, 2 ) - sbl( 1,2,2) ) / ( 2 * delv( 1 ) ) ) * &
                 H * D
       ENDIF

      ENDIF

     END SUBROUTINE addGasDispersionX

     SUBROUTINE  addGasDispersionY( n, pp, ip, delv, tloc, fbl, porst, satutn, &
                                   dyy, modelname )
      INTEGER*4 :: n, tp1, tp2, tp3
      INTEGER*4 :: i, j, k
      INTEGER*4, DIMENSION(:,:) :: ip
      INTEGER*4, DIMENSION(:) :: tloc
      REAL*8, DIMENSION(:) :: fbl, delv
      REAL*8, DIMENSION(:,:) :: pp
      REAL*8, DIMENSION(:,:,:) :: porst, satutn
      REAL*8 :: dyy, D, H
      REAL*8, DIMENSION(4,4,4) :: pbl, sbl
      CHARACTER (LEN=20) :: modelname

      IF( modelname == 'MacQ' ) THEN
!      IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' .OR.        &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' .OR.       &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2'  ) THEN
       IF( ip(n, 1 ) .EQ. 10 .OR. ip(n, 1) .EQ. 5 .OR. ip(n,1) .EQ. 4 ) THEN  

        DO  i = 1, 3
          DO  j = 1, 3
            DO  k = 1, 3

              tp1 = tloc(1) -2 + i
              tp2 = tloc(2) -2 + j
              tp3 = tloc(3) -2 + k
             
              pbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * porst( tp1, tp2, tp3 ) +                   &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * porst( tp1, tp2 + 1, tp3 )      +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           porst( tp1, tp2, tp3 + 1 )  +                &
                           fbl(2) * fbl(3) * porst( tp1, tp2 + 1, tp3 + 1 )
              sbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * satutn( tp1, tp2, tp3 ) +                  &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * satutn( tp1, tp2 + 1, tp3 )     +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           satutn( tp1, tp2, tp3 + 1 )  +               &
                           fbl(2) * fbl(3) * satutn( tp1, tp2 + 1, tp3 + 1 )
            ENDDO
          ENDDO
        ENDDO

!        IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' ) THEN
        IF ( ip(n, 1) .EQ. 10 ) THEN
          H = npars%N2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

!        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2' ) THEN
        ELSE IF ( ip(n, 1) .EQ. 4 ) THEN
          H = npars%O2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

!        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' ) THEN
        ELSE IF ( ip(n, 1) .EQ. 5 ) THEN
          H = npars%CO2HenryConst
          D = npars%SubstrateFreeAirDiffusionCoeff_m2_day

        ENDIF

        dyy = dyy + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 2, 3, 2 ) - pbl(2,1,2) ) / ( 2 * delv( 2 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 2, 3, 2 ) - sbl( 2,1,2) ) / ( 2 * delv( 2 ) ) ) * &
                 H * D
        ENDIF

       ELSE IF( modelname == 'MacQ1990unsat' ) THEN

       IF ( ip(n, 1) .EQ. 1. .OR. ip(n, 1) .EQ. 2 ) THEN

        DO  i = 1, 3
          DO  j = 1, 3
            DO  k = 1, 3

              tp1 = tloc(1) -2 + i
              tp2 = tloc(2) -2 + j
              tp3 = tloc(3) -2 + k
             
              pbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * porst( tp1, tp2, tp3 ) +                   &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * porst( tp1, tp2 + 1, tp3 )      +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           porst( tp1, tp2, tp3 + 1 )  +                &
                           fbl(2) * fbl(3) * porst( tp1, tp2 + 1, tp3 + 1 )
              sbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * satutn( tp1, tp2, tp3 ) +                  &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * satutn( tp1, tp2 + 1, tp3 )     +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           satutn( tp1, tp2, tp3 + 1 )  +               &
                           fbl(2) * fbl(3) * satutn( tp1, tp2 + 1, tp3 + 1 )
            ENDDO
          ENDDO
        ENDDO

        IF ( ip(n, 1) .EQ. 1 ) THEN
          H = umpars%SH
          D = umpars%SD0
        ELSE IF ( ip(n, 1) .EQ. 2 ) THEN
          H = umpars%EH
          D = umpars%ED0
        ENDIF

        dyy = dyy + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 2, 3, 2 ) - pbl(2,1,2) ) / ( 2 * delv( 2 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 2, 3, 2 ) - sbl( 2,1,2) ) / ( 2 * delv( 2 ) ) ) * &
                 H * D
      ENDIF

      ENDIF
     END SUBROUTINE addGasDispersionY

     SUBROUTINE  addGasDispersionZ( n, pp, ip, delv, tloc, fbl, porst, satutn, &
                                   dzz, modelname )
      INTEGER*4 :: n, tp1, tp2, tp3
      INTEGER*4 :: i, j, k
      INTEGER*4, DIMENSION(:,:) :: ip
      INTEGER*4, DIMENSION(:) :: tloc
      REAL*8, DIMENSION(:) :: fbl, delv
      REAL*8, DIMENSION(:,:) :: pp
      REAL*8, DIMENSION(:,:,:) :: porst, satutn
      REAL*8 :: dzz, D, H
      REAL*8, DIMENSION(4,4,4) :: pbl, sbl
      CHARACTER (LEN=20) :: modelname

      IF ( modelname == 'MacQ' ) THEN
!      IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' .OR.        &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' .OR.       &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2'  ) THEN
       IF( ip(n, 1 ) .EQ. 10 .OR. ip(n, 1) .EQ. 5 .OR. ip(n,1) .EQ. 4 ) THEN  

        DO  i = 1, 3
          DO  j = 1, 3
            DO  k = 1, 3

              tp1 = tloc(1) -2 + i
              tp2 = tloc(2) -2 + j
              tp3 = tloc(3) -2 + k
             
              pbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * porst( tp1, tp2, tp3 ) +                   &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * porst( tp1, tp2 + 1, tp3 )      +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           porst( tp1, tp2, tp3 + 1 )  +                &
                           fbl(2) * fbl(3) * porst( tp1, tp2 + 1, tp3 + 1 )
              sbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * satutn( tp1, tp2, tp3 ) +                  &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * satutn( tp1, tp2 + 1, tp3 )     +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           satutn( tp1, tp2, tp3 + 1 )  +               &
                           fbl(2) * fbl(3) * satutn( tp1, tp2 + 1, tp3 + 1 )
            ENDDO
          ENDDO
        ENDDO

!        IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' ) THEN
        IF ( ip(n, 1) .EQ. 10 ) THEN
          H = npars%N2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

        ELSE IF ( ip(n, 1) .EQ. 4 ) THEN
          H = npars%O2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

!        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' ) THEN
        ELSE IF ( ip(n, 1) .EQ. 5 ) THEN
          H = npars%CO2HenryConst
          D = npars%SubstrateFreeAirDiffusionCoeff_m2_day

        ENDIF

        dzz = dzz + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 2, 2, 3 ) - pbl(2,2,1) ) / ( 2 * delv( 3 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 2, 2, 3 ) - sbl( 2,2,1) ) / ( 2 * delv( 3 ) ) ) * &
                 H * D
       ENDIF
      ELSE IF ( modelname == 'MacQ1990unsat' ) THEN

      IF ( ip(n, 1) .EQ. 1. .OR. ip(n, 1) .EQ. 2 ) THEN

        DO  i = 1, 3
          DO  j = 1, 3
            DO  k = 1, 3

              tp1 = tloc(1) -2 + i
              tp2 = tloc(2) -2 + j
              tp3 = tloc(3) -2 + k
             
              pbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * porst( tp1, tp2, tp3 ) +                   &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * porst( tp1, tp2 + 1, tp3 )      +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           porst( tp1, tp2, tp3 + 1 )  +                &
                           fbl(2) * fbl(3) * porst( tp1, tp2 + 1, tp3 + 1 )
              sbl(i,j,k) = ( 1.0D0 - fbl(2) ) * ( 1.0D0 - fbl(3) )      &
                           * satutn( tp1, tp2, tp3 ) +                  &
                           fbl(2) * ( 1.0D0 - fbl(3) )                  &
                           * satutn( tp1, tp2 + 1, tp3 )     +          &
                           ( 1.0D0 - fbl(2) ) * fbl(3) *                &
                           satutn( tp1, tp2, tp3 + 1 )  +               &
                           fbl(2) * fbl(3) * satutn( tp1, tp2 + 1, tp3 + 1 )
            ENDDO
          ENDDO
        ENDDO

        IF ( ip(n, 1) .EQ. 1 ) THEN
          H = umpars%SH
          D = umpars%SD0
        ELSE IF ( ip(n, 1) .EQ. 2 ) THEN
          H = umpars%EH
          D = umpars%ED0
        ENDIF

        dzz = dzz + ( ( 1.0D0 -                                     &
                 satutn( tloc(1), tloc(2), tloc(3) ) )**3.33333 /       &
                 satutn( tloc(1), tloc(2), tloc(3) ) / ( 3 *            &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.66666 ) *      &
                 ( pbl( 2, 2, 3 ) - pbl(2,2,1) ) / ( 2 * delv( 3 ) ) -  &
                 porst( tloc(1), tloc(2), tloc(3) ) ** 0.33333 * (      &
                 ( 10.0D0 *                                             &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3 ) ) ) ** 2.33333 ) / &
                 ( 3 * satutn( tloc(1), tloc(2), tloc(3) ) ) +          &
                 ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 /   &
                  satutn( tloc(1), tloc(2), tloc(3) ) ** 2 )  *         &
                 ( sbl( 2, 2, 3 ) - sbl( 2,2,1) ) / ( 2 * delv( 3 ) ) ) * &
                 H * D
       ENDIF

      ENDIF
     END SUBROUTINE addGasDispersionZ

     SUBROUTINE  addGasDispersion2( n, pp, ip, tloc, z, al, at, vn, betad,  &
                                   porst, satutn, cx, cy, cz, modelname  )
      INTEGER*4 :: n 
      INTEGER*4, DIMENSION(:) :: tloc
      INTEGER*4, DIMENSION(:,:) :: ip
      REAL*8, DIMENSION(:) :: z
      REAL*8, DIMENSION(:,:) :: pp
      REAL*8, DIMENSION(:,:,:) :: porst, satutn
      REAL*8 :: al, at, cx, cy, cz, vn, betad, D, H
      CHARACTER (LEN=20) :: modelname

      IF ( modelname == 'MacQ' ) THEN
      IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' .OR.        &
           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' .OR.       &
           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2'  ) THEN

        IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' ) THEN
          H = npars%N2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2' ) THEN
          H = npars%O2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' ) THEN
          H = npars%CO2HenryConst
          D = npars%SubstrateFreeAirDiffusionCoeff_m2_day
        ENDIF

        cx = z(1)*DSQRT( (2D0 * al ) / vn +                                &
              2D0 * ( porst( tloc(1), tloc(2), tloc(3)) ** 0.333333 *      &
                   ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 2.33333 ) &
                   * H * D )
        cy = z(2)*DSQRT( (2D0 * at * vn ) +                                &
              2D0 * ( porst( tloc(1), tloc(2), tloc(3)) ** 0.333333 *      &
                   ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 2.33333 ) &
                   * H * D ) / betad
        cz = z(3)*DSQRT( (2D0 * at * vn ) +                                &
              2D0 * ( porst( tloc(1), tloc(2), tloc(3)) ** 0.333333 *      &
                   ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 2.33333 ) &
                   * H * D ) / ( betad * vn )

       ENDIF
      ENDIF
     END SUBROUTINE addGasDispersion2

     SUBROUTINE addGasDiffusion( n, ip, tloc, porst, satutn, moldiff, &
                                  newmldif, modelname  )

      INTEGER*4 :: n 
      INTEGER*4, DIMENSION(:) :: tloc
      INTEGER*4, DIMENSION(:,:) :: ip
      REAL*8, DIMENSION(:,:,:) :: porst, satutn
      REAL*8 :: D, H, moldiff, newmldif
      CHARACTER (LEN=20) :: modelname

      IF ( modelname == 'MacQ' ) THEN
!      IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' .OR.        &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' .OR.       &
!           TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2'  ) THEN
       IF( ip(n, 1 ) .EQ. 10 .OR. ip(n, 1) .EQ. 5 .OR. ip(n,1) .EQ. 4 ) THEN  

!        IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'N2' ) THEN
        IF ( ip(n, 1) .EQ. 10 ) THEN
          H = npars%N2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

!        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'O2' ) THEN
        ELSE IF ( ip(n, 1) .EQ. 4 ) THEN
          H = npars%O2HenryConst
          D = npars%EAcceptorFreeAirDiffusionCoeff_m2_day

!        ELSE IF ( TRIM( npars%SpeciesNames( ip(n, 1) ) ) .EQ. 'CO2' ) THEN
        ELSE IF ( ip(n, 1) .EQ. 5 ) THEN
          H = npars%CO2HenryConst
          D = npars%SubstrateFreeAirDiffusionCoeff_m2_day
        ENDIF

        newmldif = moldiff +                                                   &
              ( porst( tloc(1), tloc(2), tloc(3)) ** 0.333333 *               &
                   ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 ) &
                   / satutn( tloc(1), tloc(2), tloc(3) )                      &
                   * H * D
      ELSE
        newmldif = moldiff 
      ENDIF

      ELSE IF ( modelname == 'MacQ1990unsat' ) THEN

      IF ( ip(n, 1) .EQ. 1. .OR. ip(n, 1) .EQ. 2 ) THEN

        IF ( ip(n, 1) .EQ. 1 ) THEN
          H = umpars%SH
          D = umpars%SD0
        ELSE IF ( ip(n, 1) .EQ. 2 ) THEN
          H = umpars%EH
          D = umpars%ED0
        ENDIF

        newmldif = moldiff +                                                   &
              ( porst( tloc(1), tloc(2), tloc(3)) ** 0.333333 *               &
                   ( 1D0 - satutn( tloc(1), tloc(2), tloc(3) ) ) ** 3.33333 ) &
                   / satutn( tloc(1), tloc(2), tloc(3) )                      &
                   * H * D
      ELSE
        newmldif = moldiff 
      ENDIF

      ENDIF
     END SUBROUTINE addGasDiffusion

     !
     ! to prevent negative concentration. Remove some negative mass particles
     ! when the sum of conc and background concentration is less than zero
     !
     SUBROUTINE checkConc( pp, cellv, concmgl, porst, satutn, &
                          nx, ny, nz, ns, modelname )

       INTEGER*4 :: nx, ny, nz, ns, i,j,k,l, p, ns_moving
       INTEGER*4 :: cnt
!       INTEGER*4, ALLOCATABLE, TARGET :: temp(:)
       INTEGER*4, ALLOCATABLE :: temp(:)
       REAL*8, DIMENSION(:,:) :: pp
       REAL*8 :: mass, cellv, cellvLiter 
       REAL*4, DIMENSION(:,:,:,:) :: concmgl
       REAL*8, DIMENSION(:,:,:) :: porst, satutn
       CHARACTER*20 :: modelname

       REAL*4 :: bconc

       IF ( modelname == 'Chen1992' ) THEN
               ns_moving = Chen1992TotalNumberOfSpecies - 2
       ELSE IF( modelname == 'MacQ1990' ) THEN
               ns_moving = MacQ1990TotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ1990unsat' ) THEN
               ns_moving = MacQ1990unsatTotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ' ) THEN
               ns_moving = TotalNumberOfSpecies 
       ELSE IF( modelname == 'VG' ) THEN
               ns_moving = VGTotalNumberOfSpecies 
       ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
       ENDIF
       cellvLiter = cellv * 1000 !M^3 to liter
       DO l = 1, ns_moving
         DO i = 1, nx
          DO j = 1, ny
             DO k = 1, nz
               cnt = 0
               IF ( modelname == 'Chen1992' ) THEN
                 bconc =cpars%backgroundConc(l)
               ELSE IF( modelname == 'MacQ1990' ) THEN
                 bconc =mpars%backgroundConc(l)
               ELSE IF( modelname == 'MacQ1990unsat' ) THEN
                 bconc =umpars%backgroundConc(l)
               ELSE IF( modelname == 'MacQ' ) THEN
                 bconc =npars%backgroundConc(l)
               ELSE IF( modelname == 'VG' ) THEN
                 bconc =vgpars%backgroundConc(l)
               ELSE 
                  WRITE(*,*) 'ERROR: non-known model name: ', modelname
                  stop
               ENDIF
               DO WHILE ( concmgl( l, i, j, k ) + bconc     &
                     .LT. 0.0 .AND. cnt .LE. number_of_parts( l, i, j, k ) )
                cnt = cnt + 1
                mass = 0.0
                DO p = cnt + 1, number_of_parts( l, i, j, k )
                  mass = mass + pp( part_numbers(l, i,j,k)%arr( p ), 4 )
                END DO
                concmgl( l, i, j, k ) = mass/ ( cellvLiter * satutn( i, j, k ) &
                                           * porst( i, j, k ) )
                totalremoved = totalremoved + 1
                removed( totalremoved ) = part_numbers(l, i,j,k)%arr( cnt ) 
                pp( part_numbers(l, i,j,k)%arr( cnt ), 1 ) = -999
                pp( part_numbers(l, i,j,k)%arr( cnt ), 2 ) = -999
                pp( part_numbers(l, i,j,k)%arr( cnt ), 3 ) = -999
               END DO  ! DO WHILE
               !
               ! delete particles from part_numbers
               !
               IF ( cnt .GT. 0 ) THEN

                 ALLOCATE( temp( number_of_parts( l, i, j, k ) - cnt ) )

                 temp = part_numbers(l, i,j,k)%arr(                      &
                      cnt + 1 : number_of_parts( l, i, j, k  ) )

                 DEALLOCATE( part_numbers( l, i, j, k)%arr )              

                 ALLOCATE( part_numbers( l, i, j, k )%arr(               &
                                       number_of_parts(l, i, j, k ) - cnt ) )

!                 part_numbers( l, i, j, k )%arr => temp
                 part_numbers( l, i, j, k )%arr = temp

                 DEALLOCATE( temp )
                 number_of_parts( l, i, j, k ) =                         &
                            number_of_parts(l, i, j,k ) - cnt
               ENDIF
             END DO
           END DO
         END DO
       END DO

    END SUBROUTINE checkConc

     SUBROUTINE ConcToMgPerLiter( conc, nx, ny, nz, ns )
       REAL*4, DIMENSION(:,:,:,:) :: conc
       INTEGER*4 :: nx, ny, nz, ns, i,j,k,l
        conc = conc / 1000.0 !mg/L^3 to mg/l
    END SUBROUTINE ConcToMgPerLiter

    SUBROUTINE updateConc( conc, pp, ipp, np, delc, domax, sat, porosity, &
                                                      Rtard, tnext )
       REAL*4, DIMENSION(:,:,:,:) :: conc
       REAL*8, DIMENSION(:,:,:,:) :: Rtard
       REAL*8, DIMENSION(:,:,:)   :: sat, porosity
       REAL*8, DIMENSION(:,:)     :: pp
       REAL*8, DIMENSION(:)       ::    delc
       REAL*8                     :: tnext, Rstar, cellv
       INTEGER*4, DIMENSION(:,:)  :: ipp
       INTEGER*4, DIMENSION(:)    :: domax
       INTEGER*4                  :: np, i, x, y, z

       cellv = delc(1) * delc(2) * delc(3) * 1000.0 !M^3 to L 
       conc = 0.0
       DO i = 1, np
        x = IDINT( pp( i, 1 ) / delc( 1 ) ) + 1
        y = IDINT( pp( i, 2 ) / delc( 2 ) ) + 1
        z = IDINT( pp( i, 3 ) / delc( 3 ) ) + 1

        IF( ( x .GE. 1 .AND. x .LE. domax(1) ) .AND.    &
            ( y .GE. 1 .AND. y .LE. domax(2) ) .AND.    &
            ( z .GE. 1 .AND. z .LE. domax(3) ) .AND.    &
            ( pp( i, 7 ) .LE. tnext ) ) THEN

          rstar = 1.d0 + ( Rtard( ipp( i, 1 ), x, y, z )- 1.d0 )/ sat( x, y, z )
          conc( ipp(i,1), x, y, z) = conc( ipp( i, 1 ), x, y, z ) +         &
                SNGL( dble( ipp( i, 2 ) ) * pp( i, 4 )            /         &
                ( cellv * sat( x, y, z) * porosity( x, y, z) * Rstar ) )
        ENDIF 
       ENDDO
    END SUBROUTINE updateConc

     SUBROUTINE ConcAddBackground( conc, nx, ny, nz, ns, modelname )
       REAL*4, DIMENSION(:,:,:,:) :: conc
       INTEGER*4 :: nx, ny, nz, ns, i,j,k,l, ns_moving
       CHARACTER*20 :: modelname

       REAL*4 :: bconc

       IF ( modelname == 'Chen1992' ) THEN
               ns_moving = Chen1992TotalNumberOfSpecies - 2
       ELSE IF( modelname == 'MacQ1990' ) THEN
               ns_moving = MacQ1990TotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ1990unsat' ) THEN
               ns_moving = MacQ1990unsatTotalNumberOfSpecies - 1
       ELSE IF( modelname == 'MacQ' ) THEN
               ns_moving = TotalNumberOfSpecies 
       ELSE IF( modelname == 'VG' ) THEN
               ns_moving = VGTotalNumberOfSpecies 
       ELSE IF( modelname == 'noreact' ) THEN
               ns_moving = 0
       ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
       ENDIF

       DO i = 1, nx
        DO j = 1, ny
           DO k = 1, nz
             DO l = 1, ns_moving
               IF ( modelname == 'Chen1992' ) THEN
                  bconc =cpars%backgroundConc(l)
               ELSE IF( modelname == 'MacQ1990' ) THEN
                  bconc =mpars%backgroundConc(l)
               ELSE IF( modelname == 'MacQ1990unsat' ) THEN
                  bconc =umpars%backgroundConc(l)
               ELSE IF( modelname == 'MacQ' ) THEN
                  bconc =npars%backgroundConc(l)
               ELSE IF( modelname == 'VG' ) THEN
                  bconc =vgpars%backgroundConc(l)
               ELSE IF( modelname == 'noreact' ) THEN
               ELSE 
                  WRITE(*,*) 'ERROR: non-known model name: ', modelname
                  stop
               ENDIF
                conc(l, i, j, k) = conc(l, i,j,k) + bconc
                IF ( conc( l, i, j, k ) .LT. 0.0 ) conc( l, i, j, k ) = 0.0 
             ENDDO
             IF ( modelname == 'Chen1992' ) THEN
                conc(4, i, j, k) = Biomass1Conc( i, j, k)
                conc(5, i, j, k) = Biomass2Conc( i, j, k)
             ELSE IF( modelname == 'MacQ1990' ) THEN
                conc(3, i, j, k) = BiomassConc( i, j, k)
             ELSE IF( modelname == 'MacQ1990unsat' ) THEN
                conc(3, i, j, k) = UBiomassConc( i, j, k)
             ELSE IF( modelname == 'MacQ' ) THEN
             ELSE IF( modelname == 'VG' ) THEN
             ELSE IF( modelname == 'noreact' ) THEN
             ELSE 
               WRITE(*,*) 'ERROR: non-known model name: ', modelname
               stop
             ENDIF
          ENDDO
         ENDDO
       ENDDO



     END SUBROUTINE ConcAddBackground

     
END MODULE Particles
