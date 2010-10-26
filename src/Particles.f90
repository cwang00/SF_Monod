

MODULE Particles

   USE NTransport

   IMPLICIT NONE

   TYPE :: Ptr_to_array
   
      INTEGER*4, DIMENSION(:), POINTER :: arr 

   END TYPE Ptr_to_array

   TYPE ( Ptr_to_array ), DIMENSION(:,:,:,:), ALLOCATABLE :: part_numbers 
   INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: number_of_parts
!   REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmpp, tmplastprint
!   INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmpipwell
!   INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: tmpip, tmpiprp
!   INTEGER*4, DIMENSION(:), ALLOCATABLE :: tmpirp
!   INTEGER, DIMENSION(:), ALLOCATABLE :: partNumber
   INTEGER*4, DIMENSION(:), ALLOCATABLE :: removed
!   INTEGER totalremoved, rmcnter1, rmcnter2
   REAL*8, DIMENSION(TotalNumberOfSpecies) :: maxParticleMass

   CONTAINS

     SUBROUTINE countParticles( pp,ip, totalparticles, nx, ny, nz, ns, &
                                 delv )

!       TYPE( Ptr_to_array ), DIMENSION(:,:,:,:) :: part_num
!       INTEGER, DIMENSION(:,:,:,:) :: number_of_parts
       REAL*8, DIMENSION(:,:) :: pp
       INTEGER*4, DIMENSION(:,:) :: ip
       INTEGER*4 :: nx, ny, nz, ns, totalparticles, i,j,k,l, p, x, y, z
       REAL*8, DIMENSION(3) :: delv 
       INTEGER*4, DIMENSION(ns, nx, ny, nz) :: pcount

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
         IF ( ip(p,2) .EQ. 1 .AND. x .GT. 0 .AND. y .GT. 0 .AND. z .GT. 0) THEN
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
         IF ( ip(p,2) .EQ. 1 .AND. x .GT. 0 .AND. y .GT. 0 .AND. z .GT. 0) THEN
           pcount( ip(p, 1), x, y, z ) = pcount( ip(p, 1), x, y, z ) + 1
           part_numbers( ip(p, 1), x, y, z )%arr( pcount( ip(p, 1), x, y,z ))&
                                                                         = p
         ENDIF
       END DO
     END SUBROUTINE countParticles

     SUBROUTINE allocateParticles_Memory(  nx, ny, nz, ns, npmax )
        INTEGER*4 :: nx, ny, nz, ns, npmax
        ALLOCATE( number_of_parts( ns, nx, ny, nz ) )
!       ALLOCATE( tmpp( npmax, 10 ) )
!       ALLOCATE( tmplastprint( npmax, 3 ) )
!       ALLOCATE( tmpipwell( npmax, 20 ) )
!       ALLOCATE( tmpip( npmax, 10 ) )
!       ALLOCATE( tmpirp( npmax ) )
!       ALLOCATE( tmpiprp( npmax, 2 ) )
       ALLOCATE( removed( npmax ) )
!       ALLOCATE( partNumber(npmax) )
       ALLOCATE( part_numbers( ns, nx, ny, nz ) )
!       totalremoved = 0
!       rmcnter1 = 0
!       rmcnter2 = 0
     END SUBROUTINE allocateParticles_Memory

     SUBROUTINE deallocateParticles_Memory(  nx, ny, nz, ns )

        INTEGER*4 :: nx, ny, nz, ns, l, i, j, k

!         IF( ALLOCATED( number_of_parts ) ) THEN
             DEALLOCATE( number_of_parts )
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
                              porosity, sat )
       REAL*8, DIMENSION(:,:) :: pp, lastprint
       REAL*8, DIMENSION(:,:,:) :: porosity, sat
       INTEGER, DIMENSION(:,:) :: ipwell
       INTEGER*4,  DIMENSION(:,:) :: ip, iprp
       INTEGER*4,  DIMENSION(:) ::  irp
       REAL*8, DIMENSION(3) :: delv
       INTEGER :: ns, nx, ny, nz, i, j, k, l, m, totalremoved,n, idx
       REAL*4, DIMENSION(:,:,:,:) :: conc
       REAL, DIMENSION(13) :: rates
       REAL*4 :: rate, minmass,maxmass, mass, backgroundmass, particlemass
       REAL*8 :: dt_day
       INTEGER*4 ::  npmax, np
       CHARACTER*10 substratename

    totalremoved = 0
    DO i = 3, nx - 2
      DO j = 3, ny - 2
        DO k = 3, nz - 2
        
          IF ( SUM( number_of_parts(1:ns,i,j,k) ) .GT. 0 ) THEN 
            CALL VODEReactionRatesAtCell( conc, i, j, k, REAL(dt_day), rates )
            DO l = 1, ns
!              CALL partNumbersAtCellSpec( partNumber,                        &
!                                         l, i, j, k, pp, ip, np,             &
!                                         number_of_parts(l,i,j,k), delv )

!              IF( conc( l,i,j,k) .LT. 0 .AND. conc(l,i,j,k) +           &
!                 npars%backgroundConc(l) .LT. 0 ) THEN
!                 WRITE(*,*) l, i, j, k
!              ENDIF 

!              rate = ReactRateAtCell( conc, i, j, k, npars%SpeciesNames( l ) )
              IF ( TRIM(npars%SpeciesNames( l) ) .eq. 'NH4' ) THEN
                 rate = rates( 1 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'NO3' ) THEN
                 rate = rates( 2 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'CH2O' ) THEN
                 rate = rates( 3 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'O2' ) THEN
                 rate = rates( 4 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'CO2' ) THEN
                 rate = rates( 5 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'HCO3' ) THEN
                 rate = rates( 6 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'H' ) THEN
                 rate = rates( 7 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'CO3' ) THEN
                 rate = rates( 8 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'Ca2' ) THEN
                 rate = rates( 9 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'N2' ) THEN
                 rate = rates( 10 )
              ELSE IF( TRIM(npars%SpeciesNames( l )) .eq. 'OH' ) THEN
                 rate = rates( 11 )
              ELSE
                 PRINT*, 'Unknow species naem ', npars%SpeciesNames(l)

              ENDIF
              mass = rate * delv(1) * delv(2) * delv(3)                     &
                        * sat( i,j , k)  * porosity( i,j,k)
              backgroundmass = npars%backgroundConc( l )                    &
                           * delv(1) * delv(2) * delv(3)                    &
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
                  pp(idx,7) = timeCellAllSpec( i, j, k, pp, np, ip, delv ) 
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
                  pp(idx,7) = timeCellAllSpec( i, j, k, pp, np, ip, delv ) 
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
                 ELSE 
                   idx = np + 1
                   np = np + 1
                 ENDIF

                 IF ( np <= npmax ) THEN

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
                 ELSE 
                   idx = np + 1
                   np = np + 1
                 ENDIF
                 IF ( np <= npmax ) THEN
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
                  pp(idx,7) = timeCellAllSpec( i, j, k, pp, np, ip, delv ) 
                  ip(idx,1) = l
                  ip(idx,2) = 1

                 ELSE
                   PRINT*, '2 max num particles exceeded. &
                           &increase npmax parameter ',npmax
                   STOP
                 ENDIF
                ENDDO
            ENDIF ! IF( number_of_parts( l, i, j, k ) == 0 .AND. mass > 0.0)
           END DO
        ENDIF !IF ( SUM( number_of_parts(1:ns,i,j,k) ) .GT. 0 ) THEN 
        NitrifierConc( i, j, k ) = NitrifierConc( i, j, k ) + rates( 12 )
        DenitrifierConc( i, j, k ) = DenitrifierConc( i, j, k ) +  rates( 13 )

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
         PRINT*, 'NO particles are found in cell( ', x, ', ', y, ', ', z, ')' 
         PRINT*, 'for species ', s 
         PRINT*, 'when creating new particles, set time of particle to zero!' 

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

     REAL*8 FUNCTION timeCellAllSpec( x, y, z, pp, np, ip, delv )

       INTEGER ::  x, y, z, i, ix, iy, iz
       INTEGER*4 ::  np
       REAL*8, DIMENSION(:,:) :: pp
       INTEGER*4,  DIMENSION(:,:) :: ip
       REAL*8, DIMENSION(3) :: delv
       LOGICAL found
     
      found = .FALSE. 
      timeCellAllSpec = 0.0
       DO i = 1, np 
         ix = IDINT( pp( i, 1 ) / delv( 1 ) ) + 1
         iy = IDINT( pp( i, 2 ) / delv( 2 ) ) + 1
         iz = IDINT( pp( i, 3 ) / delv( 3 ) ) + 1
         IF ( ip(i,2) .EQ. 1 .AND.              &
              ix .EQ. x .AND. iy .EQ. y .AND. iz .EQ. z ) THEN
          found = .TRUE.
          timeCellAllSpec = pp( i, 7 )
          EXIT 
         ENDIF  
       END DO

       IF ( .NOT. found ) THEN
         PRINT*, 'NO particles are found in cell( ', x, ', ', y, ', ', z, ')' 
         PRINT*, 'when creating new particles, set time of particle to zero!' 

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
! 
FUNCTION ran1(idum)
INTEGER*4 idum,IA,IM,IQ,IR,NTAB,NDIV
REAL*8 ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,  &
   NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM

      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END FUNCTION

     
END MODULE Particles
