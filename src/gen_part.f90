!
! Generate particles for given concentration at indicator locations
!
! Zhengtao Cui
! Oct. 9, 2011
!
subroutine gen_part( icat, jcat, kcat, ic_cat_num, n_ic_cats, ic_cats,         &
                     ic_conc, ic_part_dens, particle, ip,                      &
                     delv, vel, num_of_parts,                                  &
                     timestep_num, constitute_num, npmax, current_conc,        &
                     time_begin, time_end, porosity, sat, bnd_cnd ) 
   USE Particles

INTEGER*4 :: ii,jj,kk, icat, jcat, kcat, n_ic_cats,  &
           num_of_parts, iii, timestep_num, ir,     &
          constitute_num, npmax, part
INTEGER*4, DIMENSION(:,:,:) :: ic_cat_num
INTEGER*4, DIMENSION(:,:) :: ip
INTEGER*4, DIMENSION(:) :: ic_cats, ic_part_dens
INTEGER :: increased_num_of_parts

real*8  :: time_begin, time_end, rantemp, part_mass, mass, current_mass
REAL*8, DIMENSION(:,:) :: particle, ic_conc
REAL*8, DIMENSION(:,:,:) :: sat, porosity
REAL*8, DIMENSION(:, :,:,:) :: vel
REAL*4, DIMENSION(:,:,:,:) ::  current_conc
REAL*8, DIMENSION(:) :: delv
CHARACTER (LEN=20) :: bnd_cnd

ir = 21
!  PRINT*,' 3rd ic: loop thru ic'
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
          DO iii = 1, n_ic_cats
           IF( ic_part_dens( iii ) .GT. 0 ) THEN
            IF( timestep_num == 1 ) THEN
               IF ( bnd_cnd == 'const_conc' ) THEN
                 mass = ( ic_conc( timestep_num, iii ) -            &
                   current_conc(constitute_num, ii, jj,kk ) ) &
                     * delv( 1 ) * delv( 2 ) * delv( 3 )      &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               ELSE           
                  mass = ic_conc( 1, iii ) * delv( 1 ) * delv( 2 ) * delv( 3 ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               ENDIF           
               IF ( mass .LT. 0.0 ) mass = 0.0
            ELSE
              IF ( bnd_cnd == 'const_flux' ) THEN
               mass = ic_conc( timestep_num, iii ) * delv( 2 ) * delv( 3 ) &
                       * ABS( vel(1, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               mass = mass +  ic_conc( timestep_num, iii ) * delv( 1 ) * delv( 3 ) &
                       * ABS( vel(2, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 *  porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               mass = mass +  ic_conc( timestep_num, iii ) * delv( 1 ) * delv( 2 ) &
                       * ABS( vel(3, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               ELSE IF ( bnd_cnd == 'const_conc' ) THEN
                 mass = ( ic_conc( timestep_num, iii ) -            &
                   current_conc(constitute_num, ii, jj,kk ) ) &
                     * delv( 1 ) * delv( 2 ) * delv( 3 )      &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 

                 IF ( mass .LT. 0.0 ) mass = 0.0

               ELSE
                WRITE(*,*) 'ERROR: non-known boundary condition: ', bnd_cnd
                stop
               ENDIF
              ENDIF
            part_mass =  mass / ic_part_dens( iii )

            IF (ic_cat_num(ii,jj,kk) == ic_cats(iii)) THEN
!               current_mass = current_conc(constitute_num, ii, jj,kk ) &
!                    * delv( 1 ) * delv( 2 ) * delv( 3 ) &
!                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 

                current_mass = 0.0
               IF ( mass .GT. current_mass .AND. mass .GE. 0   &
                    .AND. current_mass .GE.0 ) THEN

               increased_num_of_parts = ( mass - current_mass ) / part_mass

! assign particles to active locations
! check to see if mass of ic > 0
!             IF ( ABS( ic_conc(timestep_num,iii) ) >                        &
!                    ABS( current_conc( constitute_num, ii, jj, kk) ) ) THEN
                !
                ! assumes ic_mass_or_conc is mg/l, delv in meters, mass is in
                ! mg
                !
!                mass =  mass - current_conc(constitute_num, ii, jj , kk )   &
!                        * delv( 1 ) * delv( 2 ) * delv( 3 ) * 1000 *        &
!                          porosity( ii, jj,kk ) * sat(ii,jj,kk) 
!              DO  WHILE(  mass > 0.0 )
               DO I = 1, increased_num_of_parts
                IF ( totalremoved .GT. 0 ) THEN
                   part = removed( totalremoved); 
                   totalremoved = totalremoved - 1;
                ELSE
                   part = num_of_parts;
                   num_of_parts = num_of_parts + 1
                   IF (num_of_parts > npmax) THEN
                     PRINT*, 'max num particles exceeded.  increase npmax parameter ',npmax
                     WRITE(666,*) 'Max Particles Exceeded!  Increase NPMAX:',npmax
                     STOP
                   END IF ! max parts ?
                ENDIF
	        particle(part,1) = delv(1)*ran1(ir) + DBLE(ii-1)*delv(1)
                particle(part,2) = delv(2)*ran1(ir) + DBLE(jj-1)*delv(2)
                particle(part,3) = delv(3)*ran1(ir) + DBLE(kk-1)*delv(3)
                particle(part,4) = part_mass
                particle(part,5) = part
                particle(part,7) = time_begin +                &
                                           ran1(ir)*(time_end - time_begin)
                particle(part,9) = ic_cat_num(ii,jj,kk)
                ip(part,1) = constitute_num
                ip(part,2) = 1
               
                mass = mass - part_mass
                !
                ! update concentration
                !
                current_conc( constitute_num, ii, jj, kk ) =              &
                             current_conc( constitute_num, ii, jj, kk ) + &
                             part_mass / ( delv(1) * delv(2) * delv(3)    &
                             * 1000 * porosity( ii, jj, kk ) *            &
                             sat( ii, jj, kk ) )
                                  
              END DO  !DO WHILE DO I
              
      ELSE IF ( mass .LT. current_mass .AND. mass .GE. 0   &
               .AND. current_mass .GE. 0 ) THEN
        !remove particles

        DO WHILE ( current_mass .GT. mass )
          part = part_numbers( constitute_num, ii, jj, kk)%arr(    &
             number_of_parts( constitute_num, ii, jj, kk ) )
	  particle(part,1) = -999.0
          particle(part,2) = -999.0
          particle(part,3) = -999.0

          totalremoved = totalremoved + 1
          removed( totalremoved ) = part

          number_of_parts( constitute_num, ii, jj, kk ) =           &
               number_of_parts( constitute_num, ii, jj, kk ) - 1  

          current_mass = current_mass - particle( part, 4 )
                !
                ! update concentration
                !
                current_conc( constitute_num, ii, jj, kk ) =              &
                             current_conc( constitute_num, ii, jj, kk ) - &
                       particle(part, 4)/( delv(1) * delv(2) * delv(3)    &
                             * 1000 * porosity( ii, jj, kk ) *            &
                             sat( ii, jj, kk ) )
        END DO

      ELSE IF ( mass .GT. current_mass .AND. mass .LT. 0   &
               .AND. current_mass .LT. 0 ) THEN
        !remove particles

        DO WHILE ( current_mass .LT. mass )
          part = part_numbers( constitute_num, ii, jj, kk)%arr(    &
             number_of_parts( constitute_num, ii, jj, kk ) )
	  particle(part,1) = -999.0
          particle(part,2) = -999.0
          particle(part,3) = -999.0

          totalremoved = totalremoved + 1
          removed( totalremoved ) = part

          number_of_parts( constitute_num, ii, jj, kk ) =         &
               number_of_parts( constitute_num, ii, jj, kk ) - 1  

          current_mass = current_mass - particle( part, 4 )
                !
                ! update concentration
                !
                current_conc( constitute_num, ii, jj, kk ) =              &
                             current_conc( constitute_num, ii, jj, kk ) - &
                       particle(part, 4)/( delv(1) * delv(2) * delv(3)    &
                             * 1000 * porosity( ii, jj, kk ) *            &
                             sat( ii, jj, kk ) )
        END DO

      ELSE IF ( mass .LT. current_mass .AND. mass .LT. 0   &
               .AND. current_mass .LT. 0 ) THEN

               increased_num_of_parts = ( mass - current_mass ) / part_mass

               DO I = 1, increased_num_of_parts
                IF ( totalremoved .GT. 0 ) THEN
                   part = removed( totalremoved); 
                   totalremoved = totalremoved - 1;
                ELSE
                   part = num_of_parts;
                   num_of_parts = num_of_parts + 1
                   IF (num_of_parts > npmax) THEN
                     PRINT*, 'max num particles exceeded.  increase npmax parameter ',npmax
                     WRITE(666,*) 'Max Particles Exceeded!  Increase NPMAX:',npmax
                     STOP
                   END IF ! max parts ?
                ENDIF
	        particle(part,1) = delv(1)*ran1(ir) + DBLE(ii-1)*delv(1)
                particle(part,2) = delv(2)*ran1(ir) + DBLE(jj-1)*delv(2)
                particle(part,3) = delv(3)*ran1(ir) + DBLE(kk-1)*delv(3)
                particle(part,4) = part_mass
                particle(part,5) = part 
                particle(part,7) = time_begin +                &
                                           ran1(ir)*(time_end - time_begin)
                particle(part,9) = ic_cat_num(ii,jj,kk)
                ip(part,1) = constitute_num
                ip(part,2) = 1
               
                mass = mass - part_mass
                !
                ! update concentration
                !
                current_conc( constitute_num, ii, jj, kk ) =              &
                             current_conc( constitute_num, ii, jj, kk ) + &
                             part_mass / ( delv(1) * delv(2) * delv(3)    &
                             * 1000 * porosity( ii, jj, kk ) *            &
                             sat( ii, jj, kk ) )
                                  
              END DO  !DO WHILE DO I

      ELSE IF ( mass .LT. current_mass .AND. mass .LT. 0   &
               .AND. current_mass .GE. 0 ) THEN

        DO I = 1, number_of_parts( constitute_num, ii, jj, kk )
          part = part_numbers( constitute_num, ii, jj, kk)%arr( I )
	  particle(part,1) = -999.0
          particle(part,2) = -999.0
          particle(part,3) = -999.0
          totalremoved = totalremoved + 1
          removed( totalremoved ) = part
        END DO
        number_of_parts( constitute_num, ii, jj, kk ) = 0
        current_conc( constitute_num, ii, jj, kk ) = 0

        DO I = 1, ic_part_dens( iii )
                IF ( totalremoved .GT. 0 ) THEN
                   part = removed( totalremoved); 
                   totalremoved = totalremoved - 1;
                ELSE
                   part = num_of_parts;
                   num_of_parts = num_of_parts + 1
                   IF (num_of_parts > npmax) THEN
                     PRINT*, 'max num particles exceeded.  increase npmax parameter ',npmax
                     WRITE(666,*) 'Max Particles Exceeded!  Increase NPMAX:',npmax
                     STOP
                   END IF ! max parts ?
                ENDIF
	        particle(part,1) = delv(1)*ran1(ir) + DBLE(ii-1)*delv(1)
                particle(part,2) = delv(2)*ran1(ir) + DBLE(jj-1)*delv(2)
                particle(part,3) = delv(3)*ran1(ir) + DBLE(kk-1)*delv(3)
                particle(part,4) = part_mass
                particle(part,5) = part
                particle(part,7) = time_begin +                &
                                           ran1(ir)*(time_end - time_begin)
                particle(part,9) = ic_cat_num(ii,jj,kk)
                ip(part,1) = constitute_num
                ip(part,2) = 1
               
                !
                ! update concentration
                !
                current_conc( constitute_num, ii, jj, kk ) =              &
                             current_conc( constitute_num, ii, jj, kk ) + &
                             part_mass / ( delv(1) * delv(2) * delv(3)    &
                             * 1000 * porosity( ii, jj, kk ) *            &
                             sat( ii, jj, kk ) )

        END DO

      ELSE IF ( mass .GT. current_mass .AND. mass .GE. 0   &
               .AND. current_mass .LT. 0 ) THEN

        DO I = 1, number_of_parts( constitute_num, ii, jj, kk )
          part = part_numbers( constitute_num, ii, jj, kk)%arr( I )
	  particle(part,1) = -999.0
          particle(part,2) = -999.0
          particle(part,3) = -999.0
          totalremoved = totalremoved + 1
          removed( totalremoved ) = part
        END DO
        number_of_parts( constitute_num, ii, jj, kk ) = 0
        current_conc( constitute_num, ii, jj, kk ) = 0

        DO I = 1, ic_part_dens( iii )
                IF ( totalremoved .GT. 0 ) THEN
                   part = removed( totalremoved); 
                   totalremoved = totalremoved - 1;
                ELSE
                   part = num_of_parts;
                   num_of_parts = num_of_parts + 1
                   IF (num_of_parts > npmax) THEN
                     PRINT*, 'max num particles exceeded.  increase npmax parameter ',npmax
                     WRITE(666,*) 'Max Particles Exceeded!  Increase NPMAX:',npmax
                     STOP
                   END IF ! max parts ?
                ENDIF
	        particle(part,1) = delv(1)*ran1(ir) + DBLE(ii-1)*delv(1)
                particle(part,2) = delv(2)*ran1(ir) + DBLE(jj-1)*delv(2)
                particle(part,3) = delv(3)*ran1(ir) + DBLE(kk-1)*delv(3)
                particle(part,4) = part_mass
                particle(part,5) = part
                particle(part,7) = time_begin +                &
                                           ran1(ir)*(time_end - time_begin)
                particle(part,9) = ic_cat_num(ii,jj,kk)
                ip(part,1) = constitute_num
                ip(part,2) = 1
               
                !
                ! update concentration
                !
                current_conc( constitute_num, ii, jj, kk ) =              &
                             current_conc( constitute_num, ii, jj, kk ) + &
                             part_mass / ( delv(1) * delv(2) * delv(3)    &
                             * 1000 * porosity( ii, jj, kk ) *            &
                             sat( ii, jj, kk ) )

        END DO

             END IF  ! ic_conc > current_conc ? 
            END IF  ! at a ic node
           END IF 
          END DO !iii
      END DO  !ii
    END DO  !jj
  END DO  !kk

return

end subroutine gen_part 
