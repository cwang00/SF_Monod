!
! Generate particles for given concentration at indicator locations
!
! Zhengtao Cui
! Oct. 9, 2011
!
subroutine gen_part( icat, jcat, kcat, ic_cat_num, n_ic_cats, ic_cats,         &
                     ic_conc, ic_part_dens, particle, ip,                      &
                     delv, vel, total_num_of_parts,                            &
                     timestep_num, constitute_num, npmax, current_conc,        &
                     time_begin, time_end, porosity, sat, bnd_cnd ) 
   USE Particles

INTEGER*4 :: ii,jj,kk, icat, jcat, kcat, n_ic_cats,  &
           total_num_of_parts, iii, timestep_num, ir,     &
          constitute_num, npmax, part
INTEGER*4, DIMENSION(:,:,:) :: ic_cat_num
INTEGER*4, DIMENSION(:,:) :: ip
INTEGER*4, DIMENSION(:) :: ic_cats, ic_part_dens

real*8  :: time_begin, time_end, rantemp, part_mass, mass, current_mass
REAL*8, DIMENSION(:,:) :: particle, ic_conc
REAL*8, DIMENSION(:,:,:) :: sat, porosity
REAL*8, DIMENSION(:, :,:,:) :: vel
REAL*4, DIMENSION(:,:,:,:) ::  current_conc
REAL*8, DIMENSION(:) :: delv
CHARACTER (LEN=20) :: bnd_cnd

!write(*,*) bnd_cnd
ir = 21
!  PRINT*,' 3rd ic: loop thru ic'
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
          DO iii = 1, n_ic_cats
           IF( ic_part_dens( iii ) .GT. 0 ) THEN
            IF( timestep_num == 1 ) THEN
               IF( bnd_cnd == 'const_flux_septic' ) THEN
                  mass = ic_conc( 1, iii ) * delv( 1 ) * delv( 2 ) *  &
                         0.0055 * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk)
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
               ELSE IF ( bnd_cnd == 'const_flux_x' ) THEN
               mass = ic_conc( timestep_num, iii ) * delv( 2 ) * delv( 3 ) &
                       * ABS( vel(1, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               ELSE IF ( bnd_cnd == 'const_flux_y' ) THEN
               mass = mass +  ic_conc( timestep_num, iii ) * delv( 1 ) * delv( 3 ) &
                       * ABS( vel(2, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 *  porosity( ii,jj,kk ) * sat(ii,jj,kk) 

               ELSE IF ( bnd_cnd == 'const_flux_z' ) THEN
               mass = mass +  ic_conc( timestep_num, iii ) * delv( 1 ) * delv( 2 ) &
                       * ABS( vel(3, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               ELSE IF ( bnd_cnd == 'const_conc' ) THEN
                 mass = ic_conc( timestep_num, iii )              &
                         * delv( 1 ) * delv( 2 ) * delv( 3 )      &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 

                 IF ( mass .LT. 0.0 ) mass = 0.0

               ELSE IF( bnd_cnd == 'const_flux_septic' ) THEN
                  mass = ic_conc( timestep_num, iii )             &
                         * delv( 1 ) * delv( 2 ) *                &
                         0.0055 * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) &
                       * (time_end - time_begin ) 
!                write(*,*) 'poros = ', porosity( ii,jj,kk )
!                write(*,*) 'sat = ', sat( ii,jj,kk )
!                 write(*,*) 'mass = ', mass         
               ELSE
                WRITE(*,*) 'ERROR: non-known boundary condition: ', bnd_cnd
                stop
               ENDIF !  IF ( bnd_cnd == 'const_flux' ) THEN
             ENDIF  ! IF( timestep_num == 1 ) THEN
             part_mass =  mass / ic_part_dens( iii )

            IF (ic_cat_num(ii,jj,kk) == ic_cats(iii)) THEN
!               current_mass = current_conc(constitute_num, ii, jj,kk ) &
!                    * delv( 1 ) * delv( 2 ) * delv( 3 ) &
!                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
              !remove existing particles
              IF( number_of_parts(iii, ii,jj,kk) .GT. 0 ) THEN
                DO I = 1, number_of_parts(constitute_num, ii,jj,kk) 
                  totalremoved = totalremoved + 1
                  removed( totalremoved ) = &
                        part_numbers( constitute_num, ii,jj,kk )%arr(I)
                ENDDO
                number_of_parts(constitute_num, ii,jj,kk) = 0
              ENDIF

! assign particles to active locations
! check to see if mass of ic > 0
              IF ( ABS( mass ) .GT. 0.00000000001 ) THEN
               DO I = 1, ic_part_dens(iii) 
                IF ( totalremoved .GT. 0 ) THEN
                   part = removed( totalremoved); 
                   totalremoved = totalremoved - 1;
                ELSE
                   part = total_num_of_parts;
                   total_num_of_parts = total_num_of_parts + 1
                   IF (total_num_of_parts > npmax) THEN
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
	        ip(part,6) = -1  
               
                !mass = mass - part_mass
              END DO  !DO I = 1, ic_part_dens
              
                !
                ! update concentration
                !
!                current_conc( constitute_num, ii, jj, kk ) =        &
!                     current_conc( constitute_num, ii, jj, kk )  +      &
!                           ic_conc( timestep_num, iii ) * 1000
                current_conc( constitute_num, ii, jj, kk ) =        &
                           ic_conc( timestep_num, iii ) * 1000

!                number_of_parts(constitute_num, ii,jj,kk) =      &
!                   number_of_parts(constitute_num, ii,jj,kk) +      &
!                     ic_part_dens(iii)
                number_of_parts(constitute_num, ii,jj,kk) =      &
                     ic_part_dens(iii)
            END IF  ! IF ( ABS( mass ) .GT. 0.00000000001 ) THEN

           ENDIF  !IF (ic_cat_num(ii,jj,kk) == ic_cats(iii)) THEN
          END IF ! IF( ic_part_dens( iii ) .GT. 0 ) THEN
         END DO !iii
      END DO  !ii
    END DO  !jj
  END DO  !kk

return

end subroutine gen_part 
