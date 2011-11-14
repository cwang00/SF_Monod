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
                     time_begin, time_end, porosity, sat ) 

INTEGER*4 :: ii,jj,kk, icat, jcat, kcat, n_ic_cats,  &
           num_of_parts, iii, timestep_num, ir,     &
          constitute_num, npmax
INTEGER*4, DIMENSION(:,:,:) :: ic_cat_num
INTEGER*4, DIMENSION(:,:) :: ip
INTEGER*4, DIMENSION(:) :: ic_cats, ic_part_dens

real*8  :: time_begin, time_end, rantemp, part_mass, mass
REAL*8, DIMENSION(:,:) :: particle, ic_conc
REAL*8, DIMENSION(:,:,:) :: sat, porosity
REAL*8, DIMENSION(:, :,:,:) :: vel
REAL*4, DIMENSION(:,:,:,:) ::  current_conc
REAL*8, DIMENSION(:) :: delv

INTERFACE
REAL*8  FUNCTION ran1(idum)
INTEGER*4 :: idum
END FUNCTION ran1
end interface
ir = 21
!  PRINT*,' 3rd ic: loop thru ic'
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
          DO iii = 1, n_ic_cats
           IF( ic_part_dens( iii ) .GT. 0 ) THEN
            IF( timestep_num == 1 ) THEN
               mass = ic_conc( 1, iii ) * delv( 1 ) * delv( 2 ) * delv( 3 ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
            ELSE
               mass = ic_conc( 1, iii ) * delv( 2 ) * delv( 3 ) &
                       * ABS( vel(1, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               mass = mass +  ic_conc( 1, iii ) * delv( 1 ) * delv( 3 ) &
                       * ABS( vel(2, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 *  porosity( ii,jj,kk ) * sat(ii,jj,kk) 
               mass = mass +  ic_conc( 1, iii ) * delv( 1 ) * delv( 2 ) &
                       * ABS( vel(3, ii,jj,kk) ) * (time_end - time_begin ) &
                         * 1000 * porosity( ii,jj,kk ) * sat(ii,jj,kk) 
              ENDIF
            part_mass =  mass / ic_part_dens( iii )
            IF (ic_cat_num(ii,jj,kk) == ic_cats(iii)) THEN

! assign particles to active locations
! check to see if mass of ic > 0
             IF (ic_conc(timestep_num,iii) >                               &
                         current_conc( constitute_num, ii, jj, kk) ) THEN
                !
                ! assumes ic_mass_or_conc is mg/l, delv in meters, mass is in
                ! mg
                !
!                mass =  mass - current_conc(constitute_num, ii, jj , kk )   &
!                        * delv( 1 ) * delv( 2 ) * delv( 3 ) * 1000 *        &
!                          porosity( ii, jj,kk ) * sat(ii,jj,kk) 
              DO  WHILE(  mass > 0.0 )
	        particle(num_of_parts,1) = delv(1)*ran1(ir) + DBLE(ii-1)*delv(1)
                particle(num_of_parts,2) = delv(2)*ran1(ir) + DBLE(jj-1)*delv(2)
                particle(num_of_parts,3) = delv(3)*ran1(ir) + DBLE(kk-1)*delv(3)
                particle(num_of_parts,4) = part_mass
                particle(num_of_parts,5) = num_of_parts 
                particle(num_of_parts,7) = time_begin +                &
                                           ran1(ir)*(time_end - time_begin)
                particle(num_of_parts,9) = ic_cat_num(ii,jj,kk)
                ip(num_of_parts,1) = constitute_num
                ip(num_of_parts,2) = 1
                IF (num_of_parts > npmax) THEN
                  PRINT*, 'max num particles exceeded.  increase npmax parameter ',npmax
                  WRITE(666,*) 'Max Particles Exceeded!  Increase NPMAX:',npmax
                  STOP
                END IF ! max parts ?
               
                num_of_parts = num_of_parts + 1
                mass = mass - part_mass
                !
                ! update concentration
                !
                current_conc( constitute_num, ii, jj, kk ) =              &
                             current_conc( constitute_num, ii, jj, kk ) + &
                             part_mass / ( delv(1) * delv(2) * delv(3)    &
                             * 1000 * porosity( ii, jj, kk ) *            &
                             sat( ii, jj, kk ) )
                                  
                END DO  !DO WHILE
              
             END IF  ! ic_conc > current_conc ? 
            END IF  ! at a ic node
           END IF 
          END DO !iii
      END DO  !ii
    END DO  !jj
  END DO  !kk

return

end subroutine gen_part 
