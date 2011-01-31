PROGRAM react_fast
! Slim Fast - a particle-tracking contaminant-transport code in the slim family
!
! Written by:
! Reed M. Maxwell
! International GW Modeling Center
! Hydrologic Science and Engineering Program
! Department of Geology and Geologic Engineering
! Colorado School of Mines
! Golden, CO
! 80403
! rmaxwell@mines.edu
! copyright (c) 1994-2010
!
!
! version R2.0, variable retardation (by constituent), decay/ingrowth, matrix diffusion
!
! Notable features are:
!  Particle-by-particle tracking, allowing each particle to move at
!  fastest stable pace.  Similar to 'temporal sub-cycling' methods, and
!  to Shafer-Perrini and Wilson, 1992.
!  BiLinear Velocity Interpretation for Advection Correction and Dispersion
!  ala LaBolle, et al 1996
!  Other mods re-writes by RMM as needed to streamline code.
!
! Code History:
!   7/96 - 4/97 Advection only, via analytical cell method w/ Linear
!        Velocity Interpolation
!   4/97 - 5/97 Wells, via a couple different methods developed by RMM
!   5/97  Dispersion using a Bi-Linear Interpolation of Velocity
!       for mean dispersion term (often referred to as an
!       advection cerrection factor), and the random-walk
!       dispersion term.
!   6/97        Code generalized to be fully 3D.
!   1/98        More general front-end added, code documented
!	10/99		Particle Splitting and Temporal Avgg
!	6/00		Complicated, reactive mineral retardation module added
!	11/00		Code cleaned, made more std F90
!	11/00		Mulitple Constituents added, decay between constits and ingrowth
!
!  Ideas for streamlining code come from discussions with Andrew Tompson of LLNL
!  the author owes him many thanks
!
!  "Sure Like It, Man" -AFBT
!  "The faster the better" - RMM
!
!  Refs:
!
! Maxwell and Kastenberg, SERRA, 13:1-2, 1999
! Maxwell, PhD Disert., UC, Berkeley, 1998.
! Labolle, et al., WRR , 1996
! Tompson, WRR , 1993.
! Schafer-Perrini and Wilson, WRR, 1992.
! Tompson and Gelhar, WRR, 1990.
!
! Copyright 1994-2010 Reed M. Maxwell
!
! This file is part of SLIM-Fast.
!
!    SLIM-Fast is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License.
!
!    SLIM-Fast is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SLIM-Fast in /src/gpl.txt.  If not, see <http://www.gnu.org/licenses/>.

USE NTransport
USE SLIM_FAST

IMPLICIT NONE

REAL*8   vmult,   al,at,tnext, delv,  &
    moldiff, deltx, delty, deltz,    &
	  welltnext,  wells, kd_temp,	&
	  phi_const, epsi, vga_const, vgn_const, sres_sat_const

INTEGER*4   backsl, i,       &
    concprint,welltnumb,     &
    kk,  partprint,  nw,     &
    wellprint,momprint,nt,   &
    xtent,ytent,ztent,       &
	n_constituents, ii, jj,  &
	ircheck, iv_type,phi_type, &
	npmax, give_up, press
INTEGER saturated

REAL*8,allocatable::phi(:,:,:),r(:,:,:,:),rtemp(:),half_life(:),		&
    k_det(:,:,:,:),k_att(:,:,:,:),katt_temp(:),kdet_temp(:), bulk_den(:), &
    R_temp(:,:,:)


CHARACTER (LEN=100) ::  runname, slimfile, logfile,      &
    partfile,  momfile, wellover, kxfile,kyfile,kzfile,  &
	headfile, phifile,head_list_file, time_file,   &
    vgafile, vgnfile, sresfile, vtk_file


CHARACTER(LEN=100),allocatable :: wellbtcfile(:)
CHARACTER(LEN=20),allocatable :: concfile(:)
CHARACTER(LEN=100),allocatable :: min_file(:)
CHARACTER(LEN=100),allocatable :: kd_file(:)
CHARACTER(LEN=100),allocatable :: R_file(:)
CHARACTER(LEN=100),allocatable :: mat_file(:)
CHARACTER(LEN=100),allocatable :: perm_cat_file(:)

CHARACTER (LEN=100) :: velfile, inputfile, ninputfile

DIMENSION delv(3), wells(20,10)

INTERFACE  
    SUBROUTINE pf_read(x,filename,nx,ny,nz,dx2,dy2,dz2)
    real*8  :: x(:,:,:)
    character*100 :: filename
    integer*4 :: nx
    integer*4 :: ny
    integer*4 :: nz
    real*8  :: dx2
    real*8  :: dy2
    real*8  :: dz2
    END SUBROUTINE pf_read

!SUBROUTINE slimfast(xtent,ytent,ztent,delv,al,at,							 &
!    concprint,wellprint,momprint,confile, tnext,nt,partprint,vlocfile,	         &
!    partfile,nw,well,welltnext,moldiff,welltnumb, rtard,porosity,n_constituents, &
!half_life,k_att,k_det,iv_type, press,headfile,head_list_file, time_file,kxfile,kyfile,&
!kzfile,vgafile,vgnfile,sresfile,npmax,give_up,epsi,vmult,vtk_file,saturated)  	  
!INTEGER*4                :: xtent
!INTEGER*4                :: ytent
!INTEGER*4                :: ztent
!REAL*8                   :: delv(3)
!INTEGER*4                :: n_constituents
!REAL*8                   :: al
!REAL*8                   :: at
!INTEGER*4                :: concprint
!INTEGER*4                :: wellprint
!INTEGER*4                :: momprint
!CHARACTER (LEN=20)       :: confile(:)
!REAL*8                   :: half_life(:)
!REAL*8                   :: tnext
!INTEGER*4                :: nt
!INTEGER*4                :: partprint
!CHARACTER (LEN=100)      :: vlocfile
!CHARACTER (LEN=100)      :: kxfile
!CHARACTER (LEN=100)      :: kyfile
!CHARACTER (LEN=100)      :: kzfile
!CHARACTER (LEN=100)      :: vgafile
!CHARACTER (LEN=100)      :: vgnfile
!CHARACTER (LEN=100)      :: sresfile
!Integer*4                :: press
!CHARACTER (LEN=100)      :: headfile
!CHARACTER (LEN=100)      :: head_list_file
!CHARACTER (LEN=100)      :: time_file
!CHARACTER (LEN=100)      :: partfile
!INTEGER*4                :: nw
!REAL*8                   :: well(20,10)
!REAL*8                   :: welltnext
!REAL*8                   :: moldiff
!INTEGER*4                :: welltnumb
!REAL*8                   :: rtard(:,:,:,:)
!REAL*8                   :: porosity(:,:,:)
!REAL*8                   :: k_att(:,:,:,:)
!REAL*8                   :: k_det(:,:,:,:)
!INTEGER*4				 :: iv_type
!INTEGER*4				 :: npmax
!INTEGER*4				 :: give_up
!real*8                   :: epsi
!real*8                   :: vmult 
!CHARACTER (LEN=100)      :: vtk_file
!INTEGER			 :: saturated
!end subroutine slimfast

END INTERFACE

!
!    Open Input/Output Files
!
! initialize n_constits to 1
n_constituents  = 1

!CALL STIFFEX()

READ(5,97) inputfile
READ(5,97) ninputfile
OPEN(99,FILE=inputfile,STATUS='old')

!
!    Read the input file
!

!
!    Read Header
!

READ(99,*) runname
READ(99,*) logfile
READ(99,*) saturated
! velocity input type
read(99,*) iv_type
if (iv_type == 1) READ(99,*) velfile

97 FORMAT(a100)
98 FORMAT(a20)

phi_const = 1.0d0

if (iv_type == 2) then
	read(99,*) kxfile
	read(99,*) kyfile
	read(99,*) kzfile
! @RMM Added 8-6-08 
!  need to read Van Gen alpha, n and sres filenames
! for var sat
	if (saturated == 0) then 
		read(99,*) vgafile
		read(99,*) vgnfile
		read(99,*) sresfile
        else if ( saturated == 2 ) then
		read(99,*) vga_const
		read(99,*) vgn_const
		read(99,*) sres_sat_const
        else if ( saturated == 3 ) then
		read(99,*) vga_const
		read(99,*) vgn_const
		read(99,*) sres_sat_const
	end if
! press= 0 reading hydraulic head potential
! press= 1 reading pressure (and convert internally
! when we calc v)
	read(99,*) press
	read(99,*) headfile
	read(99,*) phi_type
	if (phi_type == 2) read(99,*) phifile
	if (phi_type == 1) read(99,*) phi_const
end if

! vtype =3 means we read a list of head/pressure files with times
! for a transient run
if (iv_type == 3) then
	read(99,*) kxfile
	read(99,*) kyfile
	read(99,*) kzfile
	if (saturated == 0) then
		read(99,*) vgafile
		read(99,*) vgnfile
		read(99,*) sresfile
        else if ( saturated == 2 ) then
		read(99,*) vga_const
		read(99,*) vgn_const
		read(99,*) sres_sat_const
        else if ( saturated == 3 ) then
		read(99,*) vga_const
		read(99,*) vgn_const
		read(99,*) sres_sat_const
	end if
! press= 0 reading hydraulic head potential
! press= 1 reading pressure (and convert internally
! when we calc v)
	read(99,*) press
	read(99,*) time_file
	read(99,*) head_list_file
	read(99,*) phi_type
	if (phi_type == 2) read(99,*) phifile
	if (phi_type == 1) read(99,*) phi_const
end if


OPEN(666,FILE=trim(logfile),form='formatted',status='unknown')
WRITE(666,*)
write(666,*) ' ******************************'
WRITE(666,*) ' ** Output Log for Slim-Fast **'
WRITE(666,*) '    V4 2-10'
WRITE(666,*)
WRITE(666,*) ' Run: ',trim(runname)
WRITE(666,*) ' Log Output File: ',trim(logfile)
WRITE(666,*) ' Slim-Fast input file: ',trim(slimfile)
WRITE(666,*)
WRITE(666,*)
if (iv_type == 1) then
WRITE(666,*) 'reading V.bin files from:'
WRITE(666,*)
WRITE(666,*) 'file: ',trim(velfile)
end if
if (iv_type == 2) then
WRITE(666,*) 'calculating v, reading components'
WRITE(666,*) 'file for Kx: ',trim(kxfile)
WRITE(666,*) 'file for Ky: ',trim(kyfile)
WRITE(666,*) 'file for Kz: ',trim(kzfile)
WRITE(666,*) 'file for Head: ',trim(headfile)
if (phi_type == 2) WRITE(666,*) 'file for phi: ',trim(phifile)
if (phi_type == 1) WRITE(666,'("constant phi:",f8.6)') phi_const
end if


WRITE(666,*)
WRITE(666,*)

WRITE(666,*) ' * Domain Variable List:'
WRITE(666,*)
READ(99,*)     xtent
READ(99,*)     ytent
READ(99,*)     ztent
WRITE(666,*)
WRITE(666,*) '  Domain Size:'
WRITE(666,*)
WRITE(666,'("  nx: ",i5,",   ny: ",i5,",   nz: ",i5)' ) xtent,ytent,ztent
READ(99,*) deltx
READ(99,*) delty
READ(99,*) deltz

WRITE(*,*) 'domain, nx,ny,nz, dx,dy,dz'
WRITE(666,77) deltx,delty,deltz
77 FORMAT (' delx,dely,delz [m]: ',2(f10.3,', '),f10.3)
WRITE(666,*)

! read in number of consituents
READ(99,*) n_constituents
PRINT*, 'n_constituents'
ALLOCATE(phi(xtent,ytent,ztent),r(n_constituents,xtent,ytent,ztent),                   &
	concfile(n_constituents), wellbtcfile(n_constituents),min_file(n_constituents),    &
	rtemp(n_constituents),half_life(n_constituents),kd_file(n_constituents),          &
	k_att(n_constituents,xtent,ytent,ztent),k_det(n_constituents,xtent,ytent,ztent),  &
	kdet_temp(n_constituents),katt_temp(n_constituents),mat_file(n_constituents),     &
	perm_cat_file(n_constituents),bulk_den(n_constituents),R_file(n_constituents))
	

! read in half life
DO i=1, n_constituents
  READ(99,*) half_life(i)
  print*,'half life ',half_life(i)
END DO

! read in dispersion
READ(99,*) al
WRITE(*,*) 'al ', al
READ(99,*) at
WRITE(*,*) 'at ', at
! read in retardation
DO i=1, n_constituents
  READ(99,*) rtemp(i)
  WRITE(*,*) 'Rtemp ',rtemp(i)
  IF(rtemp(i) < -10.) THEN
    READ(99,*) min_file(i)
    WRITE(*,*) 'min_file ', min_file(i)
  END IF
    IF(rtemp(i) == -1.) THEN
    READ(99,*) kd_file(i)
	read(99,*) bulk_den(i)
    WRITE(*,*) 'kd_file ', kd_file(i)
  END IF
IF(rtemp(i) == -2.) THEN
    READ(99,*) R_file(i)
    WRITE(*,*) 'R_file ', R_file(i)
  END IF
END DO
! read in attachment/detachment, matrix diffusion or kinetic sorption
DO i=1, n_constituents
  READ(99,*) katt_temp(i)
  WRITE(*,*) 'Katt_temp ',katt_temp(i)
  READ(99,*) kdet_temp(i)
  WRITE(*,*) 'Kdet_temp ',kdet_temp(i)

  IF ((katt_temp(i) < 0.).or.(kdet_temp(i) < 0.)) THEN
    READ(99,*) mat_file(i)
    WRITE(*,*) 'mat_file ', mat_file(i)
	    READ(99,*) perm_cat_file(i)
    WRITE(*,*) 'perm_cat_file ', perm_cat_file(i)
  END IF

END DO

READ(99,*) tnext
WRITE(*,*) 'tnext ', tnext
READ(99,*) nt
WRITE(*,*) 'nt ',nt
READ(99,*) wellprint
DO i=1, n_constituents
  READ(99,*) wellbtcfile(i)
  WRITE(*,*) 'well file ',wellbtcfile(i)
END DO
READ(99,*) welltnext
READ(99,*) welltnumb
READ(99,*) concprint
if (concprint == 1 ) then
DO i=1, n_constituents
  READ(99,*) concfile(i)
  WRITE(*,*) 'conc file ',concfile(i)
END DO
end if
if (concprint == 2 ) then
 read(99,*) vtk_file
DO i=1, n_constituents
  READ(99,*) concfile(i)
  WRITE(*,*) 'conc file ',concfile(i)
END DO
end if
READ(99,*) partprint
READ(99,*) partfile
WRITE(*,*) 'part file ',partfile
READ(99,*) momprint
READ(99,*) momfile
WRITE(*,*) 'moment file ',momfile
READ(99,*) backsl
READ(99,*) moldiff
!
! read in particle limit and numerics params
!
read(99,*) npmax
read(99,*) give_up
read(99,*) epsi

WRITE(666,*)
WRITE(666,*) '  * Input Variable list for Slim-Fast'
WRITE(666,*)
WRITE(666,*)

WRITE(666,'(" Number of Consituents: ",i2)') n_constituents
WRITE(666,*)

WRITE(666,'("   alpha_l [m]: ",e12.4, "     alpha_t [m]: ",e12.4,"     mol diff [m**2]: ",e12.4)') al,at,moldiff
WRITE(666,*)
DO i=1, n_constituents
  IF (rtemp(i) > 0.) THEN
    WRITE(666,'(" Linear Retardation Factor [-]: ",e12.4)') rtemp(i)
  ELSE
	if (rtemp(i) == -1.) WRITE(666,*)' Kd File',i,': ',trim(kd_file(i))
	if (rtemp(i) == -2.) WRITE(666,*)' R PFBFile',i,': ',trim(R_file(i))
    if (rtemp(i) < -10.) WRITE(666,*)' Mineralization File',i,': ',trim(min_file(i))
  END IF




  IF ((katt_temp(i) >= 0.).and.(kdet_temp(i) >= 0.)) THEN
    WRITE(666,'(" Constant attachment [1/d]: ",e12.4)') katt_temp(i)
    WRITE(666,'(" Constant dettachment [1/d]: ",e12.4)') kdet_temp(i)
  ELSE
    WRITE(*,*) 'matrix diffusion read in from file: ', trim(mat_file(i))
    WRITE(*,*) 'perm category file: ', trim(perm_cat_file(i))
  END IF

END DO
WRITE(666,*)
WRITE(666,'(" Tnext [d]: ",e12.4)') tnext
WRITE(666,'(" number of transport steps, nt : ",i5)') nt
WRITE(666,*)

IF (wellprint == 1) wellover = 'append'
IF (wellprint == 2) wellover = 'sequential'
IF (wellprint == 0) THEN
  WRITE(666,*) ' not writing well breakthrough curves'
ELSE
 DO i=1, n_constituents
  WRITE(666,*) ' writing well BTC to file: ',TRIM(wellbtcfile(i))
  OPEN(30+i,FILE=TRIM(wellbtcfile(i)),STATUS='unknown',ACCESS=wellover)
 END DO
  WRITE(666,*) ' accessing as: ',trim(wellover)
  WRITE(666,'(" calculating breakthough every ",e12.4," [d]")') welltnext
  WRITE(666,'(" for  ",i8," steps")') welltnumb
END IF

WRITE(666,*)

IF ((concprint == 1)) THEN
 DO i=1, n_constituents
  WRITE(666,*) ' writing concentration values to CNB file: ', trim(concfile(i))
 END DO
  WRITE(666,*) ' note that 001,002, etc will be appended to the file to'
  WRITE(666,*) ' correspond to each tnext.'
ELSEIF (concprint == 2) then
  WRITE(666,*) ' writing concentration values to VTK file: ', trim(vtk_file)
 DO i=1, n_constituents
  WRITE(666,*) ' using conc header: ',i, trim(concfile(i))
 END DO
  WRITE(666,*) ' note that 00000000, 00000001,00000002, etc  and .vtk will be appended to the file to'
  WRITE(666,*) ' correspond to each tnext.'
else 
  WRITE(666,*) ' not writing out concentration values'
END IF

WRITE(666,*)

IF (momprint == 0) THEN
  WRITE(666,*) ' not writing plume spatial moments'
ELSE
  WRITE(666,*) ' writing plume moments to file: ',trim(momfile)
  OPEN(41,FILE=momfile,STATUS='unknown',ACCESS='append')
END IF
WRITE(666,*)
IF (backsl == 1) WRITE(666,*) ' using reverse SL tracing'


  WRITE(666,*) 
  WRITE(666,*) ' Maximum allocated particles: ',npmax
  WRITE(666,*) ' Maximum particle steps per time loop: ',give_up
  WRITE(666,'(" Epsilon/Drop value (below this assume zero): ",e12.4)') epsi
  write(666,*)
!    Initialize Variables
!
if (phi_type == 1) then
phi = phi_const
else if(phi_type ==2) then
call pf_read(phi,phifile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
end if

DO i=1, n_constituents
  IF (rtemp(i) > 0.) THEN
    print*,i,rtemp(i)
    R(i,1:xtent,1:ytent,1:ztent) = rtemp(i)
  ELSE

	IF (rtemp(i) < -10.) THEN
		print*,' calling min_react'
		CALL react_min_lu(xtent, ytent, ztent, phi, r, min_file(i),n_constituents,i)
		print*,' back from min_react'
	end if
	IF (rtemp(i) == -1.) THEN
		print*,' reading kd file'
		ircheck = 1
		if (ircheck.eq.1) then
		open (199,file='rcheck.c1.asc.txt')
		write(199,*) xtent,ytent,ztent
		end if

		open (990,file=trim(kd_file(i)))
		read (990,*)
			do ii = 1, xtent
				do jj = 1, ytent
					do kk = 1, ztent
					read (990,*) kd_temp
					R(i,ii,jj,kk) = kd_temp*(Bulk_den(i)/Phi(ii,jj,kk)) + 1.d0 
					if (ircheck == 1) write(199,*) r(i,ii,jj,kk)
					end do
				end do
			end do
		close (990)
		if (ircheck == 1) close(199)
		print*,' done reading kd file'
	END IF
IF (rtemp(i) == -2.) THEN
		print*,' reading R from pfb file'
		!print*, 'R? ', allocated(R_temp)

	!	allocate(R_temp(xtent,ytent,ztent))
		!print*, 'R? ', allocated(R_temp)
    call pf_read(R(i,:,:,:),R_file(i),xtent,ytent,ztent,delv(1),delv(2),delv(3))		
    !R(i,:,:,:) = R_temp
    !deallocate(R_temp)
    !print*, 'R? ',allocated(R_temp)

		print*,' done reading R from file'
	END IF
  END IF
END DO

DO i=1, n_constituents
  IF ((katt_temp(i) >= 0.).and.(kdet_temp(i) >= 0.)) THEN
  print*,i,kdet_temp(i)
! katt change 8-05
    k_att(i,1:xtent,1:ytent,1:ztent) = katt_temp(i)
	  print*,i,kdet_temp(i)
    k_det(i,1:xtent,1:ytent,1:ztent) = kdet_temp(i)
  ELSE
  		print*,' calling mat diff lu'
		CALL  mat_diff(xtent, ytent, ztent, mat_file(i),perm_cat_file(i),n_constituents,k_att,k_det,i)
		print*,' back from mat diff lu'
  END IF
END DO


delv(1) = deltx
delv(2) = delty
delv(3) = deltz

!     Porosity = phi

33   CONTINUE

READ(99,*) nw

DO   kk = 1, nw
  READ(99,*) wells(kk,1),wells(kk,2),wells(kk,3),wells(kk,4), wells(kk,5)
END DO


CALL ReadNTransPars( ninputfile, xtent, ytent, ztent )

!CALL RACT() 

vmult = 1.0D0
IF (backsl == 1) vmult = -1.0D0

  
256 FORMAT (e15.7)

!
! Call Slim-Fast
!

WRITE(666,*)
WRITE(666,*) 'Calling Slim-Fast'
WRITE (*,*) 'Calling Slim-Fast'
flush(666)


!
!    call Slim-Fast
!


CALL slimfast(xtent,ytent,ztent,delv,al,at,                            &
    concprint,wellprint,momprint,concfile, tnext,nt,partprint,velfile,   &
    partfile,nw,wells,welltnext,moldiff,welltnumb, r,phi,n_constituents, &
	half_life,k_att,k_det,iv_type,press, headfile,head_list_file, time_file, &
	 kxfile,kyfile,kzfile,vgafile,vgnfile,sresfile,npmax, give_up, epsi,vmult, &
     vtk_file,saturated, vga_const, vgn_const, sres_sat_const)
	

print*, 'finished'
WRITE(666,*)
write(666,*) '*** run complete!  ***'

close(666)
END PROGRAM react_fast
