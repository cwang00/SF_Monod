MODULE SLIM_FAST

USE NTransport
USE VGTransport
USE Particles

CONTAINS

SUBROUTINE slimfast(xtent,ytent,ztent,delv,al,at, &
    concprint,wellprint,momprint,confile, tnext,nt,partprint,vlocfile,         &
    partfile,nw,well,welltnext,moldiff,welltnumb, rtard,porosity,n_constituents, &
 half_life,k_att,k_det,iv_type, press,headfile,head_list_file,time_file, kxfile,kyfile,  &
 kzfile,vgafile,vgnfile,sresfile,npmax,give_up,epsi,vmult,vtk_file, saturated, &
 vga_const, vgn_const, sres_sat_const, modelname)
! Slim-Fast Main Routine
! written by Reed M. Maxwell
! rmaxwell@mines.edu
!
! Variables
!
! P(n,prop)   particle array (real particle properties)
!  where:
!  n = particle number, ranges from 1 to np
!  prop = 1  X location
!  prop = 2  Y location
!  prop = 3  Z location
!  prop = 4  Mass
!  prop = 5  Particle Number
!  prop = 6  well number where particle is entrained
!  prop = 7  Time of Particle  
!  prop = 8  Tau, time of flight of particle
!  prop = 9  current node geotype
!  prop = 10 mag(V).
!
! iP(n,prop) integer particle properties
!  prop = 1 constituent (=1,n_constituents)
!  prop = 2 state (matrix/aqueous)
!  ** currently iP(n,2) is multiplied by any/all moves, 
!     so iP(n,2) = 1 allows particle to move and
!        iP(n,2) = 0 keeps particle fixed
!  prop = 3 whether or not particle has crossed plane
!  prop = 4 is 0 if the particle is not near a reflection boundary
!          and l = 1,2 or 3 depending upon which bddy
!
!  variables that are dependant on constituent
! 
!  Retardation
!  R(m,i,j,k)
!  m = 1, 2, 3... corresponding to iP(n,1)
!
!	Velocity (as V/R)
! V(l,i,j,k)  x,y,z velocities
!  where
!  l = 1, 2, 3, for x, y, z respectively
!  i,j,k are the gridpoint coords.
!	
! C(m,i,j,k)  concentration interpolated from particles,
!  m = 1, 2, 3... corresponding to iP(n,1)
!
! tad(l) advection time, the time to move particle from current location to
!  edge of next cell.  Always follows direction with greatest velocity.
!  l works like above, with l = 4 being the max velocity
!
! tnext the next reporting time. All particles
! allowed to act 'independently'
!  of one-another.  This time is the next 're-grouping' time.
!	for the transient case, TNEXT is specified at each velocity update
!
! ploc(l) is the integer particle location, l works as above
!
! delv(l) the velocity grid-spacing, l works as above
!
! delc(l) the concentration interpolation grid spacing, l works as above
!
! Al, At, Alt are the longitudinal, transverse and difference of the two, respectively
!
!
!
! the well array works as follows:
!
! W(wn,i)
!  where
!   wn = the well number and
!    i = 1 is the x location of the well
!    i = 2 is the y location of the well
!    i = 3 is the top z screen location of the well
!    i = 4 is the bottom screen location of the well
!    i = 5 is the pumping rate, Q, in [L^3/time]
!
! BTC(wn,time,constit) is the mass-based breakthrough
!  curve (mass accumulated over time)
!  where
!   wn = the well number and
!   time = the time of breakthrough
!   constit = the constituent breaking through
!
! WELLTAVG is the averaging time for each step of the well BT Curve
!
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



IMPLICIT NONE
INTEGER*4                :: xtent
INTEGER*4                :: ytent
INTEGER*4                :: ztent
REAL*8                   :: delv(3)
INTEGER*4                :: n_constituents
REAL*8                   :: al
REAL*8                   :: at
INTEGER*4                :: concprint, velprint
INTEGER*4                :: wellprint
INTEGER*4                :: momprint
!CHARACTER (LEN=20)       :: confile(n_constituents)
CHARACTER (LEN=20)       :: confile(:)
!REAL*8                   :: half_life(n_constituents)
REAL*8                   :: half_life(:)
REAL*8                   :: tnext
INTEGER*4                :: nt
INTEGER*4                :: partprint
CHARACTER (LEN=100)      :: vlocfile
CHARACTER (LEN=100)      :: kxfile
CHARACTER (LEN=100)      :: kyfile
CHARACTER (LEN=100)      :: kzfile
CHARACTER (LEN=100)      :: vgafile
CHARACTER (LEN=100)      :: vgnfile
CHARACTER (LEN=100)      :: sresfile
INTEGER*4                :: press
CHARACTER (LEN=100)      :: headfile
CHARACTER (LEN=100)      :: head_list_file
CHARACTER (LEN=100)      :: time_file
CHARACTER (LEN=100)      :: partfile
INTEGER*4                :: nw
REAL*8                   :: well(20,10)
REAL*8                   :: welltnext
REAL*8                   :: moldiff, mldif
INTEGER*4                :: welltnumb
REAL*8                   :: rtard(:,:,:,:)
!REAL*8                   :: rtard(n_constituents,xtent,ytent,ztent)
!REAL*8                   :: porosity(xtent,ytent,ztent)
REAL*8                   :: porosity(:,:,:)
!REAL*8                   :: k_att(n_constituents,xtent,ytent,ztent)
REAL*8                   :: k_att(:,:,:,:)
!REAL*8                   :: k_det(n_constituents,xtent,ytent,ztent)
REAL*8                   :: k_det(:,:,:,:)
INTEGER*4                :: iv_type
INTEGER*4                :: npmax
INTEGER*4                :: give_up
real*8                   :: epsi
real*8                   :: vmult 
CHARACTER (LEN=100)      :: vtk_file
INTEGER 		 :: saturated

INTEGER, PARAMETER       ::  massive_debug = 0, part_conc_write = 0, recycle_well = 0

REAL*8  time,  delc(5),vp(5), cellv, tloop,lamda,Prob,Zhl, &
    trgp,  junk, cx, cy,dx,dy,vn,vnp1,vnm1,ptemp,      &
    fact,  a(5),b(5), a1(5), stuck, pi, btp(10,8500,20), &
    del2v(5),btc(10,8500,20),welltavg(20),ivolume, &
    vbl(4,4,4,4), alt, z(6), betad, coic, &
    rxx, ryy, sxx, sxy, syy, vxz, dtw,vpw,wellq,  &
    nstepav,  dtr,  bdyx1, bdyx2, bdyy1, bdyy2, &
    bdyz1, bdyz2, npact, xlic,xuic,ylic,yuic,zlic,zuic, &
    wellx,welly,wellzt,wellzb,  &
     rzz, szz, dxx, dxy, dxz,xplaneloc(10), &
    dyz, dzz, dzx, dzy, vxx, vzz, vyy, vyz, cz, dz,vxy, &
    dyx, dyy, dlimit, ddx, ddy, ddz, timenext, &
    decayic, dectimic, dicnstep, rand, dt_lamda, &
     grad_phi_x, grad_phi_y, grad_phi_z, &
     t_att, t_prev, p_av1, p_av2,t_big_step, pulse_vol, rw_flux(2),mass_rec(50),rstar

INTEGER*4   ici(10), nplane,  wellnpart, po(4,4),planeloc, &
    ll, done, iwjb, ijunk, iwjt,movedir,welloc(20,10),iplane, &
    iwflag(13), tempavg, n_ic_cats(13),ic_part_dens(13,20), icnpart, &
    ic_cats(13, 20),ts1, iii, jjj, npcurrent,partsplit, ts, planedir(10), &
 idecsteps, ii, jj, kk, ind_ic_catagory, icat, jcat, kcat, ic_cat_count, &
   kic, n_ic_timesteps, loop2, iic, cell_sink,boundary_cond(3,2),well_r(2,10000,3), &
        rw_ind(2), n_rw_ind(2), n_rw,imax



REAL*8  tad(6), move,extemp,tdd(8),f(5),fbl(5), rand2,pdist,minpdist
REAL*4,allocatable::c(:,:,:,:)
REAL*4 cps, cmin

REAL*8,allocatable::v(:,:,:,:),mass(:,:), P(:,:),sat(:,:,:),scxyz(:,:),lastprint(:,:), ic_mass_or_conc(:, :, :), ic_time_begin(:,:), ic_time_end(:,:)
REAL*8,allocatable::hkx(:,:,:), hky(:,:,:),hkz(:,:,:),vga(:,:,:),vgn(:,:,:),sres(:,:,:)
INTEGER*4,allocatable::ic_cat_num(:,:,:,:),iP(:,:),irP(:),iprP(:,:), &
                       ic_cat_num_bg(:, :, :, :)

INTEGER*4 l,i,j,k,n,ploc(4),np,domax(4),numax,loc3p,it, &
    passes, locstuck(5),  itemp, wc,av(5), &
    tp1,tp2,tp3, tloc(5), inbounds, wellfile

integer*4 ir, icpartdiv, nr, mp

integer, allocatable :: ipwell(:,:)

CHARACTER (LEN=6) :: dotit

CHARACTER (LEN=6) :: fadd
CHARACTER (LEN=4) :: confile2
CHARACTER (LEN=20) :: confile1
CHARACTER (LEN=40) :: filename, format_desc
CHARACTER (LEN=40),allocatable ::planefile(:)
CHARACTER (LEN=100) :: ind_ic_file

REAL*8 :: vga_const, vgn_const, sres_sat_const, bnd_Xup, bnd_Xdown,    &
                                                bnd_Yup, bnd_Ydown,    &
                                                bnd_Zup, bnd_Zdown
REAL*8 :: oldloc(3), origloc(3)
INTEGER*4 :: Xup_Ref, Xdown_Ref, Yup_Ref, Ydown_Ref, Zup_Ref, Zdown_Ref
CHARACTER (LEN=20) :: modelname
interface

subroutine v_calc(v,hkx,hky,hkz,vga,vgn,sres,headfile,phi,delv,nx,ny,nz,press,sat)
real*8  :: v(:,:,:,:),hkx(:,:,:), hky(:,:,:),hkz(:,:,:),vga(:,:,:),vgn(:,:,:),sres(:,:,:)
character(100) :: headfile
real*8  :: phi(:,:,:)
real*8  :: sat(:,:,:)
real*8 :: delv(3)
integer*4 :: nx
integer*4 :: ny
integer*4 :: nz
integer*4 :: press
end subroutine v_calc

subroutine v_calc_const_sat(v,hkx,hky,hkz,vga,vgn,headfile,phi,delv,nx,ny,nz,press,sat) 
real*8  :: v(:,:,:,:),hkx(:,:,:), hky(:,:,:),hkz(:,:,:),vga(:,:,:),vgn(:,:,:)
character(100) :: headfile
real*8  :: phi(:,:,:)
real*8  :: sat(:,:,:)
real*8 :: delv(3)
integer*4 :: nx
integer*4 :: ny
integer*4 :: nz
integer*4 :: press
end subroutine v_calc_const_sat

SUBROUTINE cbin_write(x,filename,ixlim,iylim,izlim,dx,dy,dz)
REAL*4    :: x(:,:,:)
CHARACTER (LEN=40)     :: filename
INTEGER*4 :: ixlim
INTEGER*4 :: iylim
INTEGER*4 :: izlim
REAL*8                 :: dx
REAL*8                 :: dy
REAL*8                 :: dz
end subroutine cbin_write

SUBROUTINE cbin_write_real8(x,filename,ixlim,iylim,izlim,dx,dy,dz)
REAL*8    :: x(:,:,:)
CHARACTER (LEN=40)     :: filename
INTEGER*4 :: ixlim
INTEGER*4 :: iylim
INTEGER*4 :: izlim
REAL*8                 :: dx
REAL*8                 :: dy
REAL*8                 :: dz
end subroutine cbin_write_real8

SUBROUTINE vtk_write(time,x,conc_header,ixlim,iylim,izlim,dx,dy,dz,icycle,n_constituents,vtk_file)
real*8                 :: time
REAL*4    :: x(:,:,:,:)
!REAL*4    :: x(ixlim,iylim,izlim)
CHARACTER (LEN=20)     :: conc_header(:)
INTEGER*4 :: ixlim
INTEGER*4 :: iylim
INTEGER*4 :: izlim
REAL*8                 :: dx
REAL*8                 :: dy
REAL*8                 :: dz
INTEGER                :: icycle
INTEGER*4              :: n_constituents
CHARACTER (LEN=100)      :: vtk_file
end subroutine vtk_write

SUBROUTINE gnuplot_write(x,ixlim,iylim,izlim,icycle,n_constituents,bg, vtk_file)
REAL*4    :: x(:,:,:,:)
REAL, DIMENSION(:) :: bg
INTEGER*4 :: ixlim
INTEGER*4 :: iylim
INTEGER*4 :: izlim
INTEGER                :: icycle
INTEGER*4              :: n_constituents
CHARACTER (LEN=100)      :: vtk_file
end subroutine gnuplot_write

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
 
subroutine gen_part( icat, jcat, kcat, ic_cat_num, n_ic_cats, ic_cats,      &
                     ic_conc, ic_part_dens, particle, ip,                   &
                     delv, vel, num_of_parts,                             &
                     timestep_num, constitute_num, npmax, current_conc,     &
                     time_begin, time_end, porosity, sat ) 
    INTEGER*4 :: icat, jcat, kcat, n_ic_cats,       &
                 num_of_parts, timestep_num,              &
                 constitute_num, npmax
    INTEGER*4, DIMENSION(:,:,:) :: ic_cat_num
    INTEGER*4, DIMENSION(:,:) :: ip
    INTEGER*4, DIMENSION(:) :: ic_cats, ic_part_dens

    real*8  :: time_begin, time_end, cellv
    REAL*8, DIMENSION(:,:) :: particle, ic_conc
    REAL*8, DIMENSION(:,:,:) :: sat, porosity
    REAL*8, DIMENSION(:, :,:,:) :: vel
    REAL*4, DIMENSION(:,:,:,:) ::  current_conc
    REAL*8, DIMENSION(:) :: delv

END SUBROUTINE gen_part

REAL*8  FUNCTION ran1(idum)
INTEGER*4 :: idum
END FUNCTION ran1

end interface


PRINT*, 'hey made it here'
PRINT*, ' tnext ', tnext

tempavg = 0
ir = -219191
ir = 21
imax = max(xtent,ytent,ztent)

allocate(v(3,xtent,ytent,ztent), c(n_constituents,xtent,ytent,ztent), &
       ic_cat_num(13,xtent,ytent,ztent),                              &
       ic_cat_num_bg(13,xtent,ytent,ztent),                           &
       planefile(n_constituents),    &
       mass(n_constituents+1,nt+10),  &
       P(npmax,10),iP(npmax,10),ipwell(npmax,20),sat(xtent,ytent,ztent),irP(npmax)  , &
                scxyz(3,imax),hkx(xtent,ytent,ztent),hky(xtent,ytent,ztent), &
                 hkz(xtent,ytent,ztent),vga(xtent,ytent,ztent),vgn(xtent,ytent,ztent),  &
                 sres(xtent,ytent,ztent), iprP(npmax,2), lastprint(npmax,3)) 

 allocate( ic_time_begin(13, 1000  ), ic_time_end( 13, 1000 ), &
            ic_mass_or_conc( 13, 1000, 20 ) )

c = 0.d0
mass = 0.d0
sat = 1.0d0
boundary_cond = 0
!boundary_cond(2:3,:) = 1
!boundary_cond(3,1) = 1
scxyz = 1.0d0
 
! variable dx,dy,dz HARD WIRED into velocity right now, really should
! make this more permanent later...
!scxyz(1,1:11) = 75.0d0
!scxyz(1,561:domax(1)) = 75.0d0
 
!scxyz(2,1:9) = 50.0d0
!scxyz(2,159:domax(2)) = 50.0d0
 
!scxyz(3,1:10) = 25.0d0


! Open I/O files

if (iv_type == 1) then
OPEN (16,FILE=trim(vlocfile),FORM='unformatted',access='stream')
end if
if (iv_type == 2) then
t_big_step = tnext
end if
if (iv_type == 3) then
OPEN (16,FILE=trim(time_file))
OPEN (17,FILE=trim(head_list_file))
end if
ipwell = 1



fadd = ''
nr = 1
IF (nr == 1) wellfile = 31
IF (nr == 2) wellfile = 32
IF (nr == 3) wellfile = 33
IF (nr == 4) wellfile = 34

dxx = 0.d0
dxy = 0.d0
dxz = 0.d0
dx = 0.d0
ddx = 0.d0
dyy = 0.d0
dyx = 0.d0
dyz = 0.d0
dy = 0.d0
ddy = 0.d0
dzz = 0.d0
dzx = 0.d0
dzy= 0.d0
dz = 0.d0
ddz = 0.d0

ici(1) = 2
ici(2) = 7
ici(3) = 2
PRINT*, 'flag 1'
dlimit = 0.5D0
dlimit = 0.1D0

IF(partprint >= 1) THEN
  !OPEN(13,FILE=partfile,STATUS='unknown')
  OPEN (113,FILE='end_'//partfile,STATUS='unknown')
  OPEN (213,FILE='wt_'//partfile,STATUS='unknown')
  OPEN (313,FILE='ls_'//partfile,STATUS='unknown')
END IF

!! @RMM space left here for writing streamline grid information
!!
IF(partprint == 2) THEN
 ! WRITE(13,*) 'x,y,z, time,TOF,unit,flux'
END IF


OPEN (667,file='slim.warn.txt',status='unknown')

WRITE(666,*)
WRITE(666,*) 'Slim-Fast Log File'
WRITE(666,*)


WRITE(667,*)
WRITE(667,*) 'Slim-Fast Warning and Error File'
WRITE(667,*)


!
! First we read in concentrations
!
! Set up initial parameters
!
! test where we set ip=1
ip = 1
confile2='.cnb'
np = 1
! set riP() = 0 - this is the well recycle counter
irP = 0
time = 0
! set iprP() = 0 - this is a print counter that flags endprint and print for part loc
iprP = 0
! set up print-position array
lastprint = 0.d0
pdist = 0.0
minpdist = min(delv(1),delv(2),delv(3))
print*, 'min pdist', minpdist



WRITE(666,*)
WRITE(666,*) ' Well Averaging Time: ', welltnext

welltavg(1) = welltnext
welltavg(2) = welltnext
welltavg(3) = welltnext
welltavg(4) = welltnext
welltavg(5) = welltnext
welltavg(6) = welltnext
welltavg(7) = welltnext
welltavg(8) = welltnext
welltavg(9) = welltnext
welltavg(10) = welltnext
welltavg(11) = welltnext
welltavg(12) = welltnext

t_prev = 0.0D0
pi = 3.14159265358979D0

numax = 3

del2v(1) = delv(1)/2.0D0
del2v(2) = delv(2)/2.0D0
del2v(3) = delv(3)/2.0D0
delc(1) = delv(1)
delc(2) = delv(2)
delc(3) = delv(3)
cellv = delc(1)*delc(2)*delc(3)

WRITE(666,*) ' velocity nodal spacing: '
WRITE(666,'(" dx,dy,dz: ",3(f7.2,"  "))') delv(1),delv(2),delv(3)
WRITE(666,'(" vol of cell: ",e12.4)') delv(1)*delv(2)*delv(3)
WRITE(666,*)
WRITE(666,*) ' concentration nodal spacing:'
WRITE(666,'( " dx,dy,dz: ",3(f7.2,"  "))') delc(1),delc(2),delc(3)
WRITE(666,'(" vol of cell: ",e12.4)') cellv


rxx = 0.0D0
ryy = 0.0D0
rzz = 0.0D0
sxx = 0.0D0
sxy = 0.0D0
syy = 0.0D0
szz = 0.0D0

PRINT*, 'flag 2'

!
!    Clear out well BTC
!

PRINT*, 'nw= ', nw

btc = 0.
btp = 0.

domax(1) = xtent
domax(2) = ytent
domax(3) = ztent

PRINT*, 'flag 2.5'

PRINT*, 'nx ', xtent, domax(1)
PRINT*, 'ny ', ytent, domax(2)
PRINT*, 'nz ', ztent, domax(3)

DO  n = 1, npmax
  DO  l = 1, 3
    p(n,l) = 0.0D0
  END DO
  p(n,4) = 0.
  p(n,7) = 0.
END DO

!
! Set up Kron Delta, but we'll call it po(ij)
!

DO  i=1,3
  DO  j=1,3
    IF (i == j) THEN
      po(i,j) = 1
    ELSE
      po(i,j) = 0
    END IF
  END DO
END DO

!
! Set up well locations from the well file
!

WRITE(666,*)
WRITE(666,*) ' Well Info: '
DO  i=1,nw
  DO  l=1,3
    welloc(i,l) = IDINT(well(i,l)/delv(l)) + 1
  END DO
  welloc(i,4) = IDINT(well(i,4)/delv(3))  + 1
	WRITE(666,'("X: ",e12.4," Y: ",e12.4," Z1: ",e12.4," Z2: ",e12.4)') well(i,1:4) 
  WRITE(666,*) i,': i ',welloc(i,1),' j ',welloc(i,2),' k1 '		&
	,welloc(i,3),' k2',welloc(i,4)
  WRITE(666,'(i5,": Q= ",e12.4)') i,well(i,5)
END DO

! read in plane info

WRITE(666,*)
WRITE(666,*) ' Monitoring Planes '
DO i = 1,n_constituents
  READ(99,*) planefile(i)
  WRITE(666,*) ' planes written to file: ',trim(planefile(i)),' for constituent', i 
ENDDO
READ(99,*) nplane
WRITE(666,*) ' Number of Planes:' , nplane


DO i = 1,nplane
  READ(99,*) planedir(i)
  READ(99,*) xplaneloc(i)
  WRITE(666,*) i,' direction: ',planedir(i),' location '		&
	,xplaneloc(i)
  
END DO

if (iv_type == 1) then
        if( saturated == 3 ) then
	 vga=vga_const
	 vgn=vgn_const
         sat = sres_sat_const
        end if
end if ! iv_type == 1

! call v_calc early for sats for IC
! @RMM 8-6-08 - added to make SS vel consistent w/ var sat/transient
! call v_calc for steady state velocity
if (iv_type == 2) then
!t_prev = tnext
!read(16,*) tnext

! read in kx, ky, kz
call pf_read(hkx(:,:,:),kxfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
call pf_read(hky(:,:,:),kyfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
call pf_read(hkz(:,:,:),kzfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
!start sat test here	
	if(saturated == 1) then			
	vga=1.0
	 vgn=1.0
	 sres=1.0
     sat=1.0	
	end if
	
	if(saturated == 0)then		
		call pf_read(vga(:,:,:),vgafile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
		call pf_read(vgn(:,:,:),vgnfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
		call pf_read(sres(:,:,:),sresfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
	end if

	if(saturated == 2)then		
	 vga=vga_const
	 vgn=vgn_const
	 sres=sres_sat_const
	end if

        if( saturated == 3 ) then
	 vga=vga_const
	 vgn=vgn_const
         sat = sres_sat_const
        end if

hkx = hkx + 1E-15
hky = hky + 1E-15
hkz = hkz + 1E-15

if ( saturated .ne. 3 ) then
   call v_calc(v,hkx,hky,hkz,vga,vgn,sres,headfile,porosity,delv,xtent,ytent,ztent,press,sat)
else
   call v_calc_const_sat(v,hkx,hky,hkz,vga,vgn,headfile,porosity,delv,xtent,ytent,ztent,press,sat)
endif

end if ! calc'd vel ?

!
! call vcalc for transient vel simulation
if (iv_type == 3) then
!t_prev = tnext
!read(16,*) tnext
read(17,'(a100)') headfile
! read in kx, ky, kz
call pf_read(hkx(:,:,:),kxfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
call pf_read(hky(:,:,:),kyfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
call pf_read(hkz(:,:,:),kzfile,xtent,ytent,ztent,delv(1),delv(2),delv(3)) 
!start sat test here	
	if(saturated == 1)then			
	 vga=1.0
	 vgn=1.0
	 sres=1.0
     sat=1.0
	end if
	
	if(saturated == 0)then	
		call pf_read(vga(:,:,:),vgafile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
		call pf_read(vgn(:,:,:),vgnfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
		call pf_read(sres(:,:,:),sresfile,xtent,ytent,ztent,delv(1),delv(2),delv(3))
	end if
	
	if(saturated == 2)then		
	 vga=vga_const
	 vgn=vgn_const
	 sres=sres_sat_const
	end if

        if( saturated == 3 ) then
	 vga=vga_const
	 vgn=vgn_const
         sat = sres_sat_const
        end if

hkx = hkx + 1E-15
hky = hky + 1E-15
hkz = hkz + 1E-15
if ( saturated .ne. 3 ) then
   call v_calc(v,hkx,hky,hkz,vga,vgn,sres,headfile,porosity,delv,xtent,ytent,ztent,press,sat)
else
   call v_calc_const_sat(v,hkx,hky,hkz,vga,vgn,headfile,porosity,delv,xtent,ytent,ztent,press,sat)
endif
end if ! calc'd vel ?

!
!    Set up Initial Condition
!

n = 1
np = 0

DO iic = 1, n_constituents
 READ(99,*) iwflag(iic)
 IF (iwflag(iic) == 1) THEN
  READ(99,*) wellx
  WRITE(666,*) wellx
  READ(99,*) welly
  WRITE(666,*) welly
  READ(99,*) wellzb
  READ(99,*) wellzt
  READ(99,*) wellq
  READ(99,*) wellnpart
  WRITE(666,*) wellnpart
  
!
!    Loop through domain
!    Hard wire IC just for SL traces
!    RMM 7-8-98

  
  ijunk = INT((wellzt-wellzb)/delv(3))
  iwjt = INT(wellzt/delv(3))
  iwjb = INT(wellzb/delv(3))
  

  DO i=0, 360,(360/wellnpart)
    DO k=iwjb,iwjt
      p(n,1) = 0.50D0*delv(1)*cos(DBLE(i)) + wellx
      p(n,2) = 0.50D0*delv(2)*sin(DBLE(i)) + welly
      p(n,3) = delv(3)*k
      p(n,5) = n
	  iP(n,1) = iiC
	  iP(n,2) = 1
      n = n +1
    END DO
  END DO

 
 END IF

! pulse-type IC

 IF (iwflag(iic) == 2) THEN
  write(666,*)
  write(666,*) 'Pulse-type IC (type 2) for Constituent ',iic,trim(confile(iic))
  C(iiC,:,:,:) = 0.d0
  READ(99,*) xlic
  READ(99,*) xuic
  xlic = xlic !+ del2v(1)
  xuic = xuic !+ del2v(1)
  WRITE(666,'("X-lower:",e12.4," X-upper:",e12.4)') xlic , xuic
  READ(99,*) ylic
  READ(99,*) yuic
  ylic = ylic !+ del2v(2)
  yuic = yuic !+ del2v(2)
  WRITE(666,'("Y-lower:",e12.4," Y-upper:",e12.4)') ylic , yuic
  READ(99,*) zlic
  READ(99,*) zuic
  zlic = zlic !+ del2v(3)
  zuic = zuic !+ del2v(3)
  WRITE(666,'("Z-lower:",e12.4," Z-upper:",e12.4)') zlic , zuic

  READ(99,*) coic
  coic = coic * 1000.0 ! mg/l to mg/M^3
  READ(99,*) icnpart
  WRITE(666,'(" Number of Particles after this IC:",i12)') icnpart
  READ(99,*) decayic
  READ(99,*) dectimic
  READ(99,*) idecsteps
  ivolume = (xuic-xlic)*(ylic-yuic)*(zlic-zuic)
  pulse_vol = (xuic-xlic)*(yuic-ylic)*(zuic-zlic)
!  WRITE(666,'(" Total Volume [cells]:",i12)') ivolume
  write(666,'(" Total Volume [L**3]:",e12.4)') pulse_vol

!        write(666,*) coic*ivolume/dble(icnpart)

! round icnpart to nearest 3
  !icpartdiv = int(float(icnpart) ** (1./3.))
  !icnpart = icpartdiv**3
  
  
  WRITE(666,'(" C0 for Pulse Input [ppm]: ",f7.3)')  coic
  WRITE(666,'(" Initial mass of each particle [g]: ",e12.3)') (coic*pulse_vol*porosity(1,1,1))/DBLE(icnpart)
  write(666,*) 'porosity', porosity(1,1,1)
  WRITE(666,'(" pulse location x1, x2 [m]:",f8.4,",",f8.4)') xlic,xuic
  WRITE(666,'(" pulse location y1, y2 [m]:",f8.4,",",f8.4)') ylic,yuic
  WRITE(666,'(" pulse location z1, z2 [m]:",f8.4,",",f8.4)') zlic, zuic
  
  
  print*, n
  print*, icnpart
  print*, icpartdiv
  
  dicnstep = dectimic/DBLE(icnpart)
  
kk = 1
  
  DO  kk = 1, icnpart
  !do i = 1, icpartdiv
  !do j = 1, icpartdiv
  !do k = 1, icpartdiv
  
    p(n,1) = (xuic-xlic)*(ran1(ir)) + xlic
    p(n,2) = (yuic-ylic)*(ran1(ir)) + ylic
    p(n,3) = (zuic-zlic)*(ran1(ir)) + zlic

    !p(n,1) = (xuic-xlic)*(float(i)/float(icpartdiv)) + xlic
    !p(n,2) = (yuic-ylic)*(float(j)/float(icpartdiv)) + ylic
    !p(n,3) = (zuic-zlic)*(float(k)/float(icpartdiv)) + zlic

    !p(n,1) = (xuic-xlic)*(ran(ir)) + xlic
    !p(n,2) = (yuic-ylic)*(ran(ir)) + ylic
    !p(n,3) = (zuic-zlic)*(ran(ir)) + zlic
    
    ploc(1) = IDINT(p(n,1)/delv(1)) + 1
	ploc(2) = IDINT(p(n,2)/delv(2)) + 1
	ploc(3) = IDINT(p(n,3)/delv(3)) + 1

    !p(n,4) =  coic*ivolume/DBLE(icnpart)
!	p(n,4) =  coic*pulse_vol/(DBLE(icnpart)
	p(n,4) = (coic*pulse_vol*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3)))/DBLE(icnpart)

! print*, P(n,1), P(n,2), P(n,3), P(n,4) ,n

    !p(n,5) =  kk
    p(n,5) =  n
    p(n,7) = 0.0D0    !dble(n-1)*dicnstep
    iP(n,1) = iiC  
	iP(n,2) = 1  
    IF (concprint >= 1) THEN
      inbounds = 1
      DO  l = 1, numax
        ploc(l) = IDINT(p(n,l)/delv(l)) + 1
!        print*, ploc(l),l,n
! change of ploc, test
!		ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
!        IF((ploc(l) <= 1).OR.(ploc(l) >= domax(l) )) THEN
        IF((ploc(l) <= 0).OR.(ploc(l) > domax(l) )) THEN
          inbounds = 0
        END IF
      END DO  

      
      IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
        c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
            SNGL( (dble(iP(n,2))*(p(n,4)))/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )
            !SNGL( (dble(iP(n,2))*(p(n,4)))/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))) )
		!print*, c(ip(n,1),ploc(1),ploc(2),ploc(3)), ip(n,1),ploc(1),ploc(2),ploc(3)
		!print*, cellv , porosity(ploc(1),ploc(2),ploc(3)) , Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
		!print*, (dble(iP(n,2))*(p(n,4)))
      END IF
    END IF  ! printing concentrations
    n = n +1
    
   END DO !!kk
!end do !i
!end do !j
!end do !k
  
 END IF

! first Indicator IC

 IF (iwflag(iic) == 3) THEN
  
  READ(99,*) ind_ic_file
  OPEN(123,FILE=trim(ind_ic_file),STATUS='old')
  READ(99,*) ind_ic_catagory
  READ(123,*) icat, jcat, kcat
  ic_cat_count = 0
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        READ(123,*) ic_cat_num(iic, ii,jj,kk)
        IF (ic_cat_num(iic, ii,jj,kk) == ind_ic_catagory)  &
			 ic_cat_count = ic_cat_count + 1
      END DO
    END DO
  END DO
  CLOSE (123)
  PRINT*, 'number ic cells:',ic_cat_count
  
  READ(99,*) coic
  READ(99,*) icnpart
  WRITE(666,*) icnpart
  
  WRITE(666,'(" C0 for Ind File Input [ppm]: ",f7.3)')  coic
! write(666,'(" Initial mass of each particle [g]: ",f7.3)')
! 1 coic*(porosity(ploc(1),ploc(2),ploc(3))*
!    2 cellv*dble(ic_cat_count))/dble(icnpart)
  
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        IF (ic_cat_num(iic, ii,jj,kk) == ind_ic_catagory) THEN
! assign particles to active locations
          DO  kic = 1, (icnpart/ic_cat_count)
            p(n,1) = delv(1)*(ran1(ir)) + DBLE(ii-1)*delv(1) 
            p(n,2) = delv(2)*(ran1(ir)) + DBLE(jj-1)*delv(2) 
            p(n,3) = delv(3)*(ran1(ir)) + DBLE(kk-1)*delv(3) 
            p(n,4) =  coic*(sat(ii,jj,kk)*porosity(ii,jj,kk)*   &
				cellv*DBLE(ic_cat_count))/DBLE(icnpart)
            p(n,5) =  n
            p(n,7) = 0.0D0
            iP(n,1) = iiC
            iP(n,2) = 1
            IF (concprint >= 1) THEN
              inbounds = 1
              DO  l = 1, numax
                ploc(l) = IDINT(p(n,l)/delc(l)) + 1
! change in ploc test
!				ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
                
!                IF((ploc(l) <= 2).OR.(ploc(l) >= domax(l))) THEN
                IF((ploc(l) <= 0).OR.(ploc(l) > domax(l))) THEN
                  inbounds = 0
                END IF
              END DO
              
              IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
                c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
                    SNGL( dble(iP(n,2))*(p(n,4))/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )

              END IF
            END IF  ! printing concentrations
            n = n +1
          END DO
        END IF
      END DO  !ii
    END DO  !jj
  END DO  !kk
  
 END IF

! second Indicator IC

 IF (iwflag(iic) == 4 .OR. iwflag(iic) == 5 ) THEN
  PRINT*, ' ind mult'
  WRITE(666,*) ' Mulitple Indicator IC'
  WRITE(666,*)
  
  READ(99,*) ind_ic_file
  WRITE(666,*) ' reading from file:',trim(ind_ic_file)
  OPEN (124,FILE=trim(ind_ic_file),STATUS='old')

! skip comment file


  READ(124,*)
  READ(124,*)
  READ(124,*)
  READ(124,*)
  
  READ (124,*) ind_ic_file
  READ (124,*) n_ic_cats(iic)
  READ (124,*) ic_cats(iic, 1:n_ic_cats(iic))
  READ (124,*) ic_part_dens(iic, 1:n_ic_cats(iic))
  READ (124,*) n_ic_timesteps
  
  WRITE(666,*) ' Ind file:',trim(ind_ic_file)
  WRITE(666,*) ' Num Indicators:',n_ic_cats(iic)
  WRITE(666,*) ' Indicators:',ic_cats(iic, 1:n_ic_cats(iic))
  WRITE(666,*) ' Particle Density:',ic_part_dens(iic, 1:n_ic_cats(iic))
  WRITE(666,*) ' Num Timesteps:',n_ic_timesteps
  
  READ (124,*)
  DO ii = 1, n_ic_timesteps
    READ (124,*) ts1, ic_time_begin(iic,ii),ic_time_end(iic,ii),   &
                 ic_mass_or_conc(iic, ii,1:n_ic_cats(iic))
    WRITE(666,*) 'timestep:',ii
    WRITE(666,*)  'time begin, end:',ic_time_begin(iic, ii),ic_time_end(iic, ii)
    WRITE(666,*) ic_mass_or_conc(iic, ii,1:n_ic_cats(iic))
  END DO
  CLOSE (124)
  PRINT*,' read ind file'
  !OPEN(123,FILE=trim(ind_ic_file),STATUS='old',readonly)
  OPEN(123,FILE=trim(ind_ic_file),form='unformatted', access='stream',STATUS='old')

! read(99,*) ind_ic_catagory

  !READ(123,*) icat, jcat, kcat
  READ(123) icat, jcat, kcat
  ic_cat_count = 0
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        READ(123) ic_cat_num(iic, ii,jj,kk)
      END DO
    END DO
  END DO
  CLOSE (123)
 
  PRINT*,' loop thru ic'
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        DO jjj = 1, n_ic_timesteps
          DO iii = 1, n_ic_cats(iic)
            IF (ic_cat_num(iic, ii,jj,kk) == ic_cats(iic,iii)) THEN

! assign particles to active locations
! check to see if mass of ic > 0
             IF (ic_mass_or_conc(iic,jjj,iii) > 0.D0) THEN
              DO  kic = 1, ic_part_dens(iic,iii)
	            p(n,1) = delv(1)*(ran1(ir)) + DBLE(ii-1)*delv(1)
                p(n,2) = delv(2)*(ran1(ir)) + DBLE(jj-1)*delv(2)
                p(n,3) = delv(3)*(ran1(ir)) + DBLE(kk-1)*delv(3)
! hardwire particles to sit one cell below trench indicator for trench infil test
                !p(n,3) = delv(3)*(ran1(ir)) + DBLE(kk-2)*delv(3)
! hardwire particles to sit at BOTTOM of cell- only valid for trench/playa runs
                !p(n,3) =  DBLE(kk-1)*delv(3)
                p(n,4) =  ic_mass_or_conc(iic,jjj,iii)/ic_part_dens(iic,iii)
                p(n,5) =  n
                p(n,7) = ic_time_begin(iic,jjj) + ran1(ir)*(ic_time_end(iic,jjj)-ic_time_begin(iic,jjj))
                ip(n,1) = iiC
                ip(n,2) = 1
                IF (concprint >= 1) THEN
                  inbounds = 1
                  DO  l = 1, numax
                    ploc(l) = IDINT(p(n,l)/delc(l)) + 1
! change in ploc test
!					ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
!                    IF((ploc(l) <= 1).OR.(ploc(l) >= domax(l) ) ) THEN
                    IF((ploc(l) <= 0).OR.(ploc(l) > domax(l) ) ) THEN
                      inbounds = 0
                    END IF
                  END DO !l
				  IF (P(n,7) > 0.) inbounds = 0

                  IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
                    c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
                       SNGL( dble(iP(n,2))*(p(n,4))/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )
                  END IF   ! inbounds?
                END IF  ! printing concentrations?
                
                n = n +1
                IF (n > npmax) THEN
                  PRINT*, 'max num particles exceeded.  increase npmax parameter ',npmax
                  WRITE(666,*) 'Max Particles Exceeded!  Increase NPMAX:',npmax
                  STOP
                END IF ! max parts ?
               
              END DO  ! kic
             END IF  ! mass >0 ? 
            END IF  ! at a ic node
            
          END DO !iii
        END DO !jjj
      END DO  !ii
    END DO  !jj
  END DO  !kk
  
 END IF ! ic = 4 .or. ic =5

 IF ( iwflag(iic) == 5 ) THEN
  READ(99,*) ind_ic_file
  OPEN(123,FILE=trim(ind_ic_file),form='unformatted', access='stream',STATUS='old')
  READ(99,*) ind_ic_catagory
  READ(123) icat, jcat, kcat
  ic_cat_count = 0
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        READ(123) ic_cat_num_bg(iic,ii,jj,kk)
        IF (ic_cat_num_bg(iic,ii,jj,kk) == ind_ic_catagory)  &
			 ic_cat_count = ic_cat_count + 1
      END DO
    END DO
  END DO
  CLOSE (123)
  PRINT*, 'number ic cells:',ic_cat_count
  
  READ(99,*) coic
  READ(99,*) icnpart
  WRITE(666,*) icnpart
  
  WRITE(666,'(" C0 for Ind File Input [ppm]: ",f15.9)')  coic
! write(666,'(" Initial mass of each particle [g]: ",f7.3)')
! 1 coic*(porosity(ploc(1),ploc(2),ploc(3))*
!    2 cellv*dble(ic_cat_count))/dble(icnpart)
  
 IF( icnpart .GT. 0 .AND. coic .GT. 0.0 ) THEN
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        IF (ic_cat_num_bg(iic,ii,jj,kk) == ind_ic_catagory) THEN
! assign particles to active locations
          DO  kic = 1, (icnpart/ic_cat_count)
            p(n,1) = delv(1)*(ran1(ir)) + DBLE(ii-1)*delv(1) 
            p(n,2) = delv(2)*(ran1(ir)) + DBLE(jj-1)*delv(2) 
            p(n,3) = delv(3)*(ran1(ir)) + DBLE(kk-1)*delv(3) 
            p(n,4) =  coic*(sat(ii,jj,kk)*porosity(ii,jj,kk)*   &
				cellv*DBLE(ic_cat_count))/DBLE(icnpart)
            p(n,5) =  n
            p(n,7) = 0.0D0
            iP(n,1) = iiC
            iP(n,2) = 1
            IF (concprint >= 1) THEN
              inbounds = 1
              DO  l = 1, numax
                ploc(l) = IDINT(p(n,l)/delc(l)) + 1
! change in ploc test
!				ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
                
                IF((ploc(l) <= 2).OR.(ploc(l) >= domax(l))) THEN
                  inbounds = 0
                END IF
              END DO
              
              IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
                c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
                    SNGL( dble(iP(n,2))*(p(n,4))/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )

              END IF
            END IF  ! printing concentrations
            n = n +1
          END DO
        END IF
      END DO  !ii
    END DO  !jj
  END DO  !kk
 END IF
 END IF ! ic =5

! third Indicator IC, given concentratio at the boundary

 IF (iwflag(iic) == 6) THEN
  PRINT*, ' third ind mult - boundary conc.'
  WRITE(666,*) ' 3rd Mulitple Indicator IC - boundary conc'
  WRITE(666,*)
  
  READ(99,*) ind_ic_file
  WRITE(666,*) ' reading from file:',trim(ind_ic_file)
  OPEN (124,FILE=trim(ind_ic_file),STATUS='old')

! skip comment files

  READ(124,*)
  READ(124,*)
  READ(124,*)
  READ(124,*)
  
  READ (124,*) ind_ic_file
  READ (124,*) n_ic_cats(iic)
  READ (124,*) ic_cats(iic, 1:n_ic_cats(iic))
  READ (124,*) ic_part_dens(iic, 1:n_ic_cats(iic))
  READ (124,*) n_ic_timesteps
  
  WRITE(666,*) ' Ind file:',trim(ind_ic_file)
  WRITE(666,*) ' Num Indicators:',n_ic_cats(iic)
  WRITE(666,*) ' Indicators:',ic_cats(iic, 1:n_ic_cats(iic))
  WRITE(666,*) ' Particle Density:',ic_part_dens(iic, 1:n_ic_cats(iic))
  WRITE(666,*) ' Num Timesteps:',n_ic_timesteps
  
  READ (124,*)
  DO ii = 1, n_ic_timesteps
    READ (124,*) ts1, ic_time_begin(iic, ii),ic_time_end(iic, ii),ic_mass_or_conc(iic, ii,1:n_ic_cats(iic))
    WRITE(666,*) 'timestep:',ii
    WRITE(666,*)  'time begin, end:',ic_time_begin(iic, ii),ic_time_end(iic, ii)
    WRITE(666,*) ic_mass_or_conc(iic, ii,1:n_ic_cats(iic))
  END DO
  CLOSE (124)
  PRINT*,' read ind file'
  !OPEN(123,FILE=trim(ind_ic_file),STATUS='old',readonly)
  OPEN(123,FILE=trim(ind_ic_file),form='unformatted', access='stream',STATUS='old')

  READ(123) icat, jcat, kcat
  ic_cat_count = 0
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        READ(123) ic_cat_num(iic, ii,jj,kk)
      END DO
    END DO
  END DO
  CLOSE (123)

  READ(99,*) ind_ic_file
  OPEN(123,FILE=trim(ind_ic_file),form='unformatted', access='stream',STATUS='old')
  READ(99,*) ind_ic_catagory
  READ(123) icat, jcat, kcat
  ic_cat_count = 0
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        READ(123) ic_cat_num_bg(iic,ii,jj,kk)
        IF (ic_cat_num_bg(iic,ii,jj,kk) == ind_ic_catagory)  &
			 ic_cat_count = ic_cat_count + 1
      END DO
    END DO
  END DO
  CLOSE (123)
  PRINT*, 'number ic cells:',ic_cat_count
  
  READ(99,*) coic
  READ(99,*) icnpart
  WRITE(666,*) icnpart
  
  WRITE(666,'(" C0 for Ind File Input [ppm]: ",f15.9)')  coic
! write(666,'(" Initial mass of each particle [g]: ",f7.3)')
! 1 coic*(porosity(ploc(1),ploc(2),ploc(3))*
!    2 cellv*dble(ic_cat_count))/dble(icnpart)
  
 IF( icnpart .GT. 0 .AND. coic .GT. 0.0 ) THEN
  DO kk = 1, kcat
    DO jj = 1, jcat
      DO ii = 1, icat
        IF (ic_cat_num_bg(iic,ii,jj,kk) == ind_ic_catagory) THEN
! assign particles to active locations
          DO  kic = 1, (icnpart/ic_cat_count)
            p(n,1) = delv(1)*(ran1(ir)) + DBLE(ii-1)*delv(1) 
            p(n,2) = delv(2)*(ran1(ir)) + DBLE(jj-1)*delv(2) 
            p(n,3) = delv(3)*(ran1(ir)) + DBLE(kk-1)*delv(3) 
            p(n,4) =  coic*(sat(ii,jj,kk)*porosity(ii,jj,kk)*   &
				cellv*DBLE(ic_cat_count))/DBLE(icnpart)
            p(n,5) =  n
            p(n,7) = 0.0D0
            iP(n,1) = iiC
            iP(n,2) = 1
            IF (concprint >= 1) THEN
              inbounds = 1
              DO  l = 1, numax
                ploc(l) = IDINT(p(n,l)/delc(l)) + 1
! change in ploc test
!				ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
                
                IF((ploc(l) <= 2).OR.(ploc(l) >= domax(l))) THEN
                  inbounds = 0
                END IF
              END DO
              
              IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
                c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
                    SNGL( dble(iP(n,2))*(p(n,4))/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )

              END IF
            END IF  ! printing concentrations
            n = n +1
          END DO
        END IF
      END DO  !ii
    END DO  !jj
  END DO  !kk
 END IF

 END IF ! ic = 6

 END DO  ! DO iiC = 1, n_constituents

np = n-1

print*, 'de allocating'
!deallocate(ic_cat_num)

partsplit = 0
ts = 0
cmin = 0.0D0
tnext = 0.
print*, 'v call'

!@RMM moved this up to before IC so that we calc Sat and can use for IC and to keep constistancy
!if (iv_type == 2) call v_calc(v,kxfile,kyfile,kzfile,headfile,porosity,delv,xtent,ytent,ztent,press,sat)

! Check if splitting particles?

READ(99,*) partsplit
PRINT*, ' partsplit',partsplit
IF (partsplit == 1) READ(99,*) cmin
READ(99,*) tempavg
PRINT*,' tempav',tempavg
IF(tempavg == 1) WRITE (666,*) ' Temporal Averaging '

PRINT*,'ts'
READ(99,*) ts

READ(99,*) Xup_Ref
READ(99,*) bnd_Xup
READ(99,*) Xdown_Ref
READ(99,*) bnd_Xdown
READ(99,*) Yup_Ref
READ(99,*) bnd_Yup
READ(99,*) Ydown_Ref
READ(99,*) bnd_Ydown
READ(99,*) Zup_Ref
READ(99,*) bnd_Zup
READ(99,*) Zdown_Ref
READ(99,*) bnd_Zdown

if (iv_type == 1) then
WRITE(666,*)
WRITE(666,*) ' number of timesteps to skip: ',ts
nt = nt - ts
CLOSE(99)

! skip ts timesteps
do i= 1, ts
  CALL compact_vread(v,domax(1),domax(2),domax(3),  &
      timenext)

end do ! i, timestep skips
end if ! vtype

!
!    Open concentration file and print concentrations
!    (if we're asked to)
!

PRINT*, 'flag4'

IF (concprint == 1) THEN
 DO i=1, n_constituents
  dotit = '.00000.'
  filename= trim(confile(i))//dotit//confile2
  CALL cbin_write(c(i,:,:,:),filename,xtent, ytent,ztent,delv(1),delv(2),delv(3))
 END DO
  else if (concprint == 2 ) THEN

  CALL vtk_write(time, c(:,:,:,:),confile(:),xtent, ytent,ztent,delv(1),delv(2),delv(3),0,n_constituents,vtk_file)
  
END IF    ! print concentrations?
 
! set time column of mass to zero for initi condition

mass(1,1) = 0.d0

DO  n = 1, np
  mass(iP(n,1)+1,1) = mass(iP(n,1)+1,1) + P(n,4)
  rxx = rxx + p(n,1)
  ryy = ryy + p(n,2)
  rzz = rzz + p(n,3)
  
  61  FORMAT(4(e15.8,','),e15.6)
  62  FORMAT(3(f6.1,','),3(f9.4,','),f9.4)
END DO
rxx = rxx/np
ryy = ryy/np
rzz = rzz/np
DO  n = 1, np
  sxx = sxx + (p(n,1)-rxx)**2
  sxy = sxy + (p(n,1)-rxx)*(p(n,2)-ryy)
  syy = syy + (p(n,2)-ryy)**2
  szz = szz + (p(n,3)-rzz)**2
  
END DO
sxx = sxx/np
sxy = sxy/np
syy = syy/np
szz = szz/np

IF (momprint == 1) THEN
  WRITE (wellfile+10,*) 'Moment anlysis'
  WRITE (wellfile+10,*) ' time, Mass, Rx,  Ry, Rz,Sx, Sy, Sz'
  WRITE(wellfile+10,188) time,np,rxx,ryy,rzz,sxx,syy,szz
END IF    ! print moments?

PRINT*, 'flag 5'
900   FORMAT (i4,i4,e15.5,e15.5)
910   FORMAT (e10.5)

! set up well recycle arrays if we are doing this
if (recycle_well == 1) then
print*, 'recycling well'
n_rw = 2
rw_ind(1) = 28
rw_ind(2) = 41
rw_flux(1) = 0.56 
rw_flux(2) = 1.-0.56 
n_rw_ind = 0
do i = 1, xtent
 do j = 1, ytent
  do k = 1, ztent
   do ii = 1, n_rw
   if (ic_cat_num(1,i,j,k)==rw_ind(ii)) then
    n_rw_ind(ii) = n_rw_ind(ii) + 1
    well_r(ii,n_rw_ind(ii),1) = i
    well_r(ii,n_rw_ind(ii),2) = j
    well_r(ii,n_rw_ind(ii),3) = k
print*, well_r(ii,n_rw_ind(ii),1),well_r(ii,n_rw_ind(ii),2),well_r(ii,n_rw_ind(ii),3)
print*, i,j,k,ii,n_rw,n_rw_ind(ii) 
    end if
   end do
  end do
 end do
end do


end if

CALL allocateParticles_Memory(xtent, ytent, ztent, n_constituents, npmax)

!
! Big Time loop.  loop until time is up, reporting concentrations
! at the end of each minor loop.
!

PRINT*, 'flag 6'
DO  it = 1, nt

  PRINT*, 'got to big loop'

!  
! Read in Velocities for that time step from NUFT output
!
 
 
if (iv_type == 1) then
  CALL compact_vread(v,domax(1),domax(2),domax(3),  &
      timenext)
  
! convert seconds to days
  t_prev = tnext
  tnext = timenext !/86400.d0
end if
if (iv_type == 2) then
	timenext = tnext + t_big_step
  t_prev = tnext
  tnext = tnext + t_big_step
  
end if
! don't call here the first timestep, we call earlier to get sats for IC
if (iv_type == 3) then
t_prev = tnext
read(16,*) tnext
if (it > 1) then
read(17,'(a100)') headfile
call v_calc(v,hkx,hky,hkz,vga,vgn,sres,headfile,porosity,delv,xtent,ytent,ztent,press,sat)
end if ! first timestep?
end if ! calc'd vel ?

  if(massive_debug == 1) then
  print*
  print*,' reading in next velocity field time step'
  print*,' timenext',timenext
  end if

  print*,' tnext', tnext
  print*,' t_prev', t_prev

  nstepav = 0.
  WRITE(666,*)
  write(666,*)
  WRITE(666,'("*** Loop Number:",i12," ***")') it
  WRITE(666,*)
  WRITE(666,'(" v_time [s]:",e12.4)') timenext
  WRITE(666,'("Time [d]:",e12.4,", [y]:",e12.4)' ) tnext,tnext/365.24
  
  movedir = 99

DO iic = 1, n_constituents
 IF (iwflag(iic) == 6) THEN
  !
  ! generate particles for the 1st time step
  !
  np = np + 1
  CALL gen_part( icat, jcat, kcat, ic_cat_num(iic,:,:,:) , n_ic_cats(iic), &
            ic_cats(iic,:),  ic_mass_or_conc(iic, :,:), ic_part_dens(iic,:), &
             p, ip, delv, v, np, it, iic, npmax, C,                        &
             ic_time_begin( iic, it ), ic_time_end( iic, it ), porosity, sat )
  np = np - 1
 END IF ! ic = 6

 END DO  ! DO iiC = 1, n_constituents

WRITE(666,*)
WRITE(666,*) ' Number of Particles in step: ', it, ' is ', np
WRITE(666,*)

!
! Clear out old concentrations
!
   c = 0.0D0
  
 ! print*, C(1,1,1,1)

!
!reverse velocities if needed
! 
 V = V*vmult 
! mp = np
!
! Big Particle loop.  This guy is where all the particles get moved.
! This is also where SLIM_FAST can be most easily parallelized.
!
  DO  n = 1, np

! continue statement for well recycling
1999 continue

origloc( 1 ) = p(n,1)
origloc( 2 ) = p(n,2)
origloc( 3 ) = p(n,3)
    IF (partprint == 1) THEN
      IF ((p(n,1) > 0.).AND.(p(n,2) > 0.).AND.(p(n,3) > 0.))  &
      pdist = 0.0d0
      do ii = 1, 3
      pdist = pdist +  (p(n,ii)-lastprint(n,ii))**2 
      end do
      pdist = dsqrt(pdist)
      if (pdist >= minpdist) then
      !WRITE(13,61) p(n,1), p(n,2), p(n,3), p(n,5), p(n,7)
      lastprint(n,1) = p(n,1)
      lastprint(n,2) = p(n,2)
      lastprint(n,3) = p(n,3)
      end if ! min dist/ last print
    END IF
    
! if (partprint.eq.2) then
! write(13,62) P(n,1), P(n,2), P(n,3),  P(n,7)
! 1 ,P(n,8),P(n,9),P(n,10)
! end if
    
    tloop = 0
    loop2 = 0
    done = 0
    trgp = 0.
    passes = 0
    DO WHILE (done == 0)
      tad(1) = 0.0D0
      tad(2) = 0.0D0
      tad(3) = 0.0D0
      tad(4) = 0.0D0
      
      
!
! Find particle cell location
!
      
      DO  l = 1, numax
! uncomment next lines for massive debugging output
!		print*, P(n,l), delv(l), n l
!        ploc(l) = INT((p(n,l)-del2v(l))/delv(l))
! ploc test
		ploc(l) = IDINT(p(n,l)/delv(l)) + 1
!        f(l) = ((p(n,l)-del2v(l))/delv(l)) - DBLE(ploc(l))
! f test
        f(l) = (p(n,l)/delv(l)) - float(ploc(l)-1)

!		write(667,*) l,ploc(l)
        
        stuck= 1.0D0-f(l)
! if (stuck.le.epsi) then
!  f(l) = 1.0d0
!  ploc(l) = ploc(l) + 1
! end if
        
        stuck = f(l)
!        if (stuck.le.epsi) then
!  f(l) = 0.0d0
!  ploc(l) = ploc(l) -1
! end if
        
        
!        stuck = epsi*delv(l)+ del2v(l)
!        locstuck(l) = INT((p(n,l)+stuck)/delv(l))
!        IF (locstuck(l) /= ploc(l)) THEN
!          ploc(l) = locstuck(l)
!          IF (stuck > 0) THEN
!            f(l) = 0.0D0
!          ELSE
!            f(l) = 1.0D0
!          END IF
!        END IF
!        GO TO 992
        
!        IF ((f(l) > 0.0D0).AND.(f(l) < 1.0D0)) THEN
!          stuck = f(l) - epsi
!          IF(stuck < 0.d0) THEN
!            ploc(l) = ploc(l) - 1
!            f(l) = 1.0D0
!            GO TO 991
!          END IF
          
!          stuck = f(l) + epsi
!          IF (stuck > 1.0D0) THEN
!            ploc(l) = ploc(l) + 1
!            f(l) = 0.0D0
!          END IF
!          991 CONTINUE
!        END IF
!        992 CONTINUE
        
!
! Check to make sure we're still in bounds
!
! in space
!
! reset bddy condition
 iP(n,4) = 0
!        IF (ploc(l) >= (domax(l)-1)) THEN
!         IF (ploc(l) >= (domax(l))) THEN
         IF (ploc(l) > (domax(l))) THEN
		IF (Boundary_cond(l,2) == 0) THEN
          done = 1
          IF(p(n,l) >= 0.) p(n,l) = -p(n,l)
		  if(massive_debug == 1) then
!		  print*,' outside domain'
!		  print*,'+',n,ip(n,1),l,p(n,1),p(n,2),p(n,3)
		  end if

!		  print*,' outside domain'
!		  print*,'+',n,ip(n,1),l,p(n,1),p(n,2),p(n,3)

          GO TO 151
		  ELSE
		  iP(n,4) = l

		  END IF
        END IF

!        IF (ploc(l) <= 2) THEN
        IF (ploc(l) <= 0) THEN
		 IF (Boundary_cond(l,1) == 0) THEN
          IF(p(n,l) >= 0.) p(n,l) = -p(n,l)
		  if(massive_debug == 1) then
!		  print*,' outside domain'
!         print*,'-',n,ip(n,1),l,p(n,1),p(n,2),p(n,3)
		  end if

!		  print*,' outside domain'
!          print*,'-',n,ip(n,1),l,p(n,1),p(n,2),p(n,3)

          done = 1
          GO TO 151
		  ELSE
		  iP(n,4) = l
		 END IF
		END IF

       IF ( ploc( l ) >= domax( l ) ) THEN
         po( l, l ) = 0 
       ELSE
         po( l, l ) = 1 
       ENDIF

     END DO

! and in time
 
     IF (P(n,7) > tnext) THEN
          done = 1
		  if(massive_debug == 1) then
		  print*,' skipping due to time'
		  print*,' time=',time
		  end if
          GO TO 151
     END IF
 !!
 !@RMM Remove/comment below 
 !!
  ! check for inactive zone if we are not in 12 (river) or 1 (gneral subsurface) dump
! if((ic_cat_num(ploc(1),ploc(2),ploc(3)) /= 12).or.     &
!    (ic_cat_num(ploc(1),ploc(2),ploc(3)) /= 1)) then
 !if (ic_cat_num(ploc(1),ploc(2),ploc(3)) == 0) then
!done =1
!if(iprP(n,2) == 0 ) then
!iprP(n,2) = 1
!@RMM particle print fix
!if (partprint >= 1) write(313,61) P(n,1), P(n,2), P(n,3),P(n,5), P(n,7)
!end if
!go to 151
!print*, 'in inactive zone/air done',ploc(1),ploc(2),ploc(3),n
!end if
! if(sat(ploc(1),ploc(2),ploc(3)) < 0.95d0) then
!done =1
!if (iprP(n,1) == 0 )  then 
!iprP(n,1) = 1
!write(213,61) P(n,1), P(n,2), P(n,3),P(n,5), P(n,7)
!end if
!go to 151
!print*, 'low sat done',ploc(1),ploc(2),ploc(3),n
!end if
! write(*,*) 'get here?'

      passes = passes +1
      loc3p = ploc(3)
      
!
! Check to change coord frame if neg vel and cell bddy
!

      DO  l=1, numax
        IF ((f(l) < epsi).AND.(v(l,ploc(1),ploc(2),ploc(3)) < 0.0D0)) THEN
		  if(massive_debug == 1) then
		  write(667,*), ' unsticking f=',f(l)
		  write(667,*), ' l=',l,' i,j,k=',ploc(1),ploc(2),ploc(3)
		  end if
          ploc(l) = ploc(l) - 1
          f(l) = 1.0D0

        ELSE
! old way rmm 9-21		 IF ((ABS(1.D0-f(l)) < epsi).AND.(v(l,ploc(1),ploc(2),ploc(3)) > 0.0D0)) THEN
		 IF ((DABS(1.D0-f(l)) < epsi).AND.(v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) > 0.0D0)) THEN
		  if(massive_debug == 1) then
          write(667,*), ' unsticking f=',f(l)
! old way		  write(667,*), ' l=',l,' i,j,k=',ploc(1),ploc(2),ploc(3)
			write(667,*), ' l=',l,' i,j,k=',ploc(1),ploc(2),ploc(3)
		  end if
		    f(l) = 0.0D0
			ploc(l) = ploc(l) + 1
         END IF
		END IF
        
                if ( ploc( l ) > domax( l ) ) ploc( l ) = domax( l )
                if ( ploc( l ) < 1 ) ploc( l ) = 1
      END DO

! print*,'checking geotype'
! print*,ploc(1),ploc(2),ploc(3)
! print*,ic_cat_num(ploc(1),ploc(2),ploc(3))

      p(n,9) = ic_cat_num(ip(n, 1), ploc(1),ploc(2),ploc(3))

! if we are sorbed/stuck jump to calculating prob of unsticking
if (iP(n,2) == 0) go to 919

!
! Check for wells
!
      IF (nw > 0) THEN

        DO  wc=1, nw
          IF(ploc(1) == welloc(wc,1)) THEN
		    IF(ploc(2) == welloc(wc,2)) THEN

              IF((ploc(3) <= welloc(wc,3)).AND.(ploc(3) >= welloc(wc,4))) THEN
                
                IF (well(wc,5) /= 0.d0) THEN
                  
                  DO  l = 1, numax
                    IF (f(l) <= (.5)) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
                      vp(l) = v(l,ploc(1),ploc(2),ploc(3))/Rstar 
                    ELSE
                      vp(l) = v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3))/Rtard(ip(n,1),ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3))
                    END IF
                  END DO
                  
                  dtr = DSQRT(((f(1)-0.5))**2 + ((f(2)-0.5))**2)
                  IF (dtr == 0.0D0) THEN
                    WRITE(667,*) 'dtr=0!'
                    GO TO 151
                  END IF
                  vpw = (ABS(.5-f(1))/dtr)*vp(1) + (ABS(.5-f(2))/dtr)*vp(2)
                  IF (vpw == 0.0D0) THEN
                    WRITE(667,*) 'vpw=0!'
                    GO TO 151
                  END IF
                  
                  dtw = SQRT((delv(1)*(f(1)-0.5))**2 + (delv(2)*(f(2)-0.5))**2)
                  tad(4) = dtw/vpw
                  
                  
                  tad(4) = dtw/vpw
                  
                 if (recycle_well == 1) then
                 !print*, 'recycling AT well', n
                  irP(n) = irP(n) + 1
                  ! choose 1/2 trench/playa
                 z(1) = ran1(ir)
                      ii = 1
                    if(z(1) > rw_flux(1) ) ii = 2 
                    ! pick a new source block at random from either trench/playa
                    jj = dint(ran1(ir)*n_rw_ind(ii))
                    ! assign x,y to be random w/in source block, z at bottom
                    p(n,1) = delv(1)*(ran1(ir)) + DBLE(well_r(ii,jj,1)-1)*delv(1)
                    p(n,2) = delv(2)*(ran1(ir)) + DBLE(well_r(ii,jj,2)-1)*delv(2)
                    p(n,3) =  DBLE(well_r(ii,jj,3)-1)*delv(3)
                    !print*, ii, jj
                    !print*, P(n,1),P(n,2),P(n,3)
                    !print*, (well_r(ii,jj,1)-1),(well_r(ii,jj,2)-1),(well_r(ii,jj,3)-1)
              
                    !P(n,7) = P(n,7) + tad(4)
                    !print*, P(n,7), iP(n,6), n
                   IF (partprint == 1) THEN 
                     IF ((p(n,1) > 0.).AND.(p(n,2) > 0.).AND.(p(n,3) > 0.)) then 
                      IF ((p(n,1) > 0.).AND.(p(n,2) > 0.).AND.(p(n,3) > 0.)) then 
                       pdist = 0.0d0
                       do ii = 1, 3
                         pdist = pdist +  (p(n,ii)-lastprint(n,ii))**2 
                       end do
                       pdist = dsqrt(pdist)
                       if (pdist >= minpdist) then
                       !WRITE(13,61) p(n,1), p(n,2), p(n,3), p(n,5), p(n,7)
                       lastprint(n,1) = p(n,1)
                       lastprint(n,2) = p(n,2)
                       lastprint(n,3) = p(n,3)
                      end if  ! print
                      end if ! in bounds
                      end if ! in bounds
                        END IF 

                     else  ! not recycling mass from wells, just removing it from the domain
                  p(n,1) = -999
                  p(n,2) = -999
                  p(n,3) = -999
                  p(n,6) = wc
                  p(n,7) = p(n,7) +  tad(4)
                  IF (welltavg(wc) < epsi) THEN
                    WRITE(667,*) 'wt avg lt 0 !'
                    GO TO 151
                end if ! well recycling?
                  END IF
                  itemp = ABS(INT(p(n,7)/welltavg(wc))) + 1
                  IF (itemp <= 0) itemp = 1
                  IF(itemp >= welltnumb) itemp = welltnumb
                  btc(wc,itemp,ip(n,1)) = btc(wc,itemp,ip(n,1)) + p(n,4)  !*dble(iP(n,2))
                  if (recycle_well == 1) go to 1999
                    GO TO 151
                  
                ELSE    ! mon well?

					if (ipwell(n,wc) == 1) then
                  itemp = ABS(IDINT(p(n,7)/welltavg(wc))) + 1
                  IF (itemp <= 0) itemp = 1
                  IF(itemp >= welltnumb) itemp = welltnumb
                  btc(wc,itemp,ip(n,1)) = btc(wc,itemp,ip(n,1))  &
					+ p(n,4)/(cellv * dble(welloc(wc,3)-welloc(wc,4)+1) *sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3)))
!	print*, dble(welloc(wc,3)-welloc(wc,4)+1)
						ipwell(n,wc) = - 1

!                  btc(wc,itemp,ip(n,1)) = btc(wc,itemp,ip(n,1))  &
!					+ p(n,4)/( cellv* (well(wc,3)-well(wc,4))* porosity(ploc(1),ploc(2),ploc(3)) )
!						ipwell(n,wc) = - 1

					end if ! seen this well before?

                END IF ! endif or monitoring well
                
                
              END IF   !endif for z-coord in wellscreen
            END IF !endif for y-coord=wellx
          END IF  !endif for x-coord=wellx
          
        END DO
      END IF   !endif for nw=0
      
                  !  print*,' r c ',P(n,1),P(n,2),P(n,3)
! check monitoring plane
!
      IF (nplane > 0) THEN
        DO  iplane=1, nplane
          planeloc = DINT(xplaneloc(iplane)/delv(planedir(iplane))) + 1
!		planeloc = INT((xplaneloc(iplane)+del2v(l))/delv(planedir(iplane)))
!			print*,ploc(planedir(iplane)),planeloc
          IF(ploc(planedir(iplane)) == planeloc) THEN
!		  print*,ploc(planedir(iplane)),planeloc
		  IF(iP(n,4+iplane) == 1) THEN

         if(massive_debug == 1) then
           print*, ' in plane'
           print*, 'planeloc=',planeloc
		    print*, 'xplaneloc=',xplaneloc(iplane)
		    print*, 'iplane=',iplane,' planedir=',planedir(iplane)
		    print*, ploc(1),ploc(2),ploc(3)
		 end if
            itemp = ABS(INT(p(n,7)/welltavg(1)))  + 1
!			print*, itemp,p(n,7),welltavg(1)
            IF (itemp <= 0) itemp = 1
            IF(itemp >= welltnumb) itemp = welltnumb
            btp(iplane,itemp,ip(n,1)) = btp(iplane,itemp,ip(n,1)) + &
			dble(iP(n,2))*p(n,4) !/(cellv*porosity(ploc(1),ploc(2),ploc(3))*Rtard(ip(n,1),ploc(1),ploc(2),ploc(3)))
!						p(n,4) 

			ip(n,4+iplane) = -1
		  END IF ! endif for particle crossing plane already
          END IF !endif for in plane
          
          
        END DO
      END IF   !endif for nplane=0
      
      
                 !   print*,' r c ',P(n,1),P(n,2),P(n,3)
!
! Calculate linear interp velocities devided by constituent's retardation fact
! Calculate timestep to reach cell bddys
!
      DO  l = 1, numax
        
!        vp(l) =  (v(l,ploc(1),ploc(2),ploc(3))*(1.d0-f(l))					&
!			  /Rtard(ip(n,1),ploc(1),ploc(2),ploc(3)) )   +				&
!            (v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*f(l)	&
!			/Rtard(ip(n,1),ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) )

! new mod for flux calc, 2/18/03

!p_av1 = 0.5d0* ( Porosity(ploc(1)-po(l,1),ploc(2)-po(l,2),ploc(3)-po(l,3)) + Porosity(ploc(1),ploc(2),ploc(3)) )
!print*, p_av1, Porosity(ploc(1),ploc(2),ploc(3))
!p_av2 = 0.5d0* ( Porosity(ploc(1),ploc(2),ploc(3)) + Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) )
!print*, p_av2,  Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))
!        vp(l) =  (v(l,ploc(1),ploc(2),ploc(3))*(1.d0-f(l))*p_av1) +					&
!            (v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*f(l)*p_av1)
! original way

		        vp(l) =  (v(l,ploc(1),ploc(2),ploc(3))*(1.d0-f(l))					&
		  *Porosity(ploc(1),ploc(2),ploc(3))*sat(ploc(1),ploc(2),ploc(3)) )   +				&
        (v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*f(l)	&
			*Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*Sat(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) )


!		mass_bal(i,j,k) = (-V(1,i-1,j,k)*0.5*(phi(i,j,k)+phi(i-1,j,k)) + V(1,i,j,k)*0.5*(phi(i,j,k)+phi(i+1,j,k))) +    &
						  !(-V(2,i,j-1,k)*0.5*(phi(i,j,k)+phi(i,j-1,k)) + V(2,i,j,k)*0.5*(phi(i,j,k)+phi(i,j+1,k))) +    &
						  !(-V(3,i,j,k-1)*0.5*(phi(i,j,k)+phi(i,j,k-1)) + V(3,i,j,k)*0.5*(phi(i,j,k)+phi(i,j,k+1))) 


!  V's are q's so we don't need to mult v by phi to get q
!		        vp(l) =  v(l,ploc(1),ploc(2),ploc(3))*(1.d0-f(l)) +					&
!                         v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*f(l)	

! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
        vp(l) =  vp(l) / (Rstar*scxyz(l,ploc(l)))

!
! uncomment for extra debugging and a HUGE log file
! write(666,*) 'v(',l,') ',Vp(l)
! if (abs(Vp(l)).eq.0.0d0) then
! write(666,*) 'velocity effectively zero'
! write(666,*) 'part num', n
! write(666,*) 'node i j k', ploc(1),ploc(2),ploc(3)
! write(666,*) 'direction:',l
! goto 151
! endif
        
!        b(l) = (1.0d0/delv(l))*  &
!            ( (v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))/            &
!			Rtard(ip(n,1),ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) ) -  &
!            (v(l,ploc(1),ploc(2),ploc(3))/Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))) )
 
!p_av1 = 0.5d0* ( Porosity(ploc(1)-po(l,1),ploc(2)-po(l,2),ploc(3)-po(l,3)) + Porosity(ploc(1),ploc(2),ploc(3)) )
!print*, p_av1, Porosity(ploc(1),ploc(2),ploc(3))
!p_av2 = 0.5d0* ( Porosity(ploc(1),ploc(2),ploc(3)) + Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) )

! old way
        b(l) = (1.0d0/delv(l))*  &
            ( (v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*            &
			Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))           &
			*sat(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) ) -           &
            v(l,ploc(1),ploc(2),ploc(3))*Porosity(ploc(1),ploc(2),ploc(3))*Sat(ploc(1),ploc(2),ploc(3)) )
!  V's are q's so we don't need to mult v by phi to get q
!        b(l) = (1.0d0/delv(l))*  &
!            ( v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))  -  &
!          v(l,ploc(1),ploc(2),ploc(3)) )

! new way
!       b(l) = (1.0d0/delv(l))*  &
!           ( (v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*            &
!			p_av2 ) -  &
!           (v(l,ploc(1),ploc(2),ploc(3))*p_av1) )

 
!        b(l) = (1.0d0/delv(l))*  &
!            ( (v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))*            &
!			(Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))/Rtard(ip(n,1),ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))) ) -  &
!            (v(l,ploc(1),ploc(2),ploc(3))*  &
!			(Porosity(ploc(1),ploc(2),ploc(3))/Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))) ) )
       
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
 
		b(l) =  b(l) /( Rstar * Scxyz(l,ploc(l))) 
!
! Check for direction of local flow, choose location of
! velocity appropriately
!
        
        IF (vp(l) >= 0.0D0) THEN
!          a1(l) = v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))	&
!				/Rtard(ip(n,1),ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))
!          a(l) = v(l,ploc(1),ploc(2),ploc(3))		&
!				/Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))

! old way
          a1(l) = v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))	&
				*Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) &
				*Sat(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))

          a(l) = v(l,ploc(1),ploc(2),ploc(3))		&
				*Porosity(ploc(1),ploc(2),ploc(3))  &
				*Sat(ploc(1),ploc(2),ploc(3))
!  V's are q's so we don't need to mult v by phi to get q
 !         a1(l) = v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))

  !       a(l) = v(l,ploc(1),ploc(2),ploc(3))		

! new way
!p_av1 = 0.5d0* ( Porosity(ploc(1)-po(l,1),ploc(2)-po(l,2),ploc(3)-po(l,3)) + Porosity(ploc(1),ploc(2),ploc(3)) )
!p_av2 = 0.5d0* ( Porosity(ploc(1),ploc(2),ploc(3)) + Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) )

!          a1(l) = v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))	&
!				*p_av2

 !         a(l) = v(l,ploc(1),ploc(2),ploc(3))		&
!				*p_av1

!		a(l) =  a(l) / Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
!		a1(l) =  a1(l) / Rtard(ip(n,1),ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))

          fact = 1.0D0-f(l)
        ELSE
 !         a(l) = v(l,ploc(1),ploc(2),ploc(3))/Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
!          a1(l) = v(l,ploc(1),ploc(2),ploc(3))/Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
! old way
          a(l) = v(l,ploc(1),ploc(2),ploc(3))*Porosity(ploc(1),ploc(2),ploc(3))*Sat(ploc(1),ploc(2),ploc(3))
          a1(l) = v(l,ploc(1),ploc(2),ploc(3))*Porosity(ploc(1),ploc(2),ploc(3))*Sat(ploc(1),ploc(2),ploc(3))

!  V's are q's so we don't need to mult v by phi to get q
!          a(l) = v(l,ploc(1),ploc(2),ploc(3))
!          a1(l) = v(l,ploc(1),ploc(2),ploc(3))

! new way
!p_av1 = 0.5d0* ( Porosity(ploc(1)-po(l,1),ploc(2)-po(l,2),ploc(3)-po(l,3)) + Porosity(ploc(1),ploc(2),ploc(3)) )
!p_av2 = 0.5d0* ( Porosity(ploc(1),ploc(2),ploc(3)) + Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) )
!          a(l) = v(l,ploc(1),ploc(2),ploc(3))*P_av1
!          a1(l) = v(l,ploc(1),ploc(2),ploc(3))*p_av1


!   		a(l) =  a(l) / Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
!   		a1(l) =  a1(l) / Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))

          fact = f(l)
        END IF
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
		a(l) =  a(l) /(Rstar*scxyz(l,ploc(l)))   
		a1(l) =  a1(l) / (Rstar*scxyz(l,ploc(l)))
!
! Check for uniform or linear interp flow, then
! determine delt by ds node location, screening
! for "stuck" particles
!
! first we check all possible cases that would break down analytical soln, 
! if any of these cases are true, we switch over to the old way of solving 
! av(dir) is the control variable, av(l) = 0 means analytical soln, 
! av(l)=1 means std euler integration,
! av(l) = 2 means some stuck situation, so we set the delta_t large for that dir
! 		
		av(l) = 0
!		av(l) = 1

! first for uniform flow
! old way
        IF (ABS( v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))    &
		*porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) &
		*Sat(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))	 -     &
            v(l,ploc(1),ploc(2),ploc(3))*                              &
			  porosity(ploc(1),ploc(2),ploc(3)) &
			  *Sat(ploc(1),ploc(2),ploc(3)) ) < epsi)  av(l)= 1

! test RMM
!av = 1

!  V's are q's so we don't need to mult v by phi to get q
!        IF (ABS( v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))  -     &
!             v(l,ploc(1),ploc(2),ploc(3)) ) < epsi)  av(l)= 1


! new way
!p_av1 = 0.5d0* ( Porosity(ploc(1)-po(l,1),ploc(2)-po(l,2),ploc(3)-po(l,3)) + Porosity(ploc(1),ploc(2),ploc(3)) )
!p_av2 = 0.5d0* ( Porosity(ploc(1),ploc(2),ploc(3)) + Porosity(ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3)) )
!        IF (ABS( v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))    &
!		*p_av2 -     &
!              v(l,ploc(1),ploc(2),ploc(3))*                              &
!			  p_av1 ) < epsi)  av(l)= 1


! then for zero velocity        
!		IF ( DABS(Vp(l)/delv(l)) < epsi) then
!		av(l) = 2
!		tad(l) = tnext 
!		end if
		
! then for other oddities w/ analytical soln
 !       IF ( DABS(a1(l)/delv(l)) < epsi) then
  !      av(l) = 2
!		tad(l) = tnext 
!		end if

! then for zero velocity        
		IF ( ABS(Vp(l)/delv(l)) < epsi) av(l) = 1
! then for other oddities w/ analytical soln
        IF ( ABS(a1(l)/delv(l))  < epsi) av(l) = 1
!!av(l) = 1

		if (av(l) == 0 ) then
		  if ( (a1(l)/vp(l)) < 0.) then
!tad(l) = (Porosity(ploc(1),ploc(2),ploc(3))/b(l))*DLOG(a1(l)/vp(l)) 
        !print*, Porosity(ploc(1),ploc(2),ploc(3)), b(l), DLOG(a1(l)/vp(l)),tad(l)
		   av(l) = 1
		   !tad(l) = tnext 			 
		   !print*,' a1, vp diff signs'
!		   print*,' n=',n
!		   print*,ploc(1),ploc(2),ploc(3)
!		   print*, l 
!		   print*,a1(l),vp(l)
!		   print*,a1(l)/vp(l)
		  end if
		end if

	    tad(l) = 1.E15
		
		!IF ( av(l) == 1 )  tad(l) = ABS( (0.01d0*delv(l)*fact) / ( vp(l)/ Porosity(ploc(1),ploc(2),ploc(3))  ) )
!		IF ( av(l) == 1 )  tad(l) = DABS( (0.05d0*delv(l)*fact) / ( vp(l)/ Porosity(ploc(1),ploc(2),ploc(3))  ) )
		IF ( av(l) == 1 )  tad(l) = DABS( (0.1d0*delv(l) ) / &
        ( vp(l)/ (Porosity(ploc(1),ploc(2),ploc(3))*Sat(ploc(1),ploc(2),ploc(3))  ) ) )
!		IF ( av(l) == 1 )  tad(l) = ABS( (delv(l)*0.1D0)/vp(l) )
        IF ( av(l) == 0 )   tad(l) = (Sat(ploc(1),ploc(2),ploc(3))*Porosity(ploc(1),ploc(2),ploc(3))/b(l))*DLOG(a1(l)/vp(l)) 
!        print*, Porosity(ploc(1),ploc(2),ploc(3)), b(l), DLOG(a1(l)/vp(l)),tad(l) 
		!end if
!		print*, l 
!		print*,a1(l),vp(l)
!		print*,a1(l)/vp(l)
!		print*,LOG(a1(l)/vp(l))


 !       end if
   
        END DO  ! l, vx,y,z loop
        
        
                 !   print*,' r c ',P(n,1),P(n,2),P(n,3)

        p(n,10) = DSQRT((vp(1)**2)+(vp(2)**2)+(vp(3)**2))
        
        tdd(1) = 1.E15
        tdd(2) = 1.E15
        tdd(3) = 1.E15
        tdd(4) = 1.E15
        tdd(5) = 1.E15
        tdd(6) = 1.E15
        tdd(7) = 1.E15
        tdd(8) = 1.E15
          if(massive_debug == 1) then
		  print*, '________________________'
		  print*, '*Part ',n,' IT=',it
		  print*, 'tadx=',tad(1)
		  print*, 'tady=',tad(2)
		  print*, 'tadz=',tad(3)
		  print*, ' mass=',p(n,4),' time=',p(n,7)
          print*, 'Vpx=',vp(1), ' Ax=',a(1), ' A1x=',a1(1)
          print*, 'Bx=',b(1),' fx=', f(1)
		  print*, 'avx=',av(1)

		  print*, 'Vpy=',vp(2), ' Ay=',a(2), ' A1y=',a1(2)
          print*, 'By=',b(2),' fy=', f(2)
		  print*, 'avy=',av(2)

		  print*, 'Vpz=',vp(3), ' Az=',a(3), ' A1z=',a1(3)
          print*, 'Bz=',b(3),' fz=', f(3)
		  print*, 'avz=',av(3)

          print*, 'ploc1=',ploc(1), ' ploc2=',ploc(2), ' ploc3=',ploc(3)
          print*, 'x=',p(n,1), ' y=',p(n,2),' z=',p(n,3)
	      print*, ' R=',Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
		  print*, ' phi=',porosity(ploc(1),ploc(2),ploc(3))
          print*, ' vxi',v(1,ploc(1),ploc(2),ploc(3))
          print*, ' vyj',v(2,ploc(1),ploc(2),ploc(3))
		  print*, ' vyk',v(3,ploc(1),ploc(2),ploc(3))
		  l = 1
		  print*, ' po(',l,',1)=', po(l,1)
          print*, ' vxi+1',v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))
		  l = 2
		  print*, ' po(',l,',2)=', po(l,2)
          print*, ' vyj+1',v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))
		  l = 3
		  print*, ' po(',l,',3)=', po(l,3)
          print*, ' vyk+1',v(l,ploc(1)+po(l,1),ploc(2)+po(l,2),ploc(3)+po(l,3))
		  end if

        
        !IF (((al > 0).OR.(at > 0)).AND.(iP(n,4) == 0)) THEN
! hack to turn of RW if we are close to z=0 bddy
!        IF (((ploc(3)<=2).or.(al > 0).OR.(at > 0)).AND.(iP(n,4) == 0)) THEN
                IF ( ( ploc(3) > 2 .and. ploc(3) < domax(3) - 1 ) .and. ((al > 0).OR.(at > 0)).AND.(iP(n,4) == 0)) THEN
 
!
! Okay, time to do random walk dispersion and correction term
!  this is done using Bi-Linear velocity interpolation (oh boy)
!  the delta-t is taken from tad(4) above
!  the correction factor is also BLI, taken as an
!  O(delx**2) central difference
! For method and form of dispersion tensor and subsequent correction
!  factor, please see Tompson, et al (1987)
! For method of velocity interpoltaion, please see Labolle, et al (1996)
!
! set up the capital F's
!
         tloc(1) = INT(ploc(1) + f(1) - 0.5)
         tloc(2) = INT(ploc(2) + f(2) - 0.5)
         tloc(3) = INT(ploc(3) + f(3) - 0.5)

          fbl(1) = ploc(1) + f(1) - 0.5 - DBLE(tloc(1))
          fbl(2) = ploc(2) + f(2) - 0.5 - DBLE(tloc(2))
          fbl(3) = ploc(3) + f(3) - 0.5 - DBLE(tloc(3))

!@RMM changed location of blinear interp per Toru's email
!          tloc(1) = ploc(1) 
!          tloc(2) = ploc(2) 
!          tloc(3) = ploc(3) 

!          fbl(1) = f(1) 
!          fbl(2) = f(2) 
!          fbl(3) = f(3) 

          
!
! Cycle thru and get surrounding velocity matrix using bilinear interp
!
          
          DO  i = 1, 3
            DO  j = 1, 3
              DO  k = 1, 3
                
                tp1 = tloc(1) -2 + i
                tp2 = tloc(2) -2 + j
                tp3 = tloc(3) -2 + k
                
                
! if we have actual calculated velocities (not fluxes) we use the following
!
!                vbl(1,i,j,k) = (1.d0/Rtard(ip(n,1),tp1,tp2,tp3))*       &
!				    (1.0D0-fbl(2))*(1.0D0-fbl(3))*(v(1,tp1,tp2,tp3)) +  &
!                    (fbl(2))*(1.0D0-fbl(3))* (v(1,tp1,tp2+1,tp3)) +     &
!                    (1.0D0-fbl(2))*(fbl(3))*(v(1,tp1,tp2,tp3+1)) +      &
!                    (fbl(2))*(fbl(3))*(v(1,tp1,tp2+1,tp3+1))
                
!                vbl(2,i,j,k) = (1.d0/Rtard(ip(n,1),tp1,tp2,tp3))*       &
!				    (1.0D0-fbl(1))*(1.0D0-fbl(3))*(v(2,tp1,tp2,tp3)) +  &
!                    (fbl(1))*(1.0D0-fbl(3))*(v(2,tp1+1,tp2,tp3)) +      &
!                    (1.0D0-fbl(1))*(fbl(3))*(v(2,tp1,tp2,tp3+1)) +      &
!                    (fbl(1))*(fbl(3))*(v(2,tp1+1,tp2,tp3+1))
                
!                vbl(3,i,j,k) = (1.d0/Rtard(ip(n,1),tp1,tp2,tp3))*       &
!				    (1.0D0-fbl(1))*(1.0D0-fbl(2))*(v(3,tp1,tp2,tp3)) +  &
 !                   (fbl(1))*(1.0D0-fbl(2))*(v(3,tp1+1,tp2,tp3)) +      &
 !                   (1.0D0-fbl(1))*(fbl(2))*(v(3,tp1,tp2+1,tp3)) +      &
 !                   (fbl(1))*(fbl(2))*(v(3,tp1+1,tp2+1,tp3))
           

! this is some hybrid thing w/ no sat and no retard
!               
!                vbl(1,i,j,k) = (1.d0/(porosity(tloc(1),tloc(2),tloc(3))))*       &
!				    (1.0D0-fbl(2))*(1.0D0-fbl(3))*(porosity(tp1,tp2,tp3)*v(1,tp1,tp2,tp3)/Rtard(ip(n,1),tp1,tp2,tp3)) +  &
!                    (fbl(2))*(1.0D0-fbl(3))* (porosity(tp1,tp2+1,tp3)*v(1,tp1,tp2+1,tp3)/Rtard(ip(n,1),tp1,tp2+1,tp3)) +     &
!                    (1.0D0-fbl(2))*(fbl(3))*(porosity(tp1,tp2,tp3+1)*v(1,tp1,tp2,tp3+1)/Rtard(ip(n,1),tp1,tp2,tp3+1)) +      &
!                    (fbl(2))*(fbl(3))*(porosity(tp1,tp2+1,tp3+1)*v(1,tp1,tp2+1,tp3+1)/Rtard(ip(n,1),tp1,tp2+1,tp3+1))
               		        
!                vbl(1,i,j,k) = (1.0D0-fbl(2))*(1.0D0-fbl(3))*(v(1,tp1,tp2,tp3)/Rtard(ip(n,1),tp1,tp2,tp3)) +  &
!                    (fbl(2))*(1.0D0-fbl(3))* (v(1,tp1,tp2+1,tp3)/Rtard(ip(n,1),tp1,tp2+1,tp3)) +   &
!                    (1.0D0-fbl(2))*(fbl(3))*(v(1,tp1,tp2,tp3+1)/Rtard(ip(n,1),tp1,tp2,tp3+1)) +   &
!                    (fbl(2))*(fbl(3))*(v(1,tp1,tp2+1,tp3+1)/Rtard(ip(n,1),tp1,tp2+1,tp3+1))
!
!  this complicated mess is the bilinear interp 
!
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),tloc(1),tloc(2),tloc(3))-1.d0)/sat(tloc(1),tloc(2),tloc(3))

!              ! add in test for boundary valuem otherwise wire vel to zero
!             if (tp3>=1) then
                vbl(1,i,j,k) = ( 1.d0/( porosity(tloc(1),tloc(2),tloc(3))*sat(tloc(1),tloc(2),tloc(3))*Rstar ))*       &
				    ((1.0D0-fbl(2))*(1.0D0-fbl(3))*(porosity(tp1,tp2,tp3)*sat(tp1,tp2,tp3)*v(1,tp1,tp2,tp3)) +  &
                    (fbl(2))*(1.0D0-fbl(3))* (porosity(tp1,tp2+1,tp3)*sat(tp1,tp2+1,tp3)*v(1,tp1,tp2+1,tp3)) +     &
                    (1.0D0-fbl(2))*(fbl(3))*(porosity(tp1,tp2,tp3+1)*sat(tp1,tp2,tp3+1)*v(1,tp1,tp2,tp3+1)) +      &
                    (fbl(2))*(fbl(3))*(porosity(tp1,tp2+1,tp3+1)*sat(tp1,tp2+1,tp3+1)*v(1,tp1,tp2+1,tp3+1)))

                vbl(1,i,j,k) = vbl(1,i,j,k)/scxyz(1,tloc(1))

                vbl(2,i,j,k) = ( 1.d0/(porosity(tloc(1),tloc(2),tloc(3))*sat(tloc(1),tloc(2),tloc(3))*Rstar ))*       &
				    ((1.0D0-fbl(1))*(1.0D0-fbl(3))*(porosity(tp1,tp2,tp3)*sat(tp1,tp2,tp3)*v(2,tp1,tp2,tp3)) +  &
                    (fbl(1))*(1.0D0-fbl(3))*(porosity(tp1+1,tp2,tp3)*sat(tp1+1,tp2,tp3)*v(2,tp1+1,tp2,tp3)) +      &
                    (1.0D0-fbl(1))*(fbl(3))*(porosity(tp1,tp2,tp3+1)*sat(tp1,tp2,tp3+1)*v(2,tp1,tp2,tp3+1)) +      &
                    (fbl(1))*(fbl(3))*(porosity(tp1+1,tp2,tp3+1)*sat(tp1+1,tp2,tp3+1)*v(2,tp1+1,tp2,tp3+1)))
                
                vbl(2,i,j,k) = vbl(2,i,j,k)/scxyz(2,tloc(2))
                vbl(3,i,j,k) = ( 1.d0/( porosity(tloc(1),tloc(2),tloc(3))*sat(tloc(1),tloc(2),tloc(3))*Rstar ))*       &
				    ((1.0D0-fbl(1))*(1.0D0-fbl(2))*(porosity(tp1,tp2,tp3)*sat(tp1,tp2,tp3)*v(3,tp1,tp2,tp3)) +  &
                    (fbl(1))*(1.0D0-fbl(2))*(porosity(tp1+1,tp2,tp3)*sat(tp1+1,tp2,tp3)*v(3,tp1+1,tp2,tp3)) +      &
                    (1.0D0-fbl(1))*(fbl(2))*(porosity(tp1,tp2+1,tp3)*sat(tp1,tp2+1,tp3)*v(3,tp1,tp2+1,tp3)) +      &
                    (fbl(1))*(fbl(2))*(porosity(tp1+1,tp2+1,tp3)*sat(tp1+1,tp2+1,tp3)*v(3,tp1+1,tp2+1,tp3)))

                vbl(3,i,j,k) = vbl(3,i,j,k)/scxyz(3,tloc(3))
 !            else
 !               vbl(1,i,j,k) = 0.d0
  !              vbl(2,i,j,k) = 0.d0
   !             vbl(3,i,j,k) = 0.d0
    !          end if
  
                139 CONTINUE
                
                
              END DO
            END DO
          END DO
          
!
! Now that we've built up a matrix of BLI Velocities, let's crank those
! Derivatives of  the dispersion tensor
!
          
          alt = al - at
          vnp1 = DSQRT(vbl(1,3,2,2)**2 + vbl(2,3,2,2)**2 +vbl(3,3,2,2)**2)
          vnm1 = DSQRT(vbl(1,1,2,2)**2 + vbl(2,1,2,2)**2 +vbl(3,1,2,2)**2)
          
          if(vnp1*vnm1 /= 0.) then
          
          dxx = (alt/(2*delv(1) )) *( ((vbl(1,3,2,2)**2)/vnp1) -   &
           ((vbl(1,1,2,2)**2)/vnm1))  + (at/(2*delv(1)))*(vnp1 - vnm1)
          
          dxy = (alt/(2*delv(1))) *( ((vbl(1,3,2,2)*vbl(2,3,2,2))/vnp1) -  &
           ((vbl(1,1,2,2)*vbl(2,1,2,2))/vnm1))
          
          
          dxz = (alt/(2*delv(1))) *( ((vbl(1,3,2,2)*vbl(3,3,2,2))/vnp1) -	&
           ((vbl(1,1,2,2)*vbl(3,1,2,2))/vnm1))
           CALL addGasDispersionX( n, p, ip, delv, tloc, fbl, porosity, sat, &
                                   dxx, modelname )
          else
		  dxx = 0.d0
		  dxy = 0.d0
		  dxz = 0.d0
          end if

          vnp1 = DSQRT(vbl(1,2,3,2)**2 + vbl(2,2,3,2)**2+vbl(3,2,3,2)**2)
          vnm1 = DSQRT(vbl(1,2,1,2)**2 + vbl(2,2,3,2)**2+vbl(3,3,1,2)**2)
          
		  if(vnp1*vnm1 /= 0.) then
          dyy = (alt/(2*delv(2))) *( ((vbl(2,2,3,2)**2)/vnp1) -	&
           ((vbl(2,2,1,2)**2)/vnm1))  + (at/(2*delv(2)))*(vnp1 - vnm1)
          
          dyx = (alt/(2*delv(2))) *( ((vbl(1,2,3,2)*vbl(2,2,3,2))/vnp1) -	&
           ((vbl(1,2,1,2)*vbl(2,2,1,2))/vnm1))
          
          dyz = (alt/(2*delv(2))) *( ((vbl(2,2,3,2)*vbl(3,2,3,2))/vnp1) -	&
           ((vbl(2,2,1,2)*vbl(3,2,1,2))/vnm1))
           
           CALL addGasDispersionY( n, p, ip, delv, tloc, fbl, porosity, sat, &
                                   dyy, modelname )
			else
		  dyy = 0.d0
		  dyx = 0.d0
		  dyz = 0.d0
          end if

          vnp1 = DSQRT(vbl(1,2,2,3)**2 + vbl(2,2,2,3)**2+vbl(3,2,2,3)**2)
          vnm1 = DSQRT(vbl(1,2,2,1)**2 + vbl(2,2,2,1)**2+vbl(3,3,2,1)**2)
          if(vnp1*vnm1 /= 0.) then
          alt = al - at
!@RMM dzz term, delv(2) changed to delv(3) per Toru's email
          
          dzz = (alt/(2*delv(3))) *( ((vbl(3,2,2,3)**2)/vnp1) -		&
           ((vbl(3,2,2,3)**2)/vnm1))  + (at/(2*delv(2)))*(vnp1 - vnm1)
          
          dzx = (alt/(2*delv(3))) *( ((vbl(3,2,2,3)*vbl(1,2,2,3))/vnp1) -	&
           ((vbl(3,2,2,1)*vbl(1,2,2,3))/vnm1))
          
          dzy = (alt/(2*delv(3))) *( ((vbl(3,2,2,3)*vbl(2,2,2,3))/vnp1) -	&
           ((vbl(3,2,2,1)*vbl(2,2,2,1))/vnm1))
           CALL addGasDispersionZ( n, p, ip, delv, tloc, fbl, porosity, sat, &
                                   dzz, modelname )
			else
		  dzz = 0.d0
		  dzx = 0.d0
		  dzy = 0.d0
          end if
          
!
! Let's do the random walk, starting with random number gen
!
          
          z(1) = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
          z(2) = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
          z(3) = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
          
 ! @RMM changed random number gen from ran1 to ran         
          
          
          vn = DSQRT(vbl(1,2,2,2)**2 + vbl(2,2,2,2)**2 +vbl(3,2,2,2)**2)
          
          vxz = vbl(1,2,2,2)*vbl(3,2,2,2)
          vxx = vbl(1,2,2,2)**2
          vzz = vbl(3,2,2,2)**2
          vyy = vbl(2,2,2,2)**2
          vxy = vbl(1,2,2,2)*vbl(2,2,2,2)
          vyz = vbl(2,2,2,2)*vbl(3,2,2,2)
          
          
          
          betad = DSQRT(vn**2 + vbl(2,2,2,2)**2 + 2.0D0*vxz)
          
          IF (betad /= 0.d0) THEN
            cx = z(1)*DSQRT((2.0D0*al)/vn)
            cy = z(2)*DSQRT(2.0D0*at*vn)/betad
            cz = z(3)*DSQRT(2.0D0*at*vn)/(betad*vn)
          ELSE
            cx = 0.
            cy = 0.
            cz = 0.
          END IF


            cx = z(1)*DSQRT((2.0D0*al)/vn)
            cy = z(2)*DSQRT(2.0D0*at*vn)/betad
            cz = z(3)*DSQRT(2.0D0*at*vn)/(betad*vn)
  
                      
          dx = cx*vbl(1,2,2,2) - cy*vbl(2,2,2,2) - cz*(vyy+vzz+vxz)
          dy = cx*vbl(2,2,2,2) + cy*(vbl(1,2,2,2)+vbl(3,2,2,2))	&
               + cz*(vxy-vyz)
          dz = cx*vbl(3,2,2,2) - cy*vbl(2,2,2,2) + cz*(vxx+vyy+vxz)
          
          IF((dxx+dxy+dxz) /= 0.0D0) THEN
            tdd(1) = DABS(dlimit*delv(1)/(dxx+dxy+dxz))
          END IF
          IF((dyy+dyx+dyz) /= 0.0D0) THEN
            tdd(2) = DABS(dlimit*delv(2)/(dyy+dyx+dyz))
          END IF
          
         ! IF((dzz+dzx+dzy) > 0.0D0) THEN
		  IF((dzz+dzx+dzy) /= 0.0D0) THEN
            tdd(3) = DABS(dlimit*delv(3)/(dzz+dzx+dzy))
          END IF
          
          IF(dx /= 0.0D0) THEN
            tdd(4) = (dlimit*delv(1)/dx)**2
          END IF
          IF(dy /= 0.0D0) THEN
            tdd(5) = (dlimit*delv(2)/dy)**2
          END IF
          IF(dz /= 0.0D0) THEN
            tdd(6) = (dlimit*delv(3)/dz)**2
          END IF

        END IF !  IF (((ploc(3)<=2).or.(al > 0).OR.(at > 0)).AND.(iP(n,4) == 0)) THEN
        
        
        
!
! Figure out which direction has largest the velocity, or
! if it's time to re-group
!
        
        
        trgp = tnext - p(n,7)
        if(massive_debug == 1) then
		print*
		print*, ' trgp', trgp,' tnext',tnext,' p(n,7)',p(n,7)
		print*
		end if

! check for flow sink (all velocities pointing into the node)
! set tad(1-3) = trgp (we wait to see if flowfield changes later)
! print out warning message
       cell_sink = 0
!        IF( (V(1,ploc(1),ploc(2),ploc(3)) > 0.d0).and.(V(1,ploc(1)+1,ploc(2),ploc(3)) < 0.d0) )  THEN
!		  IF(  (V(2,ploc(1),ploc(2),ploc(3))> 0.d0) .and. (V(2,ploc(1),ploc(2)+1,ploc(3)) < 0.d0) )  THEN
!		    IF( (V(3,ploc(1),ploc(2),ploc(3))> 0.d0) .and. (V(3,ploc(1),ploc(2),ploc(3)+1) < 0.d0) ) THEN
!			WRITE(667,*) ' cell sink at node:',ploc(1),ploc(2),ploc(3)
!			WRITE(667,*) ' time', time,' ts',it
!			WRITE(667,*) 
!            tad(1:3) = trgp
!			cell_sink = 1
!		    END IF
!		  END IF
!		END IF



! write(*,*) trgp, tnext, P(n,7)
!     1   -(ncut-loop2)*(dt/dble(ncut))
! if (n.eq.399) write(666,*) ' trgp=',trgp
!        if (tad(1) <= epsi) tad(1) = trgp
!        if (tad(2) <= epsi) tad(2) = trgp
!        if (tad(3) <= epsi) tad(3) = trgp

!        tad(1) =  ABS((delv(1)*.1d0)/vp(1))
!        tad(2) =  ABS((delv(2)*.1d0)/vp(2))
!        tad(3) =  ABS((delv(3)*.1d0)/vp(3))
!		av = 1
		dt_lamda = 1.E15
		IF (half_life(iP(n,1)) > 0.d0) dt_lamda = 0.1d0*half_life(iP(n,1))
! hard-wire OUT atta/det to save memory; flag katt RMM 8-31-05		
		if (K_att(iP(n,1),ploc(1),ploc(2),ploc(3)) == 0.) THEN
		t_att = trgp
		ELSE
		t_att = 0.01D0/K_att(iP(n,1),ploc(1),ploc(2),ploc(3))
		END IF
!        tad(4) = DMIN1(tad(1),tad(2),tad(3),trgp,tdd(1),tdd(2),tdd(3),dt_lamda,10.) 
!        tad(4) = DMIN1(tad(1),tad(2),tad(3),trgp,tdd(1),tdd(2),tdd(3),tdd(4),tdd(5),tdd(6),t_att,dt_lamda) 
        ! change @RMM in limit for dt for dispersion
        tad(4) = DMIN1(tad(1),tad(2),tad(3),trgp,t_att,dt_lamda) 

		if(tad(4) == 0.) tad(4) = tad(4) + epsi*10.D0

        if(iP(n,2) == 0 ) then
		  if( K_det(iP(n,1),ploc(1),ploc(2),ploc(3)) == 0.D0 ) then
		    tad(4) = trgp
		  else
		    tad(4) = DMIN1(0.05D0/K_det(iP(n,1),ploc(1),ploc(2),ploc(3)),trgp,dt_lamda)
!		    tad(4) = DMIN1(trgp,tnext/100.D0)
         end if
		end if
if (massive_debug == 1) then
print*
print*,'time check'
print*, 'iP(n,2)',iP(n,2)
print*, 'tad(1)',tad(1)
print*, 'tad(2)',tad(2)
print*, 'tad(3)',tad(3)
print*, 'tad(4)',tad(4)
print*, 'trgp',trgp
print*, 'tdd(1)',tdd(1)
print*, 'tdd(2)',tdd(2)
print*, 'tdd(3)',tdd(3)
print*, 't_att', t_att
print*, 'dt_lamda',dt_lamda
print*, 'cell sink', cell_sink
print*
end if

		p(n,8) = tad(4) 
        
! if (n.eq.399) write(666,*) ' tad4=',tad(4)
! if (n.eq.399) write(666,*) ' t1-3=',tad(1), tad(2), tad(3)
        
        movedir = 99

! if (tad(4).eq.tad(1)) movedir=1
!       if (tad(4).eq.tad(2)) movedir=2
!        if (tad(4).eq.tad(3)) movedir=3
! write(666,*) 'movedir:',movedir
!
! Move them darn particles
!
        oldloc(1) = p( n, 1 )
        oldloc(2) = p( n, 2 )
        oldloc(3) = p( n, 3 )
        DO  l = 1, numax
          IF (av(l) == 0) THEN
            
            extemp = (b(l)*tad(4))/(Porosity(ploc(1),ploc(2),ploc(3))*sat(ploc(1),ploc(2),ploc(3)))
            move = (vp(l)*DEXP(extemp) - a(l))/b(l)
!			if(l <> 1) then
!			  if (move <> 0.) print*,' move<>0 ',move,extemp
!			end if
! if cell sink = 1, then move = 0
            IF( cell_sink == 1) move = 0.0
!
! to further combat roundoff error, adjust particle displacement
! in direction of limiting move
!
! if (movedir.eq.l) then
! if (Vp(l).lt.0.0) Move = 0.0d0
! if (Vp(l).gt.0.0) Move = 1.0d0*delv(l)
! write(666,*) 'movedir:',movedir,' move:',move
! end if
!
! do some error checking and write out messages in log file
!
!		if (iP(n,2) == 10) then
            IF (move < -epsi) THEN
              WRITE(667,*) 'Move is zero or less move=', move
              WRITE(667,*) 'tad=',tad(4), 'dir=',l
              WRITE(667,*) 'Vp=',vp(l), '    A=',a(l)
              WRITE(667,*) 'A1=',a1(l)
              WRITE(667,*) 'extemp=',extemp, '  dexp(extemp)=',DEXP(extemp)
              WRITE(667,*) 'Vp*e()',vp(l)*DEXP(extemp)
              WRITE(667,*) 'Vp*e() - a', (vp(l)*DEXP(extemp))-a(l)
              WRITE(667,*) 'V(i+1,j,k)=',v(l,ploc(1)+1,ploc(2),ploc(3))
              WRITE(667,*) 'B=',b(l),'  fact=',fact, ' f=',f(l)
              WRITE(667,*) 'x=',p(n,1), 'y=',p(n,2)
              WRITE(667,*) 'ploc1=',ploc(1), 'ploc2=',ploc(2)
              WRITE(667,*)  ' ploc3=',ploc(3)
              WRITE(667,*)  'n ',n
              WRITE(667,*)  'nt ', i
              WRITE(667,*)  'phi ', porosity(ploc(1),ploc(2),ploc(3))
              WRITE(667,*)  'R ', rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
              WRITE(667,*)  ' tad(1:n)' , tad
              WRITE(667,*)
			  MOVE = 0.d0
            END IF
            
            IF (move > delv(l) + epsi) THEN
              WRITE(667,*) 'Move is > delta: move=',move
              WRITE(667,*) 'tad=',tad(4),' dir=', l
              WRITE(667,*) 'Vp=',vp(l), ' A=',a(l), ' A1=',a1(l)
              WRITE(667,*) 'B=',b(l),' fact=', fact,' f=', f(l)
			  WRITE(667,*) 'av=',av
	          WRITE(667,*) 'extemp=',extemp, '  dexp(extemp)=',DEXP(extemp)
              WRITE(667,*) 'ploc1=',ploc(1), ' ploc2=',ploc(2), ' ploc3=',ploc(3)
              WRITE(667,*) 'x=',p(n,1), ' y=',p(n,2)
			  WRITE(667,*) 'n=',n, ' R=',Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
              WRITE(667,*)
			  move = delv(l)
            END IF
!		end if
            
! Update new location in particle array - note move is mulitplied by iP(n,2)
!
! for SL exit location
!            if (move <> move*dble(iP(n,2)) ) print*,move,move*dble(iP(n,2))
			ptemp = p(n,l)
!			p(n,l) = DBLE(ploc(l))*delv(l) + move*dble(iP(n,2))    &
!			 - del2v(l)
!			p(n,l) = DBLE(ploc(l))*delv(l) + move*dble(iP(n,2))    &
!			 + del2v(l)
			if (iP(n,2) == 1) p(n,l) = dfloat(ploc(l)-1)*delv(l) + move  
!            p(n,l) = DBLE(ploc(l))*delv(l) + move*dble(iP(n,2))  

!			if(l <> 1) then
!			  if (ptemp <> p(n,l)) print*,' ptemp <> p(n,l)',ptemp,p(n,l),n,l
!			end if
          ELSE
!
! for uniform cell velocity
!
            if (iP(n,2) == 1)									&
			p(n,l) = p(n,l) + tad(4)*(vp(l)/(Porosity(ploc(1),ploc(2),ploc(3))*sat(ploc(1),ploc(2),ploc(3)) ) )
!			if ( tad(4)*vp(l) <> tad(4)*vp(l)*dble(iP(n,2)) ) print*,tad(4)*vp(l),tad(4)*vp(l)*dble(iP(n,2))
!			if(l <> 1) then
!			  if (vp(l) <> 0.) print*,' v<>0 ',vp(l)
!			end if
          END IF
          
		  if(massive_debug == 1) then
		  print*, 'move ',move
		  print*,' vp(l)*tad(4)',vp(l)*tad(4)*dble(iP(n,2)) 
		  print*, 'tad=',tad(4),' dir=', l
          print*, 'Vp=',vp(l), ' A=',a(l), ' A1=',a1(l)
          print*, 'B=',b(l),' fact=', fact,' f=', f(l)
		  print*, 'av=',av
          print*, 'ploc1=',ploc(1), ' ploc2=',ploc(2), ' ploc3=',ploc(3)
          print*, 'x=',p(n,1), ' y=',p(n,2),' z=',p(n,3)
	      print*, 'n=',n, ' R=',Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))
		  print*, ' phi=',porosity(ploc(1),ploc(2),ploc(3))
		  end if
        END DO
        9194 CONTINUE
        
       CALL addGasDiffusion( n, ip, tloc, porosity, sat, moldiff, mldif,&
                  modelname  )

        IF (((al > 0.0D0).OR.(at > 0.0D0).OR.(moldiff > 0.0D0)).AND.(iP(n,4) == 0) ) THEN
        !disable dispersion at the source cell
!        IF (((al > 0.0D0).OR.(at > 0.0D0).OR.(moldiff > 0.0D0)).AND.(iP(n,4) == 0) .AND. ( p(n, 1) > 15.0 )) THEN
! Okay, now we add RW component and correction factor to our current
!  position
!
          z(4) = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
          z(5) = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
          z(6) = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
          
          ddx = z(4) * DSQRT(moldiff*2.0D0*tad(4)) *dble(iP(n,2))
          ddy = z(5) * DSQRT(moldiff*2.0D0*tad(4)) *dble(iP(n,2))
          ddz = z(6) * DSQRT(moldiff*2.0D0*tad(4)) *dble(iP(n,2))
          
         
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3)) 
          p(n,1) = p(n,1) + (dxx+dxy+dxz)*tad(4)*dble(iP(n,2)) + dx*DSQRT(tad(4))*dble(iP(n,2)) + ddx/Rstar
          p(n,2) = p(n,2) + (dyy+dyx+dyz)*tad(4)*dble(iP(n,2)) + dy*DSQRT(tad(4))*dble(iP(n,2)) + ddy/Rstar
          p(n,3) = p(n,3) + (dzz+dzx+dzy)*tad(4)*dble(iP(n,2)) + dz*DSQRT(tad(4))*dble(iP(n,2)) + ddz/Rstar
 !         p(n,1) = p(n,1)  + dx*DSQRT(tad(4))*dble(iP(n,2)) + ddx
 !         p(n,2) = p(n,2) + dy*DSQRT(tad(4))*dble(iP(n,2)) + ddy
 !         p(n,3) = p(n,3)  + dz*DSQRT(tad(4))*dble(iP(n,2)) + ddz
! 
! quickie grad_phi_term addition
!          p(n,1) = p(n,1) + (1.0D0/Porosity(ploc(1),ploc(2),ploc(3)))*(grad_phi_x)*tad(4)
!          p(n,2) = p(n,2) + (1.0D0/Porosity(ploc(1),ploc(2),ploc(3)))*(grad_phi_y)*tad(4)*dble(iP(n,2)) 
!          p(n,3) = p(n,3) + (1.0D0/Porosity(ploc(1),ploc(2),ploc(3)))*(grad_phi_z)*tad(4)*dble(iP(n,2)) 

        END IF

!       ! reflection algorithm
!       IF( p(n,3) > 4.3 ) THEN
!        p(n,3) =4.3 - ( p(n, 3) - 4.3)
!       ENDIF 
!       IF( p(n,1) < 14 ) THEN
!        p(n,1) =14.0 + 14.0 - p(n, 1) 
!       ENDIF 
       !MacQ1990 
!       IF( p(n,3) > 0.3 ) THEN
!        p(n,1) =0.3 - ( p(n, 3) - 0.3 )
!       ENDIF 

       IF ( Zup_Ref .EQ. 1 ) THEN
         IF( p(n, 3) >  bnd_Zup ) p(n,3) = bnd_Zup - ( p(n,3) - bnd_Zup ) 
       ENDIF
       IF ( Zdown_Ref .EQ. 1 ) THEN
         IF( p(n, 3) <  bnd_Zdown ) p(n,3) = bnd_Zdown - ( p(n,3) - bnd_Zdown ) 
       ENDIF
       IF ( Yup_Ref .EQ. 1 ) THEN
         IF( p(n, 2) >  bnd_Yup ) p(n,2) = bnd_Yup - ( p(n,2) - bnd_Yup ) 
       ENDIF
       IF ( Ydown_Ref .EQ. 1 ) THEN
         IF( p(n, 2) <  bnd_Ydown ) p(n,2) = bnd_Ydown - ( p(n,2) - bnd_Ydown ) 
       ENDIF
       IF ( Xup_Ref .EQ. 1 ) THEN
         IF( p(n, 1) >  bnd_Xup ) p(n,1) = bnd_Xup - ( p(n,1) - bnd_Xup ) 
       ENDIF
       IF ( Xdown_Ref .EQ. 1 ) THEN
         IF( p(n, 1) <  bnd_Xdown ) p(n,1) = bnd_Xdown - ( p(n,1) - bnd_Xdown ) 
       ENDIF
       
      
! Now we add decay/ingrowth
!

  IF (half_life(iP(n,1)) > 0.d0) THEN
   lamda = DLOG(2.D0)/half_life(iP(n,1))
   Prob =  1.0d0 -  DEXP(-lamda*tad(4))
   !Zhl = rand2(ir)
   Zhl = ran1(ir)
	!Prob =  half_life(iP(n,1))*tad(4)
		IF (Zhl.le.Prob) THEN
			iP(n,1) = iP(n,1) + 1
			if (iP(n,1) > n_constituents ) iP(n,1) = 1
			go to 399
		END IF ! decay/ingrow?
  END IF ! half life bigger than zero?
!
! Here we leave space for more stuff like this, matrix diffusion, etc
!
!
919 continue

  IF (iP(n,2) == 0) THEN
  trgp = tnext - p(n,7)
  tad(4) = dmin1(0.01D0/K_det(iP(n,1),ploc(1),ploc(2),ploc(3)),trgp)
    IF (K_det(iP(n,1),ploc(1),ploc(2),ploc(3)) > 0.d0) THEN
   Prob = K_det(iP(n,1),ploc(1),ploc(2),ploc(3))*tad(4)
!   lamda = DLOG(2.D0)/K_det(iP(n,1),ploc(1),ploc(2),ploc(3))
!   Prob =  DEXP(-lamda*tad(4))
   Zhl = ran1(ir)
		IF (Zhl.le.Prob) THEN
			iP(n,2) = 1
		END IF ! change from att-> det
	END IF ! Kdet >0
  END IF ! IP(m,2) = 1, check if we are even doing kinetics
  IF (iP(n,2) == 1) THEN
!	else
! katt change memory saver 8-05
	IF (K_att(iP(n,1),ploc(1),ploc(2),ploc(3)) > 0.d0) THEN
  Prob = K_att(iP(n,1),ploc(1),ploc(2),ploc(3))*tad(4)
 ! lamda = DLOG(2.D0)/K_att(iP(n,1),ploc(1),ploc(2),ploc(3))
 !  Prob =  DEXP(-lamda*tad(4))
   Zhl = ran1(ir)
		IF (Zhl.le.Prob) THEN
			iP(n,2) = 0
		END IF ! 
	END IF ! k_att > 0

  END IF ! IP(n,2) ==1
!

399 continue

!
! Update particle time in array and location after move
!
        p(n,7) = p(n,7) +tad(4)
		ploc(1) = IDINT(p(n,1)/delv(1)) + 1
		ploc(2) = IDINT(p(n,2)/delv(2)) + 1
		ploc(3) = IDINT(p(n,3)/delv(3)) + 1

		if(massive_debug == 1) then
        print*
		print*, 'updating time: p(n,7)',p(n,7),' tad(4)',tad(4)
        end if

! if it's time to regroup, add to current values of conc.
! if not keep looping.
!

        145 CONTINUE
        
        IF (passes > give_up) THEN
          WRITE(667,*) 'gave up ',give_up
          WRITE(667,*) 'P',n
		  WRITE(667,*) 'I',it
		  if(massive_debug == 1) then
           DO ll = 1, 3
             WRITE(667,*) 'tad=',tad(ll),' dir=', ll
             WRITE(667,*) 'Vp=',vp(ll),'  f=' , f(ll)
		   	 WRITE(667,*) 'Av=',av(ll)
           END DO
           WRITE(667,*) 'tad(4)=',tad(4)
           WRITE(667,*) 'x=',p(n,1),' y=', p(n,2), ' z=',p(n,3)
           WRITE(667,*) 'l1=',ploc(1), 'l2=',ploc(2),'l3=',ploc(3)
		   WRITE(667,*)
		  end if          
          done = 1
          GO TO 149
        END IF
        
        
        IF (tad(4) == trgp) THEN
          
! if(loop2.eq.ncut) then
        if(massive_debug == 1) then
		print*, ' DONE, moving to next particle'
		print*
		end if
          done = 1
          GO TO 151
        END IF
! even if we give up on the particle we still track concentrations    
	149 continue    
        IF (tempavg == 1) THEN
          inbounds = 1
          DO  l = 1, numax
            ploc(l) = IDINT(p(n,l)/delc(l)) + 1
! change in ploc test
!             ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
!            IF((ploc(l) <= 1).OR.(ploc(l) >= domax(l) )) THEN
            IF((ploc(l) <= 0).OR.(ploc(l) > domax(l) )) THEN
              inbounds = 0
            END IF
          END DO
          
          IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
            c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
               SNGL( dble(iP(n,2))*(p(n,4)*tad(4))/  &
               (tnext*cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )
          END IF
        END IF   ! temporal averaging
        
! loop2 = loop2 + 1
! endif
        
        IF (partprint == 1) THEN
          IF ((p(n,1) > 0.).AND.(p(n,2) > 0.).AND.(p(n,3) > 0.)) then 
                       pdist = 0.0d0
                       do ii = 1, 3
                         pdist = pdist +  (p(n,ii)-lastprint(n,ii))**2 
                       end do
                       pdist = dsqrt(pdist)
                       if (pdist >= minpdist) then
                       !WRITE(13,61) p(n,1), p(n,2), p(n,3), p(n,5), p(n,7)
                       lastprint(n,1) = p(n,1)
                       lastprint(n,2) = p(n,2)
                       lastprint(n,3) = p(n,3)
                 end if ! if pdist
          end if ! in bounds?
        END IF ! partprint ==1?
        IF (partprint == 2) THEN
         ! IF ((p(n,1) > 0.).AND.(p(n,2) > 0.).AND.(p(n,3) > 0.))	&
           !WRITE(13,62) p(n,1), p(n,2), p(n,3), p(n,7)				&
           !,p(n,8),p(n,9),p(n,10)
          
        END IF
        
        
        END DO  ! correpsonding to tnext part loop? DO WHILE (done== 0)
      
151 CONTINUE

!       IF( p(n, 1) .GE. 15.0 .AND. origloc( 1 ) < 15.0 ) THEN
!         mp = mp + 1
!         p(mp, 1 ) = origloc(1) 
!         p(mp, 2 ) = origloc(2) 
!         p(mp, 3 ) = origloc(3) 
!         p(mp,4) = p(n,4) 
!         p(mp,5) = DBLE(np)
!         p(mp,7) = p(n,7) 
!	ip(mp,1) = ip(n,1)
!	ip(mp,2) = ip(n,2)
!	ip(mp,3) = ip(n,3)
!	ip(mp,4) = ip(n,4)
!	ip(mp,5) = ip(n,5)
!	ip(mp,6) = ip(n,6)
!	ip(mp,7) = ip(n,7)
!	ip(mp,8) = ip(n,8)
!	ip(mp,9) = ip(n,9)
!	ip(mp,10) = ip(n,10)
!       ENDIF
!
!       IF( p(n, 1) .LT. 15.0 .AND. origloc( 1 ) >= 15.0 ) THEN
!               p( n, 1 ) = origloc( 1 )
!       ENDIF

!
! now we contribute to the concentration at the reporting time
!
      nstepav = nstepav + passes
      
      IF (concprint >= 1) THEN
        IF (tempavg == 1) THEN
          inbounds = 1
          DO  l = 1, numax
            ploc(l) = IDINT(p(n,l)/delc(l))  + 1
! ploc test
!			  ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
!            IF((ploc(l) <= 1).OR.(ploc(l) >= domax(l) )) THEN
            IF((ploc(l) <= 0).OR.(ploc(l) > domax(l) )) THEN
              inbounds = 0
            END IF
          END DO
          
          IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
            c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
                SNGL( dble(iP(n,2))*p(n,4)/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )
          END IF
        ELSE   ! not temporal averaging
          inbounds = 1
          DO  l = 1, numax
            ploc(l) = IDINT(p(n,l)/delc(l)) + 1
! change ploc test
!			ploc(l) = INT((p(n,l)+del2v(l))/delv(l))
!            IF((ploc(l) <= 1).OR.(ploc(l) >= domax(l) )) THEN
            IF((ploc(l) < 1).OR.(ploc(l) > domax(l) )) THEN
              inbounds = 0
            END IF
          END DO

		  if(P(n,7) > tnext) inbounds = 0
          
          IF (inbounds == 1) THEN
! we are adding a fix for unsat retardation
rstar = 1.d0 + (Rtard(ip(n,1),ploc(1),ploc(2),ploc(3))-1.d0)/sat(ploc(1),ploc(2),ploc(3))
            c(ip(n,1),ploc(1),ploc(2),ploc(3)) = c(ip(n,1),ploc(1),ploc(2),ploc(3)) +  &
                SNGL( dble(iP(n,2))*p(n,4)/(cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3))*Rstar) )
!				uncomment lines for major amounts of debugging info
!				print*,ip(n,1),ploc(1),ploc(2),ploc(3)
!               print*,p(n,4),cellv,porosity(ploc(1),ploc(2),ploc(3))
!				print*,c(ip(n,1),ploc(1),ploc(2),ploc(3))
          END IF
        END IF ! temporal averaging?
      END IF   ! print concentrations?
      IF (partprint == 1) THEN
        IF ((p(n,1) > 0.).AND.(p(n,2) > 0.).AND.(p(n,3) > 0.)) then	
                       pdist = 0.0d0
                       do ii = 1, 3
                         pdist = pdist +  (p(n,ii)-lastprint(n,ii))**2 
                       end do
                       pdist = dsqrt(pdist)
                       if (pdist >= minpdist) then
                       !WRITE(13,61) p(n,1), p(n,2), p(n,3), p(n,5), p(n,7)
                       lastprint(n,1) = p(n,1)
                       lastprint(n,2) = p(n,2)
                       lastprint(n,3) = p(n,3)
         end if ! if pdist
       end if  ! if bounds
! check if we've written out this particle yet and check if we hit land surface or WT
! if not we write, then set iprP to 2
  if (iprP(n,1) == 1) then      
  write(113,61) P(n,1), P(n,2), P(n,3),P(n,5), P(n,7)
  iprP(n,1) = 2
  end if
      END IF

      IF(partprint == 2) THEN

! check if we've written out this particle yet
! if not we write, then set iprP to 2
  if (iprP(n,1) == 1) then      
  write(113,61) P(n,1), P(n,2), P(n,3),P(n,5), P(n,7)
  iprP(n,1) = 2
  end if

        !WRITE(13,62) p(n,1), p(n,2), p(n,3), p(n,7)				&
        ! ,p(n,8),p(n,9),p(n,10)
      END IF
      
! endif
      
!      endif  ! whether or not to interrogate paticle position for conc
      
      
END DO ! the big particle loop:  DO  n = 1, np
    time =  tnext
    nstepav = nstepav/np
    WRITE(666,'(" Average number of Particle Steps:",e12.4)') nstepav
	WRITE(666,'(" Current time:",e12.4," Step:",e12.4)') time,tnext
    WRITE(666,*)
! assign current timestamp to mass array
    mass(1,it+1) =  tnext
!
! Clear out old moments
!
    rxx = 0.0D0
    ryy = 0.0D0
    sxx = 0.0D0
    syy = 0.0D0
    sxy = 0.0D0
    szz = 0.0D0
    rzz = 0.0D0

!
! now do the chemical reactions
!
  CALL countParticles( p,ip, np, xtent, ytent, ztent, n_constituents, delv, &
       tnext )

  CALL ConcToMgperLiter( c, xtent, ytent, ztent, n_constituents )

  CALL    checkConc( p, cellv, c, porosity, sat, xtent, ytent, ztent,      &
                                                  n_constituents, modelname )  

  CALL maxMassAllCellSpec( p, np, ip )

 CALL addRemoveParticles( p, ip, ipwell, irp, iprp, lastprint, &
                                    np, npmax, delv,  n_constituents, &
                                xtent, ytent, ztent, c, ( tnext - t_prev ), &
                           porosity, sat, modelname )
!  np = mp
  CALL updateConc( c, p, ip, np, delc, domax, sat, porosity, Rtard, tnext )

  CALL ConcAddBackground( c, xtent, ytent, ztent, n_constituents, modelname )
!
! now write out concentration at given time
!
    IF (concprint == 1) THEN
      WRITE (dotit,199) it
      199 FORMAT('.',i5.5)
      DO i = 1, n_constituents
      confile1 = trim(confile(i)) // dotit
      filename= trim(confile(i))//dotit//confile2
      
      CALL cbin_write(c(i,:,:,:),filename,xtent, ytent,ztent,delv(1),delv(2),delv(3))
if (part_conc_write==1) then
!	  write out each particle in each conc cell
!
!
 filename = trim(confile(i)) // dotit //'.txt'
open (1919,file=trim(filename))
DO  k=1,ztent
  DO  j=1,ytent
    DO  ii=1,xtent
      IF (c(1,ii,j,k) > 0.d0) THEN
	  write(1919,*) '+++++++++++++++++++'
	  write(1919,*) ii,j,k
	  write(1919,*) C(1,ii,j,k)
		do n = 1, np
		ploc(1) = IDINT(p(n,1)/delv(1)) + 1
		ploc(2) = IDINT(p(n,2)/delv(2)) + 1
		ploc(3) = IDINT(p(n,3)/delv(3)) + 1
		if((ploc(1) == ii).and.(ploc(2)==j).and.(ploc(3)==k)) then
		write(1919,*) ii,j,k
		write(1919,*) n,p(n,7)
		write(1919,*) p(n,4),porosity(i,j,k)
		write(1919,*) p(n,1),p(n,2),p(n,3)
		write(1919,*)
		end if
		end do
      END IF
    END DO
  END DO
END DO
close(1919)
end if ! part_conc_write
!
!
!
      END DO
        else if (concprint == 2 ) THEN
              WRITE (dotit,199) it

  !CALL vtk_write(time, c(:,:,:,:),confile(:),xtent, ytent,ztent,delv(1),delv(2),delv(3),it,n_constituents,vtk_file)
  CALL gnuplot_write(c(:,:,:,:),xtent, ytent,ztent,it,n_constituents,npars%backgroundConc,vtk_file)
      
    END IF   !print concentration

    IF( velprint .EQ. 1) THEN
    
      filename = 'velocity' // dotit //'.cnb'
      CALL cbin_write_real8(v(3,:,:,:),filename,xtent, ytent,ztent,delv(1),delv(2),delv(3))

    END IF   !print concentration

    npact =0
    bdyx1 = 2.*delv(1)
    bdyx2 = DBLE(domax(1)-1)*delv(1)
    bdyy1 = 2.*delv(2)
    bdyy2 = DBLE(domax(2)-1)*delv(2)
    bdyz1 = 2.*delv(3)
    bdyz2 = DBLE(domax(3)-1)*delv(3)
    
    
    DO  n = 1, np
!C     write(13,61) P(n,1), P(n,2), P(n,3), P(n,5), P(n,7)
!   mass(iP(n,1)+1,it+1) = mass(iP(n,1)+1,it+1) + P(n,4)

      IF (p(n,1) > bdyx1) THEN
        IF (p(n,1) < bdyx2) THEN
          IF (p(n,2) > bdyy1) THEN
            IF (p(n,2) < bdyy2) THEN
              IF(p(n,3) > bdyz1) THEN
                IF(p(n,3) < bdyz2) THEN
    mass(iP(n,1)+1,it+1) = mass(iP(n,1)+1,it+1) + P(n,4)
	              npact = npact +1
                  rxx = rxx + p(n,1)
                  ryy = ryy + p(n,2)
                  rzz = rzz + p(n,3)
                  
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
      
    END DO
    IF (npact > 0) THEN
      rxx = rxx/npact
      ryy = ryy/npact
      rzz = rzz/npact
    ELSE
      rxx = 0.
      ryy = 0.
      rzz = 0.
    END IF
    DO  n = 1, np
      IF (p(n,1) > bdyx1) THEN
        IF (p(n,1) < bdyx2) THEN
          IF (p(n,2) > bdyy1) THEN
            IF (p(n,2) < bdyy2) THEN
              IF(p(n,3) > bdyz1) THEN
                IF(p(n,3) < bdyz2) THEN
                  
                  
                  sxx = sxx + (p(n,1)-rxx)**2
                  sxy = sxy + (p(n,1)-rxx)*(p(n,2)-ryy)
                  syy = syy + (p(n,2)-ryy)**2
                  szz = szz + (p(n,3)-rzz)**2
                  
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
      
      
    END DO
    IF (npact > 0) THEN
      sxx = sxx/npact
      sxy = sxy/npact
      syy = syy/npact
      szz = szz/npact
    ELSE
      sxx = 0.
      sxy = 0.
      syy = 0.
      szz = 0.
    END IF
    
    IF (momprint == 1) THEN
      WRITE(wellfile+10,188) time,npact,rxx,ryy,rzz,sxx,syy,szz
      188  FORMAT (7(e15.5,','),e15.5)
    END IF   ! print moments?

!
! new particle splitting subroutine
!

    IF (partsplit == 1) THEN
      npcurrent = np
      DO n = 1, npcurrent
        inbounds = 1
        DO  l = 1, numax
          ploc(l) = INT(p(n,l)/delc(l)) + 1
!          IF((ploc(l) <= 1).OR.(ploc(l) >= domax(l) )) THEN
          IF((ploc(l) < 1).OR.(ploc(l) > domax(l) )) THEN
            inbounds = 0
          END IF ! in bounds?
        END DO ! in bounds?

!
! if we're in bounds, check for single particles
!

        IF (inbounds == 1) THEN
          cps = p(n,4) / (cellv*sat(ploc(1),ploc(2),ploc(3))*porosity(ploc(1),ploc(2),ploc(3)))
          IF (c(ip(n,1),ploc(1),ploc(2),ploc(3)) <= cps) THEN
            IF (c(ip(n,1),ploc(1),ploc(2),ploc(3)) >= cmin) THEN
              IF (np < npmax) THEN
                np = np + 1
                p(np,1) = p(n,1) - delv(1)/100.
                p(np,2) = p(n,2) - delv(2)/100.
                p(np,3) = p(n,3) - delv(3)/100.
                p(np,4) = p(n,4) / 2.d0
                p(np,5) = DBLE(np)
                p(np,7) = p(n,7)
				ip(np,1) = ip(n,1)
				ip(np,2) = ip(n,2)
				ip(np,3) = ip(n,3)
				ip(np,4) = ip(n,4)
				ip(np,5) = ip(n,5)
				ip(np,6) = ip(n,6)
				ip(np,7) = ip(n,7)
				ip(np,8) = ip(n,8)
				ip(np,9) = ip(n,9)
				ip(np,10) = ip(n,10)
                p(n,1) = p(n,1) + delv(1)/100.
                p(n,2) = p(n,2) + delv(2)/100.
                p(n,3) = p(n,3) + delv(3)/100.
                p(n,4) = p(n,4) / 2.d0
              END IF  ! used up all our particles?
            END IF ! bigger than min conc?
          END IF  ! need to split?
        END IF ! are we in bounds?
      END DO ! particle split
    END IF ! part split ?

!
! Next timestep
!
!@RMM
!print*, iP(1:500,6)
    WRITE(666,*) 'Current number of particles: ', np
    flush(666)

! now write out concentration at given time
!
  END DO

CALL deallocateParticles_Memory( xtent, ytent, ztent, n_constituents )

 deallocate( ic_time_begin, ic_time_end, ic_mass_or_conc )

!
! write out all constituent masses over time
!
    WRITE(666,*)
	WRITE(666,*) 'Total Mass (C+S) by Constituent'
	WRITE(666,*) 'Step	  time    ',(trim(confile(i))//'   ',i=1,n_constituents)
    DO i = 1, nt + 1
!	  WRITE(666,'(i4,1x,<n_constituents>e12.4,e12.4)') i-1,mass(1:n_constituents+1,i)
      WRITE(format_desc,'("(i4,1x,",i2,"(e10.3))")') n_constituents+1
	  WRITE(666,format_desc) i-1,mass(1:n_constituents+1,i)
!321   FORMAT(i4,',',<n_constituents>(e12.4,','),e12.4)
	END DO

! write out amount of mass recycled
if (recycle_well == 1) then
mass_rec = 0

do n =1, np
mass_rec(irP(n) + 1) = mass_rec(irP(n) + 1) + P(n,4)
end do ! n

write(666,*)
write(666,*) ' Mass of particles recycled '
write(666,*) ' no,        total mass'
do ii = 1, 50
write(666,'(i4,e12.4)') ii-1, mass_rec(ii)
end do
end if
!
! Write out Well BTC
!

  IF (wellprint >= 1) THEN
  DO k = 1, n_constituents
    DO  i = 1, nw
      DO  j = 1, welltnumb
        IF (well(i,5) /= 0) THEN
          WRITE(wellfile+k-1,930) i,DBLE(j)*welltavg(i),		&
             ((btc(i,j,k))/welltavg(i))/well(i,5)
          930   FORMAT (i4,f10.1,e15.5)
        ELSE
          WRITE(wellfile+k-1,930) i,DBLE(j)*welltavg(i),btc(i,j,k)
        END IF
      END DO ! j, num well num
    END DO ! i, num wells
  END DO ! k, num constits
  END IF   ! write out well BTC?
  !CLOSE(13)
  
!
! Write out Plane BTC
!
 
  IF (nplane >= 1) THEN
  if(massive_debug == 1) then
  print*, ' writing plane'
  print*, n_constituents
  print*, nplane
  print*, welltavg(1)
  print*, welltnumb
  end if
  DO k = 1, n_constituents
    OPEN (33,FILE=TRIM(planefile(k)))
    DO i = 1, nplane
      DO j = 1, welltnumb
        931   FORMAT (i4,f10.1,e15.5)
        WRITE(33,931) i,DBLE(j)*welltavg(1),btp(i,j,k)
      END DO  ! i, nplane
    END DO  ! j, num wells
  END DO ! k num constit
  END IF   ! write out well BTC?
  CLOSE(33)
  close (16)
  close(17)
  
  RETURN
END SUBROUTINE slimfast



SUBROUTINE compact_vread(x,nx,ny,nz,time)

INTEGER, INTENT(IN)                      :: nx
INTEGER, INTENT(IN)                      :: ny
INTEGER, INTENT(IN)                      :: nz

REAL*8, INTENT(IN OUT)                   :: x(3,nx,ny,nz)
!CHARACTER (LEN=100), INTENT(IN OUT)      :: filename

REAL*8, INTENT(IN OUT)                   :: time
!INTEGER*4, INTENT(IN OUT)                :: ts


REAL*8  timediscard, vdiscard
INTEGER*4  ii,jj, k, j, i 

!
! open file
!

!OPEN (15,FILE=trim(filename),FORM='binary',STATUS='old',readonly)

!
! skip nodes as needed
!

!DO ii = 1, ts
!  READ(15) timediscard
!  DO jj = 1,  (nx*ny*nz)
!    READ(15)  vdiscard, vdiscard, vdiscard   ! data
!  END DO !jj
!END DO !ii

READ(16) time
DO k = 1,  nz
  DO j = 1, ny
    DO i = 1, nx
      READ(16)  x(1,i,j,k),x(2,i,j,k),x(3,i,j,k)
!      x(1,i,j,k) =  x(1,i,j,k)*86400.d0
!      x(2,i,j,k) =  x(2,i,j,k)*86400.d0
!      x(3,i,j,k) =  x(3,i,j,k)*86400.d0
!	  x(1,i,j,k) =  x(1,i,j,k)/997.d0
!      x(2,i,j,k) =  x(2,i,j,k)/997.d0
!      x(3,i,j,k) =  x(3,i,j,k)/997.d0
    END DO !i
  END DO !j
END DO !k
!CLOSE(15)
!ts = ts + 1
RETURN
END SUBROUTINE compact_vread

! Moved to Particles.f90
! 
!FUNCTION ran1(idum)
!INTEGER*4 idum,IA,IM,IQ,IR,NTAB,NDIV
!REAL*8 ran1,AM,EPS,RNMX
!PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,  &
!   NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!      INTEGER j,k,iv(NTAB),iy
!      SAVE iv,iy
!      DATA iv /NTAB*0/, iy /0/
!      if (idum.le.0.or.iy.eq.0) then
!        idum=max(-idum,1)
!        do 11 j=NTAB+8,1,-1
!          k=idum/IQ
!          idum=IA*(idum-k*IQ)-IR*k
!          if (idum.lt.0) idum=idum+IM
!          if (j.le.NTAB) iv(j)=idum
!11      continue
!        iy=iv(1)
!      endif
!      k=idum/IQ
!      idum=IA*(idum-k*IQ)-IR*k
!      if (idum.lt.0) idum=idum+IM
!
!      j=1+iy/NDIV
!      iy=iv(j)
!      iv(j)=idum
!      ran1=min(AM*iy,RNMX)
!      return
!      END

!
!-----------------------------------------------------------------------
!     function to generate pseudo random numbers, uniform (0,1)
!-----------------------------------------------------------------------
!
      function rand2(iuu)
!
	real*8 rand2,rssq,randt
	integer m,l,iuu
   
!
      data m/1048576/
      data l/1027/
      n=l*iuu
      iuu=mod(n,m)
      rand2=float(iuu)/float(m)
!	rand2 = rand(0)
	iiu = iiu + 2
      return
      end function
END MODULE SLIM_FAST
