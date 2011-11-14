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

subroutine v_calc(v,hkx,hky,hkz,vga,vgn,sres,headfile,phi,delv,nx,ny,nz,press,sat)  !,vgn,vga,sres,ssat)

real*8 :: delv(3),heads, m,ssat,temp1,temp2,temp3,hstaru,hstard
integer*4 i,j,k, nx,ny,nz,press,udir
real*8  :: phi(:,:,:)
!real*8  phi(nx,ny,nz)
real*8  :: v(:,:,:,:),hkx(:,:,:),hky(:,:,:),hkz(:,:,:),sat(:,:,:),vga(:,:,:),vgn(:,:,:),sres(:,:,:)
!real*8   v(3,nx,ny,nz)
character(100) :: headfile
real*8 x0,y0,z0, kr, kr1, theta
real*8, allocatable :: head(:,:,:),scx(:),scy(:),scz(:)  !,kr(:,:,:)


interface
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

end interface


allocate(head(nx,ny,nz),scx(nx),scy(ny),scz(nz)) !,kr(nx,ny,nz)) !,  &
!phi(nx,ny,nz),v(3,nx,ny,nz))

scx = 1.0d0
scy = 1.0d0
scz = 1.0d0

! variable dx,dy,dz HARD WIRED into velocity right now, really should
! make this more permanent later...
!scx(1:11) = 75.0d0
!scx(562:nx) = 75.0d0

!scy(1:9) = 50.0d0
!scy(160:ny) = 50.0d0

!scz(1:10) = 25.0d0

!call pf_read(hkx(:,:,:),kxfile,nx,ny,nz,delv(1),delv(2),delv(3))
!call pf_read(hky(:,:,:),kyfile,nx,ny,nz,delv(1),delv(2),delv(3))
!call pf_read(hkz(:,:,:),kzfile,nx,ny,nz,delv(1),delv(2),delv(3))

call pf_read(head,headfile,nx,ny,nz,delv(1),delv(2),delv(3))


!print*, allocated(x)
!vga = 1.
!vgn = 2.
ssat = 1.
!sres = 0.5
V = 0.d0
DO k = 1,  nz
  DO j = 1, ny
    DO i = 1, nx

! we calculate Sat and relperm everywhere, later will read in
! we also pass Sat BACK to the main particle routine (but not relperm)
! right now this goofs up the initial concentrations since sat is init'd to 1 everywhere
! even though for IC's it might be something else...

                   if (Head(i,j,k) >= 0.0d0) then
                       sat(i,j,k) = 1.0d0
!					   kr(i,j,k) = 1.0d0
                   else
                        m = 1.0d0 - (1.0d0/vgn(i,j,k))
                       heads = dabs(Head(i,j,k))
                       temp1 = (vga(i,j,k)*heads)**vgn(i,j,k)
                       temp2 = (1.0d0+temp1)**m
                       sat(i,j,k) = Sres(i,j,k) + (Ssat - Sres(i,j,k))/temp2
                          !heads = dabs(Head(i+udir,j,k))
 !                      temp1 = (vga*heads)**(vgn-1)
!                       temp2 = (vga*heads)**vgn
!                       temp3 = (1.+temp2)**(m/2.)
!                       kr(i,j,k) = ( (1.- (temp1/((1+temp2)**m))**2) )/temp3 
                       end if
                       end do
         end do
end do

DO k = 1,  nz
  DO j = 1, ny
    DO i = 1, nx

         if (i > 1)       then

                 IF ( Head( i - 1, j, k ) >= 0.0d0 ) THEN
                       kr1 = 1.0d0
                 ELSE

                        m = 1.0d0 - ( 1.0d0 / vgn( i - 1, j, k ) )
                        heads = dabs(Head(i - 1, j, k))
                       temp1 = (vga(i - 1,j,k)*heads)**(vgn(i - 1,j,k)-1.0d0)
                      temp2 = (vga(i - 1,j,k)*heads)**vgn(i - 1,j,k)
                       temp3 = (1.0d0+temp2)**(m/2.0d0)
                    kr1 = ( 1.0d0- temp1/(1.0d0+temp2)**m)**2.0d0/temp3 
                 ENDIF
! do upwinding check
!                      udir = -1
!          if (Head(i,j,k) > Head(i-1,j,k) ) udir = 0 
                   if (Head(i,j,k) >= 0.0d0) then
                       kr = 1.0d0 
                   else
                        m = 1.0d0 - (1.0d0/vgn(i,j,k))
                       heads = dabs(Head(i,j,k))
                       temp1 = (vga(i,j,k)*heads)**(vgn(i,j,k)-1.0d0)
                      temp2 = (vga(i,j,k)*heads)**vgn(i,j,k)
                       temp3 = (1.0d0+temp2)**(m/2.0d0)
                    kr = ( 1.0d0- temp1/(1.0d0+temp2)**m)**2.0d0/temp3 
                   end if 
                  theta = 0.5 * ( phi( i - 1, j, k ) * sat( i - 1, j, k ) + &
                                  phi( i, j, k ) * sat( i, j, k ) )
! convert press to head if need be
Hstaru = Head(i,j,k)
Hstard = Head(i-1,j,k)
if (press == 1) then 
Hstaru = Head(i,j,k) + dfloat(k-1)*delv(3) + delv(3)/2.d0 !+z0
Hstard = Head(i-1,j,k) + dfloat(k-1)*delv(3) + delv(3)/2.d0 !+z0
end if
                 ! V(1,i,j,k) = -(Head(i,j,k) - Head(i-1,j,k))*     &
                      V(1,i,j,k) = -(Hstaru - Hstard)*     &
                (DSQRT( kr1 * kr ) /(theta*delv(1)*scx(i))) *    & 
                 (2.0d0/(1.0d0/hKx(i,j,k) + 1.0d0/hKx(i-1,j,k) ))
                         ! (hKx(i+udir,j,k)*kr(i+udir,j,k)) 
                   ! if ((i==10).and.(j==10)) print*, i,j,k,V(1,i,j,k),sat,kr,hkx(i,j+udir,k),head(i,j,k),head(i+1,j,k),udir
            end if

         if (j > 1)  then
                 IF ( Head( i, j - 1, k ) >= 0.0d0 ) THEN
                       kr1 = 1.0d0
                 ELSE

                        m = 1.0d0 - ( 1.0d0 / vgn( i, j - 1, k ) )
                        heads = dabs(Head(i, j - 1, k))
                       temp1 = (vga(i,j - 1, k)*heads)**(vgn(i,j - 1, k)-1.0d0)
                      temp2 = (vga(i,j - 1, k)*heads)**vgn(i,j - 1, k)
                       temp3 = (1.0d0+temp2)**(m/2.0d0)
                    kr1 = ( 1.0d0- temp1/(1.0d0+temp2)**m)**2.0d0/temp3 
                 ENDIF
! do upwinding check
!                      udir = -1
!          if (Head(i,j,k) > Head(i,j-1,k) ) udir = 0 
                   if (Head(i,j,k) >= 0.0d0) then
                       kr = 1.0d0
                   else
                        m = 1.0d0 - (1.0d0/vgn(i,j,k))
                       heads = dabs(Head(i,j,k))
                       temp1 = (vga(i,j,k)*heads)**(vgn(i,j,k)-1.0d0)
                       temp2 = (vga(i,j,k)*heads)**vgn(i,j,k)
                      temp3 = (1.0d0+temp2)**(m/2.0d0)
                    kr = ( 1.0d0- temp1/(1.0d0+temp2)**m)**2.0d0/temp3 
                   end if          

                  theta = 0.5 * ( phi( i, j - 1, k ) * sat( i, j - 1, k ) + &
                                  phi( i, j, k ) * sat( i, j, k ) )
! convert press to head if need be
Hstaru = Head(i,j,k)
Hstard = Head(i,j-1,k)
if (press == 1) then 
Hstaru = Head(i,j,k) + dfloat(k-1)*delv(3) + delv(3)/2.d0 !+z0
Hstard = Head(i,j-1,k) + dfloat(k-1)*delv(3) + delv(3)/2.d0 !+z0
end if


           ! V(2,i,j,k) = -(Head(i,j,k) - Head(i,j-1,k))*        &
                  V(2,i,j,k) = -(Hstaru - Hstard)*        &
                (DSQRT( kr1 * kr ) /( theta *delv(2)*scy(j)))*     &
                 (2.0d0/(1.0d0/hKy(i,j,k) + 1.0d0/hKy(i,j-1,k) ))
                   !(hKy(i,j+udir,k)*kr(i,j+udir,k))
                   
                   end if
         if (k > 1)  then

                 IF ( Head( i, j, k - 1 ) >= 0.0d0 ) THEN
                       kr1 = 1.0d0
                 ELSE

                        m = 1.0d0 - ( 1.0d0 / vgn( i, j, k - 1 ) )
                        heads = dabs(Head(i, j, k - 1))
                       temp1 = (vga(i,j, k - 1)*heads)**(vgn(i,j, k - 1)-1.0d0)
                      temp2 = (vga(i,j, k - 1)*heads)**vgn(i,j, k - 1)
                       temp3 = (1.0d0+temp2)**(m/2.0d0)
                    kr1 = ( 1.0d0- temp1/(1.0d0+temp2)**m)**2.0d0/temp3 
                 ENDIF

! do upwinding check
!                      udir = -1
!          if (Head(i,j,k) > Head(i,j,k-1) ) udir = 0 
                   if (Head(i,j,k) >= 0.0d0) then
                       kr = 1.0d0
                   else
                        m = 1.0d0 - (1.0d0/vgn(i,j,k))
                     heads = dabs(Head(i,j,k))
                     temp1 = (vga(i,j,k)*heads)**(vgn(i,j,k)-1.0d0)
                     temp2 = (vga(i,j,k)*heads)**vgn(i,j,k)
                     temp3 = (1.0d0+temp2)**(m/2.0d0)
                    kr = ( 1.0d0- temp1/(1.0d0+temp2)**m)**2.0d0/temp3 
                   end if           

                  theta = 0.5 * ( phi( i, j, k - 1 ) * sat( i, j, k - 1 ) + &
                                  phi( i, j, k ) * sat( i, j, k ) )
! convert press to head if need be
Hstaru = Head(i,j,k)
Hstard = Head(i,j,k-1)
if (press == 1) then 
Hstaru = Head(i,j,k) + dfloat(k-1)*delv(3) + delv(3)/2.d0 !+z0
Hstard = Head(i,j,k-1) + dfloat(k-2)*delv(3) + delv(3)/2.d0 !+z0
end if
         
            V(3,i,j,k) = -(Hstaru - Hstard )*        &
             (DSQRT( kr1 * kr ) /(theta * delv(3)* scz(k) ) )*   &
                 (2.0d0/(1.0d0/hKz(i,j,k) + 1.0d0/hKz(i,j,k-1) ))
                   !(hKz(i,j,k+udir)*kr(i,j,k+udir))
                   end if

    END DO !i
  END DO !j
END DO !k

deallocate(head,scx,scy,scz) !,kr)

!print*, allocated(hkx)
!print*, allocated(hky)
!print*, allocated(hkz)
!print*, allocated(head)

return

end subroutine v_calc

subroutine pf_read(x,filename,nx,ny,nz,dx2,dy2,dz2)
real*8  :: x(:,:,:)
!real*8  :: x(nx,ny,nz)
character*100 :: filename
integer*4 :: nx
integer*4 :: ny
integer*4 :: nz
real*8  :: dx2
real*8  :: dy2
real*8  :: dz2
real*8  :: x1
real*8  :: y1
real*8  :: z1
Real*8 value,  ri, rj, rk1, rk2, headsum, rsum, junk,  &
ksum, kavg,f, dx, dy, dz								
Integer*4 i,j,k, nni, nnj, nnk, ix, iy, iz,			&
ns,  rx, ry, rz, nnx, nny, nnz, is
integer*4 ijk, namelength, xtent,ytent,ztent

!allocate(x(nx,ny,nz))

       
	print*, trim(filename)

!      Open File

open(15,file=trim(filename),form='unformatted',   &
access='stream',convert='BIG_ENDIAN',status='old')

!      Read in header info

        read(15) x1
        read(15) y1
        read(15) z1

        read(15) nx
        read(15) ny
        read(15) nz

        read(15) dx
        read(15) dy
        read(15) dz

	dx2 = dx
	dy2 = dy
	dz2 = dz

        read(15) ns

print*, 'header for file:', filename
print*, 'X0, Y0, Z0: ', x1, y1, z1
print*, 'nx, ny, nz: ', nx, ny, nz
print*, 'dx, dy, dz: ', dx, dy, dz
print*, 'number of subgrids: ', ns

do is = 0, (ns-1)

  read(15) ix
  read(15) iy
  read(15) iz
  read(15) nnx
  read(15) nny
  read(15) nnz

  read(15) rx
  read(15) ry
  read(15) rz

  do  k=iz +1 , iz + nnz
    do  j=iy +1 , iy + nny
    do  i=ix +1 , ix + nnx
    read(15) value
	x(i,j,k) = value
   	end do
    end do
  end do

end do
      close(15)
      print*, 'file read in'
        return
        end
