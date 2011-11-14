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

subroutine v_calc_const_sat(v,hkx,hky,hkz,vga,vgn,headfile,phi,delv,nx,ny,nz,press,sat) 

real*8 :: delv(3),heads, m,temp1,temp2,temp3,hstaru,hstard
integer*4 i,j,k, nx,ny,nz,press,udir
real*8  :: phi(:,:,:)
!real*8  phi(nx,ny,nz)
real*8  :: v(:,:,:,:),hkx(:,:,:),hky(:,:,:),hkz(:,:,:),sat(:,:,:),vga(:,:,:),vgn(:,:,:)
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
V = 0.d0

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
                    kr1 = ( (1.0d0- (temp1/((1.0d0+temp2)**m))**2.0d0) )/temp3 
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
                      kr = ( (1.0d0- (temp1/((1.0d0+temp2)**m))**2.0d0) )/temp3 
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
                (DSQRT( kr1 * kr ) / ( theta * delv(1) * scx(i) ) ) *    & 
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
                    kr1 = ( (1.0d0- (temp1/((1.0d0+temp2)**m))**2.0d0) )/temp3 
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
                       kr = ((1.0d0-(temp1/((1.0d0+temp2)**m))**2.0d0) )/temp3 
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
                ( DSQRT( kr1 * kr ) /( theta * delv(2) * scy(j) ) )*     &
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
                    kr1 = ( (1.0d0- (temp1/((1.0d0+temp2)**m))**2.0d0) )/temp3 
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
                     kr = ((1.0d0-(temp1/((1.0d0+temp2)**m))**2.0d0) )/temp3 
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
         
            V(3,i,j,k) = -(Hstaru - Hstard + 0.01d0 )*        &
             ( DSQRT( kr1 * kr ) /( theta * delv(3) * scz(k) ) ) *   &
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

end subroutine v_calc_const_sat
