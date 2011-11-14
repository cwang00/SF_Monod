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

SUBROUTINE cbin_write_real8(x,filename,ixlim,iylim,izlim,dx,dy,dz)
REAL*8    :: x(:,:,:)
REAL*4    :: xtemp
CHARACTER (LEN=40)     :: filename
INTEGER*4 :: ixlim
INTEGER*4 :: iylim
INTEGER*4 :: izlim
REAL*8                 :: dx
REAL*8                 :: dy
REAL*8                 :: dz

!INTEGER*4 :: num_constit
!INTEGER*4 :: constit_num


REAL*4                 :: dx1
REAL*4                 :: dy1
REAL*4                 :: dz1

INTEGER*4 i,j,k, ijk, debug

!allocate (x(num_constit,ixlim,iylim,izlim))
debug = 0
!
!      Open File
!

OPEN(15,FILE=trim(filename),FORM='unformatted',  &
    access='stream', status='replace',convert='BIG_ENDIAN')

!
!      Calc domain bounds
!

!
!      Write header info
!


WRITE(15) ixlim
WRITE(15) iylim
WRITE(15) izlim

dx1 = dble(dx)
dy1 = dble(dy)
dz1 = dble(dz)

WRITE(15) dx1
WRITE(15) dy1
WRITE(15) dz1


! print*, ixlim, iylim, izlim
! print*, dx, dy, dz

DO  k=1,izlim
  DO  j=1,iylim
    DO  i=1,ixlim
      ijk = i + (j-1)*ixlim + (k-1)*ixlim*iylim
      IF (x(i,j,k) .NE. 0.d0) THEN
        xtemp = x(i,j,k)
        WRITE(15) ijk, xtemp
! print*,  i,j,k, x(constit_num,i,j,k)
      END IF
    END DO
  END DO
END DO

CLOSE(15)

RETURN
END SUBROUTINE cbin_write_real8

