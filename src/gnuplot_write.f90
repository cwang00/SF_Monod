! Copyright 2010 Reed M. Maxwell
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

SUBROUTINE gnuplot_write(x,ixlim,iylim,izlim,icycle,n_constituents,bg, vtk_file)
REAL*4    :: x(:,:,:,:)
REAL, DIMENSION(:) :: bg
INTEGER*4 :: ixlim
INTEGER*4 :: iylim
INTEGER*4 :: izlim
INTEGER                :: icycle
INTEGER*4              :: n_constituents
CHARACTER (LEN=100)      :: vtk_file



INTEGER*4 i,j,k, l
CHARACTER*8 ctime
CHARACTER*12 cfmat

!
!      Open File
!
write(cfmat, '(A,I6,5A)') '(', ixlim, 'F9.3)'
Write(ctime,'(i8.8)') icycle
OPEN(15,FILE=trim(vtk_file)//'.'//ctime//'.gnuplot')



do l = 1, n_constituents

WRITE(15,*) '#SPECIES: ', l
DO  j=1,iylim
WRITE(15,*) '#Y =  ', j
WRITE(15,*)
DO  k=izlim, 1, -1
        WRITE(15, (cfmat)) ( x(l,i,j,k)+bg(l), i = 1, ixlim )
END DO
WRITE(15,*)
END DO
end do ! l, n_constits
CLOSE(15)

RETURN
END SUBROUTINE gnuplot_write

