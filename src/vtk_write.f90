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

SUBROUTINE vtk_write(time,x,conc_header,ixlim,iylim,izlim,dx,dy,dz,icycle,n_constituents,vtk_file)
real*8                 :: time
REAL*4    :: x(:,:,:,:)
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



INTEGER*4 i,j,k, ijk, debug,l
CHARACTER*1 lf
CHARACTER*12 num1, num2, num3
CHARACTER*8 ctime
real*8 number

debug = 0
!
!      Open File
!
Write(ctime,'(i8.8)') icycle
OPEN(15,FILE=trim(vtk_file)//'.'//ctime//'.vtk',FORM='unformatted',  &
    access='stream',convert='BIG_ENDIAN')



!
!      Write header info
!
lf = char(10) ! line feed character
Write(15) "# vtk DataFile Version 3.0"//lf
Write(15) "Slim Concentration Output"//lf
Write(15) "BINARY"//lf
Write(15) "DATASET RECTILINEAR_GRID"//lf
Write(15) "FIELD FieldData 2"//lf
write(15) "TIME 1 1 double"//lf
write(15) time
write(15) lf
Write(15) "CYCLE 1 1 int"//lf
Write(15) icycle
Write(15) lf

Write(num1,'(i12)') ixlim+1
Write(num2,'(i12)') iylim+1
Write(num3,'(i12)') izlim+1

write(15) "DIMENSIONS "//num1//" "//num2//" "//num3//lf

Write(15) "X_COORDINATES "//num1//" double"//lf
do i = 0, ixlim
number = float(i)*dx 
Write(15) number
end do
write(15) lf

Write(15) "Y_COORDINATES "//num2//" double"//lf
do i = 0, iylim
number = float(i)*dy 
Write(15) number 
end do
Write(15) lf 

Write(15) "Z_COORDINATES "//num3//" double"//lf
do i = 0, izlim
number = float(i)*dz 
Write(15) number
end do
write(15) lf
Write(num1,'(i12)') ixlim*iylim*izlim

Write(15) "CELL_DATA "//num1//lf

do l = 1, n_constituents
write(15) "SCALARS "//trim(conc_header(l))//" float 1"//lf
Write(15) "LOOKUP_TABLE default"//lf

DO  k=1,izlim
  DO  j=1,iylim
    DO  i=1,ixlim
        WRITE(15) x(l,i,j,k)
    END DO
  END DO
END DO
write(15) lf
end do ! l, n_constits
CLOSE(15)

if(debug == 1) then
! asci cbin for debugging

OPEN(15,FILE=trim(vtk_file)//'.ascii.vtk')



CLOSE(15)
end if ! debug

RETURN
END SUBROUTINE vtk_write

