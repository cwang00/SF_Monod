!	Matrix Diffusion lookup
!	R. Maxwell
!	12/20/00
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
!
subroutine mat_diff(nx, ny, nz,  mat_file,perm_cat_file,n_constituents,k_att,k_det,constit_num)

integer*4 nx,ny,nz, cat_table(100), mat_diff_cat_num,  i, j, icat_temp,  finder,&
 constit_num, iic, kacheck, kdcheck,n_constituents, k

real*8 cat_value(100,30),k_att(n_constituents,nx,ny,nz), &
       k_det(n_constituents,nx,ny,nz)

integer nnx,nny,nnz
integer*4,allocatable::perm_ind(:,:,:)

character*100 mat_file, perm_cat_file
character*25 k_att_check_file(10),k_det_check_file(10)

k_att_check_file(1) = 'kacheck.c1.asc.txt'
k_att_check_file(2) = 'kacheck.c2.asc.txt'
k_att_check_file(3) = 'kacheck.c3.asc.txt'
k_att_check_file(4) = 'kacheck.c4.asc.txt'
k_att_check_file(5) = 'kacheck.c5.asc.txt'

k_det_check_file(1) = 'kdcheck.c1.asc.txt'
k_det_check_file(2) = 'kdcheck.c2.asc.txt'
k_det_check_file(3) = 'kdcheck.c3.asc.txt'
k_det_check_file(4) = 'kdcheck.c4.asc.txt'
k_det_check_file(5) = 'kdcheck.c5.asc.txt'


allocate(perm_ind(nx,ny,nz))


print*, 'nx,ny,nz',nx,ny,nz
print*, ' constit ',constit_num
open (120, file=trim(mat_file))

read(120,*)
read(120,*)
read(120,*)
read(120,*)
read(120,*) mat_diff_cat_num
read(120,*)

! write out to log file 
write(666,*) ' Matrix Diffusion Lookup'
write(666,*)
write(666,*) '	reading permeability catagories from:'
write(666,*) trim(perm_cat_file)
write(666,*) ' num perm catagories :' , mat_diff_cat_num
write(666,*)
iic = constit_num
write(666,*) ' constituent number:',iic
write(666,*)
! read in lookup table
do i = 1 , mat_diff_cat_num
	read(120,*) cat_table(i), cat_value(i,1:2)
	write(666,*) ' catagory:',cat_table(i)
	write(666,*) cat_value(i,1:2)
end do



print*, ' perm ind'
! read in perm indicators
open(130,file=trim(perm_cat_file))
	read(130,*) nnx,nny,nnz
!	print*, nnx,nny,nnz
	do k=1,nz
!	print*, k
		do j=1,ny
			do i =1,nx
				read(130,*) perm_ind(i,j,k)
			end do ! x,i
		end do !y,j
	end do !z, k
close(130)

k_att(constit_num,1:nx,1:ny,1:nz) = 0.d0
k_det(constit_num,1:nx,1:ny,1:nz) = 0.d0



do i = 1, nx
	do j = 1, ny
		do k = 1, nz
			icat_temp = finder(perm_ind(i,j,k),cat_table,mat_diff_cat_num) 
			if (icat_temp.le.0.) print*, ' icat <0! ', icat_temp
			if (icat_temp.gt.mat_diff_cat_num) print*, ' icat > max cats ', icat_temp
			k_att(constit_num,i,j,k) = cat_value(icat_temp,1)
            k_det(constit_num,i,j,k) = cat_value(icat_temp,2)
!			print*,i,j,k,icat_temp,cat_value(icat_temp,1)
		end do
	end do
end do

kacheck = 1
kdcheck = 1

if (kacheck.eq.1) then
open (199,file=trim(k_att_check_file(constit_num)))
	write(199,*) nx,ny,nz
	do k=1,nz
		do j=1,ny
			do i =1,nx
				write(199,*) k_att(iic,i,j,k)
			end do ! x,i
		end do !y,j
	end do !z, k
close(199)
end if

if (kdcheck.eq.1) then
open (199,file=trim(k_det_check_file(constit_num)))
	write(199,*) nx,ny,nz
	do k=1,nz
		do j=1,ny
			do i =1,nx
				write(199,*) k_det(iic,i,j,k)
			end do ! x,i
		end do !y,j
	end do !z, k
close(199)
end if

return
end 
