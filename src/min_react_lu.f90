!	Reactive mineral lookup
!	R. Maxwell
!	8/21/00
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
subroutine react_min_lu(nx, ny, nz, porosity, r, min_file,n_constituents,constit_num)

integer*4 nx,ny,nz, cat_table(100), perm_cat_num, min_cat_num, i, j, icat_temp, ii, finder,&
 constit_num, iic, ircheck, iphicheck,n_constituents, k, min_cat_temp

real*8 porosity(nx,ny,nz),r(n_constituents,nx,ny,nz),cat_value(100,30)

integer nnx,nny,nnz
integer*4,allocatable::perm_ind(:,:,:),min_ind(:,:,:)

character*100 min_file, min_cat_file, perm_cat_file
character*25 rcheck_file(10),phicheck_file(10)

rcheck_file(1) = 'rcheck.c1.asc.txt'
rcheck_file(2) = 'rcheck.c2.asc.txt'
rcheck_file(3) = 'rcheck.c3.asc.txt'
rcheck_file(4) = 'rcheck.c4.asc.txt'
rcheck_file(5) = 'rcheck.c5.asc.txt'

phicheck_file(1) = 'phicheck.c1.asc.txt'
phicheck_file(2) = 'phicheck.c2.asc.txt'
phicheck_file(3) = 'phicheck.c3.asc.txt'
phicheck_file(4) = 'phicheck.c4.asc.txt'
phicheck_file(5) = 'phicheck.c5.asc.txt'

allocate(perm_ind(nx,ny,nz),min_ind(nx,ny,nz))


print*, 'nx,ny,nz',nx,ny,nz
print*, ' constit ',constit_num
open (120, file=trim(min_file))

read(120,*)
read(120,*)
read(120,*)
read(120,*)
read(120,*) min_cat_file
read(120,*) perm_cat_file
read(120,*) perm_cat_num
read(120,*) min_cat_num
read(120,*)

! write out to log file 
write(666,*) ' Porosity and Chemical Retardation Lookup'
write(666,*)
write(666,*) '	reading permeability catagories from:'
write(666,*) trim(perm_cat_file)
write(666,*) '	reading mineral catagories from:'
write(666,*) trim(min_cat_file)
write(666,*) ' num perm catagories :' , perm_cat_num
write(666,*) ' num minerals :' , min_cat_num
write(666,*)
iic = constit_num
write(666,*) ' constituent number:',iic
write(666,*)
! read in lookup table
do i = 1 , perm_cat_num
	read(120,*) cat_table(i), cat_value(i,1:(2+2*min_cat_num))
	write(666,*) ' catagory:',cat_table(i)
	write(666,*) cat_value(i,1:(2+2*min_cat_num))
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

r(constit_num,1:nx,1:ny,1:nz) = 1.d0


print*, ' min ind'
! read in mineral indicators
open(140,file=trim(min_cat_file))
	read(140,*) nnx,nny,nnz
!	print*, nnx,nny,nnz
	do k=1,nz
		do j=1,ny
			do i =1,nx
				read(140,*) min_ind(i,j,k)
			end do ! x,i
		end do !y,j
	end do !z, k
close(140)

do i = 1, nx
	do j = 1, ny
		do k = 1, nz
			icat_temp = finder(perm_ind(i,j,k),cat_table,perm_cat_num) 
			if (icat_temp.le.0.) print*, ' icat <0! ', icat_temp
			if (icat_temp.gt.perm_cat_num) print*, ' icat > max cats ', icat_temp
			porosity(i,j,k) = cat_value(icat_temp,1)
!			print*,i,j,k,icat_temp,cat_value(icat_temp,1)
		end do
	end do
end do

do i = 1, nx
	do j = 1, ny
		do k = 1, nz
		icat_temp = finder(perm_ind(i,j,k),cat_table,perm_cat_num)
		r(iic,i,j,k) = 	r(iic,i,j,k) + 10.**cat_value(icat_temp,2)
		min_cat_temp = min_ind(i,j,k)
			do ii = 1,min_cat_num
!			print*,i,j,k,min_cat_temp
			 if (mod(min_cat_temp,2).eq.0.) then
				r(iic,i,j,k) = r(iic,i,j,k)  + 10**cat_value(icat_temp,2+2*ii)
				min_cat_temp = min_cat_temp/2
				else
				r(iic,i,j,k) = r(iic,i,j,k) + 10**cat_value(icat_temp,1+2*ii )
				min_cat_temp = (min_cat_temp-1)/2
				end if
			 end do !ii
		end do !k
	end do !j
end do !i

ircheck = 1
iphicheck = 1

if (ircheck.eq.1) then
open (199,file=trim(rcheck_file(constit_num)))
	write(199,*) nx,ny,nz
	do k=1,nz
		do j=1,ny
			do i =1,nx
				write(199,*) r(iic,i,j,k)
			end do ! x,i
		end do !y,j
	end do !z, k
close(199)
end if

if (iphicheck.eq.1) then
open (199,file=trim(phicheck_file(constit_num)))
	write(199,*) nx,ny,nz
	do k=1,nz
		do j=1,ny
			do i =1,nx
				write(199,*) porosity(i,j,k)
			end do ! x,i
		end do !y,j
	end do !z, k
close(199)
end if

return
end 

function finder(lookup_val,cat_table,max_cats)
integer*4 lookup_val, max_cats, i,finder
integer*4 cat_table(100)

finder = -99
!print*, max_cats, lookup_val
do i=1, max_cats
	if(cat_table(i).eq.lookup_val) finder = i
end do
return
end
		