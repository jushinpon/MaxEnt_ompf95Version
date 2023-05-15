! define all global parameters
!1.nx,ny,nz: replicated unit cell numbers in x, y, and z dimensions
!2. natom_unit: atom number in the unit cell (read input file)
!3.

MODULE INFORMATION
implicit real*8(a-h,o-z)
integer nx,ny,nz,natom_unit,nkMCshift ! replicate unit cell times in x, y, z dimensions,atom number in unit cell 
integer elemtype!total element type number, atom type
logical pbcx, pbcy, pbcz
real*8,allocatable::x(:),y(:),z(:),sx(:),sy(:),sz(:),xkeep(:),ykeep(:),zkeep(:)
integer,allocatable::atkeep(:),minatkeep(:)
real*8,allocatable::xmin(:),ymin(:),zmin(:),atomentropy(:),katomentropy(:),pairweight(:,:)
integer,allocatable::CN_No(:,:),CN_ID(:,:,:),atype(:),aid(:),tempatype(:)
integer natom,iterMC,inirand,snatom
character*30 name
real*8 rdfpeak(5),weight(5),ux(10),uy(10),uz(10),xdmin,ydmin,zdmin ! the first five rdf peaks from MS and entropy weight
real*8 xl,yl,zl,half_xl,half_yl,half_zl,rlist,confentropy,real_lattice
real*8 xlo,xhi,ylo,yhi,zlo,zhi
real*8 Entmin,kT,Entkeep,acceptratio,filterValue ! for MC
integer nadjust,naccept,keepNo ! for MC and counter for keeping ID with the large scoring values 
integer,allocatable:: keepeID(:) ! keep ID with the large scoring values

contains
!*********************************************
subroutine general_para
implicit real*8(a-h,o-z)
!real*8 shift
logical assigntype
real*8,allocatable::frac(:),MCtime(:) ! atom fraction coordinate shift
character*128 buildWay,filename,tempChar 

!call system('rm -rf output')
!call system('mkdir output')
     
open(112,file="00input.dat",status='old')
read(112,*) !skip comment (# The input file for...)
read(112,*) !skip comment (####### begin below)
read(112,*) ! the total MC iteration times,the total initial random search times
read(112,*) iterMC,inirand
allocate (MCtime(iterMC))
MCtime = 0.0
read(112,*) ! #! PBC in x, y, and z
read(112,*) pbcx,pbcy,pbcz
read(112,*)! rlist for cell list 
read(112,*)rlist
write(*,*)'rlist: ', rlist

read(112,*)! rdf peak(5)
do 9 i=1,5
	read(112,*)rdfpeak(i)
	write(*,*)"rdfpeak ",i,rdfpeak(i)
9 continue

read(112,*)!element type for HEA
read(112,*)natom,elemtype
write(*,*)'total atom number: ',natom
write(*,*)'element type number: ',elemtype

allocate (pairweight(elemtype,elemtype)) !! weighting for different pairs
write(*,*)"reading box information"
read(112,*) !ITEM: BOX BOUNDS pp pp pp
read(112,*)xlo, xhi
read(112,*)ylo, yhi
read(112,*)zlo, zhi  

write(*,*)xlo, xhi
write(*,*)ylo, yhi
write(*,*)zlo, zhi  
xl = xhi-xlo
yl = yhi-ylo
zl = zhi-zlo
half_xl = xl/2.d0    
half_yl = yl/2.d0    
half_zl = zl/2.d0

!read atom information

allocate (x(natom))
allocate (y(natom))
allocate (z(natom))

allocate (xkeep(natom))! configuration kept by MC
allocate (ykeep(natom))
allocate (zkeep(natom))

allocate (xmin(natom))
allocate (ymin(natom))
allocate (zmin(natom))

allocate (atype(natom))
allocate (atkeep(natom))
allocate (minatkeep(natom))

allocate (atomentropy(natom))
allocate (katomentropy(natom)) ! keep atom potential for the current best
allocate (keepeID(natom))

snatom = natom
allocate (aid(snatom))
allocate (tempatype(snatom))
allocate (sx(snatom))
allocate (sy(snatom))
allocate (sz(snatom))

!!! get atom information
read(112,*) !read coordinates
do 100 i=1,natom
	read(112,*)id_atom,itempatype,tempx,tempy,tempz !If id not sorted by lammps       
	x(id_atom)=tempx
	y(id_atom)=tempy
	z(id_atom)=tempz
	atype(id_atom)=itempatype
	if(i == natom)then
		write(*,*)"show the coordinates of last atom for check"
		write(*,*)id_atom,atype(id_atom),x(id_atom),y(id_atom),z(id_atom) !If id not sorted by lammps       
    endif 
100 continue
!move all atoms within the box
xdmin= minval(x)
ydmin= minval(y)
zdmin= minval(z)
x = x - xdmin
y = y - ydmin
z = z - zdmin
xl = xhi-xlo
yl = yhi-ylo
zl = zhi-zlo
half_xl = xl/2.d0    
half_yl = yl/2.d0    
half_zl = zl/2.d0

! assign the initial values for data output by lmpdata.f95
xkeep = x
ykeep = y
zkeep = z

xmin = x
ymin = y
zmin = z

atkeep =atype ! keep the atom types
	  
allocate(CN_No(natom,5)) ! consider five neighbor types
allocate(CN_ID(natom,5,50)) ! consider five neighbor types and corresponding IDs.
! assume the maximal number is 50.

!write(*,*)'1'
! B2 rdf peaks by MS for a B2 unit cell (lattice constant = 1)
! 0.87 1.0 1.41 1.65 1.73 
!!!!!!!!!!!!!!

!rdfpeak(1)=3.1
!rdfpeak(2)=4.4
!rdfpeak(3)=5.3
!rdfpeak(4)=10
!rdfpeak(5)=10
!rdfpeak(1)=0.87
!rdfpeak(2)=1.0
!rdfpeak(3)=1.4
!rdfpeak(4)=1.65
!rdfpeak(5)=1.73

!write(*,*)'3'
!!!!!!!!!!!!!!

!write(*,*)'4' 	
call lmpdata("INI",0)	
endsubroutine general_para 
	 
ENDMODULE
