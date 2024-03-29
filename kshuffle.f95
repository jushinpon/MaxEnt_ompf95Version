!subroutine kshuffle(tid,imc,natom,atype,confentropy,atomentropy,weight,CN_No,CN_ID)
subroutine kshuffle
use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift
real*8 r
integer natom1
integer atype1(natom)
integer :: i, randpos, temp

natom1 = natom
atype1 = atype

ikNo = int(keepNo*0.05)
if(ikNo .lt. 10) ikNo = 10
do 1 ikMC=1, ikNo
	call random_number(r)
	ikid = int(dint(keepNo*r)+1)
    i1 = keepeID(ikid) !site ID
	irepeat  = 0 ! parameter to count the 1112 continuoue times
	if (atype(i1) .eq. itypeNotSwap ) cycle
! pick the second one 
12345 irepeat  = irepeat + 1
if (irepeat .gt. 20) cycle ! skip this loop 
	call random_number(r)
    i2 = int(dint(natom*r)+1)
	if(i1 .ne. i2 .and. atype(i1) .ne. atype(i2) .and. atype(i2) .ne. itypeNotSwap ) then
    	ishift = atype1(i1)
    	atype1(i1) = atype1(i2)
		atype1(i2) = ishift
	else
		goto 12345		  
	endif

1 continue ! atom type swap loop
atype = atype1
return
end
