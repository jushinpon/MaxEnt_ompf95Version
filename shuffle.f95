subroutine shuffle !(natom,atype)
use information
implicit real*8(a-h,o-z)   
integer natom1
integer atype1(natom)
integer :: i, randpos, temp
real :: r
natom1 = natom
atype1 = atype
!!$OMP PARALLEL PRIVATE(i) shared(atype)
!!$OMP DO
do i = natom1, 2, -1
  call random_number(r)
  randpos = int(r * i) + 1
  if (randpos .gt. i) randpos = i 
  temp = atype1(randpos)
  atype1(randpos) = atype1(i)
  atype1(i) = temp
end do
!!$OMP END DO
!!$OMP END PARALLEL  
atype = atype1

return 
end 