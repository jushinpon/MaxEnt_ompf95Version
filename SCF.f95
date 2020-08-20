!subroutine SCF (tid,natom,atype,confentropy,atomentropy,filterValue,keepNo,weight,CN_No,CN_ID)
!subroutine SCF (natom,atype,keepNo,CN_No,CN_ID)
subroutine SCF
use information
use omp_lib
implicit real*8(a-h,o-z)
integer iscf,ishift,all_keepNo,tid
integer startid,endid,Quotient,Remainder,nthreads ! for mpi-like loop id and remainder 
integer,allocatable:: omp_counter(:),omp_keepNo(:,:)
!real*8  oldconfentropy,SCFEntmin
!real*8  newent,oldent

!integer natom,keepeID(natom),all_keepeID(natom),CN_No(natom,5),CN_ID(natom,5,50),CN_type(natom,5,10)  
integer CN_type(natom,elemtype),iswap_time  

iscf = 0	! counter for scf steps	
iswap_time = 1 ! counter for successful swap atom time
nthreads = 1 ! for serial version
!$ nthreads = omp_get_max_threads()

Quotient =  natom/nthreads
Remainder = MOD(natom,nthreads)
allocate (omp_counter(0:nthreads-1))
allocate (omp_keepNo(0:nthreads-1,0:Quotient)) !Quotient+1 for array element number

do while (iswap_time .gt. 0 .and. iscf .lt. 5)
all_keepNo = 0
 !!$ start = omp_get_wtime()  
	
!$OMP PARALLEL PRIVATE(tid,i,startid,endid,neID,ntemp,JID,itempatype,all_keepNo) &
!$OMP shared(omp_counter,omp_keepNo,CN_type,atype,iscf)
!$OMP WORKSHARE
CN_type = 0 
!$OMP END WORKSHARE
tid = 0 ! for the serial version 
!$ tid = OMP_GET_THREAD_NUM()
omp_counter(tid) = -1 ! array starts from 0
!omp_keepNo
if(tid < remainder)then ! rang from 0
    startid = 1+tid*(Quotient + 1)
    endid = (tid+1)*(Quotient + 1)
else
    startid = 1+(tid - Remainder )*Quotient + Remainder*(Quotient + 1)
    endid = (tid - Remainder + 1 )*Quotient + Remainder*(Quotient + 1)
endif
!print *,"tid: ",tid, "startid",startid, "endid", endid
do 1112 i=startid,endid
 	ntemp = CN_No(i,1)
	keepindex = 0 
	do 311 neID = 1,ntemp ! atom ID of a neighbor type
 		JID = CN_ID(i,1,neID) 
		itempatype = atype(JID)
		CN_type(i,itempatype) = CN_type(i,itempatype) + 1

		if(atype(i) .eq. atype(JID) .and. keepindex .eq. 0)then
			keepindex = 1
			omp_counter(tid) = omp_counter(tid) + 1
			all_keepNo = omp_counter(tid)
        	omp_keepNo(tid,all_keepNo) = i
 		endif 
    311 continue         
1112 continue 
!!$OMP END PARALLEL
 
!$omp barrier

!$omp single
! combine all local arrays to one
keepNo = 0
do i = 0,nthreads-1
	istart = keepNo + 1
	keepNo = keepNo + omp_counter(i) + 1 !(omp_counter(i) + 1) is real No.
	keepeID(istart:keepNo) = omp_keepNo(i,0:omp_counter(i))
!	print *,"tid ",i, "omp_counter:",omp_counter(i) 
enddo

iscf = iscf + 1 
iswap_time = 0
!$omp end single

!!$OMP PARALLEL 
!$OMP DO PRIVATE(itri) 
  do itri = 0,keepNo*(keepNo-1)/2 - 1
  !    do  jj1=ii1+1,keepNo
	 ii1 = itri/keepNo
	 jj1 = mod(itri,keepNo)
	 if (jj1<=ii1)then
     ii1 = keepNo- ii1 -2
	 jj1 = keepNo- jj1 -1
	 endif
	 ii1 = ii1 + 1
	 jj1 = jj1 + 1

        ii = keepeID(ii1)
        jj = keepeID(jj1)		
	   	itype = atype(ii) !type
	   	jtype = atype(jj)
	   if(itype .ne. jtype) then
	   ! total number of the same types at site ii and jj	
	   		ibefore_swap = CN_type(ii,itype) + CN_type(jj,jtype) 
	   		iafter_swap = CN_type(ii,jtype) + CN_type(jj,itype)	   
	   	    if (iafter_swap .lt. ibefore_swap) then ! could get a lower conentropy
				iswap_time  = 1
			ishift = atype(ii)
		        atype(ii) = atype(jj)
		        atype(jj) = ishift

				itype_org = atype(jj)
				jtype_org = atype(ii)

				itype_swap = atype(ii)
				jtype_swap = atype(jj)

				ntemp = CN_No(ii,1)! the first neighbour site number of site ii
				do neID = 1,ntemp ! atom ID of a neighbor type
				  JID = CN_ID(ii,1,neID) ! the nth atom ID of atom i's neighbor atom

				  CN_type(JID,itype_org) = CN_type(JID,itype_org) - 1 ! the nth atom ID of atom i's neighbor atom
				  CN_type(JID,itype_swap) = CN_type(JID,itype_swap) + 1 ! the nth atom ID of atom i's neighbor atom
				enddo
				ntemp = CN_No(jj,1)! the first neighbour atom number of atom ii
				do neID = 1,ntemp ! atom ID of a neighbor type
				  JID = CN_ID(jj,1,neID) ! the nth atom ID of atom i's neighbor atom
				  CN_type(JID,jtype_org) = CN_type(JID,jtype_org) - 1 ! the nth atom ID of atom i's neighbor atom
				  CN_type(JID,jtype_swap) = CN_type(JID,jtype_swap) + 1 ! the nth atom ID of atom i's neighbor atom
				enddo
          endif
       endif
	  
	  end do                    !2 continue	 

!	 end do   	 !1 continue
!$OMP END DO
!$OMP END PARALLEL	
!!$ finish = omp_get_wtime()  
!  !$      print *,"iscf:",iscf,"****maxThreads:",omp_get_max_threads(),"time consumption:",finish-start," sec."

enddo ! while loop
return
end
