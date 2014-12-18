!This fortran code is to find first solvation shell (FSS) of a defect from a MD simulation trajectory.
!It is used for visualization and selection using VMD.
!Lixin Zheng, May 2014
!Modified 20140901, add error message.
!Modified 20141014, refactor with rdf codes.
!==================================================================================================
program find_FSS
   !
   implicit none
   !
   !======Global variables======
   character(len=40)    :: filename
   integer, parameter   :: DP = selected_real_kind(14,200) !double-precision kind
   real(DP), parameter  :: pi=4.d0*atan(1.d0)
   real, parameter      :: celldm = 23.5170     !dimension of the supercell 
   real, parameter      :: convertBA= 0.529     !Transfer Bohr to Angestrom9
   real, parameter      :: r_FSS=3.5            !distance (Angestrom) between O and O for first solvation
   real, parameter      :: r_OH=1.24            !bond length of OH in a regular H2O 
   integer              :: nsp(2)               !number of each atomic species
   integer              :: cs
   integer              :: stepstart
   integer              :: stepstop
   integer              :: ncount=0
   integer              :: step 
   integer              :: ierror
   real                 :: time
   !======Variables concerning .pos======
   real,allocatable     :: rO(:,:),rH(:,:)
   !======Variables concerning .fss======
   real                 :: rO_FSS(3,15)         !atomic positions (dimension,atomic-species)
   real                 :: rH_FSS(3,30)         !atomic positions (dimension,atomic-species)
   real                 :: rIO_record(3,5)      !Ion position (dimension,numION)
   real                 :: rIH_record(3,3)      !ion positions (dimension,numH)
   real                 :: rIO(3)               !Real ion position (dimension)
   real                 :: rIH(3,3)             !Real ion position (dimension,numH)
   integer              :: new_nsp(2)           !number of each atomic species
   integer              :: iO_FSS(15)
   integer              :: iH_FSS(30)
   integer              :: iIH_record(5,4)      !index of Ion hydrogen (numH,numION)
   integer              :: iIO_record(5)        !index of O* read from Ion
   integer              :: iIO                  !index of Ion oxygen
   integer              :: iIH(3)               !index of Ion oxygen
   integer              :: numION               !Number of Ions in each timestep. 
   integer              :: numH                 !Number of Hydrogen in each ion. cs=2, numH=1; cs=3, numH=3
   integer              :: iIOp(2),iiHp(4)
   integer              :: iiH_trans            !Find the hydrogen, through which PT occurs.
   integer              :: cov_bond(4,15)       !The covalent bonded hydrogen to each oxygen
   !======Local variables======
   integer              :: e_step
   integer              :: j,k,p,q,q1
   integer              :: start=0
   integer              :: numW                 !number of water
   integer              :: error
   real                 :: d, dr(5)
!
namelist /input/ cs,filename, nsp, stepstart, stepstop
!
!initialization
call init
!
!Open files
call files1
!
!main loop
call main
!
call files2
!
!
! 
contains
!
!=========================================================================================
!=========================================================================================
!
subroutine init
!
implicit none
!read Namelist input for stdin
read(*,input)
!
!print Intro
write(*,*) '--------------------------------------------------'
write(*,*) '|     Find the Ion and First Solvation Shell      |'
write(*,*) '--------------------------------------------------'
write(*,*) 'Position file      : ', trim(filename)//'.pos'
write(*,*) 'Defect file        : ', trim(filename)//'.ion'
write(*,*) '1st solvation file : ', trim(filename)//'.fss'
write(*,*) 'Ion transfer file  : ', trim(filename)//'.trans ',trim(filename)//'.trans_no_r'
write(*,*) 'Ion error file     : ', trim(filename)//'.error'
!
!allocation
allocate(rO(3,nsp(1)))
allocate(rH(3,nsp(2)))
!Initiation
per=0
iIO=0
iIO_record=0
iIOp=0
iiHp=0
if (cs .eq. 2) numH=1
if (cs .eq. 3) numH=3
!
end subroutine init
!
!=========================================================================================
!=========================================================================================
!
subroutine files1
implicit none
open(unit=1, file=(trim(filename)//'.pos'), status='old')
open(unit=2, file=(trim(filename)//'.ion'), status='unknown')
open(unit=3, file=(trim(filename)//'.fss'), status='unknown')
open(unit=4, file=(trim(filename)//'.trans'), status='unknown')
open(unit=5, file=(trim(filename)//'.error'), status='unknown')
open(unit=10, file=(trim(filename)//'.trans_no_r'), status='unknown')
!write(5,*) "This file is to record these occasions when ion has gone out of the first solvation shell of previous step."
end subroutine files1
!
subroutine files2
implicit none
close(1)
close(2)
close(3)
close(4)
close(5)
close(10)
end subroutine files2 
!
!=========================================================================================
!=========================================================================================
!
subroutine main
!
implicit none
!
!************************************
!main loop
do 
  !**************************************************************
  !end of file
  !**************************************************************
  read(1,*,iostat=ierror) step, time 
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  ncount=ncount+1
  !**************************************************************
  !**************************************************************
  !
  !MODULE 0: READ FILE .POS
  !
  !**************************************************************
  !**************************************************************
  do i = 1, nsp(1)
    !Read position of Oxygen
    read(1,*) rO(1:3,i)
  enddo
  do i = 1,nsp(2)
    !read position of Hydrogen
    read(1,*) rH(1:3,i)
  enddo
  !**************************************************************
  !**************************************************************
  !
  !MODULE 1: FIND THE ION
  !ONLY numION, INDEX iIO(1), iiH(1:numION,1), AND COORDINATES rIO(1:3,1) NEED TO BE PASSED TO NEXT MODULE.
  !
  !**************************************************************
  !**************************************************************
  numION = 0
  numW=0
  p=0
  e_step=0
  do iO=1,nsp(1)
    q=0
    do iH=1, nsp(2)
      call get_distance(rO(1:3,iO),rH(1:3,iH),d)
      if (d*convertBA .lt. r_OH) then
        q=q+1
        iiH(q,p+1)=iH
      endif
    enddo
    !*************************
    if (q .eq. 2) then
      numW=numW+1
    else if (q .eq. numH) then
      p=p+1
      iIO(p)=iO
    end if
    !**************************
  enddo
  numION=p
  do i=1, numION
    iIO_record(i)=iIO(i)
  enddo
  !**************************
  !1. The most straight forward situation
  !**************************
  if ((numION+numW) .eq. nsp(1) .and. numION .eq. 1) then
    do j=1,3
      rIO(j,1) = rO(j,iIO(1))
    enddo
  endif
  !**************************
  !2. The different situations that would probably go wrong.
  !Just skip the judgement part, and adopt the previous ion index.
  !**************************
  if ((numION+numW) .ne. nsp(1) .or. numION .gt. 2 .or. numION .lt. 1) then 
    !if ((numION+numW) .ne. nsp(1) .or. numION .gt. 2) then
    !  write(*,*) time, numION, numION+numW
    !endif
    iIO(1)=iIOp(1)
    do i=1,numH
      iiH(i,1)=iiHp(i)
    enddo
    do j=1,3
      rIO(j,1)=rO(j,iIO(1))
    enddo
  endif
  !**************************
  !3. When numION=2, try to find the real ion
  !**************************
  if ((numION+numW) .eq. nsp(1) .and. numION .eq. 2) then
    do p=1,2
      do j=1,3
        rIO(j,p) = rO(j,iIO(p))
      enddo
    enddo
    !**************************
    !compare and change order
    !so that the total distances from O to H of iIO(1) is the shortest
    !If cs=2, OH-
    if (cs .eq. 2) then
      do p=1,numION
        call get_distance(rIO(1:3,p),rH(1:3,iiH(1,p)),dr(p))
      enddo
      do p=2,numION
        if(dr(p) .lt. dr(1)) then
         iIO(1)=iIO(p)
          dr(1)=dr(p)
          do j=1,3
            rIO(j,1)=rIO(j,p)
            iiH(1,1)=iiH(1,p)
          enddo
        endif
      enddo
    !IF cs=3, Find the mutual iiH that ion1 and ion2 has
    else if (cs .eq. 3) then
      do q=1,3
        do q1=1,3
          if (iiH(q,1) .eq. iiH(q1,2)) then
            iiH_real=iiH(q,1)
            exit
          endif
        enddo
      enddo
      do p=1,2
        call get_distance(rIO(1:3,p),rH(1:3,iiH_real),dr(p))
      enddo
      do p=2,numION
        if(dr(p) .lt. dr(1)) then
         iIO(1)=iIO(p)
          dr(1)=dr(p)
          do j=1,3
            rIO(j,1)=rIO(j,p)
            do q=1,numH
              iiH(q,1)=iiH(q,p)
            enddo
          enddo
        endif
      enddo
    endif
  endif
  !
  !
  !***
  !Modified 20140901
  !To determine if the ion has transfered too fast
  error=1
  !
  if (ncount .gt. 1) then
    do i=1,new_nsp(1)
      if (iIO(1) .eq. iO_FSS(i) ) then
        error=0
        exit
      endif
    enddo
  endif
  if (error .eq. 1 .and. ncount .gt. 1) then
    write(5,*) "Error1: The new ion has gone out of the FSS of previous ion."
    e_step=step
    write(5,*) step,time,iIO(1) 
  endif
  !***
  !
  !
  !**************************
  !Write in .ion
  !**************************
  write(2,*) step,time,numION
  write(2,*) iIO(1),iiH(1:numH,1)
  write(2,*) rIO(1:3,1)
  !
  !**************************
  !**************************
  if (iIO(1) .ne. iIOp(1)) then
    !Write in .trans
    write(4,*) step, iIO(1), iiH(1:numH,1)
    !Modified 20141014
    !We also want to write out the non-rattling-transfer index in the same time.
    if (iIO(1) .ne. iIOp(2)) then
      write(10,*) step, iIO(1)
    endif
    iIOp(2)=iIOp(1)
    iIOp(1)=iIO(1)
    !End of modification
  endif
  do i=1,numH
    iiHp(i)=iiH(i,1)
  enddo
  !
  !
  !
  !**************************************************************
  !**************************************************************
  !
  !MODULE 2: FIND THE 1ST SOLVATION SHELL 
  !ONLY numION, INDEX iIO(1), iiH(1:numION,1), AND COORDINATES rIO(1:3,1) NEED TO BE PASSED FROM LAST MODULE.
  !
  !**************************************************************
  !**************************************************************
  !
  !**************************
  !This if statement is to make sure that the 1st ion has been found;
  !If not, we will write nothing to the .fss file.
  !**************************
  if (start .eq. 0 .and. numION .eq. 0) then
    cycle
  else if (start .eq. 0 .and. numION .ne. 0) then
    start=1
  endif
  !**************************
  !
  !**************************
  !Find FSS oxygens
  !**************************
  p=0
  q=0
  cov_bond=0
  do i=1,nsp(1)
    call get_distance(rO(1:3,i),rIO(1:3,1),d)
    if (d*convertBA .lt. r_FSS) then
      p=p+1  
      iO_FSS(p)=i
      do j=1,3
        rO_FSS(j,p)=rO(j,i)
      enddo
    endif
  enddo
  new_nsp(1)=p
  !**************************
  !Find the hydrogens relates to FSS oxygens
  !**************************
  do i=1,nsp(2)
    call get_distance(rH(1:3,i),rIO(1:3,1),dr(1))
    do p=1,new_nsp(1)
      call get_distance(rH(1:3,i),rO_FSS(1:3,p),dr(2))
      if (dr(2)*convertBA .lt. r_OH .or. dr(1)*convertBA .lt. (r_FSS-r_OH)) then
        q=q+1
        iH_FSS(q)=i
        do j=1,3
          rH_FSS(j,q)=rH(j,i)
        enddo
        exit
      endif     
    enddo
  enddo
  new_nsp(2)=q
  !**************************
  !Record the link between O and H
  !**************************
  do p=1,new_nsp(1)
    k=0
    do q=1, new_nsp(2) 
      call get_distance(rH_FSS(1:3,q),rO_FSS(1:3,p),d)
      if (d*convertBA .lt. r_OH) then
        k=k+1
        cov_bond(k,p)=iH_FSS(q)
      endif
    enddo
  enddo
  !**************************
  !Write index into file.fss 
  !**************************
  if (numION .eq. 1) write(3,*) step, time, numION, new_nsp(1),new_nsp(2)
  if (numION .ne. 1) write(3,*) step, time, numION, new_nsp(1),new_nsp(2), iIO_record(1:numION)
  write(3,*) iIO(1),iiH(1:numH,1)
  do p=1,new_nsp(1)
    write(3,*) iO_FSS(p),rO_FSS(1:3,p),cov_bond(1:4,p)
  enddo
  do q=1,new_nsp(2)
    write(3,*) iH_FSS(q),rH_FSS(1:3,q)
  enddo
  !
  !
  !Modified 20140901
  !if (e_step .eq. step) then
    if (cs .eq. 2) then
      if ((new_nsp(1)*2-1) .ne. new_nsp(2) .and. numION .ne. 0) then
        write(5,*) "Error2: (nsp(1)*2-1) .ne. nsp(2)!"
        write(5,*) step, time, "numION:", numION, new_nsp(1),new_nsp(2)
        write(5,*) "iIO:",iIO_record(1:numION)
      endif
    endif
    if (cs .eq. 3) then
      if ((new_nsp(1)*2+1) .ne. new_nsp(2) .and. numION .ne. 0) then
        write(5,*) "Error2: (nsp(1)*2+1) .ne. nsp(2)!"
        write(5,*) step, time, "numION:", numION, new_nsp(1),new_nsp(2) 
        write(5,*) "iIO:",iIO_record(1:numION)
      endif
    endif
  !endif
enddo

!
end subroutine main
!
!=========================================================================================
!=========================================================================================
!
subroutine get_distance(pos1,pos2,dist)

implicit none
!
real, intent(in)    :: pos1(3),pos2(3)
real, intent(inout) :: dist
real                :: delta(3)
real, parameter     :: pcell=23.5170
integer             :: i
!
do i=1,3
  delta(i)=pos1(i)-pos2(i)
  delta(i)=delta(i)-nint(delta(i)/pcell)*pcell
enddo
dist=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
!
return
end subroutine get_distance
!
!=========================================================================================
!=========================================================================================
!
end program find_FSS
