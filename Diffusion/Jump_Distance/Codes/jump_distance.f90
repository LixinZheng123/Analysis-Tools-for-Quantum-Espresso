program jump_distance 
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: celldm = 23.5170,&   !dimension of the supercell 
                           convertBA = 0.529    !Transfer Bohr to Angestrom
character(len=40)       :: filename 
integer                 :: step,         & 
                           i,j,l1,l2,    &      !General index of atom,
                           k,k0,dk,&            !steps
                           n,           &       !For counting numbers and normalize
                           nsp(2),      &       !number of atoms in O and H
                           numION,      &       !Number of ion in each timestep
                           dstep,       &       !step difference from *.ion
                           stepstart
integer                 :: stepstop_trans
integer                 :: stepstop_trans_nr
integer,allocatable     :: step_trans(:), iO_trans(:)
integer,allocatable     :: step_trans_nr(:), iO_trans_nr(:)
integer                 :: dummyI
integer                 :: iiO,iiOp
integer                 :: ncount=0
integer                 :: ierror
real(DP)                :: jump, d     
real(DP)                :: time,rIO(3),rIOp(3)
!
namelist /input/ filename, stepstop_trans, stepstop_trans_nr
!
!initialization
call init
!
!Open files
call files1
!
!read data from *pos
call main
!
call files2
!
!
! 
contains
!******************************************************************************************
!
!******************************************************************************************
subroutine init
!
implicit none
!
!read Namelist input for stdin
read(*,input)
!
!print Intro
write(*,*) '------------------------------------------'
write(*,*) '|  Mean Square Displacememt Calculation  |'
write(*,*) '------------------------------------------'
write(*,*) 'Input file  : ', trim(filename)//'.ion', trim(filename)//'.trans'
write(*,*) 'Output file : ', trim(filename)//'.jd'
!
!
allocate(step_trans(stepstop_trans), iO_trans(stepstop_trans))
allocate(step_trans_nr(stepstop_trans_nr), iO_trans_nr(stepstop_trans_nr))
!************************************
  
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files1
        !
        implicit none
        !
        open(unit=1, file=(trim(filename)//'.ion'), status='old')
        open(unit=2, file=(trim(filename)//'.trans'), status='old')
        open(unit=3, file=(trim(filename)//'.trans_no_r'), status='old')
        open(unit=4, file=(trim(filename)//'.jd'), status='unknown')
        !
end subroutine files1

subroutine files2
        implicit none
        close(1)
        close(2)
        close(3)
        close(4)
end subroutine files2 
!
!******************************************************************************************
subroutine main
!
implicit none
!
!************************************
!
!
do k=1,stepstop_trans 
  read(2,*,iostat=ierror) step_trans(k), iO_trans(k)
enddo
do k=1,stepstop_trans_nr
  read(3,*,iostat=ierror) step_trans_nr(k), iO_trans_nr(k)
enddo
!
!***********************
!
!Read .ion file
!********
!Jump out the first few lines
!if (stepstart .gt. 1) then
!  do k=1, stepstart-1
!    do i=1,3
!      read(1,*)
!    enddo
!  enddo
!endif
!********
!
!Initialization
ncount=0
k0=1
l1=2
l2=2
jump=0
iiO=0
rIO=0
ncount=0
!
!
do 
  if (l1 .gt. stepstop_trans) exit
  !
  if (ncount .eq. 0) then
    read(1,*) step, time
    read(1,*) iiO
    read(1,*) rIO(1:3)
    ncount=ncount+1
    cycle
  endif
  !
  do j=1,3
    rIOp(j)=rIO(j)
  enddo
  iiOp=iiO
  !
  read(1,*) step, time
  read(1,*) iiO
  read(1,*) rIO(1:3)
  !
  if (step_trans(l1) .lt. step) l1=l1+1
  if (step_trans_nr(l2) .lt. step) l2=l2+1
  !
  if (step .eq. step_trans(l1)) then
    call get_distance(rIO(1:3), rIOp(1:3), d)
    if (iiO .eq. iO_trans_nr(l2-1) .or. iiO .eq. iO_trans_nr(l2)) then
      jump=jump+d
    else
      jump=jump-d
    endif
    if (jump .lt. 0) write(*,*) "Error! Jump_distance smaller than 0!"
  endif
  ncount=ncount+1
  !
  write(4,*) step, time, jump*convertBA
enddo
   
end subroutine main
!
subroutine get_distance(pos1,pos2,dist)

implicit none
!
real(DP), intent(in)    :: pos1(3),pos2(3)
real(DP), intent(inout) :: dist
real(DP)                :: delta(3)
real(DP), parameter     :: pcell=23.5170
integer                 :: i
!
do i=1,3
  delta(i)=pos1(i)-pos2(i)
  delta(i)=delta(i)-nint(delta(i)/pcell)*pcell
enddo
dist=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
!
return
end subroutine get_distance

!******************************************************
end program jump_distance
