program plot
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: celldm = 23.5170,&   !dimension of the supercell 
                           convertBA = 0.529    !Transfer Bohr to Angestrom
character(len=40)       :: filename 
integer                 :: step,stepp,         & 
                           i,j,l1,l2,    &      !General index of atom,
                           k,k0,dk,&            !steps
                           n,           &       !For counting numbers and normalize
                           nsp(2),      &       !number of atoms in O and H
                           numION,      &       !Number of ion in each timestep
                           dstep,       &       !step difference from *.ion
                           stepstart
integer                 :: stepstop_trans
integer                 :: stepstop_trans_nr
integer                 :: dummyI
integer                 :: iiO,iiOp
integer                 :: ncount=0
integer                 :: ierror
real(DP)                :: jump,jumpp
real(DP)                :: time,timep,rIO(3),rIOp(3)
!
namelist /input/ filename,stepstart
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
!************************************
  
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files1
        !
        implicit none
        !
        open(unit=1, file=(trim(filename)//'.jd'), status='unknown')
        open(unit=2, file=(trim(filename)//'.jd_p'), status='unknown')
        !
end subroutine files1

subroutine files2
        implicit none
        close(1)
        close(2)
end subroutine files2 
!
!******************************************************************************************
subroutine main
!
implicit none
!
!************************************
!
do k=1,stepstart-1
  read(1,*) stepp,timep,jumpp
enddo
!
do 
  read(1,*,iostat=ierror) step,time,jump
  if (ierror .lt. 0) exit
  write(2,*) step-stepp, time-timep, jump-jumpp
enddo
   
end subroutine main

!******************************************************
end program plot
