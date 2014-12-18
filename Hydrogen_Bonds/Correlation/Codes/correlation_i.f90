program corr
!
!This fortran file is to find first solvation shell (FSS) of a defect from a MD simulation trajectory.
!It is used for visualization and selection using VMD.
!Author: Lixin Zheng
!May 2014
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: celldm = 23.5170     !dimension of the supercell 
real(DP), parameter     :: convertBA= 0.529     !Transfer Bohr to Angestrom9
real                    :: timestep
character(len=40)       :: filename
integer                 :: i,l,l0,k,dk,last_k,last_dk
integer                 :: dstep
integer                 :: ierror,ncount=0
integer                 :: iiO0
integer                 :: last_step
integer,allocatable     :: readstep(:), iiO(:)
integer,allocatable     :: lh(:)                !Stands for h(t) (lower case of h) in equation
integer,allocatable     :: uh(:)                !Stands for H(t) (upper case of h) in equation
real(DP)                :: tot_lh
real(DP)                :: avr_lh
real(DP)                :: Ci, Cc
real(DP)                :: tau_i,tau_c
!
namelist /input/ filename, timestep, dstep, last_dk
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
write(*,*) '------------------------------------------------------'
write(*,*) '|     Proton Transfer Time Correlation Function      |'
write(*,*) '------------------------------------------------------'
write(*,*) 'Input file   : ', trim(filename)//'.trans'
write(*,*) 'Output file  : ', trim(filename)//'.corr_i'
!
!
end subroutine init
!
!=========================================================================================
!=========================================================================================
!
subroutine files1
implicit none
open(unit=4, file=(trim(filename)//'.trans'), status='unknown')
open(unit=5, file=(trim(filename)//'.corr_i'), status='unknown')
end subroutine files1
!
subroutine files2
implicit none
close(4)
close(5)
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
do 
  read(4,*,iostat=ierror) last_step
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  ncount=ncount+1
enddo
close(4)
open(unit=4, file=(trim(filename)//'.trans'), status='unknown')
!
!
last_k=last_step/dstep
write(*,'("   Total steps and time to be calculated: ",I8, F7.2)') last_k, timestep*last_step
allocate(lh(last_dk),uh(last_dk))
!allocate(lh(last_k),uh(last_k))
allocate(readstep(ncount),iiO(ncount))
!
lh=0
uh=0
tot_lh=0
!
!
do l=1,ncount
  read(4,*) readstep(l), iiO(l)
enddo
!
!
!last_k=40000
tot_lh=0
do k=1, last_k-last_dk
!do k=1, last_k
  !
  !Find iiO0
  do l=1,ncount
    if (readstep(l)/dstep .le. k) then
      cycle
    else
      exit
    endif
  enddo
  !
  l0=l-1
  iiO0=iiO(l0)
  !write(*,*) "iiO0=", iiO0
  !
  !Initialization
  l=l0
  !
  !
  do dk=0,last_dk-1
    !
    if (readstep(l0+1)/dstep .gt. (k+dk)) uh(dk+1)=uh(dk+1)+1
    !
    if (l .ge. ncount) exit
    if (iiO(l) .eq. iiO0) then
      lh(dk+1)=lh(dk+1)+1
      tot_lh=tot_lh+1
    endif
    !
    do
      if (readstep(l+1)/dstep .le. (k+dk)) then
        l=l+1
      else
        exit
      endif
    enddo
    !
  enddo
enddo
!
avr_lh=tot_lh/((last_k-last_dk)*last_dk)
!
write(*,'("  <h>=",F8.2)') avr_lh
write(*,*)
!do dk=1,50,5
!  write(*,'("  <h(0)h(",I2,")>=",F5.2)') dk, real(lh(dk))/last_k
!enddo
!
Ci=0
Cc=0
tau_c=0
do dk=0, last_dk-1
  !Ci=(real(lh(dk))/last_k)/avr_lh
  !Cc=(real(uh(dk))/last_k)/avr_lh
  Ci=(real(lh(dk))/(last_k-last_dk))
  !Cc=(real(uh(dk))/(last_k-last_dk))
  !tau_c=tau_c+Cc
  !write(5,*) (dk*dstep), timestep*(dk*dstep), Ci, Cc
  write(5,*) (dk*dstep), timestep*(dk*dstep), Ci
enddo
!tau_c=tau_c*timestep*dstep
!write(*,*) "******"
!write(*,'("  tau_exch equals to ",F6.3)') tau_c
!write(*,*) "******"
write(*,*) 
!
end subroutine main
!
!=========================================================================================
!=========================================================================================
end program corr 
