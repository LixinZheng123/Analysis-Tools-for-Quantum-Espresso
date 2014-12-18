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
integer                 :: i,l,l0,k,dk,last_k,last_dk,p,p1
integer                 :: dstep,stepstart,stepstop
integer                 :: ierror,ncount=0
integer                 :: nsp(2)
integer                 :: iiO0
integer                 :: last_step
integer                 :: dummy
integer                 :: hb0                  !First hydrogen bond index
integer,allocatable     :: readstep(:), iiO(:)
integer,allocatable     :: lh(:)                !Stands for h(t) (lower case of h) in equation
integer,allocatable     :: hb(:,:,:)    
integer,allocatable     :: tot_hb
real(DP)                :: tot_lh
real(DP)                :: avr_lh
real(DP)                :: Ci, Cc
real(DP)                :: tau_i,tau_c
real(DP)                :: Add
real(DP)                :: tau_exch=0
real(DP),allocatable    :: time(:)
real(DP),allocatable    :: uh(:)                !Stands for H(t) (upper case of h) in equation
!
namelist /input/ filename, timestep, dstep, last_dk, stepstart, stepstop, nsp
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
write(*,*) '|      Hydrogen Bond Time Correlation Function       |'
write(*,*) '------------------------------------------------------'
write(*,*) 'Input file   : ', trim(filename)//'.bond'
write(*,*) 'Output file  : ', trim(filename)//'.corr'
!
last_k=stepstop-stepstart+1
allocate(readstep(last_k),time(last_k),hb(4,nsp(1),last_k))
!
allocate(lh(last_dk),uh(last_dk))
!
end subroutine init
!
!=========================================================================================
!=========================================================================================
!
subroutine files1
implicit none
open(unit=4, file=(trim(filename)//'.bond'), status='unknown')
open(unit=5, file=(trim(filename)//'.corr'), status='unknown')
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
if (stepstart .ne. 1) then
  do k=1,(stepstart-1)*(nsp(1)+1)
    read(4,*)
  enddo
endif
! 
do k=1,last_k
  read(4,*,iostat=ierror) readstep(k), time(k)
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  ncount=ncount+1
  do i=1,nsp(1)
    read(4,*) dummy, dummy, dummy, hb(1:4,i,k)
  enddo
enddo
!
!
!
!lh=0
uh=0
tot_hb=0
!tot_lh=0
!
!last_k=40000
!do k=1, last_k-last_dk
do k=1,1
  tot_hb=0
  do i=1, nsp(1)
    do p=1,4
      if (hb(p,i,k) .eq. 0) exit
      tot_hb=tot_hb+1
    enddo
  enddo
  !
  Add=real(1)/tot_hb
  !write(*,*) "tot_hb=",tot_hb,"Add=",Add
  !
  !do i=1,nsp(1)
  !  do p=1,4
  !    if (hb(p,i,k) .eq. 0) exit
  !    hb0=hb(p,i,k)
  !    do dk=0,last_dk-1
  !      if (hb(1,i,k+dk) .eq. hb0 .or. hb(2,i,k+dk) .eq. hb0 .or. &
  !          hb(3,i,k+dk) .eq. hb0 .or. hb(4,i,k+dk) .eq. hb0) then
  !        uh(dk+1)=uh(dk+1)+Add
  !      else
  !        exit
  !      endif
  !    enddo
  !  enddo
  !enddo

  do i=1,nsp(1)
    do p=1,4
      if (hb(p,i,k) .eq. 0) exit
      hb0=hb(p,i,k)
      do dk=0,last_dk-1
        if (hb(1,i,k+dk) .eq. hb0 .or. hb(2,i,k+dk) .eq. hb0 .or. &
            hb(3,i,k+dk) .eq. hb0 .or. hb(4,i,k+dk) .eq. hb0) then
          uh(dk+1)=uh(dk+1)+Add
        else
          exit
        endif
      enddo
    enddo
  enddo

  do dk=0, last_dk-1
    write(5,*) readstep(dk+1)-readstep(1),time(dk+1)-time(1),uh(dk+1)
    tau_exch=tau_exch+uh(dk+1)
  enddo
  !
enddo
!
tau_exch=(tau_exch/(last_dk))*(time(last_dk)-time(1))
write(*,*) "******"
write(*,'("  tau_exch equals to ",F7.4)') tau_exch
write(*,*) "******"
write(*,*)

!
tau_c=0
!
!do dk=0, last_dk-1
!  Cc=uh(dk+1)/(last_k-last_dk)
!  write(*,*) uh(dk+1)
!  write(5,*) readstep(dk+1)-readstep(1),time(dk)-time(1),Cc
!  tau_c=tau_c+Cc
!  !write(5,*) (dk*dstep), timestep*(dk*dstep), Ci, Cc
!enddo
!tau_c=tau_c*timestep*dstep
!write(*,*) "******"
!write(*,'("  tau_exch equals to ",F6.3)') tau_c
!write(*,*) "******"
!write(*,*) 
!
end subroutine main
!
!=========================================================================================
!=========================================================================================
end program corr 
