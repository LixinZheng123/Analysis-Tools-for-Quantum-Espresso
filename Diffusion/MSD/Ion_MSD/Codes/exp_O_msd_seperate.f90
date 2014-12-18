program mean_square_displacement
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: celldm = 23.5170,&   !dimension of the supercell 
                           convertBA = 0.529    !Transfer Bohr to Angestrom
real(DP)                :: cutoff1=0.7*celldm   !
real(DP)                :: cutoff2=0.25*celldm  !
character(len=40)       :: filename 
integer                 :: cs,          &       !case = 1, mean square displacement of all Oxydgen
                                                !case = 2, msd of O* at OH-
                                                !case = 3, msd of O* at H3O+
                                                !case = 4, square displacement
                                                !          of one Oxygen 
                           i,ip,j,p,    &       !General index of atom,
                           k,k0,dk,&            !steps
                           n,           &       !For counting numbers and normalize
                           nsp(2),      &       !number of atoms in O and H
                           numION,      &       !Number of ion in each timestep
                           dstep,       &       !step difference from *.ion
                           stepstart
integer                 :: dummyI
integer                 :: iO,iiO,iiH
integer                 :: ncount=0
integer                 :: ierror
integer                 :: last_dk
integer                 :: average_box(3)
integer                 :: e_step(5), e_iiO(5)
real(DP)                :: sum,         &
                           msd,         &       !mean square displacement
                           tot_sd,      &       !square displacement before nomalization
                           d2,          &
                           dummyR               !a dummy variable
real(DP)                :: tmp_dr
real(DP)                :: dr(3)
real(DP)                :: amass(2), pmass(2)
real(DP),allocatable    :: time(:),     &      
                           rIO(:,:), drp(:,:)
integer,allocatable     :: box(:,:),step(:)
!
!Added Dec 2014
integer                 :: stepstop_trans
integer                 :: stepstop_trans_nr
integer,allocatable     :: step_trans(:), iO_trans(:)
integer,allocatable     :: step_trans_nr(:), iO_trans_nr(:)
REAL(DP)                :: msd_j, msd_v
!
namelist /input/ filename,  cs, stepstart, dstep, last_dk, &
                 stepstop_trans, stepstop_trans_nr         !Added Dec 2014
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
write(*,*) 'ION file : ', trim(filename)//'.ion'
write(*,*) 'MSD file : ', trim(filename)//'.msd'
!
!
!Added Dec 2014
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
        open(unit=2, file=(trim(filename)//'.trans'), status='old')      !Added Dec 2014
        open(unit=5, file=(trim(filename)//'.trans_no_r'), status='old') !Added Dec 2014
        open(unit=3, file=(trim(filename)//'.msd'), status='unknown')
        open(unit=4, file=(trim(filename)//'.msd_error'), status='unknown')
        !
end subroutine files1

subroutine files2
        implicit none
        close(1)
        close(2) !Added Dec 2014
        close(3)
        close(4)
        close(5) !Added Dec 2014
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
!Added Dec 2014
do k=1,stepstop_trans
  read(2,*,iostat=ierror) step_trans(k), iO_trans(k)
enddo
do k=1,stepstop_trans_nr
  read(5,*,iostat=ierror) step_trans_nr(k), iO_trans_nr(k)
enddo
!
!
!read from .fss to get the total ncount
ncount=0
!
do 
  !*********just want to know how many steps there is******************
  read(1,*,iostat=ierror) step, dummyR
  if (ierror .lt. 0) then
    exit
  endif
  do i=1,2
      read(1,*)
  enddo
  ncount=ncount+1
enddo
close(1)
open(unit=1, file=(trim(filename)//'.ion'), status='old')
!
allocate(time(ncount-stepstart+1),rIO(3,ncount-stepstart+1),box(3,ncount-stepstart+1))
allocate(drp(3,ncount-stepstart+1),step(ncount-stepstart+1))
!***********************
!
!
!
rIO=0
time=0
box=0
drp=0
step=0
e_step=0
!
!
!
!***********************
!Read .error file
!ncount=0
!read(2,*)
!do 
!  ncount=ncount+1
!  read(2,*,iostat=ierror) e_step(ncount), dummyR, e_iiO(ncount)
!  if (ierror .lt. 0) exit
!enddo
!
!
! 
!***********************
!Read .ion file
!
!********
!Jump out the first few lines
if (stepstart .gt. 1) then
  do k=1, stepstart-1
    do i=1,3
      read(1,*)
    enddo
  enddo
endif
!********
!
ncount=0
k0=1
do 
  read(1,*,iostat=ierror) step(ncount+1), time(ncount+1), numION
  !
  if (ierror .lt. 0) then
    !Added Dec 2014
    CALL EOF_PrintOut(ierror,ncount)
    exit
  endif
  !
  ncount=ncount+1
  read(1,*) iiO, iiH
  read(1,*) rIO(1:3,ncount) 
  !
  do j=1,3
    if (rIO(j,ncount) .gt. 0) rIO(j,ncount)=rIO(j,ncount)-int(rIO(j,ncount)/celldm)*celldm
    if (rIO(j,ncount) .lt. 0) rIO(j,ncount)=rIO(j,ncount)-(int(rIO(j,ncount)/celldm)-1)*celldm
  enddo
  !
  if ((step(ncount)-step(ncount-1)) .ne. dstep) then
    write(*,'("Find dstep change at step=",I7,"from",I6,"to",I6, "at time", F8.4, "step", I7)') & 
         ncount, dstep,(step(ncount)-step(ncount-1)), time(ncount), step(ncount)
    dstep=step(ncount)-step(ncount-1)
    k0=ncount
  endif
  !
enddo
!
write(*,'(3x,"ncount=",I7,5x,"k0=",I3)') ncount,k0
write(*,*)
!
!
!
!
write(*,'("MSD is calculated during time of:",2F8.3)') time(k0),time(ncount)
!
!
!
!
!
!************CORE CONTENTS***************
!Start MSD calculation
!do dk=1,ncount-k0
do dk=1,last_dk
  !
  n=0
  average_box=0
  tot_sd=0
  !
  do k=k0,ncount-dk,5
    !Added Dec 2014 
    l1=1
    do
      if (step_trans(l1) .le. k) then
        l1+l1+1
      else
        exit
      endif
    enddo
    !
    do j=1,3
      !
      tmp_dr=rIO(j,k+dk)-rIO(j,k)
      dr(j)=tmp_dr+box(j,k)*celldm
      !
      !
      !Adjust the parameter: box
      if (abs(dr(j)-drp(j,k)) .gt. cutoff1 ) then
        if (dr(j)-drp(j,k) .gt. 0) then
          box(j,k)=box(j,k)-1
        else 
          box(j,k)=box(j,k)+1
        endif
        dr(j)=tmp_dr+box(j,k)*celldm
      endif
      !
      !
      !Situations that can go wrong
      if (abs(dr(j)-drp(j,k)) .gt. cutoff2 ) then
        !if (                                &
        !    step(k+dk) .ne. e_step(1) .and. &
        !    step(k+dk) .ne. e_step(2) .and. &
        !    step(k+dk) .ne. e_step(3) .and. &
        !    step(k+dk) .ne. e_step(4) .and. &
        !    step(k+dk) .ne. e_step(5)       &
        !   ) then
        !  !Error! dr-drp too big!
          write(4,*) step(k+dk),time(k+dk),abs(dr(j)-drp(j,k)),j
        !endif
      endif
      !
      !
      !
      average_box(j)=average_box(j)+box(j,k)
      drp(j,k)=dr(j)
      !Added Dec 2014
      if ((k+dk) .eq. step_trans(l1)) then
        dr_tot_j(j)=dr_tot(j)
      else 
      !
    enddo
    !
    d2=(dr(1))**2+(dr(2))**2+(dr(3))**2
    !Added Dec 2014
    if ((k+dk) .eq. step_trans(l1)) then
      tot_sd_j=tot_sd_j+d2_j
      n1=n1+1
      l1=l1+1
    else
      tot_sd_v=tot_sd_v+d2_v
      n2=n2+1
    endif
    !
    tot_sd=tot_sd+d2
    n=n+1
    !
    !
    !
    !##### Checking Errors #####
    !if (mod(dk,1000) .eq. 0) write(5,'(2I8,5x,3F10.5,5x,3I3)') k, dk, rIO(1:3,k+dk),box(1:3,k)
    !##### Checking Errors Completes#####
    !
    !
    !
  enddo
  !
  msd=tot_sd*convertBA*convertBA/real(n)
  write(3,'(2F8.3)') time(dk)-time(k0), msd
  !
enddo
!
!
!
write(*,*) "Program ends!"
write(*,'("Final box adjustments are:",3F8.3)') real(average_box(1))/n,real(average_box(2))/n,real(average_box(3))/n
!
!
!
end subroutine main


  !Added Dec 2014
  SUBROUTINE EOF_PrintOut(err,k)
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in)    :: err,k
    !
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(k)
    ! 
    return
    !
  END SUBROUTINE EOF_PrintOut
  !---------------------------------------------------
!
!******************************************************
end program mean_square_displacement
