program mean_square_displacement
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: celldm = 23.5170,&   !dimension of the supercell 
                           convertBA = 0.529    !Transfer Bohr to Angestrom
character(len=40)       :: filename 
integer                 :: i,j,k,       &       
                           nsp1, nsp2,  &       !number of atoms in O and H
                           nsp(2),      &
                           stepstart,stepstop
integer                 :: readstep
integer                 :: ncount=0
integer                 :: ierror
real(DP)                :: d2,          &
                           delta(3),    &       !displacement in (x,y,z) of marked Oxygen
                           sd, msd,     &       !square displacement before nomalization
                           sd_com, msd_com=0
real(DP)                :: com(3), comi(3)
real(DP)                :: rO(3,64), rOi(3,64)
real(DP)                :: rH(3,128), rHi(3,128)
real(DP)                :: amass(2), tmass
real(DP)                :: time
!
namelist /input/ filename, stepstart, nsp1, nsp2
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
write(*,*) 'Input file   : ', trim(filename)//'.pos'
write(*,*) 'Output file  : ', trim(filename)//'.msd,  ',trim(filename)//'.cm'
!
!
!Initiation
nsp(1)=nsp1
nsp(2)=nsp2
amass(1)=15.9994_DP      ! mass of Oxygen (amu) 
amass(2)=1.00794_DP      ! mass of H (amu)
tmass=amass(1)*nsp(1)+amass(2)*nsp(2)
!************************************
  
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files1
        !
        implicit none
        !
        open(unit=1, file=(trim(filename)//'.pos'), status='old')
        open(unit=3, file=(trim(filename)//'.msd'), status='unknown')
        open(unit=4, file=(trim(filename)//'.cm'), status='unknown')
        !
end subroutine files1

subroutine files2
        implicit none
        close(1)
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
!Jump out the first few lines
if (stepstart .gt. 1) then
  do k=1, stepstart-1
    do i=1,nsp(1)+nsp(2)+1
      read(1,*)
    enddo
  enddo
endif
!***********************
do k=1,1 
  comi=0
  !
  ! Read in .pos
  !
  read(1,*) readstep, time
  !
  do i=1,nsp(1)
    read(1,*) rOi(1:3,i)
  enddo
  do i=1,nsp(2)
    read(1,*) rHi(1:3,i)
  enddo
  !
  do j=1,3
    ! Compute center of mass
    do i=1, nsp(1)
      comi(j)=comi(j)+amass(1)*rOi(j,i)
    enddo
    do i=1, nsp(2)
      comi(j)=comi(j)+amass(2)*rHi(j,i)
    enddo
    comi(j)=comi(j)/tmass
    ! Comput the coordinate under the center of mass reference
    do i=1, nsp(1)
      rOi(j,i)=rOi(j,i)-comi(j)
    enddo
    do i=1, nsp(2)
      rHi(j,i)=rHi(j,i)-comi(j)
    enddo
  enddo
  !write(4,*) readstep, msd_com, comi(1:3)
enddo
!
!
do 
  read(1,*,iostat=ierror) readstep, time
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif

  com=0
  sd=0
  sd_com=0
  !
  ! Read in .pos
  !
  !
  do i=1,nsp(1)
    read(1,*) rO(1:3,i)
  enddo
  do i=1,nsp(2)
    read(1,*) rH(1:3,i)
  enddo
  !
  do j=1,3
    ! Compute center of mass
    do i=1, nsp(1)
      com(j)=com(j)+amass(1)*rO(j,i)
    enddo
    do i=1, nsp(2)
      com(j)=com(j)+amass(2)*rH(j,i)
    enddo
    com(j)=com(j)/tmass
    ! Comput the coordinate under the center of mass reference
    do i=1, nsp(1)
      rO(j,i)=rO(j,i)-com(j)
    enddo
    do i=1, nsp(2)
      rH(j,i)=rH(j,i)-com(j)
    enddo
  enddo
  !
  do i=1, nsp(1)
    call get_distance(rO(1:3,i),rOi(1:3,i),d2)
    sd=sd+d2
  enddo
  call get_distance(com(1:3),comi(1:3),d2)
  sd_com=sd_com+d2
  !
  msd=sd*(convertBA**2)/nsp(1)
  msd_com=sd_com*(convertBA**2)
  write(3,*) time, msd
  !write(4,*) readstep, msd_com, com(1:3)
  !
  ncount=ncount+1
  !
enddo
!
end subroutine main
!
!******************************************************
!******************************************************
subroutine get_distance(pos1,pos2,dist2)

implicit none
!
real(DP), intent(in)    :: pos1(3),pos2(3)
real(DP), intent(inout) :: dist2
real(DP)                :: delta(3)
real(DP), parameter     :: pcell=23.5170
integer                 :: i
!
do i=1,3
  delta(i)=pos1(i)-pos2(i)
  delta(i)=delta(i)-nint(delta(i)/pcell)*pcell
enddo
!
dist2=delta(1)**2+delta(2)**2+delta(3)**2
!
return
end subroutine get_distance
!******************************************************************************************
end program mean_square_displacement
