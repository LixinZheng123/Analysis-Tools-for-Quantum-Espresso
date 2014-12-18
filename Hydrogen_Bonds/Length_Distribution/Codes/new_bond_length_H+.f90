program bond_length
  !**********************************
!
!This fortran file is to read from file .fss, 
!and calculate the average bond length from ion to first solvation shell atoms.
!Author: Lixin Zheng
!May 2014
!
implicit none
!
real, parameter         :: celldm = 23.5170     !dimension of the supercell 
real, parameter         :: convertBA= 0.52918     !convert Bohr into Angestrom
character(len=40)       :: filename
integer                 :: i,j,k,p,q,   &       !General index of atom,
                           iiO,iiH,     &       !the O* index
                           nsp(2),      &       !number of atoms in O and H
                           step,        &
                           stepstart, stepstop
integer                 :: step_t=0,iiO_t      !the step and oxygen index number from .trans file
integer                 :: total(4)
integer                 :: cs
integer                 :: state                !state=1: stable, ready to transfer; 
                                                !state=2: transfer; 
                                                !state=3: middle of transfer
                                                !state=4: not transfering;
integer                 :: numION,numH
integer                 :: ncount=0
integer                 :: ierror
integer                 :: iO(15),iH(30),hnum(4),hnum_t(4)
real                    :: time
real                    :: rO(3,15),rH(3,30),rIO(3),rIH(3)
real                    :: rIO_t(3),rIH_t(3)
real                    :: delta(3)
real                    :: bond(2)              !the calculated bond length
real                    :: sum(2)               !for the bond length when not transiting ions
real                    :: per                  !percentage of stable ions among all steps.
real                    :: sigma
real                    :: d(4)
real                    :: r_OH
integer                 :: iiH_H(3),iiH_t(3)
!
namelist /input/ cs, filename, stepstart, stepstop, r_OH
!
!initialization
call init
!
!Open files
call files1
!
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
!read Namelist input for stdin
read(*,input)
!
!print Intro
write(*,*) '------------------------------------------'
write(*,*) '|          Find the Bond Length           |'
write(*,*) '------------------------------------------'
write(*,*) '1st solvation file : ', trim(filename)//'.fss'
write(*,*) 'Ion transfer file  : ', trim(filename)//'.trans'
write(*,*) 'Sigma file         : ', trim(filename)//'.sigma'
!
if (cs .eq. 2) numH=1
if (cs .eq. 3) numH=3
!
sum=0
total=0
!
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files1
implicit none
open(unit=3, file=(trim(filename)//'.fss'), status='old')
open(unit=4, file=(trim(filename)//'.trans'), status='old')
open(unit=5, file=(trim(filename)//'.sigma'), status='unknown')
end subroutine files1
!
subroutine files2
implicit none
close(3)
close(4)
close(5)
end subroutine files2 
!
!******************************************************************************************
subroutine main
!
implicit none
!
!************************************
read(4,*) step_t,iiO_t
!Jump out to stepstart
if (stepstart .gt. 1) then
  do k=1, stepstart-1
    read(3,*) step,time,numION,nsp(1),nsp(2)
    read(3,*)
    if (p .gt. 0) then
      do i=1,nsp(1)+nsp(2)
        read(3,*)
      enddo
    endif
  enddo
endif
!************************************
!main loop
do
  read(3,*,iostat=ierror) step,time,numION,nsp(1),nsp(2)
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !**********************************
  if (step_t .le. step+10) read(4,*,iostat=ierror) step_t,iiO_t,iiH_t(1:3)
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !**********************************
  !
  if (nsp(1) .eq. 0) then
    write(*,*) 'Error! nsp1=0.'
    exit 
  endif
  !
  !**********************************
  !Read in .fss file
  !**********************************
  read(3,*) iiO,iiH_H(1:3)
  do i=1,nsp(1)
    read(3,*)  iO(i),rO(1:3,i)
    if (iO(i) .eq. iiO) then
      do j=1,3
        rIO(j)=rO(j,i)
      enddo 
      !write(*,*) "rIO:", rIO(1:3)
    else if (iO(i) .eq. iiO_t) then
      !Find the coordinate of following O*
      do j=1,3
        rIO_t(j)=rO(j,i)
      enddo
      !write(*,*) "rIO_t:", rIO_t(1:3)
      !Find the index of Hydrogen around the following O*
    endif
  enddo
  !
  do p=1,3
    do i=1,3
      if (iiH_t(p) .eq. iiH_H(i)) then
        iiH=iiH_H(i)
        exit
      endif
    enddo
  enddo
  !
  do i=1,nsp(2)
    read(3,*)  iH(i),rH(1:3,i)
    if (iH(i) .eq. iiH) then
      do j=1,3
        rIH_t(j)=rH(j,i)
      enddo
    endif
  enddo
  !
  !!**********************************
  !!Find the only hydrogen that would going back and forth between O* and following O*
  !!**********************************
  !d=0
  !!
  !do p=1,4
  !  if (hnum_t(p) .eq. 0) exit
  !  do i=1,nsp(2)
  !    if (iH(i) .eq. hnum_t(p)) then
  !      call get_distance(rIO(1:3),rH(j,i),d(p))
  !      !write(*,*) "when iH(i)=",iH(i),"distance=",d(p)
  !    endif
  !  enddo
  !enddo
  !do p=2,4
  !  if (d(p) .eq. 0) exit
  !  if (d(p) .lt. d(1)) then
  !    hnum_t(1)=hnum_t(p)
  !    d(1)=d(p)
  !  endif
  !enddo
  !do i=1,nsp(2)
  !  if (iH(i) .eq. hnum_t(1)) then
  !    do j=1,3
  !      rIH_t(j)=rH(j,i)
  !    enddo
  !    exit
  !  endif
  !enddo
  !!write(*,*) "The shortest distance is:",d(1),"with hydrogen coordinate",rIH_t(1:3)
  !
  !*****************CORE CONTENT 1***********************
  !find sigma
  call get_distance(rIO(1:3),rIH_t(1:3),d(1))
  call get_distance(rIO_t(1:3),rIH_t(1:3),d(2))
  sigma=convertBA*abs(d(2)-d(1))
  !write(*,*) "d(1) and d(2), and sigma equals to", d(1), d(2), sigma
  if (sigma .lt. 0.1) then
    state=4         !transfer state
  else if (sigma .ge. 0.1 .and. sigma .lt. 0.5) then
    state=3         !middle of transfer
  else if (sigma .ge. 0.5 .and. sigma .lt. r_OH) then
    state=2         !redy to transfer
  else
    state=1         !not transfering
  endif
  if (sigma .gt. 2.9) write(*,*) time, iiO, sigma, "error!"
  !  write(5,*) step,iiO,sigma, "error!"
  !  write(*,*) "rIO:", rIO(1:3)
  !  write(*,*) "rIO_t:", rIO_t(1:3)
  !  write(*,*) "rIH_t:", rIH_t(1:3)
  !  exit
  !if (state .eq. 1 .or. state .eq. 2) write(5,*) step,iiO,sigma
  write(5,*) step,time, iiO, iiH, sigma, state
  !*****************CORE CONTENT 2***********************
  !find the average bond length
  !do i=1,nsp(1)
  !  call get_distance(rO(1:3,i), rIO(1:3), d(1))
  !  sum(state)=sum(state)+d(1)
  !  norm(state)=norm(state)+1
  !enddo
  total(state)=total(state)+1
  ncount=ncount+1
enddo
!
!bond(1)=convertBA*sum(1)/norm(1)
!bond(2)=convertBA*sum(2)/norm(2)
!per=norm(1)/real((norm(1)+norm(2)))*100
!need to add percentage calculation
!write(*,*) "First Solvation Shell bond length when at stable state is:    ", bond(1)
!write(*,*) "First Solvation Shell bond length when at transition state is:", bond(2)
!write(*,*) norm(1:2)
!write(*,*) "The percent at stable state is:", per
write(*,*) "Percentage of each state (ready, transfer, middle, not transfer):"
write(*,*) real(total(1))/real(ncount),real(total(2))/real(ncount),real(total(3))/real(ncount),real(total(4))/real(ncount)

end subroutine main
!
!******************************************************
!******************************************************
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
!******************************************************************************************
end program bond_length
