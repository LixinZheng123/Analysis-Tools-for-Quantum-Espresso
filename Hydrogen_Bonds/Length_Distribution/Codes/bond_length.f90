program bond_length
  !**********************************
!
!This fortran file is to read from file .ion, .fss, 
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
                           iiO=0,       &       !the O* index
                           iFSS(60),    &       !the index of first solvation shell (FSS)
                           nsp1, nsp2,  &       !number of atoms in O and H
                           step,        &
                           stepstart, stepstop
integer                 :: norm(2)              !counter for the steps 
                                                !1: when not transiting ions; 2: when transiting ions
integer                 :: cs
integer                 :: state                !state=1: stable; state=2: transitioon
integer                 :: numION,numH
integer                 :: ncount=0
integer                 :: ierror
integer, allocatable    :: iiH(:)
real                    :: time
real                    :: r(3,60)              !the position for FSS atoms
real                    :: delta(3)
real                    :: bond(2)              !the calculated bond length
real                    :: sum(2)               !for the bond length when not transiting ions
real                    :: per                  !percentage of stable ions among all steps.
real                    :: r_ION(3)             !the ion position
real                    :: d
!
namelist /input/ cs, filename, nsp1, nsp2, stepstart, stepstop
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
write(*,*) 'First Shell file : ', trim(filename)//'.fss'
!
if (cs .eq. 2) numH=1
if (cs .eq. 3) numH=3
!allocation
allocate(iiH(numH))
!
sum(1)=0
sum(2)=0
norm(1)=0
norm(2)=0
!
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files1
implicit none
open(unit=3, file=(trim(filename)//'.fss'), status='old')
end subroutine files1
!
subroutine files2
implicit none
close(3)
end subroutine files2 
!
!******************************************************************************************
subroutine main
!
implicit none
!
!************************************
!Jump out to stepstart
if (stepstart .gt. 1) then
  do k=1, stepstart-1
    read(3,*) step,time,numION,iiO,p
    if (p .gt. 0) then
      do i=1,p
        read(3,*)
      enddo
    endif
  enddo
endif
!************************************
!main loop
!do k =stepstart, stepstop
do
  read(3,*,iostat=ierror) step,time,numION,iiO,p
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !**********************************
  !read .ion file
  !read(3,*) step,time,numION,iiO,p
  !### judge if during transition or on stable state ###
  if (numION .eq. 1) state=1
  if (numION .ne. 1) state=2
  !**********************************
  if (p .gt. 0) then
    do i=1,p
      read(3,*) iFSS(i), r(1:3,i)
      if (iFSS(i) .eq. iiO) then
        do j=1,3
          r_ION(j)=r(j,i)
        enddo
      endif
    enddo
    !*****************CORE CONTENT***********************
    !find the average bond length
    do i=1,p
      if (iFSS(i) .ge. nsp1) cycle
      call get_distance(r(1:3,i), r_ION(1:3), d)
      sum(state)=sum(state)+d
      norm(state)=norm(state)+1
    enddo
  endif
  ncount=ncount+1
enddo
bond(1)=convertBA*sum(1)/norm(1)
bond(2)=convertBA*sum(2)/norm(2)
per=norm(1)/real((norm(1)+norm(2)))*100
!need to add percentage calculation
write(*,*) "First Solvation Shell bond length when at stable state is:    ", bond(1)
write(*,*) "First Solvation Shell bond length when at transition state is:", bond(2)
write(*,*) norm(1:2)
write(*,*) "The percent at stable state is:", per
end subroutine main
!
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
