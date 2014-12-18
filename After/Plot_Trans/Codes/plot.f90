program plot
!
!This fortran file is to read from file .trans, and turn it into gnuplot-friendly format to generate graphs.
!Author: Lixin Zheng
!July 2014
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP),PARAMETER      :: dt = 8.466D-5        !3.5*2.4188843D-5 (ps)
character(len=40)       :: filename
integer                 :: i,j,k,p,q,      &    !General index of atom,
                           step
integer                 :: ierror=0
integer                 :: iiO=0, iiO_p=0
integer                 :: iiH(3), iiH_p(3)
integer                 :: ncount,cs
!
namelist /input/ filename, cs
!
!initialization
call init
!
!Open files
call files(1)
!
call main
!
call files(2)
!
!
! 
contains
!
!******************************************************************************************
subroutine init
!
implicit none
!read Namelist input for stdin
read(*,input)
!
end subroutine init
!******************************************************************************************
!
!******************************************************************************************
subroutine files(io)
!
implicit none
!
integer,intent(in) :: io
!
select case(io)
  case(1)
    open(unit=3, file=(trim(filename)//'.trans'), status='old')
    open(unit=4, file=('trans_'//trim(filename)//'.dat'), status='unknown')
  case(2)
    close(3)
    close(4)
  case default
end select

end subroutine files
!
!******************************************************************************************
!
!******************************************************************************************
subroutine main
!
implicit none
!
!main loop
do
  !======
  !Read .trans
  read(3,*,iostat=ierror) step,iiO
  if (ierror .lt. 0) then
    !write(*,*) '  End of File Reached'
    exit
  endif
  !======
  if(iiO_p .ne. 0) write(4,*) step*dt-0.00001, iiO_p
  write(4,*) step*dt, iiO
  iiO_p=iiO
  ncount=ncount+1
enddo

end subroutine main
!******************************************************************************************
!
!******************************************************************************************
end program plot
