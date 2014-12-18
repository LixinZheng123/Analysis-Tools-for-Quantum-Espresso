program hydrogen_bond
!**********************************
!
!This fortran file is to read from file .ion, .fss, 
!and calculate the average bond length from ion to first solvation shell atoms.
!Author: Lixin Zheng
!May 2014
!
implicit none
!
INTEGER, PARAMETER      :: DP = selected_real_kind(14,200) !double-precision kind
real(DP), parameter     :: pi=4.d0*atan(1.d0)
real(DP), parameter     :: celldm = 23.5170     !dimension of the supercell 
real(DP), parameter     :: convertBA= 0.52918     !convert Bohr into Angestrom
character(len=40)       :: filename
integer                 :: i,i1,i2,j,k,p,q,      &    !General index of atom,
                           step,                 &
                           stepstart, stepstop
integer                 :: nsp(2)
integer                 :: cs
integer                 :: ncount=0
integer                 :: ierror
integer                 :: num_donate,num_accept
integer                 :: AXDX(5,5,4)                    !Accept, donate, state
integer                 :: hbcase
integer                 :: h_num(2,70)
integer                 :: hb_accept(4)
real(DP)                :: time
real(DP)                :: sigma
real(DP)                :: rO(3,70), rH(3,140)
real(DP)                :: rIO(3), rIO_t(3), rIH(3,3), rIH_t(3)
real(DP)                :: delta(3,2)
real(DP)                :: d(2)
real(DP)                :: theta
real(DP)                :: rHB,r_OH
real(DP)                :: percent(5,5)
!
namelist /input/ cs, filename, rHB, r_OH, nsp, stepstart, stepstop
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
write(*,*) '|         Hydrogen Bond Analysis          |'
write(*,*) '------------------------------------------'
write(*,*) 'Input file           : ', trim(filename)//'.pos'
write(*,*) 'Output file          : ', trim(filename)//'.bond'
!
!if (cs .eq. 2) numH=1
!if (cs .eq. 3) numH=3
!
!AXDX=0
!countS=0
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
    open(unit=3, file=(trim(filename)//'.pos'), status='old')
    open(unit=9, file=(trim(filename)//'.bond'), status='unknown')
  case(2)
    close(3)
    close(9)
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
!************************************
!Jump out to stepstart
if (stepstart .gt. 1) then
  do k=1, (stepstart-1)*(nsp(1)+nsp(2)+1)
    read(3,*) 
  enddo
endif
!************************************
!main loop
do
  !********************************************
  !Read .pos
  read(3,*,iostat=ierror) step,time
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  ncount=ncount+1
  write(9,*) step,time
  !**********************************
  !Read .pos
  !Oxygen
  do i1=1,nsp(1)
    read(3,*) rO(1:3,i1)
    !if (iO(i) .eq. iiO) then
    !  do j=1,3
    !    rIO(j)=rO(j,i)
    !  enddo
    !endif
  enddo
  !Hydrogen
  do i2=1,nsp(2)
    read(3,*) rH(1:3,i2)
    !if (iH(i) .eq. iiH(q+1) .and. q .lt. numH) then
    !  q=q+1
    !  do j=1,3
    !    rIH(j,q)=rH(j,i)
    !  enddo
    !endif
  enddo
  h_num=0
  !
  !Record the link between O and H
  do i1=1,nsp(1)
    q=0
    do i2=1,nsp(2) 
      call get_distance(rO(1:3,i1),rH(1:3,i2),d(1))
      if (d(1)*convertBA .lt. r_OH) then
        q=q+1
        if (q .gt. 2) then
          write(*,'("Error! Hydrogen number around Oxygen ",I2," equals to ",I1," at step= ", I8)') &
               i1,q,step
          exit 
        endif
        h_num(q,i1)=i2
      endif
    enddo
    if (q .lt. 2) then 
      write(*,'("Error! Hydrogen number around Oxygen ",I2,"equals to ",I1,"at step= ", I8)') &
                        i1,q,step
      exit
    endif
  enddo
  !
  !********************************************
  !
  !===============count the ACCEPTING Hydrogen Bonds=================
  do i=1,nsp(1)
    p=0
    hb_accept=0
    do i1=1,nsp(1)
      if (i .ne. i1) then
        call get_distance(rO(1:3,i),rO(1:3,i1),d(1))
        if (d(1)*convertBA .lt. rHB) then
          !do q=1,4
          do q=1,2
            !if (h_num(q,i1) .eq. 0) exit
            i2=h_num(q,i1)
            call get_distance(rO(1:3,i1),rH(1:3,i2),d(2))
            call get_theta(rO(1:3,i1),rO(1:3,i),rH(1:3,i2),d(1:2),theta)
            if (theta .le. pi/6) then 
              p=p+1
              if (p .gt. 4) then
                write(*,'("Error! HB around Oxygen ",I2,"equals to ",I1,"at step= ", I8)') &
                     i, p, step
                exit
              endif
              hb_accept(p)=i2
            endif
          enddo
        endif
      endif
    enddo
    num_accept=p
    write(9,*) i,h_num(1:2,i),hb_accept(1:4)
  enddo
enddo
end subroutine main
!
!
!
!******************************************************
!Apply periodic boundary conditions to wrap 
! coordinates around the origin
!******************************************************
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

subroutine get_theta(pos1,pos2,pos3,pdist,angle)
!Try to find the angle of line 1-2, and lind 1-3.
implicit none

real(DP), intent(in)    :: pos1(3),pos2(3),pos3(3)
real(DP), intent(in)    :: pdist(2)  !pdist(1) is r(1-2),pdist(2) is r(1-3)
real(DP), intent(inout) :: angle
real(DP)                :: delta(3,2)
real(DP), parameter     :: pcell=23.5170
integer                 :: i
!
do i=1,3
  delta(i,1)=pos1(i)-pos2(i)
  delta(i,1)=delta(i,1)-nint(delta(i,1)/pcell)*pcell
  delta(i,2)=pos1(i)-pos3(i)
  delta(i,2)=delta(i,2)-nint(delta(i,2)/pcell)*pcell
enddo
angle=acos((delta(1,1)*delta(1,2)+delta(2,1)*delta(2,2)+delta(3,1)*delta(3,2))/(pdist(1)*pdist(2)))

return
end subroutine get_theta
!******************************************************************************************
end program hydrogen_bond
