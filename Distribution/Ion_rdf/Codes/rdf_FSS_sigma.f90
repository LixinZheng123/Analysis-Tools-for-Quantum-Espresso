!Radial Distribution Function g(r)O* for Quantum Espresso Output Filesres
!Calculates RDF around a certain O*.
!==================================================================================================
program rdf 
   !
   implicit none
   !
   INTEGER, PARAMETER   :: DP = selected_real_kind(14,200) !double-precision kind
   real(DP), parameter  :: pi=4.d0*atan(1.d0),&
                           celldm=23.5170*0.3,&
                           convertBA = 0.52917
   character(len=1)     :: rs,ch_atom(2)
   integer              :: binNum,        &  !number of bins for g(r)
                           nsp(2),        &  !number of each atomic species     
                           atom(2),       &
                           cs,            &  !cs=2, OH-; cs=3, H3O+
                           stepstart=1,     &  !nfi that we will start reading *.pos 
                           stepstop=-1,      &  !nfi that we will stop reading *.pos AFTER (optional)
                           ncount=0,      &  !number of steps index
                           i, j, k, ii=0, &  !general indexes
                           n,m,l,p,       &  !general index 
                           iO,iH, iiH_t,  &  !general index 
                           iiO,iiOp=0,    &  !index of O* read from Ion
                           iiH(3),        &
                           numION,        &  !Number of Ions in each timestep. If numION=0, iiO=iiOp
                           numH,          &  !Number of Hydrogen in each ion. cs=2, numH=1; cs=3, numH=3
                           step, stepp,   &  !read the nfi from *.pos 
                           nbin              !bin index
   integer              :: ierror
   integer              :: state, readstate
   real(DP), allocatable:: grt(:),        &  !unnormalized g(r)
                           temgrt(:)         !temp unnormalized g(r)
   real(DP)             :: gofr,          &  !FINAL g(r) print-out variable (not an array, print-out values)  
                           sum=0,         &
                           r,             &  !total distance (not an array, print-out value)
                           d,             &
                           cell(3,3),     &  !lattice prim vectors
                           cellinv(3,3),  &  !inverse of prim vectors, for scaled positions
                           rIO(3),        &  !Ion position
                           rIH(3,3),      &  !ion positions (dimension,numH), (dimension,numH*)
                           rO(3,15),rH(3,30),&!atomic positions (dimension,atomic-species)
                           box,           &  !length of the region (includes multiple supercells)
                           res,           &  !resolution (size) of bins  binNum/Hbox (number/length)
                           omega,         &  !volume of cell
                           dens,          &  !raw particle density
                           r2,            &  !total squared distance
                           time,          &  !timestep
                           sigma,         &
                           vshell            !volume of the infinitesimal radial shell

   character(len=30)       :: filename, filegofr
   !
   namelist /input/ filename,filegofr, atom, cs, binNum, readstate, stepstart, stepstop
   !
   !initialization
   call init
   !
   !Open files
   call files(1)
   !
   !Main Loop
   call main
   !
   !finalize the results
   call final
   !
   !close files
   call files(2)
   !
   !
   !********************************************************************************
   contains
   !********************************************************************************
      !
      !
subroutine init
!
implicit none
!
!read Namelist input for stdin
read(*,input)
!
!
write(rs,'(i1)') readstate
do i=1,2
  if (atom(i) .eq. 1) ch_atom(i)='O'
  if (atom(i) .eq. 2) ch_atom(i)='H'
enddo
!
!Print Intro
write(*,*) '------------------------------'
write(*,fmt='("|          g(r) of ",A1,"*",A1,"        |")') ch_atom(1), ch_atom(2)
write(*,*) '------------------------------'
write(*,*) '1st solvation file : ', trim(filename)//'.fss'
write(*,*) 'Sigma file         : ', trim(filename)//'.sigma'
if (readstate .eq. 0) then
  write(*,*) 'Output file        : ', trim(filegofr)
else
  write(*,*) 'Output file        : ', trim(filegofr)//'.'//trim(rs)
endif
write(*,fmt='(1X, "Bin number          : ", I8)') binNum
!
!
if (cs .eq. 2) numH=1
if (cs .eq. 3) numH=3
!
!Allocate  Files
allocate(grt(binNum), temgrt(binNum))
!
grt=0
!******************************************************
!calculate the inverse
cell=0
cell(1,1)=celldm
cell(2,2)=celldm
cell(3,3)=celldm
call invert(cell, cellinv, omega)
box=omega**(1.0d0/3.0d0)
res=binNum/box
!******************************************************
end subroutine init
!
!******************************************************
subroutine files(io)
!
implicit none
!
integer, intent(in)  ::io
!
select case(io)
  case(1)
    open(unit=1, file=(trim(filename)//'.fss'), status='old')
    open(unit=2, file=(trim(filename)//'.sigma'), status='old')
    if (readstate .eq. 0) then
      open(unit=3, file=(trim(filegofr)), status='unknown')
    else
      open(unit=3, file=(trim(filegofr)//'.'//trim(rs)), status='unknown')
    endif
    !
  case(2)
    close(1)
    close(2)
    close(3)
  case default
end select
!
end subroutine files
!
!
!
!******************************************************
!------------------------------------------------------
!Main Loop
! read-in data, calculate inverse, r_to_s,
! 
! 
!------------------------------------------------------
!******************************************************
subroutine main
!
implicit none
!
!Jump out the first few lines
if (stepstart .gt. 1) then
  do k=1, stepstart-1
    read(1,*) step, time, numION, nsp(1), nsp(2)
    do i=1,nsp(1)+nsp(2)+1
      read(1,*)
    enddo
    read(2,*)
  enddo
endif
!
do 
  !******************************************************
  read(1,*,iostat=ierror) step, time, numION, nsp(1), nsp(2)
  if (ncount .gt. stepstop .and. stepstop .ne. -1) exit
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !
  read(2,*,iostat=ierror) stepp, time, iiO, iiH_t, sigma, state
  if (ierror .lt. 0) then
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
    exit
  endif
  !
  if (stepp .ne. step) then
    write(*,*) "Error! steps not equal from .fss and .sigma!"
  endif
  if (state .ne. readstate .and. readstate .ne. 0) then
    do i=1,nsp(1)+nsp(2)+1
      read(1,*)
    enddo
    cycle
  endif
  !
  !******************************************************
  !Read in .fss
  read(1,*) iiO, iiH(1:numH)
  do i=1, nsp(1)
    read(1,*)  iO,rO(1:3,i)
    if (iO .eq. iiO) then
      do j=1,3 
        rIO(j)=rO(j,i)
      enddo
      ii=i
    endif
  enddo
  p=1
  do i=1,nsp(2)
    read(1,*)  iH,rH(1:3,i)
    if (p .le. numH) then
      if (iH .eq. iiH(p)) then 
        do j=1,3
          rIH(j,p)=rH(j,i)
        enddo
        p=p+1
      endif
    endif
  enddo
  !******************************************************
  !Zero all temp g(r)
  temgrt = 0
  !******************************************************
  !------------------Count PAIRS-----------------------
  !================first situation=====================
  if (atom(1) .eq. 1 .and. atom(2) .eq. 1) then
    !***************IMPORTANT*******************
    !because in each step new_nsp1 is different, we need to partly normalize it now.
    dens = DBLE(nsp(1)-1)/omega
    !
    do i=1,nsp(1)
      call get_distance(rIO(1:3),rO(1:3,i),d)
      !if d=0
      if (d .lt. 0.001) cycle
      !establish which bin this distance is in
      nbin = nint(d*res)
      !if we are at the end of binNum
      if (nbin .gt. binNum) nbin = binNum
      !update the temp g(r)
      temgrt(nbin) = temgrt(nbin) + 1/dens
    enddo
  endif
  !================second situation=====================
  if (atom(1) .eq. 1 .and. atom(2) .eq. 2) then
    dens = DBLE(nsp(2))/omega
    do i=1,nsp(2)
      call get_distance(rIO(1:3),rH(1:3,i),d)
      nbin = nint(d*res)
      if (nbin .gt. binNum) nbin = binNum
      temgrt(nbin) = temgrt(nbin) + 1/dens
    enddo
  endif
  !================third situation=====================
  if (atom(1) .eq. 2 .and. atom(2) .eq. 1) then
    dens = DBLE(nsp(1))/omega
    do i=1,nsp(1)
      if (ii .eq. i) cycle
      do p=1, numH
        call get_distance(rIH(1:3,p),rO(1:3,i),d)
        nbin = nint(d*res)
        if (nbin .gt. binNum) nbin = binNum
      enddo
      temgrt(nbin) = temgrt(nbin) + 1/dens
    enddo
  endif
  !================forth situation=====================
  !if (atom1 .eq. 2 .and. atom2 .eq. 2) then
  !not considering it now
  !endif
  !===================================================
  do n=1,binNum,1
    grt(n)=grt(n)+temgrt(n)
  enddo
  ncount=ncount + 1
enddo 
!
write(*,*) ' ... Main loop completed!'
!
end subroutine main
!
!
!
!******************************************************
!------------------------------------------------------
!Finalize the results, print out to the output
!
!------------------------------------------------------
!******************************************************
subroutine final
!
implicit none
!
do n=1,binNum-1,1
  !
  !construct some details of normalization and distance  
  r = DBLE(n)/res
  r2 = r**2.0d0
  vshell = 4.0*pi*r2/res
  !Calculate the final, normalized, value of g(r)! Please note that the
  !normalization constant (ncount*mic3*nsp(atom1)*dens*vshell)...
  !  ncount*mic3 = number of steps and number of extended shells
  !  nsp(atom1)*dens = number of pairs
  !  vshell = volume of the infinitesimal shell
  if (atom(1) .eq. 1) gofr = grt(n)/(DBLE(ncount*vshell))
  if (atom(1) .eq. 2) gofr = grt(n)/(DBLE(ncount*vshell*numH))
  sum=sum+gofr
  !
  write(3,*) (r*convertBA), gofr
  !
enddo
write(*,*) "sum of all gofr is :", sum
!
write(*,*) ' Program complete!'
!
end subroutine final
! 
!
!
!******************************************************
!******************************************************
!         
subroutine invert(Mi, Mo, det)

implicit none

real(8), intent(in)     :: Mi(3,3)
real(8), intent(out)    :: Mo(3,3), det
real(8)                 :: tmp(3,3), vs

integer              :: vi,vj,vk,vl  !indexes
integer              :: vn,vr     !int counters

det=0.0
vs   = 1.0
vi   = 1
vj   = 2
vk   = 3


do
  do vn=1,3,1
    det= det + vs*Mi(1,vi)*Mi(2,vj)*Mi(3,vk)
    vl=vi
    vi=vj
    vj=vk
    vk=vl
  enddo
  vi=2
  vj=1
  vk=3
  vs=-vs
  if (vs .GT. 0.0) exit
enddo

if(ABS(det) .LT. 1.0e-20) then
  write(*,*) 'Error: Singular Matrix'
  Stop
endif

vi=1
vj=2
vk=3

do vr=1,3
  tmp(vr,1) = (Mi(2,vj)*Mi(3,vk) - Mi(2,vk)*Mi(3,vj))/det
  tmp(vr,2) = (Mi(3,vj)*Mi(1,vk) - Mi(3,vk)*Mi(1,vj))/det
  tmp(vr,3) = (Mi(1,vj)*Mi(2,vk) - Mi(1,vk)*Mi(2,vj))/det
  vl=vi
  vi=vj
  vj=vk
  vk=vl
enddo

do vl=1,3
  do vk=1,3
    Mo(vk,vl)=tmp(vk,vl)
  enddo
enddo

end subroutine invert
!
!******************************************************
!******************************************************
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
!
dist=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
!
return
end subroutine get_distance
!******************************************************************************************
end PROGRAM rdf 
