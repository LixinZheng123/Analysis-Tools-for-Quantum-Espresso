!Radial Distribution Function g(r) for Quantum Espresso Output Filesres
!Calculates RDF around a certain O*.
!Modifide 20141006, add hydrogen bond classification
!Modified 20141031, refactor with ion_wfc_rdf.f90
!==================================================================================================
program rdf 
   !
   implicit none
   !
   !======Global variables======
   character(len=30)    :: filename
   integer, parameter   :: DP = selected_real_kind(14,200) !double-precision kind
   real(DP), parameter  :: pi=4.d0*atan(1.d0)
   real(DP), parameter  :: convertBA = 0.52917
   integer              :: nsp(2)            !number of each atomic species     
   integer              :: cs                !cs=2, OH-; cs=3, H3O+
   integer              :: stepstart=1       !nfi that we will start reading *.pos 
   integer              :: stepstop=-1       !nfi that we will stop reading *.pos AFTER (optional)
   integer              :: ncount=0          !number of steps index
   integer              :: step 
   integer              :: ierror
   real(DP)             :: time
   !======Variables concerning RDF======
   character(len=1)     :: surfix,ch_atom(2)
   character(len=30)    :: filegofr
   integer              :: nbin              !bin index
   integer              :: binNum            !number of bins for g(r)
   integer              :: atom(2)
   real(DP), allocatable:: grt(:)            !unnormalized g(r)
   real(DP), allocatable:: temgrt(:)         !temp unnormalized g(r)
   real(DP)             :: r                 !total distance (not an array, print-out value)
   real(DP)             :: gofr              !FINAL g(r) print-out variable (not an array, print-out values)  
   real(DP)             :: box               !length of the region (includes multiple supercells)
   real(DP)             :: res               !resolution (size) of bins  binNum/Hbox (number/length)
   real(DP)             :: omega             !volume of cell
   real(DP)             :: dens              !raw particle density
   real(DP)             :: vshell            !volume of the infinitesimal radial shell
   real(DP)             :: tot               !The add-up value of g(r)s of each r
   !======Variables concerning .fss======
   real(DP)             :: rO(3,15),rH(3,30) !atomic positions (dimension,atomic-species)
   real(DP)             :: rIO(3)            !Ion position
   real(DP)             :: rIH(3,3)          !ion positions (dimension,numH), (dimension,numH*)
   integer              :: numION            !Number of Ions in each timestep. If numION=0, iIO=iIOp
   integer              :: numH              !Number of Hydrogen in each ion. cs=2, numH=1; cs=3, numH=3
   integer              :: iO,iH             !general index 
   integer              :: iIO               !index of Ion oxygen
   integer              :: iIH(3)            !index of Ion hydrogens
   !======Variables concerning .hbcase======
   integer              :: readcase, hbnum
   !======Local variables======
   integer              :: ncount0=0         !number of total steps index
   integer              :: i, j, k, ii       !general indexes
   integer              :: n,m,l,p           !general indexes
   integer              :: stepp 
   integer              :: ctn
   real(DP)             :: d                 
   !
   !
   !
   namelist /input/ filename,filegofr, atom, cs, binNum, hbnum, stepstart, stepstop
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
write(surfix,'(i1)') hbnum
do i=1,2
  if (atom(i) .eq. 1) ch_atom(i)='O'
  if (atom(i) .eq. 2) ch_atom(i)='H'
enddo
!
!Print Intro
write(*,*) '------------------------------'
write(*,fmt='("|          g(r) of ",A1,"*",A1,"        |")') ch_atom(1), ch_atom(2)
write(*,*) '------------------------------'
write(*,*) 'Input file : ', trim(filename)//'.fss, ',trim(filename)//'.hbcase'
if (hbnum .eq. 0) then
  write(*,*) 'Output file        : ', trim(filegofr)
else
  write(*,*) 'Output file        : ', trim(filegofr)//'.'//trim(surfix)
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
!Initialization
grt=0
ii=0
box=4
omega=box*box*box
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
    open(unit=2, file=(trim(filename)//'.hbcase'), status='old')
    if (hbnum .eq. 0) then
      open(unit=3, file=(trim(filegofr)), status='unknown')
    else
      open(unit=3, file=(trim(filegofr)//'.'//trim(surfix)), status='unknown')
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
  !Read .fss
  read(1,*,iostat=ierror) step, time, numION, nsp(1), nsp(2)
  if (ncount .gt. stepstop .and. stepstop .ne. -1) exit
  if (ierror .lt. 0) then
    call EOF_PrintOut(ierror,hbnum,ncount,ncount0)
    exit
  endif
  !
  !Read .hbcase
  read(2,*,iostat=ierror) stepp, time, readcase
  if (ierror .lt. 0) then
    call EOF_PrintOut(ierror,hbnum,ncount,ncount0)
    exit
  endif
  read(2,*)
  read(2,*)
  !
  if (stepp .ne. step) then
    write(*,*) "Error! steps not equal in 2 files (.fss and .hbcase)!", step,stepp,time
    exit
  endif
  !
  ncount0=ncount0+1
  !******************************************************
  !Read .fss
  read(1,*) iIO, iIH(1:numH)
  do i=1, nsp(1)
    read(1,*)  iO,rO(1:3,i)
    if (iO .eq. iIO) then
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
      if (iH .eq. iIH(p)) then 
        do j=1,3
          rIH(j,p)=rH(j,i)
        enddo
        p=p+1
      endif
    endif
  enddo
  !******************************************************
  if (cs .eq. 2 ) then
    call Continue_Judgement(hbnum,readcase,ctn)
    if (ctn .eq. 0) cycle
  endif
  !******************************************************
  !Zero all temp g(r)
  temgrt = 0
  !******************************************************
  !------------------Count PAIRS-----------------------
  !================first situation=====================
  if (atom(1) .eq. 1 .and. atom(2) .eq. 1) then
    !***************IMPORTANT*******************
    !because in each step nsp(1) is different, we need to partly normalize it now.
    dens = DBLE(nsp(1)-1)/omega
    !
    do i=1,nsp(1)
      call get_distance(rIO(1:3),rO(1:3,i),d)
      !if d=0
      if (d .lt. 0.001) cycle
      nbin = nint(d*res)
      if (nbin .gt. binNum) nbin = binNum
      !update the temp g(r)
      temgrt(nbin) = temgrt(nbin) + DBLE(1)/dens
    enddo
  endif
  !================second situation=====================
  if (atom(1) .eq. 1 .and. atom(2) .eq. 2) then
    dens = DBLE(nsp(2))/omega
    do i=1,nsp(2)
      call get_distance(rIO(1:3),rH(1:3,i),d)
      nbin = nint(d*res)
      if (nbin .gt. binNum) nbin = binNum
      temgrt(nbin) = temgrt(nbin) + DBLE(1)/dens
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
  ncount=ncount+1
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
  !  vshell = volume of the infinitesimal shell
  vshell = 4.0*pi*(r**2)*(1/res) 
  !
  !
  !Calculate the final, normalized, value of g(r)! Please note that the
  !normalization constant (ncount*mic3*nsp(atom1)*dens*vshell)...
  !  ncount*mic3 = number of steps and number of extended shells
  !  nsp(atom1)*dens = number of pairs
  if (atom(1) .eq. 1) gofr = grt(n)/(DBLE(ncount*vshell))
  if (atom(1) .eq. 2) gofr = grt(n)/(DBLE(ncount*vshell*numH))
  tot=tot+gofr
  !
  write(3,*) (r*convertBA), gofr
  !
enddo
write(*,*) "The integrated result of curve is :", tot*(1/res)*convertBA
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
!******************************************************************************************
!
subroutine EOF_PrintOut(err,switch,k1,k2)

implicit none
!
integer, intent(in)    :: err,switch,k1,k2
!
write(*,*) '  End of File Reached'
write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(k1)
if (switch .ne. 0) write(*,*) '  Percentage calculated:', real(k1)/k2
!
return
end subroutine EOF_PrintOut
!******************************************************************************************
!******************************************************************************************
!
subroutine Continue_Judgement(hb,rd_hb,continu)

implicit none
!
integer, intent(in)    :: hb      !Hydrogen number we want to analyse
integer, intent(in)    :: rd_hb   !Hydrogen number we read
integer, intent(out)   :: continu
!
!======
!These lines below is copied from hydrogen_bond.f90 for reminder.
!if (num_accept .eq. 3 .and. num_donate .eq. 0) hbcase=1
!if (num_accept .eq. 3 .and. num_donate .eq. 1) hbcase=2
!if (num_accept .eq. 4 .and. num_donate .eq. 0) hbcase=3
!if (num_accept .eq. 4 .and. num_donate .eq. 1) hbcase=4
!Reminder ends
!======
!
continu=0
!
if (hb .eq. 3) then
  if (rd_hb .eq. 1 .or. rd_hb .eq. 2) then
    continu=1
  endif
endif
if (hb .eq. 4) then
  if (rd_hb .eq. 3 .or. rd_hb .eq. 4) then
    continu=1
  endif
endif
!
return
end subroutine Continue_Judgement
!******************************************************************************************
!******************************************************************************************
!

!
end PROGRAM rdf 
