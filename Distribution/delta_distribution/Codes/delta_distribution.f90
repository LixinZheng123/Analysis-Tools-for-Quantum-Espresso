!It is in fact a 1D rdf
!Output file: 
!delta;  percentage;  average # donating bond.
!Created 20141103, from rdf_FSS_HB_OO_Only.f90
!==================================================================================================
program delta
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
   !integer              :: atom(2)
   integer(DP), allocatable:: PairNum(:)     
   real(DP)             :: r                 !total distance (not an array, print-out value)
   real(DP)             :: gofr              !FINAL g(r) print-out variable (not an array, print-out values)  
   real(DP)             :: box               !length of the region (includes multiple supercells)
   !real(DP)             :: res               !resolution (size) of bins  binNum/Hbox (number/length)
   !real(DP)             :: omega             !volume of cell
   real(DP)             :: dens              !raw particle density
   real(DP)             :: vshell            !volume of the infinitesimal radial shell
   real(DP)             :: tot=0             !The add-up value of g(r)s of each r
   !======Variables concerning .hbcase======
   integer              :: readcase, hbnum
   !======Variables concerning donated hb distribution======
   real(DP), allocatable:: tot_dhb(:)        !unnormalized donated hydrogen bond for each delta
   !======Local variables======
   integer              :: ncount0=0         !number of total steps index
   integer              :: i, j, k, ii       !general indexes
   integer              :: n,m,l,p           !general indexes
   integer              :: stepp 
   integer              :: ctn
   real(DP)             :: d                 
   integer              :: dummy
   integer              :: c=0
   !
   !
   !
   namelist /input/ filename,filegofr, cs, binNum, hbnum, stepstart, stepstop
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
!
!Print Intro
write(*,*) '------------------------------'
write(*,fmt='("|          g(r) of delta     |")')
write(*,*) '------------------------------'
write(*,*) 'Input file : ', trim(filename)//'.sigma, ',trim(filename)//'.hbcase'
if (hbnum .eq. 0) then
  write(*,*) 'Output file        : ', trim(filegofr)
else
  write(*,*) 'Output file        : ', trim(filegofr)//'.'//trim(surfix)
endif
write(*,fmt='(1X, "Bin number         : ", I5)') binNum
!
!
!Allocate  Files
allocate(PairNum(binNum))
allocate(tot_dhb(binNum))
!
!Initialization
PairNum=0
tot_dhb=0
!Recall that delta maximu can be 2 Angestrom
box=1.2
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
    open(unit=1, file=(trim(filename)//'.sigma'), status='old')
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
    read(1,*)
    read(2,*)
    read(2,*)
    read(2,*)
  enddo
endif
!
do 
  !******************************************************
  !Read .sigma
  read(1,*,iostat=ierror) step, time, dummy, dummy, d, dummy
  !if (ncount .gt. stepstop .and. stepstop .ne. -1) exit
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
  if (cs .eq. 2 ) then
    call Continue_Judgement(hbnum,readcase,ctn)
    if (ctn .eq. 0) cycle
  endif
  !******************************************************
  !------------------Count PAIRS-----------------------
  nbin=nint(DBLE(d)/box*binNum)
  if (nbin .gt. binNum) nbin = binNum
  PairNum(nbin)=PairNum(nbin)+1
  !======
  !These lines below is copied from hydrogen_bond.f90 for reminder.
  !if (num_accept .eq. 3 .and. num_donate .eq. 0) hbcase=1
  !if (num_accept .eq. 3 .and. num_donate .eq. 1) hbcase=2
  !if (num_accept .eq. 4 .and. num_donate .eq. 0) hbcase=3
  !if (num_accept .eq. 4 .and. num_donate .eq. 1) hbcase=4
  !Reminder ends
  !======
  if (readcase .eq. 2 .or. readcase .eq. 4) tot_dhb(nbin)=tot_dhb(nbin)+1
  !
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
do n=1,binNum-1
  !
  !construct some details of normalization and distance  
  r = DBLE(n)/binNum*box
  !
  !
  gofr = DBLE(PairNum(n))/(DBLE(ncount))*box
  tot=tot+gofr
  !
  if (n .le. 22 .and. n .ge. 13) c=c+tot_dhb(n)
  if (PairNum(n) .eq. 0) then
    write(3,*) r, gofr
  else
    write(3,*) r, gofr, DBLE(tot_dhb(n))/DBLE(PairNum(n))
  endif
  !
enddo
!
write(*,*) " The integrated result of curve is :", tot
write(*,*) " The sum of DHB from 13 to 22 is:", c
!
write(*,*) ' Program completes!'
!
end subroutine final
! 
!
!
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
if (hb .eq. 0) continu=1
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
end PROGRAM delta
