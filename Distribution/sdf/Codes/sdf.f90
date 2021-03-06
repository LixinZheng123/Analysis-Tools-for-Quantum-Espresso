!Spatial Distribution Function g(r) for Quantum Espresso Output Files
!Calculates SDF around a certain O*.
!Created some time 2013
!Modified 20141103, refactored with ion_wfc_rdf.f90
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
   integer              :: nbin(3)           !bin index
   integer              :: binNum            !number of bins for g(r)
   integer              :: atom(2)
   !real(DP), allocatable:: grt(:)           !unnormalized g(r)
   !real(DP), allocatable:: temgrt(:)        !temp unnormalized g(r)
   real(DP)             :: r                 !total distance (not an array, print-out value)
   !real(DP)             :: gofr              !FINAL g(r) print-out variable (not an array, print-out values)  
   real(DP)             :: box               !length of the region (includes multiple supercells)
   !real(DP)             :: res              !resolution (size) of bins  binNum/Hbox (number/length)
   real(DP)             :: omega             !volume of cell
   !real(DP)             :: dens              !raw particle density
   real(DP)             :: vshell            !volume of the infinitesimal radial shell
   real(DP)             :: tot               !The add-up value of g(r)s of each r
   !======Variables concerning SDF======
   real(DP)             :: res(3)            !resolution (size) of bins  binNum/Hbox (number/length)
   real(DP)             :: dens(3)           !raw particle density
   real(DP), allocatable:: grt(:,:)          !unnormalized g(r)
   real(DP), allocatable:: temgrt(:,:)       !temp unnormalized g(r)
   real(DP)             :: gofr(3)           !FINAL g(r) print-out variable (not an array, print-out values)  
   !======Variables concerning .fss======
   real(DP)             :: rO(3,15),rH(3,30) !atomic positions (dimension,atomic-species)
   real(DP)             :: rIO(3)            !Ion position
   real(DP)             :: rIH(3,3)          !ion positions (dimension,numH), (dimension,numH*)
   integer              :: numION            !Number of Ions in each timestep. If numION=0, iiO=iiOp
   integer              :: numH              !Number of Hydrogen in each ion. cs=2, numH=1; cs=3, numH=3
   integer              :: iO,iH             !general index 
   integer              :: iiO               !index of O* read from Ion
   integer              :: iiH(3)        
   !======Variables concerning .hbcase======
   integer              :: readcase, hbnum
   integer              :: HBIndex(4,2)      !(num_HB, accept(1) or donate(2))
   integer              :: num_accept, num_donate
   !======Local variables======
   integer              :: ncount0=0         !number of total steps index
   integer              :: i, j, k, ii       !general indexes
   integer              :: n,m,l,p,q         !general indexes
   integer              :: stepp 
   integer              :: ctn
   real(DP)             :: d                 
   real(DP)             :: theta,phy
   real(DP)             :: sO(3,15),sH(3,30) !Spherical atomic positions
   real(DP)             :: sIO(3),sIH(3,3)   !Spherical atomic positions
   real(DP)             :: r_AHB(3,4)        !Accepted hydrogen bond by O*
   real(DP)             :: s_AHB(3,4)        !Spherical accepted hydrogen bond by O*
   real(DP)             :: rad2deg           !A function that turns angle unit from radian to degree
   !
   !
   !
   !namelist /input/ filename,filegofr, atom, cs, binNum, hbnum, stepstart, stepstop
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
!atom(1)=1
!atom(2)=1
do i=1,2
  if (atom(i) .eq. 1) ch_atom(i)='O'
  if (atom(i) .eq. 2) ch_atom(i)='H'
enddo
!
!Print Intro
write(*,*) '------------------------------'
write(*,fmt='("|      Spatial Distribution Function      |")') 
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
allocate(grt(3,binNum), temgrt(3,binNum))
!
!Initialization
grt=0
ii=0
!Recall that r_FSS=3.5 Angestrom (6.62 Bohr)
box=7
omega=box*box*box
res(1)=binNum/box
res(2)=binNum/pi
res(3)=binNum/(2*pi)
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
    open(unit=5, file=(trim(filename)//'.angle'), status='unknown')
    !
  case(2)
    close(1)
    close(2)
    close(3)
    close(5)
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
!if (stepstart .gt. 1) then
!  do k=1, stepstart-1
!    read(1,*) step, time, numION, nsp(1), nsp(2)
!    do i=1,nsp(1)+nsp(2)+1
!      read(1,*)
!    enddo
!    read(2,*)
!  enddo
!endif
!
do 
  !******************************************************
  !if (ncount0 .gt. 10000) exit
  !Read .fss
  read(1,*,iostat=ierror) step, time, numION, nsp(1), nsp(2)
  write(5,*) step, time
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
  !
  if (stepp .ne. step) then
    write(*,*) "Error! steps not equal in 2 files (.fss and .hbcase)!", step,stepp,time
    exit
  endif
  !
  !
  ncount0=ncount0+1
  !
  !******************************************************
  if (readcase .eq. 1) then
    num_accept=3
    num_donate=0
  else if (readcase .eq. 2) then
    num_accept=3
    num_donate=1
  else if (readcase .eq. 3) then
    num_accept=4
    num_donate=0
  else if (readcase .eq. 4) then
    num_accept=4
    num_donate=1
  endif
  !******************************************************
  !Read .hbcase
  read(2,*) HBIndex(1,1), HBIndex(2,1), HBIndex(3,1), HBIndex(4,1)
  read(2,*) HBIndex(1,2)
  !******************************************************
  !Read .fss
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
  q=0
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
    if (i .eq. HBIndex(1,1) .or. i .eq. HBIndex(2,1) .or. i .eq. HBIndex(3,1) .or. i .eq. HBIndex(4,1)) then
      q=q+1
      if (q .gt.4) write(*,*) "Error! AHB greater than 4!"
      do j=1,3
        !The accepted hydrogen bond by O*
        r_AHB(j,q)=rH(j,i)
      enddo
    endif
  enddo
  !
  !if (mod(ncount0,100) .ne. 0) then
  !  cycle
  !else
  !  write(*,*) ncount0
  !endif
  !******************************************************
  !CORE DIFFERENCE of SDF from rdf
  if (cs .eq. 2 ) then
    !
    !
    !=========
    !In the new coordinate system, the atoms should have coordinate like:
    !-------------------------------------------------------------------------------
    !                    |xyz                                    |sph
    !O*                  |(0,0,0)	                         |(0,0,0)
    !H*                  |(0,0,d(O*-H*))                         |(d(O*-H*),0,0)
    !1st HB donated to O*|(d(O*-H)cos(alpha),0,d(O*-H)sin(alpha))|(d(O*-H),0,alpha)
    !-------------------------------------------------------------------------------
    !(alpha is a parameter of no importance.)
    !=========
    !
    !
    !==============================
    !First ranslational transformation, so that rIO has the xyz coordinate of (0,0,0)
    !==============================
    !
    !
    !
    !
    do i=1,nsp(1)
      call TRANS(rO(1:3,i),rIO(1:3))
    enddo
    !
    call TRANS(rIH(1:3,1),rIO(1:3))
    call xyz2sph(rIH(1:3,1),sIH(1:3,1))
    !
    do i=1,num_accept
      call TRANS(r_AHB(1:3,i),rIO(1:3))
    enddo
    call xyz2sph(r_AHB(1:3,1),s_AHB(1:3,1))
    !===
    !
    !write(*,*) "========== After translation ==========="
    !write(*,fmt='("H* sph coordinate:      ",3F9.3)') sIH(1,1),sIH(2,1)/pi*180,sIH(3,1)/pi*180
    !write(*,fmt='("1st AHB* sph coordinate:",3F9.3)') s_AHB(1,1),s_AHB(2,1)/pi*180,s_AHB(3,1)/pi*180
    !write(*,*)
    !
    !
    !
    !
    !==============================
    !Then rotational transformation
    !==============================
    !
    !
    !******************************
    !1st round
    theta=sIH(2,1)
    phy=sIH(3,1)
    write(5,*) theta, phy
    !
    do i=1,nsp(1)
      call ROTATE(rO(1:3,i),theta,phy,sO(1:3,i))
    enddo
    !
    call ROTATE(rIH(1:3,1),theta,phy,sIH(1:3,1))
    !
    do i=1,num_accept
      call ROTATE(r_AHB(1:3,i),theta,phy,s_AHB(1:3,i))
    enddo
    !===
    !write(*,fmt='("Rotate theta ",F6.2)'), theta/pi*180
    !write(*,fmt='("Rotate phy   ",F6.2)'), phy/pi*180
    !write(*,*) "=========== After 1st rotation ============="
    !write(*,fmt='("H* sph coordinate:      ",3F9.3)') sIH(1,1),sIH(2,1)/pi*180,sIH(3,1)/pi*180
    !write(*,fmt='("1st AHB sph coordinate: ",3F9.3)') s_AHB(1,1),s_AHB(2,1)/pi*180,s_AHB(3,1)/pi*180
    !write(*,*)
    !
    !******************************
    !2nd round
    theta=0
    phy=s_AHB(3,1)
    write(5,*) s_AHB(3,1:num_accept)
    !
    do i=1,nsp(1)
      call ROTATE(rO(1:3,i),theta,phy,sO(1:3,i))
    enddo
    !
    call ROTATE(rIH(1:3,1),theta,phy,sIH(1:3,1))
    !
    do i=1,num_accept
      call ROTATE(r_AHB(1:3,i),theta,phy,s_AHB(1:3,i))
    enddo
    !===
    !write(*,fmt='("Rotate theta ",F6.2)'), theta/pi*180
    !write(*,fmt='("Rotate phy   ",F6.2)'), phy/pi*180
    !write(*,*) "=========== After 2nd rotation ============="
    !write(*,fmt='("H* sph coordinate:      ",3F9.3)') sIH(1,1),sIH(2,1)/pi*180,sIH(3,1)/pi*180
    !write(*,fmt='("1st AHB sph coordinate: ",3F9.3)') s_AHB(1,1),s_AHB(2,1)/pi*180,s_AHB(3,1)/pi*180
    !write(*,*)
    !write(*,fmt='("Other AHB phy:          ",3F9.3)') s_AHB(3,2:num_accept)/pi*180
    !===
    !
  endif
  !
  !******************************************************
  if (cs .eq. 2 ) then
    call Continue_Judgement(hbnum,readcase,ctn)
    if (ctn .eq. 0) cycle
  endif
  !******************************************************
  !------------------Count PAIRS-----------------------
  do i=2,hbnum
    !do j=1,3
      j=3
      nbin(j) = nint(s_AHB(j,i)*res(j))
      !if (nbin(j) .lt. 0) nbin(j)=nbin(j)+binNum
      if (nbin(j) .gt. binNum) nbin = binNum
      !write(*,*) i,nbin(j)
      grt(j,nbin(j))=grt(j,nbin(j))+1
    !enddo
  enddo
  !do n=1,binNum,1
    !do j=1,3
  !    j=3
  !    grt(j,n)=grt(j,n)+temgrt(j,n)
    !enddo
  !enddo
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
  phy=DBLE(n)/res(3)
  !  vshell = volume of the infinitesimal shell
  !vshell=4.0*pi*(r**2)*(1/res(1)) 
  !dens(1)=DBLE(hbnum)/omega
  !dens(2)=DBLE(hbnum)/pi
  dens(3)=DBLE(hbnum)/(2*pi)
  !
  !
  !Calculate the final, normalized, value of g(r)! Please note that the
  !normalization constant (ncount*mic3*nsp(atom1)*dens*vshell)...
  !  ncount*mic3 = number of steps and number of extended shells
  !  nsp(atom1)*dens = number of pairs
  !if (atom(1) .eq. 1) gofr = grt(1,n)/(DBLE(ncount*vshell))
  !if (atom(1) .eq. 2) gofr = grt(1,n)/(DBLE(ncount*vshell*numH))
  !gofr(1) = grt(1,n)/(DBLE(ncount*vshell*dens(1)))

  !gofr(2)=grt(2,n)/(DBLE(ncount)*dens(2))
  gofr(3)=grt(3,n)/(DBLE(ncount)*dens(3))
  tot=tot+gofr(3)
  !
  !write(3,*) (r*convertBA), gofr(1), gofr(2), gofr(3)
  write(3,*) , phy/pi*180, gofr(3)
  !
enddo
write(*,*) "The integrated result of curve is :", tot*(1/res(3))
!
write(*,*) ' Program complete!'
!
end subroutine final
! 
!
!
subroutine TRANS(vec1,vec2)

!This is to rotate the rx, ry, rz, so that the new z axis is towa

implicit none

real(8), intent(inout)    :: vec1(3)
real(8), intent(in)       :: vec2(3)
real(8), parameter        :: pcell=23.5170
integer                   :: j

do j=1,3
  vec1(j)=vec1(j)-vec2(j)
  vec1(j)=vec1(j)-nint(vec1(j)/pcell)*pcell
enddo

return

end subroutine TRANS
!******************************************************************************************
!******************************************************************************************
!
subroutine ROTATE(vec1,rtheta,rphy,vec2)

!This is to rotate the rx, ry, rz, so that the new z axis is towa

implicit none

real(8), intent(inout)    :: vec1(3), vec2(3)
real(8), intent(inout)    :: rtheta, rphy
real(8)                   :: temp(3)

!
!Rotation matrix of passive transformationi is R^(-1). 
!Rotates counter-clockwise about the original xy coordinate.
!(see http://en.wikipedia.org/wiki/Rotation_matrix for details)
!   
!R= | cos(gamma)  -sin(gamma) |
!   | sin(gamma)   cos(gamma) |
!
!First rotate around the old z axis, counter-clockwise, of angle phy;
!           | cos(phy) -sin(phy)     0    |
! R(phy)=   | sin(phy)  cos(phy)     0    |
!           |    0         0         1    |
!
! then rotate around the new y axis, counter-clockwise, of angle theta.
!
!           |  cos(theta)  0    sin(theta)  |
! R(theta)= |     0        1        0       |
!           | -sin(theta)  0    cos(theta)  |
!
!new_vecter=(R(theta))^(-1)*(R(phy))^(-1)
!
!            |  cos(theta)  0   -sin(theta)  |   | cos(phy)  sin(phy)   0 |
!          = |     0        1        0       | * |-sin(phy)  cos(phy)   0 |
!            |  sin(theta)  0    cos(theta)  |   |    0         0       1 |
!
!           |  cos(theta)cos(phy)   cos(theta)sin(phy) -sin(theta)  |
!         = |     -sin(phy)              cos(phy)           0       |
!           |  sin(theta)cos(phy)   sin(theta)sin(phy)   cos(theta) |
!

temp(1) = + vec1(1)*cos(rtheta)*cos(rphy)  &
          + vec1(2)*cos(rtheta)*sin(rphy)  &
          - vec1(3)*sin(rtheta)
temp(2) = - vec1(1)*sin(rphy)              &
          + vec1(2)*cos(rphy)              &
          + vec1(3)*0
temp(3) = + vec1(1)*sin(rtheta)*cos(rphy) &
          + vec1(2)*sin(rtheta)*cos(rphy) &
          + vec1(3)*cos(rtheta)


do j=1,3
  vec1(j)=temp(j)
enddo

!r=sqrt(x**2+y**2+z**2)
vec2(1)=sqrt((vec1(1))**2+(vec1(2))**2+(vec1(3))**2)
!theta=acos(z/r)
vec2(2)=acos(vec1(3)/vec2(1))
!phy=acos(x/sqrt(x**2+y**2))
!vec2(3)=atan(vec1(2)/vec1(1))

if (vec1(1) .gt. 0 .and. vec1(2) .lt. 0) then
  vec2(3)=atan(vec1(2)/vec1(1))+2*pi
!Quadrant II & III
else if (vec1(1) .lt. 0) then
  vec2(3)=atan(vec1(2)/vec1(1))+pi
!Quadrant I
else
  vec2(3)=atan(vec1(2)/vec1(1))
endif

return

end subroutine ROTATE
!******************************************************************************************
!******************************************************************************************
subroutine xyz2sph(vec1,vec2)

implicit none

real(8), intent(in)      :: vec1(3)
real(8), intent(out)     :: vec2(3)

!r=sqrt(x**2+y**2+z**2)
vec2(1)=sqrt((vec1(1))**2+(vec1(2))**2+(vec1(3))**2)
!theta=acos(z/r)
vec2(2)=acos(vec1(3)/vec2(1))
!phy=acos(x/sqrt(x**2+y**2))
!phy=atan(y/x)
!Quadrant IV
if (vec1(1) .gt. 0 .and. vec1(2) .lt. 0) then
  vec2(3)=atan(vec1(2)/vec1(1))+2*pi
!Quadrant II & III
else if (vec1(1) .lt. 0) then
  vec2(3)=atan(vec1(2)/vec1(1))+pi
!Quadrant I
else
  vec2(3)=atan(vec1(2)/vec1(1))
endif

return

end subroutine xyz2sph
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
end PROGRAM rdf 

!function rad2deg(rad) result(deg)
!  real(8), parameter  :: pi=4.d0*atan(1.d0)
!  real(8), intent(in) :: rad
!  real(8)             :: deg
!  deg=rad/pi*180
!end function rad2deg

