PROGRAM calcgOOr
!
IMPLICIT NONE
!
! Global Parameters (consistent with Quantum Espresso)
!
INTEGER, PARAMETER :: DP = selected_real_kind(14,200) !double-precision kind
REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP !pi to double-precision
REAL(DP), PARAMETER :: bohr2Ang = 0.52917720859_DP    !bohr (au) to Angstrom conversion factor
!
! Local Parameters
!
INTEGER             :: nat = 192                       !number of atoms in simulation cell
INTEGER             :: natO = 64                       !number of oxygen atoms in simulation cell
INTEGER             :: natH = 128                      !number of hydrogen atoms in simulation cell
INTEGER, PARAMETER  :: nbins = 622                     !number of histogram bins  => deltar ~ 0.01 Angstrom
REAL(DP), PARAMETER :: alat = 23.517_DP*bohr2Ang      !lattice constant for water64 (in angstrom)
REAL(DP), PARAMETER :: omega = alat**(3.0_DP)         !simulation cell volume (in Angstrom^3)
REAL(DP)            :: eta                            !number density=number of oxygen atoms/cell volume
REAL(DP), PARAMETER :: Vpi = (4.0_DP/3.0_DP)*pi       !(4/3)pi for normalization
!
! Local Variables
!
INTEGER :: ia,ib,ic,ig,nstep,nconf,nskip
REAL(DP) :: deltar,dAB,r,deltaV,time,h(3,3),ainv(3,3)
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: atxyz,gOOr
CHARACTER(len=100) :: inputfile,atom
!
! READ from standard input
READ(*,*)inputfile
READ(*,*)nconf
READ(*,*)nskip
READ(*,*)natO
READ(*,*)natH
nat=natO+natH
eta = natO/omega
!
! Open files used during program execution...
!
OPEN (UNIT=10, FILE=inputfile,ACTION='READ',FORM='FORMATTED')
OPEN (UNIT=20, FILE='gOO.dat',ACTION='WRITE',FORM='FORMATTED')
!
!Initialize cell Parameters ...
h=0.0_DP
h(1,1)=alat
h(2,2)=alat
h(3,3)=alat
ainv=0.0_DP
ainv(1,1)=1.0_DP/alat
ainv(2,2)=1.0_DP/alat
ainv(3,3)=1.0_DP/alat
!
! Initialize atxyz and gOO(r)...
! N.B.: gOOr indices are (nbins,1)=r and (nbins,2)=g(r)...
!
ALLOCATE(atxyz(3,nat)); atxyz=0.0_DP
ALLOCATE(gOOr(0:nbins,1:2)); gOOr=0.0_DP
deltar=alat/(2.0_DP*nbins) !bin/histogram width (in Angstrom)
  DO ig=0,nbins-1
    ! Compute center of each histogram bin (r)...
    r=deltar*(DBLE(ig)+0.5_DP)    
    gOOr(ig,1)=r
  END DO
!
! Skip initial configurations...
!
DO ic=1,nskip
  !
  READ(10,*)
  READ(10,*) 
  !
  ! Read in atomic coordinates (in bohr) from file...
  !
  DO ia=1,nat
    READ(10,*) atom,atxyz(1,ia),atxyz(2,ia),atxyz(3,ia)
  END DO
  !
END DO
!
! Start computimng g(r) from MD configurations after skipping initial
! configurations ...
!
DO ic=nskip+1,nconf
  !
  READ(10,*)
  READ(10,*) 
  !
  ! Read in atomic coordinates (in bohr) from file...
  !
  DO ia=1,nat
    READ(10,*) atom,atxyz(1,ia),atxyz(2,ia),atxyz(3,ia)
  END DO
  !
  ! Convert atomic coordinates from bohr into Angstroms...
  !
  DO ia=1,nat
    atxyz(1,ia)=atxyz(1,ia)*bohr2Ang
    atxyz(2,ia)=atxyz(2,ia)*bohr2Ang
    atxyz(3,ia)=atxyz(3,ia)*bohr2Ang
  END DO
  !
  ! Populate histogram bins for gOO(r)...
  !
  DO ia=1,natO
    DO ib=1,natO
      IF (ib.NE.ia) THEN
        ! Compute the distance of the closest image of atom B to atom A using minimum image convention...
        CALL CalcMinDist(h,ainv,atxyz(1,ia),atxyz(2,ia),atxyz(3,ia),atxyz(1,ib),atxyz(2,ib),atxyz(3,ib),dAB)
        ! Screen distances that are greater than alat/2 (half the box size)...
        IF (dAB.LT.alat/2.0_DP) THEN
          ig=INT(dAB/deltar)  !bin/histogram position
          gOOr(ig,2)=gOOr(ig,2)+1
        END IF
      END IF
    END DO !ib
  END DO !ia
  !
END DO !ic
!
! Final normalization of gOO(r) and printing...
!
  DO ig=0,nbins-1
    deltaV=Vpi*((ig+1)**(3.0_DP)-ig**(3.0_DP))*deltar**(3.0_DP) !bin/shell volumes...
    gOOr(ig,2)=gOOr(ig,2)/((nconf-nskip)*natO*deltaV*eta) !includes 1/nat and 1/nconf normalization factors...
    WRITE(20,'(2F13.6)') gOOr(ig,1),gOOr(ig,2)
  END DO
!
! Some housekeeping...
!
DEALLOCATE(atxyz)
DEALLOCATE(gOOr)
!
! Close files used during program execution...
!
CLOSE(10)
CLOSE(20)
!
!
CONTAINS
  !
  !
  SUBROUTINE CalcMinDist(h,ainv,xA,yA,zA,xB,yB,zB,dAB)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200) !double-precision kind
  REAL(DP) :: alat,xA,yA,zA,xB,yB,zB,dAB
  REAL(DP) :: rAB(3),rAB2(3),ainv(3,3),h(3,3)
  !Initialization of distance 
  dAB=0.0_DP
  !
  ! Compute distance between atom A and atom B (according to the minimum
  ! image convention)...
  !
  rAB(1)=xA-xB   ! r_AB = r_A - r_B   
  rAB(2)=yA-yB   ! r_AB = r_A - r_B   
  rAB(3)=zA-zB   ! r_AB = r_A - r_B   
  !
  rAB2(1)=ainv(1,1)*rAB(1)+ainv(1,2)*rAB(2)+ainv(1,3)*rAB(3)   ! s_AB =h^-1 r_AB
  rAB2(2)=ainv(2,1)*rAB(1)+ainv(2,2)*rAB(2)+ainv(2,3)*rAB(3)   ! s_AB =h^-1 r_AB
  rAB2(3)=ainv(3,1)*rAB(1)+ainv(3,2)*rAB(2)+ainv(3,3)*rAB(3)   ! s_AB =h^-1 r_AB
  !
  rAB2(1)=rAB2(1)-IDNINT(rAB2(1))   ! impose MIC on s_AB in range:[-0.5,+0.5]
  rAB2(2)=rAB2(2)-IDNINT(rAB2(2))   ! impose MIC on s_AB in range:[-0.5,+0.5]
  rAB2(3)=rAB2(3)-IDNINT(rAB2(3))   ! impose MIC on s_AB in range:[-0.5,+0.5]
  !
  rAB(1)=h(1,1)*rAB2(1)+h(1,2)*rAB2(2)+h(1,3)*rAB2(3)   ! r_AB = h s_AB(MIC)
  rAB(2)=h(2,1)*rAB2(1)+h(2,2)*rAB2(2)+h(2,3)*rAB2(3)   ! r_AB = h s_AB(MIC)
  rAB(3)=h(3,1)*rAB2(1)+h(3,2)*rAB2(2)+h(3,3)*rAB2(3)   ! r_AB = h s_AB(MIC)
  !
  dAB=DSQRT(rAB(1)*rAB(1)+rAB(2)*rAB(2)+rAB(3)*rAB(3))   ! |r_A -r_B| (MIC)
  !
  RETURN
  END SUBROUTINE CalcMinDist
  !
  !
END PROGRAM calcgOOr
