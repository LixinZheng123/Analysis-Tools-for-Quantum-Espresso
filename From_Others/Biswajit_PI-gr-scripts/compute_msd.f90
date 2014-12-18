PROGRAM calcmsd
!
IMPLICIT NONE
!
! Global Parameters (consistent with Quantum Espresso)
!
INTEGER, PARAMETER  :: DP = selected_real_kind(14,200) !double-precision kind
REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP !pi to double-precision
REAL(DP), PARAMETER :: bohr2Ang = 0.52917720859_DP    !bohr (au) to Angstrom conversion factor
REAL(DP), PARAMETER :: amu_au = 1822.888484265_DP     !bohr (au) to Angstrom conversion factor
REAL(DP), PARAMETER :: tps = 0.00048_DP               !dt in picoseconds, dt=0.48 femtoseconds * 1 step
!
! Local Parameters
!
INTEGER, PARAMETER :: nsp = 2                         !number of atom species
!
! Local variables dependent on local parameters
!
INTEGER  :: nat, na(nsp)                              !total number of atoms, number of atoms in each species
REAL(DP) :: amass(nsp),pmass(nsp),dis(nsp)            !mass amu, mass in a.u., displacement for each species 
!
! Local Variables 
!
INTEGER :: ia,ic,is,isa,i,nconf
REAL(DP) :: rdist(3),cdm(3),mean_dis,r2,time,tmass
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: tau, taui
CHARACTER(len=100) :: inputfile,atom 
!
! Provide values for number of atoms in each species ... 
!
na(1)=64                ! number of Oxygen atoms 
na(2)=128               ! number of Hydrogen/Deuterium atoms 
!
! Provide values for mass in a.m.u for each species ...
!
amass(1)=15.9994_DP      ! mass of Oxygen (amu) 
amass(2)=1.00794_DP      ! mass of H (amu)
!
! Computing total number of atoms
!
nat=SUM(na)
!
! READ from standard input
!
READ(*,*)inputfile
READ(*,*)nconf
READ(*,*)na(1)
READ(*,*)na(2)
!
nat=SUM(na)
!
! Open files used during program execution...
!
OPEN (UNIT=10, FILE=inputfile,ACTION='READ',FORM='FORMATTED')
OPEN (UNIT=20, FILE='msd.dat',ACTION='WRITE',FORM='FORMATTED')
!
! compute total mass in a.u.
!
tmass=0.0_DP
!
DO is=1,nsp
   pmass(is)=amass(is)*amu_au
   tmass=tmass+na(is)*pmass(is)
END DO
!
ALLOCATE(taui(3,nat)); taui=0.0_DP
ALLOCATE(tau(3,nat)); tau=0.0_DP
!
! Read coordinates of the first timestep and store it as initial coordinate .. in centre of mass (COM) reference frame ..
!
DO ic=1,1
  !
  READ(10,*) 
  READ(10,*) 
  !
  ! Read in atomic coordinates (in bohr) from file...
  !
  DO ia=1,nat
    READ(10,*) atom,taui(1,ia),taui(2,ia),taui(3,ia)
  END DO
  !
  ! Compute COM ... 
  !
  DO i=1,3
    !
    cdm(i)=0.0_DP
    isa=0
    DO is=1,nsp
      DO ia=1,na(is)
        isa = isa + 1
        cdm(i)=cdm(i)+taui(i,isa)*pmass(is)
      END DO
    END DO
    cdm(i)=cdm(i)/tmass
    !
  END DO
  !
  ! Transform intial coordinate in COM ref. frame .. 
  !
  DO isa=1,nat
    taui(:,isa)=taui(:,isa)-cdm(:)
  END DO
  !
END DO !ic
!
! Start loop over MD configurations...
!
DO ic=2,nconf
  !
  READ(10,*) 
  READ(10,*) 
  !
  ! Read in atomic coordinates (in bohr) from file...
  !
  DO ia=1,nat
    READ(10,*) atom,tau(1,ia),tau(2,ia),tau(3,ia)
  END DO
  !
  ! Compute COM ... 
  !
  DO i=1,3
    !
    cdm(i)=0.0_DP
    isa = 0
    DO is=1,nsp
      DO ia=1,na(is)
        isa = isa + 1
        cdm(i)=cdm(i)+tau(i,isa)*pmass(is)
      END DO
    END DO
    cdm(i)=cdm(i)/tmass
    !
  END DO
  !
  ! Compute atomic displacement in COM ref frame from initial coordiante in COM ref frame
  !
  isa=0
  !
  DO is=1, nsp
    !
    dis(is)=0.0_DP
    r2=0.0_DP
    !
    DO ia=1, na(is)
      !
      isa=isa+1
      rdist(:)=tau(:,isa)-cdm(:)
      r2=r2+SUM( ( rdist(:) - taui(:,isa) )**2 )
      !
    END DO
    !
    dis(is)=dis(is) + r2 / DBLE(na(is))
    !
  END DO
  !
  ! convert displacement to Angstrom^2
  !
  dis(:)=dis(:)*bohr2Ang*bohr2Ang
  !
  ! compute msd averaged over all atoms
  !
  mean_dis=SUM(dis(:)*na(:))/DBLE(nat)
  !
  ! compute time
  !
  time=DBLE((ic-1))*tps
  !
  WRITE(20,'(4F13.6)')time,dis(:),mean_dis  
  !
END DO !ic
!
! Some housekeeping...
!
DEALLOCATE(taui)
DEALLOCATE(tau)
!
! Close files used during program execution...
!
CLOSE(10)
CLOSE(20)
!
!
END PROGRAM calcmsd
