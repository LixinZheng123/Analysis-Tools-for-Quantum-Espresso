!==================================================================================================
!Radial Distribution Function g(r) for Quantum Espresso Output Files
!
!Author: Charles W. Swartz VI
!
!
!
!
!
!
!==================================================================================================
 PROGRAM new_rdf
   !
   USE rdf_lib
   USE omp_lib
   !
   implicit none
   !
   integer, parameter   :: namax = 20                    !max number of species   
   real, parameter      :: pi    = 4.d0*atan(1.d0)
   !
   integer              :: binNum,     &  !number of bins for g(r)
                           ntsp,       &  !number of total atomic species, limit namax
                           nsp(namax), &  !number of each atomic species     
                           nat=0,      &  !total number of atoms
                           mic,        &  !Multiple Supercells 
                                             !  mic = 1, Only the main supercell
                                             !  mic = 2, translated in the pos direction
                                             !  mic = 3, translated in both the pos and neg directions
                                             !  Also, scaled length of region
                           mic3,          &  !number of extended cells, mic**3 (line )
                           stepstart,     &  !nfi that we will start reading *.pos 
                           stepstop=-1,   &  !nfi that we will stop reading *.pos AFTER (optional)
                                             ! NOTE: a negative value will bypass and read until EOF
                           interval,      &  !interval between nfi in *pos file
                           nsteps,        &  !total number of steps
                           ncount=0,      &  !number of steps index
                           i, j, k,       &  !general indexes
                           offsetI,       &  !initial offset
                           offsetF,       &  !final offset
                           atom1,         &  !first atom type
                           atom2,         &  !first atom type
                           n,             &  !general index (extended supercell)
                           m,             &  !secondary general index
                           ns,            &  !number of species index
                           na,            &  !number of atoms index
                           readstep,      &  !read the nfi from *.pos 
                           ierror,        &  !error index
                           iADD,          &  !number of to add to pair correlation
                           nbin,          &  !bin index
                           pcount_end,    &  !first limit used in the pair counting process
                           pcount_start      !second limit used in the pair counting process
                           
   !
   real(DP), allocatable     ::   ngdr(:), &  !main pair correlation (unnormalized g(r)) array
                                 ngdrt(:)    !Temp pair correlation(unnormalized g(r)) array

   real(DP), allocatable    :: tau(:,:,:),    &  !atomic positions (dim, index, atomic-species)
                               stau(:,:,:),   &  !scaled atomic positions (dim,index, atomic-species)
                               staux(:,:,:)      !multiple scaled atomic positions (dim, index, atomic-species)
   !
   real(DP)             :: gofr,          &  !FINAL g(r) print-out variable (not an array, print-out values)  
                           r,             &  !total distance (not an array, print-out value)
                           aprim(3,3),    &  !lattice prim vectors
                           apinv(3,3),    &  !inverse of prim vectors, for scaled positions
                           box,           &  !length of the region (includes multiple supercells)
                           Hbox,          &  !length of half the region (including multiple supercells)
                                             ! Recall: we can only compute up to
                                             ! L/2 of the total region, which is
                                             ! why we need multiple super cells
                           res,           &  !resolution (size) of bins  binNum/Hbox (number/length)
                           omega,         &  !volume of cell
                           omegaT=0.0,    &  !total volume of the cell
                           dens,          &  !raw particle density
                           dx,            &  !diff in x-value
                           dy,            &  !diff in y-value
                           dz,            &  !diff in z-value
                           sdist(3),      &  !square of the components distance in scaled coordinates
                           rdist(3),      &  !square of the components distance in real coordinates
                           r2,            &  !total squared distance
                           norm=1.0,      &  !normalizing factor (used to normalize final g(r), if there are multiple atoms  
                           time,          &  !timestep
                           dmic,          &  !convert mic to real
                           vshell,        &  !volume of the infinitesimal radial shell
                           dummy,         &  !a dummy variable
                           dum(3)            !a dummy array (remove!)

   character(len=30)       :: filePos, fileCel, fileOut
   integer                 :: start_time, end_time, total_time
   !
   namelist /input/ filePos, fileCel, fileOut, binNum, ntsp, nsp, mic, &
                  stepstart, stepstop, atom1, atom2, norm, lattice
   !
   !
   !Start Time
   !$ start_time = OMP_get_wtime()   
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
   !$ end_time = OMP_get_wtime()
   !$ total_time  = end_time - start_time
   !$ write(*,'(A, I5, A)') 'Total Computing Time: ', total_time, ' s'
   !    
   !
   !********************************************************************************
   !********************************************************************************
   contains
   !********************************************************************************
   !********************************************************************************
      !
      !
      !
      !******************************************************
      !------------------------------------------------------
      !Initialization Subroutine
      !sets mic, and allocates memory
      !------------------------------------------------------
      !******************************************************
      subroutine init
         !
         implicit none
         !
         integer     :: ll
         !
         !read Namelist input for stdin
         read(*,input)
         !
         !Total number of atoms
         do ns=1,ntsp
            nat = nat + nsp(ns)
         enddo
         !
         !Print Intro
                  write(*,*) '---------------------------------------------------------'
                  write(*,*) '|          Radial Distribution Function g(r)            |'
                  write(*,*) '|               Quantum Espresso Output                 |'
                  write(*,*) '|                                                       |'
                  write(*,*) '|                 Charles W Swartz VI                   |'
                  write(*,*) '---------------------------------------------------------'
                  write(*,*)
         !
         !Check for multiple Supercells
         select case (mic)
            case(1)
               offsetI = 0
               offsetF = 0
            case(2)
               offsetI = 0
               offsetF = 1
            case(3)
               offsetI = -1
               offsetF = 1
            case default
               write(*,*) ' ERROR: incorrect mic!! Only 0-3 supported!'
               stop
         end select
         mic3 = mic**3
         !
         !Confrim Values
         write(*,fmt='(1X, " Calculating the g(r) for atom1 :" , I2," and atom2 : ", I2)' )  atom1 , atom2
         write(*,*)
         write(*,*) ' Position file          : ', filePos
         write(*,*) ' Cell file              : ', fileCel
         write(*,*) ' Output file            : ', fileOut
         write(*,fmt='(1X, " StepStart              : ", I8)') stepstart
         write(*,fmt='(1X, " StepStop               : ", I8)') stepstop
         write(*,fmt='(1X, " Bin number             : ", I8)') binNum
         write(*,*)
         write(*,fmt='(1X," Total atomic Speices   : ",4X, I3)') ntsp
         do ll=1,ntsp,1
            write(*,fmt='(1X,"  Species num ", I3, " total : ",4X, I3)')  ll,  nsp(ll) 
         enddo 
         write(*,*)
         write(*,fmt='(1X, " mic : ", I2, " offsetI : ", I2, " offsetF : ", I2)')  mic, offsetI, offsetF 
         if (norm == 1) then
            write(*,*) ' Extra normalization set to unity' 
         else
            write(*,fmt='(1X, " Warning, Extra normalization set to  ", F3.1)') norm
         endif
         write(*,*)   
         !
         !Allocate  Files
         allocate( tau(3,(mic**3*nat),ntsp)    )
         allocate( stau(3,(mic**3*nat),ntsp)   )
         allocate( staux(3,(mic**3*nat),ntsp)  )
         allocate( ngdr(binNum), ngdrt(binNum) )
         !
      end subroutine init
      !
      !
      !
      !******************************************************
      !------------------------------------------------------
      !files subroutine
      !Opens Files
      !------------------------------------------------------         
      !******************************************************
      subroutine files(io)
         !
         implicit none
         !
         integer, intent(in)  ::io
         !
         select case(io)
            case(1)
               open(unit=1, file=(trim(filePos)), status='old')
               open(unit=2, file=(trim(fileCel)), status='old')
               open(unit=3, file=(trim(fileOut)), status='unknown')
               !
               write(*,*)
               write(*,*) ' Data to be read from position file ', trim(filePos),  &
                           ' and cell-file ', trim(fileCel)
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
      ! read-in data, calculate inverse, apply mic, r_to_s,
      ! 
      ! 
      !------------------------------------------------------
      !******************************************************
      subroutine main
         !
         implicit none
         !
         write(*,*)
         write(*,*) ' Begining main loop (This may take a few moments) ... '
         !
         Main_loop: do 
            !
            !-----read *.pos file-----
            read(1,*, iostat=ierror) readstep, time
            !
            !-----Check for end of *.pos file-----
            if (ierror < 0) then
               write(*,*) '  End of File Reached'
               write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
               exit
            endif
            !
            !-----Check for stop condition in *.pos file-------
            ! Note: stepstop is less then zero for a full EOF read
            if(.not. (stepstop < 0) .and. (readstep > stepstop ) ) then
               write(*,fmt='(1X,"  Final step reached : ",4X, I7)') stepstop
               write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(ncount)
               exit
            endif
            !
            !----Check for start condition------
            if (readstep < stepstart) then
               cycle
            endif
            !
            !-----read *cel file------
            read(2, *) dummy, dummy
            do i=1,3,1
               read(2, *) (aprim(i,j),j=1,3)
            enddo
            !
            !calculate the inverse
            call invert(aprim, apinv, omega)
            omegaT = omegaT + omega
            !
            !region dimensions for THIS loop
            box  = mic*omega**(1.0d0/3.0d0)
            Hbox = box/2.0d0
            res  = binNum/Hbox
            !
            !Read in current position values
            do ns=1,ntsp,1
               do na=1,nsp(ns),1
                  read(1,*) (tau(j,na,ns),j=1,3)   
               enddo
            enddo
            !
            !
            !Convert to scaled positions 
            !write(34,'(A, I8)') 'Cycle Number: ', ncount
            do ns=1,ntsp,1
               !write(34,'(A, I8)') 'Species number : ', ns
               do na=1,nsp(ns),1
                  call r_to_s( tau(1:3,na,ns), stau(1:3,na,ns), apinv )
               enddo
            enddo
            !
            !Zero all temp g(r)
            do n=1,binNum,1
               ngdrt(n) = 0
            enddo
            !
            !Shift all values to the origin (faster)
            do ns=1,ntsp
               do na=1,nsp(ns)
                  CALL PBCS(stau(1,na, ns), stau(2,na, ns), stau(3,na, ns), 1)
               enddo
            enddo 
            !
            !Extended Multiple supercell
            do ns=1,ntsp,1
               n=0
               do i=offsetI,offsetF
                  do j=offsetI,offsetF
                     do k=offsetI,offsetF
                        do na=1,nsp(ns),1
                        n = n +1
                        staux(1,n,ns) = stau(1,na,ns) + i
                        staux(2,n,ns) = stau(2,na,ns) + j
                        staux(3,n,ns) = stau(3,na,ns) + k
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            ! 
            !---------Count PAIRS-----------
            !Atoms to be used in the pair correlation
            k = atom1
            j = atom2
            !
            !If they are the same we can simplify the process, and double count
            ! no such luck if they are different
            if ( k == j) then
               pcount_end = nsp(k) * mic3 - 1
               IADD       = 2
            else
               pcount_end = nsp(k) * mic3
               IADD       = 1
            endif
            !
            !! $OMP PARALLEL PRIVATE(dx, dy, dz, sdist, rdist, r2, nbin)
            !! $OMP DO 
            do n=1,pcount_end,1
               !
               !if we are double counting set second limit 
               if ( k == j) then
                  pcount_start = n + 1
               else
                  pcount_start = 1
               endif
               !
               do m=pcount_start, nsp(j) * mic3
                  !
                  !coordinate distance
                  dx = staux(1,n,k) - staux(1,m,j)
                  dy = staux(2,n,k) - staux(2,m,j)
                  dz = staux(3,n,k) - staux(3,m,j)
                  !
                  !Periodic adjusted distance
                  sdist(1) = dx - DNINT( dx/(DBLE(mic)) ) * DBLE(mic)
                  sdist(2) = dy - DNINT( dy/(DBLE(mic)) ) * DBLE(mic)
                  sdist(3) = dz - DNINT( dz/(DBLE(mic)) ) * DBLE(mic)
                  !
                  !convert from scaled to real
                  call s_to_r( sdist, rdist, aprim )
                  r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  !
                  !establish which bin this distance is in
                  nbin = nint(dsqrt(r2)*res)
                  !
                  !if we are at the end of binNum
                  if (nbin > binNum) nbin = binNum
                  !
                  !update the temp g(r)
                  !! $OMP CRITICAL
                  ngdrt(nbin) = ngdrt(nbin) + iAdd
                  !! $OMP END CRITICAL 
               enddo
            enddo
            !! $OMP BARRIER
            !! $OMP END PARALLEL
            !
            !update true unnormalized, raw pair correlation
            do n=1,binNum,1
               ngdr(n) = ngdr(n) + ngdrt(n)
            enddo
            !
            ncount = ncount + 1
            !
         enddo Main_loop
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
         !Total region dimensions, averaged over all loops (it should be the same as
         !individual loop region dimensions if the *.cel file is constant)
         omega = omegaT/DBLE(ncount)
         box = mic*omega**(1.0d0/3.0d0)
         Hbox = box/2.0d0
         res = binNum/ Hbox
         !
         !
         do n=1,binNum-1,1
            !
            !construct some details of normalization and distance  
            r = DBLE(n)/res
            r2 = r**2.0d0
            vshell = 4.0 * pi * r2 / res
            if(atom1 /= atom2) then
               dens = DBLE(nsp(atom1))/omega
            else
               dens = DBLE(nsp(atom1) - 1)/omega
            endif
            !

            !Calculate the final, normalized, value of g(r)! Please note that the
            !normalization constant (ncount*mic3*norm*nsp(atom1)*dens*vshell)...
            !  ncount*mic3 = number of steps and number of extended shells
            !  norm = external normalization constant (default is 1.0)
            !  nsp(atom1)*dens = number of pairs
            !  vshell = volume of the infinitesimal shell
            gofr = ngdr(n)/ &     
            &( DBLE(ncount*mic3*norm*nsp(atom1)*dens*vshell) )
            !
            write(3,*) (r*0.52917), gofr
            !
         enddo
         !
         write(*,*)
         write(*,fmt='(1X," g(r) for atom1 : ", I2, " and atom2 : ",I2, " recorded in ", A)')  &
                     atom1, atom2, trim(fileOut)
         write(*,*)
         write(*,*) ' Program complete!'
         !
      end subroutine final
      ! 
      !
      !
      !******************************************************
      !------------------------------------------------------
      !Apply periodic boundary conditions to wrap 
      ! coordinates around the origin
      !------------------------------------------------------
      !******************************************************
      subroutine PBCS(rx, ry, rz, m)

         implicit none

         real(DP), intent(inout)    :: rx, ry, rz
         integer, intent(in)        :: m
         integer                    :: mic

         mic = DBLE(m)

         rx = rx - NINT(rx/mic)*mic
         ry = ry - NINT(ry/mic)*mic
         rz = rz - NINT(rz/mic)*mic
   
         return

      end subroutine PBCS
         
 END PROGRAM new_rdf 
