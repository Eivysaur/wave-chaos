! Compile with gfortran -I/usr/local/include file.f90 -llapack95 -llapack -lopenblas -lpthread -lgfortran
! Needs OpenBLAS and OpenMP
program bucketmodel
  use pbucket
  use functions
  use mpi
  !$ use omp_lib
  implicit none

  integer(wp)                           :: hours, minutes, seconds, ch, j, k, wint, M, evenM, oddM
  real(dp)                              :: timerstart, timerstop, tgf, tgi, tli, tlf, tg, tl, time
  real(wp)                              :: sumSSdag, offdiag, ky, E
  real(wp), dimension(nx,ny)            :: bucketshape
  real(wp), allocatable                 :: neighdist(:), evenneighdist(:), oddneighdist(:)
  complex(wp), dimension(nx,ny)         :: psi
  complex(wp), dimension(N)             :: linearpsi
  complex(wp), allocatable              :: greenpotential(:,:), smatrix(:,:), eigvals(:), ssdag(:,:), minissdag(:,:)
  complex(wp), allocatable              :: evensmatrix(:,:), oddsmatrix(:,:), eveneigvals(:), oddeigvals(:), lucket(:)
  character(len=90)                     :: filename
  integer(sp), allocatable              :: seed(:)
  integer(sp)                           :: seedlen

  ! MPI variables
  integer(wp)   :: to_node, M_node
  integer(sp)   :: ich, fch, error, errorcode, nodeid, nnodes
  integer       :: status(MPI_status_size)

  ! Initialize MPI.
  call MPI_Init(error)
  ! Get the number of processes.
  call MPI_Comm_size(MPI_COMM_WORLD, nnodes, error)
  ! Get the individual process ID.
  call MPI_Comm_rank(MPI_COMM_WORLD, nodeid, error)
  if (nodeid .eq. masternode) then
     print *, "____________________________________________________________________________"
     print *, ""
     write(*,'(A,I0)') "Nodes running:         ", nnodes
     write(*,'(A,I0)') "OMP threads available: ", omp_get_max_threads()
     print *, ""
  end if
  
  ! Time the program (OpenMP is in use)
  !$ timerstart = omp_get_wtime()
  tg = 0._dp
  tl = 0._dp
  
  ! Allocate matrix with Green's function and potential
  allocate(greenpotential(N,N))

  ! Make the bucket
  call makeBucket(bucketshape,nodeid)

  ! Set discretization parameters
  dx = width/(nx-1)
  dy = height/(ny-1)
  dA = dx*dy
  if (wp .eq. sp) then
     explim = -50
  elseif (wp .eq. dp) then
     explim = -700
  end if

  if (nodeid .eq. masternode) then
     ! Prepare a files for nearest neighbor data
     open(unit = 39, file = "tnearest.txt", form='formatted')
     write(39, FMT='(1A16)') "angdist"
     open(unit = 40, file = "tnearestsplitandcombined.txt", form='formatted')
     write(40, FMT='(1A16)') "angdist"
  end if

  ! Get a seed for random number generation
  if (rndpotential .eqv. .true.) then
     call getRandomSeed(seed,seedlen,nnodes,nodeid)
  end if
     
  ! Loop over energies expressed as a wavelength when p=1
  do wint = 1, wlsteps
     ! Update wavelength
     if (wlsteps .ne. 1) then
        wl = wli + (wint-1)*(wlf-wli)/(wlsteps-1)
     else
        wl = wli
     end if

     ! Prepare files
     if (nodeid .eq. masternode) then
        print *, "____________________________________________________________________________"
        print *, ""
        write(*,FMT='(A,F0.2,A)') "wl is: ", wl*1e9, " nm"

        if (wlsteps .eq. 1) then
           ! Prepare files for sumAm data
           write(filename, '("tsumAm-", I0, "x", I0, ".txt")') nx, ny
           open(unit = 26, file = filename, form='formatted')
           write(26, FMT='(4A15)') "inch", "sumAm", "eff"

           ! Prepare a file for S matrix elements to be stored columnwise
           write(filename, '("tsmatrix-" I0, "x", I0, ".txt")') nx, ny
           open(unit = 76, file = filename, form='formatted')
           write(76, FMT='(5A15)') "x", "y", "real", "imag", "abssquare"

           ! Prepare a file for eigenvalues
           write(filename, '("teigenvalues-" I0, "x", I0, ".txt")') nx, ny
           open(unit = 123, file = filename, form='formatted')
           write(123, FMT='(2A15)') "real", "imag"

           ! Prepare a file for SSdag
           write(filename, '("tSSdag-" I0, "x", I0, ".txt")') nx, ny
           open(unit = 28, file = filename, form='formatted')
           write(28, FMT='(4A15)') "y", "x", "SSdag"

           ! Store sub SSdag (even and odd)
           write(filename, '("tSSdag-even-", I0, "x", I0, ".txt")') nx, ny
           open(unit = 129, file = filename, form='formatted')
           write(129, FMT='(4A15)') "y", "x", "SSdag"
           write(filename, '("tSSdag-odd-", I0, "x", I0, ".txt")') nx, ny
           open(unit = 130, file = filename, form='formatted')
           write(130, FMT='(4A15)') "y", "x", "SSdag"
        else
           ! Prepare files for sumAm data
           write(filename, '("tsumAm", "-wl", F0.1, ".txt")') wl*1e9
           open(unit = 26, file = filename, form='formatted')
           write(26, FMT='(4A15)') "inch", "sumAm", "eff"

           ! Prepare a file for S matrix elements to be stored columnwise
           write(filename, '("tsmatrix", "-wl", F0.1, ".txt")') wl*1e9
           open(unit = 76, file = filename, form='formatted')
           write(76, FMT='(5A15)') "x", "y", "real", "imag", "abssquare"

           ! Prepare a file for eigenvalues
           write(filename, '("teigenvalues-wl", F0.1, ".txt")') wl*1e9
           open(unit = 123, file = filename, form='formatted')
           write(123, FMT='(2A15)') "real", "imag"

           ! Prepare a file for SSdag
           write(filename, '("tSSdag-wl", F0.1, ".txt")') wl*1e9
           open(unit = 28, file = filename, form='formatted')
           write(28, FMT='(4A15)') "y", "x", "SSdag"

           ! Store sub SSdag (even and odd)
           write(filename, '("tSSdag-even-", F0.1, ".txt")') wl*1e9
           open(unit = 129, file = filename, form='formatted')
           write(129, FMT='(4A15)') "y", "x", "SSdag"
           write(filename, '("tSSdag-odd-", F0.1, ".txt")') wl*1e9
           open(unit = 130, file = filename, form='formatted')
           write(130, FMT='(4A15)') "y", "x", "SSdag"
        end if
     end if

     ! Find number of open channels based on initial wavelength and p=1
     ky = 2._wp*PI/wl
     E  = sqrt(ky**2+(PI/width)**2)
     M  = floor(E*width/PI)

     ! Allocate (every node)
     allocate(smatrix(M,M))
     allocate(ssdag(M,M))
     allocate(neighdist(M))
     allocate(eigvals(M))
     allocate(lucket(N))

     ! Divide work among the nodes
     if (nodeid .eq. masternode) then
        write(*, FMT='(A16,I0)')    "#open channels: ", M
        write(*, FMT='(A16,E11.4)') "The energy is:  ", E
        M_node = M/nnodes

        ! Check that the number of nodes is good
        if (mod(M,nnodes) .ne. 0) then
           print *, "Allocate a number of nodes that is a multiple of the number of open channels"
           call MPI_abort(MPI_comm_world,errorcode,error)
           call exit(1)
        end if

        ! Masternode splits workload for slavenodes
        do to_node = 1, nnodes-1
           ich = to_node*M_node+1
           fch = ich+M_node-1

           ! Send ich and fch to the slavenodes
           call MPI_send(ich,1,MPI_int,to_node,send_tag,MPI_comm_world,error)
           call MPI_send(fch,1,MPI_int,to_node,send_tag,MPI_comm_world,error)
        end do
        ich = 1
        fch = M_node
        write(*,'(A,I0,A,I0,A,I0)') "I'm the node-", nodeid," and handle channels ", ich, "-", fch
     else ! all except masternode
        ! Recieve ich and fch
        call MPI_recv(ich,1,MPI_int,masternode,MPI_any_tag,MPI_comm_world,status,error)
        call MPI_recv(fch,1,MPI_int,masternode,MPI_any_tag,MPI_comm_world,status,error)
        write(*,'(A,I0,A,I0,A,I0)') "I'm the node-", nodeid," and handle channels ", ich, "-", fch
     end if

     ! Fill bucket
     if (rndpotential .eqv. .true.) then
        call fillBucket(bucketshape,lucket,E,nodeid,seed)
     else
        call fillBucket(bucketshape,lucket,E,nodeid)
     end if

     ! Loop over incoming channels. The workload is shared among the nodes
     do ch = ich, fch
        ! Recalculate ky
        ky = sqrt(E**2-(ch*PI/width)**2)
        
        ! Construct
        !$ tgi = omp_get_wtime()
        call makeGreenIncoming(bucketshape,ky,greenpotential,linearpsi,ch,M,E,nodeid,lucket)
        !$ tgf = omp_get_wtime()
        tg = tg+tgf-tgi

        ! Store incoming wave and Green's function matrix (optional)
        ! call storeInGreen(linearpsi,greenpotential,ch,wl)

        ! Solve the linear system; greenpotential*psi=linearpsi, using LAPACK
        !$ tli = omp_get_wtime()
        call linSolve(greenpotential,linearpsi,nodeid,ch)
        !$ tlf = omp_get_wtime()
        tl = tl+tlf-tli

        ! Map the solution linearpsi (previously incoming wave) to original matrix shape
        call storeOut(linearpsi,psi,ch)
        
        ! Calculate S matrix
        call makeSmatrix(bucketshape,psi,smatrix,ch,M,E,ky,nodeid,nnodes)

        if (nodeid .eq. masternode) then
           print *, "__________________________________________________________________________"
        end if
     end do

     ! The master node does the S matrix statistics
     if (nodeid .eq. masternode) then
        print *, "__________________________________________________________________________"
        print *, ""
        write(*,'(A,F0.1)') "Masternode is doing S-matrix statistics for wl = ", wl*1e9
        
        close(26)

        ! Write S matrix to file
        do j = 1, M
           do k = 1, M
              write(76,*) j, k, real(smatrix(j,k)), aimag(smatrix(j,k)), abs(smatrix(k,j))**2
           end do
        end do

        close(76)

        ! Check if S-matrix is unitary
        call checkUnitarity(smatrix,M,sumSSdag,offdiag)

        ssdag = matmul(conjg(transpose(smatrix)),smatrix)
        do j = 1, M
           do k = 1, M
              write(28,*) j, k, abs(ssdag(k,j))**2
           end do
        end do
        close(28)

        ! Calculate eigenvalues of scattering matrix
        call sEigenValues(smatrix,eigvals,M)

        ! Find nearest neighbors distances between eigenvalues
        call nearestNeighborDist(eigvals,neighdist,M)

        ! Write nearest neighbors to file (file is opened at beginning)
        do j = 1, M
           write(39,*) neighdist(j)
           ! print *, neighdist(j)
        end do

        !! Do parity splitting of the S matrix
        if (mod(M,2) .eq. 0) then
           write(*,FMT='(A35,I0)') "Dimension of even parity S-matrix: ", M/2
           write(*,FMT='(A35,I0)') "Dimension of odd parity S-matrix:  ", M/2
           evenM = M/2
           oddM  = M/2
        else
           write(*,FMT='(A35,I0)') "Dimension of even parity S-matrix: ", M/2
           write(*,FMT='(A35,I0)') "Dimension of odd parity S-matrix:  ", M/2+1
           evenM = M/2
           oddM  = M/2+1
        end if
        allocate(evensmatrix(evenM,evenM))
        allocate(oddsmatrix(oddM,oddM))

        ! Split S-matrix into even and odd parity submatrices
        call splitSmatrix(smatrix,evensmatrix,oddsmatrix,evenM,oddM)

        ! Even parity
        print *, "Doing even parity..."
        allocate(evenneighdist(evenM))
        allocate(eveneigvals(evenM))
        call checkUnitarity(evensmatrix,evenM,sumSSdag,offdiag)

        minissdag = matmul(conjg(transpose(evensmatrix)),evensmatrix)
        do j = 1, evenM
           do k = 1, evenM
              write(129,*) j, k, abs(minissdag(k,j))**2
           end do
        end do
        close(129)
        
        call sEigenValues(evensmatrix,eveneigvals,evenM)
        call nearestNeighborDist(eveneigvals,evenneighdist,evenM)
        do j = 1, evenM
           write(40,*) evenneighdist(j)
           write(123,*) real(eveneigvals(j)), aimag(eveneigvals(j))
        end do

        ! Odd parity
        print *, "Doing even parity..."
        allocate(oddneighdist(oddM))
        allocate(oddeigvals(oddM))
        call checkUnitarity(oddsmatrix,oddM,sumSSdag,offdiag)

        minissdag = matmul(conjg(transpose(oddsmatrix)),oddsmatrix)
        do j = 1, oddM
           do k = 1, oddM
              write(130,*) j, k, abs(minissdag(k,j))**2
           end do
        end do
        close(130)
        
        call sEigenValues(oddsmatrix,oddeigvals,oddM)
        call nearestNeighborDist(oddeigvals,oddneighdist,evenM)
        do j = 1, oddM
           write(40,*) oddneighdist(j)
           write(123,*) real(oddeigvals(j)), aimag(oddeigvals(j))
        end do

        ! Deallocate to prepare for a new wavelength (in case sizes change)
        deallocate(evensmatrix,oddsmatrix,eveneigvals,oddeigvals,evenneighdist,oddneighdist)
        print *, "__________________________________________________________________________"
     end if
     ! Deallocate (every node)
     deallocate (lucket,ssdag,smatrix,eigvals,neighdist)
     
     write(*,'(A,I2,A)') "Node-", nodeid, " | waiting for next wavelength"
     call MPI_barrier(MPI_comm_world,error) ! Synchronize the nodes before the next wavelength
  end do
  close(39)
  close(123)
  close(40)
  deallocate(greenpotential)
  if (rndpotential .eqv. .true.) then
     deallocate(seed)
  end if

  ! Write useful info to screen
  if (nodeid .eq. masternode) then
     call cpu_time(timerstop)
     !$ timerstop = omp_get_wtime()
     time    = timerstop-timerstart
     hours   = int8(time/3600)
     minutes = int8((time-hours*3600)/60)
     seconds = int8(time-hours*3600-minutes*60)
     print *, "==============================================================="
     print *, "Run time:               ", timerstop-timerstart
     print *, "Time spent building T:  ", tg
     print *, "Time spent linSolving:  ", tl  
     print *, "#closed channels        ", clch
     print *, "==============================================================="
     write(*,FMT='(T2,A10)', advance="no") "Run time: "
     write(*,FMT='(I2,A6)') hours, " hours"
     write(*,FMT='(T12,I2,A8)') minutes, " minutes"
     write(*,FMT='(T12,I2,A8)') seconds, " seconds"
     print *, "==============================================================="

     ! Write logbook
     call logParameters(hours,minutes,seconds,tg,tl)
  end if
  call MPI_barrier(MPI_comm_world,error)
  call MPI_finalize(error)
end program bucketmodel
