module functions
  use pbucket
  use mpi
  !$ use omp_lib
  implicit none

contains

  !_______________________________________________________________________________________________
  !
  ! Chose surface and make bucket
  !_______________________________________________________________________________________________
  subroutine makeBucket(bucket,nodeid)
    implicit none
    integer(wp)                                    :: j, k, tip, hwidth, st, sl, rows
    integer(sp), intent(in)                        :: nodeid
    real(wp)                                       :: slope
    real(wp), dimension(nx,ny), intent(inout)      :: bucket
    real(wp), dimension(nx)                        :: surf, surf2
    real(wp), allocatable                          :: ps(:,:), peaks(:,:)
    integer(wp), dimension(3)                      :: vol

    surf    = 0._wp
    surf2   = 0._wp
    bucket  = 0._wp
    ! Choose surface
    select case (surface)
    case (1) ! Flat potential
       do j = 1, nx
          surf(j) = ny
       end do
       ! do j = 1, nx/4
       !    surf(j) = 0._wp
       ! end do
       ! do j = nx/4+1, 3*nx/4
       !    surf(j) = ny-1
       ! end do
       ! do j = 3*nx/4+1, nx
       !    surf(j) = 0._wp
       ! end do
       call fillLayers(surf,bucket)

    case (2) ! Sine potential
       if (amplt .gt. ny/2-1) then
          print *, "Amplt is too large"
          call exit(1)
       end if
       do j = 1, nx
          surf(j) = amplt*sin(bumps*(j-1)*PI/(nx-1))+amplt
       end do
       call fillLayers(surf,bucket)

    case (3) ! Oblique potential
       if (amplt .gt. ny) then
          print *, "amplt is too large!"
          STOP
       end if
       do j = 1, nx
          surf(j) = (j-1)*(amplt/(nx-1))
       end do
       call fillLayers(surf,bucket)

    case (4) ! House with 2 sub-houses
       if (mod(ny,4) .ne. 0 .and. mod(nx,4) .ne. 0) then
          print *, "Height and width of bucket need to be multiples of 4"
          call exit(1)
       end if

       ! House
       tip = nint(nx*ratioB)              
       do j = 1, tip
          surf(j) = (ny*heightB/tip)*(j-1)
       end do
       do j = tip+1, nx
          surf(j) = (surf(tip)/(nx-(tip+1)))*(nx-j)
       end do

       ! 1st sub-house
       tip = nint((nx/2)*ratioL)
       do j = 1, tip
          surf2(j) = (ny*heightL/tip)*(j-1)
       end do
       do j = tip+1, nx/2
          surf2(j) = (surf2(tip)/(nx/2-(tip+1)))*((nx/2)-j)
       end do
       ! 2nd sub-house
       tip = nint((nx/2)*ratioR)       
       do j = (nx/2)+1, (nx/2)+tip
          surf2(j) = (ny*heightR/tip)*((j-1)-nx/2)
       end do
       do j = (nx/2)+tip+1, nx
          surf2(j) = (surf2(tip)/(nx/2-(tip+1)))*(nx-j)
       end do
       call fillLayers(surf,bucket,surf2)

    case (5) ! Circle top with half-pipe inside
       if (mod(ny,2) .ne. 0 .and. mod(nx,2) .ne. 0) then
          print *, "Height and width of bucket needs to be even numbers"
          call exit(1)
       end if

       ! Dome
       do j = 1, nx/2
          surf(j) = sqrt(real(((ny/2))**2-(j-(nx/2)-1)**2))
       end do
       do j = nx/2+1, nx
          surf(j) = sqrt(real(((ny/2))**2-((nx/2)-j)**2))
       end do

       ! Half-pipe
       do j = 1, nx/2
          surf2(j) = real((ny/2))-sqrt(real(((ny/2))**2-(j-(nx/2)-1)**2))
       end do
       do j = nx/2, nx
          surf2(j) = real((ny/2))-sqrt(real(((ny/2))**2-((nx/2)-j)**2))
       end do
       call fillLayers(surf,bucket,surf2)

    case (6) ! Dome top with 2 dome houses
       if (mod(ny,4) .ne. 0 .and. mod(nx,4) .ne. 0) then
          print *, "Height and width of bucket need to be multiples of 4"
          call exit(1)
       end if
       ! Dome
       do j = 1, nx/2
          surf(j) = sqrt(real(((ny/2))**2-(j-(nx/2)-1)**2))
       end do
       do j = nx/2+1, nx
          surf(j) = sqrt(real(((ny/2))**2-((nx/2)-j)**2))
       end do

       ! 1st dome house
       do j = 1, nx/4
          surf2(j) = ecc1*sqrt(real(((ny/4))**2-(j-(nx/4)-1)**2))
       end do
       do j = nx/4+1, nx/2
          surf2(j) = ecc1*sqrt(real(((ny/4))**2-((nx/4)-j)**2))
       end do
       ! 2nd dome house
       do j = nx/2+1, (3*nx)/4
          surf2(j) = ecc2*sqrt(real(((ny/4))**2-((j-1-nx/2)-(nx/4))**2))
       end do
       do j = (3*nx)/4+1, nx
          surf2(j) = ecc2*sqrt(real(((ny/4))**2-((3*nx/4)-j)**2))
       end do
       call fillLayers(surf,bucket,surf2)

    case (7) ! Sine with small houses
       if (mod(ny,4) .ne. 0 .and. mod(nx,4) .ne. 0) then
          print *, "Height and width of bucket need to be multiples of 4"
          call exit(1)
       end if
       if (amplt .gt. ny/2-1) then
          print *, "Amplt is too large"
          call exit(1)
       end if
       do j = 1, nx
          surf(j) = amplt*sin(3._wp*(j-1)*PI/(nx-1))+amplt
       end do
       ! 1st sub-house
       tip = nint((nx/2)*ratioL)
       do j = 1, tip
          surf2(j) = (ny*heightL/tip)*(j-1)
       end do
       do j = tip+1, nx/2
          surf2(j) = (surf2(tip)/(nx/2-(tip+1)))*((nx/2)-j)
       end do
       ! 2nd sub-house
       tip = nint((nx/2)*ratioR)       
       do j = (nx/2)+1, (nx/2)+tip
          surf2(j) = (ny*heightR/tip)*((j-1)-nx/2)
       end do
       do j = (nx/2)+tip+1, nx
          surf2(j) = (surf2(tip)/(nx/2-(tip+1)))*(nx-j)
       end do
       call fillLayers(surf,bucket,surf2)

    case (8) ! Adjustable saw potential
       if (mod(nx,nhus) .ne. 0) then
          print *, "Width need to be a multiple of number of houses (nhus)"
          call exit(1)
       end if
       hWidth = nx/nHus
       write(*,FMT='(A33,I0)') "Number of grid points per house: ", hWidth
       do k = 1, nHus
          tip = nint((hWidth)*nRatio(k))
          if (tip .eq. 0) then
             st = (k-1)*(hWidth)+1
             sl = (hWidth)*k-tip
             write(*,FMT='(A24,I0,A5,F5.2,A8)') "Left angle of house nr. ", k," is: ",&
                  &(atan(real((height*surf(tip)/ny)/(width*tip/(hWidth)))))*180._wp/PI, " degrees"
             do j = st, sl+tip
                surf(j) = ny*nHeight(k)-ny*nHeight(k)/(hWidth)*((j-1)-(hWidth)*(k-1))! (delta y/ delta x)*j
                ! print *, j, surf(j)
             end do
          else
             st = (k-1)*(hWidth)+1
             sl = (hWidth)*k-tip
             do j = st, tip+(k-1)*(hWidth)
                surf(j) = (ny*nHeight(k)/(tip-1))*((j-1)-(hWidth)*(k-1))
                ! print *, j, surf(j)
             end do
             do j = st+tip, sl+tip
                surf(j) = (surf(tip)/((hWidth)-(tip+1)))*((hWidth)*k-j)
                ! print *, j, surf(j)
             end do
             write(*,FMT='(A24,I0,A5,F5.2,A8)') "Left angle of house nr. ", k," is: ",&
                  &(atan(real((height*surf(tip)/ny)/(width*tip/nx))))*180._wp/PI, " degrees"
          end if
       end do
       call fillLayers(surf,bucket)

    case (9) ! Chaotic dip
       if (amplt .gt. ny) then
          print *, "amplt is too large!"
          STOP
       end if
       do j = 1, nint(0.1766_wp*nx)
          !surf(j) = amplt-(j-1)*(amplt/nint(0.18_wp*nx))
          surf(j) = amplt-(amplt/(nint(0.1766_wp*nx)-1))*(j-1)!(j-1)*(amplt/nx-1)
          ! print *, surf(j)
       end do
       do j = 1, nx-nint(0.1766_wp*nx)
          surf(j+nint(0.1766_wp*nx)) = (j-1)*(amplt/(nx-nint(0.1766_wp*nx)-1))
          ! print *, surf(j+nint(0.176_wp*nx))
       end do
       call fillLayers(surf,bucket)

    case(10) ! Read tpotential1.txt
       open(unit=42,file="tpotential1-7.txt",status="old")
       call readFileLength(42,rows)
       allocate(ps(rows,2))
       allocate(peaks(rows,2))

       rewind(42)
       read(42,*)
       do j = 1, rows
          read(42,*) ps(j,:)
       end do
       close(42)

       ! Find mesh coordinates of the peaks of the potential
       ! and put these in surf
       surf = 0
       surf(1) = nint(ps(1,2)*ny)
       peaks(1,1) = 1
       peaks(1,2) = nint(ps(1,2)*ny)
       do k = 2, rows
          ! print *, nint(ps(k,1))*nx
          if (nint(ps(k,1)*nx) .eq. 0) then
             surf(1) = nint(ps(k,2)*ny)
             peaks(k,1) = 1
          else
             surf(nint(ps(k,1)*nx)) = nint(ps(k,2)*ny)
             peaks(k,1) = nint(ps(k,1)*nx)
          end if
          peaks(k,2) = nint(ps(k,2)*ny)
       end do
       ! print *, "x-coords:"
       ! do j = 1, rows
       !    print *, peaks(j,2)
       ! end do

       ! Connect the "dots"
       do k = 2, rows
          if (peaks(k,1)-peaks(k-1,1) .ne. 0) then
             slope = (peaks(k,2)-peaks(k-1,2))/(peaks(k,1)-peaks(k-1,1))
             ! print *, slope
             do j = int(peaks(k-1,1)), int(peaks(k,1))
                surf(j) = peaks(k-1,2)+slope*(j-peaks(k-1,1))
             end do
          else
             surf(int(peaks(k,1))) = peaks(k,2)
          end if
       end do

       ! Scale surf to film thickness and mesh
       surf = surf*ny/maxval(surf)
       height = ps(maxloc(peaks(:,2),1),2)*1e-6_wp
       write(*, FMT='(A22,F0.4,A3)')    "New bucket height:    ", height*1e9, " nm"

       ! Check surf
       ! do j = 1, nx
       !    print *, j, surf(j)
       ! end do

       ! Call function to fill the area under surf
       call fillRandomLayer(surf,bucket)

    case (11) ! Film with dome on top
       if ((mod(ny,2) .ne. 0) .and. (mod(nx/dupes,2) .ne. 0)) then
          print *, "Height and width/dupes of bucket need to be multiples of 2"
          call exit(1)
       end if

       ! print *, "ny needs to be at least: ", ecc1*(ny-1)/height, ecc1*(ny-1)/height+yoff*(ny-1)/height
       ! print *, ecc1*(ny-1)*dy/height, yoff*(ny-1)*dy/height
       ! call exit(1)

       ! Dome(s)
       do k = 1, dupes
          do j = 1, (nx/dupes)/2
             surf((k-1)*(nx/dupes)+j) = ecc1*(ny-1)/height*sqrt(1._wp-((j-1-(nx/dupes)/2)/dble((nx/dupes)/2))**2)+&
                  yoff*(ny-1)/height
             ! print *, "surf1: ", j, surf(j)
          end do
          do j = (nx/dupes)/2+1, (nx/dupes)
             surf((k-1)*(nx/dupes)+j) = ecc1*(ny-1)/height*sqrt(1._wp-((j-(nx/dupes)/2)/dble((nx/dupes)/2))**2)+&
                  yoff*(ny-1)/height
             ! print *, "surf1: ", j, surf(j)
          end do
       end do
       ! print *, ecc1, ecc1*dble((ny-1)/height)
       ! print *, yoff, yoff*dble((ny-1)/height)

       ! Film
       do j = 1, nx
          surf2(j) = yoff*dble((ny-1)/height)
          ! print *, "surf2: ", surf2(j)
       end do
       call fillLayers(surf,bucket,surf2)
       
    end select

    vol = 0
    ! Check volume
    do j = 1, ny
       do k = 1, nx
          if (bucket(k,j) .eq. 1) then
             vol(1) = vol(1)+1
          elseif (bucket(k,j) .eq. 2) then
             vol(2) = vol(2)+1
          else
             vol(3) = vol(3)+1
          end if
       end do
    end do

    ! Store shape
    if (nodeid .eq. masternode) then
       open(unit = 25, file = "tbucket.txt", form='formatted')
       write(25, FMT='(4A15)') "ny", "nx", "bucket"
       do j = 1, ny
          do k = 1, nx
             write(25,*) j, k, bucket(k,j)
          end do
       end do
       close(unit = 25)
    end if
  end subroutine makeBucket

  ! Fill layers
  subroutine fillLayers(surf,bucket,surf2)
    integer(wp)                                    :: j, k, level, step
    real(wp), dimension(nx,ny), intent(inout)      :: bucket
    real(wp), dimension(nx), intent(in)            :: surf
    real(wp), dimension(nx), intent(in), optional  :: surf2
    logical                                        :: inside

    ! Fill 1st layer
    level = ny - nint(maxval(surf)-minval(surf))
    do j = 1, nx
       step = nint(surf(j))
       if (step .eq. 0) then
          step = 1
       end if
       ! bucket(j,level+step) = 1._wp ! Rework any potential that requires level, it's confusing
       bucket(j,step) = 1._wp
       ! bucket(j,ny/2) = 1._wp
    end do
    ! Fill inside
    do j = 1, nx
       inside = .true.
       do k = 1, ny
          if (bucket(j,k) .eq. 0._wp .and. inside .eqv. .true.) then
             bucket(j,k) = 1._wp
          elseif (bucket(j,k) .ne. 0._wp) then
             inside = .false.
          end if
       end do
    end do

    ! Fill 2nd layer if required
    if (present(surf2)) then
       do j = 1, nx
          step = nint(surf2(j))
          bucket(j,step) = 2._wp
          ! bucket(j,step+1) = 2._wp ! only for flat methinks
       end do
       do j = 1, nx
          do k = 1, ny
             if (bucket(j,k) .ne. 2._wp) then
                bucket(j,k) = 2._wp
             else
                exit
             end if
          end do
       end do
    end if
  end subroutine fillLayers

  ! Fill random layer width simulation height specified
  subroutine fillRandomLayer(surf,bucket)
    integer(wp)                                   :: j, k, step
    real(wp), dimension(nx,ny), intent(inout)     :: bucket
    real(wp), dimension(nx), intent(in)           :: surf
    logical                                       :: inside

    ! Fill inside
    ! print *, "level is: ", level
    do j = 1, nx
       step = nint(surf(j))
       ! print *, j, step
       bucket(j,step) = 1._wp
    end do
    ! Fill inside
    do j = 1, nx
       inside = .true.
       do k = 1, ny
          if (bucket(j,k) .eq. 0._wp .and. inside .eqv. .true.) then
             bucket(j,k) = 1._wp
          elseif (bucket(j,k) .ne. 0._wp) then
             inside = .false.
          end if
       end do
    end do
  end subroutine fillRandomLayer

  !_______________________________________________________________________________________________
  !
  ! Constructs the bucket (only masternode should use it)
  !_______________________________________________________________________________________________
  subroutine fillBucket(bucketshape,lucket,E,nodeid,seed)
    integer(wp)                                 :: j, k, it
    integer(sp), intent(in), optional           :: nodeid
    real(wp)                                    :: rnd
    real(wp), intent(in)                        :: bucketshape(:,:), E
    real(wp), dimension(N)                      :: x, y
    complex(wp), intent(inout)                  :: lucket(:)
    integer(sp), intent(in), optional           :: seed(:)
    character(len=90)                           :: filename

    lucket = 0._wp
    it = 1
    if (rndpotential .eqv. .false.) then
       do j = 1, nx
          do k = 1, ny
             x(it) = dx*(j-1)
             y(it) = dy*(k-1)
             if (bucketshape(j,k) .eq. 1) then
                lucket(it) = (1-rind**2)*E**2
             elseif (bucketshape(j,k) .eq. 2) then
                lucket(it) = (1-rind2**2)*E**2
             else
                lucket(it) = 0._wp
             end if
             it = it + 1
          end do
       end do
    else ! randomized index of refraction
       call random_seed(put=seed)
       ! Store shape
       write(filename, '("tbucketpotential-node", I0, ".txt")') nodeid
       open(unit = 25, file = "tbucketpotential.txt", form='formatted')
       write(25, FMT='(4A15)') "ny", "nx", "bucketpotential"
       do j = 1, nx
          do k = 1, ny
             x(it) = dx*(j-1)
             y(it) = dy*(k-1)
             if (bucketshape(j,k) .eq. 1) then
                call random_number(rnd)
                rnd = rnd+1._wp ! random nubmer in interval [1,2]
                lucket(it) = (1-rnd**2)*E**2
                write(25,*) j, k, abs(lucket(it))**2, rnd
             elseif (bucketshape(j,k) .eq. 2) then
                call random_number(rnd)
                rnd = rnd+1._wp ! random nubmer in interval [1,2]
                lucket(it) = (1-rnd**2)*E**2
                write(25,*) j, k, abs(lucket(it))**2, rnd
             else
                lucket(it) = 0._wp
                write(25,*) j, k, abs(lucket(it))**2, 1._wp
             end if
             it = it + 1
          end do
       end do
       close(unit = 25)
    end if

  end subroutine fillBucket
  
  !_______________________________________________________________________________________________
  !
  ! Constructs Green's function matrix (by calling green) and the incoming wave
  !_______________________________________________________________________________________________
  subroutine makeGreenIncoming(bucketshape,ky,greenpotential,linearpsi,inch,M,E,nodeid,lucket)
    implicit none
    integer(wp)                                 :: j, k, it, tracker, int
    integer(wp), intent(in)                     :: inch, M
    integer(sp), intent(in)                     :: nodeid
    real(dp)                                    :: tgi, tgf
    real(wp), intent(in)                        :: ky, E
    real(wp), dimension(M)                      :: B
    real(wp), dimension(clch)                   :: A
    real(wp), dimension(N)                      :: x, y
    real(wp), intent(in)                        :: bucketshape(:,:)
    complex(wp), intent(in)                     :: lucket(:)
    complex(wp), intent(inout)                  :: linearpsi(:), greenpotential(:,:)

    ! Ensure a clean start
    greenpotential = 0._wp

    ! Generate coordinate strips
    it = 1
    do j = 1, nx
       do k = 1, ny
          x(it) = dx*(j-1)
          y(it) = dy*(k-1)
          it = it + 1
       end do
    end do

    ! Initialize some variables for green to save computation time
    do k = 1, M
       B(k) = sqrt(E**2 - k**2*PI**2/width**2)
    end do
    do k = M+1, M+clch
       A(k-M) = sqrt(k**2*PI**2/width**2 - E**2)
    end do

    ! This gives <denominator> number of progress updates
    tracker = N/prog

    ! Call green to make matrix gmat (T) and incoming wave iowave (W)
    !$ tgi = omp_get_wtime()
    !$OMP parallel ! num_threads(1)
    ! !$ write(*,FMT='(A4,I2,A13)') "Red ", omp_get_thread_num(), " standing by"
    !$OMP do schedule(dynamic) ! Divides the work in x equal parts
    do j = 1, N
       ! if (mod(j,tracker) .eq. 0) then
       !    print *, "Progress: ", (j*100)/N, "%"
       ! end if
       linearpsi(j) = (exp(-i*ky*y(j))-exp(i*ky*y(j)))*sin(inch*PI*x(j)/width) ! Sine modulated
       !if (lucket(j) .ne. 0._wp) then
       do k = j, N
          greenpotential(k,j) = -green(M,B,A,x(k),x(j),y(k),y(j))
          if (j .ne. k) then
             greenpotential(j,k) = greenpotential(k,j)*lucket(k)*dA
          end if
          greenpotential(k,j) = greenpotential(k,j)*lucket(j)*dA
       end do
       !end if
    end do
    !$OMP end do
    !$OMP end parallel
    !$ tgf = omp_get_wtime()
    write(*,'(A,I2,A,I2,A,F7.1)') "Node-", nodeid, " | inch-", inch,&
         " | seconds spent building greenpotential: ", tgf-tgi
    do j = 1, N
       greenpotential(j,j) = greenpotential(j,j)+1
    end do
  end subroutine makeGreenIncoming

  !_______________________________________________________________________________________________
  !
  ! Constructs Green's function matrix only
  !_______________________________________________________________________________________________
  subroutine makeGreenPotential(bucketshape,greenpotential,M,E)
    implicit none
    integer(wp)                                 :: j, k, it, tracker
    integer(wp), intent(in)                     :: M
    real(wp)                                    :: tgi, tgf
    real(wp), intent(in)                        :: E
    real(wp), dimension(M)                      :: B
    real(wp), dimension(clch)                   :: A
    real(wp), dimension(N)                      :: x, y
    real(wp), intent(in)                        :: bucketshape(:,:)
    complex(wp), dimension(N)                   :: lucket
    complex(wp), intent(inout)                  :: greenpotential(:,:)

    ! Ensure a clean start
    greenpotential = 0._wp
    lucket = 0._wp
  
    ! Linearize the bucket and put in the potential
    it = 1
    do j = 1, nx
       do k = 1, ny
          x(it) = dx*(j-1)
          y(it) = dy*(k-1)
          if (bucketshape(j,k) .eq. 1) then
             lucket(it) = (1-rind**2)*E**2
          elseif (bucketshape(j,k) .eq. 2) then
             lucket(it) = (1-rind2**2)*E**2
          else
             lucket(it) = 0._wp
          end if
          it = it + 1
       end do
    end do

    ! Initialize some variables for green to save computation time
    do k = 1, M
       B(k) = sqrt(E**2 - k**2*PI**2/width**2)
    end do
    do k = M+1, M+clch
       A(k-M) = sqrt(k**2*PI**2/width**2 - E**2)
    end do

    ! This gives <denominator> number of progress updates
    tracker = N/prog

    ! Call green to make matrix gmat (T) and incoming wave iowave (W)
    !$ tgi = omp_get_wtime()
    !$OMP parallel ! num_threads(1)
    !$ write(*,FMT='(A4,I2,A13)') "Red ", omp_get_thread_num(), " standing by"
    !$OMP do schedule(dynamic,16) ! Divides the work in x equal parts
    do j = 1, N
       if (mod(j,tracker) .eq. 0) then
          print *, "Progress: ", (j*100)/N, "%"
       end if
       !if (lucket(j) .ne. 0._wp) then
       do k = j, N
          greenpotential(k,j) = -green(M,B,A,x(k),x(j),y(k),y(j))
          if (j .ne. k) then
             greenpotential(j,k) = greenpotential(k,j)*lucket(k)*dA
          end if
          greenpotential(k,j) = greenpotential(k,j)*lucket(j)*dA
       end do
       !end if
    end do
    !$OMP end do
    !$OMP end parallel
    !$ tgf = omp_get_wtime()
    print *, "Build Green's function potential matrix: ", tgf-tgi
    do j = 1, N
       greenpotential(j,j) = greenpotential(j,j)+1
    end do
  end subroutine makeGreenPotential

  !_______________________________________________________________________________________________
  !
  ! Constructs the free solution wavefunction
  !_______________________________________________________________________________________________
  subroutine makeFreeSolution(linearpsi,ky,inch)
    implicit none
    integer(wp)                                 :: j, k, it, tracker
    integer(wp), intent(in)                     :: inch
    real(wp)                                    :: tgi, tgf
    real(wp), intent(in)                        :: ky
    real(wp), dimension(N)                      :: x, y
    complex(wp), intent(inout)                  :: linearpsi(:)

    ! Linearize the bucket and put in the potential
    it = 1
    do j = 1, nx
       do k = 1, ny
          x(it) = dx*(j-1)
          y(it) = dy*(k-1)
          it = it + 1
       end do
    end do

    ! This gives <denominator> number of progress updates
    tracker = N/prog

    !$ tgi = omp_get_wtime()
    !$OMP parallel ! num_threads(1)
    ! !$ write(*,FMT='(A4,I2,A13)') "Red ", omp_get_thread_num(), " standing by"
    !$OMP do
    do j = 1, N
       ! if (mod(j,tracker) .eq. 0) then
       !    print *, "Progress: ", (j*100)/N, "%"
       ! end if
       linearpsi(j) = (exp(-i*ky*y(j))-exp(i*ky*y(j)))*sin(inch*PI*x(j)/width) ! Sine modulated
    end do
    !$OMP end do
    !$OMP end parallel
    !$ tgf = omp_get_wtime()
    print *, "Build free solution: ", tgf-tgi
  end subroutine makeFreeSolution

  !_______________________________________________________________________________________________
  !
  ! Store incoming wave and Green's function matrix
  !_______________________________________________________________________________________________
  subroutine storeInGreen(linearpsi,greenpotential,inch,wl)
    implicit none
    integer(wp)                                    :: j, k, it
    integer(wp), intent(in)                        :: inch    
    real(wp), intent(in)                           :: wl
    complex(wp)                                    :: iwave(nx,ny)
    complex(wp), intent(in)                        :: linearpsi(:)
    complex(wp), intent(in)                        :: greenpotential(:,:)
    character(len=90)                              :: filename

    ! Store incoming
    it = 1
    do j = 1, nx
       do k = 1, ny
          iwave(j,k) = linearpsi(it)
          it         = it + 1
       end do
    end do
    write(filename, '("tincoming", "-inch", I0, "-wl", F0.1, ".txt")') inch, wl*1e9
    open(unit = 25, file = filename, form='formatted')
    write(25, FMT='(5A15)') "x", "y", "absincoming", "real", "imag"
    do j = 1, nx
       do k = 1, ny
          write (25,*) j, k, abs(iwave(j,k))**2, real(iwave(j,k)), aimag(iwave(j,k))
       end do
    end do
    close(unit = 25)

    ! Store Green
    ! if (N .le. 30*30) then
    !    open(unit = 25, file = "tgreen.txt", form='formatted')
    !    write(25, FMT='(4A15)') "x", "y", "green"
    !    do j = 1, N
    !       do k = 1, N
    !          write(25,*) j, k, abs(greenpotential(j,k))**2
    !       end do
    !    end do
    !    close(unit = 25)
    ! end if
  end subroutine storeInGreen

  !_______________________________________________________________________________________________
  !
  ! Solve a linear system Ax=b using ZGESV from LAPACK
  !_______________________________________________________________________________________________
  subroutine linSolve(greenpotential,linearpsi,nodeid,inch)
    implicit none
    real(dp)                                    :: tlapackf, tlapacki
    integer(wp)                                 :: INFO
    integer(wp), dimension(N)                   :: IPIV
    integer(sp), intent(in)                     :: nodeid
    integer(wp), intent(in)                     :: inch
    complex(wp), dimension(N), intent(inout)    :: linearpsi
    complex(wp), dimension(N,N), intent(inout)  :: greenpotential

    INFO = 42
    !$ tlapacki = omp_get_wtime()
    if (wp .eq. dp) then
       call ZGESV(N,1,greenpotential,N,IPIV,linearpsi,N,INFO)
    elseif (wp .eq. sp) then
       call CGESV(N,1,greenpotential,N,IPIV,linearpsi,N,INFO)
    end if
    ! write(*,FMT='(A7,I0)') "Info = ", INFO
    !$ tlapackf = omp_get_wtime()
    write(*,'(A,I2,A,I2,A,F7.1)') "Node-", nodeid, " | inch-", inch,&
         " | seconds spent linsolving: ", tlapackf-tlapacki
  end subroutine linSolve

  !_______________________________________________________________________________________________
  !
  ! Store |psi|^2
  !_______________________________________________________________________________________________
  subroutine storeOut(linearpsi,psi,inch)
    implicit none
    integer(wp)                                      :: j, k, it
    integer(wp), intent(in)                          :: inch
    complex(wp), dimension(nx,ny), intent(out)       :: psi
    complex(wp), dimension(N), intent(in)            :: linearpsi
    character(len=90)                                :: filename

    it = 1
    do j = 1, nx
       do k = 1, ny
          psi(j,k) = linearpsi(it)
          it = it + 1
       end do
    end do

    ! Store wavefunction
    write(filename, '("toutgoing", "-inch", I0, "-wl", F0.1, ".txt")') inch, wl*1e9
    open(unit = 25, file = filename, form='formatted')
    write(25, FMT='(5A15)') "y", "x", "psi", "real", "imag"
    do j = 1, ny
       do k = 1, nx
          write(25,*) j, k, abs(psi(k,j))**2, real(psi(k,j)), aimag(psi(k,j))
       end do
    end do
    close(unit = 25)
  end subroutine storeOut

  !_______________________________________________________________________________________________
  !
  ! Calculate scattering amplitudes
  !_______________________________________________________________________________________________
  subroutine makeSmatrix(bucketshape,psi,smatrix,inch,M,E,kp,nodeid,nnodes)
    implicit none
    integer(wp)                                      :: j, k, l, it, srch, from_node, error
    integer(wp), intent(in)                          :: inch, M
    integer(sp), intent(in)                          :: nodeid, nnodes
    real(wp), intent(in)                             :: E, kp
    real(wp)                                         :: km, srm, B, sumAm
    real(wp), dimension(nx)                          :: xcrd
    real(wp), dimension(ny)                          :: ycrd
    real(wp), intent(in)                             :: bucketshape(:,:)
    complex(wp), intent(in)                          :: psi(:,:)
    complex(wp), dimension(nx,ny)                    :: zbucket
    complex(wp), intent(out)                         :: smatrix(:,:)
    complex(wp), dimension(M)                        :: ams
    integer                                          :: status(MPI_status_size)

    ! Prepare bucket potential
    zbucket = 0._wp
    it = 0
    do j = 1, nx
       do k = 1, ny
          if (bucketshape(j,k) .eq. 1) then
             zbucket(j,k) = (1-rind**2)*E**2
             it = it+1
          elseif (bucketshape(j,k) .eq. 2) then
             zbucket(j,k) = (1-rind2**2)*E**2
             it = it+1
          else
             zbucket(j,k) = 0._wp
          end if
       end do
    end do

    ! Calculate coordinates
    do j = 1, nx
       xcrd(j) = dx*(j-1)
    end do
    do k = 1, ny
       ycrd(k) = dy*(k-1)
    end do

    ! Calculate Am's
    ams = 0._wp
    do l = 1, M
       B  = sqrt(E**2 - (l*PI/width)**2)
       do j = 1, nx
          do k = 1, ny
             ! if (bucket(j,k) .ne. 0) then
             ams(l) = ams(l) -sin(l*pw*xcrd(j))*sin(B*ycrd(k))*zbucket(j,k)*psi(j,k)&
                  &*2._wp*dA/(B*width)
             ! End if
          end do
       end do
    end do
    ams(inch) = ams(inch)-1._wp ! changed the delta's sign

    ! Calculate sum of |Am|^2
    sumAm = 0._wp
    do k = 1, M
       km = sqrt(E**2 - (k*PI/width)**2)
       ! print *, "sdf ", (km/kp), ams(k)!*abs(ams(k))**2
       sumAm = sumAm + (km/kp)*abs(ams(k))**2
       ! sumAm = sumAm + abs(ams(k))**2
    end do
    if (wp .eq. sp) then
       write(*,'(A,I2,A,I2,A,F10.8)') "Node-", nodeid, " | inch-", inch, " | Sum (km/kp)*|Am|^2: ", SumAm
    elseif (wp .eq. dp) then
       write(*,'(A,I2,A,I2,A,F18.16)') "Node-", nodeid, " | inch-", inch, " | Sum (km/kp)*|Am|^2: ", SumAm
    end if

    ! Add to S-matrix
    if (nodeid .ne. masternode) then
       write(*,'(A,I2,A,I2,A)') "Node-", nodeid, " | inch-", inch, " | is ready to send S-matrix elements"
       if (wp .eq. sp) then
          call MPI_send(ams,M,MPI_complex,masternode,send_tag,MPI_comm_world,error)
          call MPI_send(sumAm,1,MPI_real,masternode,send_tag,MPI_comm_world,error)
       elseif (wp .eq. dp) then
          call MPI_send(ams,M,MPI_double_complex,masternode,send_tag,MPI_comm_world,error)
          call MPI_send(sumAm,1,MPI_double_precision,masternode,send_tag,MPI_comm_world,error)
          write(*,'(A,I2,A,I2,A)') "Node-", nodeid, " | inch-", inch, " | sent S-matrix elements"
       end if
       srch = inch
       call MPI_send(srch,1,MPI_int,masternode,send_tag,MPI_comm_world,error)
       write(*,'(A,I2,A,I2,A)') "Node-", nodeid, " | inch-", inch, " | sent its channel number"
    else ! only masternode
       ! First masternode handle its own S-matrix part (ams)
       do j = 1, M
          km = sqrt(E**2-(j*PI/width)**2)
          smatrix(inch,j) = sqrt(km/kp)*ams(j)
       end do
       write(26,*) inch, sumAm, 1._wp-sumAm
       
       ! Then recieves from everyone else
       write(*,'(A,I2,A,I2,A)') "Node-", nodeid, " | inch-", inch, " | is ready to receive S-matrix elements"
       do from_node = 1, nnodes-1
          if (wp .eq. sp) then
             call MPI_recv(ams,M,MPI_complex,from_node,MPI_any_tag,MPI_comm_world,status,error)
             call MPI_recv(sumAm,1,MPI_real,from_node,MPI_any_tag,MPI_comm_world,status,error)
          elseif (wp .eq. dp) then
             call MPI_recv(ams,M,MPI_double_complex,from_node,MPI_any_tag,MPI_comm_world,status,error)
             call MPI_recv(sumAm,1,MPI_double_precision,from_node,MPI_any_tag,MPI_comm_world,status,error)
          end if
          call MPI_recv(srch,1,MPI_int,from_node,MPI_any_tag,MPI_comm_world,status,error)
          print *, "receiving from node: ", from_node
          ! print *, "error? ", E**2-(srch*PI/width)**2, E, srch, width
          srm = sqrt(E**2-(srch*PI/width)**2)
          do j = 1, M
             km  = sqrt(E**2-(j*PI/width)**2)
             smatrix(srch,j) = sqrt(km/srm)*ams(j)
          end do
          
          ! Write sum of all |Am|^2
          write(26,*) srch, sumAm, 1._wp-sumAm
          write(*,'(A,I2,A,I2,A,I2)') "Node-", nodeid, " | inch-", inch, " | received from Node-", from_node
       end do
    end if

    ! ! Write all Am's
    ! write(filename, '("tam", "-inch", I0, "-wl", I0, ".txt")') inch, nint(wl*10*1e9)
    ! print *, "filename: ", filename
    ! open(unit = 25, file = filename, form='formatted')
    ! write(25, FMT='(4A15)') "inch", "Am"
    ! do j = 1, M
    !    write (25,*) j, abs(ams(j))**2
    ! end do
    ! close(unit = 25)

    ! Check if |Am|^2 is below threshold
    ! if (sumAm .gt. hilim) then
    !    print *, "Sum |Am|^2       ", SumAm
    !    print *, "***"
    !    print *, "sum|Am|^2 is too large"
    !    print *, "***"
    !    if (chkbnds .eqv. .true.) then
    !       call exit(1)
    !    end if
    ! end if

    ! Print sum
    ! print *, "### Sum (km/kp)*|Am|^2:      ", SumAm, " ###"
    ! write(*,FMT='(A26,I0,A5,F8.4)') "The efficiency of channel ", inch, " is: ", 1._wp-SumAm
    ! write(*,FMT='(A12,F5.2)') "In percent: ", (1._wp-SumAm)*100._wp
  end subroutine makeSmatrix

  !_______________________________________________________________________________________________
  !
  ! Checks the unitarity of the S matrix
  !_______________________________________________________________________________________________
  subroutine checkUnitarity(smatrix,M,sumSSdag,offdiag)
    implicit none
    integer(wp)                              :: j, k
    integer(wp), intent(in)                  :: M
    real(wp), intent(out)                    :: sumSSdag, offdiag
    complex(wp), intent(in)                  :: smatrix(:,:)
    complex(wp), dimension(M,M)              :: SSdag

    SSdag = transpose(smatrix)
    SSdag = conjg(SSdag)
    SSdag = matmul(SSdag,smatrix) ! Should equal the identity matrix if S is unitary

    ! Off diagonal elements should be very small, sum of all below 1
    sumSSdag = 0._wp
    offdiag  = 0._wp
    do j = 1, M
       do k = 1, M
          if (j .ne. k) then
             offdiag  = offdiag + abs(smatrix(j,k))**2
             sumSSdag = sumSSdag + abs(SSdag(j,k))**2
          end if
       end do
    end do
    ! Scale according to size of S matrix and abs()**2
    offdiag  = offdiag/(M**2-M)
    sumSSdag = sumSSdag/(M**2-M)
    write(*,'(A)') "|"
    write(*, FMT='(A,F18.16)') "| Sum of all off-diagonal elements of S-matrix: ", offdiag
    write(*, FMT='(A,F18.16)') "| Sum of all off-diagonal elements of SSdag   : ", sumSSdag
    write(*,'(A)') "|"
    if (sumSSdag .ge. 1._wp) then
       print *, "Sum is too large! S is not unitary within the limits"
       if (chkbnds .eqv. .true.) then
          call exit(1)
       end if
    end if
  end subroutine checkUnitarity

  !_______________________________________________________________________________________________
  !
  ! Calculates the eigenvalues of the S-matrix after checking its correctness
  !_______________________________________________________________________________________________
  subroutine sEigenValues(smatrix,eigvals,M)
    implicit none
    integer(wp), intent(in)              :: M
    complex(wp), intent(inout)           :: smatrix(:,:)
    complex(wp), intent(out)             :: eigvals(:)
    complex(wp), dimension(M,M)          :: VL, VR, backupsmatrix
    complex(wp), allocatable             :: WORK(:)
    ! real(wp), intent(out)                :: abseig
    real(wp), dimension(2*M)             :: RWORK
    integer(wp)                          :: INFO, LWORK, j
    character(len=1)                     :: JOBVL, JOBVR
    character(len=90)                    :: filename                  

    JOBVL = 'N'
    JOBVR = 'N'
    backupsmatrix = smatrix ! Because we need to keep it for splitting later
    allocate(WORK(1))

    INFO = 42 ! For no good reason
    ! Call CGEEV/ZGEEV with workspace query
    if (wp .eq. sp) then
       call CGEEV(JOBVL,JOBVR,M,smatrix,M,eigvals,VL,M,VR,M,WORK,-1,RWORK,INFO)
    elseif (wp .eq. dp) then
       call ZGEEV(JOBVL,JOBVR,M,smatrix,M,eigvals,VL,M,VR,M,WORK,-1,RWORK,INFO)
    end if
    write(*, FMT='(A34,I0)') "Info after workspace query:       ", INFO
    LWORK = int(WORK(1))
    deallocate(WORK)
    allocate(WORK(LWORK))

    ! Call CGEEV/ZGEEV to compute eigenvalues
    if (wp .eq. sp) then
       call CGEEV(JOBVL,JOBVR,M,smatrix,M,eigvals,VL,M,VR,M,WORK,LWORK,RWORK,INFO)
    elseif (wp .eq. dp) then
       call ZGEEV(JOBVL,JOBVR,M,smatrix,M,eigvals,VL,M,VR,M,WORK,LWORK,RWORK,INFO)
    end if
    write(*, FMT='(A34,I0)') "Info after computing eigenvalues: ", INFO

    ! Print raw eigenvalues    
    print *, "The raw eigenvalues:"
    do j = 1, M
       print *, eigvals(j)
    end do

    ! Store eigenvalues
    write(filename, '("teigenvalues-wl", F0.1, ".txt")') wl*1e9
    open(unit = 25, file = filename, form='formatted')
    write(25, FMT='(2A15)') "real", "imag"
    do j = 1, M
       write(25,*) real(eigvals(j)), aimag(eigvals(j))
    end do
    close(unit = 25)

    ! Check that product of eigenvalues equal 1, a property of unitarity
    ! write(*, FMT='(A31,I0)') "The leading dimension of S is: ", M
    write(*, FMT='(A24,F19.16,A1,F19.16,A1)') "Product of eigenvalues: ", real(product(eigvals))&
         , " ", aimag(product(eigvals)), "i"
    write(*, FMT='(A42,F19.16)') "Absolute value of product of eigenvalues: ", abs(product(eigvals))

    if (abs(product(eigvals)) .gt. hilim .or. abs(product(eigvals)) .lt. lolim) then
       print *, ""
       print *, "* * * * * * * * * * * * * * * * * *"
       print *, "S is not unitary within the limits"
       print *, "* * * * * * * * * * * * * * * * * *"
       print *, ""
       if (chkbnds .eqv. .true.) then
          call exit(1)
       end if
    end if
    deallocate(WORK)

    ! Restore the old S-matrix
    smatrix = backupsmatrix
  end subroutine sEigenValues

  !_______________________________________________________________________________________________
  !
  ! Extract angles from eigenvalues and find nearest neighbor
  ! distance between numbers from  a list of real values
  !_______________________________________________________________________________________________
  subroutine nearestNeighborDist(evs,near,M)
    implicit none
    integer(wp), intent(in)                   :: M
    complex(wp), intent(in)                   :: evs(:)
    real(wp), dimension(M)                    :: angs
    real(wp), intent(out)                     :: near(:)
    integer(wp)                               :: j
    integer(sp)                               :: error, errorcode

    ! Extract eigenangles from eigenvalues on the form e^(i*angs)
    angs = 0._wp
    do j = 1, M
       angs(j) = atan2(aimag(evs(j)),real(evs(j)))
       if (angs(j) .lt. 0._wp) then
          angs(j) = angs(j)+2._wp*PI
       end if
    end do

    ! Scale eigenangles
    angs = angs*(M/(2._wp*PI))

    ! Sort eigenangles
    call quicksort(angs,int4(1),int4(M))
    
    print *, "The sorted eigenangles:"
    do j = 1, M
       print *, angs(j)
    end do

    ! Find nearest neighbors
    if (M .ge. 3) then
       near(1) = M-angs(M)+angs(1) ! Reinhold didn't include this one
       do j = 2, M
          near(j) = angs(j)-angs(j-1)
       end do
    elseif (M .lt. 3) then
       print *, "Dimension of S-matrix is too small (", M, ")"
       call MPI_abort(MPI_comm_world,errorcode,error)
       call exit(1)
    end if
  end subroutine nearestNeighborDist

  !_______________________________________________________________________________________________
  !
  ! Split S-matrix into even and odd parity submatrices
  !_______________________________________________________________________________________________
  subroutine splitSmatrix(smatrix,evensmatrix,oddsmatrix,evenM,oddM)
    implicit none
    complex(wp), intent(in)       :: smatrix(:,:)
    complex(wp), intent(inout)    :: evensmatrix(:,:), oddsmatrix(:,:)
    integer(wp), intent(in)       :: evenM, oddM
    integer(wp)                   :: j, k

    do j = 1, evenM
       do k = 1, evenM
          evensmatrix(k,j) = smatrix(k*2,j*2)
       end do
    end do

    do j = 1, oddM
       do k = 1, oddM
          oddsmatrix(k,j) = smatrix(2*(k-1)+1,2*(j-1)+1)
       end do
    end do

  end subroutine splitSmatrix

  !_______________________________________________________________________________________________
  !
  ! Quicksort algorithm.
  !_______________________________________________________________________________________________
  recursive subroutine quicksort(a, first, last)
    implicit none
    real(wp)     :: a(:), x, t
    integer(sp)  :: first, last
    integer(sp)  :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort

  !_______________________________________________________________________________________________
  !
  ! Get a seed for random number generation
  !_______________________________________________________________________________________________
  subroutine getRandomSeed(seed,seedlen,nnodes,nodeid)
    implicit none
    integer(sp), allocatable, intent(inout)    :: seed(:)
    integer(sp), intent(inout)                 :: seedlen
    integer(sp), intent(in)                    :: nnodes, nodeid
    integer(wp)                                :: to_node
    integer(sp)                                :: error
    integer                                    :: status(MPI_status_size)

    if (nodeid .eq. masternode) then
       call random_seed(size=seedlen)
       allocate(seed(seedlen))
       call random_seed(get=seed)
       
       do to_node = 1, nnodes-1
          ! Send the length of the seed to the other nodes
          call MPI_send(seedlen,1,MPI_int,to_node,send_tag,MPI_comm_world,error)
          ! Send the seed to the other nodes
          call MPI_send(seed,seedlen,MPI_int,to_node,send_tag,MPI_comm_world,error)
       end do
    else ! All other nodes
       ! Receive the length of the seed from masternode
       call MPI_recv(seedlen,1,MPI_int,masternode,MPI_any_tag,MPI_comm_world,status,error)
       allocate(seed(seedlen))
       ! Receive the seed from masternode
       call MPI_recv(seed,seedlen,MPI_int,masternode,MPI_any_tag,MPI_comm_world,status,error)
    end if
  end subroutine getRandomSeed

  !_______________________________________________________________________________________________
  !
  ! Integrates the absolute square of the wavefunction inside the bucket potential
  !_______________________________________________________________________________________________
  real(wp) function integratebucket(pot,wf)
    implicit none
    integer(wp)                               :: j, k
    real(wp), intent(in)                      :: pot(:,:)
    real(wp)                                  :: intwf
    complex(wp), dimension(nx,ny), intent(in) :: wf

    intwf = 0._wp
    do j = 1, nx
       do k = 1, ny
          if (pot(j,k) .ne. 0) then
             intwf = intwf + abs(wf(j,k))**2*dA
          end if
       end do
    end do
    integratebucket = intwf
  end function integratebucket

  !_______________________________________________________________________________________________
  !
  ! Green's function for light in a bucket with reflective boundary conditions
  !_______________________________________________________________________________________________
  complex(wp) function green(opch,B,A,x1,x2,y1,y2)
    implicit none
    integer(wp)                           :: m
    integer(wp), intent(in)               :: opch
    real(wp), intent(in)                  :: x1, x2, y1, y2
    real(wp), intent(in)                  :: B(:)
    real(wp), intent(in)                  :: A(:)
    complex(wp)                           :: opex1, opex2
    real(wp)                              :: yy, ssin, clex1, clex2

    if (        x1 .eq. 0._wp .or. x1 .eq. width &
         & .or. x2 .eq. 0._wp .or. x2 .eq. width &
         & .or. y1 .eq. 0._wp .or. y2 .eq. 0._wp ) then
       green = 0._wp
       RETURN
    end if

    green = 0._wp
    yy    = abs(y1-y2)
    do m = 1, opch
       ssin  = sin(pw*m*x1)*sin(pw*m*x2)
       opex1 = exp(i*B(m)*(y1+y2))
       opex2 = exp(i*B(m)*yy)
       green = green + (factor*i/B(m))*ssin*(opex1-opex2)
    end do
    if (opch .ne. 0) then
       do m = opch+1, opch+clch
          ssin  = sin(pw*m*x1)*sin(pw*m*x2)
          ! if (-A(m-opch)*(y1+y2) .lt. explim) then
          !    clex1 = 0._wp
          ! end if
          ! if (-A(m-opch)*yy .lt. explim) then
          !    clex2 = 0._wp
          ! end if
          ! if ((-A(m-opch)*(y1+y2) .ge. explim) .and. (-A(m-opch)*yy .ge. explim)) then
          !    clex1 = exp(-A(m-opch)*(y1+y2))
          !    clex2 = exp(-A(m-opch)*yy)
          ! elseif ((-A(m-opch)*(y1+y2) .lt. explim) .and. (-A(m-opch)*yy .lt. explim)) then
          !    exit
          ! end if

          clex1 = exp(-A(m-opch)*(y1+y2))
          clex2 = exp(-A(m-opch)*yy)
          !print *, (factor/A(m-opch)), ssin*(clex1-clex2), clex1, clex2
          green = green + (factor/A(m-opch))*ssin*(clex1-clex2)
       end do
    end if
  end function green

  !_______________________________________________________________________________________________
  !
  ! Get length of a file with header
  !_______________________________________________________________________________________________
  subroutine readFileLength(fileid,rows)
    integer(wp), intent(inout)    :: rows
    integer(sp), intent(in)       :: fileid
    integer(wp)                   :: io

    ! Count the number of lines in a file
    rows = 0
    rewind(fileid)
    do
       read(fileid,*,iostat=io)
       if (io .ne. 0) then
          ! print *, "eof reached, all good"
          exit
       end if
       rows = rows + 1
    end do
    rewind(fileid)

    ! Now rows is the total number of lines in the file
    ! not interested in the descriptors (1 line header), so:
    rows = rows-1


  end subroutine readFileLength

  !_______________________________________________________________________________________________
  !
  ! Create logbook entry
  !_______________________________________________________________________________________________
  subroutine logParameters(hours,minutes,seconds,tg,tl)
    implicit none
    integer(wp), intent(in)          :: hours, minutes, seconds
    real(dp), intent(in)             :: tg, tl
    real(dp)                         :: tsec

    tsec = hours*3600+minutes*60+seconds

    ! Create logbook
    open(unit = 25, file = "logbook.txt", form='formatted')
    write(25, FMT='(A62)') "______________________________________________________________"
    write(25, FMT='(A62)') ""
    write(25, FMT='(A15)') "Parameters list"
    write(25, FMT='(A62)') "______________________________________________________________"
    write(25, FMT='(A62)') ""
    write(25, FMT='(A26,I0,A1,I0)') "Bucket grid:              ", nx, "x", ny
    write(25, FMT='(A26,I0,A7)')    "Bucket width:             ", nint(width*1e6), " micron"
    write(25, FMT='(A26,I0,A3)')    "Bucket height:            ", nint(height*1e9), " nm"
    write(25, FMT='(A26,F7.2,A3)')  "Initial wavelength:       ", wli*1e9_wp, " nm"
    write(25, FMT='(A26,F7.2,A3)')  "Final wavelength:         ", wlf*1e9_wp, " nm"
    write(25, FMT='(A26,I0)')       "#Wavelengths used         ", wlsteps
    write(25, FMT='(A26,I0)')       "#Calls of Green:          ", N*N
    write(25, FMT='(A26,I0)')       "#Closed channels:         ", clch
    write(25, FMT='(A26,F5.2,A3)')  "RAM required by G:        ", N*N*16._wp/1e9, "GB"
    write(25, FMT='(A62)') "______________________________________________________________"
    write(25, FMT='(A62)') ""
    write(25, FMT='(A15)') "Potential structure details"
    write(25, FMT='(A62)') "______________________________________________________________"
    write(25, FMT='(A62)') ""
    select case (surface)
    case (1) ! Flat potential
       write(25, FMT='(A26,I0)')       "Flat potential            "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction:      ", real(rind), " + ", aimag(rind), "i"
    case (2) ! Sine
       write(25, FMT='(A26,I0)')       "Sine potential            "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,I0)')       "Shape parameter:          ", nint(amplt)
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction:      ", real(rind), " + ", aimag(rind), "i"
    case (3) ! Oblique potential
       write(25, FMT='(A26,I0)')       "Slant potential           "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,I0)')       "Shape parameter:          ", nint(amplt)
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction:      ", real(rind), " + ", aimag(rind), "i"
    case (4) ! House with 2 sub-houses
       write(25, FMT='(A26,I0)')       "House w/ 2 sub-houses     "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 1:     ", real(rind), " + ", aimag(rind), "i"
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 2:     ", real(rind2), " + ", aimag(rind2), "i"
       write(25, FMT='(A26,F4.1)')     "Shape parameter heightL:  ", heightL
       write(25, FMT='(A26,F4.1)')     "Shape parameter ratioL:   ", ratioL
       write(25, FMT='(A26,F4.1)')     "Shape parameter heightR:  ", heightR
       write(25, FMT='(A26,F4.1)')     "Shape parameter ratioR:   ", ratioR
    case (5) ! Circle top with half-pipe inside
       write(25, FMT='(A26,I0)')       "Slant potential           "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 1:     ", real(rind), " + ", aimag(rind), "i"
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 2:     ", real(rind2), " + ", aimag(rind2), "i"
    case (6) ! Dome with dome houses
       write(25, FMT='(A26,I0)')       "Slant potential           "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 1:     ", real(rind), " + ", aimag(rind), "i"
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 2:     ", real(rind2), " + ", aimag(rind2), "i"
       write(25, FMT='(A26,F4.1)')     "Shape parameter ecc1:     ", ecc1
       write(25, FMT='(A26,F4.1)')     "Shape parameter ecc1:     ", ecc2
    case (7) ! Sine with small houses
       write(25, FMT='(A26,I0)')       "Sine w/ 2 sub-houses     "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 1:     ", real(rind), " + ", aimag(rind), "i"
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 2:     ", real(rind2), " + ", aimag(rind2), "i"
       write(25, FMT='(A26,F4.1)')     "Shape parameter heightL:  ", heightL
       write(25, FMT='(A26,F4.1)')     "Shape parameter ratioL:   ", ratioL
       write(25, FMT='(A26,F4.1)')     "Shape parameter heightR:  ", heightR
       write(25, FMT='(A26,F4.1)')     "Shape parameter ratioR:   ", ratioR
    case (8) ! Adjustable saw potential
       write(25, FMT='(A26,I0)')       "Adjustable saw potential  "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
       write(25, FMT='(A26,F3.1,A3,F3.1,A1)') "Index of refraction 1:     ", real(rind), " + ", aimag(rind), "i"
       write(25, FMT='(A26,I0)')       "Shape parameter nHus:     ", nHus
       write(25, FMT='(A26,F3.1)')     "Shape parameter nRatio:   ", nRatio
       write(25, FMT='(A26,F3.1)')     "Shape parameter nHeight:  ", nHeight
    case (9) ! Chaotic dip potential
       write(25, FMT='(A26,I0)')       "Chaotic dip potential     "
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
    case (10) ! From file
       write(25, FMT='(A26,I0)')       "Random potential from file"
       write(25, FMT='(A26,I0)')       "Potential shape:          ", surface
    end select
    write(25, FMT='(A10)', advance="no") "Run time: "
    write(25, FMT='(I2,A6)') hours, " hours"
    write(25, FMT='(T11,I2,A8)') minutes, " minutes"
    write(25, FMT='(T11,I2,A8)') seconds, " seconds"
    write(25, FMT='(A1)')           " " !spacer
    if (tsec .ne. 0) then
       write(25, FMT='(A26,I0,A2,I0,A2)') "Time spent building T:    ", nint(tg), " (", nint(tg*100/tsec), "%)"
       write(25, FMT='(A26,I0,A2,I0,A2)') "Time spent linSolving:    ", nint(tl), " (", nint(tl*100/tsec), "%)"
    end if
    close(unit = 25)
  end subroutine logParameters

end module functions
