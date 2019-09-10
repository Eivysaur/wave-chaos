module pbucket
  implicit none

  ! Precision, choose double precision (dp) or single precision (sp)
  integer, parameter                    :: dp = kind(0.d0)
  integer, parameter                    :: sp = kind(0.0) 
  integer, parameter                    :: wp = sp ! working precision
  ! integer, parameter                    :: long_double = selected_real_kind(16)

  ! Natural constants as.
  real(wp), parameter                   :: PI     = acos(-1._wp)
  real(wp), parameter                   :: planck = 6.62607004e-34_wp
  real(wp), parameter                   :: redpla = planck/(2._wp*PI)
  real(wp), parameter                   :: me     = 9.10938356e-9_wp
  real(wp), parameter                   :: c      = 3e8_wp
  real(wp), parameter                   :: eps    = 1.0e-15_wp
  complex(wp), parameter                :: i      = (0._wp, 1._wp)

  ! Potential
  ! Surface | 1=Flat, 2=Sine, 3=Oblique, 4=House w/ small houses, 5=Submerged sphere
  !         | 6=Dome w/ small domes, 7=Sine w/ small houses, 8=Peaks, 9=Asymm, 10=random from file
  !         | 11=Film+dome
  integer(wp), parameter                :: nx      = 250
  integer(wp), parameter                :: ny      = nx
  integer(wp), parameter                :: scale   = 5 ! Scale factor
  integer(wp)                           :: clch    = 9 ! Number of closed channels
  
  ! Constants specific to shapes
  ! Slant & sine
  real(wp), parameter                   :: amplt   = nint(0.7_wp*ny) ! Amplt = ny/2 for a half-slant
  real(wp), parameter                   :: bumps   = 3._wp
  ! Small houses
  real(wp), parameter                   :: heightL = 0.25_wp
  real(wp), parameter                   :: ratioL  = 0.25_wp
  real(wp), parameter                   :: heightR = 0.25_wp
  real(wp), parameter                   :: ratioR  = 0.25_wp
  ! The Big House
  real(wp), parameter                   :: heightB = 0.5_wp
  real(wp), parameter                   :: ratioB  = 0.5_wp
  ! Small domes
  real(wp), parameter                   :: ecc1    = scale*0.692e-6_wp
  ! real(wp), parameter                   :: ecc1    = 0.6e-6_wp
  real(wp), parameter                   :: ecc2    = 2.0_wp
  real(wp), parameter                   :: yoff    = scale*400e-9_wp
  integer(wp), parameter                :: dupes   = 1
  ! Adjustable saw potential
  integer(wp), parameter                :: nHus    = 3
  real(wp), parameter, dimension(nhus)  :: nRatio  = 0.5_wp
  real(wp), parameter, dimension(nhus)  :: nHeight = 0.5_wp
  ! Random potential
  real(wp), parameter                   :: simheight = 0.3_wp ! 
  ! Select potential type
  integer(wp), parameter                :: surface = 11
  logical                               :: rndpotential = .false.
  
  ! System parameters
  real(wp), parameter                   :: wli     = 477e-9_wp
  real(wp), parameter                   :: wlf     = 497e-9_wp
  integer(wp), parameter                :: wlsteps = 41
  integer(wp), parameter                :: specch  = 1
  complex(wp), parameter                :: rind    = 1.03_wp!+i*0.05_wp ! Dome
  complex(wp), parameter                :: rind2   = 2.0_wp!+i*0.0054113_wp ! Film
  real(wp), parameter                   :: width   = scale*1000.0e-9_wp ! Width in meter
  real(wp)                              :: height  = scale*1100.0e-9_wp  ! Height in meter
  integer(wp), parameter                :: N       = nx*ny
  real(wp)                              :: dx
  real(wp)                              :: dy
  real(wp)                              :: dA
  real(wp), parameter                   :: hilim   = 1.08_wp ! Highest allowed sum|Am|^2
  real(wp), parameter                   :: lolim   = 0.92_wp ! Lowest  allowed sum|Am|^2
  logical, parameter                    :: chkbnds = .false. ! Terminate program if value>bounds
  real(wp)                              :: wl      = 1000e-9_wp ! Global wavelength

  ! Switch
  integer(wp), parameter                :: setinch = 0 ! For chaosbucket
  logical, parameter                    :: allch = .true. ! False; runs program only in channel setinch

  ! MPI
  integer(wp), parameter                :: masternode = 0
  integer(wp), parameter                :: send_tag = 2001
  integer(wp), parameter                :: return_tag = 2002

  ! Green's function
  real(wp), parameter                   :: pw     = (PI/width)
  real(wp), parameter                   :: factor = 1._wp/width
  real(wp)                              :: explim

end module pbucket
