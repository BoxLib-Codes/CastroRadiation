program testeos

  use bl_types
  use network
  use eos_module

  implicit none

  double precision G, P, C, T, dpdr, dpde
  double precision R, e
  double precision Y(2)
  integer inspec

  inspec = 2

  call network_init()
  call eos_init()

!  T = 5.887076d-1      ! 1.24218878652991360D18

  T = 5.8870763d-1     ! 1.24218883525691802D18
!  T = 5.88707635d-1    ! 1.24220609108067200D18

!  T = 5.8870764d-1     ! 1.24220609920989005D18
!  T = 5.8870765d-1     ! 1.24220611546832819D18
!  T = 5.887077d-1      ! 1.24220619676051123D18

  R = 684695765.96184516d0
  Y(1) = 1.d0
  Y(2) = 0.46117011670248287d0

  call eos_given_RTX(e, P, R, T, Y)
  print*,'e = ', e

  print *, '------------------------------------'

  T = 5.88707635d-1    ! 1.24220609108067200D18
  call eos_given_RTX(e, P, R, T, Y)
  print*,'e = ', e

  stop'done'

end program testeos

