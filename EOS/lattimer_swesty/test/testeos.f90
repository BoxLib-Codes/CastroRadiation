program testeos

  use bl_types
  use network
  use eos_module

  implicit none

  double precision G, P, C, T, dpdr, dpde
  double precision R, e
  double precision Y(2)
  integer inspec
  double precision dedt, dedy

  inspec = 2

  call network_init()
  call eos_init()

  T = 0.466094505084d0
  R = 1.56379963476d13
  Y(1) = 1.d0
  Y(2) = 0.30232939317d0

  call eos_given_RTX(e, P, R, T, Y)
  print*,'e = ', e
  call eos_get_cv(dedt, R, T, Y)
  print *, 'dedT = ', dedt
  call eos_get_dedx(dedt, dedy, R, T, Y)
  print *, 'dedT = ', dedt, ' dedy = ', dedy

  print *, '------------------------------------'

  T = 0.150000798608d0
  R = 2.54474729742d13
  Y(1) = 1.d0
!  Y(2) = 0.302919059782d0

  call eos_given_RTX(e, P, R, T, Y)
  print*,'e = ', e
  call eos_get_cv(dedt, R, T, Y)
  print *, 'dedT = ', dedt
  call eos_get_dedx(dedt, dedy, R, T, Y)
  print *, 'dedT = ', dedt, ' dedy = ', dedy

  stop'done'

end program testeos

