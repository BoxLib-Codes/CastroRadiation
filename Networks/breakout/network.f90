! network module
!
! nspec -> number of species
! naux  -> number of auxiliary variables also used in the EOS
!
! aion -> atomic number
! zion -> proton number
!

module network
 
  use bl_types
 
  implicit none
 
  integer, parameter :: nspec = 1
  integer, parameter :: naux  = 2
 
  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)

  character (len=16), save ::  aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)
 
  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)
 
  logical, save :: network_initialized = .false.
 

contains

  subroutine network_init()

    spec_names(1) = "X"
    short_spec_names(1) = "X"

    aux_names(1) = "Ye"
    short_aux_names(1) = "Ye"

    aux_names(2) = "invmu"   ! 1/mu in P = rho*R*T/mu
    short_aux_names(2) = "invmu"

    aion(:) = 1.d0
    zion(:) = 1.d0
    ebin(:) = 0.d0

    network_initialized = .true.

  end subroutine network_init

end module network
