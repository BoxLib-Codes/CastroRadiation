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
  integer, parameter :: naux  = 1
 
  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)

  character (len=16), save ::  aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)
 
  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)
 
  logical, save :: network_initialized = .false.
 

contains

  subroutine network_init()

!   integer :: ih, ihe, iox, ife

    ! integer keys -- for convenience.  In all other places, we will find
    ! these by querying based on species name using network_species_index
!   ih  = 1
!   ihe = 2
!   iox = 3
!   ife = 4
 
!   spec_names(ih)  = "hydrogen"
!   spec_names(ihe)  = "helium"
!   spec_names(iox) = "oxygen"
!   spec_names(ife) = "iron"

!   short_spec_names(ih)  = "H"
!   short_spec_names(ihe)  = "He"
!   short_spec_names(iox) = "O"
!   short_spec_names(ife) = "Fe"

!   aion(ih) = 1.d0
!   aion(ihe) = 4.d0
!   aion(iox) = 16.d0
!   aion(ife) = 56.d0

!   zion(ih) = 1.d0
!   zion(ihe) = 2.d0
!   zion(iox) = 8.d0
!   zion(ife) = 26.d0

    spec_names(1) = "X"
    short_spec_names(1) = "X"

    aux_names(1) = "Ye"
    short_aux_names(1) = "Ye"

    ebin(:) = 0.d0

    network_initialized = .true.

  end subroutine network_init

end module network
