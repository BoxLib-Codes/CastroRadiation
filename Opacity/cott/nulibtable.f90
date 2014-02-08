module nulibtable

  implicit none

  integer, save :: nulibtable_number_species
  
  type neutrino
     integer :: number_groups
     real*8, allocatable :: energies(:)
     real*8, allocatable :: inv_energies(:)
     real*8, allocatable :: ewidths(:)
     real*8, allocatable :: ebottom(:)
     real*8, allocatable :: etop(:)
     real*8, allocatable :: emissivities(:,:,:,:)
     real*8, allocatable :: absopacity(:,:,:,:)
     real*8, allocatable :: scatopacity(:,:,:,:)
  end type neutrino

  type(neutrino), target, save :: electron_neut, electron_antineut, muon_neut

  real*8, allocatable,save :: nulibtable_logrho(:)
  real*8, allocatable,save :: nulibtable_logtemp(:)
  real*8, allocatable,save :: nulibtable_ye(:)
  
  integer,save :: nulibtable_nrho
  integer,save :: nulibtable_ntemp
  integer,save :: nulibtable_nye

  real*8,save :: nulibtable_logrho_min
  real*8,save :: nulibtable_logrho_max

  real*8,save :: nulibtable_logtemp_min
  real*8,save :: nulibtable_logtemp_max

  real*8,save :: nulibtable_ye_min
  real*8,save :: nulibtable_ye_max

contains

  subroutine init_neutrino(neut, ngrp)
    type(neutrino) :: neut
    integer :: ngrp

    neut%number_groups = ngrp

    allocate(neut%energies(neut%number_groups))
    allocate(neut%inv_energies(neut%number_groups))
    allocate(neut%ewidths(neut%number_groups))
    allocate(neut%ebottom(neut%number_groups))
    allocate(neut%etop(neut%number_groups))

    allocate(neut%emissivities(nulibtable_nrho,nulibtable_ntemp, &
         nulibtable_nye,neut%number_groups))
    allocate(neut%absopacity(nulibtable_nrho,nulibtable_ntemp, &
         nulibtable_nye,neut%number_groups))
    allocate(neut%scatopacity(nulibtable_nrho,nulibtable_ntemp, &
         nulibtable_nye,neut%number_groups))
  end subroutine init_neutrino

end module nulibtable


subroutine nulibtable_single_species_single_energy(ab, sc, eta, &
     rho, ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

  use nulibtable
  implicit none

  real*8, intent(out) :: ab, sc, eta
  real*8, intent(in) :: rho, ye, temp, er
  integer, intent(in) :: inu
  logical, intent(in) :: comp_ab, comp_sc, comp_eta

  real*8 :: lgrho, lgtemp
  type(neutrino), pointer :: neutp
  integer :: ier, g
  
  if (inu .eq. 1) then
     neutp => electron_neut
  else if (inu .eq. 2) then
     neutp => electron_antineut
  else
     neutp => muon_neut
  end if

  lgrho = log10(rho)
  lgtemp = log10(temp)

  if (lgrho.lt.nulibtable_logrho_min) stop "density below nulib table minimum rho"
  if (lgrho.gt.nulibtable_logrho_max) stop "density above nulib table maximum rho"
  if (lgtemp.lt.nulibtable_logtemp_min) stop "temperature below nulib table minimum temp"
  if (lgtemp.gt.nulibtable_logtemp_max) stop "temperature above nulib table maximum temp"
  if (ye.lt.nulibtable_ye_min) stop "ye below nulib table minimum ye"
  if (ye.gt.nulibtable_ye_max) stop "ye above nulib table maximum ye"

  if (er .lt. neutp%ebottom(1)) then
     ier = 1
  else
     do ier = 1, neutp%number_groups
        if (er.lt.neutp%etop(ier)) exit
     end do
  end if

  if (comp_eta) then
     call intp3d_many_mod(lgrho,lgtemp,ye,eta, &
          neutp%emissivities(:,:,:,ier),nulibtable_nrho, &
          nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
          nulibtable_logtemp,nulibtable_ye)
     eta = 10.d0**eta / neutp%ewidths(ier)
  end if

  if (comp_ab) then
     call intp3d_many_mod(lgrho,lgtemp,ye,ab, &
          neutp%absopacity(:,:,:,ier),nulibtable_nrho, &
          nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
          nulibtable_logtemp,nulibtable_ye)
     ab = 10.d0**ab
  end if

  if (comp_sc) then
     call intp3d_many_mod(lgrho,lgtemp,ye,sc, &
          neutp%scatopacity(:,:,:,ier),nulibtable_nrho, &
          nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
          nulibtable_logtemp,nulibtable_ye)
     sc = 10.d0**sc
  end if

  nullify(neutp)

end subroutine nulibtable_single_species_single_energy
