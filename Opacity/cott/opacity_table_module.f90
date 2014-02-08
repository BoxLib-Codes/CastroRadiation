module opacity_table_module

contains

  subroutine prep_opacity(g, inu, er, der)

    use rad_params_module, only : get_ispec, nugroup, dnugroup, Hz2MeV, etafactor

    implicit none

    integer, intent(in) :: g
    integer, intent(out) :: inu
    double precision, intent(out) :: er, der
  
    inu = get_ispec(g) + 1
    er = nugroup(g) * Hz2MeV
    der = dnugroup(g) * Hz2MeV * etafactor

  end subroutine prep_opacity

!     ============================================
!
!     opacity for a given temperature, density, and ye computed
!     by an interpolation of the precalculated opacity table
!
!     this is a simplified routine with all interpolations linear
!
!     input:   er - neutrino energy (mev)
!              temp   - temperature (mev)
!              yein  - electron fraction
!              rho (in state vector) - density     (g cm^-3)
!     output:  ab  - absorptive opacity  (1/cm)
!              sc  - scattering opacity  (1/cm)
!              delta - scattering anisotropy factor
!              eta   - emissivity
!
  subroutine get_opacity_emissivity( &
       ab, sc, delta, eta,           &
       rho, ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

    implicit none

    integer, intent(in) :: inu 
    logical, intent(in) :: comp_ab, comp_sc, comp_eta
    double precision, intent(in) :: rho, ye, temp, er 
    double precision, intent(out) :: ab, sc, delta, eta

    delta = 0.d0

    call nulibtable_single_species_single_energy(ab, sc, eta, &
         rho, ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

  end subroutine get_opacity_emissivity

  subroutine get_opacities(kp, kr, rho, temp, nu, get_Planck_mean, get_Rosseland_mean)

    implicit none

    logical, intent(in) :: get_Planck_mean, get_Rosseland_mean
    double precision, intent(in) :: rho, temp, nu
    double precision, intent(out) :: kp, kr
    kp = 0.d0
    kr = 0.d0
  end subroutine get_opacities

end module opacity_table_module


subroutine init_opacity_table(iverb)
  implicit none
  integer iverb
  character(len=512), dimension(3) :: files
!
!  files(1) = "opac_e_rho100_temp50_ye40_ng35_20120425.h5"
!  files(2) = "opac_a_rho100_temp50_ye40_ng25_20120425.h5"
!  files(3) = "opac_x_rho100_temp50_ye40_ng30_20120425.h5"
!
  files(1) = "opac_e_rho50_temp50_ye35_ng35_20120425.h5"
  files(2) = "opac_a_rho50_temp50_ye35_ng25_20120425.h5"
  files(3) = "opac_x_rho50_temp50_ye35_ng30_20120425.h5"
!
  call nulibtable_reader(files)
end subroutine init_opacity_table

