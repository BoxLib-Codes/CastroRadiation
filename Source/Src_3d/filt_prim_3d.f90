
subroutine ca_filt_prim(lo, hi, &
     Stmp, Stmp_l1, Stmp_l2, Stmp_l3, Stmp_h1, Stmp_h2, Stmp_h3, &
     Snew, Snew_l1, Snew_l2, Snew_l3, Snew_h1, Snew_h2, Snew_h3, &
     mask, mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3, &
     filt_T, S, domlo,domhi, &
     delta,xlo,problo,time,level)

  use network, only : naux, nspec
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
       UFA, UFS, UFX, small_temp, small_dens, nadv
  use filter_module

  implicit none

  integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3), level
  integer, intent(in) :: filt_T, S
  double precision, intent(in) :: delta(3), xlo(3), problo(3), time
  integer, intent(in) ::   Stmp_l1,Stmp_h1,Stmp_l2,Stmp_h2,Stmp_l3,Stmp_h3
  integer, intent(in) ::   Snew_l1,Snew_h1,Snew_l2,Snew_h2,Snew_l3,Snew_h3
  integer, intent(in) ::   mask_l1,mask_h1,mask_l2,mask_h2,mask_l3,mask_h3
  double precision :: Stmp(Stmp_l1:Stmp_h1,Stmp_l2:Stmp_h2,Stmp_l3:Stmp_h3,NVAR)
  double precision :: Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,Snew_l3:Snew_h3,NVAR)
  double precision :: mask(mask_l1:mask_h1,mask_l2:mask_h2,mask_l3:mask_h3)
  ! mask has three possible values: -1.d0, 0.d0, and 1.d0.
  ! -1.d0 appears only in cells that are covered by neither this level nor the finer level.
  !       It can only appear in ghost cells. 
  !  0.d0 appears only in cells that are covered by only this level, but not the finer level.
  !  1.d0 appears only in cells that are covered by the finer level.
  !       It can appear in either valid cells or ghost cells. 

  integer :: i,j,k
  double precision :: p, e, X(nspec+naux)

end subroutine ca_filt_prim

