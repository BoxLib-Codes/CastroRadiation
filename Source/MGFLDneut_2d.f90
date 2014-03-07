
subroutine ca_accel_acoe_neut(   &
     eta1,eta1_l1,eta1_l2,eta1_h1,eta1_h2, &
     theT,theT_l1,theT_l2,theT_h1,theT_h2, &
     theY,theY_l1,theY_l2,theY_h1,theY_h2, &
     spc , spc_l1, spc_l2, spc_h1, spc_h2, &
     kap , kap_l1, kap_l2, kap_h1, kap_h2, &
     aco , aco_l1, aco_l2, aco_h1, aco_h2, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: eta1_l1,eta1_h1,eta1_l2,eta1_h2
  integer, intent(in) :: theT_l1,theT_h1,theT_l2,theT_h2
  integer, intent(in) :: theY_l1,theY_h1,theY_l2,theY_h2
  integer, intent(in) ::  spc_l1, spc_h1, spc_l2, spc_h2
  integer, intent(in) ::  kap_l1, kap_h1, kap_l2, kap_h2
  integer, intent(in) ::  aco_l1, aco_h1, aco_l2, aco_h2
  double precision, intent(in ) :: eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2)
  double precision, intent(in ) :: theT(theT_l1:theT_h1,theT_l2:theT_h2)
  double precision, intent(in ) :: theY(theY_l1:theY_h1,theY_l2:theY_h2)
  double precision, intent(in ) :: spc ( spc_l1: spc_h1, spc_l2: spc_h2,0:ngroups-1)
  double precision, intent(in ) :: kap ( kap_l1: kap_h1, kap_l2: kap_h2,0:ngroups-1)
  double precision, intent(out) :: aco ( aco_l1: aco_h1, aco_l2: aco_h2)
  double precision, intent(in) :: dt, tau

  integer :: i, j, g
  double precision :: kbar, kybar, H1, Theta, dt1, foo

  dt1 = (1.d0+tau)/dt

  do j = aco_l2, aco_h2
  do i = aco_l1, aco_h1
     kbar = 0.d0
     kybar = 0.d0
     do g=0, ngroups-1
        foo = spc(i,j,g) * kap(i,j,g)
        kbar = kbar + foo
        kybar = kybar + foo * erg2rhoYe(g)
     end do

     H1 = eta1(i,j)
     Theta = theY(i,j) - theT(i,j)

     aco(i,j) = (H1*kbar - Theta*kybar) * clight + dt1
  end do
  end do

end subroutine ca_accel_acoe_neut


subroutine ca_accel_rhs_neut(  &
     Ern , Ern_l1, Ern_l2, Ern_h1, Ern_h2,  &
     Erl , Erl_l1, Erl_l2, Erl_h1, Erl_h2,  &
     kap , kap_l1, kap_l2, kap_h1, kap_h2,  &
     etaT,etaT_l1,etaT_l2,etaT_h1,etaT_h2,  &
     etaY,etaY_l1,etaY_l2,etaY_h1,etaY_h2,  &
     theT,theT_l1,theT_l2,theT_h1,theT_h2,  &
     theY,theY_l1,theY_l2,theY_h1,theY_h2,  &
     rhs , rhs_l1, rhs_l2, rhs_h1, rhs_h2,  &
     dt)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: Ern_l1, Ern_h1, Ern_l2, Ern_h2
  integer, intent(in) :: Erl_l1, Erl_h1, Erl_l2, Erl_h2
  integer, intent(in) :: kap_l1, kap_h1, kap_l2, kap_h2
  integer, intent(in) ::etaT_l1,etaT_h1,etaT_l2,etaT_h2
  integer, intent(in) ::etaY_l1,etaY_h1,etaY_l2,etaY_h2
  integer, intent(in) ::theT_l1,theT_h1,theT_l2,theT_h2
  integer, intent(in) ::theY_l1,theY_h1,theY_l2,theY_h2
  integer, intent(in) :: rhs_l1, rhs_h1, rhs_l2, rhs_h2
  double precision, intent(in ) ::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2,0:ngroups-1)
  double precision, intent(in ) ::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2,0:ngroups-1)
  double precision, intent(in ) :: kap( kap_l1: kap_h1, kap_l2: kap_h2,0:ngroups-1)
  double precision, intent(in ) ::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2)
  double precision, intent(in ) ::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2)
  double precision, intent(in ) ::theT(theT_l1:theT_h1,theT_l2:theT_h2)
  double precision, intent(in ) ::theY(theY_l1:theY_h1,theY_l2:theY_h2)
  double precision, intent(out) :: rhs( rhs_l1: rhs_h1, rhs_l2: rhs_h2)
  double precision, intent(in) :: dt

  integer :: i, j, g
  double precision :: rt, ry, H, Theta, foo

  do j = rhs_l2, rhs_h2
  do i = rhs_l1, rhs_h1
     rt = 0.d0
     ry = 0.d0
     do g=0,ngroups-1
        foo = kap(i,j,g)*(Ern(i,j,g)-Erl(i,j,g))
        rt = rt + foo
        ry = ry + foo * erg2rhoYe(g)
     end do

     H = etaT(i,j) - etaY(i,j)
     Theta = theY(i,j) - theT(i,j)

     rhs(i,j) = clight*(H*rt + Theta*ry)
  end do
  end do

end subroutine ca_accel_rhs_neut


subroutine ca_accel_spec_neut( lo, hi, &
     Ern , Ern_l1, Ern_l2, Ern_h1, Ern_h2, &
     Erl , Erl_l1, Erl_l2, Erl_h1, Erl_h2, &
     kap , kap_l1, kap_l2, kap_h1, kap_h2, &
     etaT,etaT_l1,etaT_l2,etaT_h1,etaT_h2, &
     etaY,etaY_l1,etaY_l2,etaY_h1,etaY_h2, &
     theT,theT_l1,theT_l2,theT_h1,theT_h2, &
     theY,theY_l1,theY_l2,theY_h1,theY_h2, &
     mugT,mugT_l1,mugT_l2,mugT_h1,mugT_h2, &
     mugY,mugY_l1,mugY_l2,mugY_h1,mugY_h2, &
     spec,spec_l1,spec_l2,spec_h1,spec_h2, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2
  integer,intent(in)::etaY_l1,etaY_h1,etaY_l2,etaY_h2
  integer,intent(in)::theT_l1,theT_h1,theT_l2,theT_h2
  integer,intent(in)::theY_l1,theY_h1,theY_l2,theY_h2
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2
  integer,intent(in)::mugY_l1,mugY_h1,mugY_l2,mugY_h2
  integer,intent(in)::spec_l1,spec_h1,spec_l2,spec_h2
  double precision,intent(in )::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2,0:ngroups-1)
  double precision,intent(in )::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2,0:ngroups-1)
  double precision,intent(in )::kap ( kap_l1: kap_h1, kap_l2: kap_h2,0:ngroups-1)
  double precision,intent(in )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2)
  double precision,intent(in )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2)
  double precision,intent(in )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,0:ngroups-1)
  double precision,intent(in )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,0:ngroups-1)
  double precision,intent(out)::spec(spec_l1:spec_h1,spec_l2:spec_h2,0:ngroups-1)
  double precision,intent(in) :: dt, tau

  integer :: i, j, g
  double precision :: cdt1, rt, ry, p, q, r, s, foo, sumeps
  double precision,dimension(0:ngroups-1)::Hg, Tg, epsilon, kapt, kk

  cdt1 = 1.d0/(clight*dt)

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     rt = 0.d0
     ry = 0.d0
     do g=0,ngroups-1
        foo = kap(i,j,g)*(Ern(i,j,g)-Erl(i,j,g))
        rt = rt + foo
        ry = ry + foo * erg2rhoYe(g)
     end do

     Hg = mugT(i,j,:)*etaT(i,j) - mugY(i,j,:)*etaY(i,j)
     Tg = -mugT(i,j,:)*theT(i,j) + mugY(i,j,:)*theY(i,j)

     kapt = kap(i,j,:) + (1.d0+tau)*cdt1
     kk = kap(i,j,:) / kapt

     p = 1.d0 - sum(Hg*kk)
     q = 1.d0 - sum(Tg*erg2rhoYe*kk)
     r = sum(Hg*erg2rhoYe*kk)
     s = sum(Tg*kk)

     epsilon = ((r*Tg + q*Hg) * rt + (s*Hg + p*Tg) * ry) / kapt

     sumeps = sum(epsilon)
     if (sumeps .eq. 0.d0) then
        sumeps = 1.d-50
     end if
     
     spec(i,j,:) = epsilon / sumeps
  end do
  end do
  
end subroutine ca_accel_spec_neut


subroutine ca_check_conv_neut( lo, hi, &
     ren,ren_l1,ren_l2,ren_h1,ren_h2, &
     res,res_l1,res_l2,res_h1,res_h2, &
     re2,re2_l1,re2_l2,re2_h1,re2_h2, &
     ern,ern_l1,ern_l2,ern_h1,ern_h2, &
     Tmn,Tmn_l1,Tmn_l2,Tmn_h1,Tmn_h2, &
     Tms,Tms_l1,Tms_l2,Tms_h1,Tms_h2, &
     rYn,rYn_l1,rYn_l2,rYn_h1,rYn_h2, &
     rYs,rYs_l1,rYs_l2,rYs_h1,rYs_h2, &
     rY2,rY2_l1,rY2_l2,rY2_h1,rY2_h2, &
     rho,rho_l1,rho_l2,rho_h1,rho_h2, &
     kap,kap_l1,kap_l2,kap_h1,kap_h2, &
     jg , jg_l1, jg_l2, jg_h1, jg_h2, &
     deT,deT_l1,deT_l2,deT_h1,deT_h2, &
     deY,deY_l1,deY_l2,deY_h1,deY_h2, &
     rel_re, abs_re, &
     rel_FT, abs_FT, rel_T, abs_T, &
     rel_FY, abs_FY, rel_Y, abs_Y, &
     dt)
  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer,intent(in)::lo(2),hi(2)
  integer,intent(in)::ren_l1, ren_h1, ren_l2, ren_h2
  integer,intent(in)::res_l1, res_h1, res_l2, res_h2
  integer,intent(in)::re2_l1, re2_h1, re2_l2, re2_h2
  integer,intent(in)::ern_l1, ern_h1, ern_l2, ern_h2
  integer,intent(in)::Tmn_l1, Tmn_h1, Tmn_l2, Tmn_h2
  integer,intent(in)::Tms_l1, Tms_h1, Tms_l2, Tms_h2
  integer,intent(in)::rYn_l1, rYn_h1, rYn_l2, rYn_h2
  integer,intent(in)::rYs_l1, rYs_h1, rYs_l2, rYs_h2
  integer,intent(in)::rY2_l1, rY2_h1, rY2_l2, rY2_h2
  integer,intent(in)::rho_l1, rho_h1, rho_l2, rho_h2
  integer,intent(in)::kap_l1, kap_h1, kap_l2, kap_h2
  integer,intent(in):: jg_l1,  jg_h1,  jg_l2,  jg_h2
  integer,intent(in)::deT_l1, deT_h1, deT_l2, deT_h2
  integer,intent(in)::deY_l1, deY_h1, deY_l2, deY_h2
  double precision,intent(in   )::ren(ren_l1:ren_h1,ren_l2:ren_h2)
  double precision,intent(in   )::res(res_l1:res_h1,res_l2:res_h2)
  double precision,intent(in   )::re2(re2_l1:re2_h1,re2_l2:re2_h2)
  double precision,intent(in   )::ern(ern_l1:ern_h1,ern_l2:ern_h2,0:ngroups-1)
  double precision,intent(in   )::Tmn(Tmn_l1:Tmn_h1,Tmn_l2:Tmn_h2)
  double precision,intent(in   )::Tms(Tms_l1:Tms_h1,Tms_l2:Tms_h2)
  double precision,intent(in   )::rYn(rYn_l1:rYn_h1,rYn_l2:rYn_h2)
  double precision,intent(in   )::rYs(rYs_l1:rYs_h1,rYs_l2:rYs_h2)
  double precision,intent(in   )::rY2(rY2_l1:rY2_h1,rY2_l2:rY2_h2)
  double precision,intent(in   )::rho(rho_l1:rho_h1,rho_l2:rho_h2)
  double precision,intent(in   )::kap(kap_l1:kap_h1,kap_l2:kap_h2,0:ngroups-1)
  double precision,intent(in   ):: jg( jg_l1: jg_h1, jg_l2: jg_h2,0:ngroups-1)
  double precision,intent(in   )::deT(deT_l1:deT_h1,deT_l2:deT_h2)
  double precision,intent(in   )::deY(deY_l1:deY_h1,deY_l2:deY_h2)
  double precision,intent(inout)::rel_re, abs_re
  double precision,intent(inout)::rel_FT, abs_FT,rel_T, abs_T
  double precision,intent(inout)::rel_FY, abs_FY,rel_Y, abs_Y
  double precision,intent(in) :: dt

  integer :: i, j
  double precision :: chg, relchg, FT, FY, cdt, FTdenom, FYdenom, dTe, dYe

  cdt = clight*dt

  do j=lo(2),hi(2)
  do i=lo(1),hi(1)
     chg = abs(ren(i,j) - res(i,j))
     relchg = abs(chg/(ren(i,j)+1.d-50))
     rel_re = max(rel_re,relchg)
     abs_re = max(abs_re,chg)

     chg = abs(Tmn(i,j) - Tms(i,j))
     relchg = abs(chg/(Tmn(i,j)+1.d-50))
     rel_T = max(rel_T,relchg)
     abs_T = max(abs_T,chg)

     chg = abs(rYn(i,j) - rYs(i,j))
     relchg = abs(chg/(rYn(i,j)+1.d-50))
     rel_Y = max(rel_Y,relchg)
     abs_Y = max(abs_Y,chg/rho(i,j))

     FT = abs((ren(i,j)-re2(i,j)) - cdt*sum(kap(i,j,:)*Ern(i,j,:)-jg(i,j,:)))
     FY = abs((rYn(i,j)-rY2(i,j)) - cdt*sum(erg2rhoYe*(kap(i,j,:)*Ern(i,j,:)-jg(i,j,:))))

     dTe = Tmn(i,j)
     dYe = rYn(i,j)/rho(i,j)
     FTdenom = rho(i,j)*(abs(deT(i,j)*dTe)+abs(deY(i,j)*dYe))
     FYdenom = rYn(i,j)

     rel_FT = max(rel_FT, FT/(FTdenom+1.d-50))
     abs_FT = max(abs_FT, FT)

     rel_FY = max(rel_FY, FY/(FYdenom+1.d-50))
     abs_FY = max(abs_FY, FY)
  end do
  end do

end subroutine ca_check_conv_neut


subroutine ca_check_conv_er_neut( lo, hi,  &
     Ern , Ern_l1, Ern_l2, Ern_h1, Ern_h2, &
     Erl , Erl_l1, Erl_l2, Erl_h1, Erl_h2, &
     kap , kap_l1, kap_l2, kap_h1, kap_h2, &
     etTz,etTz_l1,etTz_l2,etTz_h1,etTz_h2, &
     etYz,etYz_l1,etYz_l2,etYz_h1,etYz_h2, &
     thTz,thTz_l1,thTz_l2,thTz_h1,thTz_h2, &
     thYz,thYz_l1,thYz_l2,thYz_h1,thYz_h2, &
     temp,temp_l1,temp_l2,temp_h1,temp_h2, &
     Ye  ,  Ye_l1,  Ye_l2,  Ye_h1,  Ye_h2, &
     rela, abso, errr, dt)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer,intent(in):: lo(2), hi(2)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2
  integer,intent(in)::temp_l1,temp_h1,temp_l2,temp_h2
  integer,intent(in)::  Ye_l1,  Ye_h1,  Ye_l2,  Ye_h2
  integer,intent(in)::etTz_l1,etTz_h1,etTz_l2,etTz_h2
  integer,intent(in)::etYz_l1,etYz_h1,etYz_l2,etYz_h2
  integer,intent(in)::thTz_l1,thTz_h1,thTz_l2,thTz_h2
  integer,intent(in)::thYz_l1,thYz_h1,thYz_l2,thYz_h2
  double precision,intent(in):: Ern( Ern_l1: Ern_h1, Ern_l2: Ern_h2,0:ngroups-1)
  double precision,intent(in):: Erl( Erl_l1: Erl_h1, Erl_l2: Erl_h2,0:ngroups-1)
  double precision,intent(in):: kap( kap_l1: kap_h1, kap_l2: kap_h2,0:ngroups-1)
  double precision,intent(in)::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2)
  double precision,intent(in)::etYz(etYz_l1:etYz_h1,etYz_l2:etYz_h2)
  double precision,intent(in)::thTz(thTz_l1:thTz_h1,thTz_l2:thTz_h2)
  double precision,intent(in)::thYz(thYz_l1:thYz_h1,thYz_l2:thYz_h2)
  double precision,intent(in)::temp(temp_l1:temp_h1,temp_l2:temp_h2)
  double precision,intent(in)::Ye  (  Ye_l1:  Ye_h1,  Ye_l2:  Ye_h2)
  double precision, intent(inout) :: rela, abso, errr
  double precision, intent(in) :: dt

  integer :: i, j, g
  double precision :: chg, tot, cdt, der, kdeT, kdeY, err_T, err_Y, err

  cdt = clight * dt
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     chg = 0.d0
     tot = 0.d0
     kdeT = 0.d0
     kdeY = 0.d0
     do g=0,ngroups-1
        der = Ern(i,j,g)-Erl(i,j,g)
        chg = chg + abs(der)
        tot = tot + abs(Ern(i,j,g))
        kdeT = kdeT + kap(i,j,g)*der
        kdeY = kdeY + erg2rhoYe(g)*kap(i,j,g)*der
     end do
     abso = max(abso, chg)
     rela = max(rela, chg / (tot + 1.d-50))

     err_T =  etTz(i,j)*kdeT - thTz(i,j)*kdeY 
     err_Y = -etYz(i,j)*kdeT + thYz(i,j)*kdeY 
     err = max(abs(err_T/(temp(i,j)+1.d-50)), abs(err_Y/(Ye(i,j)+1.d-50)))
     errr = max(errr, err)
  end do
  end do

end subroutine ca_check_conv_er_neut


subroutine ca_compute_coupty(  &
     cpt,cpt_l1,cpt_l2,cpt_h1,cpt_h2, &
     cpy,cpy_l1,cpy_l2,cpy_h1,cpy_h2, &
     kpp,kpp_l1,kpp_l2,kpp_h1,kpp_h2, &
     eg , eg_l1, eg_l2, eg_h1, eg_h2, &
     jg , jg_l1, jg_l2, jg_h1, jg_h2)
  
  use rad_params_module, only : ngroups, erg2rhoYe

  implicit none

  integer, intent(in) :: cpt_l1, cpt_h1, cpt_l2, cpt_h2
  integer, intent(in) :: cpy_l1, cpy_h1, cpy_l2, cpy_h2
  integer, intent(in) :: kpp_l1, kpp_h1, kpp_l2, kpp_h2
  integer, intent(in) ::  eg_l1,  eg_h1,  eg_l2,  eg_h2
  integer, intent(in) ::  jg_l1,  jg_h1,  jg_l2,  jg_h2
  double precision, intent(out):: cpt(cpt_l1:cpt_h1,cpt_l2:cpt_h2)
  double precision, intent(out):: cpy(cpy_l1:cpy_h1,cpy_l2:cpy_h2)
  double precision, intent(in) :: kpp(kpp_l1:kpp_h1,kpp_l2:kpp_h2,0:ngroups-1)
  double precision, intent(in) ::  eg( eg_l1: eg_h1, eg_l2: eg_h2,0:ngroups-1)
  double precision, intent(in) ::  jg( jg_l1: jg_h1, jg_l2: jg_h2,0:ngroups-1)

  integer :: i, j, g
  double precision :: foo

  cpt = 0.d0
  cpy = 0.d0

  do g=0, ngroups-1
     do j=cpt_l2, cpt_h2
     do i=cpt_l1, cpt_h1
        foo = kpp(i,j,g) * eg(i,j,g) - jg(i,j,g)
        cpt(i,j) = cpt(i,j) + foo
        cpy(i,j) = cpy(i,j) + erg2rhoYe(g) * foo
     end do
     end do
  end do

end subroutine ca_compute_coupty


subroutine ca_compute_dedx(  &
     S   ,   S_l1,   S_l2,   S_h1,   S_h2, &
     T   ,   T_l1,   T_l2,   T_h1,   T_h2, &
     Ye  ,  Ye_l1,  Ye_l2,  Ye_h1,  Ye_h2, &
     Ts  ,  Ts_l1,  Ts_l2,  Ts_h1,  Ts_h2, &
     Yes , Yes_l1, Yes_l2, Yes_h1, Yes_h2, &
     dedT,dedT_l1,dedT_l2,dedT_h1,dedT_h2, &
     dedY,dedY_l1,dedY_l2,dedY_h1,dedY_h2, &
     validStar)

  use eos_module, only : eos, eos_t, eos_input_rt, Tmin=>table_Tmin, Ymin=>table_Yemin
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS

  implicit none

  integer, intent(in) ::    S_l1,    S_h1,    S_l2,    S_h2 
  integer, intent(in) ::    T_l1,    T_h1,    T_l2,    T_h2
  integer, intent(in) ::   Ye_l1,   Ye_h1,   Ye_l2,   Ye_h2
  integer, intent(in) ::   Ts_l1,   Ts_h1,   Ts_l2,   Ts_h2
  integer, intent(in) ::  Yes_l1,  Yes_h1,  Yes_l2,  Yes_h2
  integer, intent(in) :: dedT_l1, dedT_h1, dedT_l2, dedT_h2
  integer, intent(in) :: dedY_l1, dedY_h1, dedY_l2, dedY_h2
  double precision, intent(in ) :: S   (   S_l1:   S_h1,   S_l2:   S_h2,NVAR)
  double precision, intent(in ) :: T   (   T_l1:   T_h1,   T_l2:   T_h2)
  double precision, intent(in ) :: Ye  (  Ye_l1:  Ye_h1,  Ye_l2:  Ye_h2)
  double precision, intent(in ) :: Ts  (  Ts_l1:  Ts_h1,  Ts_l2:  Ts_h2)
  double precision, intent(in ) :: Yes ( Yes_l1: Yes_h1, Yes_l2: Yes_h2)
  double precision, intent(out) :: dedT(dedT_l1:dedT_h1,dedT_l2:dedT_h2)
  double precision, intent(out) :: dedY(dedY_l1:dedY_h1,dedY_l2:dedY_h2)
  integer, intent(in) :: validStar

  integer :: i, j
  double precision :: rhoinv, e1, e2
  type(eos_t) :: eos_state
  double precision :: dT, dYe, T1, T2, Ye1, Ye2
  double precision, parameter :: fac = 0.5d0, minfrac = 1.d-8

  !$OMP PARALLEL DO PRIVATE(i,j,rhoinv,eos_state,dT,dYe,T1,T2,Ye1,Ye2,e1,e2)
  do j=dedT_l2, dedT_h2
  do i=dedT_l1, dedT_h1

     rhoinv = 1.d0/S(i,j,URHO)
     eos_state % rho = S(i,j,URHO)
     eos_state % xn  = S(i,j,UFS:UFS+nspec-1) * rhoinv
     eos_state % aux = ye(i,j)

     if (validStar > 0) then
        dT = fac*abs(Ts(i,j) - T(i,j))
        dT = max(dT, minfrac*T(i,j))
        dYe = fac*abs(Yes(i,j) - Ye(i,j))
        dYe = max(dYe, minfrac*Ye(i,j))
     else
        dT = T(i,j) * 1.d-3 + 1.d-50
        dYe = 1.d-4
     end if

     T1 = T(i,j) - dT
     if (T1 < Tmin*(1.d0+1.d-6)) then
        T1 = T(i,j)
     end if
     T2 = T(i,j) + dT

     Ye1 = Ye(i,j) - dYe
     if (Ye1 < Ymin*(1.d0+1.d-6)) then
        Ye1 = Ye(i,j)
     end if
     Ye2 = Ye(i,j) + dYe

     eos_state % T = T1
     call eos(eos_input_rt, eos_state)
     e1 = eos_state % e

     eos_state % T = T2
     call eos(eos_input_rt, eos_state)
     e2 = eos_state % e

     dedT(i,j) = (e2-e1) / (T2-T1)

     eos_state % T = T(i,j)

     eos_state % aux = Ye1
     call eos(eos_input_rt, eos_state)
     e1 = eos_state % e

     eos_state % aux = Ye2
     call eos(eos_input_rt, eos_state)
     e2 = eos_state % e

     dedY(i,j) = (e2-e1)/(Ye2-Ye1)

  end do
  end do
  !$OMP END PARALLEL DO 

end subroutine ca_compute_dedx


subroutine ca_compute_eta_the( lo, hi, &
     & etaT,etaT_l1,etaT_l2,etaT_h1,etaT_h2, &
     & etTz,etTz_l1,etTz_l2,etTz_h1,etTz_h2, &
     & etaY,etaY_l1,etaY_l2,etaY_h1,etaY_h2, &
     & etYz,etYz_l1,etYz_l2,etYz_h1,etYz_h2, &
     & eta1,eta1_l1,eta1_l2,eta1_h1,eta1_h2, &
     & theT,theT_l1,theT_l2,theT_h1,theT_h2, &
     & thTz,thTz_l1,thTz_l2,thTz_h1,thTz_h2, &
     & theY,theY_l1,theY_l2,theY_h1,theY_h2, &
     & thYz,thYz_l1,thYz_l2,thYz_h1,thYz_h2, &
     & the1,the1_l1,the1_l2,the1_h1,the1_h2, &
     & djdT,djdT_l1,djdT_l2,djdT_h1,djdT_h2, &
     & djdY,djdY_l1,djdY_l2,djdY_h1,djdY_h2, &
     & dkdT,dkdT_l1,dkdT_l2,dkdT_h1,dkdT_h2, &
     & dkdY,dkdY_l1,dkdY_l2,dkdY_h1,dkdY_h2, &
     & dedT,dedT_l1,dedT_l2,dedT_h1,dedT_h2, &
     & dedY,dedY_l1,dedY_l2,dedY_h1,dedY_h2, &
     & Ers , Ers_l1, Ers_l2, Ers_h1, Ers_h2, &
     & rho , rho_l1, rho_l2, rho_h1, rho_h2, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: etaT_l1, etaT_h1, etaT_l2, etaT_h2
  integer, intent(in) :: etTz_l1, etTz_h1, etTz_l2, etTz_h2
  integer, intent(in) :: etaY_l1, etaY_h1, etaY_l2, etaY_h2
  integer, intent(in) :: etYz_l1, etYz_h1, etYz_l2, etYz_h2
  integer, intent(in) :: eta1_l1, eta1_h1, eta1_l2, eta1_h2
  integer, intent(in) :: theT_l1, theT_h1, theT_l2, theT_h2
  integer, intent(in) :: thTz_l1, thTz_h1, thTz_l2, thTz_h2
  integer, intent(in) :: theY_l1, theY_h1, theY_l2, theY_h2
  integer, intent(in) :: thYz_l1, thYz_h1, thYz_l2, thYz_h2
  integer, intent(in) :: the1_l1, the1_h1, the1_l2, the1_h2
  integer, intent(in) :: djdT_l1, djdT_h1, djdT_l2, djdT_h2
  integer, intent(in) :: djdY_l1, djdY_h1, djdY_l2, djdY_h2
  integer, intent(in) :: dkdT_l1, dkdT_h1, dkdT_l2, dkdT_h2
  integer, intent(in) :: dkdY_l1, dkdY_h1, dkdY_l2, dkdY_h2
  integer, intent(in) :: dedT_l1, dedT_h1, dedT_l2, dedT_h2
  integer, intent(in) :: dedY_l1, dedY_h1, dedY_l2, dedY_h2
  integer, intent(in) ::  Ers_l1,  Ers_h1,  Ers_l2,  Ers_h2
  integer, intent(in) ::  rho_l1,  rho_h1,  rho_l2,  rho_h2
  double precision,intent(out)  ::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2)
  double precision,intent(out)  ::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2)
  double precision,intent(out)  ::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2)
  double precision,intent(out)  ::etYz(etYz_l1:etYz_h1,etYz_l2:etYz_h2)
  double precision,intent(out)  ::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2)
  double precision,intent(out)  ::theT(theT_l1:theT_h1,theT_l2:theT_h2)
  double precision,intent(out)  ::thTz(thTz_l1:thTz_h1,thTz_l2:thTz_h2)
  double precision,intent(out)  ::theY(theY_l1:theY_h1,theY_l2:theY_h2)
  double precision,intent(out)  ::thYz(thYz_l1:thYz_h1,thYz_l2:thYz_h2)
  double precision,intent(out)  ::the1(the1_l1:the1_h1,the1_l2:the1_h2)
  double precision,intent(inout)::djdT(djdT_l1:djdT_h1,djdT_l2:djdT_h2,0:ngroups-1)
  double precision,intent(inout)::djdY(djdY_l1:djdY_h1,djdY_l2:djdY_h2,0:ngroups-1)
  double precision,intent(in )  ::dkdT(dkdT_l1:dkdT_h1,dkdT_l2:dkdT_h2,0:ngroups-1)
  double precision,intent(in )  ::dkdY(dkdY_l1:dkdY_h1,dkdY_l2:dkdY_h2,0:ngroups-1)
  double precision,intent(in )  ::dedT(dedT_l1:dedT_h1,dedT_l2:dedT_h2)
  double precision,intent(in )  ::dedY(dedY_l1:dedY_h1,dedY_l2:dedY_h2)
  double precision,intent(in )  ::Ers ( Ers_l1: Ers_h1, Ers_l2: Ers_h2,0:ngroups-1)
  double precision,intent(in )  ::rho ( rho_l1: rho_h1, rho_l2: rho_h2)
  double precision,intent(in) :: dt, tau

  integer :: i, j
  double precision :: cdt, det, et, ey, tt, ty, sigma
  double precision :: dZdT(0:ngroups-1), dZdY(0:ngroups-1)
  double precision :: sumdZdT, sumdZdY, fooT, fooY, barT, barY

  sigma = 1.d0 + tau
  cdt = clight * dt

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     dZdT = djdT(i,j,:) - dkdT(i,j,:)*Ers(i,j,:)
     dZdY = djdY(i,j,:) - dkdY(i,j,:)*Ers(i,j,:)

     sumdZdT = sum(dZdT)
     sumdZdY = sum(dZdY)

     if (sumdZdT .eq. 0.d0) then
        sumdZdT = 1.d-50
     end if

     if (sumdZdY .eq. 0.d0) then
        sumdZdY = 1.d-50
     end if

     fooT = cdt * sumdZdT
     fooY = cdt * sumdZdY
     barT = sigma*rho(i,j)*dedT(i,j)
     barY = sigma*rho(i,j)*dedY(i,j)

     et = sigma*rho(i,j) + cdt * sum(erg2rhoYe*dZdY)
     ey = cdt * sum(erg2rhoYe*dZdT)
     tt = barY + fooY
     ty = barT + fooT
     det = ty*et - tt*ey
     etaT(i,j) = (et * fooT) / det
     etaY(i,j) = (ey * fooY) / det
     theT(i,j) = (tt * fooT) / det
     theY(i,j) = (ty * fooY) / det

     etTz(i,j) = etaT(i,j) / sumdZdT
     etYz(i,j) = etaY(i,j) / sumdZdY
     thTz(i,j) = theT(i,j) / sumdZdT
     thYz(i,j) = theY(i,j) / sumdZdY

     eta1(i,j) = (et*barT - ey*barY) / det
     the1(i,j) = (sigma*rho(i,j) * ty) / det

     djdT(i,j,:) = dZdT / sumdZdT
     djdY(i,j,:) = dZdY / sumdZdY
  end do
  end do

end subroutine ca_compute_eta_the


subroutine ca_compute_rhs_neut(  &
     rhs , rhs_l1, rhs_l2, rhs_h1, rhs_h2, &
     jg  ,  jg_l1,  jg_l2,  jg_h1,  jg_h2, &
     mugT,mugT_l1,mugT_l2,mugT_h1,mugT_h2, &
     mugY,mugY_l1,mugY_l2,mugY_h1,mugY_h2, &
     cpT , cpT_l1, cpT_l2, cpT_h1, cpT_h2, &
     cpY , cpY_l1, cpY_l2, cpY_h1, cpY_h2, &
     etaT,etaT_l1,etaT_l2,etaT_h1,etaT_h2, &
     etay,etaY_l1,etaY_l2,etaY_h1,etaY_h2, &
     theT,theT_l1,theT_l2,theT_h1,theT_h2, &
     they,theY_l1,theY_l2,theY_h1,theY_h2, &
     Er2 , Er2_l1, Er2_l2, Er2_h1, Er2_h2, &
     re2 , re2_l1, re2_l2, re2_h1, re2_h2, &
     rY2 , rY2_l1, rY2_l2, rY2_h1, rY2_h2, &
     Ers , Ers_l1, Ers_l2, Ers_h1, Ers_h2, &
     res , res_l1, res_l2, res_h1, res_h2, &
     rYs , rYs_l1, rYs_l2, rYs_h1, rYs_h2, &
     r, dt, igroup, tau)

  use rad_params_module, only : ngroups, clight

  implicit none

  integer,intent(in):: rhs_l1, rhs_h1, rhs_l2, rhs_h2
  integer,intent(in)::  jg_l1,  jg_h1,  jg_l2,  jg_h2
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2
  integer,intent(in)::mugY_l1,mugY_h1,mugY_l2,mugY_h2
  integer,intent(in):: cpT_l1, cpT_h1, cpT_l2, cpT_h2
  integer,intent(in):: cpY_l1, cpY_h1, cpY_l2, cpY_h2
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2
  integer,intent(in)::etaY_l1,etaY_h1,etaY_l2,etaY_h2
  integer,intent(in)::theT_l1,theT_h1,theT_l2,theT_h2
  integer,intent(in)::theY_l1,theY_h1,theY_l2,theY_h2
  integer,intent(in):: Er2_l1, Er2_h1, Er2_l2, Er2_h2
  integer,intent(in):: re2_l1, re2_h1, re2_l2, re2_h2
  integer,intent(in):: rY2_l1, rY2_h1, rY2_l2, rY2_h2
  integer,intent(in):: Ers_l1, Ers_h1, Ers_l2, Ers_h2
  integer,intent(in):: res_l1, res_h1, res_l2, res_h2
  integer,intent(in):: rYs_l1, rYs_h1, rYs_l2, rYs_h2
  double precision,intent(out)::rhs ( rhs_l1: rhs_h1, rhs_l2: rhs_h2)
  double precision,intent(in )::jg  (  jg_l1:  jg_h1,  jg_l2:  jg_h2,0:ngroups-1)
  double precision,intent(in )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,0:ngroups-1)
  double precision,intent(in )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,0:ngroups-1)
  double precision,intent(in )::cpT ( cpT_l1: cpT_h1, cpT_l2: cpT_h2)
  double precision,intent(in )::cpY ( cpY_l1: cpY_h1, cpY_l2: cpY_h2)
  double precision,intent(in )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2)
  double precision,intent(in )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2)
  double precision,intent(in )::Er2 ( Er2_l1: Er2_h1, Er2_l2: Er2_h2,0:ngroups-1)
  double precision,intent(in )::re2 ( re2_l1: re2_h1, re2_l2: re2_h2)
  double precision,intent(in )::rY2 ( rY2_l1: rY2_h1, rY2_l2: rY2_h2)
  double precision,intent(in )::Ers ( Ers_l1: Ers_h1, Ers_l2: Ers_h2,0:ngroups-1)
  double precision,intent(in )::res ( res_l1: res_h1, res_l2: res_h2)
  double precision,intent(in )::rYs ( rYs_l1: rYs_h1, rYs_l2: rYs_h2)
  double precision,intent(in) ::   r( rhs_l1: rhs_h1)
  double precision,intent(in) :: dt, tau             
  integer, intent(in) :: igroup

  integer :: i, j
  double precision :: Hg, thetag, dt1

  dt1 = 1.d0/dt
  do j=rhs_l2, rhs_h2
  do i=rhs_l1, rhs_h1
     Hg = mugT(i,j,igroup)*etaT(i,j) - mugY(i,j,igroup)*etaY(i,j)
     thetag = mugY(i,j,igroup)*theY(i,j) - mugT(i,j,igroup)*theT(i,j)

     rhs(i,j) = clight*(jg(i,j,igroup) + Hg*cpT(i,j) + thetag*cpY(i,j)) &
          + dt1 * (Er2(i,j,igroup) - Hg*(res(i,j)-re2(i,j))  &
          &                        - thetag*(rYs(i,j)-rY2(i,j)) &
          &        + tau*Ers(i,j,igroup))

     rhs(i,j) = r(i) * rhs(i,j)
   end do
   end do

end subroutine ca_compute_rhs_neut


subroutine ca_local_accel_neut( lo, hi,  &
     Ern , Ern_l1, Ern_l2, Ern_h1, Ern_h2, &
     Erl , Erl_l1, Erl_l2, Erl_h1, Erl_h2, &
     kap , kap_l1, kap_l2, kap_h1, kap_h2, &
     etaT,etaT_l1,etaT_l2,etaT_h1,etaT_h2, &
     etaY,etaY_l1,etaY_l2,etaY_h1,etaY_h2, &
     theT,theT_l1,theT_l2,theT_h1,theT_h2, &
     theY,theY_l1,theY_l2,theY_h1,theY_h2, &
     mugT,mugT_l1,mugT_l2,mugT_h1,mugT_h2, &
     mugY,mugY_l1,mugY_l2,mugY_h1,mugY_h2, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer,intent(in):: lo(2), hi(2)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2
  integer,intent(in)::etaY_l1,etaY_h1,etaY_l2,etaY_h2
  integer,intent(in)::theT_l1,theT_h1,theT_l2,theT_h2
  integer,intent(in)::theY_l1,theY_h1,theY_l2,theY_h2
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2
  integer,intent(in)::mugY_l1,mugY_h1,mugY_l2,mugY_h2
  double precision,intent(inout)::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2,0:ngroups-1)
  double precision,intent(in   )::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2,0:ngroups-1)
  double precision,intent(in   )::kap ( kap_l1: kap_h1, kap_l2: kap_h2,0:ngroups-1)
  double precision,intent(in   )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2)
  double precision,intent(in   )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2)
  double precision,intent(in   )::theT(theT_l1:theT_h1,theT_l2:theT_h2)
  double precision,intent(in   )::theY(theY_l1:theY_h1,theY_l2:theY_h2)
  double precision,intent(in   )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,0:ngroups-1)
  double precision,intent(in   )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,0:ngroups-1)
  double precision,intent(in) :: dt, tau

  integer :: i, j
  double precision :: cdt1, rt, ry, p, q, r, s 
  double precision,dimension(0:ngroups-1)::Hg, Tg, epsilon, kapt, kk

  cdt1 = 1.d0/(clight*dt)

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     rt = sum(kap(i,j,:)*(Ern(i,j,:)-Erl(i,j,:)))
     ry = sum(erg2rhoYe(:)*kap(i,j,:)*(Ern(i,j,:)-Erl(i,j,:)))

     Hg = mugT(i,j,:)*etaT(i,j) - mugY(i,j,:)*etaY(i,j)
     Tg = -mugT(i,j,:)*theT(i,j) + mugY(i,j,:)*theY(i,j)

     kapt = kap(i,j,:) + (1.d0+tau)*cdt1
     kk = kap(i,j,:) / kapt

     p = 1.d0-sum(Hg*kk)
     q = 1.d0-sum(Tg*erg2rhoYe*kk)
     r = sum(Hg*erg2rhoYe*kk)
     s = sum(Tg*kk)

     epsilon = ((r*Tg + q*Hg) * rt + (s*Hg + p*Tg) * ry)  &
          / (kapt*(p*q-r*s) + 1.d-50)

     Ern(i,j,:) = Ern(i,j,:) + epsilon
  end do
  end do
  
end subroutine ca_local_accel_neut


subroutine ca_opac_emis_neut( lo, hi,  &
     Snew,Snew_l1,Snew_l2,Snew_h1,Snew_h2, &
     T   ,   T_l1,   T_l2,   T_h1,   T_h2, &
     Ye  ,  Ye_l1,  Ye_l2,  Ye_h1,  Ye_h2, &
     Ts  ,  Ts_l1,  Ts_l2,  Ts_h1,  Ts_h2, &
     Yes , Yes_l1, Yes_l2, Yes_h1, Yes_h2, &
     kpp , kpp_l1, kpp_l2, kpp_h1, kpp_h2, &
     kpr , kpr_l1, kpr_l2, kpr_h1, kpr_h2, &
     jg  ,   j_l1,   j_l2,   j_h1,   j_h2, &
     djdT,djdT_l1,djdT_l2,djdT_h1,djdT_h2, &
     djdY,djdY_l1,djdY_l2,djdY_h1,djdY_h2, &
     dkdT,dkdT_l1,dkdT_l2,dkdT_h1,dkdT_h2, &
     dkdY,dkdY_l1,dkdY_l2,dkdY_h1,dkdY_h2, &
     use_dkdT, validStar, lag_opac) 

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: Snew_l1, Snew_h1, Snew_l2, Snew_h2 
  integer, intent(in) ::    T_l1,    T_h1,    T_l2,    T_h2
  integer, intent(in) ::   Ye_l1,   Ye_h1,   Ye_l2,   Ye_h2
  integer, intent(in) ::   Ts_l1,   Ts_h1,   Ts_l2,   Ts_h2
  integer, intent(in) ::  Yes_l1,  Yes_h1,  Yes_l2,  Yes_h2
  integer, intent(in) ::  kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2 
  integer, intent(in) ::  kpr_l1,  kpr_h1,  kpr_l2,  kpr_h2
  integer, intent(in) ::    j_l1,    j_h1,    j_l2,    j_h2
  integer, intent(in) :: djdT_l1, djdT_h1, djdT_l2, djdT_h2
  integer, intent(in) :: djdY_l1, djdY_h1, djdY_l2, djdY_h2
  integer, intent(in) :: dkdT_l1, dkdT_h1, dkdT_l2, dkdT_h2
  integer, intent(in) :: dkdY_l1, dkdY_h1, dkdY_l2, dkdY_h2
  double precision, intent(in ) :: Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,NVAR)
  double precision, intent(in ) :: T   (   T_l1:   T_h1,   T_l2:   T_h2)
  double precision, intent(in ) :: Ye  (  Ye_l1:  Ye_h1,  Ye_l2:  Ye_h2)
  double precision, intent(in ) :: Ts  (  Ts_l1:  Ts_h1,  Ts_l2:  Ts_h2)
  double precision, intent(in ) :: Yes ( Yes_l1: Yes_h1, Yes_l2: Yes_h2)
  double precision, intent(inout) :: kpp ( kpp_l1: kpp_h1, kpp_l2: kpp_h2,0:ngroups-1)
  double precision, intent(inout) :: kpr ( kpr_l1: kpr_h1, kpr_l2: kpr_h2,0:ngroups-1)
  double precision, intent(out  ) :: jg  (   j_l1:   j_h1,   j_l2:   j_h2,0:ngroups-1)
  double precision, intent(out  ) :: djdT(djdT_l1:djdT_h1,djdT_l2:djdT_h2,0:ngroups-1)
  double precision, intent(out  ) :: djdY(djdY_l1:djdY_h1,djdY_l2:djdY_h2,0:ngroups-1)
  double precision, intent(inout) :: dkdT(dkdT_l1:dkdT_h1,dkdT_l2:dkdT_h2,0:ngroups-1)
  double precision, intent(inout) :: dkdY(dkdY_l1:dkdY_h1,dkdY_l2:dkdY_h2,0:ngroups-1)
  integer, intent(in) :: use_dkdT, validStar, lag_opac

  integer :: i, j, g, inu
  double precision :: ab, sc, delta, eta, er, der, rho, temp
  double precision :: ab1, sc1, delta1, eta1
  double precision :: ab2, sc2, delta2, eta2 
  double precision :: dT, dYe
  double precision :: Bg, Bg1, Bg2
  logical :: comp_ab, comp_sc, comp_eta
  double precision, parameter :: fac = 0.5d0, minfrac = 1.d-8

  !$OMP PARALLEL DO PRIVATE(i,j,g,rho,temp,dT,dYe,inu,er,der,comp_ab,comp_sc,comp_eta) &
  !$OMP PRIVATE(ab ,sc ,delta ,eta ,Bg ) &
  !$OMP PRIVATE(ab1,sc1,delta1,eta1,Bg1) &
  !$OMP PRIVATE(ab2,sc2,delta2,eta2,Bg2)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     
     rho = Snew(i,j,URHO)
     temp = T(i,j)

     if (validStar > 0) then
        dT = fac*abs(Ts(i,j) - T(i,j))
        dT = max(dT, minfrac*T(i,j))
        dYe = fac*abs(Yes(i,j) - Ye(i,j))
        dYe = max(dYe, minfrac*Ye(i,j))
     else
        dT = T(i,j) * 1.d-3 + 1.d-50
        dYe = 1.d-4
     end if

     do g=0, ngroups-1
        
        call prep_opacity(g, inu, er, der)

        if (lag_opac .eq. 1) then

           dkdT(i,j,g) = 0.d0
           dkdY(i,j,g) = 0.d0              

           comp_ab = .true.
           comp_sc = .false.
           comp_eta = .true.

           call get_opacity_emissivity(ab, sc, delta, eta, &
                rho, Ye(i,j), temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg = eta * der / ab
           jg(i,j,g) = Bg * kpp(i,j,g)

           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j)-dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg1 = eta1 * der / ab1

           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j)+dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg2 = eta2 * der / ab2

           djdY(i,j,g) = (Bg2-Bg1)*kpp(i,j,g)/(2.d0*dye)

           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j), temp-dT, er, inu, comp_ab, comp_sc, comp_eta)
           Bg1 = eta1 * der / ab1
           
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j), temp+dT, er, inu, comp_ab, comp_sc, comp_eta)
           Bg2 = eta2 * der / ab2

           djdT(i,j,g) = (Bg2-Bg1)*kpp(i,j,g)/(2.d0*dT)
           
        else

           comp_ab = .true.
           comp_sc = .true.
           comp_eta = .true.
           
           call get_opacity_emissivity(ab, sc, delta, eta, &
                rho, Ye(i,j), temp, er, inu, comp_ab, comp_sc, comp_eta)
           kpp(i,j,g) = ab
           kpr(i,j,g) = ab + sc * (1.d0 - delta/3.d0)
           jg(i,j,g) = eta * der 
           
           comp_ab = .true.
           comp_sc = .false.
           comp_eta = .true.
           
           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j)-dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j)+dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           djdY(i,j,g) = (eta2-eta1)*der/(2.d0*dye)
           dkdY(i,j,g) = (ab2-ab1)/(2.d0*dye)
           if (use_dkdT .eq. 0) then
              djdY(i,j,g) = (eta2/ab2 - eta1/ab1)*der/(2.d0*dye) * ab
              dkdY(i,j,g) = 0.d0
           end if
           
           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j), temp-dT, er, inu, comp_ab, comp_sc, comp_eta)
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j), temp+dT, er, inu, comp_ab, comp_sc, comp_eta)
           djdT(i,j,g) = (eta2-eta1)*der/(2.d0*dT)
           dkdT(i,j,g) = (ab2-ab1)/(2.d0*dT)
           if (use_dkdT .eq. 0) then        
              djdT(i,j,g) = (eta2/ab2 - eta1/ab1)*der/(2.d0*dT) * ab
              dkdT(i,j,g) = 0.d0
           end if

        end if

     end do
  end do
  end do
  !$OMP END PARALLEL DO 

end subroutine ca_opac_emis_neut


subroutine ca_state_update_neut( lo, hi, &
     state,state_l1,state_l2,state_h1,state_h2, &
     rhoe,  rhoe_l1, rhoe_l2, rhoe_h1, rhoe_h2, &
     Ye  ,    Ye_l1,   Ye_l2,   Ye_h1,   Ye_h2, &
     temp,  temp_l1, temp_l2, temp_h1, temp_h2, &
     msk ,   msk_l1,  msk_l2,  msk_h1,  msk_h2, &
     derat, dTrat, dye)

  use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, UFS, UFX

  implicit none

  integer, intent(in) :: lo(2), hi(2) 
  integer, intent(in) :: state_l1, state_h1, state_l2, state_h2
  integer, intent(in) ::  rhoe_l1,  rhoe_h1,  rhoe_l2,  rhoe_h2
  integer, intent(in) ::    Ye_l1,    Ye_h1,    Ye_l2,    Ye_h2
  integer, intent(in) ::  temp_l1,  temp_h1,  temp_l2,  temp_h2
  integer, intent(in) ::   msk_l1,   msk_h1,   msk_l2,   msk_h2
  double precision, intent(in)   :: rhoe( rhoe_l1: rhoe_h1, rhoe_l2: rhoe_h2)
  double precision, intent(in)   ::   Ye(   Ye_l1:   Ye_h1,   Ye_l2:   Ye_h2)
  double precision, intent(in)   :: temp( temp_l1: temp_h1, temp_l2: temp_h2)
  double precision, intent(in)   ::  msk(  msk_l1:  msk_h1,  msk_l2:  msk_h2)
  double precision, intent(inout)::state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision, intent(inout) :: derat, dTrat, dye

  integer :: i, j
  double precision :: ei, ek, Told, Yeold

  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     ei = state(i,j,UEINT)
     derat = max(derat, abs((rhoe(i,j) - ei)*msk(i,j)/ max(ei, 1.d-50)))
     ek = state(i,j,UEDEN) - state(i,j,UEINT)
     state(i,j,UEINT) = rhoe(i,j)
     state(i,j,UEDEN) = rhoe(i,j) + ek

     Told = state(i,j,UTEMP);
     dTrat = max(dTrat, abs((temp(i,j)-Told)*msk(i,j)/ max(Told, 1.d-50)))
     state(i,j,UTEMP) = temp(i,j)

     Yeold = state(i,j,UFX) / state(i,j,URHO)
     dye = max(dye, abs((Ye(i,j)-Yeold)*msk(i,j)))
     state(i,j,UFX) = state(i,j,URHO) * Ye(i,j)
  end do
  end do

end subroutine ca_state_update_neut


subroutine ca_update_matter_neut( lo, hi,  &
     re_n,re_n_l1,re_n_l2,re_n_h1,re_n_h2, &
     rY_n,rY_n_l1,rY_n_l2,rY_n_h1,rY_n_h2, &
     Ye_n,Ye_n_l1,Ye_n_l2,Ye_n_h1,Ye_n_h2, &
     Er_n,Er_n_l1,Er_n_l2,Er_n_h1,Er_n_h2, &
     Er_l,Er_l_l1,Er_l_l2,Er_l_h1,Er_l_h2, &
     re_s,re_s_l1,re_s_l2,re_s_h1,re_s_h2, &
     rY_s,rY_s_l1,rY_s_l2,rY_s_h1,rY_s_h2, &
     re_2,re_2_l1,re_2_l2,re_2_h1,re_2_h2, &
     rY_2,rY_2_l1,rY_2_l2,rY_2_h1,rY_2_h2, &
     etaT,etaT_l1,etaT_l2,etaT_h1,etaT_h2, &
     etaY,etaY_l1,etaY_l2,etaY_h1,etaY_h2, &
     eta1,eta1_l1,eta1_l2,eta1_h1,eta1_h2, &
     theT,theT_l1,theT_l2,theT_h1,theT_h2, &
     theY,theY_l1,theY_l2,theY_h1,theY_h2, &
     the1,the1_l1,the1_l2,the1_h1,the1_h2, &
      cpt, cpt_l1, cpt_l2, cpt_h1, cpt_h2, &
      cpy, cpy_l1, cpy_l2, cpy_h1, cpy_h2, &
      kpp, kpp_l1, kpp_l2, kpp_h1, kpp_h2, &
     mugT,mugT_l1,mugT_l2,mugT_h1,mugT_h2, &
     mugY,mugY_l1,mugY_l2,mugY_h1,mugY_h2, &
     Snew,Snew_l1,Snew_l2,Snew_h1,Snew_h2, &
     dt, tau)

  use rad_params_module, only : ngroups, erg2rhoYe, clight
  use meth_params_module, only : NVAR, URHO

  implicit none

  integer,intent(in)::lo(2),hi(2)
  integer,intent(in)::re_n_l1, re_n_h1, re_n_l2, re_n_h2
  integer,intent(in)::rY_n_l1, rY_n_h1, rY_n_l2, rY_n_h2
  integer,intent(in)::Ye_n_l1, Ye_n_h1, Ye_n_l2, Ye_n_h2
  integer,intent(in)::Er_n_l1, Er_n_h1, Er_n_l2, Er_n_h2
  integer,intent(in)::Er_l_l1, Er_l_h1, Er_l_l2, Er_l_h2
  integer,intent(in)::re_s_l1, re_s_h1, re_s_l2, re_s_h2
  integer,intent(in)::rY_s_l1, rY_s_h1, rY_s_l2, rY_s_h2
  integer,intent(in)::re_2_l1, re_2_h1, re_2_l2, re_2_h2
  integer,intent(in)::rY_2_l1, rY_2_h1, rY_2_l2, rY_2_h2
  integer,intent(in)::etaT_l1, etaT_h1, etaT_l2, etaT_h2
  integer,intent(in)::etaY_l1, etaY_h1, etaY_l2, etaY_h2
  integer,intent(in)::eta1_l1, eta1_h1, eta1_l2, eta1_h2
  integer,intent(in)::theT_l1, theT_h1, theT_l2, theT_h2
  integer,intent(in)::theY_l1, theY_h1, theY_l2, theY_h2
  integer,intent(in)::the1_l1, the1_h1, the1_l2, the1_h2
  integer,intent(in):: cpt_l1,  cpt_h1,  cpt_l2,  cpt_h2
  integer,intent(in):: cpy_l1,  cpy_h1,  cpy_l2,  cpy_h2
  integer,intent(in):: kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2
  integer,intent(in)::mugT_l1, mugT_h1, mugT_l2, mugT_h2
  integer,intent(in)::mugY_l1, mugY_h1, mugY_l2, mugY_h2
  integer,intent(in)::Snew_l1, Snew_h1, Snew_l2, Snew_h2
  double precision,intent(out)::re_n(re_n_l1:re_n_h1,re_n_l2:re_n_h2)
  double precision,intent(out)::rY_n(rY_n_l1:rY_n_h1,rY_n_l2:rY_n_h2)
  double precision,intent(out)::Ye_n(Ye_n_l1:Ye_n_h1,Ye_n_l2:Ye_n_h2)
  double precision,intent(in )::Er_n(Er_n_l1:Er_n_h1,Er_n_l2:Er_n_h2,0:ngroups-1)
  double precision,intent(in )::Er_l(Er_l_l1:Er_l_h1,Er_l_l2:Er_l_h2,0:ngroups-1)
  double precision,intent(in )::re_s(re_s_l1:re_s_h1,re_s_l2:re_s_h2)
  double precision,intent(in )::rY_s(rY_s_l1:rY_s_h1,rY_s_l2:rY_s_h2)
  double precision,intent(in )::re_2(re_2_l1:re_2_h1,re_2_l2:re_2_h2)
  double precision,intent(in )::rY_2(rY_2_l1:rY_2_h1,rY_2_l2:rY_2_h2)
  double precision,intent(in )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2)
  double precision,intent(in )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2)
  double precision,intent(in )::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2)
  double precision,intent(in )::the1(the1_l1:the1_h1,the1_l2:the1_h2)
  double precision,intent(in ):: cpt( cpt_l1: cpt_h1, cpt_l2: cpt_h2)
  double precision,intent(in ):: cpy( cpy_l1: cpy_h1, cpy_l2: cpy_h2)
  double precision,intent(in ):: kpp( kpp_l1: kpp_h1, kpp_l2: kpp_h2,0:ngroups-1)
  double precision,intent(in )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,0:ngroups-1)
  double precision,intent(in )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,0:ngroups-1)
  double precision,intent(in )::Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,NVAR)
  double precision,intent(in) :: dt, tau

  integer :: i,j,g
  double precision :: cdt, H1, Theta, Thbar1, Hbar 
  double precision :: dkEE, dkEEY, foo, chg, chgY

  cdt = clight * dt
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     H1 = eta1(i,j)
     Thbar1 = the1(i,j)

     Theta = theY(i,j) - theT(i,j)
     Hbar = sum(erg2rhoYe*(mugT(i,j,:)*etaT(i,j) - mugY(i,j,:)*etaY(i,j)))

     dkEE = 0.d0
     dkEEY = 0.d0
     do g=0, ngroups-1
        foo = kpp(i,j,g)*(Er_n(i,j,g)-Er_l(i,j,g))
        dkEE = dkEE + foo
        dkEEY = dkEEY + foo * erg2rhoYe(g)
     end do

     chg = cdt*dkEE + H1*((re_2(i,j)-re_s(i,j)) + cdt*cpt(i,j)) &
          - Theta*((rY_2(i,j)-rY_s(i,j)) + cdt*cpy(i,j))

     chgY = cdt*dkEEY + Thbar1*((rY_2(i,j)-rY_s(i,j)) + cdt*cpy(i,j)) &
          - Hbar*((re_2(i,j)-re_s(i,j)) + cdt*cpt(i,j))

     re_n(i,j) = re_s(i,j) + chg
     rY_n(i,j) = rY_s(i,j) + chgY

     re_n(i,j) = (re_n(i,j) + tau*re_s(i,j)) / (1.d0+tau)
     rY_n(i,j) = (rY_n(i,j) + tau*rY_s(i,j)) / (1.d0+tau)

     Ye_n(i,j) = rY_n(i,j) / Snew(i,j,URHO)

     ! temperature will be updated after exiting this subroutine
  end do
  end do

end subroutine ca_update_matter_neut


subroutine ca_ncupdate_matter_neut( lo, hi,  &
     Tp_n, Tp_n_l1, Tp_n_l2, Tp_n_h1, Tp_n_h2,  &
     Ye_n, Ye_n_l1, Ye_n_l2, Ye_n_h1, Ye_n_h2,  &
     Er_n, Er_n_l1, Er_n_l2, Er_n_h1, Er_n_h2,  &
     re_s, re_s_l1, re_s_l2, re_s_h1, re_s_h2,  &
     rY_s, rY_s_l1, rY_s_l2, rY_s_h1, rY_s_h2,  &
     re_2, re_2_l1, re_2_l2, re_2_h1, re_2_h2,  &
     rY_2, rY_2_l1, rY_2_l2, rY_2_h1, rY_2_h2,  &
     etTz, etTz_l1, etTz_l2, etTz_h1, etTz_h2,  &
     etYz, etYz_l1, etYz_l2, etYz_h1, etYz_h2,  &
     thTz, thTz_l1, thTz_l2, thTz_h1, thTz_h2,  &
     thYz, thYz_l1, thYz_l2, thYz_h1, thYz_h2,  &
      kpp,  kpp_l1,  kpp_l2,  kpp_h1,  kpp_h2,  &
       jg,   jg_l1,   jg_l2,   jg_h1,   jg_h2,  &
     dt)

  use rad_params_module, only : ngroups, erg2rhoYe, clight
  use meth_params_module, only : NVAR, URHO

  implicit none

  integer,intent(in)::lo(2),hi(2)
  integer,intent(in)::Tp_n_l1, Tp_n_h1, Tp_n_l2, Tp_n_h2
  integer,intent(in)::Ye_n_l1, Ye_n_h1, Ye_n_l2, Ye_n_h2
  integer,intent(in)::Er_n_l1, Er_n_h1, Er_n_l2, Er_n_h2
  integer,intent(in)::re_s_l1, re_s_h1, re_s_l2, re_s_h2
  integer,intent(in)::rY_s_l1, rY_s_h1, rY_s_l2, rY_s_h2
  integer,intent(in)::re_2_l1, re_2_h1, re_2_l2, re_2_h2
  integer,intent(in)::rY_2_l1, rY_2_h1, rY_2_l2, rY_2_h2
  integer,intent(in)::etTz_l1, etTz_h1, etTz_l2, etTz_h2
  integer,intent(in)::etYz_l1, etYz_h1, etYz_l2, etYz_h2
  integer,intent(in)::thTz_l1, thTz_h1, thTz_l2, thTz_h2
  integer,intent(in)::thYz_l1, thYz_h1, thYz_l2, thYz_h2
  integer,intent(in):: kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2
  integer,intent(in)::  jg_l1,   jg_h1,   jg_l2,   jg_h2
  double precision,intent(inout)::Tp_n(Tp_n_l1:Tp_n_h1,Tp_n_l2:Tp_n_h2)
  double precision,intent(inout)::Ye_n(Ye_n_l1:Ye_n_h1,Ye_n_l2:Ye_n_h2)
  double precision,intent(in   )::Er_n(Er_n_l1:Er_n_h1,Er_n_l2:Er_n_h2,0:ngroups-1)
  double precision,intent(in   )::re_s(re_s_l1:re_s_h1,re_s_l2:re_s_h2)
  double precision,intent(in   )::rY_s(rY_s_l1:rY_s_h1,rY_s_l2:rY_s_h2)
  double precision,intent(in   )::re_2(re_2_l1:re_2_h1,re_2_l2:re_2_h2)
  double precision,intent(in   )::rY_2(rY_2_l1:rY_2_h1,rY_2_l2:rY_2_h2)
  double precision,intent(in   )::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2)
  double precision,intent(in   )::etYz(etYz_l1:etYz_h1,etYz_l2:etYz_h2)
  double precision,intent(in   )::thTz(thTz_l1:thTz_h1,thTz_l2:thTz_h2)
  double precision,intent(in   )::thYz(thYz_l1:thYz_h1,thYz_l2:thYz_h2)
  double precision,intent(in   ):: kpp( kpp_l1: kpp_h1, kpp_l2: kpp_h2,0:ngroups-1)
  double precision,intent(in   )::  jg(  jg_l1:  jg_h1,  jg_l2:  jg_h2,0:ngroups-1)
  double precision,intent(in) :: dt

   integer :: i,j,g
   double precision :: cdt1, cpT, cpY, foo, scrch_re, scrch_rY
   double precision :: dTemp, dYe
   double precision, parameter :: fac = 0.01d0

   cdt1 = 1.d0 / (clight * dt)
   do j = lo(2), hi(2)
   do i = lo(1), hi(1)

      cpT = 0.d0
      cpY = 0.d0
      do g = 0, ngroups-1
         foo = kpp(i,j,g)*Er_n(i,j,g) - jg(i,j,g)
         cpT = cpT + foo
         cpY = cpY + erg2rhoYe(g)*foo
      end do

      scrch_re = cpT - (re_s(i,j) - re_2(i,j)) * cdt1
      scrch_rY = cpY - (rY_s(i,j) - rY_2(i,j)) * cdt1

      dTemp = etTz(i,j)*scrch_re - thTz(i,j)*scrch_rY
      dYe   = -etYz(i,j)*scrch_re + thYz(i,j)*scrch_rY

      if (abs(dTemp/(Tp_n(i,j)+1.d-50)) > fac) then
         dTemp = sign(fac*Tp_n(i,j), dTemp)
      end if

      if (abs(dYe/(Ye_n(i,j)+1.d-50)) > fac) then
         dYe = sign(fac*Ye_n(i,j), dYe)
      end if

     Tp_n(i,j) = Tp_n(i,j) + dTemp
     Ye_n(i,j) = Ye_n(i,j) + dYe

  end do
  end do

end subroutine ca_ncupdate_matter_neut


subroutine ca_compute_rosseland_neut( lo, hi,  &
     kpr , kpr_l1, kpr_l2, kpr_h1, kpr_h2, &
     stat,stat_l1,stat_l2,stat_h1,stat_h2 )

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) ::  kpr_l1,  kpr_h1,  kpr_l2,  kpr_h2
  integer, intent(in) :: stat_l1, stat_h1, stat_l2, stat_h2 
  double precision, intent(out) :: kpr ( kpr_l1: kpr_h1, kpr_l2: kpr_h2,0:ngroups-1)
  double precision, intent(in ) :: stat(stat_l1:stat_h1,stat_l2:stat_h2,NVAR)

  integer :: i, j, g, inu
  double precision :: ab, sc, delta, eta, er, der, rho, temp, Ye
  logical :: comp_ab, comp_sc, comp_eta

  comp_ab = .true.
  comp_sc = .true.
  comp_eta = .false.

  !$OMP PARALLEL DO SHARED(comp_ab,comp_sc,comp_eta) &
  !$OMP PRIVATE(g,inu,er,der,i,j,rho,temp,Ye,ab,sc,delta,eta)
  do g=0, ngroups-1

     call prep_opacity(g, inu, er, der)

     do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rho = stat(i,j,URHO)
        temp = stat(i,j,UTEMP)
        Ye = stat(i,j,UFX) / rho

        call get_opacity_emissivity(ab, sc, delta, eta, &
             rho, Ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

        kpr(i,j,g) = ab + sc * (1.d0 - delta/3.d0)

     end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine ca_compute_rosseland_neut


subroutine ca_compute_planck_neut( lo, hi,  &
     kpp , kpp_l1, kpp_l2, kpp_h1, kpp_h2, &
     stat,stat_l1,stat_l2,stat_h1,stat_h2 )

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) ::  kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2
  integer, intent(in) :: stat_l1, stat_h1, stat_l2, stat_h2 
  double precision, intent(out) :: kpp ( kpp_l1: kpp_h1, kpp_l2: kpp_h2,0:ngroups-1)
  double precision, intent(in ) :: stat(stat_l1:stat_h1,stat_l2:stat_h2,NVAR)

  integer :: i, j, g, inu
  double precision :: ab, sc, delta, eta, er, der, rho, temp, Ye
  logical :: comp_ab, comp_sc, comp_eta

  comp_ab = .true.
  comp_sc = .false.
  comp_eta = .false.

  !$OMP PARALLEL DO SHARED(comp_ab,comp_sc,comp_eta) &
  !$OMP PRIVATE(g,inu,er,der,i,j,rho,temp,Ye,ab,sc,delta,eta)
  do g=0, ngroups-1

     call prep_opacity(g, inu, er, der)

     do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rho = stat(i,j,URHO)
        temp = stat(i,j,UTEMP)
        Ye = stat(i,j,UFX) / rho

        call get_opacity_emissivity(ab, sc, delta, eta, &
             rho, Ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

        kpp(i,j,g) = ab 

     end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine ca_compute_planck_neut
