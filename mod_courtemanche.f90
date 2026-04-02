! This module implements the Courtemanche Model  (M Courtemanche, RJ Ramirez, S Nattel; Am J Physiol 275:H301-H321, 1998)
! Added IKACh: formulation of Grandi (E.Grandi,S.V.Pandit,N.Voigt,A.J.Workman,D.Dobrev,J.Jalife,D.M.Bers;Circulation Research,vol.109,pp.1055-1066,2011)
! Implemented by Laura Martínez Mateu (laumarma@ci2b.upv.es)  26/02/2019
! Checked by Violeta Puche García (viopucga@ci2b.upv.es) 08/02/2023 + blocking factors due to drug action
!
!------------------------------------------------------------------------------
module mod_courtemanche
!------------------------------------------------------------------------------
  use mod_precision
  implicit none   ! Para que Fortran me obligue a declarar las variables
  include 'courtemanche_parameters.inc'
  !.
 type t_courte
    sequence
    real(rp) ::   u, v, w, d, xnO, xr, Nai, Ki, Carel, oi, ui, xnIS, xnDpIC3, &
                  f, xs, oa, ua, fCa, Cai, Caup, xnIC3, xnIC2, xnIF, xnC3,    &
                  xnC2, xnC1, xnDpIC2, xnDpIF, xnDpC3, xnDpC2, xnDpC1, xnDpO, & 
                  xnDpIS, xnDpIT, xnDIC3, xnDIC2, xnDIF, xnDC3, xnDC2, xnDC1, &
                  xnDO, xnDIS, xnDIT
  end type t_courte
  
  type t_cur
    sequence
    real(rp) ::   INa,IK1,Ito,IKur,IKr,IKs,ICaL,INaK,INaCa,IbCa,IbNa,IpCa,Irel,&
                  Itr,Iup,Ileak,IKACh
  end type t_cur
  
  !.
  public  :: get_parameter_CRN , ic_Courte, Courtemanche_A_P01
  
  private :: put_param, f_currents, f_gates, f_concentrations, safe_exp
  !.
  type, public:: t_prm
    private
    real(rp)  :: p_gNa    ! 1	
    real(rp)  :: p_gCaL   ! 2  
    real(rp)  :: p_gto    ! 3
    real(rp)  :: p_gKur   ! 4
    real(rp)  :: p_gKr    ! 5 
    real(rp)  :: p_gKs    ! 6 
    real(rp)  :: p_gK1    ! 7 
    real(rp)  :: p_gbCa   ! 8
    real(rp)  :: p_gbNa   ! 9 
    real(rp)  :: p_ACh    ! 10
    real(rp)  :: p_IKACh  ! 11
  end type t_prm
  


!------------------------------------------------------------------------------!
!  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  !
!------------------------------------------------------------------------------!
contains
!------------------------------------------------------------------------------!
pure real(rp) function safe_exp(x)
!------------------------------------------------------------------------------!
  implicit none
  real(rp), intent(in) :: x
  real(rp), parameter  :: EXP_ARG_MAX = 80.0_rp

  safe_exp = exp(max(-EXP_ARG_MAX, min(EXP_ARG_MAX, x)))
end function safe_exp
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
subroutine write_state_Courte(lu_st)
implicit none
integer(ip), intent(in)   :: lu_st
    write(lu_st,'(A)',advance='no') '% '
    write(lu_st,'(A)',advance='no') 'u '
    write(lu_st,'(A)',advance='no') 'v '
    write(lu_st,'(A)',advance='no') 'w '
    write(lu_st,'(A)',advance='no') 'd '
    write(lu_st,'(A)',advance='no') 'xnO '
    write(lu_st,'(A)',advance='no') 'xr '
    write(lu_st,'(A)',advance='no') 'Nai '
    write(lu_st,'(A)',advance='no') 'Ki '
    write(lu_st,'(A)',advance='no') 'Carel '
    write(lu_st,'(A)',advance='no') 'oi '
    write(lu_st,'(A)',advance='no') 'ui '
    write(lu_st,'(A)',advance='no') 'xnIS '
    write(lu_st,'(A)',advance='no') 'xnDpIC3 '
    write(lu_st,'(A)',advance='no') 'f '
    write(lu_st,'(A)',advance='no') 'xs '
    write(lu_st,'(A)',advance='no') 'oa '
    write(lu_st,'(A)',advance='no') 'ua '
    write(lu_st,'(A)',advance='no') 'fCa '
    write(lu_st,'(A)',advance='no') 'Cai '
    write(lu_st,'(A)',advance='no') 'Caup '   
    write(lu_st,'(A)',advance='no') 'xnIC3 '
    write(lu_st,'(A)',advance='no') 'xnIC2 '
    write(lu_st,'(A)',advance='no') 'xnIF '
    write(lu_st,'(A)',advance='no') 'xnC3 '
    write(lu_st,'(A)',advance='no') 'xnC2 '
    write(lu_st,'(A)',advance='no') 'xnC1 '
    write(lu_st,'(A)',advance='no') 'xnDpIC2 '
    write(lu_st,'(A)',advance='no') 'xnDpIF '
    write(lu_st,'(A)',advance='no') 'xnDpC3 '
    write(lu_st,'(A)',advance='no') 'xnDpC2 '
    write(lu_st,'(A)',advance='no') 'xnDpC1 '
    write(lu_st,'(A)',advance='no') 'xnDpO '
    write(lu_st,'(A)',advance='no') 'xnDpIS '
    write(lu_st,'(A)',advance='no') 'xnDpIT '
    write(lu_st,'(A)',advance='no') 'xnDIC3 '
    write(lu_st,'(A)',advance='no') 'xnDIC2 '
    write(lu_st,'(A)',advance='no') 'xnDIF '
    write(lu_st,'(A)',advance='no') 'xnDC3 '
    write(lu_st,'(A)',advance='no') 'xnDC2 '
    write(lu_st,'(A)',advance='no') 'xnDC1 '
    write(lu_st,'(A)',advance='no') 'xnDO '
    write(lu_st,'(A)',advance='no') 'xnDIS '
    write(lu_st,'(A)',advance='no') 'xnDIT '

    
end subroutine write_state_Courte
!------------------------------------------------------------------------------!
subroutine write_current_Courte(lu_cr)
implicit none
integer(ip), intent(in)   :: lu_cr

    write(lu_cr,'(A)',advance='no') '% '
    write(lu_cr,'(A)',advance='no') 'INa '
    write(lu_cr,'(A)',advance='no') 'IK1 '
    write(lu_cr,'(A)',advance='no') 'Ito '
    write(lu_cr,'(A)',advance='no') 'IKur '
    write(lu_cr,'(A)',advance='no') 'IKr '
    write(lu_cr,'(A)',advance='no') 'IKs '
    write(lu_cr,'(A)',advance='no') 'ICaL '
    write(lu_cr,'(A)',advance='no') 'INaK '
    write(lu_cr,'(A)',advance='no') 'INaCa '
    write(lu_cr,'(A)',advance='no') 'IbCa '
    write(lu_cr,'(A)',advance='no') 'IbNa '
    write(lu_cr,'(A)',advance='no') 'IpCa '
    write(lu_cr,'(A)',advance='no') 'Irel '
    write(lu_cr,'(A)',advance='no') 'Itr '
    write(lu_cr,'(A)',advance='no') 'Iup '
    write(lu_cr,'(A)',advance='no') 'Ileak '
    write(lu_cr,'(A)',advance='no') 'IKACh '
    write(lu_cr,'(A)',advance='no') 'Iion '

end subroutine write_current_Courte
!------------------------------------------------------------------------------!
function get_parameter_CRN (tcell) result (v_prm)
!------------------------------------------------------------------------------!
  implicit none
  integer(ip), intent(in) :: tcell
  real (rp)               :: v_prm(np_courte)

  select case (tcell)
  case (Courte_mod)
    v_prm = (/ p_gNa, p_gCaL, p_gto, p_gKur, p_gKr, p_gKs, p_gK1_RA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_RA/)
  case (Courte_mod_RA_PM)
    v_prm = (/ p_gNa, p_gCaL, p_gto, p_gKur, p_gKr, p_gKs, p_gK1_RA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_RA/)
  case (Courte_mod_CT_BBra)
    v_prm = (/ p_gNa, p_gCaL_CT_BBra, p_gto, p_gKur, p_gKr, p_gKs, p_gK1_RA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_RA/)
  case (Courte_mod_TVR)
    v_prm = (/ p_gNa, p_gCaL_TVR, p_gto, p_gKur, p_gKr_TVR, p_gKs, p_gK1_RA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_RA/)
  case (Courte_mod_RAA)
    v_prm = (/ p_gNa, p_gCaL_RAA, p_gto_RAA, p_gKur, p_gKr_RAA, p_gKs, p_gK1_RA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_RA/)
  case (Courte_mod_LA)
    v_prm = (/ p_gNa, p_gCaL_LA, p_gto, p_gKur, p_gKr_LA, p_gKs, p_gK1_LA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_LA/)
  case (Courte_mod_PV)
    v_prm = (/ p_gNa, p_gCaL_PV, p_gto_PV, p_gKur, p_gKr_PV, p_gKs_PV, p_gK1_PV, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_LA/)
  case (Courte_mod_LAA)
    v_prm = (/ p_gNa, p_gCaL_LAA, p_gto_LAA, p_gKur, p_gKr_LAA, p_gKs, p_gK1_LA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_LA/)   
  case (Courte_mod_MVR)
    v_prm = (/ p_gNa, p_gCaL_MVR, p_gto, p_gKur, p_gKr_MVR, p_gKs, p_gK1_LA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_LA/)
  case (Courte_mod_BBla)
    v_prm = (/ p_gNa, p_gCaL_BBla, p_gto, p_gKur, p_gKr_BBla, p_gKs, p_gK1_LA, &
               p_gbCa, p_gbNa,p_ACh,p_IKACh_LA/)
  case default
    write (*,10) tcell; stop
  end select
  !.
  return
  10 format ('###.Specified Cell type is not defined:',I3)
end function get_parameter_CRN
!------------------------------------------------------------------------------
subroutine put_param(v_prm,param)
!------------------------------------------------------------------------------
  implicit none
  real(rp), intent(in)    :: v_prm(:)
  type (t_prm)            :: param
  param%p_gNa    = v_prm(1)  ! 1	
  param%p_gCaL   = v_prm(2)  ! 2  
  param%p_gto    = v_prm(3)  ! 3
  param%p_gKur   = v_prm(4)  ! 4
  param%p_gKr    = v_prm(5)  ! 5 
  param%p_gKs    = v_prm(6)  ! 6 
  param%p_gK1    = v_prm(7)  ! 7 
  param%p_gbCa   = v_prm(8)  ! 8
  param%p_gbNa   = v_prm(9)  ! 9 
  param%p_ACh    = v_prm(10) ! 10 
  param%p_IKACh  = v_prm(11) ! 11 
  return
end subroutine put_param
!------------------------------------------------------------------------------
function ic_Courte(ict) result(courte_ic)
!------------------------------------------------------------------------------------
! This function sets the initial conditions described in courtemanche_parameters.inc 
!------------------------------------------------------------------------------------
  implicit none
  !.
  integer (ip), intent(in)  :: ict
  real(rp)                  :: courte_ic(nvar_courte)
  
  select case (ict)
  case (Courte_mod)  ! Anexo Ecuaciones diferenciales lP
    courte_ic = (/ uic, vic, wic, dic, xnOic, xric, Naiic, Kiic, Carelic, oiic, uiic,  &
                   xnISic, xnDpIC3ic, fic, xsic, oaic, uaic, fCaic, Caiic, Caupic,    &
                   xnIC3ic, xnIC2ic, xnIFic, xnC3ic, xnC2ic, xnC1ic, xnDpIC2ic, xnDpIFic, &
                   xnDpC3ic, xnDpC2ic, xnDpC1ic, xnDpOic, xnDpISic, xnDpITic, xnDIC3ic,   &
                   xnDIC2ic, xnDIFic, xnDC3ic, xnDC2ic, xnDC1ic, xnDOic, xnDISic, xnDITic/)
  case(Courte_mod_RA_PM) !.----ATRIA
    courte_ic = (/ uic_RA_PM, vic_RA_PM, wic_RA_PM, dic_RA_PM, xnOic_RA_PM, xric_RA_PM, &
                   Naiic_RA_PM, Kiic_RA_PM, Carelic_RA_PM, oiic_RA_PM, uiic_RA_PM,    &
                   xnISic_RA_PM, xnDpIC3ic_RA_PM, fic_RA_PM, xsic_RA_PM, oaic_RA_PM,  &
                   uaic_RA_PM, fCaic_RA_PM, Caiic_RA_PM, Caupic_RA_PM, xnIC3ic_RA_PM, xnIC2ic_RA_PM, &
                   xnIFic_RA_PM, xnC3ic_RA_PM, xnC2ic_RA_PM, xnC1ic_RA_PM, xnDpIC2ic_RA_PM, xnDpIFic_RA_PM, xnDpC3ic_RA_PM,     &
         	   xnDpC2ic_RA_PM, xnDpC1ic_RA_PM, xnDpOic_RA_PM, xnDpISic_RA_PM,     &
	 	   xnDpITic_RA_PM, xnDIC3ic_RA_PM, xnDIC2ic_RA_PM, xnDIFic_RA_PM,     &
		   xnDC3ic_RA_PM, xnDC2ic_RA_PM, xnDC1ic_RA_PM, xnDOic_RA_PM,         &
		   xnDISic_RA_PM, xnDITic_RA_PM/)
  case(Courte_mod_CT_BBra) !.----ATRIA
    courte_ic = (/ uic_CT_BBra, vic_CT_BBra, wic_CT_BBra, dic_CT_BBra, xnOic_CT_BBra, xric_CT_BBra, &
                   Naiic_CT_BBra, Kiic_CT_BBra, Carelic_CT_BBra, oiic_CT_BBra, uiic_CT_BBra,    &
                   xnISic_CT_BBra, xnDpIC3ic_CT_BBra, fic_CT_BBra, xsic_CT_BBra, oaic_CT_BBra,  &
                   uaic_CT_BBra, fCaic_CT_BBra, Caiic_CT_BBra, Caupic_CT_BBra, xnIC3ic_CT_BBra, xnIC2ic_CT_BBra,  &
                   xnIFic_CT_BBra, xnC3ic_CT_BBra, xnC2ic_CT_BBra, xnC1ic_CT_BBra, xnDpIC2ic_CT_BBra, xnDpIFic_CT_BBra, xnDpC3ic_CT_BBra,       &
	           xnDpC2ic_CT_BBra, xnDpC1ic_CT_BBra, xnDpOic_CT_BBra, xnDpISic_CT_BBra,       &
		   xnDpITic_CT_BBra,  xnDIC3ic_CT_BBra, xnDIC2ic_CT_BBra, xnDIFic_CT_BBra,      &
		   xnDC3ic_CT_BBra,  xnDC2ic_CT_BBra, xnDC1ic_CT_BBra, xnDOic_CT_BBra,          &
		   xnDISic_CT_BBra, xnDITic_CT_BBra/)!
  case(Courte_mod_TVR) !.----ATRIA
    courte_ic = (/ uic_TVR, vic_TVR, wic_TVR, dic_TVR, xnOic_TVR, xric_TVR, Naiic_TVR,   &
                   Kiic_TVR, Carelic_TVR, oiic_TVR, uiic_TVR, xnISic_TVR, xnDpIC3ic_TVR, &
                   fic_TVR, xsic_TVR, oaic_TVR, uaic_TVR, fCaic_TVR, Caiic_TVR,          &
                   Caupic_TVR, xnIC3ic_TVR, xnIC2ic_TVR, xnIFic_TVR, xnC3ic_TVR,         &
                   xnC2ic_TVR, xnC1ic_TVR, xnDpIC2ic_TVR, xnDpIFic_TVR, &
                   xnDpC3ic_TVR, xnDpC2ic_TVR, xnDpC1ic_TVR, xnDpOic_TVR, xnDpISic_TVR, xnDpITic_TVR, xnDIC3ic_TVR,   &
                   xnDIC2ic_TVR, xnDIFic_TVR, xnDC3ic_TVR, xnDC2ic_TVR,    &
                   xnDC1ic_TVR, xnDOic_TVR, xnDISic_TVR, xnDITic_TVR/)!
  case(Courte_mod_RAA) !.----ATRIA
    courte_ic = (/ uic_RAA, vic_RAA, wic_RAA, dic_RAA, xnOic_RAA, xric_RAA, Naiic_RAA,    &
                   Kiic_RAA, Carelic_RAA, oiic_RAA, uiic_RAA, xnISic_RAA, xnDpIC3ic_RAA,  &
                   fic_RAA, xsic_RAA, oaic_RAA, uaic_RAA, fCaic_RAA, Caiic_RAA,           &
                   Caupic_RAA, xnIC3ic_RAA, xnIC2ic_RAA, xnIFic_RAA, xnC3ic_RAA,          &
                   xnC2ic_RAA, xnC1ic_RAA, xnDpIC2ic_RAA, xnDpIFic_RAA, &
                   xnDpC3ic_RAA, xnDpC2ic_RAA, xnDpC1ic_RAA, xnDpOic_RAA, xnDpISic_RAA, xnDpITic_RAA, xnDIC3ic_RAA,   &
                   xnDIC2ic_RAA, xnDIFic_RAA, xnDC3ic_RAA, xnDC2ic_RAA,     & 
		   xnDC1ic_RAA, xnDOic_RAA, xnDISic_RAA, xnDITic_RAA/)
  case(Courte_mod_LA) !.----ATRIA
    courte_ic = (/ uic_LA, vic_LA, wic_LA, dic_LA, xnOic_LA, xric_LA, Naiic_LA, Kiic_LA,   &
                   Carelic_LA, oiic_LA, uiic_LA, xnISic_LA, xnDpIC3ic_LA, fic_LA, xsic_LA, &
                   oaic_LA, uaic_LA, fCaic_LA, Caiic_LA, Caupic_LA, xnIC3ic_LA, xnIC2ic_LA, &
                   xnIFic_LA, xnC3ic_LA, xnC2ic_LA, xnC1ic_LA, xnDpIC2ic_LA, xnDpIFic_LA,   &
                   xnDpC3ic_LA, xnDpC2ic_LA, xnDpC1ic_LA, xnDpOic_LA, xnDpISic_LA, xnDpITic_LA, xnDIC3ic_LA,   &
                   xnDIC2ic_LA, xnDIFic_LA, xnDC3ic_LA, xnDC2ic_LA, xnDC1ic_LA, xnDOic_LA, xnDISic_LA, xnDITic_LA/)
  case(Courte_mod_PV) !.----ATRIA
    courte_ic = (/ uic_PV, vic_PV, wic_PV, dic_PV, xnOic_PV, xric_PV, Naiic_PV, Kiic_PV,     &
                   Carelic_PV, oiic_PV, uiic_PV, xnISic_PV, xnDpIC3ic_PV, fic_PV, xsic_PV,   &
                   oaic_PV, uaic_PV, fCaic_PV, Caiic_PV, Caupic_PV, xnIC3ic_PV, xnIC2ic_PV,  &
                   xnIFic_PV, xnC3ic_PV, xnC2ic_PV, xnC1ic_PV, xnDpIC2ic_PV, xnDpIFic_PV,   &
                   xnDpC3ic_PV, xnDpC2ic_PV, xnDpC1ic_PV, xnDpOic_PV, xnDpISic_PV, xnDpITic_PV, xnDIC3ic_PV,   &
                   xnDIC2ic_PV, xnDIFic_PV, xnDC3ic_PV, xnDC2ic_PV, xnDC1ic_PV, xnDOic_PV, xnDISic_PV, xnDITic_PV/)
  case(Courte_mod_LAA) !.----ATRIA
    courte_ic = (/ uic_LAA, vic_LAA, wic_LAA, dic_LAA, xnOic_LAA, xric_LAA, Naiic_LAA,   &
                   Kiic_LAA, Carelic_LAA, oiic_LAA, uiic_LAA, xnISic_LAA, xnDpIC3ic_LAA, &
                   fic_LAA, xsic_LAA, oaic_LAA, uaic_LAA, fCaic_LAA, Caiic_LAA,          &
                   Caupic_LAA, xnIC3ic_LAA, xnIC2ic_LAA, xnIFic_LAA, xnC3ic_LAA, xnC2ic_LAA, & 
                   xnC1ic_LAA, xnDpIC2ic_LAA, xnDpIFic_LAA,   &
                   xnDpC3ic_LAA, xnDpC2ic_LAA, xnDpC1ic_LAA, xnDpOic_LAA, xnDpISic_LAA, xnDpITic_LAA, xnDIC3ic_LAA,   &
                   xnDIC2ic_LAA, xnDIFic_LAA, xnDC3ic_LAA, xnDC2ic_LAA, xnDC1ic_LAA, xnDOic_LAA, xnDISic_LAA, xnDITic_LAA/)
  case(Courte_mod_MVR) !.----ATRIA
    courte_ic = (/ uic_MVR, vic_MVR, wic_MVR, dic_MVR, xnOic_MVR, xric_MVR, Naiic_MVR,     &
                   Kiic_MVR, Carelic_MVR, oiic_MVR, uiic_MVR, xnISic_MVR, xnDpIC3ic_MVR,   &
                   fic_MVR, xsic_MVR, oaic_MVR, uaic_MVR, fCaic_MVR, Caiic_MVR,            &
                   Caupic_MVR, xnIC3ic_MVR, xnIC2ic_MVR, xnIFic_MVR, xnC3ic_MVR, xnC2ic_MVR, xnC1ic_MVR,  xnDpIC2ic_MVR, xnDpIFic_MVR, &
                   xnDpC3ic_MVR, xnDpC2ic_MVR, xnDpC1ic_MVR, xnDpOic_MVR, xnDpISic_MVR, xnDpITic_MVR, xnDIC3ic_MVR,   &
                   xnDIC2ic_MVR, xnDIFic_MVR, xnDC3ic_MVR, xnDC2ic_MVR, xnDC1ic_MVR, xnDOic_MVR, xnDISic_MVR, xnDITic_MVR/)
  case(Courte_mod_BBla) !.----ATRIA
    courte_ic = (/ uic_BBla, vic_BBla, wic_BBla, dic_BBla, xnOic_BBla, xric_BBla, Naiic_BBla,  &
                   Kiic_BBla, Carelic_BBla, oiic_BBla, uiic_BBla, xnISic_BBla, xnDpIC3ic_BBla, &
                   fic_BBla, xsic_BBla, oaic_BBla, uaic_BBla, fCaic_BBla, Caiic_BBla,          &
                   Caupic_BBla, xnIC3ic_BBla, xnIC2ic_BBla, xnIFic_BBla, xnC3ic_BBla,          &
                   xnC2ic_BBla, xnC1ic_BBla,  xnDpIC2ic_BBla, xnDpIFic_BBla, &
                   xnDpC3ic_BBla, xnDpC2ic_BBla, xnDpC1ic_BBla, xnDpOic_BBla, xnDpISic_BBla, xnDpITic_BBla, xnDIC3ic_BBla,   &
                   xnDIC2ic_BBla, xnDIFic_BBla, xnDC3ic_BBla, xnDC2ic_BBla, xnDC1ic_BBla, xnDOic_BBla, xnDISic_BBla, xnDITic_BBla/)
  case default
    courte_ic = (/ uic, vic, wic, dic, xnOic, xric, Naiic, Kiic, Carelic, oiic, uiic,       &
                   xnISic, xnDpIC3ic, fic, xsic, oaic, uaic, fCaic, Caiic, Caupic, xnIC3ic, &
                   xnIC2ic, xnIFic, xnC3ic, xnC2ic, xnC1ic, xnDpIC2ic, xnDpIFic, &
                   xnDpC3ic, xnDpC2ic, xnDpC1ic, xnDpOic, xnDpISic, xnDpITic, xnDIC3ic,   &
                   xnDIC2ic, xnDIFic, xnDC3ic, xnDC2ic, xnDC1ic, xnDOic, xnDISic, xnDITic/)
  end select
  !.
  return
end function ic_Courte

!------------------------------------------------------------------------------
subroutine get_me_struct(str_me,v_me)
!------------------------------------------------------------------------------
  implicit none
  type (t_courte), intent(in) :: str_me
  real(rp), intent(out)       :: v_me(nvar_courte)
  
  
  v_me(1)  =   str_me%u 
  v_me(2)  =   str_me%v 
  v_me(3)  =   str_me%w
  v_me(4)  =   str_me%d
  v_me(5)  =   str_me%xnO
  v_me(6)  =   str_me%xr
  v_me(7)  =   str_me%Nai
  v_me(8)  =   str_me%Ki
  v_me(9)  =   str_me%Carel
  v_me(10) =   str_me%oi
  v_me(11) =   str_me%ui
  v_me(12) =   str_me%xnIS
  v_me(13) =   str_me%xnDpIC3
  v_me(14) =   str_me%f
  v_me(15) =   str_me%xs
  v_me(16) =   str_me%oa
  v_me(17) =   str_me%ua
  v_me(18) =   str_me%fCa
  v_me(19) =   str_me%Cai
  v_me(20) =   str_me%Caup
  v_me(21) =   str_me%xnIC3
  v_me(22) =   str_me%xnIC2
  v_me(23) =   str_me%xnIF
  v_me(24) =   str_me%xnC3
  v_me(25) =   str_me%xnC2
  v_me(26) =   str_me%xnC1
  v_me(27) =   str_me%xnDpIC2
  v_me(28) =   str_me%xnDpIF
  v_me(29) =   str_me%xnDpC3
  v_me(30) =   str_me%xnDpC2
  v_me(31) =   str_me%xnDpC1
  v_me(32) =   str_me%xnDpO
  v_me(33) =   str_me%xnDpIS
  v_me(34) =   str_me%xnDpIT
  v_me(35) =   str_me%xnDIC3
  v_me(36) =   str_me%xnDIC2
  v_me(37) =   str_me%xnDIF
  v_me(38) =   str_me%xnDC3
  v_me(39) =   str_me%xnDC2
  v_me(40) =   str_me%xnDC1
  v_me(41) =   str_me%xnDO
  v_me(42) =   str_me%xnDIS
  v_me(43) =   str_me%xnDIT
   
  return
end subroutine get_me_struct

!------------------------------------------------------------------------------
subroutine put_me_struct(v_me,str_me)
!------------------------------------------------------------------------------
  implicit none
  real(rp), intent(in)         :: v_me(:)
  type (t_courte), intent(out) :: str_me
  str_me%u       = v_me(1)
  str_me%v       = v_me(2)
  str_me%w       = v_me(3)
  str_me%d       = v_me(4)
  str_me%xnO     = v_me(5)
  str_me%xr      = v_me(6)
  str_me%Nai     = v_me(7)
  str_me%Ki      = v_me(8)
  str_me%Carel   = v_me(9)
  str_me%oi      = v_me(10)
  str_me%ui      = v_me(11)
  str_me%xnIS    = v_me(12)
  str_me%xnDpIC3 = v_me(13)
  str_me%f       = v_me(14)
  str_me%xs      = v_me(15)
  str_me%oa      = v_me(16)
  str_me%ua      = v_me(17)
  str_me%fCa     = v_me(18)
  str_me%Cai     = v_me(19)
  str_me%Caup    = v_me(20)
  str_me%xnIC3   = v_me(21)
  str_me%xnIC2   = v_me(22)
  str_me%xnIF    = v_me(23)
  str_me%xnC3    = v_me(24)
  str_me%xnC2    = v_me(25)
  str_me%xnC1    = v_me(26)
  str_me%xnDpIC2 = v_me(27)
  str_me%xnDpIF  = v_me(28)
  str_me%xnDpC3  = v_me(29)
  str_me%xnDpC2  = v_me(30)
  str_me%xnDpC1  = v_me(31)
  str_me%xnDpO   = v_me(32)
  str_me%xnDpIS  = v_me(33)
  str_me%xnDpIT  = v_me(34)
  str_me%xnDIC3  = v_me(35)
  str_me%xnDIC2  = v_me(36)
  str_me%xnDIF   = v_me(37)
  str_me%xnDC3   = v_me(38)
  str_me%xnDC2   = v_me(39)
  str_me%xnDC1   = v_me(40)
  str_me%xnDO    = v_me(41)
  str_me%xnDIS   = v_me(42)
  str_me%xnDIT   = v_me(43)
  
   
return
end subroutine put_me_struct
! ======================================================================================
! ------------------------------------------------------------------------------------- 
! ===================================================================================== 

subroutine f_currents(Vm,prm,cur,cm_m,Iion)
! ===================================================================================== 
! ===================================================================================== 
  implicit none
  real(rp),intent(in)           :: Vm
  real(rp),intent(out)          :: Iion
  type(t_courte), intent(in)    :: cm_m
  type (t_prm)                  :: prm
  type (t_cur), intent(out)     :: cur
  real(rp)                      :: ECa,ENa,EK,sigma,fNaK,tau_tr
  real(rp)                      :: den, num, factor_ICaL, factor_IK1, factor_IKr, factor_IKs, &  ! Factors for drug simulation
								factor_IKur, factor_Ito, factor_INa, factor_INaL, factor_INaK, &
								factor_INCX, factor_IKACh
  
  
  ! Factors to PoM - Violeta 12/05/2023
  
    factor_ICaL = 1.0;
    factor_IK1 = 1.0;
    factor_IKr = 1.0;
    factor_IKs = 1.0;
    factor_IKur = 1.0;
    factor_Ito = 1.0;
    factor_INa = 1.0;
    factor_INaK = 1.0;
    factor_INCX = 1.0;
    factor_INaL = 1;
    factor_IKACh = 1.0;


    
 ! Nerst potentials
  ECa = (p_R*p_T/(2.0*p_F))*log(p_Cao/cm_m%Cai);
  ENa = (p_R*p_T/p_F)*log(p_Nao/cm_m%Nai);
  EK = (p_R*p_T/p_F)*log(p_Ko/cm_m%Ki);

  ! L type Ca channel

!  cur%ICaL = p_Cm*prm%p_gCaL*cm_m%d*cm_m%f*cm_m%fCa*(Vm-65.0);
  cur%ICaL = factor_ICaL*prm%p_gCaL*cm_m%d*cm_m%f*cm_m%fCa*(Vm-65.0);

  ! NaCa exchanger current

  den = (((p_KmNa**3)+(p_Nao**3))*(p_KmCa+p_Cao)*(1.0+p_ksat*exp(((p_gamma-1.0)*p_F*Vm)/(p_R*p_T))));
  num = exp((p_gamma*p_F*Vm)/(p_R*p_T))*(cm_m%Nai**3)*p_Cao
  num = num - exp(((p_gamma-1.0)*p_F*Vm)/(p_R*p_T))*(p_Nao**3)*cm_m%Cai
!  cur%INaCa = p_Cm*p_INaCamax*num/den
  cur%INaCa = factor_INCX*p_INaCamax*num/den

  ! Background currents

!  cur%IbNa = p_Cm*prm%p_gbNa*(Vm-ENa);
!  cur%IbCa = p_Cm*prm%p_gbCa*(Vm-ECa);
  cur%IbNa = prm%p_gbNa*(Vm-ENa);
  cur%IbCa = prm%p_gbCa*(Vm-ECa);
 
  
  
  ! Fast sodium current
  
!  cur%INa = p_Cm*prm%p_gNa*(cm_m%m**3)*cm_m%h*cm_m%j*(Vm-ENa);
 ! cur%INa = factor_INa*prm%p_gNa*(cm_m%m**3)*cm_m%h*cm_m%j*(Vm-ENa);
  cur%INa = factor_INa*prm%p_gNa*(cm_m%xnO)*(Vm-ENa);  
  ! Rapid delayed rectifier K current
  
!  cur%IKr = p_Cm*prm%p_gKr*cm_m%xr*(Vm-EK)/(1.0+exp((Vm+15.0)/22.4));
  cur%IKr = factor_IKr*prm%p_gKr*cm_m%xr*(Vm-EK)/(1.0+exp((Vm+15.0)/22.4));
  
  ! Sarcolemmal Ca pump current
  
!  cur%IpCa = p_Cm*p_IpCamax*cm_m%Cai/(0.0005+cm_m%Cai);
  cur%IpCa = p_IpCamax*cm_m%Cai/(0.0005+cm_m%Cai);
  
  ! Slow delayed rectifier K current
  
!  cur%IKs = p_Cm*prm%p_gKs*(cm_m%xs**2)*(Vm-EK);
  cur%IKs = factor_IKs*prm%p_gKs*(cm_m%xs**2)*(Vm-EK);
  
  ! Sodium potassium pump
  
  sigma = (1.0/7.0)*(exp(p_Nao/67.3)-1.0);
  fNaK = 1.0/(1.0+0.1245*exp((-0.1*p_F*Vm)/(p_R*p_T))+0.0365*sigma*exp((-p_F*Vm)/(p_R*p_T)));
!  cur%INaK = (p_Ko*p_Cm*p_INaKmax*fNaK/(1.0+((p_KmNai/cm_m%Nai)**1.5)))/(p_Ko+p_KmKo);
  cur%INaK = factor_INaK*(p_Ko*p_INaKmax*fNaK/(1.0+((p_KmNai/cm_m%Nai)**1.5)))/(p_Ko+p_KmKo);
  
  ! Time independent potassium current
  
!  cur%IK1 = p_Cm*prm%p_gK1*(Vm-EK)/(1.0+exp(0.07*(Vm+80.0)));
  cur%IK1 = factor_IK1*prm%p_gK1*(Vm-EK)/(1.0+exp(0.07*(Vm+80.0)));
  
  ! Transfer current from NSR to JSR
  
  tau_tr = 180.0;
  cur%Itr = (cm_m%Caup-cm_m%Carel)/tau_tr;
  
  ! Transient outward K current
  
!  Ito = p_Cm*prm%p_gto*(cm_m%oa**3)*cm_m%oi*(Vm-EK);
  cur%Ito = factor_Ito*prm%p_gto*(cm_m%oa**3)*cm_m%oi*(Vm-EK);
  
  ! Ultrarapid delayed rectifier K current
  
  prm%p_gKur = 0.005+0.05/(1.0+exp((Vm-15.0)/(-13.0)));
!  IKur = p_Cm*prm%p_gKur*(cm_m%ua**3)*cm_m%ui*(Vm-EK);
  cur%IKur = factor_IKur*prm%p_gKur*(cm_m%ua**3)*cm_m%ui*(Vm-EK);
  
! Inward rectifier IKACh (Grandi) 

	cur%IKACh = factor_IKACh*prm%p_IKACh*(((0.08+(0.04/(1.0+exp((Vm+91.0)/12.0))))*(Vm-EK))/(1.0+(0.03/prm%p_ACh)**2.1));

!  Iion = (INa+IK1+Ito+IKur+IKr+IKs+IbNa+IbCa+INaK+IpCa+INaCa+ICaL)/p_Cm;
  Iion = cur%INa+cur%IK1+cur%Ito+cur%IKur+cur%IKr+cur%IKs+ &
         cur%IbNa+ cur%IbCa+ cur%INaK+ cur%IpCa+ cur%INaCa+&
         cur%ICaL+cur%IKACh
  
  return
end subroutine f_currents
!-----------------------------------------------------------------------------------------------------------------

subroutine f_concentrations(dt,cur,cm_m,Istm)
! ===================================================================================== 
! ===================================================================================== 
  implicit none
  real(rp),intent(in)            :: dt,Istm
  type(t_courte), intent(inout)  :: cm_m
  type (t_cur), intent(inout)    :: cur
  real(rp)                       :: Ca_Cmdn,Ca_Trpn,Ca_Csqn,B1,B2,factor_Irel,factor_Ileak,factor_Iup

  
  ! Factors to PoM - Violeta 12/05/2023
  
    factor_Irel = 1;
    factor_Ileak = 1;
    factor_Iup = 1;
  
  
  ! Ca buffers
  
  Ca_Cmdn = (p_CMDNmax*cm_m%Cai)/(cm_m%Cai+p_KmCmdn);
  Ca_Trpn = (p_TRPNmax*cm_m%Cai)/(cm_m%Cai+p_KmTrpn);
  Ca_Csqn = (p_CSQNmax*cm_m%Carel)/(cm_m%Carel+p_KmCsqn);
  
  ! Ca leak current by the NSR
  
  cur%Ileak = factor_Ileak*p_Iupmax*cm_m%Caup/p_Caupmax;
  
  ! Ca uptake current by the NSR
  
  cur%Iup = factor_Iup*p_Iupmax/(1.0+p_Kup/cm_m%Cai);

  ! Ca release current from JSR
  if ((p_krel*(cm_m%u**2)*cm_m%v*cm_m%w*(cm_m%Carel-cm_m%Cai))<1.0E-25) then 
	cur%Irel = 0.0;
  else
  	cur%Irel = factor_Irel*p_krel*(cm_m%u**2)*cm_m%v*cm_m%w*(cm_m%Carel-cm_m%Cai);
  end if

  ! Intracellular ion concentrations
  
  cm_m%Nai = cm_m%Nai + dt*p_Cm*(-3.0*cur%INaK-(3.0*cur%INaCa+cur%IbNa+cur%INa))/(p_Vi*p_F);
  !
  cm_m%Ki = cm_m%Ki + dt*p_Cm*(2.0*cur%INaK-(Istm+cur%IK1+cur%Ito+cur%IKur+cur%IKr+cur%IKs+cur%IKACh))/(p_Vi*p_F);
  !
  B1 = p_Cm*(2.0*cur%INaCa-(cur%IpCa+cur%ICaL+cur%IbCa))/(2.0*p_Vi*p_F) + (p_Vup*(cur%Ileak-cur%Iup)+cur%Irel*p_Vrel)/p_Vi;
  B2 = 1.0 + (p_TRPNmax*p_KmTrpn)/((cm_m%Cai+p_KmTrpn)**2) + (p_CMDNmax*p_KmCmdn)/((cm_m%Cai+p_KmCmdn)**2);
  cm_m%Cai = cm_m%Cai + dt*B1/B2;
  !
  cm_m%Caup = cm_m%Caup + dt*(cur%Iup-(cur%Ileak+cur%Itr*p_Vrel/p_Vup));
  !
  cm_m%Carel = cm_m%Carel + dt*(cur%Itr-cur%Irel)/(1.0+p_CSQNmax*p_KmCsqn/((cm_m%Carel+p_KmCsqn)**2));  
  return
end subroutine f_concentrations
!-----------------------------------------------------------------------------------------------------------------

!=================================================================================

!=================================================================================
subroutine f_gates(dt, Vm, cur, cm_m)
! ===================================================================================== 
! SUBRUTINA DE COMPUERTAS:
!   - Canales de Ca y K: Hodgkin-Huxley original Courtemanche
!   - Canal de Na: Modelo Markov 23-estados con vernakalant (p_verna_conc_nM)
! ===================================================================================== 
  implicit none
  real(rp), intent(in)           :: dt, Vm
  type(t_courte), intent(inout)  :: cm_m
  type(t_cur),    intent(in)     :: cur

  ! ---------------- HH: variables auxiliares ----------------
  real(rp) :: Fn
  real(rp) :: u_inf, tau_u, v_inf, tau_v, w_inf, tau_w
  real(rp) :: d_inf, tau_d, fCa_inf, tau_fCa, f_inf, tau_f
  real(rp) :: xr_inf, alfa_xr, beta_xr, tau_xr
  real(rp) :: xs_inf, alfa_xs, beta_xs, tau_xs
  real(rp) :: oa_inf, alfa_oa, beta_oa, tau_oa
  real(rp) :: oi_inf, alfa_oi, beta_oi, tau_oi
  real(rp) :: ua_inf, alfa_ua, beta_ua, tau_ua
  real(rp) :: ui_inf, alfa_ui, beta_ui, tau_ui

  ! ---------------- MARKOV Na: parmetros ----------------
  ! Parmetros x1..x25
  real(rp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25
  real(rp) :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14
  real(rp) :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14
  real(rp) :: actshift,h1
  real(rp) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
  real(rp) :: p14_new,p15_new,p16_new,p17_new,p18,p19,p20,p21,p01
  real(rp) :: p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35

  real(rp) :: fleca_fact, verna_fact, D_verna, conc, pKa, kd_open, pH, portion
  real(rp) :: conc_dplus, conc_d, diffusion
  real(rp) :: kon,kcon,koff,kcoff,k_on,k_off,ki_on,ki_off,kc_on,kc_off
  real(rp) :: Tfactor
  real(rp) :: a11,a12,a13,b11,b12,b13
  real(rp) :: a3_ss,a3_tau,a3,b3,a2,b2,ax,bx
  real(rp) :: a13c,b13c,a13n,b13n,ax1,bx1,ax2,bx2
  real(rp) :: a22,a_22,b22,b_22,a33,b33,a_33,b_33,a44,b44,a_44,b_44

  ! Derivadas Markov (para integracin explcita con dt_sub)
  real(rp) :: dIC3,dIC2,dIF,dC3,dC2,dC1,dO_mkv,dIS
  real(rp) :: dDpIC3,dDpIC2,dDpIF,dDpC3,dDpC2,dDpC1,dDpO,dDpIS,dDpIT
  real(rp) :: dDIC3,dDIC2,dDIF,dDC3,dDC2,dDC1,dDO,dDIS,dDIT

  ! Estados previos para actualización semi-implícita
  real(rp) :: oIC3,oIC2,oIF,oC3,oC2,oC1,oO,oIS
  real(rp) :: oDpIC3,oDpIC2,oDpIF,oDpC3,oDpC2,oDpC1,oDpO,oDpIS,oDpIT
  real(rp) :: oDIC3,oDIC2,oDIF,oDC3,oDC2,oDC1,oDO,oDIS,oDIT

  ! Términos de pérdida lineal (dX = P - L*X)
  real(rp) :: lIC3,lIC2,lIF,lC3,lC2,lC1,lO,lIS
  real(rp) :: lDpIC3,lDpIC2,lDpIF,lDpC3,lDpC2,lDpC1,lDpO,lDpIS,lDpIT
  real(rp) :: lDIC3,lDIC2,lDIF,lDC3,lDC2,lDC1,lDO,lDIS,lDIT

  ! Sub-stepping
  integer(ip) :: k_sub, n_sub
  real(rp)    :: dt_sub
  real(rp)    :: max_rate


  ! Normalizacin
  real(rp) :: sum_states, norm_factor, eps

  ! =====================================================================================
  ! 1. HodgkinHuxley Ca y K (Courtemanche original)
  ! =====================================================================================

  ! Fn segn Courtemanche original:
  Fn = 1.0e3_rp*( 1.0e-15_rp*p_Vrel*cur%Irel - p_Cm*( (1.0e-15_rp/(2.0_rp*p_F))*(0.5_rp*cur%ICaL - 0.2_rp*cur%INaCa) ) )

  ! ----- u, v (Ca release gate SR) -----
  tau_u = 8.0_rp
  if (1.0_rp/(1.0_rp+exp(-(Fn-3.4175e-13_rp)/13.67e-16_rp)) < 1.0e-25_rp) then
     u_inf = 0.0_rp
  else
     u_inf = 1.0_rp/(1.0_rp+exp(-(Fn-3.4175e-13_rp)/13.67e-16_rp))
  end if
  cm_m%u = u_inf - (u_inf-cm_m%u)*exp(-dt/tau_u)

  tau_v = 1.91_rp + 2.09_rp/(1.0_rp+exp(-(Fn-3.4175e-13_rp)/13.67e-16_rp))
  v_inf = 1.0_rp - 1.0_rp/(1.0_rp+exp(-(Fn-6.835e-14_rp)/13.67e-16_rp))
  cm_m%v = v_inf - (v_inf-cm_m%v)*exp(-dt/tau_v)

  ! ----- w gate ICaL -----
  if (abs(Vm-7.9_rp) < 1.0e-10_rp) then
     tau_w = 6.0_rp*0.2_rp/1.3_rp
  else
     tau_w = 6.0_rp*(1.0_rp-exp(-(Vm-7.9_rp)/5.0_rp))/ &
             ((1.0_rp+0.3_rp*exp(-(Vm-7.9_rp)/5.0_rp))*(Vm-7.9_rp))
  end if
  w_inf = 1.0_rp - 1.0_rp/(1.0_rp+exp(-(Vm-40.0_rp)/17.0_rp))
  cm_m%w = w_inf - (w_inf-cm_m%w)*exp(-dt/tau_w)

  ! ----- d gate ICaL -----
  d_inf = 1.0_rp/(1.0_rp+exp((Vm+10.0_rp)/(-8.0_rp)))
  if (abs(Vm+10.0_rp) < 1.0e-10_rp) then
     tau_d = 4.579_rp/(1.0_rp+exp((Vm+10.0_rp)/(-6.24_rp)))
  else
     tau_d = (1.0_rp-exp((Vm+10.0_rp)/(-6.24_rp)))/ &
             (0.035_rp*(Vm+10.0_rp)*(1.0_rp+exp((Vm+10.0_rp)/(-6.24_rp))))
  end if
  cm_m%d = d_inf - (d_inf-cm_m%d)*exp(-dt/tau_d)

  ! ----- fCa -----
  fCa_inf  = 1.0_rp/(1.0_rp + cm_m%Cai/0.00035_rp)
  tau_fCa  = 2.0_rp
  cm_m%fCa = fCa_inf - (fCa_inf-cm_m%fCa)*exp(-dt/tau_fCa)

  ! ----- f -----
  f_inf = (1.0_rp + exp((Vm+28.0_rp)/6.9_rp))**(-1)
  tau_f = 9.0_rp/(0.0197_rp*exp(-(0.0337_rp**2)*(Vm+10.0_rp)**2) + 0.02_rp)
  cm_m%f = f_inf - (f_inf-cm_m%f)*exp(-dt/tau_f)

  ! =====================================================================================
  ! 2. Modelo Markov de INa con Vernakalant (23 estados)
  ! =====================================================================================

  eps = 1.0e-15_rp

  ! --- Parmetros x1..x25 (desde tu MATLAB) ---
  x1  = 1.1821516171073378_rp
  x2  = 1.8128950937709263_rp
  x3  = 0.55559223164198501_rp
  x4  = 1.2653811014166259_rp
  x5  = 0.379465955877468_rp
  x6  = 2.5485320500798077_rp
  x7  = 1.5828527326859647_rp
  x8  = 1.0237930757107652_rp
  x9  = 0.63628656115270266_rp
  x10 = 1.0802804689922176_rp
  x11 = 1.0992577563248505_rp
  x12 = 1.1001781318138173_rp
  x13 = 1.0751649122176028_rp
  x14 = 1.0400623701712255_rp*0.1_rp
  x15 = 1.0221615376561712_rp
  x16 = 1.0349852098209555_rp
  x17 = 1.0095837445899951_rp
  x18 = 1.1952465_rp
  x19 = 1.811411333_rp*1.05_rp
  x20 = 0.8_rp
  x21 = 0.9994445714285713_rp
  x22 = 1.0284021_rp
  x23 = 0.98904877_rp
  x24 = 1.01288225_rp
  x25 = 1.00647894_rp

  ! Flecainida y Vernakalant
  fleca_fact = 0.0_rp
  verna_fact = 0.0_rp   ! aqu representamos Vernakalant

  ! Factores Flecainida (no usados pero mantenidos para compatibilidad)
  f1  = 15802.517082265076_rp
  f2  = 11212.49445307003_rp
  f3  = 1.0835911457972944_rp
  f4  = 4.7213138217409707_rp
  f5  = 1.1886854610200599_rp
  f6  = 1.1293915676022106_rp
  f7  = 0.98647121091698153_rp
  f8  = 0.54542349779755606_rp
  f9  = 1.020713680006784_rp
  f10 = 2.6609951303941308_rp
  f11 = 1.1411905523331658_rp
  f12 = 5.3540204782968566_rp
  f13 = 1.028853996176295_rp
  f14 = 1.0033797797467918_rp

  ! Factores Vernakalant
  v1  = 7170.4622644564224_rp
  v2  = 15352.068978379524_rp
  v3  = 0.78329295519833497_rp
  v4  = 6.2271179964787571_rp
  v5  = 22.10184390723893_rp
  v6  = 1.5107146709941404_rp
  v7  = 1.0935642293370962_rp
  v8  = 0.65910138578597266_rp
  v9  = 2.3700350948377_rp
  v10 = 6.2071593539114822_rp
  v11 = 0.0_rp
  v12 = 0.0_rp
  v13 = 0.0_rp
  v14 = 0.0_rp

  ! Parmetros Markov base
  actshift = -15.0_rp*x19
  h1       =  2.0_rp*x20

  p1 = 8.5539_rp*x1
  p2 = 7.4392e-2_rp*x2
  p3 = 17.0_rp*x3
  p4 = 15.0_rp*x4
  p5 = 12.0_rp*x5
  p6 = 2.0373e-1_rp*x6
  p7 = 150.0_rp*x7
  p8 = 7.5215e-2_rp*x8
  p9 = 20.3_rp*x9
  p10= 2.7574_rp*x10
  p11= 5.0_rp*x11
  p12= 4.7755e-1_rp*x12
  p13= 10.0_rp*x13

  p14_new = -70.0_rp*x21
  p15_new =  3.5_rp*x22
  p16_new = 0.052_rp*2.9_rp*x23
  p17_new = 0.132_rp*1.9_rp*x24
  p18 = 13.370_rp*x14
  p19 = 43.749_rp*x15
  p20 = 3.4229e-2_rp*x16
  p21 = 1.7898e-2_rp*x17
  p01 = 41.0_rp*x25

  p22 = (fleca_fact*3.6324e-3_rp*f1)   + (verna_fact*5.6974e-03_rp*v1)
  p23 = (fleca_fact*2.6452_rp*f2)      + (verna_fact*8.4559e+01_rp*v2)
  p24 = (fleca_fact*5.7831e-5_rp*f3)   + (verna_fact*6.3992e-07_rp*v3)
  p25 = (fleca_fact*1.6689e-8_rp*f4)   + (verna_fact*1.3511e+00_rp*v4)
  p26 = (fleca_fact*2.6126e-01_rp*f5)  + (verna_fact*1.3110e-01_rp*v5)
  p27 = (fleca_fact*1.4847e3_rp*f6)    + (verna_fact*6.7067e-06_rp*v6)
  p28 = (fleca_fact*4.2385e+01_rp*f7)  + (verna_fact*1.7084e-05_rp*v7)
  p29 = (fleca_fact*1.7352e-6_rp*f8)   + (verna_fact*1.9698e-05_rp*v8)
  p30 = (fleca_fact*2.1181e+00_rp*f9)  + (verna_fact*4.8477_rp*v9)
  p31 = (fleca_fact*6.7505e-05_rp*f10) + (verna_fact*3.2976_rp*v10)
  p32 = (fleca_fact*2.4135_rp*f11)     + (verna_fact*2.4135_rp*v11)
  p33 = (fleca_fact*4.9001e-2_rp*f12)  + (verna_fact*4.9001e-2_rp*v12)
  p34 = (fleca_fact*1.0326e-03_rp*f13) + (verna_fact*1.0326e-03_rp*v13)
  p35 = (fleca_fact*2.1378e-02_rp*f14) + (verna_fact*2.1378e-02_rp*v14)

  ! --- Concentracin de vernakalant global (p_verna_conc_nM) ---
  D_verna = p_verna_conc_nM          ! nM
  conc    = D_verna * 1.0e-9_rp       ! mol/L, igual que MATLAB

  pKa    = (fleca_fact*9.3_rp) + (verna_fact*5.4_rp)
  kd_open= ((fleca_fact*11.2e-6_rp) + (verna_fact*318e-6_rp))* &
           exp(-0.7_rp*Vm*p_F/(p_R*p_T))

  pH = 7.4_rp
  portion    = 1.0_rp/(1.0_rp+10.0_rp**(pH-pKa))
  conc_dplus = portion*conc
  conc_d     = (1.0_rp-portion)*conc

  diffusion = (fleca_fact*5500.0_rp) + (verna_fact*500.0_rp)

  kon   = conc_dplus * diffusion
  kcon  = kon
  koff  = kd_open * diffusion
  kcoff = koff

  k_on  = conc_d * diffusion
  k_off = ((fleca_fact*400e-6_rp) + (verna_fact*400e-6_rp))*diffusion
  ki_on = k_on/2.0_rp
  ki_off= ((fleca_fact*5.4e-6_rp) + (verna_fact*3.4e-6_rp))*diffusion
  kc_on = k_on/2.0_rp
  kc_off= ((fleca_fact*800e-6_rp) + (verna_fact*900e-6_rp))*diffusion

  Tfactor = 1.0_rp/(3.0_rp**((37.0_rp-(p_T-273.0_rp))/10.0_rp))

  a11 = Tfactor*p1/(p2*safe_exp(-(Vm-actshift)/p3) + p6*safe_exp(-(Vm-actshift)/p7))
  a12 = Tfactor*p1/(p2*safe_exp(-(Vm-actshift)/p4) + p6*safe_exp(-(Vm-actshift)/p7))
  a13 = Tfactor*p1/(p2*safe_exp(-(Vm-actshift)/p5) + p6*safe_exp(-(Vm-actshift)/p7))
  b11 = Tfactor*p8 * safe_exp(-(Vm-actshift)/p9)
  b12 = Tfactor*p10* safe_exp(-(Vm-actshift-p11)/p9)
  b13 = Tfactor*p12* safe_exp(-(Vm-actshift-p13)/p9)

  a3_ss = 1.0_rp/(1.0_rp+safe_exp((Vm-p14_new)/p15_new))
  a3_tau= h1 + p01*safe_exp(p16_new*(Vm-p14_new))/(1.0_rp+safe_exp(p17_new*(Vm-p14_new)))
  a3    = Tfactor*a3_ss/a3_tau
  b3    = Tfactor*(1.0_rp-a3_ss)/a3_tau

  a2 = Tfactor*p18*safe_exp(Vm/p19)
  b2 = (a13*a2*a3)/(b13*b3)
  ax = p20*a2
  bx = p21*a3

  a13c = p22*a13
  if (kon > 0.0_rp) then
    b13c = (b13*kcon*koff*a13c)/(kon*kcoff*a13)
  else
    b13c = 0.0_rp
  end if

  a13n = p23*a13
  if (k_on > 0.0_rp) then
    b13n = (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on)
  else
    b13n = 0.0_rp
  end if

  ax1 = p24*ax
  bx1 = p25*bx
  ax2 = p26*ax
  if (ki_on > 0.0_rp) then
    bx2 = (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off)
  else
    bx2 = 0.0_rp
  end if

  a22  = p27*a2
  a_22 = p28*a2
  a33  = p31*a3
  b33  = p29*b3
  if (b13c > 0.0_rp) then
    b22 = (a13c*a22*a33)/(b13c*b33)
  else
    b22 = 0.0_rp
  end if

  b_33 = p30*b3
  if (ki_on > 0.0_rp) then
    a_33 = (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3)
  else
    a_33 = 0.0_rp
  end if

  if (b13n > 0.0_rp) then
    b_22 = (a_33*a13n*a_22)/(b_33*b13n)
  else
    b_22 = 0.0_rp
  end if

  a44  = p32*a2
  b44  = p33*a3
  a_44 = p34*a2
  b_44 = p35*a2

  ! -------------------------------------------------------------------------
  ! Integracin explcita con substeps (dt_sub)
  ! -------------------------------------------------------------------------
  ! Control: criterio original por Vm.
  ! Vernakalant activo: criterio por rigidez (tasa mxima) para evitar colapso numrico del Markov.
  if (D_verna <= 0.0_rp) then
    if (Vm > -65.0_rp .and. Vm < 40.0_rp) then
      n_sub = 10_ip
    else
      n_sub = 2_ip
    end if
  else
    max_rate = maxval( (/ a11,a12,a13,b11,b12,b13, a2,b2,ax,bx, a3,b3, &
                         a13c,b13c,a13n,b13n, ax1,bx1,ax2,bx2, &
                         a22,a_22,b22,b_22, a33,b33,a_33,b_33, a44,b44,a_44,b_44, &
                         kon,koff,k_on,k_off, ki_on,ki_off,kc_on,kc_off, kcon,kcoff /) )
    if (max_rate < 1.0_rp) max_rate = 1.0_rp
    n_sub = int(ceiling((dt*max_rate)/0.10_rp), kind=ip)
    if (n_sub < 6_ip) n_sub = 6_ip
    if (n_sub > 4000_ip) n_sub = 4000_ip
  end if
  dt_sub = dt/real(n_sub,rp)
  if (dt_sub <= 0.0_rp) dt_sub = dt

  do k_sub = 1, n_sub

    ! Estados previos del substep (esquema semi-implícito diagonal)
    oIC3   = cm_m%xnIC3  ; oIC2   = cm_m%xnIC2  ; oIF    = cm_m%xnIF
    oC3    = cm_m%xnC3   ; oC2    = cm_m%xnC2   ; oC1    = cm_m%xnC1
    oO     = cm_m%xnO    ; oIS    = cm_m%xnIS

    oDpIC3 = cm_m%xnDpIC3; oDpIC2 = cm_m%xnDpIC2; oDpIF  = cm_m%xnDpIF
    oDpC3  = cm_m%xnDpC3 ; oDpC2  = cm_m%xnDpC2 ; oDpC1  = cm_m%xnDpC1
    oDpO   = cm_m%xnDpO  ; oDpIS  = cm_m%xnDpIS ; oDpIT  = cm_m%xnDpIT

    oDIC3  = cm_m%xnDIC3 ; oDIC2  = cm_m%xnDIC2 ; oDIF   = cm_m%xnDIF
    oDC3   = cm_m%xnDC3  ; oDC2   = cm_m%xnDC2  ; oDC1   = cm_m%xnDC1
    oDO    = cm_m%xnDO   ; oDIS   = cm_m%xnDIS  ; oDIT   = cm_m%xnDIT

    ! ======= Bloque nativo xn* =======
    dIC3 = oIC2*b11 + oC3*b3 + ki_off*oDIC3
    lIC3 = a11 + a3 + ki_on

    dIC2 = oIC3*a11 + oIF*b12 + oC2*b3 + ki_off*oDIC2
    lIC2 = b11 + a3 + a12 + ki_on

    dIF  = oIC2*a12 + oC1*b3 + oO*a2 + ki_off*oDIF
    lIF  = b12 + a3 + b2 + ki_on

    dC3  = oIC3*a3 + oC2*b11 + oDpC3*kcoff + oDC3*kc_off
    lC3  = b3 + a11 + kcon + kc_on

    dC2  = oC3*a11 + oIC2*a3 + oC1*b12 + oDpC2*kcoff + oDC2*kc_off
    lC2  = b11 + b3 + a12 + kcon + kc_on

    dC1  = oC2*a12 + oIF*a3 + oO*b13 + oDpC1*kcoff + oDC1*kc_off
    lC1  = b12 + b3 + a13 + kcon + kc_on

    dO_mkv = oC1*a13 + oIF*b2 + oIS*bx + oDpO*koff + oDO*k_off
    lO    = b13 + a2 + ax + kon + k_on

    dIS  = oO*ax + oDIS*ki_off
    lIS  = bx + ki_on

    ! ======= Bloque Dp* (drogado atrapado) =======
    dDpIC3 = oDpIC2*b11 + oDpC3*b33
    lDpIC3 = a33 + a11

    dDpIC2 = oDpIC3*a11 + oDpIF*b12 + oDpC2*b33
    lDpIC2 = b11 + a33 + a12

    dDpIF  = oDpIC2*a12 + oDpC1*b33 + oDpO*a22 + oDpIT*b44
    lDpIF  = b12 + a33 + b22 + a44

    dDpC3  = oDpIC3*a33 + oDpC2*b11 + oC3*kcon
    lDpC3  = b33 + a11 + kcoff

    dDpC2  = oDpC3*a11 + oDpIC2*a33 + oDpC1*b12 + oC2*kcon
    lDpC2  = b11 + b33 + a12 + kcoff

    dDpC1  = oDpC2*a12 + oDpIF*a33 + oDpO*b13c + oC1*kcon
    lDpC1  = b12 + b33 + a13c + kcoff

    dDpO   = oDpC1*a13c + oDpIF*b22 + oDpIS*bx1 + oO*kon
    lDpO   = b13c + a22 + ax1 + koff

    dDpIS  = oDpO*ax1
    lDpIS  = bx1

    dDpIT  = oDpIF*a44
    lDpIT  = b44

    ! ======= Bloque D* (drogado rápido) =======
    dDIC3 = oDIC2*b11 + oDC3*b_33 + ki_on*oIC3
    lDIC3 = a_33 + a11 + ki_off

    dDIC2 = oDIC3*a11 + oDIF*b12 + oDC2*b_33 + ki_on*oIC2
    lDIC2 = b11 + a_33 + a12 + ki_off

    dDIF  = oDIC2*a12 + oDC1*b_33 + oDO*a_22 + oDIT*b_44 + ki_on*oIF
    lDIF  = b12 + a_33 + b_22 + a_44 + ki_off

    dDC3  = oDIC3*a_33 + oDC2*b11 + oC3*kc_on
    lDC3  = b_33 + a11 + kc_off

    dDC2  = oDC3*a11 + oDIC2*a_33 + oDC1*b12 + oC2*kc_on
    lDC2  = b11 + b_33 + a12 + kc_off

    dDC1  = oDC2*a12 + oDIF*a_33 + oDO*b13n + oC1*kc_on
    lDC1  = b12 + b_33 + a13n + kc_off

    dDO   = oDC1*a13n + oDIF*b_22 + oDIS*bx2 + oO*k_on
    lDO   = b13n + a_22 + ax2 + k_off

    dDIS  = oDO*ax2 + oIS*ki_on
    lDIS  = bx2 + ki_off

    dDIT  = oDIF*a_44
    lDIT  = b_44

    ! Actualización semi-implícita: X_{n+1} = (X_n + dt*P_n)/(1 + dt*L_n)
    cm_m%xnIC3   = (oIC3   + dt_sub*dIC3)   /(1.0_rp + dt_sub*lIC3)
    cm_m%xnIC2   = (oIC2   + dt_sub*dIC2)   /(1.0_rp + dt_sub*lIC2)
    cm_m%xnIF    = (oIF    + dt_sub*dIF)    /(1.0_rp + dt_sub*lIF)
    cm_m%xnC3    = (oC3    + dt_sub*dC3)    /(1.0_rp + dt_sub*lC3)
    cm_m%xnC2    = (oC2    + dt_sub*dC2)    /(1.0_rp + dt_sub*lC2)
    cm_m%xnC1    = (oC1    + dt_sub*dC1)    /(1.0_rp + dt_sub*lC1)
    cm_m%xnO     = (oO     + dt_sub*dO_mkv) /(1.0_rp + dt_sub*lO)
    cm_m%xnIS    = (oIS    + dt_sub*dIS)    /(1.0_rp + dt_sub*lIS)

    cm_m%xnDpIC3 = (oDpIC3 + dt_sub*dDpIC3) /(1.0_rp + dt_sub*lDpIC3)
    cm_m%xnDpIC2 = (oDpIC2 + dt_sub*dDpIC2) /(1.0_rp + dt_sub*lDpIC2)
    cm_m%xnDpIF  = (oDpIF  + dt_sub*dDpIF)  /(1.0_rp + dt_sub*lDpIF)
    cm_m%xnDpC3  = (oDpC3  + dt_sub*dDpC3)  /(1.0_rp + dt_sub*lDpC3)
    cm_m%xnDpC2  = (oDpC2  + dt_sub*dDpC2)  /(1.0_rp + dt_sub*lDpC2)
    cm_m%xnDpC1  = (oDpC1  + dt_sub*dDpC1)  /(1.0_rp + dt_sub*lDpC1)
    cm_m%xnDpO   = (oDpO   + dt_sub*dDpO)   /(1.0_rp + dt_sub*lDpO)
    cm_m%xnDpIS  = (oDpIS  + dt_sub*dDpIS)  /(1.0_rp + dt_sub*lDpIS)
    cm_m%xnDpIT  = (oDpIT  + dt_sub*dDpIT)  /(1.0_rp + dt_sub*lDpIT)

    cm_m%xnDIC3  = (oDIC3  + dt_sub*dDIC3)  /(1.0_rp + dt_sub*lDIC3)
    cm_m%xnDIC2  = (oDIC2  + dt_sub*dDIC2)  /(1.0_rp + dt_sub*lDIC2)
    cm_m%xnDIF   = (oDIF   + dt_sub*dDIF)   /(1.0_rp + dt_sub*lDIF)
    cm_m%xnDC3   = (oDC3   + dt_sub*dDC3)   /(1.0_rp + dt_sub*lDC3)
    cm_m%xnDC2   = (oDC2   + dt_sub*dDC2)   /(1.0_rp + dt_sub*lDC2)
    cm_m%xnDC1   = (oDC1   + dt_sub*dDC1)   /(1.0_rp + dt_sub*lDC1)
    cm_m%xnDO    = (oDO    + dt_sub*dDO)    /(1.0_rp + dt_sub*lDO)
    cm_m%xnDIS   = (oDIS   + dt_sub*dDIS)   /(1.0_rp + dt_sub*lDIS)
    cm_m%xnDIT   = (oDIT   + dt_sub*dDIT)   /(1.0_rp + dt_sub*lDIT)

  end do  ! k_sub

  if (.not.((cm_m%xnO == cm_m%xnO) .and. (cm_m%xnC1 == cm_m%xnC1) .and. (cm_m%xnIC3 == cm_m%xnIC3))) then
    cm_m%xnC3    = 1.0_rp
    cm_m%xnIC3   = 0.0_rp; cm_m%xnIC2 = 0.0_rp; cm_m%xnIF  = 0.0_rp
    cm_m%xnC2    = 0.0_rp; cm_m%xnC1  = 0.0_rp
    cm_m%xnO     = 0.0_rp; cm_m%xnIS  = 0.0_rp

    cm_m%xnDpIC3 = 0.0_rp; cm_m%xnDpIC2=0.0_rp; cm_m%xnDpIF =0.0_rp
    cm_m%xnDpC3  = 0.0_rp; cm_m%xnDpC2 =0.0_rp; cm_m%xnDpC1 =0.0_rp
    cm_m%xnDpO   = 0.0_rp; cm_m%xnDpIS =0.0_rp; cm_m%xnDpIT =0.0_rp

    cm_m%xnDIC3  = 0.0_rp; cm_m%xnDIC2=0.0_rp; cm_m%xnDIF  =0.0_rp
    cm_m%xnDC3   = 0.0_rp; cm_m%xnDC2 =0.0_rp; cm_m%xnDC1  =0.0_rp
    cm_m%xnDO    = 0.0_rp; cm_m%xnDIS =0.0_rp; cm_m%xnDIT  =0.0_rp
  end if

  ! -------------------------------------------------------------------------
  ! Clamping y normalizacin a suma = 1
  ! -------------------------------------------------------------------------
  cm_m%xnIC3   = max(0.0_rp, cm_m%xnIC3  )
  cm_m%xnIC2   = max(0.0_rp, cm_m%xnIC2  )
  cm_m%xnIF    = max(0.0_rp, cm_m%xnIF   )
  cm_m%xnC3    = max(0.0_rp, cm_m%xnC3   )
  cm_m%xnC2    = max(0.0_rp, cm_m%xnC2   )
  cm_m%xnC1    = max(0.0_rp, cm_m%xnC1   )
  cm_m%xnO     = max(0.0_rp, cm_m%xnO    )
  cm_m%xnIS    = max(0.0_rp, cm_m%xnIS   )

  cm_m%xnDpIC3 = max(0.0_rp, cm_m%xnDpIC3)
  cm_m%xnDpIC2 = max(0.0_rp, cm_m%xnDpIC2)
  cm_m%xnDpIF  = max(0.0_rp, cm_m%xnDpIF )
  cm_m%xnDpC3  = max(0.0_rp, cm_m%xnDpC3 )
  cm_m%xnDpC2  = max(0.0_rp, cm_m%xnDpC2 )
  cm_m%xnDpC1  = max(0.0_rp, cm_m%xnDpC1 )
  cm_m%xnDpO   = max(0.0_rp, cm_m%xnDpO  )
  cm_m%xnDpIS  = max(0.0_rp, cm_m%xnDpIS )
  cm_m%xnDpIT  = max(0.0_rp, cm_m%xnDpIT )

  cm_m%xnDIC3  = max(0.0_rp, cm_m%xnDIC3 )
  cm_m%xnDIC2  = max(0.0_rp, cm_m%xnDIC2 )
  cm_m%xnDIF   = max(0.0_rp, cm_m%xnDIF  )
  cm_m%xnDC3   = max(0.0_rp, cm_m%xnDC3  )
  cm_m%xnDC2   = max(0.0_rp, cm_m%xnDC2  )
  cm_m%xnDC1   = max(0.0_rp, cm_m%xnDC1  )
  cm_m%xnDO    = max(0.0_rp, cm_m%xnDO   )
  cm_m%xnDIS   = max(0.0_rp, cm_m%xnDIS  )
  cm_m%xnDIT   = max(0.0_rp, cm_m%xnDIT  )

  sum_states = cm_m%xnIC3 + cm_m%xnIC2 + cm_m%xnIF  + cm_m%xnC3  + cm_m%xnC2 + cm_m%xnC1 + &
               cm_m%xnO   + cm_m%xnIS  + &
               cm_m%xnDpIC3 + cm_m%xnDpIC2 + cm_m%xnDpIF + cm_m%xnDpC3 + cm_m%xnDpC2 + &
               cm_m%xnDpC1 + cm_m%xnDpO + cm_m%xnDpIS + cm_m%xnDpIT + &
               cm_m%xnDIC3 + cm_m%xnDIC2 + cm_m%xnDIF + cm_m%xnDC3 + cm_m%xnDC2 + &
               cm_m%xnDC1 + cm_m%xnDO + cm_m%xnDIS + cm_m%xnDIT

  if (sum_states > eps) then
    norm_factor = 1.0_rp/sum_states

    cm_m%xnIC3   = cm_m%xnIC3   * norm_factor
    cm_m%xnIC2   = cm_m%xnIC2   * norm_factor
    cm_m%xnIF    = cm_m%xnIF    * norm_factor
    cm_m%xnC3    = cm_m%xnC3    * norm_factor
    cm_m%xnC2    = cm_m%xnC2    * norm_factor
    cm_m%xnC1    = cm_m%xnC1    * norm_factor
    cm_m%xnO     = cm_m%xnO     * norm_factor
    cm_m%xnIS    = cm_m%xnIS    * norm_factor

    cm_m%xnDpIC3 = cm_m%xnDpIC3 * norm_factor
    cm_m%xnDpIC2 = cm_m%xnDpIC2 * norm_factor
    cm_m%xnDpIF  = cm_m%xnDpIF  * norm_factor
    cm_m%xnDpC3  = cm_m%xnDpC3  * norm_factor
    cm_m%xnDpC2  = cm_m%xnDpC2  * norm_factor
    cm_m%xnDpC1  = cm_m%xnDpC1  * norm_factor
    cm_m%xnDpO   = cm_m%xnDpO   * norm_factor
    cm_m%xnDpIS  = cm_m%xnDpIS  * norm_factor
    cm_m%xnDpIT  = cm_m%xnDpIT  * norm_factor

    cm_m%xnDIC3  = cm_m%xnDIC3  * norm_factor
    cm_m%xnDIC2  = cm_m%xnDIC2  * norm_factor
    cm_m%xnDIF   = cm_m%xnDIF   * norm_factor
    cm_m%xnDC3   = cm_m%xnDC3   * norm_factor
    cm_m%xnDC2   = cm_m%xnDC2   * norm_factor
    cm_m%xnDC1   = cm_m%xnDC1   * norm_factor
    cm_m%xnDO    = cm_m%xnDO    * norm_factor
    cm_m%xnDIS   = cm_m%xnDIS   * norm_factor
    cm_m%xnDIT   = cm_m%xnDIT   * norm_factor
  else
    ! Reset de seguridad
    cm_m%xnC3    = 1.0_rp
    cm_m%xnIC3   = 0.0_rp; cm_m%xnIC2 = 0.0_rp; cm_m%xnIF  = 0.0_rp
    cm_m%xnC2    = 0.0_rp; cm_m%xnC1  = 0.0_rp
    cm_m%xnO     = 0.0_rp; cm_m%xnIS  = 0.0_rp

    cm_m%xnDpIC3 = 0.0_rp; cm_m%xnDpIC2=0.0_rp; cm_m%xnDpIF =0.0_rp
    cm_m%xnDpC3  = 0.0_rp; cm_m%xnDpC2 =0.0_rp; cm_m%xnDpC1 =0.0_rp
    cm_m%xnDpO   = 0.0_rp; cm_m%xnDpIS =0.0_rp; cm_m%xnDpIT =0.0_rp

    cm_m%xnDIC3  = 0.0_rp; cm_m%xnDIC2=0.0_rp; cm_m%xnDIF  =0.0_rp
    cm_m%xnDC3   = 0.0_rp; cm_m%xnDC2 =0.0_rp; cm_m%xnDC1  =0.0_rp
    cm_m%xnDO    = 0.0_rp; cm_m%xnDIS =0.0_rp; cm_m%xnDIT  =0.0_rp
  end if

  ! =====================================================================================
  ! 3. Resto de compuertas K (Courtemanche original)
  ! =====================================================================================

  ! xr (IKr)
  if (abs(Vm+14.1_rp) < 1.0e-10_rp) then
     alfa_xr = 0.0015_rp
  else
     alfa_xr = 0.0003_rp*(Vm+14.1_rp)/(1.0_rp-exp((Vm+14.1_rp)/(-5.0_rp)))
  end if
  if (abs(Vm-3.3328_rp) < 1.0e-10_rp) then
     beta_xr = 3.7836118e-4_rp
  else
     beta_xr = 0.000073898_rp*(Vm-3.3328_rp)/(exp((Vm-3.3328_rp)/5.1237_rp)-1.0_rp)
  end if
  tau_xr = 1.0_rp/(alfa_xr+beta_xr)
  xr_inf = 1.0_rp/(1.0_rp+exp((Vm+14.1_rp)/(-6.5_rp)))
  cm_m%xr = xr_inf - (xr_inf-cm_m%xr)*exp(-dt/tau_xr)

  ! xs (IKs)
  if (abs(Vm-19.9_rp) < 1.0e-10_rp) then
     alfa_xs = 0.00068_rp
     beta_xs = 0.000315_rp
  else
     alfa_xs = 0.00004_rp*(Vm-19.9_rp)/(1.0_rp-exp((Vm-19.9_rp)/(-17.0_rp)))
     beta_xs = 0.000035_rp*(Vm-19.9_rp)/(exp((Vm-19.9_rp)/9.0_rp)-1.0_rp)
  end if
  tau_xs = 0.5_rp/(alfa_xs+beta_xs)
  xs_inf = 1.0_rp/sqrt(1.0_rp+exp((Vm-19.9_rp)/(-12.7_rp)))
  cm_m%xs = xs_inf - (xs_inf-cm_m%xs)*exp(-dt/tau_xs)

  ! Ito: oa / oi
  alfa_oa = 0.65_rp/(exp((Vm+10.0_rp)/(-8.5_rp))+exp((Vm-30.0_rp)/(-59.0_rp)))
  beta_oa = 0.65_rp/(2.5_rp+exp((Vm+82.0_rp)/17.0_rp))
  tau_oa  = 1.0_rp/((alfa_oa+beta_oa)*p_KQ10)
  oa_inf  = 1.0_rp/(1.0_rp+exp((Vm+20.47_rp)/(-17.54_rp)))
  cm_m%oa = oa_inf-(oa_inf-cm_m%oa)*exp(-dt/tau_oa)

  alfa_oi = 1.0_rp/(18.53_rp+exp((Vm+113.7_rp)/10.95_rp))
  beta_oi = 1.0_rp/(35.56_rp+exp((Vm+1.26_rp)/(-7.44_rp)))
  tau_oi  = 1.0_rp/((alfa_oi+beta_oi)*p_KQ10)
  oi_inf  = 1.0_rp/(1.0_rp+exp((Vm+43.1_rp)/5.3_rp))
  cm_m%oi = oi_inf-(oi_inf-cm_m%oi)*exp(-dt/tau_oi)

  ! IKur: ua / ui
  alfa_ua = 0.65_rp/(exp((Vm+10.0_rp)/(-8.5_rp))+exp((Vm-30.0_rp)/(-59.0_rp)))
  beta_ua = 0.65_rp/(2.5_rp+exp((Vm+82.0_rp)/17.0_rp))
  tau_ua  = 1.0_rp/((alfa_ua+beta_ua)*p_KQ10)
  ua_inf  = 1.0_rp/(1.0_rp+exp((Vm+30.3_rp)/(-9.6_rp)))
  cm_m%ua = ua_inf-(ua_inf-cm_m%ua)*exp(-dt/tau_ua)

  alfa_ui = 1.0_rp/(21.0_rp+exp((Vm-185.0_rp)/(-28.0_rp)))
  beta_ui = 1.0_rp/(exp((Vm-158.0_rp)/(-16.0_rp)))
  tau_ui  = 1.0_rp/((alfa_ui+beta_ui)*p_KQ10)
  ui_inf  = 1.0_rp/(1.0_rp+exp((Vm-99.45_rp)/27.48_rp))
  cm_m%ui = ui_inf-(ui_inf-cm_m%ui)*exp(-dt/tau_ui)

  return
end subroutine f_gates



!-------------------------------------------------------------------------------
subroutine Courtemanche_A_P01(ict, dt, Vm, Istm, Iion, v_prm, v_courte, v_cr)
  implicit none
  
  integer(ip), intent(in)    :: ict
  real(rp),    intent(in)    :: dt, Istm, Vm
  real(rp),    intent(in)    :: v_prm(:)
  real(rp),    intent(out)   :: Iion
  real(rp),    intent(inout) :: v_courte(nvar_courte)
  real(rp),    intent(out)   :: v_cr(ncur_courte)

  type(t_courte) :: courte
  type(t_prm)    :: param
  type(t_cur)    :: cur

  ! --------------------------------------------------------------
  ! Cargar parametros y estado del modelo
  ! --------------------------------------------------------------
  call put_me_struct(v_courte, courte)
  call put_param(v_prm, param)

  ! --------------------------------------------------------------
  ! ORDEN CORRECTO
  ! --------------------------------------------------------------

  ! 1) Actualizar compuertas (incluye Markov de INa)
  call f_gates(dt, Vm, cur, courte)

  ! 2) Calcular corrientes Iion
  call f_currents(Vm, param, cur, courte, Iion)

  ! 3) Actualizar concentraciones (usa cur del paso actual)
  call f_concentrations(dt, cur, courte, Istm)

  ! --------------------------------------------------------------
  ! Devolver resultados
  ! --------------------------------------------------------------
  call get_me_struct(courte, v_courte)

  v_cr(1:ncur_courte) = (/ cur%INa, cur%IK1, cur%Ito, cur%IKur, cur%IKr, cur%IKs, &
                           cur%ICaL, cur%INaK, cur%INaCa, cur%IbCa, cur%IbNa,     &
                           cur%IpCa, cur%Irel, cur%Itr, cur%Iup, cur%Ileak,       &
                           cur%IKACh, Iion /)

  return
end subroutine Courtemanche_A_P01
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 !-------------------------------------------------------------------------------
end module mod_courtemanche
!-------------------------------------------------------------------------------
