!-----------------------------------------------------------------------------------------------------------------------------
! This are the parameters for the Courtemanche Model (M Courtemanche, RJ Ramirez, S Nattel; Am J Physiol 275:H301-H321, 1998)
! Added IKACh formulation of Grandi (E.Grandi,S.V.Pandit,N.Voigt,A.J.Workman,D.Dobrev,J.Jalife,D.M.Bers;Circulation Research,vol.109,pp.1055-1066,2011)
! Implemented by Laura Martínez Mateu (laumarma@gbio.i3bh.es)  26/05/2014
! Checked  by JFR (jfrodrig@unizar.es)  31/03/2014
! Checked by Violeta Puche García (viopucga@ci2b.upv.es) 08/02/2023 
!
! p_INakmax must be increased by a 40% in order to stabilize Nai and Ki homeostasis and do not loose
! the spike and dome morphology of the AP
!
!=============
! pAF
!=============
!
! Constants
!-----------------------------------------------------------------------------------------------------------------------------
! Model Constants
! (Parameter hace que sean constantes simbólicas y que su valor no pueda cambiar a lo largo de la ejecución del programa)

! Cm=100pF

 real(rp), parameter, private :: p_R        = 8.3143;       ! J/mol·K    
 real(rp), parameter, private :: p_F        = 96.4867;      ! C/mmol  
 real(rp), parameter, private :: p_T        = 310.0;	    ! K 
 real(rp), parameter, private :: p_Cm       = 100.0;        ! pF
 
 real(rp), parameter, private :: p_Vcell    = 20100.0;      ! um^3
 real(rp), parameter, private :: p_Vi       = 13668.0;      ! um^3
 real(rp), parameter, private :: p_Vup      = 1109.52;      ! um^3
 real(rp), parameter, private :: p_Vrel     = 96.48;        ! um^3
 
 real(rp), parameter, private :: p_Ko       = 5.4;          ! mM
 real(rp), parameter, private :: p_Nao      = 140.0;        ! mM
 real(rp), parameter, private :: p_Cao      = 1.8;          ! mM
 real(rp), parameter, private :: p_verna_conc_nM = 0.0;  ! CONTROL

 !-----------------------------------------------------------------------------------------------------------
 real(rp), parameter, private :: p_gNa      = 7.8*0.58;	   		! nS/pF 
 real(rp), parameter, private :: p_gK1_RA   = 0.09;	    	! nS/pF
 real(rp), parameter, private :: p_gK1_LA   = 1.2*0.09;   	! nS/pF  - Irene 29/07/2024  remodelado pAF en LA (1.2·IK1 en todas las regiones)
 real(rp), parameter, private :: p_gto      = 0.1652;	    ! nS/pF 
 real(rp), parameter, private :: p_gKr      = 0.0294;  	    ! nS/pF 
 real(rp), parameter, private :: p_gKs      = 0.12941176;  		! nS/pF 
 real(rp), parameter, private :: p_gCaL     = 0.1238;	    ! nS/pF
 real(rp), parameter, private :: p_gbCa     = 0.00113;	    ! nS/pF
 real(rp), parameter, private :: p_gbNa     = 0.000674; 	! nS/pF
 real(rp), parameter, private :: p_gKur     = 0.005;        ! nS/pF  (valor para las CI del paper)


 ! Región CT-BBra
 real(rp), parameter, private :: p_gCaL_CT_BBra  = 1.67*0.1238;	! nS/pF  (g_CaL_CT_BBra/g_CaL=1.67)
 
 ! Región TVR
 real(rp), parameter, private :: p_gCaL_TVR = 0.8*0.1238;	    ! nS/pF  (g_CaL_TVR/g_CaL=0.8) 
 real(rp), parameter, private :: p_gKr_TVR  = 1.3*0.0294;	    ! nS/pF  (g_Kr_TVR/g_Kr=1.3) 

 ! Región RAA
 real(rp), parameter, private :: p_gto_RAA  = 0.68*0.1652;	    ! nS/pF  (gto_RAA/gto=0.68)
 real(rp), parameter, private :: p_gCaL_RAA = 0.9*0.1238;	    ! nS/pF  (g_CaL_RAA/g_CaL=0.9) 
 real(rp), parameter, private :: p_gKr_RAA  = 1.1*0.0294;	    ! nS/pF  (gKr_RAA/gKr=1.1) 
 
 ! Región LA
 real(rp), parameter, private :: p_gCaL_LA = 0.9*0.1238;	    ! nS/pF  (g_CaL_LA/g_CaL=0.9) 
 real(rp), parameter, private :: p_gKr_LA  = 2.0*0.0294;	    ! nS/pF  (gKr_LA/gKr=2.0) 
 
 ! Región BBla
 real(rp), parameter, private :: p_gKr_BBla  = 1.6*0.0294;	    ! nS/pF  (gKr_BBla/gKr=1.6) 
 real(rp), parameter, private :: p_gCaL_BBla = 1.67*0.1238;	    ! nS/pF  (g_CaL_BBla/g_CaL=1.67) 
 
 ! Región PV
 real(rp), parameter, private :: p_gto_PV  = 0.9*0.1652;	    ! nS/pF  (gto_PV/gto=0.9) 
 real(rp), parameter, private :: p_gCaL_PV  = 0.8*0.1238;	    ! nS/pF  (gCaL_PV/gCaL=0.8) 
 real(rp), parameter, private :: p_gKr_PV  = 2.5*0.0294;	    ! nS/pF  (gKr_PV/gKr=2.5) 
 real(rp), parameter, private :: p_gKs_PV  = 1.9*0.129;	    	! nS/pF  (gKs_PV/gKs=1.9) 
 real(rp), parameter, private :: p_gK1_PV  = 1.2*0.9*0.0900;	! nS/pF  (gK1_PV/gK1=0.9) ! nS/pF  - Violeta 21/06/2022  remodeling for paroxysmal AF in LA (2·IK1 in all regions of LA)

 ! Región LAA
 real(rp), parameter, private :: p_gto_LAA = 1.0*0.1652;	    ! nS/pF  (gto_LAA/gto=0.68)
 real(rp), parameter, private :: p_gCaL_LAA = 0.9*0.1238;	    ! nS/pF  (gCaL_LAA/gCaL=0.75) 
 real(rp), parameter, private :: p_gKr_LAA = 2.2*0.0294;	    ! nS/pF  (gKr_LAA/gKr=2.2) 

 ! Región MVR
 real(rp), parameter, private :: p_gCaL_MVR = 0.8*0.1238;	    ! nS/pF  (gCaL_MVR/gCaL=0.8) 
 real(rp), parameter, private :: p_gKr_MVR  = 2.6*0.0294;	    ! nS/pF  (gKr_MVR/gKr=2.6)

 !-------------------------------------------------------------------------------------------------------
 
 real(rp), parameter, private :: p_INaKmax  = 1.0*0.6;	    ! pA/pF: para hacer pruebas de mejor estabilizacion del modelo subiendo el valor
 real(rp), parameter, private :: p_INaCamax = 1.0*1600.0;   ! pA/pF
 real(rp), parameter, private :: p_IpCamax  = 0.275;	    ! pA/pF 

 real(rp), parameter, private :: p_Iupmax   = 0.005;	    ! mM/ms

 real(rp), parameter, private :: p_KQ10     = 3.0;	        ! 
 real(rp), parameter, private :: p_gamma    = 0.35;	        ! 
 
 real(rp), parameter, private :: p_KmNai    = 10.0;          ! mM
 real(rp), parameter, private :: p_KmKo     = 1.5;           ! mM
 real(rp), parameter, private :: p_KmNa     = 87.5;          ! mM
 
 real(rp), parameter, private :: p_KmCa     = 1.38;          ! 
 real(rp), parameter, private :: p_ksat     = 0.1;           !
 
 real(rp), parameter, private :: p_krel     = 30.0;          ! 1/ms 
 
 real(rp), parameter, private :: p_Kup      = 0.00092;       ! mM
 real(rp), parameter, private :: p_Caupmax  = 15.0;          ! mM
 real(rp), parameter, private :: p_CMDNmax  = 0.05;          ! mM
 real(rp), parameter, private :: p_TRPNmax  = 0.07;          ! mM
 real(rp), parameter, private :: p_CSQNmax  = 10.0;          ! mM
 real(rp), parameter, private :: p_KmCmdn   = 0.00238;       ! mM
 real(rp), parameter, private :: p_KmTrpn   = 0.0005;        ! mM 
 real(rp), parameter, private :: p_KmCsqn   = 0.8;           ! mM
  
 real(rp), parameter, private :: p_ACh      = 0.005;  	     ! uM	
 real(rp), parameter, private :: p_IKACh_RA = 1.0;           ! Rebecca 30/03/2020  remodeling for paroxysmal AF in LA (2·IKACh en LA)
 real(rp), parameter, private :: p_IKACh_LA = 2.0;	         ! Rebecca 30/03/2020  remodeling for paroxysmal AF in LA (2·IKACh en LA)
 
!------------------------------------------------------------------------------
!-- INITIAL CONDITIONS --------------------------------------------------------
!------------------------------------------------------------------------------
!
! Partiendo de las Cond.Ini. del paper
!
  real(rp), parameter :: Vi_courte    = -81.2;
  real(rp), parameter :: Naiic        = 1.12e+01;
  real(rp), parameter :: xnISic       = 0.00e+00;
  real(rp), parameter :: xnOic        = 0.00e+00;
  real(rp), parameter :: xnDpIC3ic    = 0.00e+00;
  real(rp), parameter :: Kiic         = 1.39e+02;
  real(rp), parameter :: oaic         = 3.04e-02;
  real(rp), parameter :: oiic         = 9.99e-01;
  real(rp), parameter :: uaic         = 4.96e-03;
  real(rp), parameter :: uiic         = 9.99e-01;
  real(rp), parameter :: xric         = 3.29e-05;
  real(rp), parameter :: xsic         = 1.87e-02;
  real(rp), parameter :: Caiic        = 1.02e-04;
  real(rp), parameter :: dic          = 1.37e-04;
  real(rp), parameter :: fic          = 9.99e-01;
  real(rp), parameter :: fCaic        = 7.75e-01;
  real(rp), parameter :: Carelic      = 1.49e+00;
  real(rp), parameter :: uic          = 0.00e+00;
  real(rp), parameter :: vic          = 1.00e+00;
  real(rp), parameter :: wic          = 9.99e-01;
  real(rp), parameter :: Caupic       = 1.49e+00;
  real(rp), parameter :: Fnic         = 0.0e+00; 
  real(rp), parameter :: tau_uic      = 0.0e+00;
  real(rp), parameter :: u_infic      = 0.0e+00;
  real(rp), parameter :: xnIC3ic      = 0.0e+00;
  real(rp), parameter :: xnIC2ic      = 0.0e+00;
  real(rp), parameter :: xnIFic       = 0.0e+00;
  real(rp), parameter :: xnC3ic       = 1.0e+00;
  real(rp), parameter :: xnC2ic       = 0.0e+00;
  real(rp), parameter :: xnC1ic       = 0.0e+00;
  real(rp), parameter :: xnDpIC2ic    = 0.00e+00;
  real(rp), parameter :: xnDpIFic     = 0.00e+00;
  real(rp), parameter :: xnDpC3ic     = 0.00e+00;
  real(rp), parameter :: xnDpC2ic     = 0.00e+00;
  real(rp), parameter :: xnDpC1ic     = 0.00e+00;
  real(rp), parameter :: xnDpOic      = 0.00e+00;
  real(rp), parameter :: xnDpISic     = 0.00e+00;
  real(rp), parameter :: xnDpITic     = 0.00e+00;
  real(rp), parameter :: xnDIC3ic     = 0.00e+00;
  real(rp), parameter :: xnDIC2ic     = 0.00e+00;
  real(rp), parameter :: xnDIFic      = 0.00e+00;
  real(rp), parameter :: xnDC3ic      = 0.00e+00;
  real(rp), parameter :: xnDC2ic      = 0.00e+00;
  real(rp), parameter :: xnDC1ic      = 0.00e+00;
  real(rp), parameter :: xnDOic       = 0.00e+00;
  real(rp), parameter :: xnDISic      = 0.00e+00;
  real(rp), parameter :: xnDITic      = 0.00e+00;
 ! real(rp), parameter :: alfa_hic     = 0.0e+00;
 ! real(rp), parameter :: beta_hic     = 0.0e+00;

! Condiciones iniciales para la heterogeneidad en la auricula
!

 ! RA_PM
  real(rp), parameter :: Vi_courte_RA_PM 		= -0.8017894179E+02;
  real(rp), parameter :: uic_RA_PM       		= 0.4643259802E-25;
  real(rp), parameter :: vic_RA_PM       		= 0.1000000000E+01;
  real(rp), parameter :: wic_RA_PM       		= 0.9991498677E+00;
  real(rp), parameter :: dic_RA_PM       		= 0.1549804000E-03;
  real(rp), parameter :: xnOic_RA_PM       		= 0.0000000000E+00;
  real(rp), parameter :: xric_RA_PM      		= 0.1687206099E-01;
  real(rp), parameter :: Naiic_RA_PM     		= 0.1162143384E+02;
  real(rp), parameter :: Kiic_RA_PM      		= 0.1385878017E+03;
  real(rp), parameter :: Carelic_RA_PM   		= 0.5928295288E+00;
  real(rp), parameter :: oiic_RA_PM      		= 0.9990699746E+00;
  real(rp), parameter :: uiic_RA_PM      		= 0.9886362973E+00;
  real(rp), parameter :: xnISic_RA_PM       		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_RA_PM       		= 0.0000000000E+00;
  real(rp), parameter :: fic_RA_PM       		= 0.8207060627E+00;
  real(rp), parameter :: xsic_RA_PM      		= 0.2911562135E-01;
  real(rp), parameter :: oaic_RA_PM      		= 0.3218598428E-01;
  real(rp), parameter :: uaic_RA_PM      		= 0.5516193670E-02;
  real(rp), parameter :: fCaic_RA_PM     		= 0.6897953387E+00;
  real(rp), parameter :: Caiic_RA_PM     		= 0.1569521360E-03;
  real(rp), parameter :: Caupic_RA_PM    		= 0.1622678749E+01;
  real(rp), parameter :: xnIC3ic_RA_PM                  = 0.0e+00;
  real(rp), parameter :: xnIC2ic_RA_PM                  = 0.0e+00;
  real(rp), parameter :: xnIFic_RA_PM                   = 0.0e+00;
  real(rp), parameter :: xnC3ic_RA_PM                   = 1.0e+00;
  real(rp), parameter :: xnC2ic_RA_PM                   = 0.0e+00;
  real(rp), parameter :: xnC1ic_RA_PM                   = 0.0e+00;
  real(rp), parameter :: xnDpIC2ic_RA_PM                = 0.00e+00;
  real(rp), parameter :: xnDpIFic_RA_PM                 = 0.00e+00;
  real(rp), parameter :: xnDpC3ic_RA_PM     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_RA_PM   		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_RA_PM    		= 0.00e+00;
  real(rp), parameter :: xnDpOic_RA_PM   		= 0.00e+00;
  real(rp), parameter :: xnDpISic_RA_PM                 = 0.00e+00;
  real(rp), parameter :: xnDpITic_RA_PM                 = 0.00e+00;
  real(rp), parameter :: xnDIC3ic_RA_PM                 = 0.00e+00;
  real(rp), parameter :: xnDIC2ic_RA_PM                 = 0.00e+00;  
  real(rp), parameter :: xnDIFic_RA_PM                  = 0.00E+00;
  real(rp), parameter :: xnDC3ic_RA_PM                  = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_RA_PM                  = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_RA_PM                  = 0.00e+00; 
  real(rp), parameter :: xnDOic_RA_PM                   = 0.00e+00; 
  real(rp), parameter :: xnDISic_RA_PM                  = 0.00e+00; 
  real(rp), parameter :: xnDITic_RA_PM                  = 0.00e+00; 

 
 ! CT_BBra 
  real(rp), parameter :: Vi_courte_CT_BBra    		= -0.7858070269E+02;
  real(rp), parameter :: uic_CT_BBra          		= 0.7394082440E-25;
  real(rp), parameter :: vic_CT_BBra         	 	= 0.1000000000E+01;
  real(rp), parameter :: wic_CT_BBra          		= 0.9990660868E+00;
  real(rp), parameter :: dic_CT_BBra          		= 0.1892900518E-03;
  real(rp), parameter :: xnOic_CT_BBra         		= 0.0000000000E+00;
  real(rp), parameter :: xric_CT_BBra         		= 0.4416142316E-01;
  real(rp), parameter :: Naiic_CT_BBra     	        = 0.1188539577E+02;
  real(rp), parameter :: Kiic_CT_BBra         		= 0.1381901058E+03;
  real(rp), parameter :: Carelic_CT_BBra     	 	= 0.9267038574E+00;
  real(rp), parameter :: oiic_CT_BBra         		= 0.9987231566E+00;
  real(rp), parameter :: uiic_CT_BBra         		= 0.9862434466E+00;
  real(rp), parameter :: xnISic_CT_BBra       		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_CT_BBra    		= 0.0000000000E+00;
  real(rp), parameter :: fic_CT_BBra          		= 0.7234092328E+00;
  real(rp), parameter :: xsic_CT_BBra         		= 0.3843517678E-01;
  real(rp), parameter :: oaic_CT_BBra         		= 0.3516567131E-01;
  real(rp), parameter :: uaic_CT_BBra         		= 0.6514865666E-02;
  real(rp), parameter :: fCaic_CT_BBra        		= 0.6145508829E+00;
  real(rp), parameter :: Caiic_CT_BBra        		= 0.2187829270E-03;
  real(rp), parameter :: Caupic_CT_BBra       		= 0.2241535125E+01;
  real(rp), parameter :: xnIC3ic_CT_BBra                = 0.00e+00;
  real(rp), parameter :: xnIC2ic_CT_BBra                = 0.00e+00;
  real(rp), parameter :: xnIFic_CT_BBra                 = 0.00e+00;
  real(rp), parameter :: xnC3ic_CT_BBra                 = 1.00e+00;
  real(rp), parameter :: xnC2ic_CT_BBra                 = 0.00e+00;
  real(rp), parameter :: xnC1ic_CT_BBra                 = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_CT_BBra              = 0.00e+00;
  real(rp), parameter :: xnDpIFic_CT_BBra               = 0.00e+00;
  real(rp), parameter :: xnDpC3ic_CT_BBra    		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_CT_BBra  		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_CT_BBra    		= 0.00e+00;
  real(rp), parameter :: xnDpOic_CT_BBra  		= 0.00e+00;
  real(rp), parameter :: xnDpISic_CT_BBra               = 0.00e+00;
  real(rp), parameter :: xnDpITic_CT_BBra               = 0.00e+00;
  real(rp), parameter :: xnDIC3ic_CT_BBra               = 0.00e+00;
  real(rp), parameter :: xnDIC2ic_CT_BBra               = 0.00e+00;
  real(rp), parameter :: xnDIFic_CT_BBra                = 0.00E+00;
  real(rp), parameter :: xnDC3ic_CT_BBra                = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_CT_BBra                = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_CT_BBra                = 0.00e+00; 
  real(rp), parameter :: xnDOic_CT_BBra                 = 0.00e+00; 
  real(rp), parameter :: xnDISic_CT_BBra                = 0.00e+00; 
  real(rp), parameter :: xnDITic_CT_BBra                = 0.00e+00; 



 ! TVR  
  real(rp), parameter :: Vi_courte_TVR    		= -0.8084771577E+02;
  real(rp), parameter :: uic_TVR          		= 0.3893329133E-25;
  real(rp), parameter :: vic_TVR          		= 0.1000000000E+01;
  real(rp), parameter :: wic_TVR          		= 0.9991826549E+00;
  real(rp), parameter :: dic_TVR          		= 0.1425400962E-03;
  real(rp), parameter :: xnOic_TVR          		= 0.0000000000E+00;
  real(rp), parameter :: xric_TVR         		= 0.1009298127E-01;
  real(rp), parameter :: Naiic_TVR        		= 0.1159964716E+02;
  real(rp), parameter :: Kiic_TVR         		= 0.1386668183E+03;
  real(rp), parameter :: Carelic_TVR      		= 0.4740667700E+00;
  real(rp), parameter :: oiic_TVR         		= 0.9991844535E+00;
  real(rp), parameter :: uiic_TVR         		= 0.9892974438E+00;
  real(rp), parameter :: xnISic_TVR          		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_TVR         		= 0.0000000000E+00;
  real(rp), parameter :: fic_TVR          		= 0.8723014144E+00;
  real(rp), parameter :: xsic_TVR         		= 0.2503899935E-01;
  real(rp), parameter :: oaic_TVR         		= 0.3101347320E-01;
  real(rp), parameter :: uaic_TVR         		= 0.5145110775E-02;
  real(rp), parameter :: fCaic_TVR        		= 0.7253050452E+00;
  real(rp), parameter :: Caiic_TVR        		= 0.1322202957E-03;
  real(rp), parameter :: Caupic_TVR       		= 0.1361221100E+01;
  real(rp), parameter :: xnIC3ic_TVR                    = 0.00e+00;
  real(rp), parameter :: xnIC2ic_TVR                    = 0.00e+00;
  real(rp), parameter :: xnIFic_TVR                     = 0.00e+00;
  real(rp), parameter :: xnC3ic_TVR                     = 1.00e+00;
  real(rp), parameter :: xnC2ic_TVR                     = 0.00e+00;
  real(rp), parameter :: xnC1ic_TVR                     = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_TVR                  = 0.00e+00;
  real(rp), parameter :: xnDpIFic_TVR      		= 0.00e+00;
  real(rp), parameter :: xnDpC3ic_TVR     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_TVR      		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_TVR      		= 0.00e+00;
  real(rp), parameter :: xnDpOic_TVR       	        = 0.00e+00;
  real(rp), parameter :: xnDpISic_TVR     		= 0.00e+00;
  real(rp), parameter :: xnDpITic_TVR      		= 0.00e+00;
  real(rp), parameter :: xnDIC3ic_TVR     		= 0.00e+00;
  real(rp), parameter :: xnDIC2ic_TVR                   = 0.00e+00;
  real(rp), parameter :: xnDIFic_TVR                    = 0.00E+00;
  real(rp), parameter :: xnDC3ic_TVR                    = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_TVR                    = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_TVR                    = 0.00e+00; 
  real(rp), parameter :: xnDOic_TVR                     = 0.00e+00; 
  real(rp), parameter :: xnDISic_TVR                    = 0.00e+00; 
  real(rp), parameter :: xnDITic_TVR                    = 0.00e+00; 



 ! RAA   
  real(rp), parameter :: Vi_courte_RAA   		= -0.8032425505E+02;
  real(rp), parameter :: uic_RAA         		= 0.6332796299E-25;
  real(rp), parameter :: vic_RAA         		= 0.1000000000E+01;
  real(rp), parameter :: wic_RAA         		= 0.9991571033E+00;
  real(rp), parameter :: dic_RAA         		= 0.1521870255E-03;
  real(rp), parameter :: xnOic_RAA         		= 0.0000000000E+00;
  real(rp), parameter :: xric_RAA        		= 0.1902469380E-01;
  real(rp), parameter :: Naiic_RAA       		= 0.1164611638E+02;
  real(rp), parameter :: Kiic_RAA        		= 0.1385699248E+03;
  real(rp), parameter :: Carelic_RAA     		= 0.5857070128E+00;
  real(rp), parameter :: oiic_RAA        		= 0.9990965497E+00;
  real(rp), parameter :: uiic_RAA        		= 0.9873262916E+00;
  real(rp), parameter :: xnISic_RAA         		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_RAA         		= 0.0000000000E+00;
  real(rp), parameter :: fic_RAA         		= 0.8250495611E+00;
  real(rp), parameter :: xsic_RAA        		= 0.2787387330E-01;
  real(rp), parameter :: oaic_RAA        		= 0.3192714839E-01;
  real(rp), parameter :: uaic_RAA        		= 0.5433210907E-02;
  real(rp), parameter :: fCaic_RAA       		= 0.6952162191E+00;
  real(rp), parameter :: Caiic_RAA       		= 0.1530183364E-03;
  real(rp), parameter :: Caupic_RAA      		= 0.1589463427E+01;
  real(rp), parameter :: xnIC3ic_RAA                    = 0.00e+00;
  real(rp), parameter :: xnIC2ic_RAA                    = 0.00e+00;
  real(rp), parameter :: xnIFic_RAA                     = 0.00e+00;
  real(rp), parameter :: xnC3ic_RAA                     = 1.00e+00;
  real(rp), parameter :: xnC2ic_RAA                     = 0.00e+00;
  real(rp), parameter :: xnC1ic_RAA                     = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_RAA                  = 0.00e+00;
  real(rp), parameter :: xnDpIFic_RAA     		= 0.00e+00;
  real(rp), parameter :: xnDpC3ic_RAA     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_RAA      		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_RAA      		= 0.00e+00;
  real(rp), parameter :: xnDpOic_RAA      	        = 0.00e+00;
  real(rp), parameter :: xnDpISic_RAA    		= 0.00e+00;
  real(rp), parameter :: xnDpITic_RAA      		= 0.00e+00;
  real(rp), parameter :: xnDIC3ic_RAA     		= 0.00e+00;
  real(rp), parameter :: xnDIC2ic_RAA                   = 0.00e+00;
  real(rp), parameter :: xnDIFic_RAA                    = 0.00E+00;
  real(rp), parameter :: xnDC3ic_RAA                    = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_RAA                    = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_RAA                    = 0.00e+00; 
  real(rp), parameter :: xnDOic_RAA                     = 0.00e+00; 
  real(rp), parameter :: xnDISic_RAA                    = 0.00e+00; 
  real(rp), parameter :: xnDITic_RAA                    = 0.00e+00; 



 ! LA
  real(rp), parameter :: Vi_courte_LA    		= -0.8400572389E+02;
  real(rp), parameter :: uic_LA          		= 0.5583532915E-25;
  real(rp), parameter :: vic_LA          		= 0.1000000000E+01;
  real(rp), parameter :: wic_LA          		= 0.9993211516E+00;
  real(rp), parameter :: dic_LA          		= 0.9604120565E-04;
  real(rp), parameter :: xnOic_LA          		= 0.0000000000E+00;
  real(rp), parameter :: xric_LA         		= 0.6913416041E-02;
  real(rp), parameter :: Naiic_LA        		= 0.1184482571E+02;
  real(rp), parameter :: Kiic_LA         		= 0.1384557260E+03;
  real(rp), parameter :: Carelic_LA      		= 0.4017111780E+00;
  real(rp), parameter :: oiic_LA         		= 0.9995538651E+00;
  real(rp), parameter :: uiic_LA         		= 0.9882096268E+00;
  real(rp), parameter :: xnISic_LA          		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_LA          		= 0.0000000000E+00;
  real(rp), parameter :: fic_LA          		= 0.8958224757E+00;
  real(rp), parameter :: xsic_LA         		= 0.1961180518E-01;
  real(rp), parameter :: oaic_LA         		= 0.2602897181E-01;
  real(rp), parameter :: uaic_LA         		= 0.3706142541E-02;
  real(rp), parameter :: fCaic_LA        		= 0.7487122692E+00;
  real(rp), parameter :: Caiic_LA        		= 0.1171948061E-03;
  real(rp), parameter :: Caupic_LA       		= 0.1210484599E+01;
  real(rp), parameter :: xnIC3ic_LA                     = 0.00e+00;
  real(rp), parameter :: xnIC2ic_LA                     = 0.00e+00;
  real(rp), parameter :: xnIFic_LA                      = 0.00e+00;
  real(rp), parameter :: xnC3ic_LA                      = 1.00e+00;
  real(rp), parameter :: xnC2ic_LA                      = 0.00e+00;
  real(rp), parameter :: xnC1ic_LA                      = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_LA                   = 0.00e+00;
  real(rp), parameter :: xnDpIFic_LA     		= 0.00e+00;
  real(rp), parameter :: xnDpC3ic_LA     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_LA      		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_LA      		= 0.00e+00;
  real(rp), parameter :: xnDpOic_LA      	        = 0.00e+00;
  real(rp), parameter :: xnDpISic_LA    		= 0.00e+00;
  real(rp), parameter :: xnDpITic_LA      		= 0.00e+00;
  real(rp), parameter :: xnDIC3ic_LA     		= 0.00e+00;
  real(rp), parameter :: xnDIC2ic_LA                    = 0.00e+00;
  real(rp), parameter :: xnDIFic_LA                     = 0.00E+00;
  real(rp), parameter :: xnDC3ic_LA                     = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_LA                     = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_LA                     = 0.00e+00; 
  real(rp), parameter :: xnDOic_LA                      = 0.00e+00; 
  real(rp), parameter :: xnDISic_LA                     = 0.00e+00; 
  real(rp), parameter :: xnDITic_LA                     = 0.00e+00; 


! LAA
  real(rp), parameter :: Vi_courte_LAA   		= -0.8401640587E+02;
  real(rp), parameter :: uic_LAA         		= 0.5576510930E-25;
  real(rp), parameter :: vic_LAA         		= 0.1000000000E+01;
  real(rp), parameter :: wic_LAA         		= 0.9993215779E+00;
  real(rp), parameter :: dic_LAA         		= 0.9591293562E-04;
  real(rp), parameter :: xnOic_LAA         		= 0.0000000000E+00;
  real(rp), parameter :: xric_LAA        		= 0.6767222950E-02;
  real(rp), parameter :: Naiic_LAA       		= 0.1184824114E+02;
  real(rp), parameter :: Kiic_LAA        		= 0.1384543025E+03;
  real(rp), parameter :: Carelic_LAA     		= 0.3977939214E+00;
  real(rp), parameter :: oiic_LAA        		= 0.9995547916E+00;
  real(rp), parameter :: uiic_LAA        		= 0.9882033654E+00;
  real(rp), parameter :: xnISic_LAA         		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_LAA         		= 0.0000000000E+00;
  real(rp), parameter :: fic_LAA         		= 0.8979814211E+00;
  real(rp), parameter :: xsic_LAA        		= 0.1951425375E-01;
  real(rp), parameter :: oaic_LAA        		= 0.2601346638E-01;
  real(rp), parameter :: uaic_LAA        		= 0.3702017446E-02;
  real(rp), parameter :: fCaic_LAA       		= 0.7499487528E+00;
  real(rp), parameter :: Caiic_LAA       		= 0.1164272383E-03;
  real(rp), parameter :: Caupic_LAA      		= 0.1201574598E+01;
  real(rp), parameter :: xnIC3ic_LAA                    = 0.00e+00;
  real(rp), parameter :: xnIC2ic_LAA                    = 0.00e+00;
  real(rp), parameter :: xnIFic_LAA                     = 0.00e+00;
  real(rp), parameter :: xnC3ic_LAA                     = 1.00e+00;
  real(rp), parameter :: xnC2ic_LAA                     = 0.00e+00;
  real(rp), parameter :: xnC1ic_LAA                     = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_LAA                  = 0.00e+00;
  real(rp), parameter :: xnDpIFic_LAA     		= 0.00e+00;
  real(rp), parameter :: xnDpC3ic_LAA     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_LAA      		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_LAA      		= 0.00e+00;
  real(rp), parameter :: xnDpOic_LAA      	        = 0.00e+00;
  real(rp), parameter :: xnDpISic_LAA    		= 0.00e+00;
  real(rp), parameter :: xnDpITic_LAA      		= 0.00e+00;
  real(rp), parameter :: xnDIC3ic_LAA     		= 0.00e+00;
  real(rp), parameter :: xnDIC2ic_LAA                   = 0.00e+00;
  real(rp), parameter :: xnDIFic_LAA                    = 0.00E+00;
  real(rp), parameter :: xnDC3ic_LAA                    = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_LAA                    = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_LAA                    = 0.00e+00; 
  real(rp), parameter :: xnDOic_LAA                     = 0.00e+00; 
  real(rp), parameter :: xnDISic_LAA                    = 0.00e+00; 
  real(rp), parameter :: xnDITic_LAA                    = 0.00e+00; 


  
 ! PV    
  real(rp), parameter :: Vi_courte_PV    		= -0.8380595736E+02;
  real(rp), parameter :: uic_PV          		= 0.4556871921E-25;
  real(rp), parameter :: vic_PV          		= 0.1000000000E+01;
  real(rp), parameter :: wic_PV          		= 0.9993131331E+00;
  real(rp), parameter :: dic_PV          		= 0.9846930324E-04;
  real(rp), parameter :: xnOic_PV          		= 0.0000000000E+00;
  real(rp), parameter :: xric_PV         		= 0.6254754847E-02;
  real(rp), parameter :: Naiic_PV        		= 0.1181117699E+02;
  real(rp), parameter :: Kiic_PV         		= 0.1385044636E+03;
  real(rp), parameter :: Carelic_PV      		= 0.3759856471E+00;
  real(rp), parameter :: oiic_PV         		= 0.9995367472E+00;
  real(rp), parameter :: uiic_PV         		= 0.9881987789E+00;
  real(rp), parameter :: xnISic_PV          		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_PV          		= 0.0000000000E+00;
  real(rp), parameter :: fic_PV          		= 0.9088873287E+00;
  real(rp), parameter :: xsic_PV         		= 0.1937464567E-01;
  real(rp), parameter :: oaic_PV         		= 0.2631922428E-01;
  real(rp), parameter :: uaic_PV         		= 0.3783764940E-02;
  real(rp), parameter :: fCaic_PV        		= 0.7583367168E+00;
  real(rp), parameter :: Caiic_PV        		= 0.1112858381E-03;
  real(rp), parameter :: Caupic_PV       		= 0.1141391732E+01;
  real(rp), parameter :: xnIC3ic_PV                     = 0.00e+00;
  real(rp), parameter :: xnIC2ic_PV                     = 0.00e+00;
  real(rp), parameter :: xnIFic_PV                      = 0.00e+00;
  real(rp), parameter :: xnC3ic_PV                      = 1.00e+00;
  real(rp), parameter :: xnC2ic_PV                      = 0.00e+00;
  real(rp), parameter :: xnC1ic_PV                      = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_PV                   = 0.00e+00;
  real(rp), parameter :: xnDpIFic_PV     		= 0.00e+00;
  real(rp), parameter :: xnDpC3ic_PV     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_PV      		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_PV      		= 0.00e+00;
  real(rp), parameter :: xnDpOic_PV     	        = 0.00e+00;
  real(rp), parameter :: xnDpISic_PV    		= 0.00e+00;
  real(rp), parameter :: xnDpITic_PV      		= 0.00e+00;
  real(rp), parameter :: xnDIC3ic_PV     		= 0.00e+00;
  real(rp), parameter :: xnDIC2ic_PV                    = 0.00e+00;
  real(rp), parameter :: xnDIFic_PV                     = 0.00e+00;
  real(rp), parameter :: xnDC3ic_PV                     = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_PV                     = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_PV                     = 0.00e+00; 
  real(rp), parameter :: xnDOic_PV                      = 0.00e+00; 
  real(rp), parameter :: xnDISic_PV                     = 0.00e+00; 
  real(rp), parameter :: xnDITic_PV                     = 0.00e+00; 

 
  ! MVR   
  real(rp), parameter :: Vi_courte_MVR    		= -0.8412150159E+02;
  real(rp), parameter :: uic_MVR          		= 0.4023455699E-25;
  real(rp), parameter :: vic_MVR          		= 0.1000000000E+01;
  real(rp), parameter :: wic_MVR          		= 0.9993257582E+00;
  real(rp), parameter :: dic_MVR          		= 0.9466024662E-04;
  real(rp), parameter :: xnOic_MVR          		= 0.0000000000E+00;
  real(rp), parameter :: xric_MVR         		= 0.5344114232E-02;
  real(rp), parameter :: Naiic_MVR        		= 0.1181926754E+02;
  real(rp), parameter :: Kiic_MVR         		= 0.1385054435E+03;
  real(rp), parameter :: Carelic_MVR      		= 0.3581646100E+00;
  real(rp), parameter :: oiic_MVR         		= 0.9995637572E+00;
  real(rp), parameter :: uiic_MVR         		= 0.9886117879E+00;
  real(rp), parameter :: xnISic_MVR          		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_MVR         		= 0.0000000000E+00;
  real(rp), parameter :: fic_MVR          		= 0.9167758540E+00;
  real(rp), parameter :: xsic_MVR         		= 0.1887127825E-01;
  real(rp), parameter :: oaic_MVR         		= 0.2586151028E-01;
  real(rp), parameter :: uaic_MVR         		= 0.3661706369E-02;
  real(rp), parameter :: fCaic_MVR        		= 0.7641512251E+00;
  real(rp), parameter :: Caiic_MVR        		= 0.1077868418E-03;
  real(rp), parameter :: Caupic_MVR       		= 0.1101228988E+01;
  real(rp), parameter :: xnIC3ic_MVR                    = 0.00e+00;
  real(rp), parameter :: xnIC2ic_MVR                    = 0.00e+00;
  real(rp), parameter :: xnIFic_MVR                     = 0.00e+00;
  real(rp), parameter :: xnC3ic_MVR                     = 1.00e+00;
  real(rp), parameter :: xnC2ic_MVR                     = 0.00e+00;
  real(rp), parameter :: xnC1ic_MVR                     = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_MVR                  = 0.00e+00;
  real(rp), parameter :: xnDpIFic_MVR   		= 0.00e+00;
  real(rp), parameter :: xnDpC3ic_MVR     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_MVR     		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_MVR     		= 0.00e+00;
  real(rp), parameter :: xnDpOic_MVR     	        = 0.00e+00;
  real(rp), parameter :: xnDpISic_MVR    		= 0.00e+00;
  real(rp), parameter :: xnDpITic_MVR     		= 0.00e+00;
  real(rp), parameter :: xnDIC3ic_MVR    		= 0.00e+00;
  real(rp), parameter :: xnDIC2ic_MVR                   = 0.00e+00;
  real(rp), parameter :: xnDIFic_MVR                    = 0.00e+00;
  real(rp), parameter :: xnDC3ic_MVR                    = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_MVR                    = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_MVR                    = 0.00e+00; 
  real(rp), parameter :: xnDOic_MVR                     = 0.00e+00; 
  real(rp), parameter :: xnDISic_MVR                    = 0.00e+00; 
  real(rp), parameter :: xnDITic_MVR                    = 0.00e+00; 
 
  ! BBla 
  real(rp), parameter :: Vi_courte_BBla   		= -0.8318182348E+02;
  real(rp), parameter :: uic_BBla         		= 0.1081134893E-21;
  real(rp), parameter :: vic_BBla         		= 0.1000000000E+01;
  real(rp), parameter :: wic_BBla         		= 0.9992874475E+00;
  real(rp), parameter :: dic_BBla         		= 0.1064683602E-03;
  real(rp), parameter :: xnOic_BBla         		= 0.0000000000E+00;
  real(rp), parameter :: xric_BBla        		= 0.2482783509E-01;
  real(rp), parameter :: Naiic_BBla       		= 0.1218243320E+02;
  real(rp), parameter :: Kiic_BBla        		= 0.1379577651E+03;
  real(rp), parameter :: Carelic_BBla     		= 0.7511009337E+00;
  real(rp), parameter :: oiic_BBla        		= 0.9994764614E+00;
  real(rp), parameter :: uiic_BBla        		= 0.9853091089E+00;
  real(rp), parameter :: xnISic_BBla         		= 0.0000000000E+00;
  real(rp), parameter :: xnDpIC3ic_BBla        		= 0.0000000000E+00;
  real(rp), parameter :: fic_BBla         		= 0.7784587323E+00;
  real(rp), parameter :: xsic_BBla        		= 0.2535760010E-01;
  real(rp), parameter :: oaic_BBla        		= 0.2725182252E-01;
  real(rp), parameter :: uaic_BBla        		= 0.4038346417E-02;
  real(rp), parameter :: fCaic_BBla       		= 0.6555514618E+00;
  real(rp), parameter :: Caiic_BBla       		= 0.1833468846E-03;
  real(rp), parameter :: Caupic_BBla      		= 0.1945175733E+01;
  real(rp), parameter :: xnIC3ic_BBla                   = 0.00e+00;
  real(rp), parameter :: xnIC2ic_BBla                   = 0.00e+00;
  real(rp), parameter :: xnIFic_BBla                    = 0.00e+00;
  real(rp), parameter :: xnC3ic_BBla                    = 1.00e+00;
  real(rp), parameter :: xnC2ic_BBla                    = 0.00e+00;
  real(rp), parameter :: xnC1ic_BBla                    = 0.00e+00;
  real(rp), parameter :: xnDpIC2ic_BBla                 = 0.00e+00;
  real(rp), parameter :: xnDpIFic_BBla      		= 0.00e+00;
  real(rp), parameter :: xnDpC3ic_BBla     		= 0.00e+00;
  real(rp), parameter :: xnDpC2ic_BBla       		= 0.00e+00;
  real(rp), parameter :: xnDpC1ic_BBla       		= 0.00e+00;
  real(rp), parameter :: xnDpOic_BBla      	        = 0.00e+00;
  real(rp), parameter :: xnDpISic_BBla   		= 0.00e+00;
  real(rp), parameter :: xnDpITic_BBla       		= 0.00e+00;
  real(rp), parameter :: xnDIC3ic_BBla      		= 0.00e+00;
  real(rp), parameter :: xnDIC2ic_BBla                  = 0.00e+00;
  real(rp), parameter :: xnDIFic_BBla                   = 0.00E+00;
  real(rp), parameter :: xnDC3ic_BBla                   = 0.00e+00; 
  real(rp), parameter :: xnDC2ic_BBla                   = 0.00e+00; 
  real(rp), parameter :: xnDC1ic_BBla                   = 0.00e+00; 
  real(rp), parameter :: xnDOic_BBla                    = 0.00e+00; 
  real(rp), parameter :: xnDISic_BBla                   = 0.00e+00; 
  real(rp), parameter :: xnDITic_BBla                   = 0.00e+00; 



!--------------------------------------------------------------------------
!-- STRUCTURE PARAMETERS ------------------------------------------------------
!------------------------------------------------------------------------------
  integer(ip), parameter :: Courte_mod  	= 23;  			! Courtemanche model (CRN)  
  integer(ip), parameter :: Courte_mod_RA_PM  	= 231;  		! Courtemanche model (CRN): RA/PM  
  integer(ip), parameter :: Courte_mod_CT_BBra 	= 232;  		! Courtemanche model (CRN): CT/BBra  
  integer(ip), parameter :: Courte_mod_TVR  	= 233;  		! Courtemanche model (CRN): TVR  
  integer(ip), parameter :: Courte_mod_RAA  	= 234;  		! Courtemanche model (CRN): RAA  
  integer(ip), parameter :: Courte_mod_LA  	= 235;  		! Courtemanche model (CRN): LA  
  integer(ip), parameter :: Courte_mod_LAA  	= 236;  		! Courtemanche model (CRN): LAA  
  integer(ip), parameter :: Courte_mod_PV  	= 237;  		! Courtemanche model (CRN): PV  
  integer(ip), parameter :: Courte_mod_MVR  	= 238;  		! Courtemanche model (CRN): MVR  
  integer(ip), parameter :: Courte_mod_BBla  	= 239;  		! Courtemanche model (CRN): BBla 
  integer(ip),parameter  :: nvar_courte 	= 43;  			! Number of state variables
  integer(ip),parameter  :: ncur_courte 	= 18;  			! Number of currents
  integer(ip),parameter  :: np_courte   	= 11;  			! Number of modifiable parameters  ! De momento pongo las conductancias máximas

!------------------------------------------------------------------------------
!-- ADAPTIVE TIME STEPING -----------------------------------------------------
!------------------------------------------------------------------------------
!.
  integer(ip), parameter    :: m_stpML   =    5
  real(rp),    parameter    :: m_dvdtML  =  1.0
