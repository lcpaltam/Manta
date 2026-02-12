%===============================================================================
% CellML file:   D:\Desktop\Models\courtemanche_ramirez_nattel_1998.cellml
% CellML model:  courtemanche_1998
% Date and time: 17/06/2015 at 22:53:35
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2015 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function [dY, ionCurrents] = Courtemanche_model_v2(time_instat, Y, ~,settings, ~)

%% Definition of diferential variables

dY = zeros(length(Y),1);

i_u      = 1;  i_v      = 2;  i_w      = 3;  i_d      = 4;  i_f_Ca   = 5;  i_f      = 6; 
i_h      = 7;  i_j      = 8;  i_m      = 9;  i_Ca_i   = 10; i_Ca_rel = 11; i_Ca_up  = 12; 
i_K_i    = 13; i_Na_i   = 14; i_V      = 15; i_xr     = 16; i_xs     = 17; i_oa     = 18; 
i_oi     = 19; i_ua     = 20; i_ui     = 21; 

% if Markov model for INa is included
i_IC3  = 22; i_IC2   = 23; i_IF    = 24; i_C3   = 25; i_C2   = 26; i_C1   = 27; i_O    = 28; 
i_IS   = 29; i_DpIC3 = 30; i_DpIC2 = 31; i_DpIF = 32; i_DpC3 = 33; i_DpC2 = 34; i_DpC1 = 35;
i_DpO  = 36; i_DpIS  = 37; i_DpIT  = 38; i_DIC3 = 39; i_DIC2 = 40; i_DIF  = 41; i_DC3  = 42; 
i_DC2  = 43; i_DC1   = 44; i_DO    = 45; i_DIS  = 46; i_DIT  = 47; 

%% Select drug or control
%   1    2     3    4      5   6    7     8     9 
% i_Na i_Ca_L i_to i_Kur i_Kr i_Ks i_K1 i_NaK i_NaCa

INa_Block    = 1; ICaL_Block   = 2; Ito_Block    = 3; IKur_Block   = 4; 
IKr_Block    = 5; IKs_Block    = 6; IK1_Block    = 7; INaK_Block   = 8; 
INaCa_Block  = 9; 

drug = ones(9,1);
verna_fact = 0; fleca_fact = 0;

switch settings.drug_compound
    
    case 'control'
        drug = ones(9,1);

    case 'dofetilide'
        D_dofe = settings.drug_concentration;
        
        % Block = 1/(1 + (D_dof/IC50)^nH)
        drug(IKr_Block)  = 1/(1 + ( D_dofe / settings.dofe_IKr_IC50  )^settings.dofe_IKr_nH);
        drug(INa_Block)  = 1/(1 + ( D_dofe / settings.dofe_INa_IC50  )^settings.dofe_INa_nH);
        drug(ICaL_Block) = 1/(1 + ( D_dofe / settings.dofe_ICaL_IC50 )^settings.dofe_ICaL_nH);
        drug(IKs_Block)  = 1/(1 + ( D_dofe / settings.dofe_IKs_IC50  )^settings.dofe_IKs_nH);
%         drug(Ito_Block)  = 1/(1 + ( D_dofe / settings.dofe_Ito_IC50  )^settings.dofe_Ito_nH);

    case 'flecainide'
        fleca_fact = 1;
        D_fleca = settings.drug_concentration;
        
        % Block = 1/(1 + (D_dof/IC50)^nH)
        drug(IKr_Block)  = 1/(1 + ( D_fleca / settings.fleca_IKr_IC50  )^settings.fleca_IKr_nH);
        drug(INa_Block)  = 1/(1 + ( D_fleca / settings.fleca_INa_IC50  )^settings.fleca_INa_nH);
        drug(ICaL_Block) = 1/(1 + ( D_fleca / settings.fleca_ICaL_IC50 )^settings.fleca_ICaL_nH);
        drug(IKs_Block)  = 1/(1 + ( D_fleca / settings.fleca_IKs_IC50  )^settings.fleca_IKs_nH);
        drug(Ito_Block)  = 1/(1 + ( D_fleca / settings.fleca_Ito_IC50  )^settings.fleca_Ito_nH);
        drug(IKur_Block) = 1/(1 + ( D_fleca / settings.fleca_IKur_IC50 )^settings.fleca_IKur_nH);

    case 'amiodarone'
        D_amiod = settings.drug_concentration;
        
        % Block = 1/(1 + (D_dof/IC50)^nH)
        drug(IKr_Block)  = 1/(1 + ( D_amiod / settings.amiod_IKr_IC50  )^settings.amiod_IKr_nH);
        drug(INa_Block)  = 1/(1 + ( D_amiod / settings.amiod_INa_IC50  )^settings.amiod_INa_nH);
        drug(ICaL_Block) = 1/(1 + ( D_amiod / settings.amiod_ICaL_IC50 )^settings.amiod_ICaL_nH);
        drug(IKs_Block)  = 1/(1 + ( D_amiod / settings.amiod_IKs_IC50  )^settings.amiod_IKs_nH);
        drug(Ito_Block)  = 1/(1 + ( D_amiod / settings.amiod_Ito_IC50  )^settings.amiod_Ito_nH); 

    case 'vernakalant'
        verna_fact = 1;
        D_verna = settings.drug_concentration;
        
        % Block = 1/(1 + (D_dof/IC50)^nH)
        drug(IKr_Block)  = 1/(1 + ( D_verna / settings.verna_IKr_IC50  )^settings.verna_IKr_nH);
        drug(INa_Block)  = 1/(1 + ( D_verna / settings.verna_INa_IC50  )^settings.verna_INa_nH);
        drug(ICaL_Block) = 1/(1 + ( D_verna / settings.verna_ICaL_IC50 )^settings.verna_ICaL_nH);
        drug(Ito_Block)  = 1/(1 + ( D_verna / settings.verna_Ito_IC50  )^settings.verna_Ito_nH); 
        drug(IKur_Block) = 1/(1 + ( D_verna / settings.verna_IKur_IC50 )^settings.verna_IKur_nH);

    case 'dronedarone'
        D_drone = settings.drug_concentration;
        
        % Block = 1/(1 + (D_drone/IC50)^nH)
        drug(IKr_Block)  = 1/(1 + ( D_drone / settings.drone_IKr_IC50  )^settings.drone_IKr_nH);
        drug(IKs_Block)  = 1/(1 + ( D_drone / settings.drone_IKs_IC50  )^settings.drone_IKs_nH);
        drug(INa_Block)  = 1/(1 + ( D_drone / settings.drone_INa_IC50  )^settings.drone_INa_nH);
        drug(ICaL_Block) = 1/(1 + ( D_drone / settings.drone_ICaL_IC50 )^settings.drone_ICaL_nH);

    case 'ranolazine'
        D_rano = settings.drug_concentration;
        
        % Block = 1/(1 + (D_dof/IC50)^nH)
        drug(IKr_Block)  = 1/(1 + ( D_rano / settings.rano_IKr_IC50  )^settings.rano_IKr_nH);
        drug(IKs_Block)  = 1/(1 + ( D_rano / settings.rano_IKs_IC50  )^settings.rano_IKs_nH);
        drug(INa_Block)  = 1/(1 + ( D_rano / settings.rano_INa_IC50  )^settings.rano_INa_nH);
        drug(ICaL_Block) = 1/(1 + ( D_rano / settings.rano_ICaL_IC50 )^settings.rano_ICaL_nH);
        drug(Ito_Block)  = 1/(1 + ( D_rano / settings.rano_Ito_IC50  )^settings.rano_Ito_nH);

    otherwise
        disp('\n Unknown drug compound')

end


%%
%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

% Universal constants
F = 96.4867;                % coulomb_per_millimole (in membrane)
R = 8.3143;                 % joule_per_mole_kelvin (in membrane)
T = 310.0;                  % kelvin (in membrane)  / MANTA es 295 lp

% Cell constants
Cm = 100.0;                 % picoF (in membrane) 

V_cell = 20100.0;           % micrometre_3 (in intracellular_ion_concentrations)
V_i = V_cell*0.68;          % micrometre_3 (in intracellular_ion_concentrations)
V_up = V_cell*0.0552;       % micrometre_3 (in intracellular_ion_concentrations)
V_rel = V_cell*0.0048;      % micrometre_3 (in intracellular_ion_concentrations)

% Channel constants
CMDN_max = 0.05;            % millimolar (in Ca_buffers)
CSQN_max = 10.0;            % millimolar (in Ca_buffers)
Km_CMDN = 0.00238;          % millimolar (in Ca_buffers)
Km_CSQN = 0.8;              % millimolar (in Ca_buffers)
Km_TRPN = 0.0005;           % millimolar (in Ca_buffers)
TRPN_max = 0.07;            % millimolar (in Ca_buffers)
Ca_up_max = 15.0;           % millimolar (in Ca_leak_current_by_the_NSR)
K_rel = 30.0;               % per_millisecond -1 (in Ca_release_current_from_JSR)
I_up_max = 0.005;           % millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
K_up = 0.00092;             % millimolar (in Ca_uptake_current_by_the_NSR)
I_NaCa_max = 1600.0;        % picoA_per_picoF (in Na_Ca_exchanger_current)
K_mCa = 1.38;               % millimolar (in Na_Ca_exchanger_current)
K_mNa = 87.5;               % millimolar (in Na_Ca_exchanger_current)
K_sat = 0.1;                % dimensionless (in Na_Ca_exchanger_current)
gamma = 0.35;               % dimensionless (in Na_Ca_exchanger_current)
g_B_Ca = 0.001131;          % nanoS_per_picoF (in background_currents)
g_B_K = 0.0;                % nanoS_per_picoF (in background_currents)
g_B_Na = 0.0006744375;      % nanoS_per_picoF (in background_currents)
g_Kr = .0294;               % nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
i_CaP_max = 0.275;          % picoA_per_picoF (in sarcolemmal_calcium_pump_current)
g_Ks = 0.12941176;          % nanoS_per_picoF (in slow_delayed_rectifier_K_current)
Km_K_o = 1.5;               % millimolar (in sodium_potassium_pump)
Km_Na_i = 10.0;             % millimolar (in sodium_potassium_pump)
i_NaK_max = 0.59933874;     % picoA_per_picoF (in sodium_potassium_pump)
Ca_o = 1.8;                 % millimolar (in standard_ionic_concentrations)
K_o = 5.4;                  % millimolar (in standard_ionic_concentrations)
Na_o = 140.0;               % millimolar (in standard_ionic_concentrations)
g_K1 = 0.09;                % nanoS_per_picoF (in time_independent_potassium_current)
tau_tr = 180.0;             % millisecond (in transfer_current_from_NSR_to_JSR) %equation n°70
K_Q10 = 3.0;                % dimensionless (in transient_outward_K_current)
g_to = 0.1652;              % nanoS_per_picoF (in transient_outward_K_current)
g_Ca_L = 0.12375;           % nanoS_per_picoF (in L_type_Ca_channel)


%-----------------------------------------------------------------------------
% Stimulus current (modificado Violeta)
%-----------------------------------------------------------------------------

if (time_instat < settings.stim_dur+settings.stim_offset)  && (time_instat >= settings.stim_offset) 
    i_st = settings.stim_amp;
else
    i_st = 0.0;
end


%-----------------------------------------------------------------------------
% Ion Currents 
%-----------------------------------------------------------------------------

%%% Equilibrium potentials
E_Na = R*T/F*log(Na_o/Y(i_Na_i)); %(Eq.28)
E_Ca = R*T/(2.0*F)*log(Ca_o/Y(i_Ca_i)); %(Eq.28) /2
E_K = R*T/F*log(K_o/Y(i_K_i)); %(Eq.28)

%%% Fast Na+ Current
switch settings.INa_model

    case 'INa_hh'
        
        g_Na = 7.8; % [nS/pF]
        INa = Cm*g_Na*drug(INa_Block)*Y(i_m)^3.0*Y(i_h)*Y(i_j)*(Y(i_V)-E_Na); %(Eq.29)

        %%%% m
        if (Y(i_V) == -47.13) %(Eq.30)
           alpha_m = 3.2;
        else
           alpha_m = 0.32*(Y(i_V)+47.13)/(1.0-exp(-0.1*(Y(i_V)+47.13)));
        end
        beta_m = 0.08*exp(-Y(i_V)/11.0); %(Eq.30)
        m_inf = alpha_m/(alpha_m+beta_m); %(Eq.34)
        tau_m = (alpha_m+beta_m)^(-1); %(Eq.34)

        %%%% h
        if (Y(i_V) < -40.0) %(Eq.31)
           alpha_h = 0.135*exp((Y(i_V)+80.0)/-6.8);
        else
           alpha_h = 0.0;
        end
        if (Y(i_V) < -40.0) %(Eq.31)
           beta_h = 3.56*exp(0.079*Y(i_V))+3.1e5*exp(0.35*Y(i_V));
        else
           beta_h = 1.0/(0.13*(1.0+exp((Y(i_V)+10.66)/-11.1)));
        end
        h_inf = alpha_h/(alpha_h+beta_h);
        tau_h = 1.0/(alpha_h+beta_h); %(Eq.34)

        %%%% j 
        if (Y(i_V) < -40.0) %(Eq.32)
           alpha_j = (-1.2714e5*exp(0.2444*Y(i_V))-3.474e-5*exp(-0.04391*Y(i_V)))*(Y(i_V)+37.78)/(1.0+exp(0.311*(Y(i_V)+79.23)));
        else
           alpha_j = 0.0;
        end
        if (Y(i_V) < -40.0) %(Eq.33)
           beta_j = 0.1212*exp(-0.01052*Y(i_V))/(1.0+exp(-0.1378*(Y(i_V)+40.14)));
        else
           beta_j = 0.3*exp(-2.535e-7*Y(i_V))/(1.0+exp(-0.1*(Y(i_V)+32.0)));
        end
        j_inf = alpha_j/(alpha_j+beta_j); %(Eq.34)
        tau_j = 1.0/(alpha_j+beta_j); %(Eq.34)

        dY(i_h)  = (h_inf-Y(i_h))/tau_h; %(Eq.20)
        dY(i_j)  = (j_inf-Y(i_j))/tau_j; %(Eq.20)
        dY(i_m)  = (m_inf-Y(i_m))/tau_m; %(Eq.20)


    case 'INa_mkv'
        
        % General factors (differ between courtemanche - schimdt) Los paramentros desde x1 - x25 fxx y vxx y pxx son tomados de MANTA?
        x1  = 1.1821516171073378;
        x2  = 1.8128950937709263;
        x3  = 0.55559223164198501;
        x4  = 1.2653811014166259;
        x5  = 0.379465955877468;
        x6  = 2.5485320500798077;
        x7  = 1.5828527326859647;
        x8  = 1.0237930757107652;
        x9  = 0.63628656115270266;
        x10 = 1.0802804689922176;
        x11 = 1.0992577563248505;
        x12 = 1.1001781318138173;
        x13 = 1.0751649122176028;
        x14 = 1.0400623701712255*0.1;
        x15 = 1.0221615376561712;
        x16 = 1.0349852098209555;
        x17 = 1.0095837445899951;
        x18 = 1.1952465;
        x19 = 1.811411333*1.05;
        x20 = 0.8;
        x21 = 0.9994445714285713;
        x22 = 1.0284021;
        x23 = 0.98904877;
        x24 = 1.01288225;
        x25 = 1.00647894;
        
        % Flecainide factors to optimize parameters (differ between courtemanche - schimdt)
        f1  = 15802.517082265076;
        f2  = 11212.49445307003;
        f3  = 1.0835911457972944;
        f4  = 4.7213138217409707;
        f5  = 1.1886854610200599;
        f6  = 1.1293915676022106;
        f7  = 0.98647121091698153;
        f8  = 0.54542349779755606;
        f9  = 1.020713680006784;
        f10 = 2.6609951303941308;
        f11 = 1.1411905523331658;
        f12 = 5.3540204782968566;
        f13 = 1.028853996176295;
        f14 = 1.0033797797467918;
        
        % Vernakalant factors to optimize parameters (differ between courtemanche - schimdt)
        v1  = 7170.4622644564224;
        v2  = 15352.068978379524;
        v3  = 0.78329295519833497;
        v4  = 6.2271179964787571;
        v5  = 22.10184390723893;
        v6  = 1.5107146709941404;
        v7  = 1.0935642293370962;
        v8  = 0.65910138578597266;
        v9  = 2.3700350948377;
        v10 = 6.2071593539114822;
        v11 = 0;
        v12 = 0;
        v13 = 0;
        v14 = 0;

        % Optimized Na+-channel Markov model parameters (formulation remains equal in courtemanche - schimdt)
        actshift = -15*x19;
        h1 = 2*x20;
        p1 = 8.5539*x1;
        p2 = 7.4392e-2*x2;
        p3 = 17.0*x3;
        p4 = 15.0*x4;
        p5 = 12.0*x5;
        p6 = 2.0373e-1*x6;
        p7 = 150*x7;
        p8 = 7.5215e-2*x8;
        p9 = 20.3*x9;
        p10 = 2.7574*x10;
        p11 = 5*x11;
        p12 = 4.7755e-1*x12;
        p13 = 10*x13;
        p14_new = -70*x21;
        p15_new = 3.5*x22;
        p16_new = 0.052 * 2.9*x23;
        p17_new = 0.132 * 1.9*x24;
        p18 = 13.370*x14;
        p19 = 43.749*x15;
        p20 = 3.4229e-2*x16;
        p21 = 1.7898e-2*x17;
        p01 = 41*x25;
        
        p22 = (fleca_fact*3.6324e-3*f1)   + (verna_fact*5.6974e-03*v1);
        p23 = (fleca_fact*2.6452*f2)      + (verna_fact*8.4559e+01*v2);
        p24 = (fleca_fact*5.7831e-5*f3)   + (verna_fact*6.3992e-07*v3);
        p25 = (fleca_fact*1.6689e-8*f4)   + (verna_fact*1.3511e+00*v4);
        p26 = (fleca_fact*2.6126e-01*f5)  + (verna_fact*1.3110e-01*v5);
        p27 = (fleca_fact*1.4847e3*f6)    + (verna_fact*6.7067e-06*v6);
        p28 = (fleca_fact*4.2385e+01*f7)  + (verna_fact*1.7084e-05*v7);
        p29 = (fleca_fact*1.7352e-6*f8)   + (verna_fact*1.9698e-05*v8);
        p30 = (fleca_fact*2.1181e+00*f9)  + (verna_fact*4.8477*v9);
        p31 = (fleca_fact*6.7505e-05*f10) + (verna_fact*3.2976*v10);
        p32 = (fleca_fact*2.4135*f11)     + (verna_fact*2.4135*v11);
        p33 = (fleca_fact*4.9001e-2*f12)  + (verna_fact*4.9001e-2*v12);
        p34 = (fleca_fact*1.0326e-03*f13) + (verna_fact*1.0326e-03*v13);
        p35 = (fleca_fact*2.1378e-02*f14) + (verna_fact*2.1378e-02*v14);
        
        % Equations for Markov INa formulation
        diffusion = ((fleca_fact*5500) + (verna_fact*500)); %in [1 / M ms] - % vernakalant value differ between courtemanche - schimdt (multiplicative factor in schmidt)
        %conc = 0; % concentration [mol]
        switch settings.drug_compound 
            case 'flecainide'
                conc = D_fleca*1e-9;
            case 'vernakalant'
                conc = D_verna*1e-9;
            case 'control'
                conc = 0.00;
        end
        pKa = (fleca_fact*9.3) + (verna_fact*5.4);
        kd_open = ((fleca_fact*11.2e-6) + (verna_fact*318e-6)) * exp(-0.7 * Y(i_V) * F / (R * T)); % vernakalant value differ between courtemanche - schimdt (multiplicative factor in schmidt)
        
        pH = 7.4;
        portion = 1/(1+10^(pH-pKa));
        conc_dplus = portion * conc;
        conc_d = (1-portion) * conc;
        %Equations markov rev ok lp 
        kon = conc_dplus * diffusion;
        kcon = kon;
        koff = kd_open * diffusion;
        kcoff = koff;
        
        k_on = conc_d * diffusion;
        k_off = ((fleca_fact*400e-6) + (verna_fact*400e-6)) * diffusion;    % vernakalant value differ between courtemanche - schimdt (multiplicative factor in schmidt)/ los valores son d2 schmidt lp
        ki_on = k_on / 2;
        ki_off = ((fleca_fact*5.4e-6) + (verna_fact*3.4e-6)) * diffusion;   % vernakalant value differ between courtemanche - schimdt (multiplicative factor in schmidt ) -d3
        kc_on = k_on / 2;
        kc_off = ((fleca_fact*800e-6) + (verna_fact*900e-6)) * diffusion;   % vernakalant value differ between courtemanche - schimdt (multiplicative factor in schmidt) -d4
        
        Tfactor = 1 / (3 ^ ((37 - (T - 273)) / 10.0)); % d5= 3 / Se saco de Courtemanche?? en schmidt es igual a 1 LP
        
        % Transition rates (ms-1)  
        a11 = Tfactor * p1 / (p2 * exp(-(Y(i_V) - actshift) / p3) + p6 * exp(-(Y(i_V) - actshift) / p7));
        a12 = Tfactor * p1 / (p2 * exp(-(Y(i_V) - actshift) / p4) + p6 * exp(-(Y(i_V) - actshift) / p7));
        a13 = Tfactor * p1 / (p2 * exp(-(Y(i_V) - actshift) / p5) + p6 * exp(-(Y(i_V) - actshift) / p7));
        b11 = Tfactor * p8  * exp(-(Y(i_V) - actshift) / p9);
        b12 = Tfactor * p10 * exp(-(Y(i_V) - actshift - p11) / p9);
        b13 = Tfactor * p12 * exp(-(Y(i_V) - actshift - p13) / p9);
        
        a3_ss = 1/(1+exp((Y(i_V) - p14_new) / p15_new));
        a3_tau = h1 + p01 * exp(p16_new * (Y(i_V) - p14_new)) / (1 + exp(p17_new * (Y(i_V) - p14_new)));
        a3 = Tfactor * a3_ss / a3_tau;
        b3 = Tfactor * (1 - a3_ss) / a3_tau;
        a2 = Tfactor * p18 * exp(Y(i_V) / p19);
        b2 = (a13*a2*a3) / (b13*b3);
        ax = p20 * a2;
        bx = p21 * a3;
        
        a13c = p22 * a13;
        b13c = 0; if kon>0, b13c = (b13 * kcon * koff * a13c) / (kon * kcoff * a13); end
        a13n = p23 * a13;
        b13n = 0; if k_on>0, b13n = (b13 * kc_on * a13n * k_off) / (kc_off * a13 * k_on); end
        ax1 = p24 * ax;
        bx1 = p25 * bx;
        ax2 = p26 * ax;
        bx2 = 0; if ki_on>0, bx2 = (bx * k_on * ax2 * ki_off) / (ax * ki_on * k_off); end
		a22 = p27 * a2;
        a_22 = p28 * a2;
        a33 = p31 * a3;
        b33 = p29 * b3;
        b22 = 0; if b13c>0, b22 = (a13c * a22 * a33) / (b13c * b33); end
        b_33 = p30 * b3;
        a_33 = 0; if ki_on>0, a_33 = (ki_off * a3 * kc_on * b_33) / (ki_on * kc_off * b3); end
        b_22 = 0; if b13n>0, b_22 = (a_33 * a13n * a_22) / (b_33 * b13n); end
      
       
        
        
        % DIT or DpIT ONLY in Lidocaine (Inactivated, Trapped state)
        a44 = p32 * a2;
        b44 = p33 * a3;
        a_44 = p34 * a2;
        b_44 = p35 * a2;
        
        
        % Differential equations  % Donde se describen ?¿ lp Supplemental figure 3, pp 77
        dY(i_IC3) = -Y(i_IC3) * (a11 + a3 + ki_on) + Y(i_IC2) * b11 + Y(i_C3) * b3 + ki_off * Y(i_DIC3);
        dY(i_IC2) = -Y(i_IC2) * (b11 + a3 + a12 + ki_on) + Y(i_IC3) * a11 + Y(i_IF) * b12 + Y(i_C2) * b3 + ki_off * Y(i_DIC2);
        dY(i_IF)  = -Y(i_IF)  * (b12 + a3 + b2 + ki_on) + Y(i_IC2) * a12 + Y(i_C1) * b3 + Y(i_O) * a2 + ki_off * Y(i_DIF);
        dY(i_C3)  = -Y(i_C3)  * (b3 + a11 + kcon + kc_on) + Y(i_IC3) * a3 + Y(i_C2) * b11 + Y(i_DpC3) * kcoff + Y(i_DC3) * kc_off;
        dY(i_C2)  = -Y(i_C2)  * (b11 + b3 + a12 + kcon + kc_on) + Y(i_C3) * a11 + Y(i_IC2) * a3 + Y(i_C1) * b12 + Y(i_DpC2) * kcoff + Y(i_DC2) * kc_off;
        dY(i_C1)  = -Y(i_C1)  * (b12 + b3 + a13 + kcon + kc_on) + Y(i_C2) * a12 + Y(i_IF) * a3 + Y(i_O) * b13 + Y(i_DpC1) * kcoff + Y(i_DC1) * kc_off;
        dY(i_O)   = -Y(i_O)   * (b13 + a2 + ax + kon + k_on) + Y(i_C1) * a13 + Y(i_IF) * b2 + Y(i_IS) * bx + Y(i_DpO) * koff + Y(i_DO) * k_off;
        dY(i_IS)  = -Y(i_IS)  * (bx + ki_on) + Y(i_O) * ax + Y(i_DIS) * ki_off;
        
        dY(i_DpIC3) = -Y(i_DpIC3) * (a33 + a11) + Y(i_DpIC2) * b11 + Y(i_DpC3) * b33;
        dY(i_DpIC2) = -Y(i_DpIC2) * (b11 + a33 + a12) + Y(i_DpIC3) * a11 + Y(i_DpIF) * b12 + Y(i_DpC2) * b33;
        dY(i_DpIF)  = -Y(i_DpIF)  * (b12 + a33 + b22 + a44) + Y(i_DpIC2) * a12 + Y(i_DpC1) * b33 + Y(i_DpO) * a22 + Y(i_DpIT) * b44;
        dY(i_DpC3)  = -Y(i_DpC3)  * (b33 + a11 + kcoff) + Y(i_DpIC3) * a33 + Y(i_DpC2) * b11 + Y(i_C3) * kcon;
        dY(i_DpC2)  = -Y(i_DpC2)  * (b11 + b33 + a12 + kcoff) + Y(i_DpC3) * a11 + Y(i_DpIC2) * a33 + Y(i_DpC1) * b12 + Y(i_C2) * kcon;
        dY(i_DpC1)  = -Y(i_DpC1)  * (b12 + b33 + a13c + kcoff) + Y(i_DpC2) * a12 + Y(i_DpIF) * a33 + Y(i_DpO) * b13c + Y(i_C1) * kcon;
        dY(i_DpO)   = -Y(i_DpO)   * (b13c + a22 + ax1 + koff) + Y(i_DpC1) * a13c + Y(i_DpIF) * b22 + Y(i_DpIS) * bx1 + Y(i_O) * kon;
        dY(i_DpIS)  = -Y(i_DpIS)  * (bx1) + Y(i_DpO) * ax1;
        dY(i_DpIT)  = -Y(i_DpIT)  * b44 + Y(i_DpIF) * a44;
        
        dY(i_DIC3) = -Y(i_DIC3) * (a_33 + a11 + ki_off) + Y(i_DIC2) * b11 + Y(i_DC3) * b_33 + ki_on * Y(i_IC3);
        dY(i_DIC2) = -Y(i_DIC2) * (b11 + a_33 + a12 + ki_off) + Y(i_DIC3) * a11 + Y(i_DIF) * b12 + Y(i_DC2) * b_33 + ki_on * Y(i_IC2);
        dY(i_DIF)  = -Y(i_DIF)  * (b12 + a_33 + b_22 + a_44 + ki_off) + Y(i_DIC2) * a12 + Y(i_DC1) * b_33 + Y(i_DO) * a_22 + Y(i_DIT) * b_44 + ki_on * Y(i_IF);
        dY(i_DC3)  = -Y(i_DC3)  * (b_33 + a11 + kc_off) + Y(i_DIC3) * a_33 + Y(i_DC2) * b11 + Y(i_C3) * kc_on;
        dY(i_DC2)  = -Y(i_DC2)  * (b11 + b_33 + a12 + kc_off) + Y(i_DC3) * a11 + Y(i_DIC2) * a_33  + Y(i_DC1) * b12 + Y(i_C2) * kc_on;
        dY(i_DC1)  = -Y(i_DC1)  * (b12 + b_33 + a13n + kc_off) + Y(i_DC2) * a12 + Y(i_DIF) * a_33 + Y(i_DO) * b13n + Y(i_C1) * kc_on;
        dY(i_DO)   = -Y(i_DO)   * (b13n + a_22 + ax2 + k_off) + Y(i_DC1) * a13n + Y(i_DIF) * b_22 + Y(i_DIS) * bx2 + Y(i_O) * k_on;
        dY(i_DIS)  = -Y(i_DIS)  * (bx2 + ki_off) + Y(i_DO) * ax2 + Y(i_IS) * ki_on;
        dY(i_DIT)  = -Y(i_DIT)  * b_44 + Y(i_DIF) * a_44;

        dY(i_h)  = 0.00;
        dY(i_j)  = 0.00;
        dY(i_m)  = 0.00;

        % currents
        g_Na = x18 * 23.5 * 0.31 * 0.5; % in [1/ms]
        INa = Cm * g_Na * Y(i_O) * (Y(i_V) - E_Na); % in [pA]

    otherwise
        disp('\n Unknown INa model')
end


%%% Time-Independent K+ Current
IK1 = Cm*g_K1*(Y(i_V)-E_K)/(1.0+exp(0.07*(Y(i_V)+80.0))); %(Eq.35)

%%% Transient Outward K+ Current
Ito = Cm*g_to*drug(Ito_Block)*Y(i_oa)^3.0*Y(i_oi)*(Y(i_V)-E_K); %(Eq.36)
%%%% oa
alpha_oa = 0.65*(exp((Y(i_V)+10)/-8.5)+exp((Y(i_V)-30)/-59.0))^(-1.0); %(Eq.37)
beta_oa = 0.65*(2.5+exp((Y(i_V)+82.0)/17.0))^(-1.0); %(Eq.37)
tau_oa = ((alpha_oa+beta_oa)^(-1.0))/K_Q10; %(Eq.38)
oa_infinity = (1.0+exp((Y(i_V)+20.47)/-17.54))^(-1.0); %(Eq.38)
%%%% oi
alpha_oi = (18.53+1.0*exp((Y(i_V)+113.7)/10.95))^(-1.0); %(Eq.39)
beta_oi = (35.56+1.0*exp((Y(i_V)+1.26)/-7.44))^(-1.0); %(Eq.39)
tau_oi = ((alpha_oi+beta_oi)^(-1.0))/K_Q10; %(Eq.40)
oi_infinity = (1.0+exp((Y(i_V)+43.1)/5.3))^(-1.0); %(Eq.40)

%%% Ultrarapid Delayed Rectifier K+ Current
g_Kur = (0.005+0.05/(1.0+exp((Y(i_V)-15.0)/-13.0))); %(Eq.42)
IKur = Cm*g_Kur*drug(IKur_Block)*Y(i_ua)^3.0*Y(i_ui)*(Y(i_V)-E_K); %(Eq.41)

%%%% ua
alpha_ua = 0.65*(exp((Y(i_V)+10.0)/-8.5)+exp((Y(i_V)-30.0)/-59.0))^(-1.0); %(Eq.43)
beta_ua = 0.65*(2.5+exp((Y(i_V)+82.0)/17.0))^(-1.0); %(Eq.43)
tau_ua = ((alpha_ua+beta_ua)^(-1.0))/K_Q10; %(Eq.44)
ua_infinity = (1.0+exp((Y(i_V)+30.03)/-9.6))^(-1.0); %(Eq.44)

%%%% ui
alpha_ui = (21.0+1.0*exp((Y(i_V)-185.0)/-28.0))^(-1.0); %(Eq.45)
beta_ui = exp((Y(i_V)-158.0)/(16.0));  %(Eq.45) same of 1.0/exp((Y(i_V)-158.0)/(-16.0));
tau_ui = (alpha_ui+beta_ui)^(-1.0)/K_Q10; %(Eq.46)
ui_infinity = (1.0+exp((Y(i_V)-99.45)/27.48))^(-1.0); %(Eq.46)

%%% Rapid Delayed Outward Rectifier K+ Current
IKr = Cm*g_Kr*drug(IKr_Block)*Y(i_xr)*(Y(i_V)-E_K)/(1.0+exp((Y(i_V)+15.0)/22.4)); %(Eq.47)

%%%% xr
if (abs(Y(i_V)+14.1) < 1.0e-10) %(Eq.48)
   alpha_xr = 0.0015;
else
   alpha_xr = 0.0003*(Y(i_V)+14.1)/(1.0-exp((Y(i_V)+14.1)/-5.0));
end
if (abs(Y(i_V)-3.3328) < 1.0e-10) %(Eq.48)
   beta_xr = 3.7836118e-4;
else
   beta_xr = 0.000073898*(Y(i_V)-3.3328)/(exp((Y(i_V)-3.3328)/5.1237)-1.0);
end
tau_xr = (alpha_xr+beta_xr)^(-1.0); %(Eq.49)
xr_infinity = (1.0+exp((Y(i_V)+14.1)/-6.5))^(-1.0); %(Eq.49)

%%% Slow Delayed Outward Rectifier K+ Current
IKs = Cm*g_Ks*drug(IKs_Block)*Y(i_xs)^2.0*(Y(i_V)-E_K); %(Eq.50)

%%%% xs
if (abs(Y(i_V)-19.9) < 1.0e-10) %(Eq.51)
   alpha_xs = 0.00068;
else
   alpha_xs = 0.00004*(Y(i_V)-19.9)/(1.0-exp((Y(i_V)-19.9)/-17.0));
end
if (abs(Y(i_V)-19.9) < 1.0e-10) %(Eq.51)
   beta_xs = 0.000315;
else
   beta_xs = 0.000035*(Y(i_V)-19.9)/(exp((Y(i_V)-19.9)/9.0)-1.0);
end
tau_xs = 0.5*(alpha_xs+beta_xs)^(-1.0); %(Eq.52)
xs_infinity = (1.0+exp((Y(i_V)-19.9)/-12.7))^(-0.5); %(Eq.52)

%%% I_KACh: Time-Independent, Acetylcholine-Activated K+ current (from Grandi)
% Based on experimental data from Koumi et al. 1994
I_KACh = Cm*1/(1+(0.03/settings.ACh))^2.1 * (0.08 + 0.04/(1 + exp((Y(i_V)+91)/12))) * (Y(i_V) - E_K);


%%% L-Type Ca2+ Current
ICa_L = Cm*g_Ca_L*drug(ICaL_Block)*Y(i_d)*Y(i_f)*Y(i_f_Ca)*(Y(i_V)-65.0); %(Eq.53)

%%%% d
if (abs(Y(i_V)+10.0) < 1.0e-10) %(Eq.54)
   tau_d = 4.579/(1.0+exp((Y(i_V)+10.0)/-6.24)); 
else
   tau_d = (1.0-exp((Y(i_V)+10.0)/-6.24))/(0.035*(Y(i_V)+10.0)*(1.0+exp((Y(i_V)+10.0)/-6.24)));
end
%tau_d = 0.6 + 30.0 / (1 + exp(((Y(i_V) + 40) / 15)^2.0)); % = MANTA =/ CTR (1998)

d_infinity = (1.0+exp((Y(i_V)+10.0)/(-8.0)))^(-1.0); %(Eq.54)

%%%% f
tau_f = 9.0*(0.0197*exp(-0.0337^2.0*(Y(i_V)+10.0)^2.0)+0.02)^(-1.0); %(Eq.55)
f_infinity = (1+exp((Y(i_V)+28)/6.9))^(-1); % (Eq.55) 
%f_infinity = exp(-(Y(i_V)+28.0)/6.9)/(1.0+exp(-(Y(i_V)+28.0)/6.9)); %= MANTA =/ CTR (1998)

%%%% fCa
tau_f_Ca = 2.0; %(Eq.56)
f_Ca_infinity = (1.0+Y(i_Ca_i)/0.00035)^(-1.0); %(Eq.56)

%%% Na+-K- Pump Current
sigma = 1.0/7.0*(exp(Na_o/67.3)-1.0); %(Eq.59)
f_NaK = (1.0+0.1245*exp(-0.1*F*Y(i_V)/(R*T))+0.0365*sigma*exp(-F*Y(i_V)/(R*T)))^(-1.0); %(Eq.58)
INaK = Cm*i_NaK_max*f_NaK*(1.0/(1.0+(Km_Na_i/Y(i_Na_i))^1.5))*(K_o/(K_o+Km_K_o)); %(Eq.57)

%%% Na+/Ca2+ Exchanger Current
INaCa = Cm*I_NaCa_max*(exp(gamma*F*Y(i_V)/(R*T))*Y(i_Na_i)^3.0*Ca_o-exp((gamma-1.0)*F*Y(i_V)/(R*T))*Na_o^3.0*Y(i_Ca_i))/((K_mNa^3.0+Na_o^3.0)*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*Y(i_V)*F/(R*T)))); %(Eq.60)

%%% Background Currents
Ib_Na = Cm*g_B_Na*(Y(i_V)-E_Na); %(Eq.62)
Ib_Ca = Cm*g_B_Ca*(Y(i_V)-E_Ca); %(Eq.61)
Ib_K = Cm*g_B_K*(Y(i_V)-E_K); %does not exist! g_B_K=0

%%% Ca2+ Pump Current
ICaP = Cm*i_CaP_max*Y(i_Ca_i)/(0.0005+Y(i_Ca_i));

%%% Ca2+ Release Current from JSR
Irel = K_rel*Y(i_u)^2.0*Y(i_v)*Y(i_w)*(Y(i_Ca_rel)-Y(i_Ca_i)); %(Eq.64)
Fn = 1e-12*V_rel*Irel - ((5e-13)/F)*(0.5*ICa_L-0.2*INaCa); %(Eq.68)
%Fn = 1.0e3*(1.0e-15*V_rel*Irel-1.0e-15/(2.0*F)*(0.5*ICa_L-0.2*INaCa));%(Eq.68) = MANTA =/ CTR (1998)

%%%% u
tau_u = 8.0; %(Eq.65)
u_infinity = 1/(1+exp(-(Fn-3.4175e-13)/13.67e-16));

%%%% v
tau_v = 1.91+2.09*(1.0+exp(-(Fn-3.4175e-13)/13.67e-16))^(-1.0); %(Eq.66)
v_infinity = 1.0-(1.0+exp(-(Fn-6.835e-14)/13.67e-16))^(-1.0); %(Eq.66)

%%%% w
%tau_w = 6.0*(1.0-exp(-(Y(i_V)-7.9)/5.0))/((1.0+0.3*exp(-(Y(i_V)-7.9)/5.0))*1.0*(Y(i_V)-7.9));
if (abs(Y(i_V)-7.9) < 1.0e-10) %(Eq.67)
   tau_w = 6.0*0.2/1.3; 
else
   tau_w = 6.0*(1.0-exp(-(Y(i_V)-7.9)/5.0))/((1.0+0.3*exp(-(Y(i_V)-7.9)/5.0))*1.0*(Y(i_V)-7.9));
end
w_infinity = 1.0-(1.0+exp(-(Y(i_V)-40.0)/17.0))^(-1.0); %(Eq.67)

%%% Transfer Current from NSR to JSR
Itr = (Y(i_Ca_up)-Y(i_Ca_rel))/tau_tr; %(Eq.69)

%%% Ca2+ Leak Current by the NSR
Iup_leak = I_up_max*Y(i_Ca_up)/Ca_up_max; %(Eq.72)

%%% Ca2+ Uptake Current by the NSR
Iup = I_up_max/(1.0+K_up/Y(i_Ca_i)); %(Eq.71)


%%% Total ion current (modificado Clara)
Iion=INa+IK1+Ito+IKur+IKr+IKs+Ib_Na+Ib_Ca+INaK+ICaP+INaCa+ICa_L+I_KACh;

ionCurrents = [INa, IK1, Ito, IKur, IKr, IKs, ICa_L, INaK, INaCa, ICaP,...
                Ib_Na, Ib_Ca, Ib_K, Irel, Itr, Iup_leak, Iup, i_st];


%-----------------------------------------------------------------------------
% Ca2+ Buffers
%-----------------------------------------------------------------------------

Ca_CMDN = CMDN_max*Y(i_Ca_i)/(Y(i_Ca_i)+Km_CMDN); %(Eq.73)
Ca_TRPN = TRPN_max*Y(i_Ca_i)/(Y(i_Ca_i)+Km_TRPN); %(Eq.74)
Ca_CSQN = CSQN_max*Y(i_Ca_rel)/(Y(i_Ca_rel)+Km_CSQN); %(Eq.75)


%-----------------------------------------------------------------------------
% Differential equations
%-----------------------------------------------------------------------------
%%% MYOCYTE (CNR model) =================================================


dY(i_V) = -(Iion+i_st)/Cm;

dY(i_u)  = (u_infinity-Y(i_u))/tau_u;
dY(i_v)  = (v_infinity-Y(i_v))/tau_v;
dY(i_w)  = (w_infinity-Y(i_w))/tau_w;
dY(i_d)  = (d_infinity-Y(i_d))/tau_d;
dY(i_f_Ca)  = (f_Ca_infinity-Y(i_f_Ca))/tau_f_Ca;
dY(i_f)  = (f_infinity-Y(i_f))/tau_f;

dY(i_xr) = (xr_infinity-Y(i_xr))/tau_xr;
dY(i_xs) = (xs_infinity-Y(i_xs))/tau_xs;
dY(i_oa) = (oa_infinity-Y(i_oa))/tau_oa;
dY(i_oi) = (oi_infinity-Y(i_oi))/tau_oi;
dY(i_ua) = (ua_infinity-Y(i_ua))/tau_ua;
dY(i_ui) = (ui_infinity-Y(i_ui))/tau_ui;

dY(i_Ca_up) = Iup-(Iup_leak+Itr*V_rel/V_up);
dY(i_Ca_rel) = (Itr-Irel)*(1.0+CSQN_max*Km_CSQN/(Y(i_Ca_rel)+Km_CSQN)^2.0)^(-1.0);

dY(i_Na_i) = (-3.0*INaK-(3.0*INaCa+Ib_Na+INa))/(V_i*F);
dY(i_K_i) = (2.0*INaK-(i_st+IK1+Ito+IKur+IKr+IKs+Ib_K))/(V_i*F);

B1 = (2.0*INaCa-(ICaP+ICa_L+Ib_Ca))/(2.0*V_i*F)+(V_up*(Iup_leak-Iup)+Irel*V_rel)/V_i;
B2 = 1.0+TRPN_max*Km_TRPN/(Y(i_Ca_i)+Km_TRPN)^2.0+CMDN_max*Km_CMDN/(Y(i_Ca_i)+Km_CMDN)^2.0;
dY(i_Ca_i) = B1/B2;



end
