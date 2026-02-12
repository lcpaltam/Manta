
function [t, Vm, Ca_i, StateVars, currents] = Courtemanche_main(settings)

if nargin == 0    % no initialized struct settings
    settings.marijose = 1;      % create a value to start the struct
end

settings = setDefaultSettings(settings);

% y0       = getInitialVector(settings);
y0 = importdata('y0_CRN_nSR_mkv_0p.dat');

tol = 1e-12;
options  = odeset('AbsTol', tol,'RelTol', tol);

StateVars   = [];
Ti          = [];

t_start = tic;

% display(settings)

if settings.comments == 1
    % verbose
    fprintf('Running simulation...  0%% completed');
end

% Get time (Ti) and state variables (StateVars)
for n = 1:settings.stim_num 

    if settings.comments == 1
        % verbose
        fprintf(repmat('\b', 1, length('100% completed')));
        fprintf('%3.0f%% completed', 100*(n/settings.stim_num))
    end

    % dY Integration
    [t Y] = ode15s('Courtemanche_model_v2', [settings.stim_offset (settings.stim_BCL + settings.stim_offset)], y0, options, settings);

    y0    = Y(end,:); % last row becomes initial vector for the following iteration

    Ti          = [Ti; t + (settings.stim_BCL*(n-1))];
    StateVars   = [StateVars; Y];

end

% Get currents and re-organize state variables
Cm = 100; % pF
for i = 1:length(Ti)

        [dY, IsJs] = Courtemanche_model_v2(Ti(i), StateVars(i,:), [], settings, []); 

        currents.INa(i)         = IsJs(1)./Cm;                 % [pA/pF]
        currents.IK1(i)         = IsJs(2)./Cm;                 % [pA/pF]
        currents.Ito(i)         = IsJs(3)./Cm;                 % [pA/pF]
        currents.IKur(i)        = IsJs(4)./Cm;                 % [pA/pF]
        currents.IKr(i)         = IsJs(5)./Cm;                 % [pA/pF]
        currents.IKs(i)         = IsJs(6)./Cm;                 % [pA/pF]
        currents.ICa_L(i)       = IsJs(7)./Cm;                 % [pA/pF]
        currents.INaK(i)        = IsJs(8)./Cm;                 % [pA/pF]
        currents.INaCa(i)       = IsJs(9)./Cm;                 % [pA/pF]
        currents.ICaP(i)        = IsJs(10)./Cm;                % [pA/pF]
        currents.IB_Na(i)       = IsJs(11)./Cm;                % [pA/pF]
        currents.IB_Ca(i)       = IsJs(12)./Cm;                % [pA/pF]
        currents.IB_K(i)        = IsJs(13)./Cm;                % [pA/pF]
        currents.Irel(i)        = IsJs(14)./Cm;                % [pA/pF]
        currents.Itr(i)         = IsJs(15)./Cm;                % [pA/pF]
        currents.Iup_leak(i)    = IsJs(16)./Cm;                % [pA/pF]
        currents.Iup(i)         = IsJs(17)./Cm;                % [pA/pF]
        currents.i_st(i)        = IsJs(18)./Cm;                % [pA/pF]

        % Re-organize state variables
        Vm(i)                   = StateVars(i,15);              % (mV) Membrane potential
        Ca_i(i)                 = StateVars(i,10)*1e+3;         % (µM) Ca2+ intracellular concentration

end

t = Ti';

t_elapsed = toc(t_start);
if settings.comments == 1
    fprintf('\n\t Simulation finished, took %.2f seconds \n', t_elapsed);
end

end



%%
% function y0 = getInitialVector(settings)    
% Set the initial conditions for state variables (21 state variables)

    % Initial Conditions

%     if ~isfield(settings, 'y0'), disp('Not able to set initial conditions'); end

%     y0 = importdata('y0_CRN_nSR_mkv_0p');

% end



%%
function settings = setDefaultSettings(settings)
% Set default settings for the simulation

    if nargin == 0    % no initialized struct settings
        settings.marijose = 1;      % create a value to start the struct
    end
    
    if ~isfield(settings, 'comments'),      settings.comments    = 1;                  end      % 0 if no comments, 1 if so

    % Stimulus settings from the paper [ Ref. CTR (1998) ]
	if ~isfield(settings, 'stim_num'),              settings.stim_num               = 12;                  end      % # of stimuli
	if ~isfield(settings, 'stim_offset'),           settings.stim_offset            = 0.0;                 end	    % offset for the stimuls (ms)
    if ~isfield(settings, 'stim_amp'),              settings.stim_amp               = -2000;               end      % Stimulation amplitude (pA/pF)
    if ~isfield(settings, 'stim_dur'),              settings.stim_dur               = 2;                   end      % Stimulation duration (ms)
    if ~isfield(settings, 'stim_BCL'),              settings.stim_BCL               = 1000;                end      % Basic Cycle Length (ms)
    if ~isfield(settings, 'Tsim'),                  settings.Tsim                   = 12000;               end      % Duration of simulation (ms)
           
    % Experiment configuration
    if ~isfield(settings, 'ACh'),                   settings.ACh                    = 0.0;                 end      % ACh concentration (UNIDADES)
    if ~isfield(settings, 'INa_model'),             settings.INa_model              = 'INa_hh';            end      % INa model
    if ~isfield(settings, 'drug_compound'),         settings.drug_compound          = 'control';           end      % Drug compound
    if ~isfield(settings, 'drug_concentration'),    settings.drug_concentration     = 0;                   end      % Drug concentration

    % Dofetilide parameters
    if ~isfield(settings, 'dofe_IKr_IC50'),         settings.dofe_IKr_IC50           = 0.01e+3; end %0.075e+3;            end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'dofe_IKr_nH  '),         settings.dofe_IKr_nH             = 1;                   end      % Jordi Llopis 
    if ~isfield(settings, 'dofe_INa_IC50'),         settings.dofe_INa_IC50           = 1460e+3;             end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'dofe_INa_nH  '),         settings.dofe_INa_nH             = 5.1;                 end      % Jordi Llopis 
    if ~isfield(settings, 'dofe_INaL_IC50'),        settings.dofe_INaL_IC50          = 837e+3;              end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'dofe_INaL_nH  '),        settings.dofe_INaL_nH            = 4.6;                 end      % Jordi Llopis 
    if ~isfield(settings, 'dofe_ICaL_IC50'),        settings.dofe_ICaL_IC50          = 2.3e+3;              end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'dofe_ICaL_nH '),         settings.dofe_ICaL_nH            = 5.4;                 end      % Jordi Llopis 
    if ~isfield(settings, 'dofe_IKs_IC50'),         settings.dofe_IKs_IC50           = 100e+3;              end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'dofe_IKs_nH  '),         settings.dofe_IKs_nH             = 1;                   end      % Jordi Llopis 

    % Flecainide parameters
    if ~isfield(settings, 'fleca_IKr_IC50'),        settings.fleca_IKr_IC50          = 1500;                end      % (nM) - Jordi Llopis / 1.6 manta
    if ~isfield(settings, 'fleca_IKr_nH  '),        settings.fleca_IKr_nH            = 0.88;                end      % Jordi Llopis / h igual a 1 en manta
    if ~isfield(settings, 'fleca_INa_IC50'),        settings.fleca_INa_IC50          = 5800;                end      % (nM) - Jordi Llopis /6.5 manta
    if ~isfield(settings, 'fleca_INa_nH  '),        settings.fleca_INa_nH            = 1;                   end      % Jordi Llopis 
    if ~isfield(settings, 'fleca_INaL_IC50'),       settings.fleca_INaL_IC50         = 18870;               end      % (nM) - Jordi Llopis /20 manta
    if ~isfield(settings, 'fleca_INaL_nH  '),       settings.fleca_INaL_nH           = 0.6;                 end      % Jordi Llopis / h igual a 0.6
    if ~isfield(settings, 'fleca_ICaL_IC50'),       settings.fleca_ICaL_IC50         = 26349.5;             end      % (nM) - Jordi Llopis /27.1 manta
    if ~isfield(settings, 'fleca_ICaL_nH '),        settings.fleca_ICaL_nH           = 1.185;               end      % Jordi Llopis 
    if ~isfield(settings, 'fleca_IKs_IC50'),        settings.fleca_IKs_IC50          = 20000;               end      % (nM) - Jordi Llopis / no hay en manta
    if ~isfield(settings, 'fleca_IKs_nH  '),        settings.fleca_IKs_nH            = 1;                   end      % Jordi Llopis 
    if ~isfield(settings, 'fleca_Ito_IC50 '),       settings.fleca_Ito_IC50          = 9266;                end      % (nM) - Jordi Llopis / 10 en manta
    if ~isfield(settings, 'fleca_Ito_nH   '),       settings.fleca_Ito_nH            = 0.7;                 end      % Jordi Llopis / h= 0.8 manta
    if ~isfield(settings, 'fleca_IKur_IC50 '),      settings.fleca_IKur_IC50         = 2.9e+3;              end      % (nM) - Jordi Llopis /igual manta
    if ~isfield(settings, 'fleca_IKur_nH   '),      settings.fleca_IKur_nH           = 1;                   end      % Jordi Llopis 
    
    % Amiodarone parameters
    if ~isfield(settings, 'amiod_IKr_IC50'),        settings.amiod_IKr_IC50         = 560;                 end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'amiod_IKr_nH  '),        settings.amiod_IKr_nH           = 1.318;               end      % Jordi Llopis 
    if ~isfield(settings, 'amiod_INa_IC50'),        settings.amiod_INa_IC50         = 4800;                end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'amiod_INa_nH  '),        settings.amiod_INa_nH           = 0.7;                 end      % Jordi Llopis 
    if ~isfield(settings, 'amiod_INaL_IC50'),       settings.amiod_INaL_IC50        = 9423;                end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'amiod_INaL_nH  '),       settings.amiod_INaL_nH          = 0.4;                 end      % Jordi Llopis 
    if ~isfield(settings, 'amiod_ICaL_IC50'),       settings.amiod_ICaL_IC50        = 1590.5;              end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'amiod_ICaL_nH '),        settings.amiod_ICaL_nH          = 0.645;               end      % Jordi Llopis 
    if ~isfield(settings, 'amiod_IKs_IC50'),        settings.amiod_IKs_IC50         = 1740;                end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'amiod_IKs_nH  '),        settings.amiod_IKs_nH           = 0.72;                end      % Jordi Llopis 
    if ~isfield(settings, 'amiod_Ito_IC50 '),       settings.amiod_Ito_IC50         = 3758;                end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'amiod_Ito_nH   '),       settings.amiod_Ito_nH           = 0.4;                 end      % Jordi Llopis 

    % Vernakalant parameters
    if ~isfield(settings, 'verna_IKr_IC50'),        settings.verna_IKr_IC50         = 20e+3;              end      % (nM) - MANTA
    if ~isfield(settings, 'verna_IKr_nH  '),        settings.verna_IKr_nH           = 1;                  end      % MANTA
    if ~isfield(settings, 'verna_INa_IC50'),        settings.verna_INa_IC50         = 90e+3;              end      % (nM) - MANTA
    if ~isfield(settings, 'verna_INa_nH  '),        settings.verna_INa_nH           = 1;                  end      % MANTA
    if ~isfield(settings, 'verna_INaL_IC50'),       settings.verna_INaL_IC50        = 14e+3;              end      % (nM) - MANTA
    if ~isfield(settings, 'verna_INaL_nH  '),       settings.verna_INaL_nH          = 1;                  end      % MANTA
    if ~isfield(settings, 'verna_ICaL_IC50'),       settings.verna_ICaL_IC50        = 84e+3;              end      % (nM) - MANTA
    if ~isfield(settings, 'verna_ICaL_nH '),        settings.verna_ICaL_nH          = 1;                  end      % MANTA
    if ~isfield(settings, 'verna_Ito_IC50 '),       settings.verna_Ito_IC50         = 15e+3;              end      % (nM) - MANTA
    if ~isfield(settings, 'verna_Ito_nH   '),       settings.verna_Ito_nH           = 1;                  end      % MANTA
    if ~isfield(settings, 'verna_IKur_IC50 '),      settings.verna_IKur_IC50        = 15e+3;              end      % (nM) - MANTA
    if ~isfield(settings, 'verna_IKur_nH   '),      settings.verna_IKur_nH          = 1;                  end      % MANTA

    % Ranolazine parameters
    if ~isfield(settings, 'rano_IKr_IC50'),         settings.rano_IKr_IC50          = 8.3e+3;                 end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'rano_IKr_nH  '),         settings.rano_IKr_nH            = 1;                      end      % Jordi Llopis 
    if ~isfield(settings, 'rano_INa_IC50'),         settings.rano_INa_IC50          = 83.7e+3;                end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'rano_INa_nH  '),         settings.rano_INa_nH            = 1.1;                    end      % Jordi Llopis 
    if ~isfield(settings, 'rano_INaL_IC50'),        settings.rano_INaL_IC50         = 5.95e+3;                end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'rano_INaL_nH  '),        settings.rano_INaL_nH           = 1;                      end      % Jordi Llopis 
    if ~isfield(settings, 'rano_ICaL_IC50'),        settings.rano_ICaL_IC50         = 6540e+3;                end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'rano_ICaL_nH '),         settings.rano_ICaL_nH           = 3.8;                    end      % Jordi Llopis 
    if ~isfield(settings, 'rano_IKs_IC50'),         settings.rano_IKs_IC50          = 345e+3;                 end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'rano_IKs_nH  '),         settings.rano_IKs_nH            = 1;                      end      % Jordi Llopis 
    if ~isfield(settings, 'rano_Ito_IC50 '),        settings.rano_Ito_IC50          = 36155e+3;               end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'rano_Ito_nH   '),        settings.rano_Ito_nH            = 1;                      end      % Jordi Llopis 

    % Dronedarone parameters
    if ~isfield(settings, 'drone_IKr_IC50'),         settings.drone_IKr_IC50          = 0.0508e+3;           end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'drone_IKr_nH  '),         settings.drone_IKr_nH            = 0.85;                end      % Jordi Llopis 
    if ~isfield(settings, 'drone_INa_IC50'),         settings.drone_INa_IC50          = 0.7e+3;              end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'drone_INa_nH  '),         settings.drone_INa_nH            = 1;                   end      % Jordi Llopis 
    if ~isfield(settings, 'drone_ICaL_IC50'),        settings.drone_ICaL_IC50         = 0.1778e+3;              end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'drone_ICaL_nH '),         settings.drone_ICaL_nH           = 1;                   end      % Jordi Llopis 
    if ~isfield(settings, 'drone_IKs_IC50'),         settings.drone_IKs_IC50          = 10e+3;               end      % (nM) - Jordi Llopis
    if ~isfield(settings, 'drone_IKs_nH  '),         settings.drone_IKs_nH            = 1;                   end      % Jordi Llopi
    
end
