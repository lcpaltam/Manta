
clear; clc;
%addpath(genpath('D:\UPV2025\MANTA\MANTA_MATLAB_Violeta'))

%% general config
n = 50; bcl = 1000;
settings.stim_BCL = bcl;
settings.stim_num = n;
settings.stim_amp = -1.5*1200;
settings.Tsim     = n*bcl; 

% INa hh
settings.INa_model              = 'INa_hh';

%%% control
settings.drug_compound = 'control'; % 'control' 'flecainide' 'vernakalant'
settings.drug_concentration = 0;

tic;
[t_ctrl_hh, vm_ctrl_hh, cai_ctrl_hh, stateVars_ctrl_hh, currents_ctrl_hh] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_ctrl_hh, vm_ctrl_hh, bcl*(n - 1), bcl*n);
bmrks_ctrl_hh = get_AP_biomarkers(last_t,last_Vm);
hh_t_ctrl = toc;


%%% flecainida
settings.drug_compound = 'flecainide';

% 1.5 microM
settings.drug_concentration = 1e+3*1.5; % (1.5 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_fleca_x1_hh, vm_fleca_x1_hh, cai_fleca_x1_hh, stateVars_fleca_x1_hh, currents_fleca_x1_hh] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_fleca_x1_hh, vm_fleca_x1_hh, bcl*(n - 1), bcl*n);
bmrks_fleca_x1_hh = get_AP_biomarkers(last_t,last_Vm);
hh_t_fleca_x1 = toc;

% 3 microM
settings.drug_concentration = 1e+3*3; % (10 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_fleca_x2_hh, vm_fleca_x2_hh, cai_fleca_x2_hh, stateVars_fleca_x2_hh, currents_fleca_x2_hh] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_fleca_x2_hh, vm_fleca_x2_hh, bcl*(n - 1), bcl*n);
bmrks_fleca_x2_hh = get_AP_biomarkers(last_t,last_Vm);
hh_t_fleca_x2 = toc;

%%% vernakalant 
settings.drug_compound = 'vernakalant';

% 10 microM
settings.drug_concentration = 1e+3*10; % (10 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_verna_x1_hh, vm_verna_x1_hh, cai_verna_x1_hh, stateVars_verna_x1_hh, currents_verna_x1_hh] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_verna_x1_hh, vm_verna_x1_hh, bcl*(n - 1), bcl*n);
bmrks_verna_x1_hh = get_AP_biomarkers(last_t,last_Vm);
hh_t_verna_x1 = toc;

% 20 microM
settings.drug_concentration = 1e+3*20; % (20 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_verna_x2_hh, vm_verna_x2_hh, cai_verna_x2_hh, stateVars_verna_x2_hh, currents_verna_x2_hh] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_verna_x2_hh, vm_verna_x2_hh, bcl*(n - 1), bcl*n);
bmrks_verna_x2_hh = get_AP_biomarkers(last_t,last_Vm);
hh_t_verna_x2 = toc;




%% INa markov
settings.INa_model              = 'INa_mkv';

%%% control
settings.drug_compound = 'control'; % 'control' 'flecainide' 'vernakalant'
settings.drug_concentration = 0;

tic;
[t_ctrl_mkv, vm_ctrl_mkv, cai_ctrl_mkv, stateVars_ctrl_mkv, currents_ctrl_mkv] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_ctrl_mkv, vm_ctrl_mkv, bcl*(n - 1), bcl*n);
bmrks_ctrl_mkv = get_AP_biomarkers(last_t,last_Vm);
mkv_t_ctrl = toc;


%%% flecainida
settings.drug_compound = 'flecainide';

% 1.5 microM
settings.drug_concentration = 1e+3*1.5; % (1.5 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_fleca_x1_mkv, vm_fleca_x1_mkv, cai_fleca_x1_mkv, stateVars_fleca_x1_mkv, currents_fleca_x1_mkv] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_fleca_x1_mkv, vm_fleca_x1_mkv, bcl*(n - 1), bcl*n);
bmrks_fleca_x1_mkv = get_AP_biomarkers(last_t,last_Vm);
mkv_t_fleca_x1 = toc;

% 3 microM
settings.drug_concentration = 1e+3*3; % (10 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_fleca_x2_mkv, vm_fleca_x2_mkv, cai_fleca_x2_mkv, stateVars_fleca_x2_mkv, currents_fleca_x2_mkv] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_fleca_x2_mkv, vm_fleca_x2_mkv, bcl*(n - 1), bcl*n);
bmrks_fleca_x2_mkv = get_AP_biomarkers(last_t,last_Vm);
mkv_t_fleca_x2 = toc;

%%% vernakalant 
settings.drug_compound = 'vernakalant';

% 10 microM
settings.drug_concentration = 1e+3*10; % (10 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_verna_x1_mkv, vm_verna_x1_mkv, cai_verna_x1_mkv, stateVars_verna_x1_mkv, currents_verna_x1_mkv] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_verna_x1_mkv, vm_verna_x1_mkv, bcl*(n - 1), bcl*n);
bmrks_verna_x1_mkv = get_AP_biomarkers(last_t,last_Vm);
mkv_t_verna_x1 = toc;

% 20 microM
settings.drug_concentration = 1e+3*20; % (20 microM) TIENE QUE ESTAR EN nM!!!
tic;
[t_verna_x2_mkv, vm_verna_x2_mkv, cai_verna_x2_mkv, stateVars_verna_x2_mkv, currents_verna_x2_mkv] = Courtemanche_main(settings);
[last_t, last_Vm] = get_last(t_verna_x2_mkv, vm_verna_x2_mkv, bcl*(n - 1), bcl*n);
bmrks_verna_x2_mkv = get_AP_biomarkers(last_t,last_Vm);
mkv_t_verna_x2 = toc;


%% time consumption
% tiempos = [ hh_t_ctrl     mkv_t_ctrl;
%             hh_t_fleca_x1 mkv_t_fleca_x1;
%             hh_t_fleca_x2 mkv_t_fleca_x2;
%             hh_t_verna_x1 mkv_t_verna_x1;
%             hh_t_verna_x2 mkv_t_verna_x2]';
% 
% tiempos = array2table(tiempos, ...
%     'VariableNames',{'ctrl';'fleca_x1'; 'fleca_x2'; 'verna_x1'; 'verna_x2'}, ...
%     'RowNames',{'INa_hh';'INa_mkv'});

% tiempos de computo (segundos)
%           ctrl    fleca_x1    fleca_x2    verna_x1    verna_x2
% INa_hh    19.08	21.06	    21.15	    20.63	    20.58
% INA_mkv   19.82	29.18	    29.98	    31.03	    29.30

figure; hold on; set(gca,'FontName','Aptos Narrow','FontSize',14,'TickLabelInterpreter','none');
% plot(tiempos{1,:},'k-o','LineWidth',1.5,'DisplayName','INa hh')
% plot(tiempos{2,:},'m-o','LineWidth',1.5,'DisplayName','INa mkv')
xticks(1:5); xticklabels({'ctrl';'fleca_x1'; 'fleca_x2'; 'verna_x1'; 'verna_x2'})
xlim([0.8 5.2])
ylim([15 35])
ylabel('tiempo de computacion (segundos)')
legend('Box','off')


%% VISUALIZATION APs

figure; hold on; set(gca,'FontName','Aptos Narrow','FontSize',14,'TickLabelInterpreter','none');
%t = tiledlayout(2,1)
% 
% plot(t_ctrl_hh-(n-1)*bcl, vm_ctrl_hh,'k','LineWidth',1.5,'DisplayName','INa hh - ctrl')
% plot(t_fleca_x1_hh-(n-1)*bcl, vm_fleca_x1_hh,'LineWidth',1.5,'Color','#b2e061','DisplayName','INa hh - fleca_x1')
% plot(t_fleca_x2_hh-(n-1)*bcl, vm_fleca_x2_hh,'LineWidth',1.5,'Color','#a861e0','DisplayName','INa hh - fleca_x2')
% plot(t_verna_x1_hh-(n-1)*bcl, vm_verna_x1_hh,'LineWidth',1.5,'Color','#7d9d44','DisplayName','INa hh - verna_x1')
% plot(t_verna_x2_hh-(n-1)*bcl, vm_verna_x2_hh,'LineWidth',1.5,'Color','#653a86','DisplayName','INa hh - verna_x2')

plot(t_ctrl_mkv-(n-1)*bcl, vm_ctrl_mkv,'k--','LineWidth',1.5,'DisplayName','INa mkv - ctrl')
plot(t_fleca_x1_mkv-(n-1)*bcl, vm_fleca_x1_mkv,'LineStyle','--','LineWidth',1.5,'Color','#b2e061','DisplayName','INa mkv - fleca_x1')
plot(t_fleca_x2_mkv-(n-1)*bcl, vm_fleca_x2_mkv,'LineStyle','--','LineWidth',1.5,'Color','#a861e0','DisplayName','INa mkv - fleca_x2')
plot(t_verna_x1_mkv-(n-1)*bcl, vm_verna_x1_mkv,'LineStyle','--','LineWidth',1.5,'Color','#7d9d44','DisplayName','INa mkv - verna_x1')
plot(t_verna_x2_mkv-(n-1)*bcl, vm_verna_x2_mkv,'LineStyle','--','LineWidth',1.5,'Color','#653a86','DisplayName','INa mkv - verna_x2')


xlim([-10 600])
legend('Box','off','Interpreter','none')
xlabel('tiempo (ms)')
ylabel('potencial (mV)')
