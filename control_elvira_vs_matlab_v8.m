clc; clear; close all;

%% =====================================================================
% control_elvira_vs_matlab_v5.m
% ----------------------------------------------------------------------
% Tres comparaciones CONTROL
%
%   Fila 1 - Comp1: MATLAB INa_hh   vs  MATLAB INa_mkv
%   Fila 2 - Comp2: Fortran HH      vs  Fortran Markov
%   Fila 3 - Comp3: MAT(HH+MKV)    vs  FT(HH+MKV)
%
% COLORES (consistentes entre las 3 comparaciones):
%   Comp1/MAT  HH      -> azul oscuro  (solido)
%   Comp1/MAT  Markov  -> azul claro   (discontinuo)
%   Comp2/FT   HH      -> verde oscuro (solido)
%   Comp2/FT   Markov  -> verde claro  (discontinuo)
%   Comp3: usa los mismos 4 colores anteriores
%
% SINCRONIZACION MATLAB / FORTRAN:
%   MATLAB : N pulsos estimulados desde t=0.
%            En este script se descarta el primer pulso de MATLAB.
%   Fortran: FT_SKIP_BEATS latidos espontaneos + pulsos estimulados.
%            El script descarta los primeros FT_SKIP_BEATS antes de numerar
%            y ademas elimina el ultimo pulso de Fortran.
%
% Fig A  : resumen 3x3  (APD90 | Pico INa | max dV/dt)
% Fig B1-B3 : formas de onda del pulso seleccionado
% Fig C  : comparacion de 30 pulsos completos concatenados
%
% Dependencias: get_AP_biomarkers.m , get_last.m
%% =====================================================================

%% =========================
% 1) Ajustes de usuario
%% =========================
bcl_list = [1000 500 333 ];     % BCL(s) en ms

% Pulso estimulado a graficar (mismo indice para MATLAB y Fortran)
pulse_to_plot = 30;    % <--- CAMBIA AQUI

% Numero de pulsos completos a comparar en la figura tren
n_pulses_compare = 30;

% Pulsos de estabilizacion de Fortran a descartar
% (= PULSO_OFFSET/BCL del lanzabcl.sh, tipicamente 20)
FT_SKIP_BEATS = 20;    % <--- CAMBIA AQUI si modificas lanzabcl.sh

x_axis_mode       = 'pulse';   % 'pulse' o 'time'
ina_abs_threshold = 1.0;       % pA/pF

% ---- MATLAB .mat ----
%mat_file = 'D:\UPV\MANTA\Elvira\Resultados\2025\MATLAB\rigel\CRN_custom70_local_20260218_220355.mat';
%mat_file = 'D:\UPV\MANTA\Elvira\Resultados\2025\MATLAB\1702\CRN_custom50_36416674_20260216_161145.mat';
mat_file = 'D:\UPV\MANTA\Elvira\Resultados\2025\MATLAB\CRN_batch_20260202_165610.mat';

% ---- Fortran HH ----
FT_hh.folder   = 'D:\UPV\MANTA\Elvira\Resultados\2025\vfs\controlhh_RS';
FT_hh.stat_tpl = 'validCNR_controlRShh_p20amp20__232_ctrl%dbcl_stat.dat';
FT_hh.curr_tpl = 'validCNR_controlRShh_p20amp20__232_ctrl%dbcl_curr.dat';

% ---- Fortran Markov RS----
FT_mkv.folder   = 'D:\UPV\MANTA\Elvira\Resultados\2025\vfs\controlvf\markov_1102\codex\codexrs';
FT_mkv.stat_tpl = 'validCNR_controlcodexRS_p20amp20__232_ctrl%dbcl_stat.dat';
FT_mkv.curr_tpl = 'validCNR_controlcodexRS_p20amp20__232_ctrl%dbcl_curr.dat';

% 
% FT_mkv.folder   = 'D:\UPV\MANTA\Elvira\Resultados\2025\vfs\controlvf\markov_1102\codex\codexrs\v3';
% FT_mkv.stat_tpl = 'validCNR_controlcodex3RS_p20amp20__232_ctrl%dbcl_stat.dat';
% FT_mkv.curr_tpl = 'validCNR_controlcodex3RS_p20amp20__232_ctrl%dbcl_curr.dat';

% ---- Fortran Markov pAF----
% FT_mkv.folder   = 'D:\UPV\MANTA\Elvira\Resultados\2025\vfs\controlvf\markov_1102\codex';
% FT_mkv.stat_tpl = 'validCNR_controlcodex_p20amp20__232_ctrl%dbcl_stat.dat';
% FT_mkv.curr_tpl = 'validCNR_controlcodex_p20amp20__232_ctrl%dbcl_curr.dat';

%% =========================
% Paleta de colores
%% =========================
% Comp1 / MATLAB
c_mhh  = [0.08 0.27 0.60];    % azul oscuro  - MAT HH
c_mmkv = [0.42 0.68 0.96];    % azul claro   - MAT Markov
% Comp2 / Fortran
c_fhh  = [0.05 0.45 0.15];    % verde oscuro - FT HH
c_fmkv = [0.40 0.82 0.45];    % verde claro  - FT Markov

%% =========================
% 2) Cargar MATLAB (.mat)
%% =========================
if ~isfile(mat_file)
    error(['No existe el .mat:\n  %s\n' ...
           'Asegurate de haber relanzado CRN_batch_custom_70beats.m ' ...
           'y de apuntar al fichero correcto.'], mat_file);
end
S = load(mat_file);
if ~isfield(S,'results'), error('El .mat no contiene "results".'); end
results = S.results;

% Verificar numero de pulsos en el .mat
tmp_bcl  = sprintf('BCL_%dms', bcl_list(1));
if isfield(results, tmp_bcl) && isfield(results.(tmp_bcl),'INa_hh') && ...
   isfield(results.(tmp_bcl).INa_hh,'control')
    n_mat_pulsos = results.(tmp_bcl).INa_hh.control.settings.stim_num;
    fprintf('INFO: el .mat tiene %d pulsos simulados (stim_num).\n', n_mat_pulsos);
    if n_mat_pulsos < pulse_to_plot
        warning(['El .mat solo tiene %d pulsos pero pulse_to_plot=%d.\n' ...
                 'Relanza CRN_batch_custom_70beats.m y usa el nuevo .mat.\n' ...
                 'Por ahora se usara el ultimo pulso disponible.'], ...
                 n_mat_pulsos, pulse_to_plot);
    end
end

%% =========================
% 3) Loop por BCL(s)
%% =========================
for bcl_ms = bcl_list(:).'

    fprintf('\n%s\n', repmat('=',1,80));
    fprintf('CONTROL | BCL=%d ms\n', bcl_ms);
    fprintf('%s\n', repmat('=',1,80));

    keyBCL = sprintf('BCL_%dms', bcl_ms);
    if ~isfield(results, keyBCL)
        error('No existe results.%s', keyBCL);
    end

    %% ---- MATLAB INa_hh ----
    simHH  = results.(keyBCL).INa_hh.control;
    tHH    = simHH.t(:);   VmHH  = simHH.vm(:);
    INaHH  = extract_INa(simHH.currents, numel(tHH));

    %% ---- MATLAB INa_mkv ----
    simMKV = results.(keyBCL).INa_mkv.control;
    tMKV   = simMKV.t(:);  VmMKV = simMKV.vm(:);
    INaMKV = extract_INa(simMKV.currents, numel(tMKV));

    %% ---- Fortran HH ----
    [tFhh,  VmFhh]  = read_stat(fullfile(FT_hh.folder,  sprintf(FT_hh.stat_tpl,  bcl_ms)));
    [tFhhc, INaFhh] = read_curr(fullfile(FT_hh.folder,  sprintf(FT_hh.curr_tpl,  bcl_ms)));

    %% ---- Fortran Markov ----
    [tFmkv,  VmFmkv]  = read_stat(fullfile(FT_mkv.folder, sprintf(FT_mkv.stat_tpl, bcl_ms)));
    [tFmkvc, INaFmkv] = read_curr(fullfile(FT_mkv.folder, sprintf(FT_mkv.curr_tpl, bcl_ms)));

    %% ---- Deteccion de latidos ----
    % MATLAB: segmentacion bloqueada al estimulo, descartando el primer pulso
    % Fortran: deteccion robusta por inicio de upstroke (dV/dt),
    %          descartando FT_SKIP_BEATS al inicio y el ultimo pulso al final.
    %          La ventana de Fortran empieza unos ms antes del upstroke para
    %          que APD90 se calcule bien.
    nStimHH  = get_stim_num_safe(simHH,  50);
    nStimMKV = get_stim_num_safe(simMKV, 50);

    beatsHH   = detect_beats_stimlocked(tHH,   VmHH,   bcl_ms, nStimHH,  1);
    beatsMKV  = detect_beats_stimlocked(tMKV,  VmMKV,  bcl_ms, nStimMKV, 1);

    beatsFhh  = detect_beats_upstroke(tFhh,  VmFhh,  bcl_ms, FT_SKIP_BEATS);
    beatsFmkv = detect_beats_upstroke(tFmkv, VmFmkv, bcl_ms, FT_SKIP_BEATS);

    beatsFhh  = drop_last_beat(beatsFhh);
    beatsFmkv = drop_last_beat(beatsFmkv);

    %% ---- Diagnostico de pulsos disponibles ----
    n_mhh  = beatsHH.n;
    n_mmkv = beatsMKV.n;
    n_fhh  = beatsFhh.n;
    n_fmkv = beatsFmkv.n;

    fprintf('\nPulsos estimulados detectados por fuente:\n');
    fprintf('  MAT INa_hh  : %d pulsos\n', n_mhh);
    fprintf('  MAT INa_mkv : %d pulsos\n', n_mmkv);
    fprintf('  FT  HH      : %d pulsos estimulados (latidos fisicos 1-%d descartados + ultimo eliminado)\n', n_fhh,  FT_SKIP_BEATS);
    fprintf('  FT  Markov  : %d pulsos estimulados (latidos fisicos 1-%d descartados + ultimo eliminado)\n', n_fmkv, FT_SKIP_BEATS);

    if n_mhh < pulse_to_plot || n_mmkv < pulse_to_plot
        fprintf('\n  *** ATENCION: MATLAB tiene menos pulsos (%d) que pulse_to_plot=%d ***\n', ...
                min(n_mhh,n_mmkv), pulse_to_plot);
        fprintf('  *** Relanza CRN_batch_custom_70beats.m y actualiza mat_file ***\n\n');
    end

    % Pulso solicitado clampeado por fuente
    p_mhh  = clamp_pulse(pulse_to_plot, n_mhh,  'MAT INa_hh');
    p_mmkv = clamp_pulse(pulse_to_plot, n_mmkv, 'MAT INa_mkv');
    p_fhh  = clamp_pulse(pulse_to_plot, n_fhh,  'FT  HH');
    p_fmkv = clamp_pulse(pulse_to_plot, n_fmkv, 'FT  Markov');

    fprintf('Pulso solicitado: %d  =>  MAT-HH=%d | MAT-MKV=%d | FT-HH=%d (fis.%d) | FT-MKV=%d (fis.%d)\n', ...
            pulse_to_plot, p_mhh, p_mmkv, ...
            p_fhh,  p_fhh+FT_SKIP_BEATS, ...
            p_fmkv, p_fmkv+FT_SKIP_BEATS);

    %% ---- Pulsos comunes para figura resumen ----
    n1 = min(n_mhh,  n_mmkv);
    n2 = min(n_fhh,  n_fmkv);
    n3 = min([n_mhh, n_mmkv, n_fhh, n_fmkv]);

    pm1a = min(p_mhh,  n1);
    pm1b = min(p_mmkv, n1);
    pm2a = min(p_fhh,  n2);
    pm2b = min(p_fmkv, n2);
    pm3  = min([p_mhh, p_mmkv, p_fhh, p_fmkv, n3]);

    %% ---- Metricas por pulso ----
    M1a = metrics_per_pulse(tHH,   VmHH,   tHH,    INaHH,   bcl_ms, beatsHH.s0,   beatsHH.dt,   n1);
    M1b = metrics_per_pulse(tMKV,  VmMKV,  tMKV,   INaMKV,  bcl_ms, beatsMKV.s0,  beatsMKV.dt,  n1);

    M2a = metrics_per_pulse(tFhh,  VmFhh,  tFhhc,  INaFhh,  bcl_ms, beatsFhh.s0,  beatsFhh.dt,  n2);
    M2b = metrics_per_pulse(tFmkv, VmFmkv, tFmkvc, INaFmkv, bcl_ms, beatsFmkv.s0, beatsFmkv.dt, n2);

    M3_mhh  = metrics_per_pulse(tHH,   VmHH,   tHH,    INaHH,   bcl_ms, beatsHH.s0,   beatsHH.dt,   n3);
    M3_mmkv = metrics_per_pulse(tMKV,  VmMKV,  tMKV,   INaMKV,  bcl_ms, beatsMKV.s0,  beatsMKV.dt,  n3);
    M3_fhh  = metrics_per_pulse(tFhh,  VmFhh,  tFhhc,  INaFhh,  bcl_ms, beatsFhh.s0,  beatsFhh.dt,  n3);
    M3_fmkv = metrics_per_pulse(tFmkv, VmFmkv, tFmkvc, INaFmkv, bcl_ms, beatsFmkv.s0, beatsFmkv.dt, n3);

    %% ---- Ejes X ----
    [x1, xlab] = make_xaxis(x_axis_mode, n1, beatsHH.s0,  bcl_ms);
    [x2, ~]    = make_xaxis(x_axis_mode, n2, beatsFhh.s0, bcl_ms);
    [x3, ~]    = make_xaxis(x_axis_mode, n3, beatsHH.s0,  bcl_ms);

    %% ===================================================================
    %  FIG A: Resumen 3x3
    %% ===================================================================
    ylabels  = {'APD90 (ms)', 'Pico INa (pA/pF)', 'max dV/dt (V/s)'};
    metricas = {'APD90', 'INaPeak', 'DVdtMax'};
    tit_fila = {'MATLAB HH vs MATLAB Markov', ...
                'Fortran HH vs Fortran Markov', ...
                'MATLAB (HH+MKV) vs Fortran (HH+MKV)'};

    figure('Name', sprintf('Resumen CONTROL 3x3 | BCL=%d ms', bcl_ms), ...
           'Color','w', 'Position',[40 40 1500 950]);

    for row = 1:3
        for col = 1:3
            subplot(3,3,(row-1)*3+col);
            met = metricas{col};

            switch row
                case 1
                    plot_2s(x1, M1a.(met), x1, M1b.(met), pm1a, pm1b, ...
                            c_mhh, c_mmkv, 'MAT INa\_hh', 'MAT INa\_mkv');
                case 2
                    plot_2s(x2, M2a.(met), x2, M2b.(met), pm2a, pm2b, ...
                            c_fhh, c_fmkv, 'FT HH', 'FT Markov');
                case 3
                    plot_4s(x3, M3_mhh.(met), M3_mmkv.(met), ...
                                M3_fhh.(met),  M3_fmkv.(met), pm3, ...
                            c_mhh, c_mmkv, c_fhh, c_fmkv);
            end

            if col == 2
                yline(-abs(ina_abs_threshold), 'k--', 'Umbral', ...
                      'LineWidth',1.0, 'HandleVisibility','off');
            end

            box on; grid on;
            xlabel(xlab, 'FontSize',8);
            if col == 1
                ylabel(sprintf('%s\n%s', tit_fila{row}, ylabels{col}), ...
                       'FontSize',8, 'FontWeight','bold');
            else
                ylabel(ylabels{col}, 'FontSize',9);
            end
            if row == 1
                title(metricas{col}, 'FontWeight','bold', 'FontSize',10);
            end
            legend('Location','best', 'FontSize',8);
        end
    end

    sgtitle(sprintf('Resumen CONTROL | BCL=%d ms | pulso estimulado marcado: %d', ...
                    bcl_ms, pulse_to_plot), 'FontSize',14, 'FontWeight','bold');

    %% ===================================================================
    %  FIG B1, B2, B3 - Formas de onda
    %% ===================================================================
    [tVmMHH,  VmWMHH]  = beat_window(tHH,   VmHH,   bcl_ms, beatsHH.s0,   beatsHH.dt,   p_mhh);
    [tInMHH,  InWMHH]  = beat_window(tHH,   INaHH,  bcl_ms, beatsHH.s0,   beatsHH.dt,   p_mhh);
    [tVmMMKV, VmWMMKV] = beat_window(tMKV,  VmMKV,  bcl_ms, beatsMKV.s0,  beatsMKV.dt,  p_mmkv);
    [tInMMKV, InWMMKV] = beat_window(tMKV,  INaMKV, bcl_ms, beatsMKV.s0,  beatsMKV.dt,  p_mmkv);
    [tVmFHH,  VmWFHH]  = beat_window(tFhh,  VmFhh,  bcl_ms, beatsFhh.s0,  beatsFhh.dt,  p_fhh);
    [tInFHH,  InWFHH]  = beat_window(tFhhc, INaFhh, bcl_ms, beatsFhh.s0,  beatsFhh.dt,  p_fhh);
    [tVmFMKV, VmWFMKV] = beat_window(tFmkv, VmFmkv, bcl_ms, beatsFmkv.s0, beatsFmkv.dt, p_fmkv);
    [tInFMKV, InWFMKV] = beat_window(tFmkvc,INaFmkv,bcl_ms, beatsFmkv.s0, beatsFmkv.dt, p_fmkv);

    fig_wave2( ...
        sprintf('Comp1 | MATLAB HH vs MATLAB Markov | BCL=%d ms | pulso estim. %d', bcl_ms, pulse_to_plot), ...
        tVmMHH, VmWMHH, tInMHH, InWMHH, ...
        tVmMMKV,VmWMMKV,tInMMKV,InWMMKV, ...
        'MAT INa\_hh', 'MAT INa\_mkv', c_mhh, c_mmkv, 1);

    fig_wave2( ...
        sprintf('Comp2 | Fortran HH vs Fortran Markov | BCL=%d ms | pulso estim. %d (lat.fis. %d)', ...
                bcl_ms, pulse_to_plot, pulse_to_plot+FT_SKIP_BEATS), ...
        tVmFHH, VmWFHH, tInFHH, InWFHH, ...
        tVmFMKV,VmWFMKV,tInFMKV,InWFMKV, ...
        'FT HH', 'FT Markov', c_fhh, c_fmkv, 2);

    fig_wave4( ...
        sprintf('Comp3 | MAT(HH+MKV) vs FT(HH+MKV) | BCL=%d ms | pulso estim. %d', ...
                bcl_ms, pulse_to_plot), ...
        tVmMHH, VmWMHH, tInMHH, InWMHH, ...
        tVmMMKV,VmWMMKV,tInMMKV,InWMMKV, ...
        tVmFHH, VmWFHH, tInFHH, InWFHH, ...
        tVmFMKV,VmWFMKV,tInFMKV,InWFMKV, ...
        c_mhh, c_mmkv, c_fhh, c_fmkv, 3);

    %% ===================================================================
    %  FIG C: Comparacion de 30 pulsos completos concatenados
    %% ===================================================================
    ntrain = min([n_pulses_compare, n_mhh, n_mmkv, n_fhh, n_fmkv]);

    [t30_mhh_vm,  y30_mhh_vm]  = multi_beat_window(tHH,   VmHH,   bcl_ms, beatsHH.s0,   beatsHH.dt,   1, ntrain);
    [t30_mmkv_vm, y30_mmkv_vm] = multi_beat_window(tMKV,  VmMKV,  bcl_ms, beatsMKV.s0,  beatsMKV.dt,  1, ntrain);
    [t30_fhh_vm,  y30_fhh_vm]  = multi_beat_window(tFhh,  VmFhh,  bcl_ms, beatsFhh.s0,  beatsFhh.dt,  1, ntrain);
    [t30_fmkv_vm, y30_fmkv_vm] = multi_beat_window(tFmkv, VmFmkv, bcl_ms, beatsFmkv.s0, beatsFmkv.dt, 1, ntrain);

    [t30_mhh_in,  y30_mhh_in]  = multi_beat_window(tHH,    INaHH,   bcl_ms, beatsHH.s0,   beatsHH.dt,   1, ntrain);
    [t30_mmkv_in, y30_mmkv_in] = multi_beat_window(tMKV,   INaMKV,  bcl_ms, beatsMKV.s0,  beatsMKV.dt,  1, ntrain);
    [t30_fhh_in,  y30_fhh_in]  = multi_beat_window(tFhhc,  INaFhh,  bcl_ms, beatsFhh.s0,  beatsFhh.dt,  1, ntrain);
    [t30_fmkv_in, y30_fmkv_in] = multi_beat_window(tFmkvc, INaFmkv, bcl_ms, beatsFmkv.s0, beatsFmkv.dt, 1, ntrain);

    fig_train4( ...
        sprintf('Comp4 | 30 pulsos completos | MAT vs FT | BCL=%d ms | n=%d', bcl_ms, ntrain), ...
        t30_mhh_vm,  y30_mhh_vm,  t30_mhh_in,  y30_mhh_in, ...
        t30_mmkv_vm, y30_mmkv_vm, t30_mmkv_in, y30_mmkv_in, ...
        t30_fhh_vm,  y30_fhh_vm,  t30_fhh_in,  y30_fhh_in, ...
        t30_fmkv_vm, y30_fmkv_vm, t30_fmkv_in, y30_fmkv_in, ...
        c_mhh, c_mmkv, c_fhh, c_fmkv, 4);

    %% ---- Reporte consola ----
    fprintf('\n-- Comp1: MAT-HH vs MAT-MKV --\n');
    report_ina('MAT INa_hh ', M1a.INaPeak, ina_abs_threshold);
    report_ina('MAT INa_mkv', M1b.INaPeak, ina_abs_threshold);
    fprintf('-- Comp2: FT-HH vs FT-MKV --\n');
    report_ina('FT HH      ', M2a.INaPeak, ina_abs_threshold);
    report_ina('FT Markov  ', M2b.INaPeak, ina_abs_threshold);
    fprintf('-- Comp3 --\n');
    report_ina('MAT HH     ', M3_mhh.INaPeak,  ina_abs_threshold);
    report_ina('MAT MKV    ', M3_mmkv.INaPeak, ina_abs_threshold);
    report_ina('FT  HH     ', M3_fhh.INaPeak,  ina_abs_threshold);
    report_ina('FT  MKV    ', M3_fmkv.INaPeak, ina_abs_threshold);

end % loop bcl

%% =========================================================================
%  FUNCIONES DE DIBUJO
%% =========================================================================

function plot_2s(xa,ya, xb,yb, pa,pb, cA,cB, labA,labB)
    plot(xa, ya, '-',  'LineWidth',1.8, 'Color',cA, 'DisplayName',labA); hold on;
    plot(xb, yb, '--', 'LineWidth',1.8, 'Color',cB, 'DisplayName',labB);
    mark_pt(xa,ya,pa,cA);
    mark_pt(xb,yb,pb,cB);
    hold off;
end

function plot_4s(x, ya,yb,yc,yd, p, cMHH,cMMKV,cFHH,cFMKV)
    plot(x, ya, '-',  'LineWidth',1.8, 'Color',cMHH,  'DisplayName','MAT INa\_hh');  hold on;
    plot(x, yb, '--', 'LineWidth',1.8, 'Color',cMMKV, 'DisplayName','MAT INa\_mkv');
    plot(x, yc, '-',  'LineWidth',1.8, 'Color',cFHH,  'DisplayName','FT HH');
    plot(x, yd, '--', 'LineWidth',1.8, 'Color',cFMKV, 'DisplayName','FT Markov');
    mark_pt(x,ya,p,cMHH);  mark_pt(x,yb,p,cMMKV);
    mark_pt(x,yc,p,cFHH);  mark_pt(x,yd,p,cFMKV);
    hold off;
end

function mark_pt(x,y,p,c)
    if p>=1 && p<=numel(y) && isfinite(y(p))
        plot(x(p),y(p),'o','Color',c,'LineWidth',2,'MarkerSize',7,'HandleVisibility','off');
    end
end

function fig_wave2(fig_name, tVA,VmA,tIA,InA, tVB,VmB,tIB,InB, labA,labB, cA,cB, num)
    figure('Name',fig_name,'Color','w','Position',[80+num*50,80+num*50,1000,840]);
    subplot(3,1,1);
    plot(tVA,VmA,'-', 'LineWidth',2,'Color',cA,'DisplayName',labA); hold on;
    plot(tVB,VmB,'--','LineWidth',2,'Color',cB,'DisplayName',labB);
    hold off; box on; grid on; ylabel('Vm (mV)','FontSize',11);
    title(fig_name,'FontWeight','bold','FontSize',9,'Interpreter','none');
    legend('Location','best','FontSize',10);

    subplot(3,1,2);
    plot(tIA,InA,'-', 'LineWidth',2,'Color',cA,'DisplayName',labA); hold on;
    plot(tIB,InB,'--','LineWidth',2,'Color',cB,'DisplayName',labB);
    hold off; box on; grid on; ylabel('INa (pA/pF)','FontSize',11);
    legend('Location','best','FontSize',10);

    subplot(3,1,3);
    [tdA,dvA]=dvdt_calc(tVA,VmA); [tdB,dvB]=dvdt_calc(tVB,VmB);
    plot(tdA,dvA,'-', 'LineWidth',2,'Color',cA,'DisplayName',labA); hold on;
    plot(tdB,dvB,'--','LineWidth',2,'Color',cB,'DisplayName',labB);
    hold off; box on; grid on;
    ylabel('dV/dt (V/s)','FontSize',11);
    xlabel('Tiempo dentro del latido (ms)','FontSize',11);
    legend('Location','best','FontSize',10);
end

function fig_wave4(fig_name, ...
        tVmMHH,VmMHH,tInMHH,InMHH, tVmMMKV,VmMMKV,tInMMKV,InMMKV, ...
        tVmFHH,VmFHH,tInFHH,InFHH,  tVmFMKV,VmFMKV,tInFMKV,InFMKV, ...
        cMHH,cMMKV,cFHH,cFMKV, num)
    figure('Name',fig_name,'Color','w','Position',[80+num*50,80+num*50,1000,840]);

    subplot(3,1,1);
    plot(tVmMHH, VmMHH, '-', 'LineWidth',2,'Color',cMHH, 'DisplayName','MAT INa\_hh');  hold on;
    plot(tVmMMKV,VmMMKV,'--','LineWidth',2,'Color',cMMKV,'DisplayName','MAT INa\_mkv');
    plot(tVmFHH, VmFHH, '-', 'LineWidth',2,'Color',cFHH, 'DisplayName','FT HH');
    plot(tVmFMKV,VmFMKV,'--','LineWidth',2,'Color',cFMKV,'DisplayName','FT Markov');
    hold off; box on; grid on; ylabel('Vm (mV)','FontSize',11);
    title(fig_name,'FontWeight','bold','FontSize',9,'Interpreter','none');
    legend('Location','best','FontSize',10);

    subplot(3,1,2);
    plot(tInMHH, InMHH, '-', 'LineWidth',2,'Color',cMHH, 'DisplayName','MAT INa\_hh');  hold on;
    plot(tInMMKV,InMMKV,'--','LineWidth',2,'Color',cMMKV,'DisplayName','MAT INa\_mkv');
    plot(tInFHH, InFHH, '-', 'LineWidth',2,'Color',cFHH, 'DisplayName','FT HH');
    plot(tInFMKV,InFMKV,'--','LineWidth',2,'Color',cFMKV,'DisplayName','FT Markov');
    hold off; box on; grid on; ylabel('INa (pA/pF)','FontSize',11);
    legend('Location','best','FontSize',10);

    subplot(3,1,3);
    [tdMHH, dvMHH] =dvdt_calc(tVmMHH, VmMHH);
    [tdMMKV,dvMMKV]=dvdt_calc(tVmMMKV,VmMMKV);
    [tdFHH, dvFHH] =dvdt_calc(tVmFHH, VmFHH);
    [tdFMKV,dvFMKV]=dvdt_calc(tVmFMKV,VmFMKV);
    plot(tdMHH, dvMHH, '-', 'LineWidth',2,'Color',cMHH, 'DisplayName','MAT INa\_hh');  hold on;
    plot(tdMMKV,dvMMKV,'--','LineWidth',2,'Color',cMMKV,'DisplayName','MAT INa\_mkv');
    plot(tdFHH, dvFHH, '-', 'LineWidth',2,'Color',cFHH, 'DisplayName','FT HH');
    plot(tdFMKV,dvFMKV,'--','LineWidth',2,'Color',cFMKV,'DisplayName','FT Markov');
    hold off; box on; grid on;
    ylabel('dV/dt (V/s)','FontSize',11);
    xlabel('Tiempo dentro del latido (ms)','FontSize',11);
    legend('Location','best','FontSize',10);
end

function fig_train4(fig_name, ...
        tVmMHH,VmMHH,tInMHH,InMHH, tVmMMKV,VmMMKV,tInMMKV,InMMKV, ...
        tVmFHH,VmFHH,tInFHH,InFHH,  tVmFMKV,VmFMKV,tInFMKV,InFMKV, ...
        cMHH,cMMKV,cFHH,cFMKV, num)

    figure('Name',fig_name,'Color','w','Position',[120+num*40,120+num*40,1200,860]);

    subplot(3,1,1);
    plot(tVmMHH, VmMHH, '-',  'LineWidth',1.6,'Color',cMHH, 'DisplayName','MAT INa\_hh'); hold on;
    plot(tVmMMKV,VmMMKV,'--', 'LineWidth',1.6,'Color',cMMKV,'DisplayName','MAT INa\_mkv');
    plot(tVmFHH, VmFHH, '-',  'LineWidth',1.6,'Color',cFHH, 'DisplayName','FT HH');
    plot(tVmFMKV,VmFMKV,'--', 'LineWidth',1.6,'Color',cFMKV,'DisplayName','FT Markov');
    hold off; box on; grid on;
    ylabel('Vm (mV)','FontSize',11);
    title(fig_name,'FontWeight','bold','FontSize',10,'Interpreter','none');
    legend('Location','best','FontSize',10);

    subplot(3,1,2);
    plot(tInMHH, InMHH, '-',  'LineWidth',1.6,'Color',cMHH, 'DisplayName','MAT INa\_hh'); hold on;
    plot(tInMMKV,InMMKV,'--', 'LineWidth',1.6,'Color',cMMKV,'DisplayName','MAT INa\_mkv');
    plot(tInFHH, InFHH, '-',  'LineWidth',1.6,'Color',cFHH, 'DisplayName','FT HH');
    plot(tInFMKV,InFMKV,'--', 'LineWidth',1.6,'Color',cFMKV,'DisplayName','FT Markov');
    hold off; box on; grid on;
    ylabel('INa (pA/pF)','FontSize',11);
    legend('Location','best','FontSize',10);

    subplot(3,1,3);
    [tdMHH, dvMHH]   = dvdt_calc(tVmMHH,  VmMHH);
    [tdMMKV, dvMMKV] = dvdt_calc(tVmMMKV, VmMMKV);
    [tdFHH, dvFHH]   = dvdt_calc(tVmFHH,  VmFHH);
    [tdFMKV, dvFMKV] = dvdt_calc(tVmFMKV, VmFMKV);

    plot(tdMHH,  dvMHH,  '-',  'LineWidth',1.6,'Color',cMHH, 'DisplayName','MAT INa\_hh'); hold on;
    plot(tdMMKV, dvMMKV, '--', 'LineWidth',1.6,'Color',cMMKV,'DisplayName','MAT INa\_mkv');
    plot(tdFHH,  dvFHH,  '-',  'LineWidth',1.6,'Color',cFHH, 'DisplayName','FT HH');
    plot(tdFMKV, dvFMKV, '--', 'LineWidth',1.6,'Color',cFMKV,'DisplayName','FT Markov');
    hold off; box on; grid on;
    ylabel('dV/dt (V/s)','FontSize',11);
    xlabel('Tiempo concatenado (ms)','FontSize',11);
    legend('Location','best','FontSize',10);
end

function [td,dvdt] = dvdt_calc(t,Vm)
    dv=diff(Vm); dt=diff(t);
    dvdt=dv./dt; dvdt(~isfinite(dvdt))=0;
    td=t(1:end-1);
end

%% =========================================================================
%  FUNCIONES DE DATOS
%% =========================================================================

function p = clamp_pulse(requested, n_avail, label)
    p = max(1, min(round(requested), n_avail));
    if p ~= requested
        fprintf('  AVISO: %s tiene %d pulsos; se usara %d en lugar de %d.\n', ...
                label, n_avail, p, requested);
    end
end

function [x,xlab] = make_xaxis(mode, n, s0, bcl_ms)
    if strcmpi(mode,'time')
        x    = s0+((1:n)'-1)*bcl_ms;
        xlab = 'Tiempo de estimulo (ms)';
    else
        x    = (1:n)';
        xlab = 'Pulso estimulado';
    end
end

function beats = detect_beats_stimlocked(t, Vm, bcl_ms, expected_beats, skip_beats)
% Segmentacion robusta usando el calendario de estimulos.
% Ideal para MATLAB cuando los pulsos se aplican desde t=0 cada BCL exacto.

    t  = t(:);
    Vm = Vm(:);

    dt_est = median(diff(t(1:min(end,2000))));
    if ~isfinite(dt_est) || dt_est <= 0
        error('dt invalido en detect_beats_stimlocked.');
    end

    s0 = 0;
    n_possible = floor((t(end) - s0 + 0.5*dt_est) / bcl_ms);
    n_total = min(expected_beats, n_possible);

    if n_total < 1
        error('No hay suficientes datos para segmentar pulsos MATLAB.');
    end

    valid = false(n_total,1);
    for k = 1:n_total
        t0 = s0 + (k-1)*bcl_ms;
        t1 = s0 +  k   *bcl_ms - dt_est;
        idx = (t >= t0) & (t <= t1);
        if nnz(idx) < 10
            continue;
        end
        vseg = Vm(idx);
        if (max(vseg) - min(vseg)) > 20
            valid(k) = true;
        end
    end

    first_valid = find(valid,1,'first');
    if isempty(first_valid)
        error('No se han validado pulsos en MATLAB con segmentacion por estimulo.');
    end

    valid_idx = find(valid);
    if skip_beats >= numel(valid_idx)
        error('skip_beats=%d >= pulsos validos=%d.', skip_beats, numel(valid_idx));
    end

    used_idx = valid_idx(skip_beats+1:end);

    beats = struct();
    beats.n_total = numel(valid_idx);
    beats.skip    = skip_beats;
    beats.n       = numel(used_idx);
    beats.s0      = s0 + (used_idx(1)-1)*bcl_ms;
    beats.dt      = dt_est;
    beats.k0      = used_idx(1);
    beats.kind    = 'stimlocked';
end

function beats = detect_beats_upstroke(t, Vm, bcl_ms, skip_beats)
% Deteccion robusta de latidos por inicio de upstroke usando dV/dt.
% Para APD90, la ventana del latido empieza un poco antes del upstroke.

    t  = t(:);
    Vm = Vm(:);

    dt_est = median(diff(t(1:min(end,2000))));
    dt     = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        error('dt invalido en detect_beats_upstroke.');
    end

    dVdt = gradient(Vm, t);
    dVdt(~isfinite(dVdt)) = 0;

    dvdt_thr = max(15, 0.18 * max(dVdt));
    cross_idx = find(dVdt(1:end-1) < dvdt_thr & dVdt(2:end) >= dvdt_thr) + 1;

    min_sep = round(0.60 * bcl_ms / dt);
    if isempty(cross_idx)
        error('No se detectaron upstrokes en Vm.');
    end

    keep = true(size(cross_idx));
    last = cross_idx(1);
    for i = 2:numel(cross_idx)
        if (cross_idx(i) - last) < min_sep
            keep(i) = false;
        else
            last = cross_idx(i);
        end
    end
    onsets = cross_idx(keep);

    good = false(size(onsets));
    wpre  = max(1, round(0.08*bcl_ms/dt));
    wpost = max(1, round(0.25*bcl_ms/dt));
    for i = 1:numel(onsets)
        i0 = max(1, onsets(i)-wpre);
        i1 = min(numel(Vm), onsets(i)+wpost);
        if (max(Vm(i0:i1)) - min(Vm(i0:i1))) > 20
            good(i) = true;
        end
    end
    onsets = onsets(good);

    if isempty(onsets)
        error('No quedaron latidos validos tras filtrar onsets.');
    end

    if skip_beats >= numel(onsets)
        error('skip_beats=%d >= latidos detectados=%d.', skip_beats, numel(onsets));
    end

    onsets_used = onsets(skip_beats+1:end);

    % ---- CLAVE PARA APD90: incluir linea base antes del ascenso ----
    pre_ms = min(15, 0.12*bcl_ms);
    s0 = t(onsets_used(1)) - pre_ms;
    s0 = max(s0, t(1));

    beats = struct();
    beats.n_total = numel(onsets);
    beats.skip    = skip_beats;
    beats.n       = numel(onsets_used);
    beats.s0      = s0;
    beats.dt      = dt_est;
    beats.kind    = 'upstroke';
end

function beats = drop_last_beat(beats)
% Elimina el ultimo pulso disponible sin mover s0
    if beats.n <= 1
        error('No se puede quitar el ultimo pulso: solo hay %d disponible(s).', beats.n);
    end
    beats.n = beats.n - 1;
end

function M = metrics_per_pulse(tVm,Vm,tIna,INa,bcl_ms,s0,dt_est,nB)
    APD90=nan(nB,1); DVdtMax=nan(nB,1); INaPeak=nan(nB,1);
    for k=1:nB
        t0=s0+(k-1)*bcl_ms; t1=s0+k*bcl_ms-dt_est;
        ivm=(tVm>=t0)&(tVm<=t1);
        if nnz(ivm)>=10
            bio=safe_bio(tVm(ivm)-t0, Vm(ivm));
            APD90(k)=bio.APD90; DVdtMax(k)=bio.mx_dvdt;
        end
        iina=(tIna>=t0)&(tIna<=t1);
        if nnz(iina)>=1, INaPeak(k)=min(INa(iina)); end
    end
    M=struct('APD90',APD90,'DVdtMax',DVdtMax,'INaPeak',INaPeak);
end

function [t_rel,x_seg] = beat_window(t,x,bcl_ms,s0,dt_est,k)
    t0=s0+(k-1)*bcl_ms; t1=s0+k*bcl_ms-dt_est;
    idx=(t>=t0)&(t<=t1);
    t_rel=t(idx)-t0; x_seg=x(idx);
end

function [t_rel, x_seg] = multi_beat_window(t, x, bcl_ms, s0, dt_est, k_first, nbeats)
    t0 = s0 + (k_first-1)*bcl_ms;
    t1 = s0 + (k_first+nbeats-1)*bcl_ms - dt_est;
    idx = (t >= t0) & (t <= t1);
    t_rel = t(idx) - t0;
    x_seg = x(idx);
end

function bio = safe_bio(t,v)
    bio=struct('RMP',nan,'APO',nan,'APA',nan,'mx_dvdt',nan,'APD90',nan,'APD50',nan,'APD30',nan);
    try, bio=get_AP_biomarkers(t,v); catch, end
end

function nstim = get_stim_num_safe(sim, default_n)
    nstim = default_n;
    try
        if isfield(sim,'settings') && isfield(sim.settings,'stim_num') ...
                && ~isempty(sim.settings.stim_num) && isfinite(sim.settings.stim_num)
            nstim = sim.settings.stim_num;
        end
    catch
    end
    nstim = max(1, round(nstim));
end

function INa = extract_INa(currents,N)
    if isstruct(currents)
        if     isfield(currents,'INa'), INa=currents.INa(:);
        elseif isfield(currents,'ina'), INa=currents.ina(:);
        else,                           INa=nan(N,1);
        end
    else
        INa=currents(:);
        if numel(INa)~=N, INa=nan(N,1); end
    end
end

function [t,Vm] = read_stat(fname)
    if ~isfile(fname), error('No existe stat: %s',fname); end
    A=load(fname); Vm=A(:,end-1); t=A(:,end); t=t(:); Vm=Vm(:);
end

function [t,INa] = read_curr(fname)
    if ~isfile(fname), error('No existe curr: %s',fname); end
    A=load(fname); INa=A(:,1); t=A(:,end); t=t(:); INa=INa(:);
end

function report_ina(tag,INaPeak,thr)
    lost=find(INaPeak>-abs(thr));
    if isempty(lost)
        fprintf('%s: INa NO perdida en ningun pulso.\n',tag);
    else
        fprintf('%s: perdida en pulsos %s | primero: %d\n',tag,mat2str(lost'),lost(1));
    end
end