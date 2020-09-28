%Running GABA MP experiment for a range of TEs
clc
close all
clear all

%Addpath FID-A
addpath('/Users/muhammad/Documents/MATLAB/FID-A/inputOutput')
addpath('/Users/muhammad/Documents/MATLAB/FID-A/rfPulseTools/mklassenTools')
addpath('/Users/muhammad/Documents/MATLAB/FID-A/rfPulseTools')
addpath('/Users/muhammad/Documents/MATLAB/FID-A/simulationTools/metabolites')
addpath('/Users/muhammad/Documents/MATLAB/FID-A/simulationTools')
addpath('/Users/muhammad/Documents/MATLAB/FID-A/processingTools')

%addpath depending on your folder path
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/GitHub/Simulations/Metabolite_simulations/Functions/edited_PRESS_sim')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/GitHub/Simulations/Metabolite_simulations/Functions/edited_PRESS_sim/Experim_funct')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/GitHub/Simulations/Metabolite_simulations/pulses')


% addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Run/PRESS_and_MEGA-PRESS/Area_TE_Calc/fast')
% addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Run/TE_EditON_FA_maps/MP_TE_FA_matfiles')
% addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Run/run_commands/necessary_run_codes')

% initialize metab parameters
tic
mega_or_hadam          = {'MEGA'}; %Options: HERC, HERM, MEGA
if strcmp(mega_or_hadam, 'HERC') %HERCULES
    metab                  = {'GABA','GSH'};
    TE                     = 80;
    ppm_min                = [2.8 2.75];
    ppm_max                = [3.2 3.15];
    A                      = (4.58+1.9)/2; %dual-lobe pulse
    B                      = 4.58; %single-lobe pulse
    C                      = (4.18+1.9)/2; %dual-lobe pulse
    D                      = 4.18; %single-lobe pulse
    Nmetab                 = length(metab);   %Only metabolite to be run
elseif strcmp(mega_or_hadam, 'HERM') %HERMES
    metab                  = {'GABA','GSH'};  %{'GABA','GSH','Lac'}; 
    TE                     = 80;
    ppm_min                = [2.8 2.75];% 
    ppm_max                = [3.2 3.15];% 
    A                      = (4.56+1.9)/2; %dual-lobe pulse
    B                      = 4.56; %single-lobe pulse
    C                      = 1.90; %single-lobe pulse
    D                      = 7.50; %single-lobe pulse
    Nmetab                 = length(metab);   %Only metabolite to be run
else %MEGA-PRESS
    metab                  = {'GABA'};
    TE                     = 80;
    Nmetab                 = 1;   %Only metabolite to be run
    ppm_min                = [2.80];
    ppm_max                = [3.20];
    A                      = 1.90; %single-lobe pulse
    B                      = 7.50; %single-lobe pulse
%     create a struct variable
    MRS_temp.seq           = mega_or_hadam; %MRS_temp.HERC_or_HERM    = HERC_or_HERM;
    MRS_temp.metab         = metab{Nmetab}; %{ii};
    MRS_temp.Nmetab        = Nmetab; %ii;
    MRS_temp.TEs           = num2cell(TE);
    MRS_temp.ppm_min       = ppm_min(Nmetab); %(ii);
    MRS_temp.ppm_max       = ppm_max(Nmetab); %(ii);
    MRS_temp.editON        = num2cell([A B]);
    MRS_temp               = save_param_mega_hadamard(MRS_temp); % This is the function you need to edit to change the simulation parameters (not specified above this line)
    MRS_opt                = MRS_temp; % Creating a struct variable with dimens >=1;
end

if ~strcmp(mega_or_hadam, 'MEGA')
    for ii = 1:Nmetab
        clear MRS_temp
        MRS_temp.seq             = mega_or_hadam; %MRS_temp.HERC_or_HERM    = HERC_or_HERM;
        MRS_temp.metab           = metab{ii};
        MRS_temp.Nmetab          = ii;
        MRS_temp.TEs             = num2cell(TE);
        MRS_temp.ppm_min         = ppm_min(ii); %(ii);
        MRS_temp.ppm_max         = ppm_max(ii); %(ii);
        MRS_temp.editON          = num2cell([A B C D]);
        MRS_temp                 = save_param_mega_hadamard(MRS_temp); % This is the function you need to edit to change the simulation parameters (not specified above this line)
        MRS_opt(ii)              = MRS_temp; % Creating a struct variable with dimens >=1;
    end
end

%Simulate
if ~strcmp(mega_or_hadam, 'MEGA')
    [outA, outB, outC, outD]  = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
else
    [outA, outB]              = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
end
time_sim              = datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS');
disp(time_sim)

