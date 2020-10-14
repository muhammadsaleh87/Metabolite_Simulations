%Saving parameters for running the simulation -- MGSaleh
% metabolite     - metabolite of interest
% echo time - a single or range of echo times
% editON - a single or range of editing frequencies. Only specify edit ON frequencies
% ppm_min/ppm_max -  the range of display of spectra
function MRS_opt = save_param_mega_hadamard(MRS_opt)

metabolite          = MRS_opt.metab;
spinSys             = MRS_opt.metab;

% General properties of the simulations
Npts                = 8192;     % number of spectral points
sw                  = 4000;     % spectral width [Hz]
Bfield              = 3;        % Philips magnetic field strength [Tesla]
lw                  = 2;        % linewidth of the output spectrum [Hz]
gamma               = 42577000; % gyromagnetic ratio
MRS_opt.Npts        = Npts;     % number of spectral points
MRS_opt.sw          = sw;     % spectral width [Hz]
MRS_opt.Bfield      = Bfield;        % Philips magnetic field strength [Tesla]
MRS_opt.lw          = lw;        % linewidth of the output spectrum [Hz]
MRS_opt.gamma       = gamma; % gyromagnetic ratio
MRS_opt.split_fit   = [-4.72454812433615e-39,3.49746593911910e-35,-1.17972529114192e-31,2.41257501068876e-28,-3.35334486561442e-25,3.36289930490681e-22,-2.51925567024359e-19,1.43926948019703e-16,-6.34362195286665e-14,2.16723851396576e-11,-5.73276757923450e-09,1.16606262695043e-06,-0.000179940454083095,0.0206037322492284,-1.68939288729085,93.4579199805247,-3115.57335510316,47271.5463178346];
MRS_opt.B1power_fit = [3.04275712634937e-44,-2.20316511739376e-40,7.44769386993364e-37,-1.56125144847875e-33,2.27359544767397e-30,-2.44184493594410e-27,2.00405005381190e-24,-1.28500696403995e-21,6.52644480995629e-19,-2.64612779316882e-16,8.59140538534298e-14,-2.23150800017893e-11,4.61306484949738e-09,-7.51638521375630e-07,9.50378264188365e-05,-0.00910554220880586,0.637054360250509,-30.6191846583449,901.589099210139,-12233.0077131620];

% Define the pulse waveforms here
refocWaveform       = 'univ_eddenrefo.pta';         % name of refocusing pulse waveform.
if strcmp(MRS_opt.seq, 'HERC')
    editWaveform1       = 'dl_Philips_4_58_1_90.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
    editWaveform2       = 'sl_univ_pulse.pta'; % name of 1st single editing pulse waveform. [4.58ppm]
    editWaveform3       = 'dl_Philips_4_18_1_90.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
    editWaveform4       = 'sl_univ_pulse.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
elseif strcmp(MRS_opt.seq, 'HERM')
    editWaveform1       = 'dl_Philips_univ_4_56_1_90.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
    editWaveform2       = 'sl_univ_pulse.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
    editWaveform3       = 'sl_univ_pulse.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
    editWaveform4       = 'sl_univ_pulse.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
else
    editWaveform1       = 'sl_univ_pulse.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
    editWaveform2       = 'sl_univ_pulse.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
end


% Define frequency parameters for editing targets
flipAngle               = 180;
centreFreq              = 3.0;                 % Center frequency of MR spectrum [ppm]
MRS_opt.edit_flipAngle  = 180; %flipAngle
MRS_opt.flipAngle       = flipAngle;                 % Flip angle degrees
MRS_opt.centreFreq      = centreFreq;                 % Center frequency of MR spectrum [ppm]
editOnFreq1             = MRS_opt.editON{1};%4.56; %4.58;                % Center frequency of 1st HERCULES experiment [ppm]
editOnFreq2             = MRS_opt.editON{2}; %freq_ppm(2);%7.5; %4.04;                % Center frequency of 2nd HERCULES experiment [ppm]
if ~strcmp(MRS_opt.seq, 'MEGA')
    editOnFreq3         = MRS_opt.editON{3};%4.56; %4.58;                % Center frequency of 1st HERCULES experiment [ppm]
    editOnFreq4         = MRS_opt.editON{4}; %freq_ppm(2);%7.5; %4.04;                % Center frequency of 2nd HERCULES experiment [ppm]
end
MRS_opt.editOnFreq1 = editOnFreq1;%4.56; %4.58;                % Center frequency of 1st HERCULES experiment [ppm]
MRS_opt.editOnFreq2 = editOnFreq2; %freq_ppm(2);%7.5; %4.04;                % Center frequency of 2nd HERCULES experiment [ppm]
if ~strcmp(MRS_opt.seq, 'MEGA')
    MRS_opt.editOnFreq3 = editOnFreq3;%4.56; %4.58;                % Center frequency of 1st HERCULES experiment [ppm]
    MRS_opt.editOnFreq4 = editOnFreq4; %freq_ppm(2);%7.5; %4.04;                % Center frequency of 2nd HERCULES experiment [ppm]
end

% Define pulse durations and flip angles specific for every vendor
refTp               = 7.0;              % duration of refocusing pulses [ms], set to zero for hard pulse
editTp1             = 20;               % duration of 1st editing pulse [ms]
editTp2             = 20;               % duration of 2nd editing pulse [ms]
if ~strcmp(MRS_opt.seq, 'MEGA')
    editTp3             = 20;               % duration of 2nd editing pulse [ms]
    editTp4             = 20;               % duration of 2nd editing pulse [ms]
end
editTp              = editTp1;
TE1                 = 6.55*2;           % TE1 [ms] (Use 6.96*2 for Philips Original and 6.55*2 for Universial/Siemens)
MRS_opt.refTp       = refTp;              % duration of refocusing pulses [ms], set to zero for hard pulse
MRS_opt.editTp1     = editTp1;               % duration of 1st editing pulse [ms]
MRS_opt.editTp2     = editTp2;               % duration of 2nd editing pulse [ms]
if ~strcmp(MRS_opt.seq, 'MEGA')
    MRS_opt.editTp3     = editTp3;               % duration of 2nd editing pulse [ms]
    MRS_opt.editTp4     = editTp4;               % duration of 2nd editing pulse [ms]
end
MRS_opt.TE1         = TE1;           % TE1 [ms] (Use 6.96*2 for Philips Original and 6.55*2 for Universial/Siemens)

% Define spatial resolution of simulation grid
MRS_opt.fovX        = 4.5; %size of the full simulation Field of View in the x-direction [cm]
MRS_opt.fovY        = 4.5; %size of the full simulation Field of View in the y-direction [cm]
MRS_opt.fovZ        = 4.5; %size of the full simulation Field of View in the y-direction [cm]
MRS_opt.thkX        = 3;    % slice thickness of x refocusing pulse [cm]
MRS_opt.thkY        = 3;    % slice thickness of y refocusing pulse [cm]
MRS_opt.thkZ        = 3;    % slice thickness of z excitation pulse [cm]
MRS_opt.nX          = 21;   % number of spatial points to simulate in x direction
MRS_opt.nY          = 21;   % number of spatial points to simulate in y direction
MRS_opt.nZ          = 21;   % number of spatial points to simulate in z direction
if MRS_opt.nX>1
    MRS_opt.x       = linspace(-MRS_opt.fovX/2,MRS_opt.fovX/2,MRS_opt.nX); %X positions to simulate [cm]
else
    MRS_opt.x=0;
end
if MRS_opt.nY>1
    MRS_opt.y       = linspace(-MRS_opt.fovY/2,MRS_opt.fovY/2,MRS_opt.nY); %X positions to simulate [cm]
else
    MRS_opt.y=0;
end
if MRS_opt.nZ>1
    MRS_opt.z       = linspace(-MRS_opt.fovZ/2,MRS_opt.fovZ/2,MRS_opt.nZ); %X positions to simulate [cm]
else
    MRS_opt.z=0;
end

% set up exact timing - based on 'normal' pulse set on Philips 3T - SH 07252019
tausA               = [TE1/2];               %middle 2nd EDITING to the start of readout
taus                = tausA;

% ************END of PARAMETERS**********************************


% ********************SET UP SIMULATION**********************************

% Load RF waveforms for excitation and refocusing pulses
%     excRF = io_loadRFwaveform(exciteWaveform,'exc',0); %DJL added %Come back later - SH 07252019
switch refocWaveform
    case 'gtst1203_sp.pta'
        refRF       = io_loadRFwaveform(refocWaveform,'ref',0);
        refRF.tbw   = 9.32; %BW99 (kHz) * dur (ms)
    case 'univ_eddenrefo.pta'
        refRF       = io_loadRFwaveform(refocWaveform,'ref',0);
        refRF.tbw   = 9.394; %BW99 (kHz) * dur (ms)
end
MRS_opt.refRF       = refRF;
% Load RF waveforms for editing pulses, defined up
if strcmp(MRS_opt.seq, 'HERC')
    editRF1             = io_loadRFwaveform(editWaveform4,'inv',0);
%     editRF1.tw1         = 1.8107;
    MRS_opt.editRF1     = editRF1;
    editRF2             = io_loadRFwaveform(editWaveform4,'inv',0);
    MRS_opt.editRF2     = editRF2;
    editRF3             = io_loadRFwaveform(editWaveform4,'inv',0);
%     editRF3.tw1         = 1.8107;
    MRS_opt.editRF3     = editRF3;
    editRF4             = io_loadRFwaveform(editWaveform4,'inv',0);
    MRS_opt.editRF4     = editRF4;
elseif strcmp(MRS_opt.seq, 'HERM')
    editRF1             = io_loadRFwaveform(editWaveform4,'inv',0);
%     editRF1.tw1         = 1.8107;
    MRS_opt.editRF1     = editRF1;
    editRF2             = io_loadRFwaveform(editWaveform4,'inv',0);
    MRS_opt.editRF2     = editRF2;
    editRF3             = io_loadRFwaveform(editWaveform4,'inv',0);
    MRS_opt.editRF3     = editRF3;
    editRF4             = io_loadRFwaveform(editWaveform4,'inv',0);
    MRS_opt.editRF4     = editRF4;
else
    editRF1             = io_loadRFwaveform(editWaveform1,'inv',0);
% % % %     editRF1.tw1         = editRF1.tw1*(1-0.07);
    MRS_opt.editRF1     = editRF1;
    editRF2             = io_loadRFwaveform(editWaveform2,'inv',0);
    MRS_opt.editRF2     = editRF2;
end

% Construct the editing pulses from the waveforms and defined
% frequencies
MRS_opt.editRFonA   = rf_freqshift(editRF1,editTp1,(centreFreq-editOnFreq1)*Bfield*gamma/1e6); %1.90 = GABA ON
MRS_opt.editRFonB   = rf_freqshift(editRF2,editTp2,(centreFreq-editOnFreq2)*Bfield*gamma/1e6); %7.5 = GABA OFF or MM supp
if ~strcmp(MRS_opt.seq, 'MEGA')
    MRS_opt.editRFonC   = rf_freqshift(editRF3,editTp3,(centreFreq-editOnFreq3)*Bfield*gamma/1e6); %7.5 = GABA OFF or MM supp
    MRS_opt.editRFonD   = rf_freqshift(editRF4,editTp4,(centreFreq-editOnFreq4)*Bfield*gamma/1e6); %7.5 = GABA OFF or MM supp
end
% HERCULES has the same editing pulse duration and timing for all sub-experiments. Valid for Mega-press as well:
MRS_opt.editTp      = editTp1;
MRS_opt.taus        = tausA;

% % Load the spin system definitions and pick the spin system of choice
load spinSystems
sys=eval(['sys' spinSys]);

% Set up gradients
MRS_opt.Gx          = (refRF.tbw/(refTp/1000))/(gamma*MRS_opt.thkX/10000); %[G/cm]
MRS_opt.Gy          = (refRF.tbw/(refTp/1000))/(gamma*MRS_opt.thkY/10000); %[G/cm]

for k=1:length(sys)
    sys(k).shifts   = sys(k).shifts-MRS_opt.centreFreq;
end
MRS_opt.sys=sys;
%Calculate new delays by subtracting the pulse durations from the taus
%vector;
delays         = zeros(size(taus));
delays(1)      = taus(1)-(refTp/2);
delays(2)      = 0; %taus(2)-((refTp+editTp)/2);
delays(3)      = 0; %taus(3)-((editTp+refTp)/2);
delays(4)      = 0; %taus(4)-((refTp+editTp)/2);
delays(5)      = 0; %taus(5)-(editTp/2);
MRS_opt.delays = delays;
%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian_mgs(sys,Bfield);
MRS_opt.H = H;
MRS_opt.d = d;
%Creating propagators for editing pulse
[MRS_opt.QoutONA]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonA,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
[MRS_opt.QoutONB]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonB,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
if ~strcmp(MRS_opt.seq, 'MEGA')
    [MRS_opt.QoutONC]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonC,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
    [MRS_opt.QoutOND]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonD,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
end
%Creating propagators for Refoc pulse ONLY in the X direction With an assumption x=y, i.e., the voxel is isotropic
parfor X=1:length(MRS_opt.x)  %Use this if you do have the MATLAB parallel processing toolbox
    [Qrefoc{X}]  = calc_shapedRF_propagator_refoc(MRS_opt.H,MRS_opt.refRF,MRS_opt.refTp,MRS_opt.flipAngle,    0,MRS_opt.y(X),MRS_opt.Gx);
end
MRS_opt.Qrefoc = Qrefoc; 

% %% Start the sequence
% d=sim_excite(d,MRS_opt.H,'x');                                          %EXCITE
% d=sim_apply_pfilter(d,MRS_opt.H,-1);
% d=sim_evolve(d,MRS_opt.H,MRS_opt.delays(1)/1000);
% parfor X=1:length(MRS_opt.x)
%     d_temp     = apply_propagator_refoc(d,MRS_opt.H,MRS_opt.Qrefoc{X});  %Refocusing in the X-direction
%     d_temp2{X} = d_temp;
% end
% %calculate the average density matrix (Doing this inside a separate for
% %loop because I couldn't figure out how to do this inside the parfor loop):
% d_A=struct([]);
% parfor X=1:length(MRS_opt.x)
%     d_A         =sim_dAdd(d_A,d_temp2{X});
% end
% d_A       = sim_apply_pfilter(d_A,MRS_opt.H,+1);                   %Apply p_filter
% MRS_opt.d = d_A; %From here, transfer to the fminsearch loop
% end



