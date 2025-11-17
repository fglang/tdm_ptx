%% SETUP
[fp, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(fp));

%% LOAD DEMO DATA
% demo data set for Duke voxel model contains:
%   VOPs: VOP matrices for all 16 channels, corresponds to Q^{VOP}_{full} from paper
%   VOPs_sw: VOP matrices for switching, corresponds to Q^{VOP}_{combined} from paper
%   dist: spatial coordinates of voxels withing mask, N_r x 3
%   mask: brain mask, Nx x Ny x Nz
%   s: simulated voxel- and channel-wise B1+ maps (complex, [T/V]), N_r x 16
%   stepSize: voxel size in mm
%   dB0: simulated off-resonance [Hz], N_r x 1
load("data.mat");

%% SETTINGS
nkt = 3; % number of kT points for optimization
pulsedur = 0.2e-3; % [s], duration of a single block-pulse
TR = 10e-3; % [s], repetition time of sequence, needed for duty cycle scaling of SAR
Nc = size(VOPs,1); % number of channels in "full" case
SARLIM = 10; % [W/kg], local SAR limit
algo = 'interior-point'; % fmincon optimizer, interior-point or active-set are reasonable here

targFA = 5; % [deg], target flip angle
target = get_constant_target(targFA, mask);

%% EVAL: CP MODE
p_CP = vec(repmat(getCPmode(8),1,2).'); % CP mode pulse vector
ACP = createAmat(s,dist,[0,0,0],pulsedur); % system matrix
p_CP_scaled = p_CP .* targFA / mean(abs(calcFA(p_CP, ACP)),"all","omitnan"); % scale to achieve target FA on average
[FA_CP, nrmse_CP, power_CP, SAR_VOP_CP] =...
    eval_pulse(p_CP_scaled,ACP,mask,targFA,VOPs,[1],pulsedur,TR);

%% OPTIMIZE: kT points "full16"
% simultaneous 16ch transmission on both rows  (reference case)

% initialization at 0 for optimizer
k0 = zeros(nkt,3);
p0 = [];

% run joint optimization of pulse and k-space locations under VOP SAR constraints
[p, k_opt, nrmse_opt, FAmap_opt, SAR_VOP_opt] = ...
    jointopti_kt(nkt, s, dist, VOPs, dB0, mask, target, k0, p0, SARLIM, 0, [], 0, algo, pulsedur, TR);

%% CONVERT to multi8_conv
% rearrange previously optimized 16ch pulse into 8ch multiplexing shape and evaluate
p_conv = zeros(2*nkt,Nc/2);
p_conv(1:2:end,:) = 2.*p(:,1:2:end); % upper row, double amplitude for multiplexing
p_conv(2:2:end,:) = 2.*p(:,2:2:end); % lower row, double amplitude for multiplexing

% to keep total pulse duration same as for "full16", we need to halve the sub-pulse duration
pulsedur_sw = 0.1e-3;

% switching mask: indicate state of the switch for each sub-pulse
swmask = repmat([1,2], 1, nkt); 

% re-arrange B1+ values into multiplexing form: N_r x Nc/2 x 2 (last axis according to switch state)
s_sw = permute(reshape(s, [], 2, 8), [1,3,2]); 

% switched system matrix
Asw = createAmat(s_sw, dist, repelem(k_opt,2,1), pulsedur_sw, swmask, dB0); % repeat k-space locations for multiplexing

% re-evaluate the pulse in multiplexing scenario
[FAmap_conv, nrmse_conv, power_conv, SAR_VOP_conv] = eval_pulse(p_conv, Asw, mask, targFA, VOPs_sw, swmask, pulsedur_sw, TR);

%% OPTIMIZE: kT points "multi8_opt"
% multiplexed 8ch transmission with dedicated optimization of RF pulses

% initialization at 0 for optimizer
k0 = zeros(nkt*2,3);
p0 = [];

% run joint optimization of pulse and k-space locations under VOP SAR constraints
[p_sw, k_opt_sw, nrmse_opt_sw, FAmap_opt_sw, SAR_VOP_opt_sw] = ...
    jointopti_kt(nkt, s_sw, dist, VOPs_sw, dB0, mask, target, k0, p0, SARLIM, 1, swmask, 0, algo, pulsedur_sw, TR);


%% PLOT
PLOTDATA = {{FA_CP, nrmse_CP, SAR_VOP_CP, 'CP'},...
            {FAmap_opt, nrmse_opt, SAR_VOP_opt, 'full16'},...
            {FAmap_conv, nrmse_conv, SAR_VOP_conv, 'multi8_conv'},...
            {FAmap_opt_sw, nrmse_opt_sw, SAR_VOP_opt_sw, 'multi8_opt'},...
    };

figure;
t=tiledlayout(4,4);
t.TileSpacing='tight';
ix = round(size(FA_CP)./[2,2,1.5]);
CAX = targFA.*[0.5, 1.5];
for jj=1:4
    nexttile(jj); % tra
    imagescc(abs(PLOTDATA{jj}{1}(:,:,ix(3))));
    caxis(CAX);
    title(PLOTDATA{jj}{4}, 'Interpreter','none')

    nexttile(jj+4); % cor
    imagescc(abs(PLOTDATA{jj}{1}(:,ix(2),:)));
    caxis(CAX);

    nexttile(jj+4+4); % sag
    imagescc(abs(PLOTDATA{jj}{1}(ix(1),:,:)));
    caxis(CAX);

    nexttile(jj+4+4+4); % FA histogram
    hi=histogram(abs(PLOTDATA{jj}{1}), 'Normalization','probability', 'BinWidth',0.05, 'DisplayStyle','bar');
    hi.EdgeColor="none";
    h=gca;
    h.YAxis.Visible='off';
    xlim([0,10]);
    xlabel('FA [deg]');
    title(sprintf('NRMSE=%.2f %%, SAR=%.2f W/kg', PLOTDATA{jj}{2}.*100, PLOTDATA{jj}{3}));
end
colormap turbo;
set(gcf,'color','w');