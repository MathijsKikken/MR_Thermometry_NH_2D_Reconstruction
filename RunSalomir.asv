%% Salomir's near harmonic 2D reconstruction for MR Thermometry
% 
%  HIMM: Harmonic Initialized Model-based Multi-echo  
%  Main script to extract temperature and drift fields from artificially
%  generated data
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% The basis of this algorithm was proposed by Rares Salomir:
% Salomir et al., IEEE Transactions on Medical Imaging 2011
%
% Do not reproduce, distribute, or modify without proper citation according
% to license file

clear all; close all;

% Add the software to the directory
addpath(genpath('.\Code\MRT optimization\Salomir MRT\Single-echo GRE - near harmonic 2D reconstruction'))
addpath(genpath('.\Code\MaskingTools\make_mask_fcn_v2_file_exchange'))
addpath(genpath('.\Code\ISMRM Water Fat Toolbox'))

%% Define scanner and sequence parameters

b0 = 3;                     % Tesla
frequency_system = 128*1e6; % Frequency measured by used MR system [Hz]
prc = 1;                    % Direction of precession (scanner-dependent) 
                            % +1 = heating induces negative apparent freq shift
                            % -1 = heating induces positive apparent freq shift
alpha = -0.01;              % ppm/deg C;
gamma = 42.5778;            % MHz/T
rad2degC = 1/(2*pi*b0*alpha*gamma);  % Factor to convert radians to temperature


%% Load data

% Select if we want to create artificial data or use measured data
data = 'measured'; % 'measured' or 'artificial'

if strcmp(data,'measured')

    % Patient definition
    datafile = 'Data\';
    Pt       = 1;             % subject number

    % Load dynamic (baseline/heating) data

    % Scan definition
    sequence      = 'baseline';        % 'baseline' or 'heating'
    slice         = 1;                 % usually set to 1, except for M2D acquisitions
    num_stacks    = 1;                 % Number of transverse images acquired in a single dynamic
    discard_dyns  = 0;                 % Number of dynamics that are discarded  ( 2 )
    num_avgs      = 1;                 % Number of dynamics that will be averaged ( 5 )
    
    % Specify which echoes are used for the reconstruction, and the appendix of
    % the image and info file (after '_img' or _info', e.g. '_img_shim2')
    relevant_echoes = [2:2];           % 
    echoes_drift = [];                 % Specify whether drift estimation is to be performed on a select number of echoes 2:10
                                       % Set to '[]' if all echoes are used for both drift and temperature
    img_name_append = '';  % '' for no appended info, '_shim2' for second shim | strcat('_registered_slice',num2str(slice-1))
    info_name_append = ''; % '' for no appended info, '_shim2' for second shim
    
    % Acquire the data and corresponding tes
    [scan,tes,GRE_info,~,~,~,file_name] = get_data(img_name_append, info_name_append, sequence, datafile, Pt, slice, relevant_echoes, echoes_drift);
    imgsdyn = scan;
    tes_heat = tes;

elseif strcmp(data,'artificial')
    
    % Create artificial data

    % For now I created arteficial MRT data, to show the principles of the algorithm
    % The parameters Wlib, Flib, R2starlib, and dw0lib should normally be
    % generated from imgslib using the ISMRM Fat-Water toolbox. Here, the data is
    % arteficially generated...
    dim = 60;               % dimensionality of the arteficial data
    fatPercentage = 0.5;    % percentage of fat present in regions that have both water and fat
    hssig2 = 0.011;         % sigma for the hotspot
    order = 2;              % polynomial order of the arteficially generated drift field
    dw0shift = 0.2*pi;        % background polynomial frequency shift (will be replicated to all polynomial coeffs) (rad/sec)   
    num_fat_peaks = 6;                          % Choose an x-peak model, options are:
                                                % [ 3 | 4 | 5 | 6 | 7 | 9 ]                                 
    tes_heat = [0.02];      % echo time in seconds (single-echo experiment)
    noiselevel = 1e-6;      % add in noise (0.2 for SNR of 40)
    
    % Generate the hotspot used for heating
    [x,y] = meshgrid(-dim/2:dim/2-1);
    hotspot = exp(-hssig2*(x.^2 + y.^2)); % normalized hot spot
    
    maxtemps = linspace(0,1,25); % temperature increase (in degC) as the dynamic count increases
    for dyn = 1:length(maxtemps)
        hotspot_data(:,:,dyn) = hotspot*alpha*maxtemps(dyn)*b0*2*pi*gamma;
    end
    
    % Generate the drift field
    [yc,xc] = meshgrid(linspace(-1/2,1/2,dim));
    yc = yc(:);
    xc = xc(:);
    A = [];
    for yp = 0:order
        for xp = 0:(order-yp)
            A = [A (xc.^xp).*(yc.^yp)];
        end
    end
    coefficients = zeros(size(A,2),1) + dw0shift;  
    drift = reshape(A*coefficients,[dim dim]);    % polynomial background phase shift
    
    % Generate the true water and fat images
    [x,y] = meshgrid(-dim/2:dim/2-1);
    oval_large = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % larger oval
    oval_small = x.^2 + 1.5*y.^2 <= (0.7*dim/2)^2; % smaller oval
    
    % Fat part of the image is in the outer circle
    % Water part the region within the outer circle, which constitutes fat
    F = oval_large - oval_small;
    W = oval_large;
    
    % Apply WF ratio and initialize R2* and dw0
    Wlib = W*(1-fatPercentage);
    Flib = F*fatPercentage;
    R2starlib = zeros(dim);
    dw0lib = zeros(dim)+dw0shift;
    fatmodel = get_fatmodel(num_fat_peaks);
    
    % Generate the dynamic images
    for dyn = 1:length(maxtemps)
        noise = noiselevel*randn(dim,dim,length(tes_heat));
        imgsdyn(:,:,:,dyn) = calcmeimgs(Wlib,Flib,tes_heat,b0,dw0lib+0*(dyn-1)*drift,...
            R2starlib,0*hotspot_data(:,:,dyn),prc,fatmodel)+noise;
    end
    
    % Provide a visualization
    subplot(2,2,1)
    imagesc(Wlib); colorbar; title('water'); axis('off')
    subplot(2,2,2)
    imagesc(Flib); colorbar; title('fat'); axis('off')
    subplot(2,2,3)
    imagesc(drift); colorbar; title('drift'); axis('off')
    subplot(2,2,4)
    imagesc(hotspot_data(:,:,end)*rad2degC); colorbar; title('temperature'); axis('off')
    
end

%% Masking

if strcmp(data,'measured')

    % Check if directory for body masks exists for current volunteer and create if not
    output_dir_BodyMasks = [datafile,'MRT Body\M',num2str(Pt,'%04u'),'\Masks'];
    if ~exist(output_dir_BodyMasks, 'dir')
        mkdir(output_dir_BodyMasks)
    end
    
    % Check if body mask for current dataset was already made
    file_body_mask = [output_dir_BodyMasks '\BodyMaskSalomir' sequence img_name_append '_slice' num2str(slice) '.mat'];
    if ~isfile(file_body_mask)
        % If not, use drawing tool to manually draw the mask and save it to
        % workspace as 'mask', then close the figure
        GUI = draw_mask(abs(scan(:,:,1,1))); uiwait(GUI.d.fig);
        body_mask = mask;
        save(file_body_mask,'body_mask');
    else
        % Else, we load the mask that was already generated previously
        load(file_body_mask);
    end

    % Check if directory for fat masks exists for current volunteer and create if not
    output_dir_FatMasks = [datafile,'MRT Body\M',num2str(Pt,'%04u'),'\Masks'];
    if ~exist(output_dir_FatMasks, 'dir')
        mkdir(output_dir_FatMasks)
    end
    
    % Check if mask for current dataset was already made
    file_fat_mask = [output_dir_FatMasks '\FatMaskSalomir' sequence img_name_append '_slice' num2str(slice) '.mat'];
    if ~isfile(file_fat_mask)
        % If not, use drawing tool to manually draw the mask and save it to
        % workspace as 'mask', then close the figure
        GUI = draw_mask(abs(scan(:,:,1,1))); uiwait(GUI.d.fig);
        fat_mask = mask;
        save(file_fat_mask,'fat_mask');
    else
        % Else, we load the mask that was already generated previously
        load(file_fat_mask);
    end

elseif strcmp(data,'artificial')
    
    % The drift field will be initialized using near harmonic 2D reconstruction, 
    % this algorithm requires a mask of fatty regions from which the harmonic
    % reconstruction can be initialized.
    % In addition to the fat mask, we also create masks for the body and for
    % regions with negligable signal
    % Three parameters are to be optimized for the fat, body and no-signal masks:
    %   1. th_fat    - threshold in fat/water image to assign voxels to fat mask 
    %   2. no_signal - threshold in respectively water and fat image to define low signal region
    %   3. SE        - structure element to erode thick fatty regions
    th_fat = 0.4;               % Higher value results in smaller 'fat' mask
    no_signal = [0.05, 0.05]; 	% Higher value results in larger 'no signal' mask (first corresponds to water and second value to fat)
    SE = strel('diamond', 0);   % Larger SE means thinner fat mask, but if too large we keep the original non-eroded mask
    body_factor = 0.9;          % Larger values result in larger body mask
    body_erode = 0;             % Integer value that determined how much the body mask is eroded
    
    % Compute the fat, body and no signal mask with current parameters
    [fat_mask,body_mask,signal_mask] = get_fatty_regions(Wlib(:,:,1), Flib(:,:,1), 1, [size(imgsdyn,1) size(imgsdyn,2)], th_fat, no_signal, SE, body_factor, body_erode);
    
    % Provide a visualization
    subplot(1,3,1); imagesc(fat_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('fat mask');
    subplot(1,3,2); imagesc(body_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('body mask');
    subplot(1,3,3); imagesc(signal_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('no signal mask');

end


%% Algorithm parameters

% Only algorithm parameters may be altered, which does not hold for scanner 
% and sequence since those parameters are constant

% Parameter structure: parameters to adjust the salomir fitting mask
algp.Nstd_ol = 4;               % Filter out local extreme values in the fat border
algp.neighbours = 4;            % Define local neighbourhood
algp.min_nb = 4;                % Eliminate voxels at the fat-tissue interface from the fat mask (value between 0 and 8)

% Visualize the mask that is used for the salomir fit
mask = remove_outer_border_pixels(fat_mask,abs(imgsdyn(:,:,1,1)),algp.min_nb);
imagesc(mask); axis off; axis image;


%% Perform iterative optimization to separate the temperature and drift fields
clear lib drift_recon T_recon;
disp('Optimization of multi-echo fat-suppressed MR thermometry has started..')

% Provide which dynamics are to be implemented in the optimization
show_ix = 1;   % number of dynamics that will be evaluated | size(imgsdyn,4)
init_ix = 100;                 % dynamic which will be evaluated first
step = 1;                    % step size in dynamics

tic
for ii = init_ix:step:init_ix+show_ix-1
    
    % Provide information on the duration
    disp(['dyn ' num2str(ii)]);

    % Evaluate a single dynamic
    dyn_img = imgsdyn(:,:,:,ii);

    % Define phase difference:
    % P = angle(Z) returns the phase angles, in radians, for each element of 
    % complex array Z. The angles lie between ±pi. So here we calculate the
    % differens in phase between scan 1 and scan ii.
    subtract = angle(dyn_img.*conj(imgsdyn(:,:,:,1)));

    % Get masks
    if strcmp(data,'artificial')
        [fat_mask,body_mask,signal_mask] = get_fatty_regions(Wlib(:,:,1), Flib(:,:,1), 1, [size(imgsdyn,1) size(imgsdyn,2)], th_fat, no_signal, SE, body_factor, body_erode);
    end

    % Perform near-harmonic 2D reconstruction using Salomir's algorithm
    % (initialization with the drift field estimation of the previous dynamic)
    if ii == init_ix
        drift_recon(:,:,ii) = Salomir_fit(subtract,zeros(size(dyn_img)),fat_mask,algp.Nstd_ol,algp.neighbours,algp.neighbours);
    else
        drift_recon(:,:,ii) = Salomir_fit(subtract,drift_recon(:,:,ii-step),fat_mask,algp.Nstd_ol,algp.neighbours,algp.neighbours);
    end

    % Calculate temperature (turn phase into deltaT)
    T_nocorr(:,:,ii) = subtract.*rad2degC/tes_heat;
    compphase = subtract - drift_recon(:,:,ii);
    T_recon(:,:,ii) = compphase.*rad2degC/tes_heat;

end
toc

figure; imagesc(T_recon(:,:,init_ix), [-1 1]); axis image; axis off;

%% Save data as .mat file and as NIfTI files

% Create output directory
output_dir = [datafile,'MRT Body\M',num2str(Pt,'%04u'),'\Optimization\',sequence];
if ~exist(output_dir, 'dir')    % Check if directory exists and create if not
   mkdir(output_dir)
end

% Save the optimized temperature and drift together with the body and fat mask as NIfTI files
niftiwrite(T_nocorr,[output_dir,'\T_nocorr.nii']);      	       % Reconstructed temperature
niftiwrite(T_recon,[output_dir,'\T_recon.nii']);                    % Reconstructed drift

% Save the relevant parameters as .mat file in the output directory
save([output_dir,'\optimized_data.mat'],'subtract','T_recon');


%%

x = 115;
y = 93;

figure;
subplot(1,2,1); imagesc(mean(T_recon(:,:,250:300),3), [-1 1]); hold on; axis off; axis image; colorbar; cb = colorbar; cb.FontSize = 16;
plot(x,y,'bx', 'MarkerSize', 10, 'linewidth', 3); hold off;
subplot(1,2,2); plot(squeeze(T_recon(y,x,:)), "Linewidth",5); ylim([-3 3]);


%% Compare artificially generated temperature and drift field maps to reconstructed temperature and drift field maps

dynamics_to_plot = [5,10,15,20,25];

for ii = 1:length(dynamics_to_plot)
    subplot(length(dynamics_to_plot),5,(ii-1)*5+1);
    imagesc(T_nocorr(:,:,dynamics_to_plot(ii)).*body_mask, [-1 1]); axis image; axis off; colorbar; title(['dynamic ' num2str(dynamics_to_plot(ii))]);
    subplot(length(dynamics_to_plot),5,(ii-1)*5+2);
    imagesc(T_recon(:,:,dynamics_to_plot(ii)).*body_mask, [-1 1]); axis image; axis off; colorbar;
    subplot(length(dynamics_to_plot),5,(ii-1)*5+3);
    imagesc(hotspot_data(:,:,dynamics_to_plot(ii)).*body_mask.*rad2degC, [-1 1]); axis image; axis off; colorbar;
    subplot(length(dynamics_to_plot),5,(ii-1)*5+4);
    imagesc(drift_recon(:,:,dynamics_to_plot(ii)).*body_mask); axis image; axis off; colorbar;
    subplot(length(dynamics_to_plot),5,(ii-1)*5+5);
    imagesc((dynamics_to_plot(ii)-1)*drift.*body_mask); axis image; axis off; colorbar;
end
  