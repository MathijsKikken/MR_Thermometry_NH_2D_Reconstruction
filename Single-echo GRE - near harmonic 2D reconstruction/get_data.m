function [scan_T,tes_T,GRE_info_T,scan_drift,tes_drift,GRE_info_drift,file_name] = get_data(img_name_append, info_name_append, sequence, datafile, Pt, slice, rel_echoes, rel_drift_echoes)

% Specify whether the baseline or the heated file will be loaded
if strcmp(sequence,'baseline')    % Baseline
    seq = 'MRTbase';
elseif strcmp(sequence,'heating') % Heating
    seq = 'MRTheat';
else
    disp('The parameter sequence was incorrectly specified')
end

% Find the muti-echo GRE MRT scan (and its info file) via .mat file
file_name = [datafile 'MRT Body\M',num2str(Pt,'%04u'),'\' seq '_mGRE_img' img_name_append '.mat'];
%file_name = ['D:\Thesis\Data\Processed_DIXON\M',num2str(Pt,'%04u'),'\' seq '_mGRE_img' img_name_append '.mat'];
path_file_scan = file_name;
path_file_info = [datafile 'MRT Body\M',num2str(Pt,'%04u'),'\' seq '_mGRE_info' info_name_append '.mat'];
%path_file_info = ['D:\Thesis\Data\Processed_DIXON\M',num2str(Pt,'%04u'),'\' seq '_mGRE_info' img_name_append '.mat'];

% Load data and provide a general name as it can be used in the remainder of the code
S_img = load(path_file_scan); N_img = fieldnames(S_img); eval(['GRE_scan = S_img.' N_img{1} ';']);
S_info = load(path_file_info); N_info = fieldnames(S_info); eval(['GRE_info = S_info.' N_info{1} ';']);

% The raw data has 7 dimensions, where the last 2 are GRE data
% Cast to complex and reshape the image such that there are 4 dimensions: [ nx | ny | nTE | dynamics]
if length(size(GRE_scan)) == 4
    scan = GRE_scan;
else
    scan = squeeze(GRE_scan(:,:,slice,:,:,1,1).*exp(1i*GRE_scan(:,:,slice,:,:,1,2)));
end

% Remove certain echoes: do this separately for the scan used for
% temperature measurements and the scan used for drift measurements
GRE_info_T = GRE_info;
GRE_info_drift = GRE_info;
scan_T = scan(:,:,rel_echoes,:);
GRE_info_T.imgdef.echo_time.uniq = GRE_info.imgdef.echo_time.uniq(rel_echoes);
scan_drift = scan(:,:,rel_drift_echoes,:);
GRE_info_drift.imgdef.echo_time.uniq = GRE_info.imgdef.echo_time.uniq(rel_drift_echoes);

% Acquire echo times
tes_T = GRE_info_T.imgdef.echo_time.uniq.*1e-3; 
tes_drift = GRE_info_drift.imgdef.echo_time.uniq.*1e-3;

end