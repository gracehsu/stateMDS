% catCarryingVoxel ver 0.5
% Using mask to find the voxels and extract voxel value from your data.
%
% Usage: [meanvalue, voxelvalue, voxmni, voxcor] = catCarryingVoxel(mask, data, mode, th);
%
%   Output:
%       meanvalue: the mean voxel value for each data. (mask x data matirx)
%       voxelvalue: give you every voxel value for each data. (mask x data)
%       voxmni: the MNI coordinates for each voxel in "voxelvalue".
%       voxcor: the index of 3D image matrix for each voxel.
%
%   Input:
%       mask: string or cell array for mask image. (*.nii or *.img; 
%             *_roi.mat is not working for now)
%       data: string or cell array for your data image. (*.nii or *.img)
%
%       mode: Here provide three mode to get voxel value. Please read 
%             carefully, and choose appropriate mode for your data. Input 
%             numeric number for correspoding mode.
%           mode 1: base on voxels space from "mask". (default)
%           mode 2: base on voxels space from "data".
%           mode 3: original "spm_mask" method from extract_voxel_value.
%                   The results from mode 3 should be equal to mode 2.
%           ** Rule of thumb: 
%           1. If you're interesting in mean voxel value: 
%              mask have smaller voxel size than data -> choose mode 1
%              data have smaller voxel size than mask -> choose mode 2
%           2. If you're interesting in multivoxel method which should 
%              retain every voxel value:
%              mask have smaller voxel size than data -> choose mode 2
%              data have smaller voxel size than mask -> choose mode 1
%
%       th: the threshold to your mask image. Default is 0.5 which is good
%           for binanry mask or probalistic mask.
%
% Example: 
%   [meanvalue, voxelvalue, voxmni, voxcor] = catCarryingVoxel(...
%       spm_select(Inf, 'image', 'Select mask file'), spm_select(Inf, ...
%       'image', 'Select data file'), 1, 0.5);
%
% Dependicies: SPM (tested on SPM12 revision 6906)
%
% Version log:
% Yu-Shiang Su Dec 12 2016: ver 0.1 created this script.
% Yu-Shiang Su Feb 21 2017: ver 0.2 minor modification for public version.
% Yu-Shiang Su Feb 27 2017: ver 0.3 fixed minor bug in mode 2.
% Yu-Shiang Su Feb 28 2017: ver 0.4 add feature of outputing coordinates 
%                                   in mni or voxel space.
% Yu-Shiang Su Mar 01 2017: ver 0.5 minor bug fixed

function [meanvalue, voxelvalue, voxmni, voxcor] = catCarryingVoxel(mask, data, mode, th)

if nargin < 3 || isempty(mode)
    mode = 1;
elseif nargin < 4 || isempty(th)
    th = 0.5;
end

%% Get image header
% if the data is cell array, transform to string array. It will be more
% easy to get the index from each image volume.
if iscell(data)
    data_V = spm_vol(char(data));
else
    data_V = spm_vol(data);
end

if iscell(mask)
    mask_V = spm_vol(char(mask));
else
    mask_V = spm_vol(mask);
end
mask_n = length(mask_V);

byeachdata = 0; % default
if any(ismember(cat(1,data_V.dim), data_V(1).dim, 'rows') == 0)
    fprintf('%s\n', 'Your data are not in the same space (V.dim). ROI space will calculate by each data.');
    byeachdata = 1;
end
if any(ismember(reshape(cat(1,data_V.mat)', 16, [])', reshape(data_V(1).mat', 16, [])', 'rows') == 0)
    fprintf('%s\n', 'Your data have different transformation matrix (V.mat), ROI space will calculate by each data.');
    byeachdata = 1;
end
if byeachdata == 1;
    fprintf('%s\n', 'If you have big data here. This will be dramaticaaly slower than regular processing.');
    voxmni = cell(mask_n, length(data_V));
    voxcor = cell(mask_n, length(data_V));
else
    voxmni = cell(mask_n, 1);
    voxcor = cell(mask_n, 1);
end

meanvalue = nan(mask_n, length(data_V));
voxelvalue = cell(mask_n, length(data_V));

switch mode
    %% mode 1: extract value base on voxels space from "mask"
    case 1
        for mask_counter = 1: length(mask_V)
            [mask_Y, mask_xyz] = spm_read_vols(mask_V(mask_counter));
            mask_mnixyz = mask_xyz(:, mask_Y> th);
            % add one for later matirx multiplication
            mask_mnixyz = [mask_mnixyz; ones(1,size(mask_mnixyz,2))];
            
            if byeachdata == 1
                for data_counter = 1:length(data_V)
                    mask_datavoxxyz = mask_mnixyz' * inv(data_V(data_counter).mat)';
                    rawdata = spm_get_data(data_V(data_counter), mask_datavoxxyz');
                    invalidvox = (rawdata == 0) | isnan(rawdata);
                    meanvalue(mask_counter, data_counter) = mean(rawdata(~invalidvox));
                    voxelvalue{mask_counter, data_counter} = rawdata(~invalidvox);
                    rawdata_mnixyz = (data_V(data_counter).mat*[mask_datavoxxyz(~invalidvox,1) mask_datavoxxyz(~invalidvox,2) mask_datavoxxyz(~invalidvox,3) ones(sum(~invalidvox),1)]')';
                    voxmni{mask_counter, data_counter} = rawdata_mnixyz(:,1:3);
                    voxcor{mask_counter, data_counter} = mask_datavoxxyz(~invalidvox, 1:3);
                    fprintf('%s%d%s%d%s%d%s%d%s%d%s\n', 'Mask ', mask_counter, ...
                            ', Data ', data_counter, ': Found ', ...
                            length(mask_mnixyz), ...
                            ' voxels from mask space, get ', ...
                            sum(~invalidvox), ' valid voxels and ', ...
                            sum(invalidvox), ...
                            ' invalid voxels from data space.');
                end
            else
                mask_datavoxxyz = mask_mnixyz' * inv(data_V(1).mat)';
                rawdata = spm_get_data(data_V, mask_datavoxxyz');
                invalidvox = (sum(rawdata, 1) == 0) | isnan(sum(rawdata,1));
                meanvalue(mask_counter, :) = mean(rawdata(:,~invalidvox),2);
                voxelvalue(mask_counter, :) = mat2cell(rawdata(:,~invalidvox), ones(1,size(rawdata,1)), sum(~invalidvox));
                rawdata_mnixyz = (data_V(1).mat*[mask_datavoxxyz(~invalidvox,1) mask_datavoxxyz(~invalidvox,2) mask_datavoxxyz(~invalidvox,3) ones(sum(~invalidvox),1)]')';
                voxmni{mask_counter, 1} = rawdata_mnixyz(:,1:3);
                voxcor{mask_counter, 1} = mask_datavoxxyz(~invalidvox, 1:3);
                fprintf('%s%d%s%d%s%d%s%d%s\n', 'Mask ', mask_counter, ...
                            ': Found ', ...
                            length(mask_mnixyz), ...
                            ' voxels from mask space, get ', ...
                            length(unique(rawdata(1,:))), ' valid voxels and ', ...
                            length(find(isnan(rawdata(1,:)))), ...
                            ' invalid voxels from data space.');
            end
        end
        
    %% mode 2: extract value base on voxels space from "data"
    case 2
        for mask_counter = 1:length(mask_V)
            [mask_Y, ~] = spm_read_vols(mask_V(mask_counter));
            mask_idx = find(mask_Y > th);
            mask_size = mask_V(mask_counter).dim;
            [mask_idxx, mask_idxy, mask_idxz] = ind2sub(mask_size(1:3), mask_idx);
            if byeachdata == 1
                for data_counter = 1:length(data_V)
                    [~, data_xyz] = spm_read_vols(data_V(data_counter));
                    data_mnixyz = [data_xyz; ones(1, size(data_xyz,2))];
                    data_maskvoxxyz = round(data_mnixyz' * inv(mask_V(mask_counter).mat)');
                    data_idx = find(ismember(data_maskvoxxyz(:,1:3), [mask_idxx, mask_idxy, mask_idxz], 'rows'));
                    data_size = data_V(data_counter).dim;
                    [data_voxxyzx, data_voxxyzy, data_voxxyzz]= ind2sub(data_size(1:3), data_idx);
                    data_voxxyz = [data_voxxyzx, data_voxxyzy, data_voxxyzz];
                    rawdata = spm_get_data(data_V(data_counter), data_voxxyz');
                    invalidvox = (rawdata == 0) | isnan(rawdata);
                    meanvalue(mask_counter, data_counter) = mean(rawdata(~invalidvox));
                    voxelvalue{mask_counter, data_counter} = rawdata(~invalidvox);
                    rawdata_mnixyz = (data_V(data_counter).mat*[data_voxxyz(~invalidvox,1) data_voxxyz(~invalidvox,2) data_voxxyz(~invalidvox,3) ones(sum(~invalidvox),1)]')';
                    voxmni{mask_counter, data_counter} = rawdata_mnixyz(:,1:3);
                    voxcor{mask_counter, data_counter} = data_voxxyz(~invalidvox,1);
                    fprintf('%s%d%s%d%s%d%s%d%s%d%s\n', 'Mask ', mask_counter, ...
                            ', Data ', data_counter, ': Found ', ...
                            length(mask_idx), ...
                            ' voxels from mask space, get ', ...
                            sum(~invalidvox), ' valid voxels and ', ...
                            sum(invalidvox), ...
                            ' invalid voxels from data space.');
                end
            else
                [~, data_xyz] = spm_read_vols(data_V(1));
                data_mnixyz = [data_xyz; ones(1, size(data_xyz,2))];
                data_maskvoxxyz = round(data_mnixyz' * inv(mask_V(mask_counter).mat)');
                data_idx = find(ismember(data_maskvoxxyz(:,1:3), [mask_idxx, mask_idxy, mask_idxz], 'rows'));
                data_size = data_V(1).dim;
                [data_voxxyzx, data_voxxyzy, data_voxxyzz]= ind2sub(data_size(1:3), data_idx);
                data_voxxyz = [data_voxxyzx, data_voxxyzy, data_voxxyzz];
                rawdata = spm_get_data(data_V, data_voxxyz');
                invalidvox = (sum(rawdata, 1) == 0) | isnan(sum(rawdata,1));
                meanvalue(mask_counter, :) = mean(rawdata(:,~invalidvox),2);
                voxelvalue(mask_counter, :) = mat2cell(rawdata(:,~invalidvox), ones(1,size(rawdata,1)),sum(~invalidvox));
                rawdata_mnixyz = (data_V(1).mat*[data_voxxyz(~invalidvox,1) data_voxxyz(~invalidvox,2) data_voxxyz(~invalidvox,3) ones(sum(~invalidvox),1)]')';
                voxmni{mask_counter,1} = rawdata_mnixyz(:,1:3);
                voxcor{mask_counter,1} = data_voxxyz(~invalidvox,:);
                fprintf('%s%d%s%d%s%d%s%d%s\n', 'Mask ', mask_counter, ...
                            ': Found ', ...
                            length(mask_idx), ...
                            ' voxels from mask space, get ', ...
                            sum(~invalidvox), ' valid voxels and ', ...
                            sum(invalidvox), ...
                            ' invalid voxels from data space.');
            end
        end
        
    %% mode 3: extract value using spm_mak (extract_voxel_values.m)
    case 3
        mask_idx = cell(mask_n, 1);
        [data_pth, data_nam, data_ext, data_num] = spm_fileparts(data_V(1).fname);
        if strcmp(data_ext, 'nii')
            if strcmp(data_num, '')
                data_num = num2str(data_V(1).n(1));
            end
        end
        for mask_counter = 1:mask_n
            % Create temporary mask
            spm_mask([mask_V(mask_counter).fname, ',', num2str(mask_V(mask_counter).n(1))], ...
                [data_pth, filesep, data_nam, data_ext, data_num], th)
            
            % get voxel index from temporary mask (select the voxel > 0 and non-nan)
            temp_mask_V = spm_vol([data_pth, filesep, 'm', data_nam, data_ext, data_num]);
            temp_mask_Vol = spm_read_vols(temp_mask_V);
            mask_idx{mask_counter} = find(temp_mask_Vol~=0 & ~isnan(temp_mask_Vol));
            fprintf('%s%d%s%s\n', 'Found ', length(mask_idx{mask_counter}), ...
                ' voxels in ROI', num2str(mask_counter));
            
            % delete temporary mask
            delete([data_pth, filesep, 'm', data_nam, data_ext, data_num])
            if strcmp(data_ext, '.img')
                delete([data_pth, filesep, 'm', data_nam, '.hdr'])
            end
        end
        fprintf('%s', 'Getting data ...')
        data_Vol = spm_read_vols(data_V);
        data_Vol2D = reshape(data_Vol, [], length(data_V));
        
        for mask_counter = 1:mask_n
            rawdata = data_Vol2D(mask_idx{mask_counter}, :)';
            meanvalue(mask_counter, :) = mean(rawdata, 2);
            voxelvalue(mask_counter, :) = mat2cell(rawdata(:,:), ones(1,size(rawdata,1)), size(rawdata,2));
        end
        
        fprintf('%s\n', 'Done!')
end