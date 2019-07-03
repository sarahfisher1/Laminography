% LAMINO_INTERFACE
% -------------------------------------------------------------------------
% Function that you call to run the interface. From the interface you can
% set the reconstruction parameters, corrections, algorithm, sinogram and
% start the reconstruction. This function contains the definitions of all
% the buttons, checkboxes and fields on the interface window and defines
% what happens when there is an interaction on the interface. 
%
% TESTED ENVIRONMENT:
%   Windows 10 with GPU, MATLAB version 2017a
%
% HARDWARE REQUIREMENTS:
%   This code runs on a single GPU
%
% SOFTWARE REQUIREMENTS:
%   - ASTRA version 1.8 (available from Github)
%   - CUDA driver version 8 
%
% INPUT:
%   no inputs
% OUTPUT:
%   no output
% 
% The function will save the reconstructed volume as a .mat file (MATLAB
% file) and a .vol file (file that can be read into visualisation programs
% such as Avizo). The projection geometry is saved as a .mat file. 
%
% -------------------------------------------------------------------------
% Copyright (c) 2019 S L Fisher D J Holmes J S Jørgensen P Gajjar J Behnsen W R B Lionheart P J Withers
% University of Manchester
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Contact: sarah.fisher-1@outlook.com
% -------------------------------------------------------------------------

function lamino_interface() 

fig = figure('name', 'Laminography reconstruction', 'numbertitle', 'off');
set(fig, 'Position', [500, 500, 525, 400]);

working_dir = pwd;

% Print the copyright notice 
fprintf('Laminography reconstruction software\n')
fprintf('---------------------------------\n')
fprintf('MATLAB code for laminography reconstruction using SIRT and CGLS\n');
fprintf('Copyright (c) 2019 S Fisher and D Holmes\n');
fprintf('------------------------------------------\n');
fprintf(['This program is free software: you can redistribute it and/or modify\n',...
    'it under the terms of the GNU General Public License as published by\n', ...
     'the Free Software Foundation, version 3 of the License.\n\n']);
fprintf(['This program is distributed in the hope that it will be useful\n', ...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of \n', ...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n', ...
    'GNU General Public License for more details.\n\n']);
fprintf(['You should have received a copy of the GNU General Public License \n',...
    'along with this program. If not, see <http://www.gnu.org/licenses/>.\n\n'])

% Print the dependancies
fprintf('This software requires the ASTRA Toolbox 1.8 and CUDA driver version 8.\n');
fprintf('The ASTRA toolbox folder must be in the current working directory.\n');
fprintf('This code runs on a single GPU\n');
fprintf('This software has been tested on Windows 10 and MATLAB version 2017a\n');

% Set up the workspace for a laminography reconstruction 
status = setup_recon_laminography();

% Check the driver version. The laminography tools only works using version
% 8.
driver = status.DriverVersion; 
fprintf('Using CUDA version %1.1f\n', driver);
if (driver < 8)
   error('The laminography reconstruction code can only be run using CUDA v8. Please update the driver\n');
end

% Pre allocate data structures 
parameters.Alpha = 0; 
parameters.Iterations = 0; 
parameters.COR = 0;
parameters.ROI = 0;
parameters.BinningFactor = 0; 
parameters.Algorithm = ' ';
parameters.ChopDataFactor = 1;

dimensions.Nx = 0;
dimensions.Ny = 0;
dimensions.Nz = 0;
dimensions.MountHeight = 0;
dimensions.CORPositionCorrection = 0;

file_paths.DataPath = ' ';
file_paths.AnglesFile = ' ';
file_paths.SinogramFile = ' '; 
file_paths.ScanDimensionsFile = ' ';

ctrl_variables.SaveReconstruction = false;
ctrl_variables.SaveFrequency = false;
ctrl_variables.MirrorImages = false;
ctrl_variables.LoadSinogram = false;
ctrl_variables.ChopData = false;

%% Title textbox 
%% ========================================================================

% Title textbox 
title_textbox = uicontrol('style', 'text', 'FontSize', 16, 'Position', [80 365 400 30]);
set(title_textbox, 'String', 'Laminography reconstruction code');

% 
%% Tilt angle text box 
%% ========================================================================
tilt_angle = 30; 
tilt_angle_min = 0;
tilt_angle_max = 90;

tilt_angle_field = uicontrol('Style', 'edit', 'String', tilt_angle, ...
    'Position', [100 330 70 20], 'FontSize', 11, 'Callback', @setangle);

tilt_angle_text = uicontrol('Style', 'text', 'String', 'Tilt angle', ...
    'Position', [10 330 80 20], 'FontSize', 11);

tilt_angle_unit_text = uicontrol('Style', 'text', 'String', 'deg', ...
    'Position', [180 330 30 20], 'FontSize', 11);

% Set the tilt angle whenever the user inputs a new number
function setangle(source, event)
    tilt_angle_prev = tilt_angle;
    tilt_angle_str = get(source, 'String');
    tilt_angle = str2double(tilt_angle_str);
    if(isnan(tilt_angle))
        tilt_angle = tilt_angle_prev; 
        fprintf('Tilt angle must be a number\n');
    end
    if ((tilt_angle > tilt_angle_max) || (tilt_angle < tilt_angle_min))
        fprintf('The tilt angle must be between %1.0f deg and %1.0f deg\n', tilt_angle_max, tilt_angle_min);
        tilt_angle = tilt_angle_prev; 
    end
    fprintf('The tilt angle is %1.1f�\n', tilt_angle);
end

%% Angle file name selection box
%% ========================================================================
angles_file_name = working_dir;

angles_file_path = uicontrol('Style', 'pushbutton', 'String', '...', ...
    'Position', [95 300 120 20], 'FontSize', 11, 'Callback', @setanglespath, ...
    'Visible', 'on');

angles_file_path_text = uicontrol('Style', 'text', 'String', 'Angles file ', ...
    'Position', [10 300 80 20], 'FontSize', 11, ...
    'Visible', 'on');


function setanglespath(source, event)

    angles_file_path_prev = angles_file_name;

    [file_name, folder_name] = uigetfile('*.ang', 'Select angles file');
    disp(file_name);

    angles_file_name = file_name;

    fprintf('The name of the file containing the angle data is %s\n', angles_file_name);

    set(angles_file_path, 'String', angles_file_name);

    cd(working_dir);
end

%% Reconstruction volumes 
%% ========================================================================

vol_text = uicontrol('Style', 'text', 'String', 'Reconstruction volume dimensions', ...
    'Position', [10 230 110 60], 'FontSize', 11);

%% Nx text box 
%% ------------------------------------------------------------------------
% Default & limit
nx = 500; 
nx_min = 10;

nx_field = uicontrol('Style', 'edit', 'String', nx, ...
    'Position', [160 270 50 19], 'FontSize', 11, 'Callback', @setnx);

nx_text = uicontrol('Style', 'text', 'String', 'x', ...
    'Position', [125 270 30 18], 'FontSize', 11);

function setnx(source, event)
    nx_prev = nx;
    ny_prev = ny;
    nx_str = get(source, 'String');
    nx = round(str2double(nx_str));
    ny = nx;
    if(isnan(nx))
        nx = nx_prev;
        ny = ny_prev;
        fprintf('The X dimension in the reconstruction volume must be a number\n');
    end
    if (nx <= nx_min)
        fprintf('The X dimension in the reconstruction volume must be greater than %1.0f\n', nx_min);
        nx = nx_prev;
        ny = ny_prev;
    end
    set(nx_field, 'String', nx);
    set(ny_field, 'String', ny);
end


%% Ny text box 
%% ------------------------------------------------------------------------
ny = 500; 
ny_min = 10;

ny_field = uicontrol('Style', 'edit', 'String', ny, ...
    'Position', [160 250 50 19], 'FontSize', 11, 'Callback', @setny);

ny_text = uicontrol('Style', 'text', 'String', 'y', ...
    'Position', [125 250 30 18], 'FontSize', 11);

function setny(source, event)
    ny_prev = ny;
    nx_prev = nx;
    ny_str = get(source, 'String');
    ny = round(str2double(ny_str));
    nx = ny;
    if(isnan(ny))
        ny = ny_prev;
        nx = nx_prev;
        fprintf('The Y dimension in the reconstruction volume must be a number\n');
    end
    if (ny <= ny_min)
        fprintf('The Y dimension in the reconstruction volume must be greater than %1.0f\n', ny_min);
        ny = ny_prev;
        nx = nx_prev;
    end
    set(ny_field, 'String', ny);
    set(nx_field, 'String', nx);
end

%% Nz text box 
%% ------------------------------------------------------------------------
% Defaults & limit
nz = 500;
nz_min = 10;

nz_field = uicontrol('Style', 'edit', 'String', nz, ...
    'Position', [160 230 50 19], 'FontSize', 11, 'Callback', @setnz);

nz_text = uicontrol('Style', 'text', 'String', 'z', ...
    'Position', [125 230 30 18], 'FontSize', 11);

function setnz(source, event)
    nz_prev = nz;
    nz_str = get(source, 'String');
    nz = round(str2double(nz_str));
    if(isnan(nz))
        nz = nz_prev; 
        fprintf('The Z dimension in the reconstruction volume must be a number\n');
    end
    if (nz <= nz_min)
        fprintf('The Z dimension in the reconstruction volume must be greater than %1.0f\n', nz_min);
        nz = nz_prev; 
    end
    %fprintf('The Z dimension in the reconstruction volume is %1.0f\n', nz);
    set(nz_field, 'String', nz);
end

%% Tick boxes for control variables 
%% ========================================================================

%% Center of rotation correction text box 
%% ------------------------------------------------------------------------
padsize_cor = 0; % Default

cor_field = uicontrol('Style', 'edit', 'String', padsize_cor, ...
    'Position', [100 150 70 20], 'FontSize', 11, 'Callback', @setcor);

cor_text = uicontrol('Style', 'text', 'String', 'COR', ...
    'Position', [10 150 80 20], 'FontSize', 11);

cor_unit_text = uicontrol('Style', 'text', 'String', 'pixels', ...
    'Position', [180 150 40 20], 'FontSize', 11);

% Set the tilt angle whenever the user inputs a new number
function setcor(source, event)
    cor_prev = padsize_cor;
    cor_str = get(source, 'String');
    padsize_cor = round(str2double(cor_str));
    if(isnan(padsize_cor))
        padsize_cor = cor_prev; 
        fprintf('Center of rotation corrections must be a number\n');
    end
    set(cor_field, 'String', padsize_cor);
    fprintf('The center of rotation correction is %1.0f pixels\n', padsize_cor);
end

%% Region of interest correction text box
%% ------------------------------------------------------------------------
padsize_roi = 0; %Default

roi_field = uicontrol('Style', 'edit', 'String', padsize_roi, ...
    'Position', [100 120 70 20], 'FontSize', 11, 'Callback', @setroi);

roi_text = uicontrol('Style', 'text', 'String', 'ROI', ...
    'Position', [10 120 80 20], 'FontSize', 11);

roi_unit_text = uicontrol('Style', 'text', 'String', 'pixels', ...
    'Position', [180 120 40 20], 'FontSize', 11);

% Set the tilt angle whenever the user inputs a new number
function setroi(source, event)
    roi_prev = padsize_roi;
    roi_str = get(source, 'String');
    padsize_roi = round(str2double(roi_str));
    if(isnan(padsize_roi))
        padsize_roi = roi_prev; 
        fprintf('Region of interest corrections must be a number\n');
    end
    if (padsize_roi < 0)
        padsize_roi = roi_prev; 
        fprintf('Region of interest corrections must be positive\n');
    end
    set(roi_field, 'String', padsize_roi);
    fprintf('The region of interest correction is %1.0f pixels\n', padsize_roi);
end


%% Iteration number text box 
%% ------------------------------------------------------------------------

% Defaults
iterations= 50; 
iterations_min = 1;

iterations_field = uicontrol('Style', 'edit', 'String', iterations, ...
    'Position', [100 175 70 20], 'FontSize', 11, 'Callback', @setiterations, ...
    'Visible', 'on');

iterations_text = uicontrol('Style', 'text', 'String', 'Iterations', ...
    'Position', [10 175 80 20], 'FontSize', 11, ...
    'Visible', 'on');

% Set the tilt angle whenever the user inputs a new number
function setiterations(source, event)
    iterations_prev = iterations;
    iterations_str = get(source, 'String');
    iterations = round(str2double(iterations_str));
    if(isnan(iterations))
        iterations = iterations_prev; 
        fprintf('The number of iterations must be a number\n');
    end
    if (iterations < iterations_min)
        fprintf('The number of iterations must be positive\n');
        iterations = iterations_prev; 
    end
    fprintf('The number of iterations is %1.0f \n', iterations);
    set(iterations_field, 'String', iterations);
end

%% Method drop down list 
%% ------------------------------------------------------------------------

method_str = 'CGLS';

method_button = uicontrol('Style', 'popupmenu', 'String', ...
    {'CGLS', 'SIRT'}, 'Position', [100 200 80 20], 'Callback', @setmethod);

method_text = uicontrol('Style', 'text', 'String', 'Method', ...
    'Position', [10 200 80 20], 'FontSize', 11);
% 
function setmethod(source, event)

    method_options = get(source, 'String');
    method_index = get(source, 'Value');
    method_str = method_options{method_index};
    
    disp(method_str);
    
    if (strcmpi(method_str, 'CGLS') || strcmpi(method_str, 'SIRT'))
        set(iterations_field, 'Visible', 'on');
        set(iterations_text, 'Visible', 'on');
    else
        set(iterations_field, 'Visible', 'off');
        set(iterations_text, 'Visible', 'off');
    end

end

%% Generate new sinogram tick box 
%% ------------------------------------------------------------------------

sinogram_button = uicontrol('Style', 'pushbutton', 'String', 'Generate sinogram', ...
    'Position', [250 150 250 40], 'FontSize', 16, 'Callback', @generatesinogram);

 % Control function defined further down

%% Sinogram data path selection
%% ------------------------------------------------------------------------

sinogram_data_path = working_dir;

sino_datapath_field = uicontrol('Style', 'pushbutton', 'String', '...', ...
    'Position', [350 280 150 20], 'FontSize', 11, 'Callback', @setdatapath, ...
    'Visible', 'on');

sino_datapath_text = uicontrol('Style', 'text', 'String', 'Data location', ...
    'Position', [250 280 100 20], 'FontSize', 11, ...
    'Visible', 'on');

    function setdatapath(source, event)
        
        sinogram_data_path_prev = sinogram_data_path;
        
        folder_name = uigetdir('I:\Projects\MPhys20162017\lami', 'Select folder where projection data is located');
        disp(folder_name);
        
        sinogram_data_path = folder_name;
  
        fprintf('The name of the folder containing the projection data is %s\n', sinogram_data_path);
        
        set(sino_datapath_field, 'String', sinogram_data_path);
        
        cd(working_dir);
    end

%% Shading correction data path selection
%% ------------------------------------------------------------------------

%% DARK

dark_data_path = working_dir;


shading_corr_datapath_text = uicontrol('Style', 'text', 'String', 'Shading correction', ...
    'Position', [250 255 150 20], 'FontSize', 11, ...
    'Visible', 'on');

dark_datapath_text = uicontrol('Style', 'text', 'String', 'Dark', ...
    'Position', [250 230 50 20], 'FontSize', 11, ...
    'Visible', 'on');

dark_datapath_field = uicontrol('Style', 'pushbutton', 'String', '...', ...
    'Position', [300 230 50 20], 'FontSize', 11, 'Callback', @setdarkdatapath, ...
    'Visible', 'on');

%% LIGHT 

light_data_path = working_dir;

light_datapath_text = uicontrol('Style', 'text', 'String', 'Light', ...
    'Position', [355 230 50 20], 'FontSize', 11, ...
    'Visible', 'on');

light_datapath_field = uicontrol('Style', 'pushbutton', 'String', '...', ...
    'Position', [405 230 50 20], 'FontSize', 11, 'Callback', @setlightdatapath, ...
    'Visible', 'on');



function setdarkdatapath(source, event)

    [file_name, folder_name] = uigetfile('*.tif', 'Select dark file');
    dark_data_path = strcat(folder_name, file_name);
    fprintf('The dark file is %s\n', dark_data_path);

    set(dark_datapath_field, 'String', dark_data_path);
    cd(working_dir);

end

function setlightdatapath(source, event)

    [file_name, folder_name] = uigetfile('*.tif', 'Select light file');
    light_data_path = strcat(folder_name, file_name);
    fprintf('The light file is %s\n', light_data_path);

    set(light_datapath_field, 'String', light_data_path);
    cd(working_dir);

end

%% Sinogram binning list selection box
%% ------------------------------------------------------------------------
%Defaults
binning = 1;
binning_on = false;

bin_button = uicontrol('Style', 'popupmenu', 'String', ...
    {'1', '2', '4'}, 'Position', [405 200 40 20], 'Callback', @setbinning);

bin_text = uicontrol('Style', 'text', 'String', 'Binning factor', ...
    'Position', [250 200 150 20], 'FontSize', 11);
% 

    function setbinning(source, event)
        
        binning_options = get(source, 'String');
        binning_index = get(source, 'Value');
        binning = str2double(binning_options{binning_index});
        
        if (binning ~=1)
            binning_on = true; 
            fprintf('Binning of projection images is enabled at a binning factor of %1.0f\n', binning);
        else
            binning_on = false;
            fprintf('Binning of projection images is disabled\n');
        end
        
    end

%% Load sinogram tick box 
%% ------------------------------------------------------------------------
load_sino = false; %Default
sino_file_name = 'sinogram.mat';

sino_from_file_box = uicontrol('Style', 'checkbox', 'String', 'Load sinogram from file', 'Value', 0,...
    'Position', [250 305 200 20], 'FontSize', 11, 'Callback', @sinoload);

sino_file_path = uicontrol('Style', 'edit', 'String', sino_file_name, ...
                'Position', [305 285 205 20], 'FontSize', 11, 'Callback', @sinofilepath, ...
                'Visible', 'off');

sino_file_path_text = uicontrol('Style', 'text', 'String', 'File name:', ...
                'Position', [250 285 80 20], 'FontSize', 11, 'Visible', 'off');
            

%% Generate sinogram tick box 
%% ------------------------------------------------------------------------
sino_generate_box = uicontrol('Style', 'checkbox', 'String', 'Generate sinogram', 'Value', 1,...
        'Position', [250 330 200 20], 'FontSize', 11, 'Callback', @sinogenerate);

    
    function sinofilepath(source, event)
        
        sino_file_name_prev = sino_file_name;
        sino_file_name = get(source, 'String');
        set(sino_file_path, 'String', sino_file_name);
        
    end

%% Generate sinogram tick box control functions 
%% ------------------------------------------------------------------------
generate_sino = true; % Default

    function sinogenerate(source, event)
        
        generate_sino_ans = get(source, 'Value');
        
        if (generate_sino_ans == 1)
            generate_sino = true; 
            load_sino = false;
            fprintf('Generate sinogram:\n');
            disp(generate_sino);
            fprintf('Load sinogram :\n');
            disp(load_sino);
            set(sino_from_file_box, 'Value', 0);
            set(sinogram_button, 'Visible', 'on');
            set(sino_datapath_field, 'Visible', 'on');
            set(sino_datapath_text, 'Visible', 'on');
            set(shading_corr_datapath_text, 'Visible', 'on');
            set(dark_datapath_text, 'Visible', 'on');
            set(dark_datapath_field, 'Visible', 'on');
            set(light_datapath_text, 'Visible', 'on');
            set(light_datapath_field, 'Visible', 'on');
            set(sino_file_path, 'Visible', 'off');
            set(sino_file_path_text, 'Visible', 'off');
        else
            generate_sino = false; 
            load_sino = true;
            fprintf('Generate sinogram:\n');
            disp(generate_sino);
            fprintf('Load sinogram :\n');
            disp(load_sino);
            set(sino_from_file_box, 'Value', 1);
            set(sinogram_button, 'Visible', 'off');
            set(sino_datapath_field, 'Visible', 'off');
            set(sino_datapath_text, 'Visible', 'off');
            set(shading_corr_datapath_text, 'Visible', 'off');
            set(dark_datapath_text, 'Visible', 'off');
            set(dark_datapath_field, 'Visible', 'off');
            set(light_datapath_text, 'Visible', 'off');
            set(light_datapath_field, 'Visible', 'off');
            set(sino_file_path, 'Visible', 'on');
            set(sino_file_path_text, 'Visible', 'on');
        end
    end

    function sinoload(source, event)
        load_sino_ans = get(source, 'Value');
        if (load_sino_ans == 1)
            load_sino = true;
            generate_sino = false; 
            set(sino_generate_box, 'Value', 0);
            set(sino_file_path, 'Visible', 'on');
            set(sino_file_path_text, 'Visible', 'on');
            set(sino_from_file_box, 'Value', 1);
            set(sinogram_button, 'Visible', 'off');
            set(sino_datapath_field, 'Visible', 'off');
            set(sino_datapath_text, 'Visible', 'off');
            set(shading_corr_datapath_text, 'Visible', 'off');
            set(dark_datapath_text, 'Visible', 'off');
            set(dark_datapath_field, 'Visible', 'off');
            set(light_datapath_text, 'Visible', 'off');
            set(light_datapath_field, 'Visible', 'off');
        else
            load_sino = false;
            generate_sino = true;
            set(sino_generate_box, 'Value', 1);
            set(sino_file_path, 'Visible', 'off');
            set(sino_file_path_text, 'Visible', 'off');
            set(sinogram_button, 'Visible', 'on');
            set(sino_datapath_field, 'Visible', 'on');
            set(sino_datapath_text, 'Visible', 'on');
            set(shading_corr_datapath_text, 'Visible', 'on');
            set(dark_datapath_text, 'Visible', 'on');
            set(dark_datapath_field, 'Visible', 'on');
            set(light_datapath_text, 'Visible', 'on');
            set(light_datapath_field, 'Visible', 'on');
        end
    end

function generatesinogram (source, event)
    fprintf('Generating sinogram!\n');
    % Hide sinogram button so no more sinograms can be made
    set(sinogram_button, 'Visible', 'off');
    sinogram = get_sinogram_3D(sinogram_data_path, dark_data_path, light_data_path, binning_on, binning);
    save('sinogram', 'sinogram', '-v7.3');
end

%% Scan dimensions file selection 
%% ------------------------------------------------------------------------

scan_dimensions_file_name = working_dir; % Default --> Meaningless as not pointing at a file

scan_dimensions_file_path = uicontrol('Style', 'pushbutton', 'String', '...', ...
    'Position', [95 90 120 20], 'FontSize', 11, 'Callback', @setdimensionspath, ...
    'Visible', 'on');

scan_dimensions_file_path_text = uicontrol('Style', 'text', 'String', 'Dimensions file ', ...
    'Position', [10 90 80 20], 'FontSize', 11, ...
    'Visible', 'on');

function setdimensionspath(source, event)

    scan_dimensions_file_path_prev = scan_dimensions_file_name;

    [file_name, folder_name] = uigetfile('*.xtekct', 'Select dimensions file');
    disp(file_name);

    scan_dimensions_file_name = file_name;

    fprintf('The name of the file containing the dimension data is %s\n', scan_dimensions_file_name);

    set(scan_dimensions_file_path, 'String', scan_dimensions_file_name);

    cd(working_dir);
end

%% Mount height text box 
%% ------------------------------------------------------------------------
mount_height = 180; % Default (in mm)

mount_height_field = uicontrol('Style', 'edit', 'String', mount_height, ...
                'Position', [132 60 50 20], 'FontSize', 11, 'Callback', @mountheightchanged, ...
                'Visible', 'on');
mount_height_text = uicontrol('Style', 'text', 'String', 'Mount height:', ...
                'Position', [10 60 120 20], 'FontSize', 11, 'Visible', 'on');
            
mount_height_units_text = uicontrol('Style', 'text', 'String', 'mm', ...
                'Position', [185 60 30 20], 'FontSize', 11, 'Visible', 'on');
            
function mountheightchanged(source, event)
    
    mount_height_prev = mount_height; % Save the previous source to center distance
    mount_height = str2double(get(source, 'String'));
    
    if(isnan(mount_height))
        mount_height = mount_height_prev; 
        fprintf('The mount height must be a number\n');
    end
    
    if (mount_height <=  0)
        fprintf('The mount height must be positive and non-zero.\n');
        mount_height = mount_height_prev;
    end
    
    fprintf('The mount height is %1.2f mm\n', mount_height);
    set(mount_height_field, 'String', mount_height);
    
end

%% COR position text box 
%% ------------------------------------------------------------------------
cor_pos_corr = 0; % Default

cor_pos_corr_field = uicontrol('Style', 'edit', 'String', cor_pos_corr, ...
                'Position', [132 20 50 40], 'FontSize', 11, 'Callback', @corposcorrchanged, ...
                'Visible', 'on');
cor_pos_corr_text = uicontrol('Style', 'text', 'String', 'COR position correction:', ...
                'Position', [10 20 120 40], 'FontSize', 11, 'Visible', 'on');
            
cor_pos_corr_units_text = uicontrol('Style', 'text', 'String', 'mm', ...
                'Position', [185 20 30 40], 'FontSize', 11, 'Visible', 'on');
            
function corposcorrchanged(source, event)
    
    cor_pos_corr_prev = cor_pos_corr; % Save the previous source to center distance
    cor_pos_corr = str2double(get(source, 'String'));
    
    if(isnan(cor_pos_corr))
        cor_pos_corr = cor_pos_corr_prev; 
        fprintf('The COR position correction must be a number\n');
    end
    
    if (cor_pos_corr <  0)
        fprintf('The COR position correction must be positive and non-zero.\n');
        cor_pos_corr = cor_pos_corr_prev;
    end
    
    fprintf('The COR position correction is %1.2f mm\n', cor_pos_corr);
    set(cor_pos_corr_field, 'String', cor_pos_corr);
    
end

%% Save reconstructions tick-box
%% ------------------------------------------------------------------------
% Defaults
% --------
save_freq = 1;
saving_enabled = false;

% Check-boxes
% -----------
saving_checkbox = uicontrol('Style', 'checkbox', 'Value', 0, ...
    'Position', [250 125 175 20], 'String', 'Save reconstructions', 'FontSize', 11, ...
    'Callback', @savingenabled);

saving_freq_text = uicontrol('Style', 'text', 'String', 'Save every ', ...
    'Position', [250 105 100, 20], 'FontSize', 11, 'Visible', 'off');

saving_freq_field = uicontrol('Style', 'edit', 'String', save_freq, ...
    'Position', [355 105 50 20], 'FontSize', 11, 'Callback', @savefreqchanged, ...
    'Visible', 'off');

saving_freq_units_text = uicontrol('Style', 'text', 'String', 'iterations', ...
    'Position', [410 105 70 20], 'FontSize', 11,...
    'Visible', 'off');

function savingenabled(source, event)
    
    saving_value = get(source, 'Value');
    
    if (saving_value == 1); 
        saving_enabled = true; 
        set(saving_freq_text, 'Visible', 'on');
        set(saving_freq_field, 'Visible', 'on');
        set(saving_freq_units_text, 'Visible', 'on');
        fprintf('Saving enabled\n');
    else 
        saving_enabled = false; 
        set(saving_freq_text, 'Visible', 'off');
        set(saving_freq_field, 'Visible', 'off');
        set(saving_freq_units_text, 'Visible', 'off');
        fprintf('Saving disabled\n');
    end
end

function savefreqchanged(source, event)
    
    save_freq_prev = save_freq;
    save_freq = round(str2double(get(source, 'String')));
    
    % Check for non-sensical inputs
    if(isnan(save_freq))
        fprintf('Save frequency must be a number.\n');
        save_freq = save_freq_prev; 
    end 
    
    % Check for inputs that are too large or too small
    if (save_freq <= 0)
        fprintf('Save frequency must be positive and non-zero.\n');
        save_freq = save_freq_prev; 
    elseif (save_freq >= iterations)
        fprintf('Save frequency should not be greater than or equal to the number of iterations.\n');
        save_freq = save_freq_prev; 
    end 
    
    % Check that the iteration number is divisible by the save frequency 
    if (double(int32(iterations/save_freq)) ~= (iterations/save_freq))
        fprintf('The iteration number must be divisible by the save frequency otherwise no reconstructions will be saved.\n');
        save_freq = save_freq_prev;
    end
    
    fprintf('The save frequency is %1.0f\n', save_freq);
    set(saving_freq_field, 'String', save_freq);
    
end


%% Mirror projection images button
%% ------------------------------------------------------------------------
mirror_images = false; % Default

mirror_checkbox = uicontrol('Style', 'checkbox', 'Value', 0, ...
    'Position', [250 85 110 20], 'String', 'Mirror images', 'FontSize', 11, ...
    'Callback', @mirrorenabled);

function mirrorenabled(source, event)  
    mirror_value = get(source, 'Value');
    if (mirror_value == 1) % If the box is checked, mirror the images
        mirror_images = true; 
        fprintf('Mirror projetion images enabled\n');
    else % If it is unchecked, don't mirror them
        mirror_images = false;
        fprintf('Mirror projection images disabled\n');
    end
end

%% Chop down the data-set button
%% ------------------------------------------------------------------------
% Defaults
% --------
chop_data = false; 
chop_data_factor = 2;

% GUI buttons & text fields 
% -------------------------
chop_data_checkbox_use = uicontrol('Style', 'checkbox', 'Value', 0, ...
    'Position', [365 85 150 20], 'String', 'Use 1/', 'FontSize', 11, ...
    'Callback', @chopenabled);

chop_data_text = uicontrol('Style', 'text', 'Value', 0, ...
    'Position', [465 85 60 20], 'String', ' of data', 'FontSize', 11);

chop_data_factor_text = uicontrol('Style', 'popupmenu', 'String', {2, 4, 5, 10}, ...
    'Position', [430 95 40 15], 'FontSize', 11, ...
    'Callback', @factorchanged, 'Visible', 'off');

% What to do when check-box is ticked 
% -----------------------------------
function chopenabled(source, event)  
    % If the check-box is ticked, make the drop-down box appear where you
    % can select what factor the data should be chopped down by
    % Set the variables to the correct values for this functionality
    % -------------------------------------------------------------------
    ten_value = get(source, 'Value');
    if (ten_value == 1)
        chop_data = true; 
        chop_data_factor = 2; % Set back to default value
        set(chop_data_factor_text, 'Visible', 'on');
        set(chop_data_factor_text, 'Value', 1); % Make the pop-up menu text reflect this - 
                                                % 1 as 2 is the 1st element in the array
        fprintf('Using 1/%1.0f of projection images\n', chop_data_factor);
    % If the box is unchecked then make the drop-down box disappear
    % Set the variables to the correct values for this functionality
    % -------------------------------------------------------------
    else 
        chop_data = false;
        chop_data_factor = 1;
        set(chop_data_factor_text, 'Visible', 'off');
        fprintf('Using full projection data\n');
    end
end

% What to do when the sampling factor is changed 
% ----------------------------------------------
function factorchanged(source, event)
    factor_options = get(source, 'String');
    factor_index = get(source, 'Value');
    chop_data_factor = str2double(factor_options{factor_index});
    fprintf('Using 1/%1.0f of projection images\n', chop_data_factor);
end


%% ========================================================================
%% Start reconstruction button 
%% ========================================================================

% Button
% ------
start_recon = uicontrol('Style', 'pushbutton', 'String', 'Reconstruct!', 'FontSize', 22, ...
    'Position', [250 20, 250, 60], 'Callback', @startreconstruction);

% Function that starts reconstruction
% -----------------------------------
function startreconstruction(source, event)
    
    fprintf('Starting laminography reconstruction!');
    
    % Collate all the variables into specific groups ready to parse to the reconstruction function 
    % ---------------------------------------------------------------------------------------------
    
    % Parameters
    % ----------
    parameters.Alpha = tilt_angle;
    parameters.Iterations = iterations; 
    parameters.COR = padsize_cor;
    parameters.ROI = padsize_roi;
    parameters.BinningFactor = binning;
    parameters.Algorithm = method_str;
    parameters.ChopDataFactor = chop_data_factor;

    % Dimensions
    % -----------
    dimensions.Nx = nx;
    dimensions.Ny = ny;
    dimensions.Nz = nz;
    dimensions.MountHeight = mount_height;
    dimensions.CORPositionCorrection = cor_pos_corr;

    % File paths
    % ----------
    file_paths.DataPath = working_dir;
    file_paths.AnglesFile = angles_file_name;
    file_paths.SinogramFile = sino_file_name; 
    file_paths.ScanDimensionsFile = scan_dimensions_file_name;

    % Control variables 
    % -----------------
    ctrl_variables.SaveReconstruction = saving_enabled;
    ctrl_variables.SaveFrequency = save_freq;
    ctrl_variables.MirrorImages = mirror_images;
    ctrl_variables.LoadSinogram = load_sino;

     % Start the reconstruction
     % ------------------------
     % This function returns the reconstructed volume and the projection
     % geometry
    [V, proj_geom] = laminography_reconstruct(parameters, dimensions, ctrl_variables, file_paths);
    show_slices(V);
   
    % Save the reconstructed volume and the projection geometry
    % ---------------------------------------------------------
    % The reconstruction volume is saved as a .mat file and also a .vol
    % file
    % The projection geometry is just saved as a .mat file
    date = datestr(now,'mmmm_dd_yyyy_HH_MM_AM');

    file_name = strcat('reconstruction', date);
    file_name_vol = strcat(file_name, '.vol');
    file_name_proj = strcat('proj_geom', date);

    save(file_name, 'V', '-v7.3');
    save(file_name_proj, 'proj_geom', '-v7.3');
    % write_vol(V, working_dir, file_name, 'laminography', 'double');
    fid = fopen(file_name_vol, 'w');
    fwrite(fid, V, 'double');
    fclose(fid);
    
    fprintf('Reconstruction saved as .mat file: %s.mat\n', file_name);
    fprintf('Also saved as a .vol file: %s\n', file_name_vol);
    fprintf('Projection geometry saved as a .mat file: %s.mat\n', file_name_proj);
   
end
%% ========================================================================

end