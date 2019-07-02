% LAMINOGRAPHY_RECONSTRUCTION_CODE
% -------------------------------------------------------------------------
% Function used to perform the laminography reconstruction 
% This function can be run without the interface by specifiying the data
% structures given below with the correct fields. Alternatively the
% function can be run through the laminography interface (run
% lamino_interface). This method is more user friendly. 
%
% TESTED ENVIRONMENT:
%   Windows 10, MATLAB version 2017a
% 
% HARDWARE REQUIREMENTS:
%   This code runs on a single GPU
%
% SOFTWARE REQUIREMENTS:
%   - ASTRA version 1.8 (available from Github)
%   - CUDA driver version 8 
%
% INPUT:
%   parameters - A data structure containing the scanning geometry parameters
%   The data structure must have the following fields:
%          - parameters.Alpha - the tilt angle given by the scanner
%          - parameters.Iterations - the number of iterations of the
%          - parameters.COR - the centre of rotation correction in pixels
%          - parameters.ROI - the region of interest correction in pixels 
%          - parameters.BinningFactor - the factor the projection iThis mages have been binned by 
%          (if binning has already been done by the scanner during data collection)
%          - parameters.Algorithm - the reconstruction algorithm (either
%          CGLS or SIRT)
%          - parameters.ChopDataFactor - the sampling factor if you wish to
%          use a data set of a reduced size for performance reasons or
%          during testing
%   dimensions - A data structure containing the reconstruction volume
%   dimensions. The data structure must have the following fields:
%          - dimensions.Nx - number of pixels in the x dimension in the reconstruction
%          - dimensions.Ny - number of pixels in the y dimension in the
%          reconstruction
%          - dimensions.Nz - number of pixels in the z dimension of the
%          reconstruction
%          - dimensions.MountHeight - the length of the mount used in mm
%          - dimensions.CORPositionCorrection -Depending on the machine, 
%          the centre of rotation of the manipulator stage may or may not be 
%          at the bottom of the sample mount. To adjust for this, a correction
%          to the source-to-center distance in mm should be added. A positive
%          distance shortens the source-to-center distance (sample moves closer 
%          to the source) by the amount specified, and a negative distance lengthens it. 
%   ctrl_variables - A data structure containing variables to toggle 
%   between different behaviours of the script and to control what is
%   output
%           - ctrl_variables.SaveReconstruction - True/False if
%           reconstruction should be saved during the interations. The end
%           result be saved regardless of the value of this field. This
%           should be used if you want to save during the iterations to
%           monitor how the quality of the reconstruction changes with
%           iteration number.
%           - ctrl_variables.SaveFrequency - How often you want to save the
%           iterations (depends on SaveReconstruction being True). For
%           example a value of 5 means saving the reconstruction every 5
%           iterations. 
%           - ctrl_variables.MirrorImages - Sometimes the direction in which
%           the stage rotates differs from machine to machine and you may need to 
%           reverse the order of the projection images in order to get a reconstruction
%           that makes sense. Set this to True to do this.
%           - ctrl_variables.LoadSinogram - True/False to load a pre-made
%           sinogram from a file 
%           - ctrl_variables.ChopData - True/False if you want to use a
%           sampled data set to speed up reconstruction.
%   file_paths - A data structure containing file paths to different files
%   that are needed in the reconstruction
%           - file_paths.DataPath - the file path to the projection data
%           - file_paths.AnglesFile - the file path to the projection
%           angles file 
%           - file_paths.SinogramFile - the file path to the sinogram
%           - file_paths.ScanDimensionsFile - the file path to the folder
%           containing the dimensions file output from the scanner
% 
% OUTPUT:
%   V - The reconstructed volume as a 3D Matlab matrix
%   proj_geom - The projection geometry used for reconstruction as a data
%   structure
%
% -------------------------------------------------------------------------
% Copyright (c) 2019 S Fisher and D Holmes 
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

%% Beware: Clear existing variable

function [V, proj_geom] = laminography_reconstruct(parameters, dimensions, ctrl_variables, file_paths)

% Check that ASTRA has been set-up (should be automatically set-up if
% using the interface
% -------------------------------------------------------------------------
try
    astra_clear
catch ex
    error('Remember to run setup_recon_laminography.m before running this script.\n');
end

% Assign all the parameters from the interface to variables that are used
% in the function. Read in any data from files that have been given.
% -------------------------------------------------------------------------
alpha = 90 - parameters.Alpha; % (90 - tilt angle specified by the scanner)
iterations = parameters.Iterations; % Number of iterations 
padsize_cor = parameters.COR; % Number of pixels to add for the center of rotation correction
padsize_roi = parameters.ROI; % Number of pixels to add for region of interest correction
datapath = file_paths.DataPath; % Current folder
data_load = ctrl_variables.LoadSinogram;
pixel_bin_factor = parameters.BinningFactor ; % this specifies the binning factor relative to the raw data from the XTH scanner. 
% (may be different depending on scanner being used - check the geometry files)
N_x = dimensions.Nx; % Number of pixels in the x' direction in the reconstruction volume
N_y = dimensions.Ny; % Number of pixels in the y' direction in the reconstruction volume
N_z = dimensions.Nz; % Number of pixels in the z' direction in the reconstruction volume (can be less than N as a thin object so most of reconstruction volume
% will be empty, so should not reconstruct to save time)

% Control variables
algorithm_type = parameters.Algorithm; % can be either CGLS or SIRT
save_recon = ctrl_variables.SaveReconstruction; % true if the reconstruction is to be saved
save_freq = ctrl_variables.SaveFrequency; % Save every save_freq iterations
minus_angles = ctrl_variables.MirrorImages; % specifies whether there needs to be a minus sign when using the angles.
% Necessary if the projection images have been mirrored.
chop_data_factor = parameters.ChopDataFactor;

% Scan dimensions - read in from the given .xtek file
% -------------------------------------------------------------------------
% The line numbers for the important parameters are: 
% source-origin distance - line 16 of xtek file 
% source-detector distance - line 17 of xtek file 
% detector pixel size - line 21 of xtek file 
% If your xtek file has these parameters on a different line you will need
% to update the index in the definition of these parameters below to the
% correct line numbers for your file.
fileID_dimensions = fopen(file_paths.ScanDimensionsFile, 'r');
data_from_file = textscan(fileID_dimensions, '%s %s', 'Delimiter', '=');
scan_dimensions = string(data_from_file{1, 2});
fclose(fileID_dimensions);

source_origin = double(scan_dimensions(16)); % Source to COR distance in mm
source_detector = double(scan_dimensions(17)); % Source to detector distance in mm
mount_height = dimensions.MountHeight; % Height of polystyrene mount used in mm
detector_pixel_size = double(scan_dimensions(21)); % Detector pixel size in mm given in the xtec file
corr_pos_correction = dimensions.CORPositionCorrection; % A correction to correct for the fact that the COR does not lie on the mount

% Calculate the detector pixel size in mm taking into account image binning
detector_pixel_size = detector_pixel_size * pixel_bin_factor;

% Angles - read in from given .ang file 
% -------------------------------------------------------------------------
angles_file = file_paths.AnglesFile;
fileID_angles = fopen(angles_file, 'r');
data_from_angles_file = textscan(fileID_angles, '%s %s', 'Delimiter', ':');
angles=double(string(data_from_angles_file{1, 2 }));
projection_no = length(angles)-1
angles = angles(2:projection_no+1);
fclose(fileID_angles);

% Sinogram file path 
% -------------------------------------------------------------------------
sinogram_file = file_paths.SinogramFile;

% Summary of variables - output to command line  
% -------------------------------------------------------------------------
fprintf('A laminography reconstruction will be performed using the following parameters:\n');
fprintf('alpha = %1.1f degrees\n', alpha);
fprintf('(Tilt angle = %1.1f degrees\n)', (90-alpha));
fprintf('Number of iterations = %1.0f\n', iterations);
fprintf('Center of rotation correction = %1.0f pixels\n', padsize_cor);
fprintf('Region of interest correction = %1.0f pixels\n', padsize_roi);
fprintf('Height of sample mount = %1.2f mm \n', mount_height);
fprintf('Source to center distance = %1.2f mm \n', source_origin);
fprintf('Source to detector distance = %1.2f mm \n', source_detector);
fprintf('The pixel bin factor is = %1.0f \n', pixel_bin_factor);
fprintf('Detector pixel size %1.4f mm \n', detector_pixel_size);
fprintf('Center of rotation position correction = %1.2f mm\n', corr_pos_correction);
fprintf('Reconstruction volume %1.0f x %1.0f x %1.0f\n', N_x, N_y, N_z);
fprintf('Projection angles file = %s\n', angles_file);
fprintf('Number of projection images = %1.0f\n', projection_no);
fprintf('Using 1/%1.0f of the projection data\n', chop_data_factor);
fprintf('The reconstruction will be performed using the %s algorithm\n', algorithm_type);

if (data_load == true)
    fprintf('Sinogram will be loaded from file = %s\n', sinogram_file);
end

if(save_recon == true)
    fprintf('Saving enabled every %1.0f iterations\n', save_freq);
end
fprintf('Data will be stored in = %s\n', datapath); 

fprintf('\n');
fprintf('Review the above information before the scan begins\n');
cont = input('Press any key to continue. Press c to cancel', 's');
if(strcmpi(cont, 'c'))
    V = 0;
    proj_geom = 0;
    return;
end
% -------------------------------------------------------------------------

tic % Start timing  

% Load saved sinogram data if selected
% ------------------------------------
% The sinogram that is loaded in should be of the dimensions rows x
% projection images x columns
if(data_load == true)
    disp(strcat('Loading sinogram data from file: ', sinogram_file));
    try 
        load(sinogram_file);
    catch ex
        error('No sinogram found in the workspace. You must load it from a file.\n');
    end
end

% log the sinogram
% ----------------
sinogram = -log(sinogram);

% If using binned data, sample the sinogram to the appropriate bin factor 
% -----------------------------------------------------------------------
if (chop_data_factor ~= 1)
    fprintf('Reducing sinogram size by factor %1.0f\n', chop_data_factor);
    sinogram = chop_down_sinogram(sinogram, chop_data_factor);
    angles = chop_down_angles(angles, chop_data_factor);
    projection_no = length(sinogram(1, :, 1)); % Recalculate the number of projection angles
end

save('angles', 'angles', '-v7.3');

N_x_sino = length(sinogram(:, 1, 1));
N_y_sino = length(sinogram(1, 1, :));
if (N_x_sino ~= N_y_sino) 
   error('Projection images must be square'); 
end

% Output sinogram dimension information
% -------------------------------------
fprintf('Information obtained from the sinogram file:');
fprintf('The number of projection images is %1.0f\n', projection_no);
fprintf('The image dimensions are %1.0f x %1.0f\n', N_x_sino, N_y_sino);

% Perform center of rotation correction
% -------------------------------------
if (padsize_cor ~= 0)
    disp('Performing centre of rotation corrections...');
    sinogram = pad_sinogram(sinogram, padsize_cor);
end


if (padsize_roi ~= 0)
    sinogram = pad_sinogram(sinogram, padsize_roi);
    sinogram = pad_sinogram(sinogram, (-1*padsize_roi));
end

%% Extract parameters for ASTRA

% Geometric magnification: Scaling factor mapping one object pixel in the
% center-of-rotation plane on to one detector pixel.
c = mount_height*cos(alpha*(pi/180));

% Calculate the distance from the source to the center of rotation, and the
% detector to the center of rotation

detector_origin = source_detector - source_origin;
source_origin = -1*(source_origin - c + corr_pos_correction); 
detector_origin = -1*(detector_origin + c - corr_pos_correction);
magnification = (source_origin+detector_origin)/source_origin;

% Pixel size in mm
pixel_size = detector_pixel_size / magnification;

% Detector pixel size in units of pixels 
detector_pixel_size = detector_pixel_size/pixel_size;

fprintf('Distance from source to center = %2.1f cm\n', abs(source_origin)); 
fprintf('Distance from detector to center = %2.1f cm\n', abs(detector_origin)); 
fprintf('Magnification = %2.1f\n', magnification); 

% Projections angles in radians 
% Reverse if projection images are mirrored
if (minus_angles == true)
    angles_astra = -(angles/180)*pi;

else
    angles_astra = (angles/180)*pi;
    
end

% Distance source to center-of-rotation, given in numbers of pixels.
source_center_dist = source_origin ./ pixel_size;

% Distance detector to center-of-rotation, given in numbers of pixels.
detector_center_dist = detector_origin ./ pixel_size;

%% Set up ASTRA volume geometry

disp('Assigning volume geometry...');
vol_geom = astra_create_vol_geom(N_x, N_y, N_z);
disp('Completed assigning vol_geom');

%% Set up ASTRA projection geometry

% Convert tilt angle to radians 
alpha = alpha*(pi/180);

% Preallocate the geometry specifying vectors
srcX = zeros(projection_no, 1); srcY = zeros(projection_no, 1); srcZ = zeros(projection_no, 1);
dX = zeros(projection_no, 1); dY = zeros(projection_no, 1); dZ = zeros(projection_no, 1);
uX = zeros(projection_no, 1); uY = zeros(projection_no, 1); uZ = zeros(projection_no, 1);
vX = zeros(projection_no, 1); vY = zeros(projection_no, 1); vZ = zeros(projection_no, 1);

% Create a 3D cone beam geometry specified by 3D vectors defining the
% laminography geometry 
% -------------------------------------------------------------------
for i = 1:length(angles_astra)
    
    srcX(i) = -1*source_center_dist*sin(alpha).*sin(angles_astra(i));
    srcY(i) = source_center_dist*sin(alpha).*cos(angles_astra(i));
    srcZ(i) = -1*source_center_dist*cos(alpha);

    dX(i) = detector_center_dist*sin(alpha).*sin(angles_astra(i));
    dY(i) = -1*detector_center_dist*sin(alpha).*cos(angles_astra(i));
    dZ(i) = detector_center_dist*cos(alpha);
    
    vX(i) = 1*detector_pixel_size.*cos(angles_astra(i));
    vY(i) = 1*detector_pixel_size.*sin(angles_astra(i));
    vZ(i) = 0;
   
    uX(i) = detector_pixel_size*cos(alpha).*sin((angles_astra(i)));
    uY(i) = -1*detector_pixel_size*cos(alpha).*cos(angles_astra(i));
    uZ(i) = -1*detector_pixel_size*sin(alpha);
    
end

vectors = [srcX srcY srcZ dX dY dZ uX uY uZ vX vY vZ];
disp('Assigning projection geometry...');
proj_geom = astra_create_proj_geom('cone_vec', N_y_sino+ abs(padsize_cor) + 2*padsize_roi , N_x_sino, vectors);
  
disp('Completed assigning projection geometry')

% Create a 3d object of the volume reconstruction for ASTRA
rec_id = astra_mex_data3d('create', '-vol', vol_geom);
disp('Completed mex_data3d');

%% Reconstruction algorithms

if (strcmpi(algorithm_type, 'SIRT'))
    
    % Calculate the R and C matrices for the SIRT algorithm
    fprintf('Calculating row and column sums...\n');
    % Start with a matrix of ones 
    M_R = ones(N_x, N_y, N_z);
    M_C = ones(N_x_sino, projection_no, N_y_sino+ abs(padsize_cor) + 2*padsize_roi);

    % Reduce the size of jobs sent to the GPU to stop the Windows time out
    % problem
    astra_mex('set_gpu_index', 0, 'memory', 1000e6);

    % Do a forward projection of this matrix 
    [R_id_inv, R_inv] = astra_create_sino3d_cuda(M_R, proj_geom, vol_geom);
    R = 1./R_inv;
    R(R==Inf) = 0;
    % Do a back projection of this matrix 
    [C_id, C_inv] = astra_create_backprojection3d_cuda(M_C, proj_geom, vol_geom);
    C = 1./C_inv;
    C(C==Inf) = 0;
    
    % Clean up memory
    fprintf('Cleaning up memory\n');
    astra_mex_data3d('delete', C_id);
    astra_mex_data3d('delete', R_id_inv);
    % Initial solution is a blank matrix
    V = zeros(N_x, N_y, N_z);

    fprintf('Reconstruction starting...\n');

    for i = 1:iterations;

        fprintf('Starting iteration %1.0f of %1.0f... \n', i, iterations);
        tic;

        % Reduce the size of jobs sent to the GPU to stop the Windows time out
        % problem
        astra_mex('set_gpu_index', 0, 'memory', 1000e6);

        % Perform a forward projection operation
        [sino_id, sinogram_iteration] = astra_create_sino3d_cuda(squeeze(V(:,:,:)), proj_geom, vol_geom);
        fprintf('Forward projection complete\n');

        % Find the difference between real sinogram data obtained and forward projection
        diff_raw = sinogram - sinogram_iteration;
        diff = R.*diff_raw;

        % Do a back projection of this difference
        [id, volume] = astra_create_backprojection3d_cuda(diff, proj_geom, vol_geom);
        fprintf('Back projection complete\n');

        % Get the new solution
        V(:,:,:) = V(:,:,:) + C.*volume;

        % Apply the non-negativity condition
        V(V<0) = 0; % Non-negativity

        fprintf('Iteration solution found\n');

        % Display the result of the iteration
        if (double(int32(i/50)) == i/50)
            figure('name', strcat('Reconstruction: Iteration ', num2str(i)));
            imagesc(squeeze(V(:, :, 100)));
            colormap gray
            colorbar 
        end

        % Garbage disposal 
        astra_mex_algorithm('delete', sino_id);
        astra_mex_data3d('delete', id);
        fprintf('Algorithm memory cleared\n');

        time = toc; % Record the time for each iteration for iteration timing estimation

        fprintf('Iteration %1.0f of %1.0f complete\n', i, iterations);

        if (i ~= iterations)
            fprintf('Approximate time remaining %1.0f min\n', ((iterations-i)*time)/60);
        end


        % Save the result every 10 iterations

        if (save_recon == true)
            if (double(int32(i/save_freq)) == i/save_freq)
                fprintf('Saving reconstruction...\n');
                file_name = strcat('Iteration_', algorithm_type);
                save(file_name, 'V', '-v7.3');
            end
        end

    end

elseif (strcmpi(algorithm_type, 'CGLS'))
    
    % Reduce the size of jobs sent to the GPU to stop the Windows time out
    % problem
    astra_mex('set_gpu_index', 0, 'memory', 1000e6);
        
    % Prep for iteration
    V = zeros((N_x*N_y*N_z),1);
    sinogram_vec = sinogram(:);
    [id_1, d] = astra_create_backprojection3d_cuda(sinogram, proj_geom, vol_geom);
    d = d(:);
    normr2 = d'*d;

    % Start the iteration
    for i = 1:iterations

        fprintf('Starting iteration %1.0f of %1.0f... \n', i, iterations);
        tic;

        % Reduce the size of jobs sent to the GPU to stop the Windows time out
        % problem
        astra_mex('set_gpu_index', 0, 'memory', 1000e6);

        % Appy the CGLS algorithm
        d_array = reshape(d,[N_x, N_y, N_z]);
        [id_2, Ad_array] = astra_create_sino3d_cuda(d_array, proj_geom, vol_geom);
        Ad = Ad_array(:);
        alpha = normr2./(Ad'*Ad);
        V = V + alpha.*d;
        sinogram_vec = sinogram_vec - alpha.*Ad;
        sinogram_iteration = reshape(sinogram_vec, [N_x_sino, projection_no, N_y_sino + abs(padsize_cor) + 2*padsize_roi]); 
        [id_3, s_array] = astra_create_backprojection3d_cuda(sinogram_iteration, proj_geom, vol_geom);

        s = s_array(:);
        normr2_new = s'*s;
        beta = normr2_new./normr2;
        normr2 = normr2_new;
        d = s + beta.*d;

        % Garbage disposal 
        astra_mex_algorithm('delete', id_2);
        astra_mex_data3d('delete', id_1);
        astra_mex_data3d('delete', id_3);
        fprintf('Algorithm memory cleared\n');

        fprintf('Iteration %1.0f of %1.0f complete\n', i, iterations);

        time = toc;
        if (i ~= iterations)
            fprintf('Approximate time remaining %1.0f min\n', ((iterations-i)*time)/60);
         end

        % Save the result every save_freq iterations
        if (save_recon == true)
            if (double(int32(i/save_freq)) == i/save_freq)
                fprintf('Saving reconstruction...\n');
                V_recon = reshape (V, [N_x,N_y ,N_z]); 
                V_recon(V_recon<0) = 0;
                file_name = strcat('Iteration_', algorithm_type);
               
                save(file_name, 'V_recon', '-v7.3');
            end
        end

    end

    V = reshape (V, [N_x,N_y ,N_z]); 
end

total_time = toc;

% Write parameters to file
date = datestr(now,'mmmm_dd_yyyy_HH_MM_AM');

file_id = fopen(strcat('parameters', date, '.txt'), 'w');

if (file_id == -1)
    error('An error occurred whilst writing the log file');
end

% Output summary of reconstruction
fprintf(file_id, 'Alpha (90-tilt angle) : %1.1f degrees\n', alpha*(180/pi));
fprintf(file_id, 'Center of rotation correction : %1.0f pixels\n', padsize_cor);
fprintf(file_id, 'Region of interest correction : %1.0f pixels\n', padsize_roi);
fprintf(file_id, 'Algorithm : %s\n', algorithm_type);
fprintf(file_id, 'Iterations : %1.0f\n', iterations);
fprintf(file_id, '\n');
fprintf(file_id, 'Number of projection images: %1.0f\n', projection_no);
fprintf(file_id, 'Projection images have been binned by a factor %1.0f \n', pixel_bin_factor);
fprintf(file_id, 'Projection image dimensions: %1.0fx%1.0f pixels\n', N_x_sino, N_y_sino);
fprintf(file_id, 'Reconstruction volume : %1.0fx%1.0fx%1.0f pixels\n', N_x, N_y, N_z);
fprintf(file_id, 'Source to center distance: %1.2f mm \n', source_origin);
fprintf(file_id, 'Source to detector distance: %1.2f mm \n', source_detector);
fprintf(file_id, 'Mount height : %1.2f mm \n', mount_height);
fprintf(file_id, 'Center of rotation vertical position correction %1.2f mm \n', corr_pos_correction); 
fprintf(file_id, 'Detector pixel size: %1.2f mm \n', detector_pixel_size);
fprintf(file_id, 'Voxel size: %1.4f mm \n', pixel_size);
fprintf(file_id, '\n');
fprintf(file_id, 'Data saved in file %s \n', datapath);
fprintf(file_id, 'Angle data from file %s \n', angles_file);
fprintf(file_id, 'Sinogram from file %s \n', sinogram_file);
fprintf(file_id, 'Scan dimensions from file %s \n', file_paths.ScanDimensionsFile);
fprintf(file_id, '\n');
fprintf(file_id, 'Save reconstruction : %1.0f\n', save_recon);
fprintf(file_id, 'Mirror images : %1.0f \n', minus_angles);
fprintf(file_id, '\n \n');
fprintf(file_id, 'Total time taken for reconstruction : %1.3f hours\n', (total_time/(3600)));

fclose(file_id);

end
