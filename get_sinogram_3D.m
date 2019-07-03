% GET_SINOGRAM_3D
% -------------------------------------------------------------------------
% Apply a flat field correction and any projection binning and then create 
% a 3D sinogram from 2D projection data.
% 
% TESTED ENVIRONMENT:
%   Windows 10, MATLAB version 2017a
% 
% INPUTS:
%   proj_path - file path to directory where projection data is located
%   dark file - file path to the dark file 
%   light file - file path to the light file
%   bin_data - boolean for binning data True/False
%   bin_amount - the factor for binning the projection images by. Will only
%   have an affect if the bin_data variable is set to True.
% 
% OUTPUTS:
%   sino_dbl - the constructed sinogram
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

function sino_dbl = get_sinogram_3D(proj_path, dark_file, light_file, bin_data, bin_amount)

fprintf('Generating a sinogram from projection data located at: \n %s\n', proj_path);
fprintf('Using dark file : %s \n', dark_file);
fprintf('Using light file: %s \n', light_file); 

if (bin_data == true)
    fprintf('Binning of projection images is enabled with binning factor %1.0f\n', bin_amount);
else
    fprintf('Binning of projection images is disabled\n');
end

D = dir([proj_path, '\*.tif']);
old_dir = cd(proj_path);
N = length(D(not([D.isdir])));

% Read in a test image to get image dimensions
image = imread( D(1, 1).name); % Read in the image
N_x = length(image(:, 1));
N_y = length(image(1, :));

% Read in light and dark files
dark = imread(dark_file);
light = imread(light_file);


if (bin_data == true)
    sino = zeros(N_x/bin_amount, N_y/bin_amount, N);
else
    sino = zeros(N_x, N_y, N);
end

% Read in each projection image, apply a background correction, bin the
% image if required, and add to the sinogram matrix
for j = 1:N
    
    tic
    a = D(j, 1).name; % Get the file name
    A = imread(a); % Read in the image
    cd(old_dir);
    
    % Background correction
    A = ((double(A)-double(dark))./(double(light)-double(dark)))*(65535);
    
    % Bin data
    if (bin_data == true)
        A = imresize(A, [N_x/bin_amount, N_y/bin_amount], 'nearest');
    end
    
    cd(proj_path);
    % Add to singoram
    sino(:, :, j ) = A(:, :);
    
    time = toc;
    
    % Estimate time remaining
    if ( double(int32(j/10)) == j/10)
        fprintf('Image %1.0f of %1.0f\n', j, N);
        fprintf('Approximate time remaining in this data set %1.0f min\n', (time*(N-j)/60));
    end
    
end

sino_dbl = double(sino)/65535;
sino_dbl = permute(sino_dbl, [1 3 2]); % Permute so in the right order for the reconstruction code

cd(old_dir);

