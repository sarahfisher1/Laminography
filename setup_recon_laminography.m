% SETUP_RECON_LAMINOGRAPHY
% -------------------------------------------------------------------------
% Set up the file paths and add the GPU in preparation for laminography 
% reconstruction.
% 
% TESTED ENVIRONMENT:
%   Windows 10, Matlab 2017a 
%
% INPUTS:
%   no inputs
% OUTPUTS:
%   status - the gpu status
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

function status = setup_recon_laminography()

% Add the Astra toolbox version 1.8. The Astra toolbox folder must be in the
% current working directory
addpath(fullfile('astra-1.8','mex'))
addpath(fullfile('astra-1.8','tools'))

% Add the GPU
gpu_index = 1;
status = gpuDevice(gpu_index);

end


