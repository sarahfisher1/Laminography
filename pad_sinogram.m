% PAD_SINOGRAM
% -------------------------------------------------------------------------
% Add pixels to the sinogram. Used for adding centre of rotation or region
% of interest corrections 
% 
% TESTED ENVIRONMENT:
%   Windows 10, MATLAB version 2017a
% 
% INPUTS:
%   sinogram - the sinogram to add the correction to
%   cor_padsize - the number of pixels to add to make the correction
% 
% OUTPUTS:
%   sinogram_ready - the sinogram with the correction added
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

function sinogram_ready = pad_sinogram(sinogram, cor_padsize)

N_x = length(sinogram(:, 1, 1));
N_y = length(sinogram(1, 1, :)) + abs(cor_padsize);
proj_no = length(sinogram(1, :, 1));

sinogram_ready = zeros(N_x, proj_no, N_y);

for i = 1:proj_no
    
    tic
    sinogram_1 = squeeze(sinogram(:, i, :));
    
    if cor_padsize > 0
        sinogram_cor = [repmat(sinogram_1(:,1),1,cor_padsize), sinogram_1];
    else
        sinogram_cor = [sinogram_1, repmat(sinogram_1(:,end),1,-cor_padsize)];
    end
  
    sinogram_ready(:, i, : ) = sinogram_cor(:, :);

    time = toc;
    
    if (double(int32(i/100)) == i/100)
        
        fprintf('Projection %1.0f of %1.0f. Estimated time remaining %2.1f s\n', i, proj_no, (time*(proj_no-i)));
    end
    
end

 
end