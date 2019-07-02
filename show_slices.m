% SHOW_SLICES
% -------------------------------------------------------------------------
% A function for viewing reconstructions output from the reconstruction
% algorithm as a 3D MATLAB matrix. The program allows you to view 2D
% slices through the 3D volume in the x, y and z planes. You can use the
% buttons to flick through the slices or just type what slice you want to
% view. You can adjust the grayscale and save the slices as png or tiff
% files.
% 
% TESTED ENVIRONMENT:
%   Windows 10, MATLAB version 2017a
%
% INPUTS:
%   rec - reconstruction as a 3D Matlab matrix
%
% OUTPUTS:
%   no outputs
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

function show_slices(rec)

f1 = figure('Name', 'Slice viewer');

N = length(rec(:, 1, 1));
image = zeros(N, N);

max_X = length(rec(:, 1, 1));
max_Y = length(rec(1, :, 1));
max_Z = length(rec(1, 1, :));

slice_number = max_Z/2;
resolution = 1;

i = slice_number/resolution;

cmax = max(max(rec(:, :, i*resolution)));
cmin = 0;

orientation = 'z';

if (strcmpi(orientation, 'x'))
    max_i = max_X;
elseif(strcmpi(orientation, 'y'))
    max_i = max_Y;
elseif(strcmpi(orientation, 'z'))
    max_i = max_Z;
else 
    error('Orientation must be either x, y or z');
end

        
save_btn_png = uicontrol('Style', 'pushbutton', 'String', 'Save as .png', ...
    'Position', [10 380 70 25], 'Callback', @saveimagepng);

save_btn_tiff = uicontrol('Style', 'pushbutton', 'String', 'Save as .tiff', ...
    'Position', [10 350 70 25], 'Callback', @saveimagetiff);

up_btn = uicontrol('Style', 'pushbutton', 'String', 'Up', ...
    'Position', [10 320 70 25], 'Callback', @up);

down_btn = uicontrol('Style', 'pushbutton', 'String', 'Down', ...
    'Position', [10 290 70 25], 'Callback', @down);

orientation_x = uicontrol('Style', 'pushbutton', 'String', 'x', ...
    'Position', [10 260 20 25], 'Callback', @x);

orientation_y = uicontrol('Style', 'pushbutton', 'String', 'y', ...
    'Position', [33 260 20 25], 'Callback', @y);

orientation_z = uicontrol('Style', 'pushbutton', 'String', 'z', ...
    'Position', [56 260 20 25], 'Callback', @z);

slice_field_text = uicontrol('Style','text',...
            'Position',[10 25 100 15],...
            'String', 'Type slice number: ');
        
slice_field = uicontrol('Style', 'edit', 'String', slice_number, ...
    'Position', [120 25 50 15], 'Callback', @slice);

axes_text = uicontrol('Style','text',...
            'Position',[10 200 70 20],...
            'String', 'Colour axes');

cmax_field_text = uicontrol('Style','text',...
            'Position',[25 175 40 20],...
            'String', 'Max');
        
cmax_field = uicontrol('Style', 'edit', 'String', cmax, ...
    'Position', [25 150 40 25], 'Callback', @setcmax);

auto_cmax = uicontrol('Style', 'pushbutton', 'String', 'auto', ...
    'Position', [25 135 40 15], 'Callback', @cmaxauto);

cmin_field_text = uicontrol('Style','text',...
            'Position',[25 100 40 20],...
            'String', 'Min');
        
cmin_field = uicontrol('Style', 'edit', 'String', cmin, ...
    'Position', [25 80 40 25], 'Callback', @setcmin);

auto_cmin = uicontrol('Style', 'pushbutton', 'String', 'auto', ...
    'Position', [25 65 40 15], 'Callback', @cminauto);


display_image();


    function saveimagetiff(source, event)

        image = image/max(max(max(image)));
        imwrite(image, strcat('slice', num2str(i*resolution), '_', orientation, '.tiff'));

    end

    function saveimagepng(source, event)

        image = image/max(max(max(image)));
        imwrite(image, strcat('slice', num2str(i*resolution),'_', orientation, '.png'));

    end

    function up(source, event)
        if (i <= (max_i - resolution))
            i = i + resolution;
            display_image();
        end
    end

    function down(source, event)
        if (i > resolution)
            i = i - resolution; 
            display_image();
        end
    end

    function display_image()
       
        if (strcmpi(orientation,'z'))
            image = squeeze(rec(:,:,i*resolution));
            imshow(image,[]); 
        elseif (strcmpi(orientation, 'x'))
            image =squeeze(rec(i*resolution, :, :));
            imshow(image, []); 
        elseif (strcmpi(orientation, 'y'))
            image = squeeze(rec(:, i*resolution, :));
            imshow(image, []); 
        else 
            error('Orientation must be either x, y or z');
        end
        
        slice_str = strcat('Slice', num2str(i*resolution));
        orientation_str = strcat(orientation, ' axis');
    
        slice_txt = uicontrol('Style','text',...
            'Position',[300 25 100 15],...
            'String',slice_str);
        
        orientation_txt = uicontrol('Style','text',...
             'Position',[300 10 100 15],...
             'String', orientation_str);
        
        colorbar
       
        caxis([cmin cmax]);
        
    end

    function x(source, event)
        if (i*resolution <= max_X)
           orientation = 'x'; 
           max_i = max_X;
           display_image();
        else
            fprintf('Slices must be between 1 and %1.0f\n', max_X);
        end
    end

    function y(source, event)
        if (i*resolution <= max_Y)
            orientation = 'y';
            max_i = max_Y;
            display_image();
        else
            fprintf('Slices must be between 1 and %1.0f\n', max_Y);
        end
    end

    function z(source, event)
        if(i*resolution <= max_Z)
            orientation = 'z';
            max_i = max_Z;
            display_image();
        else
            fprintf('Slices must be between 1 and %1.0f\n', max_Z);
        end
    end

    function slice(source, event)
       
        slice_string = get(source, 'String');
        try 
            slice_no = str2double(slice_string);
            if (slice_no > 0 && slice_no <= max_i)
                i = slice_no/resolution;
                slice_number = slice_no;
                display_image();
            else
                fprintf('Slices must be between 1 and %1.0f\n', max_i);
            end
        catch
            fprintf('Please enter a number as the slice number.');
        end
    end

    function setcmax(source, event)
       cmax_str = get(source, 'String');
        try
            cmax_num = str2double(cmax_str);
            cmax = cmax_num;
            cmax_field = uicontrol('Style', 'edit', 'String', cmax, ...
                    'Position', [25 150 40 25], 'Callback', @setcmax);
            try
                display_image();
            catch
                fprintf('Maximum axis value must be larger than the minimum axis value.\n');
            end
        catch
            fprintf('Maximum axis value must be a number.\n');
        end    
    end

    function setcmin(source, event)
       cmin_str = get(source, 'String');
        try
            cmin_num = str2double(cmin_str);
            cmin = cmin_num;
            cmin_field = uicontrol('Style', 'edit', 'String', cmin, ...
                    'Position', [25 80 40 25], 'Callback', @setcmin);
            try
                display_image();
            catch
                fprintf('Minimum axis value must be smaller than the maximum axis value.\n');
            end
        catch
            fprintf('Minimum axis value must be a number.\n');
        end    
    end

    function cmaxauto(source, event)
        cmax = max(max(image));
        cmax_field = uicontrol('Style', 'edit', 'String', cmax, ...
                    'Position', [25 150 40 25], 'Callback', @setcmax);
        try
            display_image();
        catch
            fprintf('Maximum axis value must be larger than the maximum axis value.\n');
        end
    end

    function cminauto(source, event)
        cmin = 0;
        cmin_field = uicontrol('Style', 'edit', 'String', cmin, ...
                'Position', [25 80 40 25], 'Callback', @setcmin);
        try
            display_image();
        catch
            fprintf('Minimum axis value must be smaller than the maximum axis value.\n');
        end
    end
end

