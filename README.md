# Laminography reconstruction code
MATLAB code for laminography reconstruction using SIRT and CGLS

This software is used for reconstruction of laminography data from a standard microCT scanner. This software was developed and used in the scanning and reconstruction process of a test object using the laminography method. This work was published in the the journal Measurement Science and Technology titled Laminography in the lab: Imaging planar objects using a conventional x-ray CT scanner (2019) (see https://iopscience.iop.org/article/10.1088/1361-6501/aafcae). 
The data-set for the original work can be downloaded from:  https://doi.org/10.5281/zenodo.2540509.

COPYRIGHT NOTICE
---------------------------------------------------------------------------------------------------------------------------------
Copyright (c) 2019 S Fisher and D Holmes 
University of Manchester
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

TESTED ENVIRONMENT
---------------------------------------------------------------------------------------------------------------------------------
  -> 3.10 GHz Intel Xeon E5-2687W v3 processor with 256 GB RAM 
  
  -> GPU (NVidia Quadro K6000)
  
  -> Operating system: Windows 10
  
  -> Matlab version 2017a
  
  -> CUDA driver version 8
  
  -> ASTRA Toolbox for MATLAB (version 1.8)
  
LAMINOGRAPHY SCAN DETAILS 
---------------------------------------------------------------------------------------------------------------------------------
The laminography scan was performed on a Nikon XTH 225 system. The scanning routine was programed using the IPC interface. The outputs of the scan were: the set of projection images in .tiff format, a light and dark image in .tiff format following calibration, a list of projection angles in .ang format, and the scanning geometry file in .xtekct format. Example .ang and .xtekct file are provided with the reconstruction code. 
  
REQUIRED SOFTWARE
---------------------------------------------------------------------------------------------------------------------------------
  -> MATLAB (Tested on version 2017a but may work on other recent versions)
  -> ASTRA Toolbox for MATLAB (version 1.8)
  -> CUDA driver version 8
  
RUNNING THIS SOFTWARE
---------------------------------------------------------------------------------------------------------------------------------
  -> Make sure the correct software is installed. The ASTRA Toolbox is available for download at: https://github.com/astra-toolbox/astra-toolbox. 
  
  -> Download the source code into your working directory.
  
  -> Put the ASTRA toolbox folder in your working directory.
  
  -> Put your .ang (projection angles file) and .xtekct (scanning geometry file) in your current working directory 
  
  -> Open MATLAB and in the MATLAB command window run the command "lamino_interface". This will bring up the laminography interface where you can set the reconstruction parameters, create a sinogram and start the laminography reconstruction. 
  
INTERFACE PARAMETER DEFINITIONS 
---------------------------------------------------------------------------------------------------------------------------------
  -> Tilt angle: This defines the inclination of the sample rotation axis relative to the vertical.
  
  -> Angles: Select the .ang file corresponding to the projection angles used in the relevant scan. The file must be in the current working directory.
  
  -> Reconstruction volume dimensions: Defines the number of pixels in the reconstruction volume in the x, y and z directions respectively .
  
  -> Method: Defines the iterative algorithm that will be used to reconstruct the data (can be SIRT aor CGLS).
  
  -> Iterations: The number of iterations of the chosen algorith that will be performed.
  
  -> COR: The number of pixels you wish to add to the detector to correct for any centre of rotation issues. Can be positive or negative. Positive values add pixels to the right hand side of the detector, while negative vaules add pixels to the left hand side.
  
  -> ROI: If performing a region of interest scan, the sinogram can be padded to remove the ring artefact that may otherwise be present. The value inputted specifies the number of pixels to pad the sinogram by. Both sides of the sinogram are padded by the amount specified.
  
  -> Dimension file: Select the relevent .xtexct file that contains the dimensions of the scan (source_detector distance, source_object dist, pixel size etc...). The file must be in the current working directory.
  
  -> Mount height: The height of the sample mount used in mm.
  
  -> COR position correction: Depending on the machine, the centre of rotation of the manipulator stage may or may not be at the bottom of the sample mount. To adjust for this, a correction to the source-to-center distance in mm should be added. A positive distance shortens the source-to-center distance (sample moves closer to the source) by the amount specified, and a negative distance lengthens it.
  
  -> Generate sinogram: Tick if you want to generate a 3D sinogram from raw data.
  
  -> Data: If generate sinogram is ticked, select the folder that 	contains the relevent raw projection data.
    
  -> Shading correction: If generate sinogram is ticked, Select 	the dark and light shading 	correction images.
    
  -> Load sinogram from file: Tick this box if you wish to load in a previously created sinogram. The sinogram must be saved in the current working directory.
  
  -> Binning factor: Select the binng factor you wish to apply to the raw data.
  
  -> Save reconstruction: Tick this box if you wish to save your reconstruction. If ticked, you must also select how frequently you wish the reconstruction to be saved ( every x iterations). The save frequency must be a multiple of the number of interations to allow regular saving intervals.
  
  -> Mirror images: Sometimes the direction in which the stage rotates differs between machines. If you are getting a nonsensical reconstruction, try ticking (or unticking) this box. This reverses the direction of rotation of the stage.
  
  -> Use 1/x of the data: If you do not wish to use the full dataset, select this and choose the factor by which you wish to reduce the number of projections by. This is useful when doing quick tests and are not so concerned about the quality of the reconstruction, to check the center of rotation correction, for example, as it drastically reduces the reconstruction time. 
  
  -> Reconstruct: Press when you have ensured that all parameters are correct and the necessary files have been selected. 

DURING THE SCAN
---------------------------------------------------------------------------------------------------------------------------------
First, you will be asked to review the parameters you have inputted on the command line in Matlab. You must press ENTER for the scan to begin once these paramters are displayed. You can press C to cancel the scan at this point (Once the reconstruction has started, the only way to cancel the scan is to force the program to stop using the Ctrl+C command. Any reconstructions that have been made up to the point that you cancel the reconstruction will be lost). The command line in MATLAB will constantly update with the reconstruction code's progress. Once iterations begin, an estimated time to completion is output after every iteration. 

SCAN OUTPUT
---------------------------------------------------------------------------------------------------------------------------------
After the scan completes, .mat and .vol files will be saved containing the resulting reconstruction, after the full number of iterations have passed. These will be named: reconstructionMonth_Day_Year_Time.mat and reconstructionMonth_Day_Year_Time.vol. 
The projection geometry will also be saved in the file proj_geomMonth_Day_Year_time.mat, as well as a .txt file containing details on all the scan parameters, saved as parametersMonth_Day_Year_Time.txt. 

VIEWING THE RECONSTRUCTION
---------------------------------------------------------------------------------------------------------------------------------
After the scan, the show_slices function will automatically run so that you can view orthogonal slices through the volume that has just been reconstructed. You can also run this after the scan on any volume (including sinograms) by running the show_slices_ function. You can also save images as .png or .tiff using this. Alternatively, you can easily load the .vol version of the reconstruction for full 3D rendering.

FUNCTION DEFINITIONS 
---------------------------------------------------------------------------------------------------------------------------------

  -> LAMINO_INTERFACE
  -------------------
Function that you call to run the interface. From the interface you can set the reconstruction parameters, corrections, algorithm, sinogram and start the reconstruction. This function contains the definitions of all the buttons, checkboxes and fields on the interface window and defines what happens when there is an interaction on the interface. 

INPUT:
   no inputs
 OUTPUT:
   no output
 
The function will save the reconstructed volume as a .mat file (MATLAB file) and a .vol file (file that can be read into visualisation programs such as Avizo). The projection geometry is saved as a .mat file. 

  -> LAMINOGRAPHY_RECONSTRUCT
  ---------------------------  
Function used to perform the laminography reconstruction. This function can be run without the interface by specifiying the data structures given below with the correct fields. Alternatively the function can be run through the laminography interface (run lamino_interface). This method is more user friendly. 

INPUT:
  parameters - A data structure containing the scanning geometry parameters
  
  The data structure must have the following fields:
  
         - parameters.Alpha - the tilt angle given by the scanner
         
         - parameters.Iterations - the number of iterations of the
         
         - parameters.COR - the centre of rotation correction in pixels
         
         - parameters.ROI - the region of interest correction in pixels 
         
         - parameters.BinningFactor - the factor the projection iThis mages have been binned by 
         (if binning has already been done by the scanner during data collection)
         
         - parameters.Algorithm - the reconstruction algorithm (either
         CGLS or SIRT)
         
         - parameters.ChopDataFactor - the sampling factor if you wish to
         use a data set of a reduced size for performance reasons or
         during testing
         
  dimensions - A data structure containing the reconstruction volume dimensions. The data structure must have the following fields:
  
         - dimensions.Nx - number of pixels in the x dimension in the reconstruction
         
         - dimensions.Ny - number of pixels in the y dimension in the
         reconstruction
         
         - dimensions.Nz - number of pixels in the z dimension of the
         reconstruction
         
         - dimensions.MountHeight - the length of the mount used in mm
         
         - dimensions.CORPositionCorrection -Depending on the machine, 
         the centre of rotation of the manipulator stage may or may not be 
         at the bottom of the sample mount. To adjust for this, a correction
         to the source-to-center distance in mm should be added. A positive
         distance shortens the source-to-center distance (sample moves closer 
         to the source) by the amount specified, and a negative distance lengthens it. 
         
  ctrl_variables - A data structure containing variables to toggle between different behaviours of the script and to control what is output
  
          - ctrl_variables.SaveReconstruction - True/False if
          reconstruction should be saved during the interations. The end
          result be saved regardless of the value of this field. This
          should be used if you want to save during the iterations to
          monitor how the quality of the reconstruction changes with
          iteration number.
          
          - ctrl_variables.SaveFrequency - How often you want to save the
          iterations (depends on SaveReconstruction being True). For
          example a value of 5 means saving the reconstruction every 5
          iterations. 
          
          - ctrl_variables.MirrorImages - Sometimes the direction in which
          the stage rotates differs from machine to machine and you may need to 
          reverse the order of the projection images in order to get a reconstruction
          that makes sense. Set this to True to do this.
          
          - ctrl_variables.LoadSinogram - True/False to load a pre-made
          sinogram from a file 
          
          - ctrl_variables.ChopData - True/False if you want to use a
          sampled data set to speed up reconstruction.
          
  file_paths - A data structure containing file paths to different files that are needed in the reconstruction
  
          - file_paths.DataPath - the file path to the projection data
          
          - file_paths.AnglesFile - the file path to the projection
          angles file 
          - file_paths.SinogramFile - the file path to the sinogram
          
          - file_paths.ScanDimensionsFile - the file path to the folder
          containing the dimensions file output from the scanner

OUTPUT:
  V - The reconstructed volume as a 3D Matlab matrix
  
  proj_geom - The projection geometry used for reconstruction as a data
  structure
  
  -> SETUP_RECON_LAMINOGRAPHY
  ---------------------------
Set up the file paths and add the GPU in preparation for laminography reconstruction.

INPUTS:
  no inputs
OUTPUTS:
   status - the gpu status 
  
  -> GET_SINOGRAM_3D
  ------------------
Apply a flat field correction and any projection binning and then create a 3D sinogram from 2D projection data.

INPUTS:
  proj_path - file path to directory where projection data is located
  
  dark file - file path to the dark file 
  
  light file - file path to the light file
  
  bin_data - boolean for binning data True/False
  
  bin_amount - the factor for binning the projection images by. Will only
  have an affect if the bin_data variable is set to True.

OUTPUTS:
  sino_dbl - the constructed sinogram

  -> CHOP_DOWN_SINOGRAM
  ---------------------
Function used to reduce the number of projection images in a sinogram by sampling. Use if you need to reduce the size of your data set for testing or performance reasons. Use in conjuction with chop_down_angles to reduce the number of projection angles in the corresponding projection angles list. 

INPUT:
  sinogram - a matrix containing the sinogram you want to sample from
  
  factor - the sampling factor 

OUTPUT:
  sinogram_chopped - the sampled sinogram
  
  -> CHOP_DOWN_ANGLES
  -------------------
Function used to reduce the number of projection angles by sampling. Use if you need to reduce the size of your data set for testing or performance reasons. Use in conjunction with chop_down_sinogram to reduce the number of projections in the corresponding sinogram.

INPUT:
  angles - An array containing the projection angles (in order)
  
  factor - The sampling factor

OUTPUT:
  angles_chopped - The sampled angles array

EXAMPLE:
  If the input angles array is [0, 1, 2, 3, 4, 5] and the sampling factor
  is 2 then the output array is [1, 3, 5].

  -> PAD_SINOGRAM
  ---------------
Add pixels to the sinogram. Used for adding centre of rotation or region of interest corrections 

INPUTS:
  sinogram - the sinogram to add the correction to
  cor_padsize - the number of pixels to add to make the correction

OUTPUTS:
  sinogram_ready - the sinogram with the correction added

  -> SHOW_SLICES
  --------------
A function for viewing reconstructions output from the reconstruction algorithm as a 3D MATLAB matrix. The program allows you to view 2D slices through the 3D volume in the x, y and z planes. You can use the buttons to flick through the slices or just type what slice you want to view. You can adjust the grayscale and save the slices as png or tiff files.

INPUTS:
  rec - reconstruction as a 3D Matlab matrix

OUTPUTS:
  no outputs
  








  





  
  
