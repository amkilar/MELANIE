%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%	M E L A N I E	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Matlab script for quantitative DNA methylation analysis in fibres
%
% Vesrsion: v2.0
%
% Author: Agata Kilar
%
% Institution: Jiří Fajkus Research Group 
% Mendel Centre for Plant Genomics and Proteomics 
% Central European Institute of Technology
% Masaryk University
%
%
%
% ImageSegmenter toolbox is required.
% 
% Results are saved in results.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Clean session
% It is recommended to run this part before each analysis

close all
clear all
clc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% VARIABLES DEFINITION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Insert name of the image with telomere
I = [''];

%%% Insert name of the image with methylation
I2 = [''];

%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define mask's properties
% The widith of the mask
max_val = 3;

%%% Define the treshold for valid signal 
% [%] of the pixel with the highest intensity within the telomere mask
perc = 0.2;

%%% Decide whether to use image enhancing
% Switch off (0) / on (1) enheacing of the image with methylation channel
% Switched off by default
enhance_meth = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%	SELF-EXECUTABLE PART 1	%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% READ TELOMERE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Read image
[X2, map2] = imread(I);

%%% Convert to RGB
X2RGB = ind2rgb(X2, map2);

%%% Convert to grey scale
X2G = rgb2gray(X2RGB);

%%%% Scale the image
prop = 1;
hight = 0.2;

%%% Check the size of the image to adjust display of additional information
image_size = size(X2RGB);
text_prop = round(image_size(1)/36);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ADJUST IMAGE ORIENTATION %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine the orientation (horizontal or vertical)
if (image_size(1) > image_size(2))
    vertical = 1;
elseif (image_size(2) / image_size(1) >= 7 )
    vertical = 0.5;
else
    vertical = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%% RUN MATLAB TOOLBOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Run imageSegmenter
% To manually extract desired telomere shape

imageSegmenter(X2G)

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%	SELF-EXECUTABLE PART 2	%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ANALYSIS OF METHYLATION COVERAGE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% FILTER OUT TELOMERE AS OBJECT %%%%%%%%%%%%%%%%%%%%%%%

%%% Filter out telomere
% Filter out the biggest object (first sort based on sie and then keep the largest object)
object = bwareafilt(BW, 1);

%%% Make skeleton
% Draw skeleton inside chosen object; filter out "branches" shorter than 10
BW2a = bwskel(object, 'MinBranchLength', 10);

%%% Remove possible loops inside the skeletone
% Create object - disc, and then use basic operation - close
SC = strel('disk', 3, 4);
methyl_fill = imclose(BW2a, SC);
BW2 = bwskel(methyl_fill);

%%%%%%%%%%%%%%%%%%% MEASURE LENGHT OF TELOMERE %%%%%%%%%%%%%%%%%%%%%%%%%

%%% Measure properties of image regions; extract area of telomere
telomere_area = cell2mat(struct2cell(regionprops(BW2, 'Area')));


%%%%%%%%%%%%%%%%%%% CREATE MASK A TELOMERE MASK %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% EXTENDED BY NOISE  %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mask is created based on the skeletone. To every skeleton's pixel there is
% added a fixed number of pixels, on the left and right side of the pixel
% (or up and down in case of horizontal orientation). 
% At the same time mask for noise estimation is created in the same manner,
% but 4 times bigger than declared siez of the telomere mask

%%% Define masks
mask = zeros(size(BW2));
mask_noise = zeros(size(BW2));

%%% Declare range for noise calculation
noise_val = 4 * max_val;

%%% Draw masks
% Drawing depends on telomere orientation (vertical or horizontal)

if (vertical == 1)

for x = (1:size(mask,1))
  
     for y = (1:size(mask,2))
                
        pix = BW2(x, y);
        
            if (pix == 1)
            
                white_row = [y-max_val:y+max_val];
                mask(x, white_row) = 1;
                
                white_row_noise = [y-noise_val:y+noise_val];
                mask_noise(x, white_row_noise) = 1;
            
            end     
    end   
end 

else
    
    for y = (1:size(mask,2))
  
    
    for x = (1:size(mask,1))
                
        pix = BW2(x, y);
        
            if (pix == 1)
            
                white_col = [x-max_val:x+max_val];
                mask(white_col, y) = 1;
                
                white_col_noise = [x-noise_val:x+noise_val];
                mask_noise(white_col_noise, y) = 1;
            end       
    end   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% CALCULATE METHYLATION COVERAGE %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% READ METHYLATION IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read methylation image
[X1, map1] = imread(I2);

%%% Convert to RGB
X1RGB = ind2rgb(X1, map1);

%%% Enhence RGB
if (enhance_meth == 1)
    X1RGB = imadjust(X1RGB, [0.1 0.1 0.1; .4 .4 .4], []);
end

%%% Convert to grey scale
X1G = rgb2gray(X1RGB);

%%%%%%%%%%%%%%%% OVERLAP MASK WITH METHYLATION IMAGE %%%%%%%%%%%%%%%%%%%

% Iterate on telomere mask, and when inside the mask take the value from
% the methylation image

% Create empty matrix
mask_area = zeros(size(X1G));

for i = (1:image_size(1))
   for j = (1:image_size(2))
       
       pixel = mask(i, j);
            
       if (pixel == 1)
          new_pixel = X1G(i, j);
       else 
           new_pixel = 0;
       end
          
       mask_area(i, j) = new_pixel;
           
   end
end

%%%%%%%%%%%%%%%% DEFINE VALID SIGNAL OF THE METHYLATION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   INSIDE THE TELOMERE MASK   %%%%%%%%%%%%%%%%%%%%%%

%%% Detect the pixel with the highest intensity within mask area
max_pixel = max(max(mask_area));

%%% Define valid singal - treshold related to the value of the pixel with 
% the highest intensity within the telomere mask area
tresh = perc * max_pixel;

%%% Treshold methylation image within the mask area
X1B = im2bw(mask_area, tresh);

%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE THE COVARAGE %%%%%%%%%%%%%%%%%%%%%%%%

% Sum up pixels with valid signal withtin mask area
valid_signal = sum(cell2mat(struct2cell(regionprops(X1B, 'area'))));

% Sum up the total area of the mask
mask_area = cell2mat(struct2cell(regionprops(mask, 'area')));

% Calculate the coverage: ratio of valid signal within mask to the area of
% mask (in percentage, rounded up to decimal)
cov = round(((valid_signal / mask_area) * 100), 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATE THE NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Substract telomere mask from noise mask (exteneded telomere mask)
noise_area = mask_noise - mask;

%%%%%%%%%%%%%%%%% OVERLAP MASK WITH METHYLATION IMAGE %%%%%%%%%%%%%%%%%%

% Iterate on noise area, and when inside the mask take the value from
% the methylation image

% Create empty matrix
outside_mask = zeros(size(X1G));

for i = (1:image_size(1))
   for j = (1:image_size(2))
       
       pixel = noise_area(i, j);
            
       if (pixel == 1)
          new_pixel = X1G(i, j);
       else 
           new_pixel = 0;
       end
          
       outside_mask(i, j) = new_pixel;
           
   end
end


%%%%%%%%%%%%%%% DEFINE VALID SIGNAL OF THE METHYLATION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   OUTSIDE THE MASK (background)   %%%%%%%%%%%%%%%%%%%

%%% Treshold methylation image within the mask area
X1B_noise = im2bw(outside_mask, tresh);

%%%%%%%%%%%%%%%%%%%%% CALCULATE THE NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum up pixels with valid signal withtin noise area
valid_signal_noise = sum(cell2mat(struct2cell(regionprops(X1B_noise, 'area'))));

% Sum up the area of the noise area
mask_area_noise = cell2mat(struct2cell(regionprops(noise_area, 'area')));

% Calculate the noise: ratio of valid signal within noise area to the area of
% mask (in percentage, rounded up to decimal)
noise = round(((valid_signal_noise / mask_area_noise) * 100), 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PLOT STAGES OF THE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------- TELOMERE LENGTH -------------------------------

%%% Prepare plot with the biggest object and its skeleton overlapped on the original image 

% The biggest identified object overlapped with the original image
telomere_length_plot_1 = imfuse(X2RGB, object, 'blend');

% Skeleton of the biggest identified object overlapped with the original image
telomere_length_plot_2 = imfuse(X2RGB, BW2, 'blend');

% Overlap the biggest identified object with its skeleton 
telomere_length_plot = imfuse(telomere_length_plot_1, telomere_length_plot_2, 'blend');

% Display calculated length of the telomere on the image
%telomere_length_plot = insertText(telomere_length_plot, [ (0.05*image_size(2)) (image_size(1) - (hight * image_size(1)) ) ], ['Telomere length: ', num2str(telomere_area), ' px'], 'BoxColor', 'white', 'FontSize', text_prop*prop );

% Save the plot in the file
% filename = ['img_', num2str(m), '_tel_length.png'];
% imwrite(telomere_length_plot, filename);


%----------------------- COVERAGE PLOT ---------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%% TELOMERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mask overlapped with original image
cov_plot_mask = imfuse(X1RGB, mask, 'blend');

% Signal detected within the mask area 
cov_plot_sig = imfuse(mask, X1B, 'blend');

% Signal detected within the mask area
cov_plot_comb = imfuse(cov_plot_mask, cov_plot_sig , 'blend');

% Mask with valid signal overlapped with original telomere image
cov_plot_mask_telomere = imfuse(X2RGB, cov_plot_sig, 'blend');

%%%%%%%%%%%%%%%%%%%%%%%%%% NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mask overlapped with original image
cov_plot_mask_noise = imfuse(X1RGB, noise_area, 'blend');

% Signal detected within the mask area 
cov_plot_sig_noise = imfuse(noise_area, X1B_noise, 'blend');

% Signal detected within the noise area
cov_plot_comb_noise = imfuse(cov_plot_mask_noise, cov_plot_sig_noise , 'blend');

% Add information about the noise
cov_plot_comb_noise_text = insertText(cov_plot_comb_noise, [ 0.05*image_size(2) (image_size(1)- hight * image_size(1)) ], ['Noise: ', num2str(noise), '%' ], 'BoxColor', 'white', 'FontSize', text_prop*prop);

%----------------------- RESULT PLOT -----------------------------------

% Mask on the telomere channel
combined_plot_1 = imfuse(X2RGB, X1B, 'blend');
combined_plot_2 = imfuse(X1RGB, mask, 'blend');
combined_plot = imfuse(combined_plot_2, combined_plot_1, 'blend');
combined_plot = insertText(combined_plot, [ 0.05*image_size(2) (image_size(1)- hight * image_size(1)) ], ['Methylation: ', num2str(cov), '%', newline, 'Telomere length: ', num2str(telomere_area), 'px' ], 'BoxColor', 'white', 'FontSize', text_prop*prop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PLOT STAGES OF THE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
if (vertical == 1)
    
    subplot(2, 4, 1)
    imshow(X2RGB)
        
    subplot(2, 4, 2)
    imshow(BW)
        
    subplot(2, 4, 3)
    imshow(object)
        
    subplot(2, 4, 4)
    imshow(telomere_length_plot)
        
    subplot(2, 4, 5)
    imshow(mask)
        
    subplot(2, 4, 6)
    imshow(cov_plot_mask_telomere)
    
    subplot(2, 4, 7)
    imshow(cov_plot_comb)
        
    subplot(2, 4, 8)
    imshow(combined_plot)

elseif (vertical == 0)
    
    subplot(3, 3, 1)
    imshow(X2RGB)
    
    subplot(3, 3, 2)
    imshow(BW)
    
    subplot(3, 3, 3)
    imshow(object)
    
    subplot(3, 3, 4)
    imshow(telomere_length_plot)
    
    subplot(3, 3, 5)
    imshow(mask)
    
    subplot(3, 3, 6)
    imshow(cov_plot_mask_telomere)
        
    subplot(3, 3, 7)
    imshow(cov_plot_comb)
    
    subplot(3, 3, 8)
    imshow(combined_plot)
    
else
    subplot(4, 2, 1)
    imshow(X2RGB)
    
    subplot(4, 2, 2)
    imshow(BW)
    
    subplot(4, 2, 3)
    imshow(object)
    
    subplot(4, 2, 4)
    imshow(telomere_length_plot)
    
    subplot(4, 2, 5)
    imshow(mask)
    
    subplot(4, 2, 6)
    imshow(cov_plot_mask_telomere)
        
    subplot(4, 2, 7)
    imshow(cov_plot_comb)
    
    subplot(4, 2, 8)
    imshow(combined_plot)
    
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% SAVE CALCULATIONS IN TXT FILE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create results table
T = readtable('results.txt');
row = {I, telomere_area, cov, noise};

%%% Write calculations to the file
T2 = [T;row];
writetable(T2,'results.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SAVE RESULTING PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filename = ['img_', num2str(m), '_1.png'];
    imwrite(BW, filename)
    
    filename = ['img_', num2str(m), '_2.png'];
    imwrite(object, filename)
    
    filename = ['img_', num2str(m), '_3.png'];
    imwrite(telomere_length_plot, filename)
    
    filename = ['img_', num2str(m), '_4.png'];
    imwrite(mask, filename)
    
    filename = ['img_', num2str(m), '_5.png'];
    imwrite(cov_plot_sig, filename)
      
    filename = ['img_', num2str(m), '_6.png'];
    imwrite(cov_plot_comb, filename)
    
    filename = ['img_', num2str(m), '_7.png'];
    imwrite(cov_plot_comb_noise_text, filename);
    
    filename = ['img_', num2str(m), '_8.png'];
    imwrite(combined_plot, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%