clc;
close all;
clear all;

% Step 1: Read the original image
img = 'Hard.jpg';      % char
og_img = imread(img);


% Step 2: Find a intensity threshold to filter out the background
figure(1);
% set a threshold for filterng background
threshold_img = og_img;
subplot(2, 1 ,1);
imshow(threshold_img, [0 255]), title('Img BEFORE thresholding');
colormap gray;
colorbar;
threshold_img(threshold_img > 110 ) = 255;


subplot(2, 1, 2);
imshow(threshold_img, [0 255]), title('Img AFTER thresholding');
colormap gray;
colorbar;

% STEP 3: convert to double data type (binary img)
binary_img = threshold_img(:,:,1);

    % Gradient of a threshold image
% [fx, fy] = imgradientxy(binary_img, 'sobel');
% figure;
% subplot(2, 1, 1);
% imshow(fx);
% title('Gradient in x- direction');
% 
% subplot(2, 1, 2);
% imshow(fy);
% title('Gradient in y-direction')
 

% STEP 4: Calculate the gradient of the binary image after thresholding 
[Gmag, Gdir] = imgradient(binary_img, 'sobel');
figure; 
imshow(Gmag), title('Gradient in both directions - Easy');


% STEP 5: Prep for -> Edge-detecting given the binary image of gradient in both directions 

[~, threshold] = edge(Gmag, 'sobel');
fudge_factor = 0.2;
edges = edge(Gmag, 'sobel', threshold * fudge_factor);

% STEP 5.1 : Dilate the iamge by classfication shapes
    % clssify type of edges
se90 = strel('line', 2, 90);
se0 = strel('line', 2, 0);
elongate_edges = imdilate(edges, [se0 se90]);

dilate_img = elongate_edges;    %fill in lines
figure;
    % remove dots - background noises
min_size = 400;
dilate_img = bwareaopen(dilate_img, min_size);
imshow(dilate_img), title('Dilated image ');

    % fill gaps 
fill_gap_img = imfill(dilate_img, 'holes');
figure;
subplot(2, 2, [1 2]),imshow(fill_gap_img), title('Filling in objects of Binary Image ')

% converting to Double 
fill_gap_img = double(fill_gap_img);
smooth_img = imgaussfilt(fill_gap_img, 1);
subplot(2,2,[3 4]), imshow(smooth_img), title('Smoothed out bianry image')


% STEP 6 - Edge -dtetecting
edge_obj = edge(smooth_img, 'canny');
figure, imshow(edge_obj), title('Edge detected Image');


%STEP 7 -  detect semicrircles
[centres, radius, metrics ] = imfindcircles(edge_obj, [40 60], 'ObjectPolarity', 'dark', 'Sensitivity', 0.9);
% Display the detected circles on the original image
figure;
imshow(dilate_img);
title('Detected Circles');
hold on;
viscircles(centres, radius, 'EdgeColor', 'r');

% STEP 8 - replacing semicircles with circles 
combined_img  = fill_gap_img;

    for i = 1:length(radius)
        % Get circle properties
        center = centres(i, :);
        rad = radius(i);

        % Create a binary mask for the circle
        [x, y] = meshgrid(1:4032, 1:3024);
        circle_mask = ((x - center(1)).^2 + (y - center(2)).^2) <= rad.^2;
       
        % Replace the semicircle with the full circle in the combined image
        combined_img(circle_mask) = 1; % Set pixels inside the circle mask to 1 (white)
    end


    % STEP 9 - Smoothing out combined image 
    se = strel('disk', 13); % Structuring element for smoothing
    smoothed_combined_image = imerode(imdilate(combined_img, se), se);
    [B, ~] = bwboundaries(smoothed_combined_image, 'noholes');
    
    % Display the results
    figure;
    subplot(2, 2, 1), imshow(fill_gap_img), title('Cleaned Binary Image');
    subplot(2, 2, 2), imshow(combined_img), title('Combined Image with Detected Circles');
    subplot(2, 2, [3 4]), imshow(smoothed_combined_image), title('Smoothed Combined Image');
    hold on;
    
    % Overlay boundaries on the smoothed combined image
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); % Overlay boundary
    end
    hold off;
    
    % Step 5: Output the number of final detected objects
    disp(['Number of objects in the final combined image: ', num2str(length(B))]);
