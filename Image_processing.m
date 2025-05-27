clc;
close all;
clear all;

img = 'Hard.jpg';      % char
og_img = imread(img);
figure;
image(og_img), title('Original image - Easy')

% % % % 2
figure;
% set a threshold for filterng background
threshold_img = og_img;
subplot(2, 1 ,1);
imshow(threshold_img, [0 255]);
colormap gray;
colorbar;
threshold_img(threshold_img > 110 ) = 255;


subplot(2, 1, 2);
imshow(threshold_img, [0 255]);
colormap gray;
colorbar;

% convert to double data type 
binary_img = threshold_img(:,:,1);

% Gradient of a threshold image
% [fx, fy] = imgradientxy(binary_img, 'sobel');

% % % % 3
% figure;
% subplot(2, 1, 1);
% imshow(fx);
% title('Gradient in x- direction');
% 
% subplot(2, 1, 2);
% imshow(fy);
% title('Gradient in y-direction')

[Gmag, Gdir] = imgradient(binary_img, 'sobel');
% % % % 4
figure;         
subplot(2, 1, 1);
imshow(Gmag);
title('Gradient in both directions');


% EDGE DETECTION 
% this is basiaclly same as finding gradient mask
[~, threshold] = edge(Gmag, 'sobel');
fudge_factor = 0.2;
edges = edge(Gmag, 'sobel', threshold * fudge_factor);
% subplot(2, 1, 2);
% imshow(edges);
% title('Binary gradient mask');


    % Step 2: Dilate the iamge 
    
    % clssify type of edges 
se90 = strel('line', 2, 90);
se0 = strel('line', 2, 45);
elongate_edges = imdilate(edges, [se0 se90]);

dilate_img = elongate_edges;    %fill in lines
figure;
imshow(dilate_img), title('Dilated image');

    % remove dots
min_size = 400;
dilate_img = bwareaopen(dilate_img, min_size);

    % fill gaps 
fill_gap_img = imfill(dilate_img, 'holes');
    % converting to Double 
    binary_circles = double(fill_gap_img);
    smooth_img = imgaussfilt(binary_circles, 1);

    all_edges = edge(smooth_img, 'canny');
    figure, imshow(all_edges), title('Edge detected Image');
    
        % detect semicrircles
        [centres, radius, metrics ] = imfindcircles(all_edges, [40 60], 'ObjectPolarity', 'dark', 'Sensitivity', 0.9);
        % Display the detected circles on the original image
        figure;
        imshow(dilate_img);
        title('Detected Circles');
        hold on;
        viscircles(centres, radius, 'EdgeColor', 'r');

        % replacing semicircles with circles 
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
        se = strel('disk', 2); % Structuring element for smoothing
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
% % % % % % % 
% % % COUNTING ROUND OBJ (complete circles only) no semis
% % % % % % % 
[labeled_img, num_objects] = bwlabel(smoothed_combined_image);

% Measure object properties
props = regionprops(labeled_img, 'Area', 'Perimeter', 'Eccentricity', 'BoundingBox', 'Circularity');
        
       
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % Counting number of screws

round_washer = [];
large_screws = [];
small_screws = [];

% pre-define the shape for LONG SCREWS 
for i = 1:num_objects
    bounding_box = props(i).BoundingBox;
    area = props(i).Area;
    perimeter = props(i).Perimeter;
    eccentricity = props(i).Eccentricity;
    h_w_ratio = bounding_box(4) / bounding_box(3);
    circularity = (4 * pi * area) / (perimeter^2);
    % disp([area]);

    % Classifications
    if circularity > 0.7 && eccentricity < 0.5
        % disp([circularity, eccentricity]);
        round_washer = [round_washer, i];
    elseif  area > 25000
        
        large_screws = [large_screws, i];
    elseif (area < 17000 && area > 4000) && eccentricity > 0.6

        small_screws = [small_screws, i];
    end
end


% % % % % FIGURE 
figure; 
imshow(fill_gap_img),title('Extreme level'),  hold on;

for i = round_washer
    rectangle('Position', props(i).BoundingBox, 'EdgeColor', 'y', 'LineWidth', 2);
end

for i = large_screws
    rectangle('Position', props(i).BoundingBox, 'EdgeColor', 'g', 'LineWidth', 2);
end

for i = small_screws
    rectangle('Position', props(i).BoundingBox, 'EdgeColor', 'b', 'LineWidth', 2);
end

hold off; 

disp(['Number of washers:', num2str(length(round_washer))]);
disp(['Number of Big screws : ', num2str(length(large_screws))]);
disp(['Number of Small screws: ', num2str(length(small_screws))]);

total_obj_count = length(round_washer) + length(small_screws) + length(large_screws);
disp(['TOTAL number of objects : ', num2str(total_obj_count)]);