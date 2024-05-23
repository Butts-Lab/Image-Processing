function CompareFeaturePlotSections
% CompareFeaturePlotSections.m
% Alex Sansom
% 16 May 2024

% Notes:
% Must remove Cell Type or Clusters UMAP folder to run functions that do
% not use it. They will often turn up as results because every pixel is
% read as being expressed.
% When creating figures from a Compare... function with an input UMAP, the
% function filters out the "most similar" UMAP to filter out the comparison
% with the input UMAP.
% To crop: drag to expand rectangle to desired size / area, then double tap
% within the rectangle.
% Optimal crop dimensions are: [190 130 1090 490].

% Call Compare... function with input UMAP jpeg image file, input Cell
% Type / cluster #, or input range of UMAPs.
    % CompareSectionScore("Nhlh2_Exp.jpeg")
    % CompareSectionCoexpress("Mki67_Exp.jpeg")
    % CompareSectionAntiexpress("WPRE_Exp.jpeg")
    % CompareSectionNetCoexpress("Sox2_Exp.jpeg")
    % CompareCellType("p")
    % CompareClusters(1)
    % CompareAll(UMAPsort(26:35))

end

function CompareAll(imgs)
% CompareAll: Compare all possible pairs of UMAPs in the input array.
% Input:
% imgs = struct array of UMAPs.
% Output:
% 2 figures with the top 8 most similar and least similar pairs. 

% Set the crop dimensions to optimal, or manually crop using an example
% UMAP.
    % img = imread("Atoh1_Exp.jpeg");
    % [img, rect] = imcrop(img);
    rect = [190 130 1090 490];

% Initialize numpairs to represent the total combinations and preallocate
% vectors to hold scores and indices of pairs.
    numimgs = length(imgs);
    numpairs = nchoosek(numimgs,2);
    scores = zeros(1,numpairs);
    indices1 = zeros(1,numpairs);
    indices2 = zeros(1,numpairs);
    count = 1;

% Record the score for every possible pair.
    for i = 1:numimgs
        img1 = imread(imgs(i).name);
        img1 = imcrop(img1,rect);
        for j = i+1:numimgs
            indices1(count) = i;
            indices2(count) = j;
            img2 = imread(imgs(j).name);
            img2 = imcrop(img2,rect);
            netnumco = UMAPnumco(img1,img2) - UMAPnumanti(img1,img2);
            % netnumco = UMAPimscore(img1,img2);
            scores(count) = netnumco;
            count = count + 1;
        end
    end

% Sort pairs in descending order.
    [scoressorted, I] = sort(scores);
    indices1sorted = indices1(I);
    indices2sorted = indices2(I);

% Display the eight least similar pairs.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(4,4,'TileSpacing','tight')
    for i = 1:4
        nexttile
        imshow(imgs(indices1sorted(i)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
        titletext = ['Net # Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);
    end
    for i = 1:4
        nexttile
        imshow(imgs(indices2sorted(i)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    end
    for i = 5:8
        nexttile
        imshow(imgs(indices1sorted(i)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
        titletext = ['Net # Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);
    end
    for i = 5:8
        nexttile
        imshow(imgs(indices2sorted(i)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    end

% Display the eight most similar pairs.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(4,4,'TileSpacing','tight')
    for i = 1:4
        nexttile
        imshow(imgs(indices1sorted(numpairs - i + 1)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
        titletext = ['Net # Pixels Coexpressing: ' num2str(scoressorted(numpairs - i))];
        title(titletext);
    end
    for i = 1:4
        nexttile
        imshow(imgs(indices2sorted(numpairs - i + 1)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    end
    for i = 5:8
        nexttile
        imshow(imgs(indices1sorted(numpairs - i + 1)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
        titletext = ['Net # Pixels Coexpressing: ' num2str(scoressorted(numpairs - i))];
        title(titletext);
    end
    for i = 5:8
        nexttile
        imshow(imgs(indices2sorted(numpairs - i + 1)).name);
        rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    end

end

function imgsnew = UMAPsort(range)
% CompareAll: Compare all possible pairs of UMAPs in the input array.
% Input:
% range = range of sorted UMAPs desired. Ex. 26:35 = 11th to 20th most
% expressing in a set of 45 genes.
% Output:
% imgsnew = struct array of desired range.

% Pull UMAPs from current folder and preallocate arrays for sorting.
    imgs = dir("*.jpeg");
    numimgs = length(imgs);
    scores = zeros(1,numimgs);
    indices = zeros(1,numimgs);

% Record the expression of each UMAP.
    for i = 1:numimgs

        img = imread(imgs(i).name);
        BW = (img(:,:,1) == img(:,:,2)) & (img(:,:,2) == img(:,:,3));
        scores(i) = sum(~BW(:));
        indices(i) = i;

    end

% Sort UMAPs and return imgsnew.
    [scoressorted, I] = sort(scores);
    imgssorted = imgs(I);
    imgsnew = imgssorted(range);

end

function CompareClusters(n)
% CompareClusters: Compare selected section of "Clusters" UMAP
% jpeg image to all other UMAP jpeg images in the current folder. Assumes
% all jpegs in the current folder are UMAPs from the relevant dataset
% and assumes all are of the same dimensions.
% Input:
% n = integer respresenting cluster number.
% Output:
% Figure with the top 3 most similar UMAPs.

% Pull input file from current folder and manually / optimally crop. 
% Save rect to crop other UMAPs.
    origimgIN = imread("UMAP_Clusters_Exp.jpeg");
    origimgIN = createMaskclusters(origimgIN,n);
    [imgIN, rect] = imcrop(origimgIN);
    % [imgIN, rect] = imcrop(origimgIN,[190 130 1090 490]);

% Pull all UMAPs from the current folder and assign # of UMAPs to numimgs.
    imgs = dir("*.jpeg");
    numimgs = length(imgs);

% Initialize pair of matrices to save scores and indices for sorting.
    scores = zeros(1,numimgs);
    indices = zeros(1,numimgs);

% Iterate through all UMAPs and assign net # coexpressing.
    for i = 1:numimgs

        img = imread(imgs(i).name);
        img = imcrop(img,rect);
        numco = UMAPnumco(imgIN,img) - (UMAPnumanti(imgIN,img));
        indices(i) = i;
        scores(i) = numco;

    end

    [scoressorted, I] = sort(scores);
    indicessorted = indices(I);

% Display input UMAP, selected area, and three most similar UMAP by
% # coexpressing. Display # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = flip((numimgs - 3):(numimgs - 1))

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

% Display input UMAP, selected area, and three least similar UMAPs by
% # coexpressing. Display # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = 1:3

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

    impixelinfo();

end

function [maskedRGBImage] = createMaskclusters(RGB,n)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 16-May-2024
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
n = (n / 16);
channel1Min = n;
channel1Max = n + 0.0625;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.000;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.000;
channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
sliderBW = ((I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max)) | ...
    (I(:,:,2) == 0);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 255;

end

function CompareCellType(type)
% CompareCellType: Compare manually selected section of "Cell Types" UMAP
% jpeg image to all other UMAP jpeg images in the current folder. Assumes
% all jpegs in the current folder are UMAPs from the relevant dataset
% and assumes all are of the same dimensions.
% Input:
% type = "p" or "i", cell type color of interest.
% Output:
% 2 figures with the top 3 most similar and least similar UMAPs.

% Pull Cell Type UMAP from current folder and manually crop. Mask to hide
% other cell type. Save rect to crop other UMAPs.
    origimgIN = imread("UMAP_Celltype_Exp.jpeg");

    if type == "p"

        BWred = (origimgIN(:,:,1) > origimgIN(:,:,3)) | ...
            (origimgIN(:,:,1) > origimgIN(:,:,2)) | ... 
            ((origimgIN(:,:,1) == origimgIN(:,:,2))) & ... 
            (origimgIN(:,:,1) == origimgIN(:,:,3));
        origimgIN(repmat(~BWred,[1 1 3])) = 255;

    end

    if type == "i"

        BWblue = ((origimgIN(:,:,1) < origimgIN(:,:,3)) & ...
            (origimgIN(:,:,1) < origimgIN(:,:,2))) | ...
            ((origimgIN(:,:,1) == origimgIN(:,:,2))) & ... 
            (origimgIN(:,:,1) == origimgIN(:,:,3));
        origimgIN(repmat(~BWblue,[1 1 3])) = 255;

    end

    [imgIN, rect] = imcrop(origimgIN);
    % [imgIN, rect] = imcrop(origimgIN,[190 130 1090 490]);


% Pull all UMAPs from the current folder and assign # of UMAPs to numimgs.
    imgs = dir("*.jpeg");
    numimgs = length(imgs);

% Initialize pair of matrices to save scores and indices for sorting.
    scores = zeros(1,numimgs);
    indices = zeros(1,numimgs);

% Iterate through all UMAPs and assign # coexpressing using imscore.
    for i = 1:numimgs

        img = imread(imgs(i).name);
        img = imcrop(img,rect);
        numco = UMAPnumco(imgIN,img) - 0.15*UMAPnumanti(imgIN,img);
        indices(i) = i;
        scores(i) = numco;

    end

    [scoressorted, I] = sort(scores);
    indicessorted = indices(I);

% Display input UMAP, selected area, and three most similar UMAPs by
% # coexpressing. Display # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = flip((numimgs - 3):(numimgs - 1))
    % for i = 1:3

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

% Display input UMAP, selected area, and three least similar UMAPs by
% # coexpressing. Display # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = 1:3

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

end

function CompareSectionCoexpress(jpeg)
% CompareSectionCoexpress: Compare manually selected section of input UMAP
% jpeg image to all other UMAP jpeg images in the current folder. Assumes
% all jpegs in the current folder are UMAPs from the relevant dataset
% and assumes all are of the same dimensions.
% Input:
% file = jpeg image file of UMAP of gene of interest.
% Output:
% 2 figures with the top 3 most similar and least similar UMAPs.

% Pull input file from current folder and manually crop. Save rect to crop
% other UMAPs.
    origimgIN = imread(jpeg);
    [imgIN, rect] = imcrop(origimgIN);
    % [imgIN, rect] = imcrop(origimgIN,[190 130 1090 490]);


% Pull all UMAPs from the current folder and assign # of UMAPs to numimgs.
    imgs = dir("*.jpeg");
    numimgs = length(imgs);

% Initialize pair of matrices to save scores and indices for sorting.
    scores = zeros(1,numimgs);
    indices = zeros(1,numimgs);

% Iterate through all UMAPs and assign # coexpressing using UMAPnumco.
    for i = 1:numimgs

        img = imread(imgs(i).name);
        img = imcrop(img,rect);
        numco = UMAPnumco(imgIN,img);
        indices(i) = i;
        scores(i) = numco;

    end

    [scoressorted, I] = sort(scores);
    indicessorted = indices(I);

% Display input UMAP, selected area, and three most similar UMAPs by
% # coexpressing. Display # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = flip((numimgs - 3):(numimgs - 1))

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

% Display input UMAP, selected area, and three least similar UMAPs by
% # coexpressing. Display # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = 1:3

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

end

function [num] = UMAPnumco(umap1,umap2)
% UMAPnumco: Compare two input UMAP images for # of pixels coexpressing.
% Input:
% umap1, umap2 = uint8 matrices representing UMAPs.
% Output:
% score = integer representing calculated # of pixels coexpressing.
    
% Create BW mask of nonexpressing pixels.
    BWumap1 = ((umap1(:,:,1) == umap1(:,:,2)) & ...
        (umap1(:,:,1) == umap1(:,:,3)));
    BWumap2 = ((umap2(:,:,1) == umap2(:,:,2)) & ...
        (umap2(:,:,1) == umap2(:,:,3)));

% Calculate # of pixels coexpressing.
    coexpress = ~BWumap1 & ~BWumap2;
    num = sum(coexpress(:));

end

function CompareSectionAntiexpress(jpeg)
% CompareSectionAntiexpress: Compare manually selected section of input UMAP
% jpeg image to all other UMAP jpeg images in the current folder. Assumes
% all jpegs in the current folder are UMAPs from the relevant dataset
% and assumes all are of the same dimensions.
% Input:
% file = jpeg image file of UMAP of gene of interest.
% Output:
% 2 figures with the top 3 most similar and least similar UMAPs.

% Pull input file from current folder and manually crop. Save rect to crop
% other UMAPs.
    origimgIN = imread(jpeg);
    [imgIN, rect] = imcrop(origimgIN);
    % [imgIN, rect] = imcrop(origimgIN,[190 130 1090 490]);

% Pull all UMAPs from the current folder and assign # of UMAPs to numimgs.
    imgs = dir("*.jpeg");
    numimgs = length(imgs);

% Initialize pair of matrices to save scores and indices for sorting.
    scores = zeros(1,numimgs);
    indices = zeros(1,numimgs);

% Iterate through all UMAPs and assign # antiexpressing using UMAPnumanti.
    for i = 1:numimgs

        img = imread(imgs(i).name);
        img = imcrop(img,rect);
        numanti = UMAPnumanti(imgIN,img);
        indices(i) = i;
        scores(i) = numanti;

    end

    [scoressorted, I] = sort(scores);
    indicessorted = indices(I);

% Display input UMAP, selected area, and three least similar UMAPs by
% # antiexpressing. Display # antiexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = flip((numimgs - 2):(numimgs))

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Antiexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

% Display input UMAP, selected area, and three most similar UMAPs by
% # antiexpressing. Display # antiexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = 2:4

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['# Pixels Antiexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

end

function [num] = UMAPnumanti(umap1,umap2)
% UMAPnumanti: Compare two input UMAP images for # of pixels antiexpressing.
% Input:
% umap1, umap2 = uint8 matrices representing UMAPs.
% Output:
% score = integer representing calculated # of pixels antiexpressing.
    
% Create BW mask of nonexpressing pixels.
    BWumap1 = ((umap1(:,:,1) == umap1(:,:,2)) & ...
        (umap1(:,:,1) == umap1(:,:,3)));
    BWumap2 = ((umap2(:,:,1) == umap2(:,:,2)) & ...
        (umap2(:,:,1) == umap2(:,:,3)));

% Calculate # of pixels antiexpressing.
    antiexpress = (BWumap1 ~= BWumap2);
    num = sum(antiexpress(:));

end

function CompareSectionNetCoexpress(jpeg)
% CompareSectionNetCoexpress: Compare manually selected section of input UMAP
% jpeg image to all other UMAP jpeg images in the current folder. Assumes
% all jpegs in the current folder are UMAPs from the relevant dataset
% and assumes all are of the same dimensions.
% Input:
% file = jpeg image file of UMAP of gene of interest.
% Output:
% 2 figures with the top 3 most similar and least similar UMAPs.

% Pull input file from current folder and manually crop. Save rect to crop
% other UMAPs.
    origimgIN = imread(jpeg);
    [imgIN, rect] = imcrop(origimgIN);
    % [imgIN, rect] = imcrop(origimgIN,[190 130 1090 490]);

% Pull all UMAPs from the current folder and assign # of UMAPs to numimgs.
    imgs = dir("*.jpeg");
    numimgs = length(imgs);

% Initialize pair of matrices to save scores and indices for sorting.
    scores = zeros(1,numimgs);
    indices = zeros(1,numimgs);

% Iterate through all UMAPs and assign # antiexpressing using UMAPnumco 
% UMAPnumanti.
    for i = 1:numimgs

        img = imread(imgs(i).name);
        img = imcrop(img,rect);
        netnumco = UMAPnumco(imgIN,img) - UMAPnumanti(imgIN,img);
        indices(i) = i;
        scores(i) = netnumco;

    end

    [scoressorted, I] = sort(scores);
    indicessorted = indices(I);

% Display input UMAP, selected area, and three most similar UMAPs by
% net # coexpressing. Display net # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = flip((numimgs - 3):(numimgs - 1))

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['Net # Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

% Display input UMAP, selected area, and three least similar UMAPs by
% net # coexpressing. Display net # coexpressing above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = 1:3

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['Net # Pixels Coexpressing: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

end

function CompareSectionScore(jpeg)
% CompareSectionScore: Compare manually selected section of input UMAP
% jpeg image to all other UMAP jpeg images in the current folder. Assumes
% all jpegs in the current folder are UMAPs from the relevant dataset
% and assumes all are of the same dimensions.
% Input:
% file = jpeg image file of UMAP of gene of interest.
% Output:
% 2 figures with the top 3 most similar and least similar UMAPs. 

% Pull input file from current folder and manually crop. Save rect to crop
% other UMAPs.
    origimgIN = imread(jpeg);
    [imgIN, rect] = imcrop(origimgIN);
    % [imgIN, rect] = imcrop(origimgIN,[190 130 1090 490]);

% Pull all UMAPs from the current folder and assign # of UMAPs to numimgs.
    imgs = dir("*.jpeg");
    numimgs = length(imgs);

% Initialize pair of matrices to save scores and indices for sorting.
    scores = zeros(1,numimgs);
    indices = zeros(1,numimgs);

% Iterate through all UMAPs and assign similarity score using imscore.
    for i = 1:numimgs

        img = imread(imgs(i).name);
        img = imcrop(img,rect);
        simscore = UMAPimscore(imgIN,img);
        indices(i) = i;
        scores(i) = simscore;

    end

    [scoressorted, I] = sort(scores);
    indicessorted = indices(I);

% Display input UMAP, selected area, and three most similar UMAPs by
% similarity score. Display similarity scores above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = flip((numimgs - 3):(numimgs - 1))

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['Similarity Score: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

% Display input UMAP, selected area, and three least similar UMAPs by
% similarity score. Display similarity scores above each result UMAP.
    figure()
    hFig2 = gcf;
    hFig2.WindowState = 'maximized';
    tiledlayout(2,3,'TileSpacing','tight')
    nexttile(1,[1 2])
    imshow(origimgIN);
    rectangle('Position',rect,'EdgeColor','r','LineWidth',3,'Curvature',0.1);
    nexttile
    imshow(imgIN);

    for i = 1:3

        nexttile
        imgOUT = imread(imgs(indicessorted(i)).name);
        imshow(imgOUT);
        titletext = ['Similarity Score: ' num2str(scoressorted(i))];
        title(titletext);

    end

    impixelinfo();

end

function [score] = UMAPimscore(umap1,umap2)
% UMAPimscore: Compare two input UMAP images for similarity.
% Input:
% umap1, umap2 = uint8 matrices representing UMAPs.
% Output:
% score = integer representing calculated similarity score.
    
% Mask all points not representing expression in both UMAPs.
    BWumap1 = ((umap1(:,:,1) == umap1(:,:,2)) & ...
        (umap1(:,:,1) == umap1(:,:,3)));
    umap1(repmat(BWumap1,[1 1 3])) = 0;
    BWumap2 = ((umap2(:,:,1) == umap2(:,:,2)) & ...
        (umap2(:,:,1) == umap2(:,:,3)));
    umap2(repmat(BWumap2,[1 1 3])) = 0;

% Calculate similarity scores.
    x = int16(umap1) - int16(umap2);
    y = abs(x);
    z = sum(sum(sum(y)));
    BW = ~(BWumap1 & BWumap2);
    total = sum(BW(:)) * 255 * 3;
    score = 1 - (z / total);

end
