% === User Inputs ===
folderPath = './Floes_bnds/figs/';   % Set the path to your image folder
frameRate = 10;                     % Frames per second for animation
saveVideo = true;                   % Set to true to save video
outputVideoName = '/dat1/makov/KobukBonded/animated_output.avi'; % Output MP4 filename

% === Get list of .jpg files ===
imageFiles = dir(fullfile(folderPath, '*.jpg'));

% Sort images based on numeric filenames
fileNames = {imageFiles.name};
numericOrder = cellfun(@(x) sscanf(x, '%d.jpg'), fileNames);
[~, sortIdx] = sort(numericOrder);
sortedFiles = fileNames(sortIdx);

% === Create figure for animation ===
hFig = figure;
axis off;

% Optional: Create MP4 video writer object
if saveVideo
    videoWriter = VideoWriter(outputVideoName, 'Motion JPEG AVI');
    videoWriter.FrameRate = frameRate;
    open(videoWriter);
end

% === Display images in a loop ===
for i = 1:length(sortedFiles)
    imgPath = fullfile(folderPath, sortedFiles{i});
    img = imread(imgPath);
    imshow(img, 'Border', 'tight');
    drawnow;

    % Write to video
    if saveVideo
        frame = getframe(hFig);
        writeVideo(videoWriter, frame);
    end
end

% === Finalize video ===
if saveVideo
    close(videoWriter);
    disp(['Video saved as ', outputVideoName]);
end
