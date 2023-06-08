clc; clearvars; close all;
% Open a UI to select the directory containing the .fig files
dirpath = uigetdir();

% Get a list of all .fig files in the directory
figfiles = dir(fullfile(dirpath,'**','*.fig'));

% Count how many sets of four files there are
numFigs = length(figfiles)/4;

% Initialize an empty object array to store the data for each figure
data = struct('filename',{},'ydata',{},'prefInc',{},'rPMLC',{},'rLCDr',{});

% Loop over each .fig file in the directory
for iter1 = 1:length(figfiles)
    % Check if the current file is the 3rd of a set of 4 files
    if mod(iter1, 4) == 3
        % Print a status message
        fprintf('Processing File %d of %d',(iter1+1)/4,numFigs)

        % Get the filename of the current .fig file
        filename = figfiles(iter1).name;

        % Define the regular expressions to extract the numbers
        regex1 = '_(\d+)pctPrefInc_';
        regex2 = 'rPMLC(\d+\.?\d*[eE]?[\+\-]?\d*)_(\w+)';
        regex3 = 'rLCDr(\d+\.?\d*[eE]?[\+\-]?\d*).(\w+)';
        
        % Use regular expressions to extract the numbers from the filename
        % The PMLC and LCDr Coefficients seem to be labeled backwards [??fix later]
        prefIncCoeff = regexp(filename, regex1, 'tokens');
        prefIncCoeff = str2double(prefIncCoeff{1}{1});
        
        rPMLCCoeff = regexp(filename, regex3, 'tokens');
        rPMLCCoeff = str2double(rPMLCCoeff{1}{1});
        
        rLCDrCoeff = regexp(filename, regex2, 'tokens');
        rLCDrCoeff = str2double(rLCDrCoeff{1}{1});
        
        % Open the .fig file
        openfig(fullfile(figfiles(iter1).folder, filename),'invisible');
        
        % Get the handle to the current figure
        fig = gcf;
        
        % Get the x and y data from the first subplot in the figure
        ax = fig.Children(1);
        ydata = ax.Children.YData;
        
        % Create a new object to store the filename and y data
        obj = struct('filename',filename,'ydata',ydata,'prefInc',prefIncCoeff,'rPMLC',rPMLCCoeff,'rLCDr',rLCDrCoeff);

        % Store the new object in the data array
        data(end+1) = obj;
        
        % Clean up
        close(fig); clc;
    end
end