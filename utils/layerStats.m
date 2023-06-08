clc; clearvars; close all;
% Load data matrix
% load('CompDataCoarse.mat');
load('CompDataAll.mat');

% Define the subranges of the data to calculate the average and standard deviation for
startIndexL1 = 1;
endIndexL1 = 19;
startIndexL2 = 21;
endIndexL2 = 24;

% Loop over each object in the data array
for iter1 = 1:length(data)
    % Extract the subrange for L1 for the current object
    subrangeL1 = data(iter1).ydata(startIndexL1:endIndexL1);
    % Calculate the average and standard deviation of L1
    [data(iter1).L1Var,data(iter1).L1Avg] = std(subrangeL1,1);
    % Calculate the average and variance of L1
%     [data(iter1).L1Var,data(iter1).L1Avg] = var(subrangeL1,1);

    % Extract the subrange for L2 for the current object
    subrangeL2 = data(iter1).ydata(startIndexL2:endIndexL2);    
    % Calculate the average and standard deviation of L2
    [data(iter1).L2Var,data(iter1).L2Avg] = std(subrangeL2,1);
    % Calculate the average and variance of L2
%     [data(iter1).L2Var,data(iter1).L2Avg] = var(subrangeL2,1);
end