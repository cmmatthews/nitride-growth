clc; clearvars; close all;
% Load data matrix
% load('CompDataAllStats3LayerWell.mat');
% load('CompDataAllStats4LayerWell.mat');
load('CompDataStats4LayerWell.mat');

% Extract data to plot
prefInc = [data(:).prefInc]'/100;
rPMLC = [data(:).rPMLC]';
rLCDr = [data(:).rLCDr]';

% Extract data for calculations
compL1 = [data(:).L1Avg]';
varL1  = [data(:).L1Var]';
compL2 = [data(:).L2Avg]';
varL2  = [data(:).L2Var]';

% Calculate % change between experimental and modelled layer compositions
compL1Exp = 0.3;
compL2Exp = 0.82;
compDiffL1 = 100*(compL1 - compL1Exp)/compL1Exp;
compDiffL2 = 100*(compL2 - compL2Exp)/compL2Exp;

% Plot
f = figure(1);
set(f,'units','normalized','outerposition',[0.16 0.1 0.68 0.8]);
subplot(2,2,1)
scatter3(rPMLC,rLCDr,prefInc,50,compDiffL1,'filled');
set(gca,'Projection','perspective');
set(gca,'XScale','log','YScale','log');
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('R_{PM-LC}','FontSize',16,'FontWeight','bold');
ylabel('R_{LC-Dr}','FontSize',16,'FontWeight','bold');
zlabel('Pref. Inc. Coeff.','FontSize',16,'FontWeight','bold');
zlim([min(prefInc) max(prefInc)])
cm = colormap(gca,'turbo');
% cm2 = cm(end-150:end,:);
% colormap(gca,cm2)
cb = colorbar(gca);
zlim([0.2 0.45])
% zlim([0.4 0.45])
clim([-25 50])
% ax.Projection = 'perspective';
% colorbar(gca,'Ticks',[]);

subplot(2,2,2)
% set(gca,'FontSize',14,'FontWeight','bold');
scatter3(rPMLC,rLCDr,prefInc,50,varL1,'filled');
set(gca,'Projection','perspective');
set(gca,'XScale','log','YScale','log');
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('R_{PM-LC}','FontSize',16,'FontWeight','bold');
ylabel('R_{LC-Dr}','FontSize',16,'FontWeight','bold');
zlabel('Pref. Inc. Coeff.','FontSize',16,'FontWeight','bold');
zlim([min(prefInc) max(prefInc)])
cm = colormap(gca,'cool');
colorbar(gca);
zlim([0.2 0.45])
% zlim([0.4 0.45])
clim([0 0.3])

subplot(2,2,3)
% set(gca,'FontSize',14,'FontWeight','bold');
scatter3(rPMLC,rLCDr,prefInc,50,compDiffL2,'filled');
set(gca,'Projection','perspective');
set(gca,'XScale','log','YScale','log');
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('R_{PM-LC}','FontSize',16,'FontWeight','bold');
ylabel('R_{LC-Dr}','FontSize',16,'FontWeight','bold');
zlabel('Pref. Inc. Coeff.','FontSize',16,'FontWeight','bold');
zlim([min(prefInc) max(prefInc)])
% colormap(gca,cm2)
cm = colormap(gca,'turbo');
colorbar(gca);
zlim([0.2 0.45])
% zlim([0.4 0.45])
clim([-25 50])

subplot(2,2,4)
scatter3(rPMLC,rLCDr,prefInc,50,varL2,'filled');
set(gca,'Projection','perspective');
% set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'XScale','log','YScale','log');
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('R_{PM-LC}','FontSize',16,'FontWeight','bold');
ylabel('R_{LC-Dr}','FontSize',16,'FontWeight','bold');
zlabel('Pref. Inc. Coeff.','FontSize',16,'FontWeight','bold');
zlim([min(prefInc) max(prefInc)])
cm = colormap(gca,'cool');
colorbar(gca);
zlim([0.2 0.45])
% zlim([0.4 0.45])
clim([0 0.3])