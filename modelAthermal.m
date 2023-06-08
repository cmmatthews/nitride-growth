%% Notes
%
%  Last Modified: 5 June 2023 by CMM

%% Startup
clc; clearvars; close all;
import mlreportgen.ppt.*;

%% User Variables
% -----Material of interest-----
global cations
% cations = "In";
% cations = "Ga";
% cations = "Al";
% cations = ["In"; "Ga"];
cations = ["Ga"; "Al"];
% cations = ["In", "Al"];

% -----RK4 parameters-----
h  = 0.01;                          % time step for numeric solver
hP = 2;                             % decimal precision of h (to avoid rounding errors)
t0 = 0;                             % initial time for numeric solver
tf = 24;                            % final time for numeric solver

% -----Growth Variables-----
% ---MME Shuttering---
tOpen        = 14;                  % duration of metal shutters open
tClosed      = 10;                  % duration of metal shutters closed
shutterDelay = 0.1;                 % transition time between open and closed states
% ---Nitrogen Flux---
fluxN        = 1.53;                % #atoms of active nitrogen/cm^2*s incident on surface [x10^15]
% fluxN        = [1.09, 1.61, 1.97];  % nitrogen flux array [!!! Not supported in crystalLayersClbk()]
fluxNProf    = 'fixed';             % profile of nitrogen flux
% ---Metal Flux---
IIIV         = 1.3;                 % III/V flux ratio
% IIIV         = [0.9, 1.5, 1.8];     % III/V flux ratio array
xA           = 0.5;                 % mole fraction of cation A
% xA           = [0.3, 0.4, 0.5]      % mole fraction array
% ---Growth Mode--- ('comp' = compositional | 'pref' = preferential incorporation)
global growthMode
growthMode   = 'pref';
% ---Temperature (currently unused)---
tSubs        = 600;                 % substrate temperature in degrees C
% tSubs        = [300, 400, 500];     % substrate temperature array

% -----Output Settings-----
% ---RK4 output options (cmd, file, both, none) (deprecated)---
outputSettings = 'none';
% ---Flag to save RK4 output---
fileOut = 0;
% fileOut = 1;
% ---PPT output (deprecated)---
% genPPT = 'on';
genPPT = 'off';

%% Physics Variables
% -----Rate Coefficients-----
global preferredIncorpA
preferredIncorpA = 0.30;                  % preferential incorporation coefficient
rPMLC            = 30;                    % LC-PM exchange rate coeff
% rPMLC            = [1, 2, 3, 5, 10];      % LC-PM exchange rate array
rLCDr            = 20;                    % droplet-LC exchange rate coeff
% rLCDr            = [1, 2, 3, 5, 10];      % droplet-LC exchange rate array

%% Initialize
% Define material of interest
global sysType
if length(cations) == 1
    sysType = 'binary';
    catStr = cations;
elseif length(cations) == 2
    sysType = 'ternary';
    catStr = sprintf('%s (Element B: %s)',cations(1),cations(2));
elseif length(cations) == 3
    sysType = 'quaternary';
    catStr = sprintf('%s (Element B: %s, Element C: %s)',cations(1),cations(2),cations(3));
else
    error('Invalid number of cations specified');
end
substrate = buildSubstrate('gan');

% Initialize physical constants
mPM     = substrate.MLCapacity;     % # of metal atoms/cm^2 in a pseudomorphic ML [x10^15]
mLC     = 4*mPM/3;                  % # of metal atoms/cm^2 in a laterally-contracted ML [x10^15]
LDIn    = 2e-6;                     % diffusion length [cm]
LDGa    = 2e-6;                     % diffusion length [cm]
LDAl    = 2e-6;                     % diffusion length [cm]

% Initialize Growth Variables
% MME Shuttering
shutterRes   = shutterDelay/h;      % time resolution for use in genFluxIII()
% Nitrogen Flux
fluxNs       = genFluxN(fluxN,tOpen,tClosed,4,0.5,1.8,fluxNProf);
% Metal Flux
fluxIIIs     = genFluxIII(IIIV,fluxN,tOpen,tClosed,tf,h,shutterDelay,shutterRes,'gradual');
x0           = xA*mPM;              % initial condition to be used if starting with any pre-existing crystal
% Temperature (currently unused)
temps = tSubs+273.15;               % substrate temperature(s) in Kelvin

% Collect above variables into a single object
physConsts = table(h,mPM,mLC,fluxNs,fluxIIIs,LDIn,LDGa,LDAl,xA,temps,rPMLC,rLCDr);
physConsts.Properties.VariableDescriptions = ...
    {'Runge-Kutta Timestep' ...
     'Num. of Atomic Sites in PM-ML' ...
     'Num. of Atomic Sites in LC-ML' ...
     'Active Nitrogen Flux' ...
     'Group-III Cation Flux' ...
     'Indium Adatom Diffusion Length' ...
     'Gallium Adatom Diffusion Length' ...
     'Aluminum Adatom Diffusion Length' ...
     'Cation-A Target Mole Fraction' ...
     'Substrate Temperature' ...
     'PM to LC Exchange Rate Coeff.' ...
     'LC to droplet Exchange Rate Coeff.'};
physConsts.Properties.VariableUnits = ...
     {'s' ...
      '10^{15} atoms cm^{-2}' ...
      '10^{15} atoms cm^{-2}' ...
      '10^{15} atoms cm^{-2} s^{-1}' ...
      '' ...
      'cm' ...
      'cm' ...
      'cm' ...
      'K' ...
      '' ...
      '' ...
      ''};

%% Defining the System
tic
global crystals currentCrystal
ics = 0;
MMEData = zeros(round(tf/h)+1,10);
system = genSystem(sysType,physConsts,ics);
crystals = cell(size(system));
currentCrystal = 1;
subLayers = [];
crystals{currentCrystal} = subLayers;

% Set up file output
if ~exist('output', 'dir')
   mkdir('output');
end
runDate = datestr(now,'yyyy-mm-dd_HH-MM-SS');
outDir = fullfile('output',runDate);
if ~exist(outDir, 'dir')
   mkdir(outDir);
end

% PPT output (deprecated)
outputSettings = 'none';
if(~strcmp(outputSettings,'none'))
    headingStr = sprintf('%6s%10s','i','t');
    labels = ["Crys", "N_PM", "N_LC", "Total"];
    for j = 1:length(labels)
        headingStr = strcat(headingStr,sprintf('%16s ',sprintf('%d:%s',j+1,labels(j))));
    end
    disp(headingStr);
end

%% Solving the System & Data Organization
% Create cell arrays for RK4 data and crystal structures
dataLabels = cell(length(fluxNs),length(fluxIIIs),length(rPMLC),length(rLCDr),length(xA),length(temps));
dataTables = dataLabels;
crystalStructs = dataLabels;

% Loop through all combinations of growth variables
for iter1 = 1:length(fluxNs)
 for iter2 = 1:length(fluxIIIs)
  for iter3 = 1:length(rPMLC)
   for iter4 = 1:length(rLCDr)
    for iter5 = 1:length(xA)
     for iter6 = 1:length(temps)
        % Solve system with RK4
        currentSys = system{iter1,iter2,iter3,iter4,iter5};
        MMEData = RK4(h,hP,t0,tf,currentSys.fns,currentSys.ics,@rk4Callback,outputSettings);
        switch sysType
            case "binary"
                varNames = {'t','VCrys','VPM','VLC','VDroplet','Accumulated','Desorbed'};
            case "ternary"
                varNames = {'t','VCrys','VPM','VLC','VDroplet',...
                    'ACrys','APM','ALC','ADroplet','Accumulated','AAcc'};
            case "quaternary"
                varNames = {'t','VCrys','VPM','VLC','VDroplet',...
                    'ACrys','BCrys','APM','BPM','ALC','BLC',...
                    'ADroplet','BDroplet','Accumulated'};
        end
        tbl = array2table(MMEData,'VariableNames',varNames);
        if strcmpi('const',fluxNProf)
            dataLabels(iter1,iter2,iter3,iter4,iter5) = {sprintf(['Nitrogen Flux = %0.2f, '...
                'III-V = %0.2f, Shutters = %d:%d, r_{PM->LC} = %0.2f, r_{LC->Dr} = %0.2f, '...
                'x = %0.2f, T = %d %cC'],fluxN,IIIV(iter2),tOpen,tClosed,...
                                      rPMLC(iter3),rLCDr(iter4),xA(iter5),round(temps(iter6)-273.15,2),176)};
        else
            dataLabels(iter1,iter2,iter3,iter4,iter5) = {sprintf(['Nitrogen Flux = %s, '...
                'III-V = %0.2f, Shutters = %d:%d, r_{PM->LC} = %0.2f, r_{LC->Dr} = %0.2f, '...
                'x = %0.2f, T = %d %cC'],func2str(fluxNs{iter1}),IIIV(iter2),tOpen,tClosed,...
                                      rPMLC(iter3),rLCDr(iter4),xA(iter5),round(temps(iter6)-273.15,2),176)};
        end
        dataTables(iter1,iter2,iter3,iter4,iter5) = {tbl};
        currentCrystal = currentCrystal + 1;

        if fileOut == 1
            mode = 'MME';
            dateFormat = 'yyyy-mm-dd_HH.MM';
            date = datestr(now,dateFormat);
            outFileName = sprintf('%s_%dpct_%0.2f_%0.1f-%0.1f_%0.2f-%0.2f_%s',mode,xA(iter5)*100,IIIV,...
                                                                              tOpen,tClosed,...
                                                                              rPMLC(iter3),rLCDr(iter4),date);
            outFileName = strrep(outFileName,'.','p');
            outFileName = [outFileName, '.txt'];
            writetable(tbl,outFileName);

            headers = sprintf(['Nitrogen Flux: %0.2f\n',...
                               'III-V Ratio: %0.2f\n',...
                               'Shutter Timing: %0.1f sec open : %0.1f sec closed with %0.2f sec delay\n',...
                               'Target Composition: %0.1f of Element A: %s\n',...
                               'Exchange Rates: %0.2f PM<->LC, %0.2f LC<->PM\n',...
                               'System Timing: %0.1f sec start : %0.1f sec end to intervals of %0.2g sec'],...
                               fluxN,...
                               IIIV(iter2),...
                               tOpen,tClosed,shutterDelay,...
                               xA(iter5)*100,catStr,...
                               rPMLC(iter3),rLCDr(iter4),...
                               t0,tf,h);
            prepend2file(headers,outFileName);

            outFileName = sprintf('%s_%dpct_%0.2f_%0.1f-%0.1f_%0.2f-%0.2f_%s_crystal',mode,xA(iter5)*100,IIIV,...
                                     tOpen,tClosed,rPMLC(iter3),rLCDr(iter4),date);
            outFileName = strrep(outFileName,'.','p');
            outFileName = [outFileName, '.txt'];
            writetable(crystal,outFileName);
        end
     end
    end
   end
  end
 end
end
toc

%% Plotting
% PPT output (deprecated)
if(strcmp(genPPT,'on'))
    ppt = Presentation('output.pptx','outputTemplate.pptx');
end

close all;
% initialize figure counter
counter = 1;
% Loop through all combinations of growth variables
for iter1 = 1:length(fluxNs)
 for iter2 = 1:length(fluxIIIs)
  for iter3 = 1:length(rPMLC)
   for iter4 = 1:length(rLCDr)
    for iter5 = 1:length(xA)
        % get data for this specific growth
        data = dataTables{iter1,iter2,iter3,iter4,iter5};
        label = dataLabels{iter1,iter2,iter3,iter4,iter5};
        crystal = crystals{iter1,iter2,iter3,iter4,iter5};

        % plot bilayer evolution
        [f2,fPM] = plotMLs(data,label,counter);
        % save bilayer figure to .fig and .tif files
        bilayerFigName = sprintf('%04dBilayer_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e',...
                                 counter,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        fprintf('%04dBilayer_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e\n',...
                                 counter,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        outFile = fullfile(outDir,bilayerFigName);
        savefig(f2,outFile);
        saveas(f2,outFile,'tiff');
        close(f2);
        % save PM figure to .fig and .tif files
        bilayerFigName = sprintf('%04dPMML_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e',...
                                 counter+1,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        fprintf('%04dPMML_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e\n',...
                                 counter+1,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        outFile = fullfile(outDir,bilayerFigName);
        savefig(fPM,outFile);
        saveas(fPM,outFile,'tiff');
        close(fPM);
        
        % add to PPT (deprecated)
        if(strcmp(genPPT,'on'))
            saveas(gcf,sprintf('tempImg%d.tif',counter));
            picSlide = add(ppt,'Content Only');
            replace(picSlide,'Content',Picture(sprintf('tempImg%d.tif',counter)));
            close(gcf);
        end
        
        % increment figure counter
        counter = counter + 2;
        
        % plot crystal composition w/o labels
        fB1 = plotCrystal(crystal,mPM,counter,false);
        % save composition figure to .fig and .tif files
        compFigName = sprintf('%04dComp_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e',...
                             counter,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        fprintf('%04dComp_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e\n',...
                                 counter,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        outFile = fullfile(outDir,compFigName);
        savefig(fB1,outFile);
        saveas(fB1,outFile,'tiff');
        close(fB1);
        counter = counter + 1;

        % plot crystal composition w/ labels
        fB2 = plotCrystal(crystal,mPM,counter,true);
        % save composition figure to .fig and .tif files
        compFigName = sprintf('%04dCompLabel_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e',...
                             counter,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        fprintf('%04dCompLabel_%dpct%s_%dpctPrefInc_rPMLC%0.0e_rLCDr%0.0e\n',...
                                 counter,xA(iter5)*100,cations(1),preferredIncorpA*100,rPMLC(iter3),rLCDr(iter4));
        outFile = fullfile(outDir,compFigName);
        savefig(fB2,outFile);
        saveas(fB2,outFile,'tiff');
        close(fB2);

        % increment figure counter
        counter = counter + 1;
    end
   end
  end
 end
end

if(strcmp(genPPT,'on'))
    close(ppt);
    rptview(ppt);
    for iter1 = 1:counter-1
        imgFile = sprintf('tempImg%d.tif',iter1);
        delete(imgFile);
    end
end

%% Misc. End of Run Commands
% Compare number of atoms in from flux and from RK4 output
f3 = fluxIIIs{1};
totalIn1 = integral(f3, t0, tf, 'ArrayValued', true);
totalIn2 = data.VCrys(end) + data.VPM(end) + data.VLC(end) + data.VDroplet(end);
ratio1 = totalIn2/totalIn1;
fprintf('ratio of III/V input to atoms in system = %0.4f \n', ratio1);
% Compare number of cation A in from flux and from RK4 output
fCatA = @(t) xA*f3(t);
catAIn1 = integral(fCatA, t0, tf, 'ArrayValued', true);
catAIn2 = data.ACrys(end) + data.APM(end) + data.ALC(end) + data.ADroplet(end);
ratio2 = catAIn2/catAIn1;
fprintf('ratio of catA input to catA in system = %0.4f \n', ratio2);

%% Main Functions
function t_w_for_plotting = RK4(h,hP,t0,tf,fns,xo,callbackFun,saveOutput)

    % Flag to enable msgbox popup on imaginary values in solution (error)
    myflag = 1;

    % Check inputs and set output flags based on saveOutput value
    if nargin < 7
        saveOutput = 'none';
    end
    if ~(strcmp(saveOutput,'cmd') || strcmp(saveOutput,'file')...
       || strcmp(saveOutput,'both') || strcmp(saveOutput,'none'))
        warning('saveOutput value not recognized. Defaulting to file output.')
        fileOut = 1;
        cmdOut = 0;
    elseif strcmp(saveOutput,'file')
        fileOut = 1;
        cmdOut = 0;
    elseif strcmp(saveOutput,'cmd')
        fileOut = 0;
        cmdOut = 1;
    elseif strcmp(saveOutput,'both')
        fileOut = 1;
        cmdOut = 1;
    else
        fileOut = 0;
        cmdOut = 0;
    end

    %Check that size of functions and initial conditions match
    if size(fns) ~= size(xo)
        error('inconsistent size for function and initial condition params');
    end
    %Check that both inputs above are vectors
    if ~isvector(fns)
        error('function input must be a vector with 1xN or Nx1 dimensions')
    end
    %Force both inputs to be row vectors to match chosen t_w format
    if iscolumn(fns)
        fns = transpose(fns);
        xo = transpose(xo);
    end

    %Initialize step size and number of steps.
    N = floor((tf-t0)/h);  %Total number of steps to be taken
    h = (tf-t0)/N;    %Actual step size so that integration ends at tfinal

    %Open an output file.
    if fileOut
        dateFormat = 'yyyy-mm-dd_HH.MM.SS.FFF';
        date = datestr(now,dateFormat);
        outFileName = sprintf('RK4Output%s.txt',date);
        OutputFile = fopen(outFileName,'w');
    end

    %Print title to the screen and to the output file.
    if cmdOut
        fprintf('\n OUTPUT FROM RK4.m \n\n');
    end
    if fileOut
        fprintf(OutputFile, '\n OUTPUT FROM RK4.m \n\n');
    end

    %Generate strings for initial conditions to write to log file
    if cmdOut || fileOut
        icVarStr = '[';
        icValStr = '[';
        for iter = 1:length(xo)
            icVarStr = strcat(icVarStr(2:end),sprintf('x%d, ',iter));
            icValStr = strcat(icValStr(2:end),sprintf('%10.6f, ',xo(iter)));
        end
        icVarStr = strcat(icVarStr(1:end-2),']');
        icValStr = strcat(icValStr(1:end-2),']');
    end

    %Print information about the method and the problem to the screen and to the output file.
    if cmdOut
        fprintf(' Using RK4 with h = %10.6f to integrate a first-order ODE IVP\n', h);
        fprintf(' from t = %10.6f to t = %10.6f with initial %s values = %s\n\n', t0, tf, icVarStr, icValStr);
    end
    if fileOut
        fprintf(OutputFile, ' Using RK4 with h = %10.6f to integrate a first-order ODE IVP\n', h);
        fprintf(OutputFile, ' from t = %10.6f to t = %10.6f with initial %s values = %s\n\n', t0, tf, icVarStr, icValStr);
    end
    
    %Array for plotting
    t_w_for_plotting = [];

    %Generate format strings for log file headings
    if cmdOut || fileOut
        headingStr = sprintf('%6s%10s','i','t');
        for iter = 1:length(xo)
            headingStr = strcat(headingStr,sprintf('%16s, ',sprintf('w%d',iter)));
        end
        sepStr = repmat('-',1,length(headingStr)+6);
    end
    %Print the column headings and separator for the results table.
    if cmdOut
        fprintf('%s\n', headingStr);
        fprintf('%s\n', sepStr);
    end
    if fileOut
        fprintf(OutputFile, '%s\n', headingStr);  
        fprintf(OutputFile, '%s\n', sepStr);
    end

    %Prepare for the main loop.
    told = t0;
    wold = xo;

    %Populate plotting array with initial conditions
    t_w_for_plotting(1,:) = [told, wold];

    %Print the information for i=0.
    if cmdOut || fileOut
        dataStr = sprintf('%6d%s%10.6f',0,repmat(' ',1,4),told);
        for iter = 1:length(wold)
            dataStr = strcat(dataStr,sprintf('%s%+1.4e',repmat(' ',1,5),wold(iter)));
        end
    end
    if cmdOut
        fprintf('%s\n', dataStr);
    end
    if fileOut
        fprintf(OutputFile, '%s\n', dataStr);
    end

    %Main loop
    for i = 0:N-1
        tnew = round(told + h,hP);

        % iterate over 4 RK increments (k)
        try
            % iterate over all functions in the system
            k = zeros(4,length(fns));
            for iter2 = 1:length(fns)
                k(1,iter2) = fns{iter2}(told, wold);    % each row contains the k_i for each function
            end
            for iter2 = 1:length(fns)
                k(2,iter2) = fns{iter2}(told+(h/2), wold+h*(k(1,:)/2));    % each row contains the k_i for each function
            end
            for iter2 = 1:length(fns)
                k(3,iter2) = fns{iter2}(told+(h/2), wold+h*(k(2,:)/2));    % each row contains the k_i for each function
            end
            for iter2 = 1:length(fns)
                k(4,iter2) = fns{iter2}(told+h, wold+h*k(3,:));    % each row contains the k_i for each function
            end
        catch ME
            fprintf('Time of Error: t = %0.2f\n',tnew);
            rethrow(ME)
        end


        %RK4 Method for calculating wnew
        wnew = wold + (h/6).*(k(1,:) + 2*k(2,:) + 2*k(3,:) + k(4,:));
        wnew = round(wnew,15); % round to 1e-15 precision (one atom precision)
        
        %Callback for external variables
        if isa(callbackFun, 'function_handle')
            callbackFun([told, wold], [tnew, wnew]);
        else
            if i == 0
                disp('RK4 Callback Function not defined')
            end
        end

        %Populate plotting array with told and wold. i+3 because A and I.C.
        %are in columns 1 and 2.
        if (~isreal(tnew) || ~isreal(wnew)) && myflag == 1
            msgfig = msgbox(sprintf('Imaginary value detected @ t = %10.6f',tnew));
            uiwait(msgfig)
            myflag = 0;
        end
        t_w_for_plotting(i+2,:) = [tnew, wnew];
        
        %Print the information for i.
        if cmdOut || fileOut
            dataStr = sprintf('%6d%s%10.6f',i+1,repmat(' ',1,4),tnew);
            for iter1 = 1:length(wnew)
                dataStr = strcat(dataStr,sprintf('%s%+1.4e',repmat(' ',1,5),wnew(iter1)));
            end
        end
        if cmdOut
            fprintf('%s\n', dataStr);
        end
        if fileOut
            fprintf(OutputFile, '%s\n', dataStr);
        end
        
        %Prepare for the next time through the loop.
        told = tnew;
        wold = wnew;
        try
            if sum([isinf(wnew), isnan(wnew)]) ~= 0
                if fileOut
                    fclose(OutputFile);
                end
                error('At least one w value is Inf of NaN at t = %0.2f',tnew);
            end
        catch
            fprintf('At least one w value is Inf of NaN at t = %0.2f',tnew);
        end
    end
    
    %Print another horizontal line & conclusion
    if cmdOut
        fprintf('%s\n', sepStr);
    end
    if fileOut
        fprintf(OutputFile, '%s\n', sepStr);
        %Close the output file.
        fclose(OutputFile);
    end
end

function rk4Callback(old, new)
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% DOES NOT HANDLE BINARY GROWTHS YET
% COMMENT OUT LINE 508 (callbackFun([told, wold], [tnew, wnew]);) TO PREVENT ERROR
% CODE WILL STILL FAIL FOR COMP PROFILE PLOTS BUT OTHERWISE RUNS FINE TO THAT POINT
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    crystalLayersClbk(old,new);
end

function sys = genSystem( type, constTbl, ics )
    if nargin == 2
        if strcmp(type,'binary')
            sys = genBinarySys(constTbl);
        elseif strcmp(type,'ternary')
            sys = genTernarySys(constTbl);
%         elseif strcmp(type,'quaternary')
%             fns = genQuatSys();
        else
            error('Invalid system type specified');
        end
    else
        if strcmp(type,'binary')
            sys = genBinarySys(constTbl,ics);
        elseif strcmp(type,'ternary')
            sys = genTernarySys(constTbl,ics);
%         elseif strcmp(type,'quaternary')
%             fns = genQuatSys();
        else
            error('Invalid system type specified');
        end
    end
end

function sys = genBinarySys( constTbl, icsIn )
    % unpack relevant physical constants from input
    h = constTbl.h;
    mPM = constTbl.mPM;
    mLC = constTbl.mLC;
    fluxNs = constTbl.fluxNs;
    fluxIIIs = constTbl.fluxIIIs;
    temps  = constTbl.temps;
    LD = constTbl.LDGa;
    
    % generate systems of equations for every combination of physical constants
    sys = cell(length(fluxNs),length(fluxIIIs));
    for iter1 = 1:length(fluxNs)
        for iter2 = 1:length(fluxIIIs)
            for iter3 = 1:length(temps)
                fluxN   = fluxNs{iter1};
                fluxIII = fluxIIIs{iter2};
                T = temps(iter3);

                if nargin == 2
                    ics = reshape(padarray(icsIn,5-length(icsIn),'post'),[1,5]);
                else
                    ics = zeros(1,5);
                end

                fns = {... % #1:Num. crystal (d/dt n_crys = 
                       @(t,xs) (growthRate(xs(2),mPM,fluxN(t),h)),...
                       ...
                       ... % #2:Num. PM
                       @(t,xs) (pmBinary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD)),...
                       ...
                       ... % #3:Num. LC
                       @(t,xs) (lcBinary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD)),...
                       ...
                       ... % #4:Num. Droplets
                       @(t,xs) (drBinary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD)),...
                       ...
                       ... % #5:Total Accumulation
                       @(t,xs) (fluxIII(t))...
                      };

                tempSys.ics = ics;
                tempSys.fns = fns;
                sys(iter1,iter2) = {tempSys};
            end
        end
    end
end

function sys = genTernarySys( constTbl, icsIn )
    % growth of A_(x)B_(1-x)N
    % written such that the atomic radius of A is larger than B
    % therefore, ternaries should be specified as InGaN, GaAlN, or InAlN
    
    % unpack relevant physical constants from input
    h      = constTbl.h;
    mPM    = constTbl.mPM;
    mLC    = constTbl.mLC;
    fluxNs = constTbl.fluxNs;
    fluxIIIs  = constTbl.fluxIIIs;
    rPMLC  = constTbl.rPMLC;
    rLCDr  = constTbl.rLCDr;
    xA     = constTbl.xA;
    temps  = constTbl.temps;
    LD     = constTbl.LDGa;
    
    % generate systems of equations for every combination of physical constants
    sys = cell(length(fluxNs),length(fluxIIIs),length(rPMLC),length(rLCDr),length(xA));
    for iter1 = 1:length(fluxNs)
        for iter2 = 1:length(fluxIIIs)
            for iter3 = 1:length(rPMLC)
                for iter4 = 1:length(rLCDr)
                    for iter5 = 1:length(xA)
                        for iter6 = 1:length(temps)
                            fluxN   = fluxNs{iter1};
                            fluxIII    = fluxIIIs{iter2};
                            r1      = rPMLC(iter3);
                            r2      = rLCDr(iter4);
                            x       = xA(iter5);
                            T       = temps(iter6);

                            N_A     = @(t) ((IIIV(t)-1)*fluxN(t));   % rate of metal accumulation on the surface, cm^-2*s^-1
                            fluxA   = @(t) (fluxIII(t)*x);
                            fluxB   = @(t) (fluxIII(t)*(1-x));

                            if nargin == 2
                                ics = reshape(padarray(icsIn,10-length(icsIn),'post'),[1,10]);
                            else
                                ics = zeros(1,10);
                            end

                            fns = {... % #1:Num. crystal
                                   @(t,xs) (growthRate(xs(2), mPM, fluxN(t), h)),...
                                   ...
                                   ... % #2:Num. PM
                                   @(t,xs) (pmBinary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD)),...
                                   ...
                                   ... % #3:Num. LC
                                   @(t,xs) (lcBinary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD)),...
                                   ...
                                   ... % #4:Num. Droplets
                                   @(t,xs) (drBinary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD)),...
                                   ...
                                   ... % #5:A crystal
                                   @(t,xs) (growthRateTernaryNew(xs(2), xs(6), mPM, concA(xs(6),xs(2)-xs(6)), fluxN(t), h)),...
                                   ...
                                   ... % #6:A PM
                                   @(t,xs) (pmTernary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD, x, r1, r2)),...
                                   ...
                                   ... % #7:A LC
                                   @(t,xs) (lcTernary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD, x, r1, r2)),...
                                   ...
                                   ... % #8:A Droplets
                                   @(t,xs) (drTernary(t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD, x, r1, r2)),...
                                   ...
                                   ... % #9:Total Accumulation
                                   @(t,xs) (fluxIII(t)),...
                                   ...
                                   ... % #10:A Accumulation
                                   @(t,xs) (fluxA(t))...
                                  };

                            tempSys.ics = ics;
                            tempSys.fns = fns;
                            sys(iter1,iter2,iter3,iter4,iter5,iter6) = {tempSys};
                        end
                    end
                end
            end
        end
    end
end

function sys = genQuatSys()

end

%% Rate Functions
function ratePM = pmBinary( t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD )
    fN = fluxN(t);
    fIII = fluxIII(t);
    
    %---------------------------------
    %-------------PM Only-------------
    %---------------------------------
%     % growth rate & flux
%     growRate = growthRate(xs(2),mPM,fN,h);
%     fluxIn = fIII;
%     % replenish
%     repFromLC = 0;

    %---------------------------------
    %-------------PM & LC-------------
    %---------------------------------
%     % growth rate & flux
%     growRate = growthRate(xs(2),mPM,fN,h);
%     fluxIn = (1-vd(xs(2),mPM,xs(3),LD,growRate,fIII,h))*fIII;
%     % replenish
%     repFromLC = replenish(xs(3),growRate,h);

    %---------------------------------
    %-------------All MLs-------------
    %---------------------------------
    % growth rate & flux
    growRate = growthRate(xs(2),mPM,fN,h);
    fluxIn = (1-vd(xs(2),mPM,xs(3),LD,growRate,fIII,h))*fIII;
    % replenish
    repFromLC = replenish(xs(3),growRate,h);
    
    
    ratePM = fluxIn + repFromLC - growRate;
end

function rateAPM = pmTernary( t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD, x, rPMLC, rLCDr )
    fN = fluxN(t);
    fIII = fluxIII(t);
    fluxA = fluxIII(t)*x;
    
    %---------------------------------
    %-------------PM Only-------------
    %---------------------------------
%     % growth & flux
%     cA = concA(xs(6), xs(2)-xs(6));
%     growRateA = growthRateTernaryNew(xs(2),xs(6),mPM,cA,fN,h);
%     fluxAIn = fluxA;
%     % replenish
%     repAFromLC = 0;
%     % VCS
%     segLCPM = 0;

    %---------------------------------
    %-------------PM & LC-------------
    %---------------------------------
%     % growth & flux
%     growRate = growthRate(xs(2),mPM,fN,h);
%     cA = concA(xs(6), xs(2)-xs(6));
%     growRateA = growthRateTernaryNew(xs(2),xs(6),mPM,cA,fN,h);
%     fluxAIn = (1-vd(xs(2),mPM,xs(3),LD,growRate,fIII,h))*fluxA;
%     % replenish
%     repFromLC = replenish(xs(3),growRate,h);
%     repAFromLC = repFromLC*concA(xs(7),xs(3)-xs(7));
%     % VCS
%     segLCPM = rPMLC*xs(6)*(xs(3)-xs(7));

    %---------------------------------
    %-------------All MLs-------------
    %---------------------------------
    % growth & flux
    growRate = growthRate(xs(2),mPM,fN,h);
    cA = concA(xs(6), xs(2)-xs(6));
    growRateA = growthRateTernaryNew(xs(2),xs(6),mPM,cA,fN,h);
    fluxAIn = (1-vd(xs(2),mPM,xs(3),LD,growRate,fIII,h))*fluxA;
    % replenish
    repFromLC = replenish(xs(3),growRate,h);
    repAFromLC = repFromLC*concA(xs(7),xs(3)-xs(7));
    % VCS
    segLCPM = rPMLC*xs(6)*(xs(3)-xs(7));
    
    
    rateAPM = fluxAIn - segLCPM + repAFromLC - growRateA;
end

function rateLC = lcBinary( t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD )
    fN = fluxN(t);
    fIII = fluxIII(t);
    
    %---------------------------------
    %-------------PM & LC-------------
    %---------------------------------
%     % growth rate & flux
%     growRate = growthRate(xs(2),mPM,fN,h);
%     fluxIn = vd(xs(2),mPM,xs(3),LD,growRate,fIII,h)*fIII;
%     % replenish
%     repFromLC = replenish(xs(3),growRate,h);
%     repFromDr = 0;

    %---------------------------------
    %-------------All MLs-------------
    %---------------------------------
    % growth rate & flux
    growRate = growthRate(xs(2),mPM,fN,h);
    fluxIn = vd(xs(2),mPM,xs(3),LD,growRate,fIII,h)*(1-vd(xs(3),mLC,xs(4),LD,growRate,fIII,h))*fIII;
    % replenish
    repFromLC = replenish(xs(3),growRate,h);
    repFromDr = replenish(xs(4),repFromLC,h);

    rateLC = fluxIn + repFromDr - repFromLC;
end

function rateALC = lcTernary( t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD, x, rPMLC, rLCDr )
    fN = fluxN(t);
    fIII = fluxIII(t);
    fluxA = fluxIII(t)*x;

    %---------------------------------
    %-------------PM & LC-------------
    %---------------------------------
%     % growth & flux
%     growRate = growthRate(xs(2),mPM,fN,h);
%     fluxAIn = vd(xs(2),mPM,xs(3),LD,growRate,fIII,h)*fluxA;
%     % replenish
%     repFromLC = replenish(xs(3),growRate,h);
%     repAFromLC = repFromLC*concA(xs(7),xs(3)-xs(7));
%     repAFromDR = 0;
%     % VCS
%     segLCPM = rPMLC*xs(6)*(xs(3)-xs(7));
%     segDrLC = 0;

    %---------------------------------
    %-------------All MLs-------------
    %---------------------------------
    % growth & flux
    growRate = growthRate(xs(2),mPM,fN,h);
    fluxAIn = vd(xs(2),mPM,xs(3),LD,growRate,fIII,h)*(1-vd(xs(3),mLC,xs(4),LD,growRate,fIII,h))*fluxA;
    % replenish
    repFromLC = replenish(xs(3),growRate,h);
    repAFromLC = repFromLC*concA(xs(7),xs(3)-xs(7));
    repFromDr = replenish(xs(4),repFromLC,h);
    repAFromDR = repFromDr*concA(xs(8),xs(4)-xs(8));
    % VCS
    segLCPM = rPMLC*xs(6)*(xs(3)-xs(7));
    segDrLC = rLCDr*xs(7)*(xs(4)-xs(8));
    
    rateALC = fluxAIn + segLCPM - repAFromLC - segDrLC + repAFromDR;
end

function rateDr = drBinary( t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD )
    fN = fluxN(t);
    fIII = fluxIII(t);

    %---------------------------------
    %-------------All MLs-------------
    %---------------------------------
    % growth & flux
    growRate = growthRate(xs(2),mPM,fN,h);
    fluxIn = vd(xs(2),mPM,xs(3),LD,growRate,fIII,h)*vd(xs(3),mLC,xs(4),LD,growRate,fIII,h)*fIII;
    % replenish
    repFromLC = replenish(xs(3),growRate,h);
    repFromDr = replenish(xs(4),repFromLC,h);
    
    rateDr = fluxIn - repFromDr;
end

function rateADr = drTernary( t, xs, fluxN, fluxIII, T, h, mPM, mLC, LD, x, rPMLC, rLCDr )
    fN = fluxN(t);
    fIII = fluxIII(t);
    fluxA = fluxIII(t)*x;

    %---------------------------------
    %-------------All MLs-------------
    %---------------------------------
    % growth & flux
    growRate = growthRate(xs(2),mPM,fN,h);
    fluxAIn = vd(xs(2),mPM,xs(3),LD,growRate,fIII,h)*vd(xs(3),mLC,xs(4),LD,growRate,fIII,h)*fluxA;    
    % replenish
    repFromLC = replenish(xs(3),growRate,h);
    repAFromLC = repFromLC*concA(xs(7),xs(3)-xs(7));
    repFromDr = replenish(xs(4),repFromLC,h);
    repAFromDR = repFromDr*concA(xs(8),xs(4)-xs(8));
    
    % VCS
    segLCPM = rPMLC*xs(6)*(xs(3)-xs(7));
    segDrLC = rLCDr*xs(7)*(xs(4)-xs(8));
    
    rateADr = fluxAIn + segDrLC - repAFromDR;
end

%% Rate Definitions
function fns = genFluxN( fluxNs, tO, tC, td, w, A, opt )
    switch opt
        case 'fixed'
            fns    = cell(length(fluxNs));
            for iter1 = 1:length(fluxNs)
                str = sprintf('@(t) (%0.2f)',fluxNs(iter1));
                fns(iter1) = {str2func(str)};
            end
        case 'sine'
            freq   = 2*pi/(w*(tO+tC));
            fns    = cell(length(fluxNs));
            for iter1 = 1:length(fluxNs)
                str = sprintf('@(t) (%0.2f + %f*sin(%f*t-%0.2f))',fluxNs(iter1),A,freq,td);
                fns(iter1) = {str2func(str)};
            end
    end
end

function fns = genFluxIII( IIIV, fluxN, tO, tC, tGrowth, h, shutterDelay, numSteps, opt, td, w, A )
    switch opt
        case 'fixed'
            fns    = cell(length(IIIV)*length(fluxN));
            for iter1 = 1:length(IIIV)
                for iter2 = 1:length(fluxN)
                    str = sprintf('@(t) (%0.2f)',IIIV(iter1)*fluxN(iter2));
                    count = iter2 + length(fluxN)*(iter1-1);
                    fns(count) = {str2func(str)};
                end
            end
        case 'abrupt'
            freq   = 2*pi/(tO+tC);
            duty   = 100*tO/(tO+tC);
            fns    = cell(length(IIIV));
            for iter1 = 1:length(IIIV)
                for iter2 = 1:length(fluxN)
                    str = sprintf('@(t) (%0.2f*(0.5*square(%f*t,%f)+0.5))',IIIV(iter1)*fluxN(iter2),freq,duty);
                    count = iter2 + length(fluxN)*(iter1-1);
                    fns(count) = {str2func(str)};
                end
            end
        case 'sine'
            freq   = 2*pi/(w*(tO+tC));
            fns    = cell(length(fluxN)*length(fluxN));
            for iter1 = 1:length(IIIV)
                for iter2 = 1:length(fluxN)
                    str = sprintf('@(t) (%0.2f + %f*sin(%f*t-%0.2f))',IIIV*fluxN(iter2),A,freq,td);
                    count = iter2 + length(fluxN)*(iter1-1);
                    fns(count) = {str2func(str)};
                end
            end
        case 'gradual'
            tCycle       = tO + tC;
            nCycles      = ceil(tGrowth/tCycle);
            fns          = cell(length(IIIV)*length(fluxN));
            for iter1 = 1:length(IIIV)
                str = '@(t) (';
                for iter2 = 1:length(fluxN)
                    for iter3 = 1:nCycles
                        open  = (iter3-1)*tCycle;
                        str = strcat(str,sprintf('myStep(t-%0.2f,%0.2f,%f,%d)+',open,IIIV(iter1)*fluxN(iter2),h,numSteps));
                        
                        clse  = open + tO;
                        str = str(1:end-1);
                        str = strcat(str,sprintf('-myStep(t-%0.2f,%0.2f,%f,%d)',clse,IIIV(iter1)*fluxN(iter2),h,numSteps));
                        if iter3 < nCycles
                            str = strcat(str,'+');
                        end
                    end
                end                
                str = strcat(str,')');
                count = iter2 + length(fluxN)*(iter1-1);
                fns(count) = {str2func(str)};
            end
    end
end

function rate = growthRateTernaryNew( PM, APM, mPM, xPM, fluxN, h )
    global preferredIncorpA growthMode;

    if strcmpi(growthMode,'pref')
        if (fluxN < 0 || PM < 0) || (APM < 0) % Error #1: if anything is negative
            ME = MException('MyException:PMUnrealistic','fluxN, PM, or APM less than 0: fluxN = %0.4g, PM = %0.4g, APM = %0.4g',...
                fluxN,PM,APM);
            throw(ME)
        elseif (PM - APM) <= -1e-15 % Error #2: if there is more catA than total (with precision of 1 atom)
            ME = MException('MyException:APMUnrealistic','APM exceeds PM: PM = %0.15g, APM = %0.15g',PM,APM);
            throw(ME)
        elseif (PM == 0 || APM == 0) % no cations trivial cases
            rate = 0;
        elseif (fluxN*h >= PM) % nitrogen rich trivial case (case #8)
            rate = APM/h;
        elseif (PM > 0 && fluxN*h < PM && (PM - APM) > -1e-15) % case #7
            growRate = growthRate(PM,mPM,fluxN,h);  % get overall growth rate
            catGrown = growRate*h;                  % find how much total is being grown
            Apref = preferredIncorpA*catGrown;      % calc how much catA would be grown if it's available
            Bpref = (1-preferredIncorpA)*catGrown;  % "    "   "    catB "     "  "     "  "    "
            BPM = PM - APM;                         % calc how much catB is in the PM ML
            if APM >= Apref     % enough catA to meet preferred incorporation ratio
                if BPM >= Bpref % and enough catB to meet preferred incorporation ratio
                    rate = Apref/h;             % grow at the preferred rate
                else            % not enough catB to meet preferred incorporation ratio
                    rate = (catGrown - BPM)/h;  % grow all available catB
                end
            else                % not enough catA to meet preferred incorporation ratio
                rate = APM/h;                   % grow all avaialble catA
            end
        else
            ME = MException('MyException:MissedCase','No cases trigerred in growthRateTernaryNew: PM = %0.4g, APM = %0.4g',PM,APM);
            throw(ME)
        end
    else
        rate = growthRate(PM,mPM,fluxN,h)*xPM;
    end
end

function rate = growthRate( PM, mPM, fluxN, h )
    pctN = 0;

    % assign growth rate based on number of atoms in the PM-ML
    if PM <=0
        rate = 0;
    elseif (PM < (pctN*(mPM - PM) + fluxN*h))
        rate = PM/h;
    else
        rate = fluxN;
    end
end

function rate = replenish( n, drivingRate, h )
    % assign replenishing rate based on number of atoms in the supplied ML
    if n <=0
        rate = 0;
    elseif (n < drivingRate*h)
        rate = n/h;
    else
        rate = drivingRate;
    end
end

%% Crystal Functions
function substrate = buildSubstrate(material)
    % import material parameters
    if cmpFileMod('genMatParams.m','materialParams.mat')
        genMatParams;
    end
    load materialParams.mat AlNProps GaNProps InNProps

    % check which material to import (default GaN)
    if (strcmpi(material,'GaN') || strcmpi(material,'AlN') || strcmpi(material,'InN'))
        material = lower(material);
    else
        material = 'gan';
    end
    
    % import material properties
    switch material
        case 'gan'
            subProps = GaNProps;
        case 'aln'
            subProps = AlNProps;
        case 'inn'
            subProps = InNProps;
    end

    % calculate substrate parameters from material properties
    ucBasalArea = (subProps.a^2)*sind(60);
    atomicLayerDensity = round(1./ucBasalArea);    % [atoms/cm^2]
    layerDensityNorm = atomicLayerDensity*1e-15;   % [10^15 atoms/cm^2]
    
    % create substrate object
    substrate = struct('cationSpecies',1,'MLCapacity',layerDensityNorm);
end

function crystalLayersClbk(old, new)
    % define relevant global variables
    global sysType crystals currentCrystal
    
    % Get change in n, catA, current time, fluxN, and h
    nCrysOld = old(2); nCrysNew = new(2);
    delN = nCrysNew - nCrysOld;
    if strcmp(sysType,"ternary")
        ACrysOld = old(6); ACrysNew = new(6);
        delA = ACrysNew - ACrysOld;
        cA = concA(ACrysOld, nCrysOld-ACrysOld);
    elseif strcmp(sysType,"quaternary")
        % incomplete
    end
    t = old(1);
    fluxNs = evalin('base','fluxNs');
    fluxN = fluxNs{1};
    fN = fluxN(t);
    h = evalin('base','h');
    
    % Initialize array to hold crystal data and set control variables
    layerCapacity = evalin('base','mPM');
    numLayers = ceil(nCrysNew/layerCapacity);                % get needed size of crystal layers array
    
    if isempty(crystals{currentCrystal})
        if strcmp(sysType,"binary")
            tempLayers = zeros(1,numLayers);
            if nCrysOld > 0
                tempLayers(1,1:end-1) = layerCapacity;
                tempLayers(1,end) = nCrysOld - sum(tempLayers(1,:));
            end
        elseif strcmp(sysType,"ternary")
            tempLayers = zeros(2,numLayers);
            if nCrysOld > 0
                tempLayers(1,1:end-1) = layerCapacity;
                tempLayers(1,end) = nCrysOld - sum(tempLayers(1,:));
                tempLayers(2,:) = (ACrysOld/nCrysOld)*tempLayers(1,:);
            end
        else % quaternary
            tempLayers = zeros(3,numLayers);
            if nCrysOld > 0
                tempLayers(1,1:end-1) = layerCapacity;
                tempLayers(1,end) = nCrysOld - sum(tempLayers(1,:));
                % add more to define quaternary
            end
        end
        crystals{currentCrystal} = tempLayers;
    end
    
    % Update crystal
    layers = crystals{currentCrystal};                         % fill with the previously generated crystal
    layerInd = size(layers,2);                                 % get index of the last layer with material in it
    if strcmp(sysType,"binary")
        if delN == 0 % no change in crystal
            return;
        elseif delN < 0 % total cations decreased
            if numLayers < length(layers) % cations decreased enough to lose a layer of crystal
                remN = abs(delN) - layers(1,layerInd); % how much is left to take away after top layer depletes
                layers = layers(:,1:end-1); % delete last layer
                layerInd = layerInd - 1;    % adjust layer index to new top layer
                layers(1,layerInd) = layers(1,layerInd) - remN; % remove remaining cations
            else % cations decreased within current top layer
                layers(1,layerInd) = layers(1,layerInd) + delN; % remove cations
            end
        elseif delN > 0 % total cations increased
            if numLayers > length(layers) % cations increased enough to overflow into a new layer
                nIntoCurrent = layerCapacity - layers(1,layerInd); % how many cations are used to fill top layer
                remN = delN - nIntoCurrent; % how much is left to add after current layer fills
                layers(1,layerInd) = layerCapacity; % fill current layer
                layers = padarray(layers,[0,1],'post');  % add a new layer
                layerInd = layerInd + 1;    % adjust layer index to new top layer
                layers(1,layerInd) = remN; % add remaining cations
            else
                layers(1,layerInd) = layers(1,layerInd) + delN; % add cations
            end
        end
    elseif strcmp(sysType,"ternary")
        if( delN == 0 && delA == 0 ) % no change in crystal
            return;
        elseif delN < 0 % total cations decreased
            if delA < 0 % and catA decreased
                if numLayers < length(layers) % cations decreased enough to lose a layer of crystal
                    remN = abs(delN) - layers(1,layerInd); % how much is left to take away after top layer depletes
                    remA = abs(delA) - layers(2,layerInd); % how much is left to take away after top layer depletes
                    layers = layers(:,1:end-1); % delete last layer
                    layerInd = layerInd - 1;    % adjust layer index to new top layer
                    layers(1,layerInd) = layers(1,layerInd) - remN; % remove remaining cations
                    layers(2,layerInd) = layers(2,layerInd) - remA; % remove remaining catA
                else % cations decreased within current top layer
                    layers(1,layerInd) = layers(1,layerInd) + delN; % remove cations
                    layers(2,layerInd) = layers(2,layerInd) + delA; % remove catA
                end
            else % and catA stayed the same or increased (error)
                ME = MException('MyException:UnrealisticCrystal','delN < 0 and delA >= 0');
                throw(ME)
            end
        elseif delN == 0 % total cations stayed constant
            if delA < 0 % and catA decreased
                if abs(delA) > layers(2,layerInd) % cations decreased enough to lose a layer of crystal
                    remA = abs(delA) - layers(2,layerInd); % how much is left to take away after top layer is depleted of catA
                    layers(2,layerInd) = 0; % remove catA from top layer
                    layerInd = layerInd - 1;    % adjust layer index to second layer from top
                    layers(2,layerInd) = layers(2,layerInd) - remA; % remove remaining catA
                else % cations decreased within current top layer
                    layers(2,layerInd) = layers(2,layerInd) + delA; % remove catA
                end
            else % and catA increased (error)
%                 ME = MException('MyException:UnrealisticCrystal','delN = 0 and delA >= 0');
%                 throw(ME)
                  warning('tnew = %0.3f :: delN = 0 and delA >= 0',new(1))
            end
        elseif delN > 0 % total cations increased
            if delA == 0 % and catA stayed constant
                if numLayers > length(layers) % cations increased enough to overflow into a new layer
                    nIntoCurrent = layerCapacity - layers(1,layerInd); % how many cations are used to fill top layer
                    remN = delN - nIntoCurrent; % how much is left to add after current layer fills
                    layers(1,layerInd) = layerCapacity; % fill current layer
                    layers = padarray(layers,[0,1],'post');  % add a new layer
                    layerInd = layerInd + 1;    % adjust layer index to new top layer
                    layers(1,layerInd) = remN; % add remaining cations
                else
                    layers(1,layerInd) = layers(1,layerInd) + delN; % add cations
                end
            elseif delA < 0
                % deplete catA from top (or top 2) layers
                if abs(delA) > layers(2,layerInd) % cations decreased enough to lose a layer of crystal
                    remA = abs(delA) - layers(2,layerInd); % how much is left to take away after top layer is depleted of catA
                    layers(2,layerInd) = 0; % remove catA from top layer
                    layerInd = layerInd - 1;    % adjust layer index to second layer from top
                    layers(2,layerInd) = layers(2,layerInd) - remA; % remove remaining catA
                else % cations decreased within current top layer
                    layers(2,layerInd) = layers(2,layerInd) + delA; % remove catA
                end
                % add new cations
                if numLayers > length(layers) % cations increased enough to overflow into a new layer
                    nIntoCurrent = layerCapacity - layers(1,layerInd); % how many cations are used to fill top layer
                    remN = delN - nIntoCurrent; % how much is left to add after current layer fills
                    layers(1,layerInd) = layerCapacity; % fill current layer
                    layers = padarray(layers,[0,1],'post');  % add a new layer
                    layerInd = layerInd + 1;    % adjust layer index to new top layer
                    layers(1,layerInd) = remN; % add remaining cations
                else
                    layers(1,layerInd) = layers(1,layerInd) + delN; % add cations
                end
            else % delA > 0
                if numLayers > length(layers) % cations increased enough to overflow into a new layer
                    nIntoCurrent = layerCapacity - layers(1,layerInd); % how many cations are used to fill top layer
                    remN = delN - nIntoCurrent; % how much is left to add after current layer fills
                    %__ needs to be adjusted for adding decomp (see OneDrive pg 30)
                    growRate = growthRate(old(3), layerCapacity, fN, h);
                    growRateA = growthRateTernaryNew(old(3), old(7), layerCapacity, concA(old(7), old(3)-old(7)), fN, h);
                    AintoCurrent = (growRateA/growRate)*(nIntoCurrent);
                    remA = delA - AintoCurrent;
                    if remA < 0 % or other error cases
%                         ME = MException('MyException:UnrealisticCrystal','RemA calc error in delN & delA > 0 case');
%                         throw(ME)
                        warning('tnew = %0.3f :: RemA calc error in delN & delA > 0 case',new(1))
                    elseif remA > remN
                        warning('tnew = %0.3f :: RemA (%0.8f) > RemN (%0.8f)',new(1),remA,remN)
                        remA = remN;
                    end
                    %__ end of above comment
                    layers(1,layerInd) = layerCapacity; % fill current layer
                    layers(2,layerInd) = layers(2,layerInd) + AintoCurrent; % add catA to current
                    layers = padarray(layers,[0,1],'post');  % add a new layer
                    layerInd = layerInd + 1;    % adjust layer index to new top layer
                    layers(1,layerInd) = remN; % add remaining cations
                    layers(2,layerInd) = remA; % add remaining catA
                else % cations only increased within current layer
                    layers(1,layerInd) = layers(1,layerInd) + delN; % add cations
                    layers(2,layerInd) = layers(2,layerInd) + delA; % add catA
                end
            end
        end
    elseif strcmp(sysType,"quaternary")
        % incomplete
    end
    
    % save results to global crystal variable
    crystals{currentCrystal} = layers;
end

%% Support Functions
function C = concA( A, B )
    if A <= 0 && B <= 0
        C = 0;
    else
        C = A/(A+B);
    end
end

function P = vd( n, m, nUpper, LD, GR, fluxM, h )
    % function to find the probability of an atom finding an open site in a ML
    % uses the average vacancy spacing from # of atoms in a finite ML
    % scale m to the correct OoM
    if (m < 1000)
        nV = (m-n)*1e15;
    else
        nV = (m-n);
    end
    
    if fluxM == 0
        P = 1;
    else
        if n >= m
            if (nUpper >= GR*h) % PM and LC MLs are full
                P = 1;
            else          % Only PM ML is full
                P = 1 + ((nUpper/h) - GR)/fluxM;  % constrain f such that all material into the layer is equal to the growth rate out
            end
        else
            LV = 1/sqrt(nV);
            P = LV/(LD+LV);
        end
    end
end

%% Plot Functions
function [f,f2] = plotMLs(data,label,counter)
    % get relevant values (global or pull in from base workspace)
    global sysType cations;
    mPM = evalin('base','mPM');
    mLC = evalin('base','mLC');
    tOpen = evalin('base','tOpen');
    tClosed = evalin('base','tClosed');
    t0 = evalin('base','t0');
    tf = evalin('base','tf');

    % calculate intervals for shading
    numCycles = ceil(tf/(tOpen + tClosed));
    tOs = zeros([numCycles,1]);
    tCs = zeros([numCycles,1]);
    tOs(1) = t0;
    tCs(1) = t0 + tOpen;
    for iter1 = 2:numCycles
        tOs(iter1) = tCs(iter1-1) + tClosed;
        tCs(iter1) = tOs(iter1) + tOpen;
    end
    if tCs(end) > tf
        tCs(end) = tf;
    end
    frameBounds = [tOs tCs];

    % generate full shading coordinates for plots
    windows = cell(numCycles);
    ys = [0 0 1 1];
    for iter1 = 1:numCycles
        windows{iter1} = [frameBounds(iter1,:) flip(frameBounds(iter1,:))];
    end

    % plot # atoms in crystal vs. time
    fa = figure(10*counter);
        hold on;
        set(gca,'FontSize',18,'FontWeight','bold');
        % plot data
        if strcmp(sysType,'ternary')
            plot(data.t,data.ACrys,'Color','#A2142F',...
                'DisplayName',sprintf('%s Atoms',cations(1)),'LineWidth',2);
            plot(data.t,data.VCrys-data.ACrys,'Color','b',...
                'DisplayName',sprintf('%s Atoms',cations(2)),'LineWidth',2);
        end
        plot(data.t,data.VCrys,'k',...
            'DisplayName','Cation Atoms','LineWidth',2);
        % plot shutter intervals
        maxY = ylim(gca);
        for iter1 = 1:numCycles
            patch(gca,windows{iter1},ys*maxY(2),'blue','FaceAlpha',0.3,'EdgeColor','none','DisplayName','none');
        end
        % legend
        lgd = legend('FontSize',14);
        lgd.String = lgd.String(1:end-1);
        lgd.Location = 'best';
        % re-order draw order
        allChildren = get(gca,'Children');
        newChildren = [allChildren(2:end); allChildren(1)];
        set(gca,'Children',newChildren);
        % labels
        yLa = ylabel('Crystal (10^{15} at. cm^{-2})','FontSize',20);
        yLa.Position(1) = yLa.Position(1) - 0.75;
        % axis limits
        ylim(gca,maxY);
    f1a = gca;

    % plot # atoms in PM-ML vs. time
    fb = figure(10*counter+1);
        hold on;
        set(gca,'FontSize',18,'FontWeight','bold');
        % plot data
        if strcmp(sysType,'ternary')
            plot(data.t,data.APM,'Color','#A2142F','LineWidth',2);
            plot(data.t,data.VPM-data.APM,'Color','b','LineWidth',2);
        end
        plot(data.t,data.VPM,'k','LineWidth',2);
        % plot shutter intervals
        maxY = ylim(gca);
        for iter1 = 1:numCycles
            patch(gca,windows{iter1},ys*maxY(2),'blue','FaceAlpha',0.3,'EdgeColor','none');
        end
        % re-order draw order
        allChildren = get(gca,'Children');
        newChildren = [allChildren(2:end); allChildren(1)];
        set(gca,'Children',newChildren);
        % axis limits
        ylim(gca,maxY);
        % labels
        lYb = yline(mPM,'r--','PM-ML Capacity','LabelVerticalAlignment','bottom','LineWidth',2);
        lYb.FontSize = 11;
        lYb.FontWeight = 'bold';
        yLb = ylabel('PM-ML (10^{15} at. cm^{-2})','FontSize',20);
        yLb.Position(1) = yLb.Position(1) - 0.75;
    f1b = gca;
    
    % plot # atoms in LC-ML vs. time
    fc = figure(10*counter+2);
        hold on;
        set(gca,'FontSize',18,'FontWeight','bold');
        % plot data
        if strcmp(sysType,'ternary')
            plot(data.t,data.ALC,'Color','#A2142F','LineWidth',2);
            plot(data.t,data.VLC-data.ALC,'Color','b','LineWidth',2);
        end
        plot(data.t,data.VLC,'k','LineWidth',2);
        % plot shutter intervals
        maxY = ylim(gca);
        for iter1 = 1:numCycles
            patch(gca,windows{iter1},ys*maxY(2),'blue','FaceAlpha',0.3,'EdgeColor','none');
        end
        % re-order draw order
        allChildren = get(gca,'Children');
        newChildren = [allChildren(2:end); allChildren(1)];
        set(gca,'Children',newChildren);
        % labels
        lYc = yline(mLC,'r--','LC-ML Capacity','LabelVerticalAlignment','bottom','LineWidth',2);
        lYc.FontSize = 11;
        lYc.FontWeight = 'bold';
        yLc = ylabel('LC ML (10^{15} at. cm^{-2})','FontSize',20);
        yLc.Position(1) = yLc.Position(1) - 0.75;
        xlabel('time (s)','FontSize',20)
        % axis limits
        ylim(gca,maxY);
    f1c = gca;
    
    % plot # atoms in droplets vs. time
    fd = figure(10*counter+3);
        hold on;
        set(gca,'FontSize',18,'FontWeight','bold');
        % plot data
        if strcmp(sysType,'ternary')
            plot(data.t,data.ADroplet,'Color','#A2142F','LineWidth',2);
            plot(data.t,data.VDroplet-data.ADroplet,'Color','b','LineWidth',2);
        end
        plot(data.t,data.VDroplet,'k','LineWidth',2);
        % plot shutter intervals
        maxY = ylim(gca);
        for iter1 = 1:numCycles
            patch(gca,windows{iter1},ys*maxY(2),'blue','FaceAlpha',0.3,'EdgeColor','none');
        end
        % re-order draw order
        allChildren = get(gca,'Children');
        newChildren = [allChildren(2:end); allChildren(1)];
        set(gca,'Children',newChildren);
        % labels
        yLd = ylabel('Droplets (10^{15} at. cm^{-2})','FontSize',20);
        yLd.Position(1) = yLd.Position(1) - 0.75;
        xlabel('time (s)','FontSize',20)
        % axis limits
        ylim(gca,maxY);
    f1d = gca;
    
    % combine previous plots into 1 figure
    f = figure(counter);
    set(f,'units','normalized','outerposition',[0.16 0.1 0.68 0.8]);
    set(gca,'Color','none','XColor','none','YColor','none')
    copies = copyobj([f1a,lgd],f);
    f1a_copy = copies(1);
    subplot(2,2,1,f1a_copy)
    f1b_copy = copyobj(f1b,f);
    subplot(2,2,2,f1b_copy)
    f1c_copy = copyobj(f1c,f);
    subplot(2,2,3,f1c_copy)
    f1d_copy = copyobj(f1d,f);
    subplot(2,2,4,f1d_copy)
    close(fa);close(fc);close(fd);

    % plot PM-ML vs time on its own
    f2 = figure(counter+1);
    set(f2,'units','normalized','outerposition',[0.16 0.1 0.68 0.8]);
    yLb.Position(1) = yLb.Position(1) + 1.75;
    copyobj(f1b,f2)
    xlabel('time (s)','FontSize',20)
    close(fb);
end

function f = plotCrystal(crys,mPM,counter,labelsTF)
    global cations
    
    f = figure(counter);
    cla;
    set(f,'units','normalized','outerposition',[0.16 0.1 0.68 0.8]);
    nL = 1:length(crys(1,:));
    b1 = bar(nL, crys(1,:)/mPM,'FaceColor','#D3D3D3');      % plot # total atoms in each layer
    hold on;
    b2 = bar(nL, crys(2,:)/mPM,'FaceColor','#A2142F');      % plot # A cations in each layer
    ylim([0,1]);
    set(gca,'FontSize',18,'FontWeight','bold');
    xL = xlabel('# of Grown Monolayers','FontSize',20);
    yL = ylabel(sprintf('%s Composition (%%)',cations(1)),'FontSize',20);
    yL.Position(1) = yL.Position(1) - 0.8;

    % add labels for the composition of each layer onto the A cation bars
    if labelsTF
        xtips1 = b2(1).XEndPoints;
        ytips1 = b2(1).YEndPoints;
        compAVals = 100*round(crys(2,:)./crys(1,:),3);
        labels1 = string(compAVals);
        labels1 = strcat(labels1,'%');
        text(xtips1,ytips1,labels1,'Rotation',90,'HorizontalAlignment','right',...
            'VerticalAlignment','middle','FontSize',16,'FontWeight','bold','Color','w');
    end
end

%% Helper Scripts
function status = appendFiles( readFile, writtenFile )
      [fr,err1] = fopen(readFile,'rt');
      [fw,err2] = fopen(writtenFile,'at');
      while feof(fr) == 0
          tline = fgetl(fr);
          fwrite(fw,sprintf('%s\n',tline));
      end
      fclose(fr);
      fclose(fw);
end

function prepend2file( string, filename, newline )
% newline:  is an optional boolean, if true will append a \n to the end 
% of the string that is sent in such that the original text starts on the 
% next line rather than the end of the line that the string is on 
% string:  a single line string 
% filename:  the file you want to prepend to 
      tempFile = tempname;
      fw = fopen(tempFile,'wt');
      if nargin < 3
          newline = true;
      end
      if newline
          fwrite(fw,sprintf('%s\n',string));
      else
          fwrite(fw,string);
      end
      fclose(fw);
      appendFiles(filename,tempFile);
      copyfile(tempFile,filename);
      delete(tempFile);
end

function tf = cmpFileMod(f1, f2)
    % checks whether f1 has been modified since f2 was last modified
    d1 = datetime(dir(f1).date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    d2 = datetime(dir(f2).date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    if (d1 >= d2)
        tf = true;
    else
        tf = false;
    end
end
