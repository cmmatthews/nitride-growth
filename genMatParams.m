function genMatParams
    clearvars; clc; close all;
    
    % GaN Lattice Parameters
    aGaN = 3.19e-8; % cm (ref. Qian et al 1996)
    cGaN = 5.18e-8; % cm (ref. Qian et al 1996)
    rhoAtomicGaN = 8.9e22; % at/cm^2
    aThermalExpansionCoeffGaN = 5.59e-6;    % K^-1 (ref. Qian et al 1996)
    cThermalExpansionCoeffGaN = 3.17e-6;    % K^-1 (ref. Qian et al 1996)
    GaNProps = struct('a',aGaN,'c',cGaN,'rho',rhoAtomicGaN,...
        'alphaA',aThermalExpansionCoeffGaN,'alphaC',cThermalExpansionCoeffGaN);
    
    % AlN Lattice Parameters
    aAlN = 3.11e-8; % cm (ref. Goldberg 2001)
    cAlN = 4.98e-8; % cm (ref. Goldberg 2001)
    rhoAtomicAlN = 9.58e22; % at/cm^2 
    aThermalExpansionCoeffAlN = 4.15e-6;    % K^-1 (ref. Sirota & Goldushko 1974 & Qian et al 1996)
    cThermalExpansionCoeffAlN = 5.27e-6;    % K^-1 (from 20 C to 800 C) (ref. Sirota & Goldushko 1974 & Qian et al 1996)
    AlNProps = struct('a',aAlN,'c',cAlN,'rho',rhoAtomicAlN,...
        'alphaA',aThermalExpansionCoeffAlN,'alphaC',cThermalExpansionCoeffAlN);
    
    % InN Lattice Parameters
    aInN = 3.53e-8; % cm (ref. Zubrilov 2001)
    cInN = 5.69e-8; % cm (ref. Zubrilov 2001)
    rhoAtomicInN = 6.4e22; % at/cm^2
    aThermalExpansionCoeffInN = 0;  % K^-1 (unknown)
    cThermalExpansionCoeffInN = 0;  % K^-1 (unknown)
    InNProps = struct('a',aInN,'c',cInN,'rho',rhoAtomicInN,...
        'alphaA',aThermalExpansionCoeffInN,'alphaC',cThermalExpansionCoeffInN);
    
    % Export to matfile
    save materialParams.mat GaNProps AlNProps InNProps
end