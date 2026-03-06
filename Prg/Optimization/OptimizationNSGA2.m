function [optim_results] = OptimizationNSGA2(data,geom,test,obs,method)

% move to folder which contains the NGSA2 algorythm
addpath('./Prg/Optimization/NSGA2');

%create optim vector
defaultopt  = nsgaopt(geom,test);

%call optim function
optim_results     = nsga2(defaultopt,geom,test,obs,method);

%go back to main folder
rmpath('./Prg/Optimization/NSGA2');
