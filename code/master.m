%% Inputs
addpath(genpath(pwd))
env = load('../input/env.mat');
par = load('../input/par.mat');
init = load('../input/init.mat');
exper = load('../input/exper.mat');
sol = load('../input/sol.mat');
figs = load('../input/figs.mat');
filename = '_';

%% Scripts
sp.exper=1; % 2 3
solver
out_robot
out_tax