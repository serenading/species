%% Script plots sample trajectories from a single recording (4x4 wells)

% author: @serenading. Jan 2021.

clear
close all

addpath('../AggScreening/')


%% Set parameters
strain = 'QX1410'; % 'N2','CB4856','MY23','QX1410','VX34','NIC58','JU1373'; 
n_subsample = 1;