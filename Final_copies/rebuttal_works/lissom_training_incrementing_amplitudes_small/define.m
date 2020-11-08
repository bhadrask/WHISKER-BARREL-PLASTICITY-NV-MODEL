function [lissom,opts] = define(neu_size)
%%  Define the LISSOM  architecture
%==============================================

x.layers = {
    struct('type', 'i') %input layer
    struct('type', 'v', 'dim', neu_size,'rf',[53 53],'exe_rad',8,'inhb_rad',32,'stride',1); %v1 layer
    %struct('type', 'v', 'dim', [50 50],'rf',[7 7],'exe_rad',3,'inhb_rad',6,'stride',1); %v2 layer
    }; 
 % Options for training.
lissom =x;
opts.p=[1];
opts.q=[10];
opts.r=[8];

opts.all=[0 0];
opts.control=0; % control=1 for feedback connection.
% opts.alphal=[0.3];
% opts.alphau=[0.6];     %parameters for the dynamics and weight update
opts.etaaff=[0.01];
opts.etaexc=[0.01];
opts.etainhib=[0.01];
opts.epoch=50;
opts.niter=15; %general settelling time will be 10 to 15 time steps:assumption

%%===========================================================================
end