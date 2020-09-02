function [lissom,opts] = define(neu_size)
%%  Define the LISSOM  architecture
%==============================================

x.layers = {
    struct('type', 'i') %input layer
    struct('type', 'v', 'dim', neu_size,'rf',[9 9],'exe_rad',3,'inhb_rad',10,'stride',1); %v1 layer
    %struct('type', 'v', 'dim', [50 50],'rf',[7 7],'exe_rad',3,'inhb_rad',6,'stride',1); %v2 layer
    }; 
 % Options for training.
 lissom =x;
% opts.p=[2];
% opts.q=[3];
% opts.r=[5];
% opts.s=[1.5];
% % opts.p=[1];
% % opts.q=[1];
% % opts.r=[1];
opts.p=[2];
opts.q=[9];
opts.r=[50];
opts.s=[1.5];

    
% opts.p2=[1];
% opts.q2=[2];
% opts.r2=[1.5];
opts.all=[0 0];
opts.control=0; % control=1 for feedback connection.

opts.epoch=10;
opts.niter=1; %general settelling time will be 10 to 15 time steps:assumption

%%===========================================================================
end