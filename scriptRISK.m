% Simple script to build the risk-informed TLP data structure.
clear;

% Define input parameters.
latG=(49:0.1:60);
lonG=(-120:0.1:-110);
latE=(49:0.2:60)+0.05;
lonE=(-120:0.2:-110)+0.05;
ML=0.1:0.2:6.0;
SAfile='temp';
DEPfile='temp';
POPfile='temp';
PlayFlag='SYN';
ReMAXn=300;
ReMAXd=30;
Nv=5;
Nt=100;
PSA_f=0;
rand_flag='random';

% Load in play boundaries.
latB=[51 57 57 52 51];
lonB=[-119 -119 -111 -112 -119];

% Define data structure for TLP map.
S=setupSTRUCT(PlayFlag,latG,lonG,latE,lonE,latB,lonB,ML,ReMAXn,ReMAXd,PSA_f);

% Load & grid all map data.
% Note that this is arbitrary data, to provide an example.  Reald data
% would need to be substituted for a user's specific purpose.
S=loadSA(S,SAfile);
S=loadDEP(S,DEPfile);
S=loadPOP(S,POPfile);

% Load in data, to continue iterations.
%load('Rmap.mat','S');

% Iteratively add vulnerability curves.
while(length(S.dVAR.dM)<Nt)
    
    % Prompt for percent done.
    100*length(S.dVAR.dM)/Nt
    
    % Create a perturbed data structure.
    S=perturbVAR(S,Nv,rand_flag);
    
    % Compute risk curves for each spatial pixel and perturbed value.
    tic; S=runRISK(S,rand_flag,PSA_f); toc;
    
    % Save data structure.
    save('Rmap_temp.mat','S');
end



