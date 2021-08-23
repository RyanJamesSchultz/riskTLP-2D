% Simple script to map the 50-50 nuisance and damage contours.
clear;

% Load in the data structure: S.
load('Rmap.mat');

% Predefine some variables.
Mw=5.25;
dGM=0.0;
Hn=10.^[-2.8974 -2.2947 -1.6121 -1.1385 -1.0078];
Hd=[0.2 0.4 0.8 1.6];
PSAf=0;

% Get some dependent variables.
latE=mean(S.MAP.latE);
lonE=mean(S.MAP.lonE);
dep=interp2(S.MAP.lonE,S.MAP.latE,S.MAP.DEP,lonE,latE,'linear');
Nx=length(S.MAP.lonG);
Ny=length(S.MAP.latG);
latG=repmat(S.MAP.latG',1,Nx); latG=latG(:);
lonG=repmat(S.MAP.lonG,Ny,1); lonG=lonG(:);
vs30=S.MAP.Vs30(:);

% Get the flag for the GMPE to use.
if(strcmpi(S.play_flag,'SYN'))
    GMPEflag='a15';
end

% Compute the earthquake ground motions.
Re=Geoid_Distance(latE,lonE,latG,lonG,'elliptical')*6371*pi()/180;
pgv=GMPE(Re,Mw,dep,vs30,dGM,  -1,GMPEflag)*0.01;
psa=GMPE(Re,Mw,dep,vs30,dGM,PSAf,GMPEflag)/980.665;

% Reshape.
Re=reshape(Re,Ny,Nx);
pgv=reshape(pgv,Ny,Nx);
psa=reshape(psa,Ny,Nx);

% Plot the results.
figure(7); clf;
contourf(S.MAP.lonG,S.MAP.latG,S.MAP.Vs30,'LineColor','none'); hold on;
h = colorbar(); ylabel(h, 'Site Amplififcation, Vs30 (m/s)');
contour(S.MAP.lonG,S.MAP.latG,pgv,Hn,'-b');
contour(S.MAP.lonG,S.MAP.latG,psa,Hd,':w');
plot(S.MAP.lonB,S.MAP.latB,'-k');
plot(lonE,latE,'wp','MarkerFaceColor','k');
xlabel('Longitude'); ylabel('Latitude'); title('Site Amplification');

