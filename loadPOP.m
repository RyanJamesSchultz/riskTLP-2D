function S=loadPOP(S,filename)
  % Simple function that loads in population data for the map-bounded area.
  
  % Predefine some values.
  Nx=length(S.MAP.lonG);
  Ny=length(S.MAP.latG);
  dx=Geoid_Distance(S.MAP.latG(round(Ny/2)),S.MAP.lonG(round(Nx/2)), S.MAP.latG(round(Ny/2)),S.MAP.lonG(round(Nx/2)+1),'elliptical')*6371*pi()/180;
  dy=Geoid_Distance(S.MAP.latG(round(Ny/2)),S.MAP.lonG(round(Nx/2)), S.MAP.latG(round(Ny/2)+1),S.MAP.lonG(round(Nx/2)),'elliptical')*6371*pi()/180;
  Nu=1*(dx*dy);
  Nr=10*(dx*dy);
  Nc=100*(dx*dy);
  
  % Make everywhere have a base 'uninhabitied' population.
  pop=Nu*ones([Ny Nx]);
  
  % Make some areas have a 'rural' population.
  latC=[54 58 58 54 54];
  lonC=[-119.5 -119.5 -115.5 -115.5 -119.5];
  LAT=repmat(S.MAP.latG',1,Nx);
  LON=repmat(S.MAP.lonG,Ny,1);
  Ir=inpolygon(LON(:),LAT(:),lonC,latC);
  pop(Ir)=Nr;
  
  % Make one area have an 'urban' population.
  latC=[54.75 55.25 55.25 54.75 54.75];
  lonC=[-116.25 -116.25 -115.75 -115.75 -116.25];
  LAT=repmat(S.MAP.latG',1,Nx);
  LON=repmat(S.MAP.lonG,Ny,1);
  Ic=inpolygon(LON(:),LAT(:),lonC,latC);
  pop(Ic)=Nc;
  
  % Stuff results into structure.
  S.MAP.POP=pop;
  
return