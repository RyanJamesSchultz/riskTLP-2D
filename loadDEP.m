function S=loadDEP(S,filename)
  % Simple function that loads in formation depth data for the map-bounded area.
  
  % Predefine some values.
  Nx=length(S.MAP.lonE);
  Ny=length(S.MAP.latE);
  
  % Just make depth a simple, spatially dependent function.
  dep=3.0*ones([Ny Nx])+repmat(linspace(1,0,Ny)',1,Nx).*repmat(linspace(0,1,Nx),Ny,1);
  
  % Stuff results into structure.
  S.MAP.DEP=dep;
  
return