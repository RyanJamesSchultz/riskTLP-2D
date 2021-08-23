function S=loadSA(S,filename)
  % Simple function that loads in site amplification data for the map-bounded area.
  
  % Predefine some values.
  Nx=length(S.MAP.lonG);
  Ny=length(S.MAP.latG);
  
  % Make an arbitrary grid of Vs30 values.
  vs30=lognrnd(0.0,0.1,[Ny Nx]).*repmat(linspace(500,800,Nx),Ny,1);
  dvs30=0.5*sqrt(vs30);
  
  % Stuff SA into output structure.
  S.MAP.Vs30 = vs30;
  S.MAP.dVs30=dvs30;
  
return