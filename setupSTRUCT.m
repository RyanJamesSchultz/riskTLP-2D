function S=setupSTRUCT(PLAYflag,latMAP,lonMAP,latEQ,lonEQ,latBOUN,lonBOUN,Mw,ReMAXn,ReMAXd,PSAf)
  % Simple function that defines the TLP risk-map data structure.
  % 
  % Writtern by Ryan Schultz.
  
  % Top level structure defintion.
  S=struct('play_flag',[],'MAP',[],'dVAR',[],'RISK',[],'ML',[],'Mw',[],'PSAf',[]);
  
  % Define ML magnitude ranges and convert to Mw.
  % Convert from ML to Mw (Luca & Mufano, 2018; doi: 10.1785/0120170303).
  %Mw=(2/3)*ML+1.14;
  %Mw(ML>4.3)=1.28*ML(ML>4.3)-1.50;
  ML=Mw;
  
  % Put information in the main structure, S.
  S.Mw=Mw;
  S.ML=ML;
  S.PSAf=PSAf;
  S.play_flag=PLAYflag;
  
  % Second level structure definitions.
  MAP=struct('latG',[],'lonG',[],'Vs30',[],'dVs30',[],'POP',[],'DEP',[],'latE',[],'lonE',[],'latB',[],'lonB',[],'Ir',[]);
  dVAR=struct('dM',[],'b',[],'dGM',[],'dF1',[],'dF2',[],'dSA',[],'dPOP',[],'dZ',[],'UPDATEflag',[]);
  RISK=struct('lat',[],'lon',[],'Nn2',[],'Nn3',[],'Nn4',[],'Nn5',[],'Nn6',[],'Nd1',[],'Nd2',[],'Nd3',[],'Nd4',[]);
  
  % Put information in the MAP structure.
  MAP.latG=latMAP;
  MAP.lonG=lonMAP;
  MAP.latE=latEQ;
  MAP.lonE=lonEQ;
  MAP.latB=latBOUN;
  MAP.lonB=lonBOUN;
  MAP.ReN_max=ReMAXn;
  MAP.ReD_max=ReMAXd;
  
  % Find the map lats/longs that are within the play boundary.
  LAT=repmat(latEQ',1,length(lonEQ));
  LON=repmat(lonEQ,length(latEQ),1);
  MAP.Ir=inpolygon(LON,LAT,lonBOUN,latBOUN);
  
  % Make a list of those lat/longs
  lat=LAT(MAP.Ir); lat=lat(:);
  lon=LON(MAP.Ir); lon=lon(:);
  
  % Loop over all in-play coords, and stuff into risk structure.
  for i=1:length(lat)
      RISK(i).lat=lat(i);
      RISK(i).lon=lon(i);
  end
  
  % Stuff second level structures into S.
  S.MAP=MAP;
  S.dVAR=dVAR;
  S.RISK=RISK;
  
return