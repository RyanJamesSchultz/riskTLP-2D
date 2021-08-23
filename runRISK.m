function S=runRISK(S,rand_flag,psa_f)
  % Compute the risk curves.
  % 
  % Written by Ryan Schultz.
  
  % Check that we're not already up to date.
  if(strcmpi(S.dVAR.UPDATEflag,'no'))
      return;
  end
  
  % Get lengths of map edges.
  Nx=length(S.MAP.lonG);
  Ny=length(S.MAP.latG);
  
  % Get lists of all of the lat/long coords of interest.
  latG=repmat(S.MAP.latG',1,Nx); latG=latG(:);
  lonG=repmat(S.MAP.lonG,Ny,1); lonG=lonG(:);
  latE=[S.RISK.lat];
  lonE=[S.RISK.lon];
  
  % Get the map properties of interest.
  DEP=interp2(S.MAP.lonE,S.MAP.latE,S.MAP.DEP, lonE,latE,'linear');
  VS30=S.MAP.Vs30(:);
  DVS30=S.MAP.dVs30(:);
  POP=S.MAP.POP(:);
  
  % Find the number of iterations needed.
  Ne=length(latE);
  Ng=length(latG);
  Nv=length(S.dVAR.dM);
  Nm=length(S.ML);
  
  % Get the flag for the GMPE to use.
  if(strcmpi(S.play_flag,'SYN'))
      GMPEflag='a15';
  end
  
  % Loop over all of the (new) perturbed values.
  ns=size(S.RISK(1).Nn2,1)+1;
  for i=ns:Nv
      
      if(strcmpi(rand_flag,'random'))
          % Get all information and perturb it.
          M=S.Mw+S.dVAR.dM(i);     % [1 Nm]
          dGM=S.dVAR.dGM(i);       % [1 1]
          dSA=S.dVAR.dSA(i);       % [1 1]
          dep=DEP+S.dVAR.dZ(i);    % [Ne 1]
          vs30=VS30.*lognrnd(0.0,log10(exp(DVS30)));       % [Ng 1]
          pop=abs(normrnd(POP,sqrt(POP)))*S.dVAR.dPOP(i);  % [Ng 1]
          dF=[S.dVAR.dF1(i) S.dVAR.dF2(i)];
      elseif(strcmpi(rand_flag,'none'))
          % Get all information and perturb it.
          M=S.Mw;                  % [1 Nm]
          dGM=S.dVAR.dGM(i);       % [1 1]
          dSA=S.dVAR.dSA(i);       % [1 1]
          dep=DEP;                 % [Ne 1]
          vs30=VS30;               % [Ng 1]
          pop=POP;                 % [Ng 1]
          dF=[S.dVAR.dF1(i) S.dVAR.dF2(i)];
      end
      
      % Reshape information into matrices [Ng Nm].
      M=repmat(M,Ng,1);
      vs30=repmat(vs30,1,Nm);
      pop=repmat(pop,1,Nm);
      
      % Loop over all of the map pixels.
      for j=1:Ne
          
          % Get distances and reshape into matrix (km).
          Re=Geoid_Distance(latE(j),lonE(j),latG,lonG,'elliptical')*6371*pi()/180; % [Ng 1]
          Re=repmat(Re,1,Nm); % [Ng Nm]
          
          % Truncate based on maximum distance.
          In=(Re(:,1)<=S.MAP.ReN_max);
          Id=(Re(:,1)<=S.MAP.ReD_max);
          
          % Compute the ground motion matrices (PGV:m/s & PGA:g).
          pgv=GMPE(Re(In,:),M(In,:),dep(j),vs30(In,:),dGM,   -1,GMPEflag)*dSA*0.01;      % [Ng(In) Nm]
          psa=GMPE(Re(Id,:),M(Id,:),dep(j),vs30(Id,:),dGM,psa_f,GMPEflag)*dSA/980.665;   % [Ng(Id) Nm]
          
          % Compute chance of nuisance observation [Ng(In) Nm].
          On2=RESPfxn(pgv,dF,'nuisance',2,psa_f);
          On3=RESPfxn(pgv,dF,'nuisance',3,psa_f);
          On4=RESPfxn(pgv,dF,'nuisance',4,psa_f);
          On5=RESPfxn(pgv,dF,'nuisance',5,psa_f);
          On6=RESPfxn(pgv,dF,'nuisance',6,psa_f);
          
          % Compute chance of damage observation [Ng(Id) Nm].
          Od1=RESPfxn(psa,dF,'damage',1,psa_f);
          Od2=RESPfxn(psa,dF,'damage',2,psa_f);
          Od3=RESPfxn(psa,dF,'damage',3,psa_f);
          Od4=RESPfxn(psa,dF,'damage',4,psa_f);
          
          % Compute expected number of impacted households (4 people/house) [1 Nm].
          Nn2=sum(On2.*pop(In,:),1)/4; Nn3=sum(On3.*pop(In,:),1)/4; Nn4=sum(On4.*pop(In,:),1)/4; Nn5=sum(On5.*pop(In,:),1)/4; Nn6=sum(On6.*pop(In,:),1)/4;
          Nd1=sum(Od1.*pop(Id,:),1)/4; Nd2=sum(Od2.*pop(Id,:),1)/4; Nd3=sum(Od3.*pop(Id,:),1)/4; Nd4=sum(Od4.*pop(Id,:),1)/4;
          
          % Stash results into the output data structure.
          S.RISK(j).Nn2=[S.RISK(j).Nn2;Nn2]; S.RISK(j).Nn3=[S.RISK(j).Nn3;Nn3]; S.RISK(j).Nn4=[S.RISK(j).Nn4;Nn4]; S.RISK(j).Nn5=[S.RISK(j).Nn5;Nn5]; S.RISK(j).Nn6=[S.RISK(j).Nn6;Nn6];
          S.RISK(j).Nd1=[S.RISK(j).Nd1;Nd1]; S.RISK(j).Nd2=[S.RISK(j).Nd2;Nd2]; S.RISK(j).Nd3=[S.RISK(j).Nd3;Nd3]; S.RISK(j).Nd4=[S.RISK(j).Nd4;Nd4];
          
      end
  end
  
  % Flip the risk routine run flag.
  S.dVAR.UPDATEflag='no';
  
return
