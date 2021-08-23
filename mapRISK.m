function [R]=mapRISK(S,TypeFlag,CDI,DS,N_c,ML_c,Pn,Pd,Ni)
  % Compute the risk maps.
  % 
  % Writtern by Ryan Schultz.
  
  % Set the interpolation type.
  interp_type='linear';
  
  % Check that we're up to date.
  if(strcmpi(S.dVAR.UPDATEflag,'yes'))
      return;
  end
  
  % Define structure and preallocate space for output maps.
  R=struct('MAPs',[]);
  MAPs=struct('Mr',[],'Nn',[],'Nd',[], 'Mavg',[]);
  
  % Predefine some important varaibles.
  J=find(S.MAP.Ir);
  M=S.Mw;
  
  % Loop over all of the map realizations.
  for i=1:length(S.dVAR.dM)
      
      Mr=NaN*zeros(size(S.MAP.Ir));
      Nn=Mr;
      Nd=Mr;
      
      % Loop over all of the (in-bounds) EQ grid.
      for j=1:length(S.RISK)
          
          % Get the nuisance level of interest.
          if(CDI==2)
              N_n=S.RISK(j).Nn2(i,:);
          elseif(CDI==3)
              N_n=S.RISK(j).Nn3(i,:);
          elseif(CDI==4)
              N_n=S.RISK(j).Nn4(i,:);
          elseif(CDI==5)
              N_n=S.RISK(j).Nn5(i,:);
          elseif(CDI==6)
              N_n=S.RISK(j).Nn6(i,:);
          end
          
          % Get the damage level of interest.
          if(DS==1)
              N_d=S.RISK(j).Nd1(i,:);
          elseif(DS==2)
              N_d=S.RISK(j).Nd2(i,:);
          elseif(DS==3)
              N_d=S.RISK(j).Nd3(i,:);
          elseif(DS==4)
              N_d=S.RISK(j).Nd4(i,:);
          end
          
          % Dealing with non-monotonic and non-unique inputs.
          dx=cumsum(ones(size(M)))*eps;
          N_n=N_n+dx;
          N_d=N_d+dx;
          
          %figure(10); plot(M,N_n);
          
          % Find intersecting values.
          if(strcmpi(TypeFlag,'nuisance'))
              Mr(J(j))=interp1(N_n,M,N_c,interp_type,'extrap');
              Mr(J(j))=min([6.0 Mr(J(j))]);
              Nn(J(j))=interp1(M,N_n,ML_c,interp_type,'extrap');
              Nd(J(j))=interp1(M,N_d,Mr(J(j)),interp_type,'extrap');
          else
              Mr(J(j))=interp1(N_d,M,N_c,interp_type,'extrap');
              Mr(J(j))=min([6.0 Mr(J(j))]);
              Nd(J(j))=interp1(M,N_d,ML_c,interp_type,'extrap');
              Nn(J(j))=interp1(M,N_n,Mr(J(j)),interp_type,'extrap');
          end
      end
      
      % Error handling.
      Nn(Nn<0)=0;
      Nd(Nd<0)=0;
      
      % Find the average mapped magnitude for this realization.
      MAPs(i).Mavg=nanmedian(Mr,'all');
      
      % Stuff map into the output structure.
      MAPs(i).Mr=Mr;
      MAPs(i).Nn=Nn;
      MAPs(i).Nd=Nd;
      
  end
  
  % Compute the averaged maps.
  if(strcmpi(TypeFlag,'nuisance'))
      Mr=prctile(cat(3,MAPs.Mr),Pn,3);
      Nn=prctile(cat(3,MAPs.Nn),Pn,3);
      Nd=prctile(cat(3,MAPs.Nd),Pd,3);
  else
      Mr=prctile(cat(3,MAPs.Mr),Pd,3);
      Nn=prctile(cat(3,MAPs.Nn),Pn,3);
      Nd=prctile(cat(3,MAPs.Nd),Pd,3);
  end
  
  % Get the equivalent impact map.
  % Loop over all of the (in-bounds) EQ grid.
  for j=1:length(S.RISK)
      
      % Get the impact of interest.
      if(strcmpi(TypeFlag,'nuisance'))
          % Get the damage impact of interest.
          if(DS==1)
              N_e=prctile([S.RISK(j).Nd1],Pd,1);
          elseif(DS==2)
              N_e=prctile([S.RISK(j).Nd2],Pd,1);
          elseif(DS==3)
              N_e=prctile([S.RISK(j).Nd3],Pd,1);
          elseif(DS==4)
              N_e=prctile([S.RISK(j).Nd4],Pd,1);
          end
          
      else
          % Get the nuisance impact of interest.
          if(CDI==2)
              N_e=prctile([S.RISK(j).Nn2],Pn,1);
          elseif(CDI==3)
              N_e=prctile([S.RISK(j).Nn3],Pn,1);
          elseif(CDI==4)
              N_e=prctile([S.RISK(j).Nn4],Pn,1);
          elseif(CDI==5)
              N_e=prctile([S.RISK(j).Nn5],Pn,1);
          elseif(CDI==6)
              N_e=prctile([S.RISK(j).Nn6],Pn,1);
          end
      end
      
      % Find intersecting values.
      if(strcmpi(TypeFlag,'nuisance'))
          Nd(J(j))=interp1(M,N_e,Mr(J(j)),interp_type,'extrap');
      else
          Nn(J(j))=interp1(M,N_e,Mr(J(j)),interp_type,'extrap');
      end
  end
  
  % Upsample the averaged maps.
  lon=linspace(min(S.MAP.lonE),max(S.MAP.lonE),Ni*length(S.MAP.lonE));
  lat=linspace(min(S.MAP.latE),max(S.MAP.latE),Ni*length(S.MAP.latE));
  [LON,LAT]=meshgrid(lon,lat);
  [LONe,LATe]=meshgrid(S.MAP.lonE,S.MAP.latE);
  Mr=scatteredInterpolant(LONe(J),LATe(J),Mr(J),'natural','linear'); Mr=Mr(LON,LAT);
  Nn=scatteredInterpolant(LONe(J),LATe(J),Nn(J),'natural','linear'); Nn=Nn(LON,LAT);
  Nd=scatteredInterpolant(LONe(J),LATe(J),Nd(J),'natural','linear'); Nd=Nd(LON,LAT);
  I=inpolygon(LON,LAT,S.MAP.lonB,S.MAP.latB);
  Mr(~I)=NaN;
  Nn(~I)=NaN;
  Nd(~I)=NaN;
  
  % Error handling.
  Nn(Nn<0)=0;
  Nd(Nd<0)=0;
  
  % Stuff results into the output structure.
  R.MAPs=MAPs;
  R.Mr=Mr;
  R.Nn=Nn;
  R.Nd=Nd;
  R.lat=lat;
  R.lon=lon;
  
return