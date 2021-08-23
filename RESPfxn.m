function P=RESPfxn(GM,dE,type_flag,level_flag,psa_f)
  % Computes the chance of observation for nuisance (PGV: m/s) or damage (PGA: g).
  %
  % References:
  % 
  % Chase, R. E., Liel, A. B., Luco, N., & Baird, B. W. (2019). Seismic loss and damage in light‐frame wood buildings from sequences of induced earthquakes. Earthquake Engineering & Structural Dynamics, 48(12), 1365-1383, doi: 10.1002/eqe.3189.
  % Schultz, Quitoriano, Wald, & Beroza (2020). Quantifying nuisance ground motion thresholds for induced earthquakes, Earthquake Spectra, doi: xx.
  % FEMA (2020). Hazus-MH 2.1. Multi-Hazard Loss Estimation Methodology Earthquake Model Technical Manual. Federal Emergency Management Agency.
  %
  % Written by Ryan Schultz.
  
  if(strcmpi(type_flag,'nuisance')) % Schultz et al., 2020 (Table 1).
      
      if(level_flag==2) % Just felt.
          B=[11.2336 3.8771];
          dB=[+6.05e-1 +4.47e-2 +1.63e-1];
      elseif(level_flag==3) % Excitement.
          B=[ 5.9376 2.5875];
          dB=[+2.48e-1 +1.84e-2 +6.68e-2];
      elseif(level_flag==4) % Somewhat frightened.
          B=[ 3.4387 2.1331];
          dB=[+8.10e-2 +4.23e-3 +1.78e-2];
      elseif(level_flag==5) % Very frigthened.
          B=[ 2.3712 2.0828];
          dB=[+2.37e-2 +6.63e-4 +3.16e-3];
      elseif(level_flag==6) % Extremely frightened.
          B=[ 2.1427 2.1260];
          dB=[+2.12e-2 +5.18e-4 +2.42e-3];
      end
      
      % Perturb the coeffcients.
      B(1)=B(1)+dE(1)*sqrt(dB(1));
      B(2)=B(2)+dE(2)*sqrt(dB(2));
      B(2)=B(2)+dE(1)*(dB(3)/sqrt(dB(1)));
      
      % Compute the probabiltiies.
      P=1./(1+exp(-(B(1)+log10(GM)*B(2))));
      
  elseif(strcmpi(type_flag,'damage'))
      
      if(psa_f==0.0) % FEMA, 2020 (Table 5.13c; low-code, average between all building types).
          if(level_flag==1) % Slight/minor.
              u= [0.20, 0.653];
              du=[0.00, 0.014];
          elseif(level_flag==2) % Moderate.
              u= [0.40, 0.672];
              du=[0.00, 0.009];
          elseif(level_flag==3) % Extensive.
              u= [0.80, 0.668];
              du=[0.00, 0.014];
          elseif(level_flag==4) % Complete.
              u= [1.60, 0.668];
              du=[0.00, 0.014];
          end
      elseif(psa_f==0.45) % Chase et al., 2019 (Figure 14).
          % Add this later?
          if(level_flag==1) % Slight/minor.
              u= [0.1580, 0.160];
              du=[0.00, 0.014];
          elseif(level_flag==2) % Moderate.
              u= [0.5475, 0.170];
              du=[0.00, 0.009];
          elseif(level_flag==3) % Extensive.
              u= [0.9390, 0.212];
              du=[0.00, 0.014];
          elseif(level_flag==4) % Complete.
              u= [1.1685, 0.213];
              du=[0.00, 0.014];
          end
      end
      
      % Perturb the coeffcients.
      u=u+dE.*du;
      
      % Compute the probabilites.
      P=normcdf(log(GM), log(u(1)), u(2));
      
      % Make the fragility function fall to zero faster?
      %P(P<1e-2)=P(P<1e-2).^(-(1/2)*log10(P(P<1e-2)));
      
  end
  
  
  
return