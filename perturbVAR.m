function S=perturbVAR(S,N,rand_flag)
  % Simple function that creates a list of perturbed variables.
  % 
  % Writtern by Ryan Schultz.
  
  % Define the b-value statistics.
  bm=1.00;
  db=0.05;
  
  % Flag for random or average behaviour.
  if(strcmpi(rand_flag,'random'))
      % Create a vector of perturbation values.
      b=normrnd(bm,db,[1 N]);
      q=rand([1 N]).*(1-exp(-(1.2.^b)))+exp(-(1.2.^b));
      dM=log10(1.2)./b-log10(-log(q))./b;
      dGM=normrnd(0.0,1.0,[1 N]);
      dF1=normrnd(0.0,1.0,[1 N]);
      dF2=normrnd(0.0,1.0,[1 N]);
      dSA=lognrnd(0.0,0.05,[1 N]);
      dPOP=abs(normrnd(mean(S.MAP.POP(:)), sqrt(mean(S.MAP.POP(:))), [1 N])/mean(S.MAP.POP(:)));
      dZ=pearsrnd(0+0.15,0.25,1.05,4,[1 N]);
  elseif(strcmpi(rand_flag,'none'))
      % Create a vector of perturbation values.
      b=bm*ones([1 N]);
      dM=zeros([1 N]);
      dGM=zeros([1 N]);
      dF1=zeros([1 N]);
      dF2=zeros([1 N]);
      dSA=ones([1 N]);
      dPOP=ones([1 N]);
      dZ=zeros([1 N]);
  end
  
  % Append new values to the structure.
  S.dVAR.b=[S.dVAR.b,b];
  S.dVAR.dM=[S.dVAR.dM,dM];
  S.dVAR.dGM=[S.dVAR.dGM,dGM];
  S.dVAR.dF1=[S.dVAR.dF1,dF1];
  S.dVAR.dF2=[S.dVAR.dF2,dF2];
  S.dVAR.dSA=[S.dVAR.dSA,dSA];
  S.dVAR.dPOP=[S.dVAR.dPOP,dPOP];
  S.dVAR.dZ=[S.dVAR.dZ,dZ];
  
  % Flag that the risk routine needs to be (re)run now.
  S.dVAR.UPDATEflag='yes';
  
return