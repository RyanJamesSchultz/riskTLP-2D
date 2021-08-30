% Simple script to plot the risk-informed TLP data structure.
clear;

% Load in the data structure: S.
load('Rmap.mat');

% Predefine some variables.
GREY=[0.85,0.85,0.85];
lat=linspace(min(S.MAP.latB),max(S.MAP.latB),7); lat(2:end-1);
lon=linspace(min(S.MAP.lonB),max(S.MAP.lonB),7); lon(2:end-1);
[~,m]=min(Geoid_Distance(lat(1),lon(1),[S.RISK.lat],[S.RISK.lon],'spherical'));
[~,i]=min(Geoid_Distance(lat(2),lon(2),[S.RISK.lat],[S.RISK.lon],'spherical'));
[~,j]=min(Geoid_Distance(lat(3),lon(3),[S.RISK.lat],[S.RISK.lon],'spherical'));
[~,k]=min(Geoid_Distance(lat(4),lon(4),[S.RISK.lat],[S.RISK.lon],'spherical'));
[~,l]=min(Geoid_Distance(lat(5),lon(5),[S.RISK.lat],[S.RISK.lon],'spherical'));
Nv=length(S.dVAR.b);
yl=[7e-2 7e+5];
ML_c=[2.0 5.0];
ML_f=3.5;
Pd=50;
Pn=50;
Ni=12;

% Get the play-specific thresholds to use.
if(strcmpi(S.play_flag,'SYN'))
    Nn3_f=1e4;
    Nd1_f=1e2;
end

% Plot all map data.
figure(1); clf;
subplot(411);
contourf(S.MAP.lonG,S.MAP.latG,S.MAP.Vs30,'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
plot(S.RISK(m).lon,S.RISK(m).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(i).lon,S.RISK(i).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(j).lon,S.RISK(j).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(k).lon,S.RISK(k).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(l).lon,S.RISK(l).lat,'wo','MarkerFaceColor','k');
xlabel('Longitude'); ylabel('Latitude'); title('Site Amplification');
h = colorbar(); colormap(gca,R_colormap('Vs30')); ylabel(h, 'Site Amplififcation, Vs30 (m/s)'); hold off;
subplot(412);
contourf(S.MAP.lonE,S.MAP.latE,S.MAP.DEP,'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
plot(S.RISK(m).lon,S.RISK(m).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(i).lon,S.RISK(i).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(j).lon,S.RISK(j).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(k).lon,S.RISK(k).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(l).lon,S.RISK(l).lat,'wo','MarkerFaceColor','k');
xlabel('Longitude'); ylabel('Latitude'); title('Formation Depth');
h = colorbar(); colormap(gca,R_colormap('Depth')); ylabel(h, 'True Vertical Depth (km)'); hold off;
subplot(413);
contourf(S.MAP.lonG,S.MAP.latG,log10(S.MAP.POP),'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-w');
plot(S.RISK(m).lon,S.RISK(m).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(i).lon,S.RISK(i).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(j).lon,S.RISK(j).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(k).lon,S.RISK(k).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(l).lon,S.RISK(l).lat,'wo','MarkerFaceColor','k');
set(gca,'Color','k');
xlabel('Longitude'); ylabel('Latitude'); title(['Population (',sprintf('%0.3g',sum(sum(S.MAP.POP))),')']);
h = colorbar(); colormap(gca,R_colormap('population')); ylabel(h, 'Population (log_{10}[people])'); hold off;
subplot(414);
contourf(S.MAP.lonE,S.MAP.latE,S.MAP.Ir,'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
plot(S.RISK(m).lon,S.RISK(m).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(i).lon,S.RISK(i).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(j).lon,S.RISK(j).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(k).lon,S.RISK(k).lat,'wo','MarkerFaceColor','k');
plot(S.RISK(l).lon,S.RISK(l).lat,'wo','MarkerFaceColor','k');
xlabel('Longitude'); ylabel('Latitude'); colorbar(); colormap(gca,'default'); title('In Bounds'); hold off;

% Plot all perturbation data.
figure(2); clf;
subplot(331); histogram(S.dVAR.dZ, round(sqrt(Nv)) );
xlabel('Depth Perturbation, dZ (km)'); ylabel('Count');
subplot(332); histogram(S.dVAR.b, round(sqrt(Nv)) );
xlabel('b-value distribution, b (-)'); ylabel('Count');
subplot(333); histogram(S.dVAR.dM, round(sqrt(Nv)) );
xlabel('Magnitude Perturbation, dM (Mw)'); ylabel('Count');
set(gca, 'YScale', 'log');
subplot(334); histogram(S.dVAR.dGM, round(sqrt(Nv)) );
xlabel('GMPE Perturbation, dGM (-)'); ylabel('Count');
subplot(335); histogram(S.dVAR.dF1, round(sqrt(Nv)) );
xlabel('RF Perturbation, dF1 (-)'); ylabel('Count');
subplot(336); histogram(S.dVAR.dF2, round(sqrt(Nv)) );
xlabel('RF Perturbation, dF2 (-)'); ylabel('Count');
subplot(337); histogram(S.dVAR.dSA, round(sqrt(Nv)) );
xlabel('Site Amp Perturbation, dSA (-)'); ylabel('Count');
subplot(338); histogram(S.dVAR.dPOP, round(sqrt(Nv)) );
xlabel('Population Perturbation Factor, dPOP (-)'); ylabel('Count');

% Plot some of the count curves.
figure(3); clf;
subplot(621);
semilogy(S.Mw,S.RISK(i).Nn2,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nn2,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nn2,1),'-b');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Nuisance: CDI 2');
ylim(yl);
subplot(622);
semilogy(S.Mw,S.RISK(i).Nd1,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nd1,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nd1,1),'-b');
line(xlim(),Nd1_f*[1 1],'Color','k','LineStyle','--');
line(ML_f*[1 1],ylim(), 'Color','k','LineStyle','--');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Damage: DS 1');
ylim(yl);
subplot(623);
semilogy(S.Mw,S.RISK(i).Nn3,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nn3,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nn3,1),'-b');
line(xlim(),Nn3_f*[1 1],'Color','k','LineStyle','--');
line(ML_f*[1 1],ylim(), 'Color','k','LineStyle','--');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Nuisance: CDI 3');
ylim(yl);
subplot(624);
semilogy(S.Mw,S.RISK(i).Nd2,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nd2,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nd2,1),'-b');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Damage: DS 2');
ylim(yl);
subplot(625);
semilogy(S.Mw,S.RISK(i).Nn4,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nn4,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nn4,1),'-b');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Nuisance: CDI 4');
ylim(yl);
subplot(626);
semilogy(S.Mw,S.RISK(i).Nd3,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nd3,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nd3,1),'-b');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Damage: DS 3');
ylim(yl);
subplot(627);
semilogy(S.Mw,S.RISK(i).Nn5,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nn6,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nn6,1),'-b');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Nuisance: CDI 6');
ylim(yl);
subplot(628);
semilogy(S.Mw,S.RISK(i).Nd4,'-','color',GREY); hold on;
semilogy(S.Mw,mean(S.RISK(i).Nd4,1),'-r');
semilogy(S.Mw,median(S.RISK(i).Nd4,1),'-b');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Damage: DS 4');
ylim(yl);
subplot(6,2,[9 11]);
semilogy(S.Mw,median(S.RISK(i).Nn2,1),'DisplayName','CDI 2'); hold on;
semilogy(S.Mw,median(S.RISK(i).Nn3,1),'DisplayName','CDI 3');
semilogy(S.Mw,median(S.RISK(i).Nn4,1),'DisplayName','CDI 4');
semilogy(S.Mw,median(S.RISK(i).Nn5,1),'DisplayName','CDI 5');
semilogy(S.Mw,median(S.RISK(i).Nn6,1),'DisplayName','CDI 6');
h=line(xlim(),Nn3_f*[1 1],'Color',GREY,'LineStyle','--'); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h=line(ML_f*[1 1],ylim(), 'Color',GREY,'LineStyle','--'); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Nuisance'); legend('Location','southeast');
%ylim(yl);
subplot(6,2,[10 12]);
semilogy(S.Mw,median(S.RISK(i).Nd1,1),'DisplayName','DS 1'); hold on;
semilogy(S.Mw,median(S.RISK(i).Nd2,1),'DisplayName','DS 2');
semilogy(S.Mw,median(S.RISK(i).Nd3,1),'DisplayName','DS 3');
semilogy(S.Mw,median(S.RISK(i).Nd4,1),'DisplayName','DS 4');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Damage'); legend('Location','southeast');
%ylim(yl);

% Plot the change in curves with proximity to the 'city.'
figure(4); clf;
subplot(121);
semilogy(S.Mw,median(S.RISK(m).Nn3,1),'-k','DisplayName','Farthest'); hold on;
semilogy(S.Mw,median(S.RISK(i).Nn3,1),'-b','DisplayName','Far');
semilogy(S.Mw,median(S.RISK(j).Nn3,1),'--b','DisplayName','Medium');
semilogy(S.Mw,median(S.RISK(k).Nn3,1),':b','DisplayName','Close');
semilogy(S.Mw,median(S.RISK(l).Nn3,1),':r','DisplayName','On-top');
h=line(xlim(),Nn3_f*[1 1],'Color',GREY,'LineStyle','--'); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h=line(ML_f*[1 1],ylim(), 'Color',GREY,'LineStyle','--'); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Nuisance (CDI 3)'); legend('Location','southeast');
subplot(122);
semilogy(S.Mw,median(S.RISK(m).Nd2,1),'-k','DisplayName','Farthest'); hold on;
semilogy(S.Mw,median(S.RISK(i).Nd2,1),'-b','DisplayName','Far');
semilogy(S.Mw,median(S.RISK(j).Nd2,1),'--b','DisplayName','Medium');
semilogy(S.Mw,median(S.RISK(k).Nd2,1),':b','DisplayName','Close');
semilogy(S.Mw,median(S.RISK(l).Nd2,1),':r','DisplayName','On-top');
xlabel('Red-Light Magnitdue (M_L)'); ylabel('Number of Households Impacted'); title('Damage (DS 1)'); legend('Location','southeast');

% Scatter plot nuisance vs damage counts.
figure(5); clf;
loglog(median(S.RISK(i).Nd1,1),median(S.RISK(i).Nn2,1),'o','DisplayName','CDI2-DS1'); hold on;
loglog(median(S.RISK(i).Nd2,1),median(S.RISK(i).Nn3,1),'o','DisplayName','CDI3-DS2');
loglog(median(S.RISK(i).Nd3,1),median(S.RISK(i).Nn4,1),'o','DisplayName','CDI4-DS3');
loglog(median(S.RISK(i).Nd4,1),median(S.RISK(i).Nn6,1),'o','DisplayName','CDI6-DS4');
ylabel('Nuisance Impact'); xlabel('Damage Impact'); legend('Location','southeast');

% Make iso-risk maps.
Rn=mapRISK(S,'nuisance',3,1,Nn3_f,ML_f,Pn,50,Ni);
Rd=mapRISK(S,'damage',  3,1,Nd1_f,ML_f,50,Pd,Ni);

%Rd.Nd(Rd.Nd>200)=200;
%Rn.Nd(Rn.Nd>200)=200;

% Make a combined map.
Mr2=min(cat(3,Rn.Mr,Rd.Mr),[],3,'omitnan');

% Get boundaries of play area.
YL=[min(S.MAP.latB)-0.1 max(S.MAP.latB)+0.1];
XL=[min(S.MAP.lonB)-0.1 max(S.MAP.lonB)+0.1];

% Plot iso-maps
figure(6); clf;
subplot(251);
contourf(Rn.lon,Rn.lat,Rn.Mr,linspace(ML_c(1),ML_c(2),35),'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title(['Iso-Nuisance Map, Nn3=',num2str(Nn3_f)]);
h = colorbar(); ylabel(h, 'Magnitude (M_L)');
colormap(gca,R_colormap('red-light')); caxis(ML_c);
xlim(XL); ylim(YL);
subplot(252);
contourf(Rn.lon,Rn.lat,Rn.Nn,'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title(['Iso-Magnitude Map, M_L=',num2str(ML_f)]);
h = colorbar(); ylabel(h, 'Impacted Household Count');
colormap(gca,R_colormap('nuisance'));
xlim(XL); ylim(YL);
subplot(253);
contourf(Rn.lon,Rn.lat,Rn.Nd,'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title('Equivalent Damage Impact Map');
h = colorbar(); ylabel(h, 'Impacted Household Count');
colormap(gca,R_colormap('damage')); %caxis([0 100]);
xlim(XL); ylim(YL);
subplot(256);
contourf(Rd.lon,Rd.lat,Rd.Mr,linspace(ML_c(1),ML_c(2),15),'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title(['Iso-Damage Map, Nd1=',num2str(Nd1_f)]);
h = colorbar(); ylabel(h, 'Magnitude (M_L)');
colormap(gca,R_colormap('red-light')); caxis(ML_c);
xlim(XL); ylim(YL);
subplot(257);
contourf(Rd.lon,Rd.lat,Rd.Nd,'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title(['Iso Magnitude Map, M_L=',num2str(ML_f)]);
h = colorbar(); ylabel(h, 'Impacted Household Count');
colormap(gca,R_colormap('damage')); %caxis([0 100]);
xlim(XL); ylim(YL);
subplot(258);
contourf(Rd.lon,Rd.lat,Rd.Nn,'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title('Equivalent Nuisance Impact Map');
h = colorbar(); ylabel(h, 'Impacted Household Count');
colormap(gca,R_colormap('nuisance'));
xlim(XL); ylim(YL);

% Plot the combination map.
%figure(7); clf;
subplot(2,5,[4 5 9 10])
contourf(Rn.lon,Rn.lat,Mr2,linspace(ML_c(1),ML_c(2),15),'LineColor','none'); hold on;
plot(S.MAP.lonB,S.MAP.latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title('Combination Map');
h = colorbar(); ylabel(h, 'Magnitude (M_L)');
colormap(gca,R_colormap('red-light')); caxis(ML_c);
xlim(XL); ylim(YL);

% Plot some histograms of the iso-magnitude impacts.
figure(7); clf;
subplot(121);
histogram(Rn.Nn(~isnan(Rn.Nn)), round(sqrt(length(Rn.Nn(~isnan(Rn.Nn))))) );
xlabel('Iso-Magnitude Nuisance Impacts'); ylabel('Counts');
set(gca, 'YScale', 'log');
subplot(122);
histogram(Rd.Nd(~isnan(Rd.Nd)), round(sqrt(length(Rd.Nd(~isnan(Rd.Nd))))) );
xlabel('Iso-Magnitude Damage Impacts'); ylabel('Counts');
set(gca, 'YScale', 'log');

% Output for DGSA.
%parameters=[S.dVAR.dM',S.dVAR.b',S.dVAR.dGM',S.dVAR.dF1',S.dVAR.dF2',S.dVAR.dSA',S.dVAR.dPOP',S.dVAR.dZ'];
%responsesN=[]; responsesD=[];
%for o=1:length(Rn.MAPs)
%    responsesN=[responsesN; Rn.MAPs(o).Mr(S.MAP.Ir)'];
%    responsesD=[responsesD; Rd.MAPs(o).Mr(S.MAP.Ir)'];
%end

