clear;
close all;
clc;

% This script is to perform an extensive parameter search for Martian
% lithosphere modelling
%
% The code is written by Bart Root
% 
% Using code by Zdenek Martinec and Nicolas Tosi: SFEC
% And GSH_package by Bart Root
% topo2crust: by Weilun Qin
% create_tomfiles: by Maxime Vincent and Bart Root

HOMEDIR = pwd;

% load dependencies
SFEC_dir = '<directory where you placed the SFEC code>';       % add your directory of the github repositories 
GSH_dir = '<directory where you placed the GSH code>/Tools/';   %        add your directory of the github repositories 

addpath(GSH_dir)
addpath('Data/')
addpath(SFEC_dir)

%% Input parameters
Model = struct();

% Gravity modelling input values
lmax        = 179;

Model.number_of_layers = 3;
Model.G     = 6.672e-11;    
Model.GM    = 4.2828372e13;    % Mars
Model.Re    = 3396000.0;       % Mars
Model.name  = 'Thin Shell'; %'Thin Shell' or 'Infinite Plate' or 'Airy'
Model.nmax  = lmax;     
Model.geoid = 'none';
Model.sfec_dir = SFEC_dir;


%% Topographic data

[topo,LonT,LatT] = gmt2matrix(load('Data/Topography_Mars_1.txt'));
topo = topo.*1e3;
base = gmt2matrix(load('Data/Basement_Mars_1.txt')).*1e3;

%% Gravity data 

V = load('Data/Mars_MRO120d.txt');

% Mapping input values
latLim   = [-89.95 89.95 .1];      %% °         %% limits of the lat lon grid
lonLim   = [0.05 359.95 .1];  

% Mapping input values 
height   = 0.0;
SHbounds = [2 90];
verbose  = 1;

V(3,3) = 0; % removing the mean gravity field

% calculating the observed grvaity field
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
lon = data.grd.lon;
lats = data.grd.lat; 

%%
tic;
runid = 1;
%variations = 2600:100:3300;
%num_var = length(variations);
%DV_mr = zeros(180,num_var);
for DC  =  60%50:5:65
for RHOC = 3050%2800:50:3200
for DRC = 500%100:100:800
for TE = 104.0042 %[ 0.010 48.2745   82.5482  104.0042  131.0371  150.0000  177.8447  208.0084  224.0702].*1e3
for DEPTH = 1000 %400:200:1600
for RADIUS = 1500 % 100:50:700
for DRM = -40 %-250:50:50
for THICK = 250 %50:50:400

% Crustal input values
Model.D_c   = DC*1e3;                   % 30-100 Average crustal thickness.
Model.rho_c = RHOC;                     % 2600-3300          crustal density [kg/m^3]
drho        = DRC;                      % 100-1000 [kg/m^3]
Model.rho_m = Model.rho_c+drho;         % mantle density [kg/m^3]
Model.core_depth = Model.Re-1.835e6;
Model.rho_i = 1220;                     % 1220 kg/m3from Zuber et al (2007)

% Flexure input values
Model.Te = TE.*1e3;                     % for Rigidity 48.2745   82.5482  104.0042  131.0371  150.0000  177.8447  208.0084  224.0702
Model.E  = 100e9;
Model.v  = 0.30;

% SFEC density model input variables
Model.mantle.shape = 'Disk'; 			 % Sphere, Disk, Pillar, (TBD: Plume)
Model.mantle.deltarho   = DRM;          %% kg/m^3    	%% set the value of the density anomaly 
Model.mantle.layer_size = 25;            %% km        	%% set the thickness of each layer 
Model.mantle.depth_anom = DEPTH;           %% km        	%% set the depth of the anomaly 
Model.mantle.center     = ([250,87]);    %% °         	%% set the center of the anomaly -- index of longitude  & co-latitude
Model.mantle.kmradius   = RADIUS;          %% km        	%% set the radius in kilometers of the anomaly 
Model.mantle.thickness   = THICK;          %% km        	%% set the thickness in kilometers of the disk 

% SFEC mantle convection input variabls

Model.sfec.min_degree  = 2;
Model.sfec.max_degree  = 10;
Model.sfec.LIV         = 1e23;
Model.sfec.UMV         = 6e20;
Model.sfec.LMV         = 1e21; 

%% SFEC Mantle Convection modelling
if Model.mantle.deltarho==0
    disp('No dynamic topography')
else
    % Construct tomography file using Create_tom.m
    disp('create tomography files for sfec...')
    create_tom_files(Model);    
    
    % outfile is outputmax.tom to be used in the sfec code
    
    %% initialise the SFEC input and modules file
    
    % SFEC.inp initialising
    disp('Create sfec input file...')
    construct_SFEC_inp(Model);
    
    % modules.f90 initialising
    disp('Construct geometry file for sfec...')
    construct_geometry_mod(Model);
    
    %% run SFEC
    disp('Compile and run sfec...')
    
    % compile SFEC       
    cd(SFEC_dir)
    [status,cmdout] = system('make');
    
    % use Matlab interface to run SFEC
    if status == 0
    [status,cmdout_sfec] = system('./sfec');
    else
    error('Compiling of SFEC went incorrect')
    end
    
    % test if file exists
    if isfile('geoid.dat')
     disp('SFEC output file exists.')
    else
     error(['Running sfec encounter errors. See: /n--------------------------/n' cmdout_sfec])
    end
    
    cd(HOMEDIR)
    
    %% load geoid file to extract dynamic topography and gravity files
    
    % load data from the sfec output
    G = load([SFEC_dir 'geoid_coeff.dat']);
    maxdeg = max(G(:,1));
    
    [VGM]  = zdenek2V(G(:,[1:2 9:10]),maxdeg);
    [VDT]  = zdenek2V(G(:,[1:2 5:6] ),maxdeg);
    [VDC]  = zdenek2V(G(:,[1:2 7:8] ),maxdeg);
    [VGTo] = zdenek2V(G(:,[1:2 3:4] ),maxdeg);
    
    %% Calculate the degree variance of the gravity signal of cmb and topo
    
    n = VGM(:,1);
    fatt = 4*pi*Model.G./(2*n+1);
    
    VGC = VDC;
    VGC(:,3) = fatt.* 1.835d6.*(1.835d6./Model.Re).^(n+1).*-3500.*VDC(:,3);
    VGC(:,4) = fatt.* 1.835d6.*(1.835d6./Model.Re).^(n+1).*-3500.*VDC(:,4);
    
    VGMC = VGM;
    VGMC(:,3) = (VGMC(:,3) + VGC(:,3)).*Model.Re/Model.GM;
    VGMC(:,4) = (VGMC(:,4) + VGC(:,4)).*Model.Re/Model.GM;
    
    % Dynamic topography
    
    sc_dyn_topo = vecml2sc(VDT(:,3),VDT(:,4),maxdeg);
    dyn_topo = GSHS(sc_dyn_topo,LonT(1,:),90-LatT(:,1),maxdeg);    
    
    % dynamic_gravity interpolate on specified grid
    
    [data_dyn_Mantle] = model_SH_synthesis(lonLim,latLim,height,[Model.sfec.min_degree Model.sfec.max_degree],VGMC,Model);
end

%% correct topography for equivalent topography
ice_thick = topo-base;
topo_equi = base + ice_thick.*Model.rho_i./Model.rho_c;

if Model.mantle.deltarho==0
    el = topo_equi;
else
    %Correct the topography for dynamic topography
    el = topo_equi - dyn_topo;   
end

%% modelling a local isostatic lithosphere

if strcmp(Model.name,'Airy')
    disp('Construct Airy model layer..')
    [airy_CM,lon_CM,lat_CM] = topo2crust(el,lmax,'Airy',Model);
    
    if verbose==1
    
       figure
       imagesc(lon_CM,lat_CM,airy_CM)
       set(gca,'YDir','normal')
    
    end
    moho = -airy_CM;

elseif strcmp(Model.name,'Infinite Plate')

    disp('Construct Infinite Plate model layer..')
    [IP_CM,lon_CM,lat_CM] = topo2crust(el,lmax,'Infinite_Plate',Model);
    
    if verbose==1
    
       figure
       imagesc(lon_CM,lat_CM,IP_CM)
       set(gca,'YDir','normal')
    
    end
    moho = -IP_CM;

elseif strcmp(Model.name,'Thin Shell')
    disp('Construct Thin Shell model layer..')
    [TS_CM,lon_CM,lat_CM] = topo2crust(el,lmax,'Thin_Shell',Model);
    moho = -TS_CM;

    if verbose==1
       %load("batlowW.mat")
       figure
       imagesc(lon_CM,lat_CM,-TS_CM./1e3);colormap(turbo);cc = colorbar;
       xlabel('Longitude [\circ]','FontSize',20)
       ylabel('Latitude [\circ]','FontSize',20)
       title('crust-mantle boundary','FontSize',20)
       ylabel(cc,'km','FontSize',20)
       set(gca,'YDir','normal','FontSize',20)

        figure
        axesm('MapProjection','ortho','origin',[-90,0]);c=colorbar;colormap(turbo) 
        framem
        %load coast;
        %plotm(lat,long,'k')
        gridm
        %[gLat,gLon] = meshgrid(lat,lon);
        surfm(data.grd.lat,data.grd.lon,-moho./1e3)
        xlabel('Longitude [\circ]','FontSize',20)
        ylabel('Latitude [\circ]','FontSize',20)
        title(['North pole'],'FontSize',20)
        ylabel(c,'mGal','FontSize',20)
        caxis([30 80])
        set(gca,'YDir','normal','FontSize',20)
    end
    
end

%% Hellas basin check and Insight value

minTopo = min(min(topo));
maxMoho = max(max(moho));

if minTopo-maxMoho < 0
   disp(['WARNING Hellas Basin has negative volume! Value is ' num2str(minTopo-maxMoho) ' km']) 
end

Finsight = scatteredInterpolant(reshape(LonT,[],1),reshape(LatT,[],1),reshape(moho,[],1));
Model.crust_insight = Finsight(224.1,4.5)./1e3;


%% Gravity modelling of the different crustal structures.

top = 25000;
bot = min(min(moho));
thick_lay = 25000;

top_layer = top:-thick_lay:bot;
bot_layer = [top_layer(2:end) bot];

for numl = 1:length(top_layer)
    disp(['Constructing layer number ' num2str(numl) '...'])

    ubound = top_layer(numl);
    lbound = bot_layer(numl);

    % LAB is shallower than 100 km
    upper_LAY = topo; 
    upper_LAY(topo>ubound) = ubound;                                                                                                                                       
    upper_LAY(topo<lbound) = lbound;

    mid_LAY = base; 
    mid_LAY(base>ubound) = ubound;                                                                                                                                       
    mid_LAY(base<lbound) = lbound; 

    low_LAY = moho; 
    low_LAY(moho>ubound) = ubound;                                                                                                                                       
    low_LAY(moho<lbound) = lbound; 

    % crust
    Model.l1.bound = upper_LAY;        
    Model.l1.dens  = Model.rho_i;
    % ice
    Model.l2.bound = mid_LAY; 
    Model.l2.dens = Model.rho_c;
    % crust
    Model.l3.bound = low_LAY; 
    Model.l3.dens = Model.rho_m;
    % mantle
    Model.l4.bound = lbound; 

    % perform spherical harmonic analyses and synthesis       
    [Vlay] = model_SH_analysis(Model);

    % add to previous coefficients           
    if numl == 1
        V_Model = Vlay;
    else
        V_Model(:,3) = V_Model(:,3) + Vlay(:,3);
        V_Model(:,4) = V_Model(:,4) + Vlay(:,4);   
    end
end 


%%

% correct observations for dynamic gravity
if Model.mantle.deltarho==0
    % do nothing
else
    disp('Correct gravity data for the deep signal...')
    for ii = 1:length(VGMC(:,3))
    
        index = (V_Model(:,1)==VGMC(ii,1)&V_Model(:,2)==VGMC(ii,2));
        V_Model(index,3) = V_Model(index,3) - VGMC(ii,3);
        V_Model(index,4) = V_Model(index,4) - VGMC(ii,4);
    end
end

V_Model(3,3) = 0;
[data_M] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_Model,Model);

%% plot the final solution of the gravity field

if verbose==1   
figure
subplot(2,2,1)
imagesc(lon(1,:),lats(:,1),data.vec.R.*1e5);c=colorbar;colormap(turbo) 
hold on
% plot(coast.long,coast.lat,'k','LineWidth',1.5);
hold off
axis([lonLim(1) lonLim(2) latLim(1) latLim(2)])
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['Observation'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)

subplot(2,2,2)
imagesc(lon(1,:),lats(:,1),data_M.vec.R.*1e5);c=colorbar;colormap(turbo) 
hold on
% plot(coast.long,coast.lat,'k','LineWidth',1.5);
hold off
axis([lonLim(1) lonLim(2) latLim(1) latLim(2)])
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['Gravity Model'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)

subplot(2,2,4)
imagesc(lon(1,:),lats(:,1),(data.vec.R-data_M.vec.R).*1e5);c=colorbar;colormap(turbo) 
hold on
% plot(coast.long,coast.lat,'k','LineWidth',1.5);
hold off
axis([lonLim(1) lonLim(2) latLim(1) latLim(2)])
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['Residual'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-500 500])
set(gca,'YDir','normal','FontSize',20)

end

%% degree variance of the observation and the model

[ng,DV_g] = degreeVariance(V);
[nm,DV_m] = degreeVariance(V_Model);

if verbose==1   
subplot(2,2,3)
loglog((ng),(DV_g),'-ro','lineWidth',2)
hold on
loglog((nm),(DV_m),'-ko','lineWidth',2)
loglog((ng),abs(DV_g-DV_m(1:length(DV_g))),'--o','Color',[0.5 0.5 0.5],'lineWidth',1)
hold off
xlabel('spherical harmonic degree')
ylabel('Degree variance of the Stokes coefficients')
axis([2 116 10e-17 10e-7])
legend('Obseration',[Model.name ' model'])
end
%% polar plots

if verbose == 1
    figure
    axesm('MapProjection','ortho','origin',[90,0]);c=colorbar;colormap(turbo) 
    framem
    %load coast;
    %plotm(lat,long,'k')
    gridm
    %[gLat,gLon] = meshgrid(lat,lon);
    surfm(data.grd.lat,data.grd.lon,-moho./1e3)
    xlabel('Longitude [\circ]','FontSize',20)
    ylabel('Latitude [\circ]','FontSize',20)
    title(['South pole'],'FontSize',20)
    ylabel(c,'mGal','FontSize',20)
    %caxis([-1000 4000])
    set(gca,'YDir','normal','FontSize',20)

end


%% gravity 

if verbose==1   
figure
subplot(2,3,1)
axesm('MapProjection','ortho','origin',[90,0]);c=colorbar;colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data.vec.R).*1e5)
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['North pole observed'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)

subplot(2,3,4)
axesm('MapProjection','ortho','origin',[-90,0]);c=colorbar;colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data.vec.R).*1e5)
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['South pole observed'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)

subplot(2,3,2)
axesm('MapProjection','ortho','origin',[90,0]);c=colorbar;colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data_M.vec.R).*1e5)
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['North pole modeled'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)

subplot(2,3,5)
axesm('MapProjection','ortho','origin',[-90,0]);c=colorbar;colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data_M.vec.R).*1e5)
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['South pole modeled'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)

subplot(2,3,3)
axesm('MapProjection','ortho','origin',[90,0]);c=colorbar;colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data.vec.R-data_M.vec.R).*1e5)
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['North pole residual'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)

subplot(2,3,6)
axesm('MapProjection','ortho','origin',[-90,0]);c=colorbar;colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data.vec.R-data_M.vec.R).*1e5)
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
title(['South pole residual'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-1000 1000])
set(gca,'YDir','normal','FontSize',20)
end

%%

figure
subplot(2,2,1)
axesm('MapProjection','ortho','origin',[90,0]);%c=colorbar;
colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data.vec.R-data_M.vec.R).*1e5)
%xlabel('Longitude [\circ]','FontSize',15)
%ylabel('Latitude [\circ]','FontSize',15)
title(['North pole'],'FontSize',15)
%ylabel(c,'mGal','FontSize',15)
caxis([-250 250])
set(gca,'YDir','normal','FontSize',15)

subplot(2,2,2)
axesm('MapProjection','ortho','origin',[-90,0]);%c=colorbar;
colormap(turbo) 
framem
%load coast;
%plotm(lat,long,'k')
gridm
%[gLat,gLon] = meshgrid(lat,lon);
surfm(data.grd.lat,data.grd.lon,(data.vec.R-data_M.vec.R).*1e5)
%xlabel('Longitude [\circ]','FontSize',15)
%ylabel('Latitude [\circ]','FontSize',15)
title(['South pole'],'FontSize',15)
%ylabel(c,'mGal','FontSize',15)
caxis([-250 250])
set(gca,'YDir','normal','FontSize',15)

subplot(2,2,3:4)
imagesc(lon(1,:),lats(:,1),(data.vec.R-data_M.vec.R).*1e5);c=colorbar;colormap(turbo) 
hold on
% plot(coast.long,coast.lat,'k','LineWidth',1.5);
hold off
axis([lonLim(1) lonLim(2) latLim(1) latLim(2)])
xlabel('Longitude [\circ]','FontSize',20)
ylabel('Latitude [\circ]','FontSize',20)
%title(['Residual'],'FontSize',20)
ylabel(c,'mGal','FontSize',20)
caxis([-250 250])
set(gca,'YDir','normal','FontSize',20)
%% fitness values

Total_Fitness_std    = std(abs(DV_g(2+1:90+1)     - DV_m(2+1:90+1))./DV_g(2+1:90+1));
Dyn_Fitness_std    = std(abs(DV_g(2+1:8+1)     - DV_m(2+1:8+1))./DV_g(2+1:8+1));
Long_Fitness_std    = std(abs(DV_g(9+1:15+1)   - DV_m(9+1:15+1))./DV_g(9+1:15+1));
Middle_Fitness_std  = std(abs(DV_g(16+1:40+1)  - DV_m(16+1:40+1))./DV_g(16+1:40+1));
Short_Fitness_std   = std(abs(DV_g(41+1:116+1) - DV_m(41+1:116+1))./DV_g(41+1:116+1));

Total_Fitness_avg    = mean(abs(DV_g(2+1:90+1)  - DV_m(2+1:90+1))./DV_g(2+1:90+1));
Dyn_Fitness_avg     = mean(abs(DV_g(2+1:8+1)    - DV_m(2+1:8+1))./DV_g(2+1:8+1));
Long_Fitness_avg    = mean(abs(DV_g(9+1:15+1)   - DV_m(9+1:15+1))./DV_g(9+1:15+1));
Middle_Fitness_avg  = mean(abs(DV_g(16+1:40+1)  - DV_m(16+1:40+1))./DV_g(16+1:40+1));
Short_Fitness_avg   = mean(abs(DV_g(41+1:116+1) - DV_m(41+1:116+1))./DV_g(41+1:116+1));

Total_Fitness_max    = max(abs(DV_g(2+1:90+1)  - DV_m(2+1:90+1))./DV_g(2+1:90+1));
Dyn_Fitness_max     = max(abs(DV_g(2+1:8+1)     - DV_m(2+1:8+1))./DV_g(2+1:8+1));
Long_Fitness_max    = max(abs(DV_g(9+1:15+1)   - DV_m(9+1:15+1))./DV_g(9+1:15+1));
Middle_Fitness_max  = max(abs(DV_g(16+1:40+1)  - DV_m(16+1:40+1))./DV_g(16+1:40+1));
Short_Fitness_max   = max(abs(DV_g(41+1:116+1) - DV_m(41+1:116+1))./DV_g(41+1:116+1));

Total_Fitness_min    = min(abs(DV_g(2+1:90+1)  - DV_m(2+1:90+1))./DV_g(2+1:90+1));
Dyn_Fitness_min     = min(abs(DV_g(2+1:8+1)     - DV_m(2+1:8+1))./DV_g(2+1:8+1));
Long_Fitness_min    = min(abs(DV_g(9+1:15+1)   - DV_m(9+1:15+1))./DV_g(9+1:15+1));
Middle_Fitness_min  = min(abs(DV_g(16+1:40+1)  - DV_m(16+1:40+1))./DV_g(16+1:40+1));
Short_Fitness_min   = min(abs(DV_g(41+1:116+1) - DV_m(41+1:116+1))./DV_g(41+1:116+1));
%%
Map_avg = mean(mean((data.vec.R-data_M.vec.R).*1e5));
Map_std = std(std((data.vec.R-data_M.vec.R).*1e5));
Map_min = min(min((data.vec.R-data_M.vec.R).*1e5));
Map_max = max(max((data.vec.R-data_M.vec.R).*1e5));

FITNESS = [Total_Fitness_avg  Total_Fitness_std  Total_Fitness_min  Total_Fitness_max;...
           Dyn_Fitness_avg    Dyn_Fitness_std    Dyn_Fitness_min    Dyn_Fitness_max;...
           Long_Fitness_avg   Long_Fitness_std   Long_Fitness_min   Long_Fitness_max;...
           Middle_Fitness_avg Middle_Fitness_std Middle_Fitness_min Middle_Fitness_max;...
           Short_Fitness_avg  Short_Fitness_std  Short_Fitness_min  Short_Fitness_max;...
           Map_avg            Map_std            Map_min            Map_max]

if verbose==1   
figure
subplot(2,2,1:2)
hold on
plot((ng),abs(DV_g-DV_m(1:length(DV_g)))./DV_g,'--o','Color',[0.5 0.5 0.5],'lineWidth',1)
hold off
xlabel('spherical harmonic degree')
ylabel('Relative difference wrt. observations')
axis([0 116 0 Dyn_Fitness_max+0.5])

[sc_obs] = vecml2sc(V(:,3),V(:,4),lmax);
[sc_mod] = vecml2sc(V_Model(:,3),V_Model(:,4),lmax);

subplot(2,2,4);
imagesc(log10(abs(sc_obs(1:91,:)-sc_mod(1:91,:))));colormap(turbo);colorbar;
toc

subplot(2,2,3)
loglog((ng),(DV_g),'-ro','lineWidth',2)
hold on
loglog((nm),(DV_m),'-ko','lineWidth',2)
loglog((ng),abs(DV_g-DV_m(1:length(DV_g))),'--o','Color',[0.5 0.5 0.5],'lineWidth',1)
hold off
xlabel('spherical harmonic degree')
ylabel('Degree variance of the Stokes coefficients')
axis([2 116 10e-17 10e-7])
legend('Obseration',[Model.name ' model'])

end
toc

%% save data

filename_output = ['Results/nodynamics/output_run_ ' num2str(runid) '_' num2str(datenum(datetime)) '.mat'];
%save(filename_output,'FITNESS','Model','V_Model');
runid = runid +1;

%% ending the loop over all parameters.
end
end
end
end
end
end
end
end


