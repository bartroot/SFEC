function [] = create_tom_files(Model)
%%%%%%%%%# author : Maxime VINCENT (University of Paris) 
%%%%%%%%%# contact: maxime_vincent@outlook.fr
%%%%%%%%%# date   : July/August 2021
%%%%%%%%%# For    : TU Delft 
%%%%%%%%%# With   : Bart Root (TU Delft) --supervisor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change record: 04-05-2022 (Bart Root) made suitable for full simulation
% Mars project, connecting to Flexure code of Weilun Qin.
% 
% 16-09-2022: Bart Root added different shapes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Shape = Model.mantle.shape; 		  	       %% Sphere, Disk, Pillar, (TBD: Plume)

deltarho=Model.mantle.deltarho;           %% kg/m^3    %% set the value of the density anomaly 
layer_size=Model.mantle.layer_size;       %% km        %% set the thickness of each layer 
depth_anom=Model.mantle.depth_anom;       %% km        %% set the depth of the anomaly 
depth_core=Model.core_depth;		  
center=Model.mantle.center;   		  %% °         %% set the center of the anomaly -- index of longitude  & co-latitude
kmradius=Model.mantle.kmradius;           %% km        %% set the radius in kilometers of the anomaly 
lmax = Model.nmax;         		  %%           %% set the maximum spherical harmonics degree of the tom file 
degree_resolution=(lmax+1)/180; 	  %% °         %% set the lat lon resolution 

%% ###### CONSTANTS AND NECESSARY DEFINITIONS    

R = Model.Re./1e3;                                  %% km        %% radius of Mars
latLim = [-89.75 89.75 degree_resolution];          %% °         %% limits of the lat lon grid
lonLim = [0.25 359.75 degree_resolution];           %% °         %% limits of the lat lon grid
lon = lonLim(1):lonLim(3):lonLim(2);                %%°          %% creation of longitude vector 
lat = -(latLim(1):latLim(3):latLim(2));             %%°          %% creation of latitude vector                   %%°          %% creation of lat/lon grid
DD=fliplr(20:layer_size:depth_core/1e3);            %% km        %% creation of radial vector -> depth wrt surface  
[~,closestIndex] = min(abs(DD-depth_anom))	       ;         %% finding the right index to have a middle layer at the wished depth 
middle_layer=closestIndex ;                                      %% index corresponding to the depth where we will place the middle of the anomaly 
tot_number_layer = length(DD);                                   %% total number of layers of model 

if strcmp(Shape,'Sphere')
    number_layer_blob=2*kmradius/layer_size;                     %% number of layers the anomaly will spread on to have a height equal to twice the radius 
elseif strcmp(Shape,'Disk')
    thickness = Model.mantle.thickness;
    number_layer_blob=2*thickness/layer_size;                    %% number of layers the anomaly will spread on to have a height equal to twice the radius 
end

%% ###### INITIALIZATION OF ARRAY'S AND VECTORS 

ind=0;                                                		 %% initialization of a var for creation of anom in shape=1
circle=zeros(length(lat),length(lon));                		 %% creation of the matrix of zeros for creation of the density anomaly 

%% ########### DENSITY FIELD INITIALIZATION

fid_read = fopen([Model.sfec_dir 'outputmax.tom'],'wt');         		 %% modify here to your own path
fprintf(fid_read,'%5d\t%5d\n',[tot_number_layer lmax]);
              
%% ########## CREATION OF THE DENSITY FIELD LAYER BY LAYER   

for ll = 1:tot_number_layer  
%disp('we are here'); disp(ll)
mat=circle;                                            		 %% we define a new var with the zeros array that will receive the density anomaly 
depth = DD(ll);
Radius = R- depth;         				         %% radius of layer                    

fprintf(fid_read,'%10.3E\n',Radius);

%% ########### SHAPE CONSTRUCTION

if strcmp(Shape,'Sphere')
    % disp('--------------Sphere model for Mantle convection is constructed-------------')
    % construct the density structure of the layer for a sphere model
    for i=1:length(lon)
    for z=1:length(lat)
    
        if ll==1                        			 %% parametrisation of 1/4 arclenght radius concatenated with the same but flipped to have a half circle 
            xxx1=(0:floor(number_layer_blob/2));
            xxx1=xxx1/max(xxx1) ;
            yyy1=sqrt(1-xxx1.^2);
            xxx2=((floor(number_layer_blob/2)+1:number_layer_blob)/number_layer_blob);
            xxx2=xxx2-min(xxx2);
            xxx2=flip(xxx2/max(xxx2));
            yyy2=sqrt(1-xxx2.^2);
            yyy=([yyy2 , yyy1]);
        end
          
        if ll>middle_layer-floor((number_layer_blob)/2) -1 && ll< middle_layer+floor((number_layer_blob)/2)
            degreeradius=rad2deg(kmradius/(R-DD(ll)));
            if i==1 && z==1
                ind=ind+1;
            end
            if sqrt((lon(i)-lon(center(1))).^2 +((lat(z)-lat(center(2)))).^2)< degreeradius*yyy(ind)        % product of radius in degree at the middle layer with the arc lenght param 
                mat(z,i)=deltarho;
            end
        end
    end
    end

elseif strcmp(Shape,'Disk')
    
    % disp('--------------Disk model for Mantle convection is constructed-------------')
    % construct the density structure of the layer for a disk model
    for i=1:length(lon)
    for z=1:length(lat)
          
        if ll>middle_layer-floor((number_layer_blob)/2) -1 && ll< middle_layer+floor((number_layer_blob)/2)
            degreeradius=rad2deg(kmradius/(R-DD(ll)));
            
            if sqrt((lon(i)-lon(center(1))).^2 +((lat(z)-lat(center(2)))).^2)< degreeradius      % product of radius in degree at the middle layer with the arc length param 
                mat(z,i)=deltarho;
            end
        end
    end
    end

elseif strcmp(Shape,'Pillar')

    % disp('--------------Pillar model for Mantle convection is constructed-------------')
    % construct the density structure of the layer for a pillar model
    degreeradius=rad2deg(kmradius/(R-DD(ll)));

    for i=1:length(lon)
    for z=1:length(lat)          
        
        if ll<middle_layer
            if sqrt((lon(i)-lon(center(1))).^2 +((lat(z)-lat(center(2)))).^2)< degreeradius      % product of radius in degree at the middle layer with the arc length param 
                mat(z,i)=deltarho;
            end
        end

    end
    end

elseif strcmp(Shape,'Plume')

    % construct the density structure of the layer for a plume model
    error(['This particular shape is not yet possible: ' Shape])

else

    error(['This particular shape is not possible: ' Shape])

end

%% ########### CONVERT COEFFICIENTS IN THE 4-PI NORMALIZATION 

cs = GSHA(mat,lmax); sc = cs2sc(cs);  

%% ########### CONVERT COEFFICIENTS TO THE ZDENEK'S NORMALIZATION.

[Clm,Slm] = sc2zdenek(sc,lmax);

Clm(1) = 0;         % no mean density, SFEC cannot work with that
Slm(1) = 0;         % no mean density, SFEC cannot work with that

%% ########### WRITE COEFFICIENT IN THE OUTPUT FILE IN THE CORRECT FORMAT 

for ii = 1:length(Slm)
	fprintf(fid_read,'%23.15E %23.15E\n',[Clm(ii) Slm(ii)]);
end

clear Clm Slm

end
 
%% #########-----------END-------------##########
