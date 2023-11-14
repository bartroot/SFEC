function [] = construct_geometry_mod(Model)
% 
% function to build sfec.inp file for the SFEC models

fid_geo = fopen([Model.sfec_dir geometry.f90'],'wt');         %% modify here to your own path

fprintf(fid_geo,'module geometry\n');
fprintf(fid_geo,'! -----------------------------------------------\n');
fprintf(fid_geo,'! constants used for defining the planet geometry\n');
fprintf(fid_geo,'! -----------------------------------------------\n\n');
% G value and finding 

[b,e] = getExp(Model.G);
fprintf(fid_geo,'real(kind=8), parameter :: gi=%.3fd%02d\n', b,e);
[b,e] = getExp(Model.GM./Model.Re.^2);
fprintf(fid_geo,'real(kind=8), parameter :: gref=%.2fd%02d\n', b,e);
fprintf(fid_geo,'real(kind=8), parameter :: gcmb=3.21d0\n');
[b,e] = getExp(Model.rho_m);   
fprintf(fid_geo,'real(kind=8), parameter :: derhot=%.1fd%02d\n', b,e);
fprintf(fid_geo,'real(kind=8), parameter :: derhoc=-3.5d3\n');
fprintf(fid_geo,'real(kind=8), parameter :: rho0=3.7d3\n');
   
[b,e] = getExp(Model.Re);
fprintf(fid_geo,'real(kind=8), parameter :: erad=%.3fd%02d\n', b,e);
fprintf(fid_geo,'real(kind=8), parameter :: radc=1.835d6\n'); % mean value from insight?
% viscosity interface
fprintf(fid_geo,'real(kind=8), parameter :: rm=2.897d6\n');
% Elastic lithosphere
[b,e] = getExp(Model.Re-Model.Te);
fprintf(fid_geo,'real(kind=8), parameter :: rlith=%.3fd%02d\n', b,e);   

fprintf(fid_geo,'\nend module geometry\n');

fclose(fid_geo);

end 

%% get function that decouples the exponential number from the value.
function [b,e] = getExp(val)

e = floor(log10(val));
b = val / 10^e;

end
