function [] = construct_SFEC_inp(Model)
% 
% function to build sfec.inp file for the SFEC models

fid_read = fopen([Model.sfec_dir 'sfec.inp'],'wt');         %% modify here to your own path

fprintf(fid_read,'# Spectral intervals and radial grid\n');
fprintf(fid_read,'Min_degree              %-5d\t no      ! Min. harmonic degree for the response\n', Model.sfec.min_degree);
fprintf(fid_read,'Max_degree              %-5d\t no      ! Max. harmonic degree for the response\n',Model.sfec.max_degree);
fprintf(fid_read,'N_node_layer            %-5d\t no      ! Number of FE-points within each tomographic layer\n',2);
fprintf(fid_read,'Theta_points            %-5d\t no      ! Number of Gauss-Legendre points\n',360);

fprintf(fid_read,'\n%s\n','# Tomographic model');
fprintf(fid_read,'Tomography_file         %-5d\t outputmax.tom  ! File containing density derived from tomography\n',999);

fprintf(fid_read,'\n%s\n','# Self gravitation');
fprintf(fid_read,'Self_gravitation        %-5d\t no      ! Selfgravitation (1=yes, 0=no)\n',1);

fprintf(fid_read,'\n%s\n','# Viscosity structure');
[b,e] = getExp(Model.sfec.UMV);
fprintf(fid_read,'Visc_UM                 %1dd%-2d\t no      ! Upper mantle viscosity\n',b,e);
[b,e] = getExp(Model.sfec.LMV);
fprintf(fid_read,'Visc_LM                 %1dd%-2d\t no      ! Lower mantle viscosity\n',b,e);
[b,e] = getExp(Model.sfec.LIV);
fprintf(fid_read,'Visc_Lith               %1dd%-2d\t no      ! Lithosphere viscosity\n',b,e);

fprintf(fid_read,'\n%s\n','# CG solver');
fprintf(fid_read,'CG_tol_test             %-5d\t no      ! Tolerance test for CG solver\n',1);
fprintf(fid_read,'CG_tol                  %1dd%-2d\t no      ! Tolerance of CG solver\n',1,-8);
fprintf(fid_read,'CG_maxit                %-5d\t no      ! Maximum number of CG iterations\n',100);

fclose(fid_read);

end

%% get function that decouples the exponential number from the value.
function [b,e] = getExp(val)

e = floor(log10(val));
b = val / 10^e;

end
