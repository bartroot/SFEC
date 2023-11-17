Mars Mantle-Lithosphere interaction on Mars

This scenario lets you study the gravitational and dynamic topogrpahy of a mantle plume underneath the Tharsis Region of Mars.

To generate a mantle plume, you need to run the `create_tom_files.m` program. Three different shapes are possible to model: Sphere, Disk, Pillar. You can specify the size, density differences and depth. When using this scenario, this is done within the `Main.m` program.

**To run this scenario** 

Please open the `Main.m` program in Matlab (version older than 2022a are tested)

Lines 18 and 19 need to be adjusted, such that the directories of SFEC and GSH are correctlly linked.

SFEC_dir = 'directory where you placed the SFEC code';       % add your directory of the github repositories 
GSH_dir = 'directory where you placed the GSH code/Tools/'; 

at line 35 you can specify what flexural theory lithosphere model you would like to use.
at line 96 you can specify what shape for the mantle mass anomaly you would like to use.

Then lines 73-80 will allow you to specify the model parameters, or even loop through a range of values for these parameters. 

Do not forget to turn on saving the results in line 643 by uncommenting that line.

Then you are ready to go! Happy modelling!


