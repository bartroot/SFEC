Mars Mantle-Lithosphere interaction on Mars

This scenario lets you study the gravitational and dynamic topogrpahy of a mantle plume underneath the Tharsis Region of Mars.

To generate a mantle plume, you need to run the `create_tom_files.m` program. Three different shapes are possible to model: Sphere, Disk, Pillar. You can specify the size, density differences and depth. When using this scenario, this is done within the `Main.m` program.

**To run this scenario** 

Please open the `Main.m` program in Matlab (version older than 2022a are tested)

Lines 18 and 19 need to be adjusted, such that the directories of SFEC and GSH are correctlly linked.

SFEC_dir = '<directory where you placed the SFEC code>';       % add your directory of the github repositories 
GSH_dir = '<directory where you placed the GSH code>/Tools/'; 

