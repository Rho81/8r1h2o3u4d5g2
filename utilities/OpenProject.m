% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% Rho Corporation @ 2016 - PhD Communications
% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *

% Clean Workspace From Previous Open Projects
bdclose all
clear
clc

% Use Simulink Project API to get the current project:
project = simulinkproject;
ProjectRootDir = project.RootFolder;

% Get Root Project Directory
%ProjectRootDir = fileparts(mfilename('fulpath'));

% Add Project Folders To Path
addpath(ProjectRootDir);
addpath(fullfile(ProjectRootDir,'data'),'-end');
addpath(fullfile(ProjectRootDir,'documents'),'-end');
addpath(fullfile(ProjectRootDir,'libraries'),'-end');
addpath(fullfile(ProjectRootDir,'models'),'-end');
addpath(fullfile(ProjectRootDir,'work'),'-end');

% Set the location of slprj to be the "work" folder of the current project:
myCacheFolder = fullfile(ProjectRootDir, 'work');
if ~exist(myCacheFolder, 'dir')
    mkdir(myCacheFolder)
end

% Move Matlab Generated Files to Work Directory
Simulink.fileGenControl('set', ... 
                        'CacheFolder', ...
                        fullfile(ProjectRootDir,'work'), ...
                        'CodeGenFolder', ...
                        fullfile(ProjectRootDir,'work') ...
                        );

% Set Project Configurations
% SDR_init();
% BST_Start;
% boardP = matlabshared.internal.BoardParameters('sdrz');
% setParam(boardP,'SDRZARMTargeting',true)
% Program_FPGA;
% Set_BST_Test_Parameters

% Change working folder to the "work" folder:
cd(myCacheFolder);