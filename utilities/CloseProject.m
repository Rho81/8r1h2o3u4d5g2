% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *
% Rho Corporation @ 2016 - PhD Communications
% -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- * -- *

% Clean up local customizations of the environment:

% Remove project folders from path
rmpath(ProjectRootDir);
rmpath(fullfile(ProjectRootDir,'data'),'-end');
rmpath(fullfile(ProjectRootDir,'documents'),'-end');
rmpath(fullfile(ProjectRootDir,'libraries'),'-end');
rmpath(fullfile(ProjectRootDir,'models'),'-end');
rmpath(fullfile(ProjectRootDir,'work'),'-end');

% Reset the location where generated code and other temporary files are
% created (slprj) to the default:
Simulink.fileGenControl('reset');