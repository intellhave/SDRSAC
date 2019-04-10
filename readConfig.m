% Configurations for the algorithsm
% Input: 
%       configName:  Name of the config, corresponding to the desired
%       dataset
function config = readConfig(configName)

    if (nargin==0); configName = 'synthetic'; end
    
    if (strcmp(configName, 'synthetic'))
        config = readConfig_synthetic(); 
        
     % To be added
%     elseif (strcmp(configName, 'redwood'))
%         config = readConfigRedWood();    
    end
    
    
end