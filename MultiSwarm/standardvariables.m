function [varMin,varMax] = standardvariables(con,n,options,varMin,varMax)
%% Standardised variables to be held constant during optimisation

% variCons = {"Variables",    "Num Of",   "Conditions"    'Transformations';...
%     "Dihedral",             "~",        "~",            '~';...
%     "Chord",                n+1,        "< Previous",   '.*AftLength';...
%     "LESweep",              n,          "~"             '~';...
%     "Semispan",             n,          "Minimum 0.5",  '~';...            
%     "Section",              n+1,        "Floor",        '~';...
%     "xOffset",              "~",        "~",            '.*AftLength';...
%     "zOffset",              "~",        "~"             '.*AftHeight/2';...
%     "UpperLength",          "~",        "~",            '~';...
%     "yUpperRad",            "~",        "~",            '~';...
%     "yBotRatio",            "~",        "~",            '~';...
%     "zUpperRad",            "~",        "~",            '~';...
%     "SideLength",           "~",        "~",            '~';...
%     "zLowerRad",            "~",        "~",            '~';...
%     "AftLength",            "~",        "~",            '~';...
%     "NoseRad",              "~",        "~",            '~';...
%     "NoseLength",           "~",        "~",            '.*NoseRad';...
%     "zNoseOffset",          "~",        "~",            '.*AftHeight/2';...
%     "ForeLength",           "~",        "~",            '~'};

stanChord = repmat(0.8,1,n+1);
stanLESweep = repmat(40,1,n);
stanSemispan = [2, repmat(1,1,n-1)];

if options.Bezier
    stanSection = [1, 0.9, 0.7, 0.5, 0.3, 0,  0, 1, 0.9, 0.7, 0.5, 0.3, 0, 0,  0, 0.035, 0.04, 0.07, 0.04, 0.05, 0,  0, -0.015, -0.02, -0.05, -0.02, -0.05, 0];
else
    stanSection = 1;
end

stanSection = repmat(stanSection,1,n+1);

if options.control
    
    stanControl = [0.4,0.7,0.7];
    % Cylindrical body
    standardVar = [0, stanChord, stanLESweep, stanSemispan, stanSection, 0,0,... % Wing
        0,1,1, 1,0,1, 8, 0.25,0.25,0, 4,... % Body
        stanControl]; % Control
    
    % X-34 Body (reverse-transformed, ie. if transformations in 
    % conditioning change these also need to)
%     standardVar = [0, stanChord, stanLESweep, stanSemispan, stanSection, 0,0,... % Wing
%         0,0.88,0.1, 0.88,0.759,0.05, 11.8956, 0.155,0.7193,-0.38, 4.4238,... % Body
%         stanControl]; % Control

else
    
    % Cylindrical body
    standardVar = [0, stanChord, stanLESweep, stanSemispan, stanSection, 0,0,... % Wing
        0,1,1, 1,0,1, 8, 0.25,0.25,0, 4]; % Body
    
    % X-34 Body (reverse-transformed, ie. if transformations in 
    % conditioning change these also need to)
%     standardVar = [0, stanChord, stanLESweep, stanSemispan, stanSection, 0,0,... % Wing
%         0,0.88,0.1, 0.88,0.759,0.05, 11.8956, 0.155,0.7193,-0.38, 4.4238]; % body
end

varMin(con) = standardVar(con);
varMax(con) = standardVar(con);