function variCons = standardvariables(variCons,options)
%% Standardised variables to be held constant during optimisation

n = options.WingPartitions;

% X-34 Body (reverse-transformed, ie. if transformations in
% conditioning change these also need to)
Definition = {...
    "Variables",        "Standard Value"};

wingDefs = {...
    "Dihedral",         5;
    "Chord",            0.8;
    "LESweep",          40;
    "Semispan",         2;
    "SectionDefinition",[];
    "xOffset",          0;
    "zOffset",          0};

aftDefs = {...
    "UpperLength",      0;
    "yUpperRad",        0.88;
    "yBotRatio",        0.1;
    "zUpperRad",        0.88;
    "SideLength",       0.759;
    "zLowerRad",        0.05;
    "AftLength",        11.8956};

% Cylindrical Aftbody
% standardAft = {...
%     "UpperLength",      0;
%     "yUpperRad",        1;
%     "yBotRatio",        1;
%     "zUpperRad",        1;
%     "SideLength",       0;
%     "zLowerRad",        1;
%     "AftLength",        8};

foreDefs = {...
    "ForeLength",       4.4238};

noseDefs = {...
    "NoseRad",          0.155;
    "NoseLength",       0.7193;
    "zNoseOffset",     -0.6};

controlDefs = {...
    "ControlSpan",     [0.4,0.7];
    "ControlChord",     0.7};

Bezier = options.Bezier;
wing = options.Wing;
aft = options.Aft;
fore = options.Fore;
nose = options.Nose;
control = options.Control;
baseline = options.Baseline;

[rows,~] = size(variCons);

if wing
    
    if Bezier
        str = "Bezier";
        vals = [1, 0.9, 0.7, 0.5, 0.3, 0,  0, 1, 0.9, 0.7, 0.5, 0.3, 0, 0,  0, 0.035, 0.04, 0.07, 0.04, 0.05, 0,  0, -0.015, -0.02, -0.05, -0.02, -0.05, 0];
    else
        str = "Section";
        vals = 1;
    end
    
    % Find which set of variables correspond to 2D aerofoil sections
    for i=size(wingDefs,1):-1:1
        
        secDef(i) = wingDefs{i,1} == "SectionDefinition";
    end
    
    % And Insert definition
    wingDefs(secDef,:) = {str,vals};
    
    Definition = [Definition; wingDefs];
end

if aft
    Definition = [Definition; aftDefs];
end

if fore
    Definition = [Definition; foreDefs];
end

if nose
    Definition = [Definition; noseDefs];
end

if control
    Definition = [Definition; controlDefs];
end

if baseline
    
    baseConfig = options.base;
    baseDefinition = baseConfig.Definition;
    baseVariables = baseConfig.Variables;
    baseVarArray = baseConfig.VarArray;
    basePartitions = baseConfig.nPartitions;
    baseSections = basePartitions + 1;
    baseBezier = baseConfig.Bezier;
    baseIsDirect = baseConfig.Direct;
    
    partitions = (1:n)';
    sections = (1:n+1)';
    
    [chordPrint,sweepPrint,spanPrint,sectionPrint] = deal(false);
    
    name = strings(rows,1);
    
    for i = 2:rows
        
        name(i) = variCons{i,1};
        
        min = variCons{i,2};
        max = variCons{i,3};
        
        minCon = isnan(min);
        maxCon = isnan(max);
        
        dim = 1:length(minCon);
        
        minArray = dim(minCon);
        maxArray = dim(maxCon);
        
        if any(minCon | maxCon)
            
            if name(i) ~= name(i-1)
                ID = 1;
            end
            
            stand = Definition{i,2};
            base = baseVariables(baseVarArray == name(i));
            
            for j = 1:length(min)
                
                % If yes, sections wish to be set as stardard. Must check that
                % optimisation variables and base variables are the same for the wing
                % (in terms of nPartitions) otherwise will not work
                if name(i) == "Chord"
                    
                    if sections(ID) <= baseSections
                        
                        min(j) = base(ID);
                        max(j) = base(ID);
                        
                    else
                        if ~chordPrint
                            
                            fprintf('Base chord used outboard until section %i \n', baseSections)
                            chordPrint = true;
                        end
                        
                        if length(stand) == 1
                            
                            min(minArray) = stand;
                            max(maxArray) = stand;
                        else
                            
                            min(minArray) = stand(minArray);
                            max(maxArray) = stand(maxArray);
                        end
                        
                        continue
                    end
                    
                    %% REPEAT ABOVE
                    
                elseif name(i) == "LESweep"
                    
                    if partitions(ID) <= basePartitions
                        
                        min(i) = base(ID);
                        max(i) = base(ID);
                        
                    else
                        if ~sweepPrint
                            
                            fprintf('Base sweep used outboard until partition %i \n', basePartitions)
                            sweepPrint = true;
                        end
                        
                        if length(stand) == 1
                            
                            min(minArray) = stand;
                            max(maxArray) = stand;
                        else
                            
                            min(minArray) = stand(minArray);
                            max(maxArray) = stand(maxArray);
                        end
                        
                        continue
                    end
                    
                elseif name(i) == "Semispan"
                    
                    if partitions(ID) <= basePartitions
                        
                        min(i) = base(ID);
                        max(i) = base(ID);
                    else
                        
                        if ~spanPrint
                            
                            fprintf('Base semispan used outboard until partition %i \n', basePartitions)
                            spanPrint = true;
                        end
                        
                        if length(stand) == 1
                            
                            min(minArray) = stand;
                            max(maxArray) = stand;
                        else
                            
                            min(minArray) = stand(minArray);
                            max(maxArray) = stand(maxArray);
                        end
                        
                        continue
                    end
                    
                    % Same as before although must also be defined same as optimisation
                    % variables in terms of sections or Bezier curves
                elseif any(name(i) == ["Section","Bezier"],2)
                    
                    ID = i - sectionArray == 0;
                    
                    if Bezier == baseBezier
                        
                        if sections(ID) <= baseSections
                            
                            baseSec = baseVariables(any(baseVarArray == ["Section","Bezier"],2));
                            baseSec = reshape(baseSec,n,[]);
                            
                            varMin = [varMin, baseSec(:,ID)];
                            varMax = [varMax, baseSec(:,ID)];
                        else
                            
                            if ~sectionPrint
                                
                                fprintf('Base sections used outboard until section %i \n', baseSections)
                                sectionPrint = true;
                            end
                            
                            varMin = [varMin, stand];
                            varMax = [varMax, stand];
                            
                        end
                        
                    else
                        
                        if ~sectionPrint
                            
                            fprintf('Standard sections cannot be set by base configuration due to \ndiffering definitions. \nUsing back-up standard variables instead \n')
                            sectionPrint = true;
                        end
                        
                        varMin = [varMin, stand];
                        varMax = [varMax, stand];
                        
                    end
                    
                else
                    
                    min(j) = base(ID);
                    max(j) = base(ID);

                end
                
                minArray(1) = [];
                maxArray(1) = [];
            
                ID = ID + 1;
                
            end
            
            variCons{i,2} = min;
            variCons{i,3} = max;

        end
    end
    
else
    
    % If no baseline specific, replace all desired variables with standard
    % ones specified above
    for i = 2:rows
        
        min = variCons{i,2};
        max = variCons{i,3};
        
        stand = Definition{i,2};
        
        minCon = isnan(min);
        maxCon = isnan(max);
        
        if length(stand) == 1
            min(minCon) = stand;
            max(maxCon) = stand;
        else
            min(minCon) = stand(minCon);
            max(maxCon) = stand(maxCon);
        end
        
        variCons{i,2} = min;
        variCons{i,3} = max;
    end
    
end

end