function [partArrays,sectionArray] = partIndexing(cond,varArray)

wing = false;
aft = false;
fore = false;
nose = false;
control = false;

[dim,~] = size(cond);
uniqueDim = size(unique(varArray));
counters = zeros(uniqueDim);
partArrays = cell(dim,3);
j = 0;

for i = 1:dim
    
    name = cond{i,1};
    array = cond{i,2};
    
    up = name == varArray;
    counters(up) = counters(up) + 1;
    
    % Use for > 1 parts (ie. wing & tail)
    % partNum = counters(up);
    
    switch name
        case {"Dihedral","Chord","LESweep","Semispan","Bezier","Section","xOffset","zOffset"}
            
            if ~wing
                j = j + 1;
                wing = j;
                partArrays{wing,1} = "Wing";
            end
            
            n = size(array);
            partArrays{wing,2} = [partArrays{wing,2} repmat(name,n)];
            partArrays{wing,3} = [partArrays{wing,3} array];
            
            if name == "Section" || name == "Bezier"
                sectionArray = array;
            end
            
        case {"NoseRad","NoseLength","zNoseOffset"}
            
            if ~nose
                j = j + 1;
                nose = j;
                partArrays{nose,1} = "Nose";
            end
            
            n = size(array);
            partArrays{nose,2} = [partArrays{nose,2} repmat(name,n)];
            partArrays{nose,3} = [partArrays{nose,3} array];
            
        case "ForeLength"
           
            if ~fore
                j = j + 1;
                fore = j;
                partArrays{fore,1} = "Forebody";
            end
            
            n = size(array);
            partArrays{fore,2} = [partArrays{fore,2} repmat(name,n)];
            partArrays{fore,3} = [partArrays{fore,3} array];
            
        case {"ControlChord","ControlSpan"}  
            
            if ~control
                j = j + 1;
                control = j;
                partArrays{control,1} = "Control";
            end
            
            n = size(array);
            partArrays{control,2} = [partArrays{control,2} repmat(name,n)];
            partArrays{control,3} = [partArrays{control,3} array];
            
        otherwise
            
            if ~aft
                j = j + 1;
                aft = j;
                partArrays{aft,1} = "Aftbody";
            end
            
            n = size(array);
            partArrays{aft,2} = [partArrays{aft,2} repmat(name,n)];
            partArrays{aft,3} = [partArrays{aft,3} array];
            
    end
    
end

partArrays = partArrays(1:j,:);