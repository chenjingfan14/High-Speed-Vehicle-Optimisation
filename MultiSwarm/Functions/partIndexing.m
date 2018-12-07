function [partArrays] = partIndexing(cond,varArray)

wing = false;
aft = false;
fore = false;
nose = false;
control = false;

[dim,~] = size(cond);
partArrays = cell(dim,3);
j = 0;

for i = 1:dim
    
    name = cond{i,1};
    
    switch name
        case {"Dihedral","Chord","LESweep","Semispan","Bezier","Section","xOffset","zOffset"}
            
            if ~wing
                j = j + 1;
                wing = j;
                partArrays{wing,1} = "Wing";
            end
            
        case {"NoseRad","NoseLength","zNoseOffset"}
            
            if ~nose
                j = j + 1;
                nose = j;
                partArrays{nose,1} = "Nose";
            end
            
        case "ForeLength"
           
            if ~fore
                j = j + 1;
                fore = j;
                partArrays{fore,1} = "Forebody";
            end
            
        case {"ControlChord","ControlSpan"}  
            
            if ~control
                j = j + 1;
                control = j;
                partArrays{control,1} = "Control";
            end
            
        otherwise
            
            if ~aft
                j = j + 1;
                aft = j;
                partArrays{aft,1} = "Aftbody";
            end            
    end
    
end

partArrays = partArrays(1:j,:);