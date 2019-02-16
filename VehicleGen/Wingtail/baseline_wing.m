function wingDefs = baseline_wing()
%% Define baseline wing

wingDefs = {...
    "Dihedral",     4;
    "Chord",       [13.85, 4.387, 1.5241];
    "TESweep",     [-5, -5];
    "Semispan",    [0.744+0.88, 2.62686];           
    "Section",     ["FSPLStrake.dat", "FSPLND.dat", "FSPLND.dat"];
    "xOffset",     -3.5;
    "zOffset",     -0.74};

[row,~] = size(wingDefs);

% Find which set of variables correspond to 2D aerofoil sections
for i = row:-1:1
    
    con(i) = wingDefs{i,1} == "Section";
end

if any(con)
    
    files = wingDefs{con,2};
    
    [~,foilName,nFoils] = getaerofoilsecdata();
    numArray = (1:nFoils)';
    numMat = repmat(numArray,1,numel(files));
    
    equals = files == foilName;
    
    wingDefs{con,2} = numMat(equals)';
end
