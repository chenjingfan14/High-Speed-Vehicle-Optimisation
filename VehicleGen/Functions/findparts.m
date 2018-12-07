function [wingDim,aftbodyDim,forebodyDim,noseDim,controlDim] = findparts(partArrays)

% Initialise existence of parts as false
[dim,~] = size(partArrays);

[wingCon,aftbodyCon,forebodyCon,noseCon,controlCon] = deal(false(dim,1));

% Loop through configuration and set parts to true at that point in the
% configuration cell when part name = "Wing" or "Nose" etc
for i = 1:dim
    wingCon(i) = strcmp(partArrays{i,1},"Wing") | strcmp(partArrays{i,1},"Tail");
    aftbodyCon(i) = strcmp(partArrays{i,1},"Aftbody");
    forebodyCon(i) = strcmp(partArrays{i,1},"Forebody");
    noseCon(i) = strcmp(partArrays{i,1},"Nose");
    controlCon(i) = strcmp(partArrays{i,1},"Control");
end

% Turn 1:dim logicals into indexes containing only numbers where said part
% type resides
array = 1:dim;
wingDim = array(wingCon);
aftbodyDim = array(aftbodyCon);
forebodyDim = array(forebodyCon);
noseDim = array(noseCon);
controlDim = array(controlCon);