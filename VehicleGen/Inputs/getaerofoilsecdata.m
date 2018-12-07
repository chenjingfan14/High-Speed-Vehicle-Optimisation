function [foilData,numDatFiles] = getaerofoilsecdata()

foilFiles = struct2cell(dir([pwd '\VehicleGen\Aerofoil_Sections']));
foilFiles = foilFiles(1,:);

j=0;
for i=1:length(foilFiles)
    if contains(foilFiles{i},["dat","DAT","Dat"])
        j=j+1;
        fid = fopen(foilFiles{i},'r');
        
        fileName = convertCharsToStrings(foilFiles{i});
        data = textscan(fid, '%f%f', 'HeaderLines', 1, 'Collect', 1);
        
        foilData(j,:) = {fileName, cell2mat(data)};
        fclose(fid);
    end
end

if j ~= 0
    numDatFiles = j;
    save('2DAerofoilSections.mat','foilData')
else
    % Try loading files
    load('2DAerofoilSections.mat')
    foilData = AerofoilSections;
    numDatFiles = numel(foilData);
end

end