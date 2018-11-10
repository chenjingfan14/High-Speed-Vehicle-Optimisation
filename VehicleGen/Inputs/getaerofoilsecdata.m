function [foilData,numDatFiles] = getaerofoilsecdata()

foilFiles = struct2cell(dir('VehicleGen\Aerofoil_Sections'));

foilFiles = foilFiles(1,:);
j=0;
for i=1:length(foilFiles)
    if contains(foilFiles{i},["dat","DAT","Dat"])
        j=j+1;
        fid = fopen(foilFiles{i},'r');
        foilData(j) = textscan(fid, '%f%f', 'HeaderLines', 1, 'Collect', 1);
        fclose(fid);
    end
end

numDatFiles = j;

end