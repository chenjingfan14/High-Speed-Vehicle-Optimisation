function foilData = getaerofoilsecdata()

foilFiles = struct2cell(dir('VehicleGen\Aerofoil_Sections'));

foilFiles = foilFiles(1,:);
j=1;
for i=3:length(foilFiles)
    fid = fopen(foilFiles{i},'r');
    foilData(j) = textscan(fid, '%f%f', 'HeaderLines', 1, 'Collect', 1);
    fclose(fid);
    j=j+1;
end