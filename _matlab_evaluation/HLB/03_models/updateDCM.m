function updateDCM(DcmPath,ParamName, ParamValue)
%UPDATEDCM Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(DcmPath);
lines =[];
foundParam = false;
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, ParamName)
        foundParam = true;
    end
    if foundParam && contains(tline, '   WERT')
        tline = ['   WERT ' num2str(ParamValue)];
        foundParam = false;
    end
    lines = [lines tline newline];
end
fclose(fid);
fid = fopen('03_models\kuki.dcm','wt');
fprintf(fid, '%s\n', lines);
fclose(fid);
end

