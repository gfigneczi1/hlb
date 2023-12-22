clear; 

folder = 'C:\git\kdp_hlb_evalframework\_temp\plots';

files = dir(fullfile(folder,'*.mat'));

for i=1:length(files)
    d = load(fullfile(files(i).folder, files(i).name));
    if (isfield(d, 'segment'))
        d = d.segment;
    end
end
