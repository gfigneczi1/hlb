function summaryPlotsLDM()
folderLDM = 'C:\git\KDP_Igneczi\publikációk\AutomotiveInnovations_2022\plots_\plots';

files = dir(fullfile(folderLDM,'Dr*Model*.mat'));

for i=1:length(files)
    load(fullfile(files(i).folder, files(i).name));
    drID = str2num(files(i).name(3:5));
    curveID = files(i).name(end-5:end-4);
    if (contains(curveID,'_'))
        curveID = curveID(2);
    end
    curveID = str2num(curveID);
    d{drID, curveID} = KPI.curves.data;
end

noDrivers = size(d,1);
noCurves = size(d,2);
for i=1:noCurves
    f = figure();
    hold on;
    for j=1:noDrivers
        [N,EDGES] = histcounts(d{j,i}.traj(:,2)-d{j,i}.ref(:,2));
        e = EDGES(1:end-1) + EDGES(2:end);
        [t] = spline(d{j,i}.traj(:,1), d{j,i}.traj(:,2), d{j,i}.ref(:,1));
        distances = t-d{j,i}.ref(:,2);
        errorbar(j, mean(distances), std(distances));
        text(j,mean(distances),num2str(mean(distances), '%1.3f'));
        ticks{j} = strcat('Dr', num2str(j));
    end
    grid on;
    ylim([-0.5,0.5]);
    xticklabels(ticks);
    xticks([1:1:15]);
    close(f);
end

