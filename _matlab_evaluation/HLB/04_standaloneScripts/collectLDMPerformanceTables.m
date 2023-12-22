clear;

folderELDM = 'C:\git\KDP_Igneczi\publikációk\ELDM\tables_ELDM_new';
folderLDM = 'C:\git\KDP_Igneczi\publikációk\ELDM\tables_LDM_new';



markers = ["o", "+", "diamond", "v", "^", "*", ".", "x", "square", "<"];

files = dir(fullfile(folderELDM,'*Model*.xlsx'));
summarySCELDM = figure();
hold on;
for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
    for j=1:size(t,1)
        z(i,j) = t(j,1);
    end
end
z(z>0.7) = 1;
imagesc(z);
grid on;
yticks(1:length(drID));
yticklabels(drID);
xlabel ('Curve ID');
title('ELDM side correctness');
summarySCELDM.Position = [100 100 450 300];
xlim([0,11]);

figure(1);
plot(min(z'));
hold on;

summarySCLDM = figure();
files = dir(fullfile(folderLDM,'*Model*.xlsx'));
hold on;
for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
%     j = max(1,round(rand(1,1)*10));    
%     p(i) = plot(t(:,1), 'marker',markers(j), 'DisplayName', drID(i), 'LineStyle', 'none');
    for j=1:size(t,1)
        z(i,j) = t(j,1);
    end
end
z(z>0.7) = 1;
xlim([0,11]); zlim([0,1]);
imagesc(z);
grid on;
yticks(1:length(drID));
yticklabels(drID);
xlabel ('Curve ID');

title('LDM side correctness');
summarySCLDM.Position = [100 100 450 300];

colorbar;

figure(1);
plot(min(z'));
avgSC = figure();

for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(4:5));
%     j = max(1,round(rand(1,1)*10));    
%     p(i) = plot(t(:,1), 'marker',markers(j), 'DisplayName', drID(i), 'LineStyle', 'none');
     z(i,1:2) = [mean(t(:,1)), std(t(:,1))];
end
disp('LDM mean side correctness:');
disp(z(:,1));

e = errorbar(z(:,1),z(:,2),'LineWidth', 2, 'DisplayName', 'LDM');
hold on;
xticks(1:length(drID));
xticklabels(drID);
ylim([0,1]);
grid on;
title('Mean side correctness');

files = dir(fullfile(folderELDM,'*Model*.xlsx'));

for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(4:5));
    z(i,1:2) = [mean(t(:,1)), std(t(:,1))];
end

e = errorbar(z(:,1),z(:,2),'LineWidth', 2, 'DisplayName', 'ELDM');
legend;

disp('ELDM mean side correctness:');
disp(t(:,1)');

%%% MEDIAN between driver and model %%%

files = dir(fullfile(folderELDM,'*Model*.xlsx'));
summaryMedELDM = figure();
grid on; hold on;
for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
%     j = max(1,round(rand(1,1)*10));    
%     p(i) = plot(t(:,1), 'marker',markers(j), 'DisplayName', drID(i), 'LineStyle', 'none');
    for j=1:size(t,1)
        z(i,j) = t(j,3);
    end
end
yticks(1:length(drID));
yticklabels(drID);
surf(z);
title('ELDM median deviation');
disp('Mean median deviation of ELDM');
disp(mean(z'));

summaryMedLDM = figure();
files = dir(fullfile(folderLDM,'*Model*.xlsx'));
grid on; hold on;
for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
    for j=1:size(t,1)
        z(i,j) = t(j,3);
    end
end
yticks(1:length(drID));
yticklabels(drID);
surf(z);
title('LDM median deviation');

disp('Mean median deviation of LDM');
disp(mean(z'));

avgMed = figure();

for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
     z(i,1:2) = [mean(t(:,3)), std(t(:,3))];
end

e = errorbar(z(:,1),z(:,2),'LineWidth', 2, 'DisplayName', 'LDM');
hold on;
xticks(1:length(drID));
xticklabels(drID);
ylim([0,1]);
grid on;
title('Mean median deviation');

files = dir(fullfile(folderELDM,'*Model*.xlsx'));

for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
    z(i,1:2) = [mean(t(:,3)), std(t(:,3))];
end

e = errorbar(z(:,1),z(:,2),'LineWidth', 2, 'DisplayName', 'ELDM');
legend;

%%% MEAN between driver and model %%%

files = dir(fullfile(folderELDM,'*Model*.xlsx'));
summaryMeanELDM = figure();
grid on; hold on;
for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
    for j=1:size(t,1)
        z(i,j) = t(j,2);
    end
end
yticks(1:length(drID));
yticklabels(drID);
surf(z);
title('ELDM mean deviation');

summaryMeanLDM = figure();
files = dir(fullfile(folderLDM,'*Model*.xlsx'));
grid on; hold on;
for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
    for j=1:size(t,1)
        z(i,j) = t(j,2);
    end
end
yticks(1:length(drID));
yticklabels(drID);
surf(z);
title('LDM mean deviation');

avgMean = figure();

for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
     z(i,1:2) = [mean(t(:,2)), std(t(:,2))];
end

e = errorbar(z(:,1),z(:,2),'LineWidth', 2, 'DisplayName', 'LDM');
hold on;
xticks(1:length(drID));
xticklabels(drID);
ylim([0,1]);
grid on;
title('Mean deviation');

files = dir(fullfile(folderELDM,'*Model*.xlsx'));

for i=1:length(files)
    t = xlsread(fullfile(files(i).folder, files(i).name));
    drID(i) = convertCharsToStrings(files(i).name(1:5));
    z(i,1:2) = [mean(t(:,2)), std(t(:,2))];
end

e = errorbar(z(:,1),z(:,2),'LineWidth', 2, 'DisplayName', 'ELDM');
legend;


