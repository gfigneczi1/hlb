close all;
clear; 

% path of optimized parameters
pathData = "C:\git\KDP_Igneczi\publikációk\LaneWandering\data\mpcOptimization\snippets_45_fminsearch_positive\successfulOptimization";
matFiles = dir(fullfile(pathData, '*mat'));

drID = 0;
parameters = [];
for i=1:length(matFiles)
    if (drID~= str2num(matFiles(i).name(23:25)))
        if (~isempty(parameters))
            driverData{drID} = parameters;
        end
        drID = str2num(matFiles(i).name(23:25));
        parameters = [];
    end    
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    parameters(end+1,1:4) = data.y(1:4)'; 
end


f = figure('Position', [100, 100, 650, 725]);
noOfDataPoints = 0;

for i=1:length(driverData)
   if (~isempty(driverData{i}))
       color{i} = rand(1,3);
   end
end

for i=1:length(driverData)
    if (~isempty(driverData{i}))
        x = linspace(noOfDataPoints+1, noOfDataPoints+size(driverData{i},1),size(driverData{i},1) );
        noOfDataPoints = x(end);
        
        for j=1:4
            subplot(4,1,j);
            if (j==2 || j==4)
                validPoints = driverData{i}(:,j) <=500;
            else 
                validPoints = driverData{i}(:,j) <=0.04;
            end
            errorbar(i,mean(driverData{i}(validPoints,j)), std(driverData{i}(validPoints,j)), 'LineWidth', 2, 'color', color{i}, 'DisplayName', strcat('Dr',num2str(i)));
            %plot(x, driverData{i}(:,j), 'Marker', 'x', 'DisplayName', num2str(i), 'LineStyle','none',  'color', color{i});
            hold on;
            plot(i, mean(driverData{i}(validPoints,j)), 'marker', 'x', 'MarkerSize', 10, 'color', color{i}, 'HandleVisibility','off')
            %plot([x(1), x(end)], [median(driverData{i}(validPoints,j)), median(driverData{i}(validPoints,j))], 'HandleVisibility','off', 'color', color{i}, 'LineWidth', 2);

%             if (j==2 || j==4)
%                 str = sprintf('%.2f', median(driverData{i}(validPoints,j)));
%                 text(x(1), median(driverData{i}(validPoints,j))+15, str);
%             else
%                 str = sprintf('%.4f', median(driverData{i}(validPoints,j)));
%                 text(x(1), median(driverData{i}(validPoints,j))+0.002, str);
%             end


            grid on;
            xticklabels({});
            xticks([]);
            
            if (j==2 ||j==4)
                ylim([0, 300]);
            else
                ylim([0,0.02]);
            end
            set(gca,'FontSize', 12);

            if (j==4)
                %xlabel('DrID');
                positionBefore = get(gca,'Position');
                legend('Orientation','horizontal', 'Location','southoutside', 'FontSize', 10);
                set(gca,'Position', positionBefore);
            end
            

            switch j
                case 1
                    title('w_{outputError}^{Drift}');
                case 2
                    title('w_{input}^{Drift}');
                case 3
                    title('w_{outputError}^{Compensation}');
                case 4
                    title('w_{input}^{Compensation}');
            end
        end
    end
end

shg;