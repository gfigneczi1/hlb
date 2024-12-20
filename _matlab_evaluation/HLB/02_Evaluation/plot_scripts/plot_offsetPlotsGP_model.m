function f = plot_offsetPlotsGP_model(offsets, path)

    relSegment = [0, 1200];
    relSegment = [1.4185e6, 1.4205e6];

    f = figure();
    f.Position = [10 10 1000 1500];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');

    for i=1:length(offsets)
        if (~isempty(offsets{i}))
            if (mean(diff(offsets{i}.X)) > 0) % 62B
                subplot(5,2,2);
                plot(offsets{i}.X, offsets{i}.offset, offsets{i}.marker,  'LineWidth', 1.5, 'DisplayName', offsets{i}.name);
                grid on; hold on;
            else % 62A
                subplot(5,2,1);
                plot(offsets{i}.X, offsets{i}.offset, offsets{i}.marker,  'LineWidth', 1.5, 'DisplayName', offsets{i}.name);
                grid on; hold on;
            end
        end
    end

    subplot(5,2,1);
    xlabel('X-UTM(m)'); ylabel("$\delta(m)$");
    yline(0, "HandleVisibility","off", 'Alpha',0.3, 'Color','k', 'LineWidth',3);

    legend("Location", "best", "Orientation","horizontal", "NumColumns",5);

    set(gca,'FontSize', 14);
    xlim(relSegment);
    ylim([-1,1]);

    subplot(5,2,2);
    xlabel('X-UTM(m)'); ylabel("$\delta(m)$");
    yline(0, "HandleVisibility","off", 'Alpha',0.3, 'Color','k', 'LineWidth',3);

    legend("Location", "best", "Orientation","horizontal", "NumColumns",5);

    set(gca,'FontSize', 14);
    xlim(relSegment);
    ylim([-1,1]);

    inputs = ["$o_t$", "$fo_t$", "$v_x(m/s)$", "$a_x(m/s^2)$", "$\omega_z(1/s)$", "$\kappa(1/m)$", "$d\kappa(1/m^2)$"];
    
    for n=1:length(offsets)
        if (~isempty(offsets{n}))
            for i=1:7
                subplot(5,2,3+(i-1));
                input = offsets{n}.U(i,:)'*offsets{n}.GP_params.s_in(i)+offsets{n}.GP_params.c_in(i);
                plot(offsets{n}.X, input, offsets{n}.marker);
                grid on; hold on;
    
                xlabel('X-UTM(m)'); ylabel(inputs(i));        
                set(gca,'FontSize', 14);
                xlim(relSegment);
            end
        end
    end

    subplot(5,2,10);
    plot(path(:,1), path(:,2), 'LineWidth',2, 'Color','k');
    hold on;
    xlim(relSegment);
    grid on;
    xlabel("X-UTM(m)"); ylabel("Y-UTM(m)");

    set(gca,'FontSize', 14);
    title("Vehicle path");
end

