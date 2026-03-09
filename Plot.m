function [] = Plot(value, angles, voltages, Data , choice, options)
arguments
    value
    angles
    voltages
    Data
    choice.Title = ''
    % Type of plot
    options.Spectrum = false
    options.Value = false
    options.Color = false
    options.Angles = false
    % Dimension of plot
    options.twoD = false
    options.threeD = false
    % Data to be plotted
    options.Transmission = false
    options.Torque = false
    options.Both = false
end
    Title = choice.Title;
    if strcmp(Title, 'Transmission')
        label = 'T(E)';
    elseif strcmp(Title, 'Current')
        label = 'I(E) [\frac{2e}{h}]';
    elseif strcmp(Title, 'Torquance')
        label = '\tau(E)';
    elseif strcmp(Title, 'Torque')
        label = '\tau(E)';
    elseif strcmp(Title, 'Angular Momentum')
        label = 'l_{z}(E)';
    elseif strcmp(Title, 'Helicality')
        label = 'h(E)';
    else
        label = 'Result (E)';
    end
    if options.twoD == true || (options.twoD == false && options.threeD == false)
    % plot the Energy/voltage dependence
        if options.Spectrum == true
            if options.Both ==true || (options.Transmission == true && options.Torque == true)
                Transmission = Data{1};
                Torque = Data{2};
                plotSpectrumBoth(value, Title, angles, voltages, Transmission, Torque)
            else
                plotSpectrum2D(value, Title, angles, voltages, Data)
            end
    % plot the Angle dependence
        elseif options.Value == true
            if options.Both ==true || (options.Transmission == true && options.Torque == true)
                Transmission = Data{1};
                Torque = Data{2};
                plotValueBoth(value, Title, angles, voltages, Transmission, Torque)
            else
                plotValue2D(value, Title, angles, voltages, Data)
            end
        end
    % plot the Data in 3D
    elseif options.threeD == true
        plot3D (value, Title, angles, voltages, Data)
    end
    % plot in Color
    if options.Color == true
        plotColor(value, Title, angles, voltages, Data)
    end
    % plot the Angles
    if options.Angles == true
        if options.twoD == true
            plotValue2D (value, Title, angles, Vals)
        elseif options.threeD == true
            plotAngles3D (value, Title, angles, Vals)
        end
    end
    %% change plot colors
    Palette = colororder();
    %colororder("gem12")
    colororder({'b', 'r', 'g', 'm', 'c', 'y', 'k'});
    %% change font size
    fontsize(12,"points")
    fontname("Helvetica")
end

function [] = resizeFig (figure)
    set(figure, 'PaperPositionMode','auto')
    Width = 14.53; %cm
    Convert = get(0, 'ScreenPixelsPerInch') / 2.54;
    % get the figure size
    PosOld = get(figure, 'Position');
    WidthOld = PosOld(3) / Convert;
    HeightOld = PosOld(4) / Convert;
    % set the figure size
    Factor = Width / WidthOld;
    WidthNew = Factor * WidthOld;
    HeightNew = Factor * HeightOld;
    PosNew = [PosOld(1), PosOld(2), WidthNew, HeightNew];
    set(figure, 'Units','centimeters', 'Position', PosNew)
end

function [] = setLabels (figure, Title, label)
    %title(Title);
    %labels = strcat('Voltage = ',cellstr(num2str(voltages.')));
    ylabel('label');
    elseif strcmp(Title, 'Torque')
        ylabel('Torque');
    end
    labels = strcat('E = ',cellstr(num2str(voltages.')));
    if false
        legend(labels, 'Location','northoutside');
    else
        legend(labels);
    end
    if true
        TickLabels = cellfun(@num2str , xticklabels, 'uniform',false);
        TickLabels = cellfun(@(x) [x,'\pi'], TickLabels, 'uniform',false);
        xticklabels(TickLabels)
    end
end

    xlabel('E');
    if strcmp(Title, 'Transmission')
        ylabel('Transmission');
        yLimits = ylim;
        ylim([0, yLimits(2)])
    elseif strcmp(Title, 'Torque')
        ylabel('Torque');
    end
    %title(Title);
    labels = strcat('N = ',cellstr(num2str(angles.')));
    if false
        labels = cellfun(@(x) [x,'\pi'], labels, 'uniform',false);
    end
    if true
        legend(labels, 'Location','northoutside');
    else
        legend(labels);
    end
    resizeFig(fig)
end

%% plotting functions: Energies (+angles)
function [fig] = plotSpectrum2D (value, Title, angles, voltages, Data)
    TransPlot = cell(1, length(angles));
    for i = 1:length(angles)
        TransPlot{i} = zeros(1, length(voltages));
        for j = 1:length(voltages)
            TransPlot{i}(j) = Data{i}(j);
        end
    end
    % plot the data
    fig = figure(value);
    hold on
    for i = 1:length(angles)
        plot(voltages, TransPlot{i}, linewidth=1);
    end
    hold off
end

%% plotting functions: Angle (+voltages)
function [fig] = plotValue2D (value, Title, angles, voltages, Data)
    TransPlot = cell(1, length(voltages));
    for i = 1:length(voltages)
        TransPlot{i} = zeros(1, length(angles));
        for j = 1:length(angles)
            TransPlot{i}(j) = Data{j}(i);
        end
    end
    % plot the data
    fig = figure(value);
    hold on
    for i = 1:length(voltages)
        plot(angles, TransPlot{i}, linewidth=1);
    end
    hold off
end

%% plot both
function [] = plotSpectrumBoth (value, ~, angles, voltages, Transmission, Torque)
    TransPlot = cell(1, length(angles));
    TorquePlot = cell(1, length(angles));
    for i = 1:length(angles)
        TransPlot{i} = zeros(1, length(voltages));
        TorquePlot{i} = zeros(1, length(voltages));
        for j = 1:length(voltages)
            TransPlot{i}(j) = Transmission{i}(j);
            TorquePlot{i}(j) = Torque{i}(j);
        end
    end
    % plot the data
    figure(value)
    % plot the Transmission
    subplot(2,1, 1);
    xlabel('Energy (units of t)');
    ylabel('Transmission (2e/h)');
    hold on
    for i = 1:length(angles)
        plot(voltages, TransPlot{i}, linewidth=1);
    end
    hold off
    title('Transmission');
    labels = strcat('Angle = ',cellstr(num2str(angles.')));
    if true
        legend(labels, 'Location','eastoutside');
    else
        legend(labels);
    end
    grid on;

    subplot(2,1, 2);
    xlabel('Energy (units of t)');
    ylabel('Torque (-1/2\pi)');
    hold on
    for i = 1:length(angles)
        plot(voltages, TorquePlot{i}, linewidth=1);
    end
    hold off
    title('Torque');
    labels = strcat('Angle = ',cellstr(num2str(angles.')));
    if true
        legend(labels, 'Location','eastoutside');
    else
        legend(labels);
    end
    grid on;
end

function [] = plotValueBoth (value, ~, angles, voltages, Transmission, Torque)
    TransPlot = cell(1, length(voltages));
    TorquePlot = cell(1, length(voltages));
    for i = 1:length(voltages)
        TransPlot{i} = zeros(1, length(angles));
        TorquePlot{i} = zeros(1, length(angles));
        for j = 1:length(angles)
            TransPlot{i}(j) = Transmission{j}(i);
            TorquePlot{i}(j) = Torque{j}(i);
        end
    end
    % plot the data
    figure(value)
    % plot the Transmission
    subplot(2,1, 1);
    xlabel('Angle (°)'); 
    ylabel('Current (a.u.)');
    hold on
    for i = 1:length(voltages)
        plot(angles, TransPlot{i}, linewidth=1);
    end
    title('Transmission');
    hold off
    labels = strcat('Voltage = ',cellstr(num2str(voltages.')));
    if true
        legend(labels, 'Location','eastoutside');
    else
        legend(labels);
    end
    grid on;

    subplot(2,1, 2);
    xlabel('Angle (°)'); 
    ylabel('Torque (a.u.)');
    hold on
    for i = 1:length(voltages)
        plot(angles, TorquePlot{i}, linewidth=1);
    end
    hold off
    title('Torque');
    labels = strcat('Voltage = ',cellstr(num2str(voltages.')));
    if true
        legend(labels, 'Location','eastoutside');
    else
        legend(labels);
    end
    grid on;
end

%% plotting functions: 3D
function [] = plot3D (value, Title, angles, voltages, Data)
    TransPlot = zeros(length(voltages), length(angles));
    for i = 1:length(Data)
        TransPlot(:, i) = Data{i}.';
    end
    figure(value)
    surf(angles, voltages, TransPlot)
    xlabel('Angle (°)');
    ylabel('Voltage (a.u.)'); 
    zlabel([Title, ' (a.u.)']);
    title(Title);
end

function [] = plotColor (value, Title, angles, voltages, Data)
    DataPlot = zeros(length(voltages), length(angles));
    for i = 1:length(Data)
        DataPlot(:, i) = Data{i}.';
    end
    figure(value)
    surf(angles, voltages, DataPlot,'EdgeColor', 'None', 'facecolor', 'interp');
    view(2);
    colorbar;
    xlabel('Angle (°)');
    ylabel('Voltage (a.u.)'); 
    zlabel([Title, ' (a.u.)']);
    title(Title);
end

%% plotting functions: Angles (Both and Difference)
function [varargout] = plotAngles2D (value, Title, angles, Vals)
    TransPlot = zeros(1, length(angles)*length(angles));
    angleDiff = zeros(1, length(angles)*length(angles));
    indices = zeros(1, length(angles)*length(angles));
    for i = 1:length(angles)
        for j = 1:length(angles)
            idx = (i-1)*length(angles) + j;
            TransPlot(idx) = Vals(i, j);
            angleDiff(idx) = angles(i) - angles(j);
            indices(idx) = idx;
        end
    end
    varargout{1} = indices;
    [angleSort, indices] = sort(angleDiff);
    TransSort = TransPlot(indices);
    % plot the data
    figure(value);
    plot(angleSort, TransSort, linewidth=1)
    title(Title);
end

function [] = plotAngles3D (value, Title, angles, Vals)
    figure(value)
    surf(angles, angles, Vals)
    title(Title);
end