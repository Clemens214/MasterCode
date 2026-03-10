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
    options.Size = false
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
        label = 'I(E) [2e/h]';
    elseif strcmp(Title, 'Torquance')
        label = '\tau (E)';
    elseif strcmp(Title, 'Torque')
        label = '\tau (E) [-1/2\pi]';
    elseif strcmp(Title, 'Angular Momentum')
        label = 'L_{z}(E)';
    elseif strcmp(Title, 'Helicity') || strcmp(Title, 'Helicality')
        label = 'h(E)';
    else
        label = 'Result (E)';
    end
    if options.twoD == true && options.threeD == false
        if min(voltages) >= 0
            options.integrate = true;
        else
            options.integrate = false;
        end
    % plot the Energy/voltage dependence
        if options.Spectrum == true
            fig = plotSpectrum2D(value, angles, voltages, Data);
            setLabels(fig, Title, label, angles, options)
    % plot the Angle dependence
        elseif options.Value == true || options.Angles == true
            fig = plotValue2D(value, angles, voltages, Data);
            setLabels(fig, Title, label, voltages, options)
        end
        %% customize plot
        resizeFig(fig)
        %% change plot colors
        Palette = colororder();
        %colororder("gem12")
        colororder({'b', 'r', 'g', 'm', 'c', 'y', 'k'});
        %% change font size
        fontsize(12,"points")
        fontname("Helvetica")
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
        if options.twoD == true && options.threeD == false
            plotValue2D (value, Title, angles, Vals)
        elseif options.threeD == true
            plotAngles3D (value, Title, angles, Vals)
        end
    end
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

function [] = setLabels (figure, Title, label, values, options)
arguments
    figure
    Title
    label
    values
    options
    %options.Spectrum = false
    %options.Value = false
    %options.Size = false
    %options.integrate = false
end
    % set the title of the plot
    if false
        title(Title);
    end
    % set the x-label of the plot
    if options.Spectrum == true && options.integrate == true
        xlabel('V');
    elseif options.Spectrum == true && options.integrate == false
        xlabel('E');
    elseif options.Value == true
        xlabel('\theta');
    end
    % set the x-axis ticks
    if options.Value == true
        TickLabels = cellfun(@num2str , xticklabels, 'uniform',false);
        TickLabels = cellfun(@(x) [x,'\pi'], TickLabels, 'uniform',false);
        xticklabels(TickLabels)
    end
    % set the y-label of the plot
    ylabel(label);
    % set the limits of the y-axis
    if strcmp(Title, 'Transmission')
        yLimits = ylim;
        ylim([0, yLimits(2)])
    end
    % set the legend of the plot
    if options.Value == true && options.integrate == true
        labels = strcat('V=',cellstr(num2str(values.')));
    elseif options.Value == true && options.integrate == false
        labels = strcat('E=',cellstr(num2str(values.')));
    elseif options.Spectrum == true
        labels = strcat('\theta=',cellstr(num2str(values.')));
        labels = cellfun(@(x) [x,'\pi'], labels, 'uniform',false);
    elseif options.Size == true
        labels = strcat('N=',cellstr(num2str(angles.')));
    end
    legend(labels, 'Location','northoutside');%, 'Interpreter','latex');
end

%% plotting functions: Energies (+angles)
function [fig] = plotSpectrum2D (value, angles, voltages, Data)
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
    grid on
end

%% plotting functions: Angle (+voltages)
function [fig] = plotValue2D (value, angles, voltages, Data)
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
    grid on
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

function [] = plotSpectrum2DBoth (value, Title, angles, voltages, Data)
    TransPlot = cell(1, length(Data));
    for i = 1:length(Data)
        TransPlot{i} = zeros(1, length(voltages));
        for j = 1:length(voltages)
            TransPlot{i}(j) = Data{i}(j);
        end
    end
    % plot the data
    figure(value)
    hold on
    %yyaxis left
    for i = 1:length(Data)%-1
        plot(voltages, TransPlot{i}, linewidth=0.5);
    end
    ylabel('Angular Momentum')
    if false
        yyaxis right
        plot(voltages, TransPlot{end}, linewidth=1);
        ylabel('Torque')
    end
    hold off
    xlabel('Energy (units of t)');
    title(Title);
    labels = strcat('Angle = ',cellstr(num2str(angles.')));
    %labels = {'Extended Molecule'; 'Semi-infinite leads'; 'Difference'};
    if true
        labels = cellfun(@(x) [x,'\pi'], labels, 'uniform',false);
    end
    legend(labels)
    fontsize(16,"points")
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