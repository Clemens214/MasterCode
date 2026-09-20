function [] = Plot(name, angles, voltages, Data , choice, options)
arguments
    name
    angles
    voltages
    Data
    choice.Title = ''
    % Type of plot
    options.Spectrum = false
    options.Value = false
    options.Size = false
    options.Color = false
    options.Angles = false
    options.integrate = false
    % Dimension of plot
    options.twoD = false
    options.threeD = false
    % Data to be plotted
    options.Transmission = false
    options.Torque = false
    options.Angular = false
    options.Helicity = false
end
    Palette = colororder();
    if options.Transmission == true
        Title = 'Transmission';
    elseif options.Torque == true
        Title = 'Torquance';
    elseif options.Angular == true
        Title = 'Angular';
    elseif options.Helicity == true
        Title = 'Helicity';
    else
        Title = choice.Title;
    end
    if options.twoD == true && options.threeD == false
    % plot the Energy/voltage dependence
        if options.Spectrum == true
            fig = plotSpectrum2D(name, angles, voltages, Data);
            setLabels(fig, Title, angles, options)
    % plot the Angle dependence
        elseif options.Value == true || options.Angles == true
            fig = plotValue2D(name, angles, voltages, Data);
            setLabels(fig, Title, voltages, options)
    % plot the Size dependence
        elseif options.Size == true
            fig = plotSpectrum2D(name, angles, voltages, Data);
            setLabels(fig, Title, angles, options)
        end
        setTicks (fig, angles, voltages, Data, options)
        resizeFig(fig)
    % plot the Data in 3D
    elseif options.twoD == false && options.threeD == true
        plot3D (name, Title, angles, voltages, Data)
    end
    % plot in Color
    if options.Color == true
        plotColor(name, Title, angles, voltages, Data)
    end
    % plot the Angles
    if options.Angles == true
        if options.twoD == true && options.threeD == false
            plotValue2D (name, Title, angles, Vals)
        elseif options.threeD == true
            plotAngles3D (name, Title, angles, Vals)
        end
    end
    % change plot colors
    %colororder("gem12")
    colororder({'b', 'r', 'g', 'm', 'c', 'y', 'k'});
    % change font size
    fontsize(12,"points")
    fontname("Helvetica")
    % export the plot
    cleanfigure;
    filename = strcat(Title, '.tex');
    matlab2tikz(filename)
end

%% Helping functions
function [] = setLabels (figure, Title, values, options)
arguments
    figure
    Title
    values
    options
end
    % set the title of the plot
    if false
        title(Title);
    end
    % set the x-label of the plot
    options.integrate = false;
    if options.Spectrum == true && options.integrate == true
        xlabel('V [t]');
    elseif options.Spectrum == true && options.integrate == false || options.Size == true
        xlabel('\omega [t]');
    elseif options.Value == true
        xlabel('\Delta\theta');
    end
    % set the y-label of the plot
    if options.Transmission == true
        label = 'T(\omega)';
        unit = '';
    elseif options.Torque == true
        label = '\tau(\omega)';
        unit = '';
    elseif options.Angular == true
        label = 'L_z(\omega)';
        unit = '[1/t]';
    elseif options.Helicity == true
        label = 'h(\omega)';
        unit = '[1/t]';
    end
    ylabel( strcat(label, ' ', unit) );
    % set the legend of the plot
    if options.Value == true && options.integrate == true
        labels = strcat('V=',cellstr(num2str(values.')));
    elseif options.Value == true && options.integrate == false
        labels = strcat('\omega=',cellstr(num2str(values.')));
    elseif options.Spectrum == true
        labels = strcat('\Delta\theta=',cellstr(num2str(values.')));
        labels = cellfun(@(x) [x,'\pi'], labels, 'uniform',false);
    elseif options.Size == true
        labels = strcat('N=',cellstr(num2str(values.')));
    end
    legend(labels, 'Location','northoutside', 'NumColumns', 2);%, 'Interpreter','latex');
end

function [] = setTicks (figure, angles, voltages, Data, options)
arguments
    figure 
    angles
    voltages
    Data
    options
end
    % set the limits of the x-axis
    xLimits = xlim;
    if options.Spectrum == true || options.Size == true
        xMin = min(voltages);
        xMax = max(voltages);
    elseif options.Value == true || options.Angles == true
        xMin = min(angles);
        xMax = max(angles);
    end
    xlim([xMin, xMax])
    % set the limits of the y-axis
    yLimits = ylim;
    yMaxima = zeros(1, length(Data));
    yMinima = zeros(1, length(Data));
    for i = 1:length(Data)
        yMaxima(i) = max(Data{i});
        yMinima(i) = min(Data{i});
    end
    yMax = 1.1*max(yMaxima);
    if options.Transmission == true
        yMin = min(yMinima);
    else
        yMin = 1.1*min(yMinima);
    end
    ylim([yMin, yMax])
    % set the x-axis ticks
    if options.Value == true
        TickLabels = cellfun(@num2str , xticklabels, 'uniform',false);
        TickLabels = cellfun(@(x) [x,'\pi'], TickLabels, 'uniform',false);
        xticklabels(TickLabels)
    end
end

function [] = resizeFig (figure)
arguments
    figure 
end
    % get the conversion factor
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
    PosNew = [PosOld(1)/Convert, PosOld(2)/Convert, WidthNew, HeightNew];
    PosNew = [PosOld(1)/Convert, PosOld(2)/Convert, WidthNew, WidthNew];
    set(figure, 'Units','centimeters', 'Position', PosNew)
end

%% plotting functions: Energies (+angles)
function [fig] = plotSpectrum2D (name, angles, voltages, Data)
    TransPlot = cell(1, length(angles));
    for i = 1:length(angles)
        TransPlot{i} = zeros(1, length(voltages));
        for j = 1:length(voltages)
            TransPlot{i}(j) = Data{i}(j);
        end
    end
    % plot the data
    fig = figure(Name=name);
    hold on
    for i = 1:length(angles)
        plot(voltages, TransPlot{i}, linewidth=1);
    end
    hold off
    grid off
end

%% plotting functions: Angle (+voltages)
function [fig] = plotValue2D (name, angles, voltages, Data)
    TransPlot = cell(1, length(voltages));
    for i = 1:length(voltages)
        TransPlot{i} = zeros(1, length(angles));
        for j = 1:length(angles)
            TransPlot{i}(j) = Data{j}(i);
        end
    end
    % plot the data
    fig = figure(Name=name);
    hold on
    for i = 1:length(voltages)
        plot(angles, TransPlot{i}, linewidth=1);
    end
    hold off
    grid off
end

%% plotting functions: 3D
function [] = plot3D (name, Title, angles, voltages, Data)
    TransPlot = zeros(length(voltages), length(angles));
    for i = 1:length(Data)
        TransPlot(:, i) = Data{i}.';
    end
    figure(Name=name);
    surf(angles, voltages, TransPlot)
    xlabel('Angle (°)');
    ylabel('Voltage (a.u.)'); 
    zlabel([Title, ' (a.u.)']);
    title(Title);
end

function [] = plotColor (name, Title, angles, voltages, Data)
    DataPlot = zeros(length(voltages), length(angles));
    for i = 1:length(Data)
        DataPlot(:, i) = Data{i}.';
    end
    figure(Name=name);
    surf(angles, voltages, DataPlot,'EdgeColor', 'None', 'facecolor', 'interp');
    view(2);
    colorbar;
    xlabel('Angle (°)');
    ylabel('Voltage (a.u.)'); 
    zlabel([Title, ' (a.u.)']);
    title(Title);
end

%% plotting functions: Angles
function [] = plotAngles3D (name, Title, angles, Vals)
    figure(Name=name);
    surf(angles, angles, Vals)
    title(Title);
end