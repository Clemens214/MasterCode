function [] = Plot(angleDiff, voltages, Data, options)
    % plots the calculated data
    arguments
        angleDiff
        voltages
        Data
        options.Title = 'Title'
        options.Names = {}
        options.Num = 1
        options.Traces = false
        options.angle = false
        options.voltage = false
        options.twoD = false
        options.threeD = false
        options.conductance = false
    end
    Num = options.Num;

    if options.angle == true
        figure(Num)
        anglePlot(angleDiff, Data)
        Num = Num+1;
    end

    if options.voltage == true
        figure(Num)
        voltagePlot(voltages, Data)
        Num = Num+1;
    end

    if options.twoD == true
        figure(Num)
        twoDPlot(angleDiff, voltages, Data)
        Num = Num+1;
    end

    if options.threeD == true
        figure(Num)
        title(options.Title)
        threeDPlot(angleDiff, voltages, Data)
        Num = Num+1;
    end

    if options.Traces == true
        figure(Num)
        totalTracePlot(angleDiff, voltages, Traces)
        Num = Num+1;
    end

    if options.conductance == true
        figure(Num)
        conductancePlot(angleDiff, voltages, Data, Title=options.Title, Names=options.Names)
    end
end

%% plotting functions
function [] = conductancePlot(angleDiff, Energies, Data, options)
    arguments
        angleDiff, 
        Energies, 
        Data
        options.Title = 'Title'
        options.Names = {}
    end
    sizeMax = ceil(sqrt(length(angleDiff)));
    sizeMin = ceil(length(angleDiff)/sizeMax);
    tcl = tiledlayout(sizeMax,sizeMin);
    colors = {'black', 'red', 'blue'};
    sgtitle(options.Title);
    for i = 1:length(angleDiff)
        nexttile(tcl)
        hold on
        Title = append('Angle: ',num2str(angleDiff(i)));
        for j = 1:length(Data)
            % define the data used in the plot
            PlotData = Data{j}{i};
            % plot the data
            if length(options.Names) == length(Data)
                plot(Energies, PlotData, colors{j}, DisplayName=options.Names{j});
            else
                plot(Energies, PlotData, colors{j});
            end
        end
        xlabel('Energy');
        ylabel('Torquance');
        title(Title) 
        hold off
    end
    Lgnd = legend('show');
    Lgnd.Layout.Tile = 'east';
end

function [] = anglePlot(angleDiff, Data)
    TransPlot = zeros(1, length(angleDiff));
    for i = 1:length(Data)
        TransPlot(i) = Data{i}(end);
    end
    plot(angleDiff, TransPlot)
end

function [] = voltagePlot(voltages, Data)
    TransPlot = Data{end};
    plot(voltages, TransPlot)
end

function [] = twoDPlot(angleDiff, voltages, Data)
    TransPlot = cell(1, length(voltages));
    disp(length(voltages))
    disp(length(angleDiff))
    for i = 1:length(voltages)
        TransPlot{i} = zeros(1, length(angleDiff));
        for j = 1:length(angleDiff)
            TransPlot{i}(j) = Data{j}(i);
        end
    end
    hold on
    for i = 1:length(voltages)
        plot(angleDiff, TransPlot{i});
    end
    hold off
end

function [] = threeDPlot(angleDiff, voltages, Data)
    TransPlot = zeros(length(voltages), length(angleDiff));
    for i = 1:length(Data)
        TransPlot(:, i) = Data{i}.';
    end
    surf(angleDiff, voltages, TransPlot)
end

function [] = totalTracePlot(angleDiff, voltages, Traces)
    woZero = voltages;
    woZero(woZero==0) = []; % removes the zero entry from woZero
    
    plotFactor = 1;
    sizeMin = floor(sqrt(length(woZero)));
    sizeMax = ceil(sqrt(length(woZero)));
    indexPlot = 0;
    for i = 1:length(voltages)
        if voltages(i) ~= 0
            indexPlot = indexPlot + 1;
            subplot(sizeMin, sizeMax, indexPlot);
            % define the data used in the plot
            TransPlot = zeros(length(Traces{i}{1}), length(angleDiff));
            for j = 1:length(angleDiff)
                TransPlot(:, j) = Traces{j}{i}.';
            end
            index = 0:(length(Traces{1}{1}) -1);
            % plot the data
            endIndex = floor(length(index)*plotFactor);
            surf(angleDiff, index(1:endIndex), TransPlot(1:endIndex,:));
        end
    end
end

function [] = tracePlot(angleDiff, voltages, Traces)
    woZero = voltages;
    woZero(woZero==0) = []; % removes the zero entry from woZero
    
    plotFactor = 1;
    sizeMin = floor(sqrt(length(woZero)));
    sizeMax = ceil(sqrt(length(woZero)));
    indexPlot = 0;
    for i = 1:length(voltages)
        if voltages(i) ~= 0
            indexPlot = indexPlot + 1;
            subplot(sizeMin, sizeMax, indexPlot);
            % define the data used in the plot
            TransPlot = zeros(length(Traces{i}{1}), length(angleDiff));
            for j = 1:length(angleDiff)
                TransPlot(:, j) = Traces{j}{i}.';
            end
            index = 0:(length(Traces{1}{1}) -1);
            % plot the data
            endIndex = floor(length(index)*plotFactor);
            surf(angleDiff, index(1:endIndex), TransPlot(1:endIndex,:));
        end
    end
end