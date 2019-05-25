%% Post-Processor
% Can be used to visualise configurations and their characteristics
% To do so it (re)runs prediction techniques. Also saves configurations and
% characteristics in parProperties, parameters, results

if ~exist('configInputs','var')
    
    load('OptimisationResults.mat')
end

%% IMPORTANT
% Can be run independently so long as: Initialiser has been run and there
% is a configuration to create, ie. configInputs is set of valid
% configuration inputs (can be direct or indirect)

% close all

%% Choose if only certain configurations to be output
select = [12];
% select = [];

%% Important Parameters
% Define what to plot, can be in single string containing various options,
% or string array
% Options: CN CA Cm Cl Cd Cl_Cd L D M Aerofoils
what = "Aerofoils";

% How many configurations & corresponding analysis graphs to plot at once
atOnce = 5;

%%
[nPop,~] = size(configInputs);

if exist('isDirect','var') && isDirect
    
    disp('Variables assumed to be actual configuration measurements')
    physicalPos = configInputs;
    
else
    % Impose conditions on particles
    [~,physicalPos] = contrans(configInputs,nPop,cond,varArray,options);
end

if isfield(options,'ChordDisc')
    
    lDisc = options.ChordDisc;
    uDisc = flip(lDisc);
else
    uDisc = options.UpperDisc;
    lDisc = options.LowerDisc;
end

line = ["k","k--","k-."];

switch aerofoilMethod
    
    case "BP3434"
        
        sections = BP3434(parPos,varArray,n,nPop,lDisc,uDisc);
        
    case "Bezier"
        
        % Create 2D aerofoil section Bezier curves from control points
        sectionArray = varArray == "Bezier";
        sectionPos = physicalPos(:,sectionArray);
        
        sections = Bezier3(sectionPos,n,options.BezierControlPoints,nPop,lDisc,uDisc);
        
    case "BezierTC"
        
        % Create 2D aerofoil section Bezier curves from control points
        sectionArray = varArray == "Bezier";
        sectionPos = physicalPos(:,sectionArray);
        
        sections = BezierTC(sectionPos,n,options.BezierControlPoints,nPop,lDisc,uDisc);
        
    case "Preloaded"
        
        % Assign 2D section matrices to particles. Foils variable = section indices
        sectionArray = varArray == "Section";
        sectionPos = physicalPos(:,sectionArray);
        
        zero = sectionPos == 0;
        sectionPos(zero) = 1;
        sections = options.foilData(sectionPos);
end

semispanArray = varArray == "Semispan";
semispans = physicalPos(:,semispanArray);

count = 1;
configCount = 1;
plotted = [];

if ~isempty(select)
    
    configs = select;
else
    configs = 1:nPop;
end

bool = isfield(options,'Base');

if bool
    
    configs(end + 1) = configs(end) + 1;
end

configs = fliplr(configs);

for i = configs
    
    if bool && i == configs(1)
        
        base = options.Base;
        
        baselinefun(variCons,options,true);
        
        configParameters = base.Parameters;
        configProperties = [];
        configResults = base.Results;
        configSections = base.Sections;
        configStr = 'Baseline Configurations';
        
        semispan = base.Parameters.Semispan;
    else
        configSections = sections(i,:);
        
        % Create configuration
        [~,configResults,configProperties,fricData,configInputs(i,:),configParameters] = particlecreator(configInputs(i,:),physicalPos(i,:),varArray,configSections,options);
        
        points = flowfinder(configProperties,options);
        
        formatSpec = 'Configuration %i';
        configStr = sprintf(formatSpec,i);
        
        semispan = semispans(i,:);
        
        plotter(points,"title",configStr)
        %         plotter(points)
    end
    
    Cl = configResults.Cl;
    Cd = configResults.Cd;
    Cm = configResults.Cm;
    CN = configResults.CN;
    CA = configResults.CA;
    rootMoment = configResults.RootMoment;
    L = configResults.Lift;
    D = configResults.Drag;
    
    %% Plot desired characteristics for all flight states
    
    flowCons = flow.FlightStates;
    
    legDim = flow.Dim(2);
    graphs = flow.Dim(3);
    firstDim = unique(flowCons(:,1));
    secondDim = unique(flowCons(:,2));
    thirdDim = unique(flowCons(:,3:end),'stable','rows');
    
    order = flow.PlotOrder;
    
    for ii = 1:numel(what)
        
        str = what(ii);
        
        for j = legDim:-1:1
            
            legStr(j,:) = string([char(order(2)) ' = ' num2str(secondDim(j))]);
        end
        
        for j = 1:graphs
            
            titleStr = configStr;
            
            for k = 1:size(thirdDim,2)
                
                titleStr = [titleStr ' ' char(order(k+2)) ' = ' num2str(thirdDim(j,k))];
            end
            
            if contains(str,"CN")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,CN(:,:,j))
                xlabel(order(1),'Interpreter','latex','FontSize',10)
                ylabel('C_N','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"CA")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,CA(:,:,j))
                xlabel(order(1),'Interpreter','latex','FontSize',10)
                ylabel('C_A','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"Cm")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,Cm(:,:,j))
                xlabel(order(1),'Interpreter','latex','FontSize',10)
                ylabel('C_m','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"Cl")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,Cl(:,:,j))
                xlabel(order(1))
                ylabel('C_L','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"Cd")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,Cd(:,:,j))
                xlabel(order(1),'Interpreter','latex','FontSize',10)
                ylabel('C_D','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"Cl_Cd")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(Cd(:,:,j),Cl(:,:,j))
                xlabel('C_D','Interpreter','latex','FontSize',10)
                ylabel('C_L','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"L")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,L(:,:,j))
                xlabel(order(1),'Interpreter','latex','FontSize',10)
                ylabel('Lift (N)','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"D")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,D(:,:,j))
                xlabel(order(1),'Interpreter','latex','FontSize',10)
                ylabel('Drag (N)','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"M")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                plot(firstDim,rootMoment(:,:,j))
                xlabel(order(1),'Interpreter','latex','FontSize',10)
                ylabel('Root Bending Moment (Nm)','Interpreter','latex','FontSize',10)
                title(titleStr,'Interpreter','latex','FontSize',10)
                legend(legStr)
                set(gca,'TickLabelInterpreter','latex','FontSize', 10);
                hold off
            end
            
            if contains(str,"Aerofoils")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                box on
                set(currentFig, 'Position', pos)
                
                yLim = [inf -inf];
                
                for k = 1:length(configSections)
                    
                    if k == 1 || semispan(k-1) > 0
                        
                        plot(configSections{k}(:,1),configSections{k}(:,2), line{k})
                        sectionLegStr(k) = sprintf("Section %i", k);
                        
                        yLim = [min(yLim(1),min(configSections{k}(:,2))) max(yLim(2),max(configSections{k}(:,2)))];
                    end
                end
                
                yLim(1) = yLim(1) - 0.1;
                yLim(2) = yLim(2) + 0.1;
                
                xlabel('x','Interpreter','latex','FontSize',14)
                ylabel('z','Interpreter','latex','FontSize',14)
                title([configStr ' Aerofoil Sections (Root Chord = Section 1)'],'Interpreter','latex','FontSize',10)
                legend(sectionLegStr,'Interpreter','latex','FontSize',14)
                set(gca,'TickLabelInterpreter','latex','FontSize', 14);
                axis('equal',[0 1 yLim]);
                
                hold off
            end
            
            clear sectionLegStr
            
            %             if contains(string,"body Cp")
            %                 % Plots Cp at every radial location of body
            %                 % Produces hundreds of figures so leave commented during simulations
            %                 for ii=1:dim
            %
            %                     meanRadLoc = mean(bodyRadLoc{ii},1);
            %
            %                     [chordPanels,spanPanels] = size(bodyCp{ii});
            %                     x = 0:1/(chordPanels-1):1;
            %                     for jj=1:spanPanels
            %
            %                         location = round(meanRadLoc(jj),1);
            %
            %                         figure
            %                         hold on
            %                         grid on
            %                         title(['Body chordwise pressure coefficient at ' num2str(location) '^o'] );
            %                         plot(x,bodyCp{ii}(:,jj))
            %                         xlabel('x')
            %                         ylabel('Cp')
            %                     end
            %                 end
            %             end
        end
    end
    
    
    %% Stops if max number of configs plotted
    % Asks for user input to continue > closes current figures
    
    plotted(count) = i;
    
    if bool && i == configs(1)
        
    elseif rem(count,atOnce) == 0 || i == configs(end)
        
        formatSpec = ['Plotted configurations ' repmat('%i, ',1,count)];
        
        if i ~= configs(end)
        
            str = sprintf([formatSpec '\nEnter to plot next set \n'],plotted);
            input(str)
            close all
        else
            fprintf([formatSpec '\n'],plotted)
        end
        
        plotted = [];
        count = 1;
    else
        count = count + 1;
    end
    
    AllProperties{configCount} = configProperties;
    AllParameters(configCount) = configParameters;
    AllResults(configCount) = configResults;
    configCount = configCount + 1;
end