%% Post-Processor
% Can be used to visualise configurations and their characteristics
% To do so it (re)runs prediction techniques. Also saves configurations and
% characteristics in parProperties, parameters, results

%% IMPORTANT
% Can be run independently so long as: Initialiser has been run and there
% is a configuration to create, ie. configInputs is set of valid
% configuration inputs (can be direct or indirect)

close all

Bezier = options.Bezier;

%% Important Parameters
% Define what to plot, can be in single string containing various options,
% or string array
what = "Aerofoils";

% How many configurations & corresponding analysis graphs to plot at once
atOnce = 5;

%%
[nPop,~] = size(configInputs);

% Impose conditions on particles
[partArrays,sectionArray] = partIndexing(cond,varArray);

[~,configVar] = size(configInputs);
[~,nVar] = size(varMin);

%% FIX: Some of the configuration corrector stuff is taking values out varMin/varMax. Thus postprocess assumes values are actual dimensions whereas some need to be transformed as they are ratios

if configVar == nVar
%     outwith = configInputs < varMin | configInputs > varMax;
    outwith = false;
else
    outwith = true;
end

%%
if (exist('isDirect','var') && isDirect) || any(outwith(:))
    disp('Variables assumed to be actual configuration measurements')
    physicalPos = configInputs;
else
    [~,physicalPos] = contrans(configInputs,nPop,cond,varArray,options);
end

sectionPos = physicalPos(:,sectionArray);

% Assign 2D section matrices to particles
if Bezier
    sections = Bezier3(sectionPos,n,foilData,nPop);
else
    sections = foilData(sectionPos);
end

count = 1;

for i = nPop:-1:1
    
    configSections = sections(i,:);
    
    % Create configuration
    [configProperties,configInputs(i,:),configParameters] = particlecreator(configInputs(i,:),physicalPos(i,:),partArrays,configSections);
    % Analyse configuration and plot
    [~,configResults] = aeroprediction(configProperties,flow,configParameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
    
    points = flowfinder(configProperties);
    
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
    
    formatSpec = 'Configuration %i';
    configStr = sprintf(formatSpec,i);
    
    plotter(points,"title",configStr)
    
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
                set(currentFig, 'Position', pos)
                plot(firstDim,CN(:,:,j))
                xlabel(order(1))
                ylabel('C_N')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"CA")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(firstDim,CA(:,:,j))
                xlabel(order(1))
                ylabel('C_A')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"Cm")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(firstDim,Cm(:,:,j))
                xlabel(order(1))
                ylabel('C_m')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"Cl")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(firstDim,Cl(:,:,j))
                xlabel(order(1))
                ylabel('C_L')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"Cd")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(firstDim,Cd(:,:,j))
                xlabel(order(1))
                ylabel('C_D')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"Cl_Cd")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(Cd(:,:,j),Cl(:,:,j))
                xlabel('C_D')
                ylabel('C_L')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"L")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(firstDim,L(:,:,j))
                xlabel(order(1))
                ylabel('Lift (N)')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"D")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(firstDim,D(:,:,j))
                xlabel(order(1))
                ylabel('Drag (N)')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"M")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                plot(firstDim,rootMoment(:,:,j))
                xlabel(order(1))
                ylabel('Root Bending Moment (Nm)')
                title(titleStr)
                legend(legStr)
                hold off
            end
            
            if contains(str,"Aerofoils")
                pos = figureposition();
                figure
                currentFig = gcf;
                hold on
                grid on
                set(currentFig, 'Position', pos)
                for k = n:-1:1
                    plot(configSections{k}(:,1),configSections{k}(:,2))
                    sectionLegStr(k,:) = ['Section ' num2str(k)];
                end
                xlabel('x')
                ylabel('z')
                title([configStr ' Aerofoil Sections (Root Chord = Section 1)'])
                legend(sectionLegStr)
                hold off
            end
            
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
    if rem(i-1,atOnce) == 0 || i == 1
        
        first = i;
        last = i + count - 1;
        
        if i ~= 1
            
            formatSpec = 'Plotted configurations %i to %i, \nEnter to plot next set \n';
            
            str = sprintf(formatSpec,first,last);
            
            input(str)
            
            close all
            
        elseif nPop == 1
            
            fprintf('Plotted configuration\n')
            
        else
            
            formatSpec = 'Plotted configurations %i to %i \n';
            fprintf(formatSpec,first,last);
            
        end
        
        count = 1;
        
    else
        count = count + 1;
        
    end
    
    AllProperties{i} = configProperties;
    AllParameters(i) = configParameters;
    AllResults(i) = configResults;
    
end