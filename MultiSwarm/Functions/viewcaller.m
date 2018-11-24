function viewcaller(parPos,cond,varArray,foilData,n,~,options)

prompt = 'Start plotting ideal configurations? (Y/N) \n';
str = input(prompt,'s');

if ~strcmp(str,'Y')
    return
end

[nPop,~] = size(parPos);

% Impose conditions on particles
[partArrays,sectionArray] = partIndexing(cond,varArray);
[~,physicalPos] = conditioning(parPos,nPop,cond,varArray,options);

% Assign 2D section matrices to particles. Foils variable = section indices
if options.Bezier
    sections = Bezier3(physicalPos(:,sectionArray),n,foilData,nPop);
else
    sections = foilData(physicalPos(:,sectionArray));
end

atOnce = 5; 
figNum = zeros(atOnce,1);
count = 1;

for i=1:nPop
    
    % Create configuration
    [parProperties,~,~] = particlecreator(parPos(i,:),physicalPos(i,:),partArrays,sections(i,:));
    
    figNum(count) = particleviewer(parProperties,i);
    
    if rem(i,atOnce) == 0 || i == nPop
        if i ~= nPop
            input(['\nPlotted optimal configurations ' num2str(i-atOnce+1) ' to '...
                num2str(i) ', \n'...
                'Enter to plot next set'])
            
             close(figNum(1:atOnce))
            
        elseif nPop == 1
            
            fprintf('\nPlotted optimal configuration')
            
        else
            
            fprintf(['\nPlotted final optimal configurations ' num2str(i-count+1) ' to '...
                num2str(i) '\n'])
        end
        
        count = 1;
        
    else
        count = count + 1;
        
    end
    
end