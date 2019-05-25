function streamline = streamlineslow(partCell,flow)

Vinf = flow.U;
model = 'barycentrictrimodel';
[timestep,newTimestep] = deal(0.001);

for ii = numel(partCell):-1:1
    
    part = partCell{ii};
    L = part.L;
    
    close all
    plotter(part,"triangle")
    
    if isfield(part,'Triangle')
       
        triangle = part.Triangle;
    else
        triangle = part;
    end
    
    x = part.Points(:,:,1);
    y = part.Points(:,:,2);
    z = part.Points(:,:,3);
    
    triID = triangle.TriID;
    
    dim = length(triID);
    
    norm = part.unitNorm;
    
    nxTri = triangle.unitNorm(:,:,1);
    nyTri = triangle.unitNorm(:,:,2);
    nzTri = triangle.unitNorm(:,:,3);
    
    skip = ~(nxTri | nyTri | nzTri);
    
    cxTri = triangle.centre(:,:,1);
    cyTri = triangle.centre(:,:,2);
    czTri = triangle.centre(:,:,3);
    
    % Plane equation
    d = nxTri .* cxTri + nyTri .* cyTri + nzTri .* czTri;
    
    %%
    % Creating x,y point matrices from 0 to [x y]-1
    % Number of points
    [rows,cols,~] = size(norm);
    
    streamRows = rows*3;
    streamCols = cols*3;
    
    [streamIDs,time] = deal(nan(streamRows,streamCols));
    nanStreamIDs = nan(streamRows,1);
    
    nant = nan(streamRows,1);
    nanStream = nan(streamRows,3);
    streams = nan(streamRows,streamCols,3);
    
    %% Centre point velocities
    
    T = crossmat(norm, permute(Vinf,[3 1 2]));
    Vcentre = crossmat(T, norm);
    
    % Translate to corner velocity
    Vcorner = L * Vcentre(:);
    Vc = reshape(Vcorner,size(part.Points));
    
    Vx = Vc(:,:,1);
    Vy = Vc(:,:,2);
    Vz = Vc(:,:,3);

    falseMat = false(dim,1);
    
    [IDArray,todo] = deal((1:dim)');
    
    IDmat = zeros(rows,cols);
    IDmat(:) = 1:numel(nxTri);

    begin = nan(streamCols,1);
    begin(1:cols) = IDmat(1,:);   
    TE = IDmat(end,:);
    
    i = 1;
    count = 1;
    maxRow = 0;
    
    while ~isempty(todo)

        ID = begin(i);
        rowTri = triID(ID,:);
        
        p1 = squeeze([cxTri(ID); cyTri(ID); czTri(ID)]);
        Vp1 = mean([Vx(rowTri); Vy(rowTri); Vz(rowTri)],2);
        
        streamCoords(1,count,1) = p1(1);
        streamCoords(1,count,2) = p1(2);
        streamCoords(1,count,3) = p1(3);
        
        tID = 0;
        j = 1;
        stop = false;
        switched = false;
        cross = true;
        [prevCross,saveCross] = deal([]);
        streamID = nanStreamIDs;
        stream = nanStream;
        t = nant;
        noCross = 0;
        
        while ~stop
            
            xp1 = p1(1);
            yp1 = p1(2);
            zp1 = p1(3);
            
            if cross
                
                % Probably need to get rid of the trailing edge condition
                if skip(ID) || any(TE == ID)
                    
                    con = todo == ID;
                    todo(con) = [];
                    ID = ID + 1;
                    
                    if any(TE == ID) || ID > dim
                        
                        break
                    else
                        rowTri = triID(ID,:);
                        continue
                    end
                end
                
                nxi = nxTri(ID);
                nyi = nyTri(ID);
                nzi = nzTri(ID);
                
                d0 = d(ID);
                
                c = [x(rowTri); y(rowTri); z(rowTri)];
                V = [Vx(rowTri); Vy(rowTri); Vz(rowTri)];
                
                diff21 = c(:,2) - c(:,1);
                diff31 = c(:,3) - c(:,1);
                
                cDiff = (diff21 + diff31)/2;
                
%                 zmajor = (abs(Vp1(3)) >= abs(Vp1(2)) & nyi ~= 0) | nzi == 0;
                zmajor = (abs(cDiff(3)) >= abs(cDiff(2)) & nyi ~= 0) | nzi == 0;
                
                streamID(j) = ID;
                stream(j,:) = p1;
                j = j + 1;
                noCross = 0;
            end
            
            if zmajor
                
                if cross
                    
                    p_2D = [c([1 3],:); V([1 3],:)];
                end
                
                pt_2D = p1([1 3]);
                Vpt_2D = Vp1([1 3]);
            else
                
                if cross
                    
                    p_2D = [c([1 2],:); V([1 2],:)];
                end
                
                pt_2D = p1([1 2]);
                Vpt_2D = Vp1([1 2]);
            end
            
            % RK4 integration
            [k1, flag(1)] = feval(model, [pt_2D; Vpt_2D] ,timestep,p_2D);
            [k2, flag(2)] = feval(model, [pt_2D; Vpt_2D] + k1/2 ,timestep,p_2D);
            [k3, flag(3)] = feval(model, [pt_2D; Vpt_2D] + k2/2 ,timestep,p_2D);
            [k4, flag(4)] = feval(model, [pt_2D; Vpt_2D] + k3 ,timestep,p_2D);
            
            if any(flag)
                
                while any(flag) && newTimestep >= 1e-10
                    
                    newTimestep = newTimestep/10;
                    
                    [k1, flag(1)] = feval(model, [pt_2D; Vpt_2D] ,newTimestep,p_2D);
                    [k2, flag(2)] = feval(model, [pt_2D; Vpt_2D] + k1/2 ,newTimestep,p_2D);
                    [k3, flag(3)] = feval(model, [pt_2D; Vpt_2D] + k2/2 ,newTimestep,p_2D);
                    [k4, flag(4)] = feval(model, [pt_2D; Vpt_2D] + k3 ,newTimestep,p_2D);
                end
                
                newTimestep = timestep;
            end
            
            delta = (k1 + 2*k2 + 2*k3 + k4)/6;
            
            pt_2D = pt_2D + delta([1 2]);
            Vpt_2D = Vpt_2D + delta([3 4]);
            
            if zmajor
                
                xp2 = pt_2D(1);
                zp2 = pt_2D(2);
                
                yp2 = (d0 - nxi * xp2 - nzi * zp2)/nyi;
                
                Vp2 = [Vpt_2D(1) Vp1(2) Vpt_2D(2)]';                 
            else
                xp2 = pt_2D(1);
                yp2 = pt_2D(2);
                
                zp2 = (d0 - nxi * xp2 - nyi * yp2)/nzi;
                
                Vp2 = [Vpt_2D([1 2]); Vp1(3)]; 
            end
            
            %% Crossings
            [con1, con2, con3] = deal(falseMat);
            
            p2 = [xp2 yp2 zp2]';
            
            if isequal(p1,p2)
                
                break
            end
            
            % Has particle crossed between panel points 1 & 2
            [crossCoords, cross1, p2, Vp2] = ifcross(p1,p2,Vp2,c(:,1),c(:,2),V(:,1),V(:,2),prevCross);
            
            if cross1
                
                con1 = sum(ismember(triID,rowTri([1 2])),2);
                saveCross = crossCoords;
            end
            
            % Has particle crossed between panel points 2 & 3
            [crossCoords, cross2, p2, Vp2] = ifcross(p1,p2,Vp2,c(:,2),c(:,3),V(:,2),V(:,3),prevCross);
            
            if cross2
                   
                con2 = sum(ismember(triID,rowTri([2 3])),2);
                saveCross = crossCoords;
            end
            
            % Has particle crossed between panel points 3 & 1
            [crossCoords, cross3, p2, Vp2] = ifcross(p1,p2,Vp2,c(:,1),c(:,3),V(:,1),V(:,3),prevCross);
            
            if cross3
                
                con3 = sum(ismember(triID,rowTri([1 3])),2);
                saveCross = crossCoords;
            end
            
            cross = cross1 || cross2 || cross3;
            
            if cross
                
                switched = false;
                
                if cross1 && cross2
                    
                    con = [con1, con2];
                    
                elseif cross2 && cross3
                    
                    con = [con2, con3];
                    
                elseif cross1 && cross3
                    
                    con = [con1, con3];
                else
                    con = [con1, con2, con3];
                    
                    pick = IDArray(any(con > 1,2));
                    
                    ID = pick(pick ~= ID);
                    rowTri = triID(ID,:);
                    
                    xp2 = p2(1);
                    yp2 = p2(2);
                    zp2 = p2(3);
                end
                
                if isempty(rowTri)
                    
                    stop = true;
                end
                
            else
                
                within = iswithintriangle(p2,c);
                if ~within
                    
                    if ~switched
                        
                        zmajor = ~zmajor;
                        switched = true;
                        continue
                    else
                        con = todo == ID;
                        todo(con) = [];
                        break
                    end
                else
                    noCross = noCross + 1;
                    
                    if noCross >= 10
                        
                        con = todo == ID;
                        todo(con) = [];
                        break
                    end
                end
            end
            
            %% Plotting every timestep
%             figure(gcf)
%             hold on
%             plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'r')
%             plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'b*')
%             plot3(cxTri(ID),cyTri(ID),czTri(ID),'b*')
%             hold off
            
            pdiff = ((xp2-xp1).^2 + (yp2-yp1).^2 + (zp2-zp1).^2).^0.5;
            
            Vavg = (Vp1+Vp2)./2;
            Vmag = (Vavg(1).^2 + Vavg(2).^2 + Vavg(3).^2).^0.5;
            
            if isnan(Vmag)
                
                break
            end
            
            tID = tID + pdiff./Vmag;
            
            if cross
                
                t(j-1) = tID;
                tID = 0;
            end
            
            Vp1 = Vp2;
            p1 = p2;
            prevCross = saveCross;
        end
        
        %% Remove IDs from todo list   
        
        % Remove NaNs first
        delete = isnan(streamID);
        streamID(delete) = [];
        
        % Find IDs that already exist in streamline matrix
%         con = ismember(streamID,streamIDs);
        
        leng = length(streamID);
        
%         % Terminating panel always remains in todo as it hasn't expanded
%         % into another panel
%         con(leng) = false;
%         revised = false(leng,1);
%         
%         if rem(leng,2)
%             
%             sepStreams = nan(leng,(leng+1)/2);
%         else
%             sepStreams = nan(leng,leng/2);
%         end
%         
%         for I = 1:leng-1
%             
%             if con(I) && con(I+1)
%                 % If sequence is repeat of one already tracked, expansion
%                 % already exists, remove initial panel
%                 revised(I) = true;
%                 
%             elseif con(I) && ~con(I+1) && ismember(streamID(con(I)),todo)
%                 % If sequence hasn't been tracked, AND panel expanding into
%                 % not yet tracked panel has not expanded into a panel
%                 % already, keep and remove from todo
%                 true
%                 
%             else
%                 
%                 
%             end
%         end
        
        if leng == 1
            
            todo(ismember(todo,streamID)) = [];
        
        elseif leng > 1
            
            todo(ismember(todo,streamID(1:end-1))) = [];
            
            rows = 1:leng;
            
            if leng > maxRow
                    
                maxRow = leng;
            end
            
            %% Failsafes to ensure matrix consistency
            if leng > streamRows
                
                need = leng - streamRows;
                streamIDs = [streamIDs; nan(need,streamCols)];
                streams = [streams; nan(need,streamCols,3)];
                streamRows = leng;
            end
            
            if count > streamCols
                
                need = count - streamCols;
                streamIDs = [streamIDs, nan(streamRows,need)];
                streams = [streams, nan(streamRows,need,3)];
                streamCols = count;
            end
            
            time(rows,count) = t(rows);
            streamIDs(rows,count) = streamID;
            streams(rows,count,1) = stream(rows,1);
            streams(rows,count,2) = stream(rows,2);
            streams(rows,count,3) = stream(rows,3);
            count = count + 1;
        end
        
        i = i + 1;
        
        if i > cols && ~isempty(todo)
            
            begin(i) = todo(1);
        end
    end
    
    figure(gcf)
    hold on
    plot3(streams(:,:,1),streams(:,:,2),streams(:,:,3),'b')
    hold off
    
    rowArray = 1:maxRow;
    colArray = 1:count-1;
    
    streamline(ii).Time = time(rowArray,colArray);
    streamline(ii).ID = streamIDs(rowArray,colArray);
end

end