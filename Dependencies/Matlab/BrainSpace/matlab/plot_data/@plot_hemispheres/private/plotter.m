function plotter(obj)
% Workhorse of plot hemispheres. Relies exclusively on properties of the
% plot_hemispheres class and should never be called directly. 
%
% For more information, please consult our <a
%   href="https://brainspace.readthedocs.io/en/latest/pages/matlab_doc/visualization/plot_hemispheres.html">ReadTheDocs</a>.

% Get data.
data = obj.plotted_data;

% Determine color limits
for ii = 1:size(data,2)
    clim(ii,:) = [min(data(~isinf(data(:,ii)),ii)); ...
        max(data(~isinf(data(:,ii)),ii))];
end

% Account for data consisting of only one value
clim = clim + (clim(:,1) == clim(:,1)) .* [-eps,eps];

% Grab views.
[views,surf_id] = obj.process_views;

% Plot the surface(s).
obj.handles.figure = figure('Color','white','Units','normalized','Position',[0.1 0.1 .9 .38]);

% modified by Zeng
row=0.92/size(obj.data,2);
col=0.92/numel(surf_id);

for jj = 1:size(obj.data,2)
    for ii = 1:numel(surf_id)
        % Get the surface and associated data. 
        if numel(obj.surface) == 1
            D = obj.plotted_data(:,jj);
            surf = obj.surface{1};
        else
            if surf_id{ii} == 1
                D = obj.plotted_data(1:end/2,jj);
                surf = obj.surface{surf_id{ii}};
            elseif surf_id{ii} == 2
                D = obj.plotted_data(end/2+1:end,jj);
                surf = obj.surface{surf_id{ii}};
            elseif surf_id{ii} == 3
                D = obj.plotted_data(:,jj);
                surf = combine_surfaces(obj.surface{1}, obj.surface{2}, 'surfstat');
            end
        end       
        
        position =[0.01+col*(ii-1), 0.93-row*jj col row]; %modified by Zeng
%         position = [-.07+ii*.11 .99-.18*jj .17 .17]; %modified by Zeng

        [obj.handles.axes(jj,ii), obj.handles.trisurf(jj,ii), obj.handles.camlight(jj,ii)] = ...
            obj.make_surface_plot(D, surf, position, views{ii}, clim(jj,:));
    end
end

% Add colorbars.
for ii = 1:size(data,2)
    obj.handles.colorbar(ii) = colorbar(obj.handles.axes(ii,end));
    obj.handles.colorbar(ii).Position = [obj.handles.axes(ii,end).Position(1:2) + [1.05*col 0.25*row] ...
        .05*col .5*row]; %modified by Zeng
    obj.handles.colorbar(ii).Ticks = clim(ii,:);
    obj.handles.colorbar(ii).FontName = 'DroidSans';
    obj.handles.colorbar(ii).FontSize = 14;
end

% Adding labels before drawing makes the process look a bit weird. Better
% to do this first.
drawnow

% Add labels
if ~isempty(obj.labeltext)
    obj.labels();
end
end