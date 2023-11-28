function [a,cb]=BladeSurfStatViewData(data,surf,title,background);

%BoSurfStatViewData is a simple viewer for surface data.
% 
% Usage: [a,cb]=BoSurfStatViewData(data, surf [,title [,background]]);
% 
% data        = 1 x v vector of data, v=#vertices
% surf.coord  = 3 x v matrix of coordinates.
% surf.tri    = 3 x t matrix of triangle indices, 1-based, t=#triangles.
% title       = any string, data name by default.
% background  = background colour, any matlab ColorSpec, such as 
%   'white' (default), 'black'=='k', 'r'==[1 0 0], [1 0.4 0.6] (pink) etc.
%   Letter and line colours are inverted if background is dark (mean<0.5).
%
% a  = vector of handles to the axes, left to right, top to bottom. 
% cb = handle to the colorbar.

if nargin<3 
    title=inputname(1);
end
if nargin<4
    background='white';
end

v=length(data);
vl=1:(v/2);
vr=vl+v/2;
t=size(surf.tri,1);
tl=1:(t/2);
tr=tl+t/2;
clim=[min(data),max(data)];
if clim(1)==clim(2)
    clim=clim(1)+[-1 0];
end

clf;
colormap(spectral(256));

h=0.9;
w=0.9;

a=axes('position',[0.05 0.05 w h]);
trisurf(surf.tri,surf.coord(1,:),surf.coord(2,:),surf.coord(3,:),...
    double(data),'EdgeColor','none','FaceColor','interp');
view(0,90); 
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting phong; material dull; shading interp;

set(a,'CLim',clim);
set(a,'Tag',['SurfStatView ' num2str(i)]);


cb=colorbar(a,'location','eastoutside'); % ,'location','eastoutside'
%set(cb,'Position',[0.9 0.1 0.02 0.1]);
set(cb,'XAxisLocation','right');
h=get(cb,'Title');
set(h,'String',title);

whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');

dcm_obj=datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@SurfStatDataCursor,'DisplayStyle','window');

return
end
