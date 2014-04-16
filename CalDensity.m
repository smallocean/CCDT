%%%%%%%%%%%%%%%%% integrate density on each triangle %%%%%%%%%%%
function d_t = CalDensity(x, delaunay_tri, img, qx, val)

t = delaunay_tri.Triangulation;
width = size(img, 2);
height = size(img, 1);

%% find pixels closed to each sampling point
density = zeros(length(x),1);   %% record densities on those closed pixels
for i=1:length(x)
    p = x(i,:);
    px = round(p(1));
    if px<=0
        px = 1;
    end
    if px>width
        px = width;
    end
    py = height-round(p(2));
    if py<=0
        py = 1;
    end
    if py>height
        py = height;
    end
    v = (256-double(img(py, px)))/256;
    density(i) = v;   
end

%% find triangles that each pixel belongs to
si = pointLocation(delaunay_tri,qx); 

%% eliminate pixels outside triangulation
no_triangle = isnan(si);
in_triangle = find(no_triangle==0);
si = si(in_triangle);
val_in = val(in_triangle);
val_idx = ones(length(in_triangle),1);

%% accumulate capacity of grayscale values inside each triangle
A = accumarray(si,val_in);      % capacities
np = accumarray(si, val_idx);   % number of pixels inside each triangle
tail = zeros(size(t,1)-size(A,1),1);
A = [A; tail];
np = [np; tail];

%% calculate average pixel density on each triangle
d_t = A./np;
no_pixel = find(np==0);   %%% For those tiny triangles that do not contain
for i=1:length(no_pixel)  %%% pixels inside, just avarage the grayscale values
    tidx = no_pixel(i);   %%% on the three triangle vertices 
    d_t(tidx) = sum(density(t(tidx,:),:))/3;
end

