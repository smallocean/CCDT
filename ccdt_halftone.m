function [embed, newt]=ccdt_halftone(nb, imgname) 

%% extracting image information
img = imread([imgname,'.jpg']);
img = rgb2gray(img);
close all;
imshow(img);
width = size(img, 2);
height = size(img, 1);

%% initialization
x = initial_input(nb, img); % density-adaptive initial sampling
dt = DelaunayTri(x(:,1),x(:,2));
t = dt.Triangulation;

%% setting up boundary information
[B,H]=findBoundary(x, t);
isbnd = zeros(length(x),1);
isbnd(B,:) = 1;

%% extracting local information
neigh_tri = cell(length(x), 1); 
for tidx = 1:length(t)
    for j = 1:3
        idx = t(tidx,j);
        neigh_tri{idx}(end+1) = tidx;
    end   
end

embed = x;
newt = t;
figure
plot(embed(:,1),embed(:,2),'.','Color','black');
axis equal

%% initial setting for density integration
qx = zeros(width*height, 2); % convert pixels to points with coordinates
idx = 1;
for j=1:height
    for i = 1:width
        qx(idx, :) = [i-0.5, j-0.5];
        idx = idx+1;
    end
end
val = zeros(width*height, 1); % record grayscale values as density
idx = 1;
for j=1:height
    for i=1:width
        val(idx) = (256-double(img(height-j+1, i)))/256;
        idx = idx+1;
    end
end


%% Iteration
E_plot = [];
dt_Ite = 0; % number of topology optimization phase
d_t = CalDensity(embed, dt, img, qx, val); % calculate density of each triangle
areas = pdetrg(embed',newt');
capacity = areas'.*d_t;
var_cap = var(capacity)*(length(newt)-1);
fprintf('%dth dt generations;variance of capacity after DT = %f \n', dt_Ite,var_cap);

ite = 0; % total passes over each vertex
while dt_Ite < 5
    E = var_cap;
    E_last = -1;
    ite_local = 0; % passes over each vertex during geometry optimization phase
    %%%%%%%%%%%%%% geometry optimization phase %%%%%%%%%%%%%%%%%%%%%%
    while ite_local<5 && abs(E_last/E-1)>1e-2
        ite_local = ite_local+1;
        E_last = E;
        ite = ite+1;
        for idx = 1:length(x)
            if ~isbnd(idx)    % only iterate at interior points
                A = zeros(2,2);
                b = zeros(2,1);
                ai = zeros(length(neigh_tri{idx}),1);
                bi = zeros(length(neigh_tri{idx}),1);
                ci = zeros(length(neigh_tri{idx}),1);
                di = zeros(length(neigh_tri{idx}),1);
                for i = 1:length(neigh_tri{idx})
                    tidx = neigh_tri{idx}(i);
                    j = find(newt(tidx,:)==idx);
                    j1 = mod(j,3)+1;
                    j2 = mod(j1,3)+1;
                    p1 = embed(newt(tidx,j1),:);
                    p2 = embed(newt(tidx,j2),:);
                    ax = p1(1); ay = p1(2); bx = p2(1); by = p2(2);
                    ai(i) = ay-by;
                    bi(i) = bx-ax;
                    ci(i) = ax*by-ay*bx;   % Ai = ai*x + bi*y + ci
                    di(i) = d_t(tidx);
                end
                A(1,1) = sum((ai.*di-ai'*di/length(neigh_tri{idx})).^2);
                A(1,2) = sum((ai.*di-ai'*di/length(neigh_tri{idx})).*(bi.*di-bi'*di/length(neigh_tri{idx})));
                b(1) = -sum((ai.*di-ai'*di/length(neigh_tri{idx})).*(ci.*di-ci'*di/length(neigh_tri{idx})));
                A(2,1) = A(1,2);
                A(2,2) = sum((bi.*di-bi'*di/length(neigh_tri{idx})).^2);
                b(2) = -sum((bi.*di-bi'*di/length(neigh_tri{idx})).*(ci.*di-ci'*di/length(neigh_tri{idx})));
                A_inv = [A(2,2),-A(1,2);-A(2,1),A(1,1)]/(A(2,2)*A(1,1)-A(1,2)*A(2,1));              
                p_center = A_inv*b;
                embed(idx,1:2) = p_center';
            end
        end
%         figure
%         plot(embed(:,1),embed(:,2),'.');
%         axis equal      
        areas = pdetrg(embed',newt');
        capacity = areas'.*d_t;
        var_cap = var(capacity)*(length(newt)-1);
        E = var_cap;
        E_plot(end+1,:) = [E,0];
        fprintf('%dth iterations; Energy = %f;\n', ite, E);
    end
 
 %%%%%%%%%%% topology optimization phase %%%%%%%%%%%%%%%%%%%%%%%%
    dt = DelaunayTri(embed(:,1),embed(:,2));
    newt = dt.Triangulation;
    [B,H] = findBoundary(embed, newt);  % reset boundary information
    isbnd = zeros(length(x),1);
    isbnd(B,:) = 1;
    neigh_tri = cell(length(embed), 1);  % reset neighboring information
    for tidx = 1:length(newt)
        for j = 1:3
            idx = newt(tidx,j);
            neigh_tri{idx}(end+1) = tidx;
        end   
    end
    d_t = CalDensity(embed, dt, img, qx, val); % update densities
    areas = pdetrg(embed',newt');
    capacity = areas'.*d_t;
    var_cap = var(capacity)*(length(newt)-1);
    dt_Ite = dt_Ite+1;
    fprintf('%dth dt generations;variance of capacity after CDT = %f \n', dt_Ite,var_cap);
    E_plot(end+1,:) = [var_cap,1];
end

figure;
plot(E_plot);  
[i,j]=find(E_plot(:,2)==1);
hold on
plot(i,E_plot(i,1), '.', 'Color', 'r');  % energy curve
figure    % sampling pattern
plot(embed(:,1),embed(:,2),'.', 'Color','black');
axis equal