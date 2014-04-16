function [embed, newt]=ccdt(nb) % nb = number of sampling points

%% initial triangulation
x = rand(nb,2)*100;
x(:,3) = 0;
dt = DelaunayTri(x(:,1), x(:,2));
t = dt.Triangulation;

%% setting up initial boundary information
[B,H] = findBoundary(x, t);  
isbnd = zeros(length(x),1);
isbnd(B,:) = 1;

%% extracting local information
neigh_tri = cell(length(x), 1); % neighboring triangles incident to each vertex
for tidx = 1:length(t)
    for j = 1:3
        idx = t(tidx,j);
        neigh_tri{idx}(end+1) = tidx;
    end   
end

embed = x;
newt = t;
close all;
drawmesh(newt, embed, B, 'white');

%% Iteration
E_plot = [];
dt_Ite = 0;  % number of topology optimization phase
d_t = ones(length(t),1);  % energy after Delaunay triangulation
while dt_Ite < 15
    ite = 0; % passes over each vertex during geometry optimization phase
    areas = pdetrg(embed', newt');
    E = var(areas')*(length(newt)-1);
    E_last = -1;  
    %%%%%%%%%%%%%% Geometry optimization phase %%%%%%%%%%%%%%%%%%
    while abs(E_last/E-1)>1e-1
        E_last = E;
        ite = ite+1;  
        for idx = 1:length(x)
            if ~isbnd(idx)     % only iterate at interior points
                A = zeros(2,2);
                b = zeros(2,1);
                ai = zeros(length(neigh_tri{idx}),1);
                bi = zeros(length(neigh_tri{idx}),1);
                ci = zeros(length(neigh_tri{idx}),1);
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
                end
          %%%%% linear 2*2 system AX=b for minimizing variance of area %%%
                A(1,1) = sum((ai-sum(ai)/length(neigh_tri{idx})).^2);
                A(1,2) = sum((ai-sum(ai)/length(neigh_tri{idx})).*(bi-sum(bi)/length(neigh_tri{idx})));
                b(1) = -sum((ai-sum(ai)/length(neigh_tri{idx})).*(ci-sum(ci)/length(neigh_tri{idx})));
                A(2,1) = A(1,2);
                A(2,2) = sum((bi-sum(bi)/length(neigh_tri{idx})).^2);
                b(2) = -sum((bi-sum(bi)/length(neigh_tri{idx})).*(ci-sum(ci)/length(neigh_tri{idx})));
                A_inv = [A(2,2),-A(1,2);-A(2,1),A(1,1)]/(A(2,2)*A(1,1)-A(1,2)*A(2,1));              
                p_opt = A_inv*b;
                embed(idx,1:2) = p_opt';
            end
        end  
%         drawmesh(newt, embed, B, 'white');
        areas = pdetrg(embed', newt');
        E = var(areas')*(length(newt)-1);
        E_plot(end+1,:) = [E,0];
        fprintf('%dth iterations; Energy = %f;\n', ite, E);
    end
 %%%%%%%%%%% Topology optimization phase %%%%%%%%%%%%%%%%%%%%%%%%
    dt = DelaunayTri(embed(:,1),embed(:,2)); 
    newt = dt.Triangulation;
    [B,H] = findBoundary(embed, newt);  % reset boundary information
    isbnd = zeros(length(x),1);
    isbnd(B,:) = 1;
%     drawmesh(newt, embed, B, 'white');
    neigh_tri = cell(length(embed), 1); % reset neighboring information
    for tidx = 1:length(newt)
        for j = 1:3
            idx = newt(tidx,j);
            neigh_tri{idx}(end+1) = tidx;
        end   
    end
    dt_Ite = dt_Ite+1;
    areas = pdetrg(embed', newt');
    E_dt = var(areas')*(length(newt)-1);
    E_plot(end+1,:)=[E_dt,1];
    fprintf('%dth dt generations; Energy = %f\n', dt_Ite, E_dt);
end

figure;
plot(E_plot(:,1));
[i,j]=find(E_plot(:,2)==1);
hold on
plot(i,E_plot(i,1), '*', 'Color', 'r');  % energy curve
drawmesh(newt, embed, B, 'white'); % final sampling triangulation
figure
plot(embed(:,1),embed(:,2),'.'); % sampling pattern
axis equal