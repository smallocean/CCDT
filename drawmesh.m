function h=drawmesh(T, X, B, color)
    if nargin<4
        color = 'None'; %interp
    end
    figure;
    if size(X,2)>3
        if range( X(:,4) )<1e-10
            X(:,4) = 0;
        end
        h=trisurf(T, X(:,1), X(:,2), X(:,3), X(:,4), 'FaceColor', color,'FaceAlpha',0);
    elseif size(X,2)>2
        h=trisurf(T, X(:,1), X(:,2), X(:,3), zeros(size(X,1),1), 'FaceColor', color, 'FaceAlpha',0 );
    else
        h=trisurf(T, X(:,1), X(:,2), 'FaceColor', color, 'FaceAlpha',0);
    end
    view(2);
    axis equal;
    hold on;

    if nargin>2 && ~isempty(B)
        tB = [B B(1)];
        plot( X(tB, 1), X(tB, 2), 'b', 'linewidth', 2);
    end
    
    axis off;

