function [B,H]=findBoundary(x, t)
    nf = size(t, 1);
    nv = size(x, 1);
    bedges=sparse(nv, nv);
    for i=1:nf
        i1 = t(i,1); i2 = t(i,2); i3 = t(i,3);
        if bedges(i1, i2)        pause; disp('inconsistent triangulation!');
        elseif bedges(i2, i1)    bedges(i2, i1) = 0;
        else                     bedges(i1, i2) = i;
        end
        
        if bedges(i2, i3)        pause; disp('inconsistent triangulation!');
        elseif bedges(i3, i2)    bedges(i3, i2) = 0;
        else                     bedges(i2, i3) = i;
        end
        
        if bedges(i3, i1)        pause; disp('inconsistent triangulation!');
        elseif bedges(i1, i3)    bedges(i1, i3) = 0;
        else                     bedges(i3, i1) = i;
        end
    end

    H = [];
    B = [];
    nh=0;
    while true
        [e1 e2] = find(bedges,1);
        if isempty(e1)
            break;
        end

        bedges(e1, e2) = 0;
        tB = [e1 e2];
        while tB(end)~=tB(1)
            e2 = find( bedges(tB(end),:), 1 );
            if ~size(e2)
                disp('incomplete boundary!');
                bedges=[];
                pause;
                break;
            end

            bedges(tB(end), e2) = 0;
            tB = [tB e2];
        end

        if length(B)<length(tB)
            if ~isempty(B)
                nh = nh+1;
                H(nh).id = B;
            end
            B = tB;
        else
            nh = nh+1;
            H(nh).id = tB;
        end
%         e1 = x(tB(1),:)-x(tB(2),:);
%         e2 = x(tB(1),:)-x(tB(3),:);
%         if e1(2)*e2(1)-e1(1)*e2(2)>0
%             B = [B tB];
%             antiClock = antiClock+1;
%         else
%             H = [H tB];
%         end
    end

%     if antiClock~=1
%         [B, H] = deal(H, B);
%     end
    B = B(1:end-1);
end