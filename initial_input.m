function x = initial_input(nv, img)  % generate initial sampling on image

height = size(img, 1);
width = size(img, 2);

x = [];
while size(x,1)<nv
    xx = rand(1,1);
    yy = rand(1,1);
    p_x = round(yy*height);
	p_y = round(xx*width);
	p_x = height-p_x;
	if (p_x == 0)
		p_x = 1;
    end
	if (p_y == 0)
		p_y = 1;
    end
	p = (256-double(img(p_x, p_y)))/256;
    
    r = rand(1,1);
    if (p >= r) 
      x(end+1,:) = [xx*width,yy*height,0];
    end
end