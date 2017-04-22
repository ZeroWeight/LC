
close all
figure;
xlim([-400, 0]);
ylim([0, 200]);
set (gcf,'Position',[200,200,1200,600], 'color','w');
hold on
axis off


for i = 1:300
    if findobj == 0 
        break
    end
    cla;
    plot_cart(0.1 * i, 0.006 * i);
    pause(0.1);
    
end


function plot_cart(theta, phi)
    R = 20;
    l = 100;
    theta = -theta;
    wheel_Y = 20;
    wheel_X = theta * R - R;
    x1 = wheel_X; y1 = wheel_Y;
    x2 = x1 + R * sin(theta); y2 = y1 + R * cos(theta);
    x3 = x1 - l * sin(phi); y3 = y1 + l * cos(phi);
    draw_circle(R, wheel_X, wheel_Y);
    plot([x1, x3], [y1, y3], 'm', 'LineWidth',8);
    plot([x1, x2], [y1, y2], 'k','LineWidth',1);
    

end


function draw_circle(radius, x, y)
	alpha = 0: pi/20: 2*pi;
	circle_X = radius * cos(alpha) + x;
	circle_Y = radius * sin(alpha) + y;
	fill(circle_X, circle_Y, 'r');


end
