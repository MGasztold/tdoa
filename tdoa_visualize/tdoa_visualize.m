% Program to plot the hyperbola
% y^2/a^2 - x^2/b^2 = 1
% The hyperbolae are open up/down, so that x is the independent variable
% for plotting. (Using the form x^2/a^2 - y^2/b^2 = 1 requires that y be
% the independent variable, which is awkward programming-wise.)
clc
close all

setenv("GNUTERM", "X11")
graphics_toolkit("gnuplot")

clear all% all variables
c=3E8; % velociy of light
figure(1), hold off % start a new figure
set(gca,'FontSize',10) % adjust fontsize
xmax = 40; ymax = 30;
x = linspace(-xmax,xmax,4001); % array of x values for plot
N=4;
x_stations = [-10, 10, 0, -10];
y_stations = [-10, -10, 10, 15];
plot(x_stations,y_stations, 'rs', 'Color', 'g','MarkerSize', 12);
hold on;
axis equal
axis([-xmax xmax -ymax ymax]) % specify axis limits
xlabel('x')
ylabel('y')
% Add axes
plot([0 0],[-ymax ymax],'k') % y axis (black line - ’k’)
plot([-xmax xmax],[0 0],'k') % x axis
%% TAG coodinates:
T = [10 1];
plot(T(1), T(2), '*', 'MarkerSize', 25, 'Color', 'm')
%% time of flight from TAG to each anchor
times_ = zeros(1,N);
for i=1:N
    times_(i) = sqrt((T(1) - x_stations(i))^2 + (T(2) - y_stations(i))^2)/c;
end
times_
%% normalized time differences of arrival due to anchor 0
dtimes_ = zeros(1,N);
for i=1:N
    if i<N
        dtimes_(i) = abs(times_(i) - times_(i+1));
    else
        dtimes_(i) = abs(times_(i) - times_(1));
    end
end
dtimes_
%% compute hyperbolas for each anchors pair and plot them
for i=1:N
    x = linspace(-xmax,xmax,1001); % array of x values for plot
        % compute middle point between anchors
        if i<N
            if x_stations(i) <= x_stations(i+1) && y_stations(i) < y_stations(i+1)
                xoffst = x_stations(i) + sqrt((x_stations(i) - x_stations(i+1))^2)/2;
                yoffst = y_stations(i) + sqrt((y_stations(i) - y_stations(i+1))^2)/2;
            elseif x_stations(i) <= x_stations(i+1) && y_stations(i) >= y_stations(i+1)
                xoffst = x_stations(i) + sqrt((x_stations(i) - x_stations(i+1))^2)/2;
                yoffst = y_stations(i+1) + sqrt((y_stations(i) - y_stations(i+1))^2)/2;
            elseif x_stations(i) > x_stations(i+1) && y_stations(i) < y_stations(i+1)
                xoffst = x_stations(i+1) + sqrt((x_stations(i) - x_stations(i+1))^2)/2
                yoffst = y_stations(i) + sqrt((y_stations(i) - y_stations(i+1))^2)/2
            elseif x_stations(i) > x_stations(i+1) && y_stations(i) >= y_stations(i+1)
                xoffst = x_stations(i+1) + sqrt((x_stations(i) - x_stations(i+1))^2)/2;
                yoffst = y_stations(i+1) + sqrt((y_stations(i) - y_stations(i+1))^2)/2;
            end
        else
            if x_stations(i) <= x_stations(1) && y_stations(i) < y_stations(1)
                xoffst = x_stations(i) + sqrt((x_stations(i) - x_stations(1))^2)/2;
                yoffst = y_stations(i) + sqrt((y_stations(i) - y_stations(1))^2)/2;
            elseif x_stations(i) <= x_stations(1) && y_stations(i) >= y_stations(1)
                xoffst = x_stations(i) + sqrt((x_stations(i) - x_stations(1))^2)/2;
                yoffst = y_stations(1) + sqrt((y_stations(i) - y_stations(1))^2)/2;
            elseif x_stations(i) > x_stations(1) && y_stations(i) < y_stations(1)
                xoffst = x_stations(1) + sqrt((x_stations(i) - x_stations(1))^2)/2
                yoffst = y_stations(i) + sqrt((y_stations(i) - y_stations(1))^2)/2
            elseif x_stations(i) > x_stations(1) && y_stations(i) >= y_stations(1)
                xoffst = x_stations(1) + sqrt((x_stations(i) - x_stations(1))^2)/2;
                yoffst = y_stations(1) + sqrt((y_stations(i) - y_stations(1))^2)/2;
            end
        end
        %% compute the angle at which hyperbolas are tilted
        theta(i) = atan(abs(xoffst-x_stations(i))/(abs(yoffst-y_stations(i))));
        if i<N
            theta(i) = 2*pi - theta(i)
        else 
            theta(i) = pi-theta(i)
        end
        %% compute hyperbola coefficients
        a = c*dtimes_(i)/2;
        if i<N
            odleglosc_czasowa_miedzy_stacjami = sqrt((x_stations(i) - x_stations(i+1))^2 + (y_stations(i) - y_stations(i+1))^2)/c;
        else
            odleglosc_czasowa_miedzy_stacjami = sqrt((x_stations(i) - x_stations(1))^2 + (y_stations(i) - y_stations(1))^2)/c;
        end
        d = c*0.5*abs(odleglosc_czasowa_miedzy_stacjami);
        b = sqrt(d^2 - a^2);
        
        %% draw hyperbola
        y = sqrt(((x.^2)./(b^2)+1).*a^2); % corresponding y values
        [x, y] = xfm1(x,y,theta(i),xoffst,yoffst);
        %if i==N
            plot(x,y)
            hold on;
        %else
            [x, y] = xfm1(x,y,0,-2*xoffst,-2*yoffst);
            plot(-x,-y, 'r') % Plot other half of hyperbola
        %end
        axis equal
        grid on
        axis([-xmax xmax -ymax ymax]) % specify axis limits
        xlabel('x')
        ylabel('y')
        set(gcf, 'color', 'white');
end

%%

disp('press return to continue') 
pause() 

