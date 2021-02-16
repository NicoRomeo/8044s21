function PhaseSpaceHO()
% Solves the harmonic oscillator, and displays the trajectory in phase
% space
%
% Uses a verlet integration scheme
%

close all

dt = .05;

q0 = 1;
p0 = +1.;

t_max = 3*2*pi+2;

t_arr= 1:dt:t_max;

q = zeros(length(t_arr),1);
p = zeros(size(q));

q(1) = q0;
p(1) = p0;


    function a = dhdq(q)
        a = sin(q);
    end
for i=1:length(t_arr(1:end-1))
    q(i+1) = q(i) + (p(i) - .5*dt*q(i))*dt;
    p(i+1) = p(i) - .5*dt*(q(i)+q(i+1));
    %p(i+1) = p(i) - .5*dt*(dhdq(q(i))+dhdq(q(i+1)));
end

figure
plot(q,p)
xlabel("q")
ylabel("p")
set(gcf,'color', "w")
daspect([1,1,1])

figure;
subplot(2,1,1)
plot(t_arr(1:floor(length(t_arr)))-1, q(1:floor(length(t_arr))))
xlabel("t")
ylabel("q")
subplot(2,1,2)
plot(t_arr(1:floor(length(t_arr)))-1, p(1:floor(length(t_arr))))
xlabel("t")
ylabel("p")
set(gcf,'color', "w")

"a"


%%%% make animation: two subplots, left marble in motion in a quadratic
%%%% potential, with velocity vector
%%%% right phase space motion, 

vid = VideoWriter(['SHO_phasespace', '.mp4'],'MPEG-4');
vid.FrameRate = 20;
vid.Quality = 100;
open(vid);
figure('units','pixels','position',[0 0 1920 1080])
set(gcf, 'color', 'w')

maxq = max(abs(q));
maxp = max(abs(p));
qrange = linspace(-1.1*maxq, 1.1*maxq,100);
for t=1:length(t_arr)
    clf
    subplot(1,2,1)
    hold on;
    plot(qrange, .5*qrange.^2, '-', 'LineWidth', 2, 'Color', 'black')
    scatter(q(t), .5*q(t) .^2,100,'red','filled')
    quiver(q(t), .5*q(t) .^2, .3*p(t), 0,'LineWidth', 2);
    
    xlim([-.3*maxp - maxq, maxq+.3*maxp])
    ylim([-.1*maxq, 1.2*.5*(1.1*maxq)^2])
    hold off
    box off
    axis off
%     yticks("","")
%     xticks("","")
    
    subplot(1,2,2)
    
    ts = max(1, t-10);
    %col = linspace(0,1,length(q(ts:t)));
    col = colormap(gray(length(q(ts:t))));
    hold on
    for i=0:(length(q(ts:t))-2)
        if length(q(ts:t)) < 5
            break
        end
        plot(q(t-i-1:t-i),p(t-i-1:t-i),'color',col(i+1,:), 'Linewidth', 2)
%     surface([q(ts:t)';q(ts:t)'], [p(ts:t)';p(ts:t)'], [col';col'],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
    end
    scatter(q(t),p(t), 100, 'red', 'filled')
    text(1.05*maxq, -.02* maxq, "q")
    text( -.05* maxq, 1.05*maxp, "p")
    text(-.05* maxq,-.05* maxq, "0")
    xline(0);
    yline(0);
    hold off;
    axis off
    box off
    daspect([1,1,1])
    xlim([-1.1*maxq, 1.1*maxq])
    ylim([-1.1*maxp,1.1*maxp])
    xlabel("q")
    ylabel("p")
    
    set(gcf, 'color', 'w')
    writeVideo(vid,getframe(gcf));
    t
end
close(vid)




end