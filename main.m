function main()

sigma= .05;
N = 20;
Lx = 1;
Ly = 1;
[traj, vels, times] = harddisk(Lx,Ly,N, sigma);
%harddisk_per(1,1,1)
    function plotTime(t)
        figure
        th = 0:pi/50:2*pi;
        hold on 
        for k=1:N
            %rectangle('Position',[traj(1,k,t)-sigma, traj(2,k,t)-sigma,sigma2,sigma2], 'Curvature',[1,1])
            xunit = sigma *cos(th) + traj(1,k,t);
            yunit = sigma * sin(th) + traj(2,k,t);
            plot(xunit, yunit, 'Color', 'b')
        end
        quiver(traj(1,:,t),traj(2,:,t), vels(1,:,t),vels(2,:,t))
        hold off
        xlim([0,Lx])
        ylim([0,Ly])
        daspect([1,1,1])
        box on
        set(gcf, 'color', 'w')
    end



    function [Xi, Vi] = interpolate_disks(X, V, T, Nq)
        % interpolates the position and velocities of disks 
        Tq = linspace(min(T), max(T), Nq);
        Xi = zeros(2,N,Nq);
        Vi = zeros(2,N,Nq);
        for m=1:N
            Xi(1,m,:) = interp1(T, squeeze(X(1,m,:)), Tq);
            Xi(2,m,:) = interp1(T, squeeze(X(2,m,:)), Tq);
            Vi(1,m,:) = interp1(T, squeeze(V(1,m,:)), Tq, 'previous');
            Vi(2,m,:) = interp1(T, squeeze(V(2,m,:)), Tq, 'previous');
        end
    end

Ntani = 1000;
[xi,vi] = interpolate_disks(traj, vels, times, Ntani);

%%%% make animation 

    vid = VideoWriter(['harddisk20', '.mp4'],'MPEG-4');
    vid.FrameRate = 20;
    vid.Quality = 100;
    open(vid);
    f1 = figure;
        th = 0:pi/50:2*pi;
        for t=1:Ntani
                clf
                hold on 
                for k=1:N
                %rectangle('Position',[traj(1,k,t)-sigma, traj(2,k,t)-sigma,sigma2,sigma2], 'Curvature',[1,1])
                    xunit = sigma *cos(th) + xi(1,k,t);
                    yunit = sigma * sin(th) + xi(2,k,t);
                    plot(xunit, yunit, 'Color', 'b')
                end
            quiver(xi(1,:,t),xi(2,:,t), vi(1,:,t),vi(2,:,t))
            hold off
            xlim([0,Lx])
            ylim([0,Ly])
            daspect([1,1,1])
            box on
            set(gcf, 'color', 'w')
            writeVideo(vid,getframe(gcf));
            
            t
        end
    close(vid)
end