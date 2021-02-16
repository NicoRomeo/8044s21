function main_single()
% Main script to run hard disk simulations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
% parameters
sigma= .1;
N = 3;
Lx = 1;
Ly = 1;
T = 100;

% density of the system
eta = N*pi*sigma^2/(Lx*Ly)

% alternative defintion: set eta, N, fix sigma
eta = .2;
sigma = sqrt(eta *Lx*Ly/(N*pi))

% make animation?
anim = 1;

% run hard disk simulation
[traj_s, vels_s, times_s, traj_d, vels_d, times_d] = harddisk_single(Lx,Ly,N, sigma, T);
%harddisk_per(1,1,1)

    function plotTime(t)
        % plotting function to plot the state at given time index.
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

Ntani = 500;
[xis,vis] = interpolate_disks(traj_s, vels_s, times_s, Ntani);

[xid,vid] = interpolate_disks(traj_d, vels_d, times_d, Ntani);

% % Compute density observable using formula from Krauth p.87, eq 2.3
% Y_dis = linspace(0, Ly, 30);
% dY = Y_dis(2) - Y_dis(1);
% rho_y = zeros(size(Y_dis));
% 
% final_t = max(times);
% 
% for t=1:Ntani  
%     for i=1:N
%         a = ceil(xi(2,i,t)/dY);
%         rho_y(a) = rho_y(a) + 1/abs(vi(2,i,t));
%     end
% end
% rho_y = rho_y / final_t;
% 
% figure
% plot(Y_dis, rho_y)
% box on
% xlabel("y")
% ylabel("\rho_y")
% set(gcf, 'color', 'w')
% 
% % speed distribution
% 
% speed = sqrt(vels(1,:,:).^2 + vels(2,:,:).^2);
% 
% vbar = mean(speed(:));
% vmax = max(speed(:));
% betam = 2/(pi*vbar^2);
% vrange = 0:(vmax/100):vmax;
% figure
% hold on
% histogram(speed, 'Normalization', 'pdf')
% plot(vrange, betam*vrange .* exp( - betam * vrange .^2/2));
% hold off
% xlabel("speed v")
% ylabel("p(v)")



%%%% make animation if anim 

if anim
    
    video = VideoWriter(['harddiskN_', num2str(N), 'singlecomp2.mp4'],'MPEG-4');
    video.FrameRate = 20;
    video.Quality = 100;
    open(video);
    f1 = figure;
        th = 0:pi/50:2*pi;
        %th = 0:pi/50:pi;
        for t=1:Ntani
                clf
                hold on 
                for k=1:N
                %rectangle('Position',[traj(1,k,t)-sigma, traj(2,k,t)-sigma,sigma2,sigma2], 'Curvature',[1,1])
                     xunit = sigma *cos(th) + xis(1,k,t);
                     yunit = sigma * sin(th) + xis(2,k,t);
                     plot(xunit, yunit, 'Color', '#0072BD', 'LineWidth', 2)
                    
                    %yunit = [sigma * sin(th) + xi(2,k,t); -sigma * sin(th) + xi(2,k,t)]';
                    %area(xunit, yunit, 'FaceColor', '#0072BD', 'EdgeColor', 'none')
                end
                quiver(xis(1,:,t),xis(2,:,t), .6*vis(1,:,t),.6*vis(2,:,t),'Color', 'black', 'LineWidth', 2 )
                for k=1:N
                %rectangle('Position',[traj(1,k,t)-sigma, traj(2,k,t)-sigma,sigma2,sigma2], 'Curvature',[1,1])
                     xunit = sigma *cos(th) + xid(1,k,t);
                     yunit = sigma * sin(th) + xid(2,k,t);
                     plot(xunit, yunit, 'Color', 'r', 'LineWidth', 2)
                    
                    %yunit = [sigma * sin(th) + xi(2,k,t); -sigma * sin(th) + xi(2,k,t)]';
                    %area(xunit, yunit, 'FaceColor', '#0072BD', 'EdgeColor', 'none')
                end
            quiver(xid(1,:,t),xid(2,:,t), .6*vid(1,:,t),.6*vid(2,:,t),'Color', 'black', 'LineWidth', 2 )
            hold off
            xlim([0,Lx])
            ylim([0,Ly])
            daspect([1,1,1])
            xticks([])
            yticks([])
            box on
            set(gcf, 'color', 'w')
            writeVideo(video,getframe(gcf));
            
            t
        end
    close(video)

%     video = VideoWriter(['harddiskN_', num2str(N), 'single.mp4'],'MPEG-4');
%     video.FrameRate = 20;
%     video.Quality = 100;
%     open(video);
%     f1 = figure;
%         th = 0:pi/50:2*pi;
%         %th = 0:pi/50:pi;
%         for t=1:Ntani
%                 clf
%                 hold on 
%                 for k=1:N
%                 %rectangle('Position',[traj(1,k,t)-sigma, traj(2,k,t)-sigma,sigma2,sigma2], 'Curvature',[1,1])
%                      xunit = sigma *cos(th) + xis(1,k,t);
%                      yunit = sigma * sin(th) + xis(2,k,t);
%                      plot(xunit, yunit, 'Color', '#0072BD', 'LineWidth', 2)
%                     
%                     %yunit = [sigma * sin(th) + xi(2,k,t); -sigma * sin(th) + xi(2,k,t)]';
%                     %area(xunit, yunit, 'FaceColor', '#0072BD', 'EdgeColor', 'none')
%                 end
%             quiver(xis(1,:,t),xis(2,:,t), .6*vis(1,:,t),.6*vis(2,:,t),'Color', 'black', 'LineWidth', 2 )
%             hold off
%             xlim([0,Lx])
%             ylim([0,Ly])
%             daspect([1,1,1])
%             xticks([])
%             yticks([])
%             box on
%             set(gcf, 'color', 'w')
%             writeVideo(video,getframe(gcf));
%             
%             t
%         end
%     close(video)
%     
%     video = VideoWriter(['harddiskN_', num2str(N), 'double.mp4'],'MPEG-4');
%     video.FrameRate = 20;
%     video.Quality = 100;
%     open(video);
%     f1 = figure;
%         th = 0:pi/50:2*pi;
%         %th = 0:pi/50:pi;
%         for t=1:Ntani
%                 clf
%                 hold on 
%                 for k=1:N
%                 %rectangle('Position',[traj(1,k,t)-sigma, traj(2,k,t)-sigma,sigma2,sigma2], 'Curvature',[1,1])
%                      xunit = sigma *cos(th) + xid(1,k,t);
%                      yunit = sigma * sin(th) + xid(2,k,t);
%                      plot(xunit, yunit, 'Color', '#0072BD', 'LineWidth', 2)
%                     
%                     %yunit = [sigma * sin(th) + xi(2,k,t); -sigma * sin(th) + xi(2,k,t)]';
%                     %area(xunit, yunit, 'FaceColor', '#0072BD', 'EdgeColor', 'none')
%                 end
%             quiver(xid(1,:,t),xid(2,:,t), .6*vid(1,:,t),.6*vid(2,:,t),'Color', 'black', 'LineWidth', 2 )
%             hold off
%             xlim([0,Lx])
%             ylim([0,Ly])
%             daspect([1,1,1])
%             xticks([])
%             yticks([])
%             box on
%             set(gcf, 'color', 'w')
%             writeVideo(video,getframe(gcf));
%             
%             t
%         end
%     close(video)
    
end
end