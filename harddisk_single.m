function [traj_s, vels_s, times_s, traj_d, vels_d, times_d] =harddisk_single(Lx, Ly, N, sigma, T)
% simulates Monte-Carlo dynamis of a hard disk gas in a periodic system.
%
% for Reference, see:
% W. Krauth, Statistical Mechanics: Algorithms and Computation, OUP
%
% outputs position of disks at each recorded timesteps

Ldim = [Lx, Ly];
%sigma = .1; %radius of the disk
sigma2 = sigma^2;

Nt = T; % number of events to simulate


%%%%%% initial condition st

% X = rand(2,N); % pick centers in the interior [sigma, L-sigma]
% X(1,:) = (Lx-2*sigma) * X(1,:) + sigma;
% X(2,:) = (Ly-2*sigma) * X(2,:) + sigma;
% count = 1;
% while min(pdist(X')) < 2*sigma % make sure spheres don't overlap
%     X = rand(2,N);
%     X(1,:) = (Lx-2*sigma) * X(1,:) + sigma;
%     X(2,:) = (Ly-2*sigma) * X(2,:) + sigma;
%     count = count + 1 ;
%     if count > 30
%         break
%     end
% end

% random deposition initial condition
X = zeros(2,N);
X(1,1) = (Lx-2*sigma) *rand() + + sigma;
X(2,1) = (Ly-2*sigma) *rand() + + sigma;
for k=2:N
    X(1,k) = (Lx-2*sigma) *rand() + sigma;
    X(2,k) = (Ly-2*sigma) *rand() + sigma;
    
    while min(sum((X(:,1:(k-1)) - X(:,k)).^2,1)) < 4*sigma2
        X(1,k) = (Lx-2*sigma) *rand() + sigma;
        X(2,k) = (Ly-2*sigma) *rand() + sigma;
    end
end

X = single(X);

theta = 2*pi*rand(1,N);
v0 = .1;
V = zeros(2,N);
V(1,:) = single(v0 * cos(theta));
V(2,:) = single(v0 * sin(theta));


% test plot initial contions
figure
th = 0:pi/50:2*pi;
hold on 
for k=1:N
    %rectangle('Position',[X(1,k)-sigma, X(2,k)-sigma,sigma2,sigma2], 'Curvature',[1,1])
    xunit = sigma *cos(th) + X(1,k);
    yunit = sigma * sin(th) + X(2,k);
    plot(xunit, yunit, 'Color', 'b')
end
quiver(X(1,:),X(2,:), V(1,:),V(2,:))
hold off
xlim([0,Lx])
ylim([0,Ly])
daspect([1,1,1])
box on
set(gcf, 'color', 'w')


%%%% function definitions
    function tp = pairtime(x, vx, y, vy)
        % computes time before next collision of 2 particles
        dx = y-x;
        dv = vy - vx;
        dxdv = dot(dx,dv);
        Tau = dxdv^2 - dot(dv,dv)*(dot(dx,dx) - 4*sigma2);
        if Tau >0 && dxdv <0
            tp = - (dxdv + sqrt(Tau))/dot(dv,dv);
        else
            tp = Inf;
        end
    end

    function tw = walltime(x,vx)
        % computes time before next wall collision
        if vx(1) > 0
            t1 = (Lx - sigma*.9 - x(1))/vx(1);
        else if vx(1) < 0
                t1 =  (sigma*.9 - x(1))/vx(1);
            else
                t1 = Inf;
            end
        end
        if vx(2) > 0
            t2 = (Ly -sigma*.9 - x(2))/vx(2);
        else if vx(2) < 0
                t2 =  (sigma*.9 - x(2))/vx(2);
            else
                t2 = Inf;
            end
        end
        tw = min(t1,t2);
    end

    function vnew = wallcollision(x, v)
        % handles wall collisions, returns new velocity vector
        vnew = v;
        if x(1) <= sigma || Lx - x(1) <= sigma
            vnew(1) = -v(1);
        end
        if x(2) <= sigma || Ly - x(2) <= sigma
            vnew(2) = -v(2);
        end
    end

    function [vn1, vn2] = paircollision(x1,x2,v1,v2)
        % handles pair collisions, returns new velocity vectors
        dx = x2 - x1;
        e_perp = dx/norm(dx);
        dv = v2 - v1;
        vn2 = v2 - e_perp * (dot(dv, e_perp));
        vn1 = v1 + e_perp * (dot(dv, e_perp));
    end


    function [Xn, Vn, tn] = update(X, V, t)
        % updates situation to next timestep
        times_coll = Inf*ones(N);
        times_wall = Inf*ones(N,1);
        for i=1:(N-1)
            for j=i+1:N
                if i~=j
                times_coll(i,j) = pairtime(X(:,i), V(:,i), X(:,j), V(:,j));
                end
            end
            times_wall(i) = walltime(X(:,i),V(:,i));
        end
        times_wall(end) = walltime(X(:,end),V(:,end));
        
        [mintw, indw] = min(times_wall);
        mintc = min(min(times_coll));
        tnext = min(mintw,mintc);
        
        Xn = X + tnext * V;
        Vn = V;
        if mintw < mintc
            Vn(:,indw) = wallcollision(Xn(:,indw), V(:,indw));
        else
            [ind1,ind2] = find(times_coll== mintc);
            [vn1,vn2] = paircollision(Xn(:,ind1),Xn(:,ind2),V(:,ind1),V(:,ind2));
            Vn(:,ind1) = vn1(:);
            Vn(:,ind2) = vn2(:);
        end
        tn = t + tnext;  
    end

    

    traj_s = zeros(2,N,Nt);
    vels_s = zeros(2,N,Nt);
    times_s = zeros(Nt,1);
    
    traj_d = zeros(2,N,Nt);
    vels_d = zeros(2,N,Nt);
    times_d = zeros(Nt,1);
    
    traj_s(:,:,1) = single(X(:,:));
    vels_s(:,:,1) = single(V(:,:));
    
    traj_d(:,:,1) = double(X(:,:));
    vels_d(:,:,1) = double(V(:,:));
    tic
    for t=2:Nt
        [Xn, Vn, tn] = update(squeeze(traj_s(:,:,t-1)), squeeze(vels_s(:,:,t-1)), times_s(t-1));
        traj_s(:,:,t) = single(Xn(:,:));
        vels_s(:,:,t) = single(Vn(:,:));
        times_s(t) = single(tn);
        
        [Xn, Vn, tn] = update(squeeze(traj_d(:,:,t-1)), squeeze(vels_d(:,:,t-1)), times_d(t-1));
        traj_d(:,:,t) = double(Xn(:,:));
        vels_d(:,:,t) = double(Vn(:,:));
        times_d(t) = double(tn);
    end
    toc
    
    
        


end