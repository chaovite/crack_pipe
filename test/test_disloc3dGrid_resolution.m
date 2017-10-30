% test the effect of grid resolution on static surface displacement.
%
% In this test, I test two pressure distribution: 2D gaussian (zero-center) or uniform.
%
source_dir = '../source';
addpath(genpath(source_dir));
%% uniform pressure.
% fault length, high resolution.
NX = [4, 8, 16, 32, 64];
sol_uni = cell(length(NX),1);
for i = 1: length(NX)
        xc = 0;
        yc = 0;
        nx = NX(i);

        zc = 1200;
        
        L   = 1500;
        W  = 1500;
        ny  = nx;
        
        dx  = W/nx;
        dy = L/ny;

        x    = [1/2:1: (nx-1/2)]*dx -1/2*dx;
        y    = [1/2:1: (ny-1/2)]*dy -1/2*dy;

        [X, Y] = meshgrid(x, y);

        mu = 10e9;
        nu  = 0.25;

        % uniform pressure.
        p    =  1e6*ones(size(X));
        % observation points.
        nq = 20;
        xq = linspace(-1,1,nq)*W ;
        yq = linspace(-1,1,nq)*W;
        [Xq, Yq] = meshgrid(xq, yq);
        Zq = zeros(nq, nq);

        % test for strike = 0, dip = 30.
        strike = 0; dip =50;
        % reference point.
        depth = zc + W/2*sind(dip);
        east   = xc + W/2*cosd(dip)*cosd(strike);
        north = yc - W/2*cosd(dip)*sind(strike);

        tic;
        dsl_grid = disloc3dGrid(xc, yc, zc, strike, dip, X, Y, mu, nu);
        toc;
        
        U    = dsl_grid.eval_disp_p(p(:), Xq(:), Yq(:));
        % compare different 
        subplot(1,length(NX),i)
        pcolor(Xq/W, Yq/W, reshape(U(3,:), nq, nq));
        xlabel('x/W');
        ylabel('y/W');
        title(sprintf('nx=%d', nx));
        shading interp
        caxis([0 0.04])
        sol_uni{i, 1} = U;
end
colorbar;
suptitle(sprintf('strike = %d, dip = %d',strike, dip));

%% gaussian pressure.
% fault length, high resolution.
NX = [4, 8, 16, 32, 64];
D   = [W/2, W/4, W/8];
climit = [22e-3, 13e-3, 4.5e-3];
sol_gauss = cell(length(NX), length(D));
for j = 1: length(D)
    for i = 1: length(NX)
        xc = 0;
        yc = 0;
        nx = NX(i);
        zc = 1200;
        
        L   = 1500;
        W  = 1500;
        ny  = nx;
        
        dx  = W/nx;
        dy = L/ny;

        x    = [1/2:1: (nx-1/2)]*dx -1/2*dx;
        y    = [1/2:1: (ny-1/2)]*dy -1/2*dy;

        [X, Y] = meshgrid(x-mean(x), y-mean(y));

        [Xgrd, Ygrd] = meshgrid(x, y);

        mu = 10e9;
        nu  = 0.25;

        % uniform pressure.
        p    =  1e6*exp(-(X.^2+Y.^2)/(2*D(j)^2));
        % observation points.
        nq = 20;
        xq = linspace(-1,1,nq)*W ;
        yq = linspace(-1,1,nq)*W;
        [Xq, Yq] = meshgrid(xq, yq);
        Zq = zeros(nq, nq);

        % test for strike = 0, dip = 30.
        strike = 0; dip = 0;
        % reference point.
        depth = zc + W/2*sind(dip);
        east   = xc + W/2*cosd(dip)*cosd(strike);
        north = yc - W/2*cosd(dip)*sind(strike);

        tic;
        dsl_grid = disloc3dGrid(xc, yc, zc, strike, dip, Xgrd, Ygrd, mu, nu);
        toc;
        
        U    = dsl_grid.eval_disp_p(p(:), Xq(:), Yq(:));
        % compare different 
        figure(1);
        subplot(length(D), length(NX), length(NX)*(j-1) + i);
        pcolor(Xq/W, Yq/W, reshape(U(3,:)*1e3, nq, nq));
        xlabel('x/W');
        ylabel('y/W');
        title(sprintf('nx=%d, D=%4.1f', nx, D(j)));
        shading interp
        caxis([0 climit(j)*1e3]);
        sol_gauss{i, j} = U;
    end
    colorbar;
    figure(2);
    subplot(length(D), 1, j);
    pcolor(X, Y, p);
    shading interp;
    drawnow;
end
suptitle(sprintf('strike = %d, dip = %d',strike, dip));

%% biharmonic distribution pressure.
% fault length, high resolution.
NX = [4, 8, 16, 32, 64];
W  = 1500;
lam   = [2*W, W, W/2];
climit = [15, 0.3, 0.11];
sol_bihormonic = cell(length(NX), length(lam));
for j = 1: length(lam)
    for i = 1: length(NX)
        xc = 0;
        yc = 0;
        nx = NX(i);
        zc = 1200;
        
        L   = 1500;
        W  = 1500;
        ny  = nx;
        
        dx  = W/nx;
        dy = L/ny;

        x    = [1/2:1: (nx-1/2)]*dx -1/2*dx;
        y    = [1/2:1: (ny-1/2)]*dy -1/2*dy;

        [X, Y] = meshgrid(x-mean(x), y-mean(y));

        [Xgrd, Ygrd] = meshgrid(x, y);

        mu = 10e9;
        nu  = 0.25;

        % uniform pressure.
        p    =  1e6*cos(2*pi/lam(j)*X).*cos(2*pi/lam(j)*Y);
        % observation points.
        nq = 20;
        xq = linspace(-1,1,nq)*W ;
        yq = linspace(-1,1,nq)*W;
        [Xq, Yq] = meshgrid(xq, yq);
        Zq = zeros(nq, nq);

        % test for strike = 0, dip = 30.
        strike = 0; dip = 0;
        % reference point.
        depth = zc + W/2*sind(dip);
        east   = xc + W/2*cosd(dip)*cosd(strike);
        north = yc - W/2*cosd(dip)*sind(strike);

        tic;
        dsl_grid = disloc3dGrid(xc, yc, zc, strike, dip, Xgrd, Ygrd, mu, nu);
        toc;
        
        U    = dsl_grid.eval_disp_p(p(:), Xq(:), Yq(:));
        % compare different 
        figure(1);
        subplot(length(lam), length(NX), length(NX)*(j-1) + i);
        pcolor(Xq/W, Yq/W, reshape(U(3,:)*1e3, nq, nq));
        xlabel('x/W');
        ylabel('y/W');
        title(sprintf('nx=%d, lam=%2.1f W', nx, lam(j)/W));
        shading interp
%         caxis([0 climit(j)*1e3]);
        sol_bihormonic{i, j} = U;
    end
    colorbar;
    figure(2);
    subplot(length(lam), 1, j);
    pcolor(X, Y, p);
    shading interp;
    drawnow;
end
suptitle(sprintf('strike = %d, dip = %d',strike, dip));


