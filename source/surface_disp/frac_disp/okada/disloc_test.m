% test the dislocation code with an horizontal open crack
%
%  The dislocation model is specified as: length, width,
%  depth, dip, strike, east, north, strike-slip, dip-slip,
%  and opening.  The coordinates (depth, east, and north)
%  specify a point at the middle of the bottom edge of
%  the fault for positive dips and the middle of the top
%  edge for negative dips.

x=linspace(-1e3,1e3,1000);
y=linspace(-1e3,1e3,1000);
[X Y]=meshgrid(x,y');
%%
stations=[X(:)';Y(:)'];
model=[100 100 1000 0 0 50 0 0 0 0.1]';
nu=0.25;
u=disloc(model,stations,nu);
%%
surf(X,Y,reshape(u(3,:),1000,1000));
shading flat;
set(gca,'fontsize',14);
xlabel('East','fontsize',18);ylabel('North','fontsize',18);
zlabel('Vertical displacement');
colorbar;
shg
%  view([0 0 1])

% quiver(X(:)',Y(:)',u(1,:),u(2,:),10);shg;xlim([-50,50]);ylim([-50,50]);
% view([0 0 1])


