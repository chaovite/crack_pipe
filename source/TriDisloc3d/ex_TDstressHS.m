% ex_TDstressHS
% triangle dislocation.
clear
close all;
%% Example: Calculate and plot the first component of displacement vector 
% on a regular grid.

tic
[X,Y,Z] = meshgrid(-2:.01:2,-0:.01:1, 0);
p1 = [0 0 0];
p2 = [0 1 0];
p3 = [0.5 0.5 -0.5];

Ps = [p1; p2; p3];

[ue,un,uv] = TDdispHS(X, Y, Z,p1,p2,p3, 1, 0, 0, .25);
%%
toc
un = reshape(un, size(X));
pcolor(X, Y, un);
colorbar();
view([0,0,1])