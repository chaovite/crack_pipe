function [L,U,p,q,B] = imex_ark4_get_lu(A,dt)
% [L,U,p,q,B] = imex_ark4_get_lu(A,dt)
% This function computes the LU decompostion needed for imex_ark4_lu.m
% Input:
% A  : A (matrix) implicit part
% dt  : time step
% Output:
% L,U : LU-decompostion of A-a*h*I (see imex_ark4_get_lu)

  a = 1/4;
  I = speye(size(A));
  B = I - a*dt*A;
  disp('LU decomposition....');
  tic
  [L,U,p,q]=lu(B,'vector');
  toc

