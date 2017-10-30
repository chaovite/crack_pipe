function [L,U,p,q,B] = imex_rk4i1_get_lu(A,dt)
% [L,U,p,q,B] = imex_rk4i1_get_lu(A,dt)
% This function computes the LU decompostion needed for imex_rk4i1_lu.m
% Input:
% A  : A (matrix) implicit part
% dt  : time step
% Output:
% L,U : LU-decompostion of A-h*I (see imex_rk4i1_get_lu)

  I = eye(size(A));
  B = I - dt*A;
  [L,U,p,q]=lu(sparse(B),'vector');

