// Gmsh project created on Thu Jan 23 15:53:43 2020
alpha = 1.0; // aspect ratio
size_cc = 0.1;
size_bc = 0.03;
//+
Point(1) = {-0.5, 0, 0, 0};
//+
Point(2) = {0.5, 0, 0, 0};
//+
Point(3) = {0, 0.5*alpha, 0, 0};
//+
Point(4) = {0, -0.5*alpha, 0, 0};
//+
Point(5) = {0, 0, 0, 0};
//+
Ellipse(1) = {2, 5, 2, 3};
//+
Ellipse(2) = {3, 5, 1, 1};
//+
Ellipse(3) = {1, 5, 1, 4};
//+
Ellipse(4) = {4, 5, 2, 2};
//+
Line Loop(1) = {2, 3, 4, 1};
//+
Characteristic Length {3, 1, 4, 2} = size_bc;
//+
Characteristic Length {5} = size_cc;
//+
Plane Surface(1) = {1};
Point{5} In Surface{1};
