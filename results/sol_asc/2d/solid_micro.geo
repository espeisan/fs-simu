// Gmsh project created on Tue May 19 13:43:01 2015
lc = 0.1;  //0.01
lb = lc;
xc = 0.5;
yc = 0.5;
R = 0.2;

Point(1) = {0, 0, 0, lb};
Point(2) = {1, 0, 0, lb};
Point(3) = {1, 2, 0, lb};
Point(4) = {0, 2, 0, lb};
Point(5) = {xc, yc, 0, lc};
Point(6) = {xc + R, yc, 0, lc};
Point(7) = {xc, yc + R, 0, lc};
Point(8) = {xc - R, yc, 0, lc};
Point(9) = {xc, yc - R, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};
Line Loop(9) = {4, 1, 2, 3};
Line Loop(10) = {6, 7, 8, 5};
Plane Surface(11) = {9, 10};
//Plane Surface(12) = {10};
Physical Line(13) = {4, 1, 2, 3};  //borde externo
Physical Line(14) = {6, 7, 8, 5};  //interface fs
Physical Surface(15) = {11};
//Physical Surface(16) = {12};

//Transfinite Line {5, 6, 7, 8} = 20;// Using Bump 0.2;
//Transfinite Line {1, 2, 3, 4} = 20;
//Transfinite Surface {12};
//Transfinite Surface {11};
//Recombine Surface {12};
//Recombine Surface {11};
