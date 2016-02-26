lw = 0.5;

//Box
Point(1) = {0,0,0,lw};
Point(2) = {5,0,0,lw};
Point(3) = {5,2,0,lw};
Point(4) = {0,2,0,lw};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};

//Physical properties
freeslip = 1;
inflow = 2;
outflow = 3;
Physical Line(freeslip) = {1,3};
Physical Line(inflow) = {4};
Physical Line(outflow) = {2};
Physical Point(freeslip) = {2,3};
Physical Surface(7) = {6};
