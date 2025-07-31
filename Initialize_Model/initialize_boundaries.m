function c2_boundary =initialize_boundaries()

%Boundaries are defined as a set of contours e.g. [x1 NaN x2 ; y1 NaN y2]

%%Kobuk Whole Lake
Lx = 3.5e4; Ly = 4e4; 
x=[-1 -1 1 1]*Lx; 
y=[-1 1 1 -1]*Ly;

c2_boundary = [x; y];

end