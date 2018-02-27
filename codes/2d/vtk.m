x = linspace(0, 10, 2000);
y = linspace(-5, 5, 3000);

[Y, X] = meshgrid(y, x);

Z = exp(-0.1 * X.^2) .* sin(Y);

surf(X, Y, Z)

nx = length(x);
ny = length(y);

tic
f = fopen('res.vtk', 'wb');
fprintf(f, '# vtk DataFile Version 3.0\n');
fprintf(f, 'Exported from MATLAB\n'); % Comment string
fprintf(f, 'BINARY\n');
fprintf(f, 'DATASET RECTILINEAR_GRID\n');
fprintf(f, 'DIMENSIONS %d %d 1\n', nx, ny);
fprintf(f, 'X_COORDINATES %d float\n', nx);
w = typecast(swapbytes(single(x)), 'uint8');
fwrite(f, w);
fprintf(f, 'Y_COORDINATES %d float\n', ny);
w = typecast(swapbytes(single(y)), 'uint8');
fwrite(f, w);
fprintf(f, 'Z_COORDINATES 1 float\n');
w = typecast(swapbytes(single(0)), 'uint8');
fwrite(f, w);
fprintf(f, 'CELL_DATA %d\n', (nx-1) * (ny-1));
% No cell data
fprintf(f, 'POINT_DATA %d\n', nx * ny);
fprintf(f, 'SCALARS z float\nLOOKUP_TABLE default\n');
w = typecast(swapbytes(single(Z(:))), 'uint8');
fwrite(f, w);
fclose(f);
toc