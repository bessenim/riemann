x = linspace(0, 10, 20);
y = linspace(-5, 5, 30);

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
fprintf(f, 'DATASET STRUCTURED_GRID\n');
fprintf(f, 'DIMENSIONS %d %d 1\n', nx, ny);
fprintf(f, 'POINTS %d float\n', nx * ny);
R = zeros(3, nx, ny);
R(1, :, :) = X;
R(2, :, :) = Y;
R(3, :, :) = Z;
w = typecast(swapbytes(single(R(:))), 'uint8');
fwrite(f, w);
fprintf(f, 'CELL_DATA %d\n', (nx-1) * (ny-1));
% No cell data
fprintf(f, 'POINT_DATA %d\n', nx * ny);
fprintf(f, 'SCALARS z float\nLOOKUP_TABLE default\n');
w = typecast(swapbytes(single(Z(:))), 'uint8');
fwrite(f, w);
fclose(f);
toc