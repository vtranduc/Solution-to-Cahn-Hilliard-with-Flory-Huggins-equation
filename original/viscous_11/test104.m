function test104

clear
clc





% Generate sample data in the rand 1.3 to 45.3
z = 3*(peaks(100)+7);
% Make colors be red above 27, and cyan below 27.
redChannel = z > 23;

size(redChannel)

greenChannel = z < 27;
blueChannel = greenChannel;
% Make the RGB image.
colors = double(cat(3, redChannel, greenChannel, blueChannel));
% Plot the surface with those colors.

size(colors)

max(max(z))

surf(z, colors*1);
sum(sum(sum(colors)))
nnz(colors)
colors
% Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

c = hadamard(2^5)


end