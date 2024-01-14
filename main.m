% ------------ MAIN function ------------

% add filename of image to apply code to
filename = 'test6.jpg';

% read image file of type jpeg -> low-light input L
L = (im2double(imread(filename)));

% set values needed for map calculations
vals_lambda = 0.15;
vals_sigma = 2;
vals_gamma = 0.8;

% apply LIME function
[I, Ti, Tf] = LIME(L, vals_lambda, vals_sigma, vals_gamma);

% plot initial image in figure window
figure(1);
imshow(L);
title('INPUT IMAGE');

% plot resulting image in figure window
figure(2);
imshow(I);
title('RESULT IMAGE');

% print result image in directory in jpeg format
imwrite(I, 'test_out.jpg');
