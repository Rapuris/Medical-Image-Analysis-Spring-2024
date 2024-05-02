%Initialize snake
c0 = 2;
initialLSF = c0*ones(size(Img));
initialLSF(10:55, 10:75) = -c0;
phi = initialLSF;

%Edge indicator function.
G = fspecial('gaussian',15,sigma);
Img_smooth = conv2(Img, G, 'same');