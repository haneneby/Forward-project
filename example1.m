%  A=[1.5 2.7 0.3 0.4 0.555 6 .7 .8 .9]
%  B=[1 2 3 4 5 6 7 8 9]
% %msize = numel(A);
% j= A(randperm(numel(A), 1))
% leng= fix(length(A)/2)
% 
% B (A==j)

%  a   = 10*rand(1) 
%  b   = a*rand(1)
%  Y= a + (0.4*rad -2*b).*rand(1)
%  if (Y > 0.4*rad)
%          X= 0.6.*rad-a +(0.2.*rad-2*a).*rand(1) % x of the center of ellipse 
%            'here'
%  else
%         X= a +(rad-2*a).*rand(1)
%  end
  P = phantom(128); 
imshow(P)
title('Original image')
R = radon(P,0:100);
size (R)
I1 = iradon(R,0:100);
I2 = iradon(R,0:100,'linear','none');

figure
subplot(1,2,1)
imshow(I1,[])
title('Filtered Backprojection')
subplot(1,2,2)
imshow(I2,[])
title('Unfiltered Backprojection')

