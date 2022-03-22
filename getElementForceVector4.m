function [fe] = getElementForceVector4(m1)
% getElementForceVector - returns element force vector for an individual element
% [fe] = getElementForceVector(func,xCoords)
% func:-function handle to RHS f(x)
% xCoords:- coordinates of 1D element [ xMin xMax ]
% fe: - element force vector
 h = m1(2)-m1(1); %element width 
 N1 = @(t) (1-t)/2; % linear hat function
 N2 = @(t) (1+t)/2;
 x = @(t) m1(1) + h*(t + 1)/2; % affine transformation
 % quadrature points and weights (2-point Gauss-legendre)
 tQ = [ -1 1 ]'/sqrt(3);
 wQ = [ 1 1 ];
 
 %Force vector matrix
 fe(1) = wQ*(((x(tQ).^3)-(2*(x(tQ).^2))-(5*(x(tQ)))+4).*N1(tQ));
 fe(2) = wQ*(((x(tQ).^3)-(2*(x(tQ).^2))-(5*(x(tQ)))+4).*N2(tQ));
 fe = fe'*h/2;