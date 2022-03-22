function l2_Error_4
h1 = [  0.01 0.05 0.1 0.2  1/3 0.5  ];  %Step Sizes
L2Error=[0  0  0  0  0  0 ]; 
for l=1:6
u = @(X) X.*(1-X).^2; %Exact Solution    

f1=0;
mesh=0:h1(l):1; % Refinements 
alpha =0;      %Boundary Conditions
beta =0;
 N1(l) = (1-0)/h1(l) 
[uhh] = compute_GLOBAL_MATRIX(h1(l),N1(l));
uhh
 Diff = @(X)(u(X) - interp1(mesh,uhh,X));
% Quadrature Points and Weights (2-point Gauss-Legendre Quadrature)
 tQ = [ -1 1 ]'/sqrt(3); %Points
 wQ = [ 1 1 ];           %Weights
for i=1:N1(l)
     x = @(X) mesh(i) + h1(l)*(X + 1)/2; % Affine Transformation 
     f1 = (wQ*(Diff(x(tQ)).^2));  %Integration
     f1 = f1*(h1(l)/2);
 L2Error(l) = L2Error(l) + f1;
 end
L2Error(l) = sqrt(L2Error(l))
end
h1
N1
L2Error
figure(1)
plot(N1,L2Error,'r');
figure(2)
plot(h1,L2Error,'b');