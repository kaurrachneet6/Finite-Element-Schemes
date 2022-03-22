%Program with dirichlet boundary conditions
% u'(t)-u''=f
%u(0,x)=cos x;
%u(t,0)=u(t,1)=0;
format  long 
clear all 
k=0;
m=0;
%N=input('N= '); 
% N denotes the number of nodal points.
   for p=1:1
   %N=20;
   N=input('Enter the nodes between 0 and 1');
   h1(p)=N; 
x=linspace(0,1,N+2);
y=x(2:N+1);
h=1/(N+1);
p=input('enter the time nodes between 0 and 1(it should be greater or equal to no. of nodes)');
tdiff=1/(p+1);
t1=0:tdiff:1;
t=t1(2:p+1);
ustart=zeros(N+2,p+2);
for i=1:N+2
ustart(i,1)=sin(pi*(x(i)));
end
for i=1:N+2
disp('The nodes are');
disp(x(i));
end
disp('No of cells are');
disp(N+1);
%y=x(2:N+1);
h=1/(N+1);
for i=1:N+1 
   elem(i,:)=[i,i+1]; 
end 
A=sparse(N+2,N+2);
K=sparse(N+2,N+2); 
F=sparse(N+2,1); 
for i=1:N+1
    elem(1,:);
     A(elem(i,:),elem(i,:))=[(2*h)/3  h/6;h/6 (2*h)/3];
     K(elem(i,:),elem(i,:))=[2/h -1/h;-1/h 2/h];
end 

% for i=1:N   
% mp=(x(i)+x(i+1))/2; 
% a=x(i);
% b=x(i+1);
% F(i,1)=(h/2)*f1_diff(mp);
% mp=(x(i+1)+x(i+2))/2;
% % For mid point rule
% %The main matrix of F
% F(i,1)=F(i,1)+(h/2)*f1_diff(mp);
% end

%uh=K\F;
disp('the A matrix is');
A
disp('The K matrix is ');
K
disp('The F matrix is');

F
 for k=2:p+2
     tdiff
     L= tdiff*K+A
L1=tdiff*F+A*ustart(1:N+2,k-1)
    %ustart(:,k)=inv(A)*ustart(:,k-1)+ tdiff*( inv(A)*(F - K*ustart(:,k-1)));
    L(1,1)=1;
    for o=2:N+2
        L(1,o)=0;
    end
     L(N+2,N+2)=1;
    for o=1:N+1
        L(N+2,o)=0;
    end
    L1(1)=0;
    L1(N+2)=0;
     ustart(1:N+2,k)= L\L1;  %(ustart(2:N+1,k-1))+ tdiff*(inv(A)*(F- K* (ustart(2:N+1,k-1))));
 end
% uh(2:N+1,:)=K\F;
for k=1:p+2
 ustart(1,k)=0;
 ustart(N+2,k)=0;
end
k=1;
ustart = ustart
% Exact solution at nodes
 for m=t1
    
for i=1:N+2 
   ue(i,k)=sin(pi*x(i))*exp((-(pi^2))*(m));   
end 
k=k+1;
 end
   end 
   ue2 =transpose(ue);
surf(x,t1,ue2);
xlabel('x');
ylabel('t');
zlabel('u_exact(x,1)');
figure;

disp('the exact solution is');
ue=ue
disp('The solution matrix  of u is');
disp('the value of u is shown at time t=0,1/6,1/3,1/2,2/3,5/6,1 ');
soln=ustart
%plot(x,ue,'-r', x,ustart,'-.b')
%plot3(x,t,ustart);
v=transpose(ustart);
%surf(y,t1,ustart)
surf(x,t1,v);
xlabel('x');
ylabel('t');
zlabel('u_calculated(x,1)');


