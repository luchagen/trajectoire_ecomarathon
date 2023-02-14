load("parametres.mat",'l','p','h')

%h=zeros(100,1);
n=100;
m= fix(length(p)/3);
c=5;
C=length(p)-2;
X=[];
Y=[];
for i=1:m
    if 4+(3*(i-1)) < c
        Px=[p(1+(3*(i-1)),1) , p(2+(3*(i-1)),1), p(3+(3*(i-1)),1),p(4+(3*(i-1)),1)];
        cx=p(4+(3*(i-1)),1);
        Py=-[p(1+(3*(i-1)),2) , p(2+(3*(i-1)),2), p(3+(3*(i-1)),2),p(4+(3*(i-1)),2)];
        cy=p(4+(3*(i-1)),2);
    elseif 4+(3*(i-1)) < C
        Px=[cx ,cx + p(2+(3*(i-1)),1), cx+p(3+(3*(i-1)),1),cx + p(4+(3*(i-1)),1)];
        cx=Px(4);
        Py=-[cy , cy+p(2+(3*(i-1)),2), cy+p(3+(3*(i-1)),2),cy+p(4+(3*(i-1)),2)];
        cy=-Py(4);
    else
        Px=[cx , p(2+(3*(i-1)),1), p(3+(3*(i-1)),1),p(4+(3*(i-1)),1)];
        cx=p(4+(3*(i-1)),1);
        Py=-[cy , p(2+(3*(i-1)),2), p(3+(3*(i-1)),2),p(4+(3*(i-1)),2)];
        cy=-p(4+(3*(i-1)),2);
    end
    [Qx,Qy]= Funct_Bezier(Px,Py,n);
    X=[X,Qx];
    Y=[Y,Qy];
end
plot(X,Y)

Rcircuit=Inf*ones(length(X),1);

for i=2:length(X)-2
    m1=-(X(i+1)-X(i-1))/(Y(i+1)-Y(i-1));
    b1=Y(i)/(m1*X(i));
    m2=-(X(i+2)-X(i))/(Y(i+2)-Y(i));
    b2=Y(i+1)/(m2*X(i+1));
    xr=(b2-b1)/(m1-m2);
    yr=m1*xr+b1;
   
    Rcircuit(i)=sqrt((xr-X(i+1))^2+(yr-Y(i+1))^2);
       
    plot(xr,yr,'.')
end


figure(92)
hcircuit=interp1(h,linspace(1,length(h),length(Rcircuit))');
theta = atan([diff(hcircuit);0]./(l/(n*m)));
theta(isnan(theta))=0;
plot(Rcircuit)
figure(93)
plot(hcircuit)
figure(94) 
plot(theta)
save("CIRCUIT.mat","Rcircuit","hcircuit","theta")

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------