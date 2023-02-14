% première méthode de calcul simple avec un bilan des puissances en
% fonction du profil de vitesses choisi 
function [Couple,Eval,cumEval,v]=voiturecalculs1()
load('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R',"angle","nubatx","nubaty","nucontrx","nucontry","nuconvx","nuconvy","numx","numy","numpx","numpy")


nubat = @(v)interp1(nubatx,nubaty,v,'linear','extrap'); %rendement batterie interpolation
nuconv= @(v)interp1(nuconvx,nuconvy,v,'linear','extrap'); %rendement conversion puissance
nucontr= @(v)interp1(nucontrx,nucontry,v,'linear','extrap'); %rendement controle moteur
num=@(v)interp1(numx,numy,v,'linear','extrap'); %rendement moteur
nump=@(v)interp1(numpx,numpy,v,'linear','extrap'); %rendement transmission
Pkin=@(v)max(0,m*v.*[diff(v);0]/(T/N)); %Puissance cinétique 
Pair=@(v)(1/2)*Cd*A*ro*v.*3; %Puissance des frottements de l'air
Proll=@(v)m*g*fr*cos(angle((T/N)*cumtrapz(v)./l)).*v; % Puissance de la résistance au roulement
Phill=@(v)m*g*sin(angle((T/N)*cumtrapz(v)./l)).*v; %Puissance dépensée pour gravir une pente

fun = @(v)trapz((1./(nubat(v).*nuconv(v).*nucontr(v).*num(v).*nump(v))).*(Pkin(v)+Pair(v)+Proll(v)+Phill(v)))*(T/N); %energie prise à la batterie
cumfun = @(v)cumtrapz((1./(nubat(v).*nuconv(v).*nucontr(v).*num(v).*nump(v))).*(Pkin(v)+Pair(v)+Proll(v)+Phill(v)))*(T/N); %energie prise à la batterie
amax= 0.3*9.81*ones(N,1);%accélération maximale
dmax=-0.3*9.81*ones(N,1);%décélération maximale

%DeltaE = deltaEcinétique + deltaEfrottements air + deltaErésistanceroulement + deltaEpotentiel

%load("couple.mat","v")
v0=zeros(N,1); 
A_LININEQ=[(diag(ones(N-1,1),1)-diag(ones(N,1)));-(diag(ones(N-1,1),1)-diag(ones(N,1)))];
b_LININEQ=[amax*(T/N);-dmax*(T/N)];
A_LINEQ= [1,zeros(1,N-2),1;ones(1,N)*(T/N)];
b_LINEQ= [0;l];

lb=zeros(N,1);
ub=10*ones(N,1);
options= optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',1000*N,'PlotFcn',{'optimplotfval';'optimplotfvalconstr';'optimplotconstrviolation'});

[v,Eval] = fmincon(fun,v0,A_LININEQ,b_LININEQ,A_LINEQ,b_LINEQ,lb,ub,@NONLCON1,options);
cumEval = cumfun(v);
Couple= (m*[diff(v);0]/(T/N)+(1/2)*Cd*A*ro*v.^2+m*g*fr*cos(angle((T/N)*cumtrapz(v)./l))+m*g*sin(angle((T/N)*cumtrapz(v)./l)))/R;
save('couple.mat',"Couple","v")


end