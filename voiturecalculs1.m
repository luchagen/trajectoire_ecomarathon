

load('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","Pkin","Phill","Proll","Pair","fun","nubat","nucontr","nuconv","num","nump")

fun = @(v)trapz((1./(nubat(v).*nuconv(v).*nucontr(v).*num(v).*nump(v))).*(Pkin(v)+Pair(v)+Proll(v)+Phill(v)))*(T/N); %energie prise à la batterie

save('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","Pkin","Phill","Proll","Pair","fun","nubat","nucontr","nuconv","num","nump")

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

Couple= (m*[diff(v);0]/(T/N)+(1/2)*Cd*A*ro*v.^2+m*g*fr)/R;
save('couple.mat',"Couple","v")