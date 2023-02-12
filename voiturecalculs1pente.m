
kv=125*pi*2/60; %constante de tension du moteur
l=3636; %longueur du circuit
T=576; %temps de parcours voulu 
N=100;
m=94.37+70; %masse de la voiture (+conducteur)
Cd=0.1290; %frontal drag coefficient
A = 1; %cross-sectional aera
ro=1.292; %densité de l'air
fr=0.0126; %coefficient de résistance au roulement
g= 9.81;
R=0.2625; %rayon d'une roue
f= 0.8; %coefficient de frottement


nubat = @(v)ones(N,1); %rendement batterie
nuconv= @(v)ones(N,1); %rendement conversion puissance
nucontr= @(v)ones(N,1); %rendement controle moteur
num=@(v)ones(N,1); %rendement moteur
nump=@(v)ones(N,1); %rendement transmission


Pkin=@(v)max(0,m*v.*[diff(v);0]/(T/N)); %Puissance cinétique 
Pair=@(v)(1/2)*Cd*A*ro*v.*3; %Puissance des frottements de l'air
Proll=@(v)m*g*fr*cos(anglesol((T/N)*cumtrapz(v),l)).*v; % Puissance de la résistance au roulement
Phill=@(v)m*g*sin(anglesol((T/N)*cumtrapz(v),l)).*v; %Puissance dépensée pour gravir une pente

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

Couple= (m*[diff(v);0]/(T/N)+(1/2)*Cd*A*ro*v.^2+m*g*fr*cos(anglesol((T/N)*cumtrapz(v),l))+m*g*sin(anglesol((T/N)*cumtrapz(v),l)))/R;
save('couple.mat',"Couple","v")

function thetac = anglesol(parcourscircuit,l)
        load("CIRCUIT.mat","theta")
        angle=@(s)interp1(theta',1+s*length(theta),'linear',0);
        
        progress = parcourscircuit./l;
        thetac = angle(progress);
        thetac(isnan(thetac))=0;
end