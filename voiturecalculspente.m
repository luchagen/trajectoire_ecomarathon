function [couple,Eval,cumEval]=voiturecalculspente(couple0,N,v0,v1,laparcourir,integration)
load("parametres.mat",'l','T','m','Cd','A','ro','fr','g','R','f')
load("CIRCUIT.mat","theta")
angle=@(s)interp1(theta',1+s*length(theta),'linear',0);
  
% interp1(Couple,linspace(1,N,l))';
% N=l;
% save('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R','f')


h=@(s)interp1(hcircuit',1+s*length(hcircuit),'linear',0); %profil en dénivelé du circuit

nubat = @(v)ones(N,1); %rendement batterie
nuconv= @(v)ones(N,1); %rendement conversion puissance
nucontr= @(v)ones(N,1); %rendement controle moteur
num=@(v)0.94*ones(N,1); %rendement moteur
nump=@(v)ones(N,1); %rendement transmission

Pkin=@(v)max(0,m*v.*[diff(v);0]/(T/N)); %Puissance cinétique 
Pair=@(v)(1/2)*Cd*A*ro*v.*3; %Puissance des frottements de l'air
Proll=@(v)m*g*fr*v; % Puissance de la résistance au roulement
Phill=@(v)m*g*[diff(h(cumtrapz(v)*(T/N)));0]; %energie dépensée pour gravir une pente

if integration == '0'
    fun = @(couple)sum((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R)*(couple).*(Ecarloc(couple))))*(T/N); %energie prise à la roue, démarrages
    cumfun = @(couple)cumsum((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R)*(couple).*(Ecarloc(couple))))*(T/N); %energie prise à la roue, démarrages
end
if integration == '0.5'   
    fun =@(couple)sum((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/4)*(couple+diag(ones(N-1,1),1)*couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N); %energie prise à la roue, coast zigzag
    cumfun = @(couple) cumsum((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/4)*(couple+diag(ones(N-1,1),1)*couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N); %energie prise à la roue, coast zigzag
end
if integration == '1'   
    fun = @(couple)sum((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/2)*(couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N) ;%energie prise à la roue, coast final
    cumfun = @(couple) cumsum((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/2)*(couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N) ;%energie prise à la roue, coast final
end
A_LININEQ=[zeros(1,N)];
b_LININEQ=[1];
A_LINEQ=[zeros(1,N)];
b_LINEQ=[0];

lb=-940*ones(N,1);
ub=854*ones(N,1);
options= optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-3,'MaxFunctionEvaluations',1000*N,'MaxIterations',1000*N,'PlotFcn',{'optimplotfval';'optimplotfvalconstr';'optimplotconstrviolation'});


[couple,Eval] = fmincon(fun,couple0,A_LININEQ,b_LININEQ,A_LINEQ,b_LINEQ,lb,ub,@NONLCONloc,options);
cumEval= cumfun(couple);

function [c,ceq]= NONLCONloc(couple)
    load("parametres.mat",'l','m','Cd','A','ro','fr','g','R','f')
    parcourscircuit=tril(ones(length(Ecarloc(couple))))*Ecarloc(couple)*(T/length(Ecarloc(couple)));
    
    load("CIRCUIT.mat","Rcircuit")
    Rc=@(s)interp1(Rcircuit,s);

    vmaxvirage=sqrt(Rc(parcourscircuit)*f*g);

    vmax=@(s)min(100*ones(length(s),1),vmaxvirage); %vitesses maximales autorisées à chaque point du circuit
    

    c=[Ecarloc(couple)-vmax(parcourscircuit),-Ecarloc(couple),(R/4)*(couple+diag(ones(N-1,1),1)*couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))-1000*(T/N)];%vitesse max circuit
    ceq=[ones(1,N)*(T/N)*Ecarloc(couple)-laparcourir,[zeros(1,N-1),1]*Ecarloc(couple)-v1];%condition de parcours du circuit + vitesse nulle à l'arrivée
end
function v= Ecarloc(couple)%renvoie la vitesse

  
    
    %load("parametres.mat",'l','m','Cd','A','ro','fr','g','R','f')
    v=zeros(length(couple),1);
    vcurrent=v0;
    
    for j=1:length(couple)
        v(j)=vcurrent;
        vcurrent=vcurrent+ (couple(j)*R-sign(vcurrent)*(1/2)*Cd*A*ro*vcurrent^2-sign(vcurrent)*m*g*fr*cos(anglesol((T/N)*trapz(vcurrent),l))-m*g*sin(anglesol((T/N)*trapz(vcurrent),l)))*(T/N)/m; %bilan des forces
        
    end
    
end
function thetac = anglesol(parcourscircuit,l)
        
        
        progress = parcourscircuit./l;
        thetac = angle(progress);
        thetac(isnan(thetac))=0;
end
end
