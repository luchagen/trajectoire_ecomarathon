function [couple,Eval,cumEval,margev]=voiturecalculs(couple0,N,v0,v1,laparcourir,integration)
load('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R','f',"Pkin","Phill","Proll","Pair","nubatx","nubaty","nucontrx","nucontry","nuconvx","nuconvy","numx","numy","numpx","numpy")


  
% interp1(Couple,linspace(1,N,l))';
% N=l;
% save('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R','f')

h=@(s)zeros(length(s),1); %profil en dénivelé du circuit


nubat = @(v)interp1(nubatx,nubaty,v,'linear','extrap'); %rendement batterie interpolation
nuconv= @(v)interp1(nuconvx,nuconvy,v,'linear','extrap'); %rendement conversion puissance
nucontr= @(v)interp1(nucontrx,nucontry,v,'linear','extrap'); %rendement controle moteur
num=@(v)interp1(numx,numy,v,'linear','extrap'); %rendement moteur
nump=@(v)interp1(numpx,numpy,v,'linear','extrap'); %rendement transmission

fun = @(v)trapz((1./(nubat(v).*nuconv(v).*nucontr(v).*num(v).*nump(v))).*(Pkin(v)+Pair(v)+Proll(v)+Phill(v)))*(T/N); %energie prise à la batterie


if integration == '0'
    fun = @(couple)trapz((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R)*(couple).*(Ecarloc(couple))))*(T/N); %energie prise à la roue, démarrages
    cumfun = @(couple)cumtrapz((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R)*(couple).*(Ecarloc(couple))))*(T/N); %energie prise à la roue, démarrages
end
if integration == '0.5'   
    fun =@(couple)trapz((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/4)*(couple+diag(ones(N-1,1),1)*couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N); %energie prise à la roue, coast zigzag
    cumfun = @(couple) cumtrapz((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/4)*(couple+diag(ones(N-1,1),1)*couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N); %energie prise à la roue, coast zigzag
end
if integration == '1'   
    fun = @(couple)trapz((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/2)*(couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N) ;%energie prise à la roue, coast final
    cumfun = @(couple) cumtrapz((1./(nubat(Ecarloc(couple)).*nuconv(Ecarloc(couple)).*nucontr(Ecarloc(couple)).*num(Ecarloc(couple)).*nump(Ecarloc(couple)))).*max(0,(R/2)*(couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))))*(T/N) ;%energie prise à la roue, coast final
end
A_LININEQ=[zeros(1,N)];
b_LININEQ=[1];
A_LINEQ=[zeros(1,N)];
b_LINEQ=[0];

lb=-940*ones(N,1);
ub=854*ones(N,1);
options= optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-1,'MaxFunctionEvaluations',1000*N,'MaxIterations',1000*N,'PlotFcn',{'optimplotfval';'optimplotfvalconstr';'optimplotconstrviolation'});


[couple,Eval] = fmincon(fun,couple0,A_LININEQ,b_LININEQ,A_LINEQ,b_LINEQ,lb,ub,@NONLCONloc,options);
cumEval= cumfun(couple);
margev= NONLCONloc(couple);

function [c,ceq]= NONLCONloc(couple)
    load("parametres.mat",'l','m','Cd','A','ro','fr','g','R','f')
    parcourscircuit=tril(ones(length(Ecarloc(couple))))*Ecarloc(couple)*(T/length(Ecarloc(couple)));
    
    load("CIRCUIT.mat","Rcircuit")
    Rc=@(s)interp1(Rcircuit,s);

    vmaxvirage=sqrt(Rc(parcourscircuit)*f*g);

    vmax=@(s)min(100*ones(length(s),1),vmaxvirage); %vitesses maximales autorisées à chaque point du circuit
    

    c=[Ecarloc(couple)-vmax(parcourscircuit),-Ecarloc(couple),(R/4)*(couple+diag(ones(N-1,1),1)*couple).*(Ecarloc(couple)+diag(ones(N-1,1),1)*Ecarloc(couple))-1000];%vitesse max circuit
    ceq=[ones(1,N)*(T/N)*Ecarloc(couple)-laparcourir,([zeros(1,N-1),1]*Ecarloc(couple)-v1)];%condition de parcours du circuit + vitesse nulle à l'arrivée
end
function v= Ecarloc(couple)%renvoie la vitesse

    load('parametres.mat','angle')
    v=zeros(length(couple),1);
    vcurrent=v0;
    
    for j=1:length(couple)
        v(j)=vcurrent;
        theta =angle((T/N)*trapz(v)./l);
        vcurrent=vcurrent+ (couple(j)*R-sign(vcurrent)*(1/2)*Cd*A*ro*vcurrent^2-sign(vcurrent)*m*g*fr*cos(theta)+m*g*sin(theta))*(T/N)/m; %bilan des forces
        
    end
end

end
