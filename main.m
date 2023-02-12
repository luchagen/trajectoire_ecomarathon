% modules nécessaires : system identification toolbox, optimisation
% toolbox, global optimisation toolbox
%% clear
clear, close all
%% parametres voiture
parametres
%% parametres circuit
Circuitsparams
%% paramfinder data

load("parametres.mat",'l','T','N','m','Cd','A','ro','fr','g','R','f','kv',"PCI","Pseudo_rendement","density_H2_Stp")
file='612_2022_06_02_14_28_raw.csv';
fid = fopen(file);
dataread = textscan(fid, '%s','delimiter',',');
fclose(fid);

datastr=dataread{1};
clear dataread
datastr(1:113)=[]; 
data=reshape(datastr,112,[]);
clear datastr
gps_speed=cellfun(@str2num,data(91,687:4305)'); %vitesse mesurée par le gps
QvH2=cellfun(@str2num,data(43,687:4305)'); % volume net d'hydrogène utilisé
jm3_sincereset=cellfun(@str2num,data(50,687:4305)'); % (temps)
T=(jm3_sincereset(end)-jm3_sincereset(1)) /1000;


equivalentJoule= QvH2 * PCI * Pseudo_rendement *density_H2_Stp; % équivalent en joules du volume mesuré
clear data
eqJoule=equivalentJoule(1);

for i=1:length(equivalentJoule)-1 % on supprime les discontinuités de mesure
    if equivalentJoule(i+1)-eqJoule(end)~=0
        eqJoule=[eqJoule;equivalentJoule(i+1)];
    end
end

eqJoule=interp1(eqJoule,linspace(1,length(eqJoule),length(equivalentJoule))');
speed=[gps_speed(1)];

for i=1:length(gps_speed)-1 % on supprime les discontinuités de mesure
    if gps_speed(i+1)-speed(end)~=0
        speed=[speed;gps_speed(i+1)];
    end
end
speed=interp1(speed,linspace(1,length(speed),length(gps_speed))');
N=length(speed);
joule=[diff(eqJoule);0];

save("parametres_identification.mat",'T','m','A','ro','g','R','speed','joule')
params=[Cd,fr];
Couple= @(params)(m*[diff(speed);0]/(T/N)+(1/2)*params(1)*A*ro*speed.^2+m*g.*params(2))/R;
energy0=@(params)(R/2)*Couple(params).*(speed+diag(ones(N-1,1),1)*speed)*(T/N);
energy1=@(params)max(0,energy0(params));

%% calcul params
figure(2)
plot(joule)
hold on
corr([energy1(params),joule])

newparams = ga(@residuals,2,-eye(2),[0;0]);

newparamscorr=corr([energy1(newparams),energy0(newparams),joule]);
X=[energy1(newparams),ones(N,1)];
[b,bint,r,rint,stats]= regress(joule,X);
plot(energy1(newparams))
legend("consommation au niveau de la pile","estimation du modèle")
title("energie consommée par unité de temps")
xlabel("correlation  : " + newparamscorr(1,2))

figure(3)
plot(eqJoule)
hold on 
plot(cumsum(energy1(newparams)))
legend( "energie totale consommée au niveau de la pile","estimation")
title("energie totale consommée")
xlabel( "consommation réelle  : "+ eqJoule(end)+ "    consommation estimée  : "+sum(energy1(newparams)))

%% optim1

voiturecalculs1
load("parametres.mat","Pkin","Phill","Proll","Pair","fun","nubat","nucontr","nuconv","num","nump")
load('couple.mat','Couple')
figure(11)
plot(Couple)
figure(12)
plot(vitesse(Couple))
figure(13)
plot(cumsum((1./(nubat(vitesse(Couple)).*nuconv(vitesse(Couple)).*nucontr(vitesse(Couple)).*num(vitesse(Couple)).*nump(vitesse(Couple)))).*(Pkin(vitesse(Couple))+Pair(vitesse(Couple))+Proll(vitesse(Couple))+Phill(vitesse(Couple)))*(T/N)))
title('energie prise à la roue')
xlabel("energie totale (J) : " + fun(vitesse(Couple))+ "  (" + fun(vitesse(Couple))/3600000+" kWh)")
%% optim2
load("parametres.mat",'l','T','N','m','Cd','A','ro','fr','g','R','f','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","Pkin","Phill","Proll","Pair")

[couple,Eval,cumEval]=voiturecalculs(Couple,N,0,0,l,'0');
figure(31)
plot(couple)
figure(32)
plot(vitesse(couple))
figure(33)
plot (cumEval)
title('energie prise à la roue - restarts')
xlabel("energie totale (J) : " + Eval + "  (" + Eval/3600000+ " kWh)")

%% optim3
load("parametres.mat",'l','T','N','m','Cd','A','ro','fr','g','R','f','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","Pkin","Phill","Proll","Pair")

[couple,Eval,cumEval]=voiturecalculs(Couple,N,0,0,l,'0.5');
figure(21)
plot(couple)
figure(22)
plot(vitesse(couple))
figure(23)
plot (cumEval)
title('energie prise à la roue - up&down')
xlabel("energie totale (J) : " + Eval + "  (" + Eval/3600000+ " kWh)")
%% optim4
load("parametres.mat",'l','T','N','m','Cd','A','ro','fr','g','R','f','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","Pkin","Phill","Proll","Pair")

[couple,Eval,cumEval]=voiturecalculs(Couple,N,0,0,l,'1');
figure(41)
plot(couple)
figure(42)
plot(vitesse(couple))
figure(43)
plot (cumEval)
title('energie prise à la roue - finalcoast')
xlabel("energie totale (J) : " + Eval + "  (" + Eval/3600000+ " kWh)")

%% optim1 pente

voiturecalculs1pente
load("parametres.mat","Pkin","Phill","Proll","Pair","fun","nubat","nucontr","nuconv","num","nump")
load('couple.mat','Couple')
figure(51)
plot(Couple)
figure(52)
plot(vitesse(Couple))
figure(53)
plot(cumsum((1./(nubat(vitesse(Couple)).*nuconv(vitesse(Couple)).*nucontr(vitesse(Couple)).*num(vitesse(Couple)).*nump(vitesse(Couple)))).*(Pkin(vitesse(Couple))+Pair(vitesse(Couple))+Proll(vitesse(Couple))+Phill(vitesse(Couple)))*(T/N)))
title('energie prise à la roue - pente')
xlabel("energie totale (J) : " + fun(vitesse(Couple))+ "  (" + fun(vitesse(Couple))/3600000+" kWh)")

%% optim3 pente

load("parametres.mat",'l','T','N','m','Cd','A','ro','fr','g','R','f','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","Pkin","Phill","Proll","Pair")

[couple,Eval,cumEval]=voiturecalculspente(Couple,N,0,0,l,'0.5');
figure(21)
plot(couple)
figure(22)
plot(vitesse(couple))
figure(23)
plot (cumEval)
title('energie prise à la roue - pente - up&down')
xlabel("energie totale (J) : " + Eval + "  (" + Eval/3600000+ " kWh)")


%% fun
function r=residuals(params)
    
    load("parametres_identification.mat",'T','m','A','ro','g','R','S','speed','joule','current')
    N=length(speed);
    Couple= @(params)(m*(diag(ones(N-1,1),1)*speed-speed)/(T/N)+(1/2)*params(1)*A*ro*speed.^2+(m*g)*params(2))/R;
    energy0=@(params)(R/2)*Couple(params).*(speed+diag(ones(N-1,1),1)*speed)*(T/N);
    energy1=@(params)max(0,energy0(params));
    X=[energy1([params(1);params(2)]),ones(N,1)];
    [~,~,r,~,~]= regress(joule,X);
    r=sum(abs(r));
end