
kv=125*pi*2/60; %constante de tension du moteur
l=3636; %longueur du circuit
T=576; %temps de parcours voulu 
N=100; % nombre de points de calcul
m=94.37+70; %masse de la voiture (+conducteur)
Cd=0.1290; %frontal drag coefficient
A = 1; %cross-sectional aera
ro=1.292; %densité de l'air
fr=0.0126; %coefficient de résistance au roulement
g= 9.81;
R=0.2625; %rayon d'une roue
f= 0.8; %coefficient de frottement
h=@(s)zeros(length(s),1); %profil en dénivelé du circuit
nubat = @(v)ones(N,1); %rendement batterie
nuconv= @(v)ones(N,1); %rendement conversion puissance
nucontr= @(v)ones(N,1); %rendement controle moteur
num=@(v)ones(N,1); %rendement moteur
nump=@(v)ones(N,1); %rendement transmission
Pkin=@(v)max(0,m*v.*[diff(v);0]/(T/N)); %Puissance cinétique 
Pair=@(v)(1/2)*Cd*A*ro*v.*3; %Puissance des frottements de l'air
Proll=@(v)m*g*fr*v; % Puissance de la résistance au roulement
Phill=@(v)m*g*[diff(h(cumtrapz(v)*(T/N)));0]; %energie dépensée pour gravir une pente
PCI=119930; %J/g batterie
Pseudo_rendement = 0.5; %batterie
density_H2_Stp=2/(0.0821*273.15); %g/l batterie
save('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","PCI","Pseudo_rendement","density_H2_Stp","Pkin","Phill","Proll","Pair","nubat","nucontr","nuconv","num","nump")