%circuit (courbes de béziers)
p= [61.304947,122.15887 ; 17.472746,82.068447; -20.045214,52.652096; -45.70309,49.444862 ;-25.657873,-3.207234 ;-64.67922,37.685003 ;-80.71539,41.694046 ;-16.03618,4.009044 ;-38.48682,-0.53454; -68.421,-39.02135; -29.93419,-38.486812; -20.98492,-62.322915; 1.31567,-82.355652; 22.3006,-20.032739; 60.11822,-48.451229; 82.28155,-34.675974; 25.223395,15.677188; 12.39111,35.347556; 2.61063,38.458319; -14.48255,4.606298; -59.21823,7.107624; -57.03355,32.3551502; 1.37596,15.9015428; 13.49277,23.1900278; 58.60224,12.7686948; 10.526014,-2.431753; 48.105094,-6.38201; 74.824536,5.930197; 22.4514131,10.345517; 115.51579,82.381087; 124.602952,89.196457; 9.087162,6.81537; 106.184292,77.50816; 120.750472,88.19894; 14.56619,10.69078; 104.91029,77.11521; 113.62198,83.45016; 73.84802,53.70065; 89.16048,59.61721; 85.05263,74.82415; -4.58283,16.96527; -29.384,11.28563; -59.77763,8.78651; -10.21464,-0.83991; -26.99451,-7.8006; -48.14833,-48.03532; -21.8562,-41.57065; -22.60999,-39.32498; -35.94775,-39.02529; -14.56497,0.32726; -17.80537,24.41163; -28.36808,49.71049; -9.59228,22.97457; 1.62293,39.71041; 13.84823,48.44064; 7.71954,5.51262; 111.64822,27.4349; 117.32345,28.09575; 6.74659,0.78559; 61.40242,2.04256; 79.47978,8.77009; 19.78805,7.36416; 13.61177,25.78989; 5.77593,40.36279; -8.30503,15.4455; -15.91874,23.2323; -62.29001,10.00247; -46.37127,-13.22985; -146.67954,-47.56055; -167.25929,-56.5141; -20.57976,-8.95353; -26.50833,-32.1411; -21.83111,-48.44454; 4.67722,-16.30347; 13.72693,-42.68424; 6.04436,-68.50207; 205.05075,275.1351; 69.89098,131.63024; 61.20472,123.74579];
%h=sin(linspace(0,10*pi,100))'; %profil altimétrique
h=zeros(100,1); %profil altimétrique
l=3636; %longueur du circuit
save('parametres.mat',"l","h","p")
Circuitsparams

%simulation
T=576; %temps de parcours voulu 
N=100; % nombre de points de calcul

%voiture
kv=125*pi*2/60; %constante de tension du moteur
m=94.37+70; %masse de la voiture (+conducteur)
Cd=0.1290; %frontal drag coefficient
A = 1; %cross-sectional aera
ro=1.292; %densité de l'air
fr=0.0126; %coefficient de résistance au roulement
g= 9.81;
R=0.2625; %rayon d'une roue
f= 0.8; %coefficient de frottement

%PAC
PCI=119930; %J/g batterie
Pseudo_rendement = 0.5; %batterie
density_H2_Stp=2/(0.0821*273.15); %g/l batterie

%rendement batterie
nubatx=linspace(0,10,10); % abscisses de mesure du rendement (vitesses)  
nubaty=ones(1,10);%ordonnée de mesure du rendement (rendement)

%rendement conversion puissance
nuconvx=linspace(0,10,10); % abscisses de mesure du rendement (vitesses)  
nuconvy=ones(1,10);%ordonnée de mesure du rendement (rendement)

%rendement controle moteur
nucontrx=linspace(0,10,10); % abscisses de mesure du rendement (vitesses)  
nucontry=ones(1,10);%ordonnée de mesure du rendement (rendement)

%rendement moteur
numx=linspace(0,10,10); % abscisses de mesure du rendement (vitesses)  
numy=ones(1,10);%ordonnée de mesure du rendement (rendement)

%rendement transmission 
numpx=linspace(0,10,10); % abscisses de mesure du rendement (vitesses)  
numpy=ones(1,10);%ordonnée de mesure du rendement (rendement)


%pente
load("CIRCUIT.mat","theta")
angle=@(parcours)interp1(theta',parcours*length(theta),'nearest','extrap');



save('parametres.mat','l','T','N','m','Cd','A','ro','fr','g','R','f',"kv","PCI","Pseudo_rendement","density_H2_Stp","nubatx","nubaty","nucontrx","nucontry","nuconvx","nuconvy","numx","numy","numpx","numpy","angle",'p','h')