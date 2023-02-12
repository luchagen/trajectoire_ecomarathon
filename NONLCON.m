function [c,ceq]= NONLCON(couple,N)
    load("parametres.mat",'l','T','m','Cd','A','ro','fr','g','R','f')
    parcourscircuit=tril(ones(length(vitesse(couple))))*vitesse(couple)*(T/length(vitesse(couple)))./l;
    

    load("CIRCUIT.mat","Rcircuit")
    Rc=@(s)interp1(Rcircuit,1+s*(length(Rcircuit)));

    vmaxvirage=sqrt(Rc(parcourscircuit)*f*g);

    vmax=@(s)min(100*ones(length(s),1),vmaxvirage); %vitesses maximales autorisées à chaque point du circuit
    

    c=vitesse(couple)-vmax(parcourscircuit);%vitesse max circuit
    ceq=[ones(1,N)*(T/N)*vitesse(couple)-l,[zeros(1,N-1),1]*vitesse(couple)];%condition de parcours du circuit + vitesse nulle à l'arrivée
end