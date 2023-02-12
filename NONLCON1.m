
function [c,ceq]= NONLCON1(v)
    
    load("parametres.mat",'l','T','N','m','Cd','A','ro','fr','g','R','f')
    parcourscircuit=tril(ones(length(v)))*v*(T/length(v))./l;
    load("CIRCUIT.mat","Rcircuit")
    Rc=@(s)interp1(Rcircuit,1+s*(length(Rcircuit)));

    vmaxvirage=sqrt(Rc(parcourscircuit)*f*g);

    vmax=@(s)min(100*ones(length(s),1),vmaxvirage); %vitesses maximales autorisées à chaque point du circuit
    

    c=v-vmax(parcourscircuit);%limite de vitesse
    ceq=0;
end
