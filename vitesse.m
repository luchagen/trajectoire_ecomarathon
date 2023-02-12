function v= vitesse(couple)%renvoie la vitesse

  
    
    load("parametres.mat",'l','T','m','Cd','A','ro','fr','g','R','f')
    v=zeros(length(couple),1);
    v0=0;
% 
%     Ecar=zeros(length(couple),1);
%     Ecartot=1;
%     deltaEcar = zeros(length(couple),1);
    
    for i=1:length(couple)
        v(i)=v0;
        v0=v0+ (couple(i)*R-sign(v0)*(1/2)*Cd*A*ro*v0^2-sign(v0)*m*g*fr*cos(anglesol((T/length(v))*cumtrapz(v0),l))-m*g*sin(anglesol((T/length(v))*cumtrapz(v0),l)))*(T/length(v))/m; %bilan des forces
        
    end
    function theta = anglesol(parcourscircuit,l)
        load("CIRCUIT.mat","hcircuit")
        h=@(s)interp1(hcircuit',1+s*length(hcircuit),'linear',0); %profil en dénivelé du circuit
        progress = parcourscircuit./l;
        variation = [diff(h(progress));0];
        theta = atan(variation./[diff(parcourscircuit);0]);
        theta(isnan(theta))=0;
    end
end

