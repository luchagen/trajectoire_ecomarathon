function vmaxvirage =vmaxvirages(s)
%paramètres du circuit
CIRCUIT= readmatrix("CIRCUIT.txt");

vmaxvirage=100*ones(length(s),1);
for i=1:length(s)
    atteint=0;
    for j=1:length(CIRCUIT)
        if s(i)>CIRCUIT(1,j) && atteint==0
            atteint=1;
            vmaxvirage(i)=sqrt(CIRCUIT(2,j)*f*g);
        end
    end
end