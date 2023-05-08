
figure;
i_plotit(0,1)
i_plotit(0.25,2)
i_plotit(0.5,3)
i_plotit(0.75,4)
i_plotit(1,5)

function i_plotit(a,id)    
    f2=0:0.15:pi+0.15;
    f1=-1.6:0.2:1.6;
    [X,Y]=meshgrid(f1,f2);
    Z=zeros(length(f2),length(f1));
    for k=1:length(f2)
        for l=1:length(f1)
            gates = [ryGate(1,a*pi); ryGate(2,f2(k)); ...
                     cryGate(1,2,f1(l))];
            c = quantumCircuit(gates);
            s = simulate(c);
            Z(k,l)=sqrt(probability(s,2,"1"));
        end
    end
    subplot(2,3,id+1)
    surface(X,Y,Z);
    fs=25;
    t=sprintf('Angle of 1st qubit, \\phi_{1} = %g\\pi',a);
    title(t,'FontSize',fs,'FontWeight','normal',...
        'Interpreter','tex');    
    xlabel('Rotation angle of CRY gate \theta','FontSize',fs);
    ylabel('Angle of 2nd qubit, \phi_{2}','FontSize',fs);
    xlim([min(X(:)) max(X(:))])
    ylim([min(Y(:)) max(Y(:))])
    c=colorbar;
    c.Label.String = 'Amplitude of 2nd qubit |1>';
    c.Label.FontSize=fs;
end
