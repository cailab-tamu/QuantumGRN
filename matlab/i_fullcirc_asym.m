function [po,f1,statev]=i_fullcirc_asym(theta0,C0,n)
    C=C0;
    %n=size(layer1,1);
    %a=nchoosek(1:n,2);
    %layer2=[];
    %c=1;
    %for k=1:size(a,1)
    %    layer2=[layer2; cryGate(a(k,1),a(k,2),theta0(c))];
    %    c=c+1;
    %    layer2=[layer2; cryGate(a(k,2),a(k,1),theta0(c))];
    %    c=c+1;
    %end
    %C = quantumCircuit([layer1; layer2]);
    if nargin<3
        n=size(C.Gates,1)-size(theta0,1);
    end
    n1=n+1;
    n2=size(C.Gates,1);
    for k=n1:n2
        C.Gates(k).Angles=theta0(k-n);
    end

    S = simulate(C);
    [statev,po] = querystates(S);
    if nargout>1        
        f1=zeros(1,n);
        for k=1:n
            f1(k)=probability(S,k,"1");
        end
    end
end