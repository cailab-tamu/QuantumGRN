function [po,f1,statev]=i_fullcirc_symm(theta0,C)
    np=length(theta0);
    n=(1+sqrt(1+8*np))/2;
    assert(nchoosek(n,2)==np)

    a=nchoosek(1:n,2);
    c=1;
    for k=1:size(a,1)
        C.Gates(n+c) = cryGate(a(k,1),a(k,2),theta0(c));
        C.Gates(n+c) = cryGate(a(k,2),a(k,1),theta0(c));
        c=c+1;
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


% https://www.mathworks.com/matlabcentral/answers/1987423-change-of-parameters-of-gates-in-quantumcircuit-does-not-take-effect?s_tid=mlc_ans_email_view#answer_1265263