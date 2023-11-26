% C:\ProgramData\MATLAB\SupportPackages\R2023a\toolbox\matlab\quantum\
s1_readdataset;
%f0 = [0.2339    0.6698    0.3565    0.2833    0.5844    0.3170];
%pt =


n=length(f0);
layer1=[];
for k=1:n, layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))]; end
np=2*nchoosek(n,2);
x0=i_randpmpi(np,1);
theta0=x0;

    a=nchoosek(1:n,2);
    layer2=[];
    c=1;
    for k=1:size(a,1)
       layer2=[layer2; cryGate(a(k,1),a(k,2),theta0(c))];
       c=c+1;
       layer2=[layer2; cryGate(a(k,2),a(k,1),theta0(c))];
       c=c+1;
    end
    C = quantumCircuit([layer1; layer2]);


%%

methodid=4;

bestxa=0;
bestf1a=0;
bestkl=1;

for K=1:1
    K
    switch methodid
        case 4 
            options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon);
            lb=-0.5*pi*ones(np,1);
            ub=0.5*pi*ones(np,1);
            [xa,fval] = particleswarm(@i_obj,length(x0),lb,ub,options,pt,C,f0);            
        case 1        
            options = optimoptions('fmincon','Display','iter');
            lb=-0.5*pi*ones(np,1);
            ub=0.5*pi*ones(np,1);
            [xa,fval] = fmincon(@i_obj,x0,[],[],[],[], ...
                lb,ub,[],options,pt,C,f0);
        case 2
            options = optimset('Display','iter');
            [xa,fval] = fminsearch(@i_obj,x0,options,pt,C,f0);
        case 3
            options = optimset('Display','iter');
            [xa,fval] = fminunc(@i_obj,x0,options,pt,C,f0);
    end
    [poa,f1a,statev]=i_fullcirc_asym(xa,C);
    
    
    if i_kldiverg(f0,f1a) < bestkl
        bestf1a=f1a;
        bestxa=xa;
    end

end




function [y]=i_obj(x,pt,C,f0)
    [po,f1]=i_fullcirc_asym(x,C);
    y1=i_kldiverg(pt,po,true);
    %y2=i_kldiverg(f0,f1,false);
    y=y1; %+y2;
end

