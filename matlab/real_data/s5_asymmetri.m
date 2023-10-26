% C:\ProgramData\MATLAB\SupportPackages\R2023a\toolbox\matlab\quantum\

showfig=true;

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



if showfig
    subplot(2,2,4)
    bar([pt poa])
    set(gca,'XTick',1:length(statev));
    set(gca,'XTickLabel',statev);
    ylabel('# of cells');
    xlabel('Expression pattern');
    legend({'Target','Observed'})
    title(sprintf('kl=%f',i_kldiverg(pt,poa,true)));
end


At=zeros(n);   % qGRN outupt
%a_old=a;
%a=nchoosek(1:n,2);
c=1;
for k=1:size(a,1)
    At(a(k,1),a(k,2))=bestxa(c);
    c=c+1;
    At(a(k,2),a(k,1))=bestxa(c);
    c=c+1;
end


if showfig
    figure;
    subplot(2,2,1)
        At2=At;
        At2(At<0)=0;
        plot(digraph(At2,genes,'omitselfloops'));   
        %imagesc(abs(A)+abs(A'))
        %imagesc(A); title('Ground truth')
    subplot(2,2,2)
        %imagesc(abs(At)+abs(At'))
        imagesc(At); title('qGRN')
    subplot(2,2,3)
        imagesc(sc_pcnet(X)); title('PCR')
    subplot(2,2,4)
        bar([f0' f1a'])
        %set(gca,'XTick',1:n);
        %set(gca,'XTickLabel',statev);
        %ylabel('Gene');
        xlabel('Activation frequency');
        legend({'Target','Observed'})
        title(sprintf('kl=%f',i_kldiverg(f0,f1a)));
end


At
i_kldiverg(f0,bestf1a) 
[bestpo,bestf1]=i_fullcirc_asym(bestxa,C);

isequal(bestf1a,bestf1)

i_kldiverg(bestpo,pt) 


function [y]=i_obj(x,pt,C,f0)
    [po,f1]=i_fullcirc_asym(x,C);
    y1=i_kldiverg(pt,po,true);
    %y2=i_kldiverg(f0,f1,false);
    y=y1; %+y2;
end

