% C:\ProgramData\MATLAB\SupportPackages\R2023a\toolbox\matlab\quantum\

n=length(f0);

layer1=[];
for k=1:n, layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))]; end


methodid=3;
np=2*nchoosek(n,2);
x0=i_randpmpi(np,1);

switch methodid
    case 1        
        options = optimoptions('fmincon','Display','iter');
        lb=-pi*ones(np,1);
        ub=pi*ones(np,1);
        [xa,fval] = fmincon(@i_obj,x0,[],[],[],[], ...
            lb,ub,[],options,pt,layer1,f0);
    case 2
        options = optimset('Display','iter');
        [xa,fval] = fminsearch(@i_obj,x0,options,pt,layer1,f0);
    case 3
        options = optimset('Display','iter');
        [xa,fval] = fminunc(@i_obj,x0,options,pt,layer1,f0);
end

[poa,f1a]=i_fullcirc(xa,layer1);

subplot(2,2,4)
bar([pt poa])
set(gca,'XTick',1:length(states));
set(gca,'XTickLabel',states);
ylabel('# of cells');
xlabel('Expression pattern');
legend({'Target','Observed'})
title(sprintf('kl=%f',i_kldiverg(pt,poa,true)));

At=zeros(n);
c=1;
a=nchoosek(1:n,2);
for k=1:size(a,1)
    At(a(k,1),a(k,2))=xa(c);
    c=c+1;
    At(a(k,2),a(k,1))=xa(c);
    c=c+1;
end


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
    %set(gca,'XTickLabel',states);
    %ylabel('Gene');
    xlabel('Activation frequency');
    legend({'Target','Observed'})
    title(sprintf('kl=%f',i_kldiverg(f0,f1a)));


function [y]=i_obj(x,pt,layer1,f0)
    [po,f1]=i_fullcirc(x,layer1);
    y1=i_kldiverg(pt,po,true);
    y2=i_kldiverg(f0,f1,false);
    y=y1+y2;
end

