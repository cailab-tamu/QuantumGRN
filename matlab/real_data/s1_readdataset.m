addpath('..');

% Pearson-residual normalised expression, as used for the published 6-gene
% result. Binarising at 0 then means "above the expected value". The raw
% matrix is NOT interchangeable here: on raw counts >0 means "any count
% detected", which gives quite different activation ratios (IRF4 0.670,
% REL 0.584 vs 0.349, 0.345) and hence a different encoder and a different
% target distribution.
% X=readmatrix('../../dataset/expr_matrix_pearsonresidual_7.txt');
T=readtable('../../dataset/expr_matrix_pearsonresidual_7.txt', ...
            'CollectOutput',true,'ReadVariableNames',false);
genes=string(T.Var1);
X=T.(T.Properties.VariableNames{2});

[M]=dec2bin(0:2^size(X,1)-1, size(X,1))-'0';   % all 2^n binary states, gene 1 = most significant bit
y=0+(X>0).';
[~,idx]=ismember(y,M,'rows');
t=tabulate(idx);
n=zeros(size(M,1),1);
n(t(:,1))=t(:,2);

figure;
subplot(2,2,1)
bar(n)
for k=1:size(M,1)
    txt{k}=sprintf('%d',M(k,:));
end
txt=string(txt);
set(gca,'XTick',1:size(M,1));
set(gca,'XTickLabel',txt);
ylabel('# of cells');
xlabel('Expression pattern');


% ---------------------------------------------------------
f0=sum(y>0)./size(X,2);   % per gene initial activate freq.
pt=n./size(X,2);          % target cell state freq.    
% ---------------------------------------------------------
