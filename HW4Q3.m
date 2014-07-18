% Robert Ackerman


% Homework 4
% Problem 3
% November 23, 2013

% Gregory, Allan W. and Veall, Michael R. (1985) "Formulating Wald Tests of
% Nonlinear Restrictions". Econometrica Vol. 53, No. 6(Nov., 1985), pp.
% 1465-1468).

%% Step 1: Prelminary Settings
clear all; clc;
tic;
reps=1000;
nparam = 3;
nvec = (20:10:400);
empsizeg = zeros(size(nvec));   
empsizeh = zeros(size(nvec));
%% (a) n=20 & (b) Looping for Different n 

for k=1:39
    n=nvec(k);
% Fix seed for v1.
defaultStream = RandStream.getGlobalStream;
defaultStream.reset(0);

% DGP [X1,t X2,t]'= [0.6 0.3: 0.3 0.6][X1,t-1 X2,t-1]'+[V1,t V2,t]' 
v1 = randn(n,1);

% Fix seed for v2.
defaultStream = RandStream.getGlobalStream;
defaultStream.reset(1);
v2 = randn(n,1);

% X
x = zeros(n,3);
x(:,1) = 1;
x(1,2) = 0 + v1(1);
x(1,3) = 0 + v2(1);
for i=2:n
    x(i,2) = x(i-1,2)*0.6 + x(i-1,3)*0.3 + v1(i);
    x(i,3) = x(i-1,2)*0.3 + x(i-1,3)*0.6 + v2(i);
end    
   
    
% Model & Wald Test
wg = zeros(reps,1);
wh = zeros(reps,1);
sizeg = zeros(reps,1);
sizeh = zeros(reps,1);
pvalg = zeros(reps,1);
pvalh = zeros(reps,1);
wald = zeros(reps,2);
for i=1:1000
    eps = randn(n,1);   
    y = 1 + 10*x(:,2) + 0.1*x(:,3) + eps;
    b = ((x'*x) \ eye(3)) *x'*y;
    resid = y-x*b;
    b_var = (x'*x\eye(3))*(x'*((1/(n-nparam)).*resid'*resid)*x)*(x'*x \eye(3));
    g = (b(2) - 1/b(3));
    G = [0 1 1/(b(3)^2)];
    wg(i) = (b(2) - 1/b(3))^2/((b_var(2,2) - 2*b(3)^-2*b_var(2,3) +  + b(3)^-4*b_var(3,3)));
    pvalg(i) = 1 - chi2cdf(wg(i),1);
        if pvalg(i) < 0.05 
        sizeg(i) = 1; 
        else
        sizeg(i) = 0;
        end
    h = (b(2)*b(3)-1);
    H = [0 b(3) b(2)];
    wh(i) = ((b(2)*b(3)-1)^2)/(b(3)^2*b_var(2,2) + 2*b(2)*b(3)*b_var(2,3) + b(2)^2*b_var(3,3));
    pvalh(i) = 1 - chi2cdf(wh(i),1);
        if pvalh(i) < 0.05 
        sizeh(i) = 1; 
        else
        sizeh(i) = 0;
        end
   wald(i,:) = [wg(i), wh(i)];     
end
empsizeg(k) = mean(sizeg);
empsizeh(k) = mean(sizeh);
wald_size_n = zeros(length(nvec),2);
crit = chi2inv(1-0.05,1);
wald_size(1) = sum(wald(:,1)>crit)/reps;
wald_size(2) = sum(wald(:,2)>crit)/reps;
wald_size_n(n,:) = wald_size;
end

% Plot
figure(1)
q = plot(nvec,empsizeg, 'LineWidth', 1.5);
ylim([0 0.35]);
set(q,'Color',[40/255 163/255 151/255]);
hold on
p = plot(nvec, empsizeh, 'LineWidth', 1.5);
hold on
xlabel('Sample Size');
ylabel('Empirical Size');
title('Empirical Size Under Hg & Hh, Nominal Size = 0.05')
legend('H0: Hg', 'H0: Hh')
set(p,'Color',[245/255 106/255 67/255])


second = toc;
%% (c) Bootstrap 

% Settings
run = 1;

if run == 1;
B = 999;
eps_bs = zeros(n,B);
y_bs = zeros(n,B);
b_bs = zeros(nparam,B);
resid_bs = zeros(n,B);
SE_bs = zeros(nparam,B);
t_bs = zeros(1,B);
t_q = zeros(1,1);
size_bs = zeros(1,2);
for null = 1:2

data = [y x];
b0_0=1;
b1_0=10;
b2_0=0.1;
b_0 = [b0_0 b1_0 b2_0];
Obj = @(b_0)NLSObjHW4(b_0,data);
% for hg/hh nulls
    if null == 1
        con = @(b_0)NLS_Cons_g(b_0);
    else
        con = @(b_0)NLS_Cons_h(b_0);
    end
% fmincon for nls estimates
options = optimoptions('fmincon','Display','off');
b_nls = (fmincon(Obj,b_0,[],[],[],[],[],[],con,options))';
resid_nls = y - x*b_nls;  

% generate bs samples
    for s = 1:B
            eps_bs(:,s) = datasample(resid_nls,n);
            y_bs(:,s) = x*b_nls + eps_bs(:,s);
           
            b_bs(:,s) = x \ y_bs(:,s);
            resid_bs(:,s) = y_bs(:,s) - x*b_bs(:,s);
            b_bs_var = (x'*x \ eye(3))*(x'*((1/(n-nparam)).*resid_bs(:,s)'*resid_bs(:,s))*x)*(x'*x \ eye(3));
            t_bs(s) = (b_bs(:,s) - b_nls)'*(b_bs_var \ eye(3))*(b_bs(:,s) - b_nls);
    end
       
        % sort 
        t_bs = sort(t_bs);
        t_q(1,1) = quantile(t_bs,0.95);
       
        % test original t-statistic from part a using t_quant
        size_bs(1,null) = sum(wald(:,null)>t_q(1,1))/reps;
end
disp(['Empirical Size = ' num2str(wald_size)]);
disp(['Empirical Size with Bootstrap = ' num2str(size_bs)]);
    


end



