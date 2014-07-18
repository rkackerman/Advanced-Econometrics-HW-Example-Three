% Robert Ackerman


% Homework 4
% Problem 2
% November 24, 2013

% Cameron and Trivedi. Microeconometrics: Methods and Applications. Problem
% 9-3 page 335.

%% Step 1: Load Data and Preliminary Settings
clear all; clc;

%load('HW4Q2data.mat');
lnmed=csvread('870hw4q2.csv');
lnmed = sort(lnmed);
med = exp(lnmed);
h1 = 100;
h2 = 0.25;
n=5006;
lnmedhat = zeros(n,1);
medhat = zeros(n,1);

%% Step 2: Create Estimates

for i = 1:n
     uvec = (lnmed - lnmed(i)) / h2;
     uvec2 = (med - med(i)) / h1;
     den = sum(normpdf(uvec));
     den2 = sum(normpdf(uvec2));
     lnmedhat(i) = 1/(n*h2)*sum(normpdf(uvec));
     medhat(i) = 1/(n*h1)*sum(normpdf(uvec2));
     
end;


% for i=1:n
%     uveclevel = (levelexp - levelexp(i)) / h;
%     fhatlevel(i) = 1/(n*h)*sum(normpdf(uveclevel));
%     uveclog = (logexp - logexp(i)) / h;
%     fhatlog(i) = 1/(n*h)*sum(normpdf(uveclog));
% end
%% Step 3: Plots

figure(1)
a = plot(med,medhat,'LineWidth',1.5);
xlabel('Medical Expenditures'); xlim([med(1) med(4950)]);
title('Kernel Density Estimate of Medical Expenditures h=100')
set(a,'Color',[40/255 163/255 151/255]);
legend('Kernel Estimate')

figure(2)
b = plot(lnmed,lnmedhat,'LineWidth',1.5);
xlabel('Log Medical Expenditures');
xlim([lnmed(1) lnmed(5000)]);
title('Kernel Density Estimate of Log Medical Expenditures h=0.25')
set(b,'Color',[40/255 163/255 151/255]);
legend('Kernel Estimate')

figure(3)
hist(lnmed, 50);
xlabel('Log Medical Expenditures');
title('Histogram of Log Medical Expenditures')
c = findobj(gca,'Type','patch');
set(c,'FaceColor',[40/255 163/255 151/255],'EdgeColor','w');


figure(4)
mu = 6.2;
sd = 1.5;
ix = -3*sd+2:1e-3:3*sd+8; 
iy = pdf('normal', ix, mu, sd);
d = plot(ix,iy,'k--','LineWidth',1.5); 
set(d,'Color',[245/255 106/255 67/255])
hold on
e = plot(lnmed,lnmedhat,'b','LineWidth',1.5);
set(e,'Color',[40/255 163/255 151/255]);
title('Kernel Estimate and Normal Density h=0.25') 
xlabel('Log Medical Expenditures');
legend('Normal Density', 'Kernel Estimate')




