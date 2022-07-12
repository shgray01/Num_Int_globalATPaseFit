%%Dbp5 and Dbp5Gle1 binding mADP Two Step
close all
clearvars
[filename1,path1] = uigetfile; % f[Gle1] data
datapath1 = strcat(path1,filename1); 
[filename2,path2] = uigetfile; % f[mADP] data
datapath2 = strcat(path2,filename2); 
nit = 1000;
img = imread('Dbp5_Dbp5Gle1_bindmADP_2step.png');
imshow(img);
title('Scheme being modeled')
%% Get Time Course Data for fGle1
ds1 = readtable(datapath1);
Dbp5_input_concentration = 1; %input('What is the protein concentration?');
mADP_input_conc = 100; %input('What is the concentration of mADP?');
%% Get Time Course Data for fmADP
ds2 = readtable(datapath2);
Dbp5_input_concentration_2 = 0.5; %input('What is the protein concentration?');
Gle1_input_conc = 10; %input('What is the concentration of Gle1?');
%% Data for Irriversible Dissociation
Gle1_conc_ds2_lam1 = [1,2,3,5,10, 30, 44];
Gle1_conc_ds2_lam2 = [0,1,2,3,5 10, 30, 44];
lambda1 = [304.12892,120.21526,81.2142,34.25997,18.7,22.6,22.8];
lambda2 = [2.3,1.67,1.62,1.6,1.5,1.1,1.2,1.2];
relative_weigt_lambda2 = 1;
ds3ID = 3.*ones(length(Gle1_conc_ds2_lam2).*relative_weigt_lambda2,1); %Based on relative error between lambda1 and lambda2
dsID = [ones(length(ds1.Time),1); 2.*ones(length(Gle1_conc_ds2_lam1),1);ds3ID ;4.*ones(length(ds2.Time),1)];
%% Set-up Initial concentration matrix
titconc = unique(ds1(:,1));
titconcRange = table2array(titconc);
titconc2 = unique(ds2(:,1));
titconcRange2 = table2array(titconc2);
%% Rate Constants for Dif Eq
k21 = 400;
Kd_step1 = 102; %Kd for first step of Dbp5 bind mD
k23 = 120; %isomerization of HmD to HmD'
k32 = 0.4;
q1 = 1; %quantum yield HmD'
q2 = 1.5; %quantum yield HGmD and/or HGmD'
q1_2 = q1;
q2_2 = q2;
for l = 1:nit
%% Rate Constants for Dif Eq
a1 = 0;
b1 = 20;
a3 = 50;
b3 = 500;
a2 = 0;
b2 = 1;
k34 = (b3-1).*rand(1,1) + a1; %binding of apoDbp5 to Gle1
KdGle1_1 = 0.5;
k54 = 60; %unbinding of mD by Dbp5Gle1
k65 = 20; %reverse isomerization of HGmD' to HGmD
k67 = 900; %binding of Gle1 to HmD
KdGle1_2 = 0.5;
k78 = 900; %binding of Gle1 to HmD'
KdGle1_3 = 0.5;
offset_initial = ones(1,length(titconcRange)).*0;
offset_initial2 = ones(1,length(titconcRange2)).*0;
ks=[k21,Kd_step1,k23,k32,k34,KdGle1_1,k54,k65,k67,KdGle1_2,k78,KdGle1_3,q1,q2,offset_initial,q1_2,q1_2,offset_initial2];
%ks=[k21,Kd_step1,k23,k32,k34,k43,k54,k65,k67,k76,k78,k87,q1,q2,offset_initial];
%% Fitting, etc.
lsqoptions = optimoptions(@lsqcurvefit, 'UseParallel', true, 'Display','iter','FinDiffRelStep',1e-3,'MaxFunctionEvaluations',1.900000e+04);
%lsqoptions = optimoptions(@lsqcurvefit, 'UseParallel', false, 'Display','final');
ft = @(ks,t) fitkinetics(ks, t, ds1.titconc,Dbp5_input_concentration,mADP_input_conc,titconcRange,ds2.titconc,Dbp5_input_concentration_2,Gle1_input_conc,titconcRange2,dsID);
B1(l,:) = ks;
%ks=[k12,Kd_step1,k23,k32,k34,KdGle1_1,k54,k65,k67,KdGle1_2,k78,KdGle1_3,q1,q2,offset_initial];
lb = [350,80,83,0.05,0,0,0,0,800,0,800,0,0,0,zeros(1,length(offset_initial)),0,0,0.*ones(1,length(offset_initial2))];
up = [1000,125,120,2.6,1000,1,1000,1000,1000,1,1000,1,1000,1000,10.*ones(1,length(offset_initial)),1000,1000,10.*ones(1,length(offset_initial2))];
[B3(l,:),Rsdnrm3,Rsd3,ExFlg3,OptmInfo3,Lmda3,Jmat3] = lsqcurvefit(ft,B1(l,:),[ds1.Time; Gle1_conc_ds2_lam1'; repmat(Gle1_conc_ds2_lam2',relative_weigt_lambda2,1);ds2.Time],[ds1.signal ; lambda1'; repmat(lambda2',relative_weigt_lambda2,1); ds2.signal],lb,up,lsqoptions);
ci{l} = nlparci(B3(l,:),Rsd3,'jacobian',Jmat3);
end
%% Histogram Plot

%% Fit Function
function S = fitkinetics(B, tobs1, titrant,Dbp5,mD,titconcRange,titrant2,Dbp5_2,Gle1_2,mD_tit,datasetID)
tobs = tobs1(datasetID == 1);
tobs2 = tobs1(datasetID == 4);
Gle1_lam1 = tobs1(datasetID == 2);
Gle1_lam2 = tobs1(datasetID == 3);
%ks=[k12,Kd_step1,k23,k32,k34,KdGle1_1,k54,k65,k67,KdGle1_2,k78,KdGle1_3,q1,q2,offset_initial];
k21 = B(1);
Kd_step1 = B(2); %Kd for first step of Dbp5 bind mD
k12 = k21./Kd_step1; %mD binding to Dbp5Gle1 1st step
k23 = B(3); %isomerization of HmD to HmD'
k32 = B(4);
k34 = B(5); %binding of apoDbp5 to Gle1
KdGle1_1 = B(6);
k43 = KdGle1_1.*k34;
k54 = B(7); %unbinding of mD by Dbp5Gle1
k65 = B(8); %reverse isomerization of HGmD' to HGmD
k67 = B(9); %binding of Gle1 to HmD
KdGle1_2 = B(10);
k76 = KdGle1_2.*k67;
k78 = B(11); %binding of Gle1 to HmD'
KdGle1_3 = B(12);
k87 = KdGle1_3.*k78;
q1 = B(13); %quantum yield HmD'
q2 = B(14); %quantum yield HGmD and/or HGmD'
off_sets = B(15:15+length(titconcRange)-1);
q1_2 = B(15+length(titconcRange));
q2_2 = B(15+length(titconcRange)+1);
off_sets2 = B(15+length(titconcRange)+2:end);
k45 = (k43.*k54.*k67.*k12)./(k34.*k76.*k21); %constrained from detailed balance
k56 = (k23.*k76.*k65.*k78)./(k32.*k67.*k87); %constrained from detailed balance
ks = [k12,k21,k23,k32,k34,k43,k45,k54,k56,k65,k67,k76,k78,k87];
KdG = k43./k34;
%% Lambda functions
lam1_calc = (-0.5.*(-k54-k56-k65+((KdGle1_2.*(-k21-k23+k54+k56))./(Gle1_lam1+KdGle1_2))+((KdGle1_3.*(-k32+k65))./(Gle1_lam1+KdGle1_3))-(((((Gle1_lam1.^2).*(k54+k56+k65)+Gle1_lam1.*KdGle1_2.*(k21+k23+k65)+Gle1_lam1.*KdGle1_3.*(k32+k54+k56)+KdGle1_2.*KdGle1_3.*(k21+k23+k32)).^2)-(4.*(Gle1_lam1+KdGle1_2).*(Gle1_lam1.*k54+k21.*KdGle1_2).*(Gle1_lam1+KdGle1_3).*(Gle1_lam1.*k65+k32.*KdGle1_3)))./((Gle1_lam1+KdGle1_2).^2.*(Gle1_lam1+KdGle1_3).^2)).^0.5));
lam2_calc = (-0.5.*(-k54-k56-k65+((KdGle1_2.*(-k21-k23+k54+k56))./(Gle1_lam2+KdGle1_2))+((KdGle1_3.*(-k32+k65))./(Gle1_lam2+KdGle1_3))+(((((Gle1_lam2.^2).*(k54+k56+k65)+Gle1_lam2.*KdGle1_2.*(k21+k23+k65)+Gle1_lam2.*KdGle1_3.*(k32+k54+k56)+KdGle1_2.*KdGle1_3.*(k21+k23+k32)).^2)-(4.*(Gle1_lam2+KdGle1_2).*(Gle1_lam2.*k54+k21.*KdGle1_2).*(Gle1_lam2+KdGle1_3).*(Gle1_lam2.*k65+k32.*KdGle1_3)))./((Gle1_lam2+KdGle1_2).^2.*(Gle1_lam2+KdGle1_3).^2)).^0.5));
%% Time course 1 - fGle1
for i = 1:length(titconcRange)
    HG(i) = ((Dbp5 + titconcRange(i) + KdG)-sqrt((Dbp5 + titconcRange(i) + KdG).^2-4.*Dbp5.*titconcRange(i)))./2;
    H(i) = Dbp5 - HG(i);
    G(i) = titconcRange(i) - HG(i);
end
mADP = mD.*ones(1, length(titconcRange));
HmD = zeros(1, length(titconcRange));
HGmD = zeros(1, length(titconcRange));
HmD2 = zeros(1, length(titconcRange));
HGmD2 = zeros(1, length(titconcRange));
initial = [H',mADP',G',HmD',HG',HGmD',HmD2',HGmD2'];
T = cell([1 length(titconcRange)]);
Sv = cell([1 length(titconcRange)]);
totalsig1 = cell([1 length(titconcRange)]);
m = zeros(1,length(titconcRange));
ft1 = @(tobs , initial ) DifEq(tobs , initial, ks);
for i = 1:length(titconcRange)
    %options2 = odeset('NonNegative',1:length(initial(1,:)));
    options2 = odeset();
    try
        [T{i}, Sv{i}] = ode15s(ft1,tobs(titrant == titconcRange(i)),initial(i,1:end),options2);
    if ~all(size(T{i},1)==size(tobs(titrant == titconcRange(i)),1))
        for k = 1:length(titconcRange)
        Sv{k} = 1e10.*ones(size(tobs(titrant == titconcRange(k)),1),size(Sv{1},2));
        totalsig1{k} = q1.*Sv{1,k}(:,7)+0.*q2.*Sv{1,k}(:,6)+q2.*Sv{1,k}(:,8)+off_sets(k); %assumes last Gle1 states is fluorescent
        m(k) = length(totalsig1{1,k}(:,1));
        i
        end
        break
    end
    catch ME
        ME.identifier
        ME.message
    end
    totalsig1{i} = q1.*Sv{1,i}(:,7)+0.*q2.*Sv{1,i}(:,6)+q2.*Sv{1,i}(:,8)+off_sets(i);
    m(i) = length(totalsig1{1,i}(:,1));
end
j = 1;
cellindex = zeros(1,length(titconcRange));
for i = 1:length(titconcRange)
    if i == 1
    cellindex(i)=1;
    else
    cellindex(i) = m(j)+cellindex(i-1);
    j = j+1;
    end
end
for i = 1:length(totalsig1)
    if i < length(totalsig1)
    totalsig(cellindex(i):(cellindex(i+1)-1),1) = totalsig1{1,i}(:,1);
    else
    totalsig(cellindex(i):(((cellindex(i)+length(totalsig1{1,i}(:,1))))-1),1) = totalsig1{1,i}(:,1);    
    end
end
%% Time course 2 -  fmADP
HG_2 = ones(1, length(mD_tit)).*((Dbp5_2 + Gle1_2 + KdG)-sqrt((Dbp5_2 + Gle1_2 + KdG).^2-4.*Dbp5_2.*Gle1_2))./2;
H_2 = (Dbp5_2 - HG_2).*ones(1, length(mD_tit));
G_2 = (Gle1_2 - HG_2).*ones(1, length(mD_tit));
mADP_2 = mD_tit;
HmD_2 = zeros(1, length(mD_tit));
HGmD_2 = zeros(1, length(mD_tit));
HmD2_2 = zeros(1, length(mD_tit));
HGmD2_2 = zeros(1, length(mD_tit));
initial2 = [H_2',mADP_2,G_2',HmD_2',HG_2',HGmD_2',HmD2_2',HGmD2_2'];
T2 = cell([1 length(mD_tit)]);
Sv2 = cell([1 length(mD_tit)]);
totalsig2 = cell([1 length(mD_tit)]);
m2 = zeros(1,length(mD_tit));
for i = 1:length(mD_tit)
    %options2 = odeset('NonNegative',1:length(initial(1,:)));
    options2 = odeset();
    try
        [T2{i}, Sv2{i}] = ode15s(ft1,tobs2(titrant2 == mD_tit(i)),initial2(i,1:end),options2);
    if ~all(size(T2{i},1)==size(tobs2(titrant2 == mD_tit(i)),1))
        for k = 1:length(mD_tit)
        Sv2{k} = 1e10.*ones(size(tobs2(titrant2 == mD_tit(k)),1),size(Sv2{1},2));
        totalsig2{k} = q1_2.*Sv2{1,k}(:,7)+0.*q2_2.*Sv2{1,k}(:,6)+q2_2.*Sv2{1,k}(:,8)+off_sets2(k); %assumes both Gle1 states are equiliy fluorescent
        m2(k) = length(totalsig2{1,k}(:,1));
        i
        end
        break
    end
    catch ME
        ME.identifier
        ME.message
    end
    totalsig2{i} = q1_2.*Sv2{1,i}(:,7)+0.*q2_2.*Sv2{1,i}(:,6)+q2_2.*Sv2{1,i}(:,8)+off_sets2(i);
    m2(i) = length(totalsig2{1,i}(:,1));
end
j2 = 1;
cellindex2 = zeros(1,length(mD_tit));
for i = 1:length(mD_tit)
    if i == 1
    cellindex2(i)=1;
    else
    cellindex2(i) = m2(j2)+cellindex2(i-1);
    j2 = j2+1;
    end
end
for i = 1:length(totalsig2)
    if i < length(totalsig2)
    totalsig_2(cellindex2(i):(cellindex2(i+1)-1),1) = totalsig2{1,i}(:,1);
    else
    totalsig_2(cellindex2(i):(((cellindex2(i)+length(totalsig2{1,i}(:,1))))-1),1) = totalsig2{1,i}(:,1);    
    end
end
S = [totalsig; lam1_calc; lam2_calc; totalsig_2];
end
%% Diff. Equation
function df = DifEq(tspan, initial, ks)
    k12 = ks(1);
    k21 = ks(2);
    k23 = ks(3);
    k32 = ks(4);
    k34 = ks(5);
    k43 = ks(6);
    k45 = ks(7);
    k54 = ks(8);
    k56 = ks(9);
    k65 = ks(10);
    k67 = ks(11);
    k76 = ks(12);
    k78 = ks(13);
    k87 = ks(14);
    
    Dbp5 = initial(1);
    mADP = initial(2);
    Gle1 = initial(3);
    Dbp5mADP = initial(4);
    Dbp5Gle1 = initial(5);
    Dbp5Gle1mADP = initial(6);
    Dbp5mADP2 = initial(7);
    Dbp5Gle1mADP2 = initial(8);
    
    dDbp5 = -k12.*Dbp5.*mADP + k21.*Dbp5mADP - k34.*Dbp5.*Gle1 + k43.*Dbp5Gle1;
    dmADP = -k12.*Dbp5.*mADP + k21.*Dbp5mADP - k45.*Dbp5Gle1.*mADP + k54.*Dbp5Gle1mADP;
    dGle1 = -k34.*Gle1.*Dbp5 + k43.*Dbp5Gle1 - k67.*Dbp5mADP.*Gle1 + k76.*Dbp5Gle1mADP - k78.*Dbp5mADP2.*Gle1 + k87.*Dbp5Gle1mADP2;
    dDbp5Gle1 = k34.*Dbp5.*Gle1 - k43.*Dbp5Gle1 - k45.*Dbp5Gle1.*mADP + k54.*Dbp5Gle1mADP;
    dDbp5mADP = k12.*Dbp5.*mADP - k21.*Dbp5mADP - k23.*Dbp5mADP + k32.*Dbp5mADP2 - k67.*Dbp5mADP.*Gle1 + k76.*Dbp5Gle1mADP;
    dDbp5mADP2 = k23.*Dbp5mADP - k32.*Dbp5mADP2 - k78.*Dbp5mADP2.*Gle1 + k87.*Dbp5Gle1mADP2;
    dDbp5Gle1mADP = k45.*Dbp5Gle1.*mADP - k54.*Dbp5Gle1mADP - k56.*Dbp5Gle1mADP + k65.*Dbp5Gle1mADP2 - k76.*Dbp5Gle1mADP + k67.*Dbp5mADP.*Gle1;
    dDbp5Gle1mADP2 = k78.*Dbp5mADP2.*Gle1 - k87.*Dbp5Gle1mADP2 + k56.*Dbp5Gle1mADP - k65.*Dbp5Gle1mADP2;
    
    df = [dDbp5;dmADP;dGle1;dDbp5mADP;dDbp5Gle1;dDbp5Gle1mADP;dDbp5mADP2;dDbp5Gle1mADP2];
end