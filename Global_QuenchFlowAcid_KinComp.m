%% Quench Flow Acid Quench Num. Fit
close all
clearvars
[filename1,path1] = uigetfile('','Select Quench-Flow Data');
datapath1 = strcat(path1,filename1); %QuenchFlow Data
[filename2,path2] = uigetfile('','Select Kinetic Comp Data');
datapath2 = strcat(path2,filename2); %KinCompData
nit = 100;
img = imread('QuenchFlowScheme.png');
imshow(img);
title('Scheme being modeled')
%% Input Variables - QuenchFlow
weight = 10; %weight for Quench Flow data
E_tot = 18; %initial enzyme concentration
ds1 = readtable(datapath1);
weight_1 = length(ds1.Time);
ds1 = repmat(ds1,weight,1);
titconc = unique(ds1(:,1));
titconcRange = table2array(titconc);
S_tot = titconcRange;
%k12 = 0.1;
%k21 = 4.1;
%% Input Variables - KinComp
ds2 = readtable(datapath2);
Dbp5_input_concentration = 1;%input('What is the protein concentration?');
mADP_input_conc = 20;%input('What is the concentration of mADP?');
%% Set-up Initial concentration matrix
titconc2 = unique(ds2(:,1));
titconcRange2 = table2array(titconc2);
dsID = [ones(length(ds1.Time),1); 2.*ones(length(ds2.Time),1)];
for l = 1:nit
%% Rate Constants for Dif Eq
aa1 = 1;
bb1 = 0;
k12 = (bb1-aa1).*rand(1,1) + aa1;
%k12 = 0.1;
aa2 = 0 ;
bb2 = 100;
k21 = (bb2-aa2).*rand(1,1) + aa2;
%k21 = 5
if k12 < 0
    k12 = -1.*k12;
end
if k21 < 0
    k21 = -1.*k12;
end
a1 = 0;
b1 = 100;
k23 = (bb2-aa2).*rand(1,1) + aa2;
k32 = (bb2-aa2).*rand(1,1) + aa2;
k34 = (bb2-aa2).*rand(1,1) + aa2;
k43 = (bb2-aa2).*rand(1,1) + aa2;
offset_init = zeros(1,length(titconcRange2));
ks = [k12,k21,k23,k32,k34,k43,1,offset_init];
%% Fitting, etc.
lsqoptions = optimoptions(@lsqcurvefit, 'UseParallel', true, 'Display','iter','FinDiffRelStep',1e-3,'MaxFunctionEvaluations',5000);
%lsqoptions = optimoptions(@lsqcurvefit, 'UseParallel', false, 'Display','final');
ft = @(ks,t) fitkinetics(ks, t, ds1.titconc,weight_1,weight,E_tot,titconcRange,ds2.titconc,Dbp5_input_concentration,mADP_input_conc,titconcRange2,dsID);
B1(l,:) = ks;
%ks=[k12,Kd_step1,k23,k32,k34,KdGle1_1,k54,k65,k67,KdGle1_2,k78,KdGle1_3,q1,q2,offset_initial];
lb = [0.05,2,0,0,0,0,0,-10.*ones(1,length(titconcRange2))];
up = [0.3,10,100,2000,100,100,100,10.*ones(1,length(titconcRange2))];
[B3(l,:),Rsdnrm3,Rsd3,ExFlg3,OptmInfo3,Lmda3,Jmat3] = lsqcurvefit(ft,B1(l,:),[ds1.Time ; ds2.Time],[ds1.signal ; ds2.signal],lb,up,lsqoptions);
ci{l} = nlparci(B3(l,:),Rsd3,'jacobian',Jmat3);
end
%% Fit Function 
function S = fitkinetics(B, tobs1, titrant,weight_1,weight_2,E_tot,titconcRange,titrant2,Dbp5_input,mD_input,titconcRange2,datasetID)
%ks=[k12,Kd_step1,k23,k32,k34,KdGle1_1,k54,k65,k67,KdGle1_2,k78,KdGle1_3,q1,q2,offset_initial];
tobsv1 = tobs1(datasetID == 1);
tobs = tobsv1(1:weight_1);
tobs2 = tobs1(datasetID == 2);
kmD12 = 1.2;
kmD21 = 64;
kmD23 = 26;
kmD32 = 5.5;
k12 = B(1);
k21 = B(2);
k23 = B(3);
k32 = B(4);
k34 = B(5);
k43 = B(6);
q1 = B(7);
offsets = B(8:end);
ks = [k12,k21,k23,k32,k34,k43];
ks2 = [kmD12,kmD21,kmD23,kmD32,k12,k21,k23,k32,k34,k43];
%% Initial Conc. Matrix - Quench Flow
E_0 = E_tot.*ones(1,length(titconcRange));
S_0 = titconcRange;
ES_0 = zeros(1,length(titconcRange));
EP_0 = zeros(1,length(titconcRange));
P_0 = zeros(1,length(titconcRange));
initial = [E_0',S_0,ES_0',EP_0',P_0'];
titrant = titrant(1:weight_1);
%% Quench Flow
T = cell([1 length(titconcRange)]);
Sv = cell([1 length(titconcRange)]);
totalsig1 = cell([1 length(titconcRange)]);
m = zeros(1,length(titconcRange));
for i = 1:length(titconcRange)
    ft1 = @(tobs , initial ) DifEq(tobs , initial, ks);
    %options2 = odeset('NonNegative',1:length(initial(1,:)));
    options2 = odeset('NonNegative',length(initial(1,:)));
    try
        [T{i}, Sv{i}] = ode15s(ft1,tobs(titrant == titconcRange(i)),initial(i,1:end),options2);
    if ~all(size(T{i},1)==size(tobs(titrant == titconcRange(i)),1))
        for k = 1:length(titconcRange)
        Sv{k} = 1e10*ones(size(tobs(titrant == titconcRange),1),size(Sv{1},2));
        totalsig1{k} = Sv{1,i}(:,4)+Sv{1,i}(:,5);
        m(k) = length(totalsig1{1,k}(:,1));
        end
        break
    end
    catch ME
        ME.identifier
        ME.message
    end
    totalsig1{i} = Sv{1,i}(:,4)+Sv{1,i}(:,5);
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

%% Kin Comp
% Initial Conc. Matrix - KinComp
Dbp5Gle1 = Dbp5_input*ones(1,length(titconcRange2));
mD = mD_input*ones(1,length(titconcRange2));
Dbp5Gle1mD=zeros(1,length(titconcRange2));
Dbp5Gle1mD2=zeros(1,length(titconcRange2));
Dbp5Gle1titconc=zeros(1,length(titconcRange2));
Dbp5Gle1Prod = zeros(1,length(titconcRange2));
Prod = zeros(1,length(titconcRange2));
initial2 = [Dbp5Gle1',mD',titconcRange2,Dbp5Gle1mD',Dbp5Gle1mD2',Dbp5Gle1titconc',Dbp5Gle1Prod',Prod'];

T2 = cell([1 length(titconcRange2)]);
Sv2 = cell([1 length(titconcRange2)]);
totalsig2 = cell([1 length(titconcRange2)]);
m2 = zeros(1,length(titconcRange2));
for i = 1:length(titconcRange2)
    ft2 = @(tobs2 , initial2 ) dDbp5Gle1kincomp(tobs2 , initial2, ks2);
    %options3 = odeset('NonNegative',1:length(initial2(1,:)));
    options3 = odeset('NonNegative',length(initial2(1,:)));
    try
        [T2{i}, Sv2{i}] = ode15s(ft2,tobs2(titrant2 == titconcRange2(i)),initial2(i,1:end),options3);
    if ~all(size(T2{i},1)==size(tobs2(titrant2 == titconcRange2(i)),1))
        for k = 1:length(titconcRange2)
        Sv2{k} = 1e10*ones(size(tobs2(titrant2 == titconcRange2),1),size(Sv2{1},2));
        totalsig2{k} = q1.*Sv2{1,i}(:,4)+q1.*Sv2{1,i}(:,5)+offsets(i);
        m2(k) = length(totalsig2{1,k}(:,1));
        end
        break
    end
    catch ME
        ME.identifier
        ME.message
    end
    totalsig2{i} = q1.*Sv2{1,i}(:,4)+q1.*Sv2{1,i}(:,5) + offsets(i);
    m2(i) = length(totalsig2{1,i}(:,1));
end
j2 = 1;
cellindex2 = zeros(1,length(titconcRange2));
for i = 1:length(titconcRange2)
    if i == 1
    cellindex2(i)=1;
    else
    cellindex2(i) = m2(j2)+cellindex2(i-1);
    j2 = j2+1;
    end
end
for i = 1:length(totalsig2)
    if i < length(totalsig2)
    totalsigv2(cellindex2(i):(cellindex2(i+1)-1),1) = totalsig2{1,i}(:,1);
    else
    totalsigv2(cellindex2(i):(((cellindex2(i)+length(totalsig2{1,i}(:,1))))-1),1) = totalsig2{1,i}(:,1);    
    end
end
S = [repmat(totalsig,weight_2,1);totalsigv2];
end
%% Diff. Equation
function df = DifEq(tspan, initial, ks)
    k12 = ks(1);
    k21 = ks(2);
    k23 = ks(3);
    k32 = ks(4);
    k34 = ks(5);
    k43 = ks(6);
    
    E = initial(1);
    S = initial(2);
    ES = initial(3);
    EP = initial(4);
    P = initial(5);

    dE = -k12.*E.*S + k21.*ES + k34.*EP - k43.*E.*P;
    dS = -k12.*E.*S + k21.*ES;
    dES = k12.*E.*S - k21.*ES - k23.*ES + k32.*EP;
    dEP = k23.*ES - k32.*EP - k34.*EP + k43.*E.*P;
    dP = k34.*EP - k43.*E.*P; 
    
    df = [dE;dS;dES;dEP;dP];
end

%% Differential Equation Function
function df = dDbp5Gle1kincomp(tspan,initial,ks)
    k12 = ks(1); %mD binding to Dbp5Gle1 1st step
    k21 = ks(2); %mD unbinding from Dbp5Gle1 1st step
    k23 = ks(3); %forward isomerization of Dbp5Gle1mD
    k32 = ks(4); %reverse isomerization of Dbp5Gle1mD
    k34 = ks(5); %binding of ATP to Dbp5Gle1
    k43 = ks(6); %unbinding of ATP to Dbp5Gle1
    k45 = ks(7); %hydrolysis of ATP
    k54 = ks(8); %reverse hydrolysis of ATP
    k56 = ks(9); % Pi release
    k65 = ks(10);
    
    Dbp5Gle1 = initial(1);
    mD = initial(2);
    ATP = initial(3);
    Dbp5Gle1mD = initial(4);
    Dbp5Gle1mD2 = initial(5);
    Dbp5Gle1ATP = initial(6);
    Dbp5Gle1Prod = initial(7);
    Prod = initial(8);
    
    dDbp5Gle1 = -k12.*Dbp5Gle1.*mD + k21.*Dbp5Gle1mD - k34.*Dbp5Gle1.*ATP + k43.*Dbp5Gle1ATP + k56.*Dbp5Gle1Prod - k65.*Dbp5Gle1.*Prod;
    dmD = -k12.*Dbp5Gle1.*mD + k21.*Dbp5Gle1mD;
    dATP = -k34.*Dbp5Gle1.*ATP + k43.*Dbp5Gle1ATP;
    dDbp5Gle1mD = k12.*Dbp5Gle1.*mD - k21.*Dbp5Gle1mD - k23.*Dbp5Gle1mD + k32.*Dbp5Gle1mD2;
    dDbp5Gle1mD2 = k23.*Dbp5Gle1mD - k32.*Dbp5Gle1mD2;
    dDbp5Gle1ATP = k34.*Dbp5Gle1.*ATP - k43.*Dbp5Gle1ATP - k45.*Dbp5Gle1ATP + k54.*Dbp5Gle1Prod;
    dDbp5Gle1Prod = k45.*Dbp5Gle1ATP - k54.*Dbp5Gle1Prod - k56.*Dbp5Gle1Prod + k65.*Dbp5Gle1.*Prod;
    dProd = k56.*Dbp5Gle1Prod - k65.*Dbp5Gle1.*Prod;
  


df = [dDbp5Gle1;dmD;dATP;dDbp5Gle1mD;dDbp5Gle1mD2;dDbp5Gle1ATP;dDbp5Gle1Prod;dProd;];
end