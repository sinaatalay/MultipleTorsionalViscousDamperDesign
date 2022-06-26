%% *Multiple Torsional Viscous Damper Design - ME 425 - Spring 2022*
% *NOTE:* The MATLAB code is divided into sections and does not depend on each 
% other so that you can run the code section by section.
%% Part 1
clc
clear
close
tic

% User Inputs:
n=5;

% *Program 1:*
Inertia = 100/n;
Stiffness = 25*n;

% *Configuring the matrices:*
I=eye(n)*Inertia;
K=zeros(n);
for i=1:n
    if(i==1)
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    elseif(i==n)
        K(i,i-1)=-Stiffness;
        K(i,i)=Stiffness;
    else
        K(i,i-1)=-Stiffness;
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    end
end

% *Modal Analysis:*
IToTheMinusHalf=I^(-0.5);
KTilde=IToTheMinusHalf*K*IToTheMinusHalf;
[P,Lambda]=eig(KTilde);
ModeShapes=P;
NaturalFrequencies=sqrt(Lambda*ones(size(I,1),1));

% *Printing Results:*
str="";
for i=1:n
    str=str+"Modeshape " +int2str(i) + " ("+num2str(NaturalFrequencies(i),'%.4f')+" rad/s)"+" is ";
    for j=1:n
        if(j==n)
            str=str+num2str(NaturalFrequencies(i),'%+.4f')+".\n";
        else
            str=str+num2str(NaturalFrequencies(i),'%+.4f')+", ";
        end
    end
end
str=str+"The elapsed time of Program 1 is "+num2str(toc,'%.4f')+" seconds.";
fprintf(str);

% *Plots:*
figure;
%Plot Settings:
colors={'r','g','b','c','m'};
xlim([0,n+1]);
ylim([0, n+1]);
xticks(1:n);
yticks(1:n);
xlabels= cell(1,n);
ylabels= cell(1,n);
for i=1:n
    xlabels{i}="\omega_{n"+int2str(i)+"}="+sprintf('%.3f',NaturalFrequencies(i))+" rad/s";
    ylabels{i}="Disk "+int2str(i);
end
xticklabels(xlabels);
yticklabels(ylabels);
xlabel("Angular Position");
hold on;

ModeShapes=ModeShapes/1.2;
for i=1:n
    xline(i,'--','color',colors{i}); % The equilibrium position
    plot(i+ModeShapes(:,i),1:n,'o','color',colors{i});
    plot(i+ModeShapes(:,i),1:n,'color',colors{i});
end
hold off;
%% Part 2
clc
clear
close
tic

% User Inputs:
n=5;
ObservedDisk=5;

% *Program 2:*
Inertia = 100/n;
Stiffness = 25*n;

% *Configuring the matrices:*
I=eye(n)*Inertia;

K=zeros(n);
for i=1:n
    if(i==1)
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    elseif(i==n)
        K(i,i-1)=-Stiffness;
        K(i,i)=Stiffness;
    else
        K(i,i-1)=-Stiffness;
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    end
end

F=zeros(size(I,1),1);
F(1)=Stiffness; % Base excitation

% *Modal Analysis:*
IToTheMinusHalf=I^(-0.5);
KTilde=IToTheMinusHalf*K*IToTheMinusHalf;
[P,Lambda]=eig(KTilde);
S=IToTheMinusHalf*P;
NaturalFrequencies=sqrt(Lambda*ones(size(I,1),1));

ModalF=P'*IToTheMinusHalf*F;
ModalT=@(omega) ModalF./(NaturalFrequencies.^2-omega.^2);
T=@(omega) S*ModalT(omega);

% *Printing Results:*
fprintf("The elapsed time of Program 2 is "+num2str(toc,'%.4f')+" seconds.");

% *Plot:*
omega=linspace(0,1.5*max(NaturalFrequencies),1001);
T=abs(T(omega));

loglog(omega,T(ObservedDisk,:));

ylim([min(T(ObservedDisk,:)), 1000]);
xlim([10^-2, max(omega)])
xlabel("Angular Base Excitation Frequency (rad/s)");
ylabel("Transmissibility (\Theta_n/\phi)");
%% Part 3
clc
clear
close
tic

% User Inputs:
n=5;
mu=0.3;
ObservedDisk=5;
FirstDampedDisk=1;
SecondDampedDisk=5;

% *Program 3:*
Inertia = 100/n;
Stiffness = 25*n;

% *Configuring the matrices:*
K=zeros(n);
for i=1:n
    if(i==1)
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    elseif(i==n)
        K(i,i-1)=-Stiffness;
        K(i,i)=Stiffness;
    else
        K(i,i-1)=-Stiffness;
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    end
end

I=eye(n+2)*Inertia;
I(n+1,n+1)=mu/2;
I(n+2,n+2)=mu/2;
Iinv=I^-1;

K=[K zeros(n,2)];
K=[K;zeros(2,n+2)];

C=zeros(n+2);

% *Optimization:*

objectiveFunction=@(x) GetPeakT_partC(x(1), x(2), FirstDampedDisk, SecondDampedDisk, I, Iinv, C, K, Stiffness, ObservedDisk);
x0=[0.15,0.15];
lowerBound=[1e-5,1e-5];
options = optimset('fminimax');
options.Display = "off";
[x, ~, fval] = fminimax(objectiveFunction,x0,[],[],[],[],lowerBound,[],[],options);

% *Printing Results:*
str="Optimized ca1 value is "+num2str(x(1),'%.4f')+" N*m*s/rad.\n";
str=str+"Optimized ca2 value is "+num2str(x(2),'%.4f')+" N*m*s/rad.\n";
str=str+"Peak transmissibility is "+num2str(fval,'%.4f')+".\n";
str=str+"The elapsed time of Program 3 is "+num2str(toc,'%.4f')+" seconds.";
fprintf(str);

% *Plot:*
[~, I, C, DampedNaturalFrequencies]=objectiveFunction(x);
omega=linspace(min(DampedNaturalFrequencies)/1.5,max(DampedNaturalFrequencies)*1.2,1001);
omega=[omega DampedNaturalFrequencies'];
omega=sort(omega);

T=zeros(size(omega));
for i=1:length(omega)
    T(i)=GetT(omega(i), I, C, K, Stiffness, ObservedDisk);
end

loglog(omega,T);

ylim([min(T), max(T)*1.2]);
xlim([min(omega), max(omega)])
xlabel("Angular Base Excitation Frequency (rad/s)");
ylabel("Transmissibility (\Theta_n/\phi)");
%% Part 4
clc
clear
close
tic

% User Inputs:
n=5;
mu=0.3;
ObservedDisk=5;

% *Program 4:*
Inertia = 100/n;
Stiffness = 25*n;

% *Configuring the matrices:*
K=zeros(n);
for i=1:n
    if(i==1)
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    elseif(i==n)
        K(i,i-1)=-Stiffness;
        K(i,i)=Stiffness;
    else
        K(i,i-1)=-Stiffness;
        K(i,i)=2*Stiffness;
        K(i,i+1)=-Stiffness;
    end
end

I=eye(n+2)*Inertia;

K=[K zeros(n,2)];
K=[K;zeros(2,n+2)];

C=zeros(n+2);

F=zeros(n+2,1);
F(1)=Stiffness;

% *Optimization:*
combinations=nchoosek(1:n,2);
fvalmin=99999;

options = optimset('fminimax');
options.Display = "off";

for i=1:size(combinations,1)
    FirstDampedDisk=combinations(i,1);
    SecondDampedDisk=combinations(i,2);
    objectiveFunction=@(x) GetPeakT_partD(x(1),x(2),x(3),x(4), FirstDampedDisk, SecondDampedDisk, I, C, K, Stiffness, ObservedDisk);
    x0=[0.5,0.5,mu/2,mu/2];
    A=[-1,-1,0,0];
    b=-10^(-6);
    Aeq=[0,0,1,1];
    beq=mu;
    lowerBound=[1e-5,1e-5,1e-5,1e-5];
    [x, ~, fval] = fminimax(objectiveFunction,x0,A,b,Aeq,beq,lowerBound,[],[], options);
    %Note the better result:
    if(fval<fvalmin)
        result=[x, FirstDampedDisk, SecondDampedDisk];
        fvalmin=fval;
    end
end

% *Printing Results:*
str="Optimized ca1 value is "+num2str(result(1),'%.4f')+" N*m*s/rad.\n";
str=str+"Optimized ca2 value is "+num2str(result(2),'%.4f')+" N*m*s/rad.\n";
str=str+"Optimized Ia1 value is "+num2str(result(3),'%.4f')+" kg*m^2.\n";
str=str+"Optimized Ia2 value is "+num2str(result(4),'%.4f')+" kg*m^2.\n";
str=str+"Optimized location of the first damped disk is "+int2str(result(5))+".\n";
str=str+"Optimized location of the second damped disk is "+int2str(result(6))+".\n";
str=str+"Peak transmissibility is "+num2str(max(fvalmin),'%.4f')+".\n";
str=str+"The elapsed time of Program 4 is "+num2str(toc,'%.4f')+" seconds.";
fprintf(str);

% *Plot:*
[~, I, C, DampedNaturalFrequencies]=objectiveFunction(x);
omega=linspace(min(DampedNaturalFrequencies)/1.5,max(DampedNaturalFrequencies)*1.2,1001);
omega=[omega DampedNaturalFrequencies'];
omega=sort(omega);

T=zeros(size(omega));
for i=1:length(omega)
    T(i)=GetT(omega(i), I, C, K, Stiffness, ObservedDisk);
end

loglog(omega,T);

ylim([min(T), max(T)*1.2]);
xlim([min(omega), max(omega)])
xlabel("Angular Base Excitation Frequency (rad/s)");
ylabel("Transmissibility (\Theta_n/\phi)");