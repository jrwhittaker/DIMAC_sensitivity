% dimac_sensitivity_analysis.m
%
% Sensitivity analysis of spoiled GRE sensitivity to pulsatile changes in
% CBV and CBFV in large intracranial arteries. 
% 
% Requires in path:
% arterial_spoiledgre_signalmodel.m
% pulsatile_scaling_function.m
%
% Author: JR Whittaker, 4th June 2020.

clear
clc

%% Function handles
sigmodel    = @arterial_spoiledgre_signalmodel; % DIMAC signal model
psf         = @pulsatile_scaling_function; % Simulate pulsatile waveform
sumsq       = @(x) sum((x-mean(x)).^2);

%% Basic sensitivity analysis
% MR simulation parameters
FA = rad2deg(pi./2.^[1:4]);
TR = [10:5:100]./1000;
TE = 0; % Only consider longitudinal magnetisation


% Physiological baseline and pulsatile ranges
CBFV_base   = 0;
CBFV_var    = 100;

CBV_base    = 0.5;
CBV_var     = 0.05*CBV_base;

% Sensitivity analysis
niter=1000; % number of iterations
sens_mat    = zeros(length(FA),length(TR));
Sdiff=nan(length(FA),length(TR));

for n = 1:length(FA)
    for m = 1:length(TR)
        
        if (m == 1)
        % Sensitivity analysis parameters
        CBFV_iter   = CBFV_base + CBFV_var.*rand(niter,1);
        CBV_iter    = CBV_base + CBV_var.*rand(niter,1);
        end
        
        % Reference (baseline) signal
        Sref = sigmodel(CBV_base,CBFV_base,FA(n),TR(m),TE);
        
        for ni = 1:niter
            
            % Signal for nth iteration
            Sn = sigmodel(CBV_iter(ni),CBFV_iter(ni),FA(n),TR(m),TE);
            
            % Deviation from reference
            Sdiff(ni,m) = (Sn - Sref);
            
        end
        
        % Desgin matrix for linear regression
        X(:,1) = CBV_iter-min(CBV_iter);
        X(:,1) = X(:,1)./max(X(:,1));
        X(:,2) = CBFV_iter-min(CBFV_iter);
        X(:,2) = X(:,2)./max(X(:,2));
        
        % Linear regression
        beta = regress(Sdiff(:,m),X);
        
        % Relative sensativity
        sens_mat(n,m) = beta(2)/beta(1);
        
    end
end

figure
subplot(1,2,1)
plot(1000.*TR,log10(sens_mat)','LineWidth',2);
axis('square',[10 100 -0.1 2]);
set(gca,'xtick',[20 40 60 80 100],'xticklabel',{20,40,60,80,100});
set(gca,'ytick',log10([1 10:10:100]),'yticklabel',...
    {'10^{0}','10^{1}','','','','','','','','','10^{2}'});
xlabel('TR (ms)','Interpreter','latex','FontSize',28);
ylabel('Relative sensitivity $\frac{{\Delta}S}{{\Delta}CBFV} / \frac{{\Delta}S}{{\Delta}CBV}$','Interpreter','latex', ...
    'FontSize',28);
hAxes.TickLabelInterpreter = 'latex';
lgd=legend({'$90^{\circ}$','$45^{\circ}$','$23^{\circ}$','$11^{\circ}$'},'Interpreter','latex','Box','off');
lgd.Location='northeastoutside';
lgd.FontSize=16;
set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',16);
grid on;

%% Effect of basal partial volume on sensitivity

% MR simulation parameters
FA = 90;
TR = 15./1000;
TE = 0; % Only consider longitudinal magnetisation

% Physiological baseline and pulsatile ranges
CBFV_base   = 0;
CBFV_var    = 100; %7*CBFV_base;

CBV_base    = [0.1:0.2:0.9];
CBV_var     = 0.05.*CBV_base;

% sensitivity analysis
niter=1000; % number of iterations
Sdiff=nan(niter,1);

for n=1:length(CBV_base)
    
    % Sensitivity analysis parameters
    CBFV_iter   = CBFV_base + CBFV_var.*rand(niter,1);
    CBV_iter    = CBV_base(n) + CBV_var(n).*rand(niter,1);
    
    % Reference (baseline) signal
    Sref = sigmodel(CBV_base(n),CBFV_base,FA,TR,TE);
    
    for ni = 1:niter
        
        % Signal for nth iteration
        Sn = sigmodel(CBV_iter(ni),CBFV_iter(ni),FA,TR,TE);
        
        % Deviation from reference
        Sdiff(ni,:) = (Sn - Sref);
        
    end
    
    subplot(1,2,2),np(n)=plot(CBFV_iter,Sdiff+Sref,'.'); hold on;
    
end

hold off
axis('square',[10 80 0 1]);
xlabel('CBFV (cm $s^{-1}$)','Interpreter','latex','FontSize',28);
ylabel('S (a.u)','Interpreter','latex', 'FontSize',28);
hAxes.TickLabelInterpreter = 'latex';
set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',16);
lgd=legend(np,{'0.1','0.3','0.5','0.7','0.9'},'Interpreter','latex','Box','off');
lgd.Location='northeastoutside';
lgd.FontSize=16;
set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',16);
grid on

   

        
%% Effect of SNR and partial volume

SNR=[25 50 75 100 125];
PV=[0.1 0.5 0.75 1];

% Physiological baseline and pulsatile ranges
CBFV_base   = 30;
CBFV_delta    = 1;

% pulsatile cbfv
CBFV_pulse = psf(CBFV_base,CBFV_delta);
S_pulse = nan(64,1);

trace_mat=nan(4,5*256+4);

for n=1:length(PV)
    for m=1:length(SNR)
        
        S_reps=[];
        for reps=1:4
            for t=1:64
                S_pulse(t,:) = sigmodel(PV(n),CBFV_pulse(t),90,0.015,0.0075,'SNR',SNR(m));
            end
            S_reps=[S_reps; S_pulse];
        end
        
        trace_mat(n,(m-1)*256+[1:256]+(m-1))=S_reps;
        
    end
end
       
figure
for p=1:4
    subplot(4,1,p)
    for n=1:4
        line([n*(256+1) n*(256+1)],[0 1]),hold on;
    end
    plot(trace_mat(p,:),'k','LineWidth',2),axis([1 5*256-5 0 1]);
    set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',12,'Box','On');
    ylabel(num2str(PV(p)),'Interpreter','latex','FontSize',16);
    if (p == 4)
        xt=[];
        for xi=1:5
            xt=[xt (xi-1)*(256+1)+128];
        end
        set(gca,'xtick',xt,'xticklabel',{'25','50','75','100','125'});
        xlabel('SNR','Interpreter','latex','FontSize',16);
    else
        set(gca,'xtick',[]);
    end
end

SNR=[25 50 75];
NAVES=[1 5 10 100];
R2mu=cell(1,3);
R2se=cell(1,3);
for s=1:3
    R2mu{s}=nan(length(PV),length(NAVES));
    R2se{s}=nan(length(PV),length(NAVES));
end

for s=1:length(SNR)
    for n=1:length(PV)
        
        for ave=1:length(NAVES)
            R2vals=nan(10,1);
            for ni=1:10
                S_reps=[];
                for na=1:NAVES(ave)
                    for t=1:64
                        S_pulse(t,:) = sigmodel(PV(n),CBFV_pulse(t),90,0.015,0.0075,'SNR',SNR(s));
                    end
                    S_reps=[S_reps S_pulse];
                end
                S_ave=mean(S_reps,2);
                [~,~,err]=regress(S_ave,CBFV_pulse);
                R2vals(ni,:)=1-sumsq(err)/sumsq(S_ave);
            end
            R2mu{s}(n,ave)=mean(R2vals);
            R2se{s}(n,ave)=std(R2vals)/sqrt(10);
        end
    end
end
    
figure
for n=1:4
    for m=1:3
    subplot(1,4,n),errorbar(1:4,R2mu{m}(n,:),R2se{m}(n,:),'LineWidth',2),hold on;
    end
    set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',12,'Box','On');
    axis([0.5 4.5 0 1.1]);
    set(gca,'xtick',[1 2 3 4],'xticklabel',{1,5,10,100});
    xlabel('Number of averages','Interpreter','latex','FontSize',16);
    set(gca,'ytick',[0 0.5 0.75 0.95 1],'yticklabel',{0 0.5 0.75 0.95 1});
    ylabel('$R^{2}$','Interpreter','latex','FontSize',16);
    title(['$PV = ' num2str(PV(n)) '$'],'Interpreter','latex','FontSize',16);
    lgd=legend({'25','50','75'},'Interpreter','latex','Box','off');
    lgd.Location='SouthEast';
    lgd.FontSize=14;
    grid on
end
    

    

        
       
        
        
        
        
        
        
        
        
        
        









        
        
        
        