
function [r_fir,n_spike,output_V,cISI1] = trial(D,V_DS,gg,g,tao,delay,matrix)
rng(1);

N_node=2;
step=0.01; v_rev_E=2.0;lamda=3;sita=-0.25;r=1;k=8;
%% Par
T_tot =2500; T0=500;
p=0.14;gc=2;
E_Na=50; E_K=-100; E_SL=-70; E_DL=-70;
g_Na=20; g_K=20; g_SL=2; g_DL=2;
Cm=2;beta_m=-1.2;beta_w=0;gama_m=18;gama_w=10;fhi_w=0.15;
Rand=rand(6,2);
V_S(1:N_node)=0; V_D(1:N_node)=0; V_S0(1:N_node)=0; V_D0(1:N_node)=0; w(1:N_node)=0; w0(1:N_node)=0;IsnyE(1:N_node)=0;

n_spike(1:N_node)=0.0;  flag_r(1:N_node)=0.0; firingtime(1:N_node)=0.0;ft_pre(1:N_node)=0.0; ISI(1:N_node)=0.0;
Nt=round(T_tot/step);
time=step:step:T_tot;
out_firtime(1:N_node,1:round(Nt))=0; output_V(1:N_node,1:round(Nt))=0; v2(1,2)=0;


%% time loop
for iii  = 1:Nt
       
    a=139.0;b=400;
    IsnyE=-matrix.*g.*((1./(1+exp(-k*(v2-sita))))'*(V_S0-v_rev_E));   

    Isny=sum(IsnyE,1);

    m_inf=0.5*(1+tanh((V_S-beta_m)./gama_m));
    w_inf=0.5*(1+tanh((V_S-beta_w)./gama_w));
    tao_w=1./cosh((V_S-beta_w)./(2*gama_w));

    I_DS=gc*(V_D0-V_S0+V_DS);
    I_ion_S=-g_Na.*m_inf.*(V_S0-E_Na)-g_K.*w.*(V_S0-E_K)-g_SL.*(V_S0-E_SL);
    I_ion_D=-g_DL.*(V_D0-E_DL);

    Rand=rand(2,N_node);
    g_w=(-4*D*step*log(Rand(1,:))).^0.5.*(cos(2*pi*Rand(2,:)));

    V_S=V_S0+step*(I_DS/p+I_ion_S)/Cm; %%+g_w
    V_D=V_D0+step*(Isny/(1-p)-I_DS/(1-p)+I_ion_D)/Cm;%%Isny/(1-p)
    w=w0+(fhi_w.*(w_inf-w0)./(tao_w))*step;
    if time(iii)>T0 %&& time(iii)<100
    n_spike(V_S0<0 & V_S>=0)=n_spike(V_S0<0 & V_S>=0)+1;
    end

    output_V(:,iii)=V_S;
    output_w(:,iii)=w;
    
    if time(iii)>delay+10
        delayNt=round(iii-(delay/step)); 
        v2(1,:)=output_V(1,delayNt);
    end
    
    
    g=0;
    if time(iii)>delay+10 
        g=gg;
    end
    
    if time(iii)>delay+10
        ft_pre(flag_r==0 & output_V(:,iii-round(delay/step))'>=0)=firingtime(1);
        firingtime(flag_r==0 & output_V(:,iii-round(delay/step))'>=0)=time(iii);
        flag_r(output_V(:,iii-round(delay/step))<0)=0; flag_r(output_V(:,iii-round(delay/step))>=0)=1;
        
    end
       
    if time(iii)>500
        if firingtime(1)==time(iii) 
            ISI(r)=firingtime(1)-ft_pre(1);
            r=r+1;
        end   
    end
    out_firtime(:,iii)=firingtime;

    n_pre=1; n_post=1;
    Isny1(iii)=IsnyE(n_post,n_pre);
    Isny2(iii)=Isny(1);
    output2(:,iii)=V_S;
    
    
    V_S0=V_S; V_D0=V_D; w0=w;
end

cISI={ISI(1:end)};

    aaa=delay;
    bbb=cell2mat(cISI);
    clear ISI
    ISI(2,:)=bbb;
    ISI(1,:)=aaa;     
    output_V=output_V(:,1:10:end);
    output_w=output_w(:,1:10:end);


cISI1={ISI(2,1:end)};
r_fir=n_spike*1000/(T_tot-T0);
result=length(find (r_fir(:)>0))/100;

 

