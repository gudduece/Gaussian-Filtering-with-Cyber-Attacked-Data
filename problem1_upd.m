clc
clear all;    
% close all;

%% Parameter Intialization
    rand('seed',1)
  randn('seed',15)

n=6;
step=400;
sigmaf=sqrt(25*0.001);
sigmaa=sqrt(0.8*0.001); 
Q=1*diag([sigmaf^2 sigmaf^2 sigmaf^2 sigmaa^2 sigmaa^2  sigmaa^2 ]);
sigmar=sqrt(0.9);

R=1*diag([sigmar^2 sigmar^2]);
x0=[200 500 1000 3 4 3];
% X=zeros(n,step);
P0=1*diag([20 20 20 0.5 0.5 0.5]);
T=0.25*10^-3;
delay=3  ;


 m=2*n;
chai=eye(n)*sqrt(n);
chai=[chai -chai];
wt=(1/m)*ones(1,m);

count=7  %Count 2:No filtering, 3:Gemetric delay filtering  1: Binomial delay filtering 2-delay 4: False data atack proposed 5: False Attack Abhinoy sir 6: MLCKF 7: MISSIng
f=1;          % f=0 :truth Measurment Delay data,  Z=Y. 1: No truth Measurment Delay data Z~=Y


mu=1*ones(1,2);
sigma_fd=4*eye(2);


 pna=0.2;
 
 
 pfa=0.8/3;
 pda=pfa;
  Pb= pda
 Pg=0.1;
 p=Pb;
 B_mean=1-pda;
 
for mc=1:500
    
    mc
   
    X(:,1)=x0;
for k=1:step
     X(:,k+1)=(eye(n)*X(:,k))+mvnrnd(zeros(1,n),Q)';  % True state
    a=0;b=0;
    for j=1:3
        a=a+X(j+3,k+1)*cos(2*pi*X(j,k+1)*(k+1)*T);
        b=b+X(j+3,k+1)*sin(2*pi*X(j,k+1)*(k+1)*T);
    end
   Z(:,k+1)=[a b]'+mvnrnd(zeros(1,2),R)';              % True measurment state
     
end


   

 
 
    for k=1:step
        
           Fda=mvnrnd(mu,sigma_fd)'; 
        
        
        
        G=[];Gi=[];
        
         pnad = makedist('Binomial','N',1,'p', pna);
            pfad = makedist('Binomial','N',1,'p',pfa);
              pdad = makedist('Binomial','N',1,'p',pda);
        
               alpha = random(pnad,1);
                beta = random(pfad,1);
                gamma = random(pdad,1);
         
         pd = makedist('Binomial','N',1,'p',Pg);
    
     
         
        G = random(pd,k,1);
        sumz=0;
        for i=1:k
            if (G(i)==1);
          sumz=sumz+G(i)*Z(:,k+1-i);
          Gi(i)=G(i);
           break;
            else
              sumz=sumz+G(i)*Z(:,k+1-i);
  
            end
         end
        
        Y(:,k+1)=(1-alpha*beta)*Z(:,k+1)+alpha*(1-beta)*Fda +alpha*beta*(1-gamma)*sumz;
        
    end
    
  
%  
    if (f==0)
    Z=Y;
    end
    
   
Pest=P0;
x1=[190 900 2100 4.5 4.5 2.5]';
Xest(:,1)=mvnrnd(x0,P0); % Intial Estimated state
Px=zeros(n,2,1);
Py=zeros(2,2,1);

for k=1:step
    
 if count==7
     F=eye(n);
    X_pred(:,k+1)=F*Xest(:,k);
    sum_a=0;
   sum_f=0;
   for b=1:3
       sum_a=sum_a+X_pred(b+3,k+1)*cos(2*pi*X_pred(b,k+1)*(k+1)*T);
       sum_f=sum_f+X_pred(b+3,k+1)*sin(2*pi*X_pred(b,k+1)*(k+1)*T);
   end
   
   Y_pred(:,k+1)=B_mean*[sum_a;sum_f];
   omet(k+1)=2*pi*(k+1)*T;
 
  H=[-X_pred(4,k+1)*omet(k+1)*sin(omet(k+1)*X_pred(1,k+1)) -X_pred(5,k+1)*omet(k+1)*sin(omet(k+1)*X_pred(2,k+1)) -X_pred(6,k+1)*omet(k+1)*sin(omet(k+1)*X_pred(3,k+1)) cos(omet(k+1)*X_pred(1,k+1)) cos(omet(k+1)*X_pred(2,k+1)) cos(omet(k+1)*X_pred(3,k+1));
      X_pred(4,k+1)*omet(k+1)*cos(omet(k+1)*X_pred(1,k+1)) X_pred(5,k+1)*omet(k+1)*cos(omet(k+1)*X_pred(2,k+1)) X_pred(6,k+1)*omet(k+1)*cos(omet(k+1)*X_pred(3,k+1)) sin(omet(k+1)*X_pred(1,k+1)) sin(omet(k+1)*X_pred(2,k+1)) sin(omet(k+1)*X_pred(3,k+1))];
      
P_pred=F*Pest*F'+Q;
  K=P_pred*H'*B_mean*(B_mean*H*P_pred*H'+R*B_mean+(B_mean-(B_mean)^2)*[sum_a;sum_f]*([sum_a;sum_f])')^(-1);
  Xest(:,k+1)=X_pred(:,k+1)+K*(Y(:,k+1)-Y_pred(:,k+1));
  Pest=(eye(n)-K*B_mean*H)*P_pred;
  
 else
   
    
    
    
    
    
    
    
 
% Time Update
S = chol((Pest),'lower');
xmean=0;
Pzz=0;
Pxz=0;
Pup=0;

for i=1:m
Xin(:,i)=S*chai(:,i)+Xest(:,k);
    Xinp(:,i)=(eye(n)*Xin(:,i));
    xmean=xmean+Xinp(:,i)*wt(i);
end
X_up=xmean;

for i=1:m
Pup=Pup+wt(i)*((Xinp(:,i)*Xinp(:,i)')-(xmean*(xmean)'));
end
Pup=Pup+Q;

%Measurment Update
Sm=chol(Pup,'lower');
zmean=0;


for i=1:m
Xin1(:,i)=Sm*chai(:,i)+X_up;

a1=0;b1=0;
 for j=1:n/2
        a1=a1+Xin1(j+n/2,i)*cos(2*pi*Xin1(j,i)*(k+1)*T);
        b1=b1+Xin1(j+n/2,i)*sin(2*pi*Xin1(j,i)*(k+1)*T);
 end
z(:,i)=[a1 b1]';
zmean=zmean+z(:,i)*wt(i);
end
 
 Z_m(:,k+1)=zmean;

for i=1:m
Pzz=Pzz+wt(i)*((z(:,i)*z(:,i)')-(zmean*(zmean)'));
Pxz=Pxz+wt(i)*((Xin1(:,i)*z(:,i)')-(X_up*zmean'));

end

Pzz=Pzz+R;


if (count==1)

Pzzi(:,:,k)=Pzz;
Pxzi(:,:,k)=Pxz;

y2(:,1)=Z_m(:,1);
Py(:,:,1)=Pzzi(:,:,1);
Px(:,:,1)=Pxzi(:,:,1);

    a=0;
    for i=0:delay-1
        if(i<k)
        a=a+p^i*Z_m(:,k+1-i);
        end
    end
    y2(:,k+1)=(1-p)*a+p^(delay)*y2(:,k);
   
     
    b=0;c=0;
    for i=0:delay-1
         if(i<k)
        b=b+p^i*Pzzi(:,:,k-i);
        c=c+(p^i*(1-p^i*(1-p)))*(Z_m(:,k+1-i)*Z_m(:,k+1-i)');
         end
    end
    Py(:,:,k+1)=(1-p)*b+(1-p)*c+p^(delay)*Py(:,:,k);
   
   
     Pyy= Py(:,:,k+1);
     
   d=0;
 for i=0:delay-1
      if(i<k)
        d=d+p^i*Pxzi(:,:,k-i);
      end
 end
  Px(:,:,k+1)=(1-p)*d+p^(delay)* Px(:,:,k);
    
 Pxy= Px(:,:,k+1);
 
% Kalman Gain
 Kg1=Pxy*(inv(Pyy));

% Estimated state

Xest(:,k+1)=X_up+Kg1*(Y(:,k+1)-y2(:,k+1));

%Estimated Error Covariance
Pest=(Pup-(Kg1*Pyy*Kg1')); 


end


if (count==5)

Pzzi(:,:,k)=Pzz;
Pxzi(:,:,k)=Pxz;


    ym(:,k+1)= Z_m(:,k+1)+pfa*mu';
   
     
   
    Py(:,:,k+1)=Pzzi(:,:,k)+pfa*sigma_fd+pfa*(1-pfa)*mu'*mu;
   
   
     Pyy= Py(:,:,k+1);
     
  
    
 Pxy= Pxzi(:,:,k);
 
% Kalman Gain
 Kg1=Pxy*(inv(Pyy));

% Estimated state

Xest(:,k+1)=X_up+Kg1*(Y(:,k+1)-ym(:,k+1));

%Estimated Error Covariance
Pest=(Pup-(Kg1*Pyy*Kg1')); 


end


if (count==6)
    
 Sm1=chol(Pup,'lower');

for j=1:m
Xin1y(:,j)=Sm1*chai(:,j)+X_up;
end
    
for i=1:m 
a1=0;b1=0;

for j=1:n/2
        a1=a1+Xin1y(j+n/2,i)*cos(2*pi*Xin1y(j,i)*(k+1)*T);
        b1=b1+Xin1y(j+n/2,i)*sin(2*pi*Xin1y(j,i)*(k+1)*T);
end
 
y2(:,i,k)=[a1 b1]';
for i1=0:delay
        if(i1<k)
        y_hat(:,i,i1+1)=y2(:,i,k-i1);
        else 
           y_hat(:,i,i1+1)= y2(:,i,k);
        end
end

end

  
  for i=0:delay
     ymean=0;
     
for j=1:m
ymean=ymean+y_hat(:,j,i+1)*wt(j);
end
  
 y_hats_upd(:,i+1)=ymean;
  end
  

  
  
for i=0:delay
    
     Pyy=0;Pxy=0;
     
for j=1:m
Pyy=Pyy+wt(j)*((y_hat(:,j,i+1)*y_hat(:,j,i+1)')-(y_hats_upd(:,i+1)*(y_hats_upd(:,i+1))'));
Pxy=Pxy+wt(j)*((Xin1y(:,j)*y_hat(:,j,i+1)')-(X_up*y_hats_upd(:,i+1)'));
end

Pyys(:,:,i+1)=Pyy+R;
Pxys(:,:,i+1)=Pxy;
KG(:,:,i+1)=Pxys(:,:,i+1)*inv(Pyys(:,:,i+1));
end
  
  

  
  p_lam=0;
 tau=delay;
for i=0:delay
     p_lam=p_lam+tau^i;
end
 
for i=0:delay
     mud(i+1)=tau^i/p_lam;
end
     
  del=0;
for i=0:delay
   del=del+mud(i+1)*mvnpdf(Y(:,k+1),y_hats_upd(:,i+1),Pyys(:,:,i+1));
end
 

for i=0:delay
   X_hat(:,i+1)=X_up+KG(:,:,i+1)*(Y(:,k+1)-y_hats_upd(:,i+1));
   alpha(i+1)=mud(i+1)*mvnpdf(Y(:,k+1),y_hats_upd(:,i+1),Pyys(:,:,i+1))/del;
end
 
 sumx=0;
for i=0:delay
   sumx=sumx+alpha(i+1)*KG(:,:,i+1)*(Y(:,k+1)-y_hats_upd(:,i+1));
end
 
Xest(:,k+1)=X_up+sumx;

sump=0;
for i=0:delay
sump=sump+alpha(i+1)*(-KG(:,:,i+1)*Pyys(:,:,i+1)*KG(:,:,i+1)'+(Xest(:,k+1)-X_hat(:,i+1))*(Xest(:,k+1)-X_hat(:,i+1))');
end


%Estimated Error Covariance
Pest=Pup+sump; 

end







if (count==2)
% Kalman Gain
Kg=Pxz*(inv(Pzz));

% Estimated state

Xest(:,k+1)=X_up+Kg*(Z(:,k+1)-Z_m(:,k+1));

%Estimated Error Covariance
Pest=(Pup-(Kg*Pzz*Kg')); 
end


if(count==3)
    
Pzzi(:,:,k+1)=Pzz;
Pxzi(:,:,k+1)=Pxz;
sumcov=0;sumb=0;summ=0;sumd=0;


for i=1:min(k,10)
 
summ=summ+(1-Pg)^(i-1)*Pg*Z_m(:,k+1-i);    
    
sumcov=sumcov+(1-Pg)^(i-1)*Pzzi(:,:,k+1-i);

% sumd=sumd+Pb*Pg*(1-Pb)*(1-Pg)^(i-1)*Pzzi(:,:,k+1-i);
  
  
end

y1(:,k+1)=Pb*Z_m(:,k+1)+(1-Pb)*summ;

for i=1:min(k,10)
 
for j=1:min(k,10)
   
    if (i==j)
    sumb=sumb+((Pb*(1-Pg)^(i-1)*Pg)-(1-Pb)*(1-Pg)^(2*i-2)*(Pg^2)*(3*Pb-1))* Z_m(:,k+1-i)* Z_m(:,k+1-i)';
    else
    sumb=sumb+(Pb*Pg^2*(1-Pg)^(i+j-2)-2*Pb*(1-Pb)*Pg^2*(1-Pg)^(i+j-2)+(1-Pb)^2*(1-Pg)^(i+j-2)*Pg^2)*Z_m(:,k+1-i)* Z_m(:,k+1-j)';  
    end

end
end


Py(:,:,k+1)=Pb*Pzzi(:,:,k+1)+Pb*(1-Pb)*(Z_m(:,k+1)* Z_m(:,k+1)')+(1-Pb)*Pg*sumcov+sumb+2*sumd;
Pyy= Py(:,:,k+1);

sumc=0;

for i=1:min(k,10)
  sumc=sumc+Pb*(1-Pg)^(i-1)*Pg*Pxzi(:,:,k+1-i);
 end

Px(:,:,k+1)=Pb*Pxzi(:,:,k+1)+sumc;

    Pxy= Px(:,:,k+1);
    
    % Kalman Gain
 Kg1=Pxy*(inv(Pyy));

% Estimated state

Xest(:,k+1)=X_up+Kg1*(Y(:,k+1)-y1(:,k+1));

%Estimated Error Covariance
Pest=(Pup-(Kg1*Pyy*Kg1')); 
    
end   



if(count==4)
    
Pzzi(:,:,k+1)=Pzz;
Pxzi(:,:,k+1)=Pxz;
sumcov=0;sumb=0;summ=0;sumd=0;


for i=1:min(k,10)
   

summ=summ+(1-Pg)^(i-1)*Pg*Z_m(:,k+1-i);

sumcov=sumcov+(1-Pg)^(i-1)*Pzzi(:,:,k+1-i);
  
end

y1(:,k+1)= (1-pna*pfa)*Z_m(:,k+1)+pna*(1-pfa)*mu'+pna*pfa*(1- pda)*summ;

for i=1:min(k,10)
 
for j=1:min(k,10)
   
    if (i==j)
    sumb=sumb+(pna*pfa*(1-pda)*Pg^(i-1)*Pg)*(1-(pna*pfa*(1-pda)*Pg^(i-1)*Pg))* Z_m(:,k+1-i)* Z_m(:,k+1-i)';
    else
    sumb=sumb+(pna*pfa*(1-pda)*Pg^(i+j-2)*Pg^2)*(1-(pna*pfa*(1-pda)*Pg^(i+j-2)*Pg^2))*Z_m(:,k+1-i)* Z_m(:,k+1-j)';  
    end

end
end


Py(:,:,k+1)=(1-pna*pfa)*Pzzi(:,:,k+1)+(1-pna*pfa)*(1-(1-pna*pfa))*(Z_m(:,k+1)* Z_m(:,k+1)')+pna*(1-pfa)*sigma_fd+(1-pna*pfa)*(1-(1-pna*pfa))*mu*mu'+pna*pfa*(1-pda)*Pg*sumcov+sumb;
Pyy= Py(:,:,k+1);

sumc=0;

for i=1:min(k,10)
  sumc=sumc+(1-Pg)^(i-1)*Pxzi(:,:,k+1-i);
 end

Px(:,:,k+1)=(1-pna*pfa)*Pxzi(:,:,k+1)+pna*pfa*(1-pda)*Pg*sumc;

    Pxy= Px(:,:,k+1);
    
    % Kalman Gain
 Kg1=Pxy*(inv(Pyy));

% Estimated state

Xest(:,k+1)=X_up+Kg1*(Y(:,k+1)-y1(:,k+1));

%Estimated Error Covariance
Pest=(Pup-(Kg1*Pyy*Kg1')); 
    
end   
end









e1(mc,k)=(X(1,k)-Xest(1,k));
e2(mc,k)=(X(2,k)-Xest(2,k));
e3(mc,k)=(X(3,k)-Xest(3,k));
e4(mc,k)=(X(4,k)-Xest(4,k));
e5(mc,k)=(X(5,k)-Xest(5,k));
e6(mc,k)=(X(6,k)-Xest(6,k));
    
end
 
end


  if (count==1||count==4||count==3||count==5||count==6||count==7)
 
 for i=1:step
    rms=0;
for k1=1:mc
    rms=rms+e1(k1,i)*e1(k1,i);
end
rms_value_1f(i)=sqrt(rms/mc);
 end
 rmf1=mean(rms_value_1f);

 for i=1:step
    rms=0;
for k1=1:mc
    rms=rms+e2(k1,i)*e2(k1,i);
end
rms_value_2f(i)=sqrt(rms/mc);
end
rmf2=mean(rms_value_2f);
for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e3(k,i)*e3(k,i);
end
rms_value_3f(i)=sqrt(rms/mc);
end
rmf3=mean(rms_value_3f);

for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e4(k,i)*e4(k,i);
end
rms_value_4f(i)=sqrt(rms/mc);
end
rmf4=mean(rms_value_4f);

for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e5(k,i)*e5(k,i);
end
rms_value_5f(i)=sqrt(rms/mc);
end
rmf5=mean(rms_value_5f);
for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e6(k,i)*e6(k,i);
end
rms_value_6f(i)=sqrt(rms/mc);
end
rmf6=mean(rms_value_6f);

rmsef=(rms_value_1f+rms_value_2f+rms_value_3f)/3     
rmsea=(rms_value_4f+rms_value_5f+rms_value_6f)/3



ARMSE_freq= rmf1+ rmf2+ rmf3/3
ARMSE_amp= rmf4+ rmf5+ rmf6/3




 else
     
for i=1:step
    rms=0;
for k1=1:mc
    rms=rms+e1(k1,i)*e1(k1,i);
end
rms_value_1(i)=sqrt(rms/mc);
end
rm1=mean(rms_value_1);
 for i=1:step
    rms=0;
for k1=1:mc
    rms=rms+e2(k1,i)*e2(k1,i);
end
rms_value_2(i)=sqrt(rms/mc);
end
rm2=mean(rms_value_2);

for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e3(k,i)*e3(k,i);
end
rms_value_3(i)=sqrt(rms/mc);
end
rm3=mean(rms_value_3);

for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e4(k,i)*e4(k,i);
end
rms_value_4(i)=sqrt(rms/mc);
end
rm4=mean(rms_value_4);

for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e5(k,i)*e5(k,i);
end
rms_value_5(i)=sqrt(rms/mc);
end
rm5=mean(rms_value_5);

for i=1:step
    rms=0;
for k=1:mc
    rms=rms+e6(k,i)*e6(k,i);
end
rms_value_6(i)=sqrt(rms/mc);
end
rm6=mean(rms_value_6);
rmsef=(rms_value_1+rms_value_2+rms_value_3)/3;
rmsea=(rms_value_4+rms_value_5+rms_value_6)/3;

ARMSE_freq= rm1+ rm2+ rm3/3
ARMSE_amp= rm4+ rm5+ rm6/3

  end
 

%   
%     if (count==1)
%    rmff=(rmf1+rmf2+rmf3)/3
%    rmffa=(rmf4+rmf5+rmf6)/3
%     end
%     if(count==2)
%   rmfa=(rm1+rm2+rm3)/3
%   rmfaa=(rm4+rm5+rm6)/3
%     end
%     
%     if (count==3)
%         rmffg=(rmf1+rmf2+rmf3)/3
%    rmffag=(rmf4+rmf5+rmf6)/3
%     end
%     
%     if (count==4)
%         rmfff=(rmf1+rmf2+rmf3)/3
%    rmffaf=(rmf4+rmf5+rmf6)/3
%     end
%     
%     comput.time=toc
%   figure(17)
%    hold on
%    if (count==1)
%      plot(pb,rmff,'r--',pb,rmffa,'b--')
%    
%    end
%    if (count==2)
%    
%        plot(pb,rmfa,'m--',pb,rmfaa,'g--')
%  
%        
%    end
%     if (count==3)
% plot(pb,rmffg,'k--',pb,rmffag,'y--')
%     end
%     
%  if (count==4)
% plot(pb,rmffg,'k--',pb,rmffag,'y--')
%     end

%     
%     
    
    figure(19)
 
   if (count==1)
    
     plot(rmsef,'k')
      hold on
     plot(rmsea,'g')
     
   end
   if (count==2)
   
       plot(rmsef,'r')
      hold on
     plot(rmsea,'b')
       
   end
    if (count==3)
plot(rmsef,'m')
      hold on
     plot(rmsea,'c')
    end
    

 if (count==4)
     plot(rmsef,'c--')
      hold on
     plot(rmsea,'k--')%  
 end
  

 if (count==5)
     plot(rmsef,'r*')
      hold on
     plot(rmsea,'k*')%  
 end
    

 if (count==6)
     plot(rmsef,'b-o')
      hold on
     plot(rmsea,'y-o')%  
 end
   

 if (count==7)
     plot(rmsef,'r-o')
      hold on
     plot(rmsea,'m-o')%  
end
    
%     
%     
%     
%     
%     
    figure(18)
plot(X(1,:))
hold on
plot(Xest(1,:))
 figure(11)
plot(X(2,:))
hold on
plot(Xest(2,:))
 figure(12)
plot(X(3,:))
hold on
plot(Xest(3,:))

figure(13)
plot(X(4,:))
hold on
plot(Xest(4,:))
 figure(14)
plot(X(5,:))
hold on
plot(Xest(5,:))
 figure(15)
plot(X(6,:))
hold on
plot(Xest(6,:))



    