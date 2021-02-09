function [I,I2,I3,I4] = citymodel_3influence(T)
%研究3个城市的Plateau
Name = 'A';
Name2 = 'B';
Name3 = 'C';
Name4 = 'EE';

N = 200000;
N2 = 200000;
N3 = 200000;
N4 = 100000;

p = 0.005;
p2 = 0.005;
p3 = 0.005;
p4 = 0.0001;

I = N * p;  
S = N - I;
R = 0;
I2 = N2 * p2;  
S2 = N2 - I2;
R2 = 0;
I3 = N3 * p3;  
S3 = N3 - I3;
R3 = 0;
I4 = N4 * p4;  
S4 = N4 - I4;
R4 = 0;

r = 6;
r2 = 6;
r3 = 6;
r4 = 6;

k = 5;
k2 = 5;
k3 = 5;
k4 = 5;

d = 0;
d2 = 0;
d3 = 0;
d4 = 0;

t=0;
t2=0;
t3=0;
t4=0;

y = 0;
y2 = 0;
y3 = 0;
y4 = 0;

fp = 0.05;	%人口流动比例
fp2 = 0.05;
fp3 = 0.05;
fi = 0;
fi2 = 0;
fi3 = 0;
fi4 = 0;

beta = 0.015;                                                                   %传染概率
gama = 0.12;                                                                    %康复概率
       

for idx = 1:length(T)-1
    if (fp*N>=fp2*N2) && (fp2*N2>=fp3*N3)
        f13=N3*fp3*0.8*(fp*N/(fp*N+fp2*N2));
        f23=f13*((fp2*N2)/(fp*N));
        f12=f23*((fp*N)/(fp3*N3));
    elseif fp*N>=fp3*N3 && fp3*N3>=fp2*N2
        f23=N2*fp2*0.8*(fp3*N3/(fp*N+fp3*N3));
        f12=f23*(fp*N/(fp3*N3));
        f13=f23*(fp*N/(fp2*N2));       
    elseif fp2*N2>=fp*N && fp*N>=fp3*N3
        f13=N3*fp3*0.8*(fp*N/(fp*N+fp2*N2));
        f23=f13*(fp2*N2/(fp*N));
        f12=f23*(fp*N/(fp3*N3));        
    elseif fp2*N2>=fp3*N3 && fp3*N3>=fp*N
        f13=N*fp*0.8*(fp3*N3/(fp3*N3+fp2*N2));
        f12=f13*(fp2*N2/(fp3*N3));
        f23=f13*(fp2*N2/(fp*N));        
    elseif fp3*N3>=fp*N && fp*N>=fp2*N2
        f23=N2*fp2*0.8*(fp3*N3/(fp*N+fp3*N3));
        f12=f23*(fp*N/(fp3*N3));
        f13=f23*(fp*N/(fp2*N2));        
    elseif fp3*N3>=fp2*N2 && fp2*N2>=fp*N
        f13=N*fp*0.8*(fp3*N3/(fp3*N3+fp2*N2));
        f12=f13*(fp2*N2/(fp3*N3));
        f23=f12*(fp3*N3/(fp*N));        
    end

    S(idx+1) = S(idx) - r(idx)*beta*S(idx)*I(idx)/N + f12*(S2(idx)/N2 - S(idx)/N) + f13*(S3(idx)/N3 - S(idx)/N) + (N*fp-f12-f13)*(S4(idx)/N4 - S(idx)/N);
    I(idx+1) = I(idx) + r(idx)*beta*S(idx)*I(idx)/N - gama*I(idx) + f12*(I2(idx)/N2 - I(idx)/N) + f13*(I3(idx)/N3 - I(idx)/N) + (N*fp-f12-f13)*(I4(idx)/N4 - I(idx)/N);
    R(idx+1) = R(idx) + gama*I(idx) + f12*(R2(idx)/N2 - R(idx)/N) + f13*(R3(idx)/N3 - R(idx)/N) + (N*fp-f12-f13)*(R4(idx)/N4 - R(idx)/N);
    r(idx+1) = r(idx);
    d(idx+1) = r(idx)*beta*S(idx)*I(idx)/N;
    fi(idx+1) = f12*(I2(idx)/N2 - I(idx)/N) + f13*(I3(idx)/N3 - I(idx)/N) + (N*fp-f12-f13)*(I4(idx)/N4 - I(idx)/N);
    %有效防控力度rt
    rt(idx+1)=(gama*N)/(beta*S(idx+1));
    p(idx+1) = I(idx+1)/N;
    
    S2(idx+1) = S2(idx) - r2(idx)*beta*S2(idx)*I2(idx)/N2 + f12*(S(idx)/N - S2(idx)/N2) + f23*(S3(idx)/N3 - S2(idx)/N2) + (N2*fp2-f12-f23)*(S4(idx)/N4 - S2(idx)/N2);
    I2(idx+1) = I2(idx) + r2(idx)*beta*S2(idx)*I2(idx)/N2 - gama*I2(idx) + f12*(I(idx)/N - I2(idx)/N2) + f23*(I3(idx)/N3 - I2(idx)/N2) + (N2*fp2-f12-f23)*(I4(idx)/N4 - I2(idx)/N2);
    R2(idx+1) = R2(idx) + gama*I2(idx) + f12*(R(idx)/N - R2(idx)/N2) + f23*(R3(idx)/N3 - R2(idx)/N2) + (N2*fp2-f12-f23)*(R4(idx)/N4 - R2(idx)/N2);
    r2(idx+1) = r2(idx);
    d2(idx+1) = r2(idx)*beta*S2(idx)*I2(idx)/N2;
    fi2(idx+1) = f12*(I(idx)/N - I2(idx)/N2) + f23*(I3(idx)/N3 - I2(idx)/N2) + (N2*fp2-f12-f23)*(I4(idx)/N4 - I2(idx)/N2);
    rt2(idx+1)=(gama*N2)/(beta*S2(idx+1));    
    p2(idx+1) = I2(idx+1)/N2;
    
    S3(idx+1) = S3(idx) - r3(idx)*beta*S3(idx)*I3(idx)/N3 + f13*(S(idx)/N - S3(idx)/N3) + f23*(S2(idx)/N2 - S3(idx)/N3) + (N3*fp3-f13-f23)*(S4(idx)/N4 - S3(idx)/N3);
    I3(idx+1) = I3(idx) + r3(idx)*beta*S3(idx)*I3(idx)/N3 - gama*I3(idx) + f13*(I(idx)/N - I3(idx)/N3) + f23*(I2(idx)/N2 - I3(idx)/N3) + (N3*fp3-f13-f23)*(I4(idx)/N4 - I3(idx)/N3);
    R3(idx+1) = R3(idx) + gama*I3(idx) + f13*(R(idx)/N - R3(idx)/N3) + f23*(R2(idx)/N2 - R3(idx)/N3) + (N3*fp3-f13-f23)*(R4(idx)/N4 - R3(idx)/N3);
    r3(idx+1) = r3(idx);
    d3(idx+1) = r3(idx)*beta*S3(idx)*I3(idx)/N3;
    fi3(idx+1) = f13*(I(idx)/N - I3(idx)/N3) + f23*(I2(idx)/N2 - I3(idx)/N3) + (N3*fp3-f13-f23)*(I4(idx)/N4 - I3(idx)/N3);
    rt3(idx+1)=(gama*N3)/(beta*S3(idx+1)); 
    p3(idx+1) = I3(idx+1)/N3;
    
    S4(idx+1) = S4(idx) - r4(idx)*beta*S4(idx)*I4(idx)/N4  + (N*fp-f12-f13)*(S(idx)/N - S4(idx)/N4) + (N2*fp2-f12-f23)*(S2(idx)/N2 - S4(idx)/N4) + (N3*fp3-f13-f23)*(S3(idx)/N3 - S4(idx)/N4);
    I4(idx+1) = I4(idx) + r4(idx)*beta*S4(idx)*I4(idx)/N4 - gama*I4(idx)  + (N*fp-f12-f13)*(I(idx)/N - I4(idx)/N4) + (N2*fp2-f12-f23)*(I2(idx)/N2 - I4(idx)/N4) + (N3*fp3-f13-f23)*(I3(idx)/N3 - I4(idx)/N4);
    R4(idx+1) = R4(idx) + gama*I4(idx)  + (N*fp-f12-f13)*(R(idx)/N - R4(idx)/N4) + (N2*fp2-f12-f23)*(R2(idx)/N2 - R4(idx)/N4) + (N3*fp3-f13-f23)*(R3(idx)/N3 - R4(idx)/N4);
    r4(idx+1) = r4(idx);
    d4(idx+1) = r4(idx)*beta*S4(idx)*I4(idx)/N4;
    fi4(idx+1) = (N*fp-f12-f13)*(I(idx)/N - I4(idx)/N4) + (N2*fp2-f12-f23)*(I2(idx)/N2 - I4(idx)/N4) + (N3*fp3-f13-f23)*(I3(idx)/N3 - I4(idx)/N4);
    rt4(idx+1)=(gama*N4)/(beta*S4(idx+1)); 
    p4(idx+1) = I4(idx+1)/N4;
    

%     if d(idx+1)<40
%         t=t+1;
%         y=0;
%         if t>=14
%             if abs(r(idx+1)-10)>0.0001
%                 r(idx+1)=r(idx+1)*(10/6)^(1/k);
%             end
%             if abs(fp-0.15) > 0.0000001
%                 fp=fp*(3)^(1/k);
%             end
%         end
%     end
    if d2(idx+1)<40
        t2=t2+1;
        y2=0;
        if t2>=14
            if abs(r2(idx+1)-10)>0.0001
                r2(idx+1)=r2(idx+1)*(10/6)^(1/k2);
            end
            if abs(fp2-0.15) > 0.0000001
                fp2=fp2*(3)^(1/k2);
            end
        end
    end
    if d3(idx+1)<20
        t3=t3+1;
        y3=0;
        if t3>=14
            if abs(r3(idx+1)-10)>0.0001
                r3(idx+1)=r3(idx+1)*(10/6)^(1/k3);
            end
            if abs(fp3-0.15) > 0.0000001
                fp3=fp3*(3)^(1/k3);
            end
        end
    end
%     if d4(idx+1)<10
%         t4=t4+1;
%         y4=0;
%         if t4>=14
%             if abs(r4(idx+1)-10)>0.0001
%                 r4(idx+1)=r4(idx+1)*(10/6)^(1/k4);
%             end
%         end
%     end

    
    %防控加强
%     if d(idx+1)>=40
%         t=0;
%         y=y+1;
%         if y>=1
%             if abs(r(idx+1)-6)>0.0001
%                 r(idx+1)=r(idx+1)*(6/10)^(1/k);
%             end
%             if abs(fp-0.05) > 0.0000001
%                 fp=fp*(1/3)^(1/k);
%             end
%         end
%     end    
    if d2(idx+1)>=40
        t2=0;
        y2=y2+1;
        if y2>=1
            if abs(r2(idx+1)-6)>0.0001
                r2(idx+1)=r2(idx+1)*(6/10)^(1/k2);
            end
            if abs(fp2-0.05) > 0.0000001
                fp2=fp2*(1/3)^(1/k2);
            end
        end
    end    
    if d3(idx+1)>=20
        t3=0;
        y3=y3+1;
        if y3>=1
            if abs(r3(idx+1)-6)>0.0001
                r3(idx+1)=r3(idx+1)*(6/10)^(1/k3);
            end
            if abs(fp3-0.05) > 0.0000001
                fp3=fp3*(1/3)^(1/k3);
            end
        end
    end    
%     if d4(idx+1)>=10
%         t4=0;
%         y4=y4+1;
%         if y4>=1
%             if abs(r4(idx+1)-6)>0.0001
%                 r4(idx+1)=r4(idx+1)*(6/10)^(1/k4);
%             end
%         end
%     end    


end


plot(T,I,T,I2,T,I3,T,I4);grid on;
xlabel('day');ylabel('the number of existing infections')

legend('  A:r=6,p=5%','  B:\theta=40,r_{0}=6,r_{L}=6,r_{H}=10,p_{0}=5%,p_{L}=5%,p_{H}=15%,k=5','  C:\theta=20,r_{0}=6,r_{L}=6,r_{H}=10,p_{0}=5%,p_{L}=5%,p_{H}=15%,k=5','EE:r=6');
end

