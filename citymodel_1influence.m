function [I,I2,I3,I4] = citymodel_1influence(T)
%研究单个城市的Plateau
Name = 'A';

N = 200000;
I = N * 0.005;  
S = N - I;
R = 0;
I2 = N * 0.005;  
S2 = N - I2;
R2 = 0;
I3 = N * 0.005;  
S3 = N - I3;
R3 = 0;
I4 = N * 0.005;  
S4 = N - I4;
R4 = 0;


r = 6;
r2 = 6;
r3 = 6;
r4 = 15;

k=5;
k2=5;
k3=5;
k4=5;

d=0;
d2=0;
d3=0;
d4=0;

t=0;
t2=0;
t3=0;
t4=0;

y=0;
y2=0;
y3=0;
y4=0;

beta = 0.015;                                                                   %传染概率
gama = 0.12;                                                                    %康复概率

rt1=(gama*N)/(beta*S);
rt2=(gama*N)/(beta*S2);

for idx = 1:length(T)-1
    S(idx+1) = S(idx) - r(idx)*beta*S(idx)*I(idx)/N;
    I(idx+1) = I(idx) + r(idx)*beta*S(idx)*I(idx)/N - gama*I(idx);
    R(idx+1) = R(idx) + gama*I(idx);
    r(idx+1) = r(idx);
    d(idx+1) = r(idx)*beta*S(idx)*I(idx)/N;
    rt1(idx+1)=(gama*N)/(beta*S(idx+1));
    
    S2(idx+1) = S2(idx) - r2(idx)*beta*S2(idx)*I2(idx)/N;
    I2(idx+1) = I2(idx) + r2(idx)*beta*S2(idx)*I2(idx)/N - gama*I2(idx);
    R2(idx+1) = R2(idx) + gama*I2(idx);
    r2(idx+1) = r2(idx);
    d2(idx+1) = r2(idx)*beta*S2(idx)*I2(idx)/N;
    rt2(idx+1)=(gama*N)/(beta*S2(idx+1));
    
    S3(idx+1) = S3(idx) - r3(idx)*beta*S3(idx)*I3(idx)/N;
    I3(idx+1) = I3(idx) + r3(idx)*beta*S3(idx)*I3(idx)/N - gama*I3(idx);
    R3(idx+1) = R3(idx) + gama*I3(idx);
    r3(idx+1) = r3(idx);
    d3(idx+1) = r3(idx)*beta*S3(idx)*I3(idx)/N;
    
    S4(idx+1) = S4(idx) - r4(idx)*beta*S4(idx)*I4(idx)/N;
    I4(idx+1) = I4(idx) + r4(idx)*beta*S4(idx)*I4(idx)/N - gama*I4(idx);
    R4(idx+1) = R4(idx) + gama*I4(idx);
    r4(idx+1) = r4(idx);
    d4(idx+1) = r4(idx)*beta*S4(idx)*I4(idx)/N;
    
    %防控降低
    if d(idx+1)<20
        t=t+1;
        if t>=14
            if abs(r(idx+1)-8)>0.0001
                r(idx+1)=r(idx+1)*(8/6)^(1/k);
            end
        end
    end
    if d2(idx+1)<20
        t2=t2+1;
        if t2>=14
            if abs(r2(idx+1)-10)>0.0001
                r2(idx+1)=r2(idx+1)*(10/6)^(1/k2);
            end
        end
    end
    if d3(idx+1)<40
        t3=t3+1;
        y3=0;
        if t3>=14
            if abs(r3(idx+1)-10)>0.0001
                r3(idx+1)=r3(idx+1)*(10/6)^(1/k3);
            end
        end
    end
    if d4(idx+1)<20
        t4=t4+1;
        y4=0;
        if t4>=14
            if abs(r4(idx+1)-15)>0.0001
                r4(idx+1)=r4(idx+1)*(15/6)^(1/k4);
            end
        end
    end
    
    %防控加强
    if d(idx+1)>=20
        t=0;
        if abs(r(idx+1)-6)>0.0001
            r(idx+1)=r(idx+1)*(6/8)^(1/k);
        end
    end    
    if d2(idx+1)>=20
        t2=0;
        if abs(r2(idx+1)-6)>0.0001
            r2(idx+1)=r2(idx+1)*(6/10)^(1/k2);
        end
    end    
    if d3(idx+1)>=40
        t3=0;
        y3=y3+1;
        if y3>=1
            if abs(r3(idx+1)-6)>0.0001
                r3(idx+1)=r3(idx+1)*(6/10)^(1/k3);
            end
        end
    end    
    if d4(idx+1)>=20
        t4=0;
        y4=y4+1;
        if y4>=1
            if abs(r4(idx+1)-6)>0.0001
                r4(idx+1)=r4(idx+1)*(6/15)^(1/k4);
            end
        end
    end    

end



plot(T,I,T,I2,T,I3,T,I4);grid on;
xlabel('day');ylabel('the number of existing infections')
legend('\theta=20,r_{0}=6,r_{L}=6,r_{H}=8,k=5','\theta=20,r_{0}=6,r_{L}=6,r_{H}=10,k=5','\theta=40,r_{0}=6,r_{L}=6,r_{H}=10,k=5','\theta=20,r_{0}=15,r_{L}=6,r_{H}=15,k=5')

end

