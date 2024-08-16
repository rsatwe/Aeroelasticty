clc 
clear all
h =0.0001;
tf=1;
t =0:h:tf;          % time scale

%% IC-1
initial_h1           = 0.1;  
initial_dhdt1        = 0.2;
initial_alpha1       = 0.3;
initial_dalphadt1    = 0.4;
[t1,y1]=ode45( @aee, t,[initial_h1 initial_alpha1 initial_dhdt1 initial_dalphadt1] );

%% IC-2
initial_h2           = 0.5;  
initial_dhdt2        = 0.6;
initial_alpha2       = 0.7;
initial_dalphadt2    = 0.8;
[t2,y2]=ode45( @aee, t,[initial_h2 initial_alpha2 initial_dhdt2 initial_dalphadt2] );

%% Figures and Plots
figure(1)
subplot (4,1,1)
plot(t1,y1(:,1),'k',t2,y2(:,1),'b');
title('Responses at flutter airspeed')
xlabel('t(sec)'); ylabel('\kappa (rad)');
legend('IC-1','IC-2')
subplot (4,1,2)
plot(t1,y1(:,2),'k',t2,y2(:,2),'b');
xlabel('t(sec)'); ylabel('\alpha(rad)');
legend('IC-1','IC-2')
subplot (4,1,3)
plot(t1,y1(:,3),'k',t2,y2(:,3),'b');
xlabel('t(sec)'); ylabel('d\kappa/dt(rad/s)');
legend('IC-1','IC-2')
subplot (4,1,4)
plot(t1,y1(:,4),'k',t2,y2(:,4),'b');
xlabel('t(sec)'); ylabel('d\alpha/dt(rad/s)');
legend('IC-1','IC-2')

figure(2)
subplot(2,1,1)
plot(y1(:,1),y1(:,3),y2(:,1),y2(:,3))
hold on 
plot(initial_h1,initial_dhdt1,'o')
plot(initial_h2,initial_dhdt2,'o')
title('IVP Phase Portrait')
axis tight
xlabel('\kappa')
ylabel('d\kappa/dt')
legend('IC-1','IC-2','IC 1 coordinates','IC 2 coordinates')

subplot(2,1,2)
plot(y1(:,2),y1(:,4),y2(:,2),y2(:,4))
hold on 
plot(initial_alpha1,initial_dalphadt1,'o')
plot(initial_alpha2,initial_dalphadt2,'o')
axis tight
xlabel('\alpha')
ylabel('d\alpha/dt')
legend('IC-1','IC-2','IC 1 coordinates','IC 2 coordinates')

%% Time Period Calculation of Periodic Solution
xinv1=flip(y1(:,1));
TF1=islocalmax(xinv1);
k1=find(TF1==1);
k_T1=(k1(2)-k1(1))*h

xinv2=flip(y1(:,2));
TF2=islocalmax(xinv2);
k2=find(TF2==1);
k_T2=(k2(2)-k2(1))*h

xinv3=flip(y2(:,1));
TF3=islocalmax(xinv3);
k3=find(TF3==1);
k_T3=(k3(2)-k3(1))*h

xinv4=flip(y2(:,2));
TF4=islocalmax(xinv4);
k4=find(TF4==1);
k_T4=(k4(2)-k4(1))*h

%% Function

    function dydt=aee(t,y)
     % System parameters
        s           = 7.5           ;% semi span
        c           = 2             ;% chord
        m           = 100           ;% unit mass / area of wing
        kappa_freq  = 5             ;% flapping freq in Hz
        theta_freq  = 10            ;% pitch freq in Hz
        xcm         = 0.5*c         ;% position of centre of mass from nose
        xf          = 0.48*c        ;% position of flexural axis from nose
        e           = xf/c - 0.25   ;% eccentricity between flexural axis and aero centre (1/4 chord)
        a           = 2*pi          ;% 2D lift curve slope
        rho         = 1.225         ;% air density
        Mthetadot   = -1.2          ;% unsteady aero damping term
        
        
      % INPUT
       
      U           = 182.4768        ;% FLUTTER SPEED
      %U           = 190 ;
      
      % Set up system matrices
      % Inertia matrix
        m11=(m*s^3*c)/3;                        % I kappa
        m22= m*s*(c^3/3 - c*c*xf + xf*xf*c);    % I theta
        m12 = m*s*s/2*(c*c/2 - c*xf);           % I kappa theta
        m21 = m12;
        M=[m11 m12; m21 m22];
        T=M^(-1);
        
      % Structural stiffness matrix (heave stiffness varying)
        k1  = (kappa_freq*pi*2)^2*m11; % k kappa heave stiffness
        k2  = (theta_freq*pi*2)^2*m22; % k theta pitch stiffness
        K   = [k1 0;0 k2] ;     
        
      % Aero stiffness matrix (multiplied by U^2)
        ka12 = rho*c*s^2*a/4;
        ka22 = -rho*c^2*s*e*a/6;
        Ka   = [0 ka12; 0 ka22] ;
        
      % Aero Damping (multiplied by U)
        c11 = rho*c*s^3*a/6;
        c21 = -rho*c^2*s^2*e*a/4;
        c22 = -rho*c^3*s*Mthetadot/8;
        C   = [c11 0; c21 c22];
        
        X   = -T*Ka;
        Y   = -T*K;
        Z   = -T*C;
        
        dydt_1 = y(3);
        dydt_2 = y(4);
        
        dydt_3 =    (X(1,1)*(U^2)+Y(1,1))*y(1)+...
                    (X(1,2)*(U^2)+Y(1,2))*y(2)+...
                    Z(1,1)*U*y(3)+Z(1,2)*U*y(4);
        dydt_4 =    (X(2,1)*(U^2)+Y(2,1))*y(1)+...
                    (X(2,2)*(U^2)+Y(2,2))*y(2)+...
                    Z(2,1)*U*y(3)+Z(2,2)*U*y(4);
        
        dydt=[dydt_1; dydt_2; dydt_3; dydt_4];
    end