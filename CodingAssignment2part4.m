clc 
clear all
h =0.0001;
tf=6;
t =0:h:tf;          % time scale

global V            % air velocity
                    

initial_h           = 0.1;  
initial_alpha       = 0.1;
initial_dhdt        = 0;
initial_dalphadt    = 0;

V_in  = 180;
V_fin = 185;
dV    = 1;

error = 0.000000001;

while V_fin-V_in >= error

    for V = V_in:dV:V_fin


        [t,y]=ode45( @aee, t,[initial_h initial_alpha initial_dhdt initial_dalphadt] );

        j = 0;
        for i=1:tf/h-1
            if((y(i+1,1) > y(i+2,1)) && (y(i+1,1) > y(i,1)))
                j = j+1;
                L_max(j) = y(i+1,1);
            end

        end

%         if (L_max(j-1) > L_max(j))
%             output=[double(V) "Stable"];
%             disp(output)
%         else
%             output=[double(V) "Unstable"];
%             disp(output)
%         end

        if (L_max(j-1) < L_max(j))
            V_in  = V-dV;
            V_fin = V;
            dV    = dV*0.1;
            break;
        end

    end

end
output=["Flutter Velocity" double(V-dV)];
disp(output)

    function dydt=aee(t,y)
        global V;
        % System parameters
        s = 7.5;            % semi span
        c = 2;              % chord
        m = 100;            % unit mass / area of wing
        kappa_freq = 5;     % flapping freq in Hz
        theta_freq = 10;    % pitch freq in Hz
        xcm = 0.5*c;        % position of centre of mass from nose
        xf = 0.48*c;        % position of flexural axis from nose
        e = xf/c - 0.25;    % eccentricity between flexural axis and aero centre (1/4 chord)
        a = 2*pi;           % 2D lift curve slope
        rho = 1.225;        % air density
        Mthetadot = -1.2;   % unsteady aero damping term

        % Set up system matrices
        % Inertia matrix
        a11=(m*s^3*c)/3;                        % I kappa
        a22= m*s*(c^3/3 - c*c*xf + xf*xf*c);    % I theta
        a12 = m*s*s/2*(c*c/2 - c*xf);           % I kappa theta
        a21 = a12;
        A=[a11,a12;a21,a22];
        

        % Aero/Structural stiffness matrix
        k1 = (kappa_freq*pi*2)^2*a11; % k kappa heave stiffness
        k2 = (theta_freq*pi*2)^2*a22; % k theta pitch stiffness
        E = [k1 0; 0 k2];
        K = (rho*V^2*[0,c*s^2*a/4; 0,-c^2*s*e*a/6])+E;

        % Aero Damping
        C = rho*V*[c*s^3*a/6,0;-c^2*s^2*e*a/4,-c^3*s*Mthetadot/8];


        X=-(A^(-1))*K;
        Y=-(A^(-1))*C;

        dydt_1 = y(3);
        dydt_2 = y(4);
        dydt_3 = X(1,1)*y(1)+X(1,2)*y(2)+Y(1,1)*y(3)+Y(1,2)*y(4);
        dydt_4 = X(2,1)*y(1)+X(2,2)*y(2)+Y(2,1)*y(3)+Y(2,2)*y(4);

        dydt=[dydt_1; dydt_2; dydt_3; dydt_4];
    end


    

