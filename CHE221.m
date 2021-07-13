% Substance = Toluene, Equation of State = Peng-Robinson
% All analysis is done for 1 mole of Toluene (N = 1).
clc
clear

% Specifics of Peng-Robinson equation:
global sigma;
sigma = 1 + sqrt(2);
global epsilon;
epsilon = 1 - sqrt(2);
global omega;
omega = 0.07780;
global psi; 
psi = 0.45724;
global R 
R = 8.314;  %Units = Pa.m^3/(mol.K) => SI units
global R_new_bar;
R_new_bar = 0.0832; %L.bar/mol.K  (0.0821(in L.atm/mol.K)*1.01325(bar/atm))


% Specifics of Toluene:

global Tc;
Tc = 591.8 ;    %In Kelvins 
global Pc;
Pc = 41.06e5 ;    %In Pascals (Pc = 41.06 bar)
global Pc_bar;
Pc_bar = Pc/1e5;  %Pc in bars.
global w;
w = 0.262 ;       % Omega term in alpha(Tr, w)

%Note that a, being dependent on temperature has a separate function for
%it, and can't be defined as a single value. We can do this for b though.
global b;   %The parameter of Toluene
b = omega*R*Tc/Pc;
global b_new;
b_new = b*1e3;  %SI units of m^3/mol to L/mol.

%--------------------------------------------------------------------------
% Characterisation of the system over. Now begins the main code.

P_sat_vec = []; %Vector for storing pressures
V_vec1 = []; %Vector for storing volumes of liq phase
V_vec2 = []; %Vector for storing volumes of vap phase

xlim([0,7]);
ylim([0, 40]);
xlabel("ʋ = Molar volume (in Litres)");
ylabel("P = Pressure (in bars)");
title("Pressure (P) v/s Molar volume (ʋ) graph - for Toluene by Peng-Robinson EOS");
hold on;

%PLOTTING THE DOME SHAPED CURVE!
for t = 400: 1: 600
         
    %Method 1 of finding P_sat- Antoine's eqn.
%     p_sat = Antoine_eqn(t); %The simple case of calculating by Antoine's eqn.
%     p_sat = p_sat/1e5;
    
%      p_sat = find_Psat(t, 0.5);    %THE TOLERANCE VALUE CAN BE CHANGED FROM HERE!!
    
     [p_sat, vl, vg] = P_by_fsolve(t);
    
    %These following commented lines are used in case of using the find_Psat(t, 1) and Antoine's eqn line.
    
%     alpha_dome = find_alpha(t);
%     a_dome = psi*alpha_dome*(R_new_bar*Tc)^2/(Pc_bar);    
%     eqn = [p_sat, p_sat*b_new*(sigma+epsilon-1)-R_new_bar*t, b_new*(b_new*sigma*epsilon-(sigma+epsilon)*(b_new+ R_new_bar*t)), a_dome-sigma*epsilon*(b_new^2)*(p_sat*b_new+R_new_bar*t)];
%     v = roots(eqn);  %In Litres
%     vg = max(v);
%     vl = min(v(v>0));
            
%    Check if isreal() can be of any help!
%    display(t);
    
    if ~isreal(vl)      %Break from the for loop when vl (and similarly vg) come out to be imaginary
        break
    end
    
%     display(vl);
%     display(vg);
        
    V_vec1 = [V_vec1, vg];
    V_vec2 = [V_vec2, vl];
    P_sat_vec = [P_sat_vec, p_sat];
       
end

plot(V_vec1, P_sat_vec, 'r', 'Linewidth', 4);
plot(V_vec2, P_sat_vec, 'r', 'Linewidth', 4);

%Easily seen from the graph that v(c) = 1.267 L, P(c) = 23.87 bars.
plot(1.368, 20.31, 'o', 'Markersize', 10, 'MarkerFaceColor', 'k'); 


%PLOTTING THE ISOTHERMS! (Func. definition at the end)
for T = 450: 30: 700
    plot_isotherm(T, 'blue');
end


% This procedure is the same as my function, just that I wanted to
% emphasize the critical curve more by it's linewidth and extra volume range. Hence, used the
% function here, with a very slight modification.

T = 524; %Observed from the graph

a = find_a(T);
hold on;

V_vec = [];
P_vec = [];

% for v = 0.4: 0.05: 6.5
for v = 0.04: 0.01: 6.5
    P = R_new_bar*T/(v - b_new) - a/((v+sigma*b_new)*(v+epsilon*b_new));
    P_vec = [P_vec, P];
    V_vec = [V_vec, v];
end

plot(V_vec, P_vec, 'color', 'k', 'Linewidth', 2);
hold off

%--------------------------------------------------------------------------
% All Functions below

function Treduced = Tred(T)  %To simply calculate the reduced temperature from T, Tc
    global Tc;
    Treduced = T/Tc;
end

function alpha = find_alpha(T)
    Tr = Tred(T);
    global w;
    alpha = (1 + (0.37464 + 1.54226*w - 0.26992*(w^2))*(1 - sqrt(Tr)))^2 ;
end

function a = find_a(T)
    alpha_Tr = find_alpha(T);
    global Tc Pc_bar R_new_bar psi;
    a = psi*alpha_Tr*(R_new_bar^2)*(Tc^2)/(Pc_bar);
end

function P = Antoine_eqn(T)
%The following values correspond to P in bars, and T in K.
    if (T < 308)   %Lower limit is 273K, but don't need to put that. 
%We won't go that low, and even if we do, the Antoine eqns. will still be approximately valid
        A = 4.14157;
        B = 1377.578;
        C = -50.507;
    elseif (T >= 308 &&  T < 384)
        A = 4.07827;
        B = 1343.943;
        C = -53.773;
    elseif (T >= 384)  % We use these parameters to get the guess value of P_sat for T > Tc too.
        A = 4.54436;
        B = 1738.123;
        C = 0.394; % A big change from the -ve values...
    end
    
    P = 10^(A - B/(T + C)); % Gives P in bars
    P = P*1e5; % Gives P in Pa
 
end

function mu = find_mu(T, v)
    global epsilon sigma b_new R_new_bar;
    a = find_a(T);
    
    Z = v/(v-b_new) - a*v/(R_new_bar*T*(v + epsilon*b_new)*(v + sigma*b_new)) ;
    %If P is given, then Z = Pv/RT.
    q = a/(b_new*R_new_bar*T);
    mu = -v*log(1 - b_new/v) - (q/sigma-epsilon)*log((v+sigma*b_new)/(v+epsilon*b_new)) + (Z-1) - log(Z);
    %This actually calculates mu wrt the ideal gas at same conditions of
    %T, P (i.e. by the departure function).
    mu = mu*R_new_bar*T;

end

% Calculating P_sat by equating chemical potentials of liq and vap phases
% as taught to us in DH and in Assg 3, Q.4 : abs(mu_l - mu_g) < tolerance.
% Works perfectly well, is just a little slow to compute owing to the large
% amounts of computations.
function P_sat = find_Psat(T, tolerance)
    
    global R_new_bar b_new sigma epsilon
    P = Antoine_eqn(T);  %This is our initial guess.
    P = P/1e5;
    P_ini = P;
    a = find_a(T);
    eqn = [P, P*b_new*(sigma+epsilon-1)-R_new_bar*T, b_new*(b_new*sigma*epsilon-(sigma+epsilon)*(b_new+ R_new_bar*T)), a-sigma*epsilon*(b_new^2)*(P*b_new+R_new_bar*T)];
    v = roots(eqn);
    v_l = min(v(v>0));
    v_g = max(v);
    mu_l = find_mu(T, v_l);
    mu_g = find_mu(T, v_g);
    diff = abs(mu_g - mu_l);
    
    while(diff > tolerance)  
        
        P1 = P_ini + 0.05; %Computing at a slightly higher pressure. 0.05 can be changed for convergence
        eqn = [P1, P1*b_new*(sigma+epsilon-1)-R_new_bar*T, b_new*(b_new*sigma*epsilon-(sigma+epsilon)*(b_new+ R_new_bar*T)), a-sigma*epsilon*(b_new^2)*(P1*b_new+R_new_bar*T)];
        v = roots(eqn);
        v_l1 = min(v(v>0));
        v_g1 = max(v);
        mu_l1 = find_mu(T, v_l1);
        mu_g1 = find_mu(T, v_g1);
        diff1 = abs(mu_g1 - mu_l1);
        
        P2 = P_ini - 0.05;  %Computing at a slightly lower pressure. 0.05 can be changed for convergence
        eqn = [P2, P2*b_new*(sigma+epsilon-1)-R_new_bar*T, b_new*(b_new*sigma*epsilon-(sigma+epsilon)*(b_new+ R_new_bar*T)), a-sigma*epsilon*(b_new^2)*(P2*b_new+R_new_bar*T)];
        v = roots(eqn);
        v_l2 = min(v(v>0));
        v_g2 = max(v);
        mu_l2 = find_mu(T, v_l2);
        mu_g2 = find_mu(T, v_g2);
        diff2 = abs(mu_g2 - mu_l2);
        
        if(diff1 < diff2)
            diff = diff1;
            P_ini = P1;
        else
            diff = diff2;
            P_ini = P2;
        end
             
    end
    P_sat = P_ini; %This P_ini has actually changed to the correct final P_sat.
end

function [P, vl, vg] = P_by_fsolve(T)
    global b_new R_new_bar sigma epsilon
    
    a = find_a(T);
    q = a/(b_new*R_new_bar*T);
    P_guess = Antoine_eqn(T);
    P_guess = P_guess/1e5;
        
    eqn = [P_guess, P_guess*b_new*(sigma+epsilon-1)-R_new_bar*T, b_new*(b_new*sigma*epsilon-(sigma+epsilon)*(b_new+ R_new_bar*T)), a-sigma*epsilon*(b_new^2)*(P_guess*b_new+R_new_bar*T)];
    v = roots(eqn);
    vl = min(v(v>0));
    vg = max(v);
        
        %F = @(x) [log((x(2)/x(3))*(x(3)- b_new)/(x(2)-b_new)) - x(1)*(x(3)-x(2))/(R_new_bar*T) + (q/sigma-epsilon)*log(((x(3)+sigma*b_new)/(x(2)+sigma*b_new))*((x(2)+epsilon*b_new)/(x(3)+epsilon*b_new))) + log(x(3)/x(2)); (x(2)^3)*x(1) + (x(2)^2)*(x(1)*b_new*(sigma+epsilon-1) - R_new_bar*T) + x(2)*(b_new*(b_new*sigma*epsilon - (sigma+epsilon)*(b_new + R_new_bar*T))) + a-sigma*epsilon*(b_new^2)*(b_new*x(1) + R_new_bar*T); (x(3)^3)*x(1) + (x(3)^2)*(x(1)*b_new*(sigma+epsilon-1) - R_new_bar*T) + x(3)*(b_new*(b_new*sigma*epsilon - (sigma+epsilon)*(b_new + R_new_bar*T))) + a-sigma*epsilon*(b_new^2)*(b_new*x(1) + R_new_bar*T)];
    F = @(x) [log(abs((x(2)/x(3))*((x(3)- b_new)/(x(2)-b_new)))) - ((R_new_bar*T/(x(2)-b_new)) - (a/((x(2)+sigma*b_new)*(x(2)+epsilon*b_new))))*(x(3)-x(2))/(R_new_bar*T) + (q/sigma-epsilon)*log(abs(((x(3)+sigma*b_new)/(x(2)+sigma*b_new))*((x(2)+epsilon*b_new)/(x(3)+epsilon*b_new)))) + log(abs(x(3)/x(2)));
        (x(2)^3)*x(1) + (x(2)^2)*(x(1)*b_new*(sigma+epsilon-1) - R_new_bar*T) + x(2)*(b_new*(b_new*sigma*epsilon - (sigma+epsilon)*(b_new + R_new_bar*T))) + a-sigma*epsilon*(b_new^2)*(b_new*x(1) + R_new_bar*T);
        (x(3)^3)*x(1) + (x(3)^2)*(x(1)*b_new*(sigma+epsilon-1) - R_new_bar*T) + x(3)*(b_new*(b_new*sigma*epsilon - (sigma+epsilon)*(b_new + R_new_bar*T))) + a-sigma*epsilon*(b_new^2)*(b_new*x(1) + R_new_bar*T)];
        
    ini_guesses = [P_guess, vl, vg];
    [x_ans, func_val_at_roots] = fsolve(F, ini_guesses);
    
    
    %IMP!!!!!!!!
    % x, x_ans is a vector with elements: P_sat, vl, vg.
     
    P = x_ans(1);
    vl = x_ans(2);
    vg = x_ans(3);
end

function plot_isotherm(T, col)
    
    global R_new_bar b_new sigma epsilon
    a = find_a(T);
    hold on;
    
    V_vec = [];
    P_vec = [];
        for v = 0.05: 0.01: 6    
        P = R_new_bar*T/(v - b_new) - a/((v+sigma*b_new)*(v+epsilon*b_new));
        P_vec = [P_vec, P];
        V_vec = [V_vec, v];
    end
    
    plot(V_vec, P_vec, 'color', col);
end

% For our system, by the code, vc obtained = 1.368 L/mol, Pc = 20.31 bar, Tc = 524 K.