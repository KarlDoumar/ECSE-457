Sets
            m                       Index of energy storage use cases 1 to M    /1/
            s                       Index of energy storage subsystems 1 to n to S /1*3/
            Li(s)                   Dynamic set of lithium-ion batteries /2*3/
            Fly(s)                  Dynamic set of flywheels /1/
            t                       Index of time periods 0 to T /1*288/
*            t0(m)                   Set of time periods during which use cas m is deployed
            Alias(s,s_prime);
    
Parameters
*alpha(t) has been defined in TOU
            alpha(t)                Price of energy in time period t ($ per kWh)
            C_eff(s)                Chargin efficiency of storage subsystem s (p.u.) /set.Li 0.99, set.Fly 0.93/
            D_eff(s)                Discharging efficiency of storage subsystem s (p.u.) /set.Li 0.85, set.Fly 0.97/
            D_eff_prime(s_prime)    Discharging efficiency of storage subsystem s_prime (p.u.) /set.Li 0.85, set.Fly 0.97/
            delta                   Time step duration (h) /0.0833333/
************ Do three trials one where assumption is batteries start at full charge, another at no charge amd finally at mid charge?
            e_0(s)                  Inital state of charge of storage subsystem s (kWh)/set.Li 100, set.Fly 25/ 
            g_a(t)                  Actual renewable energy generation in time period t (kW)
*Flat out energy (by doing an integral)--> solar and/or wind
            g_t(t)                  Target renewable energy generation in time period t (kW)
************ Energy and power capital costs of 2018 used if we switch to 2025 it is estimated to be for li-ion:30.33 and 121.32 respectively, see which values to use! 
            beta(s)                 Cost of power capacity for storage subsystem technology s ($ per kW) /set.Li 173.96, set.Fly 283.80/
            gam(s)                  Cost of energy storage capacity for storage subsystem technology s ($ per kWh) /set.Li 43.49, set.Fly 2936.04/ 
            rho(s)                  Minimum energy-to-power ratio of energy storage subsystem s (kWh per kW) /set.Li 1, set.Fly 4/  
            epsylon(s)              Coefficient relating the maximum and minimum state of charge of storage subsystem s (p.u.) /set.Li 0.15, set.Fly 0.2/
            Fvv(s)                  Coefficient relating the maximum and minimum power of storage subsystem s (p.u.)/set.Li 0.8, set.Fly 0.5/
*            theta_max               Investment budget ($) /15000000/
            e_max(s)                Max energy per subsystem
            /1 25
            2*3 100/
            e_min(s)
            /1  5
            2*3 15
            /
            p_max(s)
            /1 6.25
            2*3 100
            /
            p_min(s)
            /1  -5
            2*3 -80
            /

            lambda(t,m)             Value of the use case m at time t
*                /1*288.1   0.00
/                1*5.1    0.5
                 6*122.1  0.06
                 123*194.1 0.01
                 195*218.1 0.06
                 219*288.1 0.03
*/                1*72.1    0.09
*                 73*122.1  0.18
*                 123*194.1 0.03
*                 195*218.1 0.18
*                 219*288.1 0.09
                 
                /
            
;

*$include "TOU.gms";
$include "TOU_day_winter.gms";
*$include "TOU_day_summer.gms";
$include "Solar_Day_winter.gms"
*$include "Solar_Day_summer.gms"
$include "Target_Solar_Day_winter.gms"
*$include "Target_Solar_Day_summer.gms"

**Something about options??? TO CHECK OUT WHAT THAT IS

Variables
            p(s,t,m)                Net power of storage subsystem s during time period t allocated to storage use case m (kW)
            p_c(s,s_prime,t)        Charging power of storage subsystem s from storage subsystem s_prime during time period t (kW)
            p_d(s,s_prime,t)        Discharging power of storage subsystem s from storage subsystem s_prime during time period t (kW)
            e(s,t)                  State of charge of storage subsystem s in period t (kWh)
*            e_max(s)                Maximum state of charge of storage subsystem s (kWh)
*            e_min(s)                Minimum state of charge of storage subsystem s (kWh)
*            p_max(s)                Maximum power of storage subsystem s (kWh)
*            p_min(s)                Minimum power of storage subsystem s (kWh)
            
            F_obj                   Variable for the objective F
            G_obj                   Variable for the objective G
            B(s,t,m)                Benefits for each use case
            Cost(s)                 Linear increasing costs
            y(s,t,m)                For linear purposes
            u(s,t)                  1 is for charging and 0 is for discharging

*Variable types
positive variables p, p_c, p_d, e;



*Initial conditions 

e.fx('1','1') = e_0('1');
e.fx('2','1') = e_0('2');
e.fx('3','1') = e_0('3');
*e.fx('4','1') = e_0('4');
*e.fx('5','1') = e_0('5');
*e.fx('6','1') = e_0('6');

*e.fx('7','1') = e_0('7');
*e.fx('8','1') = e_0('8');
*e.fx('9','1') = e_0('9');
*e.fx('10','1') = e_0('10');

**********************************************Check and see if this is okay *******************************************
p.fx('1','1','1') = 0.1*(1/rho('1')*e_0('1'));
p.fx('2','1','1') = 0.4*(1/rho('2')*e_0('2'));
p.fx('3','1','1') = 0.1*(1/rho('3')*e_0('3'));
*p.fx('4','1','1') = 0*(1/rho('1')*e_0('4'));
*p.fx('5','1','1') = 0.6*(1/rho('2')*e_0('5'));
*p.fx('6','1','1') = 0.1*(1/rho('3')*e_0('6'));
*p.fx('1','1','2') = 0.1*(1/rho('1')*e_0('1'));
*p.fx('2','1','2') = 0.4*(1/rho('2')*e_0('2'));
*p.fx('3','1','2') = 0.1*(1/rho('3')*e_0('3'));
*p.fx('4','1','2') = 0*(1/rho('1')*e_0('4'));
*p.fx('5','1','2') = 0.6*(1/rho('2')*e_0('5'));
*p.fx('6','1','2') = 0.1*(1/rho('3')*e_0('6'));
*p.fx('7','1','1') = 0.1*(1/rho('1')*e_0('7'));
*p.fx('8','1','1') = 0.1*(1/rho('2')*e_0('8'));
*p.fx('9','1','1') = 0.1*(1/rho('3')*e_0('9'));
*p.fx('10','1','1') = 0.1*(1/rho('3')*e_0('10'));
*p.fx('1','1','3') = 0.1*(1/rho('1')*e_0('1'));
*p.fx('2','1','3') = 0.4*(1/rho('2')*e_0('2'));
*p.fx('3','1','3') = 0.1*(1/rho('3')*e_0('3'));

p_c.fx(s,s_prime,'1') = 0;

p_d.fx('1',s_prime,'1') = 2;
p_d.fx('2',s_prime,'1') = 5;
p_d.fx('3',s_prime,'1') = 70;
*p_d.fx('4',s_prime,'1') = 70;
*p_d.fx('5',s_prime,'1') = 20;
*p_d.fx('6',s_prime,'1') = 60;
*p_d.fx('7',s_prime,'1') = 80;
*p_d.fx('8',s_prime,'1') = 2;
*p_d.fx('9',s_prime,'1') = 70;
*p_d.fx('10',s_prime,'1') = 70;

*p.fx('4','1',m) = 1/rho('4')*e_0('4');
*p.fx('5','1',m) = 1/rho('5')*e_0('5');
*p.fx('6','1',m) = 1/rho('6')*e_0('6');
*p.fx('7','1',m) = 1/rho('7')*e_0('7');
*p.fx('8','1',m) = 1/rho('8')*e_0('8');
*p.fx('9','1',m) = 1/rho('9')*e_0('9');
*p.fx('10','1',m) =1/rho('10')* e_0('10');

u.fx('1','1') = 0; u.fx('2','1') = 0; u.fx('3','1') = 0
************************************************************************************************************************



Equations
*Constraints
            e_balance_lo(s,t)       Lower bound of the energy balance equation
            e_balance_up(s,t)       Upper bound of the energy balance
            p_balance_lo(s,t)       Lower bound of the power contribution
            p_balance_up(s,t)       Upper bound of the power contribution
            e_min_lo(s)             Lower bound of the minimum state of charge
            e_min_up(s)             Upper bound of the minimum state of charge
            e_max_lo(s)             Lower bound of the maximum state of charge
            p_max_lo(s)             Lower bound of the maximum power
            p_min_up(s)             Upper bound of the minimum power
            p_d_lo(s,s_prime,t)     Lower bound of the discharging power of storage subsystem s
            p_d_up(s,s_prime,t)     Upper bound of the discharging power of storage subsystem s
            p_c_lo(s,s_prime,t)     Lower bound of the charging power of storage subsystem s
            p_c_up(s,s_prime,t)     Upper bound of the charging power of storage subsystem s
            abs_const(s,t,m)        Constraint for absolute function
            abs_pos(s,t,m)          Function for absolute function
            abs_neg(s,t,m)          Function for absolute function
            p_min_lo(s)
            
*Equations

            e_balance_2(s,t)        Energy balance for each storage subsystem
            p_balance(s,t)          Power contribution for each use case and storage subsystem
*            e_min_max(s)            Relationship between min and max state of charge
            Use_case1(s,t,m)        Use cases m is 1
*            Use_case2(s,t,m)        Use cases m is 2 (absolute value)
*            Use_case3(s,t,m)        Use cases m is 3
            Budget(s)               Cost of the whole operation
         

*Objective functions
*To be maximized
            Objective_F
*To be minimized
            Objective_G
;

*Constraints
e_balance_lo(s,t)..             e(s,t) =g= e_min(s);

e_balance_up(s,t)..             e(s,t) =l= e_max(s);

p_balance_lo(s,t)..             sum(m, p(s,t,m)) =g= p_min(s);

p_balance_up(s,t)..             sum(m, p(s,t,m)) =l= p_max(s);
        
p_d_lo(s,s_prime,t)..           p_d(s,s_prime,t) =g= 0;
 
p_d_up(s,s_prime,t)..           p_d(s,s_prime,t) =l= p_min(s)*(-1);
                  
*p_d_up(s,s_prime,t)..           p_d(s,s_prime,t) =l= p_min(s)*(u(s,t)-1);
  
p_c_lo(s,s_prime,t)..           p_c(s,s_prime,t) =g= 0;

p_c_up(s,s_prime,t)..           p_c(s,s_prime,t) =l= p_max(s);

*p_c_up(s,s_prime,t)..           p_c(s,s_prime,t) =l= p_max(s)*u(s,t);

***recharge the battery at the end of the day                         
*Equations

e_balance_2(s,t)$(ord(t) NE 1)..                e(s,t) =e= e(s,t-1) + delta*(C_eff(s)*sum(s_prime,p_c(s,s_prime,t))- sum(s_prime,p_d(s,s_prime,t)));
     
p_balance(s,t)$(ord(t) NE 1)..                  sum(m,p(s,t,m)) =e= sum(s_prime,D_eff_prime(s_prime)*p_d(s,s_prime,t)-p_c(s,s_prime,t));

*e_min_max(s)..                    e_min(s) =e= epsylon(s)*e_max(s);                   

Budget(s)..                       Cost(s) =e= beta(s)*p_max(s)+ gam(s)*e_max(s);

*Objective functions

Objective_F..                     F_obj =e= sum((t,s,m),B(s,t,m));


Objective_G..                     G_obj =e= sum(s,Cost(s))- sum((t,s,m), B(s,t,m));

*Different equations of B depending on the use case

*Peak Shaving (m = 1)

Use_case1(s,t,m)$ (ord(m) EQ 1)..          B(s,t,'1') =e= lambda(t,'1')*p(s,t,'1')+alpha(t)*p(s,t,'1');    


*Balancing (m = 2)

            

*Use_case3(s,t,m)$(ord(m) EQ 3)..           B(s,t,'3') =g= (-lambda(t,'3'))*y(s,t,'3') + alpha(t)*p(s,t,'3');

*abs_pos(s,t,m)$(ord(m) EQ 3)..             y(s,t,'3') =g= (g_a(t) + p(s,t,'1') - g_t(t));

*abs_neg(s,t,m)$(ord(m) EQ 3)..             y(s,t,'3') =g= -1*(g_a(t) + p(s,t,'1') - g_t(t));

*abs_const(s,t,m)$(ord(m) EQ 3)..           y(s,t,'3') =g= 0;



*Price arbitrage (m = 3)

*Use_case2(s,t,m)$(ord(m) EQ 2)..          B(s,t,'2') =e= alpha(t)*p(s,t,'2');

*Will change the all later!
*model benefits /Objective_F, e_balance_lo, e_balance_up, abs_pos, abs_neg, abs_const,p_balance_lo, p_balance_up, p_d_lo, p_d_up, p_c_lo, p_c_up, e_balance_2, p_balance, Use_case1, Use_case2,Use_case3/;

*model costs /Objective_G, e_balance_lo, e_balance_up, p_balance_lo, p_balance_up, p_d_lo, p_d_up, p_c_lo, p_c_up, e_balance_2, p_balance, Use_case1, Use_case2,Use_case3, Budget/;

*model benefits /Objective_F, e_balance_lo, e_balance_up, abs_pos, abs_neg, abs_const,p_balance_lo, p_balance_up, p_d_lo, p_d_up, p_c_lo, p_c_up, e_balance_2, p_balance, Use_case1, Use_case2,Use_case3/;

*model costs /Objective_G, e_balance_lo, e_balance_up, p_balance_lo, p_balance_up, p_d_lo, p_d_up, p_c_lo, p_c_up, e_balance_2, p_balance, Use_case1, Use_case2,Use_case3, Budget/;

model benefits /Objective_F, e_balance_lo, e_balance_up, p_balance_lo, p_balance_up, p_d_lo, p_d_up, p_c_lo, p_c_up, e_balance_2, p_balance, Use_case1/;

model costs /Objective_G, e_balance_lo, e_balance_up, p_balance_lo, p_balance_up, p_d_lo, p_d_up, p_c_lo, p_c_up, e_balance_2, p_balance, Use_case1 Budget/;
*Export results to gds, or export them into mathlab
*Do a grid for the index of each subsystem and then change place in the grad for each iteration
*for(count = 3 to (card(s)+1)
    
   
solve benefits maximizing F_obj using mip;

solve costs minimizing G_obj using mip;
        


 


