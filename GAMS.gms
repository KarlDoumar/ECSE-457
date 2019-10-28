Sets
            m                       Index of energy storage use cases 1 to M    /1, 2, 3/
            s                       Index of energy storage subsystems 1 to S
            t                       Index of time periods 0 to T
            t0(m)                   Set of time periods during which use cas m is deployed            
            s_prime                 Index of energy storage subsystem that charges  another subsystem or who is being charge by another subsystem
    
Parameters
            alpha(t)                Price of energy in time period t ($ per kWh)
            C_eff(s)                Chargin efficiency of storage subsystem s (p.u.)
            D_eff(s)                Discharging efficiency of storage subsystem s (p.u.)
            delta                   Time step duration (h)
            lambda(t,m)             Value of energy storage use case m in time period t ($ per kW)
            e_0(s)                  Inital state of charge of storage subsystem s (kWh)
            g_a(t)                  Actual renewable energy generation in time period t
            g_t(t)                  Target renewable energy generation in time period t
            beta(s)                 Cost of power capacity for storage subsystem technology s ($ per kWh)
            gamma(s)                Cost of energy storage capacity for storage subsystem technology s ($ per kW)
            rho(s)                  Minimum energy-to-power ratio of energy storage subsystem s (kWh per kW)
            epsylon(s)              Coefficient relating the maximum and minimum state of charge of storage subsystem s (p.u.)
            Fvv(s)                  Coefficient relating the maximum and minimum power of storage subsystem s (p.u.)
            theta_max               Investment budget ($)

**Something about options??? TO CHECK OUT WHAT THAT IS

**If we decide to do two different files must do : $Inlcude: "NameOfTheFile";

Variables
            p(s,t,m)                Net power of storage subsystem s during time period t allocated to storage use casem m (kW)
            p_c(s,s_prime,t)        Charging power of storage subsystem s from storage subsystem s_prime during time period t (kW)
            p_d(s,s_prime,t)        Discharging power of storage subsystem s from storage subsystem s_prime during time period t (kW)
            e(s,t)                  State of charge of storage subsystem s in period t (kWh)
            e_max(s)                Maximum state of charge of storage subsystem s (kWh)
            e_min(s)                Minimum state of charge of storage subsystem s (kWh)
            p_max(s)                Maximum power of storage subsystem s (kWh)
            p_min(s)                Minimum power of storage subsystem s (kWh)
            u(s,t)                  Indicator for subsystem is charging or discharging 
            
            F_obj                   Variable for the objective function F
            G_obj                   Variable for the objective function G
            B(s,t,m)                Benefits for each use case

*Variable types
positive variables p, p_c, p_d, e, e_max, e_min, p_max;
negative variables p_min;
binary variables u;


Equations
*Constraints
            e_balance_lo(s,t)       Lower bound of the energy balance equation
            e_balance_up(s,t)       Upper bound of the energy balance
            p_balance_lo(s,t)     Lower bound of the power contribution
            p_balance_up(s,t)     Upper bound of the power contribution
            e_min_lo(s)             Lower bound of the minimum state of charge
            e_min_up(s)             Upper bound of the minimum state of charge
            e_max_lo(s)             Lower bound of the maximum state of charge
            p_max_lo(s)             Lower bound of the maximum power
            p_min_up(s)             Upper bound of the minimum power
            p_d_lo(s,s_prime,t)     Lower bound of the discharging power of storage subsystem s
            p_d_up(s,s_prime,t)     Upper bound of the discharging power of storage subsystem s
            p_c_lo(s,s_prime,t)     Lower bound of the charging power of storage subsystem s
            p_c_up(s,s_prime,t)     Upper bound of the charging power of storage subsystem s
            beta_lo(s)              Lower bound of the parameter beta(s)
            gamma_lo(s)             Lower bound of the parameter gamma(s)
            budget_lo               Lower bound of the budget
            epsylon_lo(s)           Lower bound of epsylon
            epsylon_up(s)           Upper bound of epsylon
            Fvv_lo(s)               Lower bound of Fvv
            Fvv_up(s)               Upper bound of Fvv
            lambda_lo(t,m)          Lower bound of lambda
            lambda_eq(t,m)          Upper bound of lambda
*Equations

            e_balance(s,t)          Energy balance for each storage subsystem
            p_balance(s,t)        Power contribution for each use case and storage subsystem
            e_min_max(s)            Relationship between min and max state of charge
            p_min_max(s)            Relationship between min and max power
            Use_case1(s,t,m)        Use cases m = 1
            Use_case2(s,t,m)        Use cases m = 2
            Use_case3(s,t,m)        Use cases m = 3
*            B(p(s,t,m))             Benefits for each use case
            Cost(s)                 Linear increasing costs
*            Cost(p_max(s),e_max(s))    Linear increasing costs            

*Objective functions
*To be maximized
            Objective_F
*To be minimized
            Objective_G
;

*Constraints
e_balance_lo(s,t)..             e_min(s) =l= e(s,t);

e_balance_up(s,t)..             e(s,t) =l= e_max(s);             

p_balance_lo(s,t)..             p_min(s) =l= sum(m, p(s,t,m));

p_balance_up(s,t)..             sum(m, p(s,t,m)) =l= p_max(s);

e_min_lo(s)..                   0 =l= e_min(s);

e_min_up(s)..                   e_min(s) < e_max(s);

e_max_lo(s)..                   rho(s)*p_max(s) =l= e_max(s);

p_max_lo(s)..                   10**(-16) =l= p_max(s);

*Check to see if we keep this, because alreadys tated that p_min will be a negative number!
p_min_up(s)..                   p_min(s) =l= 10**(-16);

p_d_lo(s,s_prime,t)..           0 =l= p_d(s,s_prime,t);                      

p_d_up(s,s_prime,t)..           p_d(s,s_prime,t) =l= p_min(s)*(u(s,t)-1);
  
p_c_lo(s,s_prime,t)..           0 =l= p_c(s,s_prime,t);

p_c_up(s,s_prime,t)..           p_c(s,s_prime,t) =l= p_max(s)*u(s,t);

beta_lo(s)..                    0 =l= beta(s);

gamma_lo(s)..                   0 =l= gamma(s);

epsylon_lo(s)..                 0 =l= epsylon(s);

epsylon_up(s)..                 epsylon(s) =l= 0.9999999999999999;

Fvv_lo(s)..                     10**(-16) =g= Fvv(s);

Fvv_up(s)..                     Fvv(s) =l= 1;

budget_lo..                     sum(s,Cost(s)) =l= theta_max;
*budget_lo..                     sum(s,C(p_max(s),e_max(s))) =l= theta_max;
               
*Equations

e_balance(s,t)..                e(s,t) =e= e(s,t-1) + delta*(C_eff(s)*sum(s_prime,p_c(s,s_prime,t))- sum(s_prime,p_d(s,s_prime,t)));
     
p_balance(s,t)..                sum(m,p(s,t,m)) =e= sum(s_prime,D_eff(s_prime)*p_d(s,s_prime,t)-p_c(s,s_prime,t));
    
e_min_max(s)..                  e_min(s) =e= epsylon(s)*e_max(s);                   

p_min_max(s)..                  p_min(s) =e= (-Fvv(s))*p_max(s);


*Objective functions

Objective_F..                   F_obj = sum(t,sum(s, sum(m, B(s,t,m))));
*Objective_F..                   F = sum(t,sum(s, sum(m, B(p(s,t,m)))));

Objective_G..                   G_obj = sum(s,Cost(s)- sum(t,sum(s, sum(m, B(s,t,m))));
*Objective_G..                   G = sum(s,C(p_max(s),e_max(s)))- sum(t,sum(s, sum(m, B(p(s,t,m)))));

*Different equations of B depending on the use case

*Peak Shaving (m = 1)

Use_case1(s,t,m)$ (ord(m) EQ 1)..       B(s,t,m) =e= lambda(t,'1')*p(s,t,1)+alpha(t)*p(s,t,1);    
*Use_cases(p(s,t,m))$ (ord(m) EQ 1)..    B(s,t,m) =e= lambda(t,1)*p(s,t,1)+alpha(t)*p(s,t,1);   

*Balancing (m = 2)

Use_case2(s,t,m)$ (ord(m) EQ 2)..       B(s,t,m) =e= (-lambda(t,'2'))*abs(g_a(t) + p(s,t,'2') - g_t(t)) + alpha(t)*p(s,t,'2');
*Use_cases(p(s,t,m))$ (ord(m) EQ 2)..    B(s,t,m) =e= (-lambda(t,2))*abs(g_a(t) + p(s,t,2) - g_t(t)) + alpha(t)*p(s,t,2);

*Price arbitrage (m = 3)

Use_case3(s,t,m)$ (ord(m) EQ 3)..       B(s,t,m) =e= alpha(t)*p(s,t,'3');
*Use_cases(p(s,t,m))$ (ord(m) EQ 3)..    B(s,t,m) =e= alpha(t)*p(s,t,3);

lambda_lo(t,m)$(ord(m) EQ 1)..          lambda(t,m) =g= 10**(-16);

lambda_lo(t,m)$(ord(m) EQ 2)..          lambda(t,m) =g= 0;

lambda_eq(t,m)$(ord(m) EQ 3)..          lambda(t,m) =e= 3;

*Will change the all later!
model benefits /all/;

model costs /all/;

solve benefits using lp maximizing F_obj;

solve costs using lp minimizing G_obj;
           



