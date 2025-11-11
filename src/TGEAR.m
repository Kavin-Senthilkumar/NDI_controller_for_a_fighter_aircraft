% TGEAR: Computes the power command (TGEAR) based on the throttle setting (THTL).
% Purpose: Translates the pilot's throttle input (a normalized value between 0 and 1) 
%          into a power command that the engine control system can use. This function 
%          models the nonlinear relationship between throttle position and engine power 
%          demand, derived from empirical data for the F-16 engine.
% Why Used: It provides a realistic interface between the pilot's control input and 
%           the engine's power dynamics, ensuring accurate power commands (P1) for 
%           subsequent dynamic calculations. The piecewise linear mapping accounts 
%           for different engine operating regimes (e.g., low power vs. high power).
% Input: THTL - Throttle setting (normalized, 0 to 1)
% Output: TGEAR - Power command value (in arbitrary units, e.g., % or scaled power)



function TGEAR = TGEAR(THTL)
    
    if THTL <= 0.77
        TGEAR = 64.94*THTL;
    else
        TGEAR = 217.38*THTL - 117.38;
    end
end


