% PDOT: Computes the rate of change of actual power (PDOT) based on the power command (P1) 
%       and current actual power (P3).
% Purpose: Models the dynamic response of the engine by calculating how quickly the actual 
%          power (P3) should adjust toward the commanded power (P1). This introduces a 
%          first-order lag to simulate real engine behavior, preventing abrupt power changes.
% Why Used: It ensures a smooth transition of engine power in response to throttle inputs, 
%           which is critical for flight stability, especially in a fighter aircraft like 
%           the F-16 operating near stall/post-stall conditions. The threshold (50.0) 
%           differentiates low-power and high-power regimes with different time constants.
% Input: P3 - Actual power (in arbitrary units, e.g., % or kW)
%        P1 - Power command (from TGEAR, same units as P3)
% Output: PDOT - Rate of change of power (units: power/time, e.g., %/s or kW/s)




function PDOT = PDOT(P3,P1)

if P1>=50
    if P3>=50
       TF = 5.0; % TF = 1/TAU_T Inverse Time factor
       P2 = P1;
    else
        P2 = 60;
        TF = TAU(P2-P3);
    end
else
    if P3>=50
        TF = 5.0;
        P2 = 40;
    else
        P2 = P1;
        TF = TAU(P2-P3);
    end
end

PDOT = TF * (P2-P3);

end

function TAU = TAU(DP)

if DP <= 25
    TAU = 1.0;
elseif DP >= 50
    TAU = 0.1;
else
    TAU = 1.9 - 0.36 * DP;    
end

end
