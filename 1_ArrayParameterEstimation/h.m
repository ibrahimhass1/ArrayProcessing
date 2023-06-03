function y = h(t)
    if t >= 0 && t < 0.25
        y = 1;
    elseif t >= 0.25 && t < 0.5
        y = -1;
    elseif t >= 0.5 && t < 0.75
        y = 1;
    elseif t >= 0.75 && t < 1
        y = -1;
    else
        y = 0; % default value if t is outside the defined intervals
    end
end
