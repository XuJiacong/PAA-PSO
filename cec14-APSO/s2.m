function state = s2(f)

if f <= 0.2
    state = 0;
elseif f <= 0.3
    state = 10*f - 2;
elseif f <= 0.4
    state = 1;
elseif f <= 0.6
    state = -5*f + 3;
else
    state = 0;
end

end