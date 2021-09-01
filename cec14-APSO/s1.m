function state = s1(f)

if f <= 0.4
    state = 0;
elseif f <= 0.6
    state = 5*f - 2;
elseif f <= 0.7
    state = 1;
elseif f <= 0.8
    state = -10*f + 8;
else
    state = 0;
end

end