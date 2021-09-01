function state = s4(f)

if f <= 0.7
    state = 0;
elseif f <= 0.9
    state = 5*f - 3.5;
else
    state = 1;
end

end