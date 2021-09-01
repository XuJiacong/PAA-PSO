function state = s3(f)

if f <= 0.1
    state = 1;
elseif f <= 0.3
    state = -5*f + 1.5;
else
    state = 0;
end

end