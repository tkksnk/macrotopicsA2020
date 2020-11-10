function util = CRRA(cons, gamma)

if gamma ~= 1
%     util = cons.^(1-gamma)./(1-gamma);
    util = (cons.^(1-gamma)-1)./(1-gamma);
else
    util = log(cons);
end