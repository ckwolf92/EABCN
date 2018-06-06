% specify function

test_fun = @(x) x^2;

% specify evaluation point

x0 = 10;
x0 = myAD(x0);

% play around with myAD

getvalues(x0)
getderivs(x0)

getvalues(test_fun(x0))
getderivs(test_fun(x0))