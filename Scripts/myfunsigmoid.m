function F = myfunsigmoid(a,xdata)
%a = x(1);
%b = x(2);
%t0 = x(3);
%c = x(4);
%F = a + b*atan(xdata - t0) + c*xdata;

%a = x(1);
%b = x(2);
%F = a ./ (1 + exp(-b*xdata) );

F = tanh(a(1)*xdata).*(xdata>=0) + tanh(a(2)*xdata).*(xdata<0);

end