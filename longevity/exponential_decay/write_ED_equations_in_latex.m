P=get_parameters;
figure
set(gca, P.axes_properties{:})
text(0.1,1, '$N=N_{1}e^{-(t-1)/\tau}+\varepsilon$', 'interpreter', 'latex', 'fontsize', P.font_size)

str = ['$N_{1}=\beta_{0}^{N_{1}} + \beta_{AP}^{N_{1}} + \beta_{DV}^{N_{1}}' ...
       '+ \beta_{ML}^{N_{1}} + \beta_{SP}^{N_{1}} + \beta_{SO}^{N_{1}}$'];
text(0.1,0.8, str, 'interpreter', 'latex', 'fontsize', P.font_size)

str = ['$\tau=\beta_{0}^{\tau} + \beta_{AP}^{\tau} + \beta_{DV}^{\tau}' ...
       '+ \beta_{ML}^{\tau} + \beta_{SP}^{\tau} + \beta_{SO}^{\tau}$'];
text(0.1,0.6, str, 'interpreter', 'latex', 'fontsize', P.font_size)