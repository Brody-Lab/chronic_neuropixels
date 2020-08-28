P=get_parameters;
figure('pos', [100,100, 1000, 600])
set(gca, P.axes_properties{:}, ...
         'xtick', [], ...
         'ytick', [])

yoffset = 0.1;
y = 1;

str = '$N = Poisson(\lambda)$';
text(0.1,y, str, 'interpreter', 'latex', 'fontsize', P.font_size)

str = '$\lambda=Ae^{k_{fast}(t-1)}+Be^{k_{slow}(t-1)}$';
y = y - yoffset;
text(0.1,y, str, 'interpreter', 'latex', 'fontsize', P.font_size)

str = ['$A = \beta^A_0 + \beta^A_{AP} \cdot AP + \beta^A_{DV>-2} \cdot I_{DV>-2} + \beta^A_{ML} \cdot ML + \beta^A_{age} \cdot age + \beta^A_{use} \cdot use$'];
y = y - yoffset;
text(0.1,y, str, 'interpreter', 'latex', 'fontsize', P.font_size)

str = ['$B = \beta^B_0 + \beta^B_{AP} \cdot AP + \beta^B_{DV>-2} \cdot I_{DV>-2} + \beta^B_{ML} \cdot ML + \beta^B_{age} \cdot age + \beta^B_{use} \cdot use$'];
y = y - yoffset;
text(0.1,y, str, 'interpreter', 'latex', 'fontsize', P.font_size)

str = ['$k_{fast} = \beta^{k}_{fast} + k$'];
y = y - yoffset;
text(0.1,y, str, 'interpreter', 'latex', 'fontsize', P.font_size)

str = [' $k_{slow} = \beta^{k}_{slow} + k$'];
y = y - yoffset;
text(0.1,y, str, 'interpreter', 'latex', 'fontsize', P.font_size)

str = ['$k = \beta^{k}_{AP} \cdot AP +\beta^k_{DV>0} \cdot I_{DV>0} + \beta^{k}_{ML} \cdot ML +  \beta^{k}_{SP} \cdot SP +  \beta^{k}_{SO} \cdot SO + \beta^{k}_{age} \cdot age + \beta^{k}_{use} \cdot use$'];
y = y - yoffset;
text(0.1,y, str, 'interpreter', 'latex', 'fontsize', P.font_size)

for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'elaborated_mdl_eqs'], P.figure_image_format{i})
end