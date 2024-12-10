function add_pvalue(p_values, x1, x2, text_heights, line_heights)
for i = 1:numel(p_values)
    if p_values(i) < 0.001
        stars = '***';
    elseif p_values(i) < 0.01
        stars = '**';
    elseif p_values(i) < 0.05
        stars = '*';
    else
        stars = 'n.s';
    end
    line([x1(i), x2(i)], [line_heights(i),line_heights(i)], 'Color', 'black', 'LineWidth', 1);
    text((x1(i) + x2(i)) / 2, line_heights(i) + text_heights, stars, 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 12);
end
end

