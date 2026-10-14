function x = project_box_local(x, lb, ub)
x = min(max(x, lb), ub);
end

