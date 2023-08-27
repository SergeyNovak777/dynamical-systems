function get_fixed_point(sys, jac_sys, params, dim, box)
    ds = CoupledODEs(sys, zeros(dim), params)
    fp, _, _ = fixedpoints(ds, box, jac_sys)
end
function get_matrix(fixedpoint, normalize)
    
end