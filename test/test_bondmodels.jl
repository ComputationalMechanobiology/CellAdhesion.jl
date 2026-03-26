
println("===============================================")
println("Testing bondmodels.jl")
println("===============================================")


function _check_SlipBondModel()

    k_on_params = (k_on_0 = 0.2,)
    k_off_params = (k_off_0 = 0.8, f_1e = 1.0)

    model = SlipBondModel(k_on_params, k_off_params)

    typeof(model) == BondModel{Slip}

end

@test _check_SlipBondModel()



function _check_propertynames()

    k_on_params = (k_on_0 = 0.2,)
    k_off_params = (k_off_0s = 0.8, f_1es = 1.0, k_off_0c = 0.8, f_1ec = 1.0)

    model = CatchBondModel(k_on_params, k_off_params)


    Base.propertynames(model) == keys(getfield(model,:p))

end

@test _check_SlipBondModel()


