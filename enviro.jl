using Pkg
function activate_env(path_to_repo, env_name)
    
    path = path_to_repo * "dynamical-systems\\env\\" * env_name
    Pkg.activate(path)

end
function activate_env(env_name)
    
    path = pwd() * "\\env\\"*env_name
    Pkg.activate(path)

end