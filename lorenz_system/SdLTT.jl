using Plots

# ==========================
# Sistema de Lorenz
# ==========================
@kwdef mutable struct Lorenz
    dt::Float64 = 0.02
    σ::Float64 = 10.0
    ρ::Float64 = 28.0
    β::Float64 = 8/3
    x::Float64 = 2.0
    y::Float64 = 1.0
    z::Float64 = 1.0
end

# ==========================
# Campo vetorial
# ==========================
function f(l::Lorenz, x, y, z)
    dx = l.σ*(y-x)
    dy = x*(l.ρ-z) - y
    dz = x*y - l.β*z
    return dx,dy,dz
end

# ==========================
# RK4
# ==========================
function step_RK4!(l::Lorenz)
    dt = l.dt

    k1x,k1y,k1z = f(l,l.x,l.y,l.z)
    k2x,k2y,k2z = f(l,
        l.x + 0.5dt*k1x,
        l.y + 0.5dt*k1y,
        l.z + 0.5dt*k1z
    )
    k3x,k3y,k3z = f(l,
        l.x + 0.5dt*k2x,
        l.y + 0.5dt*k2y,
        l.z + 0.5dt*k2z
    )
    k4x,k4y,k4z = f(l,
        l.x + dt*k3x,
        l.y + dt*k3y,
        l.z + dt*k3z
    )
    l.x += dt/6*(k1x+2k2x+2k3x+k4x)
    l.y += dt/6*(k1y+2k2y+2k3y+k4y)
    l.z += dt/6*(k1z+2k2z+2k3z+k4z)
end

# ==========================
# CONFIGURAÇÃO
# ==========================
ϵ = 1e-10
Ntotal = 3000
Nplot  = 1500

save_path = "C:/Users/3viei/OneDrive/Documentos/LNCC/Estrutura_de_dados/Implementacoes/"

# ==========================
# SIMULAÇÃO E PLOT
# ==========================
function simular_plotar(traj_base::Lorenz, traj_pert::Lorenz,
                         titulo::String, nome_arquivo::String)

    xb = Float64[]; yb = Float64[]; zb = Float64[]
    xp = Float64[]; yp = Float64[]; zp = Float64[]

    for i ∈ 1:Ntotal

        step_RK4!(traj_base)
        step_RK4!(traj_pert)

        push!(xb,traj_base.x); push!(yb,traj_base.y); push!(zb,traj_base.z)
        push!(xp,traj_pert.x); push!(yp,traj_pert.y); push!(zp,traj_pert.z)

    end

    # --------- manter apenas os últimos pontos ---------
    xb = xb[end-Nplot+1:end]
    yb = yb[end-Nplot+1:end]
    zb = zb[end-Nplot+1:end]

    xp = xp[end-Nplot+1:end]
    yp = yp[end-Nplot+1:end]
    zp = zp[end-Nplot+1:end]

    plt = plot3d(
        xb,yb,zb,
        label="Sem perturbação",
        linewidth=2,
        xlim=(-30,30), ylim=(-30,30), zlim=(0,60),
        title=titulo
    )

    plot!(
        plt,
        xp,yp,zp,
        label="Perturbado",
        linewidth=2
    )

    savefig(plt, save_path * nome_arquivo)

    println("✅ Figura salva: ", nome_arquivo)
end

# ==========================
# EXECUÇÃO DOS 3 CASOS
# ==========================

# Perturbação em X
simular_plotar(
    Lorenz(),
    Lorenz(x = 2.0 + ϵ),
    "Perturbação em x — Lorenz (RK4)",
    "lorenz_perturb_x.png"
)

# Perturbação em Y
simular_plotar(
    Lorenz(),
    Lorenz(y = 1.0 + ϵ),
    "Perturbação em y — Lorenz (RK4)",
    "lorenz_perturb_y.png"
)

# Perturbação em Z
simular_plotar(
    Lorenz(),
    Lorenz(z = 1.0 + ϵ),
    "Perturbação em z — Lorenz (RK4)",
    "lorenz_perturb_z.png"
)
