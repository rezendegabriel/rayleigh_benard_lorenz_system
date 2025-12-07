using Plots

# ==========================
# Sistema de Lorenz
# ==========================
@kwdef mutable struct Lorenz
    dt::Float64 = 0.01
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
    dx = l.σ * (y - x)
    dy = x*(l.ρ - z) - y
    dz = x*y - l.β*z
    return dx, dy, dz
end

# ==========================
# RK4
# ==========================

function step_RK4!(l::Lorenz)
    dt = l.dt

    k1x,k1y,k1z = f(l, l.x, l.y, l.z)

    k2x,k2y,k2z = f(l,
        l.x + 0.5dt*k1x,
        l.y + 0.5dt*k1y,
        l.z + 0.5dt*k1z)

    k3x,k3y,k3z = f(l,
        l.x + 0.5dt*k2x,
        l.y + 0.5dt*k2y,
        l.z + 0.5dt*k2z)

    k4x,k4y,k4z = f(l,
        l.x + dt*k3x,
        l.y + dt*k3y,
        l.z + dt*k3z)

    l.x += dt/6 * (k1x + 2k2x + 2k3x + k4x)
    l.y += dt/6 * (k1y + 2k2y + 2k3y + k4y)
    l.z += dt/6 * (k1z + 2k2z + 2k3z + k4z)
end

# ==========================
# Diretório de saída
# ==========================
save_path = "C:/Users/3viei/OneDrive/Documentos/LNCC/Estrutura_de_dados/Implementacoes/"

# ==========================
# Condições iniciais
# ==========================
initial_conditions = [
    (2.0, 1.0, 1.0),
    (3.0, 2.0, 1.0),
    (1.0, 3.0, 2.0),
    (-2.0, -1.5, 4.0),
    (5.0, 2.5, 1.5),
    (-4.0, 3.0, 2.0)
]

# ==========================
# Simulação
# ==========================

N = 60000

for (i, (x0,y0,z0)) in enumerate(initial_conditions)

    traj = Lorenz(x=x0, y=y0, z=z0)

    xs = Float64[]
    ys = Float64[]
    zs = Float64[]

    # ========= integração =========
    for k ∈ 1:N
        step_RK4!(traj)
        push!(xs, traj.x)
        push!(ys, traj.y)
        push!(zs, traj.z)
    end

    # ========= FIGURA 3D =========
    p3d = plot3d(
        xs, ys, zs,
        linewidth = 1.8,
        xlim=(-30,30),
        ylim=(-30,30),
        zlim=(0,60),
        legend=false,
        title="Atrator de Lorenz – Trajetória $(i)"
    )

    savefig(p3d, save_path * "lorenz_3D_traj_$(i).png")

    # ========= PROJEÇÕES 2D =========

    # x × y
    p_xy = plot(
        xs, ys,
        linewidth=1.8,
        xlabel="x",
        ylabel="y",
        title="Projeção XY – Trajetória $(i)",
        legend=false
    )

    savefig(p_xy, save_path * "lorenz_XY_traj_$(i).png")

    # x × z
    p_xz = plot(
        xs, zs,
        linewidth=1.8,
        xlabel="x",
        ylabel="z",
        title="Projeção XZ – Trajetória $(i)",
        legend=false
    )

    savefig(p_xz, save_path * "lorenz_XZ_traj_$(i).png")

    # y × z
    p_yz = plot(
        ys, zs,
        linewidth=1.8,
        xlabel="y",
        ylabel="z",
        title="Projeção YZ – Trajetória $(i)",
        legend=false
    )

    savefig(p_yz, save_path * "lorenz_YZ_traj_$(i).png")

    println("Trajetória $(i) finalizada.")

end

println("✅ Todas as figuras foram geradas com sucesso!")
