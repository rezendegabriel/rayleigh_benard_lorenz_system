using Plots

# ==========================
# Estrutura do sistema
# ==========================
@kwdef mutable struct Lorenz
    dt::Float64 = 0.01
    σ::Float64 = 10.0
    ρ::Float64 = 28
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
    dy = x * (l.ρ - z) - y
    dz = x * y - l.β * z
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
        l.z + 0.5dt*k1z)

    k3x,k3y,k3z = f(l,
        l.x + 0.5dt*k2x,
        l.y + 0.5dt*k2y,
        l.z + 0.5dt*k2z)

    k4x,k4y,k4z = f(l,
        l.x + dt*k3x,
        l.y + dt*k3y,
        l.z + dt*k3z)

    l.x += dt/6*(k1x+2k2x+2k3x+k4x)
    l.y += dt/6*(k1y+2k2y+2k3y+k4y)
    l.z += dt/6*(k1z+2k2z+2k3z+k4z)
end

# ==========================
# Duas trajetórias
# ==========================
ϵ = 1e-15

traj1 = Lorenz()
traj2 = Lorenz(x = 2.0 + ϵ)   # perturbacao em x

# ==========================
# Armazenamento
# ==========================
N = 6000

t = collect(0:traj1.dt:(N-1)*traj1.dt)

x1 = zeros(N); y1 = zeros(N); z1 = zeros(N)
x2 = zeros(N); y2 = zeros(N); z2 = zeros(N)

# ==========================
# Simulação
# ==========================
for i in 1:N

    step_RK4!(traj1)
    step_RK4!(traj2)

    x1[i] = traj1.x
    y1[i] = traj1.y
    z1[i] = traj1.z

    x2[i] = traj2.x
    y2[i] = traj2.y
    z2[i] = traj2.z

end

# ==========================
# Diretório de salvamento
# ==========================
save_path = "C:/Users/3viei/OneDrive/Documentos/LNCC/Estrutura_de_dados/Implementacoes/"

# ==========================
# Gráfico x(t)
# ==========================
px = plot(
    t, x1,
    label = "Sem perturbação",
    xlabel = "Tempo",
    ylabel = "x(t)",
    title = "Comparação de x(t)",
    linewidth = 2
)

plot!(
    px,
    t, x2,
    label = "Com perturbação",
    linewidth = 2
)

savefig(px, save_path * "xt_comparacao.png")

# ==========================
# Gráfico y(t)
# ==========================
py = plot(
    t, y1,
    label = "Sem perturbação",
    xlabel = "Tempo",
    ylabel = "y(t)",
    title = "Comparação de y(t)",
    linewidth = 2
)

plot!(
    py,
    t, y2,
    label = "Com perturbação",
    linewidth = 2
)

savefig(py, save_path * "yt_comparacao.png")

# ==========================
# Gráfico z(t)
# ==========================
pz = plot(
    t, z1,
    label = "Sem perturbação",
    xlabel = "Tempo",
    ylabel = "z(t)",
    title = "Comparação de z(t)",
    linewidth = 2
)

plot!(
    pz,
    t, z2,
    label = "Com perturbação",
    linewidth = 2
)

savefig(pz, save_path * "zt_comparacao.png")

println("✅ Gráficos xt, yt e zt gerados com sucesso!")
