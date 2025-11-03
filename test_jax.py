import numpy as np
import os
import jax
import jax.numpy as jnp
from jax import lax
import time
from tqdm import tqdm
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

def read_conditions_fput_parallel(filename, num_condizioni, N):
    """
    Legge num_condizioni random dal file binario simultaneamente usando memmap.

    Args:
        filename: str, path al file binario
        num_condizioni: int, numero di condizioni da leggere
        N: int, numero di masse mobili

    Returns:
        x: jnp.array, shape (num_condizioni, N+2)
        p: jnp.array, shape (num_condizioni, N+2)
        xiL: jnp.array, shape (num_condizioni,)
        xiR: jnp.array, shape (num_condizioni,)
    """
    neq = 2*(N+2)+2
    cond_size_bytes = neq * 8  # double = 8 byte

    # dimensione totale file
    total_size = os.path.getsize(filename)
    total_cond = total_size // cond_size_bytes

    if num_condizioni > total_cond:
        raise ValueError(f"Requested {num_condizioni} conditions, but file has only {total_cond}")

    # seleziona indici casuali
    indices = np.arange(total_cond)
    #np.random.shuffle(indices)
    selected_indices = indices[:num_condizioni//2]

    # crea memmap (file letto solo quando necessario)
    mm = np.memmap(filename, dtype=np.float64, mode='r', shape=(total_cond, neq))

    # leggi le condizioni richieste in batch
    data = mm[selected_indices]

    # separa x, p, xiL, xiR
    x = jnp.array(data[:, 0:-2:2])
    p = jnp.array(data[:, 1:-2:2])
    xiL = jnp.array(data[:, -2])
    xiR = jnp.array(data[:, -1])

    x=jnp.concatenate([x,x],axis=0)
    p=jnp.concatenate([p,-p],axis=0)
    xiL=jnp.concatenate([xiL,-xiL],axis=0)
    xiR=jnp.concatenate([xiR,-xiR],axis=0)
    return x, p, xiL, xiR


# ------------------------------------------------------
# RHS DEL SISTEMA (FPUT + TERMOSTATI)
# ------------------------------------------------------
def rhs_fput(x, p, xi_L, xi_R):
    """Equazioni del moto per tutte le catene."""
    # Derivate
    dxdt = p / m
    dpdt = jnp.zeros_like(p)

    # differenze tra posizioni adiacenti
    r_right = x[:, 2:] - x[:, 1:-1] - a
    r_left  = x[:, 1:-1] - x[:, :-2] - a

    # forza FPUT interna
    F_int = chi * (x[:, 2:] + x[:, :-2] - 2*x[:, 1:-1]) \
            + alpha * (r_right**2 - r_left**2) \
            + beta * (r_right**3 - r_left**3)

    # assegna alle particelle interne
    dpdt = dpdt.at[:, 1:-1].set(F_int)

    # accoppiamento con termostati (Nosé–Hoover)
    dpdt = dpdt.at[:, 1].add(-xi_L * p[:, 1])
    dpdt = dpdt.at[:, -2].add(-xi_R * p[:, -2])

    # derivate dei termostati
    kin_L = (p[:, 1]**2) / m
    kin_R = (p[:, -2]**2) / m

    dxi_L = (kin_L / Tl - 1.0) / thetaL**2
    dxi_R = (kin_R / Tr - 1.0) / thetaR**2

    # estremi fissi
    dxdt = dxdt.at[:, 0].set(0.0)
    dxdt = dxdt.at[:, -1].set(0.0)
    dpdt = dpdt.at[:, 0].set(0.0)
    dpdt = dpdt.at[:, -1].set(0.0)

    return dxdt, dpdt, dxi_L, dxi_R


# ------------------------------------------------------
# PASSO RK4
# ------------------------------------------------------
def rk4_step(x, p, xi_L, xi_R, dt):
    def f(x, p, xi_L, xi_R):
        return rhs_fput(x, p, xi_L, xi_R)

    k1x, k1p, k1xiL, k1xiR = f(x, p, xi_L, xi_R)
    k2x, k2p, k2xiL, k2xiR = f(x + 0.5*dt*k1x,
                               p + 0.5*dt*k1p,
                               xi_L + 0.5*dt*k1xiL,
                               xi_R + 0.5*dt*k1xiR)
    k3x, k3p, k3xiL, k3xiR = f(x + 0.5*dt*k2x,
                               p + 0.5*dt*k2p,
                               xi_L + 0.5*dt*k2xiL,
                               xi_R + 0.5*dt*k2xiR)
    k4x, k4p, k4xiL, k4xiR = f(x + dt*k3x,
                               p + dt*k3p,
                               xi_L + dt*k3xiL,
                               xi_R + dt*k3xiR)
    x_new = x + dt*(k1x + 2*k2x + 2*k3x + k4x)/6.0
    p_new = p + dt*(k1p + 2*k2p + 2*k3p + k4p)/6.0
    xi_L_new = xi_L + dt*(k1xiL + 2*k2xiL + 2*k3xiL + k4xiL)/6.0
    xi_R_new = xi_R + dt*(k1xiR + 2*k2xiR + 2*k3xiR + k4xiR)/6.0
    return x_new, p_new, xi_L_new, xi_R_new


# ------------------------------------------------------
# ESECUZIONE
# ------------------------------------------------------
# ------------------------------------------------------
# PARAMETRI DEL SISTEMA FPUT CON TERMOSTATI
# ------------------------------------------------------
N = 110               # numero di masse mobili
m = 1.0
a = 1.0               # distanza di equilibrio
chi = 1.0
alpha = 1.0
beta = 1.0
grad_T=0.1

Tl = 1.0    # temperature dei termostati sinistro e destro
Tr=Tl+N*grad_T

thetaL, thetaR = 1.0, 1.0

dt = 0.01
t_steps = 200000
save_every = 1

n_chains = 1000000     # simulazioni parallele

filename = f"../condizioni_{N}.bin"
x0, p0, xiL0, xiR0 = read_conditions_fput_parallel(filename, n_chains, N)

def omega0_fn(T):
    return xiR0*p0[:,-2]**2*(1/Tr-1/T)+xiL0*p0[:,1]**2*(1/Tl-1/T)

omega0=omega0_fn(1.)

def observable_bulk(x, p):    
    bd_paticle = int(N * 0.15)
    segment = x[:, bd_paticle: N - bd_paticle + 1]     
    r = jnp.diff(segment, axis=1) - a                   
    flux = (chi*r + alpha*r**2 + beta*r**3) * p[:, bd_paticle : N - bd_paticle] / m
    return flux.mean(axis=1)                          


# ------------------------------------------------------
# SIMULAZIONE MULTI-CATENA
# ------------------------------------------------------
def simulate_multi(x0, p0, xiL0, xiR0, dt, t_steps, save_every):
    """
    Esegue t_steps passi di integrazione con RK4 e ritorna solo le osservabili salvate ogni save_every passi.
    """
    n_save = t_steps // save_every

    def body(carry, step):
        x, p, xiL, xiR, store = carry
        x, p, xiL, xiR = rk4_step(x, p, xiL, xiR, dt)
        obs_val = (omega0 * observable_bulk(x, p)).mean()
        store = jax.lax.cond(
            (step % save_every) == 0,
            lambda s: s.at[step // save_every].set(obs_val),
            lambda s: s,
            store
        )
        return (x, p, xiL, xiR, store), None

    store = jnp.zeros(n_save)
    init = (x0, p0, xiL0, xiR0, store)
    (x_final, p_final, xiL_final, xiR_final, store), _ = lax.scan(
        body, init, jnp.arange(t_steps)
    )
    return x_final, p_final, xiL_final, xiR_final, store


# ------------------------------------------------------
# ESECUZIONE A BLOCCHI
# ------------------------------------------------------
simulate_multi_jit = jax.jit(simulate_multi, static_argnames=("t_steps", "save_every"))

results_dir = "results"
os.makedirs(results_dir, exist_ok=True)


block_size = 1000

n_blocks = t_steps // block_size
all_store=[]
start = time.time()
print(f"Inizio simulazione ({n_blocks} blocchi da {block_size} passi)...")

x, p, xiL, xiR = x0, p0, xiL0, xiR0

#for b in tqdm(range(n_blocks), desc="Simulazione a blocchi"):
x, p, xiL, xiR, store_block = simulate_multi_jit(x, p, xiL, xiR, dt, block_size, save_every)
    #all_store.append(store_block)
end = time.time()
print(f"Simulazione completata in {end - start:.2f} s")
start = time.time()
x, p, xiL, xiR, store_block = simulate_multi_jit(x, p, xiL, xiR, dt, block_size, save_every)

end = time.time()
print(f"Simulazione completata in {end - start:.2f} s")
#final_store = np.array(jnp.concatenate(all_store).block_until_ready(), dtype=np.float32)

print(f"Simulazione completata in {end - start:.2f} s")
#np.save(os.path.join(results_dir, "store_full.npy"), final_store)

# print(f"Risultati salvati in: {results_dir}/store_block_XXXXX.npy")

# ------------------------------------------------------
# UNIONE (opzionale, dopo la simulazione)
# ------------------------------------------------------
# def concat_results(results_dir):
#     files = sorted(f for f in os.listdir(results_dir) if f.startswith("store_block_") and f.endswith(".npy"))
#     arrays = [np.load(os.path.join(results_dir, f)) for f in files]
#     full_store = np.concatenate(arrays)

#     out_name = os.path.join(results_dir, "store_full.npy")
#     np.save(out_name, full_store)
#     print(f"File completo salvato in: {out_name}")

#     # elimina i file dei blocchi
#     for f in files:
#         os.remove(os.path.join(results_dir, f))
#     print("File dei blocchi intermedi eliminati.")

#     return full_store


#concat_results(results_dir)
