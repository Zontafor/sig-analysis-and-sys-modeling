import numpy as np


def hL(t):
    t = np.asarray(t)
    return np.where(t >= 0.0, np.exp(-t), 0.0)


def fft_convolve(a, b, dt):
    n = len(a) + len(b) - 1
    nfft = 1 << (n - 1).bit_length()
    A = np.fft.rfft(a, nfft)
    B = np.fft.rfft(b, nfft)
    y = np.fft.irfft(A * B, nfft)[:n]
    return y * dt


def main():
    t_max = 20.0
    dt = 1e-5
    t = np.arange(0.0, t_max + dt, dt)
    u = np.ones_like(t)

    y_NL_num = fft_convolve(hL(t), u**2, dt)[: len(t)]
    y_NL_true = 1.0 - np.exp(-t)

    v_LN = fft_convolve(hL(t), u, dt)[: len(t)]
    y_LN_num = v_LN**2
    y_LN_true = (1.0 - np.exp(-t)) ** 2

    print("Problem 2 numeric backcheck")
    print(f"  max|y_NL_num - y_NL_true| = {np.max(np.abs(y_NL_num - y_NL_true)):.6e}")
    print(f"  max|y_LN_num - y_LN_true| = {np.max(np.abs(y_LN_num - y_LN_true)):.6e}")

    n_samples = t.size
    print(f"  dt = {dt:.3e}")
    print(f"  t_max = {t_max:.3e}")
    print(f"  n_samples = {n_samples:d}")


if __name__ == "__main__":
    main()