#!/usr/bin/env python3
# p1.py
# Numeric backcheck for Problem 1 (K=2).

import numpy as np


def h(t):
    t = np.asarray(t)
    return np.where(t >= 0.0, (np.exp(-t) - np.exp(-6.0 * t)) / 5.0, 0.0)


def y_step_closed(t):
    t = np.asarray(t)
    return np.where(
        t >= 0.0,
        1.0 / 6.0 - (1.0 / 5.0) * np.exp(-t) + (1.0 / 30.0) * np.exp(-6.0 * t),
        0.0,
    )


def y_cos_closed(t):
    t = np.asarray(t)
    return np.where(
        t >= 0.0,
        (1.0 / 100.0) * np.cos(2.0 * t)
        + (7.0 / 100.0) * np.sin(2.0 * t)
        - (1.0 / 25.0) * np.exp(-t)
        + (3.0 / 100.0) * np.exp(-6.0 * t),
        0.0,
    )


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
    ht = h(t)

    y_step_num = fft_convolve(ht, u, dt)[: len(t)]
    y_step_true = y_step_closed(t)

    x_cos = np.cos(2.0 * t)
    y_cos_num = fft_convolve(ht, x_cos, dt)[: len(t)]
    y_cos_true = y_cos_closed(t)

    print("Problem 1 numeric backcheck (K=2)")
    print(f"  max|y_step_num - y_step_true| = {np.max(np.abs(y_step_num - y_step_true)):.6e}")
    print(f"  max|y_cos_num  - y_cos_true | = {np.max(np.abs(y_cos_num  - y_cos_true )):.6e}")

    n_samples = t.size
    print(f"  dt = {dt:.3e}")
    print(f"  t_max = {t_max:.3e}")
    print(f"  n_samples = {n_samples:d}")

if __name__ == "__main__":
    main()