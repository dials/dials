def work_dqe():
  from parallax import dqe, derive_absorption_coefficient_Si
  t0 = 0.032
  for j in range(30, 201, 1):
    energy_kev = 0.1 * j
    mu = derive_absorption_coefficient_Si(energy_kev)
    print energy_kev, dqe(0.032, 0.0, mu), dqe(0.05, 0.0, mu), dqe(0.1, 0.0, mu)

if __name__ == '__main__':
  work_dqe()
