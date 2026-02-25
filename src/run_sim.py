import covasim as cv
import covasim.data as cvdata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import epyestim.covid19 as covid19
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout

def define_sim_parameters(pop_size, pop_type,
                          n_days, location,
                          pop_infected, n_imports):
    pars = dict(
    pop_size = pop_size,
    pop_type = pop_type,
    n_days = n_days,
    location = location,
    pop_infected = pop_infected,   # no automatic initial infections
    n_imports = n_imports,      # no imported infections
    )
    return pars


def check_age_household_dist(location):
    try:
        age_data = cvdata.get_age_distribution(location)
        print(f"Age distribution for {location}:")
        print(age_data)
    except ValueError as e:
        print(f"No age data for {location}, using default. ({e})")

    try: 
        household_size = cvdata.get_household_size(location) 
        print(f"Household size distribution for {location}:") 
        print(household_size) 
    except ValueError as e: 
        print(f"No household data for {location}, using default. ({e})")


def assign_people(sim, n_regions, pop_size):
    # assign each person to a region
    region_ids = np.repeat(np.arange(n_regions), pop_size // n_regions)
    np.random.shuffle(region_ids)
    sim.people.region = region_ids
    return sim


def initiate_infection(sim, row_index, col_index, n_infected):
    sim.people.infect(inds=np.where(sim.people.region == row_index)[col_index][:n_infected])
    return sim


def calculate_new_infections(sim, n_regions):
    new_cases = np.zeros((len(sim.results['new_infections']), n_regions))
    for t in range(len(sim.results['new_infections'])):
        newly_infected = np.where(sim.people.date_infectious == t)[0]
        for r in range(n_regions):
            new_cases[t, r] = np.sum(sim.people.region[newly_infected] == r)
    return new_cases


def eclipse_model(y, t, b, k, delta, p, mu, c):
    T, I1, I2, Vi, Vni = y
    dydt = [-b*Vi*T, b*Vi*T - k*I1, k*I1 - delta*I2, p*mu*I2 - c*Vi, p*(1.-mu)*I2 - c*Vni]
    return dydt


def solve_eclipse():
    b = 5e-5
    c_clear = 10
    k = 6
    mu = 1e-4
    p = 1e5
    delta = 0.5
    t = np.linspace(0, 40, 41)
    y0 = [1.33e5, 0, 1/30, 0, 0]
    # solve ode
    sol = odeint(eclipse_model, y0, t, args=(b, k, delta, p, mu, c_clear))
    c = sol[:,3:].sum(axis=1)
    c = c/c.sum()
    return c


def main():
    ### STEP 1 ####
    pars = define_sim_parameters(100, "hybrid",
                                 180, "Zambia",
                                 0, 0)
    check_age_household_dist("Zambia")
    # initialize the simulation with the parameters
    sim = cv.Sim(pars)
    sim.initialize()
    # Define regions for the simulation
    n_rows, n_cols = 2, 2
    n_regions = n_rows * n_cols
    # Assign people to each region
    sim = assign_people(sim, n_regions, 180)
    # initiate infection
    # 20 initial infections in Region_0_0
    sim = initiate_infection(sim, 0, 0, 20)
    sim.run()
    #### STEP 2 ####
    # calculate number of new cases in each region
    new_cases = calculate_new_infections(sim, n_regions)
    ans = solve_eclipse()
    # create an array of zeros the same simze of new cases
    ww_signal = np.zeros_like(new_cases)
    # this calculates viral shedding in each day per region
    for i in range(new_cases.shape[1]):
        ww_signal[:, i] = np.convolve(new_cases[:, i], ans[::-1], mode='same')


if __name__ == "__main__":
    main()