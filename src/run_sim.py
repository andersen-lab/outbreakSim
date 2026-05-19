from os import sep
from pathlib import Path
import random
from Bio import SeqIO
import covasim as cv
import covasim.data as cvdata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import epyestim.covid19 as covid19
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
import argparse

from molecular_clock import molecular_clock_evolve
from substitution_model import JukesCantor

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_REFERENCE = _REPO_ROOT / "assets" / "NC_045512_Hu-1.fasta"


def load_fasta_sequence(path: Path) -> str:
    """Load the first FASTA record as a single uppercase nucleotide string."""
    record = next(SeqIO.parse(path, "fasta"))
    return str(record.seq).upper()


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


def initiate_infection(sim, region_id):
    sim.people.infect(inds=np.where(sim.people.region == region_id))
    return sim


def calculate_new_infections(sim, n_regions):
    new_cases = np.zeros((len(sim.results['new_infections']), n_regions))
    for t in range(len(sim.results['new_infections'])):
        # infected in each time point
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


def viral_shedding_simple(new_cases):
    ans = solve_eclipse()
    # create an array of zeros the same simze of new cases
    ww_signal = np.zeros_like(new_cases)
    # this calculates viral shedding in each day per region
    for i in range(new_cases.shape[1]):
        ww_signal[:, i] = np.convolve(new_cases[:, i], ans[::-1], mode='same')
    return ww_signal


def get_viral_loads(sim, t_start, t_end):
    # Code credit to https://github.com/amath-idm/covasim_methods_paper/blob/main/figs/viral_load_plot.py
    viralLoad = np.zeros((len(sim.people.date_exposed),t_end-t_start))
    par1 = sim.pars['beta_dist']['par1']
    par2 = sim.pars['beta_dist']['par2']
    mean  = np.log(par1**2 / np.sqrt(par2 + par1**2)) # Computes the mean of the underlying normal distribution
    sigma = np.sqrt(np.log(par2/par1**2 + 1)) # Computes sigma for the underlying normal distribution
    rel_trans = np.random.lognormal(mean=mean, sigma=sigma, size=len(sim.people.date_exposed))

    for i in np.arange(t_end-t_start):
        viralLoad[:,i] = cv.utils.compute_viral_load(t_start+i, sim.people.date_infectious,\
                 sim.people.date_recovered, sim.people.date_dead,\
                 sim.pars['viral_dist']['frac_time'], sim.pars['viral_dist']['load_ratio'],\
                 sim.pars['viral_dist']['high_cap'])*rel_trans
    return viralLoad.T

  
def viral_shedding_covasim(sim,start,end):
    viral_load_matrix = get_viral_loads(sim,start,end)
    n_days, n_people = viral_load_matrix.shape
    n_regions = sim.people.region.max() + 1
    regional_viral_load = np.zeros((n_days, n_regions))
    for t in range(n_days):
        for p in range(n_people):
            r = sim.people.region[p]  # Get the region for this specific person
            regional_viral_load[t, r] += viral_load_matrix[t, p]
    return regional_viral_load

def assign_sequences(sim, reference):
    sub_model = JukesCantor()
    seq_rng = random.Random(42)

    df = pd.DataFrame(sim.people.infection_log)
    df.to_csv('infection_log.tsv', sep='\t')
    # Track each person's currently carried lineage so reinfections can replace it.
    person_current_seq = {}
    person_current_start_date = {}

    # Store sequence by transmission event (source -> target).
    event_seq_by_row = {}

    # Process in temporal order so parent lineage state is up to date.
    for idx, row in df.sort_values("date", kind="stable").iterrows():
        target = int(row["target"])
        event_date = float(row["date"])

        if pd.isna(row["source"]):
            event_seq_by_row[idx] = reference
            person_current_seq[target] = reference
            person_current_start_date[target] = event_date
            continue

        source = int(row["source"])
        parent_seq = person_current_seq.get(source, reference)
        parent_start_date = person_current_start_date.get(source, event_date)
        branch_time = max(0.0, event_date - parent_start_date)

        child_seq = molecular_clock_evolve(
            parent_seq,
            branch_time,
            rate=1e-7,  # mutations per site per day
            model=sub_model,
            rng=seq_rng,
        )

        event_seq_by_row[idx] = child_seq
        person_current_seq[target] = child_seq
        person_current_start_date[target] = event_date

    df["sequence"] = pd.Series(event_seq_by_row)

    return df
    


def plot_shedding(data,name):
    fig, ax = plt.subplots(figsize=(6,10))
    im = ax.imshow(data, aspect='auto', cmap='viridis')

    # Discrete x-axis labels
    x_labels = ['1', '2', '3', '4']   # replace with your actual category names
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels)

    ax.set_xlabel("regions")
    ax.set_ylabel("days")

    fig.colorbar(im, ax=ax, label="viral load")
    plt.savefig(f'{name}.pdf')


def main():
    parser = argparse.ArgumentParser(
        description="Define simulation parameters"
    )

    parser.add_argument("--pop_size", type=int, required=True,
                        help="Population size")

    parser.add_argument("--pop_type", type=str, required=True,
                        help="Population type (e.g. random, hybrid, etc.)")

    parser.add_argument("--n_days", type=int, required=True,
                        help="Number of simulation days")

    parser.add_argument("--location", type=str, required=True,
                        help="Location name")

    parser.add_argument("--pop_infected", type=int, default=0,
                        help="Initial infected population (default: 0)")
    parser.add_argument("--n_imports", type=int, default=0,
                        help="Number of imported infections (default: 0)")
    parser.add_argument("--n_rows", type=int, default=2,
                        help="Number of rows in the region (default: 2)")
    parser.add_argument("--n_cols", type=int, default=2,
                        help="Number of columns in the region (default: 2)")
    parser.add_argument("--n_init_inf", type=int, default=20,
                        help="Number of initial infections in the region (default: 20)")
    parser.add_argument("--r_init_inf", type=int, default=0,
                        help="Row index for where the infection initiated (default: 0)")
    parser.add_argument("--c_init_inf", type=int, default=0,
                        help="column index for where the infection initiated (default: 0)")
    parser.add_argument("--reference", type=Path, default=_DEFAULT_REFERENCE,
                        help="Path to reference genome FASTA (default: NC_045512_Hu-1.fasta)")

    args = parser.parse_args()
    ### STEP 1 ####
    pars = define_sim_parameters(
        pop_size=args.pop_size,
        pop_type=args.pop_type,
        n_days=args.n_days,
        location=args.location,
        pop_infected=args.pop_infected,
        n_imports=args.n_imports,
    )
    check_age_household_dist(args.location)
    # initialize the simulation with the parameters
    sim = cv.Sim(pars)
    sim.initialize()
    # Define regions for the simulation
    n_rows, n_cols = args.n_rows, args.n_cols
    n_regions = n_rows * n_cols
    # Assign people to each region
    sim = assign_people(sim, n_regions, args.pop_size)
    # initiate infection
    # 20 initial infections in Region_0_0
    sim = initiate_infection(sim, 0)
    sim.run()
    #### STEP 2 ####
    # calculate number of new cases in each region
    new_cases = calculate_new_infections(sim, n_regions)
    # calculate wastewater shedding per region per time point using
    # basic shedding model
    # shedding_simple = viral_shedding_simple(new_cases)
    # # covasim viral load model (adjust zero if start date is not day 0)
    # shedding_covasim = viral_shedding_covasim(sim,0,args.n_days)
    # plot values
    # plot_shedding(shedding_simple,"simple")
    # plot_shedding(shedding_covasim,"covasim")

    reference_seq = load_fasta_sequence(args.reference)
    
    transmission_df = assign_sequences(sim, reference_seq)

    print(transmission_df['sequence'].nunique())
    
    transmission_df.to_csv("transmission.tsv", sep='\t')

if __name__ == "__main__":
    main()