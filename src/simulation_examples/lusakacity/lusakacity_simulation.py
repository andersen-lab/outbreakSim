import collections, dataclasses
import glob
import os
import baltic as bt
import covasim as cv
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from PIL import Image

# Set working directory to script location so relative data paths resolve correctly
os.chdir(os.path.dirname(os.path.abspath(__file__)))

plot_dir = 'outputs'  # All plots are saved here as PDFs instead of being shown
os.makedirs(plot_dir, exist_ok=True)

# ==========================================
# 1. Viral shedding per agent heatmap
# ==========================================
sim = cv.Sim(
    pop_size=20_000,
    pop_type='rank',
    density_raster='data/lusaka_city_raster.tif',
    pop_infected=20,
    location='zambia',
    n_days=100,
    n_imports=0,
    boundary_shapefile='data/lusaka_city_adm4.geojson',
    sewershed_file='data/lusaka_city_sewershed.geojson',
)
sim.run()

viral_shedding = sim.results['viral_shedding_hist']
im = plt.imshow(
    viral_shedding[:100],  # first 100 agents only
    aspect='auto',  # prevents squashing
    origin='lower',  # agent 0 at bottom
    cmap='viridis',
)
plt.colorbar(im, label='Viral shedding')
plt.xlabel('time')
plt.ylabel('agent')
plt.title('Viral shedding by agent and time (first 100 agents)')
plt.tight_layout()
plt.savefig(f'{plot_dir}/01_shedding_heatmap.pdf')
plt.close()


# ==========================================
# 2. Viral shedding values over time by sewershed
# ==========================================
sewersheds = np.asarray(sim.people.sewershed)
plt.figure(figsize=(12, 6))
for sid in np.unique(sewersheds):
  mask = sewersheds == sid
  name = sim.sewershed_labels.get(
      sid, 'Outside all sewersheds'
  )  # IDs of sewersheds in the file -> names; 0 and -1 are not in the file
  plt.plot(
      viral_shedding[mask].sum(axis=0),
      label=f'{name} ({mask.sum()} people)',
  )
plt.plot(
    viral_shedding.sum(axis=0),
    linewidth=2,
    color='red',
    label='Total',
)
plt.xlabel('Day')
plt.ylabel('Total viral shedding')
plt.title('Total viral shedding by sewershed')
plt.legend()
plt.tight_layout()
plt.savefig(f'{plot_dir}/02_shedding_by_sewershed.pdf')
plt.close()


# ==========================================
# 3. Multi-variant simulation & fitness tracking
# ==========================================
ref_fasta = 'data/sars-2-data/NC_045512_Hu-1.fasta'
fitness_data_path = 'data/sars-2-data/nt_fitness.csv'
ay4_fasta = 'data/sars-2-data/variants/AY.4_consensus.fasta'
ba2_fasta = 'data/sars-2-data/variants/BA.2.9_consensus.fasta'
xbb15_fasta = 'data/sars-2-data/variants/XBB.1.5.18_consensus.fasta'
artic_primers = 'data/sars-2-data/ARTIC_V4-1.bed'

# Variants — rel_beta is derived from founding mutations via the Bloom lab fitness landscape
ay4 = cv.variant(
    label='AY.4',
    days=30,  # introduce AY.4 on day 30
    n_imports=50,
    founding_fasta=ay4_fasta,
)
ba2 = cv.variant(
    label='BA.2',
    days=60,  # introduce BA.2 on day 60
    n_imports=50,
    founding_fasta=ba2_fasta,
)
xbb15 = cv.variant(
    label='XBB.1.5',
    days=100,  # introduce XBB.1.5 on day 100
    n_imports=50,
    founding_fasta=xbb15_fasta,
)

# WastewaterSampler captures each unique evolved haplotype
_haplotype_days = list(range(110, 151, 5))  # days 110, 115, 120 … 150
ww = cv.WastewaterSampler(
    days=_haplotype_days,
    label='WastewaterHaplotype',
)

sim = cv.Sim(
    dict(
        pop_size=20_000,
        location='zambia',
        pop_type='rank',
        pop_infected=20,  # initial infections
        start_day='2022-10-01',
        end_day='2023-04-01',
        evo_pars=dict(
            enable=True,
            reference=ref_fasta,
            mol_clock_rate=1e-6,
            fitness_model='bloom_nt',
            fitness_data_path=fitness_data_path,
            fitness_scale=0.05,
        ),
    ),
    variants=[ay4, ba2, xbb15],
    analyzers=[ww],
    verbose=0,
)
sim.run()

# Show rel_beta multiplier derived from each variant's founding mutations
tracker = sim.people.sequence_tracker
fitness_model = tracker.fitness_model

print(f'{"variant":12s}  {"founding SNPs":>13s}  {"rel_beta (×wt)":>15s}')
print('-' * 46)
for label, muts in [('wild', frozenset())] + list(
    tracker.variant_founding_mutations.items()
):
  mult = fitness_model.compute_fitness(muts)
  print(f'{label:12s}  {len(muts):>13d}  {mult:>15.4f}')

variant_names = list(
    sim['variant_map'].values()
)  # e.g. ['wild', 'AY.4', 'BA.2', 'XBB.1.5']
n_steps = sim.npts

n_inf = np.array(sim.results['variant']['n_infectious_by_variant'].values)
total = n_inf.sum(axis=0)
total[total == 0] = 1  # avoid division by zero

loads = {v: n_inf[i] for i, v in enumerate(variant_names)}
props = {v: n_inf[i] / total for i, v in enumerate(variant_names)}

# Color palette
palette = {
    'wild': '#4C72B0',
    'AY.4': '#DD8452',
    'BA.2': '#C44E52',
    'XBB.1.5': '#55A868',
}
_defaults = plt.rcParams['axes.prop_cycle'].by_key()['color']
for i, v in enumerate(variant_names):
  palette.setdefault(v, _defaults[i % len(_defaults)])

tvec = np.arange(n_steps)
intro_days = {'AY.4': 30, 'BA.2': 60, 'XBB.1.5': 100}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# Absolute infectious counts
for vname in variant_names:
  ax1.plot(
      tvec, loads[vname], label=vname, color=palette[vname], linewidth=1.8
  )
for vname, day in intro_days.items():
  ax1.axvline(day, color=palette[vname], linestyle=':', linewidth=1.2, alpha=0.7)
ax1.set_ylabel('Infectious agents')
ax1.set_title('Infectious count by variant — absolute')
ax1.legend(framealpha=0.9)

# Scaled to 1 (stacked area)
stacked = np.array([props[v] for v in variant_names])
ax2.stackplot(
    tvec,
    stacked,
    labels=variant_names,
    colors=[palette[v] for v in variant_names],
    alpha=0.85,
)
for vname, day in intro_days.items():
  ax2.axvline(day, color='white', linestyle=':', linewidth=1.5, alpha=0.9)
ax2.set_xlabel('Simulation day')
ax2.set_ylabel('Proportion of infectious agents')
ax2.set_title('Variant proportion over time — scaled to 1')
ax2.set_ylim(0, 1)
ax2.legend(loc='upper left', framealpha=0.9)

plt.tight_layout()
plt.savefig(f'{plot_dir}/03_variants_counts_and_proportions.pdf')
plt.close()


# ==========================================
# 4. Haplotype trimming & Bygul read simulation
# ==========================================
ww = sim.get_analyzer('WastewaterHaplotype')
mixed_day = next(
    (
        t
        for t in sorted(ww.samples)
        if ww.samples[t] is not None
        and all(props[v][t] > 0.01 for v in ['wild', 'AY.4', 'BA.2', 'XBB.1.5'])
    ),
    120,
)

print(f'Sampling day {mixed_day}  ({sim.date(mixed_day)})')
print('  Variant-level proportions at that day:')
for v in ['wild', 'AY.4', 'BA.2', 'XBB.1.5']:
  print(f'    {v:10s}  {props[v][mixed_day]:.1%}')

sample = ww.samples.get(mixed_day)
if sample is None:
  print('\nNo infectious agents were present on that day.')
else:
  print(f'\n  Unique haplotypes in pool : {sample.n_genotypes}')

  order = sorted(
      range(sample.n_genotypes),
      key=lambda i: sample.proportions[i],
      reverse=True,
  )[:5]
  top_seqs = [sample.mutation_sets[i] for i in order]
  top_loads = [sample.raw_loads[i] for i in order]
  top_props = [sample.proportions[i] for i in order]
  top_total = sum(top_props)
  top_props = [p / top_total for p in top_props]
  top_props[-1] = 1.0 - sum(top_props[:-1])  # absorb floating-point residual

  print(
      '  Top-5 proportions (renormalized):'
      f' {[round(p, 4) for p in top_props]}'
  )
  print(f'  These cover {top_total:.1%} of the full pool by viral shedding')

  ww.samples[mixed_day] = dataclasses.replace(
      sample,
      mutation_sets=top_seqs,
      proportions=top_props,
      raw_loads=top_loads,
      n_genotypes=len(top_seqs),
  )

outdir = './amplicon_reads/'
try:
  ww.simulate_sample(
      day=mixed_day,
      primers=artic_primers,
      reference=ref_fasta,
      outdir=outdir,
      readcnt=300,
      redo=True,
      simulation_mode='amplicon',
      seed=42,
  )
  fastq_files = (
      glob.glob(os.path.join(outdir, '**/*.fastq*'), recursive=True)
      + glob.glob(os.path.join(outdir, '*.fastq*'))
  )
  print(f'\nBygul wrote {len(fastq_files)} FASTQ file(s) to {outdir!r}:')
  for f in sorted(set(fastq_files))[:12]:
    print(f'  {os.path.relpath(f)}  ({os.path.getsize(f)/1024:.1f} kB)')
except ImportError as exc:
  print('\n[Bygul not installed — skipping read simulation]')
  print(str(exc))
  print('Install: pip install bygul   (https://github.com/andersen-lab/Bygul)')


# ==========================================
# 5. Dual Clinical & Wastewater Simulation
# ==========================================
sample_days = list(range(40, 100, 10))
dual = cv.Sim(
    dict(
        pop_size=20_000,
        pop_infected=20,
        n_days=100,
        pop_type='rank',
        location='zambia',
        density_raster='data/lusaka_city_raster.tif',
        evo_pars=dict(
            enable=True,
            reference=ref_fasta,
            mol_clock_rate=1e-6,
            fitness_model='bloom_nt',
            fitness_data_path=fitness_data_path,
            fitness_scale=0.05,
        ),
        boundary_shapefile='data/lusaka_city_adm4.geojson',
        sewershed_file='data/lusaka_city_sewershed.geojson',
    ),
    variants=[ay4],
    interventions=[cv.test_prob(symp_prob=0.2, start_day=0)],
    analyzers=[
        cv.ClinicalSequencer(days=sample_days, n_samples=50, label='clinical')
    ]
    + [
        cv.WastewaterSampler(days=sample_days, sid=s, label=f'ww_sewershed_{s}')
        for s in (1, 2, 3)
    ],
    verbose=0,
)
dual.run()

clinical = dual.get_analyzer('clinical')
for day in sample_days:
  row = f'day {day}: clinical n={len(clinical.samples.get(day, [])):3d}'
  for s in (1, 2, 3):
    smp = dual.get_analyzer(f'ww_sewershed_{s}').samples.get(day)
    row += (
        f' | sewershed {s}: '
        + (f'{smp.n_genotypes} haplotypes' if smp else 'no shedders')
    )
  print(row)

sewershed_gdf = gpd.read_file('data/lusaka_city_sewershed.geojson')
sewershed_gdf['sid_num'] = sewershed_gdf['sid'].map(
    {sid: i for i, sid in dual.sewershed_labels.items()}
)
region_gdf = gpd.read_file('data/lusaka_city_adm4.geojson')
region_gdf['region_num'] = region_gdf['region_code'].astype(int)

outline = gpd.GeoSeries([region_gdf.union_all()], crs=region_gdf.crs)


def add_unmonitored_outline(ax):
  '''Hatch the whole area beneath sewersheds for unmonitored zones.'''
  outline.plot(
      ax=ax,
      facecolor='white',
      edgecolor='0.6',
      hatch='///',
      linewidth=0.8,
      zorder=0,
  )
  ax.legend(
      handles=[
          Patch(
              facecolor='white',
              edgecolor='0.6',
              hatch='///',
              label='Not in a sewershed',
          )
      ],
      loc='upper left',
      bbox_to_anchor=(0, 0.02),
      frameon=False,
      fontsize=8,
  )


people = dual.people
shedding_by_agent = dual.results['viral_shedding_hist'].sum(axis=1)
diagnosed = ~np.isnan(people.date_diagnosed)

for sid in sewershed_gdf['sid_num']:
  sewershed_gdf.loc[
      sewershed_gdf['sid_num'] == sid, 'shedding_per_capita'
  ] = shedding_by_agent[people.sewershed == sid].mean()
sewershed_gdf['label'] = sewershed_gdf['sid']

for code in region_gdf['region_num']:
  in_region = people.region == code
  row = region_gdf['region_num'] == code
  region_gdf.loc[row, 'cases'] = diagnosed[in_region].sum()
  region_gdf.loc[row, 'cases_per_1000'] = (
      1000 * diagnosed[in_region].sum() / in_region.sum()
      if in_region.any()
      else np.nan
  )

fig, (ax_cl, ax_ww) = plt.subplots(1, 2, figsize=(14, 7))
panels = [
    (
        ax_cl,
        region_gdf,
        'cases_per_1000',
        'Cumulative diagnosed cases per 1,000 residents',
        'Clinical: by admin level 4',
    ),
    (
        ax_ww,
        sewershed_gdf,
        'shedding_per_capita',
        'Cumulative Wastewater viral load (a.u.)',
        'Wastewater: by sewershed',
    ),
]
for ax, gdf, column, label, title in panels:
  gdf.plot(
      column=column,
      cmap='YlOrRd',
      vmin=0,
      edgecolor='0.3',
      legend=True,
      legend_kwds={
          'label': label,
          'orientation': 'horizontal',
          'shrink': 0.7,
      },
      ax=ax,
  )
  if 'label' in gdf.columns:
    for _, row in gdf.iterrows():
      c = row.geometry.centroid
      dark = row[column] > 0.5 * gdf[column].max()
      ax.annotate(
          row['label'],
          xy=(c.x, c.y),
          ha='center',
          fontsize=7,
          fontweight='bold',
          color='white' if dark else 'black',
      )
  ax.set_title(title)
  ax.set_axis_off()
add_unmonitored_outline(ax_ww)
plt.tight_layout()
plt.savefig(f'{plot_dir}/04_clinical_vs_wastewater_map.pdf')
plt.close()


# ==========================================
# 6. Small simulation & time-series animation
# ==========================================
class RecordViralLoad(cv.Analyzer):
  '''Save every agent's viral load on every day'''

  def initialize(self, sim):
    super().initialize(sim)
    self.viral_load = np.zeros((sim['pop_size'], sim.npts))

  def apply(self, sim):
    self.viral_load[:, sim.t] = sim.people.viral_load


small = cv.Sim(
    pop_size=1000,
    pop_infected=20,
    n_days=100,
    pop_type='rank',
    location='zambia',
    density_raster='data/lusaka_city_raster.tif',
    boundary_shapefile='data/lusaka_city_adm4.geojson',
    sewershed_file='data/lusaka_city_sewershed.geojson',
    interventions=[cv.test_prob(symp_prob=0.2, start_day=0)],
    analyzers=[RecordViralLoad(label='viral_load')],
    verbose=0,
)
small.run()
viral_load_hist = small.get_analyzer('viral_load').viral_load
small_people = small.people
small_diagnosed = ~np.isnan(small_people.date_diagnosed)
n_days = small.npts
shedding_hist = small.results['viral_shedding_hist']

clinical_ts = np.zeros((len(region_gdf), n_days))
wastewater_ts = np.zeros((len(sewershed_gdf), n_days))
for k, code in enumerate(region_gdf['region_num']):
  in_region = small_people.region == code
  diag_days = small_people.date_diagnosed[
      in_region & small_diagnosed
  ].astype(int)
  clinical_ts[k] = (
      1000
      * np.bincount(diag_days[diag_days < n_days], minlength=n_days)
      / in_region.sum()
      if in_region.any()
      else np.nan
  )
for k, sid in enumerate(sewershed_gdf['sid_num']):
  wastewater_ts[k] = shedding_hist[small_people.sewershed == sid].mean(axis=0)


def smooth(ts):
  return (
      pd.DataFrame(ts.T).rolling(7, min_periods=1).mean().to_numpy().T
  )  # 7-day trailing mean


clinical_ts, wastewater_ts = smooth(clinical_ts), smooth(wastewater_ts)

agent_cmap = plt.get_cmap('RdPu')
agent_norm = Normalize(vmin=0, vmax=viral_load_hist.max())
agent_xy = np.column_stack((small_people.x, small_people.y))


def agents_on_day(day):
  load = viral_load_hist[:, day]
  infectious = load > 0
  return agent_xy[~infectious], agent_xy[infectious], load[infectious]


fig, (ax_agents, ax_clinical, ax_wastewater) = plt.subplots(
    1, 3, figsize=(15, 5.5)
)
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.02, top=0.9, wspace=0.02)
region_gdf.boundary.plot(ax=ax_agents, color='0.5', linewidth=1, zorder=2)
not_infectious, infectious_xy, infectious_load = agents_on_day(0)
grey_scatter = ax_agents.scatter(
    not_infectious[:, 0],
    not_infectious[:, 1],
    color='lightgrey',
    s=12,
    zorder=1,
)
agent_scatter = ax_agents.scatter(
    infectious_xy[:, 0],
    infectious_xy[:, 1],
    c=infectious_load,
    cmap=agent_cmap,
    norm=agent_norm,
    s=12,
    zorder=3,
)
fig.colorbar(
    agent_scatter,
    ax=ax_agents,
    orientation='horizontal',
    shrink=0.7,
    label='Viral load (relative transmissibility)',
)
ax_agents.set_title('Agents')
ax_agents.set_axis_off()

panels = [
    (
        ax_clinical,
        region_gdf,
        clinical_ts,
        'Clinical (by admin level 4)',
        'New diagnoses per 1,000 residents',
    ),
    (
        ax_wastewater,
        sewershed_gdf,
        wastewater_ts,
        'Wastewater (by sewershed)',
        'Wastewater viral load (a.u.)',
    ),
]
for ax, gdf, ts, title, label in panels:
  gdf['value'] = np.nan_to_num(ts[:, 0])
  gdf.plot(
      column='value',
      cmap='YlOrRd',
      vmin=0,
      vmax=np.nanmax(ts),
      edgecolor='0.3',
      ax=ax,
  )
  fig.colorbar(
      plt.cm.ScalarMappable(
          norm=plt.Normalize(0, np.nanmax(ts)), cmap='YlOrRd'
      ),
      ax=ax,
      orientation='horizontal',
      shrink=0.7,
      label=label,
  )
  ax.set_title(title)
  ax.set_axis_off()

add_unmonitored_outline(ax_wastewater)
xmin, ymin, xmax, ymax = region_gdf.total_bounds
for ax in (ax_agents, ax_clinical, ax_wastewater):
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)


def update(day):
  not_infectious, infectious_xy, infectious_load = agents_on_day(day)
  grey_scatter.set_offsets(not_infectious)
  agent_scatter.set_offsets(infectious_xy)
  agent_scatter.set_array(infectious_load)
  for ax, gdf, ts, title, label in panels:
    ax.collections[0].set_array(ts[:, day])
  fig.suptitle(f'Day {day}')


def save_gif(fig, update, n_frames, path, duration=100):
  frames = []
  for i in range(n_frames):
    update(i)
    fig.canvas.draw()
    frames.append(
        Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
    )
  plt.close(fig)

  sample = frames[::5]
  stack = Image.new('RGB', (sample[0].width, sample[0].height * len(sample)))
  for i, frame in enumerate(sample):
    stack.paste(frame, (0, i * frame.height))
  palette = stack.quantize(
      colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE
  )
  frames = [
      frame.quantize(palette=palette, dither=Image.Dither.NONE)
      for frame in frames
  ]
  frames[0].save(
      path, save_all=True, append_images=frames[1:], duration=duration, loop=0
  )


save_gif(fig, update, n_days, f'{plot_dir}/05_clinical_vs_wastewater_over_time.gif')


# ==========================================
# 7. Transmission Phylogeny
# ==========================================
ay4_day = 30
tree_ay4 = cv.variant(
    label='AY.4', days=ay4_day, n_imports=3, founding_fasta=ay4_fasta
)
outbreak = cv.Sim(
    pop_size=1000,
    pop_type='rank',
    pop_infected=3,
    location='zambia',
    n_days=100,
    rand_seed=3,
    evo_pars=dict(
        enable=True,
        reference=ref_fasta,
        mol_clock_rate=1e-6,
        fitness_model='bloom_nt',
        fitness_data_path=fitness_data_path,
        fitness_scale=0.05,
    ),
    variants=[tree_ay4],
    verbose=0,
)
outbreak.run()
log = outbreak.people.infection_log
n_days = outbreak.npts

infections = collections.defaultdict(list)
for k, entry in enumerate(log):
  infections[int(entry['target'])].append((entry['date'], k))


def active_infection(agent, day):
  return [k for date, k in infections[agent] if date <= day][-1]


parent = [
    (
        None
        if e['source'] is None
        else active_infection(int(e['source']), e['date'])
    )
    for e in log
]
date = np.array([e['date'] for e in log])
n_muts = np.zeros(len(log))
for k, entry in enumerate(log):
  n_muts[k] = (
      n_muts[parent[k]] if parent[k] is not None else 0
  ) + entry['n_mutations']

children = collections.defaultdict(list)
for k, p in enumerate(parent):
  if p is not None:
    children[p].append(k)
roots = [k for k in range(len(log)) if parent[k] is None]

ROOT_LEN = 0.001


def to_newick(day, measure):
  at = lambda k: date[k] if measure == 'days' else n_muts[k]

  def subtree(k, parent_at):
    kids = [c for c in children[k] if date[c] <= day]
    if not kids:
      return f'tip{k}:{at(k) - parent_at}'
    parts = [subtree(c, at(k)) for c in kids] + [f'tip{k}:0']
    return f'({",".join(parts)}):{at(k) - parent_at}'

  return (
      '('
      + ','.join(subtree(r, -ROOT_LEN) for r in roots if date[r] <= day)
      + '):0;'
  )


tip_order = []


def collect_tips(k):
  for c in children[k]:
    collect_tips(c)
  tip_order.append(f'tip{k}')


for r in roots:
  collect_tips(r)
tip_y = {name: y for y, name in enumerate(tip_order)}
tip_variant = {f'tip{k}': e['variant'] for k, e in enumerate(log)}


def branch_colour(branch):
  variants = (
      {tip_variant[name] for name in branch.leaves}
      if branch.is_node()
      else {tip_variant[branch.name]}
  )
  return palette[variants.pop()] if len(variants) == 1 else '0.5'


def draw_tree(ax, day, measure, point_size):
  tree = bt.make_tree(to_newick(day, measure), 'divergence')
  tree.traverse_tree()
  y_cache = {}

  def y_of(branch):
    if id(branch) not in y_cache:
      y_cache[id(branch)] = (
          tip_y[branch.name]
          if branch.is_leaf()
          else np.mean([y_of(c) for c in branch.children])
      )
    return y_cache[id(branch)]

  tree.plot_tree(
      ax,
      colourFxn=branch_colour,
      width=0.5,
      autoSort=False,
      yCoordinateFxn=y_of,
  )
  tree.plot_points(
      ax,
      targetFxn=lambda k: k.is_leaf(),
      colourFxn=branch_colour,
      yCoordinateFxn=y_of,
      pointSize=point_size,
      outline=False,
  )
  ax.set_ylim(-1, len(tip_order))
  ax.set_yticks([])
  for side in ('top', 'right', 'left'):
    ax.spines[side].set_visible(False)
  return tree


legend_handles = [Line2D([], [], color=palette[v], label=v) for v in ('wild', 'AY.4')] + [
    Line2D([], [], color='0.5', label='mixed')
]

fig, (ax_days, ax_muts) = plt.subplots(1, 2, figsize=(14, 12))
for ax, measure, title, xlabel in [
    (ax_days, 'days', 'Time-scaled', 'Simulation day'),
    (ax_muts, 'mutations', 'Mutations', 'Mutations from the seed infection'),
]:
  newick = to_newick(n_days - 1, measure)
  with open(f'{plot_dir}/06_transmission_phylogeny_{measure}.nwk', 'w') as f:
    f.write(newick)
  tree = draw_tree(ax, n_days - 1, measure, point_size=3)
  ax.set_title(f'{title} ({len(tree.get_external())} infected agents)')
  ax.set_xlabel(xlabel)
plt.tight_layout(rect=(0, 0.03, 1, 1))
fig.legend(handles=legend_handles, loc='lower center', ncol=3, frameon=False)
plt.savefig(f'{plot_dir}/06_transmission_phylogeny.pdf')
plt.close()

n_infectious = np.array(
    outbreak.results['variant']['n_infectious_by_variant'].values
)
fig, (ax_tree, ax_counts) = plt.subplots(
    1, 2, figsize=(14, 7), gridspec_kw=dict(width_ratios=[2, 1])
)
for v, name in outbreak['variant_map'].items():
  ax_counts.plot(
      np.arange(n_days), n_infectious[v], color=palette[name], label=name
  )
ax_counts.set_xlim(0, n_days - 1)
ax_counts.set_xlabel('Simulation day')
ax_counts.set_ylabel('Infectious agents')
ax_counts.legend(frameon=False)
ax_counts.axvline(
    ay4_day, color=palette['AY.4'], linestyle=':', linewidth=1.2
)
today = ax_counts.axvline(0, color='0.3', linewidth=1)
for side in ('top', 'right'):
  ax_counts.spines[side].set_visible(False)


def update_tree_anim(day):
  ax_tree.clear()
  draw_tree(ax_tree, day, 'days', point_size=6)
  ax_tree.set_xlim(0, n_days - 1)
  ax_tree.set_xlabel('Simulation day')
  ax_tree.set_title(f'Infected so far: {sum(date <= day)}')
  today.set_xdata([day, day])
  fig.suptitle(f'Day {day}')


save_gif(
    fig,
    update_tree_anim,
    n_days,
    f'{plot_dir}/07_transmission_phylogeny_over_time.gif',
)