"""API for the MVP runtime project.

Notes
-----
- simple code to help process performance from run time jobs
"""
from pythonimports import *

import MVP_summary_functions as mvp  # http://github.com/ModelValidationProgram/MVP-offsets

plt.rcParams.update({'font.family' : 'serif'})


# navigation
resdir = '/work/lotterhos/brandon/runtime'

dirs = {
    'ind' : {
        'run_20220919_0-225': [
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_0-225/00500',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_0-225/05000',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_0-225/10000',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_0-225/20000'
        ],
        'run_20220919_225-450': [
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_225-450/00500',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_225-450/05000',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_225-450/10000',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_225-450/20000'
        ],
        'run_20220919_450-675': [
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_450-675/00500',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_450-675/05000',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_450-675/10000',
            '/work/lotterhos/brandon/ind_runtimes/run_20220919_450-675/20000'
        ]
    },
    
    'pooled' : {
        'run_20220919_0-225': [
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_0-225/00500',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_0-225/05000',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_0-225/10000',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_0-225/20000'
        ],
        'run_20220919_225-450': [
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_225-450/00500',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_225-450/05000',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_225-450/10000',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_225-450/20000'
        ],
        'run_20220919_450-675': [
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_450-675/00500',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_450-675/05000',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_450-675/10000',
            '/work/lotterhos/brandon/pooled_runtimes/run_20220919_450-675/20000'
        ]
    }
}

# figure making API - update from MVP Offset project
mvp.hue_order.update(
    {
        'num_loci' : ['500', '5000', '10000', '20000'],
        'final_la_bin' : ['0.27 < LA ≤ 0.42', '0.42 < LA ≤ 0.58'],
        'source' : ['ind', 'pooled']
    }
)
hue_order = mvp.hue_order.copy()

mvp.boxplot_kwargs['palette'].update(
    {
        # sns.cubehelix_palette(start=.5, rot=-.75, n_colors=7)
        '500' : (0.8423298817793848, 0.8737404427964184, 0.7524954030731037),  # 0
        '5000' : (0.6486603420129703, 0.797769603957374, 0.6186454636675479),  # 1
        '10000' : (0.32562863725703667, 0.5824294714811111, 0.551260440725878),  # 3
        '20000': (0.22630233856123372, 0.25904677946860183, 0.4176219861426439),  # -2
        '0.27 < LA ≤ 0.42': (0.8892638312853967, 0.8490264305563623, 0.7570511784894085),
        '0.42 < LA ≤ 0.58': (0.6779472567428826, 0.4089021118923688, 0.5211323732841375),
        'ind' : 'cornflowerblue',
        'pooled' : 'navy'
    }
)
boxplot_kwargs = mvp.boxplot_kwargs.copy()

mvp.factor_names['num_loci'] = 'Number of loci'
for factor in [500, 5000, 10000, 20000]:
#     mvp.factor_names[factor] = factor
    mvp.factor_names[str(factor)] = str(factor)
mvp.factor_names['0.27 < LA ≤ 0.42'] = '0.27 < LA ≤ 0.42'
mvp.factor_names['0.42 < LA ≤ 0.58'] = '0.42 < LA ≤ 0.58'
mvp.factor_names['final_la_bin'] = 'Local Adaptation (ΔSA)'
mvp.factor_names['ind'] = 'Individual-level'
mvp.factor_names['pooled'] = 'Population-level'
mvp.factor_names['source'] = 'Genetic Source'
factor_names = mvp.factor_names

hline_kwargs = dict(linestyle='--', color='gainsboro', linewidth=1, zorder=0)

perf_label = "Performance (Kendall's $\\tau$)"

ylim = (0.4, -1)


# functions
add_legend = mvp.add_legend

def process_performance_pkl(pkl, source=None):
    """Format saved performance pkl.
    
    Notes
    -----
    - only used within:
        - 02_pooled_runs/03_gather_pooled_scores.ipynb
        - 01_individual_runs/03_gather_individual_scores.ipynb
    """
    from pythonimports import pklload
    import pandas as pd
    from os import path as op
    import MVP_summary_functions as mvp
    
    # read in performance scores
    scores = pklload(pkl)['garden_performance'][source]['all']
    
    # reformat scores to dataframe
    scores.name = 'score'
    score_df = pd.DataFrame(scores).reset_index()
    score_df.columns = score_df.columns.str.replace('index', 'garden')
    
    # add metadata
    seed = op.basename(pkl).split("_")[0]
    for param in params.columns[-9:]:
        score_df[param] = params.loc[seed, param]
    score_df['simulation_level'] = score_df[score_df.columns[3:]].apply(lambda x: '_'.join(x.astype(str)), axis=1)
    score_df['rep'] = pkl.split('/')[5].split("_")[-1]
    score_df['num_loci'] = pkl.split("/")[6]
    score_df['seed'] = seed
    score_df['demography'] = score_df['popsize'] + '_' + score_df['migration']

    return score_df


@timer
def load_results():
    """Read in performance data created in ind and pooled notebooks.
    
    Notes
    -----
    used in:
        01_individual_runs/03_gather_individual_scores.ipynb
        02_pooled_runs/03_gather_pooled_scores.ipynb
    """
    combo_cols = ['simulation_level', 'rep', 'num_loci', 'garden']
    
    df_ind = pd.read_table(f'{resdir}/ind_performance.txt')
    df_ind['source'] = 'ind'
    df_ind['offset_level'] = df_ind[combo_cols].apply(lambda x: '_'.join(x.astype(str)), axis=1)
    
    df_pooled = pd.read_table(f'{resdir}/pooled_performance.txt')
    df_pooled['source'] = 'pooled'
    df_pooled['offset_level'] = df_pooled[combo_cols].apply(lambda x: '_'.join(x.astype(str)), axis=1)
    
    print('ind shape = ', df_ind.shape)
    print('pooled shaped = ', df_pooled.shape)
    
    results = pd.concat([df_ind, df_pooled])
    
    results['seed_garden'] = results.seed.astype(str) + '_' + results.garden.astype(str)
    
    results['final_la_bin'] = results.final_LA.apply(lambda x: '0.27 < LA ≤ 0.42' if x <= 0.42 else '0.42 < LA ≤ 0.58')
    
    return results


def latest_commit():
    """Print today's date, author info, and commit hashes of pythonimports and MVP_offsets."""
    import pythonimports as pyimp

    pyimp_info = pyimp._git_pretty(pyimp._find_pythonimports())
    mvp_info = pyimp._git_pretty('/home/b.lind/code/MVP-offsets')
    rt_info = pyimp._git_pretty('/work/lotterhos/brandon/code/06_run_time_project')
    current_datetime = "Today:\t" + time.strftime("%B %d, %Y - %H:%M:%S %Z") + "\n"
    version = "python version: " + sys.version.split()[0] + "\n"

    width = max([len(x) for x in flatten([pyimp_info.split('\n'), mvp_info.split('\n'), current_datetime])])
    hashes = '#########################################################\n'
    
    try:
        env = 'conda env: %s\n' % os.environ['CONDA_DEFAULT_ENV']
    except KeyError as e:
        env = ''

    print(
        hashes
        + current_datetime
        + version + f'{env}\n'
        + f"Current commit of %s:\n" % ColorText("pythonimports").bold()
        + pyimp_info + '\n'
        + "Current commit of %s:\n" % ColorText('MVP_offsets').bold().blue()
        + mvp_info + '\n'
        + "Current commit of %s:\n" % ColorText('MVP_runtime').bold().custom('purple')
        + rt_info
        + hashes
    )
    
    pass
