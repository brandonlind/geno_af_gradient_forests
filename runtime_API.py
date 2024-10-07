"""API for the MVP runtime project.

Notes
-----
- simple code to help process performance from run time jobs
"""
from pythonimports import *  # http://github.com/brandonlind/pythonimports

import MVP_summary_functions as mvp  # http://github.com/ModelValidationProgram/MVP-offsets

plt.rcParams.update({'font.family' : 'serif', 'mathtext.default': 'regular'})


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
        'num_loci' : ['500', '5000', '10000'],
        'final_la_bin' : ['0.27 < LA ≤ 0.42', '0.42 < LA ≤ 0.58'],
        'source' : ['ind', 'pooled'],
        'model' : ['geno-model', 'af-model']
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
        'pooled' : 'navy',
        'geno-model' : 'cornflowerblue',
        'af-model' : 'navy'
    }
)
boxplot_kwargs = mvp.boxplot_kwargs.copy()

mvp.factor_names['num_loci'] = 'Number of loci'
for factor in [500, 5000, 10000, 20000]:
    mvp.factor_names[str(factor)] = str(factor)
mvp.factor_names['0.27 < LA ≤ 0.42'] = '0.27 < LA ≤ 0.42'
mvp.factor_names['0.42 < LA ≤ 0.58'] = '0.42 < LA ≤ 0.58'
mvp.factor_names['final_la_bin'] = 'Local Adaptation (ΔSA)'
mvp.factor_names['ind'] = '$\it{GO}_{geno, ind}$' # 'Individual-level'
mvp.factor_names['pooled'] = '$\it{GO}_{AF, pop}$' # 'Population-level'
mvp.factor_names['ind-avg'] = '$\it{GO}_{geno, pop}$'
mvp.factor_names['af-ind'] = '$\it{GO}_{AF, ind}$'
mvp.factor_names['geno-model'] = '$\it{GF}$' + '$_{geno}$'
mvp.factor_names['af-model'] = '$\it{GF}$' + '$_{AF}$'
mvp.factor_names['source'] = 'Workflow'
mvp.factor_names['model'] = 'Workflow'
factor_names = mvp.factor_names

hline_kwargs = dict(linestyle='--', color='gainsboro', linewidth=1, zorder=0)

perf_label = "Performance (Kendall's $\\tau$)"

ylim = (0.4, -1)


# functions
add_legend = mvp.add_legend

def get_rt_seeds():
    """Get the seed IDs for first 3 replicates of 2-trait sims."""
    import MVP_10_train_lfmm2_offset as mvp10

    # load params this way to avoid progress bar from mvp.read_params_file (I don't need annotation either)
    params = mvp10.read_params_file('/home/b.lind/offsets/run_20220919_0-225/slimdir')
    
    reps = params.iloc[:675]  # three reps from 225 levels

    seeds = reps.seed[reps.arch.str.contains('2-trait')].astype(str).tolist()  # just the 2-trait levels

    assert len(seeds) == 540
    
    return seeds

seeds = get_rt_seeds()  # add so I can access as rt.seeds

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
    param_cols = ['final_LA', 'glevel', 'plevel', 'pleio', 'slevel', 'landscape', 'popsize', 'migration']
    seed = op.basename(pkl).split("_")[0]
    for param in param_cols:
        score_df[param] = params.loc[seed, param]
    score_df['simulation_level'] = score_df[param_cols[1:]].apply(lambda x: '_'.join(x.astype(str)), axis=1)
    score_df['rep'] = pkl.split('/')[5].split("_")[-1]
    score_df['num_loci'] = pkl.split("/")[6]
    score_df['seed'] = seed
    score_df['demography'] = score_df['popsize'] + '_' + score_df['migration']

    return score_df


@timer
def load_results(source=None, ignore_20k=False):
    """Read in performance data created in ind and pooled notebooks.
    
    Parameters
    ----------
    source : [None, str]
        retrieve performance data from both 'ind' and 'pooled' data (if `source is None`)
            otherwise, retrieve `source` (which is either 'ind' or 'pooled')
    ignore_20k : bool
        whether to remove records from models trained using 20k loci (not all ind-level
            jobs finished, so comparing ind- to pop- would exclude some simulation levels

    Notes
    -----
    used in:
        01_individual_runs/03_gather_individual_scores.ipynb
        02_pooled_runs/03_gather_pooled_scores.ipynb
    """
    combo_cols = ['simulation_level', 'rep', 'num_loci', 'garden']

    if source is None:
        sources = ['ind', 'pooled']
    else:
        assert source in ['ind', 'pooled']
        sources = [source]

    if ignore_20k is True:
        print(ColorText('removing records for models using 20k loci').warn())
    elif ignore_20k is False:
        print(ColorText('keeping records for models using 20k loci').warn())

    dfs = []
    for source in sources:
        df = pd.read_table(f'{resdir}/{source}_performance.txt').astype({'num_loci' : int})
        df['source'] = source
        df['offset_level'] = df[combo_cols].apply(lambda x: '_'.join(x.astype(str)), axis=1)

        if ignore_20k is True:
            df = df[df.num_loci != 20_000]

        print(f'{source} shape = {df.shape}')

        dfs.append(df)

    results = pd.concat(dfs).reset_index(drop=True)

    results['seed_garden'] = results.seed.astype(str) + '_' + results.garden.astype(str)

    results['final_la_bin'] = results.final_LA.apply(lambda la: '0.27 < LA ≤ 0.42' if la <= 0.42 else '0.42 < LA ≤ 0.58')

    results['num_loci'] = results.num_loci.astype(int).astype(str)
    
    results.index = results[['seed', 'garden', 'num_loci']].apply(lambda x: '_'.join(x.astype(str)), axis=1)

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
