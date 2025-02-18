import matplotlib.pyplot as plt

def normal_text():
    # Define Standard Units
    fsize = 18
    tsize = 18
    tdir = 'in'
    major = 5.0
    minor = 3.0
    style = 'default'
    set_params(style, fsize, tsize, tdir, major, minor)


def bigger_text():
    # Define Standard Units
    fsize = 24
    tsize = 24
    tdir = 'in'
    major = 5.0
    minor = 3.0
    style = 'default'
    plt.rcParams['xtick.major.size'] = 5
    set_params(style, fsize, tsize, tdir, major, minor)

    
def large_text():
    # Define Standard Units
    fsize = 28
    tsize = 28
    tdir = 'in'
    major = 5.0
    minor = 3.0
    style = 'default'
    plt.rcParams['xtick.major.size'] = 5
    set_params(style, fsize, tsize, tdir, major, minor)




def set_params(style='default', fsize=18, tsize=18, tdir='in', major=5.0, minor=3.0):
    # Set all parameters for the plot
    plt.style.use(style)
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = fsize
    plt.rcParams['legend.fontsize'] = tsize-2
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'Helvetica'
    plt.rcParams['xtick.direction'] = tdir
    plt.rcParams['ytick.direction'] = tdir
    plt.rcParams['xtick.major.size'] = major
    plt.rcParams['xtick.minor.size'] = minor
    plt.rcParams['ytick.major.size'] = major
    plt.rcParams['ytick.minor.size'] = minor
    plt.rcParams['xtick.labelsize'] = fsize-2
    plt.rcParams['ytick.labelsize'] = fsize-2
    plt.rcParams['font.stretch'] = 'expanded'

    # LaTeX preamble to force sans-serif for numbers and avoid math mode
    plt.rcParams['text.latex.preamble'] = r'''
        \usepackage{helvet}  % Use Helvetica for text
        \renewcommand{\rmdefault}{phv}  % Set default font family to Helvetica
        \renewcommand{\sfdefault}{phv}  % Set sans-serif font
        \usepackage{sfmath}  % Make math font sans-serif
        \everymath={\displaystyle}  % Ensure display math is formatted properly
    '''

    # Additional Customizations
    plt.rcParams['lines.linewidth'] = 3.0
    plt.rcParams['lines.markersize'] = 6
    plt.rcParams['lines.markeredgewidth'] = 3.0
    plt.rcParams['grid.color'] = 'gray'
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = 0.8
    plt.rcParams['grid.alpha'] = 0.7
    plt.rcParams['axes.grid'] = True
    plt.rcParams['figure.figsize'] = [8, 6]
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab10.colors)
    plt.rcParams['figure.autolayout'] = True
    plt.rcParams['axes.titlesize'] = tsize + 2
    plt.rcParams['axes.titleweight'] = 'normal'
    plt.rcParams['axes.labelweight'] = 'normal'
    plt.rcParams['axes.linewidth'] = 2.0 
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.left'] = True
    plt.rcParams['axes.spines.bottom'] = True
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['xtick.minor.width'] = 0.75
    plt.rcParams['ytick.minor.width'] = 0.75
    plt.rcParams['xtick.major.pad'] = 7
    plt.rcParams['ytick.major.pad'] = 7
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['legend.framealpha'] = 0.5  # Semi-transparent background
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['savefig.format'] = 'png'

    #add error bar caps
    plt.rcParams['errorbar.capsize'] = 4
