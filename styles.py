from matplotlib import rc, cm
from cycler import cycler

font   = {"family":"serif", "serif":"Times New Roman", "size":16.0}
text   = {"usetex":True, "latex.preamble":r"\usepackage{helvet},\usepackage{amsmath},\usepackage{sfmath},\renewcommand{\familydefault}{\sfdefault},\boldmath"}

rc("font", **font)
rc("text", **text)
rc("axes", linewidth=0.5, labelsize="small", titlesize="large", prop_cycle=cycler("color", cm.rainbow(range(10))))
rc("xtick.major", width=0.3, pad=2)
rc("xtick", labelsize="x-small")
rc("ytick.major", width=0.3, pad=2)
rc("ytick", labelsize="x-small")
rc("lines", linewidth=1.0, markeredgewidth=0.0, markersize=7)
rc("patch", linewidth=1.0)
rc("legend", numpoints=1, fontsize="xx-small", frameon=False)
rc("savefig", format="pdf", dpi=92, bbox="tight")

