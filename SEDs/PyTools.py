from numpy import ones, array, digitize, average, lexsort, linspace, median, mean, std, argmin, interp, trapz
from scipy.stats.mstats import scoreatpercentile, mode
from struct import calcsize, unpack
from fortranfile import FortranFile

def nbins(sample, range_ = None) :
  sample_ = array(sample)
  IQR = lambda x    : scoreatpercentile(x, 75.0) - scoreatpercentile(x, 25.0)
  if range_ is None : mn, mx = sample_.min(), sample_.max()
  else              : mn, mx = range_

  mask    = (sample_ >= mn) & (sample_ <= mx)
  binsize = (2 * IQR(sample_[mask]) / mask.sum() ** (1. / 3))

  return (mx - mn) / binsize, mn, mx, binsize

def percent_labels(ax, ticks, labels):
  x = ax.get_xlim()
  y = ax.get_ylim()
  ax1 = ax.twinx()
  ax2 = ax.twiny()
  ax1.set_ylim(y)
  ax2.set_xlim(x)
  ax1.set_yticks(ticks)
  ax2.set_xticks(ticks)
  ax1.set_yticklabels([r"-" + l for l in labels])
  ax2.set_xticklabels(labels)

def count_nearest(vec, v_i):
  return argmin(abs(vec-v_i))

def err_slope(percent):
  n_percent = percent*10/100.
  m1 = (10 + n_percent)/(10 - n_percent)
  m2 = 1./m1
  return m1, m2

def selection(chi1, chi2, par1, par2):
  """
  Compute the best parameter using minimum chi square criteria
  """
  result = [None]*len(par1)
  for i in range(len(par1)):
    if chi1[i] <= chi2[i]: result[i] = par1[i]
    else: result[i] = par2[i]
  return array(result)

def binner(x, y, w_sta, nbins = None, rang = None, ebar = False, per = None):
	'''
	This function reduces (x, y) scatter by binning in the x range and computing a given statistic for
	y values in each bin.
	
	input parameters:
	----------------
	*  x, y  : arraylikes containing the data to be reduced.
	*  w_sta : which statistic is to be computed in order to reduce the data. Allowed values are:
						 'mean', 'median' and 'mode'.
	*  nbins : number of bins to be computed in x. Defaults to 10.
	*  rang  : range where the data will be binned. Defaults to (x.min(), x.max())
	*  ebar  : True to return a typical deviation from the central value in each bin. Defaults to False.
	*  per   : a tuple with de percentiles to be computed in order to return deviations from the
	           central value in each bin. Ignored if ebar = False.
	           
	output parameters:
	-----------------
	*  x_cen : an array containing the center of each x bin.
	*  y_sta : an array containing the central value of y in each bin.
  * [yer]  : if ebar = True, an array of shape = (2, nbins) containing the typical deviation from the
             central value in each bin.
	'''
	ind    = lexsort((y, x))
	xs, ys = x[ind], y[ind]

	if rang is None : mn, mx = min(xs), max(xs)
	else            : mn, mx = rang

	if nbins is None : nbins = 10
    #samp  = (xs >= mn) & (xs <= mx)
    #n     = samp.sum()
    #IQR   = score(xs[samp], 75.0) - score(xs[samp], 25.0)
    #nbins = (mx - mn) / (2 * IQR * n ** (1. / 3))  # Freedman-Diaconis rule for selecting bin size

	bins  = linspace(mn, mx, nbins + 1)
	x_cen = (bins[: - 1] + bins[1:])*0.5
	bins  = linspace(mn, mx, nbins)
	ibins = digitize(xs, bins)

	if w_sta   == "median" : y_sta = array([median(ys[ibins == i]) for i in range(1, bins.size + 1)])
	elif w_sta == "mean"   : y_sta = array([average(ys[ibins == i]) for i in range(1, bins.size + 1)])
	elif w_sta == "mode"   : y_sta = array([mode(ys[ibins == i])[0] for i in range(1, bins.size + 1)])

	if ebar   == False                : return x_cen, y_sta
	elif ebar == True and per == None :
		myer = abs(array([scoreatpercentile(ys[ibins == i], 16.0) for i in range(1, bins.size + 1)]) - y_sta)
		pyer = abs(array([scoreatpercentile(ys[ibins == i], 84.0) for i in range(1, bins.size + 1)]) - y_sta)
		yer  = array([myer, pyer])
		return x_cen, y_sta, yer

	elif ebar == True and per != None :
		myer = abs(array([scoreatpercentile(ys[ibins == i], per[0]) for i in range(1, bins.size + 1)]) - y_sta)
		pyer = abs(array([scoreatpercentile(ys[ibins == i], per[1]) for i in range(1, bins.size + 1)]) - y_sta)
		yer = array([myer, pyer])
		return x_cen, y_sta, yer

def err(xr, x, relative = True) :
  '''
  Compute relative error between theoretical quantity x and real xr.
  '''

  if relative     : return (x - xr)/xr
  if not relative : return x - xr

def vec_rep(vec, nrep):
  '''
  Repeats every element in vec nrep times, preserving same order. Resulting vector will have
  size = vec.size*nrep
  '''

  output = ones(nrep*vec.size)
  for i in range(vec.size): output[i*nrep:(i + 1)*nrep] *= vec[i]
  return output

def read_mybin(fname, to_return = ["fluxes"], kind = "theo", which_sed = []):
  '''
  This function reads a Fortran binary file written by programs
  LIB_CONSTRUCTOR, THEO_PHOTOMETRY and THEO_SPECTROSCOPY.

  Variables
  ---------

  fname:      Name of the file to be readed.

  to_return:  A list of the columns that suppose to be in the file and
          you want to get. Available options are:

          * lambda
          * fluxes
          * sigmas
          * sfr       (t_form, A_burts, t_burst, tau)
          * physical  (mass, <log(t)>_m, <log(t)>_fr,
                      10**(<log(t)>_m), 10**(<log(t)>_fr))
          * indexes   (B-V, g-r, H_delta, Dn(4000))

          If a item in  to_return is no in file, then this will be
          ignored. By default to_return = ["fluxes"]

  kind:       Format of the file to be readed. There are two options:

          * theo: Only  lambda and fluxes  could be returned (even
                  if you  insist in  request  other columns).
          * empi: All parameters  availables for  to_return can be
                  found in this type of files.

          By default kind = "theo"

  Return
  ------

  A dictionary  which keys are  to_return  depending on  the kind, and
  values are the lists of the columns found in the file.
  '''

  i = calcsize("i")
  r = calcsize("f")

  IO_msg = "Error reading record from data file"

  wlen   = []
  fluxes = []
  sigmas = []
  sfr        = []
  physic = []
  indexs = []

  returnable = {
                "lambda"   : wlen,
                "fluxes"   : fluxes,
                "sigmas"   : sigmas,
                "sfr"      : sfr,
                "physical" : physic,
                "indexes"  : indexs
               }

  rkeys = returnable.viewkeys()
  if rkeys & set(to_return) == set(): return None

  f = FortranFile(fname)

  reg_size = unpack("i",f._read_exactly(i))[0]
  size     = unpack("i",f._read_exactly(i))[0]

  wlen.append(unpack("f"*size,f._read_exactly(r*size)))
  if reg_size != f._read_check(): raise IOError(IO_msg)

  if kind == "theo":
    fmt = "i" + "f"*size
    count = 1
    while True:
            try:
                    if which_sed != []:
                            if count in which_sed: fluxes.append(unpack(fmt,f.readRecord())[1:])
                            else: f.readRecord()
                            count += 1
                    else: fluxes.append(unpack(fmt,f.readRecord())[1:])
            except struct.error:
                    name = fname.split("/")[-1]
                    exit("ERROR: Check if 'kind' has the value corresponding to file '"+name+"'.")
            except IOError: break

    dict_to_return = dict()
    for key in to_return:
            values = returnable[key]
            if values != []: dict_to_return.setdefault(key,values)

    return dict_to_return

  elif kind == "empi":
    nsfr = 2*size + 4
    nphy = nsfr + 5
    ncin = nphy + 4

    fmt   = "i" + "f"*size*2 + "f"*4 + "f"*5 + "f"*4
    matr  = []
    count = 1
    while True:
      try:
        if which_sed != []:
          if count in which_sed: matr.append(unpack(fmt,f.readRecord())[1:])
          else: f.readRecord()
          count += 1
        else: matr.append(unpack(fmt,f.readRecord())[1:])
      except struct.error:
        name = fname.split("/")[-1]
        exit("ERROR: Check if 'kind' has the value corresponding to file '"+name+"'.")
      except IOError: break

    matr = array(matr)
    fluxes.append(matr[:,:size])
    sigmas.append(matr[:,size:2*size])
    sfr.append(matr[:,2*size:nsfr])
    physic.append(matr[:,nsfr:nphy])
    indexs.append(matr[:,nphy:ncin])

    dict_to_return = dict()
    for key in to_return:
      values = returnable[key]
      if values != []: dict_to_return.setdefault(key,values)

  return dict_to_return

def binned_stat(table, bin_size, stat = "mean") :
	'''
	This function computes the given statistic in a table (row sense)
	by bins of a given size.
	
	input parameters:
	----------------
	*  table    : an arraylike table to be binned_averaged.
	*  bin_size : the number of elements in each bin. It's mandatory
	              that bin_size * int(# of bins) = table.shape[0].
  *  stat     : character string with the name of the statistic to compute. Available options are:
                - mean (default)
                - median
                - stdev
	'''
	N_bins = table.shape[0] / bin_size
	if N_bins * bin_size == table.shape[0] :
		if stat == "mean" :
			return array([mean(table[i*bin_size:(i+1)*bin_size], 0) for i in xrange(N_bins)])
		elif stat == "median" :
			return array([median(table[i*bin_size:(i+1)*bin_size], 0) for i in xrange(N_bins)])
		elif stat == "stdev" :
			return array([std(table[i*bin_size:(i+1)*bin_size], 0) for i in xrange(N_bins)])
	else :
		raise ValueError, "bin_size must fulfill: bin_size * int(# of bins) = table.shape[0]."

def flux(ABmag, weff) :
  return 10 ** (-0.4 * (ABmag - 22.5)) * 3.631e-6 * 1e-23 * (3e18 / weff ** 2)

def integrated_flux(SED, passband):
  mask = (passband[0, 0] <= SED[:, 0])&(SED[:, 0] <= passband[-1, 0])
  
  ipassband = interp(SED[mask, 0], passband[:, 0], passband[:, 1])

  return trapz(SED[mask, 0]*SED[mask, 1]*ipassband, SED[mask, 0])/trapz(SED[mask, 0]*ipassband, SED[mask, 0])

"""reduccion de los datos con HISTOGRAM2D"""
#H, xe, ye = histogram2d(x, y, bins = 100, weights = etm, normed = 1)
#H_m = zeros(H.shape)
#for xx in range(xe.size - 1):
        #for yy in range(ye.size - 1):
                #if xe[xx] != xe[-1]: xbin = (x >= xe[xx])&(x < xe[xx + 1])
                #else: xbin = (x >= xe[xx])&(x <= xe[xx + 1])
                #if ye[yy] != ye[-1]: ybin = (y >= ye[yy])&(y < ye[yy + 1])
                #else: ybin = (y >= ye[yy])&(y <= ye[yy + 1])
                #H_m[xx, yy] = np.median(etm[xbin&ybin])
#ex = [ye[0], ye[-1], xe[-1], xe[0]]
#hb = imshow(H_m.T, extent = ex, origin = "lower", interpolation = "nearest", aspect = "auto",
                                                #cmap = cm.bwr, vmin = -0.5, vmax = 0.5)
