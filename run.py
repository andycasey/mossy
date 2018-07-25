import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from glob import glob

from pyraf import iraf

from  utils  import (fit_continuum, radial_velocity_correction, cross_correlate,
                     load_gemini_spectrum, _get_closest)

iraf.fitsutil()
iraf.gemini()
iraf.gmos()

# For cleaning purposes.
before = glob("*")

association = Table.read("summary.fits")


flat_kwds = dict(fl_vardq=True, fl_fulldq=True, fl_detec=True, 
                 fl_oversize=False, fl_inter=False, order="13,11,28")

gswavelength_kwds = dict(fwidth=11, cradius=5, minsep=5, order=9, match=-5,
                         fitcxord=4, fitcyord=3, thresh=1000, nsum=3, 
                         step=1, fl_inter="no", trace="yes", 
                         fl_addfeat=False, coordlist="gmos$data/CuAr_GMOS.dat")

extract_kwds = dict(
    fl_inter=False, find=True, back="fit", bfunc="chebyshev", border=1,
    tfunct="spline3", torder=5, tnsum=20, tstep=50, refimage="",
    apwidth=1.3, recent=True, trace=True, fl_vardq=True, Stdout=1,
    weights="variance")

normalization_kwds = dict(function="spline",
                          high_sigma_clip=0.5,
                          knot_spacing=100.0,
                          low_sigma_clip=3.0,
                          max_iterations=5,
                          order=5)


# Create master bias.
biases = association[association["obstype"] == "BIAS"]
with open("bias.list", "w") as fp:
    fp.write("\n".join(["{}".format(os.path.basename(b)) \
                        for b in biases["path"]]))

iraf.unlearn("gireduce")
iraf.unlearn("gbias")
iraf.unlearn("gemextn")

iraf.gbias("@bias.list", "master_bias.fits", 
           log="bias.log", fl_vardq=True, verbose="yes")


# Create master flat.
flats = association[association["obstype"] == "FLAT"]
with open("flat.list", "w") as fp:
    fp.write("\n".join(["{}".format(os.path.basename(b)) \
                        for b in flats["path"]]))

iraf.unlearn("gireduce")
iraf.unlearn("gflat")
iraf.unlearn("gemextn")

iraf.gsflat("@flat.list",
            "master_flat.fits",
            bias="master_bias.fits",
            logfile="flat.log",
            **flat_kwds)

# Reduce arcs.
arcs = association[association["obstype"] == "ARC"]

for arc in arcs:

    path = os.path.basename(arc["path"])

    iraf.unlearn("gireduce")
    iraf.unlearn("gsflat")
    iraf.unlearn("gemextn")   
    
    iraf.gsreduce(path, bias="master_bias.fits",
                  fl_fixpix=False, fl_flat=False, fl_oversize=False)

    # Do wavelength calibration on arc frames.
    iraf.unlearn("gswavelength")
    iraf.gswavelength("gs{}".format(path), **gswavelength_kwds)



# Reduce all science frames.
science_frames = association[(association["obsclass"] == "science") \
    + ((association["obstype"] == "OBJECT") * (association["obsclass"] == "partnerCal"))]

# Load template spectrum.
template_disp, template_flux, _ = np.loadtxt("Atlas.Arcturus.372_926nm.txt",
                                             skiprows=1).T
template_disp *= 10.0

for frame in science_frames:

    path = os.path.basename(frame["path"])

    iraf.unlearn("gireduce")
    iraf.unlearn("gsreduce")
    iraf.unlearn("gemextn")

    # Reduce.
    iraf.gsreduce(path, 
                  bias="master_bias.fits", flatim="master_flat.fits",
                  fl_fixpix=False, fl_oversize=False, fl_vardq=True, 
                  fl_fulldq=True)#fl_gscrrej=True)

    # Get the closest arc.
    arc = _get_closest(frame, arcs, match_keys=("grating", ))
    arc_path = "gs{}".format(os.path.basename(arc["path"]))

    iraf.unlearn("gstransform")
    iraf.gstransform("gs{}".format(path),
                     wavtraname=arc_path,
                     fl_vardq=True)

    # Sky subtraction.
    iraf.unlearn("gsskysub")
    iraf.gsskysub("tgs{}".format(path), fl_vardq=True, fl_inter=False)

    # Extraction.
    iraf.unlearn("gsextract")
    iraf.gsextract("stgs{}".format(path), **extract_kwds)

    # Load and perform continuum normalisation.
    disp, flux, ivar, dq, meta = load_gemini_spectrum("estgs{}".format(path))

    # split into channels and normalize.
    split_points = [4460, 5520]
    split_indices = np.hstack([
        0, 
        np.repeat(disp.searchsorted(split_points), 2), 
        disp.size
    ]).reshape(-1, 2)

    continuum = np.nan * np.ones(disp.size)
    for si, ei in split_indices:
        try:
            continuum[si:ei] = fit_continuum(disp[si:ei], flux[si:ei], 
                                             ivar[si:ei], **normalization_kwds)
        except:
            continue

    # Measure RV with some template?
    try:
        v_rad, ccf = cross_correlate(disp[si:ei], (flux/continuum)[si:ei],
                                     template_disp, template_flux)
    except:
        v_rad = np.nan

    # Barycentric correction.
    v_bary = radial_velocity_correction(**meta)

    image = fits.open("estgs{}".format(path))
    image[0].header["V_BARY"] = v_bary.value
    image[0].header["V_RAD"] = v_rad.value

    hdu = fits.ImageHDU()
    hdu.data = continuum
    hdu.header["EXTNAME"] = "CONTINUUM"
    image.append(hdu)

    image.writeto("final_{}".format(path), overwrite=True)


# Clean up intermediate products.
after = glob("*.fits")
intermediate_paths = [ea for ea in list(set(after).difference(before)) \
                      if not ea.startswith(("estgs", "final_"))]

for path in intermediate_paths:
    os.unlink(path)

