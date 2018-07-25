import numpy as np
import logging
from scipy import (interpolate, polyfit, optimize as op)
from astropy import units as u
from astropy.io import fits

from astropy.coordinates  import SkyCoord, EarthLocation
from astropy.time import Time


def fit_continuum(dispersion, flux, ivar, knot_spacing=200, low_sigma_clip=1.0, 
    high_sigma_clip=0.2, max_iterations=3, order=3, exclude=None, include=None, 
    additional_points=None, function='spline', scale=1.0, **kwargs):
    """
    Fits the continuum for a given `Spectrum1D` spectrum.
    
    Parameters
    ----
    knot_spacing : float or None, optional
        The knot spacing for the continuum spline function in Angstroms. Optional.
        If not provided then the knot spacing will be determined automatically.
    
    sigma_clip : a tuple of two floats, optional
        This is the lower and upper sigma clipping level respectively. Optional.
        
    max_iterations : int, optional
        Maximum number of spline-fitting operations.
        
    order : int, optional
        The order of the spline function to fit.
        
    exclude : list of tuple-types containing floats, optional
        A list of wavelength regions to always exclude when determining the
        continuum. Example:
        
        >> exclude = [
        >>    (3890.0, 4110.0),
        >>    (4310.0, 4340.0)
        >>  ]
        
        In the example above the regions between 3890 A and 4110 A, as well as
        4310 A to 4340 A will always be excluded when determining the continuum
        regions.

    function: only 'spline' or 'poly'

    scale : float
        A scaling factor to apply to the normalised flux levels.
        
    include : list of tuple-types containing floats, optional
        A list of wavelength regions to always include when determining the
        continuum.

    :param rv:
        A radial velocity correction (in km/s) to apply to the spectrum.
    """


    exclusions = []
    continuum_indices = range(len(flux))

    # Snip left and right
    finite_positive_flux = np.isfinite(flux) * (ivar > 0)

    if np.sum(finite_positive_flux) == 0:
        # No valid continuum points, return nans
        no_continuum = np.nan * np.ones_like(dispersion)

        raise a

    function = str(function).lower()
    left = np.where(finite_positive_flux)[0][0]
    right = np.where(finite_positive_flux)[0][-1]

    # See if there are any regions we need to exclude
    if exclude is not None and len(exclude) > 0:
        exclude_indices = []
        
        if isinstance(exclude[0], float) and len(exclude) == 2:
            # Only two floats given, so we only have one region to exclude
            exclude_indices.extend(range(*np.searchsorted(dispersion, exclude)))
            
        else:
            # Multiple regions provided
            for exclude_region in exclude:
                exclude_indices.extend(
                    range(*np.searchsorted(dispersion, exclude_region)))
    
        continuum_indices = np.sort(list(set(continuum_indices).difference(
            np.sort(exclude_indices))))
        
    # See if there are any regions we should always include
    if include is not None and len(include) > 0:
        include_indices = []
        
        if isinstance(include[0], float) and len(include) == 2:
            # Only two floats given, so we can only have one region to include
            include_indices.extend(range(*np.searchsorted(dispersion, include)))
            
        else:
            # Multiple regions provided
            for include_region in include:
                include_indices.extend(range(*np.searchsorted(
                    dispersion, include_region)))
    
    # We should exclude non-finite values from the fit.
    non_finite_indices = np.where(~np.isfinite(flux * ivar))[0]
    continuum_indices = np.sort(list(set(continuum_indices).difference(
        non_finite_indices)))

    # We should also exclude zero or negative flux points from the fit
    zero_flux_indices = np.where(0 >= flux)[0]
    continuum_indices = np.sort(list(set(continuum_indices).difference(
        zero_flux_indices)))

    if 1 > continuum_indices.size:
        no_continuum = np.nan * np.ones_like(dispersion)
        raise a


    original_continuum_indices = continuum_indices.copy()

    if knot_spacing is None or knot_spacing == 0:
        knots = []

    else:
        knot_spacing = abs(knot_spacing)
        
        end_spacing = ((dispersion[-1] - dispersion[0]) % knot_spacing) /2.
        if knot_spacing/2. > end_spacing: end_spacing += knot_spacing/2.
            
        knots = np.arange(dispersion[0] + end_spacing, 
            dispersion[-1] - end_spacing + knot_spacing, 
            knot_spacing)

        if len(knots) > 0 and knots[-1] > dispersion[continuum_indices][-1]:
            knots = knots[:knots.searchsorted(dispersion[continuum_indices][-1])]
            
        if len(knots) > 0 and knots[0] < dispersion[continuum_indices][0]:
            knots = knots[knots.searchsorted(dispersion[continuum_indices][0]):]

    # TODO: Use inverse variance array when fitting polynomial/spline.
    for iteration in range(max_iterations):
        
        if 1 > continuum_indices.size:

            no_continuum = np.nan * np.ones_like(dispersion)
            raise a
            """            
            failed_spectrum = self.__class__(dispersion=dispersion,
                flux=no_continuum, ivar=no_continuum, metadata=self.metadata)

            if kwargs.get("full_output", False):
                return (failed_spectrum, no_continuum, 0, dispersion.size - 1)

            return failed_spectrum
            """

        splrep_disp = dispersion[continuum_indices]
        splrep_flux = flux[continuum_indices]
        splrep_weights = ivar[continuum_indices]**0.5

        median_weight = np.nanmedian(splrep_weights)

        # We need to add in additional points at the last minute here
        if additional_points is not None and len(additional_points) > 0:

            for point, flux, weight in additional_points:

                # Get the index of the fit
                insert_index = int(np.searchsorted(splrep_disp, point))
                
                # Insert the values
                splrep_disp = np.insert(splrep_disp, insert_index, point)
                splrep_flux = np.insert(splrep_flux, insert_index, flux)
                splrep_weights = np.insert(splrep_weights, insert_index, 
                    median_weight * weight)

        if function == 'spline':
            if order > 5:
                logging.warn("Splines can only have a maximum order of 5. "
                             "Limiting order value to 5.")
                order = 5

            try:
                tck = interpolate.splrep(splrep_disp, splrep_flux,
                    k=order, task=-1, t=knots, w=splrep_weights)

            except:
                logging.exception("Exception in fitting continuum:")
                continuum = np.nan * np.ones_like(dispersion)

            else:
                continuum = interpolate.splev(dispersion, tck)

        elif function in ("poly", "polynomial"):

            coeffs = polyfit(splrep_disp, splrep_flux, order)

            popt, pcov = op.curve_fit(lambda x, *c: np.polyval(c, x), 
                splrep_disp, splrep_flux, coeffs, 
                sigma=ivar[continuum_indices], absolute_sigma=False)
            continuum = np.polyval(popt, dispersion)

        else:
            raise ValueError("Unknown function type: only spline or poly "\
                "available ({} given)".format(function))
        
        difference = continuum - flux
        sigma_difference = difference / np.std(difference[np.isfinite(flux)])

        # Clipping
        upper_exclude = np.where(sigma_difference > high_sigma_clip)[0]
        lower_exclude = np.where(sigma_difference < -low_sigma_clip)[0]
        
        exclude_indices = list(upper_exclude)
        exclude_indices.extend(lower_exclude)
        exclude_indices = np.array(exclude_indices)
        
        if len(exclude_indices) is 0: break
        
        exclusions.extend(exclude_indices)
        
        # Before excluding anything, we must check to see if there are regions
        # which we should never exclude
        if include is not None:
            exclude_indices \
                = set(exclude_indices).difference(include_indices)
        
        # Remove regions that have been excluded
        continuum_indices = np.sort(list(set(continuum_indices).difference(
            exclude_indices)))
    
    # Snip the edges based on exclude regions
    if exclude is not None and len(exclude) > 0:

        # If there are exclusion regions that extend past the left/right,
        # then we will need to adjust left/right accordingly

        left = np.max([left, np.min(original_continuum_indices)])
        right = np.min([right, np.max(original_continuum_indices)])
    
    # Apply flux scaling
    continuum *= scale

    continuum[:left] = np.nan
    continuum[right:] = np.nan

    if kwargs.get("full_output", False):
        return (continuum, left, right)
    return continuum





def cross_correlate(obs_disp, obs_flux, template_disp, template_flux, apodize=0,
    dispersion_range=None, resample="template"):
    """
    Cross-correlate the observed spectrum against a rest-frame template spectrum
    and measure the radial velocity of the source.
    """

    if dispersion_range is None:
        # Use the common ranges.
        dispersion_range = (
            np.max([
                obs_disp[0],
                template_disp[0]
            ]),
            np.min([
                obs_disp[-1],
                template_disp[-1]
            ])
        )

    if not isinstance(dispersion_range, (tuple, list, np.ndarray)) \
    or len(dispersion_range) != 2:
        raise TypeError("wavelength region must be a two length list-type")

    if apodize != 0:
        raise NotImplementedError("apodization not implemented yet")
        
    resample = resample.lower()
    if resample == "template":
        idx = np.searchsorted(obs_disp, dispersion_range)
        finite = np.isfinite(obs_flux[idx[0]:idx[1]])

        dispersion = obs_disp[idx[0]:idx[1]][finite]
        observed_flux = obs_flux[idx[0]:idx[1]][finite]

        func = interpolate.interp1d(
            template_disp, template_flux,
            bounds_error=False, fill_value=0.0)
        template_flux = func(dispersion)

    elif resample == "observed":
        raise NotImplementedError("why would you do this?")

    else:
        raise ValueError("resample must be 'template' or 'observed'")


    # Perform the cross-correlation
    padding = observed_flux.size + template_flux.size
    # Is this necessary?: # TODO
    x_norm = observed_flux - np.mean(observed_flux[np.isfinite(observed_flux)])
    y_norm = template_flux - np.mean(template_flux[np.isfinite(template_flux)])

    Fx = np.fft.fft(x_norm, padding, )
    Fy = np.fft.fft(y_norm, padding, )
    iFxy = np.fft.ifft(Fx.conj() * Fy).real
    varxy = np.sqrt(np.inner(x_norm, x_norm) * np.inner(y_norm, y_norm))

    fft_result = iFxy/varxy

    # Put around symmetry axis.
    num = len(fft_result) - 1 if len(fft_result) % 2 else len(fft_result)

    fft_y = np.zeros(num)
    fft_y[:num//2] = fft_result[num//2:num]
    fft_y[num//2:] = fft_result[:num//2]

    fft_x = np.arange(num) - num/2

    # Get initial guess of peak.
    p0 = np.array([fft_x[np.argmax(fft_y)], np.max(fft_y), 10])

    gaussian = lambda p, x: p[1] * np.exp(-(x - p[0])**2 / (2.0 * p[2]**2))
    errfunc = lambda p, x, y: y - gaussian(p, x)

    try:
        p1, ier = op.leastsq(errfunc, p0.copy(), args=(fft_x, fft_y))

    except:
        raise


    # Create functions for interpolating back onto the dispersion map
    fft_points = (0, p1[0])
    interp_x = np.arange(num/2) - num/4

    wl_points = []
    for point in fft_points:
        idx = np.searchsorted(interp_x, point)
        f = interpolate.interp1d(interp_x[idx-3:idx+3], dispersion[idx-3:idx+3],
            bounds_error=True, kind='cubic')
        wl_points.append(f(point))

    # Calculate velocity 
    c = 299792458e-3 # km/s
    f, g = wl_points
    rv = c * (1 - g/f) 

    # Create a CCF spectrum.
    ccf = np.array([fft_x * (rv/p1[0]), fft_y])

    return (rv * u.km/u.s, ccf)



def load_gemini_spectrum(path):
    with fits.open(path) as image:
        dispersion = image[2].header["CRVAL1"] \
                   + (np.arange(1, 1 + image[2].header["NAXIS1"]) - image[2].header["CRPIX1"]) \
                   * image[2].header["CD1_1"]
                
        flux = image[2].data.flatten()
        ivar = 1.0/image[3].data.flatten()
        dq = image[4].data.flatten()

        ivar[dq == 1] = 0

        metadata = dict(ra=image[0].header["RA"],
                        dec=image[0].header["DEC"],
                        obs_time="{} {}".format(image[0].header["DATE-OBS"], 
                                                image[0].header["TIME-OBS"]),
                        telescope=image[0].header["TELESCOP"],
                        path=path)

    return (dispersion, flux, ivar, dq, metadata)


def _get_closest(needle, haystack, mask=None, match_keys=None):

    N = len(haystack)

    if mask is None:
        mask = np.ones(N, dtype=bool)

    if match_keys is not None:
        for key in match_keys:
            mask *= (haystack[key] == needle[key])

    idx = np.arange(N)
    dt = np.abs(needle["obsepoch"] - haystack["obsepoch"])
    
    return haystack[idx[mask][np.argmin(dt[mask])]]




def radial_velocity_correction(ra, dec, obs_time, telescope, kind="barycentric",
    **kwargs):

    obstime = Time(obs_time, scale="utc")

    # From https://www.gemini.edu/sciops/telescopes-and-sites/locations
    telescopes = dict([
        ("Gemini-South", EarthLocation.from_geodetic(lon="-70:44:12.096", 
                                                     lat="-30:14:26.700", 
                                                     height=2722 * u.m,
                                                     ellipsoid="WGS84")),
        ("Gemini-North", EarthLocation.from_geodetic(lon="-155:28:08.616", 
                                                     lat="19:49:25.7016", 
                                                     height=4213 * u.m,
                                                     ellipsoid="WGS84")),
    ])

    location = telescopes[telescope]

    return SkyCoord(ra=ra * u.deg, dec=dec * u.deg).radial_velocity_correction(
        kind=kind, obstime=obstime, location=location).to(u.km/u.s)

