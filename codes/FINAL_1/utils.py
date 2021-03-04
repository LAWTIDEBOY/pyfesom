def lag_linregress_3D(x, y, lagx=0, lagy=0):

    # Input: Two xr.Datarrays of any dimensions with the first dim being time. 
    # Thus the input data could be a 1D time series, or for example, have three 
    # dimensions (time,lat,lon). 
    # Datasets can be provided in any order, but note that the regression slope 
    # and intercept will be calculated for y with respect to x.
    # Output: Covariance, correlation, regression slope and intercept, p-value, 
    # and standard error on regression between the two datasets along their 
    # aligned time dimension.  
    # Lag values can be assigned to either of the data, with lagx shifting x, and
    # lagy shifting y, with the specified lag amount. 
    
    #1. Ensure that the data are properly alinged to each other. 
    x,y = xr.align(x,y)

    #2. Add lag information if any, and shift the data accordingly
    if lagx!=0:

        # If x lags y by 1, x must be shifted 1 step backwards. 
        # But as the 'zero-th' value is nonexistant, xr assigns it as invalid 
        # (nan). Hence it needs to be dropped
        x   = x.shift(time = -lagx).dropna(dim='time')

        # Next important step is to re-align the two datasets so that y adjusts
        # to the changed coordinates of x
        x,y = xr.align(x,y)

    if lagy!=0:
        y   = y.shift(time = -lagy).dropna(dim='time')
        x,y = xr.align(x,y)

    #3. Compute data length, mean and standard deviation along time axis: 
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)

    #4. Compute covariance along time axis
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

    #5. Compute correlation along time axis
    cor   = cov/(xstd*ystd)

    #6. Compute regression slope and intercept:
    slope     = cov/(xstd**2)
    intercept = ymean - xmean*slope  

    #7. Compute P-value and standard error
    #Compute t-statistics
    tstats = cor*np.sqrt(n-2)/np.sqrt(1-cor**2)
    stderr = slope/tstats

    from scipy.stats import t
    pval   = t.sf(tstats, n-2)*2
    pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

    return cov,cor,slope,intercept,pval,stderr

def nans(shape, dtype=float):
    a = numpy.empty(shape, dtype)
    a.fill(numpy.nan)
    return a

def linregsum(x,y):
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    import statsmodels.api as sm
    X2 = sm.add_constant(x)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    print(est2.summary())
    return est2

def linreg(x,y):
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    import statsmodels.api as sm
    X2 = sm.add_constant(x)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    return est2

def ols(da, param='slope'):
    """
    Apply statsmodel's OLS regression to a 1d DataArray.
    Handles datetimes (in that case, regression is against days).

    Parameters
    ----------
    data : xarray DataArray
    str, one of ['slope', 'intercept', 'slope_pvalue',
        'intercept_pvalue', 'slope_se', 'intercept_se',]
        sought-after regression parameter

    Notes
    -----
    https://www.statsmodels.org/dev/generated/statsmodels.regression.linear_model.OLS.html

    License
    -------
    GNU-GPLv3, (C) A. Randelhoff
    (https://github.com/poplarShift/python-data-science-utils)
    """
    import statsmodels.api as sm

    data = da.dropna(dim=da.dims[0])

    # specify function with which to retrieve sought-after
    # parameter from statsmodels RegressionResultsWrapper
    res_fn = {
        'intercept': lambda res: res.params[0],
        'slope': lambda res: res.params[1],
        'intercept_pvalue': lambda res: res.pvalues[0],
        'slope_pvalue': lambda res: res.pvalues[1],
        'intercept_se': lambda res: res.bse[0],
        'slope_se': lambda res: res.bse[1],
    }

    if len(data)>=2:
        y = data.values
        xdata = data[data.dims[0]]
        if xdata.dtype.kind in ['M']:
            x = xdata.astype(float).values/1e9/86400.
        else:
            x = xdata.values
        ols = sm.OLS(y, sm.add_constant(x))
        res = ols.fit()
        return res_fn[param](res)
    else:
        return np.nan