import numpy as np
import scipy.interpolate as interpol
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import multiprocessing as mp
import cosmolopy.distance as cd
import sys, logging, os, re

"""
Overview of SATMC

Required Python modules: Scipy, Numpy, Matplotlib, Cosmolopy

Given an input observed SED and N SED templates, SATMC will 
determine the best fit model parameters along with their 68% 
confidence intervals using Monte Carlo Markov Chains and Bayesian 
statistics.  By benefit of Bayes' theorem, additional contraints 
may be imposed  through the use of parameter/likelihood priors. 

Seth Johnson 05/27/2013
"""

np.seterr(divide='ignore') #ignore any divide by 0 errors

#set cosmology and some constants
cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)

c=2.9979e14           #speed of light (microns/s)
h=6.626e-27
dist=50.              #Mpc
z50=70.*50./3.e5      #Redshift at D=50 Mpc

def cal_lum(sed):
    """
    Simple function to compute luminosity (L_sun) of a set of SED 
    templates.

    Input:
         sed - A NumPy array of Template objects
    Output:
         sed - The input sed now containing luminosity information
    """
  
    #assume D(z)=50 Mpc
    for i in xrange(len(sed)):
        lnu=sed[i].fluxsim*dist**2./(1.+z50)/(1.13e-8)/c
        lum=abs(np.trapz(lnu,c/sed[i].wavesim))
        sed[i].lum=lambda: None
        setattr(sed[i].lum,'lum',lum)
    return sed

def value_locate(refx,x):
    """
    Python version of IDL's useful VALUE_LOCATE procedure.
    """
    
    refx=np.array(refx)
    x=np.atleast_1d(x)
    loc=np.zeros(len(x),dtype='int')

    for i in xrange(len(x)):
        ix=x[i]
        ind=((refx-ix) <= 0).nonzero()[0]
        if len(ind) == 0: loc[i]=-1
        else: loc[i]=ind[-1]

    return loc

def find_nearest_grid(grid,p):
    """
    Return indicies of grid bracketing a set of parameters

    Input:
         grid - NumPy array of shape (M,N) containing the set of N
                parameters for M models
         p - NumPy array of shape (N,) containing a set of N parameters
    Output:
         neargrid - A NumPy array of shape (2**N,) containing the locations
                    in grid that bracket p
         intp - A NumPy array of shape (N,) containing the distance
                between p and its nearest grid location
         status - Boolean that returns False if p is outside grid and 
                  True if p is bound by grid
    """
    nparm=len(p)
    ngrid=2**nparm
    status=True
    missing=0
    max_grid=np.empty(nparm,dtype=float)
    min_grid=np.empty(nparm,dtype=float)
    neargrid=np.empty([ngrid],dtype=int)
    t_grid=grid

    #find grid points that bracket all points in P
    for i in xrange(nparm):
        ind=value_locate(grid[:,i],p[i])
        min_grid[i]=grid[max(ind,0),i]
        max_grid[i]=grid[min(ind+1,len(grid[:,0])-1),i]
    n_grid=np.column_stack((min_grid,max_grid))
    intp=(p-min_grid)/(max_grid-min_grid)
    
    #check if interpolation points lie within the grid
    if (np.isnan(intp)).any():
        intp[np.isnan(intp)]=0.0 #point lies on the grid, no problem
    if (np.isinf(intp)).any():
        status=False                #point lies outside the grid
    if ((intp < 0) | (intp > 1)).any():
        status=False                #irregular gridding prevents bracketing
    if not status: return {'neargrid':neargrid,'intp':intp,'status':status}

    #get indices of nearest grid points and mark their location
    if len(grid[:,0]) > 2:
        for i in xrange(ngrid):
            ind=(np.where(grid[:,0] == n_grid[0,i % 2]))[0]
            grid=grid[ind,:]
            for j in xrange(1,nparm):
                g_ind=i/2**j % 2
                sel=grid[:,j] == n_grid[j,g_ind]
                if np.invert(sel).all():
                    missing=1
                    break
                else:
                    grid=grid[sel,:]
                    ind=ind[sel]
            if missing == 1:
                neargrid[i]=-1
                missing=0
            else: neargrid[i]=ind
            grid=t_grid
        grid=t_grid
        
    if (neargrid[:] == -1).any(): status=False
    return (neargrid,intp,status)


def verify_parameters(params,params_old,params_std,param_name,seds,param_min,
                      param_max,**kwargs):
    """
    Confirm that parameter step is bound by parameter grid and priors
    
    Input:
        params - NumPy array of shape (N,) containing current parameter values
        params_old - NumPy array of shape (N,) containing previous parameter 
                     values
        params_std - NumPy array of shape (N,N) containing the parameter 
                     covariance matrix 
        param_name - NumPy array of shape (N,) parameter string names
        seds - NumPy array containing SED objects
        param_min - NumPy array of shape (N,) containing the minimum 
                    parameter values
        param_max - NumPy array of shape (N,) containing the maximum 
                    parameter values
        
        **kwargs - Additional keyword arguments.  Arguments that are 
                   checked for are:
                   priors - List of prior string names
                   prior_mode - List of prior modes, a value of 1
                                indicating a parameter prior

    Output:
        Returns a Boolean False if params is not bound by parameter grid or
        if params does not satisfy given parameter priors.  Returns True
        otherwise

    """
    
    sel=((params >= param_min) & (params <= param_max))
    if not sel.all():
        inv_sel=np.invert(sel)
        params[inv_sel]=(np.random.multivariate_normal(params_old,
                                                       params_std))[inv_sel]
        return False

    #check against parameter grid
    for i in xrange(len(seds)):
        sel=np.array([p in seds[i].param_name for p in param_name])
        if sel.any():
            null,intp,status=find_nearest_grid(seds[i].grid,params[sel])
            if not status:
                if not np.isfinite(intp).all():
                    sel[sel]=np.isfinite(intp)
                params[sel]=(np.random.multivariate_normal(params_old,
                                                           params_std))[sel]
                return False
        
    #check against parameter priors
    if 'priors' in kwargs:
        priors=np.atleast_1d(kwargs.get('priors'))
        p_prior=(np.atleast_1d(kwargs.get('prior_mode')) == 1)
        if p_prior.any():
            priors=priors[p_prior]
            for i in xrange(len(priors)):
                prior_list=re.split('\(|\)',priors[i])
                prior_name=prior_list.pop(0)
                if prior_list[-1] == '' : null=prior_list.pop(-1)
                prior_args=prior_list
                prior_kwargs={'params':params,'param_name':param_name}
                exec("from "+prior_name+" import "+prior_name+" as prior_mod")
                status=prior_mod(*prior_args,**prior_kwargs)
                if not status:
                    params=(np.random.multivariate_normal(params_old,
                                                          params_std))
                    return False

    return True

def shift_sed(wave,flux,redshift):
    """
    Shift a rest-frame model SED to observed frame

    Input:
        wave - NumPy array containing model wavelength values
        flux - NumPy array containing model flux values
        redshift - Float value containing source redshift
    Output:
        s_wave - Observed frame wavelengths
        s_flux - Observed frame fluxes
    """
    owave=wave*(1+redshift)
    ldist=cd.luminosity_distance(redshift,**cosmo)
    oflux=flux*(dist/ldist)**2.*(1+redshift)
    return owave,oflux

def get_model_observed(wave,flux,goodwavel,wavelim,**kwargs):
    """
    Extract model observed fluxes from redshifted SED using observed 
    wavelengths or filters

    Input:
        wave - Numpy array of model wavelenghts
        flux - Numpy array of model fluxes
        goodwavel - Numpy array of observed wavelength bands
        wavelim - Numpy array of bands containing upper limits
        
        **kwargs - Additional keyword arguments.  Arguments checked for:
              filter_list - Numpy array containing filter information

    Output:
        mgood - Model fluxes at goodwavel/filter wavelengths
        mlimit - Model fluxes at wavelim/filter wavelengths
    
    """
    sort=np.argsort(wave)
    wave=wave[sort]
    flux=flux[sort]
    whmod=(np.where((goodwavel < max(wave)) & (goodwavel > min(wave))))[0]
    mgood=np.zeros(len(goodwavel))
    if wavelim != None: mlimit=wavelim*0.
    else: mlimit=None

    if 'filter_list' in kwargs:
        filter_list=kwargs.get('filter_list')
        good=filter_list[:]['good']
        for i in xrange(len(whmod)):
            fwave=filter_list[good][i]['filter']['wave']
            fresp=filter_list[good][i]['filter']['response']
            norm=filter_list[good][i]['filter']['norm']
            fsel=np.interp(fwave,wave,flux,left=0,right=0)
            sel=(fwave !=0 ) & (fsel > 0)
            mgood[whmod[i]]=np.trapz(fsel[sel]*fresp[sel]/(h*c/fwave[sel]),
                                     c/fwave[sel])/norm
        if wavelim != None:
            limits=np.invert(filter_list[:]['good'])
            whmod=(wavelim < max(wave)) & (wavelim > min(flux))
            for i in xrange(len(whmod)):
                fwave=filter_list[limits][i]['filter']['wave']
                fresp=filter_list[limits][i]['filter']['response']
                norm=filter_list[limits][i]['filter']['norm']
                fsel=np.interp(fwave,wave,flux,left=0,right=0)
                sel=(fwave !=0 ) & (fsel > 0)
                mlimit[whmod[i]]=np.trapz(fsel[sel]*fresp[sel]/(h*c/fwave[sel]),
                                          c/fwave[sel])/norm
    else:
        mgood[whmod]=np.interp(goodwavel[whmod],wave,flux)
        if wavelim != None:
            whmod=(wavelim < max(wave)) & (wavelim > min(flux))
            mlimit[whmod]=np.interp(wavelim[whmod],wave,flux)

    return (mgood,mlimit)

def os_gauss(x,s,upperlim):
    """
    Return the likelihood function as either a step function 
    or one-sided Gaussian. 

    Input:
        x - Location of where to compute the likelihood function
        s - Flux value of the upper limit
        upperlim - String containing the type of upper limit prescription
                   to use and its defining parameters
                   1: Step function ('1').  Returns 1 if x < s and 
                      0 if x > s.
                   2: One-sided Gaussian ('2,mu,sigma').  Returns 1 
                      if x < mu and exp(-(x-mu)**2/(2*sigma**2)) 
                      if x > mu.  If mu and/or sigma are not
                      specified, assumes default values of mu=s/1.5
                      and sigma=3*(s-break)/(5*s)
    """
    npt=len(upperlim)
    pdf=np.empty(npt,dtype='longdouble')
    for i in xrange(npt):
        uparm=upperlim[i].split(',')
        if int(uparm[0]) == 1:
            if x[i] < s[i]: pdf[i]=1 
            else: pdf[i]=np.exp(np.longdouble(-5000))
        if int(uparm[0]) == 2:
            if len(uparm) < 2: level=5
            elif uparm[1]=='': level=5
            else: level=float(uparm[1])
            if len(uparm) < 3: mu=s[i]*(3./level)
            elif uparm[2] == '': mu=s[i]*(3./level)
            else: mu=float(uparm[2])
            if len(uparm) < 4: sigma=(s[i]-mu)/level
            elif uparm[3] == '' or mu == 0: sigma=(s[i]-mu)/level
            else: sigma=float(uparm[3])
            if x[i] < mu: pdf[i]=1 
            else: pdf[i]=max(np.exp(np.longdouble(-(x[i]-mu)**2/(2*sigma**2))),
                             np.exp(np.longdouble(-5000)))
    return pdf

def run_chain(queue,chain_id,nstep,params_init,lnl_init,seds,goodwavel,
              goodflux,goodfluxerr,wavelim,fluxlim,upperlim,param_name,
              norm_arr,param_min,param_max,params_std,dr,temp,**kwargs):
    """
    Run a block of steps from a Monte Carlo Markov Chain

    Input:
       queue - Multiprocessing Queue object given from parent process
       chain_id - Integer identifying the chain
       nstep - Integer number of MCMC steps to take
       params_init - Numpy array of parameter values to start from
       lnl_init - Longdouble of the log-likelihood value to start from
       seds - Numpy array of SED objects
       goodwavel - Numpy array of observed wavelengths
       goodflux - Numpy array of observed fluxes
       goodfluxerr - Numpy array of observed flux errors
       wavelim - Numpy array of upper limit wavelengths
       fluxlim - Numpy array of flux upper limits
       upperlim - Numpy array of upper limit information
       param_name - Numpy array of parameter names
       norm_arr - Numpy boolean array indicating if normalizations are to 
                  be applied
       param_min - Numpy array of minimum parameter values
       param_max - Numpy array of maximum parameter values
       params_std - Parameter-parameter covariance matrix
       dr - Numpy boolean array if large dynamic ranges are present
       temp - Chain statistical 'temperature' for parallel tempering

       **kwargs - Additional keyword arguments.  Arguments checked for:
             no_interp - Boolean disables interpolation on the model grid
             synthesis_model - String indicating an SED synthesis model
                               is to be used
             priors - Numpy string array containing priors to use
             prior_mode - Numpy array indicating if priors are
                          parameter priors (1) or likelihood priors (2)

    Output:
        chain_result - Dict containing chain results put into queue and
                       saved in a Numpy save file

    """

    param_arr=np.zeros([len(params_init),nstep])
    lnl_arr=np.zeros(nstep,dtype='longdouble')
    tot_lum_arr=np.zeros(nstep,dtype='longdouble')
    sed=np.zeros(nstep)
    nsed=len(seds)
    lum_arr=np.zeros(nsed)
    if not 'no_interp' in kwargs: method='linear'
    else: method='nearest'

    #set log file, primarily for debug purposes
    log=open('chain_'+str(chain_id)+'.log',"a")
    sys.stdout=log

    #start chain segment
    for n in xrange(nstep):
        mod_obs=[]
        
        #make proposal step
        if n == 0 and lnl_init == 0:
            params_new=params_init
            params_old=params_init
        else:
            if n == 0: params_old=params_init
            else: params_old=param_arr[:,n-1]
            params_new=np.random.multivariate_normal(params_old,params_std)
        while True:
            status=verify_parameters(params_new,params_old,params_std,
                                     param_name,seds,param_min,param_max,
                                     **kwargs)
            if status: break
        
        param_arr[:,n]=params_new

        for l in xrange(len(seds)):
            if (l == 0 and 'synthesis_model' in kwargs):
                #build a dict to contain parameter values
                synth_tags=seds[0].synth_param
                synth_val=[(params_new[np.array(param_name) == 
                                       synth_tags[i]])[0] 
                           if (synth_tags[i] in param_name) 
                           else seds[0].synth_val[i] for i in 
                           xrange(len(synth_tags))]
                synth_kwargs=dict(zip(synth_tags,synth_val))
                synthesis_model=kwargs.get('synthesis_model')

                #load in and call SED synthesis routine
                exec("from "+synthesis_model+" import "+synthesis_model+
                     " as synth_mod")
                wave,flux,err=synth_mod('model'+str(chain_id),**synth_kwargs)
                
                #calculate bolometric luminosity
                lnu=flux*dist**2/(1+z50)*(3.124e-7)
                lum_arr[0]=abs(np.trapz(lnu,c/wave))

                #correct SED for redshift
                s_wave,s_flux=shift_sed(wave,flux,params_new[-1])
                synth_sed={'wave_synth':s_wave,'flux_synth':s_flux}
                good,limits=get_model_observed(s_wave,s_flux,goodwavel,wavelim,
                                               **kwargs)
                mod_obs.append([{'loc':0,'good':good,'limits':limits}])
            else:
                sel=np.array([p in seds[l].param_name for p in param_name])
                if sel.any():
                    int_grid,intp,status=find_nearest_grid(seds[l].grid,
                                                params_new[sel])
                else: int_grid=[0]
                
                #extract SED information from models
                wave=[s['wavesim'] for s in seds[l].Template[int_grid]]
                flux=[s['fluxsim'] for s in seds[l].Template[int_grid]]
                lum=[s['lum'] for s in seds[l].Template[int_grid]]
                if norm_arr[l]:
                    norm=np.array(param_name) == 'norm'+str(l+1)
                    if dr[norm[0:len(norm)-1]]:
                        flux=flux*10**params_new[norm]
                        lum=lum*10**params_new[norm]
                    else:
                        flux=flux*params_new[norm]
                        lum=lum*params_new[norm]
                else: pass
                pgrid=seds[l].grid[int_grid,:]
                lum_arr[l]=interpol.griddata(pgrid,np.array(lum),(params_new[sel]).reshape((1,sum(sel))),method=method)
                
                t_obs=np.empty(0,dtype='object')
                for i in xrange(len(int_grid)):
                    s_wave,s_flux=shift_sed(wave[i],flux[i],params_new[-1])
                    good,limits=get_model_observed(s_wave,s_flux,goodwavel,
                                                   wavelim,**kwargs)
                    t_obs=np.append(t_obs,{'loc':int_grid[i],'good':good,
                                           'limits':limits})
                mod_obs.append(t_obs)

        #combine fluxes from SED sets
        sel=np.array([('norm' not in name) and ('redshift' not in name) 
                      for name in param_name])
        ndim=sum(sel)
        if wavelim == None:
            mod_tot=np.zeros(2**ndim,dtype=[('params',float,(ndim,)),
                                            ('good',float,len(goodwavel))])
        
        else:
            mod_tot=np.zeros(2**ndim,dtype=[('params',float,(ndim,)),
                                            ('good',float,len(goodwavel)),
                                            ('limits',float,len(wavelim))])
        pdim=1
        p0=0
        for l in xrange(nsed):
            nmod=len(mod_obs[l])
            nparm=seds[l].nparm
            for i in xrange(2**ndim):
                ind=(i/pdim % nmod)
                mod_tot[i]['params'][p0:p0+nparm]= \
                    seds[l].grid[mod_obs[l][ind]['loc'],:]
                mod_tot[i]['good']+=mod_obs[l][ind]['good']
                if wavelim != None:
                    mod_tot[i]['limits']+=mod_obs[l][ind]['limits']
            pdim*=nmod
            p0+=nparm

        #exclude wavelengths outside SEDs
        whgood=(mod_tot[0]['good'] != 0)
        gf=goodflux[whgood]
        gw=goodwavel[whgood]
        gfe=goodfluxerr[whgood]
        if wavelim != None:
            whlimit=(mod_tot[0]['limits'] != 0)
            if sum(whlimit) > 0:
                wl=wavelim[whlimit]
                fl=fluxlim[whlimit]
            else: wavelim=None

        #determine likelihood at each combination of models
        lnL=np.zeros(2**ndim,dtype='longdouble')
        for i in xrange(2**ndim):
            lnL[i]=sum(-(gf-mod_tot[i]['good'][whgood])**2./(2*gfe**2))
            if (not 'ig_lim' in kwargs) and (wavelim != None):
                sel=np.array([int(u[0]) for u in upperlim]) != 0
                lnL[i]+=sum(np.log(os_gauss(mod_tot[i]['limits'][whlimit], 
                                            fl,upperlim[sel[whlimit]])))
            #verify step exists, only happens for synthesized SEDs
            if (mod_tot[i]['good'] == 0).all(): lnL[i]=-10**np.longdouble(3000)
        
        #interpolate likelihoods over model grid
        
        if ndim > 0:
            l=lnL
            maxl=max(l)
            l-=maxl
            l=np.exp(l)
            l=interpol.griddata(mod_tot[:]['params'],l,
                                (params_new[sel]).reshape((1,ndim)),
                                method=method)
            lnL=np.log(max(l,np.exp(np.longdouble(-5000))))+maxl
        
        #check versus likelihood priors
        if 'priors' in kwargs:
            priors=np.atleast_1d(kwargs.get('priors'))
            l_prior=(np.atleast_1d(kwargs.get('prior_mode')) == 2)
            if l_prior.any():
                priors=priors[l_prior]
                for i in xrange(len(priors)):
                    prior_list=re.split('\(|\)',priors[i])
                    prior_name=prior_list.pop(0)
                    if prior_list[-1] == '' : null=prior_list.pop(-1)
                    prior_args=prior_list
                    prior_kwargs={'params':params_new,'param_name':param_name,
                                  'lum_arr':lum_arr,'seds':seds}
                    exec("from "+prior_name+" import "+prior_name+
                         " as prior_mod")
                    lnL+=prior_mod(*prior_args,**prior_kwargs)
        
        lnl_arr[n]=lnL[0]
        tot_lum_arr[n]=sum(lum_arr)
        
        #finally, decide whether or not to take the step
        if not (n == 0 and lnl_init == 0):
            if n == 0:
                alpha=min(1,np.exp((lnl_arr[n]-np.longdouble(lnl_init))/(1+temp)))
            else:
                alpha=min(1,np.exp((lnl_arr[n]-lnl_arr[n-1])/(1+temp)))
            urand=np.random.random()
            if urand < alpha: param_arr[:,n]=params_new 
            else:
                param_arr[:,n]=params_old
                if n == 0: lnl_arr[n]=lnl_init 
                else: lnl_arr[n]=lnl_arr[n-1]
        else: param_arr[:,n]=params_new
    
    #report results of the chain and save to file
    chain_result={'param_arr':param_arr,'lnl_arr':lnl_arr,
                  'lum_arr':tot_lum_arr,'chain_id':chain_id}
    if 'synthesis_model' in kwargs: chain_result['sed']=synth_sed
    queue.put(chain_result)
    np.save('chain_'+str(chain_id),chain_result)


def cal_lgrid(seds,wave,flux,fluxerr,upperlim,norm_arr,param_max,
              param_min,param_name,nparm,dr,minz,maxz,**kwargs):
    """
    Pre-compute likelihoods for grid-based templates.  Likelihood grid
    will be used to determine initial photo-z
    """

    seddim=np.empty(0)
    nsed=len(seds)
    for i in xrange(nsed):
        if norm_arr[i]: seddim=np.append(seddim,31)
        seddim=np.append(seddim,len(seds[i].Template))
    nz=np.rint((maxz-minz)/0.1)
    zgrid=np.linspace(minz,maxz,nz+1)
    seddim=np.append(seddim,nz+1)
    ngrid=np.prod(seddim)

    grid=np.empty((ngrid,nparm+1),dtype=float)
    lgrid=np.empty(ngrid,dtype=np.longdouble)

    j=0
    indx=np.unravel_index(np.arange(0,ngrid,dtype=int),seddim)
    for i in xrange(nsed):
        if norm_arr[i]:
            norm=np.array(param_name) == 'norm'+str(i+1)
            normgrid=np.linspace(param_min[j],param_max[j],31)
            grid[:,norm]=np.reshape(normgrid[indx[j]],(ngrid,1))
            j+=1
        sel=np.array([p in seds[i].param_name for p in param_name])
        grid[:,sel]=seds[i].grid[indx[j],:]
        j+=1
    zp=zgrid[indx[j]]
    grid[:,nparm]=zp

    if not 'flux_grid' in kwargs:
        ldist=cd.luminosity_distance(zgrid,**cosmo)
        flux_grid=np.empty((len(wave),ngrid))

        for i in xrange(int(ngrid)):
            indx=np.unravel_index(i,seddim)
            j=0
            for k in xrange(nsed):
                mw=seds[k].Template[0].wavesim
                if norm_arr[k]:
                    sel=np.array(param_name) == 'norm'+str(k+1)
                    norm=grid[i,sel]
                    if dr[sel[:-1]]: norm=10**norm
                    mf=seds[k].Template[indx[j+1]].fluxsim*norm
                    j+=2
                else:
                    mf=seds[k].Template[indx[j]].fluxsim
                    j+=1
                mws=mw*(1+grid[i,nparm])
                mfs=mf*(dist/ldist[indx[j]])**2*(1+grid[i,nparm])

                flux_grid[:,i]+=np.interp(wave,mws,mfs)

    good=np.array([int(a[0]) for a in upperlim]) == 0
    bad=np.invert(good)
    n_good=sum(good)
    n_bad=sum(bad)
    goodflux=(np.resize(flux[good],(ngrid,n_good))).T
    goodfluxerr=(np.resize(fluxerr[good],(ngrid,n_good))).T
    if n_bad != 0: fluxlim=(np.resize(flux[bad],(ngrid,nbad))).T
    
    lgrid=np.sum(((goodflux-flux_grid[good,:])**2/(2*goodfluxerr**2)),0)
    if not 'ig_lim' in kwargs and n_bad != 0:
        if n_bad == 1:
            lgrid+=np.log(os_gauss(flux_grid[bad,:],fluxlim,upperlim[bad]))
        else: lgrid+=np.sum(np.log(os_gauss(flux_grid[bad,:],fluxlim,
                                            upperlim[bad])),0)

    #check with priors
    if 'priors' in kwargs:
        priors=kwargs.get('priors')
        l_priors=kwarg.get('prior_mode') == 2
        if l_prior.any():
            priors=priors[l_prior]
            for i in xrange(len(priors)):
                prior_list=re.split('\(|\)',priors[i])
                prior_name=prior_list.pop(0)
                if prior_list[-1] == '' : null=prior_list.pop(-1)
                prior_args=prior_list
                prior_kwargs={'params':grid,'param_name':param_name,
                              'lnl':lgrid,'lum_arr':lum_arr,'seds':seds}
                exec("from "+prior_name+" import "+prior_name+" as prior_mod")
                lgrid+=prior_mod(*prior_args,**prior_kwargs)

    return (grid,lgrid)
    
#Python object for containing SED templates
class sed_object:
    def __init__(self,Template,synthesis_model=None):
        self.Template=Template
        #Get parameter names and construct parameter grid
        tags=[x for x in Template.dtype.names]
        req_tags=['modelname','wavesim','fluxsim','lum']
        if synthesis_model == None:
            param_name = [s for s in tags if s.lower() not in req_tags ]
            grid=np.empty([len(Template),len(param_name)])
            for i in xrange(len(param_name)):
                grid[:,i]=Template[param_name[i]]
            self.grid=grid
            self.ngrid=[len(np.unique(grid[:,i])) 
                           for i in xrange(len(param_name))]
        else:
            synth_param=[s for s in tags]
            self.synth_param=synth_param
            self.synth_val=Template[synth_param][0]
            param_name=[]
            grid=[]
            for i in xrange(len(synth_param)):
                if len(np.atleast_1d(self.synth_val[i])) == 2:
                    param_name.append(synth_param[i])
                    grid.append(self.synth_val[i])
            self.grid=np.array(grid).T
            self.ngrid=[2 for i in xrange(len(param_name))]
        self.param_name=param_name
        self.nparm=len(param_name)
                   
def satmc(filename,*args,**kwargs):
    """
    Python based SATMC

    Inputs:
         filename - String name containing the observations of current
                    source.  The source name will be taken from filename
                    
         *args - Additional arguments.  Each argument is interpreted 
                 as a Numpy array of Template objects to be used as a
                 separate SED template set

         **kwargs - Additional keyword arguments.  Arguments checked for:
              norm* - Model normlization(s)
              redshift - Source redshift
              findz - Boolean to turn on photometric redshift estimation
              axis - Set the ranges for plot axes
              flam - Boolean to plot the SED and best fit model in 
                     lambdaF_lambda units
              no_plot - Boolean to disable plotting
              nchain - Number of chains to run
              nstep - Total number of MCMC steps to keep per chain
              step_int - Number of MCMC steps in each chain segment.
                         Chain results will be returned to parent level
                         and saved after every step_int steps.
              param_min - Numpy array of user defined parameter minimums
              param_max - Numpy array of user defined parameter maximums
              params_std - User defined parameter covariance matrix
              no_temp - Boolean to disable parallel tempering
              ig_lim - Boolen to ignore upper limits
              adaptive_step - Boolean to disable processing the burn-in
              no_interp - Boolean disables interpolation on the model grid
              synthesis_model - String name of the SED sythnesis model 
                                to use in place of a SED template
              filters - List or Numpy array of string names for the
                        photometric filters
              priors - Numpy string array containing priors to use
              prior_mode - Numpy array indicating if priors are
                          parameter priors (1) or likelihood priors (2)
              restart - Integer value specifying the point at which to 
                        begin processing (restoring from a previous point)
                         1 - start from the burn-in period
                         2 - start after the burn-in
              outfile - Root string name of the output text file and Numpy
                        save file
              psfile - String name of the output (encapsulated) post-script
                       to save the plot

    Output:
        Best fit parameter values, luminosity, redshift and log-likelihood
        for the data given the models.  Information is printed to screen
        and saved in a text file given by outfile.

        Plot of the data and best fit models.  Plot is saved to a file
        given by psfile.
              
    """


    #read in acceptable keyword arguments
    if len(args)==0:
        print 'At least 1 SED template required.'
        return
    else: nsed=len(args)
    redshift=kwargs.get('redshift', None)
    findz=kwargs.get('findz',False)
    if redshift==None and findz==False:
        print ('Warning: redshift and/or findz not set. '+
               'Running with redshift=2.0,findz=True.')
        redshift=2.0
        findz=True
    
    axis=kwargs.get('axis',None)       #Plot x-range
    flam=kwargs.get('flam',False)      #Set plotting in lambda F_lambda
    nchain=kwargs.get('nchain',4)      
    nstep=kwargs.get('nstep',10000)    
    step_int=kwargs.get('step_int',500)
    synthesis_model=kwargs.get('synthesis_model',None)
    filters=kwargs.get('filters',None)
    flux_grid=kwargs.get('flux_grid',None)
    if 'priors' in kwargs and not ('prior_mode' in kwargs):
        print ('Warning: Prior type not specified. '+
               ' Running with Type 2 likelihood priors.')
        kwargs['prior_mode']=np.zeros(len(priors),dtype=int)+2

    #check number of SED models passed as additional arugments
    #combine SEDs into single object
    seds=np.empty(nsed,dtype=object)
    if nsed > 0:
        for i in xrange(nsed):
            sed=(args[i])[:]
            tags=sed.dtype.names
            if (sum([s in ['wavesim','fluxsim'] for s in tags]) < 2 and 
                synthesis_model == None):
                print 'SED object must contain WAVESIM and FLUXSIM attributes.'
                return
            if 'lum' not in tags and synthesis_model == None: sed=cal_lum(sed)
            seds[i]=sed_object(args[i],synthesis_model=synthesis_model)

    #read in file and units
    with open(filename) as obsfile:
        units=obsfile.readline().splitlines()
    units=units[0].split()
    data=np.genfromtxt(filename,skiprows=1,dtype=None)
    wavel,flux,fluxerr,upperlim=[data[i] for i in data.dtype.names]
    if units[0].lower()[0] == 'l':
        wavel=(c/10**wavel)
    elif units[0].lower()[0] == 'a':
        wavel*=1.e-4
    elif units[0].lower()[0] != 'm':
        wavel=(c/10**np.append(units[0],wavel))
    if units[1].lower()[0] == 'm':
        flux*=1.e-3
    elif units[1].lower()[0] == 'u':
        flux*=1.e-6
    elif units[1].lower()[0] != 'j':
        flux=np.append(units[1],flux)*1.e-3
    if units[2].lower()[0] == 'm':
        fluxerr*=1.e-3
    elif units[2].lower()[0] == 'u':
        fluxerr*=1.e-6
    elif units[2].lower()[0] != 'j':
        fluxerr=np.append(units[2],fluxerr)*1.e-3
    if len(units) == 4:
        if units[3].isdigit():
            upperlim=np.append(units[3],upperlim)

    #check for compatible filters
    if filters != None:
        filters=np.array(filters)
        if len(filters) != len(wavel):
            print ('Number of photometric filters must match'+
                   ' number of observed wavelengths.')
            return

    #remove 0 flux elements
    sel=np.nonzero(flux)
    flux=flux[sel]
    fluxerr=fluxerr[sel]
    wavel=wavel[sel]
    upperlim=upperlim[sel]
    if filters != None: filters=filters[sel]
    if flux_grid != None: flux_grid=flux_grid[sel,:]

    #sort by wavelength
    sort=np.argsort(wavel)
    wavel=wavel[sort]
    flux=flux[sort]
    fluxerr=fluxerr[sort]
    if filters != None: filters=filters[sort]
    if flux_grid != None: flux_grid=flux_grid[sort,:]

    #check for upper limits
    good=np.array([int(a) for a in upperlim]) == 0
    bad=np.invert(good)
    n_good=sum(good)
    n_bad=sum(bad)
    if n_bad != 0:
        wavelim=wavel[bad]
        fluxlim=max(flux[bad],fluxerr[bad])
    else:
        wavelim=None
        fluxlim=None
    goodwavel=wavel[good]
    goodflux=flux[good]
    goodfluxerr=fluxerr[good]

    #check there is enough points for a reasonable fit
    print 'Number of useful data points:',n_good
    if n_bad != 0: print 'Number of useful upper limits:',n_bad
    if n_good <= 2:
        print 'Not enough data points, quitting.'
        return

    #plot observed
    if not 'no_plot' in kwargs:
        plt.ion()
        plt.clf()
        if axis == None:
            axis=np.zeros(4,dtype=float)
        if axis[0] == 0:
            axis[0]=min(10.**(np.floor(np.log10(goodwavel))))*0.1
        if axis[1] == 0:
            axis[1]=max(10.**(np.ceil(np.log10(goodwavel))))*10.
        if axis[2] == 0:
            axis[2]=min(goodflux)*0.1
        if axis[3] == 0:
            axis[3]=max(goodflux)*10.
        if flam == False:
            plt.loglog(goodwavel,goodflux,'ks',mfc='none')
            plt.errorbar(goodwavel,goodflux,goodfluxerr,fmt=None,ecolor='black')
            if n_bad > 0:
                plt.errorbar(wavelim,fluxlim,fluxlim*.2,lolims=True,
                             ecolor='black')
            plt.xlabel(r'Observed Wavelength ($\mu$m)')
            plt.ylabel('Observed Flux (Jy)')
        else:
            lfl=goodflux*3.e-9/goodwavel
            err=goodfluxerr*3.e-9/goodwavel
            plt.loglog(goodwavel,lfl,'ks',mfc='none')
            plt.errorbar(goodwavel,lfl,err,fmt=None,ecolor='black')
            if n_bad > 0:
                lfl=fluxlim*3.e-9/wavelim
                plt.errorbar(wavelim,fluxlim,fluxlim*.2,lolims=True,
                             ecolor='black')
            plt.xlabel(r'Observed Wavelength ($\mu$m)')
            plt.ylabel(r'$\lambda$F$_{\lambda}$ (ergs cm$^{-2}$ s$^{-1}$)')
            axis[2]=min(lfl)*0.1
            axis[3]=max(lfl)*10.
        title=filename.rstrip('.dat').split('/')
        title=np.array(title)     #Ensures a single string is properly indexed
        plt.title(title[-1])
        plt.axis(axis)
        
    #---------------The MCMC parts--------------
    #Setup parameter arrays
        
    #Boolean array to turn on model normalization
    norm_arr=np.array(['norm'+str(i+1) in kwargs for i in xrange(nsed)])
    if synthesis_model != None: norm_arr[0]=False
    #setup filters
    if filters != None:
        filter_list=np.empty(len(filters),dtype=[('good',bool),
                                                 ('filter',object)])
        filter_list[:]['good']=good
        for i in xrange(len(filters)):
            fwave,fresponse=np.loadtxt(filters[i],unpack=True)
            fnorm=abs(np.trapz(c/fwave,fresponse/(h*c/fwave)))
            filter_list[i]['filter']={'wave':fwave,'response':fresponse,
                                      'norm':fnorm}
        kwargs['filter_list']=filter_list

    #Get parameter grid
    param_name=[]
    pmin=np.empty(1)
    pmax=np.empty(1)
    npgrid=[]
    nparm=0
    for i in xrange(nsed):
        if norm_arr[i]:
            param_name+=(['norm'+str(i+1)])
            norm=np.atleast_1d(kwargs.get('norm'+str(i+1)))
            if len(norm) == 2:
                pmin=np.append(pmin,[norm[0]])
                pmax=np.append(pmax,[norm[1]])
            else:
                pmax=np.append(pmax,norm[0])
                pmin=np.append(pmin,1.0)
            pmax=np.append(pmax,np.amax(seds[i].grid,axis=0))
            pmin=np.append(pmin,np.amin(seds[i].grid,axis=0))
            nparm+=1
            npgrid+=[2]
        else:
            pmax=np.append(pmax,np.amax(seds[i].grid,axis=0))
            pmin=np.append(pmin,np.amin(seds[i].grid,axis=0))
        param_name+=seds[i].param_name
        nparm+=seds[i].nparm
        npgrid+=seds[i].ngrid
    pmin=np.delete(pmin,0)
    pmax=np.delete(pmax,0)
        
    if 'param_max' in kwargs:
        param_max=np.minimum(kwargs.get('param_max'),pmax)
    else: param_max=pmax
    if 'param_min' in kwargs:
        param_min=np.maximum(kwargs.get('param_min'),pmin)
    else: param_min=pmin
    
    #Check for large dynamical xranges
    r=abs(np.log10(param_max)-np.log10(param_min))
    dr=(((r >= 2) & (param_min != 0)) | 
        ((np.log10(param_max) >= 2) & (param_min == 0)))

    if dr.any():
        param_max[dr]=np.log10(param_max[dr])
        param_min[dr]=np.log10(param_min[dr])
        np.place(param_min,np.isinf(param_min),0)
        for i in xrange(nsed):
            if seds[i].nparm > 1:
                gs=np.array([p in seds[i].param_name for p in param_name])
                seds[i].grid[:,dr[gs]]=np.log10(seds[i].grid[:,dr[gs]])

    #The MCMC array that will contain all the steps
    chain=np.zeros(nchain,dtype=[('lnl_iter',np.longdouble,(nstep,)),
                                 ('lum_iter',np.longdouble,(nstep,)),
                                 ('params_iter',float,(nparm+1,nstep)),
                                 ('params_std',float,(nparm+1,nparm+1))])
    
    #set initial guesses
    if 'params_init' not in kwargs:
        for i in xrange(nchain):
            chain['params_iter'][i,0:nparm,0]= \
                np.random.uniform(param_min,param_max)
    else:
        params_init=np.array(kwargs.get('params_init'))
        #check user supplied initial guess is within parameter space
        if len(dr) != 0:
            np.put(params_init,dr,np.log10(params_init[dr]))
        if ((params_init < param_min) | (params_init > param_max)).any():
            print 'Initial parameter guess outside parameter bounds.'
            return
        else:
            chain['params_iter'][:,0:nparm,0]= \
                np.resize(params_init,[nchain,len(params_init)])
    
    #Determine if redshift a free or fixed parameter
    z=np.empty([nchain,nstep])
    if findz:
        maxz=kwargs.get('maxz',6.0)
        minz=kwargs.get('minz',1.e-3)
        delz=kwargs.get('delz',0.01)
        if redshift != None:
            z[0,:]=redshift
        else:
            z[0,:]=np.random.uniform(minz,maxz)
    else:
        z[:]=redshift
        minz=redshift
        maxz=redshift
        delz=0.
    param_name+=['redshift']
    chain['params_iter'][:,nparm,:]=z
    param_max=np.append(param_max,maxz)
    param_min=np.append(param_min,minz)

    #Multi-normal covariance matrix to draw steps
    if 'params_std' not in kwargs:
        params_std=(0.01*(param_max-param_min)/(nparm+1))**2.
        for i in xrange(len(npgrid)):
            if npgrid[i] > 2:
                params_std[i]=((param_max[i]-param_min[i])/npgrid[i])**2.
        params_std[nparm]=delz**2.
        chain['params_std']=np.resize(np.diag(params_std),
                                      [nchain,nparm+1,nparm+1])
    else:
        params_std=np.array(kwargs.get('params_std'))
        if (params_std.shape)[0] != nparm+1:
            print ('Input covariance matrix must have dimensions'+
                   ' equal to the number of available parameters.')
            return
        chain['params_std']=np.resize(params_std,[nchain,nparm+1,nparm+1])
    
    #pre-compute likelihood grid for photo-z estimation
    if 'cal_lgrid' in kwargs and 'findz' in kwargs:
        pgrid,lgrid=cal_lgrid(seds,wavel,flux,fluxerr,upperlim,norm_arr,
                              param_max,param_min,param_name,nparm,dr,minz,
                              maxz,**kwargs)
        if redshift != None:
            sort=(np.argsort(lgrid))[:len(lgrid)-nchain-1:-1]
            chain['params_iter'][:,:,0]=pgrid[sort,:]
        sort=(np.argsort(lgrid))[:len(lgrid)-max(int(len(lgrid)*0.001),100):-1]
        cov=np.cov(np.float32(pgrid[sort,:].T))*1.e-5
        chain['params_std']=np.resize(cov,(nchain,nparm+1,nparm+1))
    
    #do the MCMC chain
    thres=2./step_int         #threshold value for acceptance rate
    pchange=0                 #rate of accepted steps
    rold=0
    pairold=0
    lastswap=0
    t_period=0
    
    #set chain temperatures
    temp=np.array([10**(3.5*i/nchain)-1 for i in xrange(nchain)])
    if 'no_temp' in kwargs: temp[:]=0
    
    #run MCMC chains in blocks of length step_init, first test for 
    #convergence and acceptance rate (i.e. the burn-in period)
    if not 'adaptive_step' in kwargs:
        if kwargs.get('restart',None) == 1:
            #restore chains
            chain=np.load('parent_burnin.npy')

        if kwargs.get('restart',None) <= 1:
            while (pchange <= 0.2 or pchange >= 0.26):
                #if doing parallel tempering, make swaps between 
                #fudical cold and warm chains
                if not (temp == 0).all() and t_period == 2:
                    t_period=0
                    psel=np.argmax(chain['lnl_iter'][:,0:step_int*2],axis=1)
                    ind=np.indices(psel.shape)
                    maxl=((chain['lnl_iter'][:,0:step_int*2])[ind,psel])[0]
                    maxsel=((maxl == np.max(maxl)).nonzero())[0]
                    
                    if maxsel[0] != 0:
                        chain['params_iter'][0,:,0]= \
                            chain['params_iter'][maxsel[0],:,psel[maxsel]]
                        chain['params_iter'][maxsel[0],:,0]= \
                            chain['params_iter'][0,:,psel[0]]
                        chain['lnl_iter'][0,0]= \
                            chain['lnl_iter'][maxsel[0],psel[maxsel]]
                        chain['lnl_iter'][maxsel[0],0]= \
                            chain['lnl_iter'][0,psel[0]]
                else: t_period+=1
    
                #run the chain
                procs=[]
                chain_queue=mp.Queue()

                for i in xrange(nchain):
                    p=mp.Process(target=run_chain,
                                 args=(chain_queue,i,step_int*2,
                                       chain['params_iter'][i,:,0],
                                       chain['lnl_iter'][i,0],seds,
                                       goodwavel,goodflux,goodfluxerr,
                                       wavelim,fluxlim,upperlim,param_name,
                                       norm_arr,param_min,param_max,
                                       chain['params_std'][i,:,:],dr,temp[i]),
                                 kwargs=kwargs)
                    procs.append(p)
                    p.start()
        
                #get chain results
                for p in procs:
                    chain_result=chain_queue.get()
                    ind=chain_result['chain_id']
                    if synthesis_model != None:
                        sed=chain_result['sed']
                        if not 'sed' in chain.dtype.names:
                            dt=chain.dtype.descr+[('sed',float,sed.shape)]
                            chain=np.empty(chain.shape,dtype=dt)
                        chain['sed'][ind]=sed
                    chain['params_iter'][ind,:,0:step_int*2]= \
                        chain_result['param_arr']
                    chain['lnl_iter'][ind,0:step_int*2]=chain_result['lnl_arr']
                    chain['lum_iter'][ind,0:step_int*2]=chain_result['lum_arr']

                #wait for processes to finish
                for p in procs:
                    p.join()

                #check acceptance rate of chains
                params_std_old=chain['params_std']
                pchange_ch=np.empty(nchain)

                for i in xrange(nchain):
                    dlnl_iter=chain['lnl_iter'][i,step_int-1:step_int*2-1]- \
                        chain['lnl_iter'][i,step_int:step_int*2]
                    whchange=abs(dlnl_iter[0:step_int-2]) > 0
                    pchange_ch[i]=np.float(sum(whchange))/(step_int-1)
                
                #tune the multi-normal covariance to reach optimal acceptance
                if 'no_temp' in kwargs: pchange=np.average(pchange_ch)
                else: pchange=pchange_ch[0]
                if (pchange == 0).any(): pchange=0.
                if (pchange < 0.20) or (pchange > 0.26):
                    for i in xrange(nchain):
                        cov=chain['params_std'][i,:,:]
                        if pchange_ch[i] < 0.05: cov=cov*0.05
                        elif pchange_ch[i] > 0.26:
                            diag_cov=np.diag(cov)
                            if np.array_equal(cov,np.diag(diag_cov)):
                                cov=np.cov(np.float32(chain['params_iter'][i,:,step_int:step_int*2]))
                            else:
                                frac=pchange_ch[i]/0.23
                                cov=cov*frac
                        else:
                            cov=np.cov(np.float32(chain['params_iter'][i,:,step_int:step_int*2]))
                            
                        sel=((np.sqrt(np.diag(cov)) < 
                              np.sqrt(params_std)*1.e-6) |
                             (np.sqrt(np.diag(cov)) > 
                              np.sqrt(params_std)*1.e4))

                        if (np.array_equal(cov,np.zeros(cov.shape)) or 
                            sel.any()): 
                            cov=np.diag(params_std)
                        chain['params_std'][i,:,:]=cov

                #check convergence using the Geweke dianogstic
                if (pchange >= 0.20) and (pchange <= 0.26):
                    if not 'no_temp' in kwargs: ra=xrange(1)
                    else: ra=xrange(nchain)
                    for i in ra:
                        av1=np.average(chain['lnl_iter'][i,0:int(0.2*step_int)])
                        av2=np.average(chain['lnl_iter'][i,step_int:2*step_int])
                        s=np.std(chain['lnl_iter'][i,step_int:2*step_int])
                        if s != 0: 
                            if abs(av1-av2)/s > 2.0: pchange=0
                            
                #check for bad chains, chains that are either unable to
                #move about (flat likelihood space) or too far from 
                #current maximum likelihood
                maxl_ch=np.empty(nchain)
                for i in xrange(nchain):
                    lnl_chain=chain['lnl_iter'][i,:]
                    maxl_ch[i]=np.amax(lnl_chain[np.nonzero(lnl_chain)])
                diffl=maxl_ch-np.amax(maxl_ch)

                if (pchange >= 0.2) and (pchange <= 0.26):
                    badchain=(pchange_ch < thres) | (abs(diffl) > 1.e2)
                    goodchain=np.invert(badchain)
                    badchain=(badchain.nonzero())[0]
                else:
                    badchain=(pchange_ch < thres) & (diffl != 0)
                    goodchain=np.invert(badchain)
                    badchain=(badchain.nonzero())[0]
                if len(badchain) != 0:
                    if len(badchain) == nchain:
                        for i in badchain:
                            chain['params_iter'][i,:,0]= \
                                np.random.uniform(param_min,param_max)
                            chain['lnl_iter'][i,0]=0
                            chain['params_std'][i,:,:]=params_std_old[i,:,:]
                    else:
                        for i in badchain:
                            rind=np.rint(np.random.random()*(step_int-1))
                            if (pchange >= 0.2) and (pchange <= 0.26):
                                rchain=0
                            else:
                                rchain=np.random.choice((goodchain.nonzero())[0])
                            
                            chain['params_iter'][i,:,0]= \
                                chain['params_iter'][rchain,:,rind]
                            chain['lnl_iter'][i,0]= \
                                chain['lnl_iter'][rchain,rind]
                            chain['params_std'][i,:,:]= \
                                params_std_old[rchain,:,:]
                        good=goodchain.nonzero()
                        chain['params_iter'][good,:,0]= \
                            chain['params_iter'][good,:,2*step_int-1]
                        chain['lnl_iter'][good,0]= \
                            chain['lnl_iter'][good,2*step_int-1]
                    pchange=0
                else:
                    chain['params_iter'][:,:,0]= \
                        chain['params_iter'][:,:,2*step_int-1]
                    chain['lnl_iter'][:,0]=chain['lnl_iter'][:,2*step_int-1]

                #create/update savefile
                np.save('parent_burnin',chain)
            params_std=chain['params_std'][0,:,:]
        period=2
    else: period=0
    #With burn-in complete, run the full chains
    #separate the chains into blocks of STEP_INT length for restart points
    if kwargs.get('restart',None) == 2:
        #restore chains if necessary
        chain=np.load('parent_chain.npy')

    if kwargs.get('restart',None) <= 2:
        swap=0
        n=period
        while n < nstep/step_int:
            #perform a Metropolis-Hastings swap for a pair of adjacent chains
            if not 'no_temp' in kwargs:
                while True:
                    rindex=np.random.choice(xrange(nchain))
                    while True:
                        if np.random.random() > 0.5:
                            pair=min(max(0,rindex+1),nchain-1)
                        else: pair=min(max(0,rindex-1),nchain-1)
                        if pair != rindex: break
                    if not (((rold == rindex) or (rold == pair)) and
                            ((pairold == rindex) or (pairold == pair))): break
                
                tpair=10**(3.5*pair/nchain)-1
                trindex=10**(3.5*rindex/nchain)-1
                sel=np.array(xrange((n-1)*step_int,n*step_int))
                maxpair=np.amax(chain['lnl_iter'][pair,sel])
                maxrindex=np.amax(chain['lnl_iter'][rindex,sel])
                psel=np.argmax(chain['lnl_iter'][pair,sel])
                rsel=np.argmax(chain['lnl_iter'][rindex,sel])
                alpha=min(1,np.exp(maxrindex-maxpair)** \
                              (1/(1+tpair)-(1/(1+trindex))))
                urand=np.random.random()
                if ((abs(maxpair-maxrindex) > 1.0) and (urand < alpha) and 
                    (swap == 0) and (n > period)):
                    chain['params_iter'][rindex,:,n*step_int-1]= \
                        chain['params_iter'][pair,:,sel[psel]]
                    chain['params_iter'][pair,:,n*step_int-1]= \
                        chain['params_iter'][pair,:,sel[rsel]]
                    chain['lnl_iter'][rindex,n*step_int-1]= \
                        chain['lnl_iter'][pair,sel[psel]]
                    chain['lnl_iter'][pair,n*step_int-1]= \
                        chain['lnl_iter'][pair,sel[rsel]]
                    if (pair == 0) or (rindex == 0):
                        swap=1
                        lastswap=n
                    else: swap=0
                    rold=rindex
                    pairold=pair

            #run the chains
            procs=[]
            chain_queue=mp.Queue()
    
            for i in xrange(nchain):
                p=mp.Process(target=run_chain,
                             args=(chain_queue,i,step_int,
                                   chain['params_iter'][i,:,n*step_int-1],
                                   chain['lnl_iter'][i,n*step_int-1],seds,
                                   goodwavel,goodflux,goodfluxerr,
                                   wavelim,fluxlim,upperlim,param_name,
                                   norm_arr,param_min,param_max,
                                   chain['params_std'][i,:,:],dr,temp[i]),
                             kwargs=kwargs)
                procs.append(p)
                p.start()
        
            #get chain results
            for p in procs:
                chain_result=chain_queue.get()
                ind=chain_result['chain_id']
                if synthesis_model != None:
                    sed=chain_result['sed']
                    if not 'sed' in chain.dtype.names:
                        dt=chain.dtype.descr+[('sed',float,sed.shape)]
                        chain=np.empty(chain.shape,dtype=dt)
                    chain['sed'][ind]=sed
                chain['params_iter'][ind,:,n*step_int:(n+1)*step_int]= \
                    chain_result['param_arr']
                chain['lnl_iter'][ind,n*step_int:(n+1)*step_int]= \
                    chain_result['lnl_arr']
                chain['lum_iter'][ind,n*step_int:(n+1)*step_int]= \
                    chain_result['lum_arr']
            
            #wait for processes to finish
            for p in procs:
                p.join()

            #if a swap was made with chain 0 (the fiducial 'cold' chain)
            #check accpetance/convergence and update if necessary
            if not 'no_temp' in kwargs:
                if (swap == 1) and (n > period):
                    dlnl_iter=chain['lnl_iter'][0,step_int-1:step_int*2-1]- \
                        chain['lnl_iter'][0,step_int:step_int*2]
                    whchange=abs(dlnl_iter[0:step_int-2]) > 0
                    pchange=np.float(sum(whchange))/(step_int-1)

                    if (pchange < 0.20) or (pchange > 0.26):
                        cov=chain['params_std'][0,:,:]
                        if pchange < 0.05: cov=cov*0.05
                        elif pchange > 0.26:
                            frac=pchange/0.23
                            cov=cov*frac
                        else:
                            cov=np.cov(np.float32(chain['params_iter'][0,:,n*step_int:(n+1)*step_int]))
                            
                        sel=((np.sqrt(np.diag(cov)) < 
                              np.sqrt(np.diag(params_std))*1.e-6) |
                             (np.sqrt(np.diag(cov)) > 
                              np.sqrt(np.diag(params_std))*1.e4))

                        if (np.array_equal(cov,np.zeros(cov.shape)) or 
                            sel.any()): 
                            cov=params_std
                        chain['params_std'][0,:,:]=cov
                        chain['params_iter'][0:1,:,n*step_int-1]= \
                            chain['params_iter'][0:1,:,(n+1)*step_int-1]
                        chain['lnl_iter'][0:1,n*step_int-1]= \
                            chain['lnl_iter'][0:1,(n+1)*step_int-1]
                        n-=1
                    else:
                        swap=0
                        if (n >= nstep/(step_int*2)):
                            chain['params_iter'][0:1,:,period*step_int-1]= \
                                chain['params_iter'][0:1,:,(n+1)*step_int-1]
                            chain['lnl_iter'][0:1,period*step_int-1]= \
                                chain['lnl_iter'][0:1,(n+1)*step_int-1]
                            n=period-1
                            lastswap=0
                            pchange=0

                    #ensure convergence is still maintained
                    if (pchange > 0.2) and (pchange < 0.26):
                        av1=np.average(chain['lnl_iter'][0,n*step_int: 
                                                         (n+0.1)*step_int])
                        av2=np.average(chain['lnl_iter'][0,(n+0.5)*step_int: 
                                                         (n+1)*step_int])
                        s=np.std(chain['lnl_iter'][0,(n+0.5)*step_int: 
                                                   (n+1)*step_int-1])
                        if s != 0:
                            if abs(av1-av2)/s > 2.0:
                                chain['params_iter'][:,:,n*step_int-1]= \
                                    chain['params_iter'][:,:,(n+1)*step_int-1]
                                chain['lnl_iter'][:,n*step_int-1]= \
                                    chain['lnl_iter'][:,(n+1)*step_int-1]
                                n-=1
                                swap=1

            #after running the temperd chains for a while and no swaps
            #are made to the 'cold' chain, reset all to 'cold' and 
            #run to completion.  This provides more samples around the
            #posterior
            if (not 'no_temp' in kwargs) and (n-lastswap == (nstep/step_int/2)):
                kwargs['no_temp']=True
                temp[:]=0
                for i in xrange(nchain):
                    maxindx=(n+1)*step_int-1
                    minindx=max(period,lastswap)*step_int
                    indx=np.random.choice(xrange(minindx,maxindx))

                    chain['params_iter'][i,:,period*step_int-1]= \
                        chain['params_iter'][0,:,indx]
                    chain['lnl_iter'][i,period*step_int-1]= \
                        chain['lnl_iter'][0,indx]
                    chain['params_std'][i,:,:]=chain['params_std'][0,:,:]
                n=period-1
            n+=1
            
            #create/update savefile
            np.save('parent_chain',chain)

    #keep track of only accepted steps
    lnl_good=np.empty(0,dtype=np.longdouble)
    params_good=np.empty((0,nparm+1))
    lum_good=np.empty(0)
    if synthesis_model != None: sed_good=[]

    for i in xrange(nchain):
        dlnl_iter=chain['lnl_iter'][i,1:]-chain['lnl_iter'][i,0:nstep-1]
        whchange=((abs(dlnl_iter) > 0) & (xrange(nstep) > step_int*2))
        lnl_good=np.concatenate((lnl_good,chain['lnl_iter'][i,whchange]))
        params_good=np.concatenate((params_good,chain['params_iter'][i,:,whchange]))
        lum_good=np.concatenate((lum_good,chain['lum_iter'][i,whchange]))
        if synthesis_model != None:
            sed_good.extend(((chain['sed'][i,whchange.nonzero()])[0]).tolist())

    dt=[('lnl_good',np.longdouble),('lum_good',np.longdouble),
        ('params_good',float,(nparm+1,))]
    if 'sed' in chain.dtype.names: dt+=[('sed_good',float,sed.shape)]
    chain_good=np.zeros(len(lnl_good),dtype=dt)
    chain_good['lnl_good']=lnl_good
    chain_good['lum_good']=lum_good
    if synthesis_model != None: chain_good['sed_good']=np.array(sed_good)

    #restore paramter space from compressed log-space
    if dr.any():
        params_good[:,dr]=10**(params_good[:,dr])
        param_max[dr]=10**param_max[dr]
        param_min[dr]=10**param_min[dr]
        for i in xrange(nsed):
            if seds[i].nparm > 1:
                gs=np.array([p in seds[i].param_name for p in param_name[:-1]])
                seds[i].grid[:,dr[gs]]=10**(seds[i].grid[:,dr[gs]])
    chain_good['params_good']=params_good
    
    #determine 68% confidence intervals
    maxl=np.amax(lnl_good)
    whmx=np.argmax(lnl_good)
    params_max=params_good[whmx,:]
    lum_max=lum_good[whmx]
    params_1sig=np.empty((nparm+1,2))
    sort=(np.sort(lnl_good))[::-1]
    ind68=np.floor(.682*len(lnl_good))
    lnL68=sort[ind68]
    wh1sig=lnl_good >= lnL68
    params_low=np.amin(params_good[wh1sig,:],axis=0)
    params_high=np.amax(params_good[wh1sig,:],axis=0)
    params_1sig[:,0]=params_high-params_max
    params_1sig[:,1]=params_max-params_low
    lum_1sig=[max(lum_good[wh1sig])-lum_max,lum_max-min(lum_good[wh1sig])]
    
    #plot the best fit model
    #for grid based, find the model closest to the maximum
    best_mods=np.empty(nsed,dtype='object')
    for i in xrange(nsed):
        if (i == 0) and synthesis_model != None:
            #grab best fit from synthesized SEDs
            sed=chain_good['sed_good'][whmx,:,:]
            if norm_arr[0]: sed[:,1]*=params_max[0]
            lnu=sed[:,1]**dist**2./(1.+z50)/(1.13e-8)/c
            best_mod={'lum':abs(np.trapz(lnu,c/sed[:,0])),'sed':sed}
        else:
            sel=np.array([p in seds[i].param_name for p in param_name])
            if sel.any():
                mod_best,ind,status=find_nearest_grid(seds[i].grid,
                                                      params_max[sel])
                ndim=sum(sel)
                grid=np.array([np.array(range(2**ndim))/2**j % 2 
                               for j in xrange(ndim)])
                mod_best=interpol.griddata(grid.T,mod_best,
                                           ind.reshape((1,ndim)),
                                           method='nearest')
            else: mod_best=np.atleast_1d(0)
            best_wave=seds[i].Template[mod_best[0]]['wavesim']
            best_flux=seds[i].Template[mod_best[0]]['fluxsim']
            best_lum=seds[i].Template[mod_best[0]]['lum']
            if norm_arr[i]:
                norm=np.array(param_name) == 'norm'+str(i+1)
                best_flux*=params_max[norm]
                best_lum*=params_max[norm]
                
            s_wave,s_flux=shift_sed(best_wave,best_flux,params_max[-1])
            best_mod={'lum':best_lum,'sed':np.column_stack((s_wave,s_flux))}
            
        best_mods[i]=best_mod

    if not 'no_plot' in kwargs:
        #set colors for SEDs
        rainbow=plt.get_cmap('gist_rainbow')
        cnorm=colors.Normalize(vmin=0,vmax=nsed)
        scalarmap=cmx.ScalarMappable(norm=cnorm,cmap=rainbow)

        #co-add SEDs
        flux_best=best_mods[0]['sed'][:,1]
        wave_best=best_mods[0]['sed'][:,0]
        
        for i in xrange(1,nsed):
            #determine overlap region and interpolate to corsest binning
            maxwave=min(np.amax(wave_best),np.amax(best_mods[i]['sed'][:,0]))
            minwave=max(np.amin(wave_best),np.amin(best_mods[i]['sed'][:,0]))
            sel=(wave_best <= maxwave) & (wave_best >= minwave)
            sels=((best_mods[i]['sed'][:,0] <= maxwave) & 
                  (best_mods[i]['sed'][:,0] >= minwave))

            if sum(sel) > sum(sels):
                mod_int=np.interp(best_mods[i]['sed'][sels,0],
                                  wave_best[sel],flux_best[sel])
                flux_int=best_mods[i]['sed'][sels,1]+mod_int
                wave_int=best_mods[i]['sed'][sels,0]
            else:
                mod_int=np.interp(wave_best[sel],best_mods[i]['sed'][sels,0],
                                  best_mods[i]['sed'][sels,1])
                flux_int=flux_best[sel]+mod_int
                wave_int=wave_best[sel]

            #check if there are points outside the interpolation region
            #and add them back to the co-added SED
            sel=wave_best < minwave
            if sel.any():
                wave_int=np.concatenate((wave_int,wave_best[sel]))
                flux_int=np.concatenate((flux_int,flux_best[sel]))
            sel=best_mods[i]['sed'][:,0] < minwave
            if sel.any():
                wave_int=np.concatenate((wave_int,best_mods[i]['sed'][sel,0]))
                flux_int=np.concatenate((flux_int,best_mods[i]['sed'][sel,1]))
            sel=wave_best > maxwave
            if sel.any():
                wave_int=np.concatenate((wave_best[sel],wave_int))
                flux_int=np.concatenate((flux_best[sel],flux_int))
            sel=best_mods[i]['sed'][:,0] > maxwave
            if sel.any():
                wave_int=np.concatenate((best_mods[i]['sed'][sel,0],wave_int))
                flux_int=np.concatenate((best_mods[i]['sed'][sel,1],flux_int))

            wave_best=wave_int
            flux_best=flux_int

        #plot co-added and component SEDs
        if flam == False:
            plt.plot(wave_best,flux_best,'k-',mfc='none')
            for i in xrange(nsed):
                colorval=scalarmap.to_rgba(i)
                plt.plot(best_mods[i]['sed'][:,0],best_mods[i]['sed'][:,1],
                         color=colorval,label='SED '+str(i+1),
                         linestyle='dotted')
            plt.legend(shadow=True, fancybox=True)
        else:
            lfl=flux_best*3.e-9/wave_best
            plt.plot(wave_best,lfl,'k-',mfc='none')
            for i in xrange(nsed):
                colorval=scalarmap.to_rgba(i)
                lfl=best_mods[i]['sed'][:,1]*3.e-9/best_mods[i]['sed'][:,0]
                plt.plot(best_mods[i]['sed'][:,0],lfl,color=colorval,
                         label='SED '+str(i+1),linestyle='dotted')
            plt.legend(shadow=True, fancybox=True)
        if 'psfile' in kwargs:
            plt.savefig(kwargs.get('psfile'), format='eps')

    #print results to screen (and outfile if set)
    logging.basicConfig(level=logging.INFO,format='%(message)s')
    if 'outfile' in kwargs:
        fh=logging.FileHandler(kwargs.get('outfile')+'.dat',mode='w')
        logging.getLogger('').addHandler(fh)
        
    logging.info('Best fit model parameters:')
    for i in xrange(nparm):
        logging.info('%s: %.4f+%.4f-%.4f' % \
                         (param_name[i],params_max[i],params_1sig[i,0],
                          params_1sig[i,1]))
    if findz:
        logging.info('%s: %.4f+%.4f-%.4f' % \
                         (param_name[-1],params_max[-1],params_1sig[-1,0],
                          params_1sig[-1,1]))
    else: logging.info('%s: %.4f' % (param_name[-1],params_max[-1]))
    logging.info('%s: %.3e+%.3e-%.3e' % \
                     ('Total luminosity (L_sun)',lum_max,lum_1sig[0],
                      lum_1sig[1]))
    logging.info('%s: %.3f' % ('Maximum Likelihood value',maxl))

    if 'outfile' in kwargs:
        savez_kwargs={'chain':chain,'chain_good':chain_good,'lum_max':lum_max,
                      'maxl':maxl,'params_max':params_max,'flux':flux,
                      'param_name':param_name,'step_int':step_int,
                      'param_max':param_max,'param_min':param_min,
                      'best_mods':best_mods,'wavel':wavel,'fluxerr':fluxerr,
                      'upperlim':upperlim}
        np.savez(kwargs.get('outfile'),**savez_kwargs)

    #clean up mid-level files
    for i in xrange(nchain):
        os.remove('chain_'+str(i)+'.log')
        os.remove('chain_'+str(i)+'.npy')
    os.remove('parent_burnin.npy')
    os.remove('parent_chain.npy')
   
               
