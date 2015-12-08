function satmc_random, d1,binomial=binomial,double=double,gamma=gamma, $
                       normal=normal,poisson=poisson,uniform=uniform, $
                       long=long
;Calls the IDL function RANDOMU without requiring a seed.  This 'should'
;alieviate reproducing the same 'random' numbers when re-using the same
;seed, particularly the IDL generic seed.  See description of RANDOMU
;for details on kewords
;

common satmc_seed, satmc_seed
if not keyword_set(satmc_seed) then satmc_seed=systime(/seconds)/!dpi

;see if it's just a scalar
if not keyword_set(d1) then d1 = 1
;return scalar if just want scalar, not array of 1
if d1 eq 1 then begin
   d1 = 1
   temp = randomu(satmc_seed,d1,binomial=binomial,double=double, $
                  gamma=gamma,normal=normal,poisson=poisson, $
                  uniform=uniform,long=long)
   return, temp[0]
endif else begin
   return, randomu(satmc_seed,d1,binomial=binomial,double=double, $
                   gamma=gamma,normal=normal,poisson=poisson, $
                   uniform=uniform,long=long)
endelse
end

function satmc_mrandom,cov,dl,double=double
;Based on the IDL function MRANDOMN, this function draws randomly
;distributed samples from a multivariate normal distribution.  The 
;difference here (aside from calls to SATMC_RANDOM instead of RANDOMN)
;is that this function uses SVD instead of Cholesky decomposition on 
;the covariance matrix COV.  Though SVD is slower, it does not run
;into the problems Cholesky does when there are 0 value elements 
;along the diagonal.
;
;The theory behind reconstructing a multivariate normal distribution 
;from its covariance matrix follows from geometric interpretation and
;reconstructing a standard Gaussian.  For a normal Gaussian N(mu,sigma)
;which has mean mu and standard deviation sigma, the expression
; mu+sqrt(sigma)N(0,1) returns the same distribution as N(mu,sigma),
;where N(0,1) is a Gaussian with 0 mean and stddev. 1.  Similarly,
;the multivariate normal distribution N(mu_vec,Sigma) with mean
;mu_vec=[mu_1,mu_2,...,mu_n] and covariance matrix Sigma can be
;reproduced by the expression mu_vec+sqrt(Sigma)*N(0,I) where I is the
;identiy matrix.  Here, we use SVD to construct sqrt(Sigma) since
;(by SVD) Sigma=U*w*V^T thus sqrt(Sigma)=U*sqrt(w).  This is true 
;since by definition the covariance matrix is positive semi-definite
;and the eigenvalue decomposition Sigma=U*w*U^T is also a SVD. 
;
; INPUTS:
; 
; COV - covariance matrix of the multivariate distribution
; DL - number of random variables to be generated from the 
;      multivariate distribution
;

if n_elements(dl) eq 0 then dl=1

;check that covariance matrix COV is a square matrix
if size(cov,/n_dim) ne 2 then begin
   message,/cont,'COV must be a matrix'
   return,0
endif
np=(size(cov))[1]
if (size(cov))[2] ne np then begin
   message,/cont,'COV must be a square matrix'
   return,0
endif

;perform the SVD
svdc,cov,w,u,v,double=double
a=diag_matrix(w)

;generate random variates
rand=make_array(dl,np,double=double)
for i=0,np-1 do rand[*,i]=satmc_random(dl,/normal,double=double)
rand=u##sqrt(a)##rand

return,rand
end

function shift_sed,model,redshift
;Simple routine to apply redshift correction to an observed SED
  redshift=redshift(0)
  dist=50.
  owave=model[*,0]*(1+redshift)
  ldist=lumdist(redshift,Lambda0=0.7,Omega_M=0.3,/silent)
  oflux=model[*,1]*(dist/ldist)^2.*(1+redshift)
  error=check_math(mask=32)     ;suppress floating underflow errors
  omodel=[[owave],[oflux]]
  return,omodel
end


function find_nearest_grid,grid,p,intp=intp,status=status
;locate given parameter values on a parameter grid
; INPUTS:
;  GRID - 2d array returned by GET_SED_PARAM containing parameter values
;         at each grid point
;  P - array of parameter values from current step
;
; OPTIONAL OUTPUTS:
;  INTP - array containing the location of each element of P in
;         respect to the bracketing grid points, used for interpolation
;  STATUS - returns -1 if P is not bound by GRID
  np=n_elements(p)
  ngrid=2^np
  status=0
  missing=0
  max_grid=fltarr(np)
  min_grid=fltarr(np)
  neargrid=fltarr(np+1,ngrid)
  grid=float(grid)
  t_grid=float(grid)

  ;for each p value, find the grid points in the respective grid dimension
  ;that bracket p
  for j=0,np-1 do begin
     ind=value_locate(grid[*,j],p[j])
     min_grid[j]=grid[(ind > 0),j]
     max_grid[j]=grid[((ind+1) < (n_elements(grid[*,0])-1)),j]
  endfor
  n_grid=[[min_grid],[max_grid]]
  intp=(p-min_grid)/(max_grid-min_grid)
  error=check_math()            ;ignore math errors

  ;check that the interpolation points lie within the grid
  whnan=where(finite(intp,/nan) eq 1,nwhnan)
  if nwhnan ne 0 then intp[whnan]=0.0 ;point lies on the grid, no problem
  whfin=where(finite(intp,/inf) eq 1,nwhinf)
  if nwhinf ne 0 then status=-1 ;point is outside the grid
  whinvalid=where(intp lt 0 or intp gt 1,nwhinv)
  if nwhinv ne 0 then status=-1 ;irregular gridding prevents bracketing
  if status eq -1 then return,neargrid
  
  ;get indicies of nearest grid points and mark their location
  if n_elements(grid[*,0]) gt 2 then begin
     for i=0,ngrid-1 do begin
        ind=where(grid[*,0] eq n_grid[0,(i mod 2)])
        grid=grid[ind,*]
        neargrid[1,i]=(i mod 2)
        for j=1,np-1 do begin
           g_ind=i/2^j mod 2
           neargrid[j+1,i]=g_ind
           sel=where(grid[*,j] eq n_grid[j,g_ind],nsel)
           if nsel eq 0 then begin
              missing=1
              break
           endif else begin
              grid=grid[sel,*]
              ind=ind[sel]
           endelse
        endfor
        if missing eq 1 then begin
           neargrid[0,i]=-1 
           missing=0
        endif else neargrid[0,i]=ind
        grid=t_grid
     endfor
     grid=t_grid
  endif
  
  whmissing=where(neargrid[0,*] eq -1,nmissing)
  if nmissing ne 0 then status=-1
  return,neargrid
end

function os_gauss,x,s,upperlim
;    Return the likelihood function as either a step function 
;    or one-sided Gaussian. 
;
;    x - Location of where to compute the likelihood function
;    s - Flux value of the upper limit
;    upperlim - String containing the type of upper limit prescription
;               to use and its defining parameters
;               1: Step function ('1').  Returns 1 if x < s and 
;                  0 if x > s.
;               2: One-sided Gaussian ('2,mu,sigma').  Returns 1 
;                  if x < mu and exp(-(x-mu)**2/(2*sigma**2)) 
;                  if x > mu.  If mu and/or sigma are not
;                  specified, assumes default values of mu=s/1.5
;                  and sigma=3*(s-break)/(5*s)

  npt=n_elements(x)
  pdf=x*0.
  for i=0L,npt-1 do begin
     uparm=strsplit(upperlim[i],',',/extract,/preserve_null)
     if fix(uparm[0]) eq 1 then begin
        if x[i] lt s[i] then pdf[i]=1 else pdf[i]=exp(-745d)
     endif
     if fix(uparm[0]) eq 2 then begin
        if n_elements(uparm) lt 2 then begin
           level=5
        endif else begin
           if uparm[1] eq '' then level=5 else level=float(uparm[1])
        endelse
        if n_elements(uparm) lt 3 then begin
           mu=s[i]*(3./level)
        endif else begin
           if uparm[2] eq '' then mu=s[i]*(3./level) else mu=float(uparm[2])
        endelse 
        if n_elements(uparm) lt 4 then begin
           sigma=(s[i]-mu)/level
        endif else begin
           if (uparm[3] eq '' or mu eq 0) then begin
              sigma=(s[i]-mu)/level 
           endif else sigma=float(uparm[3])
        endelse
        if x[i] lt mu then pdf[i]=1 else $ 
           pdf[i]=gaussian(x[i],[1,mu,sigma]) > exp(-745d)
     endif
  endfor
  return,pdf
end

pro get_sed_param,sed,grid,param_name,nparams,norm=norm,ngrid=ngrid, $
                  synthesis_model=synthesis_model,no_interp=no_interp
;Get SED information including number of parameters and parameter names
  tags=tag_names(sed)
  if not keyword_set(synthesis_model) then begin
      params=where(tags ne 'MODELNAME' and tags ne 'WAVESIM' and $
                   tags ne 'FLUXSIM' and tags ne 'LUM',nparams,ncomp=ndeftags)
      if ndeftags lt 3 then begin
          message,/cont,'Warning: SED has invalid structure format. Verify'+$
            ' tag names WAVESIM, FLUXSIM, (LUM) are present.'
          return
      endif
  endif else begin
     nparams=0
     params=[-999]
     for i=0,n_tags(sed)-1 do begin
        if n_elements(sed.(i)) gt 1 then begin
           nparams++
           params=[params,i]
        endif
     endfor
     params=params[1:*]
  endelse

  ;define parameter grid
  if nparams gt 0 then begin
     if not keyword_set(synthesis_model) then begin
        grid=sed.(params[0])
        for i=1,nparams-1 do grid=[[grid],[sed.(params[i])]]
     endif else begin
        grid=fltarr(2,nparams)
        for i=0,nparams-1 do grid[*,i]=sed.(params[i])
     endelse
     if keyword_set(norm) then begin
        param_name=['Norm',tags(params)]
        nparams++
     endif else param_name=tags(params)
  endif else begin
     grid=[0]
     if keyword_set(norm) then begin
        param_name=['Norm']
        nparams++
     endif else param_name=['']
     if keyword_set(no_interp) and n_elements(sed) gt 1 then begin
        param_name=[param_name,'Model']
        nparams++
        grid=indgen(n_elements(sed))
     endif
  endelse
  ngrid=[-999]
  for i=0,n_elements(grid[0,*])-1 do ngrid=[ngrid,n_elements(uniq(grid[sort(grid[*,i]),i]))]
  ngrid=ngrid[1:*]
  return
end

pro cal_lum,struct,xrange=xrange
  ;calculate observed bolometric luminosity given input model SED
  ;assume flux in units of Jy and wavelength in microns
  c=2.9979e14                     ;(microns/s)  
  lum=dblarr(n_elements(struct))
  tmp=create_struct(struct[0],'lum',0.0d)
  tmp=replicate(tmp,n_elements(struct))
  ;assume D(z)=50 Mpc
  z50=70.*50/3.e5
  for i=0L,n_elements(struct)-1 do begin
      lnu=struct[i].fluxsim*50.^2./(1.+z50)/(1.13e-8)/c
      lum[i]=abs(integral(c/struct[i].wavesim,lnu,xrange=xrange))
  endfor
  struct_assign,struct,tmp
  tmp.lum=lum
  struct=tmp
  return
end

function get_model_observed,model,goodwavel,wavelim, $
                            filter_struct=filter_struct,mbad=mbad
;Extract model observed fluxes from redshifted SED using either
;observed wavelengths or filters
;
  c=2.9979e14                   ;microns/s
  whmod=where(goodwavel ge min(model[*,0]) and goodwavel le max(model[*,0]))
  mgood=goodwavel*0.-1.
  if n_elements(wavelim) ne 0 then mbad=wavelim*0.
        
  if n_elements(filter_struct) gt 1 then begin
     ;use filters for extracting fluxes
     good=where(filter_struct.flag eq 1)
     fwave=filter_struct[good].wave
     fresponse=filter_struct[good].response
     norm=filter_struct[good].norm
     fsel=interpol(model[*,1],model[*,0],fwave) > 0
     for i=0,n_elements(whmod)-1 do begin
        if not array_equal(fwave[*,i],0) then begin
           sel=where(fwave[*,i] ne 0 and fsel[*,i] gt 0)
           wsel=fwave[sel,i]
           rsel=fresponse[sel,i]
           h=6.626e-27
           mgood[whmod[i]]=integral(c/wsel,fsel[sel,i]* $
                                    rsel/(h*c/wsel))/norm[i]
        endif else begin
           mgood[whmod[i]]=interpol(model[*,1],model[*,0],goodwavel[whmod[i]])
        endelse
     endfor
  endif else begin
     ;extract fluxes directly at observed wavelengths
     mgood[whmod]=interpol(model[*,1],model[*,0],goodwavel[whmod])
  endelse
  if n_elements(wavelim) ne 0 then begin
     whmod=where(wavelim ge min(model[*,0]) and wavelim le max(model[*,0]),nmod)
     if nmod ne 0 then begin
        if n_elements(filter_struct) gt 1 then begin
           ;use filters for extracting fluxes
           bad=where(filter_struct.flag eq 0)
           fwave=filter_struct[bad].wave
           fresponse=filter_struct[bad].response
           norm=filter_struct[bad].norm
           fsel=interpol(model[*,1],model[*,0],fwave) > 0
           for i=0,n_elements(whmod)-1 do begin
              if not array_equal(fwave[*,i],0) then begin
                 sel=where(fwave[*,i] ne 0 and fsel[*,i] gt 0)
                 wsel=fwave[sel,i]
                 rsel=fresponse[sel,i]
                 h=6.626e-27
                 mbad[whmod[i]]=integral(c/wsel,fsel[sel,i]*$
                                         rsel/(h*c/wsel))/norm[i]
              endif else begin
                 mbad[whmod[i]]=interpol(model[*,1],model[*,0], $
                                         wavelim[whmod[i]])
              endelse
           endfor
        endif else begin
           ;extract fluxes directly at observed wavelengths
           mbad[whmod]=interpol(model[*,1],model[*,0],wavelim[whmod])
        endelse
     endif
  endif
  return,mgood
end


pro cal_lgrid,sed_struct,grid_struct,grid,lgrid,wave,flux,$
              fluxerr,upperlim,norm_arr,nparm_arr,param_max,$
              param_min,param_name=param_name,minz=minz,maxz=maxz,$
              dr=dr,flux_grid=flux_grid,priors=priors,prior_mode=prior_mode
;Pre-compute likelihoods for grid-based templates.  In the case of
;photo-z estimate and/or normalizations, a corse grid is assumed.
;The resulting grid is used to determine initial points for photo-z
;estimates
;

  c=2.9979e14
 
  ;setup grids to hold likelihoods and corresponding parameters
  seddim=[1]
  nsed=n_tags(sed_struct)
  nparams=total(nparm_arr)
  for l=0,nsed-1 do begin
     if norm_arr[l] eq 1 then seddim=[seddim,31]
     seddim=[seddim,n_elements(sed_struct.(l))]
  endfor
  seddim=seddim[1:*]
  nz=round((maxz-minz)/0.1); < 30
  zgrid=dindgen(nz+1)/nz*(maxz-minz)+minz
  seddim=[seddim,nz+1]
  ngrid=product(seddim)
  
  grid=dblarr(ngrid,nparams+1)
  lgrid=dblarr(ngrid)
  lumgrid=dblarr(ngrid)

  ;fill in parameter grid
  indx=array_indices(seddim,lindgen(ngrid),/dim)
  if n_elements(seddim) gt 1 then indx=transpose(indx)
  j=0
  for l=0,nsed-1 do begin
     p0=nparams-total(nparm_arr[l:*])
     if norm_arr[l] eq 1 then begin
        normgrid=dindgen(31)/30.*(param_max[p0]-param_min[p0])+ $
                 param_min[p0]
        grid[*,p0]=normgrid[indx[*,j]]
        grid[*,p0+1:p0+nparm_arr[l]-1]=grid_struct.(l)[indx[*,j+1],*]
        lumgrid[*]+=normgrid[indx[*,j]]*sed_struct.(l)[indx[*,j+1],*].lum
        j+=2
     endif else begin
        grid[*,p0:p0+nparm_arr[l]-1]=grid_struct.(l)[indx[*,j],*]
        lumgrid[*]+=sed_struct.(l)[indx[*,j],*].lum
        j++
     endelse
  endfor
  zp=zgrid[indx[*,j]]
  grid[*,nparams]=zp

  ;get flux_grid if needed
  if not keyword_set(flux_grid) then begin
     ldist=lumdist(zgrid,Lambda0=0.7,Omega_M=0.3,/silent)
     flux_grid=dblarr(n_elements(wave),ngrid)

     for k=0L,ngrid-1 do begin
        indx=array_indices(seddim,k,/dim)
        j=0
        
        for l=0,nsed-1 do begin
           p0=nparams-total(nparm_arr[l:*])
           mw=sed_struct.(l)[0].wavesim
           if norm_arr[l] eq 1 then begin
              if where(dr eq p0) ne -1 then begin
                 mf=sed_struct.(l)[indx[j+1]].fluxsim*10^(grid[k,p0])
                 err=check_math()
              endif else begin
                 mf=sed_struct.(l)[indx[j+1]].fluxsim*grid[k,p0]
              endelse
              j+=2
           endif else begin
              mf=sed_struct.(l)[indx[j]].fluxsim
              j++
           endelse
           zind=indx[j]
           mws=mw*(1+grid[k,nparams])
           mfs=mf*(50./ldist[zind])^2.*(1+grid[k,nparams])

           v=value_locate(mws,wave)
           int=(wave-mws[v])/(mws[v+1]-mws[v])
           flux_grid[*,k]+=interpolate(mfs,v+int)
        endfor
     endfor
  endif

  good=where(upperlim eq 0,n_good,complement=bad,ncomplement=n_bad)
  goodflux=rebin(flux[good],n_good,ngrid)
  goodfluxerr=rebin(fluxerr[good],n_good,ngrid)
  if n_bad ne 0 then fluxlim=rebin(flux[bad],n_bad,ngrid)

  lgrid=total(-(goodflux-flux_grid[good,*])^2./(2.*goodfluxerr^2.),1)
  if not keyword_set(ig_lim) and n_bad ne 0 then begin
     if n_bad eq 1 then begin
        lgrid+=alog(os_gauss(flux_grid[bad,*],fluxlim,upperlim[bad]))
     endif else lgrid+=total(alog(os_gauss(flux_grid[bad,*],fluxlim, $
                                           upperlim[bad])),1)
  endif
  ;check with priors
  if keyword_set(priors) then begin
     l_priors=where(prior_mode eq 2,nsel)
     if nsel ge 1 then begin
        for i=0,nsel-1 do begin
           comm=priors[l_priors[i]]+',sed=sed_struct,params='+ $
                'grid,param_name=param_name,lnl=lgrid,lum=lum_arr'
           null=execute(comm)
        endfor
     endif
  endif

  return
end

function md_interpolate,array,loc,intp
; return multi-linear interpolation of array based on grid locations loc 
; and distances intp

  nint=n_elements(intp)

  val=array
  for j=0,nint-1 do val*=abs(1-intp[j]-loc[j,*])
  return,total(val)
end

function verify_parameters,params,params_old,params_std,param_name,grid, $
                           param_min,param_max,priors=priors, $
                           prior_mode=prior_mode
;check input parameters are bound by parameter grid and additional
;priors

;first bound by parameter grid
sel=where(params lt param_min or params gt param_max,nsel)
if nsel ne 0 then begin
   params[sel]=(float(satmc_mrandom(params_std,/double)))[sel]+ $
               params_old[sel]
   status=-1
   goto,return
endif

p0=0
nsed=n_elements(tag_names(grid))
for l=0,nsed-1 do begin
   sed_name=gettok(param_name[p0],' ')
   norm_pos=where(strupcase(param_name) eq sed_name+' NORM')
   sel=where(strpos(param_name,sed_name) eq 0)
   if norm_pos ne -1 then sel=sel[1:*]
   null=find_nearest_grid(grid.(l),params[sel],status=status,intp=intp)
   if status eq -1 then begin
      bsel=where(finite(intp) eq 0,nbad)
      if nbad ne 0 then begin
         whbad=sel[bsel]
         params[whbad]=(float(satmc_mrandom(params_std,/double)))[whbad]+ $
                       params_old[whbad]
      endif else begin
         params[sel]=(float(satmc_mrandom(params_std,/double)))[sel]+ $
                     params_old[sel]
      endelse
      goto,return
   endif
   p0=sel[n_elements(sel)-1]+1
endfor

;check against parameter priors
if keyword_set(priors) then begin
   p_prior=where(prior_mode eq 1,nsel)
   if nsel ne 0 then begin
      for i=0,nsel-1 do begin
         comm=priors[p_prior[i]]+',params=params,param_name=param_name,'+$
              'status=status'
         null=execute(comm)
         if status eq -1 then begin
            params=float(satmc_mrandom(params_std,/double))+params_old
            goto,return
         endif
      endfor
   endif
endif

return:
return,status
end

pro run_mcmc_chain,nstep,params_init,lnl_init,sed_struct,grid_struct, $
                   goodwavel,goodflux,goodfluxerr,wavelim,fluxlim, $
                   upperlim,param_name,norm_arr,nparm_arr,param_min,param_max, $
                   params_std,param_arr,lnl_arr,tot_lum_arr, $
                   sed,models,synthesis_model=synthesis_model, $
                   outfile=outfile,temp=temp,no_interp=no_interp, $
                   filter_struct=filter_struct,dr=dr,ig_lim=ig_lim, $
                   priors=priors,prior_mode=prior_mode
;Take NSTEP number of steps in a single Monte Carlo Markov Chain.  
;For each step,compute likelihoods, luminosities and other relevant 
;information then write current chain information to file.

param_step=fltarr(n_elements(params_init),nstep)
param_arr=fltarr(n_elements(params_init),nstep)
lnl_step=dblarr(nstep)
lnl_arr=dblarr(nstep)
tot_lum_arr=dblarr(nstep)
sed=dblarr(nstep)

sed_arr=tag_names(sed_struct)
nsed=n_tags(sed_struct)
nparams=total(nparm_arr)
n_good=n_elements(goodwavel)
n_bad=n_elements(wavelim)
lum_arr=dblarr(nsed)
fir_lum_arr=dblarr(nsed)
c=2.9979e14                     ;(microns/s)

for n=0,nstep-1 do begin
   mod_intp=[-999]
   param0=0
    
   ;Make a proposal step
   if n eq 0 and lnl_init eq 0 then begin
      params_new=params_init 
      params_old=params_new
   endif else begin
      if n eq 0 then params_old=params_init else $
         params_old=param_arr[*,n-1]
      params_new=float(satmc_mrandom(params_std,/double))+params_old
   endelse
   ;check that parameters are within their given limits
   ;if not, pick a different step and check it again
   status=-1
   while status eq -1 do begin
      status=verify_parameters(params_new,params_old,params_std, $
                               param_name,grid_struct,param_min,param_max, $
                               priors=priors,prior_mod=prior_mode)
      ;print,params_new,status
   endwhile
   
   param_step[*,n]=params_new
   
   for l=0,nsed-1 do begin
      if l eq 0 and keyword_set(synthesis_model) then begin
         ;get tag names for SED synthesis routine
         synth_tag_name=tag_names(sed_struct.(0))
         ;set respective parameters in current IDL level
         for i=0,n_elements(synth_tag_name)-1 do begin
            psel=strpos(param_name,synth_tag_name[i]) 
            sel=where(psel ne -1,nsel)
            if nsel ne 0 then begin
               if where(dr eq sel[0]) ne -1 then begin
                  (scope_varfetch(synth_tag_name[i],/enter))= $
                     10^(params_new[sel])
               endif else begin
                  (scope_varfetch(synth_tag_name[i],/enter))=params_new[sel]
               endelse
            endif else begin
               (scope_varfetch(synth_tag_name[i],/enter))= $
                  sed_struct.(0).(i)
            endelse
         endfor
            
         ;make call to SED synthesis routine, passing along
         ;required keywords
         comm=synthesis_model+','+(scope_varname(models,level=-1))[0]
         for i=0,n_elements(synth_tag_name)-1 do $
            comm+=','+synth_tag_name[i]+'='+synth_tag_name[i]
         print,comm
         null=execute(comm)
            
         mod_intp=[mod_intp,0]  ;for consistency

         ;calculate luminosities
         models=(scope_varfetch(scope_varname(models,level=-1)))
         wave=models.wavelength
         flux=models.total
         lnu=flux*50.^2./(1.0118)*(3.124e-7)
            
         lum_arr[0]=abs(integral(c/wave,lnu)) ;bolometric
         ;40-120 micron
         fir_lum_arr[0]=abs(integral(c/wave,lnu,xrange=[c/120.,c/40.]))

         ;correct SED for redshift
         shift=shift_sed([[wave],[flux]],params_new[nparams])
         int=create_struct('loc',0,'wavel',wave,'flux',flux,'good', $
                              dblarr(n_good),'selg',intarr(n_good))
         if n_bad ne 0 then begin
            int=create_struct(int,'bad',dblarr(n_bad),'selb',intarr(n_bad))
         endif
                                
         int.good=get_model_observed(shift,goodwavel,wavelim, $
                                     filter_struct=filter_struct,mbad=mbad)
         
         sel=where(int.good ge 0,nsel)
         if nsel gt 0 then int.selg[sel]=1 else int.selg[*]=1
         if n_bad ne 0 then begin
            int.bad=mbad
            sel=where(mbad ge 0,nsel)
            if nsel gt 0 then int.selb[sel]=1 else int.selb[*]=1
         endif
         
         param0=nparm_arr[0]
         if n_elements(wave) gt 2 then begin
            if size(sed,/n_dim) eq 1 then begin
               sed=rebin(shift,[size(shift,/dim),nstep])
               sed[*]=0
            endif
            sed[*,*,n]=shift
         endif
         
         int_struct=create_struct(sed_arr[0],int)
         delvarx,models,wave,flux,lnu,shift
      endif else begin
         if norm_arr[l] then begin
            if nparm_arr[l] gt 1 then begin
               int_grid=find_nearest_grid(grid_struct.(l),params_new[param0+1:param0+nparm_arr[l]-1],intp=intp,status=status)
            endif else begin
               int_grid=[0,0]
               intp=0
               status=1
            endelse
         endif else begin
            int_grid=find_nearest_grid(grid_struct.(l),params_new[param0:param0+nparm_arr[l]-1],intp=intp,status=status)
         endelse
         n_intp=n_elements(intp)
         if keyword_set(no_interp) then begin
           for i=0,n_intp-1 do $
               int_grid=int_grid[*,where(int_grid[i+1,*] eq round(intp[i]))]
            intp=0
         endif
         mod_intp=[mod_intp,intp]
         ;re-check that the current step lies in model grid space
         ;a fail-safe stop
         if status eq -1 then begin
            print,status
            stop
         endif
         npt=n_elements(sed_struct.(l)[0].wavesim)
         
         int=create_struct('loc',dblarr(n_intp),'wavel', $
                               dblarr(npt),'flux',dblarr(npt),'good', $
                               dblarr(n_good),'selg',intarr(n_good))
	 if n_bad ne 0 then begin
             int=create_struct(int,'bad',dblarr(n_bad),'selb',intarr(n_bad))
         endif
         int=replicate(int,n_elements(int_grid[0,*]))
                   
         ;extract SED information from models
         pos=reform(int_grid[0,*])
         int.wavel=sed_struct.(l)[pos].wavesim
         if norm_arr[l] then begin
            ;apply any normalization factor
            if where(dr eq param0) ne -1 then begin
               int.flux=sed_struct.(l)[pos].fluxsim*10^(params_new[param0])
               lum=sed_struct.(l)[pos].lum*10^(params_new[param0])
            endif else begin
               int.flux=sed_struct.(l)[pos].fluxsim*params_new[param0]
               lum=sed_struct.(l)[pos].lum*params_new[param0]
            endelse
         endif else begin
            int.flux=sed_struct.(l)[pos].fluxsim
            lum=sed_struct.(l)[pos].lum
         endelse
         if keyword_set(no_interp) then begin
            int.loc=int_grid[1:*]
         endif else int.loc=int_grid[1:*,*]
         lum_arr[l]=md_interpolate(lum,int.loc,intp)
         
         for j=0,n_elements(int_grid[0,*])-1 do begin
            shift=shift_sed([[int[j].wavel],[int[j].flux]], $
                            params_new[nparams])
            int[j].good=get_model_observed(shift,goodwavel,wavelim, $
                                           filter_struct=filter_struct,$
                                           mbad=mbad)
            sel=where(int.good ge 0,nsel)
            if nsel ne 0 then int[j].selg[sel]=1
            if n_bad ne 0 then begin
               sel=where(mbad ge 0,nsel)
               if nsel ne 0 then begin
                  int[j].bad=mbad
                  int[j].selb[sel]=1
               endif
            endif
         endfor
         param0+=nparm_arr[l]
         
         if l eq 0 then int_struct=create_struct(sed_arr[l],int) else $
            int_struct=create_struct(int_struct,sed_arr[l],int)
      endelse
   endfor
   mod_intp=mod_intp[1:*]
   
   sel=where(strpos(param_name,'Norm') ne -1 or param_name eq 'z',ncomp=ndim)
   if keyword_set(synthesis_model) then ndim-=nparm_arr[0]
   mod_tot=create_struct('loc',dblarr(ndim > 1),'good',dblarr(n_good),'selg', $
                         intarr(n_good))
   if n_bad ne 0 then begin
      mod_tot=create_struct(mod_tot,'bad',dblarr(n_bad),'selb',intarr(n_bad))
   endif
   mod_tot=replicate(mod_tot,2^ndim)
   
   ;combine flux measurements for different models
   pdim=1
   p0=0
   for k=0,nsed-1 do begin
      nmod=n_elements(mod_tot)
      n_grid=n_elements(int_struct.(k))
      nparm=n_elements(int_struct.(k)[0].loc)
      for i=0,nmod-1 do begin
         mod_tot[i].loc[p0:p0+nparm-1]=int_struct.(k)[i/pdim mod n_grid].loc
         mod_tot[i].good=mod_tot[i].good+ $
                         int_struct.(k)[i/pdim mod n_grid].good
         mod_tot[i].selg=mod_tot[i].selg+ $
                         int_struct.(k)[i/pdim mod n_grid].selg
         if n_bad ne 0 then begin
            mod_tot[i].bad=mod_tot[i].bad+ $
                           int_struct.(k)[i/pdim mod n_grid].bad
            mod_tot[i].selb=mod_tot[i].selb+ $
                            int_struct.(k)[i/pdim mod n_grid].selb
         endif
      endfor
      pdim=pdim*n_grid
      p0+=nparm
   endfor

   ;now determine likelihood at each possible combination of models
   lnL=dblarr(n_elements(mod_tot))
   for i=0,n_elements(lnL)-1 do begin
      ;exclude observed data points outside SED models
      whgood=where(mod_tot[0].selg gt 0)
      gf=goodflux[whgood]
      gw=goodwavel[whgood]
      gfe=goodfluxerr[whgood]
      if n_bad ne 0 then begin
        whbad=where(mod_tot[0].selb gt 0,nwhbad)
        if nwhbad ne 0 then begin
            wl=wavelim[whbad]
            fl=fluxlim[whbad]
        endif else ig_lim=1
      endif

      lnL[i]=total(-(gf-mod_tot[i].good[whgood])^2./(2.*gfe^2.))
      ;include upper limits?
      if not keyword_set(ig_lim) and n_bad ne 0 then begin
         sel=where(upperlim ne 0)
         lnL[i]+=total(alog(os_gauss(mod_tot[i].bad[whbad],fl, $
                                     upperlim[sel[whbad]])))
      endif
      ;check if potential step exists, only happens for synthesized SEDs
      if array_equal(mod_tot.good,0) then lnL[i]=-1.d100
   endfor
   
   ;preform N-d interpolation over model grid to get likelihood of 
   ;current step
   if ndim gt 0 then begin
      l=lnL
      maxl=max(l)
      l-=maxl
      l=exp(l)
      l=md_interpolate(l,mod_tot.loc,mod_intp)
      lnL=alog(l)+maxl
   endif

   if keyword_set(priors) then begin
      l_priors=where(prior_mode eq 2,nsel)
      if nsel ge 1 then begin
         for i=0,nsel-1 do begin
            if keyword_set(synthesis_model) then begin
               comm=priors[l_priors[i]]+',sed=sed[*,*,n],params='+ $
                    'params_new,param_name=param_name,lnl=lnl,lum=lum_arr'
            endif else begin
               comm=priors[l_priors[i]]+',sed=sed_struct,params='+ $
                    'params_new,param_name=param_name,lnl=lnl,lum=lum_arr'
            endelse
            print,comm
            null=execute(comm)
         endfor
      endif
   endif

   lnl_arr[n]=lnL
   lnl_step[n]=lnL
   tot_lum_arr[n]=total(lum_arr)

   ;decide whether or not to take the step
   if not (n eq 0 and lnl_init eq 0) then begin
      if n eq 0 then $
         alpha=min([1,exp((lnl_arr[n]-lnL_init)/(1.+temp))]) else $
            alpha=min([1,exp((lnl_arr[n]-lnl_arr[n-1])/(1.+temp))])
      error=check_math()        ;suppress math errors
      urand=satmc_random(1,/uniform)
      if urand le alpha then begin
         param_arr[*,n]=params_new
      endif else begin
         param_arr[*,n]=params_old
         if n eq 0 then lnl_arr[n]=lnl_init else $
            lnl_arr[n]=lnl_arr[n-1]
      endelse
   endif else param_arr[*,n]=params_new

endfor
message,/cont,"Chain complete."

;create save file to store chain results
save,lnl_arr,param_arr,tot_lum_arr,lnl_step,param_step,filename=outfile

return
end

function create_children,nchain,savefile,synthesis_model=synthesis_model, $
                         wavelim=wavelim,fluxlim=fluxlim,upperlim=upperlim, $
                         dir=dir,no_temp=no_temp,no_interp=no_interp, $
                         ig_lim=ig_lim,priors=priors,prior_mode=prior_mode
;Create child IDL processes using IDL_IDLBRIDGE for each chain to
;enable use of multiple processors
  child=objarr(nchain)
  for m=0,nchain-1 do begin
     child[m]=Obj_New('IDL_IDLBridge',output=dir+'child'+strtrim(m,2)+'.log')
     child[m]->Execute,'@'+PREF_GET('IDL_STARTUP')
     child[m]->Execute,"cd,'"+dir+"'"
     ;pass currently set keywords to child IDL processes
     if keyword_set(synthesis_model) then  $
       child[m]->SetVar,'synthesis_model',synthesis_model
     if keyword_set(ig_lim) then child[m]->SetVar,'ig_lim',ig_lim
     if keyword_set(no_temp) then child[m]->SetVar,'temp',0 else $
        child[m]->SetVar,'temp',10.^(3.5*m/nchain)-1.
     if keyword_set(no_interp) then child[m]->SetVar,'no_interp',no_interp
     if n_elements(wavelim) ne 0 then begin
        child[m]->SetVar,'wavelim',wavelim
        child[m]->SetVar,'fluxlim',fluxlim
        child[m]->SetVar,'upperlim',upperlim
     endif
     if keyword_set(priors) then child[m]->SetVar,'priors',priors
     if keyword_set(prior_mode) then child[m]->SetVar,'prior_mode',prior_mode
          
     child[m]->execute,"resolve_routine,'satmc',/compile_full"
     child[m]->Execute,"restore,'"+savefile
  endfor
  return,child
end

pro satmc,fname,sed1,sed2,sed3,norm1=norm1,norm2=norm2,norm3=norm3, $
          redshift=redshift,findz=findz,outfile=outfile,psfile=psfile, $
          no_plot=no_plot,xrange=xrange,ig_lim=ig_lim, $
          nstep=nstep,step_int=step_int,nchain=nchain, $
          params_init=params_init,params_std=params_std,delz=delz, $
          param_max=param_max,param_min=param_min,minz=minz,maxz=maxz, $
          adaptive_step=adaptive_step,synthesis_model=synthesis_model, $
          no_async=no_async,restart=restart,no_temp=no_temp, $
          no_interp=no_interp,filters=filters,use_xvfb=use_xvfb, $
          flux_grid=flux_grid,cal_lgrid=cal_lgrid,priors=priors, $
          prior_mode=prior_mode,axes=axes,save_all=save_all
;    Using user supplied SED template(s), this procedure takes the 
;    observed SED given for a source and determines the set of best-fit 
;    parameters with '68%' confidence levels.  Additional constraints may
;    be imposed by including parameter/likelihood priors.  Descriptions of 
;    required and optional keywords follows. 
;
;WARNING: 
;    The default behavior requires a completely up-to-date install of IDL 
;    (tested using IDL 6.*/7.*).  If you encounter errors or suspect 
;    something may be wrong, check your IDL version and associated linked 
;    libraries.
;
;REQUIRED KEYWORDS:
;  FNAME - [string] the filename containing the fluxes in different
;          bands for a single source.  The file should be formated 
;          with the following columns:
;          1. wavelength of observing bands (available units: 
;             log10(Hz) (default), micron, Angstrom)
;          2. observed flux density (Jy, mJy (default), uJy)
;          3. error on the flux density measurements (Jy, mJy (default), uJy)
;          4. quality flag (0 if flux density is not an upper limit, 1
;             to treat the upper limit as a step function, 2 to apply a 
;             one-sided Gaussian, see documentation for more details)
;          If a format other than the default is used, the appropriate format
;          should be specified on the first line.
;  SED1 - [structure] IDL structure containing a set of SEDs to fit to 
;         the data.  The structure should have the following tags:
;         (MODELNAME) - [optional] String name of the current SED (not used)
;         PARAMS - Each model parameter should be given its own tag with
;                  values of the parameter corresponding to each SED listed
;         WAVESIM - n-element array containing the wavelength grid of the SED
;         FLUXSIM - n-element array containing the flux grid corresponding to 
;                   WAVESIM
;         (LUM) - [optional] bolometric luminosity (in L_sun) for current model
;         For examples of how to format a library of SEDs to conform with the 
;         required settings, see the README file.
;
;OPTIONAL KEYWORDS:
;
;  SED2,SED3 - [structure] IDL structures containting additional sets of 
;              SED models that will be combined to SED1 to create a
;              composite SED set.  Should following the same formating
;              as listed above.
;  NORM1,NORM2,NORM3 - [float array] if set, normalization is assumed
;                      to be an additional free parameter with
;                      NORM*[0] and NORM*[1] setting the lower and
;                      upper limits respecitvely.  If NORM* is a
;                      single element array, NORM* is assumed to be
;                      the upper limit with a lower limit of 1.0. 
;  REDSHIFT - [float] redshift of the source.  If the keyword FINDZ
;             is set, REDSHIFT is assumed to be the initial redshift
;             guess.  If REDSHIFT is not set, then FINDZ=True.
;  PSFILE - [boolean] if set, produces a postscript file of name
;           (psfile).ps of the observed and best-fit SED
;  OUTFILE - [string] output file name that will contain the results
;            of SED fitting.  Produces an ascii text file and IDL save
;            file of the name (outfile).fit and (outfile).sav, respectively.
;  NO_PLOT - [boolean] if set, suppresses plotting of data
;  AXES - [integer] Integar value specifying format of axes on output SED 
;         plot.  Acceptable values are:
;         1 - S_v (Jy) vs wavelength (micron) (default)
;         2 - lambda*F_lambda vs wavelength (micron) (note this is the 
;             same as nu*F_nu)
;  XRANGE - [array] Upper and lower limits for the output plot, same
;           as XRANGE in PLOT.
;  IG_LIM - [boolean] if set, ignore any flux upper limits in the 
;           input SED file
;  FINDZ - [boolean] if set, will search over a range of
;          redshifts.  If REDSHIFT is not specified, initial guesses will
;          be chosen at random.
;  NCHAIN - [integer] number of markov chains to run (def =5)
;  NSTEP - [long] number of iterations to keep for each markov chain 
;             (def =10000L)
;  STEP_INT - [long] each NSTEP chain will be separated into smaller,
;             more managable segments of length STEP_INT.  The length of 
;             a single burn-in iteration is set by 2*STEP_INT.  Save files
;             will be updated after every STEP_INT steps.
;  PARAMS_INIT - [double array] array of initial parameter guesses. 
;                Defualt is to randomly determine PARAMS_INIT from
;                SED parameter range.
;  PARAMS_STD - [NxN double array] covariance matrix of the proposal
;               parameter distribution to draw potential steps from.  If 
;               not supplied, will assume a default diagonal covariance
;               matrix based on parameter limits. 
;  PARAMS_MIN - [double array] array of parameter lower limits.  Will
;               change the default behavior and only include models
;               with parameter values above PARAMS_MIN
;  PARAMS_MAX - [double array] array of parameter upper limits similar
;               to PARAMS_MIN
;  MINZ,MAXZ,DELZ - [double] redshift lower limit, upper limit and
;                   step size to be used when performing the MC
;                   chains.  If the keyword FINDZ is not set, these
;                   are ignored.  Def=0.0,6.0,0.5 
;  ADAPTIVE_STEP - [boolean] sets whether or not to adaptively change
;                  the parameter step size when performing the MC
;                  chains. Defualt=true
;  SYNTHESIS_MODEL - [string] string name of SED synthesis routine.
;                    If set, SED1 is taken to be an input structure 
;                    with tag names corresponding to inputs to the
;                    SED synthesis routine.  The values assocated with 
;                    each tag should either be in [min,max] or [val] 
;                    format.  Models will be generated at potential
;                    steps given the synthesis routine and parameter
;                    ranges/values set by SED1. For an example, see
;                    GRASIL.pro for more information.  Note that NORM
;                    should NOT be set in combination with SYNTHESIS_MODEL.
;  NO_ASYNC - [boolean] if set, modifies the default behavior of creating 
;             child IDL processes for each chain.  Instead, each step of
;             each chain is handled one at a time.  Use this if your IDL 
;             setup does not support IDL_IDLBridge objects.
;  USE_XVFB - [boolean] if set, child processes under asyncronous operation
;             will be created under a virtual X session.  Useful if running
;             in the background or through virtual terminals (i.e. screen).
;             Note that Xvfb must be installed for this option to work.
;  RESTART - [integar] if set, the MCMC will restart at the given position, 
;            restoring the save files generated from a previous run and
;            continue. Accepted values are:
;                1- restart during burn-in period and iterations over 
;                   step size
;                2- restart during full run of the MCMC
;  NO_TEMP - [boolean] flag to turn off parallel tempering of chains, 
;           requires at least 2 chains (default =false).   
;  NO_INTERP - [boolean] flag to turn off interpolation of likelihoods, 
;              to be used when fitting with classes of templates that
;              are not linked by physical parameters
;  FILTERS - [string array] string names of filters used to determine
;            fluxes at observed wavelengths, default is to perform 
;            interpolation on the template SEDs. Should contain the same 
;            number of elements as observed wavelengths in the same order
;            given FNAME.  If no filter is desired for a particular
;            wavelength but still want to include for others, specify
;            by leaving the corresponding entry blank, i.e. [''].
;  FLUX_GRID - [double array] Array of previously interpolated flux values,
;              used when pre-computing the likelihood grid when performing
;              photo-z estimates.  Should contain the same number of elements
;              observed wavelenghts in the same order given in FNAME.  If
;              FLUX_GRID is initially empty, it will return the 
;              interpolated flux grid used for the current fit.
;  CAL_LGRID - [boolean] if set, the full likelihood grid is computed
;              before the MCMC process.  'Best guess' initial points
;              for the chains will be determined from the likelihood
;              grid.  This should only be used in cases of determining
;              photo-z's where you believe redshift and other
;              parameters to be highly degenerate.
;  PRIORS - [string array] list of priors to consider during the fitting.  
;           Priors should be defined by their own IDL function and called
;           here using their basic syntax. See documentation for more
;           details.
;  PRIOR_MODE - [integar array] Array defining the type of prior(s) used
;               Accepted values are:
;                  1) Parameter coupling
;                  2) Additional likelihood modifiers
;               See documentation for more details.
;  SAVE_ALL - [boolean] flag for determining if all steps are to be saved
;             or just the final sampling.  Default is to have this option
;             turned off to conserve disk space.
;


   !PATH = !PATH + ':~/IDL/SATMC/'   

if n_params() le 1 then begin
    print,'CALLING SEQUENCE - sed_fit_mcmc,fname,sed1,sed2,sed3,'+$
      'norm1=norm1,norm2=norm2,norm3=norm3,redshift=redshift,'+$
      'findz=findz,outfile=outfile,psfile=psfile,no_plot=no_plot,'+$
      'xrange=xrange,ig_lim=ig_lim,axes=axes'+$
      'nstep=nstep,step_int=step_int,nchain=nchain,'+$
      'params_init=params_init,params_std=params_std,delz=delz,'+$
      'param_max=param_max,param_min=param_min,minz=minz,maxz=maxz,'+$
      'adaptive_step=adaptive_step,'+$
      'synthesis_model=synthesis_model,no_async=no_async,restart=restart,'+$
      'no_temp=no_temp,no_interp=no_interp,filters=filters,'+ $
      'use_xvfb=use_xvfb,priors=priors,prior_mode=prior_mode,save_all=save_all'
    return
endif

if n_elements(redshift) eq 0 then begin
   if not keyword_set(findz) then begin
      message,/cont,'Warning: REDSHIFT nor FINDZ is set.  Assuming FINDZ=True.'
      findz=1
   endif
endif else redshift=float(redshift)
if n_elements(nchain) eq 0 then nchain=5
if nchain eq 1 then begin
   ;turn off multiple thread techniques for a single chain MCMC
   no_temp=1
   no_async=1
endif
if n_elements(nstep) eq 0 then nstep=10000L else nstep=long(nstep)
if n_elements(step_int) eq 0 then step_int=500L else step_int=long(step_int)
nsed=n_params()-1
if not keyword_set(restart) then restart=0
if n_elements(axes) eq 0 then axes=1
if keyword_set(priors) and not keyword_set(prior_mode) then begin
   message,/cont,'Type of priors not specified, assuming modifier of'+$
           ' likelihoods'
   prior_mode=intarr(n_elements(priors))+2
endif

;place SED models into a single structure (helps simplify matters when
;combining SEDs)
sed_arr='sed'+strtrim(indgen(nsed)+1,2)
for i=0,nsed-1 do begin
   if i eq 0 then begin
      ;check for LUM tag, calculate and add if missing
      ;ignore if synthesized models are to be used
      if not keyword_set(synthesis_model) then begin
         if tag_exist(scope_varfetch(sed_arr[0]),'lum') eq 0 then $
            cal_lum,scope_varfetch(sed_arr[0])
      endif
      sed_struct=create_struct(scope_varname(scope_varfetch(sed_arr[0]),level=-1),scope_varfetch(sed_arr[0]))
   endif else begin
      if tag_exist(scope_varfetch(sed_arr[i]),'lum') eq 0 then $
         cal_lum,scope_varfetch(sed_arr[i])
      sed_struct=create_struct(sed_struct,scope_varname(scope_varfetch(sed_arr[i]),level=-1),scope_varfetch(sed_arr[i]))
   endelse
endfor

;setup useful constants
c=2.9979e14                     ;(microns/s)

;read in format
openr,lun,fname,/get
text=''
readf,lun,text
free_lun,lun
xf=gettok(text,' ')
yf=gettok(text,' ')
zf=gettok(text,' ')
ulf=gettok(text,' ')
readcol,fname,x,y,z,ul,/silent,skipline=1,format='D,D,D,A',delim=' '
case strlowcase(strmid(xf,0,1)) of
   'l': wavel=(c/10^x) 
   'm': wavel=x
   'a': wavel=x*1.e-4
   else: wavel=(c/10^[xf,x])
endcase
case strlowcase(strmid(yf,0,1)) of
   'm': flux=y*1.e-3 
   'j': flux=y
   'u': flux=y*1.e-6
   else: flux=[yf,y]*1.e-3
endcase

case strlowcase(strmid(zf,0,1)) of
   'm': fluxerr=z*1.e-3 
   'j': fluxerr=z
   'u': fluxerr=z*1.e-6
   else: fluxerr=[zf,z]*1.e-3
endcase

if not valid_num(strmid(ulf,0,1)) then begin
   upperlim=ul
endif else upperlim=[ulf,ul]
;endcase

;check that filters is compatible if set
if keyword_set(filters) then begin
   if n_elements(filters) ne n_elements(wavel) then begin
      message,'FILTERS must have same number of elements as input SED.'
      return
   endif
endif

;remove 0 flux value elements
sel=where(flux ne 0)
flux=flux[sel]
fluxerr=fluxerr[sel]
wavel=wavel[sel]
upperlim=upperlim[sel]
if keyword_set(filters) then filter=filters[sel]
if keyword_set(flux_grid) then flux_grid=flux_grid[sel,*]

;sort by wavelength
sortindex = reverse(sort(wavel))
wavel=wavel[sortindex]
flux=flux[sortindex]
fluxerr=fluxerr[sortindex]
upperlim=upperlim[sortindex]
if keyword_set(filters) then filter=filter[sortindex]
if keyword_set(flux_grid) then flux_grid=flux_grid[sortindex,*]

;check if there are any upper limits
good=where(upperlim eq 0,n_good,complement=bad,ncomplement=n_bad)
if n_bad ne 0 then begin 
   wavelim=wavel(bad)                   ;wavelength at upper limits
   fluxlim=(flux(bad)) > (fluxerr(bad)) ;flux upper limit
endif 

;check if there is enough points for reasonable fit
message,/cont,'Number of useful data points: '+string(n_good)
if n_bad ne 0 then message,/cont,'Number of useful upper limits: '+string(n_bad)
if n_good le 2 then begin
   message,/cont,'Not enough data points, quitting.'
   return
endif

;isolate only the good measurements
goodwavel=wavel(good)
goodflux=flux(good)
goodfluxerr=fluxerr(good)

if not keyword_set(no_plot) then begin
    if keyword_set(psfile) then begin
       set_plot,'ps'
       device,filename=psfile+'.ps',bits=8,color=1,/encap
       color=0
    endif else begin
       device,decomposed=0
       window,0,xsize=500,ysize=500,retain=2
       color=1
    endelse
    usersym,[-1,1,0,0,-1,0,1],[0,0,0,-2,-1,-2,-1]
    newerrpos=goodflux+goodfluxerr ; error bars
    newerrneg=goodflux-goodfluxerr
    sname=gettok(fname,'.dat')
    while strpos(sname,'/') ne -1 do null=gettok(sname,'/')
    if not keyword_set(xrange) then begin
       minwave=10.^(floor(alog10(goodwavel)))
       maxwave=10.^(ceil(alog10(goodwavel)))
       xrange=[min(minwave)*0.1,max(maxwave)*10]
    endif
    case axes of
       1: begin
          plot,goodwavel,goodflux,tit=sname, $
               xtit='Observed Wavelength (micron)',ytit='Observed Flux (Jy)', $
               /ylog,/xlog,psym=6,xrange=xrange, $
               yrange=[min(goodflux)*0.1,max(goodflux)*10.],subtit='' 
          errplot,goodwavel,newerrneg,newerrpos ;plot flux error
          if bad[0] ne -1 then oplot,wavelim,fluxlim,psym=8,symsize=2
       end
       2: begin
          lfl=goodflux*3.e-9/goodwavel
          errneg=newerrneg*3.e-9/goodwavel
          errpos=newerrpos*3.e-9/goodwavel
          plot,goodwavel,lfl,tit=sname,xtit='Observed Wavelength (micron)', $
               ytit=textoidl('Observed \lambda F_{\lambda}'),/ylog,/xlog, $
               psym=6,xrange=xrange, $
               yrange=[min(goodflux)*0.1,max(goodflux)*10.],subtit='' 
          errplot,goodwavel,errneg,errpos ;plot flux error
          if bad[0] ne -1 then begin
             lfl=fluxlim*3.e-9/wavelim
             oplot,wavelim,lfl,psym=8,symsize=2
          endif
       end
       else: begin
          message,/cont,'Current plot mode not specified.  Quitting.'
          return
       end
    endcase
endif

; **********  MCMC  **********************************
;extract sed parameter information and setup parameter arrays
nparams=0
param_name=['']
p_max=[-999.]
p_min=[-999.]
npgrid=[-999]

nparm_arr=intarr(nsed)          ;int array to hold number of params in each SED
norm_arr=intarr(nsed)           ;boolean array for model normalization
for i=0,nsed-1 do begin
    if keyword_set(scope_varfetch('norm'+strtrim(string(i+1),2),level=0)) $
      then norm_arr[i]=1
    ;make sure NORM is false if using synthesized models
    if i eq 0 and keyword_set(synthesis_model) then norm_arr[i]=0
endfor

;setup filters
if keyword_set(filters) then begin
   maxpt=max(file_lines(filter))
   filter_struct=create_struct('wave',fltarr(maxpt),'response',fltarr(maxpt),$
                              'norm',0.0d,'flag',0)
   filter_struct=replicate(filter_struct,n_elements(filter))

   filter_struct[good].flag=1
   for i=0,n_elements(filter)-1 do begin
      if keyword_set(filter[i]) then begin
         rdfloat,filter[i],fwave,fresponse,/silent
         h=6.626e-27
         filter_struct[i].wave=fwave
         filter_struct[i].response=fresponse
         filter_struct[i].norm=integral(c/fwave,fresponse/(h*c/fwave))
      endif
   endfor
endif else filter_struct=[0]

for i=0,nsed-1 do begin
   get_sed_param,sed_struct.(i),sedgrid,sed_param_name,nparm, $
                 synthesis_model=synthesis_model,norm=norm_arr[i], $
                 ngrid=ngrid,no_interp=no_interp                 
   if i eq 0 then grid_struct=create_struct(sed_arr[i],sedgrid) else $
      grid_struct=create_struct(grid_struct,sed_arr[i],sedgrid)
   param_name=[param_name,replicate(scope_varname(scope_varfetch(sed_arr[i]), $
                                                  level=-1),nparm)+' '+ $
               sed_param_name]
   nparm_arr[i]=nparm
   if norm_arr[i] then begin
      norm=scope_varfetch('norm'+strtrim(string(i+1),2),level=0)
      if n_elements(norm) eq 1 then begin
         max=norm[0]
         min=[1.0]
      endif else begin
         max=norm[1]
         min=norm[0]
      endelse
      if nparm gt 1 then begin
         max=[max,max(sedgrid,dim=1)]
         min=[min,min(sedgrid,dim=1)]
         ngrid=[2,ngrid]
      endif else ngrid=[2]
   endif else begin
      ;check that the fit is meaningful, i.e. at least one parameter
      ;to search over (including normalization)
      if nparm eq 0 then begin
         message,'Warning: No parameters in SED'+strtrim(string(i+1),2)+ $
                 ' and keyword Norm not set. Quitting.'
      endif else max=max(sedgrid,dim=1,min=min)
   endelse
   p_max=[p_max,max]
   p_min=[p_min,min]
   npgrid=[npgrid,ngrid]
endfor
nparams=total(nparm_arr)
param_name=param_name[1:*]
npgrid=npgrid[1:*]
if n_elements(param_max) eq 0 then param_max=p_max[1:*] else begin
   param_max=min([[param_max],[p_max[1:*]]],dim=2)
endelse
if n_elements(param_min) eq 0 then param_min=p_min[1:*] else begin
   param_min=max([[param_min],[p_min[1:*]]],dim=2)
endelse

;check for large dynamical ranges, likely to cause bad fits
dr=where((abs(alog10(param_max)-alog10(param_min)) ge 2 and param_min ne 0) or $
         (alog10(param_max) ge 2 and param_min eq 0),ndr)
if ndr ne 0 then begin
   param_max[dr]=alog10(param_max[dr])
   param_min[dr]=alog10(param_min[dr])
   sel=where(finite(param_min[dr],/inf),nsel)
   if nsel ne 0 then param_min[dr[sel]]=0.

   p0=0
   for l=0,nsed-1 do begin
      if norm_arr[l] then begin
         if nparm_arr[l] gt 1 then begin
            gs=where(dr ge p0+1 and dr le p0+nparm_arr[l]-1,nsel)
            if nsel gt 0 then grid_struct.(l)[*,dr[gs]-1]=alog10(grid_struct.(l)[*,dr[gs]-1])
         endif
      endif else begin
         gs=where(dr ge p0 and dr le p0+nparm_arr[l]-1,nsel)
         if nsel gt 0 then grid_struct.(l)[*,dr[gs]]=alog10(grid_struct.(l)[*,dr[gs]])
      endelse
      p0+=nparm_arr[l]
   endfor
endif
err=check_math()                ;ignore math errors

;setup mcmc arrays
template=create_struct('lnl_iter',dblarr(nstep,nchain),'lum_iter', $
                       dblarr(nstep,nchain),'params_iter',$
                       fltarr(nparams+1,nstep,nchain),'params_std', $
                       fltarr(nparams+1,nparams+1,nchain))
shmmap,template=template,destroy_segment=1,get_name=chainmmap
chain_struct=shmvar(chainmmap)

if n_elements(params_init) eq 0 then begin
   ;set initial guesses
   for i=0,nchain-1 do begin
      chain_struct.params_iter[0:nparams-1,0,i]= $
         satmc_random(n_elements(param_max),/uniform)* $
         (param_max-param_min)+param_min
   endfor
endif else begin
   ;if user supplied initial parameter guess, make sure it falls
   ;within parameter space
   if ndr ne 0 then params_init[dr]=alog10(params_init[dr])
   sel=where(params_init lt param_min or params_init gt param_max,nout)
   if nout ne 0 then begin
      message,/cont,'Initial parameter guess outside parameter bounds.'
      shmunmap,chainmmap
      return
   endif else begin
      chain_struct.params_iter[0:nparams-1,0,*]= $
         rebin(params_init,n_elements(params_init),nchain)
   endelse
endelse 

;determine whether the MCMC runs over redshift
z=fltarr(1,nstep,nchain)
if keyword_set(findz) then begin
   if not keyword_set(maxz) then maxz=6.0d
   if not keyword_set(minz) then minz=1.d-3
   if not keyword_set(delz) then delz=0.01
   if keyword_set(redshift) then z[0,0,*]=redshift else $
      z[0,0,*]=satmc_random(nchain,/uniform)*(maxz-minz)+minz
endif else begin
   z[0,*,*]=redshift
   maxz=redshift
   minz=redshift
   delz=0.0
endelse
param_name=[param_name,'z']
chain_struct.params_iter[nparams,*,*]=z
param_max=[param_max,maxz]
param_min=[param_min,minz]

if n_elements(params_std) eq 0 then begin
   params_std=(0.01*(param_max-param_min)/(nparams+1))^2.
   if ndr ne 0 then begin
      params_std[dr]=(0.01*(param_max[dr]-param_min[dr])/(nparams+1))^2.
   endif
   for i=0,n_elements(npgrid)-1 do begin
      if npgrid[i] gt 2 then $
         params_std[i]=((param_max[i]-param_min[i])/npgrid[i])^2.
   endfor
   params_std[nparams]=delz^2.
   chain_struct.params_std[*]=rebin(diag_matrix(params_std),nparams+1,nparams+1,nchain)
   params_std=chain_struct.params_std[*,*,0]
endif else begin
   if (size(params_std))[1] ne nparams+1 then begin
      message,/cont,'Input covariance matrix must have dimensions '+$
              'equal to the number of available parameters.'
      return
   endif
   chain_struct.params_std[*]=rebin(params_std,nparams+1,nparams+1,nchain)
endelse

;if doing photo-z's, precompute likelihood grid for grid-based 
;templates and pick initial points/covariance matrix from likelihood grid
;avoids problems of local maxima when redshift space is widely distributed
if keyword_set(cal_lgrid) and keyword_set(findz) then begin
   cal_lgrid,sed_struct,grid_struct,grid,lgrid,wavel,flux,fluxerr,upperlim, $
             norm_arr,nparm_arr,param_max,param_min,param_name=param_name, $
             minz=minz,maxz=maxz,dr=dr,flux_grid=flux_grid,priors=priors,$
             prior_mode=prior_mode
   if not keyword_set(redshift) then begin
      sel=(reverse(sort(lgrid)))[0:nchain-1]
      chain_struct.params_iter[*,0,*]=transpose(grid[sel,*])
    endif
   sel=(reverse(sort(lgrid)))[0:n_elements(lgrid)*0.001 > 99]
   cov=correlate(transpose(grid[sel,*]),/cov)*1.e-5
   chain_struct.params_std[*]=rebin(cov,nparams+1,nparams+1,nchain)
endif

;do the mcmc chain
thres=float(2)/step_int         ;threshold value for acceptance rate
maxchain=100                    ;maximum number of MCMC chains to attempt
nretry=0
pchange=0
rold=0
pairold=0
lastswap=0
t_period=0

if keyword_set(save_all) then begin
   params_all=chain_struct.params_iter[*,0,*]*0.d
   lnl_all=chain_struct.lnl_iter[0,*]*0.d
endif

if not keyword_set(no_async) then begin
   ;create virtual X buffer if requested
   ;allows the code run in the background and continue without having
   ;to stay logged into an X-session
   if keyword_set(use_xvfb) then begin
      ;first check for an available display and save current settings
      current_display=getenv('DISPLAY')
      display=1
      spawn,'ps -aef|grep Xvfb',ps_xvfb
      ps_sel=where(strmatch(ps_xvfb,'*Xvfb :*') eq 1,nsel)
      if nsel ne 0 then begin
         display_in_use=[-999]
         for i=0,nsel-1 do begin
            xvfb_dis=strmid(ps_xvfb[ps_sel[i]], $
                            strpos(ps_xvfb[ps_sel[i]],'Xvfb')+6,1)
            display_in_use=[display_in_use,xvfb_dis]
         endfor
         sel=where(display eq display_in_use,nsel)
         while nsel ne 0 do begin
            display++
            sel=where(display eq display_in_use,nsel)
         endwhile
      endif
      spawn,['Xvfb',':'+strtrim(display,2)+' -screen 0 1024x748x24'], $
            pid=xvfb_pid,/noshell,unit=unit
      setenv,'DISPLAY=:'+strtrim(display,2)+'.0'
   endif

   ;create a single IDL child process that will control the 'chain children' 
   ;necessary to avoid leaking memory when performing asynchronous processes
   parent=obj_new('IDL_IDLBridge',output='parent.log')
      
   if keyword_set(lgrid) then begin
      save,sed_struct,goodwavel,goodflux,goodfluxerr,param_name, $
           grid_struct,norm_arr,nparm_arr,param_max,param_min, $
           filter_struct,dr,filename='vars.sav'
   endif else begin
      save,sed_struct,goodwavel,goodflux,goodfluxerr,param_name, $
           grid_struct,norm_arr,nparm_arr,param_max,param_min, $
           filter_struct,dr,filename='vars.sav'
   endelse

   cd,current=current           ;necessary so children look in correct dir
   savefile=current+'/vars.sav'
   
   if keyword_set(synthesis_model) then  $
     parent->SetVar,'synthesis_model',synthesis_model
   if keyword_set(no_temp) then parent->SetVar,'no_temp',no_temp
   if keyword_set(no_interp) then parent->SetVar,'no_interp',no_interp
   if keyword_set(ig_lim) then parent->SetVar,'ig_lim',ig_lim
   if keyword_set(priors) then parent->SetVar,'priors',priors
   if keyword_set(prior_mode) then parent->SetVar,'prior_mode',prior_mode
   if n_elements(wavelim) ne 0 then begin
      parent->SetVar,'wavelim',wavelim
      parent->SetVar,'fluxlim',fluxlim
      parent->SetVar,'upperlim',upperlim
   endif
        
   ;establish links between parent and main routine variables
   parent->Execute,"template=create_struct('lnl_iter',dblarr("+ $
                   strtrim(nstep,2)+","+strtrim(nchain,2)+ $
                   "),'lum_iter',dblarr("+strtrim(nstep,2)+","+ $
                   strtrim(nchain,2)+"),'params_iter',fltarr("+ $
                   strtrim(nparams+1,2)+","+strtrim(nstep,2)+","+ $
                   strtrim(nchain,2)+"),'params_std',fltarr("+ $
                   strtrim(nparams+1,2)+","+strtrim(nparams+1,2)+","+ $
                   strtrim(nchain,2)+"))"
   parent->Execute,"shmmap,'"+chainmmap+"',template=template"
   parent->Execute,"chain_struct=shmvar('"+chainmmap+"')"
   
   ;now have parent IDL process create a child for each chain
   parent->SetVar,'savefile',savefile
   parent->Execute,'@'+PREF_GET('IDL_STARTUP')
   parent->Execute,'@~/IDL/SATMC/start_satmc.pro'
   parent->execute,"resolve_routine,'satmc',/compile_full"
   parent->Execute,'child=create_children('+strtrim(nchain)+','+$
                   'savefile,synthesis_model=synthesis_model,'+$
                   'wavelim=wavelim,fluxlim=fluxlim,upperlim=upperlim,'+$
                   'dir="'+current+'/",no_temp=no_temp,no_interp=no_interp,'+$
                   'ig_lim=ig_lim,priors=priors,prior_mode=prior_mode)'

   ;remove setup savefiles
   spawn,'rm '+savefile
   
   ;have parent start the MCMC chains, first to test for convergence
   ;and acceptance rate
   if not keyword_set(adaptive_step) then begin
      if restart eq 1 then begin
         restore,'parent_burnin.sav'
         t_chain_struct=chain_struct
         ;verify all tags exist between saved run and current
         nstags=n_tags(t_chain_struct)
         stname=tag_names(t_chain_struct)
         chain_struct=shmvar(chainmmap)
         ntags=n_tags(chain_struct)
         tname=tag_names(chain_struct)
         if ntags ne nstags then begin
            ;tags are missing, will have to destory chainmmap and remake
            template=t_chain_struct
            parent->execute,"template=chain_struct"
            for i=0,nstags-1 do begin
               if where(stname[i] eq tname) eq -1 then begin
                  parent->setvar,'mtname',stname[i]
                  parent->setvar,'mtval',t_chain_struct[0].(i)
                  parent->execute,"template=create_struct(template,"+$
                                  "mtname,mtval)"
               endif
            endfor
            parent->Execute,"shmunmap,'"+chainmmap+"'"
            shmunmap,chainmmap
            shmmap,template=template,destroy_segment=1,get_name=chainmmap
            chain_struct=shmvar(chainmmap)
            parent->Execute,"shmmap,'"+chainmmap+"',template=template"
            parent->Execute,"chain_struct=shmvar('"+chainmmap+"')"
         endif
         chain_struct[*]=temporary(t_chain_struct)
         pchange=0
      endif
      
      if restart le 1 then begin
         while (pchange le 0.20 or pchange ge 0.26) do begin
            for m=0,nchain-1 do begin
               parent->Execute,"child["+strtrim(m,2)+"]->SetVar,"+$
                               "'param_init',chain_struct.params_iter[*,0,"+ $
                               strtrim(m,2)+"]"
               parent->Execute,"child["+strtrim(m,2)+"]->SetVar,"+$
                                "'param_std',chain_struct.params_std[*,*,"+ $
                               strtrim(m,2)+"]"
               parent->Execute,"outfile='"+current+"/child"+ $
                               strtrim(m,2)+".sav'"
               parent->Execute,"child["+strtrim(m,2)+"]->"+$
                               "Setvar,'outfile',outfile"
               parent->Execute,"child["+strtrim(m,2)+"]->"+$
                               "Setvar,'lnl_init',chain_struct.lnl_iter[0,"+ $
                               strtrim(m,2)+"]"
            endfor
               
            ;if doing parallel tempering, make potential swaps between
            ;the fudical cold and warm chains
            if not keyword_set(no_temp) and t_period eq 2 ne 0 then begin
               t_period=0
               maxl=max(chain_struct.lnl_iter[0:2*step_int-1,*],dim=1,psel)
               msel=(where(maxl eq max(maxl)))[0]
               psel=(array_indices([2*step_int,nchain],psel,/dim))[0,*]

               if msel[0] ne 0 then begin
                  parent->Execute,"child["+strtrim(msel[0],2)+"]->SetVar,"+$
                                  "'param_init',chain_struct.params_iter[*,"+ $
                                  strtrim(psel[msel],2)+",0]"
                  parent->Execute,"child[0]->SetVar,"+$
                                  "'param_init',chain_struct.params_iter[*,"+ $
                                  strtrim(psel[msel],2)+","+strtrim(msel,2)+"]"
                  parent->Execute,"child["+strtrim(msel,2)+"]->"+$
                                  "Setvar,'lnl_init',chain_struct.lnl_iter["+ $
                                  strtrim(psel[0],2)+",0]"
                  parent->Execute,"child[0]->Setvar,'lnl_init',"+$
                                  "chain_struct.lnl_iter["+ $
                                  strtrim(psel[msel],2)+","+strtrim(msel,2)+"]"
                  parent->Execute,"child[0]->SetVar,'param_std',"+$
                                  "chain_struct.params_std[*,*,0]*0.001"
               endif
            endif else t_period++
            
            for m=0,nchain-1 do begin
               parent->Execute,"child["+strtrim(m,2)+"]->Execute,"+$
                 "'run_mcmc_chain,"+strtrim(step_int*2,2)+",param_init,"+$
                 "lnl_init,sed_struct,grid_struct,goodwavel,goodflux,"+$
                 "goodfluxerr,wavelim,fluxlim,upperlim,param_name,norm_arr,"+$
                 "nparm_arr,param_min,param_max,param_std,param_arr,"+$
                 "lnl_arr,lum_arr,sed,models"+strtrim(m,2)+$
                 ",synthesis_model=synthesis_model,outfile=outfile,"+$
                 "temp=temp,no_interp=no_interp,filter_struct=filter_struct,"+$
                 "dr=dr,ig_lim=ig_lim,priors=priors,prior_mode=prior_mode',"+$
                 "/nowait"
            endfor
            ;wait for chains to finish...
            for m=0,nchain-1 do begin
               while parent-> $
                  GetVar('child['+strtrim(m,2)+']->status()') do wait,0.01
            endfor
           
            ;...and continue
            ;combine child results
            for m=0,nchain-1 do begin
               if keyword_set(synthesis_model) then begin
                  if not tag_exist(chain_struct,'SED') then begin
                     ;in order to add the SED information for each
                     ;parameter step, we'll need to destory and
                     ;recreate CHAIN_STRUCT
                     old_chain_struct=chain_struct
                     parent->Execute,"old_chain_struct=chain_struct"
                     parent->Execute,"sed=child[0]->GetVar('sed')"
                     i=0
                     while n_elements(size(parent->GetVar('sed'),/dim)) eq 1 $
                     and i lt nchain-1 do begin
                        i++
                        parent->Execute,"sed=child["+strtrim(i,2)+ $
                                        "]->GetVar('sed')"
                     endwhile
                     if n_elements(size(parent->GetVar('sed'),/dim)) ne 1 then begin
                        parent->Execute,"template_sed=rebin(sed,["+ $
                                        "(size(sed,/dim))[0:1],"+$
                                        strtrim(nstep,2)+","+$
                                        strtrim(nchain,2)+"])"
                        parent->Execute,"template=create_struct("+$
                                        "old_chain_struct,'SED',template_sed)"
                        template_sed=parent->GetVar('template_sed')
                        template=create_struct(old_chain_struct,'SED', $
                                               template_sed)
                        parent->Execute,"shmunmap,'"+chainmmap+"'"
                        shmunmap,chainmmap
                        
                        shmmap,template=template,destroy_segment=1, $
                               get_name=chainmmap
                        chain_struct=shmvar(chainmmap)
                        parent->Execute,"shmmap,'"+chainmmap+"',template=template"
                        parent->Execute,"chain_struct=shmvar('"+chainmmap+"')"
                        chain_struct[*]=template
                     endif
                  endif 
                  parent->Execute,"sed=child["+strtrim(m,2)+"]->GetVar('sed')"
                  if n_elements(size(parent->GetVar('sed'),/dim)) ne 1 then begin
                     parent->Execute,"chain_struct.SED[*,*,0:"+ $
                                     strtrim(step_int*2-1,2)+","+$
                                     strtrim(m,2)+"]=sed"
                  endif
               endif
               parent->Execute,"chain_struct.params_iter[*,0:"+ $
                               strtrim(step_int*2-1,2)+","+ $
                               strtrim(m,2)+"]=child["+strtrim(m,2)+$
                               "]->GetVar('param_arr')"
               parent->Execute,"chain_struct.lnL_iter[0:"+ $
                               strtrim(step_int*2-1,2)+","+ $
                               strtrim(m,2)+"]=child["+strtrim(m,2)+$
                               "]->GetVar('lnl_arr')"
               parent->Execute,"chain_struct.lum_iter[0:"+ $
                               strtrim(step_int*2-1,2)+","+ $
                               strtrim(m,2)+"]=child["+strtrim(m,2)+$
                               "]->GetVar('lum_arr')"
            endfor

            if keyword_set(save_all) then begin
               params_all=[[params_all], $
                           [chain_struct.params_iter[*,0:step_int*2-1,*]]]
               lnl_all=[lnl_all,chain_struct.lnl_iter[0:step_int*2-1,*]]
            endif

            ;check acceptance rate of mcmc steps
            params_std_old=chain_struct.params_std
            pchange_ch=dblarr(nchain)
            
            for m=0,nchain-1 do begin
               dlnL_iter=chain_struct.lnL_iter[step_int:2*step_int-1,m]- $
                         chain_struct.lnL_iter[step_int+1:2*step_int,m]
               whchange=where(abs(dlnL_iter[0:step_int-2]) gt 0.0,nchange)
               pchange_ch[m]=float(nchange)/float(step_int-1)
            endfor

            if keyword_set(no_temp) then begin
               pchange=avg(pchange_ch)
               sel=where(pchange_ch eq 0,nsel)
               if nsel ne 0 then pchange=0.
               if pchange le 0.20 or pchange ge 0.26 then begin
                  for m=0,nchain-1 do begin
                     cov=chain_struct.params_std[*,*,m]
                     if pchange_ch[m] gt 0.26 then begin
                        diag_pstd=diag_matrix(chain_struct.params_std[*,*,m])
                        pstd=chain_struct.params_std[*,*,m]
                        if array_equal(diag_matrix(diag_pstd),pstd) then begin
                           cov=correlate(chain_struct.params_iter[*,step_int:$
                                                                  2*step_int-1,m],/cov,/double)
                        endif else begin
                           frac=pchange_ch[m]/0.23
                           cov=chain_struct.params_std[*,*,m]*frac
                        endelse
                     endif
                     if pchange_ch[m] lt 0.20 and pchange_ch[m] gt 0.05 then begin
                        cov=correlate(chain_struct.params_iter[*,step_int: $
                                                               2*step_int-1, $
                                                               m],/cov,/double)
                     endif
                     if pchange_ch[m] lt 0.05 then begin
                        cov=chain_struct.params_std[*,*,m]*0.05
                     endif
                     if array_equal(cov,0) then cov=params_std
                     diag=diag_matrix(cov)
                     diag_pstd=diag_matrix(params_std)
                     sel=where(sqrt(diag) lt sqrt(diag_pstd)*1.e-6 or $
                               sqrt(diag) gt sqrt(diag_pstd)*1.e4,nsel)
                     if nsel ne 0 then cov=params_std
                     chain_struct.params_std[*,*,m]=cov
                  endfor                    
               endif
            endif else begin
               pchange=pchange_ch[0]
               sel=where(pchange_ch eq 0,nsel)
               if nsel ne 0 then pchange=0.
               if pchange le 0.20 or pchange ge 0.26 then begin
                  for m=0,nchain-1 do begin
                     cov=chain_struct.params_std[*,*,m]
                     if pchange_ch[m] le 0.20 and $
                        pchange_ch[m] gt 0.05 then begin
                        cov=correlate(chain_struct.params_iter[*,step_int: $
                                                               2*step_int-1,$
                                                               m],/cov,/double)
                     endif 
                     if pchange_ch[m] gt 0.26 then begin
                        diag_pstd=diag_matrix(chain_struct.params_std[*,*,m])
                        pstd=chain_struct.params_std[*,*,m]
                        if array_equal(diag_matrix(diag_pstd),pstd) then begin
                           cov=correlate(chain_struct.params_iter[*,step_int:$
                                                                  2*step_int-1,m],/cov,/double)
                        endif else begin
                           frac=pchange_ch[m]/0.23
                           cov=chain_struct.params_std[*,*,m]*frac
                        endelse
                     endif
                     if pchange_ch[m] lt 0.05 then begin
                        cov=chain_struct.params_std[*,*,m]*0.05
                     endif
                     if array_equal(cov,0) then cov=params_std
                     diag=diag_matrix(cov)
                     diag_pstd=diag_matrix(params_std)
                     sel=where(sqrt(diag) lt sqrt(diag_pstd)*1.e-6 or $
                               sqrt(diag) gt sqrt(diag_pstd)*1.e4,nsel)
                     if nsel ne 0 then cov=params_std
                     chain_struct.params_std[*,*,m]=cov
                  endfor
               endif
            endelse

            ;determine current maximum likelihood value
            maxL_ch=dblarr(nchain)
            for m=0,nchain-1 do $
               maxL_ch[m]=max(chain_struct.lnl_iter[where(chain_struct.lnl_iter[*,m] ne 0.0),m])
            diffl=(maxL_ch-max(maxL_ch))
            
            if pchange ge 0.20 and pchange le 0.26 then begin
               ;Geweke dianogstic
                if not keyword_set(no_temp) then begin
                  av1=avg(chain_struct.lnL_iter[0:0.2*step_int-1,0])
                  av2=avg(chain_struct.lnL_iter[step_int:2*step_int-1,0])
                  s=stddev(chain_struct.lnL_iter[step_int:2*step_int-1,0])
                  if s ne 0 then begin
                     if abs(av1-av2)/s gt 2.0 then pchange=0.
                  endif
               endif else begin
                  for m=0,nchain-1 do begin
                     av1=avg(chain_struct.lnL_iter[0:0.2*step_int-1,m])
                     av2=avg(chain_struct.lnL_iter[step_int:2*step_int-1,m])
                     s=stddev(chain_struct.lnL_iter[step_int:2*step_int-1,m])
                     if s ne 0 then begin
                        if abs(av1-av2)/s gt 2.0 then pchange=0.
                     endif
                  endfor
               endelse
            endif

            ;check for bad chains, chains that either are unable
            ;move about (flat likelihood space) or too far from
            ;current maximum likelihood value
            if pchange ge 0.20 and pchange le 0.26 then begin
               badchain=where(pchange_ch lt thres or $
                              abs(diffl) gt 1.e2,nbc,com=goodchain)
            endif else begin
               badchain=where(pchange_ch lt thres and diffl ne 0, $
                              nbc,com=goodchain)
            endelse
            if nbc ne 0 then begin
               if nbc eq nchain then begin
                  for i=0,nbc-1 do begin
                     chain_struct.params_iter[*,0,badchain[i]]= $
                        satmc_random(n_elements(param_max),/uniform)* $
                            (param_max-param_min)+param_min

                     chain_struct.lnl_iter[0,badchain[i]]=0
                     chain_struct.params_std[*,*,badchain[i]]= $
                        params_std_old[*,*,badchain[i]]
                  endfor
                  nretry++
               endif else begin
                  ;use the best estimate from other chains to determine
                  ;new initial step
                  params_good=chain_struct.params_iter[*,step_int: $
                                                       2*step_int-1, $
                                                       goodchain]
                  if n_elements(goodchain) eq 1 then begin
                     for i=0,nbc-1 do begin
                        rind=round(satmc_random(/uniform)*(step_int-1))
                        chain_struct.params_iter[*,0,badchain[i]]= $
                           params_good[*,rind]
                        chain_struct.lnl_iter[0,badchain[i]]= $
                           chain_struct.lnl_iter[rind,goodchain]
                        chain_struct.params_std[*,*,badchain[i]]= $
                           params_std_old[*,*,goodchain]
                     endfor
                  endif else begin
                     for i=0,nbc-1 do begin
                        rind=round(satmc_random(/uniform)*(step_int-1))
                        if pchange ge 0.20 and pchange le 0.26 then begin
                           rchain=0
                        endif else begin
                           ng=n_elements(goodchain)
                           rchain=round(satmc_random(/uniform)*(ng-1))
                        endelse
                        chain_struct.params_iter[*,0,badchain[i]]= $
                           params_good[*,rind,rchain]
                        chain_struct.lnl_iter[0,badchain[i]]= $
                           chain_struct.lnl_iter[rind,goodchain[rchain]]
                        chain_struct.params_std[*,*,badchain[i]]= $
                           params_std_old[*,*,goodchain[rchain]]
                     endfor
                  endelse
                  chain_struct.params_iter[*,0,goodchain]= $
                     chain_struct.params_iter[*,2*step_int-1,goodchain]
                  chain_struct.lnl_iter[0,goodchain]= $
                     chain_struct.lnl_iter[2*step_int-1,goodchain]
               endelse
               pchange=0.
            endif else begin
               chain_struct.params_iter[*,0,*]= $
                  chain_struct.params_iter[*,2*step_int-1,*]
               chain_struct.lnl_iter[0,*]= $
                  chain_struct.lnl_iter[2*step_int-1,*]
            endelse

            ;create/update savefile
            save,chain_struct,pchange,filename='parent_burnin.sav'
         endwhile
         params_std=chain_struct.params_std[*,*,0]
      endif
      period=2
   endif else period=0
   ;Now do the full chains, seperate into blocks of STEP_INT length
   ;to set a restart point in case of chain failure.
   ;Save the combined results as well as individual chains.
   if restart le 2 then begin
      swap=0
      if restart eq 2 then begin
         restore,'parent_chain.sav'
         t_chain_struct=chain_struct
         chain_struct=shmvar(chainmmap)
         ;verify all tags exist between saved run and current
         nstags=n_tags(t_chain_struct)
         stname=tag_names(t_chain_struct)
         chain_struct=shmvar(chainmmap)
         ntags=n_tags(chain_struct)
         tname=tag_names(chain_struct)
         if ntags ne nstags then begin
            ;tags are missing, will have to destory chainmmap and remake
            template=t_chain_struct
            parent->execute,"template=chain_struct"
            for i=0,nstags-1 do begin
               if where(stname[i] eq tname) eq -1 then begin
                  parent->setvar,'mtname',stname[i]
                  parent->setvar,'mtval',t_chain_struct[0].(i)
                  parent->execute,"template=create_struct(template,"+$
                                  "mtname,mtval)"
               endif
            endfor
            parent->Execute,"shmunmap,'"+chainmmap+"'"
            shmunmap,chainmmap
            
            shmmap,template=template,destroy_segment=1,get_name=chainmmap
            chain_struct=shmvar(chainmmap)
            parent->Execute,"shmmap,'"+chainmmap+"',template=template"
            parent->Execute,"chain_struct=shmvar('"+chainmmap+"')"
         endif
         chain_struct[*]=temporary(t_chain_struct)
      endif
      for n=period,(nstep/step_int)-1 do begin
         for m=0,nchain-1 do begin
            parent->Execute,"child["+strtrim(m,2)+"]->SetVar,"+$
                            "'param_init',chain_struct.params_iter[*,"+ $
                            strtrim(n*step_int-1,2)+","+strtrim(m,2)+"]"
            parent->Execute,"child["+strtrim(m,2)+"]->SetVar,"+$
                            "'param_std',chain_struct.params_std[*,*,"+ $
                            strtrim(m,2)+"]"
            parent->Execute,"outfile='"+current+"/child"+ $
                            strtrim(m,2)+".sav'"
            parent->Execute,"child["+strtrim(m,2)+"]->"+$
                            "Setvar,'outfile',outfile"
            parent->Execute,"child["+strtrim(m,2)+"]->"+$
                            "Setvar,'lnl_init',chain_struct.lnl_iter["+ $
                            strtrim(n*step_int-1,2)+","+strtrim(m,2)+"]"
         endfor
         ;perform a Metropolis-Hastings swap for a pair of adjacent chains
         ;this allows chains to talk to each other and is the basis
         ;for parallel tempering
         if not keyword_set(no_temp) then begin
            rindex=round(satmc_random(/uniform)*(nchain-1))
            repeat begin 
               if satmc_random(/uniform) gt 0.5 then $
                  pair= 0 > (rindex+1) < (nchain-1) else  $
                     pair= 0> (rindex-1) < (nchain-1)
            endrep until pair ne rindex
            ;make sure we're not stuck switching the same chains with
            ;each other
            while (rold eq rindex or rold eq pair) and $
               (pairold eq rindex or pairold eq pair) do begin
               rindex=round(satmc_random(/uniform)*(nchain-1))
               repeat begin 
                  if satmc_random(/uniform) gt 0.5 then $
                     pair= 0 > (rindex+1) < (nchain-1) else  $
                        pair= 0> (rindex-1) < (nchain-1)
               endrep until pair ne rindex
            endwhile
            tpair=10.^(3.5*pair/nchain)-1
            trindex=10.^(3.5*rindex/nchain)-1
            sel=lindgen(step_int)+(n-1)*step_int
            maxpair=max(chain_struct.lnl_iter[sel,pair],psel)
            maxrindex=max(chain_struct.lnl_iter[sel,rindex],rsel)
            alpha=min([1,exp(maxrindex-maxpair)^ $
                       (1./(1+tpair)-(1./(1+trindex)))])
            urand=satmc_random(1,/uniform)
            if abs(maxpair-maxrindex) gt 1.0 and urand le alpha and $
                        swap eq 0 and n gt period then begin
               parent->Execute,"child["+strtrim(rindex,2)+"]->SetVar,"+$
                               "'param_init',chain_struct.params_iter[*,"+ $
                               strtrim(sel[psel],2)+","+ $
                               strtrim(pair,2)+"]"
               parent->Execute,"child["+strtrim(pair,2)+"]->SetVar,"+$
                               "'param_init',chain_struct.params_iter[*,"+ $
                               strtrim(sel[rsel],2)+","+ $
                               strtrim(rindex,2)+"]"
               parent->Execute,"child["+strtrim(rindex,2)+"]->"+$
                               "Setvar,'lnl_init',chain_struct.lnl_iter["+ $
                               strtrim(sel[psel],2)+","+ $
                               strtrim(pair,2)+"]"
               parent->Execute,"child["+strtrim(pair,2)+"]->"+$
                               "Setvar,'lnl_init',chain_struct.lnl_iter["+ $
                               strtrim(sel[rsel],2)+","+ $
                               strtrim(rindex,2)+"]"
               if (pair eq 0 or rindex eq 0) then begin
                  swap=1
                  lastswap=n
               endif else swap=0
               rold=rindex
               pairold=pair
            endif
         endif

         for m=0,nchain-1 do begin
            parent->Execute,"child["+strtrim(m,2)+"]->Execute,"+$
              "'run_mcmc_chain,"+strtrim(step_int,2)+",param_init,"+$
              "lnl_init,sed_struct,grid_struct,goodwavel,goodflux,"+$
              "goodfluxerr,wavelim,fluxlim,upperlim,param_name,norm_arr,"+$
              "nparm_arr,param_min,param_max,param_std,"+$
              "param_arr,lnl_arr,lum_arr,sed,models"+strtrim(m,2)+$
              ",synthesis_model=synthesis_model,outfile=outfile,temp=temp,"+$
              "no_interp=no_interp,filter_struct=filter_struct,dr=dr,"+$
              "ig_lim=ig_lim,priors=priors,prior_mode=prior_mode',/nowait"
         endfor
         ;wait for chains to finish...
         for m=0,nchain-1 do begin
            while parent-> $
               GetVar('child['+strtrim(m,2)+']->status()') eq 1 do wait,0.01
         endfor

         ;...and continue
         ;combine child results
         for m=0,nchain-1 do begin
            if keyword_set(synthesis_model) then begin
               if not tag_exist(chain_struct,'SED') then begin
                  ;in order to add the SED information for each
                  ;parameter step, we'll need to destory and
                  ;recreate CHAIN_STRUCT
                  old_chain_struct=chain_struct
                  parent->Execute,"old_chain_struct=chain_struct"
                  parent->Execute,"sed=child[0]->GetVar('sed')"
                  parent->Execute,"template_sed=rebin(sed,["+ $
                                  "(size(sed,/dim))[0:1],"+$
                                  strtrim(nstep,2)+","+$
                                  strtrim(nchain,2)+"])"
                  parent->Execute,"template=create_struct("+$
                                  "old_chain_struct,'SED',template_sed)"
                  template_sed=parent->GetVar('template_sed')
                  template=create_struct(old_chain_struct,'SED', $
                                         template_sed)
                  parent->Execute,"shmunmap,'"+chainmmap+"'"
                  shmunmap,chainmmap

                  shmmap,template=template,destroy_segment=1, $
                         get_name=chainmmap
                  chain_struct=shmvar(chainmmap)
                  parent->Execute,"shmmap,'"+chainmmap+"',template=template"
                  parent->Execute,"chain_struct=shmvar('"+chainmmap+"')"
               endif 
               parent->Execute,"sed=child["+strtrim(m,2)+"]->GetVar('sed')"
               if n_elements(size(parent->GetVar('sed'),/dim)) ne 1 then begin
                  parent->Execute,"chain_struct.SED[*,*,"+ $
                                  strtrim(n*step_int,2)+":"+$
                                  strtrim((n+1)*step_int-1,2)+","+$
                                  strtrim(m,2)+"]=sed"
               endif

              
            endif
            parent->Execute,"chain_struct.params_iter[*,"+ $
                            strtrim(n*step_int,2)+":"+$
                            strtrim((n+1)*step_int-1,2)+","+ $
                            strtrim(m,2)+"]=child["+strtrim(m,2)+ $
                            "]->GetVar('param_arr')"
            parent->Execute,"chain_struct.lnL_iter["+ $
                            strtrim(n*step_int,2)+":"+$
                            strtrim((n+1)*step_int-1,2)+","+ $
                            strtrim(m,2)+"]=child["+strtrim(m,2)+ $
                            "]->GetVar('lnl_arr')"
            parent->Execute,"chain_struct.lum_iter["+ $
                            strtrim(n*step_int,2)+":"+$
                            strtrim((n+1)*step_int-1,2)+","+ $
                            strtrim(m,2)+"]=child["+strtrim(m,2)+ $
                            "]->GetVar('lum_arr')"
         endfor

         if keyword_set(save_all) then begin
            params_all=[[params_all],[chain_struct.params_iter[*,n*step_int:(n+1)*step_int-1,*]]]
            lnl_all=[lnl_all,chain_struct.lnl_iter[n*step_int:(n+1)*step_int-1,*]]
         endif

         ;if a swap was made with chain 0 (the fiducial 'cold' chain)
         ;check its acceptance rate and alter step size if necessary
         if not keyword_set(no_temp) then begin
            if swap eq 1 and n gt period then begin
               dlnL_iter=chain_struct.lnL_iter[n*step_int: $
                                               (n+1)*step_int-2,0]- $
                         chain_struct.lnL_iter[n*step_int+1: $
                                               (n+1)*step_int-1,0]
               whchange=where(abs(dlnL_iter[0:step_int-2]) gt 0.0,nchange)
               pchange=float(nchange)/float(step_int-1)

               if pchange lt 0.20 or pchange gt 0.26 then begin
                  if pchange le 0.20 then begin
                     cov=correlate(chain_struct.params_iter[*,n*step_int: $
                                                            (n+1)*step_int-1,$
                                                            0],/cov,/double)
                  endif 
                  if pchange gt 0.26 then begin
                     frac=pchange/0.23
                     cov=chain_struct.params_std[*,*,0]*frac
                  endif
                  if pchange lt 0.05 then begin
                     cov=chain_struct.params_std[*,*,0]*0.05
                  endif
                  if array_equal(cov,0) then cov=params_std
                  diag=diag_matrix(cov)
                  diag_pstd=diag_matrix(params_std)
                  sel=where(sqrt(diag) lt sqrt(diag_pstd)*1.e-6 or $
                            sqrt(diag) gt sqrt(diag_pstd)*1.e4,nsel)
                  if nsel ne 0 then cov=params_std
                  chain_struct.params_std[*,*,0]=cov
                  chain_struct.params_iter[*,n*step_int-1,0:1]= $
                     chain_struct.params_iter[*,(n+1)*step_int-1,0:1]
                  chain_struct.lnl_iter[n*step_int-1,0:1]= $
                     chain_struct.lnl_iter[(n+1)*step_int-1,0:1]
                  n--
               endif else begin
                  swap=0
                  if n ge nstep/(step_int*2) then begin 
                     chain_struct.params_iter[*,period*step_int-1,0:1]= $
                        chain_struct.params_iter[*,(n+1)*step_int-1,0:1]
                     chain_struct.lnl_iter[period*step_int-1,0:1]= $
                        chain_struct.lnl_iter[(n+1)*step_int-1,0:1]
                     n=period-1
                     lastswap=0
                     pchange=0.
                  endif
               endelse
                     
               ;ensure convergence is still maintained
               if pchange gt 0.20 and pchange lt 0.26 then begin
                  av1=avg(chain_struct.lnL_iter[n*step_int: $
                                                (n+0.1)*step_int-1,0])
                  av2=avg(chain_struct.lnL_iter[(n+0.5)*step_int: $
                                                (n+1)*step_int-1,0])
                  s=stddev(chain_struct.lnL_iter[(n+0.5)*step_int: $
                                                 (n+1)*step_int-1,0])
                  if s ne 0 then begin
                     if abs(av1-av2)/s gt 2.0 then begin
                        chain_struct.params_iter[*,n*step_int-1,*]= $
                           chain_struct.params_iter[*,(n+1)*step_int-1,*]
                        chain_struct.lnl_iter[n*step_int-1,*]= $
                           chain_struct.lnl_iter[(n+1)*step_int-1,*]
                        n--
                        swap=1
                     endif
                  endif
               endif
            endif
         endif 

         ;After running the tempered chains for a while and no swaps
         ;are being made to the 'cold' chain, reset all chains to 
         ;'cold' chains and re-run.  This provides more samples for 
         ;determining confidence intervals.
         if not keyword_set(no_temp) and $
            n-lastswap eq (nstep/step_int/2) then begin
            no_temp=1
            for m=0,nchain-1 do begin
               parent->Execute,"child["+strtrim(m,2)+"]->Execute,'temp=0'"
               maxindx=(n+1)*step_int-1
               minindx=(period > lastswap)*step_int
               indx=floor(satmc_random(1,/uniform)*(maxindx-minindx))

               chain_struct.params_iter[*,period*step_int-1,m]= $
                  chain_struct.params_iter[*,minindx+indx[0],0]
               chain_struct.lnl_iter[period*step_int-1,m]= $
                  chain_struct.lnl_iter[minindx+indx[0],0]
               chain_struct.params_std[*,*,m]=chain_struct.params_std[*,*,0]
            endfor
            n=period-1
         endif
         
         ;create/update save file
         save,chain_struct,n,filename='parent_chain.sav'
      endfor
   endif
endif else begin
   ;for running without asyncronous chains
   if not keyword_set(adaptive_step) then begin
      if restart le 1 then begin
         n=0
         if restart eq 1 then begin
            restore,'parent_burnin.sav'
            t_chain_struct=chain_struct
            chain_struct=shmvar(chainmmap)
            ;verify all tags exist between saved run and current
            nstags=n_tags(t_chain_struct)
            stname=tag_names(t_chain_struct)
            ntags=n_tags(chain_struct)
            tname=tag_names(chain_struct)
            if ntags ne nstags then begin
               ;tags are missing, will have to destory chainmmap and remake
               template=t_chain_struct
               shmunmap,chainmmap
               shmmap,template=template,destroy_segment=1,get_name=chainmmap
               chain_struct=shmvar(chainmmap)
            endif
            chain_struct[*]=temporary(t_chain_struct)
         endif

         while (pchange le 0.20 or pchange ge 0.26) do begin
            ;if doing parallel tempering, make potential swaps between
            ;the fudical cold and warm chains
            if not keyword_set(no_temp) and t_period eq 2 then begin
               t_period=0
               rindex=round(satmc_random(/uniform)*(nchain-1))
               while rindex eq 0 do $
                  rindex=round(satmc_random(/uniform)*(nchain-1))
               maxpair=max(chain_struct.lnl_iter[0:2*step_int-1,0],psel)
               maxrindex=max(chain_struct.lnl_iter[0:2*step_int-1,rindex], $
                             rsel)
               trindex=10.^(3.5*rindex/nchain)-1
               alpha=min([1,exp(maxrindex-maxpair)^(1.-1./(1+trindex))])
    
               urand=satmc_random(1,/uniform)
               if urand le alpha then begin
                  params_old=chain_struct.params_iter[*,rsel,rindex]
                  lnl_old=chain_struct.lnl_iter[rsel,rindex]
                  chain_struct.params_iter[*,0,rindex]= $
                     chain_struct.params_iter[*,psel,0]
                  chain_struct.params_iter[*,0,0]=params_old
                  chain_struct.lnl_iter[0,rindex]=chain_struct.lnl_iter[psel,0]
                  chain_struct.lnl_iter[0,0]=lnl_old
                  chain_struct.params_std[*,*,0]*=0.001
               endif
            endif else t_period++
         
            for m=0,nchain-1 do begin
               param_init=chain_struct.params_iter[*,0,m] 
               lnl_init=chain_struct.lnL_iter[0,m]

               if keyword_set(no_temp) then temp=0 else $
                  temp=10.^(3.5*m/nchain)-1

               run_mcmc_chain,step_int*2,param_init,lnl_init,sed_struct, $
                 grid_struct,goodwavel,goodflux,goodfluxerr,wavelim,fluxlim, $
                 upperlim,param_name,norm_arr,nparm_arr,param_min,param_max, $
                 chain_struct.params_std[*,*,m],param_arr,lnl_arr,lum_arr, $
                 sed,models,synthesis_model=synthesis_model, $
                 outfile='mcmc_chain.sav',temp=temp,no_interp=no_interp, $
                 filter_struct=filter_struct,dr=dr,ig_lim=ig_lim, $
                 priors=priors,prior_mode=prior_mode

               if keyword_set(synthesis_model) then begin
                  if not tag_exist(chain_struct,'SED') then begin
                     s=size(sed,/dim)
                     chain_struct=create_struct(chain_struct,'SED', $
                                                rebin(sed,[s,nchain]))
                  endif 
                  chain_struct.SED[*,*,0:step_int*2-1,m]=sed
               endif
               chain_struct.lnL_iter[0:step_int*2-1,m]=lnl_arr
               chain_struct.params_iter[*,0:2*step_int-1,m]=param_arr
               chain_struct.lum_iter[0:2*step_int-1,m]=lum_arr
            endfor

            if keyword_set(save_all) then begin
               params_all=[[params_all], $
                           [chain_struct.params_iter[*,0:2*step_int-1,*]]]
               lnl_all=[lnl_all,chain_struct.lnl_iter[0:2*step_int-1,*]]
            endif
            
            ;check acceptance rate of MCMC steps
            params_std_old=chain_struct.params_std
            pchange_ch=dblarr(nchain)
            for m=0,nchain-1 do begin
               dlnL_iter=chain_struct.lnL_iter[step_int:2*step_int-1,m]- $
                 chain_struct.lnL_iter[step_int+1:2*step_int,m]
               whchange=where(abs(dlnL_iter[0:step_int-2]) gt 0.0,nchange)
               pchange_ch[m]=float(nchange)/float(step_int-1)
            endfor

            if keyword_set(no_temp) then begin
               pchange=avg(pchange_ch)
               if pchange le 0.20 or pchange ge 0.26 then begin
                  for m=0,nchain-1 do begin
                     cov=chain_struct.params_std[*,*,m]
                     if pchange_ch[m] gt 0.26 then begin
                        diag_pstd=diag_matrix(chain_struct.params_std[*,*,m])
                        pstd=chain_struct.params_std[*,*,m]
                        if array_equal(diag_matrix(diag_pstd),pstd) then begin
                           cov=correlate(chain_struct.params_iter[*,step_int:$
                                                                  2*step_int-1,m],/cov,/double)
                        endif else begin
                           frac=pchange_ch[m]/0.23
                           cov=chain_struct.params_std[*,*,m]*frac
                        endelse
                     endif
                     if pchange_ch[m] lt 0.20 and pchange_ch[m] gt 0.05 then begin
                        cov=correlate(chain_struct.params_iter[*,step_int: $
                                                               2*step_int-1, $
                                                               m],/cov,/double)
                     endif
                     if pchange_ch[m] lt 0.05 then begin
                        cov=chain_struct.params_std[*,*,m]*0.05
                     endif
                     if array_equal(cov,0) then cov=params_std
                     diag=diag_matrix(cov)
                     diag_pstd=diag_matrix(params_std)
                     sel=where(sqrt(diag) lt sqrt(diag_pstd)*1.e-6 or $
                               sqrt(diag) gt sqrt(diag_pstd)*1.e4,nsel)
                     if nsel ne 0 then cov=params_std
                     chain_struct.params_std[*,*,m]=cov
                  endfor                    
               endif
            endif else begin
               pchange=pchange_ch[0]
               if pchange le 0.20 or pchange ge 0.26 then begin
                  for m=0,nchain-1 do begin 
                     cov=chain_struct.params_std[*,*,m]
                     if pchange_ch[m] le 0.20 and pchange_ch[m] gt 0.05 then begin
                        cov=correlate(chain_struct.params_iter[*,step_int: $
                                                  2*step_int-1,m], $
                                      /cov,/double)
                     endif 
                     if pchange_ch[m] lt 0.05 then begin
                        cov=chain_struct.params_std[*,*,m]*0.05
                     endif
                     if pchange_ch[m] gt 0.26 then begin
                        diag_pstd=diag_matrix(chain_struct.params_std[*,*,m])
                        pstd=chain_struct.params_std[*,*,m]
                        if array_equal(diag_matrix(diag_pstd),pstd) then begin
                           cov=correlate(chain_struct.params_iter[*,step_int:$
                                                                  2*step_int-1,m],/cov,/double)
                        endif else begin
                           frac=pchange_ch[m]/0.23
                           cov=chain_struct.params_std[*,*,m]*frac
                        endelse
                     endif
                     if array_equal(cov,0) then cov=params_std
                     diag=diag_matrix(cov)
                     diag_pstd=diag_matrix(params_std)
                     sel=where(sqrt(diag) lt sqrt(diag_pstd)*1.e-6 or $
                               sqrt(diag) gt sqrt(diag_pstd)*1.e4,nsel)
                     if nsel ne 0 then cov=params_std
                     chain_struct.params_std[*,*,m]=cov
                  endfor
               endif
            endelse

            sel=where(pchange_ch eq 0,nsel)
            if nsel ne 0 then pchange=0.

            ;determine current maximum likelihood value
            maxL_ch=dblarr(nchain)
            for m=0,nchain-1 do $
               maxL_ch[m]=max(chain_struct.lnl_iter[where(chain_struct.lnl_iter[*,m] ne 0.0),m])
            diffl=(maxL_ch-max(maxL_ch))

            ;test for convergence
            if pchange ge 0.20 and pchange le 0.26 then begin
               ;Geweke dianogstic
               if not keyword_set(no_temp) then begin
                  av1=avg(chain_struct.lnL_iter[0:0.2*step_int-1,0])
                  av2=avg(chain_struct.lnL_iter[step_int:2*step_int-1,0])
                  s=stddev(chain_struct.lnL_iter[step_int:2*step_int-1,0])
                  if s ne 0 then begin
                     if abs(av1-av2)/s gt 2.0 then pchange=0.
                  endif
               endif else begin
                  for m=0,nchain-1 do begin
                     av1=avg(chain_struct.lnL_iter[0:0.2*step_int-1,m])
                     av2=avg(chain_struct.lnL_iter[step_int:2*step_int-1,m])
                     s=stddev(chain_struct.lnL_iter[step_int: $
                                                    2*step_int-1,m])
                     if s ne 0 then begin
                        if abs(av1-av2)/s gt 2.0 then pchange=0.
                     endif
                  endfor
               endelse
            endif

            ;check for bad chains, chains that either are unable
            ;move about (flat likelihood space) or too far from
            ;currently maximum likelihood value
            if pchange ge 0.20 and pchange le 0.26 then begin
               badchain=where(pchange_ch lt thres or $
                              abs(diffl) gt 1.e2,nbc,com=goodchain)
            endif else begin
               badchain=where(pchange_ch lt thres and diffl ne 0,nbc, $
                              com=goodchain)
            endelse
            
            if nbc ne 0 then begin
               pchange=0.
               if nbc eq nchain then begin
                  for i=0,nbc-1 do begin
                     chain_struct.params_iter[*,0,badchain[i]]= $
                        satmc_random(n_elements(param_max),/uniform)* $
                            (param_max-param_min)+param_min

                     chain_struct.lnl_iter[0,badchain[i]]=0
                     chain_struct.params_std[*,*,badchain[i]]= $
                        params_std_old[*,*,badchain[i]]
                  endfor
               endif else begin
                  ;use the best estimate from other chains to determine
                  ;new initial step
                  params_good=chain_struct.params_iter[*,step_int: $
                                                       2*step_int-1, $
                                                       goodchain]
                  if n_elements(goodchain) eq 1 then begin
                     for i=0,nbc-1 do begin
                        rind=round(satmc_random(/uniform)*(step_int-1))
                        chain_struct.params_iter[*,0,badchain[i]]= $
                           params_good[*,rind]
                        chain_struct.lnl_iter[0,badchain[i]]= $
                           chain_struct.lnl_iter[rind,goodchain]
                        chain_struct.params_std[*,*,badchain[i]]= $
                           chain_struct.params_std[*,*,goodchain]
                     endfor
                  endif else begin
                     for i=0,nbc-1 do begin
                        rind=round(satmc_random(/uniform)*(step_int-1))
                        rchain=round(satmc_random(/uniform)* $
                                     (n_elements(goodchain)-1))
                        chain_struct.params_iter[*,0,badchain[i]]= $
                           params_good[*,rind,rchain]
                        chain_struct.params_std[*,*,badchain[i]]= $
                           chain_struct.params_std[*,*,goodchain[rchain]]
                        chain_struct.lnl_iter[0,badchain[i]]= $
                          chain_struct. lnl_iter[rind,goodchain[rchain]]
                     endfor
                  endelse
                  chain_struct.params_iter[*,0,goodchain]= $
                     chain_struct.params_iter[*,2*step_int-1,goodchain]
                  chain_struct.lnl_iter[0,goodchain]= $
                     chain_struct.lnl_iter[2*step_int-1,goodchain]
               endelse
            endif else begin
               chain_struct.params_iter[*,0,*]= $
                  chain_struct.params_iter[*,2*step_int-1,*]
               chain_struct.lnl_iter[0,*]= $
                  chain_struct.lnl_iter[2*step_int-1,*]
            endelse
            ;create/update savefile
           save,chain_struct,pchange,filename='parent_burnin.sav'
           n++
         endwhile
      endif
      period=2
   endif else period=0
    ;Now do the full chains, seperate into blocks of STEP_INT length
    ;to set a restart point in case of chain failure.
    ;Save the combined results as well as individual chains.
   if restart le 2 then begin
      swap=0
      i=1
      if restart eq 2 then begin
         restore,'parent_chain.sav'
         t_chain_struct=chain_struct
         chain_struct=shmvar(chainmmap)
         ;verify all tags exist between saved run and current
         nstags=n_tags(t_chain_struct)
         stname=tag_names(t_chain_struct)
         ntags=n_tags(chain_struct)
         tname=tag_names(chain_struct)
         if ntags ne nstags then begin
            ;tags are missing, will have to destory chainmmap and remake
            template=t_chain_struct
            shmunmap,chainmmap
            shmmap,template=template,destroy_segment=1,get_name=chainmmap
            chain_struct=shmvar(chainmmap)
            endif
         chain_struct[*]=temporary(t_chain_struct)
      endif
      for n=period,(nstep/step_int)-1 do begin
         ;perform a Metropolis-Hastings swap for a pair of adjacent chains
         ;this allows chains to talk to each other and is the basis
         ;for parallel tempering
         if not keyword_set(no_temp) then begin
            rindex=round(satmc_random(/uniform)*(nchain-1))
            
            repeat begin 
               if satmc_random(/uniform) gt 0.5 then $
                  pair= 0 > (rindex+1) < (nchain-1) else  $
                     pair= 0> (rindex-1) < (nchain-1)
            endrep until pair ne rindex
            ;make sure we're not stuck switching the same chains with
            ;each other
            while (rold eq rindex or rold eq pair) and $
               (pairold eq rindex or pairold eq pair) do begin
               rindex=round(satmc_random(/uniform)*(nchain-1))
               repeat begin 
                  if satmc_random(/uniform) gt 0.5 then $
                     pair= 0 > (rindex+1) < (nchain-1) else  $
                        pair= 0> (rindex-1) < (nchain-1)
               endrep until pair ne rindex
            endwhile
            tpair=10.^(3.5*pair/nchain)-1
            trindex=10.^(3.5*rindex/nchain)-1
            sel=lindgen(step_int)+(n-1)*step_int
            maxpair=max(chain_struct.lnl_iter[sel,pair],psel)
            maxrindex=max(chain_struct.lnl_iter[sel,rindex],rsel)
            alpha=min([1,exp(maxrindex-maxpair)^(1./(1+tpair)-1./(1+trindex))])

            urand=satmc_random(1,/uniform)
            if urand le alpha and swap eq 0 and n gt period then begin
               params_old=chain_struct.params_iter[*,sel[rsel],rindex]
               lnl_old=chain_struct.lnl_iter[sel[rsel],rindex]
               chain_struct.params_iter[*,n*step_int-1,rindex]= $
                  chain_struct.params_iter[*,sel[psel],pair]
               chain_struct.params_iter[*,n*step_int-1,pair]=params_old
               chain_struct.lnl_iter[n*step_int-1,rindex]= $
                  chain_struct.lnl_iter[sel[psel],pair]
               chain_struct.lnl_iter[n*step_int-1,pair]=lnl_old
               if (pair eq 0 or rindex eq 0) then begin
                  swap=1 
                  lastswap=n
               endif else swap=0
               rold=rindex
               pairold=pair
            endif
         endif

         for m=0,nchain-1 do begin
            param_init=chain_struct.params_iter[*,(n*step_int)-1,m]
            lnl_init=chain_struct.lnl_iter[(n*step_int)-1,m]

            if keyword_set(no_temp) then temp=0 else temp=10.^(3.5*m/nchain)-1
            
            run_mcmc_chain,step_int,param_init,lnl_init,sed_struct, $
              grid_struct,goodwavel,goodflux,goodfluxerr,wavelim,fluxlim, $
              upperlim,param_name,norm_arr,nparm_arr,param_min,param_max, $
              chain_struct.params_std[*,*,m],param_arr,lnl_arr,lum_arr, $
              sed,models,synthesis_model=synthesis_model, $
              outfile='mcmc_chain.sav',temp=temp,no_interp=no_interp, $
              filter_struct=filter_struct,dr=dr,ig_lim=ig_lim, $
              priors=priors,prior_mode=prior_mode

            if keyword_set(synthesis_model) then begin
               if not tag_exist(chain_struct,'SED') then begin
                  s=size(sed,/dim)
                  chain_struct=create_struct(chain_struct,'SED', $
                                             rebin(sed,[s,nchain]))
               endif 
               chain_struct.SED[*,*,n*step_int:(n+1)*step_int-1,m]=sed
            endif
            chain_struct.lnL_iter[n*step_int:(n+1)*step_int-1,m]=lnl_arr
            chain_struct.params_iter[*,n*step_int:(n+1)*step_int-1,m]=param_arr
            chain_struct.lum_iter[n*step_int:(n+1)*step_int-1,m]=lum_arr
         endfor

         if keyword_set(save_all) then begin
            params_all=[[params_all],[chain_struct.params_iter[*,n*step_int:(n+1)*step_int-1,*]]]
            lnl_all=[lnl_all,chain_struct.lnl_iter[n*step_int:(n+1)*step_int-1,*]]
         endif

         ;if a swap was made with chain 0 (the fiducial 'cold' chain)
         ;check its acceptance rate and alter step size if necessary
         if not keyword_set(no_temp) then begin
            if swap eq 1 and n gt period then begin
               dlnL_iter=chain_struct.lnL_iter[n*step_int: $
                                               (n+1)*step_int-2,0]- $
                         chain_struct.lnL_iter[n*step_int+1: $
                                               (n+1)*step_int-1,0]
               whchange=where(abs(dlnL_iter[0:step_int-2]) gt 0.0,nchange)
               pchange=float(nchange)/float(step_int-1)
               if pchange lt 0.20 or pchange gt 0.26 then begin
                  if pchange lt 0.20 then begin
                     cov=correlate(chain_struct.params_iter[*,n*step_int: $
                                                            (n+1)*step_int-1,$
                                                            0],/cov,/double)
                  endif
                  if pchange lt 0.05 then begin
                     cov=chain_struct.params_std[*,*,0]*0.05
                  endif
                  if pchange gt 0.26 then begin
                     frac=pchange/0.23
                     cov=chain_struct.params_std[*,*,0]*frac
                  endif
                  if array_equal(cov,0) then cov=params_std
                  diag=diag_matrix(cov)
                  diag_pstd=diag_matrix(params_std)
                  for i=0,n_elements(diag)-1 do begin
                     if sqrt(diag[i]) lt sqrt(diag_pstd[i])*1.e-6 and $
                        sqrt(diag[i]) lt sqrt(diag_pstd[i])*1.e4 then begin
                        cov[i,i]=diag_pstd[i]
                     endif
                  endfor
                  chain_struct.params_std[*,*,0]=cov 
                  
                  chain_struct.params_iter[*,n*step_int-1,0:1]= $
                     chain_struct.params_iter[*,(n+1)*step_int-1,0:1]
                  chain_struct.lnl_iter[n*step_int-1,0:1]= $
                     chain_struct.lnl_iter[(n+1)*step_int-1,0:1]
                  n--
               endif else begin
                  swap=0
                  if n gt nstep/(step_int*2) then begin 
                     chain_struct.params_iter[*,period*step_int-1,0:1]= $
                        chain_struct.params_iter[*,(n+1)*step_int-1,0:1]
                     chain_struct.lnl_iter[period*step_int-1,0:1]= $
                        chain_struct.lnl_iter[(n+1)*step_int-1,0:1]
                     n=period-1
                     pchange=0.
                     lastswap=0
                  endif
               endelse
            endif
         
            ;ensure convergence is still maintained
            if pchange gt 0.20 and pchange lt 0.26 then begin
               av1=avg(chain_struct.lnL_iter[n*step_int: $
                                             (n+0.1)*step_int-1,0])
               av2=avg(chain_struct.lnL_iter[(n+0.5)*step_int: $
                                             (n+1)*step_int-1,0])
               s=stddev(chain_struct.lnL_iter[(n+0.5)*step_int: $
                                              (n+1)*step_int-1,0])
               if s ne 0 then begin
                  if abs(av1-av2)/s gt 2.0 then begin
                     chain_struct.params_iter[*,n*step_int-1,*]= $
                        chain_struct.params_iter[*,(n+1)*step_int-1,*]
                     chain_struct.lnl_iter[n*step_int-1,*]= $
                        chain_struct.lnl_iter[(n+1)*step_int-1,*]
                     n--
                     swap=1
                  endif
               endif
            endif
            i++
         endif
               
         ;After running the tempered chains for a while and no swaps
         ;are being made to the 'cold' chain, reset all chains to 
         ;'cold' chains and re-run.  This provides more samples for 
         ;determining confidence intervals.
         if not keyword_set(no_temp) and $
            n-lastswap eq (nstep/step_int/2) then begin
            no_temp=1
            for m=1,nchain-1 do begin
               indx=floor(satmc_random(/uniform)*(n-lastswap)*step_int)+$
                    lastswap*step_int-1
               chain_struct.params_iter[*,period*step_int-1,m]= $
                  chain_struct.params_iter[*,indx,0]
               chain_struct.lnl_iter[period*step_int-1,m]= $
                  chain_struct.lnl_iter[indx,0]
               chain_struct.params_std[*,*,m]=chain_struct.params_std[*,*,0]
            endfor
            n=period-1
         endif

        ;create/update save file
         save,chain_struct,n,filename='parent_chain.sav'
      endfor
   endif
endelse
;clean up MCMC products
if keyword_set(synthesis_model) then begin
   if not keyword_set(no_async) then begin
      for m=0,nchain-1 do begin
         str='rm -f MODELS'+strtrim(m,2)+'.*'
         parent->SetVar,'str',str
         parent->Execute,"child["+strtrim(m,2)+"]->SetVar,'str',str"
         parent->Execute,"child["+strtrim(m,2)+"]->Execute,'spawn,str'"
      endfor
      ;kill child processes
      parent->Execute,'obj_destroy,child'
      obj_destroy,parent
   endif else spawn,'rm -f MODELS*'
endif else begin
   if not keyword_set(no_async) then begin
     ;kill child processes
      parent->Execute,'obj_destroy,child'
      obj_destroy,parent
   endif
endelse

;keep track of accepted steps from each chain
lnL_good=[-999.]
params_good=replicate(-999.,nparams+1)
lum_good=[-999.]
if tag_exist(chain_struct,'sed') then sed_good=chain_struct.sed[*,*,0,0]

dlnl_iter=dblarr(nstep)

for m=0,nchain-1 do begin
   dlnL_iter=[0,chain_struct.lnL_iter[1:*,m]- $
              chain_struct.lnL_iter[0:nstep-1,m]]
   whchange=where(abs(dlnL_iter) gt 0.0 and $
                  lindgen(nstep) gt step_int*2,nchange)
   if nchange ne 0 then begin
      lnL_good=[lnL_good,chain_struct.lnL_iter[whchange,m]]
      params_good=[[params_good],[chain_struct.params_iter[*,whchange,m]]]
      lum_good=[lum_good,chain_struct.lum_iter[whchange,m]]
      if tag_exist(chain_struct,'sed') then $
         sed_good=[[[sed_good]],[[chain_struct.sed[*,*,whchange,m]]]]
   endif
endfor
lnL_good=lnL_good[1:*]
params_good=params_good[*,1:*]
if keyword_set(save_all) then begin
   lnl_all=lnl_all[1:*,*]
   params_all=params_all[*,1:*,*]
endif
;restore parameter space from compressed log-space
if ndr ne 0 then begin
   params_good[dr,*]=10^(params_good[dr,*])
   if keyword_set(save_all) then params_all[dr,1:*,*]=10^(params_all[dr,1:*,*])
   param_max[dr]=10^(param_max[dr])
   param_min[dr]=10^(param_min[dr])
   p0=0
   for l=0,nsed-1 do begin
      if norm_arr[l] then begin
         if nparm_arr[l] gt 1 then begin
            gs=where(dr ge p0+1 and dr le nparm_arr[l]-1,nsel)
            if nsel gt 0 then grid_struct.(l)[*,dr[gs]-1]=10^(grid_struct.(l)[*,dr[gs]-1])
         endif
      endif else begin
         gsel=where(dr ge p0 and dr le nparm_arr[l]-1,nsel)
         if nsel gt 0 then grid_struct.(l)[*,dr[gsel]]=10^(grid_struct.(l)[*,dr[gsel]])
      endelse
      p0+=nparm_arr[l]
   endfor
endif
lum_good=lum_good[1:*]

chain_struct_good=create_struct('lnl_good',lnl_good,'params_good', $
                                params_good,'lum_good',lum_good)
if tag_exist(chain_struct,'sed') then begin
   sed_good=sed_good[*,*,1:*]
   chain_struct_good=create_struct(chain_struct_good,'sed_good',sed_good)   
endif

;determine 68% (1-sigma) confidence intervals for model parameters
maxL=max(lnL_good,whmx)
params_max=params_good[*,whmx]
lum_max=lum_good[whmx]
params_1sig=fltarr(nparams+1,2)
ss=reverse(sort(lnL_good))

ind68=floor(erf(1./sqrt(2.))*n_elements(lnL_good))
lnL68=lnL_good[ss[ind68]]
wh1sig=where(lnL_good[ss] ge lnL68,nwh1sig)
params_low=min(params_good[*,ss[wh1sig]],dim=2)
params_high=max(params_good[*,ss[wh1sig]],dim=2)
params_1sig[*,0]=params_high-params_max
params_1sig[*,1]=params_max-params_low
lum_1sig=[max(lum_good[ss[wh1sig]])-lum_max,lum_max-min(lum_good[ss[wh1sig]])]

;plot the best fit model, find the closest grid model to maximum 
;likelihood parameters
param0=0
for i=0,nsed-1 do begin
    if i eq 0 and keyword_set(synthesis_model) then begin
        ;grab best fit and 1-sigma confidence region from synthesized models
        wave=chain_struct_good.sed_good[*,0,whmx]
        mflux=chain_struct_good.sed_good[*,1,whmx]

        npt=n_elements(wave)
        best_mod=create_struct('lum',0.0d,'mm',0.0d,'sed',dblarr(npt,2), $
                               'sed_low',dblarr(npt,2),'sed_high', $
                               dblarr(npt,2))
        if norm_arr[0] then flux*=params_max[param0]
        best_mod.sed=[[wave],[mflux]]
        
        lnu=mflux*50.^2./(1.0118)/(1.13e-8)/c
        best_mod.lum=abs(integral(c/wave,lnu))
        best_mod.mm=interpol(best_mod.sed[*,1],best_mod.sed[*,0],1100)
    endif else begin
       if norm_arr[i] then begin
          if nparm_arr[i] gt 1 then begin
             mod_best=find_nearest_grid(grid_struct.(i), $
                      params_max[param0+1:param0+nparm_arr[i]-1],intp=ind)
             mod_best=md_interpolate(mod_best[0,*],mod_best[1:*,*], $
                                     round(ind))
          endif else mod_best=[0]
          best_wave=sed_struct.(i)[mod_best].wavesim
          best_flux=sed_struct.(i)[mod_best].fluxsim*params_max[param0]
          best_lum=sed_struct.(i)[mod_best].lum*params_max[param0]
       endif else begin
          mod_best=find_nearest_grid(grid_struct.(i), $
                   params_max[param0:param0+nparm_arr[i]-1],intp=ind)
          mod_best=md_interpolate(mod_best[0,*],mod_best[1:*,*], $
                                  round(ind))
          
          best_wave=sed_struct.(i)[mod_best].wavesim
          best_flux=sed_struct.(i)[mod_best].fluxsim
          best_lum=sed_struct.(i)[mod_best].lum
       endelse
       npt=n_elements(best_wave)
       best_mod=create_struct('lum',0.0d,'mm',0.0d,'sed',dblarr(npt,2))
       best_mod.sed=shift_sed([[best_wave],[best_flux]],params_max[nparams])
       best_mod.mm=interpol(best_mod.sed[*,1],best_mod.sed[*,0],1100)
       best_mod.lum=best_lum
    endelse
    
    if i eq 0 then begin
        best_mods=create_struct(sed_arr[i],best_mod)
    endif else begin
        best_mods=create_struct(best_mods,sed_arr[i],best_mod)
    endelse
    param0=total(nparm_arr[0:i])
endfor

if not keyword_set(no_plot) then begin
   colorlegend=['black','white','blue','red','green']
   tvlct,255*[0,1,0,1,0],255*[0,1,0,0,1],255*[0,1,1,0,0]

   ;co-add SEDS
   flux_best=best_mods.(0).sed[*,1]
   wave_best=best_mods.(0).sed[*,0]
   
   for i=1,nsed-1 do begin
       ;determine overlap region and interpolate to corsest binning
       max=min([max(wave_best),max(best_mods.(i).sed[*,0])])
       min=max([min(wave_best),min(best_mods.(i).sed[*,0])])
       sel=where(wave_best le max and wave_best ge min,nsel)
       sels=where(best_mods.(i).sed[*,0] le max and $
                  best_mods.(i).sed[*,0] ge min,nsels)

       if nsel gt nsels then begin
           mod_int=interpol(flux_best[sel],wave_best[sel], $
                            best_mods.(i).sed[sels,0])
           flux_int=best_mods.(i).sed[sels,1]+mod_int
           wave_int=best_mods.(i).sed[sels,0]
       endif else begin
           mod_int=interpol(best_mods.(i).sed[sels,1], $
                            best_mods.(i).sed[sels,0],wave_best[sel])
           flux_int=flux_best[sel]+mod_int
           wave_int=wave_best[sel]
       endelse
       
       ;check if there are points outside the interpolation region
       ;and add them back to co-added SED
       sel=where(wave_best lt min,nsel)
       if nsel gt 0 then begin
           wave_int=[wave_int,wave_best[sel]]
           flux_int=[flux_int,flux_best[sel]]
       endif
       sel=where(best_mods.(i).sed[*,0] lt min,nsel)
       if nsel gt 0 then begin
           wave_int=[wave_int,best_mods.(i).sed[sel,0]]
           flux_int=[flux_int,best_mods.(i).sed[sel,1]]
       endif
       sel=where(wave_best gt max,nsel)
       if nsel gt 0 then begin
           wave_int=[wave_best[sel],wave_int]
           flux_int=[flux_best[sel],flux_int]
       endif
       sel=where(best_mods.(i).sed[*,0] gt max,nsel)
       if nsel gt 0 then begin
           wave_int=[best_mods.(i).sed[sel,0],wave_int]
           flux_int=[best_mods.(i).sed[sel,1],flux_int]
       endif

       wave_best=wave_int
       flux_best=flux_int
   endfor
   
   ;plot co-added and component SEDs
   case axes of
      1: begin
         oplot,wave_best,flux_best,linestyle=0
         for i=0,nsed-1 do oplot,best_mods.(i).sed[*,0], $
                                 best_mods.(i).sed[*,1],color=i+2,linestyle=i+1
      end
      2: begin
         lfl_best=flux_best*3.e-9/wave_best
         oplot,wave_best,lfl_best,linestyle=0
         for i=0,nsed-1 do begin
            wave_mods=best_mods.(i).sed[*,0]
            f_mods=best_mods.(i).sed[*,1]
            lfl_mods=f_mods*3.e-9/wave_mods
            oplot,wave_mods,lfl_mods,color=i+2,linestyle=i+1   
         endfor
      end
   endcase
   if keyword_set(psfile) then begin
      device,/close
      set_plot,'x'
   endif
endif

mm_tot=0.
for i=0,nsed-1 do mm_tot=mm_tot+best_mods.(i).mm

;print results to screen (and outfile if set)
message,/cont,'Best fit model parameters:'
for i=0,nparams-1 do begin
   message,/cont,param_name[i]+': '+strtrim(string(params_max[i]),2)+$
           '+'+strtrim(string(params_1sig[i,0]),2)+$
           '-'+strtrim(string(params_1sig[i,1]),2)
endfor
if keyword_set(findz) then begin
   message,/cont,'Redshift:'+string(params_max[nparams],'(F6.4)')+$
           '+'+string(params_1sig[nparams,0],'(F6.4)')+$
           '-'+string(params_1sig[nparams,1],'(F6.4)')
endif else message,/cont,'Redshift: '+string(params_max[nparams],'(F6.4)')
log_lum=floor(alog10(lum_max))
message,/cont,'Total luminosity (L_sun): '+ $
        string(lum_max*10.^(-log_lum),'(F5.2)')+'+'+ $
        string(lum_1sig[0]*10.^(-log_lum),'(F5.2)')+'-'+ $
        string(lum_1sig[1]*10.^(-log_lum),'(F5.2)')+'e'+ $
        strtrim(log_lum,2)
if nsed ge 2 then begin
   message,/cont,'SED1 contribution (bol,1.1mm): '+$
           string(best_mods.(0).lum/lum_max,'(F4.2)')+$
           ','+string(best_mods.(0).mm/(mm_tot),'(F4.2)')
endif
message,/cont,'Maximum Likelihood value:'+strtrim(maxL,2)

if keyword_set(outfile) then begin
   openw,lun,outfile+'.fit',/get
   printf,lun,'Best fit model parameters:'
   for i=0,nparams-1 do begin
      printf,lun,param_name[i]+': '+strtrim(string(params_max[i]),2)+$
             '+'+strtrim(string(params_1sig[i,0]),2)+$
             '-'+strtrim(string(params_1sig[i,1]),2)
   endfor
   if keyword_set(findz) then begin
      printf,lun,'Redshift:'+string(params_max[nparams],'(F6.4)')+$
             '+'+string(params_1sig[nparams,0],'(F6.4)')+$
             '-'+string(params_1sig[nparams,1],'(F6.4)')
   endif else printf,lun,'Redshift:'+string(params_max[nparams],'(F6.4)')
   log_lum=floor(alog10(lum_max))
   printf,lun,'Total luminosity (L_sun): '+ $
          string(lum_max*10.^(-log_lum),'(F5.2)')+'+'+ $
          string(lum_1sig[0]*10.^(-log_lum),'(F5.2)')+'-'+ $
          string(lum_1sig[1]*10.^(-log_lum),'(F5.2)')+'e'+ $
          strtrim(log_lum,2)
   if n_elements(stype) ne 0 then begin
      log_lx=floor(alog10(lx_max))
      printf,lun,'Estimated X-ray luminosity: '+$
             string(lx_max*10.^(-log_lx),'(F5.2)')+'+'+ $
             string(lx_1sig[0]*10.^(-log_lx),'(F5.2)')+'-'+ $
             string(lx_1sig[1]*10.^(-log_lx),'(F5.2)')+'e'+ $
             strtrim(log_lx,2)
   endif
   if nsed ge 2 then begin
      printf,lun,'SED1 contribution (bol,1.1mm): '+$
             string(best_mods.(0).lum/lum_max,'(F4.2)')+$
             ','+string(best_mods.(0).mm/(mm_tot),'(F4.2)')
   endif
   printf,lun,'Maximum Likelihood value:'+strtrim(maxL,2)
   free_lun,lun

   ;save chain results
   save,filename=outfile+'.sav',chain_struct,chain_struct_good, $
     lum_max,maxL,params_max,param_name,step_int, $
     param_max,param_min,best_mods,wavel,flux,fluxerr,upperlim
   if keyword_set(save_all) then $
      save,filename=outfile+'fullchains.sav',params_all,lnl_all
endif
;finish clean up
shmunmap,chainmmap
spawn,'rm parent_burnin.sav parent_chain.sav'
if keyword_set(no_async) then spawn,'rm mcmc_chain.sav' else begin
    for m=0,nchain-1 do spawn,'rm child'+strtrim(m,2)+'.*'
    spawn,'rm parent.log'
    ;stop virtual X-server if one was made
    if keyword_set(use_xvfb) then begin
       setenv,'DISPLAY='+current_display
       spawn,'kill -9 '+strtrim(xvfb_pid)
       free_lun,unit
    endif
endelse

return
end
