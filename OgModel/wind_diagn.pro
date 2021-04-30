pro calcola_equil_3,zeta,ioneq_game,height_in,loro

; Calculates the equilibirum ioneq once T is specified,
; using the same method as in MAKE_IONEQ_ALL.PRO
;
; interpolates over heights rather than temperature

common rate,rec_all,ion_all,ttt,hhh

ioneq_min=1.0d-99
rate_limite=1.0d-100

nions=zeta+1

ioneq_game=dblarr(nions)
ion_rate=dblarr(nions)
rec_rate=dblarr(nions)

for i=0,nions-1 do begin

   if min(ion_all(loro,i)) gt rate_limite then $
;      g1=spline(alog10(height(loro)),alog10(ion_all(loro,i)),alog10(height_in)) else g1=-300.0d0
      g1=spline(alog10(hhh(loro)),alog10(ion_all(loro,i)),alog10(height_in)) else g1=-300.0d0
   ion_rate(i)=g1

   if min(rec_all(loro,i)) gt rate_limite then $
;      g1=spline(alog10(height(loro)),alog10(rec_all(loro,i)),alog10(height_in)) else g1=-300.0d0
      g1=spline(alog10(hhh(loro)),alog10(rec_all(loro,i)),alog10(height_in)) else g1=-300.0d0
   rec_rate(i)=g1

endfor

ioneq_game(0)=1.0d0

for i=1,nions-1 do begin

   peo=total(ion_rate(0:i-1))-total(rec_rate(1:i))
   ioneq_game(i)=10^peo

endfor

ioneq_game=ioneq_game/total(ioneq_game)
questi=where(ioneq_game lt ioneq_min)
if questi(0) ge 0 then ioneq_game(questi)=0.0d0
ioneq_game=ioneq_game/total(ioneq_game)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wind_diagn,zeta,smooth=smooth,remote=remote,quiet=quiet,shift=shift,t_start=t_start,$
                    t_init=t_init,v_start=v_start,dens=dens,tempi=tempi,no_plot=no_plot,$
                    no_save=no_save,file_model=file_model,serial_calc=serial_calc,istart=istart,$
                    non_maxw=non_maxw,nthe=nthe,ntht=ntht,here=here,double_dens=double_dens,$
                    check_equil=check_equil,photoion=photoion,serial_phot=serial_phot,$
                    begin_ioneq=begin_ioneq,h_start=h_start,correct_phot=correct_phot,$
                    flare=flare,fileflare=fileflare,dist_flare=dist_flare,day_photoion=day_photoion,$
                    source_flare_dist=source_flare_dist,vel_factor=vel_factor,same_flux=same_flux,$
                    fileout=fileout

;  Note: The model input file must include the following quantities
;
;       mass:        in g
;       height:      in solar radii
;       velocity:    in km/s
;       temp:        in K
;
;  INPUTS
;
;  zeta          Atomic number:  1 (H), 2 (He) etc'
;
;  KEYWORDS
;
;  begin_ioneq   Uses an user-defined initial charge state distribution. Needs a distribution 
;                as an input
;
;  check_equil   Checks whether time scales are fast enough to have equil
;
;  correct_phot  Uses a user-defined correction factor to multiply the incident photoionizing radiation
;                then used for the calculation
;
;  day_photoion  Day when we want to calculate the photoion rate (asks directly if this
;                keyword is not given). Format: 26-03-2007'
;
;  dens          Asks for the method to use to calculate Ne. By default it assumes and H-only plasma
;
;  dist_flare    Distance from the limb when the flare goes off
;
;  double_dens   Doubles the value of the electron density everywhere, to test the effects
;
;  file_model    Interactively selects the model file (keyword depends on /remote and /no_plot)
;
;  flare         Includes flare photoionization. Asks as an input: r_flare (in solar radii),
;                files with radiance in f(time) (output of converti_spettro.pro). These can also
;                be given using the keywords /dist_flare and /fileflare. Assumes /photoion.
;
;  fileflare     Flare emission file
;
;  h_start       Sets the initial height where to start the calculation
;
;  here          Read files in the directory where you are working
;
;  istart        From which model in the series does the program start to calculate (starts from zero)
;
;  no_plot       Writes results on file, but does not do the plot
;
;  no_save       Saves results every 500 points
;
;  non_maxw      Adds second maxwellian with T2=ntht including nthe % electrons
;
;  nthe          % of electrons belonging to the second Maxwellian
;
;  ntht          T of second Maxwellian (in K)'
;
;  photoion      Includes photoionization as an ionization process. Asks for an input file 
;                for the EUV irradiance. May be used with /day_photoion.
;
;  quiet         Does not provide any detail of how the calculation is going
;
;  remote        Does not write and does not plot anything (to run from remote)
;
;  same_flux     (together with vel_factor) decreases/increases the density when the velocity is
;                changed by vel_factor, in order to keep the mass flux constant.
;
;  serial_calc   Asks for a save file with many models in an IDL structure, and runs for each of them
;
;  serial_phot   Carries out the calculation for a series of photoion values, once an input file
;                of type TIMED/SEE is given. Must be used together with /photoion
;
;  shift         Shifts velocity curve to make it start at the first point of hhh where the
;                winds starts with v=0.5 km/s
;
;  smooth        Smooths the temperature profile to remove numerical noise in the original model
;
;  source_flare_ Parameter controlling how distant is the flare from the wind source region.
;     dist       Needed to calculate the flare dilution factor
;
;  tempi         Tells how long it took to carry out each task
;
;  t_start       Starting T
;
;  t_init        Temperature to be used to calculates the initial plasma charge state composition
;                (note that it is different from t_start, which is the T from which the calculation
;                is started)
;
;  v_start       Initial velocity
;
;  vel_factor    multiplies the velocity by that factor, to speed or slow down the model
;
; HISTORY
;
;     V.2  - 29-Jan-2013 - Cambiato il modo di calcolare i rec/ion rates
;
;     V.3  -  9-Jan-2014 - Aggiunta la possibilita' di fare calcoli seriali
;
;     V.4  - 10-Jan-2014 - Messo un tetto limite al numero di punti della soluzione, 
;                          che, una volta raggiunto, forza il programma a tenere a 
;                          mente solo meta' dei punti calcolati
;
;     V.5  - 23-Jan-2014 - Levato svariati piccoli bug; messo che rifiuta tutte le
;                          altezze la cui differenza dallla precedente e' minore dello
;                          0.005%.
;
;     V.6  - 28-Jan-2014 - Messo il criterio di Smith & Hughes 2010 (usando Huges & Helmand 
;                          1985 per il calcolo) per determinare se il plasma e' in equilibrio
;                          o no a una certa altezza
;
;     V.7  - 11-Feb-2014 - Messo la possibilita' di leggere files nella directory in cui
;                          si lavora, quando si usa senza la widget
;
;     V.8  - Jan/Feb-2014- Messa la possibilita' di includere una seconda Maxwelliana
;
;     V.9  -  6-Nov-2014 - Messo la possibilita' di raddoppiare la densita', come test
;
;     V.10 - 16-Jan-2015 - Aggiunta la photoionization alla ionization totale
;
;     V.11 -  5-Feb-2015 - Aggiunte le keyowrds di poter iniziare da una T o H a scelta,
;                          e di usare una distribuzione iniziale arbitraria da dargli in
;                          input
;
;     V.12 - Early 2016  - Aggiunto flare photoionization, con un solo flare
;
;     V.13 - 11-Jun-2016 - Salvato il tempo di percorrenza del vento
;
;     V.14 - 17-Jun-2016 - Messa la distanza wind source-flare location come input opzionale,
;                          senza input e' 0.2 solar radii (0.8 nel programma)
;
;     V.15 - 27-Jun-2017 - Messa la keyword vel_factor, che moltiplica la velocita' per il 
;                          fattore vel_factor (guarda cosa cambia se il vento e' piu lento
;                          del modello, o piu' veloce; aggiunto anche same_flux x la densita'

common rate,rec_all,ion_all,ttt,hhh

month=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
elem=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S',$
      'Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
stage=['I     ','II    ','III   ','IV     ','V     ','VI    ','VII   ','VIII  ','IX    ',$
       'X     ','XI    ','XII   ','XIII   ','XIV   ','XV    ','XVI   ','XVII  ','XVIII ',$
       'XIX   ','XX    ','XXI   ','XXII   ','XXIII ','XXIV  ','XXV   ','XXVI  ','XXVII ',$
       'XXVIII','XXIX  ','XXX   ','XXXI   ']

; 1 - Checks the input

; 1.1 - Arbitrary parameters

rsun=6.97d5

ioneq_min=1.0d-07        ; Minimum accepted ioneq at any time except first step
ioneq_min_start=1.0d-08  ; Minimum accepted ioneq at first step
delta_min=0.1            ; Maximum accepted fractional variation
ioneq_threshold=1.0d-6   ; Minimum ioneq for which the DELTA_MIN limit is valid
delta_ioneq_min=1.0d-4   ; Minimum difference beyond which evolution is assumed to have stopped
size_result_limit=30000.0 ; Size limit for results, to avoid that big calculations kill memory

if not keyword_set(flare) then sfd=0.0d0
if keyword_set(flare) then begin

   if not keyword_set(source_flare_dist) then sfd=0.8d0
   if keyword_set(source_flare_dist) then sfd=double(1.0-source_flare_dist)

endif

; 1.2 - Gets very upset if the element is not given as input

if n_params() lt 1 then begin

   print,''
   print,'No, You have to tell me AT LEAST which element you want to do!'
   print,''
   print,'     IDL> wind_diagn,zeta(,/smooth,/remote,/quiet,/shift,/density,t_start=t_start,t_init=t_init,$'
   print,'                           v_start=v_start,/dens,/tempi,/no_plot,/no_save,/file_model,/non_maxw,'   
   print,'                           nthe=nthe,ntht=ntht,/here,/double_dens,/check_equil,/photoion,/serial_phot,$'
   print,'                           correct_phot=correct_phot,begin_ioneq=begin_ioneq,h_start=h_start,$'
   print,'                           flare=flare,source_flare_dist=source_flare_dist,vel_factor=vel_factor,same_flux=same_flux)'
   print,'                           fileout=fileout)'
   print,''
   print,' See code header for keyword meaning'
   return

endif

; 1.3 - In case of a second Maxwellian, checks that all parameters are given

if keyword_set(non_maxw) then begin

   if not keyword_set(nthe) or not keyword_set(ntht) then begin

      print,''
      print,'No: if you want to use a second Maxwellian you need to give me the following parameters:'
      print,''
      print,'   ntht          Second Maxwellian temperature (in K)'
      print,'   nthe          Percentage of total electrons belonging to second Maxwellian'
      print,''
      return

   endif

endif

; 1.4 - Selects element and defines number of ions

nions=zeta+1

if not keyword_set(remote) then print,'Faccio l`elemento : '+$
      strcompress(string(elem(zeta-1)),/remove_all)+' che ha '+$
      strcompress(string(nions),/remove_all)+' ioni'

elem_name=strlowcase(elem(zeta-1))
zion2name,zeta,1,meo
peo=str2arr(meo,'_')
peo=peo(0)
if not keyword_set(fileout) then filesave='results_'+strcompress(peo,/remove_all)+'_3.save'
if keyword_set(fileout) then filesave = fileout
filewrite='ioneq_start_'+strcompress(peo,/remove_all)+'_3.dat'

; 1.5 - Selects theoretical model

if not keyword_set(serial_calc) then begin

   if not keyword_set(file_model) then begin
   
      if keyword_set(remote) then begin

         file_model=''
         read,'Give me the wind model (tipo model_cranmer_2007/cranmer_model_ch_2007.save): ',file_model
         if not keyword_set(here) then file_model='~/giuda/lepri/jacob/fast_wind/'+file_model

      endif
   
      if not keyword_set(remote) then begin

         if keyword_set(no_plot) then begin

            file_model=''
            read,'Give me the wind model (tipo model_cranmer_2007/cranmer_model_ch_2007.save): ',file_model
            if not keyword_set(here) then file_model='~/giuda/lepri/jacob/fast_wind/'+file_model

         endif
         if not keyword_set(no_plot) then $
            file_model=dialog_pickfile(filter='*.save',title='Choose the wind model (e.g slow_wind.save)')

      endif

   endif

   if not keyword_set(remote) then begin 

      print,''
      print,'File scelto : ',file_model
      print,''

   endif

   restore,file_model

   height=double(height)
   temp=double(temp)
   velocity=double(velocity)
   mass=double(mass)

endif

if keyword_set(vel_factor) then begin

   velocity=velocity*vel_factor

endif

if keyword_set(serial_calc) then begin

   filein=''
   print,'The input file needs to be an IDL save file with a structure named MODEL'
   read,'Tell me the model file I need to read: ',filein
   restore,filein
   nmodels=n_elements(model)

   npoints_res=size_result_limit

   str={   hhh_used: fltarr(npoints_res), $
           ttt_used: fltarr(npoints_res), $
           ddd_used: fltarr(npoints_res), $
           vvv_used: fltarr(npoints_res), $
           ioneq_evol: fltarr(nions,npoints_res), $
           ioneq_equil: fltarr(nions,npoints_res), $
           total_time: fltarr(npoints_res), $
           time_travel: fltarr(npoints_res)};, $
;           photoion_used: strarr(nions,npoints_res), $
;           recomb_used: strarr(nions,npoints_res), $
;           ioniz_used: strarr(nions,npoints_res), $
;           avg_ioneq: strarr(nions,npoints_res)  }
   results=replicate(str,nmodels)

   tagli=fltarr(nmodels)

   if not keyword_set(istart) then istart=0.0

   riparti:

   height=reform(model[istart].height)
   temp=reform(model[istart].temp)
   velocity=reform(model[istart].velocity)
   mass=reform(model[istart].mass)

   height=double(height)
   temp=double(temp)
   velocity=double(velocity)
   mass=double(mass)

   taglia=0.0

   loro_start=where(height gt 0)

   height=height(loro_start)
   temp=temp(loro_start)
   velocity=velocity(loro_start)
   mass=mass(loro_start)

endif

if keyword_set(serial_phot) then begin

   filein=''
   print,'Note: the input file needs to be an IDL save file of the TIMED/SEE type (~/giuda/see/*.save)'
   read,'Give me the input irradiance file (full address): ',filein
   restore,filein
   nmodels=n_elements(year_see)

   npoints_res=size_result_limit

   str={   hhh_used: fltarr(npoints_res), $
           ttt_used: fltarr(npoints_res), $
           ddd_used: fltarr(npoints_res), $
           vvv_used: fltarr(npoints_res), $
           ioneq_evol: fltarr(nions,npoints_res), $
           ioneq_equil: fltarr(nions,npoints_res), $
           total_time: fltarr(npoints_res), $
           time_travel: fltarr(npoints_res)};, $
;           photoion_used: strarr(nions,npoints_res), $
;           recomb_used: strarr(nions,npoints_res), $
;           ioniz_used: strarr(nions,npoints_res), $
;           avg_ioneq: strarr(nions,npoints_res)  }

   results=replicate(str,nmodels)

   tagli=fltarr(nmodels)

endif

; 1.6 - Calculates electron doensity, assuming that the average mass is 0.61 m_p even if it is not fully ionized.

m_prot=1.6726d-24
m_avg=0.61             ; Not entirely correct, since it assumes full ionization

if not keyword_set(dens) then begin

   ratio=proton_dens(alog10(temp),/hydrogen,abund_file='/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2011_caffau.abund',ioneq_file='/usr/local/ssw/packages/chianti/dbase/ioneq/chianti.ioneq')
   density=mass/m_prot/ratio

endif

if keyword_set(dens) then begin

   print,'How do you want to calculate Ne: '
   print,''
   print,'    1 - Using the mean particle mass'
   print,'    2 - Using like Cranmer 2005 a fully ionized plasma with 5% helium'
   print,''
   read,'Your choice: ',scelta

   if scelta eq 1 then density=mass/m_prot/m_avg
   if scelta eq 2 then density=mass*0.916/m_prot

endif

if keyword_set(double_dens) then density=2.0d0*density

if keyword_set(vel_factor) then begin

   if keyword_set(same_flux) then density=density/vel_factor

endif

; 2 - Selects starting point

; Assume that everything with: Ne > d_limit and T < t_limit and H < h_limit is in ionization equilibrium
; Also, start by default at 10,000 K, to avoid problems with charge exchange

if keyword_set(v_start) then v_limit=v_start
if not keyword_set(v_start) then v_limit=min(velocity)

;d_limit=1.0d11 Original value
d_limit=1.0d12
if not keyword_set(h_start) then h_limit=1.002
if keyword_set(h_start) then begin

   if h_start le min(height) then h_limit=1.002
   if h_start gt min(height) then h_limit=h_start
   

endif

if not keyword_set(t_start) then t_limit=1e4

if keyword_set(t_start) then begin

   if t_start lt 1e4 then t_limit=1e4
   if t_start ge 1e4 then t_limit=t_start

endif

loro=where(height ge h_limit and temp gt t_limit and density lt d_limit and velocity ge v_limit)

; 2.1 - Defines the T, Ne, V to be used

nmodel=n_elements(loro)

ddd=double(density(loro))
ttt=double(temp(loro))
hhh=double(height(loro))

if not keyword_set(shift) then vvv=double(velocity(loro))
if keyword_set(shift) then begin

   v1=double(velocity(loro))
   queste=where(velocity gt 0.5)
   vvv=double(velocity(queste))
   
endif

; 2.2 - Removes points at the same or very similar height

diff_check=abs(hhh-shift(hhh,1))/hhh*100.0
levami_h=where(diff_check le 5e-3)
if levami_h(0) ge 0 then remove,levami_h,hhh,ddd,ttt,vvv

; 2.3 - Continues definition of T, Ne, V to be used

if keyword_set(smooth) then begin 
   ttt_orig=ttt
   ttt=smooth(ttt,3)
   if not keyword_set(remote) then begin
      if not keyword_set(no_plot) then begin
         plot,hhh,ttt,psym=2,xr=[min(hhh),1.02],charsize=1.5
         oplot,hhh,ttt
      endif
   endif
endif

; 2.4 - Gets ion_rates e rec_rates

ntemp=n_elements(ttt)
rec_all=dblarr(ntemp,nions)
ion_all=dblarr(ntemp,nions)

for j=0,nions-1 do begin

   zion2name,zeta,j+1,gname
   if j gt 0 then begin
      rpeo=recomb_rate(gname,ttt)
      if keyword_set(non_maxw) then begin
         rpeo1=recomb_rate(gname,ntht)
         rpeo2=(1-nthe/100.0d0)*rpeo+nthe/100.0d0*rpeo1(0)
         rpeo=rpeo2
      endif
      rec_all(*,j)=rpeo(*)
   endif
   if j lt nions-1 then begin
      ipeo=ioniz_rate(gname,ttt)
      if keyword_set(non_maxw) then begin
         ipeo1=ioniz_rate(gname,ntht)
         ipeo2=(1-nthe/100.0d0)*ipeo+nthe/100.0d0*ipeo1(0)
         ipeo=ipeo2
      endif
      ion_all(*,j)=ipeo(*)
   endif

endfor

; 2.4.1 Reads the photoionization rates and input spectrum in case one is interested

if keyword_set(photoion) then begin

   ; Reads in the photoionization cross section parameters and utilizes them to calculate the cross section

   stronzata=''
   photoion_params=dblarr(30,31,9)
   peo=dblarr(9)

   openr,1,'photoion_params.txt'
   readf,1,stronzata

   while not eof(1) do begin

      readf,1,a1,a2,peo,format='(i3,i4,9e10.3)'
      photoion_params(a1-1,a2,*)=peo(*)

   endwhile

   close,1

   ; Chooses energies in order to always have Wvl_min=1.0 A 

   peo=reform(photoion_params(zeta-1,*,0))
   peo=peo(where(peo gt 0))
   
   n_energies=30000.0                          ; Old energy definition

   sigma_photoion=dblarr(nions-1,n_energies)
   energia_sigma_photoion=dblarr(nions-1,n_energies)

   for i=0,nions-2 do begin

      params=reform(photoion_params(zeta-1,i,*))

      ethr=params(0)
      e0=params(2)
      y0=params(7)
      y1=params(8)
      yw=params(6)
      ppp=params(5)
      ya=params(4)
      sigma0=params(3)

      deltae_needed=12398.0d0-ethr
      step_needed=deltae_needed/n_energies
      energia=ethr+step_needed*findgen(n_energies);/10.0   ; Old energy definition

      xxx=energia/e0-y0
      yyy=sqrt(xxx^2+y1^2)
      fff=((xxx-1)^2+yw^2)*yyy^(0.5*ppp-5.5)*((1+sqrt(yyy/ya))^(-ppp))
      sigma_photoion(i,*)=sigma0*fff           ; Old energy definition
      energia_sigma_photoion(i,*)=energia      ; Old energy definition

   endfor
   
   wvl_photoion=12398.0d0/energia_sigma_photoion
   dwvl_photoion=0.*wvl_photoion
   sigma_photoion=sigma_photoion*1.0d-18    ; in cm^2

   ; Reads in the TIMED/SEE irradiance spectrum

   restore,'timed_see_daily_irradiance_all.save'
   flux_see=flux_see*0.1     ; converts nm-1 to A-1
   units_see_new='erg cm-2 s-1 sr-1 A-1'
   wvl_see_down=wvl_see-5.0d0
   wvl_see_up=wvl_see+5.0d0
   date_limits=minmax(year_see)

   if not keyword_set(serial_phot) then begin

      if not keyword_set(serial_calc) then begin

         ; Asks for a date for the irradiance spectrum and checks it

         input_date:

         date_see=''
         if not keyword_set(day_photoion) then read,'Which date do you want to look at (e.g. 26-03-2007, must be between 1-2-2002 and 30-11-2017): ',date_see
         if keyword_set(day_photoion) then date_see=day_photoion

         pepe=str2arr(date_see,'-')
         meme='00'+pepe(1)+pepe(0)
         date2doy,meme,doy_see
         year_see_sel=double(pepe(2))+double(doy_see)/365.25

         if year_see_sel lt date_limits(0) or year_see_sel gt date_limits(1) then begin
      
            print,''
            print,'I do not have this date in the TIMED SEE data set'
            print,''
            goto,input_date
      
         endif

      endif

      if keyword_set(serial_calc) then begin

         if istart eq 0 then begin

            ; Asks for a date for the irradiance spectrum and checks it

            input_date_2:

            date_see=''
            if not keyword_set(day_photoion) then read,'Which date do you want to look at (e.g. 26-03-2007): ',date_see
            if keyword_set(day_photoion) then date_see=day_photoion

            pepe=str2arr(date_see,'-')
            meme='00'+pepe(1)+pepe(0)
            date2doy,meme,doy_see
            year_see_sel=double(pepe(2))+double(doy_see)/365.25

            if year_see_sel lt date_limits(0) or year_see_sel gt date_limits(1) then begin

               print,''
               print,'I do not have this date in the TIMED SEE data set'
               print,''
               goto,input_date_2

            endif

         endif

      endif

   endif

   if keyword_set(serial_phot) then begin

      if not keyword_set(istart) then istart=0.0

      riparti_phot:

      year_see_sel=year_see(istart)

   endif

   ; Identifies the date

   pepe=abs(year_see-year_see_sel)
   questo_see=where(pepe eq min(pepe))

   ; Checks if a flare needs to be included - if yes, the flux of the given date is considered bkg

   flux_see_inc=dblarr(nions-1,n_energies)
   r_flare=-1.0d0
      
   if keyword_set(flare) then begin

      if not keyword_set(fileflare) then $
         fileflare=dialog_pickfile(filter='~/giuda/lepri/flare_phot/*.save',tit='Flare radiance file')

      restore,fileflare
      flux_flare=flux
      wvl_flare=wvl
      time_flare=time_final*60.0d0            ; Time in seconds
      radius_flare=radius_source              ; In solar radii

      tmin=min(temp_final)                    ; Eliminates any emission before/after flare
      these_ones=where(temp_final gt tmin)
      flux_flare=flux_flare(*,these_ones)
      time_flare=time_flare(these_ones)
      time_flare=time_flare-time_flare(0)

      ntimes_flare=n_elements(time_flare)
      flux_see_inc_flare=dblarr(nions-1,n_energies,ntimes_flare)

      if keyword_set(dist_flare) then r_flare=dist_flare
      if not keyword_set(dist_flare) then read,'At what distance from the limb (at 1.0) does the flare erupt (in solar radii): ',r_flare

   endif

   ; Calculates the actual photionization rate

   for j=0,nions-2 do begin

      dwvl_photoion(j,*)=abs(wvl_photoion(j,*)-shift(wvl_photoion(j,*),1))
      dwvl_photoion(j,0)=0.0

      for i=0,n_energies-1 do begin

         lei=where(wvl_see_down le wvl_photoion(j,i) and wvl_see_up gt wvl_photoion(j,i))
         perepeo=flux_see(lei,questo_see(0))
         if min(perepeo) ge 0.0 then flux_see_inc(j,i)=flux_see(lei,questo_see)
         if min(perepeo) lt 0.0 then flux_see_inc(j,i)=0.0d0

      endfor
 
      if keyword_set(flare) then begin

         for i=0,n_energies-2 do begin

            lei=where(wvl_flare le wvl_photoion(j,i) and wvl_flare gt wvl_photoion(j,i+1))
            if lei(0) ge 0 then perepeo=total(flux_flare(lei,*),1)
            if min(perepeo) ge 0.0 then flux_see_inc_flare(j,i,*)=perepeo(*)

         endfor

      endif

   endfor

   if keyword_set(correct_phot) then flux_see_inc=flux_see_inc*correct_phot
   if keyword_set(correct_phot) and keyword_set(flare) then flux_see_inc_flare=flux_see_inc_flare*correct_phot

   integrand=flux_see_inc*sigma_photoion*wvl_photoion*dwvl_photoion*6.326d8  ; 6.326d8=4*!pi*1e-8/hc
   photoion_rate=total(integrand,2)     ; Photoionization rate at Sun's surface

   if keyword_set(flare) then begin
   
      integrand_flare=0.d0*flux_see_inc_flare
      for k=0,ntimes_flare-1 do integrand_flare(*,*,k)=flux_see_inc_flare(*,*,k)*sigma_photoion*wvl_photoion*dwvl_photoion*6.326d8  ; 6.326d8=4*!pi*1e-8/hc
      photoion_rate_flare=total(integrand_flare,2)     ; Flare photoionization rate at Sun's surface

   endif

   if not keyword_set(remote) and not keyword_set(flare) then print,year_see_sel,photoion_rate
   if not keyword_set(remote) and keyword_set(flare) then begin
      print,year_see_sel,photoion_rate
      print,year_see_sel,photoion_rate_flare(*,0)
      print,year_see_sel,photoion_rate_flare(*,0)/photoion_rate
   endif

endif

; 2.5 - Outputs the data for the starting point 

if not keyword_set(remote) then begin

   print,'Starting point for the calculation: '
   print,''
   print,'   H  = '+string(hhh(0),'(f12.3)')+' R_sun'
   print,'   V  = '+string(vvv(0),'(f12.3)')+' Km/s' 
   print,'   Ne = '+string(ddd(0),'(e12.3)')+' cm-3' 
   print,'   T  = '+string(ttt(0),'(f12.1)')+' K'
   print,''

endif

; 2.6 - Calculates initial charge state distribution, using the method in the MAKE_IONEQ_ALL.PRO

if not keyword_set(begin_ioneq) then begin

   if not keyword_set(t_init) then begin
   
      temp_in=ttt(0)
      height_in=hhh(0)
      diff=abs(hhh-hhh(0))
      lui=where(diff eq min(diff))
      zeo=minmax(lui)
      loro=[zeo(0)-1-findgen(2),lui,zeo(1)+1+findgen(2)]
      loro=loro(sort(loro))
      loro=loro(where(loro ge 0 and loro le ntemp-1))
      help, loro
      calcola_equil_3,zeta,ioneq_game,height_in,loro
   
   endif

   if keyword_set(t_init) then begin

      temp_in=t_init
      height_in=hhh(0)
      diff=abs(hhh-hhh(0))
      lui=where(diff eq min(diff))
      zeo=minmax(lui)
      loro=[zeo(0)-1-findgen(2),lui,zeo(1)+1+findgen(2)]
      loro=loro(sort(loro))
      loro=loro(where(loro ge 0 and loro le ntemp-1))
      help, loro
      calcola_equil_3,zeta,ioneq_game,height_in,loro

   endif

endif

if keyword_set(begin_ioneq) then begin

   n_roba=n_elements(begin_ioneq)
   if n_roba eq nions then ioneq_game=begin_ioneq/total(begin_ioneq)
   if n_roba ne nions then begin

      print,'' 
      print,'You need to give me an input with the exact number of ions allowed by the element!'
      print,''
      return

   endif

endif

ioneq_start=ioneq_game

openw,1,filewrite

printf,1,'Starting point for the calculation:'
printf,1,''
printf,1,'   H  = '+string(hhh(0),'(f12.3)')+' R_sun'
printf,1,'   V  = '+string(vvv(0),'(f12.3)')+' Km/s'
printf,1,'   Ne = '+string(ddd(0),'(e12.3)')+' cm-3'
if not keyword_set(t_init) then printf,1,'   T  = '+string(ttt(0),'(f12.1)')+' K'
if keyword_set(t_init) then printf,1,'   T  = '+string(t_init,'(f12.1)')+' K'
printf,1,''

for j=0,nions-1 do printf,1,j,strupcase(elem_name)+' '+stage(j),ioneq_start(j),format='(i3,a15,e12.3)'


close,1

; 3 - Now the real calculation begins

tempo_inizio=systime(1)

; 3.1 - Final variables

ioneq_evol=dblarr(nions)
dioneq_evol=dblarr(nions)
ioneq_equil=dblarr(nions)
new_ioneq=dblarr(nions)
vvv_used=0.0d0
hhh_used=0.0d0
ddd_used=0.0d0
ttt_used=0.0d0
total_time=0.0d0
time_travel=0.0d0
photoion_used=dblarr(nions)
recomb_used=dblarr(nions)
ioniz_used=dblarr(nions)
recomb_used=0.0d0
ioniz_used=0.0d0
n_test=0.0d0
ion_stage=1.0+findgen(nions)
avg_ioneq=0.0d0
taglia=0.0
equil=0.0           
equilibrium=0.0      ; 0 is out of equilibrium, 1 is in equilibrium

; 3.2 - Sets the starting point for the Runge-Kutta iteration

i=0L

old_ioneq=ioneq_start
avg_ioneq0=total(ion_stage)*ioneq_start

teo0=ttt(0)
deo0=ddd(0)
veo0=vvv(0)
heo0=hhh(0)
delta_h=10.0d0    ; in km

passare=0.0
negativo=0.0
grosso=0.0
quale_neg=0.0

flare_trigger=0.0      ; 0: flare not started yet     1: flare has started

; 3.3 - Starts Runge Kutta

rifare:

equilibrium=0.0

; 3.4 - Corrects the step if necessary

if passare eq 1 then begin

   delta_h=delta_h/2.0d0
   time_giro_1=systime(1)

endif

if passare eq 0 then begin

   time_start=systime(1)
   time_giro_1=systime(1)
   negativo=0.0
   grosso=0.0
   quale_neg=0.0

endif

; 3.5 - Calculates plasma parameters

time=delta_h/veo0

heo1=heo0+(veo0*time/2.0)/rsun
if heo1 gt max(hhh) then begin
   passare=1.0
   goto,rifare
endif
loro=where(hhh gt heo0-(heo1-1.0)*2.0 and hhh lt heo0+(heo1-1.0)*2.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo1-1.0)*10.0 and hhh lt heo0+(heo1-1.0)*10.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo1-1.0)*100.0 and hhh lt heo0+(heo1-1.0)*100.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo1-1.0)*1000.0 and hhh lt heo0+(heo1-1.0)*1000.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo1-1.0)*10000.0 and hhh lt heo0+(heo1-1.0)*10000.0)
if n_elements(loro) lt 3 then begin
   delta_h=delta_h*10d0
   goto,rifare
endif
veo1=10^(spline(alog10(hhh(loro)),alog10(vvv(loro)),alog10(heo1)))
deo1=10^(spline(alog10(hhh(loro)),alog10(ddd(loro)),alog10(heo1)))
teo1=10^(spline(alog10(hhh(loro)),alog10(ttt(loro)),alog10(heo1)))

heo2=heo0+(veo0*time)/rsun
if heo2 gt max(hhh) then begin
   passare=1.0
   goto,rifare
endif
loro=where(hhh gt heo0-(heo2-1.0)*2.0 and hhh lt heo0+(heo2-1.0)*2.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo2-1.0)*10.0 and hhh lt heo0+(heo2-1.0)*10.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo2-1.0)*100.0 and hhh lt heo0+(heo2-1.0)*100.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo2-1.0)*1000.0 and hhh lt heo0+(heo2-1.0)*1000.0)
if loro(0) lt 0 then loro=where(hhh gt heo0-(heo2-1.0)*10000.0 and hhh lt heo0+(heo2-1.0)*10000.0)
veo2=10^(spline(alog10(hhh(loro)),alog10(vvv(loro)),alog10(heo2)))
deo2=10^(spline(alog10(hhh(loro)),alog10(ddd(loro)),alog10(heo2)))
teo2=10^(spline(alog10(hhh(loro)),alog10(ttt(loro)),alog10(heo2)))

if strcompress(string(deo2),/remove_all) eq 'NaN' then begin
   print,'Careful: I have a model with two points at the same height'
   stop
endif
if strcompress(string(deo2),/remove_all) eq '-NaN' then begin
   print,'Careful: I have a model with two points at the same height'
   stop
endif

; 3.6 - Refuses the steps that lead to temperature changes larger than 10%

if abs((teo2-teo0)/teo0) gt 0.1 then begin

   passare=1.0
   goto,rifare     

endif

; 3.7 - Defines ion and rec rates (de-facto initializing them)

rec_rate0=dblarr(nions)
rec_rate1=dblarr(nions)
rec_rate2=dblarr(nions)
ion_rate0=dblarr(nions)
ion_rate1=dblarr(nions)
ion_rate2=dblarr(nions)

; 3.8 - Remove ionization stages with too low abundances, which can only further decrease

if i eq 0 then begin    ; Uses only the charge states larger than ioneq_min_start for the lower ions

   loro=where(old_ioneq gt ioneq_min_start)
   stage_min=min(loro)
   stage_max=nions-1
   npoints_result=1

endif

if i eq 1 then begin

   loro_buoni=where(ioneq_evol(*) gt ioneq_min)
   stage_min=min(loro_buoni)
   stage_max=nions-1

endif

if i gt 1 then begin   ; Discriminates between rising and falling T

   deltat=ttt_used(npoints_result-1)-ttt_used(npoints_result-2)
   loro_buoni=where(ioneq_evol(*,npoints_result-2) gt ioneq_min)
   scarto=where(ioneq_evol(*,npoints_result-2) le ioneq_min and dioneq_evol(*,npoints_result-2) lt 0)

   if scarto(0) ge 0 and scarto(0) lt min(loro_buoni) then begin

       if deltat ge 0 then begin

          stage_min=min(loro_buoni)
          stage_max=float(nions-1)

       endif
     
       if stage_min ge stage_max then stage_min=min(loro_buoni)

   endif

endif

; 3.9 - Determines the points where to calculation ionization and recombination rates for this step

t_stage_1=systime(1)

diff0=abs(hhh-heo0)
lui0=where(diff0 eq min(diff0))
zeo0=min(lui0)
diff1=abs(hhh-heo2)
lui1=where(diff1 eq min(diff1))
zeo1=max(lui1)
nzeo=zeo1-zeo0+11
if nzeo lt 15 then nzeo=15.0          ; Ensures we have at least 3 points for the following interpolation
loro=zeo0-1.0+5.0*findgen(nzeo/5.0)   ; Chooses a sparser version of the interpolation range, to increase speed

if min(loro) lt 0 then begin

   questesi=where(loro gt 0)
   questeno=where(loro le 0)
   loro=[0,loro(questesi)]
   loro=loro(sort(loro)) 

endif
if max(loro) ge n_elements(hhh) then goto,festa  ; Has reached the limit

; 3.10 - Determines the flare photoionizing flux

if keyword_set(photoion) then begin

   if not keyword_set(flare) then photoion_rate_used=photoion_rate      ; There is no flare

   if keyword_set(flare) then begin                                     ; There is a flare

      if heo0 lt r_flare then begin                                     ; Flare hasn't started yet

         photoion_rate_used=photoion_rate
         ph_flare=dblarr(nions)

      endif
         
      if heo0 ge r_flare then begin

         if flare_trigger eq 0 then begin                               ; Flare is starting NOW

            photoion_rate_used=photoion_rate+reform(photoion_rate_flare(*,0))
            flare_trigger=flare_trigger+1.0
            flare_time=0.0
            ph_flare=dblarr(nions)

         endif

         if flare_trigger eq 1 then begin                               ; Flare has already started

            flare_time=flare_time+time

            if flare_time le max(time_flare) then begin                 ; Flare is ongoing

               diff_time=abs(flare_time-time_flare)
               loro_b=where(diff_time lt 1000)          ; Chooses near times only for interpolation
               ph_flare=dblarr(nions)
               for iflare=0,nions-2 do $
                      ph_flare(iflare)=10^spline((time_flare(loro_b)),alog10(reform(photoion_rate_flare(iflare,loro_b))),(flare_time))
               photoion_rate_used=photoion_rate+ph_flare


            endif

            if flare_time gt max(time_flare) then begin                 ; Flare is over

               photoion_rate_used=photoion_rate 
               ph_flare=dblarr(nions)

            endif

         endif

      endif

   endif

endif

; 3.11 - Calculates ionization and recombination rates for that step

rate_limite=1.0d-100

for j=stage_min,stage_max do begin

   if keyword_set(photoion) then begin
      dilution=0.5*(1-sqrt(1-1/[heo0,heo1,heo2]^2))
      pho_rate0=photoion_rate*dilution(0) 
      pho_rate1=photoion_rate*dilution(1) 
      pho_rate2=photoion_rate*dilution(2) 
      if keyword_set(flare) then begin
         dilution_flare=0.5*(1-sqrt(1-radius_flare^2/([heo0,heo1,heo2]-sfd)^2))
         pho_rate0=photoion_rate*dilution(0)+ph_flare*dilution_flare(0)
         pho_rate1=photoion_rate*dilution(1)+ph_flare*dilution_flare(1)
         pho_rate2=photoion_rate*dilution(2)+ph_flare*dilution_flare(2)
      endif
   endif else pho_rate0=0.0d0
   

   if j eq 0 then begin

      if min(ion_all(loro,j)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(ion_all(loro,j)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      if not keyword_set(photoion) then begin
         ion_rate0(j)=peo(0)*deo0
         ion_rate1(j)=peo(1)*deo1
         ion_rate2(j)=peo(2)*deo2
      endif
      if keyword_set(photoion) then begin
         ion_rate0(j)=peo(0)*deo0+pho_rate0(j)
         ion_rate1(j)=peo(1)*deo1+pho_rate1(j)
         ion_rate2(j)=peo(2)*deo2+pho_rate2(j)
      endif

      if min(rec_all(loro,j+1)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(rec_all(loro,j+1)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      rec_rate0(j+1)=peo(0)*deo0
      rec_rate1(j+1)=peo(1)*deo1
      rec_rate2(j+1)=peo(2)*deo2

   endif

   if j gt 0 and j lt nions-1 then begin

      if min(ion_all(loro,j)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(ion_all(loro,j)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      if not keyword_set(photoion) then begin
         ion_rate0(j)=peo(0)*deo0
         ion_rate1(j)=peo(1)*deo1
         ion_rate2(j)=peo(2)*deo2
      endif
      if keyword_set(photoion) then begin
         ion_rate0(j)=peo(0)*deo0+pho_rate0(j)
         ion_rate1(j)=peo(1)*deo1+pho_rate1(j)
         ion_rate2(j)=peo(2)*deo2+pho_rate2(j)
      endif
      if min(rec_all(loro,j)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(rec_all(loro,j)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      rec_rate0(j)=peo(0)*deo0
      rec_rate1(j)=peo(1)*deo1
      rec_rate2(j)=peo(2)*deo2

      if min(rec_all(loro,j+1)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(rec_all(loro,j+1)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      rec_rate0(j+1)=peo(0)*deo0
      rec_rate1(j+1)=peo(1)*deo1
      rec_rate2(j+1)=peo(2)*deo2

      if min(ion_all(loro,j-1)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(ion_all(loro,j-1)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      if not keyword_set(photoion) then begin
         ion_rate0(j-1)=peo(0)*deo0
         ion_rate1(j-1)=peo(1)*deo1
         ion_rate2(j-1)=peo(2)*deo2
      endif
      if keyword_set(photoion) then begin
         ion_rate0(j-1)=peo(0)*deo0+pho_rate0(j-1)
         ion_rate1(j-1)=peo(1)*deo1+pho_rate1(j-1)
         ion_rate2(j-1)=peo(2)*deo2+pho_rate2(j-1)
      endif

   endif

   if j eq nions-1 then begin

      if min(ion_all(loro,j-1)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(ion_all(loro,j-1)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      if not keyword_set(photoion) then begin
         ion_rate0(j-1)=peo(0)*deo0
         ion_rate1(j-1)=peo(1)*deo1
         ion_rate2(j-1)=peo(2)*deo2
      endif
      if keyword_set(photoion) then begin
         ion_rate0(j-1)=peo(0)*deo0+pho_rate0(j-1)
         ion_rate1(j-1)=peo(1)*deo1+pho_rate1(j-1)
         ion_rate2(j-1)=peo(2)*deo2+pho_rate2(j-1)
      endif
      if min(rec_all(loro,j)) gt rate_limite then $
         peo=10^spline(alog10(hhh(loro)),alog10(rec_all(loro,j)),alog10([heo0,heo1,heo2])) else peo=dblarr(3)
      rec_rate0(j)=peo(0)*deo0
      rec_rate1(j)=peo(1)*deo1
      rec_rate2(j)=peo(2)*deo2
   
   endif

endfor

t_stage_2=systime(1)
t_rate=t_stage_2-t_stage_1

; 3.12 - Calculates equilibration timescales, according to Smith & Hughes 2010

if keyword_set(check_equil) then begin

   t_test_1=systime(1)
   
   a_matrix=dblarr(nions,nions)
   
   for j=0,nions-1 do begin
   
      if j eq 1 then begin       ; Normalizes following Hughes & Helmand 1985
   
         a_matrix(j,j)=-(ion_rate0(j-1)+ion_rate0(j)+rec_rate0(j))
         a_matrix(j+1,j)=(rec_rate0(j+1)-ion_rate0(j-1))
         a_matrix(j+2:nions-1,j)=-ion_rate0(j-1)
   
      endif
   
      if j ge 2 and j le nions-2 then begin

         a_matrix(j-1,j)=ion_rate0(j-1)
         a_matrix(j,j)=-rec_rate0(j)-ion_rate0(j)
         a_matrix(j+1,j)=rec_rate0(j+1)

      endif

      if j eq nions-1 then begin

         a_matrix(j-1,j)=ion_rate0(j-1)
         a_matrix(j,j)=-rec_rate0(j)

      endif

   endfor

   a_matrix=-a_matrix(1:*,1:*)

   hes=elmhes(a_matrix,/double)
   eigenval=hqr(hes,/double)

   time_lim=4.6*max(1/eigenval)     ; Max time needed to reach 1% from equilibrium

   if real_part(time_lim) lt time then equilibrium=1.0   ; Plasma is in ionization equilibrium
   
   t_test_2=systime(1)
   t_test_eq=t_test_2-t_test_1

endif else t_test_eq=0.0

; 4 - Calculates ion fractions

; 4.1 - Out of equilibrium case

if equilibrium eq 0.0 then begin

; 4.1.1 - Defines the Runge Kutta coefficients (de facto initializing them)

   k1=dblarr(nions)
   k2=dblarr(nions)
   k3=dblarr(nions)
   k4=dblarr(nions)

; 4.1.2 - Calculates k1, k2, k3, k4 coefficients

   t_coeff_1=systime(1)
   
   for j=stage_min,stage_max do begin   ; Coefficient k1

      if j eq 0 then begin

         meo_b=old_ioneq(j)
         meo_c=old_ioneq(j+1)

         k1(j)=meo_c*rec_rate0(j+1)-meo_b*(ion_rate0(j))
         k1(j)=k1(j)*time

      endif

      if j gt 0 and j lt nions-1 then begin

         meo_a=old_ioneq(j-1)
         meo_b=old_ioneq(j)
         meo_c=old_ioneq(j+1)

         k1(j)=meo_a*ion_rate0(j-1)+meo_c*rec_rate0(j+1)-meo_b*(ion_rate0(j)+rec_rate0(j))
         k1(j)=k1(j)*time

      endif

      if j eq nions-1 then begin

         meo_a=old_ioneq(j-1)
         meo_b=old_ioneq(j)

         k1(j)=meo_a*ion_rate0(j-1)-meo_b*rec_rate0(j)
         k1(j)=k1(j)*time

      endif

   endfor

   for j=stage_min,stage_max do begin   ; Coefficient k2

      if j eq 0 then begin

         meo_b=old_ioneq(j)+0.5*k1(j)
         meo_c=old_ioneq(j+1)+0.5*k1(j+1)
   
         k2(j)=meo_c*rec_rate1(j+1)-meo_b*(ion_rate1(j))
         k2(j)=k2(j)*time

      endif

      if j gt 0 and j lt nions-1 then begin
         meo_a=old_ioneq(j-1)+0.5*k1(j-1)
         meo_b=old_ioneq(j)+0.5*k1(j)
         meo_c=old_ioneq(j+1)+0.5*k1(j+1)

         k2(j)=meo_a*ion_rate1(j-1)+meo_c*rec_rate1(j+1)-meo_b*(ion_rate1(j)+rec_rate1(j))
         k2(j)=k2(j)*time

      endif

      if j eq nions-1 then begin

         meo_a=old_ioneq(j-1)+0.5*k1(j-1)
         meo_b=old_ioneq(j)+0.5*k1(j)

         k2(j)=meo_a*ion_rate1(j-1)-meo_b*rec_rate1(j)
         k2(j)=k2(j)*time

      endif

   endfor

   for j=stage_min,stage_max do begin   ; Coefficient k3

      if j eq 0 then begin

         meo_b=old_ioneq(j)+0.5*k2(j)
         meo_c=old_ioneq(j+1)+0.5*k2(j+1)
      
         k3(j)=meo_c*rec_rate1(j+1)-meo_b*(ion_rate1(j))
         k3(j)=k3(j)*time

      endif

      if j gt 0 and j lt nions-1 then begin

         meo_a=old_ioneq(j-1)+0.5*k2(j-1)
         meo_b=old_ioneq(j)+0.5*k2(j)
         meo_c=old_ioneq(j+1)+0.5*k2(j+1)

         k3(j)=meo_a*ion_rate1(j-1)+meo_c*rec_rate1(j+1)-meo_b*(ion_rate1(j)+rec_rate1(j))
         k3(j)=k3(j)*time

      endif

      if j eq nions-1 then begin
   
         meo_a=old_ioneq(j-1)+0.5*k2(j-1)
         meo_b=old_ioneq(j)+0.5*k2(j)

         k3(j)=meo_a*ion_rate1(j-1)-meo_b*rec_rate1(j)
         k3(j)=k3(j)*time

      endif

   endfor

   for j=stage_min,stage_max do begin   ; Coefficient k4

      if j eq 0 then begin

         meo_b=old_ioneq(j)+k3(j)
         meo_c=old_ioneq(j+1)+k3(j+1)

         k4(j)=meo_c*rec_rate2(j+1)-meo_b*(ion_rate2(j))
         k4(j)=k4(j)*time
   
      endif

      if j gt 0 and j lt nions-1 then begin

         meo_a=old_ioneq(j-1)+k3(j-1)
         meo_b=old_ioneq(j)+k3(j)
         meo_c=old_ioneq(j+1)+k3(j+1)

         k4(j)=meo_a*ion_rate2(j-1)+meo_c*rec_rate2(j+1)-meo_b*(ion_rate2(j)+rec_rate2(j))
         k4(j)=k4(j)*time

      endif

      if j eq nions-1 then begin

         meo_a=old_ioneq(j-1)+k3(j-1)
         meo_b=old_ioneq(j)+k3(j)

         k4(j)=meo_a*ion_rate2(j-1)-meo_b*rec_rate2(j)
         k4(j)=k4(j)*time

      endif

; 4.1.3 - Calculates the new ion fraction for each ion j of the loop

      new_ioneq(j)=old_ioneq(j)+(k1(j)+2*k2(j)+2*k3(j)+k4(j))/6.0

   endfor

   t_coeff_2=systime(1)
   t_runge=t_coeff_2-t_coeff_1 

   if stage_min gt 0 then new_ioneq(0:stage_min-1)=0.0d0

; 4.1.4 - Tests if ion abundance is < 0

   if min(new_ioneq) lt 0.0d0 then begin

      negativo=negativo+1.0
      if negativo eq 50 then stop
      peo=where(new_ioneq lt 0) ; and abs(new_ioneq) gt 1e-10)
      if min(peo) ge 0 then begin

         quale_neg=[quale_neg,peo]
         passare=1.0
         goto,rifare

      endif

   endif

; 4.1.5 - Normalizes

   tot_ioneq=total(new_ioneq)
   ntest=tot_ioneq-1.0d0
   new_ioneq=new_ioneq/total(new_ioneq)

endif

; 4.2 - If we are in equilibrium

if equilibrium eq 1 then begin

   if i gt 0 then calcola_equil_3,zeta,ioneq_game,heo0,loro
   new_ioneq=ioneq_game
   if stage_min gt 0 then new_ioneq(0:stage_min-1)=0.0d0
   tot_ioneq=total(new_ioneq)
   ntest=tot_ioneq-1.0d0
   new_ioneq=new_ioneq/total(new_ioneq)
   t_runge=0.0

endif

; 4.3 - Decides what to make of the results and whether to continue

; 4.3.1 - Calculates the ion fraction variation from the previous step

queste=where(old_ioneq gt ioneq_min and new_ioneq gt ioneq_min)
diff=abs((new_ioneq(queste)-old_ioneq(queste))/old_ioneq(queste))

; 4.3.2 - Not good: the step is too large: (diff > delta_min)

if max(diff) ge delta_min then begin

   grosso=grosso+1.0
   passare=1.0
   goto,rifare

endif

; 4.3.3 - Not good: an ion abundance is larger than 100% (program suicides)

if max(new_ioneq) gt 1.0d0 then begin

   if not keyword_set(remote) then begin
       print,''
       print,'Stop!!!! An ion fraction is ione > 100%!'
       print,''

   endif

   stop

endif

; 4.3.4 - GOOD:

;if max(diff) lt delta_ioneq_min and heo0 le 2.0 and equilibrium eq 0.0 and negativo gt 0. then begin  ; Avoids a killer loop - EL change 12 Nov 2014
if max(diff) lt delta_ioneq_min and heo0 le 2.0 and equilibrium eq 0.0 then begin  ; Avoids a killer loop

   goto,vabene

endif

if max(diff) lt delta_min and max(diff) gt delta_ioneq_min or equilibrium eq 1 then begin

   vabene:

   taglia1=0.0

   ttt_used=[ttt_used,teo0]
   ddd_used=[ddd_used,deo0]
   vvv_used=[vvv_used,veo0]
   hhh_used=[hhh_used,heo0]
   time_travel=[time_travel,time]
   total_time=[total_time,total(time_travel)]
   n_test=[n_test,ntest]
   taglia=[taglia,taglia1]
   equil=[equil,equilibrium]

   height_in=heo0
   diff1=abs(hhh-height_in)
   lui=where(diff1 eq min(diff1))
   zeo=minmax(lui)
   loro=[zeo(0)-1-findgen(2),lui,zeo(1)+1+findgen(2)]
   loro=loro(where(loro ge 0))
   loro=loro(sort(loro))
   t_equil_1=systime(1)
   calcola_equil_3,zeta,ioneq_game,height_in,loro
   t_equil_2=systime(1)
   t_equil=t_equil_2-t_equil_1

   if i gt 0 then begin 

      new_ioneq_equil=dblarr(nions,npoints_result+1)
      new_ioneq_evol=dblarr(nions,npoints_result+1)
      new_dioneq_evol=dblarr(nions,npoints_result+1)
      new_photoion_used=dblarr(nions,npoints_result+1)
      new_recomb_used=dblarr(nions,npoints_result+1)
      new_ioniz_used=dblarr(nions,npoints_result+1)

      new_ioneq_equil(*,0:npoints_result-1)=ioneq_equil(*,*)
      new_ioneq_equil(*,npoints_result)=ioneq_game(*)
      new_ioneq_evol(*,0:npoints_result-1)=ioneq_evol(*,*)
      new_ioneq_evol(*,npoints_result)=new_ioneq(*)
      new_photoion_used(*,0:npoints_result-1)=photoion_used(*,*)
      new_photoion_used(0:nions-2,npoints_result)=pho_rate0(*)
      new_recomb_used(*,0:npoints_result-1)=recomb_used(*,*)
      new_recomb_used(*,npoints_result)=rec_rate0(*)
      new_ioniz_used(*,0:npoints_result-1)=ioniz_used(*,*)
      new_ioniz_used(*,npoints_result)=ion_rate0(*)
      new_dioneq_evol(*,0:npoints_result-1)=dioneq_evol(*,*)
      new_dioneq_evol(*,npoints_result)=new_ioneq_evol(*,npoints_result)-new_ioneq_evol(*,npoints_result-1)

   endif

   if i eq 0 then begin

      new_ioneq_equil=dblarr(nions)
      new_ioneq_evol=dblarr(nions)
      new_dioneq_evol=dblarr(nions)
      new_photoion_used=dblarr(nions)
      new_recomb_used=dblarr(nions)
      new_ioniz_used=dblarr(nions)

      new_ioneq_equil(*)=ioneq_game(*)
      new_ioneq_evol(*)=new_ioneq(*)
      new_dioneq_evol(*)=new_ioneq(*)-ioneq_start(*)
      new_photoion_used(0:nions-2)=pho_rate0(*)
      new_recomb_used(*)=rec_rate0(*)
      new_ioniz_used(*)=ion_rate0(*)

      ttt_used=ttt_used(1:*)     ; Removes first zero
      hhh_used=hhh_used(1:*)
      ddd_used=ddd_used(1:*)
      vvv_used=vvv_used(1:*)
      total_time=total_time(1:*)
      time_travel=time_travel(1:*)
      n_test=n_test(1:*)
      taglia=taglia(1:*)
      equil=equil(1:*)

   endif

   npoints_result=n_elements(ttt_used)

   ioneq_equil=new_ioneq_equil
   ioneq_evol=new_ioneq_evol
   dioneq_evol=new_dioneq_evol
   avg_ioneq0=total(ion_stage*new_ioneq)
   avg_ioneq=[avg_ioneq,avg_ioneq0]
   photoion_used=new_photoion_used
   recomb_used=new_recomb_used
   ioniz_used=new_ioniz_used
   
   old_ioneq=new_ioneq

   time_end=systime(1)
   time_giro_2=systime(1)
   time_req=time_end-time_start

   if not keyword_set(quiet) then begin
      if not keyword_set(remote) then begin
         print,string(i,'(i7)')+' - OK! H ='+strcompress(string(heo0,'(f20.7)'))+$
            '   D(H) ='+strcompress(string(delta_h,'(f20.3)'))+$
            '   Log Ne ='+strcompress(string(alog10(deo0),'(f12.3)'))+$
            '   T ='+strcompress(string(teo0,'(f20.1)'))+$
            '   D(T) ='+strcompress(string(teo2-teo0,'(f20.1)'))+$
            '   V ='+strcompress(string(veo0,'(f12.1)')) +$
            '   Ntot ='+strcompress(string(ntest,'(e12.1)'))+$
            '   Time ='+strcompress(string(time_req,'(f8.1)'))+' s'+$
            '   Neg = '+strcompress(string(negativo,'(i6)'),/remove_all)+$;'/'+$
;               strcompress(string(quale_neg,'(40i3)'))+'    '+$
            '   Gr = '+strcompress(string(grosso,'(i4)'),/remove_all)+$
            '   Sta = '+strcompress(string(stage_min,'(i4)'),/remove_all)+'-'+$
                        strcompress(string(stage_max,'(i4)'),/remove_all),+$
            '   Eq = '+strcompress(string(equilibrium,'(i4)'),/remove_all)
      endif
   endif

   ; A - Removes stuff if there are more than SIZE_RESULT_LIMIT points in the result - does NOT act on TIME_TRAVEL

   if npoints_result eq size_result_limit then begin

      meo=2*findgen(size_result_limit/2.0)
      meo=[meo,npoints_result-1]               ; Keeps last point

      ttt_used=ttt_used(meo)
      hhh_used=hhh_used(meo)
      ddd_used=ddd_used(meo)
      vvv_used=vvv_used(meo)
      total_time=total_time(meo)
      time_travel=time_travel(meo)
      avg_ioneq=avg_ioneq(meo)
      dioneq_evol=dioneq_evol(*,meo)
      ioneq_evol=ioneq_evol(*,meo)
      ioneq_equil=ioneq_equil(*,meo)
      n_test=n_test(meo)
      photoion_used=photoion_used(*,meo)
      recomb_used=recomb_used(*,meo)
      ioniz_used=ioniz_used(*,meo)

      npoints_result=n_elements(ttt_used)

      taglia2=1

   endif else taglia2=0

   if taglia2 eq 1 then taglia=[taglia(0:n_elements(taglia)-2),taglia2]

   ; B - Sets the values for next step

   i=i+1.0d0

   delta_h=delta_h*3.0d0
   teo0=teo2
   deo0=deo2
   veo0=veo2
   heo0=heo2
   passare=0.0
   equilibrium=0.0

   ; C - Sends you to the next step

   if keyword_set(no_save) then begin
      salvami=5e2*findgen(1000)
      si_salvami=where(i eq salvami)
      if si_salvami(0) ge 0 then $
         save,filename=filesave,hhh_used,ddd_used,vvv_used,ttt_used,ioneq_evol,ioneq_equil,dioneq_evol,zeta,n_test,time_travel,total_time
   endif

   if not keyword_set(no_save) then $
      save,filename=filesave,hhh_used,ddd_used,vvv_used,ttt_used,ioneq_evol,ioneq_equil,dioneq_evol,zeta,n_test,time_travel,total_time

   ; D - Writes the performance times 

      time_giro=time_giro_2-time_giro_1
      t_altro=time_giro-t_equil-t_rate-t_test_eq-t_runge

   if keyword_set(tempi) then begin

      print,'This is how I spent my time:'
      print,'Equilibrium calculation: ',t_equil/time_giro*100.0,'%',format='(a23,f8.2,a1)'
      print,'Rate calculation: ',t_rate/time_giro*100.0,'%',format='(a23,f8.2,a1)'
      print,'Equilibrium check: ',t_test_eq/time_giro*100.0,'%',format='(a23,f8.2,a1)'
      print,'Runge-Kutta calculation: ',t_runge/time_giro*100.0,'%',format='(a23,f8.2,a1)'
      print,'Calculation of all the rest: ',t_altro/time_giro*100.0,'%',format='(a23,f8.2,a1)'
      print,'Calculation re-doing: ',time_req/time_giro,' volte',format='(a23,f8.2,a6)'

   endif

   goto,rifare

endif

if max(diff) lt delta_ioneq_min and heo0 le 2.0 and equilibrium eq 0.0 then begin  ; They stopped evolving, but we are too low

   delta_h=delta_h*100.0d0
   passare=2.0
   goto,rifare

endif

if max(diff) lt delta_ioneq_min and heo0 gt 2.0 then begin  ; They stopped evolving: below 2.0 Rsun does not stop

   height_in=heo0
   diff=abs(hhh-height_in)
   lui=where(diff eq min(diff))
   zeo=minmax(lui)
   loro=[zeo(0)-1-findgen(2),lui,zeo(1)+1+findgen(2)]
   loro=loro(sort(loro))
   calcola_equil_3,zeta,ioneq_game,height_in,loro

   new_ioneq_equil=dblarr(nions,npoints_result+1)
   new_ioneq_evol=dblarr(nions,npoints_result+1)

   new_ioneq_equil(*,0:npoints_result-1)=ioneq_equil(*,*)
   new_ioneq_equil(*,npoints_result)=ioneq_game(*)
   new_ioneq_evol(*,0:npoints_result-1)=ioneq_evol(*,*)
   new_ioneq_evol(*,npoints_result)=new_ioneq(*)

   ioneq_equil=new_ioneq_equil
   ioneq_evol=new_ioneq_evol

   ttt_used=[ttt_used,teo0]
   ddd_used=[ddd_used,deo0]
   vvv_used=[vvv_used,veo0]
   hhh_used=[hhh_used,heo0]
   time_travel=[time_travel,time]
   total_time=[total_time,total(time_travel)]

   npoints_result=n_elements(ttt_used)

   goto,festa

endif

; 5 - Ends and saves results

festa:

if not keyword_set(serial_calc) then begin

   save,filename=filesave,hhh_used,ddd_used,vvv_used,ttt_used,ioneq_evol,ioneq_equil,dioneq_evol,zeta,n_test,$
                          avg_ioneq,recomb_used,photoion_used,ioniz_used,time_travel,total_time

endif

if keyword_set(serial_calc) then begin

   if istart lt nmodels-1 then begin

      results[istart].hhh_used(0:npoints_result-1)=hhh_used(*)
      results[istart].vvv_used(0:npoints_result-1)=vvv_used(*)
      results[istart].ttt_used(0:npoints_result-1)=ttt_used(*)
      results[istart].ddd_used(0:npoints_result-1)=ddd_used(*)
      results[istart].total_time(0:npoints_result-1)=total_time(*)
      results[istart].time_travel(0:npoints_result-1)=time_travel(*)
      results[istart].ioneq_evol(*,0:npoints_result-1)=ioneq_evol(*,*)
      results[istart].ioneq_equil(*,0:npoints_result-1)=ioneq_equil(*,*)
;      results[istart].avg_ioneq(0:npoints_result-1)=avg_ioneq(*)
;      results[istart].photoion_used(*,0:npoints_result-1)=photoion_used(*,*)
;      results[istart].recomb_used(*,0:npoints_result-1)=recomb_used(*,*)
;      results[istart].ioniz_used(*,0:npoints_result-1)=ioniz_used(*,*)

      rapp_tot=(float(istart)+1)/float(nmodels)*100.0
      
      if not keyword_set(remote) then $
      print,'Fatto il modello '+strcompress(string(istart+1,'(i4)'),/remove_all)+' su '+$
            strcompress(string(nmodels,'(i4)'),/remove_all)+' : '+$
            strcompress(string(rapp_tot,'(f8.3)'),/remove_all)+' %'
       
      istart=istart+1.0

      if max(taglia) gt 0 then tagli(istart)=1.0

      goto,riparti

   endif

   if istart eq nmodels-1 then begin
   
      filesave1=filesave+'_final'

      save,filename=filesave1,results

   endif

endif

if keyword_set(serial_phot) then begin

   if istart lt nmodels-1 then begin

      results[istart].hhh_used(0:npoints_result-1)=hhh_used(*)
      results[istart].vvv_used(0:npoints_result-1)=vvv_used(*)
      results[istart].ttt_used(0:npoints_result-1)=ttt_used(*)
      results[istart].ddd_used(0:npoints_result-1)=ddd_used(*)
      results[istart].total_time(0:npoints_result-1)=total_time(*)
      results[istart].time_travel(0:npoints_result-1)=time_travel(*)
      results[istart].ioneq_evol(*,0:npoints_result-1)=ioneq_evol(*,*)
      results[istart].ioneq_equil(*,0:npoints_result-1)=ioneq_equil(*,*)
;      results[istart].avg_ioneq(0:npoints_result-1)=avg_ioneq(*)
;      results[istart].photoion_used(*,0:npoints_result-1)=photoion_used(*,*)
;      results[istart].recomb_used(*,0:npoints_result-1)=recomb_used(*,*)
;      results[istart].ioniz_used(*,0:npoints_result-1)=ioniz_used(*,*)

      rapp_tot=(float(istart)+1)/float(nmodels)*100.0
      
      if not keyword_set(remote) then $
      print,'Fatto il modello '+strcompress(string(istart+1,'(i4)'),/remove_all)+' su '+$
            strcompress(string(nmodels,'(i4)'),/remove_all)+' : '+$
            strcompress(string(rapp_tot,'(f8.3)'),/remove_all)+' %'

      istart=istart+1.0

      if max(taglia) gt 0 then tagli(istart)=1.0

      goto,riparti_phot

   endif

   if istart eq nmodels-1 then begin

      filesave1=filesave+'_final'

      save,filename=filesave1,results,r_flare

   endif

endif

; 6 - That's it!

tempo_fine=systime(1)

if not keyword_set(quiet) then begin
   if not keyword_set(remote) then begin
      print,'Fine!! Ci ho messo '+$
    strcompress(string(tempo_fine-tempo_inizio,'(f10.1)'),/remove_all)+'s',format='(a30)'
   endif
endif

;stop

end
