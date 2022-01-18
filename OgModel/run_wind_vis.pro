pro run_wind_vis

; Note: fast wind model: use day_photoion='26-03-2007'
;       slow wind model: use day_photoion='21-02-2002'

defsysv,'!xuvtop', '/usr/local/ssw/packages/chianti/dbase/'

wind_diagn,6,/photoion,file_model='paramvis_fast_wind.save',/remote,/no_plot,day_photoion='26-03-2007', fileout='/home/hmorenom/SSW_Files/OgModel/temp_vis/pred_c.save'

wind_diagn,8,/photoion,file_model='paramvis_fast_wind.save',/remote,/no_plot,day_photoion='26-03-2007', fileout='/home/hmorenom/SSW_Files/OgModel/temp_vis/pred_o.save'

wind_diagn,26,/photoion,file_model='paramvis_fast_wind.save',/remote,/no_plot,day_photoion='26-03-2007', fileout='/home/hmorenom/SSW_Files/OgModel/temp_vis/pred_fe.save'

;wind_diagn,8,/photoion,file_model='param_fast_wind.save',/remote,/no_plot,day_photoion='26-03-2007', fileout='/home/hmorenom/SSW_Files/OgModel/temp/pred_o.save'

end

