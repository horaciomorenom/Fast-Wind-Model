pro run_wind

; Note: fast wind model: use day_photoion='26-03-2007'
;       slow wind model: use day_photoion='21-02-2002'

defsysv,'!xuvtop', '/usr/local/ssw/packages/chianti/dbase/'

wind_diagn,6,/photoion,file_model='param_fast_wind.save',/remote,/no_plot,day_photoion='21-06-2008', fileout='/home/hmorenom/SSW_Files/OgModel/temp/pred_c.save'

wind_diagn,8,/photoion,file_model='param_fast_wind.save',/remote,/no_plot,day_photoion='21-06-2008', fileout='/home/hmorenom/SSW_Files/OgModel/temp/pred_o.save'

wind_diagn,26,/photoion,file_model='param_fast_wind.save',/remote,/no_plot,day_photoion='21-06-2008', fileout='/home/hmorenom/SSW_Files/OgModel/temp/pred_fe.save'

;wind_diagn,8,/photoion,file_model='param_fast_wind.save',/remote,/no_plot,day_photoion='26-03-2007', fileout='/home/hmorenom/SSW_Files/OgModel/temp/pred_o.save'

;03-06-2009 21-06-2008 07-04-2008 12-08-2008

end

