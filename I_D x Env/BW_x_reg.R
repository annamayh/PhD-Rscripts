
### so birth weight and day they were born also affected by region .... eh?



m1=glmmTMB(Day_seq~ Sex + MotherStatus + mum_age+mum_age_sq+Reg+BirthWt+
                           (1|BirthYear)+(1|MumCode), 
                         family=gaussian, 
                         data=surv_loc_df, 
                         na.action = na.omit,
)

summary(m1)


Day=plot(ggpredict(m1, terms = c("Reg")))



m2=glmmTMB(BirthWt~ Sex + MotherStatus + mum_age+mum_age_sq+Reg+FROH+
             (1|BirthYear)+(1|MumCode), 
           family=gaussian, 
           data=surv_loc_df, 
           na.action = na.omit,
)

summary(m2)

BWt=plot(ggpredict(m2, terms = c("Reg")))
BWt

bw_day_plot=Day+BW



ggsave(BWt,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/region_birth_wt.png",
       width = 5,
       height = 5,
       bg="white")


