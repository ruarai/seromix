

e <- new.env()
load("input_data/kucharski_2018/R_datasets/FluScapeH3_data.RData", envir = e)

modelled_years <- 1968:2009

process_individual_year <- function(data) {
  data %>% 
    t() %>%
    as_tibble() %>%
    `colnames<-`(c("year_observed", "observed_titre", "strain_year", "ix_sample", "age"))
}

obs_df <- map(1:length(e$test.list), function(i) {
  map(e$test.list[[i]], process_individual_year) %>%
    bind_rows() %>% 
    mutate(ix_subject = i)
}) %>%
  bind_rows() %>%
  
  mutate(ix_subject = as.integer(ix_subject),
         ix_t_obs = match(year_observed, modelled_years),
         ix_strain = match(strain_year, modelled_years))  %>%
  
  select(ix_subject, ix_t_obs, ix_strain, age, year_observed, strain_year, observed_titre) %>%
  drop_na(observed_titre)

tar_source()



model_data <- read_fluscape_data_neuts()

obs_A <- model_data$observations %>%
  select(ix_subject, ix_strain, ix_t_obs, age, observed_titre) %>%
  mutate(across(everything(), c)) %>%
  arrange(ix_subject, ix_strain)


obs_B <- obs_df %>%
  select(ix_subject, ix_strain, ix_t_obs, age, observed_titre) %>%
  arrange(ix_subject, ix_strain)

# If no differences, our data reading code is good
waldo::compare(obs_A, obs_B)



load("input_data/kucharski_2018/datasets/spline_fn.RData")

scalemap<-function(xx,inf_years){
  if(max(inf_years)>2012){stop("need infection range to be inside antigenic map")}
  map.range <- c(1968:2012)
  alen <- c(333.83,370.28); alenA <- (alen[2]-alen[1])/(max(xx)-min(xx)); alenB <- alen[1]-alenA*min(xx)
  s1 <- alenA*xx+alenB-alen[1]; s1*length(inf_years)/length(map.range) +alen[1]
}

outputdmatrix.fromcoord <- function(thetasigma,inf_years,anti.map.in,linearD=F){ #anti.map.in can be vector or matrix - rows give inf_years, columns give location
  
  # Check if map is 1D or 2D
  #if(length(anti.map.in)==length(inf_years)){
  if(linearD==F){
    # Exponential decay function
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
      (dmatrix=sapply(anti.map.in,function(x){exp(-thetasigma*abs(anti.map.in-x))}))
    }else{ # If spline function defined, calculate directly from input
      (dmatrix=apply(anti.map.in,1,function(x){exp(-thetasigma*sqrt(
        colSums(apply(anti.map.in,1,function(y){(y-x)^2}))
      ))}))
    }
  }else{
    # Linear decay function
    if(is.null(dim(anti.map.in))){ # check if input map is one or 2 dimensions
      (dmatrix=sapply(anti.map.in,function(x){y=abs(anti.map.in-x); y   })) 
      
    }else{ # If spline function defined, calculate directly from input
      (dmatrix=apply(anti.map.in,1,function(x){y=sqrt(
        
        colSums(apply(anti.map.in,1,function(y){(y-x)^2}))
        
      ); y # 1-1*thetasigma* //  y[y<0]=0; HAVE REMOVED BASE
      }))
    }
  }
}

xx=scalemap(e$inf_years,e$inf_years)
yy=predict(am.spl,xx)$y
antigenic.map.in=cbind(xx,yy)


dmatrix = 1 - outputdmatrix.fromcoord(-1,e$inf_years,antigenic.map.in,linearD=TRUE)

write_rds(dmatrix, "input_data/kucharski_2018/antigenic_distances_kucharski_fluscape.rds")

image(dmatrix)
image(1 - model_data$antigenic_distances)



