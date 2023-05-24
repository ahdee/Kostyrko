
# functions for survival analysis 

# create function to extract all mut.this mutants and all WT. 
# mut df must have the basic colums here 
# sample gene effect 
getmut.id = function ( mutation, mut.df){
  # grab all non silent mutations 
  
  mutation.pid <-  unique ( mut.df[!mut.df$effect=="Silent" & mut.df$gene==mutation, ]$sample )
  # for wild type we use only samples that do not have the mutation of interest.  We do not care of the mutation is.
  # 1. silent 2. even if there are x samples per patient, if that patient have at least once it will not be included in the WT.
  # the reason is because we rather err on the side of caution in case the algo did not pick up one of the sample for whatever reason. also
  # logistically we can only calculate once, censor once ie death so there is not choice 
  # to do this proper we must get all samples that have this mutation and then negatively subset it.  Do not do != mutation directly 
  # since it might include some samples from the same patient that does not contained the mutation and will thusly be put into the wt bin
  # which is not correct. 
  
  all.mut = unique ( xena.mut[xena.mut$gene==mutation, ]$sample )
  # here we cannot simply subset because samples are repeated that is for example 
  # x may have also 
  wt.pid <- unique ( xena.mut[ !xena.mut$sample %in% all.mut , ]$sample )
  
  mutation.pid = substr(mutation.pid,1,12)
  wt.pid = substr(wt.pid,1,12)
  
  return ( list (mutation.pid = mutation.pid, wt.pid=wt.pid ))
}



do_surv.hr <- function (  x, title, censor,time,xcut = 4000,breakt = 500   ){
  
  eq2 = paste( eq, collapse = "+") 
  eq2 = paste ( "s_objM ~", eq2 )
  FM=as.formula(eq2)  
  
  # create survival object and run COX multivariate
  s_objM<-Surv(
    as.numeric(x[ , time] ), 
    as.numeric(x[ , censor])
  )
  
  
  all.coxM <- coxph(FM, data =  x, method="efron")
  
  
  
  # extract and make pretty table summary 
  df = summary ( all.coxM  )
  df2 = round ( df$coefficients[, c(2,5)], digits=4)
  df2 = cbind ( df2, df$conf.int[, c(3:4)])
  colnames ( df2 ) = c("HR","p-value", "CI.95.low", "CI.95.high")
  df2 = round(df2, digits=3)
  summ = t( data.frame ( 
    concordance = df$concordance[1],
    logrank = df$sctest[3],
    wald = df$waldtest[3],
    likelihood = df$logtest[3]
  )
  )
  colnames ( summ) = "summary"
  summ = round ( summ,digits=3)
  
  # plot this 
  p = ggsurvplot(survminer::surv_fit(
    s_objM ~ x$state, data=x), 
    legend.labs = c("Normal", "High")    # change legend labels
    , risk.table = TRUE, palette = c(normcol, highcol)
    ,  xlim = c(0,xcut), # present narrower X axis, but not affect
    # survival estimates.
    break.time.by = breakt
  ) 
  
  p$plot = p$plot +xlab ( "Days") 
  p$table$labels$x = ""
  
  top = p$plot  + ggtitle ( title ) + theme(plot.title = element_text(size = 12))
  bottom = p$table #+  labs(caption = "cat(df2)") 
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = c("#cccccc","#dddddd","#eeeeee","#eeeeee"), col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="grey", fontface=4L)),
    rowhead=list(fg_params=list(col="black", fontface=3L)))
  
  p3 = grid.arrange ( top,bottom, nrow=2 )
  p4 =  tableGrob(df2, theme = tt3) 
  th <- sum(p4$heights)
  final = grid.arrange(p3,p4,nrow=2, heights = unit.c(unit(1, "null"), th))
  
  return ( final )
}


draw.surv <- function ( x, title="" ,xcut = 4000,breakt = 500   ){
  # x = ss.mut, paste ( mutation, "mutation", "gene:", test.gene)
  pfi = do_surv.hr( x, paste ( title, "Progression Free"), censor = "PFI", time = "PFI.time",xcut = 4000,breakt = breakt )
  dfi = do_surv.hr( x, paste ( title, "Disease Free"), censor = "DFI", time = "DFI.time",xcut = 4000,breakt = breakt )
  osi = do_surv.hr( x, paste ( title, "Overall Survival"), censor = "OS", time = "OS.time",xcut = 4000,breakt = breakt  )
  dss = do_surv.hr( x, paste ( title, "Disease Specific Survival"), censor = "DSS", time = "DSS.time",xcut = 4000,breakt = breakt  )
  return ( list ( pfi=pfi, dfi=dfi, osi=osi, dss=dss ))
}

# function to take coxh object and put that as a dataframe 
org.cox <- function ( df, test.gene, mutation){
  
  df = summary ( df )
  df2 = round ( df$coefficients[, c(2,5)], digits=4)
  df2 = cbind ( df2, df$conf.int[, c(3:4)])
  colnames ( df2 ) = c("HR","p-value", "CI.95.low", "CI.95.high")
  concordance = df$concordance
  logrank = df$sctest
  wald = df$waldtest
  likelihood = df$logtest
  
  ## plot cox model 
  # we need two things. 
  # a label mattrix this is the left columns giving data like name, hr, p.value
  
  # create matrix 
  # 1st lets convert to dataframe for easier extraction 
  
  df2 = data.frame ( df2, stringsAsFactors = F)
  # rename row.names so that it is more clear. 
  row.names ( df2 ) = c(
    paste ( test.gene, "state" )
    ,paste ( mutation )
    , "age", "gender", "Stage (III & IV)"
    , paste0("Interaction ", "state", ":", mutation)
  )
  # create the matrix 
  labeldf <- cbind( c( "covariate", rownames(df2) )
                    ,c("Hazard.ratio", df2$HR)
                    ,c("p.value",df2$p.value)
                    
  )
  # add summary data like ovearll pvalue and concordance
  
  labeldf <- rbind ( labeldf
                     , c("Concordance", NA,  round ( as.numeric ( df$concordance[1]), digits=3) )
                     , c("Overall Wald p.value", NA, round ( as.numeric ( wald[3] ), digits = 12)  )
                     , c("Score logrank p.value", NA, round ( as.numeric ( logrank[3] ), digits = 12)  )
  )
  
  
  # next we need to creat the confidence level for HR 
  # however we add NA to rows that are consider "summary" like title and 3 ( eg concordance) more at the end to fill in the summary data.
  meta2 = data.frame ( 
    mean = c( NA, df2$HR, NA,NA,NA )
    ,lower=c(NA, df2$CI.95.low, NA,NA,NA)
    ,upper=c(NA, df2$CI.95.high, NA,NA,NA)
  )
  
  ### 
  
  own <- fpTxtGp(
    label = gpar(fontfamily = "", cex=.9), 
    summary= gpar(fontfamily = "", cex=1.2), 
    ticks = gpar(fontfamily = "", cex=.65), 
    xlab=gpar(cex = 1.0)
  )
  
  fp= grid.grabExpr( print ( forestplot(labeldf, 
                                        meta2,
                                        new_page = TRUE,
                                        # tell it what rows are summary.  
                                        # set the rep ( FALSE, x) where x == number of rows for data, summary should be TRUE thus the name is.summary, haha
                                        is.summary=c(TRUE,rep(FALSE,nrow (df2)),TRUE,TRUE,TRUE),
                                        
                                        lineheight = "auto", 
                                        
                                        lwd.ci = 3, # sets the width of the line 
                                        # sets the box size, by each row 
                                        boxsize= rep(.3, 6),
                                        col=fpColors(box="grey",line="black", summary="#3268e5"),
                                        # define font size, see above gpar
                                        txt_gp =  own, zero = 1, # set to one and not 0 because 1 == no change in survival
                                        xlog = FALSE,
                                        #x.ticks= seq(-1.5,maxx,0.5),
                                        line.margin = .21, 
                                        xlab="(< 1) good prognostic  HR  (> 1) poor prognostic ", 
                                        
                                        ### creates lines 
                                        
                                        hrzl_lines=list( 
                                          "3" = gpar(lwd=1, lineend="butt", columns=c(1:4), col="grey") ,
                                          "4" = gpar(lwd=1, lineend="butt", columns=c(1:4), col="grey"),  
                                          "5" = gpar(lwd=1, lineend="butt", columns=c(1:4), col="grey"),
                                          "6" = gpar(lwd=1, lineend="butt", columns=c(1:4), col="grey"),
                                          "7" = gpar(lwd=1, lineend="butt", columns=c(1:4), col="grey"),
                                          "8" = gpar(lwd=3, lineend="butt", columns=c(1:4), col="#3268e5")
                                        )
  ) # end forest plot
  )) # for capturing
  #####
  
  
  
  # end plot cox model 
  
  
  
  return ( list ( cox=df2, concordance=concordance, logrank=logrank, wald=wald, likelihood=likelihood, forest = fp) )
}

# here is a wrapper for Liu that will do everything at one go 

wrapper.surv = function(test.gene, mutation, qn=.75, counts=tumor, mut.table=xena.mut, sclin=sclin, censor="OS", censor.time="OS.time", zheng=0 , groupname="" 
                        ,xcut = 4000,breakt = 1000
){
  
  
  
  # make sure tha sclin sample has only . and not - 
  sclin$sample = gsub ( "-",".", sclin$sample)
  ## tcga code, the portion 01 and the anlyte, the analyte can be safely removed here. However it should be removed in the counts if it exists.  The sclin does not contain this info!  
  # remove analyte info if it exits
  cname = colnames ( counts )
  cname = gsub ( "[A-Z]$", "", cname  )
  colnames ( counts ) = cname 
  # qn = .75 # 3 is exactly 50% 
  #test.gene = "IL6ST"
  #mutation = mut.this
  # made changes to i to make sure to retain the all important .01 or .11 to ensure we don't end up getting normal samples; however with that said its fine as well because below we remove dups and regardless of what the sample name says so long as the original df only had tumor samples it will be fine since the data are duplicated.  # zhengASC <- function (genes, count, p=.75) 
  
  if ( zheng == 0 ){
    i <- qGet(test.gene,qn, counts)
  } else if ( zheng == 1) {
    i <- zhengASC(genes = test.gene, count=counts, p=qn)
    colnames  ( i ) = c ( "count","state","sample")
    i$sample2 = stringr::str_match(i$sample, "(.*)\\.")[ , 2]
  }
  ## get the list of sample ids for all muation ( eg mut.this) active mutation and all wildtype 
  group.id = getmut.id ( mutation, mut.table)
  
  # merge clinical data 
  ss.all = merge ( i, sclin, by.x="sample", by.y="sample") # merge by sample. This is important because it assumes that filtering was already been done upstream and that some samples may have dups. 
  ss.all = ss.all[!duplicated ( ss.all$sample), ] # there should be no dups but we do it here just to be safe! 
  colnames(ss.all)[grepl("age_at_initial_pathologic_diagnosis", names(ss.all))]= "age.diag"
  
  ss.all$genome = ifelse ( ss.all$X_PATIENT %in% c(group.id$mutation.pid ), mutation, "wt")
  ss.all = ss.all[ss.all$X_PATIENT %in% c(group.id$mutation.pid, group.id$wt.pid  ),  ]
  ss.all$genome = factor ( ss.all$genome, levels=c("wt",mutation ))
  
  # let dummy code the tumor state
  ss.all$stage = ss.all$ajcc_pathologic_tumor_stage
  ss.all$stage = as.character ( sapply ( ss.all$stage, function(x) {
    st = NA
    if ( length ( x[grepl ( "Stage IV", x)] ) != 0  ){
      st = "4,3"
    }else if (length ( x[grepl ( "Stage III", x)] ) != 0 ) {
      st = "4,3"
    }else if (length ( x[grepl ( "Stage II", x)] ) != 0 ) {
      
      st = "2,1"
    }else if (length ( x[grepl ( "Stage I", x)] ) != 0 ){
      st = "2,1"}
    return (st)
  }) )
  
  # lets make sure that age is numeric 
  ss.all$age.diag = as.numeric ( ss.all$age.diag)
  # df above is ready to do analysis
  
  # startify by wt and mutation although you can also subset with ss.all; the two should be equal! 
  ss.mut = ss.all [ ss.all$X_PATIENT %in% group.id$mutation.pid, ]
  ss.wt = ss.all [ ss.all$X_PATIENT %in% group.id$wt.pid, ]
  
  # we need to first compare the expression of wildtype to mutation 
  i2 = ss.all
  i2$state = ifelse ( i2$state == 0, "Normal", "High")
  i2$state = factor ( i2$state, levels = c("Normal","High"))
  
  
  bee = ggplot( i2 , aes(y=as.numeric(count), x=as.character ( genome) ))+
    geom_quasirandom(aes(fill = as.character ( state ), color=as.character ( state) ), size = 10, shape = 21, alpha = .6) +
    theme_bw() +
    ylab(paste ( test.gene, "log2 ( CPM) ") ) +
    xlab("") +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_text(size= 12 ),
          axis.text.x = element_text(angle = 0, size= 12 ),
          axis.title.x = element_text(size=22),
          
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + scale_fill_manual(values = c(highcol, normcol ) ) + scale_colour_manual(values = c("grey", "grey" ) ) +
    ggtitle(paste( mutation, "vs.", "WT"))
  
  ###################### run cox hazard full model multivariate 
  
  
  
  s_obj<-Surv(
    as.numeric(ss.all[ , censor.time] ), 
    as.numeric(ss.all[ , censor])
  )
  
  all.cox <- coxph(s_obj ~ state*genome + age.diag + gender + stage, data =  ss.all, method="efron")
  # the model below is the same thing! 
  #all.cox <- coxph(s_obj ~ state + genome + age.diag + gender + stage + state:genome, data =  ss.all, method="efron")
  
  
  # test for hazard 
  test.ph <- cox.zph(all.cox)
  
  # plot full model 
  if ( zheng == 0 ){
    all.cox.plot = org.cox ( df=all.cox , test.gene, mutation )
  }else { all.cox.plot = NULL }
  ###
  
  
  ## plot wt vs mutation 
  ss.all$group = paste ( ss.all$state, ":", ss.all$genome )
  aa.hr.plot = ggsurvplot(survminer::surv_fit(
    s_obj ~ ss.all$group, data=ss.all) 
    , legend.labs = c(paste ( mutation, "Normal" ), "WT Normal", paste ( mutation, "High" ), 
                      "WT High"
    )    # change legend labels
    , risk.table = TRUE
    , palette = c("#ffd9b3", "#d9b3ff" , "#ff8c1a",  "#a300cc")
    ,  xlim = c(0,4000), # present narrower X axis, but not affect
    # survival estimates.
    break.time.by = breakt,     # break X axis in time intervals by 500.
    pval = FALSE,
  ) 
  
  aa.hr.plot = grid.arrange ( aa.hr.plot$plot, aa.hr.plot$table, nrow=2)
  
  ### END all hr plotting 
  if ( zheng == 1){
    test.gene = groupname 
  }
  
  
  survival.mutation = draw.surv ( ss.mut, paste ( mutation, "mutation", "gene:", test.gene ),xcut = 4000,breakt = breakt )
  plot.mut = grid.arrange(survival.mutation$osi, survival.mutation$dss, survival.mutation$dfi, survival.mutation$pfi, nrow=2)
  
  survival.wt = draw.surv ( ss.wt, paste ( mutation, "WT" , "gene:", test.gene ),xcut = 4000,breakt = breakt  )
  plot.wt = grid.arrange(survival.wt$osi,survival.wt$dss, survival.wt$dfi, survival.wt$pfi, nrow=2)
  
  survival.all = draw.surv ( x=ss.all, title=paste (  "ALL" , ":", test.gene)  ,xcut = 4000,breakt = breakt ) 
  plot.all = grid.arrange(survival.all$osi,survival.all$dss, survival.all$dfi, survival.all$pfi, nrow=2)
  
  
  
  # make a final ridge plot to summarize the distribution of ALL data and combine this with the bee plot above
  q <- quantile ( as.numeric(i2$count), probs = qn)
  q50 <- as.numeric(q[1]) 
  
  jplot = ggplot(i2, aes(x=count, y=genome, fill = genome)) +
    geom_density_ridges(
      scale = 0.9,
      #aes(point_color = "grey", point_fill = state, point_shape = "19"),
      alpha = .52, point_alpha = .1, jittered_points = TRUE, point_shape = 19
    ) +
    xlab(paste ( test.gene, "log2 ( CPM) ")) + theme_bw() +
    theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_text(size= 12 ),
          axis.text.x = element_text(angle = 90, size= 12 ),
          axis.title.x = element_text(size=22),
          
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + scale_fill_manual(values = c(normcol, highcol ) )  +
    geom_vline(xintercept = q50, linetype="dotted",
               color = "black", size=1.5) +
    geom_text(aes(x=q50, label = paste ( "\n", qn), y=1), colour="black", angle=90, text=element_text(size=15))
  
  
  
  joybee = grid.arrange( jplot, bee, nrow=1)
  
  
  return ( list ( plot.mut=plot.mut, 
                  plot.wt=plot.wt, 
                  joybee = joybee, 
                  all.forest=all.cox.plot$forest,
                  all.kp = aa.hr.plot, 
                  all.cox = all.cox, test.ph=test.ph,
                  survival.all = plot.all, 
                  df = ss.all,
                  extplots = list ( joybee = joybee, bee=bee)
  ))
  
}

# this function will create a dataframe dummy coding 0 or 1 based of gene expression. 
# the qBin == 3 is 75 percentile 
qGet <- function (gene, qBin, rnadf){
  
  values <- (rnadf[row.names(rnadf) %in% gene,])
  
  # we can either use the default quantile .25,.5,.75 so for example qBin == 4 is 75% percentile
  if ( qBin < 1 ){
    q <- quantile ( as.numeric(values), probs = qBin)
    q50 <- as.numeric(q[1]) 
  }else{
    q <- quantile(as.numeric ( values ) )
    q50 <- as.numeric(q[qBin]) 
  }
  
  
  
  state <- c()
  
  # 0 = low expression group
  for (s in values){
    
    if (s < q50){
      state <- append(state,0)
    }
    else{
      state <- append(state,1)
    }
    
  }
  
  final <- data.frame( count=as.numeric ( values) ,state,sample=names(values))
  final$sample2 <- gsub("\\.","-",substr( final$sample,1,12))
  
  # get unique only ; some sample have multiple entries this could be due to mutliple sampling? 
  # we will not be usinging these samples since it could be ambigious 
  
  usample = final[duplicated ( final$sample), ]$sample
  
  final = final [ ! final$sample %in% usample , ]
  
  row.names(final)<- final$sample
  
  return(final)
  
  
  
}


## simple ones 
## 
do_cox = function ( x, title , censor = "OS" , time = "OS.time") {
  eq2 = paste( eq, collapse = "+") 
  eq2 = paste ( "s_objM ~", eq2 )
  FM=as.formula(eq2)  
  
  # create survival object and run COX multivariate
  s_objM<-Surv(
    as.numeric(x[ , time] ), 
    as.numeric(x[ , censor])
  )
  
  
  all.coxM <- coxph(FM, data =  x, method="efron")
  
  # all.coxM <- coxph(s_objM ~ state + stage + age.diag + gender , data =  x, method="efron")
  
  
  # extract and make pretty table summary 
  df = summary ( all.coxM  )
  df2 = round ( df$coefficients[, c(2,5)], digits=4)
  df2 = cbind ( df2, df$conf.int[, c(3:4)])
  colnames ( df2 ) = c("HR","p-value", "CI.95.low", "CI.95.high")
  df2 = round(df2, digits=3)
  summ = t( data.frame ( 
    concordance = df$concordance[1],
    logrank = df$sctest[3],
    wald = df$waldtest[3],
    likelihood = df$logtest[3]
  )
  )
  colnames ( summ) = "summary"
  summ = round ( summ,digits=3)
  df2 = data.frame ( df2 )
  
  df2$nszie = df$n
  df3 = df2["state", , drop=F]
  df3$censor = censor 
  
  
  return ( df3 )
  
}



dosimple = function ( test.gene, mutation, qn, counts, mut.table, sclin   ){
  mut = mutation
  # make sure tha sclin sample has only . and not - 
  sclin$sample = gsub ( "-",".", sclin$sample)
  ## tcga code, the portion 01 and the anlyte, the analyte can be safely removed here. However it should be removed in the counts if it exists.  The sclin does not contain this info!  
  # remove analyte info if it exits
  cname = colnames ( counts )
  cname = gsub ( "[A-Z]$", "", cname  )
  colnames ( counts ) = cname 
  # qn = .75 # 3 is exactly 50% 
  #test.gene = "IL6ST"
  #mutation = mut.this
  # made changes to i to make sure to retain the all important .01 or .11 to ensure we don't end up getting normal samples; however with that said its fine as well because below we remove dups and regardless of what the sample name says so long as the original df only had tumor samples it will be fine since the data are duplicated.  
  i <- qGet(test.gene,qn, counts)
  
  
  
  
  
  ## get the list of sample ids for all muation ( eg mut.this) active mutation and all wildtype 
  group.id = getmut.id ( mutation, mut.table)
  
  # merge clinical data 
  ss.all = merge ( i, sclin, by.x="sample", by.y="sample") # merge by sample. This is important because it assumes that filtering was already been done upstream and that some samples may have dups. 
  ss.all = ss.all[!duplicated ( ss.all$sample), ] # there should be no dups but we do it here just to be safe! 
  colnames(ss.all)[grepl("age_at_initial_pathologic_diagnosis", names(ss.all))]= "age.diag"
  
  ss.all$genome = ifelse ( ss.all$X_PATIENT %in% c(group.id$mutation.pid ), mutation, "wt")
  ss.all = ss.all[ss.all$X_PATIENT %in% c(group.id$mutation.pid, group.id$wt.pid  ),  ]
  ss.all$genome = factor ( ss.all$genome, levels=c("wt",mutation ))
  
  # let dummy code the tumor state
  ss.all$stage = ss.all$ajcc_pathologic_tumor_stage
  ss.all$stage = as.character ( sapply ( ss.all$stage, function(x) {
    st = NA
    if ( length ( x[grepl ( "Stage IV", x)] ) != 0  ){
      st = "4,3"
    }else if (length ( x[grepl ( "Stage III", x)] ) != 0 ) {
      st = "4,3"
    }else if (length ( x[grepl ( "Stage II", x)] ) != 0 ) {
      
      st = "2,1"
    }else if (length ( x[grepl ( "Stage I", x)] ) != 0 ){
      st = "2,1"}
    return (st)
  }) )
  
  # lets make sure that age is numeric 
  ss.all$age.diag = as.numeric ( ss.all$age.diag)
  # df above is ready to do analysis
  
  # startify by wt and mutation although you can also subset with ss.all; the two should be equal! 
  ss.mut = ss.all [ ss.all$X_PATIENT %in% group.id$mutation.pid, ]
  ss.wt = ss.all [ ss.all$X_PATIENT %in% group.id$wt.pid, ]
  
  # 
  combine = data.frame ()
  for ( type in c("mut","wt","all")){
    
    if ( type == "mut"){
      x=ss.mut 
      type = mut 
    }  else if ( type=="wt"){
      x=ss.wt 
    } else {
      x = ss.all 
    }
    
    pfi = do_cox( x, paste ( title, "Progression Free"), censor = "PFI", time = "PFI.time")
    dfi = do_cox( x, paste ( title, "Disease Free"), censor = "DFI", time = "DFI.time")
    osi = do_cox( x, paste ( title, "Overall Survival"), censor = "OS", time = "OS.time" )
    dss = do_cox( x, paste ( title, "Overall Specific Survival"), censor = "DSS", time = "DSS.time" )
    
    df = rbind ( pfi, dfi, osi, dss )
    df$gene = gene.this 
    df$type = type 
    combine = rbind ( combine, df )
  }
  
  return ( combine )
  
  
}



zhengASC <- function (genes, count, p=.75){
  
  test1 <- count[row.names(count) %in% genes,] # get all the values for each genes
  sum_test <- apply(test1,2,sum)
  mean_test <- mean(sum_test)
  
  
  q <- quantile(as.numeric ( sum_test), prob = p )
  q50 <- as.numeric(q[1]) 
  mean_test <- q50
  # one can also use quantile using different quartiles, eg 50%
  sum_test <- as.vector(sum_test)
  state <- c()
  
  # 0 = low expression group
  for (s in sum_test){
    
    if (s < mean_test){
      state <- append(state,0)
    }
    else{
      state <- append(state,1)
    }
  }
  
  final <- data.frame(sum_test,state)
  row.names(final) <- colnames(test1)
  
  
  final$barcode <- row.names ( final )
  
  return(final )
}