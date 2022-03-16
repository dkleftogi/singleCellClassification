#BEGIN COPYRIGHT NOTICE

#utilityFuncs -- (c) 2022 Dimitrios Kleftogiannis -- UiB and CCBIO
#Copyright 2022 University of Bergen (UiB), Norway.
#This Program is free software licensed under the MIT License.
#You may only use the source code in this repository in compliance with the license provided in this repository. 
#For more details, please refer to the file named "LICENSE.md".
#This Program is distributed as a service to the research community and is experimental in nature and may have hazardous properties. 
#The Program is distributed WITHOUT ANY WARRANTY, express or implied. In particular all warranties as to SATISFACTORY QUALITY or FITNESS FOR A PARTICULAR PURPOSE are excluded. 
#Published reports of research using this code (or a modified version) should cite the relevant article of this code
#Comments and bug reports are welcome.
#Email to dimitrios.kleftogiannis@uib.no
#I would also appreciate hearing about how you used this code, improvements that you have made to it.
#You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

#END COPYRIGHT NOTICE

#UTILITY
#Different utility functions required by generateFeatures code.

#function to normalise outliers based on quantiles, defautl values are 0.02 and 0.98 but change if needed
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(0.02, 0.98), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

#find appropriate cutoffs for DREMI estimation, adopted from the original DREMI implementation in MATLAB (c)
find_data_cutoffs <- function(m,cutoff,num_slices){
  
  X <- m[,1]
  Y <- m[,2]
  Xmin <- min(X,na.rm = TRUE)
  Xmax <- max(X, na.rm = TRUE)
  Ymin <- min(Y,na.rm = TRUE)
  Ymax <- max(Y,na.rm=TRUE)
  Z <- my_kde2D(m,n=num_slices+1,limits=c(Xmin,Xmax,Ymin,Ymax))
  bandwidth <- Z$bandwidth
  density <- Z$density
  Grid_X <- Z$X
  Grid_Y <- Z$Y
  #round small values of density
  density[density<0.00001] <- 0
  row=0;
  total_cells_so_far = 0; 
  while(total_cells_so_far<cutoff & row < (num_slices-1)){
    row = row+1;
    total_cells_so_far = length(find(Y>Grid_Y[num_slices-row,1]));
  }
  maxy = Grid_Y[num_slices-row,1];
  total_cells_so_far = 0;
  start_row = 0;
  while(total_cells_so_far<cutoff){
    start_row <- start_row + 1;
    total_cells_so_far <- length(find(Y<Grid_Y[start_row,1]));
  }
  miny <- Grid_Y[start_row,1];
  
  total_cells_so_far = 0;
  col=0;
  breakCounter <- 1
  while(total_cells_so_far<cutoff & col < (num_slices/2-1)){
    col <- col+1;
    total_cells_so_far <- length(find(X>Grid_X[1,num_slices-col]));
  }
  maxx = Grid_X[1,num_slices-col];
  total_cells_so_far = 0;
  start_col=0;
  while(total_cells_so_far<cutoff & start_col <= (num_slices/2)){
    start_col = start_col+1;
    total_cells_so_far = length(find(X<Grid_X[1,start_col]));
  }
  minx = Grid_X[1,start_col];
  out <- list(minx,miny,maxx,maxy)
  return(out)
}

#use this version of kde - based on the original DREMI publication (c)
my_kde2D <- function (data, n = 2^8, limits = NULL) {
  ndhist <- function(data, M) {
    bins <- zeros(nrow(data), ncol(data))
    for (i in 1:ncol(data)) {
      bins[, i] <- histc(data[, i], seq(0, 1, 1/M))$bin
      bins[, i] <- apply(as.matrix(bins[, i]), 1, function(y) min(y, 
                                                                  M))
    }
    #i added this part of code to fix the bug 
    a <- which(bins[,1]>0 & bins[,2]>0)
    bins <- bins[a,]
    out <- accumarray(bins[all(bins > 0), ], rep((1/nrow(data)), 
                                                 nrow(data)), M * (ones(1, ncol(data))))
    return(out)
  }
  dct2d <- function(data) {
    dct1d <- function(x) {
      x <- rbind(x[seq(1, nrow(x), 2), ], x[seq(nrow(x), 
                                                2, -2), ])
      out <- Re(weights * mvfft(x))
      return(out)
    }
    w <- as.matrix(c(1, 2 * (exp(-(0+1i) * (1:(nrow(data) - 
                                                 1)) * pi/(2 * nrow(data))))))
    weights <- w[, ones(1, ncol(data))]
    data <- t(dct1d(t(dct1d(data))))
    return(data)
  }
  idct2d <- function(data) {
    idct1d <- function(x) {
      y <- Re(ifft_mat(weights * x))
      out <- zeros(nrow(data), ncol(data))
      out[seq(1, nrow(data), 2), ] <- y[1:(nrow(data)/2), 
      ]
      out[seq(2, nrow(data), 2), ] <- y[seq(nrow(data), 
                                            nrow(data)/2 + 1, -1), ]
      return(out)
    }
    ifft_mat <- function(x) {
      x <- apply(x, 2, function(y) ifft(y))
      return(x)
    }
    w <- as.matrix(exp((0+1i) * (0:(nrow(data) - 1)) * pi/(2 * 
                                                             nrow(data))))
    weights <- w[, ones(1, ncol(data))]
    out <- idct1d(t(idct1d(data)))
    return(out)
  }
  K <- function(s) {
    if (s == 0) {
      out <- (-1)^s/sqrt(2 * pi)
    }
    else {
      out <- (-1)^s * prod(seq(from = 1, to = 2 * s - 
                                 1, by = 2))/sqrt(2 * pi)
    }
    return(out)
  }
  psi <- function(s, time) {
    w <- c((exp(-I * pi^2 * time))) * c(cbind(1, 0.5 * ones(1, 
                                                            length(I) - 1)))
    wx <- t(w * (I^s[1]))
    wy <- t(w * (I^s[2]))
    out <- (-1)^sum(s) * (wy %*% A2 %*% t(wx)) * pi^(2 * 
                                                       sum(s))
    return(out)
  }
  func <- function(s, t) {
    if (sum(s) <= 4) {
      sum_func <- func(c(s[1] + 1, s[2]), t) + func(c(s[1], 
                                                      s[2] + 1), t)
      const <- (1 + 1/2^(sum(s) + 1))/3
      time <- c((-2 * const * K(s[1]) * K(s[2])/N/sum_func)^(1/(2 + 
                                                                  sum(s))))
      out <- psi(s, time)
    }
    else {
      out <- psi(s, t)
    }
    return(out)
  }
  evolve <- function(t) {
    sum_func <- func(c(0, 2), t) + func(c(2, 0), t) + 2 * 
      func(c(1, 1), t)
    time <- (2 * pi * N * sum_func)^(-1/3)
    out <- (t - time)/time
    return(c(out, time))
  }
  subtract_evolve <- function(t) {
    return(t - evolve(t))
  }
  N <- nrow(data)
  n <- 2^ceiling(log2(n))
  if (is.null(limits)) {
    max <- c(max(data[, 1]), max(data[, 2]))
    min <- c(min(data[, 1]), min(data[, 2]))
    range <- max - min
    max_xy <- max + range/4
    min_xy <- min - range/4
  }
  else {
    max_xy <- c(limits[2], limits[4])
    min_xy <- c(limits[1], limits[3])
  }
  scaling <- max_xy - min_xy
  data_trans <- (data - repmat(min_xy, N, 1))/repmat(scaling, 
                                                     N, 1)
  data_init <- ndhist(data_trans, n)
  data_DCT <- dct2d(data_init)
  I <- (0:(n - 1))^2
  A2 <- data_DCT^2
  t_star <- fzero(subtract_evolve, c(0, 0.1))[[1]]
  p_02 <- func(c(0, 2), t_star)
  p_20 <- func(c(2, 0), t_star)
  p_11 <- func(c(1, 1), t_star)
  t_y <- (p_02^(3/4)/(4 * pi * N * p_20^(3/4) * (p_11 + sqrt(p_20 * 
                                                               p_02))))^(1/3)
  t_x <- (p_20^(3/4)/(4 * pi * N * p_02^(3/4) * (p_11 + sqrt(p_20 * 
                                                               p_02))))^(1/3)
  data_DCT_smoothed <- (t(t(exp(-(0:(n - 1))^2 * pi^2 * c(t_x)/2))) %*% 
                          t(exp(-(0:(n - 1))^2 * pi^2 * c(t_y)/2))) * data_DCT
  density <- idct2d(data_DCT_smoothed) * (numel(data_DCT_smoothed)/prod(scaling))
  grid <- meshgrid(seq(min_xy[1], max_xy[1], by = scaling[1]/(n - 
                                                                1)), seq(min_xy[2], max_xy[2], by = scaling[2]/(n - 
                                                                                                                  1)))
  bw <- sqrt(cbind(t_x, t_y)) * scaling
  return(list(bandwidth = bw, density = density, X = grid[[1]], 
              Y = grid[[2]]))
}

#funtion that estimates the DREMI score - adopted from the original DREMI implementation in MATLAB (c)
compute_dremi <- function(m, markerX,markerY,noise_threshold,set_maxy,maxy,makePlot,myLabel){
  #start with the basic parameters from the original script
  compute_pvalue <- 0;
  set_maxy  <- set_maxy;
  num_permutations <- 0; 
  max_yval <- 0; 
  num_slices <- 8; 
  non_kde_style <- 0;
  min_yval <- 0;
  maxy_val <- maxy
  use_min_y <- max(c(0, min(m[,2],na.rm = TRUE)));
  
  if(non_kde_style==0){
    
    if(set_maxy==0){
      #by default
      num_slices <- 256
      myMethod <- 'non_averaged_pts' 
      fixy <- 0
      minyval <- 0
      maxyval <- max(m[,2],na.rm = TRUE)
      fix_limits <- 0
      Limits <- c()
      Z <- pairwise_visualise(m,markerX,markerY,myMethod, 
                              noise_threshold,
                              num_slices,
                              fixy,minyval,maxyval,
                              fix_limits,Limits,
                              makePlot,myLabel);
    }
    else{
      #pairwise_visualize(channel1_name, channel2_name,'non_averaged_pts', noise_threshold,'no_plot', 'MinMaxY', use_min_y, maxy_val);
      #by default
      num_slices <- 256
      myMethod <- 'non_averaged_pts' 
      fixy <- 1
      minyval <- 0
      maxyval <- max(m[,2],na.rm = TRUE)
      fix_limits <- 0
      Limits <- c()
      minyval <- use_min_y
      maxyval <- maxy
      Z <- pairwise_visualise(m,markerX,markerY,
                              myMethod, 
                              noise_threshold,
                              num_slices,
                              fixy,minyval,maxyval,
                              fix_limits,Limits,makePlot,myLabel);
      
    }
    points_x <- Z[[1]]
    points_y <- Z[[2]]
    point_weights <- Z[[3]]
    
    total_slice_samples <- sum(point_weights) * 1000;
    samples_x <- c();
    samples_y <- c();
    
    for(i in 1:length(points_x)){
      num_repeats <- floor(point_weights[i] * 1000);
      new_samples_x <- matrix(1,nrow = num_repeats,ncol=1)*points_x[i];
      new_samples_y <- matrix(1,nrow = num_repeats,ncol=1)*points_y[i];
      
      samples_x = rbind(samples_x, new_samples_x);
      samples_y = rbind(samples_y, new_samples_y);
    }
    data <- cbind(samples_x, samples_y);
    
    minx <- min(data[,1])
    miny <- min(data[,2])
    maxx <- max(data[,1])
    maxy <- max(data[,2])
    
  }else{
    
    res <- find_data_cutoffs(X,Y, 50, 255);
    minx <- res[[1]]
    miny <- res[[3]]
    maxx <- res[[2]]
    maxy <- res[[4]]
    
  }
  
  if(set_maxy==1){
    maxy <- maxy_val;
  }
  
  if(non_kde_style == 0){
    num_slices <- 32
    out <- delta_entropyweight_rawdata(data, minx, miny, maxx, maxy,num_slices,num_slices)
  }else{
    #dremi = delta_entropyreweight(obj, channel1_name, channel2_name, num_slices, num_slices, minx, miny, maxx, maxy)
    
  }
  pvalue <- 0;
  if(num_permutations == 0){
    compute_pvalue <- 0; 
  }
  if(compute_pvalue ==1){
    #compute p value after permutations
    #TODO once finish with delta_entropy
  }
  return(out)
  
}


#argument m is a 2D matrix with X,and Y markers after removing outliers
#markerX and markerY are the names of markers, just to help visualisation labels
#Limits is an array in the format (minx,miny,maxx,maxy)
pairwise_visualise <- function(m,markerX,markerY,
                               myMethod,
                               noise_threshold,
                               num_slices,
                               fixy,minyval,maxyval,
                               fix_limits,Limits,
                               makePlot,myLabel)
  
  
  
  
{
  mySample <- myLabel
  #str <- paste('Processing sample:',mySample,' ',markerX,' vs. ',markerY,sep='')
  #print(str)
  
  X <- m[,1]
  Y <- m[,2]
  
  if(myMethod=='non_averaged_pts'){
    avg_pts <- 0
    avg_pts_threshold <- noise_threshold;
  }else{
    avg_pts <- 1
    avg_pts_threshold <- 0.9;
  }
  
  #utility variables from the original Matlab code
  num_slices <- num_slices;
  minxval <- 0;
  minyval <- minyval;
  cutoff <- 50;
  draw_contour <- 0;
  show_density <- 0;
  #make the plot the first time
  draw_plot <- makePlot;
  avg_pts <- avg_pts;
  fix_limits <- fix_limits;
  maxyval <- maxyval ;
  maxxval <- max(m[,1],na.rm = TRUE);
  fixy <- fixy;
  visual_threshold <- 0; 
  axes_specified <- 0; 
  axes_handle <- 0;
  
  
  if (fix_limits == 0){
    Xmin = minxval;
    Xmax = maxxval;
    Ymin = minyval;
    Ymax = maxyval;
  }else{
    minx <- Limits[1]
    miny <- Limits[2]
    maxx <- Limits[3]
    maxy <- Limits[4]
    
    Xmin = minx;
    Xmax = maxx;
    Ymin = miny;
    Ymax = maxy;
  }
  
  if(fixy == 1){
    Ymin = minyval;
    Ymax = maxyval;
  }
  
  #compute the 2D density,
  if(fix_limits == 0){
    Z <- my_kde2D(m,n=num_slices,limits=c(Xmin,Xmax,Ymin,Ymax))
  }else{
    #here we make a change because the code does not work with reduced box
    Z <- my_kde2D(m,n=num_slices,limits=c(minx,maxx,miny,maxy))
  }
  
  #get the variables we need
  bandwidth <- Z$bandwidth
  density <- Z$density
  Grid_X <- Z$X
  Grid_Y <- Z$Y
  
  #round small values of density
  density[density<0.00001] <- 0
  
  #weeding out spare ends
  row <- 0;
  total_cells_so_far <- 0; 
  while(total_cells_so_far<cutoff & row < (num_slices-1)){
    row <- row+1
    total_cells_so_far <- length(find(Y>Grid_Y[num_slices-row,1]))
  }
  maxy <- Grid_Y[num_slices-row,1]
  
  total_cells_so_far <- 0;
  start_row <- 0
  while(total_cells_so_far<cutoff){
    start_row <- start_row + 1
    total_cells_so_far <- length(find(Y<Grid_Y[start_row,1]))
  }
  miny <- Grid_Y[start_row,1]
  
  
  total_cells_so_far <- 0;
  col <- 0 
  while(total_cells_so_far<cutoff & col < (num_slices/2-1)){
    col <- col+1
    total_cells_so_far <- length(find(X>Grid_X[1,num_slices-col]))
  }
  maxx = Grid_X[1,num_slices-col];
  
  total_cells_so_far <- 0;
  start_col <- 0;
  while(total_cells_so_far<cutoff & start_col <= (num_slices/2)){
    start_col = start_col+1;
    total_cells_so_far = length(which(X<Grid_X[1,start_col]))
  }
  minx = Grid_X[1,start_col]
  
  if(fix_limits == 1){
    start_row <- 1
    start_col <- 1
    row <- 0
    col <- 0
  }
  
  if(fixy==1){
    start_row <- 1
    row <- 0
  }
  #filter density
  density <- density[start_row:(num_slices-row),start_col:(num_slices-col)]
  num_cols <- ncol(density)
  num_rows <- nrow(density)
  xaxis = Grid_X[1,start_col:(num_slices-col)]
  
  yaxis = Grid_Y[start_row:(num_slices-row),1]
  
  normalized_density <- matrix(0L,nrow = num_rows,ncol = num_cols)
  prob_normalized_density <- matrix(0L,nrow = num_rows,ncol = num_cols)
  #normalized by column for plotting the data 
  for(i in 1:num_cols){
    #normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
    normalized_density[,i] <- density[,i]/max(density[,i],na.rm=TRUE)
    prob_normalized_density[,i] <- density[,i]/Norm(density[,i])
  }
  
  #now create the side bars
  colsum <- colSums(density,na.rm = TRUE);
  normalized_colsum <- colsum/max(colsum,na.rm = TRUE);
  rowsum <- rowSums(density,na.rm = TRUE);
  normalized_rowsum <- rowsum/max(rowsum,na.rm = TRUE);
  
  #the corner is a fudge 
  #blueval = min(normalized_colsum);
  blueval <- 0;
  corner <- matrix(1,nrow = 11,ncol = 11)*blueval;
  #make the top bar
  #yaxis_increment = abs(yaxis(2)-yaxis(1,1));
  yaxis_increment <- 0.01;
  yaxis_top_bar <- c();
  top_bar <- c();
  zero_vector <- matrix(0L,nrow = 1,ncol = length(normalized_colsum));
  for(i in 1:1){
    top_bar <- rbind(top_bar, zero_vector)
    yaxis_top_bar <- rbind(yaxis_top_bar, max(yaxis)+(yaxis_increment*i))
  }
  for(i in 1:10){
    top_bar <- rbind(top_bar, normalized_colsum);
    yaxis_top_bar <- rbind(yaxis_top_bar, max(yaxis)+(yaxis_increment*i));  
  }
  
  #make the side bar
  #xaxis_increment = abs(xaxis(2)-xaxis(1));
  xaxis_increment <- 0.01
  xaxis_side_bar <- c()
  side_bar <- c()
  zero_vector <- matrix(0L,nrow=length(normalized_rowsum),ncol=1)
  for(i in 1:1){
    side_bar <- cbind(side_bar, zero_vector); 
    xaxis_side_bar <- cbind(xaxis_side_bar, max(xaxis)+(xaxis_increment*i))    
  }
  for(i in 1:10){
    side_bar <- cbind(side_bar, normalized_rowsum)
    xaxis_side_bar <- cbind(xaxis_side_bar, max(xaxis)+(xaxis_increment*i))
  }
  
  #find the trace through the peak regions for the return value
  points_x <- c()
  points_y <- c()
  point_weights <- c()
  if(avg_pts==1){
    for(i in 1:num_cols){
      max_indices <- which(normalized_density[,i]>= avg_pts_threshold)
      # points_y = [points_y mean(Grid_Y(max_indices,i))];
      points_x <- cbind(points_x, xaxis[i])
      #new_point_y = dot(Grid_Y(start_row+max_indices,start_col+i),normalized_density(max_indices,i));
      points_y <- cbind(points_y, mean(yaxis[max_indices]) )
      #points_y = [points_y new_point_y];
    }
    point_weights <- matrix(1,nrow = 1,ncol = length(points_y))
  }else{
    #this part of the code implements the non average points method
    #it is more complicated and it will be skipped for now, since we are working with default parameters
    for(i in 1:num_cols){
      max_indices <- find(normalized_density[,i]>= avg_pts_threshold)
      new_points <- matrix(1,nrow=1,ncol=length(max_indices))*xaxis[i];
      new_point_weights <- t(normalized_density[max_indices,i])
      new_point_weights <- new_point_weights / (sum(new_point_weights));
      points_x <- cbind(points_x, new_points);
      y_indices <- max_indices;
      new_points_y <- t(yaxis[y_indices]);
      points_y <- cbind(points_y, new_points_y);
      point_weights <- cbind(point_weights, new_point_weights);
    }
  }
  
  smoothed_normalized_density <- matrix(0L,nrow = num_rows,ncol = num_cols)
  for(i in 1:num_rows){
    for(j in 2:(num_cols-1)){
      smoothed_normalized_density[i,j] = (normalized_density[i,j-1]+normalized_density[i,j]+normalized_density[i,j+1])/3; 
    }
  }
  
  if(visual_threshold>0){
    smoothed_normalized_density <- (smoothed_normalized_density>visual_threshold)*smoothed_normalized_density;
  }
  
  matrix_to_plot <- cbind(smoothed_normalized_density, side_bar);
  top_bar <- cbind(top_bar, corner)
  matrix_to_plot <- rbind(matrix_to_plot,top_bar)
  #xaxis_side_bar[1] <- xaxis_side_bar[1]-0.01
  #yaxis_side_bar[1] <- yaxis_side_bar[1]-0.01
  xaxis_to_plot <- c(xaxis, xaxis_side_bar)
  yaxis_to_plot <- c(yaxis, yaxis_top_bar)
  
  #a1 <- sample(1:length(xaxis_to_plot),108)
  #a1 <- xaxis_to_plot[a1]
  
  #a2 <- sample(1:length(yaxis_to_plot),108)
  #a2 <- yaxis_to_plot[a2]
  
  #FROM NOW ON WE MAKE THE PLOTS
  #make the plot
  if(draw_plot==1){
    
    if(axes_specified==1){
      #here this part of the code is quite similar to the next....but we work with defaults
    }else{
      #this is the matlab plot
      #imagesc(matrix_to_plot,col=jet.colors(12),xlab='marker1',ylab='marker2');
      
      #save the plot
      myfile<-paste(plot_dir,myLabel,'.pdf',sep='')
      pdf(myfile,onefile = TRUE)
      image(x=c(1:ncol(matrix_to_plot)),y=c(1:nrow(matrix_to_plot)),t(matrix_to_plot),
            axes=FALSE,xlab=markerX,ylab=markerY,
            col = hcl.colors(12, "Spectral", rev = TRUE))
      d <- c(1:ncol(smoothed_normalized_density))
      d <- quantile(d)
      #1,50,100,150,200,ncol(smoothed_normalized_density)
      a <- c(d[1],d[2],d[3],d[4],d[5])
      a <- round(a)
      axis(1,at=a,labels = round(xaxis_to_plot[a],2))
      #a <- c(1,50,100,150,200,nrow(smoothed_normalized_density))
      d <- c(1:nrow(smoothed_normalized_density))
      d <- quantile(d)
      a <- c(d[1],d[2],d[3],d[4],d[5])
      a <- round(a)
      axis(2,at=a,labels = round(yaxis_to_plot[a],2))
      dev.off()
      
    }
    #additional plot with two colors not implemented
    #density_filtered_v1 <- matrix_to_plot>0.6
    #matrix_to_plot_v1 <- matrix_to_plot*density_filtered_v1;
    #imagesc(matrix_to_plot)
  }
  #draw contours
  if(draw_contour==1){
    Z <- kde2D(m,num_slices+1,limits=c(minx, maxx,miny, maxy))
    rdensity <- Z$density
    rGrid_X <- Z$X
    rGrid_Y <- Z$Y
    
  }
  #this can be imporeved using ggplot - TODO
  if(show_density==1){
    f <- kde(m[,1]);
    Fx <- approx(f[1,],f[2,],xout=points_x)
    Fx <- Fx$y
    Fx[Fx<0.00001] <- 0
    plot(points_x,Fx,col='black',cex=1);
  }
  
  out <- list(points_x,points_y,point_weights,density,normalized_density,
              xaxis_to_plot,yaxis_to_plot,top_bar,side_bar)
  return(out)
  
}

# adopted from the original DREMI implementation in MATLAB
delta_entropyweight_rawdata <- function(data, minx, miny, maxx, maxy, num_slicesx,  num_slicesy){
  
  fixed_partition <- 0;
  
  #skip the part with extra arguments
  
  if(fixed_partition == 0){
    xincrement <- (maxx-minx)/num_slicesx;
    yincrement <- (maxy-miny)/num_slicesy;
    edges_x <- c(minx);
    
    for(i in 1:(num_slicesx-1)){
      edges_x <- cbind(edges_x, edges_x[i]+xincrement);
    }
    edges_x <- cbind(edges_x, maxx);
    edges_y <- c(miny);
    
    for(i in 1:(num_slicesy-1)){
      edges_y <- cbind(edges_y, edges_y[i]+yincrement);
    }
    edges_y <- cbind(edges_y, maxy);
    edges <- list()
    edges[[1]] <- edges_x;
    edges[[2]] <- edges_y;
  }
  
  #this part of the code replaces the hist3 implemented in matlab
  
  nbins <- length(edges_x)
  x.bin <- edges_x
  y.bin <- edges_y
  freq <-  as.data.frame(table(findInterval(data[,1], x.bin,left.open=TRUE,rightmost.closed=TRUE),
                               findInterval(data[,2], y.bin,left.open=TRUE,rightmost.closed=TRUE)))
  freq[,1] <- as.numeric(freq[,1])
  freq[,2] <- as.numeric(freq[,2])
  freq2D <- diag(nbins)*0
  freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
  #need to be a bit careful here with the transpose operation
  hist_bins <- t(freq2D)
  hist_bin_totals = colSums(hist_bins);
  bin_totals <- t(hist_bin_totals);
  
  num_rows <- nbins
  num_cols <- nbins
  
  #here comes the entropy computations
  
  entropy_y <- 0;
  column_weights <-  matrix(1,nrow=1,ncol=num_cols);
  slice_entropies <- matrix(0,nrow=1,ncol=num_cols);
  total <- sum(column_weights[1,]);
  
  
  for(i in 1:num_rows){
    p_y = 0; 
    for(j in 1:num_cols){
      column_total = sum(hist_bins[,j]);
      if(column_total > 0){
        p_y = p_y + (hist_bins[i,j]*(column_weights[j]/column_total));
      }
    }
    p_y <- p_y/total;
    if(p_y==0){
      p_y <- 0;
      log_p_y <- 0;
    }else{
      log_p_y <- log2(p_y);
    }
    entropy_y <- entropy_y+(p_y * log_p_y);
  }#end of for loop
  entropy_y <- -1*entropy_y;
  cond_entropy_y  <- 0
  
  
  for(j in 1:num_cols){
    
    slice_entropy_y <- 0;
    col_sum <- sum(hist_bins[,j]);
    if(col_sum == 0){
      next;
    }
    norm_col_hist <- hist_bins[,j]/col_sum;
    
    
    for(i in 1:num_rows){
      slice_p_y <- norm_col_hist[i];
      if(slice_p_y==0){
        slice_p_y <- 0;
        log_slice_p_y <- 0;
      }else{
        slice_p_y;
        log_slice_p_y <- log2(slice_p_y);
      }
      slice_entropy_y = slice_entropy_y + (slice_p_y * log_slice_p_y);
    }
    slice_entropy_y <- -1*slice_entropy_y;
    slice_entropies[j] <- slice_entropy_y;
    col_prob <- column_weights[1,j]/sum(column_weights[1,]);
    cond_entropy_y <- cond_entropy_y + (col_prob * slice_entropy_y);
    #print(cond_entropy_y)
  }
  
  #      slice_entropies
  #      entropy_y
  delta_before_division <- entropy_y - cond_entropy_y;
  cond_entropy_y_before_division <- cond_entropy_y;
  #normalize this by the size of the grid 
  smaller_dim <- min(num_rows, num_cols);
  log_smaller_dim <- log2(smaller_dim);
  #this is supposed to be the maximum achievable value on a grid this
  #size
  
  delta <- delta_before_division/log_smaller_dim
  cond_entropy_y <- cond_entropy_y/log_smaller_dim
  
  out <- list(delta, entropy_y, cond_entropy_y, hist_bins, data, edges)
  return(out)
}   



