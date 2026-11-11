// parallel section over each motif for this n,L combination, count can be quite large
#pragma omp parallel for num_threads(N_threads)
for(int i=0; i<count; i++) {
  // here we have a motif, part of a necklace, first one should be lexicographically lowest.
  struct motif* m = get_motif(M[n][L],i);
  // loop through the motif rotations
  while(m) {
    struct steptypes* steps = motif_to_steptypes(m);
    int steps_length = steptypes_length(steps);
    // prepare the output buffer
    //printf("allocating %d sized buffer\n",steps_length);
    struct output_buffer* output_buffer = create_output_buffer(steps_length);
    // initialize for averages
    double num = 0.0;
    double den = 0.0;
    
    // average over the rotations of the steptypes:
    for(int j=0; j<steps_length; j++) {
      struct array* C_eff_stats;
      if(self_avoiding) {
        C_eff_stats = get_C_eff_both_dir_sap(m,steps,rna_lib,N_samples,bin_size,d_min,helix_lib,sap_report);
      } else {
        C_eff_stats = get_C_eff_both_dir(steps,rna_lib,N_samples);
      }
      // unpack the stats
      double C_eff_mean = array_get(C_eff_stats,0);
      double C_eff_error = array_get(C_eff_stats,1);
      
      // write stats to the output buffer
      char buffer[500]; // allocated buffer to print the stats for this record
	    int offset = sprint_motif(buffer,m,'\t');
	    offset += sprint_steptypes(buffer+offset,steps,'\t');
	    offset += sprintf(buffer+offset,"%.8f\t%.8f\t%.8f\t%.8f\n",C_eff_mean,C_eff_error,dist_sq_mean,theta_sq_mean);
	    //printf("%d\t%d\t%s\n",j,offset,buffer);
	    output_buffer_append(output_buffer,buffer);
            
	    // compute relative error 
	    double C_eff_relerr = C_eff_error / C_eff_mean;	    
	    // logarithmic, inverse-relerr weighted average
	    num += log(C_eff_mean)/(C_eff_relerr*C_eff_relerr);
	    den += 1.0/(C_eff_relerr*C_eff_relerr);
	    // memory management:
	    free_array(C_eff_stats);
	    // rotate the steptypes and add to running average
	    steps = rotate_steptypes(steps);
    }
    // weighted average over steptypes rotations:
    double full_C_eff_mean = exp(num/den);
    double full_C_eff_err = sqrt(1.0/den)*full_C_eff_mean;
    
    // should print out each motif's steptype rotations and overall average:
    #pragma omp critical
    {
      // print all motif steptype rotations, result stored in output_buffer
      for(int j=0; j<steps_length; j++) {
        fprintf(OUT_ALL, "%s", output_buffer_get(output_buffer,j));
      }
      fflush(OUT_ALL);
      free_output_buffer(output_buffer);
      // print average of all rotations
      fprint_motif(OUT, m, '\t');
      fprintf(OUT,"%.8f\t%.8f\t%.8f\t%.8f\n",full_C_eff_mean,full_C_eff_err,full_dist_sq_mean,full_theta_sq_mean);    
      fflush(OUT);
    }
    m = get_next_motif(m); // rotate the motif to get an equivalent one in different rotation
    free_steptypes(steps);
  }
 }
}






for(int n=N_min; n<=N_max; n++) {
  for(int L=L_min; L<=L_max; L++) {
    // update standard output log
    int count = get_motif_group_count(M[n][L]);
    printf("# N_pred(n=%d,L=%d)=%d, N_obs=%d\n",n,L,N[n][L],count);
    fflush(stdout);
    
    // parallel section over each motif for this n,L combination, count can be quite large
    #pragma omp parallel for num_threads(N_threads)
    for(int i=0; i<count; i++) {
      // here we have a motif, part of a necklace, first one should be lexicographically lowest.
      struct motif* m = get_motif(M[n][L],i);
      // loop through the motif rotations
      while(m) {
        struct steptypes* steps = motif_to_steptypes(m);
        int steps_length = steptypes_length(steps);
        // prepare the output buffer
        //printf("allocating %d sized buffer\n",steps_length);
        struct output_buffer* output_buffer = create_output_buffer(steps_length);
        // initialize for averages
        double num = 0.0;
        double den = 0.0;
        double full_dist_sq_mean = 0.0;
        double full_theta_sq_mean = 0.0;
        
        // average over the rotations of the steptypes:
        for(int j=0; j<steps_length; j++) {
          struct array* C_eff_stats;
          if(self_avoiding) {
            C_eff_stats = get_C_eff_both_dir_sap(m,steps,rna_lib,N_samples,bin_size,d_min,helix_lib,sap_report);
          } else {
            C_eff_stats = get_C_eff_both_dir(steps,rna_lib,N_samples);
          }
          // unpack the stats
          double C_eff_mean = array_get(C_eff_stats,0);
          double C_eff_error = array_get(C_eff_stats,1);
          double dist_sq_mean = array_get(C_eff_stats,2);
          double theta_sq_mean = array_get(C_eff_stats,3);
          
          // write stats to the output buffer
          char buffer[500]; // allocated buffer to print the stats for this record
          int offset = sprint_motif(buffer,m,'\t');
          offset += sprint_steptypes(buffer+offset,steps,'\t');
          offset += sprintf(buffer+offset,"%.8f\t%.8f\t%.8f\t%.8f\n",C_eff_mean,C_eff_error,dist_sq_mean,theta_sq_mean);
          //printf("%d\t%d\t%s\n",j,offset,buffer);
          output_buffer_append(output_buffer,buffer);
          
          // compute relative error 
          double C_eff_relerr = C_eff_error / C_eff_mean;	    
          // logarithmic, inverse-relerr weighted average
          num += log(C_eff_mean)/(C_eff_relerr*C_eff_relerr);
          den += 1.0/(C_eff_relerr*C_eff_relerr);
          // averages over step-type rotations
          full_dist_sq_mean += dist_sq_mean;
          full_theta_sq_mean += theta_sq_mean;
          // memory management:
          free_array(C_eff_stats);
          // rotate the steptypes and add to running average
          steps = rotate_steptypes(steps);
        }
        // weighted average over steptypes rotations:
        double full_C_eff_mean = exp(num/den);
        double full_C_eff_err = sqrt(1.0/den)*full_C_eff_mean;
        full_dist_sq_mean /= steps_length;
        full_theta_sq_mean /= steps_length;
        
        // should print out each motif's steptype rotations and overall average:
        #pragma omp critical
        {
          // print all motif steptype rotations, result stored in output_buffer
          for(int j=0; j<steps_length; j++) {
            fprintf(OUT_ALL, "%s", output_buffer_get(output_buffer,j));
          }
          fflush(OUT_ALL);
          free_output_buffer(output_buffer);
          // print average of all rotations
          fprint_motif(OUT, m, '\t');
          fprintf(OUT,"%.8f\t%.8f\t%.8f\t%.8f\n",full_C_eff_mean,full_C_eff_err,full_dist_sq_mean,full_theta_sq_mean);    
          fflush(OUT);
        }
        m = get_next_motif(m); // rotate the motif to get an equivalent one in different rotation
        free_steptypes(steps);
      }
    }
  }
}
