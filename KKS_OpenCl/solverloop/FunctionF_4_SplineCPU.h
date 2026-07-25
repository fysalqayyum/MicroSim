static double f4_cubic4_host(double p0, double p1, double p2, double p3,
                             double coordinate) {
  return -p0*(coordinate-1.0)*(coordinate-2.0)*(coordinate-3.0)/6.0
         +p1*coordinate*(coordinate-2.0)*(coordinate-3.0)/2.0
         -p2*coordinate*(coordinate-1.0)*(coordinate-3.0)/2.0
         +p3*coordinate*(coordinate-1.0)*(coordinate-2.0)/6.0;
}

void FunctionF_4_SplineCPU(struct propmatf4spline *table) {
  long it;
  int ip, is, is1, is2;
  double sample_temperature;
  double reference_temperature = T;
  double max_relative_A = 0.0;
  double max_relative_B = 0.0;
  double max_relative_C = 0.0;

  if (f4_table_count < 2 || f4_table_dT <= 0.0) {
    fprintf(stderr,
            "Invalid Function_F=4 lookup table: count=%ld dT=%le\n",
            f4_table_count, f4_table_dT);
    exit(1);
  }

  for (it=0; it<f4_table_count; it++) {
    sample_temperature = f4_table_Tmin + it*f4_table_dT;
    function_A(sample_temperature, ceq);

    for (ip=0; ip<npha; ip++) {
      for (is1=0; is1<nsol; is1++) {
        for (is2=0; is2<nsol; is2++) {
          table[it].A[ip][is1][is2] = A[ip][is1][is2];
        }
      }
      for (is=0; is<nsol; is++) {
        table[it].B[ip][is] = function_B(sample_temperature, is, ip);
      }
      table[it].C[ip] = function_C(sample_temperature, ip);
    }
  }

  /* The legacy second buffer is retained for ABI compatibility. */
  memcpy(propf4spline1, table,
         f4_table_count*sizeof(struct propmatf4spline));

  /*
   * Audit the device's piecewise-linear interpolation against the original
   * GSL cubic splines at every table-cell midpoint.
   */
  for (it=0; it<f4_table_count-1; it++) {
    double interp_temperature =
        f4_table_Tmin + (it+0.5)*f4_table_dT;
    long start = it-1;
    double coordinate;
    if (start < 0) start = 0;
    if (start > f4_table_count-4) start = f4_table_count-4;
    coordinate = (double)(it-start)+0.5;
    function_A(interp_temperature, ceq);
    for (ip=0; ip<npha; ip++) {
      for (is1=0; is1<nsol; is1++) {
        for (is2=0; is2<nsol; is2++) {
          double linear = f4_cubic4_host(
              table[start].A[ip][is1][is2],
              table[start+1].A[ip][is1][is2],
              table[start+2].A[ip][is1][is2],
              table[start+3].A[ip][is1][is2], coordinate);
          double scale = fmax(fabs(A[ip][is1][is2]), 1.0);
          double error = fabs(linear-A[ip][is1][is2])/scale;
          if (error > max_relative_A) max_relative_A = error;
        }
        {
          double direct = function_B(interp_temperature, is1, ip);
          double linear = f4_cubic4_host(
              table[start].B[ip][is1], table[start+1].B[ip][is1],
              table[start+2].B[ip][is1], table[start+3].B[ip][is1],
              coordinate);
          double error = fabs(linear-direct)/fmax(fabs(direct), 1.0);
          if (error > max_relative_B) max_relative_B = error;
        }
      }
      {
        double direct = function_C(interp_temperature, ip);
        double linear = f4_cubic4_host(
            table[start].C[ip], table[start+1].C[ip],
            table[start+2].C[ip], table[start+3].C[ip], coordinate);
        double error = fabs(linear-direct)/fmax(fabs(direct), 1.0);
        if (error > max_relative_C) max_relative_C = error;
      }
    }
  }

  if (rank == MASTER) {
    long mid = f4_table_count/2;
    double ref_coordinate =
        (reference_temperature-f4_table_Tmin)/f4_table_dT;
    long ref_lower = (long)floor(ref_coordinate);
    double ref_weight;
    if (ref_lower < 0) ref_lower = 0;
    if (ref_lower > f4_table_count-2) ref_lower = f4_table_count-2;
    ref_weight = ref_coordinate-ref_lower;
    if (ref_weight < 0.0) ref_weight = 0.0;
    if (ref_weight > 1.0) ref_weight = 1.0;
    printf("F4_THERMO_TABLE count=%ld Tmin=%.12g Tmax=%.12g dT=%.12g\n",
           f4_table_count, f4_table_Tmin, f4_table_Tmax, f4_table_dT);
    printf("F4_THERMO_INTERP max_relative_A=%.12e "
           "max_relative_B=%.12e max_relative_C=%.12e\n",
           max_relative_A, max_relative_B, max_relative_C);
    for (ip=0; ip<npha; ip++) {
      printf("F4_THERMO_PHASE phase=%d "
             "A00={%.12e,%.12e,%.12e} "
             "B0={%.12e,%.12e,%.12e} "
             "C={%.12e,%.12e,%.12e}\n",
             ip,
             table[0].A[ip][0][0],
             table[mid].A[ip][0][0],
             table[f4_table_count-1].A[ip][0][0],
             table[0].B[ip][0],
             table[mid].B[ip][0],
             table[f4_table_count-1].B[ip][0],
             table[0].C[ip],
             table[mid].C[ip],
             table[f4_table_count-1].C[ip]);
    }
    function_A(reference_temperature, ceq);
    for (ip=0; ip<npha; ip++) {
      double direct_B = function_B(reference_temperature, 0, ip);
      double direct_C = function_C(reference_temperature, ip);
      long ref_start = ref_lower-1;
      double ref_coordinate;
      if (ref_start < 0) ref_start = 0;
      if (ref_start > f4_table_count-4) ref_start = f4_table_count-4;
      ref_coordinate = (double)(ref_lower-ref_start)+ref_weight;
      double linear_A = f4_cubic4_host(
          table[ref_start].A[ip][0][0], table[ref_start+1].A[ip][0][0],
          table[ref_start+2].A[ip][0][0], table[ref_start+3].A[ip][0][0],
          ref_coordinate);
      double linear_B = f4_cubic4_host(
          table[ref_start].B[ip][0], table[ref_start+1].B[ip][0],
          table[ref_start+2].B[ip][0], table[ref_start+3].B[ip][0],
          ref_coordinate);
      double linear_C = f4_cubic4_host(
          table[ref_start].C[ip], table[ref_start+1].C[ip],
          table[ref_start+2].C[ip], table[ref_start+3].C[ip],
          ref_coordinate);
      printf("F4_THERMO_REFERENCE phase=%d T=%.12g "
             "A={%.16e,%.16e,%.3e} "
             "B={%.16e,%.16e,%.3e} "
             "C={%.16e,%.16e,%.3e}\n",
             ip, reference_temperature,
             A[ip][0][0], linear_A, linear_A-A[ip][0][0],
             direct_B, linear_B, linear_B-direct_B,
             direct_C, linear_C, linear_C-direct_C);
    }
  }

  /* Restore the host-side reference-temperature matrices. */
  function_A(reference_temperature, ceq);
  for (ip=0; ip<npha; ip++) {
    for (is=0; is<nsol; is++) {
      B[ip][is] = function_B(reference_temperature, is, ip);
    }
    C[ip] = function_C(reference_temperature, ip);
  }
}
