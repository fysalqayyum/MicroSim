void CL_Shift() {
  long local_interface_position;
  long global_interface_position;
  long shift_cells;
  long x, y, z, idx;
  long local_refill_cells, global_refill_cells;
  double local_before = 0.0, local_discarded = 0.0;
  double local_after = 0.0, local_refill_phase_max = 0.0;
  double local_refill_xsi_max = 0.0, local_refill_temp_max = 0.0;
  double global_before, global_discarded, global_after;
  double global_refill_phase_max, global_refill_xsi_max;
  double global_refill_temp_max, expected_after, residual;

  CL_DeviceToHost();
  ret = clEnqueueReadBuffer(cmdQ, d_cscl, CL_TRUE, 0,
                            nxnynz*sizeof(struct csle), cscl,
                            0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read cscl for moving window: %d\n", ret);
    exit(1);
  }

  local_interface_position = check_SHIFT_local(gridinfomN);
  MPI_Allreduce(&local_interface_position, &global_interface_position,
                1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);

  if (global_interface_position <= shiftj) {
    return;
  }

  shift_cells = global_interface_position - shiftj;
  local_refill_cells =
      (mpiparam.rows_x-2)*(mpiparam.rows_z-2)*shift_cells;

  for (x = 1; x < mpiparam.rows_x-1; x++) {
    for (z = 1; z < mpiparam.rows_z-1; z++) {
      for (y = 1; y < mpiparam.rows_y-1; y++) {
        idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        local_before += gridinfomN[idx].composition[0];
        if (y <= shift_cells)
          local_discarded += gridinfomN[idx].composition[0];
      }
    }
  }

  shift_OFFSET += shift_cells;
  /*
   * The refill temperature must use the GLOBAL cumulative shift.
   * shift_OFFSET is reset to 0 on every restart (see the TEMPGRADY block in
   * microsim_kks_opencl.c), so on a restarted run it counts only the shifts of
   * the current leg, while the timestep passed below is absolute.  Mixing the
   * two left the refilled rows colder than the imposed gradient by exactly the
   * shift accumulated in all previous legs.  shift_position holds that
   * pre-restart total (0 on a cold start), so shift_position + shift_OFFSET is
   * the global shift at all times -- the same identity farfield_guard.h uses to
   * reconstruct the absolute front position.
   */
  apply_shiftY_local(gridinfomN, cscl, shift_cells,
                     shift_position + shift_OFFSET,
                     t + STARTTIME);

  for (x = 1; x < mpiparam.rows_x-1; x++) {
    for (z = 1; z < mpiparam.rows_z-1; z++) {
      for (y = 1; y < mpiparam.rows_y-1; y++) {
        double expected_temperature;
        double error;
        idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        local_after += gridinfomN[idx].composition[0];
        if (y > mpiparam.rows_y-2-shift_cells) {
          /*
           * Same global-shift identity as the refill above.  This audit used to
           * repeat the refill's per-leg shift_OFFSET, so the two agreed with
           * each other and reported a refill temperature error of exactly zero
           * on every restarted leg -- which is why this check never caught the
           * stale-row defect.  Both sides must use shift_position + shift_OFFSET
           * for the audit to be independent of the bug it is meant to detect.
           */
          expected_temperature =
              temperature_gradientY.base_temp
              + (temperature_gradientY.DeltaT/temperature_gradientY.Distance)
                  * (((y-1) + shift_position + shift_OFFSET)*deltay
                     - temperature_gradientY.velocity
                           * (t + STARTTIME)*deltat);
          if (gridinfomN[idx].phia[0] > local_refill_phase_max)
            local_refill_phase_max = gridinfomN[idx].phia[0];
          error = fabs(gridinfomN[idx].composition[0]
                       - cfill[NUMPHASES-1][NUMPHASES-1][0]);
          if (error > local_refill_xsi_max)
            local_refill_xsi_max = error;
          error = fabs(gridinfomN[idx].temperature-expected_temperature);
          if (error > local_refill_temp_max)
            local_refill_temp_max = error;
        }
      }
    }
  }

  MPI_Reduce(&local_before, &global_before, 1, MPI_DOUBLE, MPI_SUM,
             MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&local_discarded, &global_discarded, 1, MPI_DOUBLE, MPI_SUM,
             MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&local_after, &global_after, 1, MPI_DOUBLE, MPI_SUM,
             MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&local_refill_phase_max, &global_refill_phase_max, 1,
             MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&local_refill_xsi_max, &global_refill_xsi_max, 1,
             MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&local_refill_temp_max, &global_refill_temp_max, 1,
             MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&local_refill_cells, &global_refill_cells, 1,
             MPI_LONG, MPI_SUM, MASTER, MPI_COMM_WORLD);

  if (rank == MASTER) {
    FILE *audit;
    expected_after =
        global_before - global_discarded
        + global_refill_cells
            * cfill[NUMPHASES-1][NUMPHASES-1][0];
    residual = global_after - expected_after;
    audit = fopen("DATA/moving_window_events.csv", "a+");
    if (audit == NULL) {
      perror("DATA/moving_window_events.csv");
      exit(1);
    }
    fseek(audit, 0, SEEK_END);
    if (ftell(audit) == 0)
      fprintf(audit,
              "step,shift_cells,cumulative_shift_cells,"
              "inventory_before,discarded_inventory,expected_after,"
              "inventory_after,inventory_residual,"
              "inventory_relative_residual,refill_FCC_A1_max,"
              "refill_xSi_max_abs_error,refill_temperature_max_abs_error_K\n");
    fprintf(audit,
            "%ld,%ld,%ld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g\n",
            t + STARTTIME, shift_cells, shift_position + shift_OFFSET,
            global_before, global_discarded, expected_after, global_after,
            residual, residual/fmax(fabs(expected_after), 1.0e-30),
            global_refill_phase_max, global_refill_xsi_max,
            global_refill_temp_max);
    fclose(audit);
  }

  /*
   * Deliberately the PER-LEG shift_OFFSET, not the global one.  On restart the
   * device reconstructs the global shift through TPosOffset, which is derived
   * from gradient_OFFSET after that has absorbed -shift_position*deltay.
   * Adding shift_position here as well would double-count it.
   */
  pfmdat.shift_OFFSET = shift_OFFSET;

  ret = clEnqueueWriteBuffer(cmdQ, d_pfmdat, CL_TRUE, 0,
                             sizeof(struct pfmval), &pfmdat,
                             0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to upload moving-window parameters: %d\n", ret);
    exit(1);
  }
  ret = clEnqueueWriteBuffer(cmdQ, d_gridinfomN, CL_TRUE, 0,
                             nxnynz*sizeof(struct fields), gridinfomN,
                             0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to upload shifted fields: %d\n", ret);
    exit(1);
  }
  ret = clEnqueueWriteBuffer(cmdQ, d_cscl, CL_TRUE, 0,
                             nxnynz*sizeof(struct csle), cscl,
                             0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to upload shifted phase compositions: %d\n", ret);
    exit(1);
  }
}
