void CL_Shift() {
  long local_interface_position;
  long global_interface_position;
  long shift_cells;

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
  shift_OFFSET += shift_cells;
  apply_shiftY_local(gridinfomN, cscl, shift_cells, shift_OFFSET,
                     t + STARTTIME);

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
