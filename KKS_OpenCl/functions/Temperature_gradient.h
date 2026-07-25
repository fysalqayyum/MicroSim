#ifndef TEMPERATURE_GRADIENT_H_
#define TEMPERATURE_GRADIENT_H_

void apply_temperature_gradientY(struct fields* gridinfo, long shift_OFFSET, long t) {
  long x, y, z;
  long gidy;
//   double BASE_POS=0;
//   double GRADIENT;
//   double temp_bottom;
//    
//   BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET + ((temperature_gradientY.velocity/deltay)*(t*deltat));
//   GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
//   temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT + (workers_mpi.offset[Y]-workers_mpi.offset_y)*GRADIENT;
  
  /*
   * gridinfo is the MPI-local, ghosted array after mpi_xy().  The original
   * implementation used the global, unghosted rows_* dimensions and
   * layer_size, leaving part of the local array uninitialised.  The first
   * physical y cell is local y=1; initialise the ghost cells as the linear
   * continuation of the imposed gradient so every value copied to OpenCL is
   * finite before boundary kernels run.
   */
  for(x=0; x < mpiparam.rows_x; x++) {
   for(z=0; z < mpiparam.rows_z; z++) {
    for(y=0; y < mpiparam.rows_y; y++) {
      gidy = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
      gridinfo[gidy].temperature = temp_bottom + (y - 1)*GRADIENT;
//       gridinfo_w[gidy].temperature  = temperature_gradientY.base_temp;
//       if (gridinfo1[gidy].temperature > 1.0) {
//          gridinfo1[gidy].temperature = 1.0;
//       } else if (gridinfo1[gidy].temperature < temperature_gradientY.base_temp) {
//          gridinfo[gidy].temperature = temperature_gradientY.base_temp;
//       }
      }
    }
  }
}
#endif
