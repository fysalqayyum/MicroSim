#ifndef APPLY_SHIFTY_H
#define APPLY_SHIFTY_H

long check_SHIFT_local(struct fields* local_grid) {
  long x, y, z, center;
  long interface_position = 0;

  for (x = 1; x < mpiparam.rows_x-1; x++) {
    for (z = 1; z < mpiparam.rows_z-1; z++) {
      for (y = 2; y < mpiparam.rows_y-1; y++) {
        center = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        if ((local_grid[center-1].phia[NUMPHASES-1] < 0.5) &&
            (local_grid[center].phia[NUMPHASES-1] >= 0.5) &&
            ((y-1) > interface_position)) {
          interface_position = y-1;
        }
      }
    }
  }
  return interface_position;
}

void apply_shiftY_local(struct fields* local_grid, struct csle* local_cscl,
                        long shift_cells, long new_shift_offset,
                        long absolute_timestep) {
  long x, y, z, src, dst, ip, is;
  const long first_y = 1;
  const long last_y = mpiparam.rows_y - 2;
  double c[NUMCOMPONENTS-1];
  double local_temperature;
  double imposed_gradient =
      temperature_gradientY.DeltaT/temperature_gradientY.Distance;

  if (shift_cells <= 0 || shift_cells >= mpiparam.rows_y-2) {
    return;
  }

  for (x = 0; x < mpiparam.rows_x; x++) {
    for (z = 0; z < mpiparam.rows_z; z++) {
      /*
       * Shift only physical rows.  rows_y includes one ghost row at each
       * boundary; treating the upper ghost as physical used to refill one
       * fewer physical row than was removed.
       */
      for (y = first_y; y <= last_y-shift_cells; y++) {
        dst = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        src = (y + shift_cells)
              + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        local_grid[dst] = local_grid[src];
        local_cscl[dst] = local_cscl[src];
      }

      for (y = last_y-shift_cells+1; y <= last_y; y++) {
        dst = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        local_temperature =
            temperature_gradientY.base_temp
            + imposed_gradient
                * (((y-1) + new_shift_offset)*deltay
                   - temperature_gradientY.velocity
                         * absolute_timestep*deltat);

        for (ip = 0; ip < NUMPHASES-1; ip++) {
          local_grid[dst].phia[ip] = 0.0;
        }
        local_grid[dst].phia[NUMPHASES-1] = 1.0;
        for (is = 0; is < NUMCOMPONENTS-1; is++) {
          c[is] = cfill[NUMPHASES-1][NUMPHASES-1][is];
          local_grid[dst].composition[is] = c[is];
        }
        Mu(c, local_temperature, NUMPHASES-1, local_grid[dst].compi);
        local_grid[dst].temperature = local_temperature;

        for (ip = 0; ip < NUMPHASES; ip++) {
          for (is = 0; is < NUMCOMPONENTS-1; is++) {
            local_cscl[dst].comie[ip][is] = ceq[ip][ip][is];
          }
        }
      }

      /* Restore no-flux y ghosts immediately after the host-side shift. */
      local_grid[mpiparam.rows_y*(z + mpiparam.rows_z*x)] =
          local_grid[first_y + mpiparam.rows_y*(z + mpiparam.rows_z*x)];
      local_cscl[mpiparam.rows_y*(z + mpiparam.rows_z*x)] =
          local_cscl[first_y + mpiparam.rows_y*(z + mpiparam.rows_z*x)];
      local_grid[(mpiparam.rows_y-1)
                 + mpiparam.rows_y*(z + mpiparam.rows_z*x)] =
          local_grid[last_y + mpiparam.rows_y*(z + mpiparam.rows_z*x)];
      local_cscl[(mpiparam.rows_y-1)
                 + mpiparam.rows_y*(z + mpiparam.rows_z*x)] =
          local_cscl[last_y + mpiparam.rows_y*(z + mpiparam.rows_z*x)];
    }
  }
}

void apply_shiftY(struct fields* gridinfo, long INTERFACE_POS_GLOBAL) {
  //Shift by one cell in the negative y-direction
  long x, y, z;
  long gidy;
  double chemical_potential;
  double c[NUMCOMPONENTS-1];
  
  for(x=0; x < rows_x; x++) {
    for(z=0; z < rows_z; z++) {
      for (y=0; y <= (rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y++) {
        gidy = x*layer_size + z*rows_y + y;
        for (b=0; b < NUMPHASES; b++) {
          gridinfo[gidy].phia[b] = gridinfo[gidy+(INTERFACE_POS_GLOBAL-shiftj)].phia[b];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          gridinfo[gidy].compi[k] = gridinfo[gidy+(INTERFACE_POS_GLOBAL-shiftj)].compi[k];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          gridinfo[gidy].composition[k] = gridinfo[gidy+(INTERFACE_POS_GLOBAL-shiftj)].composition[k];
        }
        gridinfo[gidy].temperature = gridinfo[gidy + (INTERFACE_POS_GLOBAL-shiftj)].temperature;
      }
//       if (workers_mpi.lasty==1) {
      for (y=(rows_y-(INTERFACE_POS_GLOBAL-shiftj)); y<=(rows_y-1); y++) {
        gidy = x*layer_size + z*rows_y + y;
        for (b=0; b < NUMPHASES-1; b++) {
          gridinfo[gidy].phia[b] = 0.0;
        }
        gridinfo[gidy].phia[NUMPHASES-1] = 1.0;
        for (k=0; k < NUMCOMPONENTS-1; k++) {
//          c[k] = ceq[NUMPHASES-1][NUMPHASES-1][k];
          //c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
          c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
        }
        Mu(c, Teq, NUMPHASES-1, gridinfo[gidy].compi); 
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          //chemical_potential         = Mu(c, Tfill, NUMPHASES-1, k);
          gridinfo[gidy].composition[k]    = c[k]; //chemical_potential;
        }
      }
//       }
    }
  }
//   init_propertymatrices(T);
}
long check_SHIFT(long x) {
  long center;
  long y, z;
  long INTERFACE_POS_MAX = 0;
  for (z=0; z < rows_z; z++) {
    for (y=1; y <=(rows_y-1); y++) {
  //     center =  gidy   + (x)*numy[levels];
      center = x*layer_size + z*rows_y + y; 
  //     printf("center=%ld\n",center);
      if ((gridinfo[center-1].phia[NUMPHASES-1]-(1.0-gridinfo[center-1].phia[NUMPHASES-1]) < 0.0) 
        && (gridinfo[center].phia[NUMPHASES-1]-(1.0-gridinfo[center].phia[NUMPHASES-1]) > 0.0) ) {
        if (y > INTERFACE_POS_MAX) {
          INTERFACE_POS_MAX = y;
        }
      }
    }
  }
  if (INTERFACE_POS_MAX > 0) {
    return INTERFACE_POS_MAX;
  } else {
    return 0;
  }
}

void apply_shiftY_cscl(struct csle *cscl, long INTERFACE_POS_GLOBAL) {
  //Shift by one cell in the negative y-direction
  long x, y, z;
  long gidy;
  double chemical_potential;
  double c[NUMCOMPONENTS-1];
  
  //for(x=0; x < rows_x; x++) {
  //  for(z=0; z < rows_z; z++) {
  //    for (y=0; y <= (rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y++) {
  //      gidy = x*layer_size + z*rows_y + y;
  //      for (b=0; b < NUMPHASES; b++) {
  //        for (k=0; k < NUMCOMPONENTS-1; k++) {
  //          cscl[gidy].comie[b][k] = cscl[gidy+(INTERFACE_POS_GLOBAL-shiftj)].comie[b][k];
  //        }
  //      }
  //    }
  //    
  //    for (y=(rows_y-(INTERFACE_POS_GLOBAL-shiftj)); y<=(rows_y-1); y++) {
  //      gidy = x*layer_size + z*rows_y + y;
  //      for (b=0; b < NUMPHASES; b++) {
  //        for (k=0; k < NUMCOMPONENTS-1; k++) {
  //          cscl[gidy].comie[b][k] = ceq[b][b][k];
  //        }
  //      }
  //    }
  //  }
  //}
}

#endif
